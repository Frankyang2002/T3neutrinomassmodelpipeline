"""Batch 2: read-only, provenance-checked analysis of Chat A all-16 output.

No new physical inputs, numerical UV evaluations, RG running, or EFT matching.
All mass claims remain conditional on the frozen Weyl dictionary.
"""
from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import statistics
from typing import Any

import numpy as np

from Numerical.lo_certification.direct_uv.native_amplitude_contract import model_dictionary
from Numerical.lo_certification.direct_uv.native_uv_loop_reference import decode

ROOT_OUTPUT = Path('output/general_observable/lo_pattern_study16')
SOURCE_INPUT = ROOT_OUTPUT / 'provisional_01'
TARGET_OUTPUT = ROOT_OUTPUT / 'analysis_batch02_01'
T_VALUES = (-1., 0., .5, 1., 2.)
HEAVY_VALUES = (.5, 1., 2.)
MODES = ('y1','y2','lambdaT3','heavy_mass')
PAIR_LABELS = ('y1_HC x y1_HC','y1_HC x y2_HC','y2_HC x y2_HC')


class AnalysisError(ValueError):
    """Input/output violation, model mismatch, or failed identity."""


def _require(ok: bool, msg: str) -> None:
    if not ok:
        raise AnalysisError(msg)


def _load(path: Path) -> dict:
    _require(path.is_file(), 'Missing source result: '+str(path))
    return json.loads(path.read_text(encoding='utf-8-sig'))


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _num(value: Any, label: str) -> float:
    try:
        v=float(value)
    except (ValueError,TypeError,OverflowError) as exc:
        raise AnalysisError(label+': non-numeric value') from exc
    _require(math.isfinite(v), label+': non-finite value')
    return v


def _matrix(raw: Any, label: str) -> np.ndarray:
    return decode(raw,label,(3,3))


def _norm(mat: np.ndarray) -> float:
    return float(np.linalg.norm(mat))


def _relative(a: np.ndarray, b: np.ndarray, *, conditioning: float = 1.) -> float:
    # The conditioning parameter allows a relative comparison to remain
    # meaningful when the observable is cancelled well below individual terms.
    return _norm(a-b)/max(_norm(a),_norm(b),1e-12*conditioning,1e-100)


def _rows_csv(path:Path) -> list[dict[str,str]]:
    _require(path.is_file(), 'Missing source CSV '+str(path))
    with path.open('r',encoding='utf-8-sig',newline='') as fp:
        reader=csv.DictReader(fp)
        _require(bool(reader.fieldnames), 'CSV has no header: '+str(path))
        return list(reader)


def _write_json(path:Path, data:dict) -> None:
    path.write_text(json.dumps(data,allow_nan=False,indent=2)+'\n',encoding='utf-8')


def _write_csv(path:Path, fields:list[str], data:list[dict]) -> None:
    with path.open('w',encoding='utf-8',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(data)


def _source_hash_guard(root:Path, source:dict, summary:dict, model:str) -> None:
    """Verify frozen benchmark/component bytes AND original guarded executable files."""
    bpath=root/source['benchmark']
    _require(bpath.is_file(),model+': frozen benchmark missing')
    _require(_hash(bpath)==source['benchmark_sha256'],model+': frozen benchmark SHA256 drift')
    # Existing result records absolute input paths on Windows. Relative paths
    # from a custom provider are resolved against the project root.
    raw=Path(source['input'])
    inp=raw if raw.is_absolute() else root/raw
    _require(inp.is_file(),model+': candidate Weyl input missing '+str(inp))
    _require(_hash(inp)==source['input_sha256'],model+': Weyl input SHA256 drift')
    provenance=summary['provenance']
    for key,meta in provenance['source_audit']['source_files'].items():
        path=root/meta['path']
        _require(path.is_file() and _hash(path)==meta['sha256'],
                 'Executable source guard SHA256 drift: '+key+' '+str(path))
    audit_dir=root/'output/general_observable/lo_direct_one_loop/source_audits'
    for name,dig in provenance['source_audit_sha256'].items():
        p=audit_dir/name
        if dig is not None:
            _require(p.is_file() and _hash(p)==dig,'Source audit SHA256 drift: '+name)


def _validate_scans(raw:list[dict[str,str]], expected:set[str]) -> dict[tuple[str,str,float],dict]:
    scans={}
    for r in raw:
        model=r.get('model'); mode=r.get('mode')
        _require(model in expected and mode in MODES,'Unknown scan model/mode: '+str((model,mode)))
        t=_num(r.get('t'),f'{model}/{mode} scan t')
        allowed=T_VALUES if mode!='heavy_mass' else HEAVY_VALUES
        _require(t in allowed,f'{model}/{mode}: unrecognized t={t}')
        key=(model,mode,t)
        _require(key not in scans,'Duplicate scan '+str(key))
        _require(r.get('status')=='conditional_uv_scanned',
                 'Blocked/incomplete scan '+str(key)+': '+str(r.get('reason')))
        row={'model':model,'mode':mode,'t':t,'status':r['status']}
        for field in ('m1_ev','m2_ev','m3_ev','mass_norm_gev','cancellation_factor'):
            row[field]=_num(r.get(field),str(key)+'/'+field)
        _require(0<=row['m1_ev']<=row['m2_ev']<=row['m3_ev'],
                 str(key)+': unordered or negative Takagi masses')
        _require(row['mass_norm_gev']>=0 and row['cancellation_factor']>=0,
                 str(key)+': invalid mass norm / cancellation factor')
        residual=r.get('matrix_scaling_residual')
        row['matrix_scaling_residual']=(None if residual in (None,'') else _num(residual,str(key)+'/matrix residual'))
        scans[key]=row
    want={(m,mode,t) for m in expected for mode in MODES
          for t in (HEAVY_VALUES if mode=='heavy_mass' else T_VALUES)}
    _require(set(scans)==want,
             f'Incomplete scan grid: expected {len(want)}, got {len(scans)}; missing={sorted(want-set(scans))[:6]}')
    return scans


def _validate_model(model:str, detail:dict, row:dict, source_summary:dict, root:Path,
                    verify_hashes:bool) -> dict:
    _require(detail.get('model')==model and detail.get('status')=='conditional_uv_evaluated',
             model+': per-model result label/status changed')
    _require(detail.get('physical_normalisation_certified') is False,
             model+': false certification promotion in per-model JSON')
    source=detail['source']
    _require(row.get('benchmark_sha256')==source['benchmark_sha256'] and
             row.get('component_input_sha256')==source['input_sha256'],
             model+': summary/detail SHA256 conflict')
    if verify_hashes:
        _source_hash_guard(root,source,source_summary,model)
    m=_matrix(detail['mnu_gev'],model+' mnu GeV')
    m_ev=_matrix(detail['mnu_ev'],model+' mnu eV')
    _require(_relative(m_ev,1e9*m)<1e-12,model+': GeV/eV matrix conversion mismatch')
    _require(_norm(m-m.T)<=1e-10*max(_norm(m),1e-100),model+': mnu not symmetric')
    ms=[_num(v,model+' Takagi mass') for v in detail['takagi_masses_ev_sorted']]
    _require(len(ms)==3 and 0<=ms[0]<=ms[1]<=ms[2],model+': invalid Takagi spectrum')
    eigen=np.sort(np.linalg.svd(m,compute_uv=False)*1e9)
    _require(max(abs(eigen-np.array(ms)))<=2e-6*max(ms[2],1e-100),
             model+': Takagi/SVD disagreement beyond numerical zero-mode tolerance')
    _require(max(abs(float(row[f'm{i+1}_ev'])-ms[i]) for i in range(3))<=
             1e-10*max(ms[2],1e-100),model+': summary/model mass mismatch')
    pairs=detail['source_vertex_pair_matrices_gev']
    _require(set(pairs)==set(PAIR_LABELS),model+': source pair labels changed')
    a,b,c=(_matrix(pairs[x],model+' '+x) for x in PAIR_LABELS)
    natural=sum(_norm(x) for x in (a,b,c))
    _require(_relative(m,a+b+c,conditioning=natural)<2e-5,
             model+': vertex-pair sum cannot reconstruct mnu')
    partition=detail['source_charge_and_group']['charge_partition_direct_exact']['charge_channel_matrices']
    qmat={lab:_matrix(rec['matrix_gev'],model+' '+lab) for lab,rec in partition.items()}
    _require(len(qmat)>0,model+': empty exact charge decomposition')
    qsum=sum(qmat.values(),np.zeros((3,3),complex))
    _require(_relative(m,qsum,conditioning=sum(_norm(q) for q in qmat.values()))<2e-5,
             model+': charge partition cannot reconstruct mnu')
    chs=detail['source_charge_and_group']['native_group_audit_channels']
    predicted={float(q['charge_numeric']) for q in chs if q['nonzero']}
    _require(len(predicted)==int(row['candidate_topology_charge_channel_count']),
             model+': candidate charge-channel count drift')
    for q in predicted:
        _require(f'{q:g}|{q:g}' in qmat,model+': CG channel absent from exact charge partition')
    # Charge partition can contain Yukawa-only null topology channels and
    # cross-charge interference. Do not equate all of these with CG channels.
    active={f'{q:g}|{q:g}':qmat[f'{q:g}|{q:g}'] for q in sorted(predicted)}
    qnorm=sum(_norm(x) for x in active.values())
    # Include the remainder explicitly: it may contain same-vertex pieces or
    # cross-charge interference in arbitrary transformed input dictionaries.
    remainder=m-sum(active.values(),np.zeros((3,3),complex))
    qnorm+=_norm(remainder)
    coherent=_norm(m)/max(qnorm,1e-100)
    _require(coherent<=1+2e-5,model+': incoherent/coherent charge norm paradox')
    overlaps=[]
    charges=sorted(predicted)
    for i,q in enumerate(charges):
        for q2 in charges[i+1:]:
            x=active[f'{q:g}|{q:g}'];y=active[f'{q2:g}|{q2:g}']
            cos=float(np.vdot(x,y).real/max(_norm(x)*_norm(y),1e-100))
            _require(-1-1e-8<=cos<=1+1e-8,model+': invalid channel overlap')
            overlaps.append({'charge_1':q,'charge_2':q2,'cosine_frobenius':max(-1.,min(1.,cos)),
                             'complex_overlap_real':float(np.vdot(x,y).real),
                             'complex_overlap_imag':float(np.vdot(x,y).imag)})
    return {'matrix_gev':m,'pair_matrices':(a,b,c),'masses_ev':ms,'charge_matrices':qmat,
            'cg_charges':charges,'candidate_charge_matrices':active,
            'unattributed_charge_remainder_norm_gev':_norm(remainder),
            'charge_coherence_ratio':coherent,'channel_overlaps':overlaps,
            'cg_weights':chs,'input_source':source,
            'cancellation_factor':_num(row['cancellation_factor'],model+'/cancellation factor'),
            'source_reference_relative_error':_num(row['source_reference_relative_error'],model+'/reference error')}


def _scan_stats(model:str, rec:dict, scans:dict) -> dict:
    m=rec['matrix_gev'];mnorm=_norm(m)
    ms=rec['masses_ev'];a,b,c=rec['pair_matrices']
    stats={'yukawa_polynomial_max_mass_residual_relative_to_m3':{},
           'yukawa_linear_max_matrix_residual':{},
           'lambda_sign_flip_matrix_residual':None,
           'lambda_zero_m3_ev':None,
           'lambda_t2_linear_matrix_residual':None,
           'lambda_t05_linear_matrix_residual':None,
           'lambda_nonlinearity_norm_ratio_t2':None,
           'heavy_mass_ratio_of_m3_t2_to_t05':None,
           'heavy_mass_log_slope_m3_t05_to_t2':None}
    for mode in MODES:
        original=scans[(model,mode,1.)]
        _require(max(abs(original[f'm{i+1}_ev']-ms[i]) for i in range(3))<=
                 1e-8*max(ms[2],1e-100),model+'/'+mode+': t=1 baseline mismatch')
    for mode in ('y1','y2'):
        worstpoly=0.;worstlin=0.
        for t in T_VALUES:
            row=scans[(model,mode,t)]
            pred=(t*t*a+t*b+c) if mode=='y1' else (a+t*b+t*t*c)
            pm=np.sort(np.linalg.svd(pred,compute_uv=False)*1e9)
            deviation=float(np.max(np.abs(pm-np.array([row[f'm{i+1}_ev'] for i in range(3)]))))
            worstpoly=max(worstpoly,deviation/max(ms[2],1e-100))
            if row['matrix_scaling_residual'] is not None:
                worstlin=max(worstlin,row['matrix_scaling_residual'])
        # This quadratic Yukawa identity is exact for any fixed spectra,
        # unlike the simpler linearity claim which requires null self-pairs.
        _require(worstpoly<4e-5,model+'/'+mode+': exact Yukawa quadratic identity failed')
        stats['yukawa_polynomial_max_mass_residual_relative_to_m3'][mode]=worstpoly
        stats['yukawa_linear_max_matrix_residual'][mode]=worstlin
    for t,label in ((-1.,'lambda_sign_flip_matrix_residual'),
                    (2.,'lambda_t2_linear_matrix_residual'),
                    (.5,'lambda_t05_linear_matrix_residual')):
        stats[label]=scans[(model,'lambdaT3',t)]['matrix_scaling_residual']
    stats['lambda_zero_m3_ev']=scans[(model,'lambdaT3',0.)]['m3_ev']
    stats['lambda_nonlinearity_norm_ratio_t2']=(scans[(model,'lambdaT3',2.)]['mass_norm_gev'] /
                                                      max(2*mnorm,1e-100))
    hlow=scans[(model,'heavy_mass',.5)]['m3_ev'];hhi=scans[(model,'heavy_mass',2.)]['m3_ev']
    stats['heavy_mass_ratio_of_m3_t2_to_t05']=hhi/max(hlow,1e-100)
    stats['heavy_mass_log_slope_m3_t05_to_t2']=(math.log(hhi/hlow)/math.log(4.)
        if hlow>0 and hhi>0 else None)
    # Exact lambda sign flip is expected in the restricted B=0 two-scalar
    # block structure; report numerical residual rather than hiding it.
    return stats


def _plots(out:Path, rows:list[dict], scan_rows:list[dict]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return ['matplotlib unavailable; no figures produced']
    models=[r['model'] for r in rows]
    mass=np.array([[r[f'm{i}_ev'] for i in (1,2,3)] for r in rows])
    fig,ax=plt.subplots(figsize=(11,4.5))
    xx=np.arange(len(models))
    for i in range(3):
        ax.plot(xx,np.maximum(mass[:,i],1e-18),marker='o',markersize=3,label=f'sorted m{i+1}')
    ax.set_yscale('log');ax.set_xticks(xx,models,rotation=60,ha='right')
    ax.set_ylabel('Takagi singular values [eV], provisional')
    ax.set_title('Frozen restricted T3 benchmarks: full 3-mass spectrum')
    ax.legend();fig.tight_layout();fig.savefig(out/'batch02_all_three_masses.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for family in 'ABCDE':
        sub=[r for r in rows if r['family']==family]
        ax.plot([r['alpha'] for r in sub],[r['m3_ev'] for r in sub],marker='o',label=family)
    ax.set_xlabel(r'Hypercharge assignment $\alpha$')
    ax.set_ylabel(r'Largest sorted Takagi mass $m_3$ [eV], provisional')
    ax.set_title('Hypercharge scans across distinct frozen benchmarks')
    ax.legend();fig.tight_layout();fig.savefig(out/'batch02_m3_by_alpha.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,4))
    ax.bar(xx,[r['charge_coherence_ratio'] for r in rows])
    ax.set_ylim(0,1.05);ax.set_xticks(xx,models,rotation=60,ha='right')
    ax.set_ylabel('||sum charge matrices|| / sum ||charge matrices||')
    ax.set_title('Charge-channel coherence, including residual partition')
    fig.tight_layout();fig.savefig(out/'batch02_charge_coherence.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for family in 'ABCDE':
        vals=[r for r in rows if r['family']==family]
        ax.scatter([r['lambda_t2_linear_matrix_residual'] for r in vals],
                   [r['m3_ev'] for r in vals],label=family)
    ax.set_xlabel(r'$||m_\nu(2\lambda)-2m_\nu(\lambda)||/||m_\nu(\lambda)||$')
    ax.set_ylabel('Frozen largest mass [eV]')
    ax.set_title('Finite-mixing nonlinearity versus frozen mass (conditional)')
    ax.legend();fig.tight_layout();fig.savefig(out/'batch02_lambda_nonlinearity.png',dpi=150);plt.close(fig)
    return []


def _markdown(report:dict) -> str:
    rows=report['models'];fams=report['families']
    lines=['# Chat A — Batch 2 numerical pattern analysis', '',
           'Scope: the **existing** all-16 restricted direct-UV output at fixed '
           '`mu0=1e7 GeV`, not a new numerical evaluation. No EFT or RGE.',
           '', '**Certification: 0/16 physically certified.** All masses and '
           'cross-family rankings are conditional on the unverified source-to-Weyl dictionary.', '',
           '## Frozen benchmark comparison', '',
           '| Model | alpha | m1 (eV) | m2 (eV) | m3 (eV) | channels | charge coherence | lambda t=2 residual |',
           '|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['model']} | {r['alpha']:+d} | {r['m1_ev']:.5g} | "
                     f"{r['m2_ev']:.5g} | {r['m3_ev']:.6g} | "
                     f"{r['candidate_channels']} | {r['charge_coherence_ratio']:.5g} | "
                     f"{r['lambda_t2_linear_matrix_residual']:.3g} |")
    lines+=['','`m1,m2,m3` are sorted singular values, **not** experimental normal/inverted ordering labels.',
            '', '## Family statistics of frozen m3', '',
            '| Family | count | mean m3 (eV) | min | max |',
            '|---|---:|---:|---:|---:|']
    for family in 'ABCDE':
        f=fams[family]
        lines.append(f"| {family} | {f['count']} | {f['mean_m3_ev']:.6g} | "
                     f"{f['min_m3_ev']:.6g} | {f['max_m3_ev']:.6g} |")
    lines += ['', '## Mechanism and robust numerical interpretation','',
              '1. **Yukawa scan:** with fixed spectra, `m(t*y1,y2)=t^2 M11+t M12+M22` and '
              '`m(y1,t*y2)=M11+t M12+t^2 M22`. The report checks all three sorted '
              'singular values against these exact quadratic identities. It separately '
              'reports whether the special topology-only linear law is numerically accurate.',
              '2. **Scalar mixing:** in the restricted input the off-diagonal S1–S2 block '
              'is scaled. With B=0 and only topology cross-terms the full mass is odd '
              'under a mixing-sign flip, but linear scaling at finite mixing is approximate. '
              'The t=2 residual quantifies departures without interpreting them as source bugs.',
              '3. **Heavy mass:** `M -> t M` changes both mass insertions and loop functions. '
              'The fitted t=0.5 to t=2 log slope is empirical and is not a universal -1 law.',
              '4. **Charge channels:** `||sum M_Q||/sum ||M_Q||` is a coherent/interfering '
              'Frobenius diagnostic, NOT a probabilistic branching fraction. Any remainder '
              'outside the source-CG topology charges is included in the denominator. '
              'Negative pairwise real overlap identifies destructive matrix interference.',
              '5. **Individually sorted masses do not add:** the mass matrix, not its Takagi '
              'spectrum, is additive across internal charge sectors.',
              '', '## E-family two-channel overlap', '']
    for r in rows:
        if r['family'] != 'E':continue
        pair=r['channel_overlaps']
        if pair:
            p=pair[0]
            lines.append(f"- `{r['model']}`: Q={p['charge_1']:g},{p['charge_2']:g}; "
                         f"normalised real Frobenius overlap = {p['cosine_frobenius']:+.6g}; "
                         f"charge coherence = {r['charge_coherence_ratio']:.6g}.")
        else:
            lines.append(f"- `{r['model']}`: no two-channel overlap available.")
    lines+=['', '## Source and audit limits', '',
            '- No cross-package Matchete/RGBeta spinor, portal, or CG normalisation is promoted.',
            '- All old files are read-only. Original file SHA256 values are recorded in '
            '`SOURCE_MANIFEST.json`, and the existing source guard and frozen inputs are checked.',
            '- Cross-model differences can reflect different Yukawa matrices, heavy masses, '
            'portal masses, scalar mixing, interference, and any still-unverified source mapping.',
            '- No EFT, threshold matching, running, new fitting, or altered physics input is performed.', '']
    return '\n'.join(lines)


def analyse(root:Path, *, study_dir:Path|None=None, output_dir:Path|None=None,
            figures:bool=True) -> dict:
    """Produce a separate Batch-2 output directory only after complete preflight.

    Input is the previous runner's persisted output; no direct_one_loop.evaluate call.
    """
    root=root.resolve()
    owned=(root/ROOT_OUTPUT).resolve()
    source=(study_dir or root/SOURCE_INPUT).resolve()
    out=(output_dir or root/TARGET_OUTPUT).resolve()
    _require(source.is_relative_to(owned) and source!=owned,'Input must be within Chat A output tree')
    _require(out.is_relative_to(owned) and out!=owned and out!=source,
             'Analysis output must be a separate NEW directory within the Chat A tree')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Refusing to replace existing analysis: '+str(out))
    report=_load(source/'all16_summary.json')
    _require(report.get('format')=='t3-conditional-pattern-study-all16-v1',
             'Unexpected source study schema')
    _require(report.get('physical_normalisation_certified') is False and
             report['counts']['physical_source_certified']==0,
             'False/promoted physical certification detected in source study')
    _require(report['provenance']['source_audit']['status']=='source_forms_verified' and
             report['provenance']['source_audits_verified_for_this_input'] is True,
             'Source-run executable/audit guard was not verified')
    _require(report['counts']['catalogue']==16 and
             report['counts']['conditional_evaluated']==16 and report['counts']['blocked']==0,
             'Source study is not a complete provisional all-16 run')
    _require(report['counts']['scan_passes']==288 and report['counts']['scan_blocked']==0,
             'Source study scans are incomplete or blocked')
    names=[x['model'] for x in model_dictionary()]
    expected=set(names)
    _require(len(expected)==16,'Source catalogue unexpectedly changed')
    rows=report['models']
    _require(len(rows)==16 and [r['model'] for r in rows]==names,
             'Source study model inventory or order differs from live catalogue')
    summary_csv=_rows_csv(source/'all16_summary.csv')
    _require([r['model'] for r in summary_csv]==names,'CSV/JSON catalogue mismatch')
    scans=_validate_scans(_rows_csv(source/'scans.csv'),expected)
    weights=_load(source/'channel_weights.json')
    _require(weights.get('status')=='native_group_CG_candidate_only_not_physical_normalization' and
             set(weights['models'])==expected,'Source channel weight report missing/mutated')
    common=_load(source/'common_kernel_group_only.json')
    _require(common.get('physical_normalisation_certified') is False and len(common['rows'])==16,
             'Formal group comparison input missing / incorrectly certified')
    diagnostics=[];charge_rows=[];out_rows=[]
    manifest={'original_directory':str(source),'original_file_sha256':{},
              'current_checkout_source_hashes_checked':True,'physical_source_certified':False}
    required_paths=[source/x for x in ('all16_summary.json','all16_summary.csv','scans.csv',
                                     'channel_weights.json','common_kernel_group_only.json')]
    required_paths += [source/'models'/f'{m}.json' for m in names]
    for path in required_paths:
        _require(path.is_file(),'Missing prior result '+str(path))
        manifest['original_file_sha256'][str(path.relative_to(source))]=_hash(path)
    for i,model in enumerate(names):
        row=rows[i];csvrow=summary_csv[i]
        _require(row['status']=='conditional_uv_evaluated' and
                 row.get('physical_normalisation_certified') is False and
                 csvrow['status']==row['status'] and
                 str(csvrow['physical_normalisation_certified']).lower()=='false',
                 model+': source summary status/certification conflict')
        detail=_load(source/'models'/f'{model}.json')
        rec=_validate_model(model,detail,row,report,root,verify_hashes=True)
        _require(weights['models'][model].get('relative_physical_normalization_certified') is False,
                 model+': source CG weight falsely certified')
        _require(weights['models'][model]['channels']==rec['cg_weights'],
                 model+': channel_weights / per-model CG audit conflict')
        stats=_scan_stats(model,rec,scans)
        mf=rec['masses_ev']
        result={'model':model,'family':model[0],'alpha':int(model[2:]),
                'm1_ev':mf[0],'m2_ev':mf[1],'m3_ev':mf[2],
                'sum_masses_ev':sum(mf),'mass_frobenius_ev':_norm(rec['matrix_gev'])*1e9,
                'candidate_channels':len(rec['cg_charges']),
                'candidate_charges':rec['cg_charges'],
                'charge_coherence_ratio':rec['charge_coherence_ratio'],
                'unattributed_charge_remainder_norm_gev':rec['unattributed_charge_remainder_norm_gev'],
                'channel_overlaps':rec['channel_overlaps'],
                'source_reference_relative_error':rec['source_reference_relative_error'],
                'cancellation_factor':rec['cancellation_factor'],
                'benchmark_sha256':row['benchmark_sha256'],
                'component_input_sha256':row['component_input_sha256'],
                'physical_normalisation_certified':False,**stats}
        out_rows.append(result)
        for lab, mat in rec['charge_matrices'].items():
            charge_rows.append({'model':model,'charge_pair':lab,
                                'frob_norm_gev':_norm(mat),
                                'matrix_fraction_of_total_frob':_norm(mat)/max(_norm(rec['matrix_gev']),1e-100),
                                'in_source_topology_candidate':lab in rec['candidate_charge_matrices'],
                                'physical_normalisation_certified':False})
        for mode in ('y1','y2'):
            diagnostics.append({'model':model,'test':mode+'_exact_quadratic_mass_identity',
                                'relative_residual':stats['yukawa_polynomial_max_mass_residual_relative_to_m3'][mode],
                                'status':'passed_conditional'})
        diagnostics.append({'model':model,'test':'lambda_sign_oddness',
                            'relative_residual':stats['lambda_sign_flip_matrix_residual'],
                            'status':'measured_conditional'})
    family_stats={}
    for family in 'ABCDE':
        vals=[r['m3_ev'] for r in out_rows if r['family']==family]
        family_stats[family]={'count':len(vals),'mean_m3_ev':statistics.mean(vals),
                              'min_m3_ev':min(vals),'max_m3_ev':max(vals)}
    all_m3=[r['m3_ev'] for r in out_rows]
    out_report={'format':'t3-chatA-readonly-batch02-numerical-patterns-v1',
                'input_run_sha256':manifest['original_file_sha256']['all16_summary.json'],
                'counts':{'models':16,'source_scan_rows':len(scans),'physical_source_certified':0,
                          'yukawa_polynomial_checks_passed':32},
                'global_m3_statistics':{'min_ev':min(all_m3),'max_ev':max(all_m3),
                                       'ratio_max_over_min':max(all_m3)/min(all_m3),
                                       'mean_ev':statistics.mean(all_m3),'median_ev':statistics.median(all_m3)},
                'models':out_rows,'families':family_stats,
                'verification':'Persisted run, frozen source bytes and original executable guards checked; no UV recomputation',
                'scope':'restricted fixed-scale direct UV conditional, no EFT',
                'physical_normalisation_certified':False}
    # Nothing is written before all model and source checks have passed.
    out.mkdir(parents=True,exist_ok=True)
    _write_json(out/'batch02_analysis.json',out_report)
    _write_json(out/'SOURCE_MANIFEST.json',manifest)
    _write_csv(out/'batch02_model_summary.csv',
               ['model','family','alpha','m1_ev','m2_ev','m3_ev','sum_masses_ev',
                'candidate_channels','charge_coherence_ratio','source_reference_relative_error',
                'cancellation_factor','lambda_t2_linear_matrix_residual',
                'lambda_sign_flip_matrix_residual','heavy_mass_log_slope_m3_t05_to_t2',
                'unattributed_charge_remainder_norm_gev','physical_normalisation_certified'],out_rows)
    _write_csv(out/'batch02_charge_matrix_norms.csv',
               ['model','charge_pair','frob_norm_gev','matrix_fraction_of_total_frob',
                'in_source_topology_candidate','physical_normalisation_certified'],charge_rows)
    _write_json(out/'batch02_numerical_checks.json',{'checks':diagnostics,
                 'caution':'Numerical identities do not certify physical spinor/CG mapping'})
    (out/'BATCH02_REPORT.md').write_text(_markdown(out_report),encoding='utf-8')
    figure_warnings=_plots(out,out_rows,list(scans.values())) if figures else []
    handoff={'format':'t3-chatA-batch02-handoff-v1',
             'input_study_sha256':manifest['original_file_sha256']['all16_summary.json'],
             'source_local_revision':report['provenance']['revision'],
             'model_count':16,'scanned_points':len(scans),
             'numerical_source_checks_passed':True,
             'physical_source_certified':False,
             'safe_to_integrate_for_physical_mass':False,
             'safe_to_integrate_for_conditional_numerical_patterns':True,
             'unresolved':['B1 native four-spinor to Weyl physical dictionary',
                           'B2 native SU2/portal/Matchete-RGBeta coupling dictionary'],
             'figure_warnings':figure_warnings}
    _write_json(out/'HANDOFF.json',handoff)
    (out/'HANDOFF.md').write_text(
        '# Chat A — Batch 2 handoff\n\n'
        f'- Input fixed-source run: `{source}`, SHA256 `{handoff["input_study_sha256"]}`.\n'
        f'- Model evaluations: 16/16; full-grid scan points: {len(scans)}.\n'
        '- Full spectrum, exact Yukawa quadratic identities, charge interference, '
        'and finite-mixing/heavy-mass diagnostics were derived from persisted results.\n'
        '- All underlying source guard, audit, benchmark and component-input hashes matched.\n'
        '- No direct UV source masses were recomputed, no input changed, no previous result overwritten.\n'
        '- Physical cross-package normalization certification remains **0/16**; '
        'do not integrate as a physical-mass certificate.\n',encoding='utf-8')
    return out_report
