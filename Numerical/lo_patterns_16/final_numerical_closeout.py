"""Final read-only closeout of the 16 source-conditional T3 neutrino-mass results.

Consumes existing frozen study, Batch 5 and Batch 6 files. Recomputes no
physical amplitudes; produces only a concise numerical summary and CSVs.
Crucially, it cannot certify the Matchete/RGBeta/Weyl action dictionary.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

BASE = Path('output/general_observable/lo_pattern_study16')
STUDY = BASE / 'provisional_01'
B5 = BASE / 'analysis_batch05_01'
B6 = BASE / 'analysis_batch06_01'
DEST = BASE / 'final_numerical_closeout_01'
MODEL_NAMES = tuple(f'{f}_{a:+d}' for f, alphas in (
    ('A', (-4,-2,0)), ('B', (-3,-1,1)), ('C', (-3,-1,1)),
    ('D', (-2,0,2)), ('E', (-4,-2,0,2))) for a in alphas)
TEXTURES = ('A_-2','C_-1','C_+1','E_-2')
SCENARIOS = ('baseline','scalar_common_shift','heavy_generation_tilt',
             'large_scalar_mixing','charge_split_diagnostic')
SCAN_GRID = {'y1':(-1.0,0.0,0.5,1.0,2.0),
             'y2':(-1.0,0.0,0.5,1.0,2.0),
             'lambdaT3':(-1.0,0.0,0.5,1.0,2.0),
             'heavy_mass':(0.5,1.0,2.0)}


class CloseoutBlocked(ValueError):
    """Refuse a closeout rather than publishing inconsistent provisional data."""


def require(ok: Any, explanation: str) -> None:
    if not ok:
        raise CloseoutBlocked(explanation)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'Missing input: '+str(path))
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), 'Expected JSON object: '+str(path))
    return obj


def read_csv(path: Path) -> list[dict[str,str]]:
    require(path.is_file(), 'Missing input: '+str(path))
    with path.open('r', encoding='utf-8-sig', newline='') as fh:
        rows = list(csv.DictReader(fh))
    require(bool(rows), 'Empty CSV: '+str(path))
    return rows


def checked_float(value: Any, name: str, *, nonnegative: bool = True) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as ex:
        raise CloseoutBlocked(f'{name}: invalid number {value!r}') from ex
    require(math.isfinite(result) and (not nonnegative or result >= 0),
            f'{name}: nonfinite/negative numerical result')
    return result


def close(a: Any, b: Any, *, rel: float = 2e-7, abs_tol: float = 1e-10) -> bool:
    return math.isclose(checked_float(a,'numeric comparison'),
                        checked_float(b,'numeric comparison'),rel_tol=rel,abs_tol=abs_tol)


def index_unique(rows: list[dict], key: str, label: str) -> dict:
    require(all(key in r for r in rows),label+': missing '+key)
    indexed = {r[key]:r for r in rows}
    require(len(indexed) == len(rows),label+': duplicate '+key)
    return indexed


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_digest(path: Path, expected: str, label: str) -> str:
    """Accept exact saved bytes, or an otherwise identical CRLF checkout.

    This is a *source-file content check*, not a current-HEAD equivalence proof.
    """
    require(isinstance(expected,str) and len(expected)==64 and
            all(c in '0123456789abcdef' for c in expected),
            label+': invalid recorded SHA256')
    require(path.is_file(),label+': missing '+str(path))
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()==expected:
        return 'raw-bytes'
    normalized=raw.replace(b'\r\n',b'\n')
    require(normalized != raw and hashlib.sha256(normalized).hexdigest()==expected,
            label+': SHA256 differs beyond checkout line endings; '+str(path))
    return 'CRLF-to-LF-only'


def _csv_table(rows: list[dict], fields: tuple[str,...]) -> str:
    buf=io.StringIO(newline='')
    writer=csv.DictWriter(buf,fieldnames=fields,lineterminator='\n',extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    return buf.getvalue()


def _fmt(x: Any, precision: int = 6) -> str:
    return format(float(x),f'.{precision}g')


def compute(root: Path) -> tuple[dict,dict[str,str]]:
    """Validate all inputs and return only in-memory results, with no writes."""
    root=root.resolve()
    p_study=root/STUDY/'all16_summary.json'
    p_b5=root/B5/'batch05_synthesis.json'
    p_b5manifest=root/B5/'SOURCE_MANIFEST.json'
    p_b6=root/B6/'batch06_factorization.json'
    p_b6manifest=root/B6/'SOURCE_MANIFEST.json'
    study=read_json(p_study)
    h=read_json(root/STUDY/'HANDOFF.json')
    b5=read_json(p_b5)
    m5=read_json(p_b5manifest)
    b6=read_json(p_b6)
    m6=read_json(p_b6manifest)

    require(study.get('format')=='t3-conditional-pattern-study-all16-v1',
            'Study schema mismatch')
    require(b5.get('format')=='t3-chatA-batch05-analytical-yukawa-texture-v1',
            'Batch 5 schema mismatch')
    require(b6.get('format')=='t3-chatA-batch06-heavy-kernel-factorization-v1',
            'Batch 6 schema mismatch')
    require(m5.get('format')=='t3-chatA-batch05-source-manifest-v1' and
            m6.get('format')=='t3-chatA-batch06-manifest-v1',
            'Source manifest schema mismatch')
    for label, obj in [('study',study),('B5',b5),('B5 manifest',m5),
                       ('B6',b6),('B6 manifest',m6)]:
        require(obj.get('physical_normalisation_certified') is False,
                label+': physical-normalisation flag has changed')
    for label,obj in [('B5',b5),('B6',b6)]:
        require(obj.get('source_conditional_only') is True and
                obj.get('safe_for_physical_source_certification') is False,
                label+': conditional scope guard failed')
    require(b5.get('model_count')==b6.get('model_count')==16,
            'Wrong B5/B6 model count')
    require(b5.get('reference_texture_models')==list(TEXTURES) and
            close(b5.get('common_scalar_mixing_abs'),0.1),
            'B5 texture-control interpretation has changed')
    require(b6.get('reference_model')=='A_-2' and
            b6.get('scenarios')==list(SCENARIOS),
            'B6 reference/scenario list has changed')
    counts=study['counts']
    require(counts == {'catalogue':16,'conditional_evaluated':16,
                       'blocked':0,'physical_source_certified':0,
                       'scan_passes':288,'scan_blocked':0},
            'Study evidence counts do not equal 16/16, 288/288, 0/16')
    require(h.get('physical_source_certified') is False and
            h.get('checks')==counts and
            h.get('safe_to_integrate_for_physical_mass') is False,
            'Original frozen-study handoff inconsistent')

    # Source-output chain: this is not a claim the current physics source is certified.
    require(m5['study_sha256']==b5['source_manifest_sha256']['study'],
            'Batch 5 study lineage mismatch')
    require(m6['source_sha256']==b6['source_sha256'],
            'Batch 6 source lineage mismatch')
    require(m6['source_sha256']['study']==m5['study_sha256'],
            'Batch 5/6 frozen-study SHA mismatch')
    digests={
        'provisional_01/all16_summary.json':verify_digest(p_study,m5['study_sha256'],'study'),
        'analysis_batch05_01/batch05_synthesis.json':verify_digest(
            p_b5,b6['source_sha256']['batch05'],'B5'),
        'analysis_batch05_01/SOURCE_MANIFEST.json':verify_digest(
            p_b5manifest,b6['source_sha256']['batch05_manifest'],'B5 manifest'),
    }

    sr=study['models']
    five=b5['models']
    require([x['model'] for x in sr]==list(MODEL_NAMES) and
            [x['model'] for x in five]==list(MODEL_NAMES) and
            b6.get('models')==list(MODEL_NAMES),
            'Frozen 16-model catalogue/order mismatch')
    csvrows=read_csv(root/STUDY/'all16_summary.csv')
    cr=index_unique(csvrows,'model','study CSV')
    require(set(cr)==set(MODEL_NAMES) and len(csvrows)==16,
            'Study CSV model catalogue mismatch')
    b5r=index_unique(five,'model','B5 models')
    require(set(m5['models'])==set(MODEL_NAMES) and
            set(m6['models'])==set(MODEL_NAMES),
            'Source manifests do not cover all 16 models')
    results=[]
    for m in sr:
        name=m['model']; old=cr[name]; b=b5r[name]; r5=m5['models'][name];r6=m6['models'][name]
        require(m['status']=='conditional_uv_evaluated' and
                m.get('physical_normalisation_certified') is False and
                old['status']=='conditional_uv_evaluated' and
                old['physical_normalisation_certified']=='False' and
                b['physical_normalisation_certified'] is False and b['conditional_only'] is True,
                name+': unexpected conditional status')
        require(m['family']==old['family']==b['family']==name[0] and
                int(m['alpha'])==int(old['alpha'])==int(b['alpha']) and
                m['fermion_type']==old['fermion_type']==b['fermion_type'],
                name+': inconsistent model identity')
        require((m['benchmark_sha256']==old['benchmark_sha256']==
                 r5['benchmark_sha256']==r6['benchmark_sha256']) and
                (m['component_input_sha256']==old['component_input_sha256']==
                 r5['candidate_weyl_input_sha256']==r6['input_sha256']) and
                r5['study_model_sha256']==r6['study_model_sha256'],
                name+': benchmark/input lineage mismatch')
        for k in ('m1_ev','m2_ev','m3_ev'):
            require(close(m[k],old[k],rel=2e-12,abs_tol=1e-15) and
                    close(m[k],b['observed'][k],rel=2e-7,abs_tol=1e-10),
                    name+': original/B5 mass mismatch '+k)
        m1,m2,m3=(checked_float(m[k],name+' '+k) for k in ('m1_ev','m2_ev','m3_ev'))
        require(0<m1<=m2<=m3,'Wrong positive Takagi mass order: '+name)
        require(close(m1/m3,b['observed']['m1_over_m3']),
                name+': B5 hierarchy ratio mismatch')
        a=checked_float(b['observed']['m1_geometry_A'],name+' geometry')
        coh=checked_float(b['observed']['m1_coherence_C'],name+' coherence')
        require(coh<=1+1e-10 and a<=1+1e-10,
                name+': B5 geometry/coherence outside [0,1]')
        err=checked_float(m['source_reference_relative_error'],name+' reference error')
        require(err<1e-7,name+': source/reference numerical agreement failed')
        results.append({
            'model':name,'family':name[0],'alpha':m['alpha'],
            'fermion_type':m['fermion_type'],
            'charge_channels':m['candidate_topology_charge_channel_count'],
            'm1_ev':m1,'m2_ev':m2,'m3_ev':m3,
            'm1_over_m3':m1/m3,'m2_over_m3':m2/m3,
            'm1_geometry_A':a,'m1_coherence_C':coh,
            'm1_mechanism':b['observed']['m1_mechanism'],
            'source_reference_relative_error':err,
            'physical_normalisation_certified':'False',
        })

    # Independent scans recorded in the frozen 288-point baseline.
    scans=read_csv(root/STUDY/'scans.csv')
    expected={(m,mode,t) for m in MODEL_NAMES for mode,ts in SCAN_GRID.items() for t in ts}
    found={}
    scan_by_model=Counter()
    zero_checks=0
    baseline_checks=0
    max_linear_mass_residual=0.0
    for r in scans:
        key=(r['model'],r['mode'],checked_float(r['t'],'scan t',nonnegative=False))
        require(key in expected and key not in found,'Scan missing/duplicate/unexpected: '+str(key))
        found[key]=r;scan_by_model[key[0]]+=1
        require(r['status']=='conditional_uv_scanned','Blocked scan '+str(key))
        ms=[checked_float(r[k],str(key)+' '+k) for k in ('m1_ev','m2_ev','m3_ev')]
        baseline=next(b for b in results if b['model']==key[0])
        if key[2]==1.0:
            require(all(close(x,baseline[k],rel=2e-7,abs_tol=1e-10)
                        for x,k in zip(ms,('m1_ev','m2_ev','m3_ev'))),
                    'Scan t=1 differs from frozen baseline: '+str(key))
            baseline_checks+=1
        if key[2]==0.0:
            require(max(ms)<1e-10,'Zero Yukawa/quartic scan not zero: '+str(key))
            zero_checks+=1
        if key[1] in ('y1','y2') and key[2] in (0.5,2.0):
            errs=[abs(x-key[2]*baseline[k])/max(key[2]*baseline[k],1e-30)
                  for x,k in zip(ms,('m1_ev','m2_ev','m3_ev'))]
            max_linear_mass_residual=max(max_linear_mass_residual,max(errs))
            require(max(errs)<1e-7,'Yukawa scan not linear: '+str(key))
    require(len(found)==len(scans)==len(expected)==288 and
            set(found)==expected and all(scan_by_model[n]==18 for n in MODEL_NAMES),
            'Scan grid not complete (16 x 18)')
    require(baseline_checks==64 and zero_checks==48,
            'Wrong scan baseline or zero-control coverage')

    controlled=[]
    for name in MODEL_NAMES:
        obs=b5r[name]['common_mixing_controlled']
        require(set(obs)==set(TEXTURES),'B5 texture control missing '+name)
        for tex in TEXTURES:
            row=obs[tex]
            masses=[checked_float(row[k],name+' controlled '+k)
                    for k in ('m1_ev','m2_ev','m3_ev')]
            require(0<masses[0]<=masses[1]<=masses[2],
                    name+' controlled masses not ordered')
            require(close(row['m1_over_m3'],masses[0]/masses[2],rel=1e-9,abs_tol=1e-12),
                    name+' controlled ratio incorrect')
            controlled.append({'model':name,'reference_texture':tex,
                               'lambdaT3_abs':0.1,'m1_ev':masses[0],
                               'm2_ev':masses[1],'m3_ev':masses[2],
                               'm1_over_m3':masses[0]/masses[2],
                               'm1_geometry_A':row['m1_geometry_A'],
                               'm1_coherence_C':row['m1_coherence_C'],
                               'counterfactual_not_fitted':'True'})
    require(len(controlled)==64,'Wrong controlled-texture coverage')

    coeff=b6['coefficient_comparisons']
    coeff_csv=read_csv(root/B6/'batch06_coefficient_tests.csv')
    require(len(coeff)==len(coeff_csv)==80,'Wrong Batch 6 coefficient coverage')
    by_coeff={}
    for x in coeff:
        key=(x['scenario'],x['model'])
        require(key[0] in SCENARIOS and key[1] in MODEL_NAMES and key not in by_coeff,
                'Duplicate/missing B6 scenario model')
        by_coeff[key]=x
        checked_float(x['coefficient_residual'],'B6 coefficient error')
        checked_float(x['max_group_fractional_spread'],'B6 group spread')
    for x in coeff_csv:
        key=(x['scenario'],x['model'])
        require(key in by_coeff,'B6 coefficient CSV key mismatch')
        rec=by_coeff[key]
        require(x['reference_model']=='A_-2' and
                close(x['coefficient_residual'],rec['coefficient_residual'],rel=1e-12,abs_tol=1e-15) and
                close(x['max_group_fractional_spread'],rec['max_group_fractional_spread'],rel=1e-12,abs_tol=1e-15),
                'B6 coefficient JSON/CSV inconsistent '+str(key))
    require(set(by_coeff)=={(s,n) for s in SCENARIOS for n in MODEL_NAMES},
            'B6 scenario grid not complete')
    reconstruction=b6['tensor_reconstruction_checks']
    require(len(reconstruction)==80,'Wrong B6 reconstruction grid')
    rec_keys=set()
    max_rec_error=0.0
    for x in reconstruction:
        key=(x['scenario'],x['model'])
        require(key in by_coeff and key not in rec_keys,'B6 reconstruction key mismatch')
        rec_keys.add(key)
        err=checked_float(x['exact_tensor_reconstruction_error'],'B6 reconstruction')
        require(err<1e-8,'B6 reconstruction no longer numerically exact')
        max_rec_error=max(max_rec_error,err)
    require(rec_keys==set(by_coeff),'B6 reconstruction coverage mismatch')
    for item in b6['model_summary']:
        name=item['model']
        require(name in MODEL_NAMES and
                close(item['baseline_coefficient_residual'],
                      by_coeff[('baseline',name)]['coefficient_residual'],
                      rel=1e-10,abs_tol=1e-15),
                'B6 model_summary not consistent')
    require(len(b6['model_summary'])==16,'B6 model summary incomplete')
    scenario_rows=[]
    for s in SCENARIOS:
        data=[by_coeff[(s,n)] for n in MODEL_NAMES]
        scenario_rows.append({'scenario':s,'models':16,
                              'max_coefficient_residual':max(checked_float(x['coefficient_residual'],'scenario') for x in data),
                              'max_group_fractional_spread':max(checked_float(x['max_group_fractional_spread'],'group') for x in data),
                              'restricted_UV_spectrum':str(s!='charge_split_diagnostic'),
                              'source_conditional_only':'True'})
    require(scenario_rows[0]['max_coefficient_residual']<1e-8,
            'Common-mixing B6 baseline coefficient comparison unexpectedly poor')
    baseline_err={n:by_coeff[('baseline',n)]['coefficient_residual'] for n in MODEL_NAMES}
    for res in results:
        res['b6_baseline_coefficient_residual']=baseline_err[res['model']]
    digest_info={
        'study_sha256':m5['study_sha256'],
        'B5_sha256':b6['source_sha256']['batch05'],
        'B5_manifest_sha256':b6['source_sha256']['batch05_manifest'],
        'B6_sha256':sha256(p_b6),
        'B6_manifest_sha256':sha256(p_b6manifest),
        'B6_csv_sha256':sha256(root/B6/'batch06_coefficient_tests.csv'),
        'study_csv_sha256':sha256(root/STUDY/'all16_summary.csv'),
        'scans_csv_sha256':sha256(root/STUDY/'scans.csv'),
        'hash_verification_modes':digests,
        'frozen_revision_recorded':study['provenance']['revision']['local_head'],
    }
    report_data={'results':results,'controls':controlled,'scenarios':scenario_rows,
                 'scan_by_model':scan_by_model,'scan_count':288,
                 'zero_checks':zero_checks,'baseline_checks':baseline_checks,
                 'max_linear_mass_residual':max_linear_mass_residual,
                 'max_tensor_reconstruction_error':max_rec_error,
                 'digests':digest_info}
    return report_data, render_files(report_data)


def render_files(d: dict) -> dict[str,str]:
    rows=d['results']
    resfields=('model','family','alpha','fermion_type','charge_channels',
               'm1_ev','m2_ev','m3_ev','m1_over_m3','m2_over_m3',
               'm1_geometry_A','m1_coherence_C','m1_mechanism',
               'source_reference_relative_error',
               'b6_baseline_coefficient_residual','physical_normalisation_certified')
    confields=('model','reference_texture','lambdaT3_abs','m1_ev','m2_ev',
               'm3_ev','m1_over_m3','m1_geometry_A','m1_coherence_C',
               'counterfactual_not_fitted')
    scenfields=('scenario','models','max_coefficient_residual',
                'max_group_fractional_spread','restricted_UV_spectrum',
                'source_conditional_only')
    contents={'all16_provisional_results.csv':_csv_table(rows,resfields),
              'controlled_texture_comparisons_64.csv':_csv_table(d['controls'],confields),
              'coefficient_scenario_summary_5.csv':_csv_table(d['scenarios'],scenfields)}
    md=['# Final numerical closeout — T3, 16 provisional one-loop benchmarks',
        '', '**Numerical phase: CLOSED FOR SOURCE-CONDITIONAL RESULTS ONLY.**',
        '', '**Not a physical neutrino-mass certificate:** the Matchete/RGBeta →',
        'left-Weyl signed action normalisation and basis/phase conversion remain',
        'unresolved. Physical certification is **0/16**. No EFT, RG running,',
        'threshold matching, oscillation fit or experimentally validated masses',
        'are inferred by this closeout.',
        '', '## Frozen results (all masses in eV; source-conditional)',
        '', '| Model | F | m1 | m2 | m3 | m1/m3 | A1 | C1 | Mechanism | B6 baseline eps |',
        '|---|---|---:|---:|---:|---:|---:|---:|---|---:|']
    for x in rows:
        md.append('| '+ ' | '.join((x['model'],x['fermion_type'],
             _fmt(x['m1_ev']),_fmt(x['m2_ev']),_fmt(x['m3_ev']),
             _fmt(x['m1_over_m3']),_fmt(x['m1_geometry_A']),
             _fmt(x['m1_coherence_C']),x['m1_mechanism'],
             _fmt(x['b6_baseline_coefficient_residual'],3)))+' |')
    low=min(rows,key=lambda x:x['m1_over_m3'])
    high=max(rows,key=lambda x:x['m1_over_m3'])
    example=next(x for x in rows if x['model']=='C_+1')
    md += ['', 'Values are the **original frozen** candidate left-Weyl',
           'one-loop results. `A1`, `C1` and the mechanism labels are the',
           'Batch 5 source-conditional mode diagnostics; `eps` is the',
           'Batch 6 Yukawa-independent coefficient residual at a common',
           '|lambdaT3|=0.10 comparison (not the frozen-model physical error).',
           '', '## Findings justified by the stored numerics',
           '', f'- **16/16** frozen candidate mass matrices evaluated;',
           f'  min m1/m3 = {_fmt(low["m1_over_m3"])} ({low["model"]});',
           f'  max m1/m3 = {_fmt(high["m1_over_m3"])} ({high["model"]}).',
           f'- **C_+1:** m1/m3 = {_fmt(example["m1_over_m3"])};',
           f'  heavy-eigenspace coherence C1 = {_fmt(example["m1_coherence_C"])}.',
           '  The saved phasor analysis attributes its small lightest',
           '  singular value to destructive heavy-eigenspace interference',
           '  within the candidate Weyl convention.',
           f'- **288/288** recorded scan points passed: 16 models ×',
           '  [5 y1 + 5 y2 + 5 lambdaT3 + 3 heavy-mass points].',
           f'  All {d["zero_checks"]} y1/y2/lambdaT3 zero points vanish;',
           f'  all {d["baseline_checks"]} t=1 points reproduce the frozen',
           '  values. Max positive Yukawa-rescaling relative mass error:',
           f'  {_fmt(d["max_linear_mass_residual"],3)}.',
           '- **64** common-mixing Yukawa-texture comparisons are preserved',
           '  separately in `controlled_texture_comparisons_64.csv`.',
           '  These are *counterfactual* substitutions, not new fitted models;',
           '  they isolate which texture effects persist for fixed candidate',
           '  representation/fermion constructions.',
           f'- **80** Batch 6 coefficient comparisons and 80 independent',
           f'  tensor-reconstruction checks passed the recorded grid.',
           '  Maximum reconstruction residual:',
           f'  {_fmt(d["max_tensor_reconstruction_error"],3)}.',
           '', '## Batch 6 coefficient-comparison scenarios',
           '', '| Scenario | Max eps | Max heavy-group spread | Interpretation |',
           '|---|---:|---:|---|']
    for s in d['scenarios']:
        label=('Diagnostic only; SU(2)-split masses outside restricted ansatz'
               if s['scenario']=='charge_split_diagnostic'
               else 'Conditional fixed-source comparison')
        md.append('| '+ ' | '.join((s['scenario'],
              _fmt(s['max_coefficient_residual'],4),
              _fmt(s['max_group_fractional_spread'],4),label))+' |')
    md += ['', 'Near-zero coefficient residual at the chosen masses establishes',
           'numerical proportionality of the stored candidate coefficient',
           'tensors **at those parameters**; it is not proof for arbitrary',
           'masses or a physical cross-package normalisation.',
           '', '## Source and evidence boundaries',
           '', 'The run compares the frozen summary JSON/CSV, all 288 scan CSV',
           'rows, Batch 5 mass/texture diagnostics, and Batch 6 coefficient',
           'JSON/CSV and reconstruction grid. It checks the recorded',
           'study→Batch 5→Batch 6 SHA256 lineage, including Git CRLF-only',
           'checkout normalization when applicable. It does **not** remeasure',
           'native Matchete/RGBeta source actions or claim current-HEAD',
           'equivalence to the frozen source revision.',
           '', '| Input/evidence | SHA256 |', '|---|---|']
    dig=d['digests']
    for label,key in [('Frozen study JSON','study_sha256'),
                      ('Batch 5 JSON','B5_sha256'),
                      ('Batch 5 manifest','B5_manifest_sha256'),
                      ('Batch 6 JSON','B6_sha256'),
                      ('Batch 6 manifest','B6_manifest_sha256'),
                      ('Frozen scan CSV','scans_csv_sha256')]:
        md.append(f'| {label} | `{dig[key]}` |')
    md += ['', f'Original frozen-study recorded commit: `{dig["frozen_revision_recorded"]}`.',
           'Its earlier revision is retained as provenance, not imposed as the',
           'current local checkout HEAD.',
           '', '## Stop rule',
           '', '**The requested rough numerical milestone is complete** if',
           'this closeout validates against the Windows checkout. Do not',
           'continue SU(2) phase, physical normalisation or EFT work as part',
           'of this milestone. Reopen those only if independently required',
           'by the thesis scope or by a specific physical claim.', '']
    contents['FINAL_NUMERICAL_CLOSEOUT.md']='\n'.join(md)
    return contents


def write(root: Path, output: Path) -> tuple[dict,list[str]]:
    root=root.resolve()
    target=(output if output.is_absolute() else root/output).resolve()
    allowed=(root/BASE).resolve()
    require(target != allowed and allowed in target.parents,
            'Output must be a new directory beneath '+str(BASE))
    require(not target.exists(), 'Refusing to overwrite existing output: '+str(target))
    data,contents=compute(root)
    require(not target.exists(), 'Output appeared during validation: '+str(target))
    target.mkdir(parents=True,exist_ok=False)
    for name,body in contents.items():
        (target/name).write_text(body,encoding='utf-8',newline='\n')
    return data,list(contents)
