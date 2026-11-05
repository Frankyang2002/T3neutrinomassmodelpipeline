"""D(-2) and A(-4) fixed-Yukawa transfer for the uniform-0.01 T3 benchmark.

Use NEW accepted uniform-0.01 reference Yukawas Y1,Y2 as frozen UV matrices,
transfer each pair to all 16 representations without refitting, and re-evolve
using the existing UV, EFT1, conditional one-scale C5 and low-energy evaluators.

This is NOT the same as the separately refitted 16 reference models, and is
NOT a certified sequential UV -> all fermions -> all scalars -> SMEFT match.
No existing physics files, original fits or GitHub data are modified.

--frozen-symbolic is explicitly opt-in for reproducibility of the EXISTING,
independently refitted 0.01 dataset against its identical saved canonical
beta/C5 files, even when the repository's live symbolic-source fingerprints
have since changed. The unmodified canonical source-fingerprint gate stays
intact; this mode validates cached symbolic structure, accepted campaign
source SHA256s, and a fresh baseline physics evaluation. It is NOT a claim
that cached betas correspond to today's edited symbolic source files.

Run from project root:
 python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --preflight
 python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --run-id fixed_yukawa_001_v1
 python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --run-id fixed_yukawa_001_v1 --resume --frozen-symbolic
 python -X utf8 -m scripts.BuildT3Uniform001FixedYukawaTables --run-id fixed_yukawa_001_v1 --render-only
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path, PureWindowsPath
import re
import shutil
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001/fixed_yukawa_transfer'
REPORT_ROOT = ROOT / 'Reports/output/thesis_general_5_1/uniform001_fixed_yukawa_tables'
REF_ROOT = (ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001/campaigns'
            / 'uniform001_all16_reference_v1')
SOURCE_REFS = ('D_-2', 'A_-4')
MODELS = tuple(f'{family}_{i:+d}' for family, vals in {
    'A':(-4,-2,0),'B':(-3,-1,1),'C':(-3,-1,1),
    'D':(-2,0,2),'E':(-4,-2,0,2)}.items() for i in vals)
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
       'sin2_theta13','sin2_theta23')
SCHEMA = 't3-uniform001-fixed-yukawa-two-reference-transfer-v1'


def tag(model: str) -> str:
    return model.replace('-', 'm').replace('+','p')


def title(model: str) -> str:
    family, alpha = model.split('_')
    return f'{family}({int(alpha):+d})'.replace('+0','0')


def load(path: Path) -> dict:
    value=json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(value,dict):
        raise ValueError(f'Expected JSON object: {path}')
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def c5_from_json(value: dict) -> np.ndarray:
    if not isinstance(value,dict) or 'real' not in value or 'imag' not in value:
        raise ValueError('C5 must be a complex real/imag matrix')
    real=np.asarray(value['real'],dtype=float)
    imag=np.asarray(value['imag'],dtype=float)
    if real.shape!=(3,3) or imag.shape!=(3,3):
        raise ValueError('C5 must be 3x3')
    c=real+1j*imag
    if not np.isfinite(c).all():
        raise ValueError('C5 not finite')
    if not np.allclose(c,c.T,rtol=1e-7,atol=1e-24):
        raise ValueError('C5 not symmetric in the saved flavour basis')
    return c


def low_c5(details: dict) -> dict:
    record=details.get('c5_at_low_scale',details.get('C5_low_gev_inv'))
    if record is None:
        raise ValueError('Missing low-scale C5 in full-general evaluator')
    c5_from_json(record)
    return deepcopy(record)


def proportion(target: np.ndarray, reference: np.ndarray) -> tuple[complex,float]:
    """Least-squares complex k and normalized residual in fixed flavour basis."""
    if target.shape!=(3,3) or reference.shape!=(3,3):
        raise ValueError('Expected two 3x3 low-scale C5 matrices')
    den=float(np.vdot(reference,reference).real)
    target_norm=float(np.linalg.norm(target))
    if not math.isfinite(den) or den<=0 or target_norm<=0:
        raise ValueError('Zero/nonfinite C5 matrix cannot be compared')
    k=complex(np.vdot(reference,target)/den)
    eps=float(np.linalg.norm(target-k*reference)/target_norm)
    if not all(math.isfinite(x) for x in (k.real,k.imag,eps)):
        raise ValueError('Nonfinite C5 proportionality')
    if eps<-1e-13 or eps>1+1e-10:
        raise ValueError('C5 residual outside [0,1]')
    return k,max(0.0,min(eps,1.0))


def result_c5(result: dict) -> np.ndarray:
    return c5_from_json(result['C5_low_gev_inv'])


def signature(model: str, ref: str, source_paths: dict, reference_files: dict,
              lambda_t3: float) -> dict:
    return {'model':model,'yukawa_reference':ref,'lambdaT3':lambda_t3,
            'original_source_sha256':{k:sha(v) for k,v in sorted(source_paths.items())},
            'fitted_yukawa_reference_sha256':{k:sha(v) for k,v in sorted(reference_files.items())}}


def reference_files(reference_root: Path, ref: str) -> dict[str,Path]:
    folder=reference_root/'reference'/tag(ref)
    return {'config':folder/'fitted_config.json','result':folder/'result.json'}


def source_inputs(model: str, canonical_run: str, selection, zee):
    ns=argparse.Namespace(fit_run_id='merged', canonical_run_id=canonical_run)
    return zee.sources(model,ns,selection[model])


# The cached Matchete input-form paths are absolute paths from the machine
# that exported the canonical snapshot (for example C:\\thesistemp\\...).
# In frozen mode, rebase only those three references for structural validation.
# The original cached files and their recorded SHA256 hashes are immutable.
PINNED_INPUTFORMS = {
    'interaction_inputform_file': 'lagrangian.json.interaction.inputform.txt',
    'uv_lagrangian_inputform_file': 'lagrangian.json.uv.inputform.txt',
    'bsm_lagrangian_inputform_file': 'lagrangian.json.bsm.inputform.txt',
}
PINNED_STATUS_SHA256 = {
    'lagrangian.json': 'lagrangian_sha256',
    'lagrangian.json.interaction.inputform.txt': 'lagrangian_interaction_sha256',
    'lagrangian.json.uv.inputform.txt': 'lagrangian_uv_sha256',
    'uv_rgbeta.json': 'uv_rgbeta_sha256',
    'eft1_rgbeta.json': 'eft1_rgbeta_sha256',
}


def validate_relocated_pinned_cache(canonical, model: str, folder: Path,
                                    entry: dict, kind: str, pinned: dict) -> None:
    """Validate immutable saved export with relocated absolute INPUTFORM paths.

    Temporary JSON metadata is for the EXISTING validator only. This does not
    rewrite the saved canonical files, skip its structural/RGE checks, or claim
    its source snapshot matches today's symbolic Lagrangian.
    """
    lag_path = folder / 'lagrangian.json'
    lag = load(lag_path)
    status = load(folder / 'status.json')
    hashes = status.get('sources')
    if not isinstance(hashes, dict):
        raise ValueError(f'{model}: canonical status lacks source SHA256s')
    for filename, key in PINNED_STATUS_SHA256.items():
        f = folder / filename
        if not f.is_file() or f.stat().st_size == 0:
            raise FileNotFoundError(f'{model}: missing pinned symbolic file {f}')
        expected = hashes.get(key)
        if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-f]{64}', expected):
            raise ValueError(f'{model}: missing/invalid pinned SHA256 for {filename}')
        if sha(f) != expected:
            raise ValueError(f'{model}: pinned symbolic SHA256 mismatch for {filename}')
    relocated = deepcopy(lag)
    for key, expected_name in PINNED_INPUTFORMS.items():
        recorded = lag.get(key)
        if not isinstance(recorded, str) or not recorded:
            raise ValueError(f'{model}: missing symbolic inputform pointer {key}')
        # PureWindowsPath handles stored Windows paths on Windows and Linux.
        if PureWindowsPath(recorded).name != expected_name:
            raise ValueError(f'{model}: unexpected pinned inputform filename in {key}')
        local = folder / expected_name
        if not local.is_file() or local.stat().st_size == 0:
            raise FileNotFoundError(f'{model}: missing pinned symbolic inputform {local}')
        relocated[key] = str(local.resolve())
    # No change to any saved provenance, beta file, Lagrangian, or status JSON.
    # Validate a temporary COPY of the original metadata using the existing
    # canonical structural validator against the original pinned source hashes.
    with tempfile.TemporaryDirectory(prefix='t3-pinned-symbolic-check-') as tmp:
        shadow = Path(tmp)
        for name in ('provenance.json', 'status.json', 'uv_rgbeta.json', 'eft1_rgbeta.json'):
            shutil.copy2(folder/name, shadow/name)
        (shadow/'lagrangian.json').write_text(
            json.dumps(relocated, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        canonical.validate_cached(model, shadow, entry, kind, pinned)


class PinnedCanonicalNonMixedFactory:
    """Use saved, reference-hash-matched beta/C5 cache; never relabel live sources.

    Canonical symbolic structural checks are still called on the saved export,
    using its own SOURCE snapshot stamp. The current-source fingerprint check
    is intentionally not claimed; this opt-in mode is a frozen-cache
    reproducibility study, and every target baseline is re-evaluated below.
    """
    def __init__(self, model: str, original: dict, paths: dict, target,
                 canonical_run: str):
        from scripts import RunT3GeneralSymbolicStudy as canonical
        from scripts import RunT3CompleteGeneralScalarSymbolic as scalar
        from scripts import RunT3GeneralAll16FromReference as general
        from scripts import RunT3UniformScalarBenchmark as uniform
        from scripts.AuditGeneralT3QuadraticBasis import audit
        from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
            load_manifest, quadratic_kind)
        from Numerical.matching.GeneralT3FlavourC5 import GeneralT3FlavourC5

        folder = canonical.DEFAULT_DATA_ROOT / canonical_run / canonical.model_tag(model)
        stamp = canonical.load_json(folder / 'provenance.json')
        pinned = stamp.get('source_sha256')
        if stamp.get('model') != model or not isinstance(pinned, dict) or not pinned:
            raise ValueError(f'{model}: invalid pinned canonical source provenance')
        current = canonical.fingerprint()
        self.changed_symbolic_sources = sorted(
            name for name in set(pinned) | set(current)
            if pinned.get(name) != current.get(name))
        if self.changed_symbolic_sources:
            print(f'  {model}: opt-in frozen canonical export; live symbolic '
                  f'sources differ ({len(self.changed_symbolic_sources)} files)',flush=True)
        registry = canonical.load_json(scalar.REGISTRY)
        manifest = load_manifest(scalar.MANIFEST)
        entry = registry['models'][model]
        kind = quadratic_kind(manifest,model)
        if model not in audit()['models']:
            raise ValueError(f'{model}: missing quadratic model audit')
        validate_relocated_pinned_cache(canonical,model,folder,entry,kind,pinned)
        # The exact UV/EFT1 and matching bytes are independently checked by
        # campaign.validated_reference for EACH target model before use.
        self.uv = canonical.load_json(paths['uv_beta'])
        self.eft = canonical.load_json(paths['eft1_beta'])
        extra = canonical.expected_yukawa(model)
        actual = set(self.uv['report_betas']) & {'y1Conj','y2Conj'}
        if actual != ({extra} if extra else set()):
            raise ValueError(f'{model}: cached UV extra Yukawa inventory mismatch')
        if set(self.eft['report_betas']) & {'y1','y2','MF','y1Conj','y2Conj'}:
            raise ValueError(f'{model}: cached EFT1 still contains heavy fermions')
        matched = general.ensure_general_matching(model,generate=False,wolframscript='wolframscript')
        if matched.resolve() != paths['matching'].resolve():
            raise ValueError(f'{model}: cached Matchete coefficient mismatch')
        kernel = GeneralT3FlavourC5(model,matched)
        if extra is not None and extra not in kernel.required_yukawas:
            raise ValueError(f'{model}: cached C5 omits symmetry-allowed {extra}')
        self.model = model
        self.original = original
        self.target = target
        self.kernel = kernel
        self.frozen, self.all_quartics, self.all_conjugates = uniform.uniform_overrides(
            canonical.load_json(paths['overrides']))

    def build(self, lam: float):
        from scripts import RunT3GeneralAll16FromReference as general
        from scripts import RunT3UniformScalarBenchmark as uniform
        changed, ordinary_quartics = uniform.uniform_config(self.original,lam)
        fixed = deepcopy(self.frozen)
        fixed['lambdaT3'] = {'real':float(lam),'imag':0.0}
        evaluator = general.NonMixedEvaluator(self.model,self.uv,self.eft,
                                               self.kernel,fixed,self.target)
        boundary = {'ordinary_quartics':ordinary_quartics,
                    'all_quartics_in_uv_layout':self.all_quartics,
                    'conjugate_yukawas_zero':self.all_conjugates,
                    'extra_quartic_value':uniform.QUARTIC,
                    'matching':'full-general-matchete-one-scale-channels',
                    'symbolic_provenance_mode':'pinned-canonical-export',
                    'live_symbolic_source_fingerprint_matches':
                        not bool(self.changed_symbolic_sources)}
        return evaluator,changed,fixed,boundary


def verify_pinned_baseline(factory, accepted: dict, cfg: dict,
                           *, relative_tolerance: float=1e-5) -> None:
    """Reproduce target model's saved 0.01 reference before fixed-Y transfer.

    This catches stale/misaligned beta, matching or evaluator behaviour rather
    than accepting a cached file merely because its JSON exists.
    """
    from scripts import RunT3QualitativeScalarRefits as base
    value=base.real_lambda(cfg)
    evaluator,_,_,_=factory.build(value)
    x=base.fitted_vector(cfg)
    fit,details=evaluator.evaluate(base.parameter_dict(x),detail=True)
    observed=base.fit_prediction(fit)
    expected=accepted['prediction']
    for name in OBS:
        a,b=float(observed[name]),float(expected[name])
        if not math.isfinite(a) or abs(a-b)>max(1e-10,relative_tolerance*abs(b)):
            raise ValueError(f'{factory.model}: frozen symbolic reference baseline differs '
                             f'for {name}: recomputed={a:g}, saved={b:g}')
    new_c5=c5_from_json(low_c5(details))
    old_c5=c5_from_json(low_c5(accepted['numerical_details']))
    error=float(np.linalg.norm(new_c5-old_c5)/np.linalg.norm(old_c5))
    if not math.isfinite(error) or error>relative_tolerance:
        raise ValueError(f'{factory.model}: frozen symbolic reference C5 baseline '
                         f'differs from accepted 0.01 result (relative={error:.3g})')
    print(f'  {factory.model}: frozen canonical baseline reproduced '
          f'(C5 relative deviation {error:.3g})',flush=True)


def reference_inputs(ref: str, root: Path, originals: dict, zee, campaign,
                     chi2_cutoff: float):
    orig, _, paths=originals[ref]
    row,cfg=campaign.validated_reference(
        root/'reference'/tag(ref),model=ref,source_paths=paths,
        manifest_path=zee.MERGED_FILE,original_config=orig,chi2_cutoff=chi2_cutoff)
    return row,cfg


def reuse_diagonal(saved: dict, model: str, ref: str, sig: dict) -> dict:
    """The self-transfer is identical to the independently fitted new reference."""
    details=saved['numerical_details']
    c5=low_c5(details)
    return {'schema':SCHEMA,'model':model,'yukawa_reference':ref,
            'status':'ReusedReference','fit_type':'fixed-Yukawa, no new minimisation',
            'lambdaT3':sig['lambdaT3'],'chi2':float(saved['chi2']),
            'prediction':dict(saved['prediction']),
            'neutrino_masses_ev':list(saved['neutrino_masses_ev']),
            'C5_low_gev_inv':c5,'C5_scale_gev':100.0,
            'numerical_details_source':'accepted-new-reference-fit',
            'matching_scope':'conditional-one-scale-not-certified-sequential',
            'yukawa_conjugate_uv_zero':True,'extra_scalar_quartics_uv':0.01,
            'signature':sig}


def evaluate_fixed(factory, model: str, ref: str, ref_cfg: dict, sig: dict):
    from scripts import RunT3QualitativeScalarRefits as base
    lam=sig['lambdaT3']
    evaluator,changed,_,boundary=factory.build(lam)
    x=base.fitted_vector(ref_cfg)  # exactly the source model's 18 fitted real Yukawas
    frozen=base.yukawa_matrices(x)
    changed['base_state']['ordinary'].update(deepcopy(frozen))
    fit,detail=evaluator.evaluate(base.parameter_dict(x),detail=True)
    p=base.fit_prediction(fit)
    m=base.low_masses(detail,changed,str(fit.ordering))
    if len(m)!=3 or not np.isfinite(m).all():
        raise ValueError('Invalid low-scale masses')
    c5=low_c5(detail)
    return {'schema':SCHEMA,'model':model,'yukawa_reference':ref,
            'status':'Success','fit_type':'fixed-Yukawa, no new minimisation',
            'lambdaT3':lam,'chi2':float(fit.chi2),'prediction':p,
            'neutrino_masses_ev':m,'C5_low_gev_inv':c5,'C5_scale_gev':100.0,
            'uv_transferred_yukawas':frozen,'scalar_boundary':boundary,
            'matching_scope':'conditional-one-scale-not-certified-sequential',
            'yukawa_conjugate_uv_zero':True,'extra_scalar_quartics_uv':0.01,
            'signature':sig}


def fmt(v: float, precision: int=6) -> str:
    return f'{float(v):.{precision}g}'


def render_markdown(headers: list[str], body: list[list[str]], caption: str) -> str:
    return (f'### {caption}\n\n'
            +'| '+' | '.join(headers)+' |\n'
            +'| '+' | '.join(['---']*len(headers))+' |\n'
            +''.join('| '+' | '.join(x)+' |\n' for x in body))


def tex_escape(x: str) -> str:
    return (x.replace('\\','\\textbackslash{}').replace('_','\\_')
            .replace('%','\\%').replace('&','\\&'))


def render_tex(headers: list[str], body: list[list[str]], caption: str,
               label: str) -> str:
    n=len(headers)
    out=['\\begin{table}[htbp]','\\centering',
         '\\small',f'\\caption{{{tex_escape(caption)}}}',
         f'\\label{{tab:{label}}}',f'\\begin{{tabular}}{{{"l"+"r"*(n-1)}}}',
         '\\hline', ' & '.join(tex_escape(x) for x in headers)+r' \\',
         '\\hline']
    out.extend(' & '.join(tex_escape(v) for v in row)+r' \\' for row in body)
    out+=['\\hline','\\end{tabular}','\\end{table}']
    return '\n'.join(out)+'\n'


def write_table(dest: Path, name: str, rows: list[dict],
                fields: list[str], headers: list[str], caption: str):
    dest.mkdir(parents=True,exist_ok=True)
    with (dest/(name+'.csv')).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows(rows)
    body=[[fmt(row[k]) if isinstance(row.get(k),float) else str(row.get(k,''))
           for k in fields] for row in rows]
    (dest/(name+'.md')).write_text(render_markdown(headers,body,caption)+'\n',encoding='utf-8')
    (dest/(name+'.tex')).write_text(render_tex(headers,body,caption,name),encoding='utf-8')


def finished_rows(data_dir: Path, models=MODELS, refs=SOURCE_REFS) -> dict:
    out={}
    for model in models:
        for ref in refs:
            path=data_dir/tag(model)/tag(ref)/'result.json'
            if not path.is_file():
                raise FileNotFoundError(f'Missing physical result: {path}')
            record=load(path)
            if record.get('schema')!=SCHEMA or record.get('model')!=model or \
                    record.get('yukawa_reference')!=ref or \
                    record.get('status') not in ('Success','ReusedReference'):
                raise ValueError(f'Incomplete or wrong fixed-Yukawa result: {path}')
            if (record.get('extra_scalar_quartics_uv')!=0.01 or
                record.get('yukawa_conjugate_uv_zero') is not True):
                raise ValueError(f'Wrong benchmark in {path}')
            if not set(OBS).issubset(record.get('prediction',{})):
                raise ValueError(f'Missing mixing observables in {path}')
            if not len(record.get('neutrino_masses_ev',[]))==3:
                raise ValueError(f'Missing masses in {path}')
            result_c5(record)
            out[(model,ref)]=record
    return out


def build_tables(data_dir: Path, report_dir: Path, models=MODELS,
                 refs=SOURCE_REFS) -> dict:
    if tuple(refs)!=SOURCE_REFS:
        raise ValueError('Reference order must be D_-2, A_-4')
    got=finished_rows(data_dir,models,refs)
    c5refs={r:result_c5(got[(r,r)]) for r in refs}
    mass_fields=['model','m1_ev','m2_ev','m3_ev','delta_m21_1e5_ev2',
                 'delta_m3l_1e3_ev2','chi2_5obs','epsilon_C5']
    mass_heads=['Model','m1 (eV)','m2 (eV)','m3 (eV)',
                'Dm21 (1e-5 eV2)','Dm3l (1e-3 eV2)','chi2','epsilon C5']
    pms=[];coeff=[];entries=[]
    for model in models:
        cr={'model':title(model)}
        for ref in refs:
            row=got[(model,ref)]
            c5=result_c5(row)
            k,eps=proportion(c5,c5refs[ref])
            r={'model':title(model),'source_model':model,'yukawa_reference':ref,
               'm1_ev':float(row['neutrino_masses_ev'][0]),
               'm2_ev':float(row['neutrino_masses_ev'][1]),
               'm3_ev':float(row['neutrino_masses_ev'][2]),
               'delta_m21_1e5_ev2':float(row['prediction']['delta_m21_sq_ev2'])*1e5,
               'delta_m3l_1e3_ev2':float(row['prediction']['delta_m3l_sq_ev2'])*1e3,
               'chi2_5obs':float(row['chi2']),'epsilon_C5':eps}
            label=tag(ref).replace('_','')
            cr['k_'+label+'_real']=k.real
            cr['k_'+label+'_imag']=k.imag
            cr['epsilon_'+label]=eps
            pms.append({'model':title(model),'yukawa_reference':title(ref),
                        'sin2_theta12':float(row['prediction']['sin2_theta12']),
                        'sin2_theta13':float(row['prediction']['sin2_theta13']),
                        'sin2_theta23':float(row['prediction']['sin2_theta23'])})
            for i in range(3):
                for j in range(3):
                    entries.append({'model':model,'yukawa_reference':ref,'i':i+1,'j':j+1,
                                    'C5_real_gev_inv':float(c5[i,j].real),
                                    'C5_imag_gev_inv':float(c5[i,j].imag)})
        coeff.append(cr)
    for ref in refs:
        data=[]
        for model in models:
            row=got[(model,ref)]
            k,eps=proportion(result_c5(row),c5refs[ref])
            data.append({'model':title(model),
                   'm1_ev':float(row['neutrino_masses_ev'][0]),
                   'm2_ev':float(row['neutrino_masses_ev'][1]),
                   'm3_ev':float(row['neutrino_masses_ev'][2]),
                   'delta_m21_1e5_ev2':float(row['prediction']['delta_m21_sq_ev2'])*1e5,
                   'delta_m3l_1e3_ev2':float(row['prediction']['delta_m3l_sq_ev2'])*1e3,
                   'chi2_5obs':float(row['chi2']), 'epsilon_C5':eps})
        name='table_C1_Dm2' if ref=='D_-2' else 'table_C2_Am4'
        write_table(report_dir,name,data,mass_fields,mass_heads,
                    f'Uniform 0.01 general-T3, fixed {title(ref)} Yukawas; low scale 100 GeV')
    pm_fields=['model','yukawa_reference','sin2_theta12','sin2_theta13','sin2_theta23']
    write_table(report_dir,'table_C3_PMNS',pms,pm_fields,
                ['Model','Fixed Yukawa reference','sin2 theta12','sin2 theta13','sin2 theta23'],
                'PMNS mixing observables from 32 fixed-Yukawa transfers')
    k_fields=['model','k_Dm2_real','k_Dm2_imag','epsilon_Dm2',
              'k_Am4_real','k_Am4_imag','epsilon_Am4']
    write_table(report_dir,'table_C4_Weinberg_proportionality',coeff,k_fields,
                ['Model','Re k(D-2)','Im k(D-2)','epsilon(D-2)',
                 'Re k(A-4)','Im k(A-4)','epsilon(A-4)'],
                'Low-scale Weinberg C5 proportionality in saved flavour basis at 100 GeV')
    with (report_dir/'table_C5_Weinberg_matrix_entries.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(entries[0]))
        w.writeheader();w.writerows(entries)
    provenance={'schema':SCHEMA,'models':list(models),'yukawa_references':list(refs),
        'points':len(got),'uniform_scalar_quartic_uv':0.01,
        'conjugate_yukawa_uv_zero':True,
        'C5_evaluation_scale_gev':100.0,
        'k_definition':'Tr(C5_ref^dagger C5_model)/Tr(C5_ref^dagger C5_ref)',
        'epsilon_definition':'FrobeniusNorm(C5_model-k*C5_ref)/FrobeniusNorm(C5_model)',
        'complex_k':True,
        'mass_squared_columns':{'delta_m21_1e5_ev2':'10^-5 eV^2',
                                'delta_m3l_1e3_ev2':'10^-3 eV^2'},
        'matching_scope':'existing conditional one-scale C5; not certified sequential EFT',
        'inputs':{f'{m}:{ref}':sha(data_dir/tag(m)/tag(ref)/'result.json')
                  for m in models for ref in refs}}
    save(report_dir/'provenance.json',provenance)
    return provenance


def export_reference_table(reference_root: Path, report_dir: Path) -> None:
    src=reference_root/'reference_comparison.csv'
    if not src.is_file():
        raise FileNotFoundError(src)
    with src.open(newline='',encoding='utf-8-sig') as f:
        raw=list(csv.DictReader(f))
    by={r['model']:r for r in raw}
    if set(by)!=set(MODELS) or len(raw)!=16:
        raise ValueError('Reference comparison must contain all 16 models')
    records=[]
    for m in MODELS:
        r=by[m]
        if r['status'] not in ('AcceptedApproxFit','AcceptedRescaledSeed') or float(r['chi2'])>5.86:
            raise ValueError(f'Unaccepted source model {m}')
        records.append({'model':title(m),
            'm1_ev':float(r['m1_ev']),'m2_ev':float(r['m2_ev']),
            'm3_ev':float(r['m3_ev']),
            'delta_m21_1e5_ev2':float(r['delta_m21_sq_ev2'])*1e5,
            'delta_m3l_1e3_ev2':float(r['delta_m3l_sq_ev2'])*1e3,
            'sin2_theta12':float(r['sin2_theta12']),
            'sin2_theta13':float(r['sin2_theta13']),
            'sin2_theta23':float(r['sin2_theta23']),
            'chi2_5obs':float(r['chi2'])})
    cols=['model','m1_ev','m2_ev','m3_ev','delta_m21_1e5_ev2',
          'delta_m3l_1e3_ev2','sin2_theta12','sin2_theta13','sin2_theta23','chi2_5obs']
    write_table(report_dir,'table_R0_individually_refitted',records,cols,
                ['Model','m1(eV)','m2(eV)','m3(eV)','Dm21(1e-5)','Dm3l(1e-3)',
                 'sin2 t12','sin2 t13','sin2 t23','chi2'],
                'Saved new uniform-0.01 INDIVIDUALLY REFITTED references (NOT fixed-Yukawa transfers)')


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='uniform001_Dm2_Am4_fixed_v1')
    p.add_argument('--reference-run-id',default='uniform001_all16_reference_v1')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--models',default='all',help='all, or comma-separated model IDs')
    p.add_argument('--preflight',action='store_true',help='Validate references without physics')
    p.add_argument('--reference-table-only',action='store_true',
                   help='Render existing 16 independently refitted models, no transfers')
    p.add_argument('--render-only',action='store_true',
                   help='Produce C1-C4 from already computed 32 fixed transfers')
    p.add_argument('--resume',action='store_true',help='Resume checkpointed transfer results')
    p.add_argument('--frozen-symbolic',action='store_true',
                   help='Opt in to original validated beta/C5 snapshot when live symbolic fingerprints have changed; baseline reproduced per target')
    p.add_argument('--chi2-cutoff',type=float,default=5.86)
    args=p.parse_args(argv)
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id):
        p.error('Invalid run ID')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.reference_run_id):
        p.error('Invalid reference run ID')
    if not (math.isfinite(args.chi2_cutoff) and args.chi2_cutoff>0):
        p.error('Invalid chi2 cutoff')
    wanted=MODELS if args.models=='all' else tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not wanted or len(set(wanted))!=len(wanted) or set(wanted)-set(MODELS):
        p.error('Invalid model selection')
    wanted=tuple(m for m in MODELS if m in wanted)
    root=REF_ROOT.parent/args.reference_run_id
    out=DATA_ROOT/args.run_id
    report=REPORT_ROOT/args.run_id
    if args.reference_table_only:
        export_reference_table(root,report)
        print(f'Wrote individually refitted reference table: {report}')
        return 0
    if args.render_only:
        export_reference_table(root,report)
        result=build_tables(out,report,wanted)
        print(f'Rendered fixed-Yukawa report from {result["points"]} results: {report}')
        return 0
    from scripts import RunT3GeneralZeeSensitivity as zee
    from scripts import RunT3QualitativeScalarRefits as base
    from scripts import RunT3UniformScalarBenchmark as uniform
    from scripts import RunT3UniformScalarCampaign as campaign
    if uniform.QUARTIC!=0.01:
        raise RuntimeError('The uniform scalar benchmark no longer uses quartic 0.01')
    selections=zee.fit_selections(list(set(wanted)|set(SOURCE_REFS)),'merged')
    source={m:source_inputs(m,args.canonical_run_id,selections,zee)
            for m in MODELS if m in set(wanted)|set(SOURCE_REFS)}
    refs={m:reference_inputs(m,root,source,zee,campaign,args.chi2_cutoff)
          for m in SOURCE_REFS}
    if args.preflight:
        for ref,(_,cfg) in refs.items():
            x=base.fitted_vector(cfg)
            print(f'{ref}: accepted, chi2={refs[ref][0]["chi2"]:.8g}, '
                  f'lambdaT3={base.real_lambda(cfg):g}, 18 finite Yukawa entries; '
                  f'Ymax={max(abs(x)):.6g}')
        for m in wanted:
            print(f'{m}: canonical source, UV/EFT1 betas and Matchete inputs present')
        print(f'Ready: {len(wanted)} models x 2 fixed Yukawa references; no physics run')
        return 0
    if out.exists() and not args.resume:
        p.error(f'Output exists: {out}. Use --resume or a new --run-id')
    out.mkdir(parents=True,exist_ok=True)
    target=zee.target_object()
    failed=0
    for model in wanted:
        original,fit,paths=source[model]
        factory=None
        lam=base.real_lambda(original)
        for ref in SOURCE_REFS:
            record_path=out/tag(model)/tag(ref)/'result.json'
            fpaths=reference_files(root,ref)
            sig=signature(model,ref,paths,fpaths,lam)
            if record_path.is_file():
                previous=load(record_path)
                if previous.get('signature') != sig:
                    raise ValueError(f'Cannot resume different source files: {record_path}')
                if previous.get('status') in ('Success','ReusedReference'):
                    print(f'{model} fixed {ref}: already saved',flush=True)
                    continue
            print(f'{model} fixed {ref}: running',flush=True)
            try:
                if model==ref:
                    row=reuse_diagonal(refs[ref][0],model,ref,sig)
                else:
                    if factory is None:
                        if args.frozen_symbolic and model not in base.MIXED:
                            # Require that this target model's own NEW 0.01
                            # reference actually used the identical current
                            # beta/C5 bytes. This verifies the source snapshot
                            # rather than trusting stale canonical provenance.
                            target_reference,target_config=campaign.validated_reference(
                                root/'reference'/tag(model),model=model,
                                source_paths=paths,manifest_path=zee.MERGED_FILE,
                                original_config=original,chi2_cutoff=args.chi2_cutoff)
                            factory=PinnedCanonicalNonMixedFactory(
                                model,original,paths,target,args.canonical_run_id)
                            verify_pinned_baseline(factory,target_reference,target_config)
                        else:
                            factory=uniform.ModelFactory(model,original,fit,paths,
                                                         target,args.canonical_run_id)
                    row=evaluate_fixed(factory,model,ref,refs[ref][1],sig)
                    row['symbolic_provenance_mode']=(
                        'pinned-canonical-export-baseline-reproduced'
                        if args.frozen_symbolic and model not in base.MIXED
                        else 'live-canonical-source-validated')
                save(record_path,row)
                print(f'  {row["status"]}: chi2={row["chi2"]:.7g}',flush=True)
            except Exception as exc:
                failed+=1
                error={'schema':SCHEMA,'model':model,'yukawa_reference':ref,
                       'status':'Failed','signature':sig,
                       'error':f'{type(exc).__name__}: {exc}'[:1200]}
                save(record_path,error)
                print(f'  FAILED: {error["error"]}',flush=True)
    if failed:
        print(f'Blocked: {failed} point(s) failed. Inspect saved result.json; '
              'run --resume after investigating, no fabricated table rows.')
        return 2
    if wanted != MODELS:
        print('Partial pilot saved. Run --models all --resume to complete 16.')
        return 0
    export_reference_table(root,report)
    result=build_tables(out,report)
    print(f'Complete {result["points"]}/32 fixed-Yukawa outputs, tables in: {report}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
