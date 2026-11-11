"""Read-only all-16 LLSS -> Weinberg contraction strength diagnostic.

Project path: scripts/DiagnoseT3LLSSContractionAll16.py

Inputs (from existing completed diagnostics):
  corrected_llss16_v4/<model>/cases.json, case 1111
  Reports/.../fermion_threshold_all16/<model>/summary.json

No RGEs, fitting, threshold reconstruction or physics reruns occur.
The norm ratio kappa_T3 = ||Delta C5_log(lambdaT3)||_F / ||C_LLSS||_F
compares different tensor spaces; it is not an anomalous-dimension eigenvalue.
Saved matching conventions and caveats, especially B_-1/C_-1, are retained.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODELS = tuple(f'{f}_{a:+d}' for f, alphas in {
    'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
    'D': (-2, 0, 2), 'E': (-4, -2, 0, 2)
}.items() for a in alphas)
CAMPAIGNS = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
REPORTS = ROOT / 'Reports/output/thesis_qualitative/eft1_llss_contributions'


def tag(model: str) -> str:
    family, value = model.split('_')
    alpha = int(value)
    return f'{family}_{"m" if alpha < 0 else "p"}{abs(alpha)}'


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mat(raw: dict, label: str) -> np.ndarray:
    if not isinstance(raw, dict) or not {'real', 'imag'} <= set(raw):
        raise ValueError(f'{label}: missing complex matrix')
    a = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError(f'{label}: invalid 3x3 complex matrix')
    return a


def norm(a: np.ndarray) -> float:
    return float(np.linalg.norm(a))


def residual(a: np.ndarray, b: np.ndarray) -> float:
    return norm(a-b) / max(norm(b), 1e-40)


def projection(a: np.ndarray, total: np.ndarray) -> float:
    return float(np.vdot(total, a).real / norm(total)**2)


def write_csv(path: Path, rows: list[dict], columns: list[str]):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def analyse(model: str, case_path: Path, threshold_path: Path):
    saved = load(case_path)
    threshold = load(threshold_path)
    if saved.get('model') != model or saved.get('status') != 'Computed':
        raise ValueError(f'{model}: corrected campaign status/model mismatch')
    if threshold.get('model') != model:
        raise ValueError(f'{model}: threshold diagnostic model mismatch')
    if threshold.get('reference_case') != '1111':
        raise ValueError(f'{model}: threshold not compared to case 1111')
    if threshold.get('reference_campaign') is not None and threshold['reference_campaign'] != case_path.parent.parent.name:
        raise ValueError(f'{model}: threshold and campaign run IDs differ')
    saved_sha = threshold.get('reference_cases_sha256')
    if saved_sha and saved_sha != digest(case_path):
        raise ValueError(f'{model}: saved cases SHA256 differs from threshold reconstruction')
    c = saved['cases']['1111']
    total = mat(c['c5_induced_wilson_log'], f'{model} total log')
    finite = mat(c['c5_joint_scalar_bubble'], f'{model} finite bubble')
    full = mat(c['c5_at_scalar_threshold'], f'{model} C5(muS)')
    if residual(total + finite, full) > 1e-8:
        raise ValueError(f'{model}: total C5 != log + bubble')
    llss = float(threshold['C_LLSS_formal_frobenius_gev_inverse'])
    if not math.isfinite(llss) or llss <= 0:
        raise ValueError(f'{model}: nonpositive LLSS norm')
    if abs(llss*float(threshold['log_over_formal_LLSS_norm'])-norm(total)) > max(1e-7*norm(total), 1e-28):
        raise ValueError(f'{model}: threshold comparison log norm disagrees with cases')
    channels = c.get('per_general_interaction')
    if not isinstance(channels, list):
        raise ValueError(f'{model}: missing channel decomposition')
    quartics = [e for e in channels if e['category'] == 'quartic_decomposition']
    sources = [e for e in channels if e['category'] == 'fermion_tree_decomposition']
    if not quartics or not sources:
        raise ValueError(f'{model}: missing source or quartic category')
    native = [x for x in quartics if x['name'] == 'native:lambdaT3']
    if len(native) != 1:
        raise ValueError(f'{model}: expected exactly one lambdaT3 quartic channel')
    t3 = mat(native[0]['c5_leading_log'], f'{model} lambdaT3 log')
    other = total - t3
    summed = sum((mat(e['c5_leading_log'], e['name']) for e in quartics), np.zeros((3,3), complex))
    source_sum = sum((mat(e['c5_leading_log'], e['name']) for e in sources), np.zeros((3,3), complex))
    residual_quartic = residual(summed, total)
    residual_sources = residual(source_sum, total)
    if max(residual_quartic, residual_sources) > 1e-7:
        raise ValueError(f'{model}: saved source/quartic leading-log additivity failed')

    log_scale = math.log(1e5/1e3)/(16*math.pi**2)
    denom = norm(total)**2
    kappa_total = norm(total)/llss
    kappa_t3 = norm(t3)/llss
    kappa_other = norm(other)/llss
    # Signed projection sums to unity, unlike naive norm fractions.
    p_t3 = projection(t3,total)
    p_other = projection(other,total)
    # Coherence: 1 indicates parallel, -1 antiparallel, 0 orthogonal.
    cosine = float(np.vdot(t3,other).real/(norm(t3)*norm(other))) if norm(other)>0 else None
    provenance = c.get('quartic_provenance') or {}
    detail_map = {v['name']:v for v in provenance.get('quartic_channels',[]) if 'name' in v}
    aliases = provenance.get('assumed_native_coupling_aliases') or {}
    chrows = []
    for item in quartics:
        name = item['name']
        matrix = mat(item['c5_leading_log'], name)
        key = name.split(':',1)[1]
        metadata = detail_map.get(key, {})
        channel_lambda = None
        channel_imag = None
        if 'real' in metadata and 'imag' in metadata:
            channel_lambda = float(metadata['real'])
            channel_imag = float(metadata['imag'])
        chrows.append({
            'model':model, 'channel':name,
            'log_norm_GeV_inv':norm(matrix),
            'signed_projection_percent':100*projection(matrix,total),
            'norm_over_total':norm(matrix)/norm(total),
            'kappa_channel':norm(matrix)/llss,
            'coupling_real_at_muS_if_saved':channel_lambda,
            'coupling_imag_at_muS_if_saved':channel_imag,
            'coupling_name':aliases.get(key,key),
            'source_origin':metadata.get('origin',''),
            'conditional_log_normalization_factor':metadata.get('conditional_log_normalization_factor',1),
        })
    rows = {
        'model':model,
        'C_LLSS_norm_GeV_inv':llss,
        'C5_log_total_norm_GeV_inv':norm(total),
        'C5_log_lambdaT3_norm_GeV_inv':norm(t3),
        'C5_log_other_norm_GeV_inv':norm(other),
        'kappa_total':kappa_total,
        'kappa_lambdaT3':kappa_t3,
        'kappa_other':kappa_other,
        'kappa_lambdaT3_over_log_prefactor':kappa_t3/log_scale,
        'lambdaT3_signed_projection_percent':100*p_t3,
        'other_signed_projection_percent':100*p_other,
        'lambdaT3_norm_over_total':norm(t3)/norm(total),
        'other_norm_over_total':norm(other)/norm(total),
        'lambdaT3_other_cosine':cosine,
        'quartic_additivity_residual':residual_quartic,
        'fermion_source_additivity_residual':residual_sources,
        'lambdaT3_muS_coupling_value':'not_saved_in_cases',
        'fermion_source_count':len(sources),
        'quartic_channel_count':len(quartics),
        'mixed_scalar_muS_basis_assumption':model in ('B_-1','C_-1'),
    }
    return rows, chrows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='corrected_llss16_v4')
    p.add_argument('--models',default=','.join(MODELS))
    p.add_argument('--campaign',type=Path,default=None)
    p.add_argument('--threshold-dir',type=Path,default=None)
    p.add_argument('--output',type=Path,default=None)
    args=p.parse_args(argv)
    requested=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not requested or len(set(requested)) != len(requested) or set(requested)-set(MODELS):
        p.error('Invalid model IDs')
    if not args.run_id.replace('_','').replace('-','').isalnum():
        p.error('Invalid run ID')
    campaign=args.campaign or CAMPAIGNS/args.run_id
    threshold=args.threshold_dir or REPORTS/args.run_id/'fermion_threshold_all16'
    output=args.output or REPORTS/args.run_id/'lambdaT3_contraction_all16'
    output=output.resolve()
    if output==campaign.resolve() or output==threshold.resolve() or campaign.resolve() in output.parents or threshold.resolve() in output.parents:
        raise ValueError('Output must not overwrite input directories')
    rows=[]; channels=[]; errors=[]
    for model in requested:
        try:
            case=campaign/tag(model)/'cases.json'
            summary=threshold/tag(model)/'summary.json'
            row,ch=analyse(model,case,summary)
            rows.append(row);channels.extend(ch)
            print(f'{model}: kappa_T3={row["kappa_lambdaT3"]:.6e}; '
                  f'kappa_other={row["kappa_other"]:.6e}; '
                  f'T3 projection={row["lambdaT3_signed_projection_percent"]:.3f}%; '
                  f'other={row["other_signed_projection_percent"]:.3f}%',flush=True)
        except (FileNotFoundError,ValueError,KeyError,TypeError,ZeroDivisionError) as exc:
            errors.append({'model':model,'error':f'{type(exc).__name__}: {exc}'})
            print(f'{model}: FAILED {type(exc).__name__}: {exc}',flush=True)
    output.mkdir(parents=True,exist_ok=True)
    summary_columns=[
        'model','C_LLSS_norm_GeV_inv','C5_log_total_norm_GeV_inv',
        'C5_log_lambdaT3_norm_GeV_inv','C5_log_other_norm_GeV_inv',
        'kappa_total','kappa_lambdaT3','kappa_other','kappa_lambdaT3_over_log_prefactor',
        'lambdaT3_signed_projection_percent','other_signed_projection_percent',
        'lambdaT3_norm_over_total','other_norm_over_total','lambdaT3_other_cosine',
        'quartic_additivity_residual','fermion_source_additivity_residual',
        'lambdaT3_muS_coupling_value','fermion_source_count','quartic_channel_count',
        'mixed_scalar_muS_basis_assumption'
    ]
    channel_columns=['model','channel','log_norm_GeV_inv','signed_projection_percent',
                     'norm_over_total','kappa_channel','coupling_real_at_muS_if_saved',
                     'coupling_imag_at_muS_if_saved','coupling_name','source_origin',
                     'conditional_log_normalization_factor']
    write_csv(output/'all16_contraction_summary.csv',rows,summary_columns)
    write_csv(output/'all16_quartic_channels.csv',channels,channel_columns)
    manifest={
        'source_campaign':str(campaign),'threshold_diagnostic':str(threshold),
        'models_requested':list(requested),'models_computed':[x['model'] for x in rows],
        'errors':errors,'reference_case':'1111',
        'loop_prefactor_log_muF_over_muS_over_16pi2':math.log(100)/ (16*math.pi**2),
        'definitions':{
            'kappa_lambdaT3':'Frobenius norm of lambdaT3 C5 leading log divided by formal six-index LLSS norm',
            'kappa_other':'Frobenius norm of sum of all non-lambdaT3 C5 leading logs divided by formal LLSS norm',
            'signed_projection_percent':'100 Re Tr[C5_total^dagger C5_channel] / ||C5_total||_F^2',
            'lambdaT3_other_cosine':'Re Tr[C5_T3^dagger C5_other]/(||C5_T3||_F ||C5_other||_F)',
            'lambdaT3_muS_coupling_value':'NOT saved in corrected campaign cases; no guessed substitute',
        },
        'warnings':[
            'Kappa is a ratio of norms in different tensor spaces; not an anomalous dimension or an individual coupling.',
            'Kappa_other is a norm of the SUM of all other quartic contributions; signed projections are additive, norms are not.',
            'No RGE evolution or refitting; use existing saved fitted results only.',
            'Reference gammaCC-to-Weyl normalization, extra HHSS basis, and B/C conditional factors are uncertified.',
            'B_-1/C_-1 fermion-source tensor uses muS scalar mixing; not strictly a same-scale muF tensor.',
            'No complete one-loop hard matching at muF or complete EFT1 dimension-five anomalous dimensions.'
        ],
        'input_sha256':{
            model:{
                'cases':digest(campaign/tag(model)/'cases.json'),
                'threshold_summary':digest(threshold/tag(model)/'summary.json')
            } for model in requested if model in {row['model'] for row in rows}
        }
    }
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Output:',output)
    return 2 if errors else 0


if __name__=='__main__':
    raise SystemExit(main())
