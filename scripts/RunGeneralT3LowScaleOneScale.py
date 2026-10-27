"""SM+Weinberg endpoint of already-calculated frozen-scalar one-scale C5.

This is a diagnostic, NOT fermion-first T3 EFT1 matching: SM-only RGEs
approximate the renormalisable background from the UV benchmark to the scalar
matching scale, then the existing SM+Weinberg RGE evolves C5 to the low scale.
The one-scale Matchete coefficient is imported unchanged at the scalar scale.

    python -m scripts.RunGeneralT3LowScaleOneScale --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / 'Numerical/benchmarks/fitted_reference'
SOURCE = ROOT / 'output/general_t3_numerical_c5_frozen_scalar_scenario'
OUTPUT = ROOT / 'output/general_t3_low_scale_one_scale'
MODELS = {
    'B_-1': 'T3_dS1_2_dS2_2_dF_1_alpha_m1',
    'C_-1': 'T3_dS1_2_dS2_2_dF_3_alpha_m1',
}
STATUS = 'DirectOneScaleC5SMWeinbergEndpointNotSequentialEFT1'
SOURCE_STATUS = 'FrozenScalarDiagonalDirectOneScaleC5NotPhysicalEFT1'


def fail(message: str) -> None:
    raise ValueError('T3 direct one-scale SM endpoint fail-closed: '+message)


def _read_json(path: Path) -> dict:
    if not path.is_file(): fail(f'missing {path}')
    payload = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(payload, dict): fail(f'expected JSON object: {path}')
    return payload


def _finite_positive(value, label: str) -> float:
    try: v = float(value)
    except (TypeError, ValueError) as e: raise ValueError(f'{label}: not numeric') from e
    if not np.isfinite(v) or v <= 0: fail(f'{label} must be finite and positive')
    return v


def _complex_matrix(value, label: str) -> np.ndarray:
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        fail(f'{label}: expected real/imag complex matrix')
    try:
        re = np.asarray(value['real'], dtype=float)
        im = np.asarray(value['imag'], dtype=float)
    except (ValueError, TypeError) as e: raise ValueError(f'{label}: nonnumeric matrix') from e
    if re.shape != (3, 3) or im.shape != (3, 3) or not np.isfinite(re).all() or not np.isfinite(im).all():
        fail(f'{label}: expected finite 3x3 matrix')
    return re + 1j * im


def _complex_entry_matrix(data, label: str) -> np.ndarray:
    if not isinstance(data, list) or len(data) != 3 or any(not isinstance(row,list) or len(row)!=3 for row in data):
        fail(f'{label}: expected 3x3 complex-entry list')
    try:
        out = np.array([[complex(float(z['real']),float(z['imag'])) for z in row]
                        for row in data],dtype=np.complex128)
    except (ValueError, KeyError, TypeError) as e: raise ValueError(f'{label}: malformed complex entries') from e
    if not np.isfinite(out.real).all() or not np.isfinite(out.imag).all(): fail(f'{label}: nonfinite complex matrix')
    return out


def _entries(matrix: np.ndarray) -> list[list[dict]]:
    return [[{'real':float(z.real),'imag':float(z.imag)} for z in row] for row in np.asarray(matrix)]


def _masses(c5: np.ndarray, vev: float) -> list[float]:
    return np.sort(np.linalg.svd((-0.5*vev**2*c5)*1e9,compute_uv=False))[::-1].tolist()


def _check_inputs(model: str, result: dict, frozen: dict, manifest: dict,
                  frozen_path: Path) -> tuple[np.ndarray, dict, float, float, float, float]:
    if model not in MODELS: fail(f'unsupported model {model}')
    expected=MODELS[model]
    rows=manifest.get('models')
    if not isinstance(rows,list) or len(rows)!=16: fail('expected complete 16-model benchmark manifest')
    match=[r for r in rows if r.get('model')==expected]
    if len(match)!=1: fail(f'{model}: frozen manifest entry missing or duplicated')
    r=match[0]
    expected_path=f'Numerical/benchmarks/fitted_reference/{expected}.json'
    if r.get('config_path')!=expected_path or r.get('sha256')!=sha256(frozen_path.read_bytes()).hexdigest():
        fail(f'{model}: frozen reference SHA256/path mismatch')
    dims=(2,2,1,-1) if model=='B_-1' else (2,2,3,-1)
    rep=frozen.get('representation',{})
    if tuple(rep.get(k) for k in ('d_s1','d_s2','d_f','alpha'))!=dims or rep.get('shared_scalar') is not False:
        fail(f'{model}: frozen representation mismatch')
    if (result.get('model')!=model or result.get('status')!=SOURCE_STATUS or
        result.get('frozen_scalar_aligned_one_scale_scenario') is not True or
        result.get('mass_basis_compatible_with_frozen_reference') is not True or
        result.get('source_mass_parameters_used_as_numeric_scalar_input') is not False or
        result.get('n_fermion_flavours')!=3 or
        result.get('sequential_EFT1_performed') is not False or
        result.get('physical_c5_certified') is not False):
        fail(f'{model}: source is not the uncertified frozen-scalar-aligned direct C5 result')
    cfg=result.get('benchmark_config',{})
    p=cfg.get('reference_provenance',{})
    if (cfg.get('format')!='t3-frozen-scalar-diagonal-direct-C5-input-v1' or
        cfg.get('physical_input_certified') is not False or
        cfg.get('sequential_EFT1') is not False or
        p.get('frozen_config_path')!=expected_path or
        p.get('frozen_sha256')!=r['sha256'] or
        p.get('physical_EFT1_matching_performed') is not False):
        fail(f'{model}: source frozen-reference provenance mismatch')
    if p.get('quadratic_mixing_source') not in {
        'zero_by_explicit_default_scenario_not_in_fitted_reference',
        'user_supplied_command_line_not_a_fitted_parameter',
    }:
        fail(f'{model}: source quadratic-mixing scenario provenance missing')
    c5=_complex_entry_matrix(result.get('symmetrized_C5_gev_inv'), 'symmetrized_C5_gev_inv')
    if not np.allclose(c5,c5.T,rtol=1e-10,atol=1e-17):
        fail(f'{model}: source C5 is not symmetric')
    v=_finite_positive(frozen.get('numerical',{}).get('vev_gev'),'vev_gev')
    if not np.isclose(v,_finite_positive(cfg.get('vev_gev'),'source vev'),rtol=0,atol=1e-10):
        fail(f'{model}: source and frozen vev differ')
    scales=frozen.get('scales',{})
    mf=_finite_positive(scales.get('mu_fermion_threshold_gev'),'fermion matching scale')
    ms=_finite_positive(scales.get('mu_scalar_threshold_gev'),'scalar matching scale')
    ml=_finite_positive(scales.get('mu_low_gev'),'low scale')
    mu_uv=_finite_positive(frozen.get('base_state',{}).get('mu_gev'),'UV scale')
    if not mu_uv > mf > ms > ml: fail(f'{model}: unexpected threshold ordering')
    base=frozen.get('base_state',{}).get('sm',{})
    if not isinstance(base,dict): fail('missing frozen SM state')
    if (frozen.get('benchmark_search',{}).get('enabled') is not False or
        frozen.get('sensitivity',{}).get('enabled') is not False):
        fail('frozen reference unexpectedly enables fitting/scan')
    return c5,base,v,mu_uv,ms,ml


def _state_from_base(base, K, state_cls):
    return state_cls(
        gY=float(base['gY']),g2=float(base['g2']),g3=float(base['g3']),
        lambdaH=float(base['lambdaH']),
        ye=_complex_matrix(base['ye'],'ye'),
        yu=_complex_matrix(base['yu'],'yu'),
        yd=_complex_matrix(base['yd'],'yd'),
        K=np.asarray(K, dtype=np.complex128),
    ).validated()


def evaluate(model: str, source_file: Path, frozen_path: Path,
             manifest_path: Path, *, evolution=None, state_cls=None) -> dict:
    source=_read_json(source_file)
    frozen=_read_json(frozen_path)
    manifest=_read_json(manifest_path)
    c5,base,v,mu_uv,mu_s,mu_low=_check_inputs(model,source,frozen,manifest,frozen_path)
    if evolution is None or state_cls is None:
        from Numerical.running.SMWeinbergEvolution import (
            SMWeinbergInitialConditions, evolve_sm_weinberg,
        )
        if evolution is None: evolution=evolve_sm_weinberg
        if state_cls is None: state_cls=SMWeinbergInitialConditions
    numerical=frozen.get('numerical',{})
    try:
        rtol=_finite_positive(numerical.get('rtol'),'rtol')
        atol=_finite_positive(numerical.get('atol'),'atol')
        # Pure SM background. No BSM, no EFT1 and no C5 above the scalar scale.
        uv=_state_from_base(base, np.zeros((3,3),dtype=complex), state_cls)
        sm_background=evolution(uv,mu_uv,mu_s,rtol=rtol,atol=atol)
        if not sm_background.solver_success:
            fail('SM-only background evolution was unsuccessful')
        if np.max(np.abs(sm_background.K))>1e-15:
            fail('SM-only background unexpectedly generated C5')
        final_initial=state_cls(
            gY=sm_background.gY,g2=sm_background.g2,g3=sm_background.g3,
            lambdaH=sm_background.lambdaH,ye=sm_background.ye,
            yu=sm_background.yu,yd=sm_background.yd,K=c5,
        ).validated()
        low=evolution(final_initial,mu_s,mu_low,rtol=rtol,atol=atol)
        if not low.solver_success: fail('SM+Weinberg endpoint evolution unsuccessful')
    except (FloatingPointError,OverflowError) as e:
        raise ValueError('RGE evolution encountered a nonfinite arithmetic result') from e
    c5_low=np.asarray(low.K,dtype=complex)
    if c5_low.shape!=(3,3) or not np.isfinite(c5_low.real).all() or not np.isfinite(c5_low.imag).all():
        fail('low-scale C5 matrix is invalid')
    if not np.allclose(c5_low,c5_low.T,rtol=1e-10,atol=1e-17):
        fail('low-scale C5 is not symmetric')
    m_hi=_masses(c5,v);m_low=_masses(c5_low,v)
    if not np.isclose(m_hi[0],source['largest_takagi_mass_eV'],rtol=1e-9,atol=1e-12):
        fail('source C5/Takagi value is internally inconsistent')
    return {
        'status':STATUS,
        'model':model,
        'source_result_sha256':sha256(source_file.read_bytes()).hexdigest(),
        'frozen_reference_sha256':sha256(frozen_path.read_bytes()).hexdigest(),
        'source_scalar_masses_gev':source.get('scalar_masses_gev'),
        'quadratic_mixing_scenario':source['benchmark_config']['reference_provenance']['quadratic_mixing_mu2_gev2'],
        'mu_UV_gev':mu_uv,'mu_scalar_boundary_gev':mu_s,'mu_low_gev':mu_low,
        'vev_gev_fixed':v,
        'masses_eV_at_scalar_boundary_descending':m_hi,
        'masses_eV_at_low_scale_descending':m_low,
        'relative_largest_mass_change':m_low[0]/m_hi[0]-1 if m_hi[0]>0 else None,
        'C5_scalar_boundary_gev_inv':_entries(c5),
        'C5_low_gev_inv':_entries(c5_low),
        'SM_only_background_evolved_from_UV_to_scalar_threshold':True,
        'SM_plus_Weinberg_evolved_from_scalar_threshold_to_low':True,
        'renormalizable_BSM_RGE_included_above_scalar_threshold':False,
        'fermion_first_EFT1_matching_performed':False,
        'scalar_threshold_physical_matching_performed':False,
        'Matchete_RGBeta_nonzero_quartic_map_certified':False,
        'physical_C5_prediction_certified':False,
        'reproduces_original_oscillation_fit':False,
        'limitations':'One-scale all-heavy Matchete C5 inserted at scalar boundary; SM-only above scalar boundary; not a fermion-first EFT1 calculation.',
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=tuple(MODELS))
    parser.add_argument('--source-dir',type=Path,default=SOURCE)
    parser.add_argument('--frozen-dir',type=Path,default=FROZEN)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    args=parser.parse_args()
    results=[]
    for model in (tuple(MODELS) if args.all else (args.model,)):
        slug=MODELS[model]
        src=args.source_dir/(model.replace('+','p').replace('-','m')+'_final_numerical_c5.json')
        frozen=args.frozen_dir/(slug+'.json')
        result=evaluate(model,src,frozen,args.frozen_dir/'benchmark_manifest.json')
        args.output_dir.mkdir(parents=True,exist_ok=True)
        destination=args.output_dir/(model.replace('+','p').replace('-','m')+'_final_low_scale_one_scale.json')
        destination.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        results.append({'model':model,'file':str(destination),
            'mass_at_boundary_eV':result['masses_eV_at_scalar_boundary_descending'],
            'mass_at_low_scale_eV':result['masses_eV_at_low_scale_descending'],
            'relative_largest_mass_change':result['relative_largest_mass_change']})
        print(f'{model}: SM-only endpoint, m_max({result["mu_low_gev"]:g} GeV)='
              f'{result["masses_eV_at_low_scale_descending"][0]:.9g} eV '
              f'(boundary {result["masses_eV_at_scalar_boundary_descending"][0]:.9g} eV); '
              'NOT fermion-first EFT1 matching',flush=True)
    summary=args.output_dir/'summary_final.json'
    summary.write_text(json.dumps({'status':STATUS,'models':results,
        'fermion_first_EFT1_matching_performed':False,'physical_C5_prediction_certified':False},
        indent=2)+'\n',encoding='utf-8')
    print('Summary:',summary)

if __name__=='__main__': main()
