"""B(-1)/C(-1) candidate mixed-mass C5 convention audit, read-only inputs.

Consumes NEW Batch 9 paired pilot results and historical fitted inputs.
Does not rerun RG, Matchete, fit, or matching. Evaluates the SAVED one-scale
ordered Matchete coefficient with the candidate RGBeta mixed-scalar endpoint.
Both identity and Stage23/26 CONDITIONAL MINUS maps are shown; neither is a
certified cross-package physical map. No sequential EFT1 C5 matching.

Usage:
  python -m scripts.AnalyzeT3BCMatchingConventions --all \
    --source-run-id pilot_q0p1_c5_unmixed --run-id matching_conventions_v1
  python -m scripts.AnalyzeT3BCMatchingConventions --all \
    --source-run-id pilot_q0p1_c5_unmixed \
    --mixing-run-id pilot_q0p1_mu0p01 --run-id matching_with_mu_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
from scripts import RunT3GeneralScalarBCPilot as b9

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'output/studies/general_scalar_extended'
EVIDENCE23 = ROOT / 'output/general_t3_stage22/stage23_quartic_prefactor_candidates.json'
EVIDENCE26 = ROOT / 'output/general_t3_stage22/stage26_conditional_quartic_action.json'
FORMAT = 't3-bc-conditional-convention-mixed-one-scale-preview-v1'
MAPS = {'identity_prior_preview': +1, 'stage23_conditional_minus': -1}
VEV_GEV = 246.22


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError('B/C convention audit: ' + reason)


def finite_float(value: Any, name: str) -> float:
    try:
        out = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f'B/C convention audit: {name} is not numeric') from exc
    require(math.isfinite(out), name + ' nonfinite')
    return out


def check_convention_evidence(stage23: dict, stage26: dict, *, models: tuple[str, ...]) -> dict:
    """Read exact signed 1D direct-C5 blocks; never promote conditional to physical."""
    require(stage23.get('format') == 'stage23-exact-source-prefactors-conditional-signed-maps-v1',
            'unexpected Stage23 report')
    require(stage23.get('source_prefactors_verified') is True and
            stage23.get('global_su2_component_basis_verified') is False and
            stage23.get('signed_physical_coupling_map_certified') is False and
            stage23.get('totals', {}).get('slots') == 50 and
            stage23.get('totals', {}).get('blocks') == 46,
            'Stage23 source/certification contract changed')
    require(stage26.get('format') == 'stage26-conditional-quartic-action-polynomials-v1' and
            stage26.get('source_relative_action_sign_supported') is True and
            stage26.get('signed_polynomial_identity_conditional_verified') is True and
            stage26.get('physical_coupling_mapping_certified') is False and
            stage26.get('models_checked') == 16 and
            stage26.get('slots_checked') == 50,
            'Stage26 conditional action/source contract changed')
    blocks = stage23.get('blocks')
    require(isinstance(blocks, list), 'Stage23 blocks missing')
    selected = {}
    for model in models:
        require(model in b9.MODELS, 'unsupported model')
        selected[model] = {}
        for rgb, match in b7.DIRECT.items():
            hits = [r for r in blocks if r.get('model') == model and
                    r.get('rgbeta_couplings') == [rgb] and
                    r.get('matchete_couplings') == [match]]
            require(len(hits) == 1, f'{model}: missing/ambiguous direct-C5 map for {rgb}')
            record = hits[0]
            require(record.get('block') in ('020020', '020200') and
                    record.get('r_in_m') == [['1']] and
                    record.get('lambda_rgbeta_from_lambda_matchete_if_assumptions_hold') == [['-1']] and
                    record.get('assumed_relative_potential_sign') == '-1' and
                    record.get('component_prefactor_one_ordered_orientation') == '1/4' and
                    record.get('physical_mapping_certified') is False,
                    f'{model}: conditional geometry/sign changed for {rgb}')
            selected[model][rgb] = {'matchete': match, 'block': record['block'],
                                    'geometry_class': record['geometry_class'],
                                    'r_in_m': '1', 'matchete_from_rgbeta_conditional': '-1'}
    return {'signed_map_is_conditional': True,
            'physical_mapping_certified': False, 'direct_blocks': selected}


def relative(reference: np.ndarray, shifted: np.ndarray) -> float:
    a = np.asarray(reference, dtype=complex)
    b = np.asarray(shifted, dtype=complex)
    require(a.shape == b.shape == (3,3), 'invalid C5 shapes')
    n = float(np.linalg.norm(a,'fro'))
    require(math.isfinite(n) and n > 0, 'zero/nonfinite reference C5')
    return float(np.linalg.norm(b-a,'fro')/n)


def decompose_effects(reference: np.ndarray, diagonal: np.ndarray,
                      mixed: np.ndarray, quartic_unmixed: np.ndarray,
                      joint: np.ndarray) -> dict:
    """Exact telescoping decomposition at *one* mapping sign."""
    ref = np.asarray(reference,dtype=complex)
    diag = np.asarray(diagonal,dtype=complex)
    mix = np.asarray(mixed,dtype=complex)
    q0 = np.asarray(quartic_unmixed,dtype=complex)
    combined = np.asarray(joint,dtype=complex)
    require(all(z.shape==(3,3) and np.isfinite(z.real).all() and
                np.isfinite(z.imag).all() and
                np.allclose(z,z.T,rtol=1e-10,atol=1e-17)
                for z in (ref,diag,mix,q0,combined)), 'invalid symmetric C5 matrix')
    scale = float(np.linalg.norm(ref,'fro'))
    require(math.isfinite(scale) and scale>0, 'zero C5 baseline')
    parts = {
        'diagonal_mass_and_lambdaT3_background':diag-ref,
        'generated_mixed_quadratic_at_zero_direct_quartics':mix-diag,
        'direct_quartics_at_unmixed_general_background':q0-diag,
        'quartic_mass_mixing_interference':combined-mix-q0+diag,
    }
    residual = combined-ref-sum(parts.values(),np.zeros((3,3),dtype=complex))
    residual_norm=float(np.linalg.norm(residual,'fro')/scale)
    require(residual_norm<2e-10, 'C5 decomposition closure failed')
    return {
        'relative_total':float(np.linalg.norm(combined-ref,'fro')/scale),
        'relative_background_diagonal':float(np.linalg.norm(parts['diagonal_mass_and_lambdaT3_background'],'fro')/scale),
        'relative_mixing_at_zero_direct_quartics':float(np.linalg.norm(parts['generated_mixed_quadratic_at_zero_direct_quartics'],'fro')/scale),
        'relative_direct_quartics_unmixed':float(np.linalg.norm(parts['direct_quartics_at_unmixed_general_background'],'fro')/scale),
        'relative_quartic_mass_mixing_interference':float(np.linalg.norm(parts['quartic_mass_mixing_interference'],'fro')/scale),
        'relative_direct_quartics_at_actual_mixing':float(np.linalg.norm(combined-mix,'fro')/scale),
        'relative_mixing_at_nonzero_direct_quartics':float(np.linalg.norm(combined-q0,'fro')/scale),
        'telescoping_residual_relative':residual_norm,
        'relative_terms_are_norms_not_additive':True,
    }


def endpoint(scene: dict) -> dict:
    e = scene.get('eft1_end')
    require(isinstance(e,dict) and e.get('mu_gev') == 1000.0,
            'missing 1 TeV EFT1 endpoint')
    m1,m2=finite_float(e.get('mS1Sq_gev2'),'mS1Sq'),finite_float(e.get('mS2Sq_gev2'),'mS2Sq')
    require(min(m1,m2)>0, 'nonpositive diagonal m²')
    mu=b7.complex_of(e.get('quadratic_mu2_gev2'))
    lam=b7.complex_of(e.get('lambdaT3'))
    require(math.isfinite(abs(mu)) and math.isfinite(abs(lam)), 'nonfinite quadratic/topology scalar')
    extra=e.get('reused_quartics')
    require(isinstance(extra,dict) and all(k in extra for k in b7.DIRECT),
            'direct quartics missing at EFT1 endpoint')
    q={k:b7.complex_of(extra[k]) for k in b7.DIRECT}
    return {'m1':m1,'m2':m2,'mu':mu,'lambdaT3':lam,'quartics':q}


def check_pilot_input(report:dict, case:str, model:str, source:dict, *,
                      require_preview:bool) -> dict:
    require(report.get('format') == b9.REPORT_FORMAT and
            report.get('case') == case and report.get('model') == model,
            'Batch9 model/case/status differs')
    require(report.get('historical_fitted_config_sha256') == b7.digest(source['cfg_path']) and
            report.get('historical_fit_record_sha256') == b7.digest(source['fit_path']) and
            report.get('fitted_yukawas_refitted') is False and
            report.get('historical_outputs_overwritten') is False and
            report.get('nonzero_quadratic_c5_computed') is False and
            report.get('matching_with_full_general_scalar_potential_completed') is False and
            report.get('candidate_basis') == 'Native RGBeta, physical Matchete map uncertified',
            'Batch9 provenance or scope flags changed')
    require(report.get('uv_on') == source['uv_on'] and
            report.get('passes_original_fit_screen') == source['fit_accepted'] and
            math.isclose(float(report.get('chi2_fit')),source['chi2'],rel_tol=1e-10,abs_tol=1e-10),
            'historical fit/UV mode changed')
    scenarios=report.get('scenarios',{})
    require(set(('restricted_zero','general_quartics')).issubset(scenarios),
            'Batch9 required paired scenarios absent')
    for name in ('restricted_zero','general_quartics'):
        require(scenarios[name].get('scenario')==name and
                scenarios[name].get('uv_integrated')==source['uv_on'],
                'Batch9 RG mode/provenance changed')
    for k in ('UV','EFT1'):
        f=b7.beta_path(model,k)
        require(report.get('beta_source_sha256',{}).get(b7.relpath(f))==b7.digest(f),
                'Batch9 beta file changed')
    require(report.get('registry_sha256') == b7.digest(b7.REGISTRY),
            'Batch9 scalar registry changed')
    if require_preview:
        preview=report.get('C5_preview')
        require(isinstance(preview,dict) and
                preview.get('status')=='UncertifiedOneScaleDirectQuarticSensitivity' and
                preview.get('c5_comparison_keeps_mu2_zero_by_projection') is True,
                'Batch9 validation preview unavailable')
    return scenarios


def historical_yukawa_fermion(source: dict) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    conf=source['cfg']['base_state']['ordinary']
    m=b9._three_by_three_complex(conf['MF'],'MF')
    masses=np.real(np.diag(m))
    require(np.all(masses>0) and np.max(np.abs(m-np.diag(masses)))<1e-8,
            'historical MF is not positive real diagonal')
    y1=b9._three_by_three_complex(conf['y1'],'y1')
    y2=b9._three_by_three_complex(conf['y2'],'y2')
    return masses,y1,y2


def make_matchete_evaluator(model: str, source: dict, coefficient: Path) -> Callable[...,np.ndarray]:
    """Evaluate the SAVED candidate Matchete expression with actual scalar U."""
    import mpmath as mp
    from scripts.RunGeneralT3NumericalC5 import compile_coefficient
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
    require(coefficient.is_file(), 'missing saved Matchete ordered source: '+str(coefficient))
    compiled=compile_coefficient(coefficient.read_text(encoding='utf-8'))
    masses,y1,y2=historical_yukawa_fermion(source)
    mp.mp.dps=55
    heavy=[mp.mpf(str(x)) for x in masses]
    fy1=[[mp.mpc(str(z.real),str(z.imag)) for z in row] for row in y1]
    fy2=[[mp.mpc(str(z.real),str(z.imag)) for z in row] for row in y2]

    def evaluate(scene:dict, quartics:dict[str,complex], *, sign:int,
                 unmixed:bool)->np.ndarray:
        require(sign in (-1,+1), 'map sign must be +/-1')
        mu=0j if unmixed else complex(scene['mu'])
        basis=MixedDoubletMassBasis.build(model,scene['m1'],scene['m2'],mu)
        require(min(basis.eigenvalues_sq)>0 and not basis.degenerate,
                'nonpositive or degenerate mixed scalar mass²; C5 blocked')
        U=basis.U
        env={f'T3MEU{i+1}{j+1}':mp.mpc(str(U[i,j].real),str(U[i,j].imag))
             for i in range(2) for j in range(2)}
        env.update({'MSX':mp.mpf(str(math.sqrt(basis.eigenvalues_sq[0]))),
                    'MSY':mp.mpf(str(math.sqrt(basis.eigenvalues_sq[1]))),
                    'lambdaT3':mp.mpc(str(scene['lambdaT3'].real),str(scene['lambdaT3'].imag)),
                    'hbar':1/(16*mp.pi**2), 'conj':mp.conj, 'sqrt':mp.sqrt, 'log':mp.log})
        for rgb,mat in b7.DIRECT.items():
            z=sign*complex(quartics[rgb])
            env[mat]=mp.mpc(str(z.real),str(z.imag))
        a=[[mp.mpc(0) for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                for k,MF in enumerate(heavy):
                    local=dict(env,MF=MF,y1a=fy1[i][k],y1b=fy1[j][k],
                                     y2a=fy2[i][k],y2b=fy2[j][k])
                    val=eval(compiled,{'__builtins__':{}},local)
                    require(bool(mp.isfinite(val.real) and mp.isfinite(val.imag)),
                            'nonfinite ordered Matchete coefficient')
                    a[i][j]+=val
        C=np.array([[complex(a[i][j]+a[j][i]) for j in range(3)]
                    for i in range(3)],dtype=complex)
        require(np.isfinite(C.real).all() and np.isfinite(C.imag).all() and
                np.allclose(C,C.T,rtol=1e-10,atol=1e-17),
                'nonfinite/non-symmetric evaluated C5')
        return C
    return evaluate


def matrix_json(a:np.ndarray)->dict:
    z=np.asarray(a,dtype=complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist()}


def model_output(model: str, report:dict, source:dict, evaluate:Callable,
                 *, mixing_scenario:dict|None = None) -> dict:
    scenarios=report['scenarios']
    zero=endpoint(scenarios['restricted_zero'])
    general=endpoint(scenarios['general_quartics'])
    require(zero['mu']==0j,'restricted mixed quadratic nonzero; baseline decomposition blocked')
    no_q={k:0j for k in b7.DIRECT}
    C00=evaluate(zero,no_q,sign=+1,unmixed=True)
    Cdiag=evaluate(general,no_q,sign=+1,unmixed=True)
    Cmixed=evaluate(general,no_q,sign=+1,unmixed=False)
    mappings={}
    for name,sign in MAPS.items():
        Cq_ref=evaluate(zero,general['quartics'],sign=sign,unmixed=True)
        Cq_unmixed=evaluate(general,general['quartics'],sign=sign,unmixed=True)
        Cjoint=evaluate(general,general['quartics'],sign=sign,unmixed=False)
        effects=decompose_effects(C00,Cdiag,Cmixed,Cq_unmixed,Cjoint)
        masses=np.sort(np.linalg.svd(-(VEV_GEV**2)/2*Cjoint,
                                      compute_uv=False))*1.0e9
        mappings[name]={
            'matchete_equals_rgbeta_times':sign,
            'relative_prior_Batch9_unmixed_projection':relative(C00,Cq_ref),
            'effects':effects,
            'C5_one_scale_candidate_gev_minus1':matrix_json(Cjoint),
            'C5_reference_one_scale_gev_minus1':matrix_json(C00),
            'neutrino_takagi_singular_values_eV_from_one_scale_candidate':masses.tolist(),
            'mass_values_not_refitted_or_physically_certified':True,
        }
    previous=report.get('C5_preview',{}).get('epsilon_C5_frobenius_uncertified')
    if previous is not None:
        require(math.isclose(mappings['identity_prior_preview']['relative_prior_Batch9_unmixed_projection'],
                             float(previous),rel_tol=2e-7,abs_tol=2e-9),
                'Batch9 identity unmixed C5 preview did not reproduce')
    additional=None
    if mixing_scenario is not None:
        mixed=endpoint(mixing_scenario)
        additional={
            'candidate_input_mixing_run_endpoint_mu2_over_min_diagonal_m2':float(abs(mixed['mu'])/min(mixed['m1'],mixed['m2'])),
            'mapping_results':{},
        }
        for name,sign in MAPS.items():
            c=evaluate(mixed,mixed['quartics'],sign=sign,unmixed=False)
            comparison=mappings[name]['C5_one_scale_candidate_gev_minus1']
            cprev=np.asarray(comparison['real'])+1j*np.asarray(comparison['imag'])
            additional['mapping_results'][name]={
                'C5_one_scale_candidate_gev_minus1':matrix_json(c),
                'relative_change_from_quartic_only_endpoint':relative(cprev,c),
            }
    return {
        'model':model, 'historical_fit_chi2':source['chi2'],
        'historical_fit_passes_5p86':source['fit_accepted'],
        'initial_scalar_boundary_scale_gev':1e7 if source['uv_on'] else 1e5,
        'zero_endpoint_mu2_over_min_m2':float(abs(zero['mu'])/min(zero['m1'],zero['m2'])),
        'general_endpoint_mu2_over_min_m2':float(abs(general['mu'])/min(general['m1'],general['m2'])),
        'general_endpoint_mu2_gev2':b7.complex_json(general['mu']),
        'actual_candidate_mixed_mass_matrix_used_in_one_scale_expression':True,
        'quartics_used_in_matchete_C5_only':list(b7.DIRECT),
        'mappings':mappings,
        'explicit_mixing_scenario':additional,
        'physical_mixed_C5_certified':False,
        'sequential_eft1_general_C5_performed':False,
        'yukawas_and_MF_frozen_historical_boundary':True,
    }


def check_output(case: str, run_id: str, *, root:Path = ROOT)->Path:
    require(case in b9.CASES, 'unknown historical case')
    tag=b9.parse_run_id(run_id)
    out=(Path(root)/'output/studies/general_scalar_extended'/case/tag)
    require(out.resolve().is_relative_to((Path(root)/'output/studies/general_scalar_extended').resolve()),
            'output path outside isolated general-scalar tree')
    require(not out.exists(), 'output run already exists (choose a NEW --run-id): '+str(out))
    return out


def pilot_path(case:str, source_run_id:str, model:str, *,root:Path=ROOT)->Path:
    b9.parse_run_id(source_run_id)
    return (Path(root)/'output/studies/general_scalar_extended'/case/source_run_id/
            f'{b9.safe_tag(model)}_paired_pilot.json')


def format_summary(case:str,run_id:str,rows:dict,*,mixed_run:bool)->str:
    lines=['# B/C candidate mass-mixed C5 and matching-convention diagnostic','',
           f'Historical study: `{case}`; new run: `{run_id}`','',
           '**Not a physical or sequential EFT1 C5 prediction.** Uses the saved Matchete one-scale ordered expression with candidate RGBeta scalar endpoint values and historical frozen Yukawas/MF.',
           'Two convention hypotheses: prior identity m=r, and Stage23/26 conditional signed map m=-r for the two direct quartic operators. Neither physically certified.','',
           '| Model | Historical fit passes | abs(mu²)/min(diagonal m²) | Identity joint ΔC5 / C5_0 | Conditional-minus joint ΔC5 / C5_0 | Conditional-minus mixing-only |',
           '|---|---|---:|---:|---:|---:|']
    for model,r in rows.items():
        identity=r['mappings']['identity_prior_preview']['effects']
        signed=r['mappings']['stage23_conditional_minus']['effects']
        lines.append(f"| {model} | {'yes' if r['historical_fit_passes_5p86'] else 'NO'} | "
                     f"{r['general_endpoint_mu2_over_min_m2']:.6g} | "
                     f"{identity['relative_total']:.6g} | {signed['relative_total']:.6g} | "
                     f"{signed['relative_mixing_at_zero_direct_quartics']:.6g} |")
    lines+=['','## Per-model decomposition (all norms divided by ||C5_ref||_F)','',
            '| Model | Sign | Background diagonal | Mixing at q=0 | Direct quartics, unmixed | Quartic×mixing interference | Total |',
            '|---|---|---:|---:|---:|---:|---:|']
    for m,r in rows.items():
        for name,h in r['mappings'].items():
            e=h['effects']
            lines.append(f"| {m} | {name} | {e['relative_background_diagonal']:.6g} | "
                         f"{e['relative_mixing_at_zero_direct_quartics']:.6g} | "
                         f"{e['relative_direct_quartics_unmixed']:.6g} | "
                         f"{e['relative_quartic_mass_mixing_interference']:.6g} | "
                         f"{e['relative_total']:.6g} |")
    if mixed_run:
        lines+=['','## Explicit quadratic-mixing scenario', '',
                'Additional input-mixing run is compared to the same signed/identity quartic-only endpoint; its RG-evolved quartics and masses may also differ.',
                '| Model | Map | abs(mu²)/min(m²) | Relative C5 shift vs quartics-only |',
                '|---|---|---:|---:|']
        for m,r in rows.items():
            s=r['explicit_mixing_scenario'];assert s is not None
            for mode,entry in s['mapping_results'].items():
                lines.append(f"| {m} | {mode} | {s['candidate_input_mixing_run_endpoint_mu2_over_min_diagonal_m2']:.6g} | {entry['relative_change_from_quartic_only_endpoint']:.6g} |")
    lines+=['','## Interpretation boundaries','',
            '- The reference C5 omits radiatively regenerated quartics in the zero-input branch, consistent with the Batch9 fixed-order preview. Their endpoint values remain saved upstream.',
            '- Background term includes diagonal scalar masses and lambdaT3 changes. Yukawas/MF remain at saved historical fitted boundary, not RG endpoints.',
            '- The mixed-mass U and eigenmasses are inserted in an **existing** one-scale Matchete coefficient, not a freshly matched general-potential EFT.',
            '- All four decomposition components add as complex matrices, not as their displayed Frobenius norms.',
            '- 1D signed maps are conditionally source-supported; a physical global Matchete/RGBeta field identification and nonzero mixed threshold matching remain unverified.',
            '- Additional five scalar quartics can affect RG evolution but are not directly inserted into this saved one-scale C5 expression.',
            '- Historical results and preceding Batch9 runs remain untouched.', '']
    return '\n'.join(lines)


def execute_case(case:str,models:tuple[str,...],*, source_run_id:str,
                 mixing_run_id:str|None,run_id:str,root:Path=ROOT)->dict:
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    dest=check_output(case,run_id,root=root)
    evidence23=b7.read_json(Path(root)/EVIDENCE23.relative_to(ROOT))
    evidence26=b7.read_json(Path(root)/EVIDENCE26.relative_to(ROOT))
    maps=check_convention_evidence(evidence23,evidence26,models=models)
    registry=load_registry(Path(root)/b7.REGISTRY.relative_to(ROOT))
    rows={}
    input_shas={}
    for model in models:
        source=b9.study_source(case,model,registry,root=root)
        path=pilot_path(case,source_run_id,model,root=root)
        report=b7.read_json(path)
        scenarios=check_pilot_input(report,case,model,source,require_preview=True)
        other_scene=None
        if mixing_run_id is not None:
            mpath=pilot_path(case,mixing_run_id,model,root=root)
            mrep=b7.read_json(mpath)
            check_pilot_input(mrep,case,model,source,require_preview=False)
            require('general_quartics_and_mixing' in mrep['scenarios'] and
                    mrep['scenarios']['general_quartics_and_mixing'].get('scenario') == 'general_quartics_and_mixing' and
                    mrep['scenarios']['general_quartics_and_mixing'].get('uv_integrated') == source['uv_on'] and
                    mrep.get('mass_mixing_input_fraction_of_min_m2') is not None and
                    0 < float(mrep['mass_mixing_input_fraction_of_min_m2']) <= 0.20 and
                    mrep.get('quartic_input_fraction_of_abs_lambdaT3') == report.get('quartic_input_fraction_of_abs_lambdaT3'),
                    'explicit mixing run missing, invalid, or uses different quartic input')
            other_scene=mrep['scenarios']['general_quartics_and_mixing']
            input_shas[f'{model}_explicit_mixing_pilot']=b7.digest(mpath)
        name={'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model]
        coefficient=DEFAULT_SOURCE/(name+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        evaluator=make_matchete_evaluator(model,source,coefficient)
        row=model_output(model,report,source,evaluator,mixing_scenario=other_scene)
        row['pilot_source_json_sha256']=b7.digest(path)
        row['one_scale_ordered_coefficient_sha256']=b7.digest(coefficient)
        row['historical_fitted_config_sha256']=b7.digest(source['cfg_path'])
        row['original_fit_record_sha256']=b7.digest(source['fit_path'])
        rows[model]=row
        input_shas[f'{model}_paired_pilot']=b7.digest(path)
        input_shas[f'{model}_one_scale_ordered_coefficient']=b7.digest(coefficient)
    result={'format':FORMAT,'case':case,'run_id':run_id,'models':list(models),
            'source_run_id':source_run_id,'mixing_run_id':mixing_run_id,
            'status':'ConditionalConventionsOnlyNotPhysicalC5',
            'source_stage23_sha256':b7.digest(Path(root)/EVIDENCE23.relative_to(ROOT)),
            'source_stage26_sha256':b7.digest(Path(root)/EVIDENCE26.relative_to(ROOT)),
            'checked_direct_maps':maps,'source_files_sha256':input_shas,
            'model_results':rows,'original_historical_results_modified':False,
            'upstream_Batch9_results_modified':False,
            'physical_rgbeta_matchete_map_certified':False,
            'physical_mixed_eft_C5_matched':False}
    # Write only after every requested model and consistency check succeeds.
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (dest/'summary.md').write_text(format_summary(case,run_id,rows,
                           mixed_run=mixing_run_id is not None),encoding='utf-8')
    return result


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='All three historical cases')
    parser.add_argument('--case',choices=tuple(b9.CASES))
    parser.add_argument('--model',choices=b9.MODELS)
    parser.add_argument('--source-run-id',default='pilot_q0p1_c5_unmixed',
                        help='Existing Batch9 quartic-only paired pilot run ID')
    parser.add_argument('--mixing-run-id',help='Optional existing Batch9 explicit-mixing run ID')
    parser.add_argument('--run-id',default='matching_conventions_v1',
                        help='New immutable run ID; must not exist')
    args=parser.parse_args()
    require(not(args.case and args.all),'Choose --all or --case, not both')
    b9.parse_run_id(args.source_run_id)
    if args.mixing_run_id:
        b9.parse_run_id(args.mixing_run_id)
    cases=tuple(b9.CASES) if args.all or not args.case else (args.case,)
    models=(args.model,) if args.model else b9.MODELS
    # Refuse an already-used run ID in ANY selected study before writing the first.
    for case in cases:
        check_output(case,args.run_id)
    for case in cases:
        result=execute_case(case,models,source_run_id=args.source_run_id,
                            mixing_run_id=args.mixing_run_id,run_id=args.run_id)
        for model,row in result['model_results'].items():
            ri=row['mappings']['identity_prior_preview']['effects']['relative_total']
            rs=row['mappings']['stage23_conditional_minus']['effects']['relative_total']
            print(f'{case}: {model} joint delta_C5/||C5_0|| identity={ri:.6g}, '
                  f'conditional_minus={rs:.6g}; physical status=UNCERTIFIED',flush=True)
        print('Saved:',OUTPUT/case/args.run_id/'summary.md',flush=True)

if __name__=='__main__':
    main()
