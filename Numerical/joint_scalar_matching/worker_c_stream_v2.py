"""Worker C, additive numerical scalar-threshold handoff and 16-model audit.

UV -> all fermions removed at muF -> EFT1(SM+S1+S2) -> all scalars at muS
-> SMEFT.  Source releases lacking signed spinor normalization stay CONDITIONAL.
No production solver changes; no independent hard scalar remainder is fabricated.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

MU_F = 100000.0
MU_S = 1000.0
MODELS = (
    'A_-4','A_-2','A_+0', 'B_-3','B_-1','B_+1',
    'C_-3','C_-1','C_+1', 'D_-2','D_+0','D_+2',
    'E_-4','E_-2','E_+0','E_+2',
)
FORMAT = 't3-worker-c-eft1-smeft-batch-v2'
RESTRICTED = 'Y1_Y2_only_no_independent_conjugate'


class HandoffError(ValueError):
    pass


def require(condition: Any, message: str) -> None:
    if not bool(condition):
        raise HandoffError(message)


def read_doc(path: Path) -> dict:
    value = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), f'Expected object: {path}')
    return value


def fingerprint(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def slug(model: str) -> str:
    require(model in MODELS, f'Unknown T3 model: {model}')
    return model.replace('-', 'm').replace('+', 'p')


def _scale(value: Any, expected: float, label: str) -> float:
    try:
        actual = float(value)
    except (TypeError, ValueError) as exc:
        raise HandoffError(f'{label} missing or nonnumeric') from exc
    require(math.isfinite(actual) and actual > 0 and
            math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10),
            f'{label} must equal {expected:g} GeV (got {actual:g})')
    return actual


def _encode(a: Any) -> dict:
    z = np.asarray(a, dtype=np.complex128)
    require(z.shape == (3,3) and np.isfinite(z).all(), 'Invalid 3x3 complex C5')
    return {'real': z.real.tolist(), 'imag': z.imag.tolist()}


def _decode(obj: Any) -> np.ndarray:
    if isinstance(obj, dict) and {'real','imag'} <= obj.keys():
        z = np.asarray(obj['real'],float) + 1j * np.asarray(obj['imag'],float)
    else:
        z = np.asarray(obj,complex)
    require(z.shape == (3,3) and np.isfinite(z).all(), 'Bad 3x3 complex matrix')
    return z


def _norm(a: Any) -> float:
    return float(np.linalg.norm(_decode(a)))


def _symmetry(a: Any) -> float:
    z = _decode(a)
    return float(np.linalg.norm(z - z.T) / max(float(np.linalg.norm(z)), 1e-100))


def _json_dump(path: Path, value: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive create: never overwrite historical or same-run results.
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def _input_npz(meta: Path) -> Path:
    doc = read_doc(meta)
    name = doc.get('tensor_npz')
    require(isinstance(name,str) and name not in ('','.', '..') and
            Path(name).name == name, 'Tensor NPZ must be adjacent basename')
    path = Path(meta).parent/name
    require(path.is_file(), f'Missing scalar NPZ: {path}')
    return path


def validate_endpoint_contract(meta_path: Path) -> dict:
    """Audit an *actual muS evolved* Wilson+background handoff, not an old toy.

    The canonical array shapes and convention are validated by the established
    adapter when run_endpoint is called; this function validates ownership and
    bookkeeping flags not present in that general adapter.
    """
    path = Path(meta_path)
    doc = read_doc(path)
    from Numerical.joint_scalar_matching.adapter import CONVENTION, normalize_model
    model = normalize_model(doc.get('model',''))
    _scale(doc.get('mu_f_gev'), MU_F, 'mu_f_gev')
    _scale(doc.get('mu_s_gev'), MU_S, 'mu_s_gev')
    require(doc.get('tensor_convention') == CONVENTION, 'Unknown tensor convention')
    require(doc.get('renormalization_scheme') == 'MSbar', 'Requires MSbar')
    require(doc.get('synthetic_fixture') is not True, 'Synthetic inputs excluded from endpoint')
    require(doc.get('scalar_mass_origin') in ('unbroken_electroweak','EWSB_resummed'),
            'Endpoint scalar mass origin missing')
    require(doc.get('yukawa_sector') == RESTRICTED,
            'Restricted Yukawa source must be explicitly documented; no y2Conj assumption')
    for key in ('all_fermions_integrated_at_mu_f','all_scalars_joint_at_mu_s',
                'C_LLSS_evolved_to_mu_s','C5_EFT1_evolved_to_mu_s',
                'no_double_counted_eft1_log','canonical_real_matching_basis_verified'):
        require(doc.get(key) is True, f'Endpoint declaration missing: {key}')
    if model in ('B_-1','C_-1'):
        require(doc.get('mixed_scalar_mass_basis_verified') is True,
                'B/C mixed-scalar source-to-eigenbasis mapping not verified')
    require(isinstance(doc.get('provenance'),dict) and
            isinstance(doc['provenance'].get('source'),str) and doc['provenance']['source'],
            'Endpoint source provenance required')
    npz = _input_npz(path)
    flags = {
        'spinor_normalization_verified': doc.get('spinor_normalization_verified') is True,
        'scalar_quartic_and_mass_provenance_verified': doc.get('scalar_quartic_and_mass_provenance_verified') is True,
        'full_scalar_one_loop_operator_inventory_verified': doc.get('full_scalar_one_loop_operator_inventory_verified') is True,
    }
    # Even all three being true does not imply we computed omitted nonbubble terms.
    return {'model':model, 'meta':str(path), 'npz':str(npz),
            'meta_sha256':fingerprint(path),'npz_sha256':fingerprint(npz),
            'verification_flags':flags,
            'all_physical_checks_passed':all(flags.values()),
            'safe_to_claim_full_one_loop':False}


def run_endpoint(meta_path: Path, *, tolerance: float=1e-8) -> dict:
    """Use the actual evolved Wilson tensors, add a *single* scalar bubble.

    Data contract: adapter.load_input v1 JSON/NPZ plus explicit endpoint fields
    checked above.  Does not manufacture C5(MS) from a muF tree boundary.
    """
    audit = validate_endpoint_contract(meta_path)
    from Numerical.joint_scalar_matching.adapter import load_input, scalar_groups
    from Numerical.joint_scalar_matching.worker_c_joint import match_joint_scalars
    inp = load_input(Path(meta_path))
    require(inp.model == audit['model'], 'Adapter model changed during handoff')
    wilsons = {
        'C5_EFT1':inp.c5_eft1, 'C_LLSS':inp.c_llss,
        'mu_s_gev':inp.mu_s_gev, 'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
        'renormalization_scheme':'MSbar',
        'provenance':{'source':inp.provenance['source'],
                      'endpoint':str(meta_path), 'sha256':audit['meta_sha256']},
    }
    background = {
        'mass_squared_gev2':inp.mass_squared,
        'lambda_HHSS':inp.lambda_hhss,
        'scalar_groups':scalar_groups(inp.model),
        'scalar_mass_origin':inp.scalar_mass_origin,
        'mu_s_gev':inp.mu_s_gev,
        'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
        'renormalization_scheme':'MSbar',
        'provenance':{'source':inp.provenance['source'],
                      'npz_sha256':audit['npz_sha256']},
    }
    if inp.quartic_channels is not None:
        background['lambda_HHSS_channels'] = inp.quartic_channels
        background['quartic_channel_names'] = list(inp.quartic_channel_names)
    matched = match_joint_scalars(wilsons,background,MU_S)
    inherited = np.asarray(matched['c5_inherited'],complex)
    bubble = np.asarray(matched['c5_finite_scalar'],complex)
    total = np.asarray(matched['c5_total'],complex)
    require(np.allclose(total,inherited+bubble,rtol=1e-11,atol=1e-28),
            'Inherited coefficient/bubble added incorrectly')
    residual = float(matched['weinberg_projection_residual'])
    require(residual <= tolerance, f'Weinberg Ward projection failed: {residual:g}')
    pair = matched['contribution_breakdown']['eigenstate_pairs']
    channels = matched['contribution_breakdown']['quartic_channels']
    check_pairs = sum((np.asarray(z,complex) for z in pair.values()),np.zeros((3,3),complex))
    require(np.allclose(bubble,check_pairs,rtol=1e-8,atol=1e-27),'Eigenpair sum mismatch')
    return {
        'format':FORMAT,'status':'EvolvedEFT1_OneJointScalarBubble_Evaluated',
        'model':inp.model,'mu_f_gev':MU_F,'mu_s_gev':MU_S,
        'c5_inherited_gev_inv':_encode(inherited),
        'c5_finite_scalar_LLSS_HHSS_gev_inv':_encode(bubble),
        'c5_total_one_bubble_only_gev_inv':_encode(total),
        'c5_norms_gev_inv':{'inherited':float(np.linalg.norm(inherited)),
                            'bubble':float(np.linalg.norm(bubble)),
                            'total':float(np.linalg.norm(total))},
        'eigenpair_contributions_gev_inv':{k:_encode(v) for k,v in pair.items()},
        'quartic_channel_contributions_gev_inv':{k:_encode(v['c5_finite_scalar']) for k,v in channels.items()},
        'scalar_eigenvalues_sq_gev2':np.asarray(matched['scalar_eigenvalues_sq_gev2'],float).tolist(),
        'weinberg_projection_residual':residual,
        'source_audit':audit,'restricted_Y1_Y2_declared':True,
        'nonbubble_finite_hard_contribution_included':False,
        'complete_scalar_matching_certified':False,
        'safe_to_add_as_second_scalar_bubble':False,
        'remarks':'Evolved EFT1 inherited once, joint LLSS-HHSS finite scalar bubble once; no additional hard scalar diagrams invented.',
    }


def compare_archived_conditional(current: dict, archive_path: Path) -> dict:
    """Regression comparison ONLY for identical source/background SHA inputs.

    An old D(-2) 3 TeV benchmark or a different scalar basis is not comparable.
    Any source mismatch is disclosed, never silently rescaled or subtracted.
    """
    path = Path(archive_path)
    if not path.is_file():
        return {'status':'NoArchivedReference'}
    archived=read_doc(path)
    old=archived.get('response')
    now=current['conditional_response_per_unit_unknown_normalization']
    require(isinstance(old,dict) and archived.get('model')==current['model'],
            'Malformed archived conditional response')
    keys=('model_source_boundary_sha256','background_json_sha256',
          'background_npz_sha256','mu_f_gev','mu_s_gev')
    mismatch=[key for key in keys if old.get(key)!=now.get(key)]
    if mismatch:
        return {'status':'IncomparableDifferentInputs','mismatched_fields':mismatch,
                'archive_sha256':fingerprint(path)}
    matrices=('C5_joint_scalar_finite_bubble_per_N_gev_inverse',
              'C5_EFT1_induced_leading_log_per_N_gev_inverse',
              'C5_LLSS_total_one_loop_per_N_gev_inverse')
    relative={}
    for key in matrices:
        lhs=_decode(now[key]);rhs=_decode(old[key])
        relative[key]=float(np.linalg.norm(lhs-rhs)/max(np.linalg.norm(lhs),np.linalg.norm(rhs),1e-100))
    largest=max(relative.values())
    require(largest<2e-9,f'Archived same-input numerical regression: {largest:g}')
    return {'status':'IdenticalInputsNumericalRegressionPassed',
            'archive_sha256':fingerprint(path),
            'max_relative_matrix_difference':largest,
            'individual_relative_differences':relative}


def survey_one(source_path: Path, background_path: Path, expected_model: str) -> dict:
    """Recompute all SU(2) and eigenpairs using existing project kernel."""
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.large_batch import load_background,evaluate_conditional
    require(expected_model in MODELS,'Unknown model')
    source=load_source(source_path,expected_model=expected_model)
    background=load_background(background_path,source)
    _scale(source.mu_f_gev,MU_F,'source.mu_f_gev')
    _scale(background.mu_f,MU_F,'background.mu_f_gev')
    _scale(background.mu_s,MU_S,'background.mu_s_gev')
    _scale(background.reference_scale,MU_S,'background.reference_scale_gev')
    data=evaluate_conditional(source,background)
    require(data.get('model')==expected_model, 'Kernel model mismatch')
    bubble=_decode(data['C5_joint_scalar_finite_bubble_per_N_gev_inverse'])
    running=_decode(data['C5_EFT1_induced_leading_log_per_N_gev_inverse'])
    net=_decode(data['C5_LLSS_total_one_loop_per_N_gev_inverse'])
    require(np.allclose(net,bubble+running,rtol=1e-9,atol=1e-27),
            'Conditional LLSS log + bubble do not close')
    require(_symmetry(bubble) < 1e-8 and _symmetry(net)<1e-8,
            'Weinberg flavour symmetry failed')
    r=float(data['relative_matching_scale_cancellation_residual'])
    require(r<=2e-9, f'Scalar matching-scale cancellation failed: {r:g}')
    npz=_input_npz(background_path)
    return {
        'model':expected_model,
        'status':'ConditionalSourceNormalisation_NotPhysicalMatching',
        'mu_f_gev':MU_F,'mu_s_gev':MU_S,
        'source_boundary_sha256':fingerprint(source_path),
        'scalar_background_json_sha256':fingerprint(background_path),
        'scalar_background_npz_sha256':fingerprint(npz),
        'source_signed_GammaCC_Weyl_normalization_certified':False,
        'restricted_Y1_Y2_source_certified':False,
        'C5_from_evolved_EFT1_Wilsons':False,
        'conditional_response_per_unit_unknown_normalization':data,
        'diagnostics':{
            'bubble_norm_gev_inv_per_unit_N':float(np.linalg.norm(bubble)),
            'leading_log_norm_gev_inv_per_unit_N':float(np.linalg.norm(running)),
            'total_norm_gev_inv_per_unit_N':float(np.linalg.norm(net)),
            'matching_scale_cancellation_relative':r,
            'flavour_symmetry_relative':_symmetry(net),
        },
        'production_C5_not_modified':True,
        'full_one_loop_scalar_matching_certified':False,
    }


def survey_all(root: Path, source_dir: Path, background_dir: Path, *, models=MODELS,
               output_dir: Path, keep_going: bool=True, archive_dir: Path | None=None) -> dict:
    """Fresh auditable all-16 reference calculation; default is fail-closed.

    No entry is called physical without independently evolved Wilson inputs.
    Failed models are retained with reason and no placeholder matrices.
    """
    root=Path(root).resolve()
    source_dir=Path(source_dir)
    background_dir=Path(background_dir)
    output_dir=Path(output_dir)
    require(not output_dir.exists(),f'Output already exists: {output_dir}')
    ordered=tuple(models)
    require(len(ordered)==len(set(ordered)) and all(m in MODELS for m in ordered),
            'Invalid requested model subset')
    output_dir.mkdir(parents=True,exist_ok=False)
    results=[]
    for model in ordered:
        source=source_dir/slug(model)/'boundary.json'
        background=background_dir/(slug(model)+'.json')
        try:
            require(source.is_file(),f'Missing fermion release: {source}')
            require(background.is_file(),f'Missing scalar background: {background}')
            result=survey_one(source,background,model)
            if archive_dir is not None:
                result['archived_comparison']=compare_archived_conditional(
                    result,Path(archive_dir)/(slug(model)+'.json'))
            status='PASS_CONDITIONAL'
        except (ValueError,OSError,KeyError,TypeError,ImportError) as exc:
            result={'model':model,'status':'BLOCKED','reason':str(exc),
                    'physical_coefficients_fabricated':False}
            status='BLOCKED'
        _json_dump(output_dir/(slug(model)+'.json'), result)
        results.append({'model':model,'status':status,
                        'report':slug(model)+'.json',
                        **({'reason':result['reason']} if status=='BLOCKED' else
                           {'bubble_norm_per_unit_N':result['diagnostics']['bubble_norm_gev_inv_per_unit_N']})})
        if status=='BLOCKED' and not keep_going:
            break
    passed=sum(row['status']=='PASS_CONDITIONAL' for row in results)
    report={'format':FORMAT,'mode':'source_resolved_16_model_conditional_survey',
            'status':'ConditionalSurveyCompleted' if passed==len(ordered) else 'Incomplete_BlockedModels',
            'required_mu_f_gev':MU_F,'required_mu_s_gev':MU_S,
            'requested_model_count':len(ordered),'models_processed':len(results),
            'conditionally_evaluated':passed,
            'blocked':len(results)-passed,
            'no_conjugate_yukawa_assumption_made':True,
            'absolute_physical_C5_available':False,
            'evolved_EFT1_C5_available_from_these_inputs':False,
            'full_scalar_one_loop_matching_completed':False,
            'source_directory':str(source_dir),'background_directory':str(background_dir),
            'archive_directory':str(archive_dir) if archive_dir else None,
            'rows':results,
            'remaining_for_full_numerical_pipeline':remaining_work(),
            'historical_results_modified':False}
    _json_dump(output_dir/'summary.json',report)
    return report


def remaining_work() -> list[dict]:
    """Conservative cross-worker status, not claims about unseen local code."""
    return [
        {'owner':'Worker A','item':'Finite renormalized UV -> EFT1 hard fermion matching in restricted Y1/Y2',
         'status':'NOT_ESTABLISHED_IN_THIS_PACKAGE'},
        {'owner':'Worker B','item':'All-16 canonical running scalar mass/quartic tensors including special bilinears',
         'status':'REQUIRES_LOCAL_WORKER_B_HANDOFF_AND_PROVENANCE'},
        {'owner':'Worker D','item':'Actual evolved C_LLSS(muS) and inherited C5_EFT1(muS), including 7-term running, exported in canonical endpoint NPZ',
         'status':'REQUIRES_LOCAL_EFT1_ENDPOINT_EXPORT'},
        {'owner':'Worker A/D','item':'Uniform explicit Y1,Y2-only source in all 16 models with y2Conj absent',
         'status':'REQUIRES_SOURCE_REGENERATION_OR_VERIFIED_RESTRICTION'},
        {'owner':'Worker C','item':'Additional independent finite scalar topologies/counterterms and UV-minus-EFT1 soft subtraction',
         'status':'NOT_COMPUTED; DO_NOT_ADD_UNSUBTRACTED_MATchete_FINITE_PART'},
        {'owner':'Workers A/B/C','item':'Signed GammaCC/Weyl operator and CG normalization alignment',
         'status':'CONDITIONAL_IN_EXISTING_REFERENCE_SOURCE'},
        {'owner':'Worker D','item':'Run 16 model numerical benchmarks and all 256 stage freezes with new certified handoffs',
         'status':'AFTER_ENDPOINTS_AVAILABLE; PRODUCTION_SOLVER_UNTOUCHED'},
    ]


def endpoint_contract_template(model: str) -> dict:
    """Explicit metadata template only: NO fictitious arrays or certifications."""
    require(model in MODELS,'Unknown model')
    return {'schema_version':1,'model':model,
            'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
            'renormalization_scheme':'MSbar',
            'mu_f_gev':MU_F,'mu_s_gev':MU_S,'reference_scale_gev':MU_S,
            'scalar_mass_origin':None, 'tensor_npz':slug(model)+'.npz',
            'synthetic_fixture':False,
            'yukawa_sector':RESTRICTED,
            'all_fermions_integrated_at_mu_f':False,
            'all_scalars_joint_at_mu_s':False,
            'C_LLSS_evolved_to_mu_s':False,
            'C5_EFT1_evolved_to_mu_s':False,
            'no_double_counted_eft1_log':False,
            'canonical_real_matching_basis_verified':False,
            'mixed_scalar_mass_basis_verified':False,
            'spinor_normalization_verified':False,
            'scalar_quartic_and_mass_provenance_verified':False,
            'full_scalar_one_loop_operator_inventory_verified':False,
            'provenance':{'source':'REPLACE_WITH_ACTUAL_RUNNING_EFT1_ENDPOINT'},
            'quartic_channel_names':[],
            'REQUIRED_NPZ_ARRAYS':{
                'mass_squared_gev2':'(n,n) full real scalar mass Hessian at muS',
                'C_LLSS':'(3,3,2,2,n,n) evolved to muS in canonical matching basis',
                'lambda_HHSS':'(4,4,n,n) at muS, complete source HHSS',
                'C5_EFT1':'(3,3) complex symmetric, already RG evolved to muS',
                'lambda_HHSS_channels':'optional (k,4,4,n,n), sum equals lambda_HHSS'},
            'WARNING':'INCOMPLETE CONTRACT. Do not fill truth flags without checks.'}
