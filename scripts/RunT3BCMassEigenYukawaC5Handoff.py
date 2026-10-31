"""Read-only B_-1/C_-1 mass-eigenstate Yukawa -> conditional one-scale C5 handoff.

Uses *existing* successful Matchete heavy-mass-basis/Yukawa certificates,
registered 110101 CG evidence, a saved Matchete one-loop ordered C5 expression,
and the saved conditional scalar EFT1 endpoint. It does not perform matching.

Run from repository root:
  python -m scripts.RunT3BCMassEigenYukawaC5Handoff --all --run-id bc_yukawa_c5_v1
  python -m scripts.RunT3BCMassEigenYukawaC5Handoff --all --run-id bc_yukawa_c5_replay_v1 --recompute-c5

The optional replay calls the existing one-scale numerical evaluator on the
saved Matchete expression; it does not compute a new one-loop diagram.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

ROOT = Path(__file__).resolve().parents[1]
MODELS = ('B_-1', 'C_-1')
OUT = ROOT / 'output/audits/t3_bc_mass_eigen_yukawa_c5'
CG_SUMMARY = ROOT / ('output/audits/t3_bc_matchete_registered_cg/'
                     'bc_live_registered_cg_v1/summary.json')
FORMAT = 't3-bc-native-yukawa-mass-eigen-one-scale-handoff-v1'
STATUS = 'ExistingNativeYukawaMassBasisAndConditionalOneScaleC5Checked'
EXPECTED_Y2_CG = {'B_-1': -1, 'C_-1': 1}


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError('B/C mass-eigen Yukawa handoff: ' + reason)


def read_json(path: Path) -> dict:
    require(path.is_file(), 'missing file: ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'not a JSON object: ' + str(path))
    return data


def digest(path: Path) -> str:
    require(path.is_file(), 'missing hash input: ' + str(path))
    return sha256(path.read_bytes()).hexdigest()


def pair(raw: Any, field: str) -> complex:
    require(isinstance(raw, dict) and set(raw) == {'real', 'imag'}, field + ': missing complex pair')
    try:
        z = complex(float(raw['real']), float(raw['imag']))
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(field + ': invalid complex number') from exc
    require(np.isfinite([z.real, z.imag]).all(), field + ': nonfinite complex number')
    return z


def complex_matrix(raw: Any, name: str, *, shape: tuple[int, ...] | None = None) -> np.ndarray:
    require(isinstance(raw, list) and bool(raw) and all(isinstance(row, list) for row in raw),
            name + ': malformed complex matrix')
    require(len(set(map(len, raw))) == 1 and len(raw[0]) > 0, name + ': ragged matrix')
    arr = np.asarray([[pair(x, name) for x in row] for row in raw], dtype=np.complex128)
    require(arr.ndim == 2 and (shape is None or arr.shape == shape), name + ': unexpected shape')
    return arr


def split_matrix(raw: Any, name: str, n: int) -> np.ndarray:
    require(isinstance(raw, dict) and set(raw) == {'real', 'imag'}, name + ': no real/imag matrix')
    try:
        x = np.asarray(raw['real'], dtype=float)
        y = np.asarray(raw['imag'], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(name + ': nonnumeric matrix') from exc
    require(x.shape == (n, n) and y.shape == (n, n) and
            np.isfinite(x).all() and np.isfinite(y).all(), name + ': invalid matrix')
    return x + 1j * y


def pairing_matrix(a: np.ndarray) -> list[list[dict[str, float]]]:
    return [[{'real':float(complex(v).real),'imag':float(complex(v).imag)} for v in row]
            for row in np.asarray(a, dtype=np.complex128)]


def yukawa_weights(model: str, U: np.ndarray, matchete: dict) -> dict:
    """Exact source native-CG factors; y2=-alpha2*U[1] (epsilon sign).

    Does NOT merge y1/y2 spinor structures: y1 uses CConj[F], y2 uses F.
    """
    require(model in MODELS and U.shape == (2, 2), 'unsupported model or U')
    require(np.allclose(U.conj().T @ U, np.eye(2), atol=1e-12, rtol=0), 'U not unitary')
    y1 = matchete.get('native_yukawa1_tensor')
    y2 = matchete.get('native_yukawa2_epsilon_tensor')
    require(isinstance(y1, dict) and isinstance(y2, dict), 'missing native Yukawa tensors')
    require(y1.get('ExactTensorProjectionVerified') is True and
            y2.get('ExactTensorProjectionVerified') is True, 'native tensor projection unverified')
    require(y1.get('OldInvariantCount') == y1.get('NewInvariantCount') == 1 and
            y2.get('OldInvariantCount') == y2.get('NewInvariantCount') == 1,
            'incorrect native Yukawa invariant count')
    require(y1.get('NativeCG') == 'T3MENativeY1CG' and
            y2.get('NativeCG') == 'T3MENativeY2CG', 'unknown native Yukawa CG')
    require(y1.get('PhysicalCouplingMapCertified') is False and
            y2.get('PhysicalCouplingMapCertified') is False,
            'unexpected upstream coupling certification')
    dF = 1 if model == 'B_-1' else 3
    require(y1.get('FermionSU2Dimension') == dF and
            y2.get('FermionSU2Dimension') == dF, 'Yukawa fermion dimension changed')
    # Source-derived *verified for these two models only*. These are exact
    # native tensor factors, not freely selected overall physical phases.
    require(y1.get('OldToNewCGFactor') == '1', 'native y1 CG factor changed')
    require(y2.get('OldTimesEpsilonCGFactor') == str(EXPECTED_Y2_CG[model]),
            'native y2 epsilon CG factor changed')
    w1 = U[0].copy()
    w2 = -EXPECTED_Y2_CG[model] * U[1].copy()
    require(abs(np.vdot(w1, w2)) < 3e-12, 'rotated source scalar row orthogonality failed')
    require(np.isclose(np.vdot(w1, w1), 1, atol=1e-12) and
            np.isclose(np.vdot(w2, w2), 1, atol=1e-12),
            'source Yukawa coefficient norm not preserved')
    return {'y1':w1, 'y2':w2, 'y2_epsilon_CG_factor':EXPECTED_Y2_CG[model],
            'y2_source_field_substitution_minus_sign':-1}


def check_live_cg(record: dict) -> None:
    require(record.get('format') == 't3-bc-live-registered-cg-110101-cross-tool-identity-v1' and
            record.get('status') == 'LiveMatcheteRegisteredCGExactWithGatedRGBetaIdentityNoEFTMatching',
            'live registered CG report incompatible')
    require(record.get('live_registered_matchete_components_compared') == 64 and
            record.get('signed_direct_and_conjugate_coefficients_compared') == 128 and
            record.get('live_matchete_registered_cg_component_identity_verified') is True and
            record.get('cross_tool_signed_component_identity_verified_in_exported_su2_basis') is True,
            'live CG / cross-tool coefficient gate not passed')
    for key in ('physical_coupling_mapping_certified','mass_basis_matching_certified',
                'sequential_eft_matching_certified','direct_live_matchete_action_subtraction_performed'):
        require(record.get(key) is False, 'live CG overclaims certification: '+key)
    rows = record.get('models')
    require(isinstance(rows, list) and len(rows) == 2 and
            {x.get('model') for x in rows} == set(MODELS), 'missing live CG models')
    for row in rows:
        require(row.get('live_registered_components_equal_native_stage22_exact') == 32 and
                row.get('live_matchete_registered_cg_direct_and_hc_coefficients_checked') == 64 and
                row.get('exact_registered_component_disagreements') == 0,
                'registered Matchete CG model check failed')


def check_one_scale(endpoint: dict, model: str) -> dict:
    require(endpoint.get('status') == 'ConditionalEFT1ScalarEndpointOneScaleC5NotPhysicalSequential' and
            endpoint.get('model') == model and endpoint.get('endpoint') == 'final',
            'wrong candidate one-scale record')
    for key in ('physical_matchete_rgbeta_mapping_certified',
                'scalar_threshold_matching_performed',
                'sequential_EFT1_wilson_running_performed',
                'yukawas_evolved_with_candidate_scalar_background',
                'physical_neutrino_mass_prediction'):
        require(endpoint.get(key) is False, 'one-scale status unexpectedly certified: '+key)
    candidate = complex_matrix(endpoint.get('candidate_one_scale_mnu_ev'), 'candidate mnu', shape=(3,3))
    baseline = complex_matrix(endpoint.get('baseline_one_scale_mnu_ev'), 'baseline mnu', shape=(3,3))
    for name,m in [('candidate',candidate),('baseline',baseline)]:
        require(np.allclose(m,m.T,atol=1e-14,rtol=1e-12), name + ' mnu not symmetric')
    eigs = sorted(np.linalg.svd(candidate,compute_uv=False).tolist(),reverse=True)
    old = sorted(np.linalg.svd(baseline,compute_uv=False).tolist(),reverse=True)
    for label,got in [('candidate',eigs),('baseline',old)]:
        src = endpoint.get(label+'_takagi_masses_ev_descending')
        require(isinstance(src,list) and len(src)==3 and
                np.allclose(got,src,rtol=2e-10,atol=1e-12),
                label + ': Takagi mass spectrum mismatch')
    diff = float(np.linalg.norm(candidate-baseline))
    relative = diff / float(np.linalg.norm(baseline))
    require(np.isfinite(relative) and
            np.isclose(diff,endpoint.get('delta_mnu_frobenius_ev'),rtol=2e-9,atol=1e-15) and
            np.isclose(relative,endpoint.get('relative_delta_mnu_frobenius'),rtol=2e-9,atol=1e-12),
            'C5 response norm mismatch')
    return {'candidate':candidate,'baseline':baseline,'masses':eigs,
            'relative':relative,'baseline_masses':old}


def paths(root: Path, model: str) -> dict[str, Path]:
    require(model in MODELS, 'unsupported model')
    tag = model.replace('-','m')
    base = root/'output/general_t3_mass_eigen_matchete'
    from scripts.RunGeneralT3ReferenceC5 import REFERENCE_ID
    return {
        'matchete':base/f'{tag}_final_matchete.json',
        'parameters':base/f'{tag}_final_parameters.json',
        'one_loop':base/'one_loop_probe'/f'{tag}_final_one_loop_probe.json',
        'one_loop_expression':base/'one_loop_probe'/f'{tag}_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt',
        'endpoint':root/'output/studies/general_scalar_extended/quadratic_endpoint_c5'/f'{tag}_final_endpoint_one_scale.json',
        'frozen':root/'Numerical/benchmarks/fitted_reference'/f'{REFERENCE_ID[model][0]}.json',
        'frozen_manifest':root/'Numerical/benchmarks/fitted_reference/benchmark_manifest.json',
        'mass_basis':root/'output/general_t3_mass_basis'/f'{tag}_final_mass_basis.json',
    }


def prepare_one(model: str, root: Path, cg: dict, *, recompute: bool=False) -> dict:
    p = paths(root,model)
    mat = read_json(p['matchete'])
    params = read_json(p['parameters'])
    one_loop = read_json(p['one_loop'])
    ep = read_json(p['endpoint'])
    require(mat.get('status') == 'MatcheteMassEigenHeavyBasisPassedPhysicalMatchingNotPerformed' and
            mat.get('model') == model and
            mat.get('mass_basis_ordinary_check_passed') is True and
            mat.get('heavy_mass_basis_detailed_check_passed') is True and
            mat.get('both_yukawa_interactions_present') is True and
            mat.get('all_old_scalar_fields_removed') is True and
            mat.get('all_17_quartics_present') is True and
            mat.get('physical_matching_performed') is False,
            'native heavy-mass-basis Matchete gate failed')
    require(params.get('format') == 't3-experimental-mass-eigen-matchete-parameters-v1' and
            params.get('model') == model and params.get('endpoint') == 'final' and
            params.get('quartic_source_count') == 17 and params.get('yukawa_source_count') == 2,
            'prepared native mass-eigen inputs incompatible')
    require(one_loop.get('model') == model and
            one_loop.get('status') == 'OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending' and
            one_loop.get('one_loop_match_stage_completed') is True and
            one_loop.get('physical_matching_performed') is False and
            one_loop.get('production_sequential_matching_performed') is False,
            'one-loop Matchete expression not validated')
    require(digest(p['one_loop_expression']) == ep.get('one_scale_expression_sha256'),
            'saved Matchete one-loop expression differs from candidate numerical run')
    mass_basis = p['mass_basis']
    require(digest(mass_basis) == params.get('source_sha256'),
            'prepared Matchete mass-basis source changed')
    U = split_matrix(params.get('unitary_U'), 'mass eigen U', 2)
    A = split_matrix(params.get('original_mass_squared_A_gev2'), 'mass matrix A', 2)
    masses = np.asarray(params.get('heavy_masses_gev'), dtype=float)
    require(masses.shape == (2,) and np.isfinite(masses).all() and np.all(masses > 0) and
            masses[0] < masses[1], 'nonpositive or unordered scalar eigenmasses')
    basis = MixedDoubletMassBasis.build(model,float(A[0,0].real),float(A[1,1].real),complex(A[0,1]))
    require(np.allclose(U.conj().T@A@U,np.diag(masses**2),atol=1e-5,rtol=1e-10) and
            np.allclose(masses**2,basis.eigenvalues_sq,atol=1e-5,rtol=1e-10),
            'mass diagonalization disagrees with source eigenmasses')
    for i in (1,2):
        for j in (1,2):
            z = pair(mat.get('rotation_spurion_values_from_numerical_U',{}).get(f'T3MEU{i}{j}'),
                     f'T3MEU{i}{j}')
            require(abs(z-U[i-1,j-1]) < 1e-12, 'saved Matchete spurion U differs')
    weights = yukawa_weights(model,U,mat)
    require(np.allclose(masses, ep.get('scalar_masses_gev'),rtol=1e-10,atol=1e-6),
            'one-scale endpoint uses different scalar masses')
    result = check_one_scale(ep,model)
    require(digest(p['frozen']) == ep.get('frozen_yukawa_source_sha256'),
            'frozen Yukawa data changed since numerical C5')
    from scripts.RunGeneralT3ReferenceC5 import make_transfer
    config = make_transfer(model,read_json(p['frozen_manifest']),p['frozen'],p['parameters'])
    yy1 = complex_matrix(config['y1'], 'frozen y1',shape=None)
    yy2 = complex_matrix(config['y2'], 'frozen y2',shape=None)
    require(yy1.shape==yy2.shape and yy1.shape[0]==3,
            'frozen Yukawa matrices have incompatible flavour indices')
    # y1 and y2 carry DIFFERENT heavy-fermion conjugations; preserve distinct
    # vertex matrices. Never add the two matrix structures as one vertex.
    vertices = {
        'X':{'y1_CConjF':pairing_matrix(weights['y1'][0]*yy1),
             'y2_F':pairing_matrix(weights['y2'][0]*yy2)},
        'Y':{'y1_CConjF':pairing_matrix(weights['y1'][1]*yy1),
             'y2_F':pairing_matrix(weights['y2'][1]*yy2)},
    }
    recomputed = False
    error = 0.0
    if recompute:
        from scripts.RunT3QuadraticEndpointC5 import numerical_c5, QUARTIC_MAP
        source = p['one_loop_expression'].read_text(encoding='utf-8')
        mapped = {name:pair(val,name) for name,val in ep['matchete_quartics_used'].items()}
        require(set(mapped) == {'lambdaT3',*QUARTIC_MAP},
                'numerical Matchete quartic input set changed')
        zeros = dict(mapped)
        for name in QUARTIC_MAP: zeros[name] = 0j
        zero_c5,zero_mnu,zero_masses = numerical_c5(source,config,basis,zeros)
        cand_c5,cand_mnu,cand_masses = numerical_c5(source,config,basis,mapped)
        error = float(max(np.max(np.abs(cand_mnu-result['candidate'])),
                          np.max(np.abs(zero_mnu-result['baseline']))))
        scale = max(1.,float(np.max(abs(result['candidate']))),float(np.max(abs(result['baseline']))))
        require(error <= scale*5e-8 and
                np.allclose(cand_masses,result['masses'],rtol=2e-8,atol=1e-10) and
                np.allclose(zero_masses,result['baseline_masses'],rtol=2e-8,atol=1e-10),
                'Matchete one-loop numerical replay differs from saved C5')
        recomputed = True
    # This reconstructs C5 solely from the saved mnu and vev. When --recompute
    # is enabled, that saved matrix is independently re-evaluated above.
    vev = float(config['vev_gev'])
    require(np.isfinite(vev) and vev > 0,'invalid reference electroweak vev')
    c5_from_saved = -2 * result['candidate']/(vev**2*1e9)
    return {
        'format':FORMAT, 'status':STATUS, 'model':model,
        'source':'existing native Matchete mass eigenstate Lagrangian and saved one-scale expression',
        'heavy_scalar_masses_gev':masses.tolist(),
        'scalar_U':{'real':U.real.tolist(), 'imag':U.imag.tolist()},
        'yukawa_scalar_weights':{k:[{'real':float(v.real),'imag':float(v.imag)} for v in weights[k]]
                                 for k in ('y1','y2')},
        'native_CG_factor_y1':1,
        'native_epsilon_projection_factor_y2':weights['y2_epsilon_CG_factor'],
        'y2_field_substitution_sign':weights['y2_source_field_substitution_minus_sign'],
        'fermion_spinor_structures_distinct':True,
        'frozen_flavour_yukawa_matrices_mass_eigenbasis':vertices,
        'source_U_unitarity_and_scalar_mass_diagonalization_verified':True,
        'native_yukawa_CG_factors_checked_against_live_Matchete_certificate':True,
        'original_two_yukawa_scalar_rows_norms_preserved':True,
        'one_scale_candidate_mnu_ev':pairing_matrix(result['candidate']),
        'one_scale_candidate_C5_gev_inverse_from_saved_mnu':pairing_matrix(c5_from_saved),
        'one_scale_candidate_takagi_masses_ev_descending':result['masses'],
        'one_scale_baseline_takagi_masses_ev_descending':result['baseline_masses'],
        'relative_extra_quartic_mnu_response':result['relative'],
        'Matchete_one_loop_coefficient_numerically_recomputed':recomputed,
        'max_recomputed_one_scale_mnu_absolute_error_ev':error if recomputed else None,
        'original_Matchete_mass_eigenbasis_validated':True,
        'live_registered_110101_CG_and_RGBeta_identity_previously_verified':True,
        'input_sha256':{key:digest(path) for key,path in p.items()},
        'source_eft1_scalar_endpoint_scale_gev':float(params['source_eft1_scale_gev']),
        'Yukawa_running_on_candidate_scalar_background_performed':False,
        'fermion_first_EFT1_threshold_matching_performed':False,
        'individual_scalar_mass_threshold_matching_performed':False,
        'full_physical_neutrino_mass_prediction_certified':False,
        'production_modified':False,
    }


def run(root: Path, out: Path, models: tuple[str, ...], *, recompute: bool=False) -> dict:
    cg_path = root/'output/audits/t3_bc_matchete_registered_cg/bc_live_registered_cg_v1/summary.json'
    check_live_cg(read_json(cg_path))
    payloads = [prepare_one(model,root,read_json(cg_path),recompute=recompute)
                for model in models]
    require(not out.exists(), 'refusing to overwrite an existing run folder: '+str(out))
    out.mkdir(parents=True,exist_ok=False)
    for x in payloads:
        name=x['model'].replace('-','m')+'_yukawa_c5.json'
        (out/name).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary={'format':FORMAT,'status':STATUS,'models_checked':len(payloads),
             'model_results':[{'model':x['model'],
                'scalar_masses_gev':x['heavy_scalar_masses_gev'],
                'yukawa2_CG_factor':x['native_epsilon_projection_factor_y2'],
                'takagi_masses_ev_descending':x['one_scale_candidate_takagi_masses_ev_descending'],
                'relative_extra_quartic_mnu_response':x['relative_extra_quartic_mnu_response'],
                'Matchete_one_loop_coefficient_numerically_recomputed':x['Matchete_one_loop_coefficient_numerically_recomputed']}
                 for x in payloads],
             'native_eigenstate_yukawa_coefficients_checked':True,
             'live_registered_quartic_CG_gate_passed':True,
             'full_physical_neutrino_mass_prediction_certified':False,
             'sequential_matching_performed':False,
             'output_directory':str(out)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return summary


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--all',action='store_true')
    mode.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--recompute-c5',action='store_true',
                        help='Re-evaluate the existing large Matchete expression (slower)')
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output-root',type=Path,default=None)
    args=parser.parse_args()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,
            'unsafe run ID')
    root=args.root.resolve()
    outroot=(args.output_root or root/'output/audits/t3_bc_mass_eigen_yukawa_c5')
    report=run(root,outroot/args.run_id,MODELS if args.all else (args.model,),
               recompute=args.recompute_c5)
    for row in report['model_results']:
        print(f"{row['model']}: y2 native CG={row['yukawa2_CG_factor']:+d}; "
              f"mX,mY={row['scalar_masses_gev']} GeV; "
              f"conditional mnu={row['takagi_masses_ev_descending']} eV; "
              f"C5 replay={row['Matchete_one_loop_coefficient_numerically_recomputed']}",flush=True)
    print('Physical sequential matching: NOT CERTIFIED')
    print('Summary:',outroot/args.run_id/'summary.json')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
