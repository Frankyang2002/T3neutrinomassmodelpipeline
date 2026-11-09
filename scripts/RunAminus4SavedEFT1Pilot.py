"""A(-4) *provisional* EFT1 Wilson transport from one SAVED sequential benchmark.

This does not rerun the UV fit or assume a Matchete/RGBeta basis identification.
It uses one existing fixed-order sequential run for the F-tree source and
for the *same run's* EFT1 boundary, then RE-runs the archived native EFT1
renormalisable RGEs to make intermediate tensor snapshots. Refuses an
inconsistent run/beta/source bundle or a mismatched EFT1 endpoint.

Typical run, from repository root:
    python -X utf8 -m scripts.RunAminus4SavedEFT1Pilot --run-id batch16_v1

Expected existing inputs (do not overwrite):
    output/general_sequential_pilot/all16_dynamic_scalar_threshold/\
        batch16_v1/A_m4/{result.json,dynamic_matching.npz}
    output/general_observable/fermion_first_all16_v1/A_m4/data/\
        eft1_after_F_wilson_seed.json
    output/studies/general_scalar_extended/complete_scalar_symbolic/\
        full16_general_v1/A_m4/eft1_rgbeta.json

Outputs remain explicitly FORMAL. No scalar-threshold matching or SMEFT
running is calculated. The old fixed-order log must never be added to the
new evolved C5: the two are alternative EFT treatments, not contributions.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEQUENTIAL = ROOT / 'output/general_sequential_pilot/all16_dynamic_scalar_threshold'
SOURCE_SEED = (ROOT / 'output/general_observable/fermion_first_all16_v1/'
               'A_m4/data/eft1_after_F_wilson_seed.json')
BETA = (ROOT / 'output/studies/general_scalar_extended/complete_scalar_symbolic/'
        'full16_general_v1/A_m4/eft1_rgbeta.json')
DEST = ROOT / 'output/eft1_stageB_saved_a_m4'
SM_KEYS = ('gY','g2','g3','lambdaH','yu','yd','ye')


class SavedPilotError(ValueError):
    pass


def _check(ok: bool, reason: str) -> None:
    if not ok:
        raise SavedPilotError(reason)


def read_json(path: Path) -> dict:
    item = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    _check(isinstance(item, dict), f'Expected JSON object: {path}')
    return item


def decode_complex(value: Any) -> np.ndarray:
    """Accept supported saved matrix representations without losing phases."""
    if isinstance(value, Mapping) and {'real','imag'} <= value.keys():
        return (np.asarray(value['real'], dtype=float)
                + 1j * np.asarray(value['imag'], dtype=float))
    if isinstance(value, list):
        return np.asarray([decode_complex(v) for v in value], dtype=complex)
    return np.asarray(value, dtype=complex)


def same_numbers(actual: Any, expected: Any, *, rtol: float,
                 atol: float, name: str) -> float:
    a, b = decode_complex(actual), decode_complex(expected)
    _check(a.shape == b.shape, f'{name}: incompatible shapes {a.shape}/{b.shape}')
    _check(np.isfinite(a).all() and np.isfinite(b).all(), f'{name}: nonfinite values')
    scale = atol + rtol * np.maximum(np.abs(a), np.abs(b))
    max_units = float(np.max(np.abs(a-b) / scale))
    _check(max_units <= 1., f'{name}: inconsistent saved/transport values, '
                               f'max tolerance units={max_units:.6g}')
    return float(np.max(np.abs(a-b)))


def repackage_saved_source(saved_path: Path, result: dict,
                           original_seed: Path, export_npz: Path) -> dict:
    """Read *existing* native dynamic matching output; verify its source hash."""
    _check(result.get('status') == 'ConditionalDynamicFToBothScalarsToSMEFTCompleted'
           and result.get('model') == 'A_-4',
           'Input is not an existing completed A(-4) sequential benchmark')
    flags = result.get('flags', {})
    for k in ('UV_evolved_fresh','all_heavy_fermions_integrated_out_at_F',
              'LLSS_F_boundary_evaluated_from_current_Yukawas'):
        _check(flags.get(k) is True, f'Saved dynamic source is missing {k}')
    with np.load(saved_path, allow_pickle=False) as data:
        _check(str(data['model'].item()) == 'A_-4', 'NPZ model is not A(-4)')
        C = np.array(data['C_LLSS_muF'], dtype=np.complex128, copy=True)
        by = np.array(data['C_LLSS_by_source_muF'], dtype=np.complex128, copy=True)
        mu_f = float(data['muF_gev'].item())
        mu_s = float(data['muS_gev'].item())
        old_log = np.array(data['C5_direct_log'], dtype=np.complex128, copy=True)
    _check(C.shape == (3,3,2,2,8,8), f'Unexpected A(-4) F source shape {C.shape}')
    _check(by.ndim == 7 and by.shape[1:] == C.shape, 'Source decomposition shape mismatch')
    _check(np.isfinite(C).all() and np.isfinite(by).all(), 'Nonfinite F source')
    _check(np.allclose(by.sum(axis=0), C, rtol=1e-12, atol=1e-27),
           'Stored per-source components fail to reconstruct full LLSS')
    _check(np.allclose(C, C.transpose(1,0,3,2,4,5), rtol=1e-11, atol=1e-27)
           and np.allclose(C, C.swapaxes(4,5), rtol=1e-11, atol=1e-27),
           'Source field interchange symmetry fails')
    ft = result['fermion_tree']
    _check(ft.get('shape') == list(C.shape) and
           ft.get('source_count') == len(by) and
           np.isclose(float(ft.get('norm_gev_inverse')), np.linalg.norm(C),
                      rtol=1e-11, atol=1e-27),
           'Saved fermion-tree metadata differs from NPZ')
    _check(mu_f == result['scales_gev']['F'] and
           mu_s == result['scales_gev']['S'] and mu_f > mu_s > 0,
           'Saved NPZ thresholds differ from result.json')
    _check(old_log.shape == (3,3) and np.isfinite(old_log).all(),
           'Invalid old fixed-order C5')
    digest = sha256(original_seed.read_bytes() + C.tobytes()).hexdigest()
    _check(digest == ft.get('source_hash'),
           'Matchete tree seed or saved source differs from result.json hash; '
           'do NOT combine different revisions')
    indices = tuple(int(row['source_index']) for row in
                    ft['source_geometry']['sources'])
    _check(len(indices) == len(by), 'Per-source index metadata mismatch')
    labels = ['NewScalar1']*2 + ['NewScalar2']*6
    # Identical layout and flags to the repo's write_provisional_source.
    meta = {
        'format': 't3-joint-scalar-source-bridge-v1',
        'model':'A_-4',
        'status':'FormalLLSSExported_PhysicalNormalizationUncertified',
        'source_boundary_path':str(original_seed),
        'source_boundary_sha256':digest,
        'mu_f_gev':mu_f, 'n_heavy_real':8,
        'source_scalar_real_labels':labels,
        'PL_source_indices':list(indices),
        'gamma_cc_weyl_normalization_certified':False,
        'eft1_muS_running_included':False,
        'source_field_order':'1-based Higgs 1..4; heavy real 5..12',
        'source_tensor_shape':'(flavor_p,flavor_q,lepton_su2_i,lepton_su2_j,heavy_a,heavy_b)',
        'full_physical_C5_computed':False,
        'note':'Verified original saved Matchete F-tree source; spinor/basis map not certified.',
    }
    np.savez_compressed(export_npz, C_LLSS_formal_muF=C,
                        C_LLSS_by_PL_source_formal_muF=by,
                        matchete_PL_source_indices=np.asarray(indices, dtype=np.int64))
    export_npz.with_suffix('.json').write_text(
        json.dumps(meta,indent=2)+'\n', encoding='utf-8')
    return {'source_hash':digest,'mu_F_gev':mu_f,'mu_S_gev':mu_s,
            'source_norm_gev_inverse':float(np.linalg.norm(C)),
            'old_fixed_order_direct_log_frobenius_gev_inverse':float(np.linalg.norm(old_log)),
            'original_saved_input':str(saved_path)}


def make_eft1_initial_document(result: dict) -> dict:
    """Obtain the exact same-run EFT1 state at F, NEVER UV-scale values."""
    from Numerical.running.Aminus4EFT1Couplings import BETA_NAMES
    at_f = result['rge']['renormalisable_parameters_at_F']
    eft1 = {k:at_f[k] for k in BETA_NAMES}
    _check(set(eft1) == set(BETA_NAMES), 'Incomplete 17-key EFT1 initial point')
    return {'mu_gev':float(result['scales_gev']['F']),
            'sm':{k:eft1.pop(k) for k in SM_KEYS}, **eft1}


def check_eft1_endpoint(numerical: Mapping, saved: Mapping,
                         *, rtol: float = 3e-4, atol: float = 1e-8) -> dict:
    """Require native re-integration to reproduce saved S-scale endpoint."""
    from Numerical.running.Aminus4EFT1Couplings import BETA_NAMES
    errors = {}
    _check(set(BETA_NAMES) <= set(numerical) and set(BETA_NAMES) <= set(saved),
           'Missing full EFT1 endpoints')
    for k in BETA_NAMES:
        errors[k] = same_numbers(numerical[k], saved[k],
                                 rtol=rtol, atol=atol, name=f'EFT1@S/{k}')
    return {'match':True,'rtol':rtol,'atol':atol,
            'maximum_absolute_residual_by_coupling':errors}


def save_same_run_couplings(initial_doc: dict, beta_path: Path,
                            mu_s: float, samples: int, output_dir: Path,
                            saved_at_s: dict) -> dict:
    """Regenerate full EFT1 RGBeta reference tensors, not borrowed snapshots."""
    from scripts.ExportAminus4EFT1Couplings import make_state
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.Aminus4EFT1Couplings import (
        BETA_NAMES, tensor_point, state_point, quartic_columns, source_contract,
    )
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    initial = make_state(initial_doc)
    beta = read_json(beta_path)
    _check(beta.get('status') == 'SuccessCandidateMassBeta'
           and beta.get('metadata',{}).get('GeneralScalarModel') == 'A_-4'
           and beta.get('metadata',{}).get('GeneralScalarStage') == 'EFT1'
           and set(beta.get('report_betas',{})) == set(BETA_NAMES),
           'Native archived A(-4) EFT1 beta metadata/keys mismatch')
    _check(not ({'MF','y1','y2','y1Conj','y2Conj'} & set(beta['report_betas'])),
           'Heavy fermion remains in EFT1 running')
    source = source_contract(ROOT)  # Audit native model without treating basis as certified.
    save_mu = np.geomspace(float(initial.mu_gev), mu_s, samples)
    running = run_intermediate_scalar_segment(
        initial, beta, mu_s, save_scales_gev=save_mu,
        rtol=1e-7, atol=1e-10)
    _check(running.solver_success and running.n_points == samples,
           'EFT1 numerical solve did not complete requested snapshots')
    end = _state_values(running.final_state)
    end_check = check_eft1_endpoint(end, saved_at_s)
    columns = quartic_columns()
    snapshots = [tensor_point(float(mu), state_point(running.state_at_index(i)),
                              columns=columns)
                 for i,mu in enumerate(running.mu_gev)]
    np.savez_compressed(output_dir/'A_m4_eft1_tensors.npz',
        mu_gev=np.asarray(running.mu_gev,dtype=float),
        Y=np.stack([p.Y for p in snapshots]),
        lambda4=np.stack([p.lambda4 for p in snapshots]),
        **{f'generator_{kind}_{rep}':np.stack(snapshots[0].gauge[kind][rep])
           for kind in ('SU2','U1','SU3') for rep in ('fermion','scalar')})
    def encode(v):
        if isinstance(v, np.ndarray):
            return {'real':v.real.tolist(),'imag':v.imag.tolist()}
        if isinstance(v, complex):
            return {'real':float(v.real),'imag':float(v.imag)}
        return float(v)
    info = snapshots[0].to_machine()
    info.update({
        'sampled_scales_gev':[float(x) for x in running.mu_gev],
        'both_heavy_scalars_active_until_mu_s':True,
        'running_couplings':[{k:encode(v) for k,v in point.couplings.items()}
                             for point in snapshots],
        'rge_payload_sha256':sha256(beta_path.read_bytes()).hexdigest(),
        'source_contract':source,
        'saved_benchmark_eft1_endpoint_verified':True,
        'saved_endpoint_comparison':end_check,
        'matchete_scalar_and_Wilson_basis_map_certified':False,
        'physical_full_Wilson_running_ready':False,
        'why_not_fully_certified':'Spinor normalization, native Matchete-to-RGBeta map, operator closure.',
    })
    (output_dir/'A_m4_eft1_tensor_basis.json').write_text(
        json.dumps(info,indent=2,allow_nan=False)+'\n', encoding='utf-8')
    return {'source':str(beta_path),'beta_sha256':info['rge_payload_sha256'],
            'n_points':samples,'solver_nfev':int(running.nfev),
            'mu_S_gev':mu_s,'endpoint_validation':end_check,
            'native_source_fingerprints':source.get('source_sha256',{})}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id',default='batch16_v1',
                        help='Exact existing sequential batch directory (default: batch16_v1)')
    parser.add_argument('--saved-folder',type=Path,
                        help='Override saved A_m4 folder with result.json + dynamic_matching.npz')
    parser.add_argument('--seed',type=Path,default=SOURCE_SEED)
    parser.add_argument('--beta',type=Path,default=BETA)
    parser.add_argument('--output',type=Path,
                        help='New additive output directory (default: output/eft1_stageB_saved_a_m4/<run-id>)')
    parser.add_argument('--samples',type=int,default=7)
    parser.add_argument('--prepare-only',action='store_true',
                        help='Verify saved boundary and rebuild EFT1 tensors without Wilson ODE')
    parser.add_argument('--rtol',type=float,default=1e-6,
                        help='Wilson ODE relative tolerance; not the EFT1 beta tolerance')
    args=parser.parse_args(argv)
    _check(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,
           'Unsafe/invalid run ID')
    _check(args.samples >= 3, 'At least 3 native running points required')
    _check(0 < args.rtol < 1, 'Invalid Wilson ODE relative tolerance')
    if args.output is None:
        args.output = DEST/args.run_id
    inp = (Path(args.saved_folder) if args.saved_folder else
           SEQUENTIAL/args.run_id/'A_m4')
    _check(not args.output.exists(), f'Refusing to overwrite output: {args.output}')
    saved = read_json(inp/'result.json')
    _check(saved.get('model')=='A_-4' and saved.get('status')==
           'ConditionalDynamicFToBothScalarsToSMEFTCompleted',
           'Must select completed A(-4) dynamic sequential benchmark')
    mu_f=float(saved['scales_gev']['F'])
    mu_s=float(saved['scales_gev']['S'])
    _check(0 < mu_s < mu_f and
           saved['scales_gev']['UV']>mu_f>mu_s>saved['scales_gev']['low']>0,
           'Invalid original UV/F/S/low hierarchy')
    # Strict general-Lagrangian gate: the matched Matchete EFT1 must retain
    # BOTH independent Higgs-S2 portal invariants before this transport can
    # be advertised as deriving from the general UV/EFT1 operator inventory.
    # The archived source currently misses lambdaH2Inv1, so STOP here.
    from scripts.AuditAminus4GeneralLagrangian import audit_source
    general_audit = audit_source(args.seed.parent.parent, args.beta)
    if not general_audit['full_general_renormalizable_eft1_source_ready']:
        raise SavedPilotError(
            'General-Lagrangian matched-source gate failed; missing quartics: '
            + ', '.join(general_audit['missing_matched_tree_quartic_symbols'])
            + '. Saved FORMAL results are preserved; regenerate and reconcile '
              'the Matchete full-general EFT1 source before further running. '
              'Run: python -X utf8 -m scripts.AuditAminus4GeneralLagrangian'
        )
    # Preserve all prior numerical outputs. This is an additive run directory.
    args.output.mkdir(parents=True,exist_ok=False)
    source_path=args.output/'A_m4_formal_F_source.npz'
    source_summary=repackage_saved_source(inp/'dynamic_matching.npz',saved,
                                          args.seed,source_path)
    initial_doc=make_eft1_initial_document(saved)
    coupling_summary=save_same_run_couplings(
        initial_doc, args.beta, mu_s, args.samples, args.output,
        saved['rge']['renormalisable_parameters_at_S'])
    report = {'model':'A_-4','general_lagrangian_source_audit': general_audit,
              'source':source_summary,
              'couplings':coupling_summary,
              'reproduces_saved_sequential_benchmark':True,
              'physical_Wilson_normalization_certified':False,
              'matchete_RGBeta_basis_map_certified':False,
              'operator_class_closure_certified':False,
              'source_sequential_path':str(inp),
              'generated_utc':datetime.now(timezone.utc).isoformat(),
              'old_fixed_order_direct_log_is_alternative_not_additive':True}
    (args.output/'input_validation.json').write_text(
        json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Native A(-4) F-tree source hash VERIFIED')
    print(f'Actual same-run EFT1 endpoint VERIFIED at muS={mu_s:g} GeV')
    if args.prepare_only:
        print('Saved source/couplings only; --prepare-only requested')
        return 0
    from Numerical.running.T3EFT1WilsonODETransport import transport
    summary,C = transport(
        source_path, args.output/'A_m4_eft1_tensors.npz',
        couplings_json=args.output/'A_m4_eft1_tensor_basis.json',
        allow_formal_identity=True, rtol=args.rtol)
    summary['input_validation_file']='input_validation.json'
    summary['original_sequential_fixed_order_direct_log_not_added']=True
    (args.output/'pilot_summary.json').write_text(
        json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    np.savez_compressed(args.output/'formal_master_C_at_muS.npz',C=C)
    print('A(-4) FORMAL full-tensor ODE completed; NOT a physical C5 prediction.')
    print('C5 formal 3x3 matrix (GeV^-1):')
    print(np.asarray(summary['C5_LLHH_formal_gev_inverse']['real']) +
          1j*np.asarray(summary['C5_LLHH_formal_gev_inverse']['imag']))
    print('Summary:',args.output/'pilot_summary.json')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
