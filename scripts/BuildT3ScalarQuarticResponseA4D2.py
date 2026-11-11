"""Build full 12^4 real-scalar quartic response at MF for A_-4/D_-2.

Only uses previously implemented Matchete/RGBeta reference tensor functions.
Checks direct reconstruction after non-unit multi-coupling perturbations;
outputs a response NPZ and report JSON. DOES NOT certify the physical
native/RGBeta map, running mass Hessian, or scalar threshold matching.
No production solver or prior output is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from Numerical.running.T3ScalarQuarticLinearResponse import (
    ScalarQuarticLinearResponse, symmetry_residual, scalar,
)

BASE = ('lambdaH', 'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12',
        'lambdaS1Adj', 'lambdaS2Adj', 'lambdaH1Adj', 'lambdaH2Adj', 'lambda12Adj')
D_EXTRA = ('lambdaT3GeneralQ000004I1', 'lambdaT3GeneralQ000013I1',
           'lambdaT3GeneralQ001102I1', 'lambdaT3GeneralQ020101I1',
           'lambdaT3GeneralQ110002I1')


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str):
    if not condition:
        raise ValueError(message)


def build(root: Path, model: str, output: Path) -> dict:
    from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
        scalar_quartic_from_matchete, add_d2_extra_hhss, add_d2_extra_scalar_only,
    )
    from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import heavy_source_from_reference
    from Numerical.running.GeneralT3DynamicFermionWilson import uv_from_json, nonmixed_boundary
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis

    require(model in ('A_-4', 'D_-2'), 'Only A_-4 and D_-2 supported')
    registry_path = root / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    registry = json.loads(registry_path.read_text(encoding='utf-8-sig'))
    extras = registry['models'][model]['new_rgbeta_symbols']
    names_extra = {x['name'] for x in extras}
    require(names_extra == (set() if model == 'A_-4' else set(D_EXTRA)),
            f'{model}: registry extra-coupling inventory changed')
    snap = load_snapshot(root, model)
    seed_path = snap.quartic_seed
    seed = json.loads(seed_path.read_text(encoding='utf-8-sig'))
    fstate = snap.trajectory.at_fermion_threshold()
    uv = uv_from_json(model, snap.uv_values, mu_f=snap.running.mu_f)
    boundary, _ = nonmixed_boundary(model, uv, repo_root=root, mu_f=snap.running.mu_f)
    # SourceBoundary exposes `path`, not `seed_path` (verified against current class).
    fermion_seed_path = Path(boundary.path)
    require(fermion_seed_path.is_file(),
            f'{model}: fermion Matchete tree seed missing: {fermion_seed_path}')
    source_seed = json.loads(fermion_seed_path.read_text(encoding='utf-8-sig'))
    native_R, _ = heavy_source_from_reference(model, source_seed)
    mass_R = np.eye(8)
    if model == 'D_-2':
        mass_R = np.asarray(build_special_basis(model, snap.running, boundary)[
            'R_original_from_boundary'], float)
    require(native_R.shape == mass_R.shape == (8, 8), 'Unexpected heavy rotation dimensions')
    for name, R in [('native_to_reference', native_R), ('original_to_mass', mass_R)]:
        require(np.allclose(R.T @ R, np.eye(8), rtol=1e-10, atol=1e-10),
                'Invalid orthogonal rotation ' + name)

    # Depend only on actual present couplings; never invent a missing potential
    # term. The source constructor itself checks all required base couplings.
    relevant = [name for name in BASE if name in fstate]
    relevant.append('lambdaT3')
    relevant.extend(D_EXTRA if model == 'D_-2' else ())
    require(len(relevant) == len(set(relevant)) and all(k in fstate for k in relevant),
            'Missing running scalar coupling from reference snapshot')
    complex_names = {'lambdaT3'} | (set(D_EXTRA) if model == 'D_-2' else set())
    # The five D directions may be real in the numerical state; probe their
    # complex orientations regardless, since the potential is Hermitian.

    def direct(values):
        q, _ = scalar_quartic_from_matchete(seed, values, model, mass_R,
                                            native_source_from_reference=native_R)
        if model == 'D_-2':
            add_d2_extra_hhss(q, root, values, mass_R)
            add_d2_extra_scalar_only(q, values, mass_R)
        return q

    response = ScalarQuarticLinearResponse.construct(
        fstate, direct, relevant, complex_names, tolerance=2e-8)
    # Verify against independent fixed-state variations, including complex
    # phases, rather than comparing only a decomposition against itself.
    rng = np.random.default_rng(1326)
    tests = []
    for scale in (0.03, 0.11):
        trial = dict(fstate)
        for name in relevant:
            base = scalar(fstate[name])
            delta = scale * (rng.standard_normal() + (1j*rng.standard_normal() if name in complex_names else 0))
            trial[name] = {'real': (base + delta).real, 'imag': (base + delta).imag} if name in complex_names else float((base+delta).real)
        rel = response.compare(trial)
        tests.append(rel)
        require(rel < 2e-8, f'{model}: direct tensor reproduction failed: {rel}')

    output = output.resolve()
    require(not output.exists(), 'Refusing to overwrite existing output: ' + str(output))
    output.mkdir(parents=True, exist_ok=False)
    npz = output / 'quartic_response.npz'
    np.savez_compressed(npz, reference_tensor=response.reference_tensor,
                        response_columns=response.columns,
                        native_rotation=native_R, heavy_mass_rotation=mass_R)
    report = {
        'schema': 't3-batch3b-full-real-quartic-response-v1',
        'model': model, 'status': 'FullReferenceQuarticAffineResponseChecked_PhysicalNormalizationOpen',
        'mu_f_gev': snap.running.mu_f,
        'tensor': response.metadata(),
        'independent_perturbation_residuals': tests,
        'basis_fixed_during_running_required': True,
        'input_sha256': {'quartic_seed': sha256(seed_path),
                         'fermion_tree_seed': sha256(fermion_seed_path),
                         'rgbeta_registry': sha256(registry_path)},
        'output_npz_sha256': sha256(npz),
        'D2_extra_quartic_count': len(extras),
        'D2_missing_quadratic': registry['models'][model]['quadratics_not_implemented'],
        'normalization_warning': 'Source-native topology and RGBeta reference invariants not independently identified; do not use as physically certified potential',
        'complete_mass_hessian_certified': False,
        'joint_scalar_threshold_certified': False,
        'production_solver_modified': False,
    }
    (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return report


def self_test() -> None:
    from Numerical.running.T3ScalarQuarticLinearResponse import ScalarQuarticLinearResponse
    def direct(data):
        a = float(data['a'])
        z = scalar(data['z'])
        x = np.zeros((12,12,12,12),float)
        x[0,0,0,0] = 2*a + 3*z.real
        x[1,1,1,1] = -4*z.imag
        return x
    ref = {'a': 1.25, 'z': {'real': 0.2, 'imag': -0.7}}
    res = ScalarQuarticLinearResponse.construct(ref,direct,['a','z'],{'z'})
    trial = {'a':0.6,'z':{'real':-0.13,'imag':0.21}}
    assert res.compare(trial) < 1e-14
    assert len(res.directions) == 3
    print('PASS: affine real/complex directions, full tensor symmetry, non-unit verification')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--model',choices=['A_-4','D_-2','both'],default='both')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test();return
    if args.project is None or args.output is None:
        parser.error('--project and --output required')
    root=args.project.resolve()
    require(root.is_dir(),'Missing project root')
    mods = ('A_-4','D_-2') if args.model == 'both' else (args.model,)
    for model in mods:
        report=build(root, model, args.output/model)
        print(model,report['status'],
              'directions=',len(report['tensor']['directions']),
              'max_relative_test=',max(report['independent_perturbation_residuals']),flush=True)
    print('No existing result or solver modified. Physical normalisation remains open.')


if __name__=='__main__':
    main()
