#!/usr/bin/env python3
"""Evaluate A_+0 at p=0 from a LIVE Matchete scalar-action/CG witness.

No imports from the provisional UV core or previous independent-loop script.
The only external physics inputs are the native witness and source benchmark.
Select --portal-assignment benchmark-identity explicitly when comparing the
legacy RGBeta-style lambdaH1,2 numbers to Matchete's native portal couplings.

Project path: scripts/ReviewPhysicalA0NativeScalarQuadratic.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = Path('Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json')
FROZEN = Path('output/general_observable/lo_direct_one_loop/provisional/A_+0_direct_one_loop.json')
NATIVE_STATUS = 'NativeA0QuadraticActionAndRegisteredCGExtracted_IndependentLoopPending'


def require(condition: bool, detail: str) -> None:
    if not condition:
        raise ValueError(detail)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing source: {path}')
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), f'Expected JSON object: {path}')
    return value


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_complex(value: dict) -> np.ndarray:
    return np.asarray(value['real'], dtype=float) + 1j*np.asarray(value['imag'], dtype=float)


def cpayload(value: dict) -> np.ndarray:
    shape = tuple(value['shape'])
    re = np.asarray(value['real_flat'], dtype=float)
    im = np.asarray(value['imag_flat'], dtype=float)
    require(re.size == math.prod(shape) == im.size, 'Native tensor flat dimension mismatch')
    v = (re + 1j*im).reshape(shape)
    require(np.isfinite(v).all(), 'Nonfinite native tensor')
    return v


def encode(value: np.ndarray) -> dict:
    return {'real': value.real.tolist(), 'imag': value.imag.tolist()}


def loop_b0(mf2: float, ms2: float, mu2: float) -> float:
    """Exact finite p=0 scalar/fermion B0 via its Feynman-parameter primitive.

    B0_fin = - integral_0^1 dx ln[(x mf2+(1-x) ms2)/mu2].
    For nearby masses use a stable log1p expansion. No UV core imported.
    """
    require(min(mf2, ms2, mu2) > 0, 'Non-positive squared mass in B0')
    x = (ms2 - mf2)/mf2
    if abs(x) < 1e-5:
        ratio = 1 + x*(0.5 + x*(-1/6 + x*(1/12 + x*(-1/20))))
    else:
        ratio = (1+x)*math.log1p(x)/x
    return 1 - math.log(mf2/mu2) - ratio


def validate_native(native: dict) -> dict[str, np.ndarray]:
    require(native.get('format') == 't3-b1-a0-live-native-scalar-quadratic-v1'
            and native.get('status') == NATIVE_STATUS
            and native.get('model') == 'A_+0'
            and native.get('native_matchete_checked_lagrangian') is True
            and native.get('physical_source_certified') is False,
            'Native Matchete action witness missing or unexpectedly certified')
    require(native.get('source_determined_operator_coefficients') is True,
            'Native operator prefactors not measured')
    flags = native.get('restricted_environment', {}).get('effective', {})
    require(len(flags) == 4 and all(x == '0' for x in flags.values()),
            'Native action did not enforce original restricted topology')
    checks = native.get('native_component_checks', {})
    require(len(checks) == 7 and all(x is True for x in checks.values()),
            'Native source CG-component validation failed')
    expected = {'scalar1_mass': ['-1'], 'scalar2_mass': ['-1'],
                'scalar1_portal': ['1'], 'scalar2_portal_inv1': ['1/2', '1/2'],
                'scalar2_portal_inv2': ['1/2', '1/2'],
                'mixing_original': ['1'], 'mixing_hc': ['1']}
    measured = native['assembled_native_operator_terms']
    for name, coefficients in expected.items():
        require(measured[name]['coefficients'] == coefficients,
                f'Unexpected measured native coefficient: {name}')
        require(len(measured[name]['native_terms_inputform']) == len(coefficients),
                f'Missing source action terms: {name}')
    require(native['scalar_free_kinetic_coefficients'] == ['1', '1'],
            'Canonical scalar kinetic normalisation missing')
    shapes = {'T3Y1CG': (2, 2), 'T3Y2CG': (2, 2, 3), 'T3MixCG': (2, 2, 3),
              'T3HiggsPortal1CGInv1': (2, 2),
              'T3HiggsPortal2CGInv1': (2, 2, 3, 3),
              'T3HiggsPortal2CGInv2': (2, 2, 3, 3)}
    cg = {name: cpayload(native['native_tensors'][name]) for name in shapes}
    require(all(cg[name].shape == size for name, size in shapes.items()),
            'Native tensor shape mismatch')
    require(np.allclose(cg['T3Y1CG'][0, 1], 1, rtol=0, atol=1e-13),
            'Neutral y1 Clebsch is not 1')
    # Reconstruct projections afresh rather than importing prior/conditional matrices.
    h1 = cg['T3HiggsPortal1CGInv1'][1, 1]
    p21 = cg['T3HiggsPortal2CGInv1'][1, 1]
    p22 = cg['T3HiggsPortal2CGInv2'][1, 1]
    p21 = (p21 + p21.conj().T)/2
    p22 = (p22 + p22.conj().T)/2
    require(abs(h1-1) < 1e-13 and np.linalg.norm(p21) < 1e-13
            and np.linalg.norm(p22-np.eye(3)) < 1e-13,
            'Native Higgs portal projections violate measured identity/zero result')
    cm = cg['T3MixCG'][1, 1]
    cy2 = cg['T3Y2CG'][0, 1]
    require(np.linalg.norm(cm+cy2) < 1e-13,
            'Source y2/mixing axes or relative phases disagree')
    return cg


def calculate(native: dict, benchmark: dict, portal: str) -> dict:
    cg = validate_native(native)
    rep = benchmark['representation']
    require(tuple(rep[x] for x in ('d_s1', 'd_s2', 'd_f', 'alpha')) == (1, 3, 2, 0),
            'Source benchmark is not A_+0')
    ordinary = benchmark['base_state']['ordinary']
    y1, y2, mf = (source_complex(ordinary[x]) for x in ('y1', 'y2', 'MF'))
    require(y1.shape == y2.shape == mf.shape == (3, 3), 'Unexpected source flavour shape')
    require(np.linalg.norm(mf-np.diag(np.diag(mf))) < 1e-10 and
            np.max(np.abs(np.diag(mf).imag)) < 1e-12 and
            np.min(np.diag(mf).real) > 0,
            'Only the source positive diagonal flavour MF case is measured')
    mfs = mf.diagonal().real
    v = float(benchmark['numerical']['vev_gev'])
    mu = float(benchmark['base_state']['mu_gev'])
    require(v > 0 and mu > 0, 'Invalid source benchmark vev/scale')
    lam = source_complex(ordinary['lambdaT3']).item()
    require(portal in ('zero', 'benchmark-identity'), 'Unknown portal assignment')
    if portal == 'zero':
        h1 = h21 = h22 = 0.0
    else:
        # This is an EXPLICIT IDENTIFICATION of coupling values, not a
        # certified RGBeta -> Matchete mapping. Matchete has two S2 couplings.
        h1 = float(ordinary['lambdaH1'])
        h21 = 0.0
        h22 = float(ordinary['lambdaH2'])
    m1sq = float(ordinary['mS1Sq'])
    m2sq = float(ordinary['mS2Sq'])
    p21 = cg['T3HiggsPortal2CGInv1'][1, 1]
    p22 = cg['T3HiggsPortal2CGInv2'][1, 1]
    p21 = (p21+p21.conj().T)/2
    p22 = (p22+p22.conj().T)/2
    A = np.zeros((4, 4), dtype=complex)
    A[0, 0] = m1sq - v*v/2*h1*cg['T3HiggsPortal1CGInv1'][1, 1]
    A[1:, 1:] = m2sq*np.eye(3) - v*v/2*(h21*p21 + h22*p22)
    A[1:, 0] = -lam*v*v/2*cg['T3MixCG'][1, 1, :]
    A[0, 1:] = A[1:, 0].conj()
    require(np.linalg.norm(A-A.conj().T) < 1e-10, 'Scalar mass Hessian not Hermitian')
    ms2, U = np.linalg.eigh(A)
    require(np.min(ms2) > 0 and np.isfinite(ms2).all(), 'Unphysical scalar squared mass')
    mu2 = mu*mu
    mnu = np.zeros((3, 3), dtype=complex)
    pole = np.zeros_like(mnu)
    # Exact Dirac chirality-flip propagator <xi chi> ~ MF.
    # L_Y = -nu [Y2 xi S2 + Z1 chi S1*] + H.c. in the declared
    # two-component convention, from the ORIGINAL four-component y1,y2 H.c.
    # Diagonal complex scalar S = U Phi, Phi = (R+i I)/sqrt(2).
    # Summing R and I cancels YY,ZZ and gives YZ+ZY once (no 1/2).
    for r, M in enumerate(mfs):
        Y2 = np.zeros((3, 4), dtype=complex)
        Z1 = np.zeros((3, 4), dtype=complex)
        Y2[:, 1:] = -np.conj(y2[:, r, None]*cg['T3Y2CG'][0, 1, :][None, :])
        Z1[:, 0] = -np.conj(y1[:, r]*cg['T3Y1CG'][0, 1])
        gY = Y2 @ U
        gZ = Z1 @ U.conj()
        for k, mass2 in enumerate(ms2):
            pair = (np.outer(gY[:, k], gZ[:, k])+
                    np.outer(gZ[:, k], gY[:, k]))
            pole += M/(16*math.pi**2)*pair
            mnu += M/(16*math.pi**2)*loop_b0(M*M, float(mass2), mu2)*pair
    # This cancellation is physical: no local tree neutrino mass counterterm.
    pole_residual = float(np.linalg.norm(pole))
    require(pole_residual < 1e-15, f'UV divergence did not cancel: {pole_residual}')
    require(np.linalg.norm(mnu-mnu.T) < 1e-15, 'Neutrino mass matrix not symmetric')
    mass_ev = 1e9*mnu
    singular = np.linalg.svd(mass_ev, compute_uv=False)
    return {
        'model': 'A_+0',
        'status': 'independent_A0_loop_evaluated_from_live_native_scalar_action_not_certified',
        'portal_assignment': portal,
        'portal_coupling_values': {'lambdaH1Inv1': h1, 'lambdaH2Inv1': h21,
                                   'lambdaH2Inv2': h22},
        'benchmark_portal_equivalence_physically_certified': False,
        'normalisation': 'Dirac mass chirality flip once, complex scalar real-mode sum, 1/(16*pi^2)',
        'scalar_native_cartesian_mass_squared_gev2': encode(A),
        'scalar_eigenmass_squared_gev2': ms2.tolist(),
        'scalar_eigenmass_gev': np.sqrt(ms2).tolist(),
        'neutrino_majorana_mass_matrix_ev': encode(mass_ev),
        'neutrino_takagi_singular_values_ev_descending': singular.tolist(),
        'uv_pole_absolute_frobenius_gev': pole_residual,
        'physical_source_certified': False,
        'no_eft_or_running': True,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', type=Path, default=ROOT)
    ap.add_argument('--native-json', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--portal-assignment', choices=('zero', 'benchmark-identity'), default='zero')
    ap.add_argument('--compare-frozen', action='store_true')
    args = ap.parse_args()
    root = args.repo_root.resolve()
    def path(p: Path) -> Path:
        return p if p.is_absolute() else root / p
    out = path(args.output_dir)
    if out.exists():
        print(f'BLOCKED: refusing to overwrite existing evidence directory {out}', file=sys.stderr)
        return 2
    try:
        require(not(args.compare_frozen and args.portal_assignment != 'benchmark-identity'),
                '--compare-frozen needs --portal-assignment benchmark-identity')
        npath = path(args.native_json)
        bpath = root / BENCHMARK
        native = read_json(npath)
        record = calculate(native, read_json(bpath), args.portal_assignment)
        record['evidence_sources_sha256'] = {
            str(npath): sha(npath), str(BENCHMARK): sha(bpath)}
        if args.compare_frozen:
            old = read_json(root / FROZEN)
            old_matrix = source_complex(old['total']['complex_symmetric_matrix_gev'])*1e9
            new_matrix = source_complex(record['neutrino_majorana_mass_matrix_ev'])
            relative = float(np.linalg.norm(new_matrix-old_matrix)/np.linalg.norm(old_matrix))
            record['conditional_frozen_comparison'] = {
                'relative_frobenius_residual': relative,
                'not_a_source_certification': True,
                'frozen_path': str(FROZEN)}
        out.mkdir(parents=True, exist_ok=False)
        target = out / 'A_+0_native_quadratic_independent_loop.json'
        target.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    except (ValueError, KeyError, TypeError, OSError, np.linalg.LinAlgError) as exc:
        print(f'BLOCKED: {exc}', file=sys.stderr)
        return 3
    print('A_+0 native Cartesian scalar A (GeV^2):')
    print(source_complex(record['scalar_native_cartesian_mass_squared_gev2']))
    print('Independent A_+0 mass matrix (eV):')
    print(source_complex(record['neutrino_majorana_mass_matrix_ev']))
    print('Takagi singular values (eV), descending:',
          record['neutrino_takagi_singular_values_ev_descending'])
    if args.compare_frozen:
        print('CONDITIONAL frozen relative residual:',
              record['conditional_frozen_comparison']['relative_frobenius_residual'])
    print('Evidence saved:', target)
    print('Physical source-certified: False (benchmark portal identification and installed spinor phase pending)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
