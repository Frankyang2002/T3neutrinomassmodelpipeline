#!/usr/bin/env python3
"""Independent restricted A_+0 zero-momentum UV one-loop mass.

Inputs are the fitted-reference *source couplings*, native Matchete CG
component audit, native SU(2) generator eigenvectors, and previous native
action coefficient witness. No production loop/core or provisional tensor
input is imported. The optional frozen-result comparison happens LAST.

The 'zero_portals' point is physically specified by setting both portal
couplings to zero. The 'legacy_portal_assumption' result assumes the old
benchmark's unverified portal coefficient mapping and is NOT certified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

CASE = 'T3_dS1_1_dS2_3_dF_2_alpha_p0'
BENCH = Path('Numerical/benchmarks/fitted_reference') / (CASE + '.json')
CG_AUDIT = Path('output/general_observable/lo_direct_one_loop/source_audits/component_support_audit.json')
GEN_AUDIT = Path('output/general_observable/lo_direct_one_loop/source_audits/matchete_generator_analysis.json')
ACTION = Path('output/general_observable/physical_normalisation/b1_action_contract_20261103_103545/native_action_coefficients.json')
FROZEN = Path('output/general_observable/lo_direct_one_loop/provisional/A_+0_direct_one_loop.json')


def read_json(path: Path) -> dict:
    with path.open('r', encoding='utf-8') as fh:
        return json.load(fh)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar(value: str | int | float) -> complex:
    return complex(sp.N(sp.sympify(value, locals={'I': sp.I, 'sqrt': sp.sqrt}), 35))


def complex_number(raw: object) -> complex:
    if isinstance(raw, dict) and 'real' in raw:
        return complex(raw['real'], raw['imag'])
    return complex(raw)


def matrix(raw: dict) -> np.ndarray:
    return np.array(raw['real'], dtype=float) + 1j * np.array(raw['imag'], dtype=float)


def encode(m: np.ndarray) -> dict:
    return {'real': np.asarray(m.real).tolist(), 'imag': np.asarray(m.imag).tolist()}


def native_coefficients(root: Path) -> tuple[complex, complex, complex, dict]:
    audit = read_json(root / CG_AUDIT)
    cases = [c for c in audit['cases'] if c['model'] ==
             {'dS1': 1, 'dS2': 3, 'dF': 2, 'alpha': 0}]
    if len(cases) != 1:
        raise ValueError('Expected exactly one native source component case for A_+0')
    c = cases[0]['vertices']

    def get(vertex: str, indices: list[int]) -> complex:
        entries = [z for z in c[vertex]['components'] if z['indices'] == indices]
        if len(entries) != 1:
            raise ValueError(f'Missing native CG {vertex}{indices}')
        return scalar(entries[0]['cg'])

    gen = read_json(root / GEN_AUDIT)
    triplets = [r for r in gen['representations'] if r['dimension'] == 3]
    doubles = [r for r in gen['representations'] if r['dimension'] == 2]
    if len(triplets) != 1 or len(doubles) != 1:
        raise ValueError('Missing doublet/triplet native generator witnesses')
    U3 = np.array([[scalar(s) for s in row]
                    for row in triplets[0]['U_columns_are_weight_eigenvectors']], complex)
    U2 = np.array([[scalar(s) for s in row]
                    for row in doubles[0]['U_columns_are_weight_eigenvectors']], complex)
    if np.linalg.norm(U2 - np.eye(2)) > 1e-11:
        raise ValueError('Doublet weight alignment must be rederived for this source')
    if np.linalg.norm(U3.conj().T @ U3 - np.eye(3)) > 1e-11:
        raise ValueError('Triplet charge eigenvectors are not orthonormal')
    if triplets[0]['cartan_spectrum_descending'] != ['1', '0', '-1']:
        raise ValueError('Triplet weight ordering differs from verified audit')

    # Neutrino is L component 1; neutral Dirac F is component 2;
    # S2 has Y=+1 and so T3=-1 (last triplet weight) is neutral.
    u_neutral = U3[:, 2]
    g1 = get('T3Y1CG', [1, 2, 1])
    g2_cart = np.array([get('T3Y2CG', [1, 2, s]) if any(
        x['indices'] == [1, 2, s] for x in c['T3Y2CG']['components']) else 0
        for s in (1, 2, 3)], complex)
    mix_cart = np.array([get('T3MixCG', [2, 2, 1, s]) if any(
        x['indices'] == [2, 2, 1, s] for x in c['T3MixCG']['components']) else 0
        for s in (1, 2, 3)], complex)
    # Each coefficient multiplies Bar[S2_cart^s]; S2_cart=U3 S2_weight.
    # Thus Bar[S2_cart^s]=conj(U3[s,w]) Bar[S2_weight^w].
    g2 = np.dot(g2_cart, np.conj(u_neutral))
    mix = np.dot(mix_cart, np.conj(u_neutral))
    if abs(g1) < 1e-14 or abs(g2) < 1e-14 or abs(mix) < 1e-14:
        raise ValueError('A_+0 neutral loop CG unexpectedly vanishes')
    if abs(g2 + mix) > 1e-11 * max(abs(g2), abs(mix)):
        raise ValueError('A_+0 CG relative sign differs from native record')

    contract = read_json(root / ACTION)
    a = contract['native_action_witnesses']['A_p0']['couplings_and_chiral_projectors']
    if (a['mass']['coefficient'] != '-1/1' or
        any(a[yn][h]['coefficient'] != '1/1'
            for yn in ('y1', 'y2') for h in ('original', 'hc'))):
        raise ValueError('Previously measured native mass/Yukawa coefficients changed')
    return g1, g2, mix, {
        'triplet_cartesian_to_weight_U': encode(U3),
        'neutral_triplet_column_zero_based': 2,
        'native_y2_cartesian_coefficients': encode(g2_cart),
        'native_mix_cartesian_coefficients': encode(mix_cart),
        'higgs_neutral_component_index_one_based': 2,
        'mass_and_yukawa_coefficients': 'native action: -1 Dirac MF, +1 y1,y2 and each H.c.'
    }


def divided_b0_by_feynman_parameter(mf2: float, s_lo: float, s_hi: float) -> mp.mpf:
    """[B0(mf2,s_hi)-B0(mf2,s_lo)]/(s_hi-s_lo), independent integral.

    B0(0)=-integral_0^1 dx log[(x mf2+(1-x)s)/mu2]+divergence.
    The difference needs neither subtraction scale nor cancellation of B0's.
    """
    mass = mp.mpf(str(mf2))
    lo, hi = mp.mpf(str(s_lo)), mp.mpf(str(s_hi))
    d = hi - lo
    if min(mass, lo, hi) <= 0:
        raise ValueError('Non-positive physical squared mass')
    if d == 0:
        return -mp.quad(lambda x: (1-x)/(x*mass+(1-x)*lo),
                        [0, mp.mpf('0.0001'), mp.mpf('0.001'), mp.mpf('0.01'), 1])
    return -mp.quad(lambda x: mp.log1p((1-x)*d / (x*mass+(1-x)*lo))/d,
                    [0, mp.mpf('0.0001'), mp.mpf('0.001'), mp.mpf('0.01'),
                     mp.mpf('0.1'), 1])


def calculate(y1: np.ndarray, y2: np.ndarray, masses: np.ndarray,
              lam: complex, vev: float, g1: complex, g2: complex,
              cg_mix: complex, a: float, b: float) -> dict:
    # The native mixing term is +lambda H H S1 Bar[S2] + H.c.; V=-L.
    # H^0=v/sqrt(2) means A[2,1] = -lambda*v**2*CGmix/2.
    off = -lam * (vev * vev/2.) * cg_mix
    A = np.array([[a, off.conjugate()], [off, b]], complex)
    eig = np.linalg.eigvalsh(A)
    if min(eig) <= 0:
        raise ValueError('Tachyonic physical scalar mass')

    # Chosen four-to-two-component Weyl convention from native-action witness:
    # L_Y = -nu chi [ -conj(y1*CG1) S1* ]
    #       -nu xi  [ -conj(y2*CG2) S2  ] + H.c.
    # L_M = -M chi xi + H.c., with no extra symmetry factor.
    v1 = -np.conjugate(y1 * g1)
    v2 = -np.conjugate(y2 * g2)
    mp.mp.dps = 55
    matrix_gev = np.zeros((3, 3), complex)
    per_heavy = []
    for r, m in enumerate(masses):
        L = float(divided_b0_by_feynman_parameter(float(m*m), float(eig[0]), float(eig[1])))
        # Scalar propagator <S1* S2> = A21 * divided-B0.
        coefficient = (m / (16.*math.pi**2)) * off * L
        contribution = coefficient * (np.outer(v1[:, r], v2[:, r]) +
                                      np.outer(v2[:, r], v1[:, r]))
        matrix_gev += contribution
        per_heavy.append({'MF_gev': float(m),
                          'divided_B0_gev_minus2': L,
                          'mass_insertion_coefficient_gev': encode(np.asarray(coefficient)),
                          'matrix_contribution_ev': encode(contribution*1e9)})

    # A separate small-mixing limit (remove |off|^2 in scalar eigenvalues):
    # source-propagator insertion instead of physical-mass eigenstates.
    original = sorted((a, b))
    insertion = np.zeros((3, 3), complex)
    for r, m in enumerate(masses):
        L0 = float(divided_b0_by_feynman_parameter(float(m*m), *original))
        insertion += (m/(16.*math.pi**2)*off*L0) * (
            np.outer(v1[:, r], v2[:, r]) + np.outer(v2[:, r], v1[:, r]))
    return {
        'diagonal_scalar_mass_inputs_gev2': [a, b],
        'physical_scalar_mass_squared_gev2': eig.tolist(),
        'physical_scalar_masses_gev': np.sqrt(eig).tolist(),
        'scalar_A21_gev2': encode(np.asarray(off)),
        'y1_neutral_vertex_coefficient': encode(np.asarray(g1)),
        'y2_neutral_vertex_coefficient': encode(np.asarray(g2)),
        'mix_neutral_coefficient': encode(np.asarray(cg_mix)),
        'mass_matrix_ev': encode(matrix_gev*1e9),
        'takagi_singular_values_ev_descending':
            np.linalg.svd(matrix_gev, compute_uv=False).__mul__(1e9).tolist(),
        'per_heavy': per_heavy,
        'small_mixing_mass_insertion_matrix_ev': encode(insertion*1e9),
        'small_mixing_relative_residual': float(np.linalg.norm(matrix_gev-insertion) /
            max(np.linalg.norm(matrix_gev), 1e-100)),
        'norm_mass_matrix_ev': float(np.linalg.norm(matrix_gev)*1e9),
        'nu_mass_symmetric': bool(np.linalg.norm(matrix_gev-matrix_gev.T) < 1e-18),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument('--output-dir', type=Path, required=True,
                    help='Fresh directory only. Existing evidence is never overwritten.')
    ap.add_argument('--compare-frozen', action='store_true',
                    help='Compare ONLY the legacy-portal conditional result, AFTER computing.')
    args = ap.parse_args()
    root = args.repo_root.resolve()
    out = args.output_dir.resolve()
    if out.exists():
        ap.error(f'Output directory already exists, refusing to overwrite: {out}')
    source_paths = [BENCH, CG_AUDIT, GEN_AUDIT, ACTION]
    for p in source_paths:
        if not (root/p).is_file():
            ap.error(f'Missing source witness: {root/p}')

    g1, g2, gm, witness = native_coefficients(root)
    bench = read_json(root/BENCH)
    if bench['representation'] != dict(d_s1=1, d_s2=3, d_f=2, alpha=0,
                                      shared_scalar=False):
        raise ValueError('Benchmark is not restricted A_+0')
    o = bench['base_state']['ordinary']
    y1, y2, MF = matrix(o['y1']), matrix(o['y2']), matrix(o['MF'])
    if y1.shape != (3,3) or y2.shape != (3,3) or MF.shape != (3,3):
        raise ValueError('Unexpected flavour array shape')
    if np.linalg.norm(MF-np.diag(np.diag(MF))) > 1e-9 or \
       np.max(np.abs(np.diag(MF).imag)) > 1e-9 or \
       np.min(np.diag(MF).real) <= 0:
        raise ValueError('Standalone witness expects positive real diagonal Dirac MF')
    masses = np.diag(MF).real
    lam = complex_number(o['lambdaT3'])
    vev = float(bench['numerical']['vev_gev'])
    m1, m2 = float(o['mS1Sq']), float(o['mS2Sq'])

    # Full independent source-to-loop evaluation finishes before any frozen result is read.
    zero = calculate(y1, y2, masses, lam, vev, g1, g2, gm, m1, m2)
    legacy_a = m1 - float(o['lambdaH1'])*vev*vev/2
    legacy_b = m2 - float(o['lambdaH2'])*vev*vev/2
    legacy = calculate(y1, y2, masses, lam, vev, g1, g2, gm, legacy_a, legacy_b)
    report = {
        'schema': 't3-A0-independent-charge-selected-loop-v1',
        'physical_source_certified': False,
        'reason_not_certified': ('Legacy portal-to-native-CG normalization not independently '
            'measured; installed absolute C matrix phase is convention-dependent.'),
        'scope': 'restricted A_+0 Dirac UV, p=0, 1-loop; zero extra operators',
        'not_used_in_independent_calculation': [
            'Numerical/direct_one_loop/core.py',
            'scripts/BuildT3DirectOneLoopInputs.py',
            str(FROZEN)],
        'source_sha256': {str(p): sha(root/p) for p in source_paths},
        'git_HEAD_local': subprocess.run(['git','rev-parse','HEAD'], cwd=root,
            check=False, capture_output=True, text=True).stdout.strip(),
        'native_CG_charge_witness': witness,
        'source_quartic_lambdaT3': [lam.real, lam.imag],
        'zero_portals_source_point': zero,
        'legacy_portal_assumption_conditional': legacy,
    }
    if args.compare_frozen:
        frozen = root/FROZEN
        if not frozen.is_file():
            ap.error(f'Comparison was requested but frozen report is missing: {frozen}')
        z = read_json(frozen)['total']['complex_symmetric_matrix_gev']
        reference = matrix(z)*1e9
        candidate = matrix(legacy['mass_matrix_ev'])
        diff = candidate-reference
        report['frozen_provisional_comparison_conditional'] = {
            'path': str(FROZEN), 'sha256': sha(frozen),
            'relative_frobenius_residual': float(np.linalg.norm(diff)/
                max(np.linalg.norm(reference),1e-90)),
            'absolute_frobenius_difference_ev': float(np.linalg.norm(diff)),
            'max_abs_entry_difference_ev': float(np.max(np.abs(diff))),
            'elementwise_complex_difference_ev': encode(diff),
            'comparison_cannot_certify_legacy_portal_mapping': True,
        }

    out.mkdir(parents=True, exist_ok=False)
    destination = out/'A_+0_independent_source_loop.json'
    destination.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n',encoding='utf-8')
    print('Independent source-point matrix (eV):')
    print(np.array2string(matrix(zero['mass_matrix_ev']), precision=10))
    print('Source zero-portal Takagi singular values (eV), descending:',
          zero['takagi_singular_values_ev_descending'])
    print('Conditional legacy-portal matrix (eV):')
    print(np.array2string(matrix(legacy['mass_matrix_ev']), precision=10))
    if 'frozen_provisional_comparison_conditional' in report:
        print('Legacy-only conditional relative residual:',
              report['frozen_provisional_comparison_conditional']['relative_frobenius_residual'])
    print('Evidence saved:', destination)
    print('Physical source-certified: False; independent loop: evaluated')


if __name__ == '__main__':
    main()
