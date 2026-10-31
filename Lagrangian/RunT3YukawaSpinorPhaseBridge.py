#!/usr/bin/env python3
"""Read-only A(0) native charge-conjugation + explicit Weyl phase witness.

This is not another CG/normalization audit. Native Matchete tests its actual
charge-conjugation and H.c. algebra; independent exact 4x4 matrices prove the
four-component-to-two-component identity in ONE explicitly chosen basis:

  gamma^0 = [[0, I], [I, 0]], gamma^i = [[0, sigma_i], [-sigma_i, 0]]
  C = i gamma^2 gamma^0 = diag(epsilon, -epsilon), epsilon_12 = +1.

The overall alternative C -> -C and a simultaneous field rephasing are
conventions, and RGBeta's abstract eps has not been independently assigned
numeric fundamental components. Those freedoms are reported explicitly;
the script does not claim a basis-independent absolute coupling certificate.

Usage (repository root, with Matchete available through wolframscript):
  python Lagrangian/RunT3YukawaSpinorPhaseBridge.py --self-test
  python Lagrangian/RunT3YukawaSpinorPhaseBridge.py --all

Writes:
  Reports/output/workstream_a/yukawa_spinor_phase/native_matchete/A_p0.json
  Reports/output/workstream_a/yukawa_spinor_phase/native_matchete/A_p0.log
  Reports/output/workstream_a/yukawa_spinor_phase/a0_spinor_phase.json
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import TypeAlias

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('Reports/output/workstream_a/yukawa_spinor_phase')
NATIVE = OUT / 'native_matchete/A_p0.json'
DRIVER = Path('Lagrangian/RunT3YukawaSpinorPhaseNative.wl')
PRIOR = Path('Reports/output/workstream_a/yukawa_a0_action_bridge/a0_action_bridge.json')
STATUS = 'NativeChargeConjugationCheckedCanonicalWeylPhaseFixedByChoicePhysicalDictionaryPending'
NATIVE_STATUS = 'NativeMatcheteAbstractChargeConjugationAndHcPassed'
NATIVE_KEYS = (
    'gamma_cc_square_minus_one',
    'gamma_cc_transpose_antisymmetric',
    'cconj_field_definition_in_bilinear',
    'cconj_flips_left_projector',
    'hermitian_conjugate_matches_explicit_cconj',
    'distinct_native_y1_and_y1conj_present',
)


class PhaseError(ValueError):
    pass


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise PhaseError(reason)


# Explicit Grassmann (exterior) product; coefficient symbols are commutative,
# spinor components anticommute. Monomial basis uses fixed ordering.
Monomial: TypeAlias = tuple[str, ...]
Poly: TypeAlias = dict[Monomial, sp.Expr]
NAMES = (
    'l0', 'l1', 'chi0', 'chi1',
    'lb0', 'lb1', 'chib0', 'chib1', 'xi0', 'xi1', 'xib0', 'xib1',
)
ORDER = {name: idx for idx, name in enumerate(NAMES)}
DAGGER = {
    'l0': 'lb0', 'l1': 'lb1',
    'chi0': 'chib0', 'chi1': 'chib1',
    'lb0': 'l0', 'lb1': 'l1',
    'chib0': 'chi0', 'chib1': 'chi1',
    'xi0': 'xib0', 'xi1': 'xib1',
    'xib0': 'xi0', 'xib1': 'xi1',
}


def variable(name: str) -> Poly:
    require(name in ORDER, 'Unknown Grassmann component ' + name)
    return {(name,): sp.Integer(1)}


def product(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            raw = ma + mb
            if len(set(raw)) != len(raw):
                continue
            swaps = sum(ORDER[raw[i]] > ORDER[raw[j]]
                        for i in range(len(raw)) for j in range(i + 1, len(raw)))
            key = tuple(sorted(raw, key=ORDER.__getitem__))
            out[key] = sp.expand(out.get(key, sp.Integer(0)) + (-1)**swaps * ca * cb)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def add(*items: Poly) -> Poly:
    out: Poly = {}
    for item in items:
        for key, val in item.items():
            out[key] = sp.expand(out.get(key, sp.Integer(0)) + val)
    return {k: v for k, v in out.items() if sp.simplify(v) != 0}


def scale(coef: sp.Expr, p: Poly) -> Poly:
    return {k: sp.expand(coef * v) for k, v in p.items()
            if sp.simplify(coef * v) != 0}


def hc(p: Poly) -> Poly:
    # Hermitian conjugation reverses the *order* of the Grassmann factors.
    result: Poly = {}
    for key, coef in p.items():
        term: Poly = {(): sp.conjugate(coef)}
        for component in reversed(key):
            term = product(term, variable(DAGGER[component]))
        result = add(result, term)
    return result


def components_bilinear(left: list[Poly], matrix: sp.Matrix,
                        right: list[Poly]) -> Poly:
    require(len(left) == len(right) == matrix.rows == matrix.cols == 4,
            'Expected 4-component Dirac bilinear')
    terms = [scale(matrix[i, j], product(left[i], right[j]))
             for i in range(4) for j in range(4) if matrix[i, j] != 0]
    return add(*terms)


def derive_weyl_identity() -> dict:
    eye = sp.eye(2)
    z = sp.zeros(2)
    sigma = [sp.Matrix([[0, 1], [1, 0]]),
             sp.Matrix([[0, -sp.I], [sp.I, 0]]),
             sp.diag(1, -1)]

    def block(a: sp.Matrix, b: sp.Matrix,
              c: sp.Matrix, d: sp.Matrix) -> sp.Matrix:
        return sp.BlockMatrix([[a, b], [c, d]]).as_explicit()

    gammas = [block(z, eye, eye, z)] + [block(z, s, -s, z) for s in sigma]
    metric = sp.diag(1, -1, -1, -1)
    ident4 = sp.eye(4)
    for mu in range(4):
        for nu in range(4):
            require(gammas[mu]*gammas[nu] + gammas[nu]*gammas[mu]
                    == 2 * metric[mu, nu] * ident4,
                    'Clifford-algebra verification failed')
    c = sp.simplify(sp.I*gammas[2]*gammas[0])
    eps = sp.Matrix([[0, 1], [-1, 0]])
    require(c == sp.diag(1, 1, 1, 1) * block(eps, z, z, -eps),
            'Weyl charge-conjugation matrix not in declared convention')
    require(c.T == -c and c*c == -ident4,
            'Charge-conjugation transpose/square failed')
    for g in gammas:
        require(sp.simplify(c * g.T * c.inv()) == -g,
                'Charge-conjugation gamma intertwiner failed')

    empty: Poly = {}
    left_bar = [empty, empty, variable('lb0'), variable('lb1')]
    # F = (chi_alpha, bar(xi)^dot(alpha)); bar F = (xi^alpha, bar chi_dot(alpha)).
    f_bar_t = [variable('xi0'), variable('xi1'),
               variable('chib0'), variable('chib1')]
    matchete = components_bilinear(left_bar, c, f_bar_t)
    lambda_chi = add(product(variable('l0'), variable('chi1')),
                     scale(sp.Integer(-1), product(variable('l1'), variable('chi0'))))
    weyl_hc = hc(lambda_chi)
    require(matchete == weyl_hc,
            'Canonical Weyl identity NCM[Bar L,CConj F] = (l eps chi)^dagger failed')
    require(matchete and len(matchete) == 2,
            'Unexpected zero or incomplete charge-conjugate bilinear')
    require(components_bilinear(left_bar, -c, f_bar_t) != weyl_hc,
            'Negative control failed to distinguish C versus -C')
    require(hc(weyl_hc) == lambda_chi,
            'Hermitian conjugation involution failed')

    return {
        'clifford_identity_verified': True,
        'charge_conjugation_matrix_4_by_4': [[str(v) for v in row] for row in c.tolist()],
        'cc_antisymmetry_and_square_verified': True,
        'cc_gamma_intertwiner_verified': True,
        'weyl_spinor_epsilon': [[int(v) for v in row] for row in eps.tolist()],
        'four_component_to_weyl_bilinear_exact': True,
        'identity_in_declared_weyl_basis':
            'Bar[L] C (Bar[F])^T = (l_0 chi_1 - l_1 chi_0)^dagger',
        'bilinear_grassmann_monomials': len(matchete),
        'explicit_cc_sign_choice': 'C = +i gamma^2 gamma^0',
        'opposite_cc_sign_negative_control_rejected': True,
        'weyl_bilinear_eta_in_chosen_basis': '+1',
        'matchete_abstract_cc_sign_is_native_fixed': False,
        'rgbeta_abstract_eps_12_plus_one_is_native_fixed': False,
        'conditional_y1R_over_y1M_in_chosen_basis_with_rho_plus_one': '-1',
        'conditional_y1ConjR_over_y1ConjM_in_chosen_basis_with_rho_plus_one': '-1',
    }


def check_prior(root: Path) -> None:
    p = root / PRIOR
    require(p.is_file(), 'Missing passed A(0) action witness: ' + str(p))
    r = json.loads(p.read_text(encoding='utf-8'))
    require(r.get('status') ==
            'A0ExactFormalWeylActionIdentityBothIndependentYukawasPhysicalPhasePending'
            and r.get('formal_action_identity_under_declared_assumptions') is True
            and r.get('absolute_matchete_rgbeta_yukawa_coupling_mapping_certified') is False
            and r.get('model') == 'A_p0',
            'Existing A(0) conditional action witness missing or changed')
    candidates = r.get('conditional_candidate_maps', [])
    require(isinstance(candidates, list) and len(candidates) == 2 and
            {i.get('epsilon_relative_orientation_rho') for i in candidates} == {-1, 1},
            'Previous A(0) epsilon alternatives missing')


def check_native(p: Path) -> dict:
    require(p.is_file(), 'Native Matchete report missing: ' + str(p))
    r = json.loads(p.read_text(encoding='utf-8'))
    require(r.get('status') == NATIVE_STATUS and r.get('model') == 'A_p0'
            and r.get('full_A0_matchete_build_passed') is True,
            'Matchete native spinor convention test did not pass')
    c = r.get('checks', {})
    require(isinstance(c, dict) and all(c.get(k) is True for k in NATIVE_KEYS),
            'One or more native CConj/H.c. algebra checks failed')
    require(r.get('native_concrete_four_by_four_gamma_cc_sign_measured') is False
            and r.get('rgbeta_epsilon_component_orientation_measured') is False
            and r.get('physical_cross_package_coupling_map_certified') is False,
            'Native JSON overclaims physical phase or gauge-basis certification')
    return {'native_checks_passed': len(NATIVE_KEYS),
            'native_charge_conjugate_bilinear':
                r.get('matchete_charge_conjugated_bilinear_inputform', '')}


def run(root: Path) -> dict:
    check_prior(root)
    result = derive_weyl_identity()
    output = root / OUT
    output.mkdir(parents=True, exist_ok=True)
    native_path = root / NATIVE
    native_path.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update({'T3_EXPERIMENTAL_GENERAL_FERMIONS': '1',
                'T3_EXPERIMENTAL_FLAVOR_MASSES': '1',
                'T3_EXPERIMENTAL_GENERAL_QUARTICS': '1',
                'T3_EXPERIMENTAL_GENERAL_QUADRATICS': '0'})
    # Prevent a stale success report from a previous run being accepted.
    if native_path.is_file():
        native_path.unlink()
    cmd = ['wolframscript', '-file', str(root/DRIVER), str(native_path)]
    try:
        proc = subprocess.run(cmd, cwd=root, env=env, capture_output=True,
                              text=True, errors='replace', check=False)
    except OSError as exc:
        raise PhaseError('Cannot start wolframscript: ' + str(exc)) from exc
    (output / 'native_matchete/A_p0.log').write_text(
        proc.stdout + '\nSTDERR:\n' + proc.stderr, encoding='utf-8')
    if proc.returncode != 0:
        diagnostic = ''
        if native_path.is_file():
            try:
                obj = json.loads(native_path.read_text(encoding='utf-8'))
                diagnostic = str(obj.get('failed_step', '')) + ': ' + str(obj.get('reason', ''))
            except (OSError, ValueError):
                pass
        raise PhaseError(f'Native Matchete spinor probe failed (exit {proc.returncode})'
                         + (': '+diagnostic if diagnostic else '')
                         + '; see ' + str(output / 'native_matchete/A_p0.log'))
    proof = check_native(native_path)
    report = {
        'format': 't3-workstream-a-a0-spinor-phase-witness-v1',
        'status': STATUS,
        'model': 'A_p0',
        'matchete_native_abstract_cc_algebra_verified': True,
        'matchete_native': proof,
        'exact_weyl_component_witness': result,
        'matched_existing_a0_action_dictionary': True,
        'coupling_dictionary_in_explicitly_chosen_CC_and_epsilon_basis': {
            'y1_R': '-y1_M',
            'y1Conj_R': '-y1Conj_M',
            'depends_on': ['C = +i gamma2 gamma0', 'eps_12 = +1',
                           'identical phase choice for the left-Weyl heavy field']
        },
        'matchete_installed_concrete_CC_orientation_proved': False,
        'rgbeta_installed_numeric_epsilon_orientation_proved': False,
        'global_dirac_mass_and_y2_phase_transport_proved': False,
        'physical_cross_package_yukawa_dictionary_certified': False,
        'production_files_modified': False,
        'next_gate': 'Transport one common heavy-fermion/Weyl phase through y2 and MF and fix RGBeta eps component orientation',
    }
    dest = output / 'a0_spinor_phase.json'
    dest.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return {'report': str(dest), 'result': report}


def self_test() -> None:
    exact = derive_weyl_identity()
    require(exact['weyl_bilinear_eta_in_chosen_basis'] == '+1',
            'Incorrect canonical Weyl eta')
    # Direct anticommution and H.c. are covered by the exact matrix witness.
    require(product(variable('l0'), variable('chi0')) ==
            scale(-1, product(variable('chi0'), variable('l0'))),
            'Grassmann interchange did not change sign')
    require(product(variable('l0'), variable('l0')) == {},
            'Grassmann variable did not square to zero')
    print('PASS: exact gamma algebra, C = +i gamma2 gamma0, Grassmann H.c., eta=+1 in chosen Weyl basis')


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--repo-root', type=Path, default=ROOT)
    a = ap.parse_args()
    try:
        if a.self_test:
            self_test()
            return 0
        require(a.all, 'Specify --self-test or --all')
        r = run(a.repo_root.resolve())
        print(r['result']['status'])
        print('Native Matchete abstract spinor checks: 6/6 PASS')
        print('Explicit Weyl-basis C choice: eta=+1; y1_R=-y1_M, y1Conj_R=-y1Conj_M')
        print('Installed absolute C sign / RGBeta eps axes / MF and y2 phase: NOT CERTIFIED')
        print('Report: ' + r['report'])
        return 0
    except (PhaseError, OSError, ValueError, TypeError, KeyError) as exc:
        print('A(0) native spinor phase bridge BLOCKED: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
