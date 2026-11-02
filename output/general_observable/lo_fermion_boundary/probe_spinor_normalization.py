"""Independent four-component/Weyl algebra check; never certifies unknown Matchete GammaCC.

Install as output/general_observable/lo_fermion_boundary/probe_spinor_normalization.py
Does not modify the 16-model boundary or any source repository file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

FORMAT = 't3-gammacc-weyl-matrix-probe-v1'


def matrices() -> dict:
    I = sp.I
    two = sp.eye(2)
    zero = sp.zeros(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -I], [I, 0]])
    sz = sp.diag(1, -1)
    def blocks(a, b, c, d):
        return sp.BlockMatrix([[a, b], [c, d]]).as_explicit()
    g0 = blocks(zero, two, two, zero)
    gspatial = [blocks(zero, s, -s, zero) for s in (sx, sy, sz)]
    g5 = sp.diag(-1, -1, 1, 1)
    PL = (sp.eye(4) - g5) / 2
    PR = (sp.eye(4) + g5) / 2
    C = sp.I * gspatial[1] * g0
    eps = sp.Matrix([[0, 1], [-1, 0]])
    return dict(gamma=[g0] + gspatial, gamma5=g5, PL=PL, PR=PR,
                C=C, C_inverse=C.inv(), epsilon=eps)


def check_exact_algebra() -> dict:
    m = matrices()
    eye = sp.eye(4)
    gamma = m['gamma']
    eta = [1, -1, -1, -1]
    for mu in range(4):
        for nu in range(4):
            target = 2 * eta[mu] * eye if mu == nu else sp.zeros(4)
            if gamma[mu] * gamma[nu] + gamma[nu] * gamma[mu] != target:
                raise AssertionError('Chiral Dirac matrices violate Clifford algebra')
    C = m['C']
    if C.T != -C or C.inv() * C != eye or C.inv() != -C:
        raise AssertionError('Charge-conjugation identities failed')
    for g in gamma:
        if C.inv() * g * C != -g.T:
            raise AssertionError('Charge-conjugation transformation failed')
    eps = m['epsilon']
    if (C * m['PL'])[:2, :2] != eps or (C * m['PR'])[:2, :2] != sp.zeros(2):
        raise AssertionError('Expected left-handed charge-conjugation block failed')
    chi1, chi2, eta1, eta2 = sp.symbols('chi1 chi2 eta1 eta2')
    left = sp.Matrix([chi1, chi2, 0, 0])
    right = sp.Matrix([eta1, eta2, 0, 0])
    two_component = (sp.Matrix([chi1, chi2]).T * eps * sp.Matrix([eta1, eta2]))[0]
    in_C = (left.T * C * m['PL'] * right)[0]
    in_inverse = (left.T * C.inv() * m['PL'] * right)[0]
    if sp.simplify(in_C - two_component) != 0 or sp.simplify(in_inverse + two_component) != 0:
        raise AssertionError('Four-component to two-component bilinear reduction failed')
    if sp.simplify(two_component - (chi1 * eta2 - chi2 * eta1)) != 0:
        raise AssertionError('Two-component epsilon convention changed')
    return {
        'basis': 'chiral/Weyl gamma matrices, metric diag(+---)',
        'chosen_two_component_epsilon': [[0, 1], [-1, 0]],
        'candidate_if_GammaCC_is_C_eq_i_gamma2_gamma0': '+1',
        'candidate_if_GammaCC_is_C_inverse': '-1',
        'candidate_if_GammaCC_is_minus_C': '-1',
        'independent_matrix_identities_pass': True,
        'caution': ('These are mathematical conditional identities. They do NOT establish which '
                    'matrix Matchete GammaCC actually represents or its ordering/Grassmann convention.'),
    }


def check_master_orbit_normalization() -> dict:
    # The source master uses 1/4 times the fourth formal derivative, including
    # the two fermionic and two scalar derivatives, with commuting placeholders.
    ps1, ps2, ph1, ph2, k = sp.symbols('ps1 ps2 ph1 ph2 k')
    cases = [(k*ps1**2*ph1**2, ps1, ps1, ph1, ph1, k, 1),
             (k*ps1**2*ph1*ph2, ps1, ps1, ph1, ph2, k/2, 2),
             (k*ps1*ps2*ph1*ph2, ps1, ps2, ph1, ph2, k/4, 4)]
    for p, a, b, c, d, expected, multiplicity in cases:
        factor = sp.diff(p, a, b, c, d) / 4
        if sp.simplify(factor-expected) != 0 or sp.simplify(factor*multiplicity-k) != 0:
            raise AssertionError('Master ordered-coordinate orbit convention failed')
    return {'derivative_rule': 'C_ijab = (1/4) d_i d_j d_a d_b P',
            'distinct_orbit_sizes_verified': [1, 2, 4],
            'no_additional_orbit_multiplication': True,
            'caution': 'Formal polynomial normalization, not an external Matchete spinor convention.'}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report(repo: Path, evidence: Path | None) -> dict:
    file_names = ('Lagrangian/RunMatching.wl',
                  'RGE/matching/WeinbergTensorAdapter.py',
                  'RGE/general/WilsonTensorRGE.py')
    sources = {}
    for name in file_names:
        f = repo / name
        if not f.is_file():
            raise FileNotFoundError(f'Missing local source: {f}')
        sources[name] = sha(f)
    matching = (repo / file_names[0]).read_text(encoding='utf-8-sig')
    adapter = (repo / file_names[1]).read_text(encoding='utf-8-sig')
    if '"Proj[-1]"], "PL"' not in matching or '"Proj[1]"], "PR"' not in matching:
        raise ValueError('Project Matchete PL/PR source convention changed; re-audit')
    if 'sp.diff(' not in adapter or 'derivative)' not in adapter:
        raise ValueError('Master real-basis adapter changed; re-audit')
    e = None
    if evidence is not None:
        if not evidence.is_file():
            raise FileNotFoundError(f'Missing Wolfram probe evidence: {evidence}')
        e = {'file': str(evidence.resolve()), 'sha256': sha(evidence),
             'bytes': evidence.stat().st_size,
             'head': evidence.read_text(encoding='utf-8-sig', errors='replace')[:12000]}
    return {
        'format': FORMAT,
        'status': 'ExactAlgebraConditionalOnly_ActualMatcheteGammaCCUndefined',
        'repo_root': str(repo), 'local_repository_sha256': sources,
        'matchete_definition_evidence': e,
        'matrix_algebra': check_exact_algebra(),
        'master_polynomial_algebra': check_master_orbit_normalization(),
        'physical_absolute_normalization': {
            'certified': False, 'value': None,
            'reason': 'No independently verified Matchete GammaCC matrix AND Weyl epsilon phase convention',
        },
        'ready_for_physical_eft1_running': False,
        'tree_f_matching_only': True,
    }


def save_new(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, allow_nan=False)
        f.write('\n')


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[3])
    p.add_argument('--matchete-evidence', type=Path,
                   help='Wolfram console evidence captured using probe_matchete_gamma_cc.wl')
    p.add_argument('--check-only', action='store_true')
    args = p.parse_args(argv)
    repo = args.repo_root.resolve()
    try:
        result = report(repo, args.matchete_evidence)
        if args.check_only:
            print(json.dumps(result, indent=2))
            return 0
        base = repo / 'output/general_observable/lo_fermion_boundary/normalization_probes'
        key = hashlib.sha256(json.dumps({'sources': result['local_repository_sha256'],
              'evidence': result['matchete_definition_evidence'] and
                  result['matchete_definition_evidence']['sha256']}, sort_keys=True).encode()).hexdigest()[:16]
        path = base / f'probe_{key}.json'
        save_new(path, result)
        print(f'CREATED: {path}')
        print('STATUS: conditional Clifford/Weyl reduction only; physical sign/factor UNRESOLVED')
        return 0
    except (OSError, KeyError, ValueError, AssertionError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2

if __name__ == '__main__':
    raise SystemExit(main())
