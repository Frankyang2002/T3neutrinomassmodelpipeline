"""Exact Grassmann p=0 heavy-fermion exchange vs LLSS contact diagnostic.

Project path: output/general_observable/lo_fermion_boundary/signed_zero_momentum_check.py

The UV action is an explicitly DECLARED two-component convention, not a
Matchete-to-Weyl identity. Cross-check selected *native* source term signs and
1/2 factors using the current verified release and SHA-bound local inputs.
Never promotes a physical Wilson normalization.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

# Self-contained: script-local import, not a change to production/RGE code.
import build_boundary as bound

PATH = Path(__file__).resolve()
ROOT = PATH.parents[3]
OUT = PATH.parent
RELEASE = 'source_021dcc1e8fe8ba0f'
FORMAT = 't3-signed-zero-momentum-grassmann-diagnostic-v1'


class G:
    """Finite exterior algebra over Q, canonical ordered Grassmann generators."""
    def __init__(self, monomials=None):
        self.c = {tuple(k): Q(v) for k, v in (monomials or {}).items() if v}
        for k in self.c:
            if tuple(sorted(k)) != k or len(set(k)) != len(k):
                raise ValueError('Noncanonical exterior monomial')

    @staticmethod
    def var(name):
        return G({(name,): 1})

    @staticmethod
    def scalar(n):
        return G({(): Q(n)})

    def __add__(self, other):
        if not isinstance(other, G):
            other = G.scalar(other)
        d = self.c.copy()
        for k, v in other.c.items():
            d[k] = d.get(k, Q(0)) + v
        return G(d)

    __radd__ = __add__

    def __neg__(self):
        return G({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + -other if isinstance(other, G) else self + -G.scalar(other)

    def __rsub__(self, other):
        return G.scalar(other) - self

    def __mul__(self, other):
        if not isinstance(other, G):
            other = G.scalar(other)
        d = {}
        for a, av in self.c.items():
            for b, bv in other.c.items():
                if set(a).intersection(b):
                    continue
                inversions = sum(x > y for x in a for y in b)
                key = tuple(sorted(a + b))
                d[key] = d.get(key, Q(0)) + av * bv * (-1 if inversions % 2 else 1)
        return G(d)

    __rmul__ = __mul__

    def __truediv__(self, n):
        return self * (Q(1) / Q(n))

    def __eq__(self, other):
        if not isinstance(other, G):
            other = G.scalar(other)
        return self.c == other.c

    def left_derivative(self, v):
        """∂^L_v(v w) = w; derivative sign follows Grassmann ordering."""
        d = {}
        for key, coeff in self.c.items():
            if v in key:
                j = key.index(v)
                rest = key[:j] + key[j + 1:]
                d[rest] = d.get(rest, Q(0)) + coeff * (-1 if j % 2 else 1)
        return G(d)

    def substitute(self, repl):
        out = G()
        for key, coeff in self.c.items():
            monomial = G.scalar(coeff)
            for v in key:
                monomial = monomial * repl.get(v, G.var(v))
            out = out + monomial
        return out


def bil(u, v):
    """u^T epsilon v, epsilon_12=+1; u,v contain anticommuting entries."""
    return u[0] * v[1] - u[1] * v[0]


def linear_solve(mat, rhs):
    """Exact rational Gaussian elimination; rhs entries are Grassmann odd polynomials."""
    n = len(mat)
    a = [[Q(x) for x in row] + [rhs[i]] for i, row in enumerate(mat)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            raise ValueError('Singular heavy mass kernel')
        a[k], a[pivot] = a[pivot], a[k]
        p = a[k][k]
        for j in range(k, n):
            a[k][j] /= p
        a[k][n] /= p
        for i in range(n):
            if i == k:
                continue
            f = a[i][k]
            if f:
                for j in range(k, n):
                    a[i][j] -= f * a[k][j]
                a[i][n] -= f * a[k][n]
    return [a[i][n] for i in range(n)]


def integrate_out(lagrangian, heavy):
    """Solve exact left-Euler equations and substitute back into UV action."""
    eoms = [lagrangian.left_derivative(v) for v in heavy]
    all_zero = {h: G() for h in heavy}
    source = [e.substitute(all_zero) for e in eoms]
    mat = []
    for e in eoms:
        row = []
        for h in heavy:
            coefficient = e.c.get((h,), Q(0))
            row.append(coefficient)
        residual = e - source[len(mat)] - sum((row[j] * G.var(h) for j, h in enumerate(heavy)), G())
        if residual != G():
            raise ValueError('UV EOM is not linear in the heavy fields')
        mat.append(row)
    sol = linear_solve(mat, [-x for x in source])
    replacements = dict(zip(heavy, sol))
    if any(e.substitute(replacements) != G() for e in eoms):
        raise AssertionError('Heavy equations of motion not satisfied')
    return lagrangian.substitute(replacements), replacements


def amplitude_benchmarks():
    """Tree-level zero-external-momentum UV action; exact signed contact tests."""
    a = [G.var('a1'), G.var('a2')]
    b = [G.var('b1'), G.var('b2')]
    h = [G.var('h1'), G.var('h2')]
    u = [G.var('u1'), G.var('u2')]
    v = [G.var('v1'), G.var('v2')]
    m, ya, yb = Q(7), Q(3), Q(-2)
    # Form of the *assumed* UV mass/Yukawa terms follows the explicitly
    # declared two-component epsilon convention. This is not an extraction
    # of Matchete's field-phase or global GammaCC normalization.
    uv_maj = -m * bil(h, h) / 2 + ya * bil(h, a) + yb * bil(h, b)
    eff_maj, repl_maj = integrate_out(uv_maj, ['h1', 'h2'])
    exp_maj = (ya * ya * bil(a, a) + 2 * ya * yb * bil(a, b) + yb * yb * bil(b, b)) / (2 * m)
    assert eff_maj == exp_maj

    uv_dir = -m * bil(u, v) + ya * bil(u, a) + yb * bil(v, b)
    eff_dir, repl_dir = integrate_out(uv_dir, ['u1', 'u2', 'v1', 'v2'])
    exp_dir = ya * yb * bil(a, b) / m
    assert eff_dir == exp_dir
    assert bil(a, b) == bil(b, a)
    assert bil(a, a) == 2 * G.var('a1') * G.var('a2')
    assert bil(a, b).c[('a1', 'b2')] == 1
    assert bil(a, b).c[('a2', 'b1')] == -1

    return {
        'convention': {
            'epsilon_12': 1, 'epsilon_21': -1,
            'uv_majorana_mass': '-(M/2) h^T epsilon h',
            'uv_dirac_mass': '-M u^T epsilon v',
            'uv_yukawa': '+y_a heavy^T epsilon a + y_b heavy^T epsilon b',
            'Lagrangian_not_amplitude_i_factor': True,
        },
        'benchmark_inputs': {'M': int(m), 'y_a': int(ya), 'y_b': int(yb)},
        'majorana_contact': '+(y_a^2 aa + 2 y_a y_b ab + y_b^2 bb)/(2M)',
        'dirac_contact': '+y_a y_b ab/M',
        'zero_momentum_Grassmann_saddle_exact': True,
        'majorana_same_channel_symmetry_factor': '1/2',
        'majorana_cross_channel_symmetry_factor': '1',
        'dirac_cross_channel_factor': '1',
        'spinor_exchange_signs_verified': True,
        'fully_physical_matchete_to_master_normalization_certified': False,
    }


def _require(term, *, index, factor_names, ending, negative=False, chirality='PL'):
    actual = term['source_term_inputform']
    assert term['source_index'] == index
    assert [t['name'] for t in term['factors']] == factor_names
    assert actual.rstrip().endswith(ending), (index, actual[-150:])
    assert actual.lstrip().startswith('-') == negative
    assert bound.chirality(actual) == chirality
    return {'source_index': index, 'factor_names': factor_names, 'negative_CG_prefactor': negative,
            'ending_verified': ending, 'projector': chirality}


def source_benchmarks(root, release):
    base = root/'output/general_observable'
    rows = {}
    source_sha = {}
    summary = bound.read(release/'summary.json')
    assert summary['n_models_verified'] == 16
    for model in ('A_-4','B_-1','C_-1'):
        oroot, proot, aroot = [base/name for name in
            ('ordered_wilson_tensor_v4','general_wilson_polynomial_batch_v2',
             'general_wilson_normalization_audit_v1')]
        regenerated, _ = bound.compile_model(model, root, oroot, proot, aroot)
        output = bound.read(release/bound.tag(model)/'boundary.json')
        assert output['checksums_sha256'] == regenerated['checksums_sha256']
        ordered = bound.read(oroot/bound.tag(model)/'ordered_tensor.json')
        native_path = bound.source_file(model, root, oroot, ordered)
        native = bound.read(native_path)
        assert bound.digest(native_path) == output['checksums_sha256']['native_source']
        source_sha[model] = bound.digest(native_path)
        terms = {t['source_index']:t for t in native['operator_terms']}
        rows[model] = terms
    a,b,c=rows['A_-4'],rows['B_-1'],rows['C_-1']
    tests = {
       'A_-4_Dirac_cross': _require(a[2],index=2,factor_names=['y1','y2'],
                     ending='/Coupling[MF, {Index[d$$1, NFlavor]}, 0]'),
       'B_-1_Majorana_diagonal': _require(b[12],index=12,factor_names=['y1','y1'],
                     ending='/(2*Coupling[MF, {}, 0])'),
       'B_-1_Majorana_cross': _require(b[14],index=14,factor_names=['y1','y2'],
                     ending='/Coupling[MF, {}, 0]'),
       'C_-1_Majorana_diagonal': _require(c[13],index=13,factor_names=['y1','y1'],
                     ending='/(2*Sqrt[3]*Coupling[MF, {}, 0])'),
       'C_-1_Majorana_cross_CG_minus': _require(c[14],index=14,factor_names=['y1','y2'],
                     ending='/(Sqrt[3]*Coupling[MF, {}, 0])',negative=True),
    }
    assert sum(bound.chirality(x['source_term_inputform'])=='PL' for x in b.values()) == 10
    assert sum(bound.chirality(x['source_term_inputform'])=='PR' for x in b.values()) == 10
    assert sum(bound.chirality(x['source_term_inputform'])=='PL' for x in c.values()) == 12
    assert sum(bound.chirality(x['source_term_inputform'])=='PR' for x in c.values()) == 12
    return {'source_checks': tests, 'source_sha256': source_sha,
            'matchete_source_checks_passed': len(tests),
            'native_CG_signs_are_not_spinor_signs': True,
            'mixed_PR_conjugates_excluded': True}


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--offline',action='store_true',help='Exact Grassmann tests; no repository input')
    ap.add_argument('--repo-root',type=Path,default=ROOT)
    ap.add_argument('--release',type=Path,help='Previously generated source_<hash> release directory')
    ap.add_argument('--output-dir',type=Path,help='New diagnostic directory (not a release)')
    args=ap.parse_args(argv)
    try:
        result={'format':FORMAT,'status':'ConditionalSignedContactOnly_NoPhysicalCertification',
                'grassmann_zero_momentum':amplitude_benchmarks(),
                'physical_master_C_ijab_normalization': None,
                'assumed_GammaCC_phase_certified': False,
                'momentum_dependent_UV_Feynman_amplitude_computed': False,
                'native_matchete_UV_yukawa_mass_action_independently_matched': False,
                'does_not_mutate_previous_boundary': True}
        if args.offline:
            print('PASS: exact signed Dirac/Majorana p=0 Grassmann UV-to-contact algebra')
            print('NO PHYSICAL CERTIFICATION: Matchete GammaCC phase and master Weyl identification open')
            return 0
        root=args.repo_root.resolve()
        release=(args.release or (root/'output/general_observable/lo_fermion_boundary/releases'/RELEASE)).resolve()
        assert release.is_dir(),f'Missing 16-model boundary release: {release}'
        result['local_native_source_checks'] = source_benchmarks(root,release)
        result['input_release_summary_sha256'] = bound.digest(release/'summary.json')
        encoded=json.dumps(result,sort_keys=True,indent=2)+'\n'
        hashkey=hashlib.sha256(encoded.encode()).hexdigest()[:16]
        dest=(args.output_dir or (OUT/'normalization_probes')).resolve()
        assert dest != release and release not in dest.parents
        target=dest/f'signed_zero_momentum_{hashkey}.json'
        if target.exists():
            print(f'UNCHANGED: {target} (same source snapshot)')
        else:
            dest.mkdir(parents=True,exist_ok=True)
            with target.open('x',encoding='utf-8') as f:
                f.write(encoded)
            print('CREATED:',target)
        print('PASS: exact Grassmann benchmark + 5 native-source sign/symmetry-factor checks')
        print('NOT CERTIFIED: absolute GammaCC/Weyl phase; independent native Feynman amplitude still required')
        return 0
    except (AssertionError,ValueError,KeyError,TypeError,OSError) as e:
        print(f'BLOCKED: {type(e).__name__}: {e}',file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
