"""Signed zero-momentum Dirac-F exchange check for native Matchete A_-4.

Project path: output/general_observable/lo_fermion_boundary/native_zero_momentum_a_m4.py

The UV algebra uses the actual native A_-4 interaction ordering:
  + y1* L^T C P_L F  S1*  + y2* Fbar P_L L S2  - M Fbar F.
The heavy F,Fbar Wick contraction, both vertex time orderings and both
external L orderings are computed with a finite exterior (Grassmann) algebra.
The native EFT vertex is separately authenticated from an existing on-disk
FeynmanRules report; the CG tensor contraction and the absolute master-Weyl
phase are NOT computed or certified here.

This program never edits the certified 16-model boundary.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
DEFAULT_REPORT = HERE / 'normalization_probes' / (
    'native_feynman_a_m4_45f6bbfe96e54fc79eb25362cb3e9d92.json'
)


class Exterior:
    """Small exact Grassmann polynomial with ordered, distinct generators."""

    def __init__(self, terms=None):
        self.terms = {tuple(k): sp.expand(v) for k, v in (terms or {}).items()
                      if sp.expand(v) != 0}

    @staticmethod
    def generator(index):
        return Exterior({(index,): sp.Integer(1)})

    def __add__(self, other):
        if not isinstance(other, Exterior):
            return NotImplemented
        terms = dict(self.terms)
        for k, v in other.terms.items():
            terms[k] = sp.expand(terms.get(k, 0) + v)
        return Exterior(terms)

    def __radd__(self, other):
        if other == 0:
            return self
        return NotImplemented

    def __neg__(self):
        return Exterior({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        if not isinstance(other, Exterior):
            return self.scale(other)
        terms = {}
        for left, a in self.terms.items():
            for right, b in other.terms.items():
                if len(set(left) | set(right)) != len(left) + len(right):
                    continue
                swaps = sum(i > j for i in left for j in right)
                seq = tuple(sorted(left + right))
                terms[seq] = sp.expand(terms.get(seq, 0) + (-1)**swaps * a*b)
        return Exterior(terms)

    def __rmul__(self, other):
        return self.scale(other)

    def scale(self, scalar):
        return Exterior({k: sp.expand(v*scalar) for k, v in self.terms.items()})

    def equal(self, other):
        return isinstance(other, Exterior) and not (self-other).terms

    def coefficients(self):
        return {','.join(map(str, k)): str(v) for k, v in sorted(self.terms.items())}


E = sp.Matrix([[0, 1], [-1, 0]])  # epsilon = i sigma_2
C = sp.diag(E, -E)
PL = sp.diag(1, 1, 0, 0)
M = sp.Symbol('M', positive=True)
Y1I, Y2J, Y1J, Y2I = sp.symbols('Y1i Y2j Y1j Y2i')
# Generator ordering: F[0..3], Fbar[4..7], Li[8..11], Lj[12..15].
F = [Exterior.generator(a) for a in range(4)]
FB = [Exterior.generator(a + 4) for a in range(4)]
LI = [Exterior.generator(a + 8) for a in range(4)]
LJ = [Exterior.generator(a + 12) for a in range(4)]


def spin_bilinear(left, right, c=C):
    mat = c * PL
    return sum((left[a]*right[b]).scale(mat[a, b])
               for a in range(4) for b in range(4) if mat[a, b] != 0)


def yukawa_one(lepton, coupling, c=C):
    """The native Bar[y1] L^T C P_L F orientation, not the EFT operator."""
    mat = c*PL
    return sum((lepton[a]*F[b]).scale(coupling*mat[a,b])
               for a in range(4) for b in range(4) if mat[a,b] != 0)


def yukawa_two(lepton, coupling):
    """The native Bar[y2] Fbar P_L L orientation."""
    return sum((FB[a]*lepton[b]).scale(coupling*PL[a,b])
               for a in range(4) for b in range(4) if PL[a,b] != 0)


def contract_free_dirac_at_zero(expr):
    """Compute <F_a Fbar_b> = delta_ab/M via the normalized Gaussian.

    With Berezin order d(Fbar)dF and action -M Fbar F:
       <F Fbar> = [integral F Fbar(1-M Fbar F)]/(-M) = +1/M.
    The generator sorting already includes all Wick crossing signs.
    Only terms with one F and one Fbar are accepted.
    """
    terms = {}
    for monomial, coeff in expr.terms.items():
        fs = [a for a in monomial if 0 <= a < 4]
        bars = [a for a in monomial if 4 <= a < 8]
        if len(fs) != 1 or len(bars) != 1:
            raise ValueError('Heavy Gaussian contraction expects exactly F and Fbar')
        if fs[0] != bars[0] - 4:
            continue
        remaining = tuple(a for a in monomial if a >= 8)
        terms[remaining] = sp.expand(terms.get(remaining, 0) + coeff/M)
    return Exterior(terms)


def connected_exchange(c=C):
    """Second order in S_int: 1/2(V1 V2 + V2 V1), including swaps."""
    v1i, v2j = yukawa_one(LI, Y1I, c), yukawa_two(LJ, Y2J)
    v1j, v2i = yukawa_one(LJ, Y1J, c), yukawa_two(LI, Y2I)
    ij = contract_free_dirac_at_zero((v1i*v2j + v2j*v1i).scale(sp.Rational(1,2)))
    ji = contract_free_dirac_at_zero((v1j*v2i + v2i*v1j).scale(sp.Rational(1,2)))
    return ij, ji


def exact_check():
    assert C.T == -C and C*PL == PL*C
    assert E.T == -E and E*E == -sp.eye(2)
    propagator_zero = sp.simplify(sp.I*M*sp.eye(4)/(0-M**2))
    assert propagator_zero == -sp.I/M*sp.eye(4)
    assert sp.I**2*(-sp.I) == sp.I  # i^2 vertices, -i/M propagator
    ij, ji = connected_exchange()
    expected_ij = spin_bilinear(LI,LJ).scale(Y1I*Y2J/M)
    expected_ji = spin_bilinear(LJ,LI).scale(Y1J*Y2I/M)
    assert ij.equal(expected_ij), 'UV ij exchange disagrees with contact sign'
    assert ji.equal(expected_ji), 'UV ji exchange disagrees with contact sign'
    assert spin_bilinear(LI,LJ).equal(spin_bilinear(LJ,LI))
    # A phase change C -> -C changes UV and native contact together: the
    # relative amplitude passes, but this cannot fix the physical Weyl phase.
    reverse_ij, reverse_ji = connected_exchange(-C)
    assert reverse_ij.equal(-ij) and reverse_ji.equal(-ji)
    return {
        'uv_ij_equals_native_contact_ij_exact': True,
        'uv_ji_equals_native_contact_ji_exact': True,
        'external_fermion_antisymmetrization_exact': True,
        'two_vertex_time_orderings_included': True,
        'exchange_kernel_uv_over_native_contact': '1',
        'independent_gamma_cc_phase_sign_still_unresolved': True,
        'spinor_basis': 'C=diag(i sigma2,-i sigma2), epsilon12=+1, PL=diag(1,1,0,0)',
        'connected_kernel_coefficients': (ij+ji).coefficients(),
        'heavy_dirac_propagator_at_p_zero': '-i/M times identity_4',
        'minkowski_vertex_and_propagator_factor': '(i)^2 * (-i) = +i, before full external fermion and CG factors',
        'amplitude_i_factor': '+i times connected local kernel',
    }


def required_report(rec):
    if rec.get('format') != 't3-native-matchete-fr-evidence-v1' or rec.get('model') != 'A_-4':
        raise ValueError('Wrong A_-4 report schema/model')
    if rec.get('status') != 'NativeUV3pointAndEFT1ContactVerticesExtracted_AmplitudeOpen':
        raise ValueError('Native report is not a successful vertex extraction')
    for name in ('uv_three_point_rules', 'uv_yukawa_lagrangian',
                 'eft1_four_point_rules', 'uv_mass_lagrangian_source',
                 'original_holomorphic_contact'):
        val = rec.get(name)
        if not isinstance(val, dict) or val.get('inputform_complete') is not True:
            raise ValueError(f'Native {name} inputform missing or truncated')
        s = val.get('inputform')
        if not isinstance(s, str) or hashlib.sha256(s.encode('utf8')).hexdigest() != val.get('inputform_sha256'):
            raise ValueError(f'Native {name} SHA256 failed')
        if len(s) != val.get('inputform_characters'):
            raise ValueError(f'Native {name} length mismatch')
    src = rec['original_holomorphic_contact']['inputform']
    if src != rec.get('native_source_inputform'):
        raise ValueError('Native matching source changed since FR extraction')
    if not all(x in src for x in ['NCM[Transp[Field[l,', 'GammaCC, Proj[-1]', 'Coupling[MF,']):
        raise ValueError('Holomorphic source chirality/mass inconsistent')
    mass = rec['uv_mass_lagrangian_source']['inputform']
    if not ('-(Coupling[MF,' in mass and
            'NCM[Bar[Field[NewFermion' in mass and
            'Field[NewFermion' in mass):
        raise ValueError('Heavy Dirac mass -M Fbar F not established')
    uv = rec['uv_yukawa_lagrangian']['inputform']
    if not all(x in uv for x in [
        'Bar[Coupling[y1,', 'Bar[Coupling[y2,',
        'NCM[Transp[Field[l,', 'DiracProduct[GammaCC, Proj[-1]]',
        'NCM[Bar[Field[NewFermion', 'DiracProduct[Proj[-1]]',
    ]):
        raise ValueError('UV Lagrangian cannot support required selected two-vertex Wick contraction')
    vertices = rec['uv_three_point_rules']['inputform']
    if ('{Conj[NewScalar1], NewFermion, l} -> I*Bar[Coupling[y1,' not in vertices or
        '{NewScalar2, l, Conj[NewFermion]} -> (-I)*Bar[Coupling[y2,' not in vertices):
        raise ValueError('UV Feynman-rule relative signs/flow differ from pinned native conventions')
    eft = rec['eft1_four_point_rules']['inputform']
    marker = '{Conj[NewScalar1], NewScalar2, l, l} ->'
    if marker not in eft:
        raise ValueError('Holomorphic contact vertex missing')
    holomorphic = eft.split(marker,1)[1].split('|>',1)[0]
    if not holomorphic.count('I*Bar[Coupling[y1,') == 2:
        raise ValueError('Both native EFT lepton orderings with +i are not identifiable')
    if 'NCM[Transp[Field[Ext[l,' not in holomorphic or 'GammaCC, Proj[-1]' not in holomorphic:
        raise ValueError('Native EFT spinor flow/PL projector invalid')
    if any(rec.get(name) is not False for name in (
            'physical_master_wilson_normalization_certified',
            'independent_connected_uv_four_point_amplitude_computed',
            'matchete_gamma_cc_to_weyl_phase_certified')):
        raise ValueError('Input report has unexpected physical certification')
    return {
        'native_source_index': rec['native_source_index'],
        'uv_vertex_signs': {'y1_bar': '+i', 'y2_bar': '-i'},
        'holomorphic_native_eft1_vertex_orderings': 2,
        'native_eft1_spinor_rule_i_factor': '+i',
        'checked_native_report_sha256': None,
        'native_su2_cg_product_checked': False,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report', type=Path, default=DEFAULT_REPORT)
    ap.add_argument('--check-only', action='store_true')
    args = ap.parse_args(argv)
    try:
        source_bytes = args.report.read_bytes()
        original = json.loads(source_bytes.decode('utf-8-sig'))
        native = required_report(original)
        calculated = exact_check()
        native['checked_native_report_sha256'] = hashlib.sha256(source_bytes).hexdigest()
        result = {
            'format': 't3-a-m4-native-zero-momentum-dirac-wick-v1',
            'model': 'A_-4',
            'status': 'SignedUVToNativeEFT1SpinorKernelVerified_CGAndMasterPhaseOpen',
            'source': str(args.report),
            'native_report': native,
            'exact_uv_exchange': calculated,
            'connected_heavy_fermion_exchange_evaluated_at_p_zero': True,
            'native_su2_cg_tensor_identity_certified': False,
            'master_weyl_spinor_phase_certified': False,
            'full_physical_16_model_normalization_certified': False,
            'source_resolved_16_model_release_changed': False,
            'notes': [
                'The signed Gaussian/Wick contraction uses exact native UV ordering and mass term; it is not a re-run of Matchete Match.',
                'The native EFT1 Feynman-rule i and fermion leg orderings are independently checked from the stored native report.',
                'The full SU(2) CG tensor product was not computed from native Yukawa CG matrices.',
                'Changing C to -C flips UV and contact together: their agreement cannot fix the master Weyl phase.',
                'Only A_-4 Dirac exchange is tested; Majorana B_-1/C_-1 need separate native checks.',
            ]
        }
        dest = HERE / 'normalization_probes' / ('native_zero_momentum_a_m4_' +
            hashlib.sha256(source_bytes).hexdigest()[:16] + '.json')
        print('PASS: signed p=0 Dirac F exchange and both external-lepton orderings')
        print('PASS: native A_-4 UV Feynman-rule signs and EFT1 spinor vertex structure')
        print('OPEN: SU(2) CG tensor contraction; master GammaCC/Weyl phase')
        if args.check_only:
            print('CHECK ONLY: no output written')
            return 0
        if dest.exists():
            raise FileExistsError(f'Refusing overwrite of independent report: {dest}')
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n',encoding='utf8')
        print('CREATED:', dest)
        return 0
    except (OSError, ValueError, KeyError, AssertionError, TypeError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
