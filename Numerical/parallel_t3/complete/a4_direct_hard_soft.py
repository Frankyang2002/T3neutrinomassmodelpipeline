"""Source-verified A(-4) direct LLHH hard/scalar-soft Laurent bridge.

Reads the F and S expressions *actually emitted* by the three-generation,
one-F-threshold/two-scalar-joint-threshold Matchete calculation.  Checks the
six explicit direct first-threshold LLHH structures and the scalar-soft
Laurent polynomial; computes their algebraic pole pairing and a finite
*ordered raw-operator* numerical diagnostic for the common-coupling limit.

NOT a renormalised EFT1 Wilson boundary, physical C5, or a replacement for
full operator mixing, counterterms, GammaCC/Weyl normalisation, or RG running.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import sympy as sp

from Numerical.parallel_t3.background.general_batch1_api import decode_native, encode_native


class DirectBridgeError(ValueError):
    """Input does not obey the explicitly verified source contraction."""


@dataclass(frozen=True)
class DirectSourceAudit:
    total_first_threshold_terms: int
    direct_llhh_terms: int
    unique_raw_hh_structures: int
    direct_source_indices: dict[str, dict[str, int]]
    soft_residue: str
    soft_finite: str
    formal_ordered_pole_cancelled: bool


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise DirectBridgeError(message)


def digest_file(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def split_matchete_sum(expression: str) -> list[str]:
    """Split only top-level + and - in Mathematica InputForm.

    Retains leading signs and never splits nested coupling, CG or NCM tensors.
    This is a bounded lexical splitter, *not* a general Mathematica evaluator.
    """
    s = expression.strip().lstrip('\ufeff')
    _require(bool(s), 'Empty Matchete expression')
    levels = {'(': 0, '[': 0, '{': 0}
    opens = '([{'
    closes = {')': '(', ']': '[', '}': '{'}
    quoted = escaped = False
    start = 0
    terms: list[str] = []
    for i, ch in enumerate(s):
        if quoted:
            if escaped:
                escaped = False
            elif ch == '\\':
                escaped = True
            elif ch == '"':
                quoted = False
            continue
        if ch == '"':
            quoted = True
            continue
        if ch in opens:
            levels[ch] += 1
        elif ch in closes:
            other = closes[ch]
            levels[other] -= 1
            _require(levels[other] >= 0, 'Unbalanced Matchete bracket')
        elif (ch in '+-' and i > 0 and s[i-1].isspace() and
              i+1 < len(s) and s[i+1].isspace() and
              all(depth == 0 for depth in levels.values())):
            term = s[start:i].strip()
            _require(bool(term), 'Empty Matchete summand')
            terms.append(term)
            start = i
    _require(not quoted and not any(levels.values()), 'Unbalanced input')
    last = s[start:].strip()
    _require(bool(last), 'Empty trailing Matchete summand')
    terms.append(last)
    return terms


_LOG_PATTERN = re.compile(
    r'\*Log\[μbar2/Coupling\[MF, \{Index\[[^,\]]+, NFlavor\]\}, 0\]\^2\]'
)
_YUK_PATTERN = {
    name: re.compile(r'Coupling\[' + name +
        r', \{Index\[([^,\]]+), Flavor\], Index\[([^,\]]+), NFlavor\]\}, 0\]')
    for name in ('y1', 'y2')
}
_MF_PATTERN = re.compile(
    r'Coupling\[MF, \{Index\[([^,\]]+), NFlavor\]\}, 0\]'
)


def _hard_signature(term: str) -> tuple[str, str, str]:
    _require(term.startswith('(2*hbar*'), 'Changed direct LLHH prefactor or sign')
    _require(term.count('hbar') == 1 and term.count('GammaCC') == 1,
             'Malformed direct LLHH loop/spinor structure')
    _require(term.count('Field[H, Scalar') == 2 and
             term.count('Field[l, Fermion') == 2,
             'Direct LLHH term does not contain exactly HH and LL')
    _require('NewScalar' not in term, 'Unintegrated scalar in direct LLHH')
    _require(term.count('CG[') == 2, 'Direct LLHH CG tensor inventory changed')
    _require(term.count('Coupling[lambdaT3,') == 1,
             'Direct LLHH no longer has exactly one mixing quartic')
    _require('y2Conj' not in term, 'Nonbaseline Yukawa in A(-4)')
    heavy = []
    for name, pattern in _YUK_PATTERN.items():
        matches = pattern.findall(term)
        _require(len(matches) == 1, 'Changed Yukawa inventory for ' + name)
        heavy.append(matches[0][1])
    mf = _MF_PATTERN.findall(term)
    _require(len(mf) in (1, 2) and set(mf + heavy) == {heavy[0]},
             'Heavy flavour indices do not contract within one mass eigenstate')
    _require(term.count('Sqrt[3]') == 1,
             'Changed 2/sqrt(3) CG prefactor')
    if 'Proj[-1]' in term and 'Proj[1]' not in term:
        chirality = 'holomorphic_HH'
        _require('Bar[Coupling[y1,' in term and
                 'Bar[Coupling[y2,' in term and
                 'Bar[Coupling[lambdaT3,' not in term and
                 term.count('CG[Bar[eps[SU2L]]') == 2 and
                 'Bar[Field[H, Scalar' not in term,
                 'Changed holomorphic CG/coupling conventions')
    elif 'Proj[1]' in term and 'Proj[-1]' not in term:
        chirality = 'anti_holomorphic_HdagHdag'
        _require('Bar[Coupling[y1,' not in term and
                 'Bar[Coupling[y2,' not in term and
                 'Bar[Coupling[lambdaT3,' in term and
                 term.count('CG[eps[SU2L]') == 2 and
                 term.count('Bar[Field[H, Scalar') == 2,
                 'Changed antiholomorphic CG/coupling conventions')
    else:
        raise DirectBridgeError('Ambiguous LLHH projector')
    has_pole = 'ϵ' in term or 'ε' in term
    has_log = 'Log[' in term
    _require(not (has_pole and has_log), 'Unexpected mixed pole/log summand')
    kind = 'pole' if has_pole else ('log' if has_log else 'constant')
    if kind == 'pole':
        _require('/(Sqrt[3]*ϵ*Coupling[MF,' in term,
                 'Different pole normalization')
    else:
        _require('/(Sqrt[3]*Coupling[MF,' in term,
                 'Different regular normalization')
    if kind == 'log':
        _require(len(_LOG_PATTERN.findall(term)) == 1,
                 'Unknown direct LLHH mass logarithm')
    else:
        _require('Log[' not in term, 'Unknown source logarithm')
    base = term.replace('ϵ*', '').replace('ε*', '')
    if kind == 'log':
        base, count = _LOG_PATTERN.subn('', base)
        _require(count == 1, 'Could not remove canonical mass log')
    return chirality, kind, base


def audit_hard_first_threshold(source: str) -> tuple[int, dict[str, dict[str, int]]]:
    terms = split_matchete_sum(source)
    direct = [(source_position, term)
              for source_position, term in enumerate(terms, start=1)
              if 'GammaCC' in term and
                 term.count('Field[H, Scalar') == 2 and
                 term.count('Field[l, Fermion') == 2 and
                 'NewScalar' not in term]
    _require(len(direct) == 6,
             f'Expected exactly six fresh A(-4) direct LLHH terms; found {len(direct)}')
    groups: dict[str, dict[str, tuple[str, int]]] = {}
    for source_position, term in direct:
        # The additive separator is retained by the balanced splitter.
        # Remove ONLY an explicit positive separator; negative amplitudes
        # must still fail the exact normalization check.
        if term.startswith('+ '):
            term = term[2:].strip()
        chi, kind, base = _hard_signature(term)
        row = groups.setdefault(chi, {})
        _require(kind not in row, 'Duplicate direct LLHH Laurent monomial')
        row[kind] = (base, source_position)
    _require(set(groups) == {'holomorphic_HH', 'anti_holomorphic_HdagHdag'},
             'Missing direct LLHH chirality')
    for chi, row in groups.items():
        _require(set(row) == {'constant', 'pole', 'log'},
                 'Incomplete direct Laurent triplet for ' + chi)
        _require(len({base for base, _ in row.values()}) == 1,
                 'Pole/log/constant have different operator tensors for ' + chi)
    return len(terms), {chi: {k: i for k, (_, i) in row.items()}
                        for chi, row in groups.items()}


_SOFT_W = re.compile(
    r'Bar\[Coupling\[T3EFT1W12, \{Index\[[^,\]]+, Flavor\], '
    r'Index\[[^,\]]+, Flavor\]\}, 0\]\]'
)
_SOFT_NAMES = (
    ('Coupling[lambdaT3, {}, 0]', 'L'),
    ('Coupling[MS1, {}, 0]', 'M1'),
    ('Coupling[MS2, {}, 0]', 'M2'),
    ('μbar2', 'Q'),
    ('ϵ', 'eps'),
    ('ε', 'eps'),
    ('hbar', 'HB'),
)


def audit_scalar_soft(source: str) -> tuple[str, str]:
    """Symbolically check the complete native scalar-soft expression.

    Accepts any algebraically equivalent sum/product ordering. No guessed CG
    normalization: the numerical coefficient is derived from the actual text.
    """
    s = source.strip().lstrip('\ufeff')
    s, count = _SOFT_W.subn('WB', s)
    _require(count >= 2, 'No conjugated three-flavour W12 in soft expression')
    for old, new in _SOFT_NAMES:
        s = s.replace(old, new)
    _require(not re.search(r'(Coupling\[|Bar\[|Field\[|CG\[|Index\[)', s),
             'Unexpected scalar-soft couplings or operator structures')
    s = s.replace('Sqrt[3]', 'sqrt(3)').replace('Log[', 'log(')
    # Mathematica functions have only the validated Sqrt and Log heads now.
    s = s.replace('[', '(').replace(']', ')').replace('^', '**')
    _require(re.fullmatch(r'[A-Za-z0-9_+*/().\s-]+', s) is not None,
             'Unexpected token in symbolic scalar-soft source')
    symbols = sp.symbols('WB L M1 M2 Q eps HB', nonzero=True)
    WB, L, M1, M2, Q, eps, HB = symbols
    try:
        expr = sp.sympify(s, locals={**{str(x): x for x in symbols},
                                    'sqrt': sp.sqrt, 'log': sp.log})
        expr = sp.cancel(expr.subs(HB, 1))
        residue = sp.simplify(sp.limit(eps*expr, eps, 0))
        finite = sp.simplify(sp.limit(expr-residue/eps, eps, 0))
        wanted_residue = -2*WB*L/sp.sqrt(3)
        wanted_finite = wanted_residue * (
            1 + (M1**2*sp.log(Q/M1**2)-M2**2*sp.log(Q/M2**2)) /
            (M1**2-M2**2))
        _require(sp.simplify(residue-wanted_residue) == 0,
                 'Native scalar-soft residue is not -2/sqrt(3) Wbar lambda')
        _require(sp.simplify(finite-wanted_finite) == 0,
                 'Native scalar-soft finite bubble differs from source expectation')
    except (TypeError, ValueError, SyntaxError) as exc:
        raise DirectBridgeError('Could not parse native scalar-soft Laurent polynomial') from exc
    return sp.sstr(residue), sp.sstr(finite)


def audit_direct_sources(fresh_hard: str, fresh_soft: str) -> DirectSourceAudit:
    count, terms = audit_hard_first_threshold(fresh_hard)
    soft_residue, soft_finite = audit_scalar_soft(fresh_soft)
    # At this strictly ordered, scalar/CG-stripped level the source hard
    # residue is +2/sqrt(3) lambda sum_r (bar y1_pr bar y2_qr/Mr)
    # while the source scalar residue is -2/sqrt(3) lambda bar W_pq.
    # The saved F tree explicitly defines W_pq=sum_r y1_pr y2_qr/Mr.
    return DirectSourceAudit(
        total_first_threshold_terms=count,
        direct_llhh_terms=6,
        unique_raw_hh_structures=2,
        direct_source_indices=terms,
        soft_residue=soft_residue,
        soft_finite=soft_finite,
        formal_ordered_pole_cancelled=True,
    )


def _mass_log_divided_difference(a: float, b: float, mass_sq: np.ndarray) -> np.ndarray:
    """(a log(a/m^2)-b log(b/m^2))/(a-b), stable for a=b."""
    _require(a > 0 and b > 0 and np.isfinite([a,b]).all(),
             'Scalar masses must be finite and positive')
    mass_sq = np.asarray(mass_sq, float)
    _require(np.all(mass_sq > 0) and np.isfinite(mass_sq).all(),
             'Heavy fermion masses must be finite and positive')
    t = (a-b)/b
    if abs(t) < 1e-5:
        # (1+t) log(1+t)/t = 1+t/2-t^2/6+t^3/12-t^4/20+t^5/30 + ...
        f = 1+t*(.5+t*(-1/6+t*(1/12+t*(-1/20+t/30))))
    else:
        f = (1+t)*math.log1p(t)/t
    return np.log(b/mass_sq)+f


def formal_numeric_direct_remainder(f_matching: dict[str, Any],
                                    s_matching: dict[str, Any], *,
                                    mu_f: float = 1e5) -> dict[str, Any]:
    """Fixed-coupling F-hard + S-soft finite ordered coefficient, diagnostic.

    Uses all 3 heavy masses, y1 and y2 flavour matrices. The F-boundary
    coupling is deliberately held fixed across the one-loop pairing; a
    genuine EFT1 Wilson RGE and CG/Weyl dictionary are not supplied here.
    """
    _require(f_matching.get('model') == 'A_-4' and
             s_matching.get('model') == 'A_-4', 'Model is not A(-4)')
    _require(f_matching.get('scale_GeV') == mu_f,
             'F matching scale must equal requested F threshold')
    _require(s_matching.get('scale_GeV') == 1e3,
             'Scalar endpoint must be 1e3 GeV')
    _require(f_matching.get('all_heavy_fermions_to_integrate_out_together') is True and
             s_matching.get('both_heavy_scalar_multiplets_to_integrate_out_together') is True,
             'F or S multiplets not integrated at common thresholds')
    values = f_matching['native_uv_rgbeta_couplings_at_muF']
    mf = np.asarray(decode_native(values['MF']), complex)
    y1 = np.asarray(decode_native(values['y1']), complex)
    y2 = np.asarray(decode_native(values['y2']), complex)
    _require(all(x.shape == (3,3) and np.isfinite(x).all() for x in (mf,y1,y2)),
             'Expected finite full 3x3 heavy MF/y1/y2')
    _require(np.linalg.norm(mf) > 0, 'Zero heavy mass matrix')
    offdiag = float(np.linalg.norm(mf-np.diag(np.diag(mf))) / np.linalg.norm(mf))
    _require(offdiag <= 1e-8,
             f'Unverified general heavy-flavour chiral rotation (offdiag={offdiag:g})')
    masses = np.diag(mf)
    _require(np.min(masses.real) > 0 and
             np.max(abs(masses.imag)) <= 1e-8*np.max(abs(masses.real)),
             'Complex or nonpositive heavy eigenmasses need chiral basis rotation')
    masses = masses.real
    lam = complex(np.asarray(decode_native(values['lambdaT3'])))
    _require(np.isfinite(lam), 'Invalid lambdaT3')
    sm = s_matching['scalar_quadratic_terms_GeV2']
    a, b = float(sm['mS1Sq']), float(sm['mS2Sq'])
    loop = _mass_log_divided_difference(a, b, masses**2)
    w = (y1.conj()/masses[None,:]) @ y2.conj().T
    common = 2/(math.sqrt(3)*16*math.pi**2)
    hard_residue = common*lam*w
    soft_residue = -hard_residue
    # Analytic form after cancellation of constant and 1/eps source pieces.
    formal_finite = common*lam*((y1.conj()*loop[None,:]/masses[None,:]) @ y2.conj().T)
    _require(np.isfinite(formal_finite).all(), 'Nonfinite direct finite remainder')
    return {
        'status': 'SourceVerifiedA4DirectHardSoftFiniteRemainderOnly_NOTPhysicalC5',
        'model':'A_-4',
        'integrated_heavy_fields_in_order':[['F'],['S1','S2']],
        'common_coupling_approximation':'F-boundary lambdaT3, y1, y2 held fixed in one-loop direct channel',
        'MF_mass_basis':'native nearly diagonal only; off-diagonal relative Frobenius <=1e-8',
        'MF_offdiagonal_relative_Frobenius':offdiag,
        'all_three_heavy_masses_GeV':masses.tolist(),
        'scalar_mass_squares_from_EFT1_renormalisable_RG_GeV2':[a,b],
        'scalar_hessian_to_native_Matchete_mass_mapping_certified':False,
        'direct_channel_scale_independence_analytic':True,
        'ordered_hard_pole_GeVinv':encode_native(hard_residue),
        'ordered_soft_pole_GeVinv':encode_native(soft_residue),
        'ordered_formal_finite_remainder_GeVinv':encode_native(formal_finite),
        'ordered_formal_finite_frobenius_GeVinv':float(np.linalg.norm(formal_finite)),
        'full_operator_basis_pole_cancellation_certified':False,
        'fermion_hard_minus_eft1_soft_plus_counterterms_performed':False,
        'complete_eft1_wilson_anomalous_dimension_running_performed':False,
        'gamma_cc_to_weyl_normalization_certified':False,
        'flavour_symmetrization_to_physical_C5_performed':False,
        'physical_C5_computed':False,
        'physical_neutrino_masses_computed':False,
        'caveat':'Only the directly shared GammaCC LLHH channel pairs algebraically after stripping common operator structure; no proof of full-basis subtraction, field/Clebsch dictionary or counterterms.',
    }
