"""Certified raw-operator UV-pole cancellation at two T3 matching thresholds.

Uses the already-produced fermion-first (F heavy) and restricted scalar-first
(X,Y heavy) Matchete InputForm outputs.  This is a symbolic/lexical matching
audit, NOT an independent EFT1 RGE, physical Weinberg projection, or MSbar
renormalization.  Source inputs are never modified.

    python -m scripts.RunGeneralT3TwoThresholdPoleAssembly --all
    python -m scripts.RunGeneralT3TwoThresholdPoleAssembly --model B_-1
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
HARD_DIR = ROOT / 'output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop'
SOFT_DIR = ROOT / 'output/general_t3_mass_eigen_matchete/eft1_scalar_threshold_probe'
AUDIT_DIR = ROOT / 'output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit'
OUTPUT_DIR = ROOT / 'output/general_t3_mass_eigen_matchete/two_threshold_operator_assembly'
MODELS = ('B_-1', 'C_-1')
POLE = '\u03f5'
LOG_F = 'Log[\u03bcbar2/Coupling[MF, {}, 0]^2]'
LOG_X = 'Log[\u03bcbar2/Coupling[MSX, {}, 0]^2]'
LOG_Y = 'Log[\u03bcbar2/Coupling[MSY, {}, 0]^2]'
SCHEME = 'raw one-loop Matchete output; regulator cancelled in LLHH sector'
FORMAT = 't3-two-threshold-raw-operator-pole-assembly-v1'
STATUS = 'PoleCancellationVerifiedCommonScaleUnprojectedNotPhysicalC5'

# The epsilon tensor is an antisymmetric SU(2) invariant. Barred epsilon has
# the same two-index antisymmetry. A pair of exchanges changes sign twice.
EPS = re.compile(
    r'CG\[(?P<bar>Bar\[)?eps\[SU2L\]\]?, '
    r'\{(?P<i>Bar\[)?Index\[d\$\$(?P<a>\d+), SU2L\[fund\]\]\]?, '
    r'(?P<j>Bar\[)?Index\[d\$\$(?P<b>\d+), SU2L\[fund\]\]\]?\}\]'
)
SU2_INDEX = re.compile(r'Index\[d\$\$(\d+), SU2L\[fund\]\]')
SU2_DUMMY_RELABELLING = {'1': '2', '2': '1', '3': '4', '4': '3'}
FAIL = 'T3 two-threshold operator assembly fail-closed: '


def fail(reason: str) -> None:
    raise ValueError(FAIL + reason)


def stem(model: str, which: str) -> str:
    if model not in MODELS or which not in ('final', 'initial'):
        fail(f'Unsupported model/which: {model}/{which}')
    return model.replace('-', 'm').replace('+', 'p') + '_' + which


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        fail(f'Missing certificate {path}')
    try:
        value = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError) as exc:
        fail(f'Cannot read certificate {path}: {exc}')
    if not isinstance(value, dict):
        fail(f'Invalid certificate object {path}')
    return value


def load_text(path: Path) -> str:
    if not path.is_file():
        fail(f'Missing source expression {path}')
    text = path.read_text(encoding='utf-8').strip()
    if not text or text == '0':
        fail(f'Empty/zero source expression {path}')
    return text


def split_terms(raw: str) -> list[tuple[int, str]]:
    """Lossless top-level InputForm Plus splitting; do not parse SU(2) tensors."""
    raw = raw.strip()
    if not raw or raw == '0' or raw.endswith((' +', ' -')):
        fail('No nonzero/complete expression')
    stack = [0, 0, 0]
    begin, sign = 0, 1
    parts: list[tuple[int, str]] = []
    for i, ch in enumerate(raw):
        if ch in '([{':
            stack['([{'.index(ch)] += 1
        elif ch in ')]}':
            stack[')]}'.index(ch)] -= 1
        if min(stack) < 0:
            fail(f'Unbalanced source expression at character {i}')
        if (ch in '+-' and stack == [0, 0, 0]
                and i > 0 and raw[i-1] == ' '
                and i+1 < len(raw) and raw[i+1] == ' '):
            term = raw[begin:i].strip()
            if not term:
                fail('Empty top-level term')
            parts.append((sign, term))
            sign = 1 if ch == '+' else -1
            begin = i + 1
    if stack != [0, 0, 0]:
        fail('Unbalanced source expression at end')
    term = raw[begin:].strip()
    if not term:
        fail('Empty final term')
    parts.append((sign, term))
    if render_terms(parts) != raw:
        fail('Top-level decomposition is not lossless')
    return parts


def render_terms(terms: list[tuple[int, str]]) -> str:
    if not terms:
        return '0'
    out = ''
    for i, (sign, value) in enumerate(terms):
        if sign not in (-1, 1) or not value:
            fail('Invalid signed term')
        if i == 0:
            out = ('-' if sign < 0 else '') + value
        else:
            out += (' + ' if sign > 0 else ' - ') + value
    return out


def strip_outer_unary_minus(sign: int, term: str) -> tuple[int, str]:
    """Matchete may write -((x)/(y)) as a single first term."""
    while term.startswith('-(') and term.endswith(')'):
        # Only strip when the opening parenthesis encloses the whole term.
        depth, closes_at = 0, -1
        for i, c in enumerate(term[1:], start=1):
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
                if depth == 0:
                    closes_at = i
                    break
        if closes_at != len(term) - 1:
            break
        sign, term = -sign, term[2:-1]
    return sign, term


def normalized_subset(terms: list[tuple[int, str]]) -> Counter[tuple[int, str]]:
    return Counter(strip_outer_unary_minus(s, t) for s, t in terms)


def kind(term: str, *, hard: bool) -> str:
    has_pole = POLE in term
    has_log = 'Log[' in term
    if has_pole and has_log:
        fail('Pole and logarithm inside the same expanded term')
    if has_pole:
        if term.count(POLE) != 1 or POLE + '*' not in term:
            fail('Unknown regulator monomial')
        return 'pole'
    if has_log:
        allowed = (LOG_F,) if hard else (LOG_X, LOG_Y)
        if sum(term.count(x) for x in allowed) != 1 or term.count('Log[') != 1:
            fail('Unexpected matching-scale logarithm')
        return 'log'
    return 'finite'


def chirality(term: str) -> str:
    if term.count('Proj[-1]') == 1 and 'Proj[1]' not in term:
        return 'PL'
    if term.count('Proj[1]') == 1 and 'Proj[-1]' not in term:
        return 'PR'
    fail('Missing/multiple Majorana chirality projectors')


def verify_llhh_term(term: str, *, hard: bool) -> None:
    if (term.count('GammaCC') != 1 or term.count('hbar') != 1
            or term.count('Field[H, Scalar') != 2
            or term.count('Field[l, Fermion') != 2
            or term.count('Coupling[MF, {}, 0]') < 1
            or any(f'Field[{x},' in term for x in ('NewScalarX','NewScalarY','NewFermion'))):
        fail('Input is not an explicit one-loop LLHH dimension-five operator')
    chirality(term)
    kind(term, hard=hard)


def pole_residue(sign: int, term: str) -> tuple[int, str]:
    sign, term = strip_outer_unary_minus(sign, term)
    if kind(term, hard=True) != 'pole':
        fail('Expected one simple pole')
    return sign, term.replace(POLE + '*', '', 1)


def relabel_su2_dummies(term: str) -> str:
    return SU2_INDEX.sub(
        lambda match: 'Index[d$$' + SU2_DUMMY_RELABELLING.get(
            match.group(1), match.group(1)) + ', SU2L[fund]]', term
    )


def canonicalize_su2_epsilon(sign: int, term: str) -> tuple[int, str]:
    """Canonicalise antisymmetric epsilon arguments and their exact signs."""
    def exchange(match: re.Match[str]) -> str:
        nonlocal sign
        a, b = int(match['a']), int(match['b'])
        if (match['bar'] is None) != (match['i'] is None):
            fail('Unsupported epsilon index/bar orientation')
        if match['i'] != match['j']:
            fail('Mismatched epsilon-index conjugations')
        if a <= b:
            return match.group(0)
        sign *= -1
        old = match.group(0)
        aa = f'Index[d$${a}, SU2L[fund]]'
        bb = f'Index[d$${b}, SU2L[fund]]'
        if match['i']:
            aa, bb = f'Bar[{aa}]', f'Bar[{bb}]'
        return old.replace(aa, '<T3SWAP>').replace(bb, aa).replace('<T3SWAP>', bb)
    replacement, count = EPS.subn(exchange, term)
    if count != 2 or replacement.count('CG[eps[SU2L]') + replacement.count('CG[Bar[eps[SU2L]]') != 2:
        fail('Expected exactly two SU(2) antisymmetric epsilon CG tensors')
    # Every SU(2) index is a closed dummy (two appearances).  No external
    # indexed leg is silently permuted by the C_-1 PR dummy relabelling.
    count_map = Counter(SU2_INDEX.findall(replacement))
    if count_map != {'1': 2, '2': 2, '3': 2, '4': 2}:
        fail(f'Unexpected SU(2) index contraction: {count_map!r}')
    return sign, replacement


def canonical_pole(sign: int, term: str, *, model: str, soft: bool) -> tuple[int, str]:
    sign, term = pole_residue(sign, term)
    if soft and model == 'C_-1' and chirality(term) == 'PL':
        # PL conjugate structure is represented with a different dummy-index
        # ordering after the scalar Match.  Relabel 1<->2 and 3<->4.
        term = relabel_su2_dummies(term)
    return canonicalize_su2_epsilon(sign, term)


def as_sectors(raw: str, *, hard: bool) -> dict[str, list[tuple[int, str]]]:
    sectors: dict[str, list[tuple[int, str]]] = {'pole': [], 'log': [], 'finite': []}
    for sign, term in split_terms(raw):
        verify_llhh_term(term, hard=hard)
        sectors[kind(term, hard=hard)].append((sign, term))
    return sectors


def _read_matching_inputs(model: str, which: str, hard_dir: Path, soft_dir: Path,
                          audit_dir: Path) -> tuple[dict, dict, dict, dict]:
    short = stem(model, which)
    hbase = short + '_eft1_fermion_first_oneloop'
    sbase = short + '_eft1_scalar_threshold_probe_c5_relevant'
    hcert = load_json(hard_dir / (hbase + '.json'))
    scert = load_json(soft_dir / (sbase + '.json'))
    acert = load_json(audit_dir / (hbase + '_hard_pole_audit.json'))
    h_expected = {
        'model': model, 'status': 'FermionFirstOneLoopRawSectorsExtractedNotRenormalized',
        'raw_hard_llhh_candidate_term_count': 204,
        'scalar_threshold_matching_performed': False,
        'physical_C5_computed': False,
    }
    s_expected = {
        'model': model,
        'status': 'ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running',
        'matching_scope': 'LeadingC5OneLoopSectorOnly',
        'excluded_four_heavy_scalar_quartic_terms': 122,
        'full_eft_validated': False,
        'tree_matching_completed': True,
        'one_loop_matching_completed': True,
        'hermiticity_check_passed': True,
        'scalar_heavy_mass_basis_detailed_check_passed': True,
        'scalar_free_lagrangian_rebased': True,
        'scalar_one_loop_correction_hbar_present': True,
        'raw_scalar_llhh_term_count': 232,
        'scalar_loop_difference_term_count': 758,
        'no_scalar_fields_in_matched_eft': True,
        'eft1_running_performed': False,
        'physical_sequential_C5_computed': False,
    }
    for cert, checks in ((hcert, h_expected), (scert, s_expected)):
        for key, value in checks.items():
            if type(cert.get(key)) is not type(value) or cert.get(key) != value:
                fail(f'Invalid matching certificate {model}: {key}: {cert.get(key)!r}')
    if acert.get('model') != model or acert.get('pole_log_finite_coefficients_identical_termwise') is not True:
        fail('Missing/invalid independently completed first-stage hard audit')
    hpath = hard_dir / (hbase + '_raw_hard_llhh_sector.inputform.txt')
    hdiff = hard_dir / (hbase + '_eft1_one_loop_correction.inputform.txt')
    spath = soft_dir / (sbase + '_raw_scalar_llhh.inputform.txt')
    sdiff = soft_dir / (sbase + '_scalar_one_loop_correction.inputform.txt')
    hraw, hfull, sraw, sfull = map(load_text, (hpath, hdiff, spath, sdiff))
    source_files = (hpath, hdiff, spath, sdiff)
    if sha(hpath) != acert.get('source_sha256', {}).get(hpath.name):
        fail('First-stage hard expression changed since independent audit')
    if sha(hdiff) != acert.get('source_sha256', {}).get(hdiff.name):
        fail('First-stage full loop expression changed since independent audit')
    if len(split_terms(hfull)) != hcert.get('loop_difference_term_count'):
        fail('Hard one-loop term-count mismatch')
    if len(split_terms(sfull)) != scert.get('scalar_loop_difference_term_count'):
        fail('Scalar one-loop term-count mismatch')
    if normalized_subset(split_terms(hraw)) - normalized_subset(split_terms(hfull)):
        fail('Hard LLHH not an exact subset of first-stage correction')
    if normalized_subset(split_terms(sraw)) - normalized_subset(split_terms(sfull)):
        fail('Scalar LLHH not an exact subset of scalar-stage correction')
    return ({'raw': hraw, 'cert': hcert}, {'raw': sraw, 'cert': scert},
            {'hard_audit': acert}, {'files': source_files})


def _at_mu_f(sign: int, term: str) -> tuple[int, str] | None:
    if LOG_F in term:
        if term.count(LOG_F) != 1:
            fail('Multiple hard logs in term')
        return None  # The log is zero when mu_bar^2 = MF^2.
    term = term.replace(LOG_X, 'Log[Coupling[MF, {}, 0]^2/Coupling[MSX, {}, 0]^2]')
    term = term.replace(LOG_Y, 'Log[Coupling[MF, {}, 0]^2/Coupling[MSY, {}, 0]^2]')
    if '\u03bcbar2' in term:
        fail('Unexpected matching scale in fixed-muF export')
    return sign, term


def assemble(model: str, *, which: str = 'final', hard_dir: Path = HARD_DIR,
             soft_dir: Path = SOFT_DIR, audit_dir: Path = AUDIT_DIR,
             output_dir: Path = OUTPUT_DIR) -> dict[str, Any]:
    hard, soft, ref, provenance = _read_matching_inputs(model, which, hard_dir, soft_dir, audit_dir)
    hs = as_sectors(hard['raw'], hard=True)
    ss = as_sectors(soft['raw'], hard=False)
    hc = {k: len(v) for k, v in hs.items()}
    sc = {k: len(v) for k, v in ss.items()}
    if hc != {'pole': 68, 'log': 68, 'finite': 68}:
        fail(f'Unrecognized first-stage LLHH inventory {hc}')
    if sc != {'pole': 68, 'log': 100, 'finite': 64}:
        fail(f'Unrecognized scalar-stage LLHH inventory {sc}')
    xlogs = sum(LOG_X in t for _, t in ss['log'])
    ylogs = sum(LOG_Y in t for _, t in ss['log'])
    if (xlogs, ylogs) != (50, 50):
        fail(f'Unrecognized two-mass logarithm distribution ({xlogs},{ylogs})')
    h_poles = Counter(canonical_pole(s, t, model=model, soft=False) for s, t in hs['pole'])
    s_poles = Counter((-s2, op) for s, t in ss['pole']
                      for s2, op in [canonical_pole(s, t, model=model, soft=True)])
    if h_poles != s_poles or len(h_poles) != 68 or any(v != 1 for v in h_poles.values()):
        fail(f'Two-threshold UV-pole residues did not cancel: hard-only={sum((h_poles-s_poles).values())}, scalar-only={sum((s_poles-h_poles).values())}')
    survivors = hs['finite'] + hs['log'] + ss['finite'] + ss['log']
    if len(survivors) != 300:
        fail('Wrong number of finite and logarithmic terms in raw composite')
    at_f = [v for pair in survivors if (v := _at_mu_f(*pair)) is not None]
    if len(at_f) != 232:  # hard log (68) has vanished
        fail('Unexpected term count at mu_bar^2=MF^2')
    for _, t in survivors + at_f:
        if POLE in t:
            fail('Uncancelled UV regulator in assembled one-loop sector')
    operator_base = stem(model, which) + '_two_threshold_raw_llhh'
    output_dir.mkdir(parents=True, exist_ok=True)
    names = {
        'hard_pole_residue_canonical': operator_base + '_hard_pole_residue_canonical.inputform.txt',
        'scalar_negative_pole_residue_canonical': operator_base + '_scalar_negative_pole_residue_canonical.inputform.txt',
        'fixed_common_mu_unprojected': operator_base + '_fixed_common_mu_unprojected.inputform.txt',
        'fixed_muF_unprojected': operator_base + '_fixed_muF_unprojected.inputform.txt',
        'hard_finite': operator_base + '_hard_finite.inputform.txt',
        'soft_finite': operator_base + '_soft_finite.inputform.txt',
        'soft_log': operator_base + '_soft_logs.inputform.txt',
    }
    # Sorted canonical residues guarantee deterministic comparison of the
    # two signed operator multisets, without evaluating a large Mathematica sum.
    poles = sorted(h_poles.elements(), key=lambda p: (p[1], p[0]))
    texts = {
        'hard_pole_residue_canonical': render_terms(poles),
        'scalar_negative_pole_residue_canonical': render_terms(sorted(s_poles.elements(), key=lambda p: (p[1], p[0]))),
        'fixed_common_mu_unprojected': render_terms(survivors),
        'fixed_muF_unprojected': render_terms(at_f),
        'hard_finite': render_terms(hs['finite']),
        'soft_finite': render_terms(ss['finite']),
        'soft_log': render_terms(ss['log']),
    }
    for label, name in names.items():
        (output_dir / name).write_text(texts[label] + '\n', encoding='utf-8')
    report = {
        'format': FORMAT, 'status': STATUS, 'model': model, 'which': which,
        'input_provenance_sha256': {p.name: sha(p) for p in provenance['files']},
        'hard_matching_scope': 'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'scalar_matching_scope': 'LeadingC5OneLoopSectorOnly',
        'excluded_scalar_self_quartic_terms': 122,
        'source_one_loop_term_counts': {'fermion_first': hard['cert']['loop_difference_term_count'],
                                        'scalar_threshold': soft['cert']['scalar_loop_difference_term_count']},
        'raw_hard_sector_counts': hc, 'raw_scalar_sector_counts': sc,
        'scalar_logarithm_terms_by_mass': {'MSX': xlogs, 'MSY': ylogs},
        'hard_and_scalar_pole_residues_opposite': True,
        'matched_canonical_pole_term_count': len(h_poles),
        'pole_cancellation_basis': ('exact signed operator matching, SU(2) epsilon antisymmetry; '
           'C_-1 PL additionally uses closed-dummy relabelling 1<->2,3<->4'),
        'independent_eft1_rge_consistency_verified': False,
        'complete_scalar_potential_hermiticity_verified': False,
        'fixed_order_common_scale_operator_assembly_only': True,
        'fixed_order_common_mu_nonpole_terms': len(survivors),
        'fixed_muF_nonpole_terms': len(at_f),
        'fixed_muF_choice': 'μbar2=MF^2; two scalar logarithms become Log[MF^2/MSX^2] and Log[MF^2/MSY^2]',
        'finite_hard_operator_is_renormalized_MSbar_C5': False,
        'weinberg_flavour_su2_projection_performed': False,
        'numerical_neutrino_masses_computed': False,
        'eft1_running_performed': False,
        'physical_sequential_C5_computed': False,
        'hbar_convention': 'symbolic hbar retained; if numerical use 1/(16 Pi^2)',
        'all_outputs_are_full_inputform_LLHH_operators_not_C5_matrices': True,
        'outputs': names,
    }
    report_name = stem(model, which) + '_two_threshold_raw_llhh_audit.json'
    (output_dir / report_name).write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    option = p.add_mutually_exclusive_group(required=True)
    option.add_argument('--all', action='store_true')
    option.add_argument('--model', choices=MODELS)
    p.add_argument('--which', choices=('final', 'initial'), default='final')
    p.add_argument('--hard-dir', type=Path, default=HARD_DIR)
    p.add_argument('--soft-dir', type=Path, default=SOFT_DIR)
    p.add_argument('--hard-audit-dir', type=Path, default=AUDIT_DIR)
    p.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = p.parse_args()
    reports = []
    for model in MODELS if args.all else (args.model,):
        rep = assemble(model, which=args.which, hard_dir=args.hard_dir,
                       soft_dir=args.soft_dir, audit_dir=args.hard_audit_dir,
                       output_dir=args.output_dir)
        print(f'{model}: hard+scalar UV poles cancel 68/68; scalar logs X/Y=50/50; '
              f'nonpole LLHH at common μ=300 terms, at μ=MF=232 terms; '
              'physical C5 NOT computed', flush=True)
        reports.append({'model': model, 'status': rep['status'],
                        'pole_pairs': rep['matched_canonical_pole_term_count'],
                        'unprojected_operator_file': rep['outputs']['fixed_muF_unprojected']})
    summary = {'format': FORMAT, 'status': STATUS, 'models': reports,
               'independent_rge_certified': False, 'physical_C5_computed': False}
    out = args.output_dir / ('summary_' + args.which + '.json')
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print('Summary:', out)


if __name__ == '__main__':
    main()
