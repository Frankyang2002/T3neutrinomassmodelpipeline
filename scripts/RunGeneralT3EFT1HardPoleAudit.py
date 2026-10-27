"""Audit B_-1/C_-1 F-first EFT1 Matchete raw LLHH poles and matching logs.

Reads existing *successful* Matchete exports without rerunning Matchete.
No field-contraction reduction, MSbar certification, or physical C5 claim.

    python -m scripts.RunGeneralT3EFT1HardPoleAudit --all
    python -m scripts.RunGeneralT3EFT1HardPoleAudit --model B_-1
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop'
OUTPUT = ROOT / 'output/general_t3_mass_eigen_matchete/eft1_hard_pole_audit'
MODELS = ('B_-1', 'C_-1')
INPUT_FORMAT = 't3-mass-eigen-fermion-first-eft1-oneloop-v1'
INPUT_STATUS = 'FermionFirstOneLoopRawSectorsExtractedNotRenormalized'
OUTPUT_FORMAT = 't3-fermion-first-hard-llhh-pole-log-audit-v1'
OUTPUT_STATUS = 'RawHardLLHHOperatorPoleLogFactorizationVerifiedNotPhysicalC5'
REGULATOR = '\u03f5'  # Exact Unicode symbol in the uploaded Matchete InputForm.
MATCH_LOG = 'Log[\u03bcbar2/Coupling[MF, {}, 0]^2]'
ZERO_ERROR = 'T3 EFT1 hard-pole audit fail-closed: '
QUARTICS = (
    'lambdaT3',
    'lambdaT3ExtraQ020020I1',
    'lambdaT3ExtraQ020200I1',
)
COUPLING = re.compile(r'(?:Bar\[)?Coupling\[([A-Za-z][A-Za-z0-9]*),')


def fail(message: str) -> None:
    raise ValueError(ZERO_ERROR + message)


def model_stem(model: str, which: str = 'final') -> str:
    if model not in MODELS or which not in ('final', 'initial'):
        fail(f'Unsupported model/which: {model!r}/{which!r}')
    return model.replace('-', 'm').replace('+', 'p') + '_' + which + '_eft1_fermion_first_oneloop'


def _balanced_terms(raw: str) -> list[tuple[int, str]]:
    """Losslessly split *top-level* InputForm Plus into signed terms.

    This is a lexical decomposition, not Mathematica evaluation or an
    unsound attempt to parse or simplify internal tensor contractions.
    """
    raw = raw.strip()
    if not raw or raw == '0':
        fail('Raw hard LLHH expression is absent or zero')
    par = square = curly = 0
    start = 0
    sign = 1
    parts: list[tuple[int, str]] = []
    for index, char in enumerate(raw):
        if char == '(':
            par += 1
        elif char == ')':
            par -= 1
        elif char == '[':
            square += 1
        elif char == ']':
            square -= 1
        elif char == '{':
            curly += 1
        elif char == '}':
            curly -= 1
        if min(par, square, curly) < 0:
            fail(f'Unbalanced Mathematica delimiters at offset {index}')
        if (
            char in '+-' and par == square == curly == 0
            and index > 0 and raw[index - 1] == ' '
            and index + 1 < len(raw) and raw[index + 1] == ' '
        ):
            part = raw[start:index].strip()
            if not part:
                fail('Empty top-level term')
            parts.append((sign, part))
            sign = 1 if char == '+' else -1
            start = index + 1
    if (par, square, curly) != (0, 0, 0):
        fail('Unbalanced Mathematica delimiters at end of source')
    part = raw[start:].strip()
    if not part:
        fail('Empty final term')
    parts.append((sign, part))
    if render_terms(parts) != raw:
        fail('Top-level term decomposition is not byte-reconstructible')
    return parts


def render_terms(terms: list[tuple[int, str]]) -> str:
    if not terms:
        return '0'
    result = ''
    for i, (sign, term) in enumerate(terms):
        if sign not in (-1, 1) or not term:
            fail('Malformed signed term')
        if i == 0:
            result = ('- ' if sign < 0 else '') + term
        else:
            result += (' + ' if sign > 0 else ' - ') + term
    return result


def _sector(term: str) -> str:
    pole = REGULATOR in term
    logged = 'Log[' in term
    if pole and logged:
        fail('Term mixes pole and matching-scale log; cannot isolate')
    if pole:
        if term.count(REGULATOR) != 1 or REGULATOR + '*' not in term:
            fail('Unknown regulator representation or multiple poles')
        return 'pole'
    if logged:
        if term.count(MATCH_LOG) != 1 or term.count('Log[') != 1:
            fail('Unknown logarithm or multiple logarithms')
        return 'log'
    return 'finite'


def _operator_coefficient_term(term: str, sector: str) -> str:
    """Remove only known *multiplicative* regulator or logarithm tokens.

    Pole and log residues are still FULL LLHH OPERATORS (with CGs/H/L).
    Nothing is stripped from SU(2), Lorentz or flavour tensor factors.
    """
    if sector == 'pole':
        result = term.replace(REGULATOR + '*', '')
    elif sector == 'log':
        result = term.replace('*' + MATCH_LOG, '')
    else:
        result = term
    # Original B expression stores 1/(ϵ*MF), finite counterpart 1/MF.
    # Parentheses around the *sole* MF denominator do not change content.
    result = result.replace('/(Coupling[MF, {}, 0])', '/Coupling[MF, {}, 0]')
    if REGULATOR in result or 'Log[' in result:
        fail('Residual unprocessed regulator/logarithm after isolation')
    if result == term and sector != 'finite':
        fail(f'Failed to extract the {sector} coefficient')
    return result


def _chirality(term: str) -> str:
    if term.count('Proj[-1]') == 1 and 'Proj[1]' not in term:
        return 'holomorphic_PL'
    if term.count('Proj[1]') == 1 and 'Proj[-1]' not in term:
        return 'conjugate_PR'
    fail('Unrecognized or mixed Majorana projector structure')


def _quartic(term: str) -> str:
    names = [x for x in COUPLING.findall(term) if x in QUARTICS]
    if len(names) != 1:
        fail(f'Expected exactly one quartic source coupling, got {names!r}')
    return names[0]


def _validate_operator(term: str) -> None:
    if (term.count('hbar*') != 1 or 'hbar^' in term or
            term.count('GammaCC') != 1):
        fail('Wrong loop order or Majorana bilinear content')
    if term.count('Field[H, Scalar') != 2 or term.count('Field[l, Fermion') != 2:
        fail('Not exactly two Higgs and two lepton fields')
    if any(f'Field[{s},' in term for s in ('NewFermion', 'NewScalarX', 'NewScalarY')):
        fail('Unexpected heavy/light BSM field in hard LLHH sector')
    if term.count('Coupling[MF, {}, 0]') < 1:
        fail('Missing fermion mass in loop-order dimension-five term')
    _chirality(term)
    _quartic(term)
    if 'ε' in term or '\\[Epsilon]' in term or 'epsilon' in term.lower():
        fail('Unrecognized epsilon-like regulator')


def factorize_hard_sector(raw: str, *, expected_count: int = 204) -> dict[str, Any]:
    terms = _balanced_terms(raw)
    if len(terms) != expected_count:
        fail(f'Expected {expected_count} raw LLHH terms; found {len(terms)}')
    buckets: dict[tuple[int, str], dict[str, list[tuple[int, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    by_sector: Counter[str] = Counter()
    by_chirality: Counter[str] = Counter()
    by_quartic: Counter[str] = Counter()
    source_sectors: dict[str, list[tuple[int, str]]] = {
        'finite': [], 'pole': [], 'log': []
    }
    coefficient_sectors: dict[str, list[tuple[int, str]]] = {
        'finite': [], 'pole': [], 'log': []
    }
    chirality_sectors: dict[str, list[tuple[int, str]]] = {
        'holomorphic_PL': [], 'conjugate_PR': []
    }
    for sign, term in terms:
        _validate_operator(term)
        sector = _sector(term)
        coeff = _operator_coefficient_term(term, sector)
        chirality = _chirality(term)
        quartic = _quartic(term)
        key = (sign, coeff)
        buckets[key][sector].append((sign, term))
        source_sectors[sector].append((sign, term))
        coefficient_sectors[sector].append((sign, coeff))
        by_sector[sector] += 1
        by_chirality[f'{chirality}.{sector}'] += 1
        by_quartic[f'{quartic}.{sector}'] += 1
        if sector == 'finite':
            chirality_sectors[chirality].append((sign, coeff))
    if set(by_sector) != {'finite', 'pole', 'log'}:
        fail(f'Unexpected sector inventory: {dict(by_sector)}')
    if len(buckets) * 3 != len(terms):
        fail('Operator structures are not grouped into 3-way triplets')
    for key, sectors in buckets.items():
        if {k: len(v) for k, v in sectors.items()} != {
            'finite': 1, 'pole': 1, 'log': 1
        }:
            fail('Pole/log/finite coefficients do not match termwise')
    n = len(buckets)
    if any(by_sector[k] != n for k in ('finite', 'pole', 'log')):
        fail('Unequal pole/log/finite term counts')
    if by_chirality['holomorphic_PL.finite'] != n // 2 or n % 2:
        fail('Unexpected holomorphic/projector distribution')
    if by_chirality['conjugate_PR.finite'] != n // 2:
        fail('Unexpected conjugate/projector distribution')
    finite_key = Counter(coefficient_sectors['finite'])
    for k in ('pole', 'log'):
        if Counter(coefficient_sectors[k]) != finite_key:
            fail('Pole/log coefficient not equal to finite coefficient')
    source_strings = {k: render_terms(v) for k, v in source_sectors.items()}
    coefficient_strings = {k: render_terms(v) for k, v in coefficient_sectors.items()}
    if any(REGULATOR in v or 'Log[' in v for v in coefficient_strings.values()):
        fail('Non-finite factor in extracted coefficient sector')
    return {
        'term_count': len(terms),
        'independent_signed_operator_structures': n,
        'sector_counts': dict(by_sector),
        'chirality_counts': dict(by_chirality),
        'quartic_source_counts': dict(by_quartic),
        'pole_log_finite_coefficients_identical_termwise': True,
        'original_sector_strings': source_strings,
        'coefficient_sector_strings': coefficient_strings,
        'holomorphic_finite_operator': render_terms(chirality_sectors['holomorphic_PL']),
        'conjugate_finite_operator': render_terms(chirality_sectors['conjugate_PR']),
    }


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _read_report(path: Path, model: str, stem: str) -> dict[str, Any]:
    if not path.is_file():
        fail(f'Input certificate not found: {path}')
    try:
        cert = json.loads(path.read_text(encoding='utf-8'))
    except (ValueError, OSError) as exc:
        fail(f'Unreadable matching certificate {path}: {exc}')
    if not isinstance(cert, Mapping):
        fail('Matchete matching certificate must be a JSON object')
    expected = {
        'format': INPUT_FORMAT, 'status': INPUT_STATUS, 'model': model,
        'matching_approach': 'FermionFirstTreeAndOneLoopScalarsRemainLight',
        'matching_EFT_order': 5, 'matching_loop_order': 1,
        'source_quartic_count': 17, 'source_yukawa_count': 2,
        'fermion_heavy_mass_basis_detailed_check_passed': True,
        'tree_matching_completed': True, 'one_loop_matching_completed': True,
        'raw_hard_llhh_no_fermion_field': True,
        'loop_difference_fermion_absent': True,
        'raw_hard_llhh_scheme_renormalized': False,
        'raw_hard_llhh_pole_rge_certified': False,
        'eft1_rge_performed': False,
        'scalar_threshold_matching_performed': False,
        'one_loop_hard_fermion_matching_performed': False,
        'physical_C5_computed': False,
        'sequential_EFT1_completed': False,
        'physical_matching_performed': False,
    }
    for key, target in expected.items():
        if cert.get(key) != target or type(cert.get(key)) is not type(target):
            fail(f'Invalid input matching certificate {key}: {cert.get(key)!r}')
    required_files = {
        'raw_hard_llhh_file': stem + '_raw_hard_llhh_sector.inputform.txt',
        'one_loop_difference_file': stem + '_eft1_one_loop_correction.inputform.txt',
        'llss_inputform_file': stem + '_llss_tree.inputform.txt',
    }
    for key, expected_name in required_files.items():
        if cert.get(key) != expected_name:
            fail(f'Invalid source filename for {key}')
    if cert.get('raw_hard_llhh_candidate_term_count') != 204:
        fail('Unexpected raw hard term count in input certificate')
    return cert


def audit_model(model: str, *, input_dir: Path = INPUT, output_dir: Path = OUTPUT,
                which: str = 'final') -> dict[str, Any]:
    stem = model_stem(model, which)
    source_path = input_dir / (stem + '.json')
    cert = _read_report(source_path, model, stem)
    raw_path = input_dir / cert['raw_hard_llhh_file']
    difference_path = input_dir / cert['one_loop_difference_file']
    llss_path = input_dir / cert['llss_inputform_file']
    for path in (raw_path, difference_path, llss_path):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f'Missing/empty source artifact: {path}')
    raw = raw_path.read_text(encoding='utf-8').strip()
    full_difference = difference_path.read_text(encoding='utf-8').strip()
    llss = llss_path.read_text(encoding='utf-8').strip()
    if not all(k in llss for k in ('GammaCC', 'Field[NewScalarX', 'Field[NewScalarY')):
        fail('Input LLSS tree expression is incomplete')
    if 'Field[NewFermion' in full_difference or 'hbar' not in full_difference:
        fail('One-loop EFT1 correction violates source certificate')
    source_term_counter = Counter(_balanced_terms(full_difference))
    if len(_balanced_terms(full_difference)) != cert['loop_difference_term_count']:
        fail('One-loop correction term count disagrees with certificate')
    hard_terms = Counter(_balanced_terms(raw))
    if hard_terms - source_term_counter:
        fail('Raw hard sector is not an exact signed subset of the one-loop EFT')
    analyzed = factorize_hard_sector(raw, expected_count=204)
    name = stem + '_hard_pole_audit'
    texts = {
        'raw_pole_sector': analyzed['original_sector_strings']['pole'],
        'raw_log_sector': analyzed['original_sector_strings']['log'],
        'raw_finite_sector': analyzed['original_sector_strings']['finite'],
        'pole_residue_llhh': analyzed['coefficient_sector_strings']['pole'],
        'matching_log_coefficient_llhh': analyzed['coefficient_sector_strings']['log'],
        'finite_at_muF_llhh': analyzed['coefficient_sector_strings']['finite'],
        'holomorphic_finite_at_muF_llhh': analyzed['holomorphic_finite_operator'],
        'conjugate_finite_at_muF_llhh': analyzed['conjugate_finite_operator'],
    }
    files = {k: name + '_' + k + '.inputform.txt' for k in texts}
    report = {
        'format': OUTPUT_FORMAT,
        'status': OUTPUT_STATUS,
        'model': model,
        'source_matching_json': source_path.name,
        'source_sha256': {p.name: _sha(p) for p in
                          (source_path, raw_path, difference_path, llss_path)},
        'raw_hard_term_count': analyzed['term_count'],
        'distinct_signed_llhh_structures': analyzed['independent_signed_operator_structures'],
        'sector_counts': analyzed['sector_counts'],
        'chirality_counts': analyzed['chirality_counts'],
        'quartic_source_counts': analyzed['quartic_source_counts'],
        'pole_log_finite_coefficients_identical_termwise': True,
        'pole_log_coefficient_equality_is_independent_rge_test': False,
        'raw_operator_factorization':
            'A_LLHH * (1 + 1/ϵ + Log[μbar2/MF^2])',
        'A_LLHH_contains_hbar': True,
        'hbar_value_to_use_if_numerical': '1/(16 Pi^2)',
        'matching_scale_scenario': 'μbar2 = MF^2, Log[μbar2/MF^2] = 0',
        'operator_hard_finite_at_muF_is_scheme_certified': False,
        'independent_counterterm_pole_check_performed': False,
        'independent_eft1_llss_to_c5_rge_check_performed': False,
        'flavour_su2_c5_operator_projection_performed': False,
        'physical_majorana_c5_computed': False,
        'eft1_running_performed': False,
        'scalar_threshold_matching_performed': False,
        'sequential_matching_completed': False,
        'output_expression_files': files,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for kind, filename in files.items():
        (output_dir / filename).write_text(texts[kind] + '\n', encoding='utf-8')
    report_path = output_dir / (name + '.json')
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n',
                           encoding='utf-8')
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--model', choices=MODELS)
    group.add_argument('--all', action='store_true')
    parser.add_argument('--which', choices=('initial', 'final'), default='final')
    parser.add_argument('--input-dir', type=Path, default=INPUT)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    rows = []
    for model in MODELS if args.all else (args.model,):
        report = audit_model(model, input_dir=args.input_dir,
                             output_dir=args.output_dir, which=args.which)
        rows.append({'model': model, 'status': report['status'],
                     'counts': report['sector_counts'],
                     'termwise_identity_verified': True})
        print(f'{model}: raw LLHH=204, pole=68, log=68, finite=68; '
              'pole/log/finite operator coefficients identical termwise; '
              'physical C5 NOT computed', flush=True)
    summary = {
        'format': OUTPUT_FORMAT, 'status': OUTPUT_STATUS,
        'models': rows, 'independent_pole_rge_certified': False,
        'physical_c5_computed': False,
        'sequential_eft1_matching_performed': False,
    }
    summary_path = args.output_dir / ('summary_' + args.which + '.json')
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n',
                            encoding='utf-8')
    print('Summary:', summary_path)


if __name__ == '__main__':
    main()
