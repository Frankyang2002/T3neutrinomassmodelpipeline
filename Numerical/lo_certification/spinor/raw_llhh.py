"""Audit native Matchete raw LLHH candidate sectors without normalizing C5.

Destination: Numerical/lo_certification/spinor/raw_llhh.py
The *matched* F-only loop expression is input; its terms are neither assumed
pure hbar nor pole-subtracted. Everything is source-resolved, with exact strings.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re

from .boundary import BoundaryError, MODELS, require, tag

FORMAT = 't3-f-first-raw-llhh-term-inventory-v1'
BATCH_FORMAT = 't3-f-only-raw-one-loop-native-evidence-v1'
MIXED = {'B_-1', 'C_-1'}
FIELD_RE = re.compile(r'(?<![A-Za-z0-9_])Field\[\s*([A-Za-z][A-Za-z0-9_]*)\s*,\s*(Fermion|Scalar|Vector)\s*,')
COUPLING_RE = re.compile(r'(?<![A-Za-z0-9_])Coupling\[\s*([A-Za-z][A-Za-z0-9_]*)\s*,')
PL_RE = re.compile(r'DiracProduct\[\s*GammaCC\s*,\s*Proj\[\s*-1\s*\]\s*\]')
PR_RE = re.compile(r'DiracProduct\[\s*GammaCC\s*,\s*Proj\[\s*1\s*\]\s*\]')
POLE_RE = re.compile(r'(?<![A-Za-z0-9_])(?:ϵ|ε|\\\[Epsilon\]|\\\[CurlyEpsilon\])(?![A-Za-z0-9_])')


def sha_bytes(b: bytes) -> str:
    return sha256(b).hexdigest()


def _get_json(path: Path) -> dict:
    require(path.is_file(), f'Missing native JSON: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Wrong native JSON type: {path}')
    return obj


def split_top_level(text: str) -> list[str]:
    """Lossless split on infix top-level +/-. Nested heads are left untouched.

    Matchete's actual InputForm is infix and uses spaces surrounding top-level
    signs. If it changes representation, refuse to claim an exhaustive split.
    """
    require(isinstance(text, str) and bool(text.strip()), 'Empty native InputForm')
    s = text.strip()
    if s == '0':
        return []
    stack: list[str] = []
    starts = [0]
    in_quote, escaped = False, False
    pairs = {']': '[', ')': '(', '}': '{'}
    for i, c in enumerate(s):
        if in_quote:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == '"':
                in_quote = False
            continue
        if c == '"':
            in_quote = True
            continue
        if c in '[({':
            stack.append(c)
        elif c in '])}':
            require(bool(stack) and stack[-1] == pairs[c],
                    f'Unbalanced InputForm at character {i}')
            stack.pop()
        elif not stack and c in '+-' and i > 0 and i + 1 < len(s) and \
                s[i - 1].isspace() and s[i + 1].isspace():
            starts.append(i)
    require(not stack and not in_quote, 'Unclosed Matchete expression/quote')
    pieces = [s[a:b].strip() for a, b in zip(starts, starts[1:] + [len(s)])]
    require(all(p and p not in ('+', '-') for p in pieces), 'Empty additive term')
    # Splitting retains the sign on every later term, never losing negatives.
    require(' '.join(pieces) == s, 'InputForm split is not lossless in normalized spacing')
    return pieces


def _close_square(s: str, start: int) -> int:
    """Closing square bracket for a Field[...], respecting nested brackets."""
    depth, quote, escaped = 0, False, False
    for i in range(start, len(s)):
        c = s[i]
        if quote:
            if escaped:
                escaped = False
            elif c == '\\':
                escaped = True
            elif c == '"':
                quote = False
            continue
        if c == '"':
            quote = True
        elif c == '[':
            depth += 1
        elif c == ']':
            depth -= 1
            if depth == 0:
                return i
    raise BoundaryError('Unterminated Field[...] in candidate term')


def fields_with_powers(s: str) -> list[tuple[str, str, int]]:
    out = []
    for m in FIELD_RE.finditer(s):
        end = _close_square(s, s.index('[', m.start()))
        tail = s[end + 1:]
        exponent = re.match(r'\s*\^\s*(\d+)', tail)
        power = int(exponent.group(1)) if exponent else 1
        require(1 <= power <= 16, 'Unsupported field exponent')
        out.append((m.group(1), m.group(2), power))
    return out


def classify(term: str) -> dict:
    fields = fields_with_powers(term)
    counter: Counter[str] = Counter()
    for name, _, power in fields:
        counter[name] += power
    pl, pr = bool(PL_RE.search(term)), bool(PR_RE.search(term))
    # Strict external species: no new light scalars, SM charged leptons,
    # vectors, etc., even when their name does not begin NewScalar.
    candidate = counter == Counter({'l': 2, 'H': 2}) and (pl != pr)
    if candidate:
        require(term.count('GammaCC') == 1, 'Ambiguous LLHH charge-conjugation structure')
    chirality = 'PL' if candidate and pl else 'PR' if candidate and pr else None
    return {
        'is_LLHH': candidate,
        'is_one_loop_LLHH': candidate and 'hbar' in term,
        'chirality': chirality,
        'fields': dict(sorted(counter.items())),
        'hbar': 'hbar' in term,
        'explicit_epsilon_pole': bool(POLE_RE.search(term)),
        'contains_log': 'Log[' in term,
        'couplings': sorted(set(COUPLING_RE.findall(term))),
    }


def analyze_text(model: str, text: str) -> dict:
    require(model in MODELS, f'Unknown model {model}')
    pieces = split_top_level(text)
    selected = []
    nonloop_llhh = []
    counts = Counter()
    for index, term in enumerate(pieces, start=1):
        entry = classify(term)
        if entry['is_LLHH']:
            entry = {**entry, 'source_additive_term_index': index,
                     'inputform': term,
                     'term_sha256': sha_bytes(term.encode('utf8'))}
            if entry['is_one_loop_LLHH']:
                selected.append(entry)
                counts[entry['chirality']] += 1
                counts['pole' if entry['explicit_epsilon_pole'] else 'no_explicit_pole'] += 1
                counts['with_Log' if entry['contains_log'] else 'no_Log'] += 1
            else:
                # Retain baseline/tree or counterterm pieces separately:
                # their presence cannot be interpreted as a hard loop term.
                nonloop_llhh.append(entry)
                counts['nonloop_LLHH'] += 1
        if not entry['hbar']:
            counts['non_hbar_terms_in_raw_export'] += 1
    require(sum(counts[ch] for ch in ('PL', 'PR')) == len(selected), 'Chirality tally invalid')
    require(counts['PL'] == counts['PR'],
            f'{model}: PL/PR terms disagree; cannot exclude HC safely')
    return {
        'format': FORMAT, 'model': model,
        'classification': 'Raw loop-sector LLHH candidates, NOT renormalized C5',
        'native_source_additive_terms': len(pieces),
        'selected_llhh_term_count': len(selected),
        'PL_holomorphic_terms': counts['PL'],
        'PR_conjugate_terms_excluded_from_PL': counts['PR'],
        'explicit_epsilon_pole_terms': counts['pole'],
        'nonpole_terms_not_yet_MSbar_subtracted': counts['no_explicit_pole'],
        'terms_with_Log': counts['with_Log'],
        'non_hbar_terms_in_raw_export': counts['non_hbar_terms_in_raw_export'],
        'nonloop_LLHH_term_count_retained_separately': counts['nonloop_LLHH'],
        'nonloop_LLHH_terms_not_counted_as_hard_matching': nonloop_llhh,
        'all_scalar_multiplets_active_at_fermion_threshold': True,
        'native_gamma_cc_to_master_weyl_factor': None,
        'MSbar_subtraction_certified': False,
        'physical_C5': None,
        'selected_terms': selected,
    }


def _secure_file(base: Path, rel: str) -> Path:
    require(isinstance(rel, str) and rel and '\\' not in rel and
            not rel.startswith('/') and Path(rel).parts and
            all(x not in ('.', '..') for x in Path(rel).parts),
            'Unsafe native relative source path')
    f = (base / rel).resolve()
    require(f.is_relative_to(base.resolve()) and f.is_file(),
            'Missing or escaped native input file')
    return f


def validate_batch_directory(directory: Path) -> tuple[dict, dict[str, dict]]:
    """Authenticate the pre-existing 16/16 native batch and its per-model SHA."""
    directory = directory.resolve()
    summary = _get_json(directory / 'summary.json')
    require(summary.get('format') == BATCH_FORMAT and
            summary.get('status') == 'RAW_F_FIRST_ONE_LOOP_EVIDENCE_COLLECTED_NOT_PHYSICAL',
            'Only the completed unnormalized F-only native batch is accepted')
    require(summary.get('n_models_with_raw_native_one_loop_evidence') == 16 and
            set(summary.get('models_requested', [])) == set(MODELS) and
            len(summary['models_requested']) == 16,
            'Native batch does not cover all 16 models')
    require(summary.get('single_simultaneous_fermion_threshold') is True and
            summary.get('both_scalar_multiplets_retained_at_F_threshold') is True and
            summary.get('SM_fermions_remain_active') is True and
            summary.get('n_models_with_MSbar_hard_sector_certified') == 0 and
            summary.get('native_GammaCC_to_master_Weyl_phase') is None and
            summary.get('physical_master_wilson_boundary') is None,
            'Unexpected native matching threshold/physical-status claims')
    require(summary.get('mu_F_gev') == 1e5, 'Unsupported mu_F: snapshot-specific native evidence')
    records = summary.get('model_records')
    require(isinstance(records, list) and len(records) == 16 and
            {x.get('model') for x in records} == set(MODELS),
            'Native batch record coverage inconsistent')
    out = {}
    for rec in records:
        model = rec['model']
        place = (directory / tag(model)).resolve()
        require(place.is_relative_to(directory) and place.is_dir(), 'Missing model directory')
        inspect_path = place / 'inspection.json'
        require(_get_json(inspect_path) == rec,
                f'{model}: native inspection disagrees with batch summary')
        log = place / 'wolfram.log'
        require(log.is_file() and sha_bytes(log.read_bytes()) == rec['native_log_sha256'],
                f'{model}: native log hash invalid')
        require(rec['native_exit_code'] == 0 and
                rec['hard_sector_renormalized'] is False and
                rec.get('physical_wilson_coefficient') is None,
                f'{model}: missing raw evidence or invalid physical claim')
        if model in MIXED:
            require(rec['status'] == 'MATCHeTE_RAW_MIXED_LOOP_LLHH_EXTRACTED_NOT_MSbar',
                    f'{model}: invalid mixed raw sector status')
            raw = _secure_file(place, rec['raw_hard_llhh_file'])
            expected_sha = rec['raw_hard_llhh_sha256']
            expected_count = rec['explicit_hard_LLHH_candidate_terms']
            report = _secure_file(place, rec['native_report_file'])
        else:
            require(rec['status'] == 'MATCHeTE_F_ONLY_LOOP_ORDER_ONE_RETURNED_RAW_SECTOR_UNRESOLVED',
                    f'{model}: invalid ordinary loop status')
            raw = _secure_file(place, rec['raw_loop_inputform_file'])
            expected_sha = rec['raw_loop_inputform_sha256']
            expected_count = None
            report = _secure_file(place, rec['native_report_file'])
        require(sha_bytes(raw.read_bytes()) == expected_sha and
                sha_bytes(report.read_bytes()) == rec['native_report_sha256'],
                f'{model}: native raw InputForm/report hash mismatch')
        out[model] = {'input_file': raw, 'native_report': report,
                      'expected_sha': expected_sha, 'expected_count': expected_count}
    return summary, out


def create_inventory(model: str, info: dict) -> dict:
    source = info['input_file']
    text = source.read_text(encoding='utf-8-sig')
    result = analyze_text(model, text)
    if info['expected_count'] is not None:
        require(result['selected_llhh_term_count'] == info['expected_count'],
                f'{model}: mixed native candidate count disagrees with independent audit')
    result['native_raw_inputform_sha256'] = sha_bytes(source.read_bytes())
    result['native_raw_inputform_bytes'] = source.stat().st_size
    result['raw_source_file'] = source.name
    return result
