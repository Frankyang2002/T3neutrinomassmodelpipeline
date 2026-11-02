"""Lossless source-resolved pole/constant/log identity audit for native LLHH.

This is an algebraic organization of raw Matchete F-only output, NOT an
MSbar subtraction or EFT1 soft-region subtraction. In particular neither
GammaCC -> master Weyl nor a physical C5 is supplied by this module.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re

from .boundary import BoundaryError, MODELS, require, tag
from .raw_llhh import FORMAT as INVENTORY_FORMAT, create_inventory, sha_bytes, validate_batch_directory

FORMAT = 't3-f-first-source-resolved-pole-log-bridge-v1'
EXPECTED_SNAPSHOT = 'e515d82d118284bab72dbad6ba593aabc529aea08f0a7e47b08cd5049361cfcd'
EPSILON = re.compile(r'(?<![A-Za-z0-9_])(?:ϵ|ε|\\\[Epsilon\]|\\\[CurlyEpsilon\])(?![A-Za-z0-9_])')
# Matchete's current exact one-loop export uses the literal Greek ϵ. Keep
# other encodings fail-closed until an independently inspected source demands them.
NATIVE_EPSILON = re.compile(r'(?<![A-Za-z0-9_])ϵ(?![A-Za-z0-9_])')
BARE_MF_DENOM = re.compile(r'/\(Coupling\[MF,[^()]+\]\)')


def digest(text: str) -> str:
    return sha256(text.encode('utf-8')).hexdigest()


def _close_square(text: str, at: int) -> int:
    require(text.startswith('Log[', at), 'Expected native Log[ factor')
    depth = 0
    for j in range(at + 3, len(text)):
        if text[j] == '[':
            depth += 1
        elif text[j] == ']':
            depth -= 1
            if depth == 0:
                return j + 1
    raise BoundaryError('Unclosed native mass logarithm')


def normalize_native_kernel(term: str) -> tuple[str, str, str]:
    """Remove exactly one optional multiplicative Log or 1/epsilon factor.

    Return (signed exact shared kernel, kind, log argument). No coefficients,
    signs, flavour indices, barred couplings or CG structures are removed.
    The first InputForm additive term has an implicit leading plus sign.
    """
    require(isinstance(term, str) and bool(term.strip()), 'Empty native LLHH term')
    value = term.strip()
    if not value.startswith(('+ ', '- ')):
        require(not value.startswith(('+', '-')), 'Unsupported initial term sign syntax')
        value = '+ ' + value

    log_argument = ''
    i = value.find('Log[')
    if i != -1:
        require(value.find('Log[', i + 4) == -1, 'More than one mass logarithm in LLHH term')
        end = _close_square(value, i)
        log_argument = value[i:end]
        require(log_argument.startswith('Log[μbar2/Coupling[MF,') and
                log_argument.endswith(']^2]'),
                'Unsupported mass log: not log(mu-bar squared / heavy MF squared)')
        require(i > 0 and value[i - 1] == '*', 'Log is not a multiplicative numerator factor')
        value = value[:i - 1] + value[end:]

    eps_match = list(NATIVE_EPSILON.finditer(value))
    require(len(eps_match) <= 1, 'Multiple epsilon factors in a native LLHH term')
    if EPSILON.search(value) and not eps_match:
        raise BoundaryError('Unsupported regulator token in native LLHH term')
    if eps_match:
        p, end = eps_match[0].span()
        if end < len(value) and value[end] == '*':
            value = value[:p] + value[end + 1:]
        elif p > 0 and value[p - 1] == '*':
            value = value[:p - 1] + value[end:]
        elif p > 0 and value[p - 1] == '/':
            value = value[:p - 1] + value[end:]
        else:
            raise BoundaryError('Epsilon is not an isolated multiplicative denominator factor')
        require(not NATIVE_EPSILON.search(value), 'Epsilon stripping incomplete')

    # Matchete prints /(eps*MF) versus /MF; stripping eps leaves /(MF).
    # Collapse parentheses only for a single heavy-mass Coupling atom.
    def unparenthesize(match: re.Match[str]) -> str:
        return '/' + match.group(0)[2:-1]
    value = BARE_MF_DENOM.sub(unparenthesize, value)

    kind = 'pole' if eps_match else 'scale_log' if log_argument else 'constant'
    require(bool(value.strip()), 'Removed all content from native term')
    return value, kind, log_argument


def validate_and_group(model: str, inventory: dict) -> dict:
    require(model in MODELS and inventory.get('model') == model and
            inventory.get('format') == INVENTORY_FORMAT,
            'Wrong inventory model/format')
    require(inventory.get('MSbar_subtraction_certified') is False and
            inventory.get('physical_C5') is None and
            inventory.get('native_gamma_cc_to_master_weyl_factor') is None,
            'Incoming inventory must remain unrenormalized and phase-uncertified')
    terms = inventory.get('selected_terms')
    require(isinstance(terms, list) and len(terms) == inventory.get('selected_llhh_term_count'),
            'Missing source-resolved LLHH term inventory')
    require(len(terms) == inventory['PL_holomorphic_terms'] +
            inventory['PR_conjugate_terms_excluded_from_PL'], 'Invalid chirality totals')

    groups = defaultdict(dict)
    indices = set()
    source_counts = Counter()
    for entry in terms:
        require(isinstance(entry, dict), 'Invalid native source term record')
        text = entry.get('inputform')
        require(isinstance(text, str) and entry.get('term_sha256') == digest(text),
                'Source term InputForm SHA256 mismatch')
        index = entry.get('source_additive_term_index')
        require(type(index) is int and index > 0 and index <= inventory['native_source_additive_terms']
                and index not in indices, 'Duplicate or invalid native source index')
        indices.add(index)
        chirality = entry.get('chirality')
        require(chirality in ('PL', 'PR') and entry.get('is_LLHH') is True and
                entry.get('is_one_loop_LLHH') is True and entry.get('hbar') is True and
                entry.get('fields') == {'H': 2, 'l': 2},
                'Not an authenticated hbar LLHH candidate')
        # Reverse-check original classification: don't trust arbitrary JSON metadata.
        from .raw_llhh import classify
        original = classify(text)
        require(original['chirality'] == chirality and original['is_one_loop_LLHH'] and
                original['explicit_epsilon_pole'] == entry['explicit_epsilon_pole'] and
                original['contains_log'] == entry['contains_log'] and
                original['couplings'] == entry['couplings'], 'Term changed since LLHH classification')
        kernel, kind, log_argument = normalize_native_kernel(text)
        require(bool(entry['explicit_epsilon_pole']) == (kind == 'pole') and
                bool(entry['contains_log']) == (kind == 'scale_log'),
                'Regulator/log classification inconsistent with original source')
        require('hbar' in kernel and 'GammaCC' in kernel and 'Coupling[MF,' in kernel,
                'Unexpected loop kernel structure')
        key = (chirality, kernel)
        require(kind not in groups[key],
                'Duplicate pole/constant/log member for same exact native kernel')
        groups[key][kind] = {'source_additive_term_index': index,
                             'term_sha256': entry['term_sha256'],
                             'original_inputform': text,
                             'log_argument': log_argument}
        source_counts[chirality] += 1
        source_counts[kind] += 1

    require(source_counts['PL'] == inventory['PL_holomorphic_terms'] and
            source_counts['PR'] == inventory['PR_conjugate_terms_excluded_from_PL'] and
            source_counts['pole'] == inventory['explicit_epsilon_pole_terms'] and
            source_counts['constant'] + source_counts['scale_log'] ==
            inventory['nonpole_terms_not_yet_MSbar_subtracted'],
            'Pole/constant/log populations disagree with native inventory')

    records = []
    for (chirality, kernel), parts in groups.items():
        require(set(parts) == {'pole', 'constant', 'scale_log'},
                f'{model}: unmatched Laurent triplet in {chirality}; '
                f'present {sorted(parts)}')
        require(parts['scale_log']['log_argument'], 'Missing log(mu^2 / M^2) argument')
        ids = {k: v['source_additive_term_index'] for k, v in parts.items()}
        records.append({
            'chirality': chirality,
            'kernel_inputform': kernel,
            'kernel_sha256': digest(kernel),
            'source_indices': ids,
            'source_term_sha256': {k: v['term_sha256'] for k, v in parts.items()},
            'mass_log_inputform': parts['scale_log']['log_argument'],
            'pole_residue_equals_constant_kernel': True,
            'log_scale_coefficient_equals_pole_residue': True,
            'coefficient_factor': '1 + 1/ϵ + Log[μbar2/MF^2]',
        })
    records.sort(key=lambda r: (r['chirality'], min(r['source_indices'].values())))
    require(len(records) * 3 == len(terms),
            'One-to-one matching of finite/pole/log source terms incomplete')
    require(len(records) * 1 == source_counts['pole'] == source_counts['scale_log'] ==
            source_counts['constant'],
            'Native Laurent triplet counts inconsistent')

    return {
        'format': FORMAT,
        'model': model,
        'status': 'EXACT_NATIVE_POLE_CONSTANT_LOG_TRIPLETS_VERIFIED_NOT_RENORMALIZED',
        'source_raw_inputform_sha256': inventory['native_raw_inputform_sha256'],
        'source_inventory_terms_sha256': None,  # caller fills with on-disk digest
        'source_term_count': len(terms),
        'exact_kernel_groups': len(records),
        'PL_groups': sum(r['chirality'] == 'PL' for r in records),
        'PR_groups': sum(r['chirality'] == 'PR' for r in records),
        'explicit_1_over_epsilon_terms': source_counts['pole'],
        'finite_constant_terms': source_counts['constant'],
        'finite_scale_log_terms': source_counts['scale_log'],
        'all_source_signs_and_complex_flavour_structures_retained': True,
        'raw_pole_residue_equals_finite_constant_equals_log_coefficient': True,
        'physical_MSbar_hard_coefficient': None,
        'master_Weyl_phase': None,
        'EFT1_soft_subtraction_performed': False,
        'counterterm_subtraction_performed': False,
        'records': records,
    }


def expressions_for(record: dict, chirality: str) -> dict[str, str]:
    """Exact native InputForm *formal* combinations, never physical C5."""
    require(chirality in ('PL', 'PR'), 'Invalid chirality')
    selected = [r for r in record['records'] if r['chirality'] == chirality]
    require(len(selected) == record[f'{chirality}_groups'], 'Invalid grouped chirality count')
    def sum_native(items: list[str]) -> str:
        return ' '.join(items) if items else '0'
    kernels = [r['kernel_inputform'] for r in selected]
    log_terms = []
    for r in selected:
        # Keep native log *inside* each exact kernel at its verified position:
        # reconstruct from original term by relying on the validated source
        # indices separately (stored in native terms.json). Here a product
        # kernel * Log is mathematically identical (scalar factor commutes).
        signed = r['kernel_inputform']
        log_terms.append(signed[:1] + ' (' + signed[2:] + ')*' + r['mass_log_inputform'])
    return {
        'pole_residue_raw': sum_native(kernels),
        'finite_constant_raw': sum_native(kernels),
        'scale_log_raw': sum_native(log_terms),
        'formal_pole_drop_NOT_MSbar': sum_native(kernels + log_terms),
    }


def audit_all(repo_root: Path, native_dir: Path, inventory_dir: Path) -> tuple[dict, dict[str, dict]]:
    """Reconstruct and compare *every* inventory against hashed native InputForm."""
    root = repo_root.resolve()
    owned = (root / 'output/general_observable/lo_certification_spinor').resolve()
    native = native_dir.resolve()
    inventories = inventory_dir.resolve()
    require(native.is_relative_to(owned) and inventories.is_relative_to(owned) and
            native.is_dir() and inventories.is_dir(), 'Input must be under owned spinor output')
    native_summary, inputs = validate_batch_directory(native)
    manifest_path = inventories / 'summary.json'
    require(manifest_path.is_file(), 'Inventory summary missing')
    inv_summary = json.loads(manifest_path.read_text('utf-8-sig'))
    require(inv_summary.get('format') == INVENTORY_FORMAT and
            inv_summary.get('status') == 'ALL16_RAW_LLHH_INVENTORIED_NOT_PHYSICAL' and
            inv_summary.get('model_count') == 16 and
            len(inv_summary.get('model_records', [])) == 16,
            'Only the released raw-LLHH inventory is supported')
    require(inv_summary.get('source_native_summary_sha256') == sha_bytes((native / 'summary.json').read_bytes()),
            'Inventory source-native summary has changed')
    require(inv_summary.get('native_gamma_cc_to_master_weyl_factor') is None and
            inv_summary.get('MSbar_subtraction_certified') is False and
            inv_summary.get('physical_C5') is None and
            inv_summary.get('ready_for_physical_EFT1_or_scalar_matching') is False,
            'Unexpected physical certification in inventory')
    require(native_summary['source_release_summary_sha256'] == EXPECTED_SNAPSHOT,
            'Unexpected source-resolved first-threshold snapshot')
    models = {}
    previous = inv_summary['model_records']
    require({x.get('model') for x in previous} == set(MODELS), 'Inventory model coverage invalid')
    for model in MODELS:
        path = inventories / tag(model) / 'terms.json'
        require(path.is_file(), f'{model}: missing terms.json')
        raw = path.read_bytes()
        observed = json.loads(raw.decode('utf-8-sig'))
        expected = create_inventory(model, inputs[model])
        require(observed == expected, f'{model}: inventory disagrees with authenticated native source')
        header = {k: v for k, v in observed.items() if k != 'selected_terms'}
        require(header == next(x for x in previous if x['model'] == model),
                f'{model}: inventory summary header mismatch')
        res = validate_and_group(model, observed)
        res['source_inventory_terms_sha256'] = sha_bytes(raw)
        models[model] = res
    totals = Counter()
    for m in MODELS:
        r = models[m]
        totals['kernels'] += r['exact_kernel_groups']
        totals['PL'] += r['PL_groups']
        totals['PR'] += r['PR_groups']
        totals['pole'] += r['explicit_1_over_epsilon_terms']
        totals['finite_constant'] += r['finite_constant_terms']
        totals['scale_log'] += r['finite_scale_log_terms']
    require(sum((totals['pole'], totals['finite_constant'], totals['scale_log'])) ==
            inv_summary['totals']['candidates'] and
            totals['PL'] * 3 == inv_summary['totals']['PL'] and
            totals['PR'] * 3 == inv_summary['totals']['PR'] and
            totals['pole'] == inv_summary['totals']['pole'] and
            totals['finite_constant'] + totals['scale_log'] == inv_summary['totals']['nonpole'],
            'All-16 inventory totals disagree with matched Laurent triplets')
    summary = {
        'format': FORMAT,
        'status': 'ALL16_EXACT_NATIVE_LAURENT_TRIPLETS_AUDITED_PHYSICAL_BLOCKED',
        'models': len(models),
        'mu_F_gev': native_summary['mu_F_gev'],
        'native_summary_sha256': sha_bytes((native / 'summary.json').read_bytes()),
        'raw_inventory_summary_sha256': sha_bytes(manifest_path.read_bytes()),
        'source_release_summary_sha256': EXPECTED_SNAPSHOT,
        'totals': dict(totals),
        'model_records': [{k: v for k, v in models[m].items() if k != 'records'} for m in MODELS],
        'gamma_cc_to_master_weyl_factor': None,
        'physical_C5': None,
        'MSbar_hard_matching_certified': False,
        'EFT1_soft_subtraction_performed': False,
        'counterterms_verified': False,
        'both_heavy_scalar_multiplets_remain_active': True,
        'SM_fermions_remain_active': True,
        'both_scalars_integrated_only_together_at_later_mu_S': True,
        'meaning_of_formal_pole_drop': 'Bare syntactic discard of 1/epsilon only: NOT a renormalized hard coefficient or SMEFT C5',
    }
    return summary, models
