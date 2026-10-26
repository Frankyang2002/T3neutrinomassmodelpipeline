"""Stage 25 read-only check of live RGBeta complex-quartic assembly.

The Wolfram companion calls RGBeta's genuine CouplingTensorMap, scalar index
map and complex scalar contraction. This verifies the *native assembly* used
inside Lam, not the sign relating a Lagrangian to a potential, and not a global
Matchete-to-RGBeta Cartesian real-field intertwiner. Fail closed.

python -m scripts.Stage25QuarticActionAudit
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from scripts.Stage23QuarticPrefactorAudit import parse_permutations, run as stage23_run

FOLDER = Path('output/general_t3_stage22')
LIVE = 'stage25_live_native_quartic_assembly.json'
DECL = 'rgbeta_native_declarations.json'
FORMAT = 'stage25-live-rgbeta-quartic-assembly-v1'
EXPECTED_STATUS = 'NativeComplexQuarticAssemblyCheckedNotPhysicalPotentialSignCertificate'
NEVER_CERTIFIED = (
    'native_complex_species_to_cartesian_real_intertwiner_verified',
    'matchete_relative_potential_sign_verified',
    'physical_coupling_mapping_certified',
    'production_modified',
)


def fail(message: str):
    raise ValueError('Stage25 fail-closed: ' + message)


def compact(s):
    return re.sub(r'\s+', '', s)


def verify(live: dict, declarations: dict, stage23: dict) -> dict:
    if live.get('format') != FORMAT or live.get('status') != EXPECTED_STATUS:
        fail('unexpected live Wolfram source or output format')
    if live.get('models_checked') != 16 or live.get('slots_checked') != 50:
        fail('incomplete live model/slot coverage')
    if any(live.get(k) is not False for k in NEVER_CERTIFIED):
        fail('native assembly evidence must not assert physical certification or production changes')
    if live.get('all_native_permutation_reconstructions_exact') is not True or \
            live.get('all_complex_conjugate_metric_two_verified') is not True:
        fail('runtime native assembly or complex scalar metric was not verified')
    dm = declarations.get('models', {})
    lm = live.get('models', {})
    if not isinstance(dm, dict) or set(dm) != set(lm) or len(dm) != 16:
        fail('live and native declaration model sets differ')
    if stage23.get('source_prefactors_verified') is not True or \
            stage23.get('signed_physical_coupling_map_certified') is not False:
        fail('Stage23 is not the conditional verified source report')
    if stage23.get('totals', {}).get('slots') != 50 or stage23.get('totals', {}).get('blocks') != 46:
        fail('Stage23 whole-model totals changed')
    if not isinstance(stage23.get('blocks'), list) or len(stage23['blocks']) != 46:
        fail('expected exactly 46 Stage23 conditional transformation blocks')
    stage23_perms = {}
    for block in stage23['blocks']:
        model = block['model']
        counts = block['native_rgbeta_permutation_counts']
        if not isinstance(counts, dict) or set(counts) != set(block['rgbeta_couplings']):
            fail('Stage23 block permutation metadata malformed')
        for symbol, count in counts.items():
            key = (model, symbol)
            if key in stage23_perms:
                fail('duplicate Stage23 RGBeta quartic')
            stage23_perms[key] = count
    validated = 0
    permutation_arrangements = 0
    across = {}
    for model in sorted(dm):
        rows = lm[model]
        native = dm[model].get('declarations', [])
        if not isinstance(rows, list) or not isinstance(native, list) or len(rows) != len(native):
            fail('live row count differs from Stage22 declarations in ' + model)
        known = {row['symbol']: row for row in native}
        if len(known) != len(native):
            fail('duplicate Stage22 declaration symbol in ' + model)
        seen = set()
        for row in rows:
            symbol = row.get('symbol')
            if symbol not in known or symbol in seen:
                fail('missing or duplicate live RGBeta symbol in ' + model)
            seen.add(symbol)
            d = known[symbol]
            n = len(parse_permutations(d['unique_arrangements_inputform']))
            if row.get('native_permutations') != n or row.get('native_assembly_rules_checked') != n:
                fail('live permutation count differs for ' + model + '/' + symbol)
            if row.get('assembly_equals_independent_permutation_reconstruction') is not True or \
                    row.get('complex_conjugate_metric_two_verified') is not True:
                fail('unverified native quartic assembly for ' + model + '/' + symbol)
            idx = row.get('field_species_positions')
            if not isinstance(idx, list) or len(idx) != 4 or any(type(x) is not int or x <= 0 for x in idx):
                fail('bad native complex species index map for ' + model + '/' + symbol)
            if compact(row.get('native_coupling_value_inputform', '')) != f'6*{symbol}':
                fail('unexpected live quartic coupling factor in ' + model + '/' + symbol)
            if stage23_perms.get((model, symbol)) != n:
                fail('Stage23/native permutation count mismatch: ' + model + '/' + symbol)
            permutation_arrangements += n
            validated += 1
        if seen != set(known):
            fail('live data missed slots in ' + model)
        across[model] = len(seen)
    if validated != 50 or len(stage23_perms) != 50:
        fail('50-slot coverage failed')
    return {
        'format': 'stage25-native-complex-quartic-assembly-verified-v1',
        'status': 'NativeAssemblyVerifiedPhysicalCouplingMapUncertified',
        'models_checked': 16,
        'slots_checked': validated,
        'blocks_covered': 46,
        'native_permutation_rules_checked': permutation_arrangements,
        'live_rgbeta_coupling_tensor_map_checked': True,
        'live_complex_field_metric_two_checked': True,
        'stage23_permutation_alignment_checked': True,
        'models': across,
        'physical_coupling_mapping_certified': False,
        'remaining_requirements': [
            'Independently establish the global Matchete/RGBeta SU2 and Cartesian-real component intertwiner, including triplets',
            'Independently establish the overall relative Lagrangian-to-potential sign in the native RGBeta action',
            'Compare the complete native quartic action, including Hermitian conjugates and identical-field factors, rather than treating CouplingTensorMap alone as physical certification',
        ],
        'production_modified': False,
    }


def run(folder=FOLDER, output=None):
    def load(name):
        p = folder / name
        if not p.is_file():
            fail(f'missing {p}; run the live Stage25 Wolfram export and regenerate Stage22 evidence')
        return json.loads(p.read_text(encoding='utf-8-sig'))
    # Regenerate Stage23 logic in memory and verify its source hashes; do not
    # trust a stale candidate file or inferred physical sign as evidence.
    stage23 = stage23_run(folder, None)
    stored = load('stage23_quartic_prefactor_candidates.json')
    if stored != stage23:
        fail('Stage23 stored report is stale against current Stage22 inputs')
    result = verify(load(LIVE), load(DECL), stage23)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=FOLDER)
    parser.add_argument('--output', type=Path,
                        default=FOLDER/'stage25_native_quartic_assembly_diagnostics.json')
    args = parser.parse_args()
    try:
        status = run(args.input_dir, args.output)
    except (ValueError, OSError, KeyError, TypeError) as e:
        parser.exit(2, str(e) + '\n')
    print('Stage25 live RGBeta assembly:', status['slots_checked'], '/ 50 slots;',
          status['blocks_covered'], '/ 46 blocks')
    print('Native permutation rules checked:', status['native_permutation_rules_checked'])
    print('Global Cartesian basis and relative action sign: UNCERTIFIED')
    print('Physical coupling mapping: UNCERTIFIED')
    print('Report:', args.output)


if __name__ == '__main__':
    main()
