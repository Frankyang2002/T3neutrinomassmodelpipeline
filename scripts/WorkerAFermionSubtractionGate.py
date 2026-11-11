"""Read-only, source-hashed UV -> EFT1 fermion-threshold subtraction gate.

Produces an operator-by-operator UV source inventory and tests whether existing
Matchete/Laurent exports are eligible for finite matching. It cannot generate
missing EFT1-soft graphs, counterterms or a physical Wilson coefficient.

Run:
 python scripts/WorkerAFermionSubtractionGate.py --sources-zip FILE --output FILE
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import zipfile

MODELS = ('A_-4', 'D_-2')
SECTORS = ('CandidateLeptonNumberViolatingLLScalarScalar',
           'CandidateLeptonNumberViolatingLLHH', 'OtherEFT1Operator')
FORBIDDEN_YUKAWAS = frozenset(('y1Conj', 'y2Conj'))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def choose(names: list[str], fragment: str, suffix: str) -> str:
    matching = [name for name in names if fragment in name and name.endswith(suffix)]
    require(len(matching) == 1, f'Expected exactly one {fragment} {suffix}, found {matching}')
    return matching[0]


def checked_json(z: zipfile.ZipFile, path: str) -> dict:
    obj = json.loads(z.read(path))
    require(isinstance(obj, dict), f'{path}: JSON root must be object')
    return obj


def source_records(z: zipfile.ZipFile, model: str) -> dict:
    names = z.namelist()
    base = f'thesis_full_general_batch2_raw_v1/{model}/'
    ledger_path = choose(names, base, 'one_loop_operator_ledger.json')
    raw_path = choose(names, base, 'native_matching/data/eft1_after_F_one_loop_inputform.txt')
    laurent_path = choose(names, f'thesis_full_general_batch2_laurent_v2/{model}/',
                          'laurent_audit.json')
    ledger, laurent = checked_json(z, ledger_path), checked_json(z, laurent_path)
    raw_hash = sha256(z.read(raw_path)).hexdigest()
    require(ledger.get('source_sha256') == raw_hash, f'{model}: ledger source SHA mismatch')
    require(laurent.get('source_sha256') == raw_hash, f'{model}: Laurent source SHA mismatch')
    require(ledger.get('status') == 'UnrenormalisedOneLoopTermsEnumerated',
            f'{model}: unexpected raw-source status')
    require(ledger.get('only_fermions_integrated_out') is True, f'{model}: not F-only')
    require(ledger.get('exact_reconstruction_by_expansion') is True,
            f'{model}: incomplete raw source reconstruction')
    terms = ledger['terms']
    require(len(terms) == ledger['term_count'], f'{model}: ledger incomplete')
    require([t['index'] for t in terms] == list(range(1, len(terms)+1)),
            f'{model}: source term positions not consecutive')
    sectors = Counter(t['sector'] for t in terms)
    require(set(sectors) <= set(SECTORS), f'{model}: unknown operator class')
    require(dict(sectors) == ledger['sector_counts'], f'{model}: sector counter mismatch')
    statuses = []
    inventory = {}
    for sector in SECTORS:
        subset = [t for t in terms if t['sector'] == sector]
        permitted, excluded = [], []
        for t in subset:
            couplings = set(t['coupling_symbols'])
            record = {'source_term_index': t['index'],
                      'source_term_sha256': sha256(t['raw_inputform'].encode('utf-8')).hexdigest(),
                      'couplings': sorted(couplings),
                      'contains_hbar': t['contains_hbar'],
                      'contains_explicit_pole_heuristic': t['contains_dim_reg_pole']}
            if couplings & FORBIDDEN_YUKAWAS:
                record['excluded_yukawas'] = sorted(couplings & FORBIDDEN_YUKAWAS)
                excluded.append(record)
            else:
                permitted.append(record)
        inventory[sector] = {'total': len(subset), 'restricted_candidate_count': len(permitted),
                             'excluded_conjugate_yukawa_count': len(excluded),
                             'restricted_candidates': permitted, 'excluded_terms': excluded}
    require(sectors['CandidateLeptonNumberViolatingLLHH'] == (6 if model == 'A_-4' else 12),
            f'{model}: direct LLHH source inventory changed')
    expected_restricted = 6
    require(len(inventory['CandidateLeptonNumberViolatingLLHH']['restricted_candidates']) == expected_restricted,
            f'{model}: restricted direct LLHH inventory changed')
    criteria = {
        'source_SHA256_consistent': True,
        'full_original_ledger_reconstructed': True,
        'all_raw_operators_retained_in_audit': ledger.get('unselected_terms_retained') is True,
        'explicit_Laurent_reconstruction_exact': laurent.get('explicit_laurent_reconstruction_exact') is True,
        'loop_function_evaluation_explicitly_verified': laurent.get('matchete_loop_function_evaluator_called') is True,
        'one_loop_hbar_linearity_exact': laurent.get('one_loop_hbar_linearity_exact') is True,
        'tree_mass_coordinate_remainder_resolved': laurent.get('unreconciled_tree_mass_coordinate_identity_verified') is True,
        'UV_minus_EFT1_soft_subtraction_performed': laurent.get('uv_hard_minus_eft1_soft_subtraction_performed') is True,
        'local_MSbar_counterterms_added': laurent.get('local_MSbar_counterterms_added') is True,
        'complete_flavour_contractions_performed': laurent.get('complete_numerical_flavour_contractions_performed') is True,
        'finite_renormalised_boundary_computed': laurent.get('renormalised_MSbar_Wilson_boundary_computed') is True,
        'matchete_to_RGBeta_dictionary_applied': ledger.get('matchete_to_rgbeta_basis_converted') is True,
    }
    status = 'BLOCKED_NO_RENORMALISED_F_THRESHOLD' if not all(criteria.values()) else 'SOURCE_METADATA_READY_FOR_INDEPENDENT_REVIEW'
    return {
        'model': model, 'status': status, 'source_sha256': raw_hash,
        'sources': {'ledger': ledger_path, 'raw_one_loop': raw_path, 'laurent_audit': laurent_path},
        'criteria': criteria, 'blocking_criteria': [k for k,v in criteria.items() if not v],
        'raw_laurent_term_counts': {'pole': laurent.get('raw_pole_residue_term_count'),
                                   'regular': laurent.get('raw_regular_term_count')},
        'operator_inventory': inventory,
        'warning': ('Restricted candidates are NOT a valid EFT1 action by themselves: '
                    'conjugate-Yukawa terms are only zero on the specified restricted '
                    'parameter surface, and counterterms/operator mixing must be recomputed accordingly.'),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources-zip', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    args = p.parse_args()
    require(args.sources_zip.is_file(), f'Missing source ZIP: {args.sources_zip}')
    with zipfile.ZipFile(args.sources_zip) as z:
        results = {model: source_records(z, model) for model in MODELS}
    report = {'schema':'t3-worker-a-fermion-subtraction-source-gate-v1',
              'matching_scale_GeV':100000,
              'integrate_out':['all heavy fermions simultaneously'],
              'retain_in_EFT1':['S1','S2'],
              'scope':'UV -> EFT1 fermion threshold; not joint scalar threshold',
              'finite_matching_arrays_written':False,
              'model_results':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Written:', args.output)
    for model, v in results.items():
        x=v['operator_inventory']
        print(f"{model}: LLSS {x[SECTORS[0]]['total']} (restricted {x[SECTORS[0]]['restricted_candidate_count']}), "
              f"LLHH {x[SECTORS[1]]['total']} (restricted {x[SECTORS[1]]['restricted_candidate_count']}), "
              f"blocked criteria {len(v['blocking_criteria'])}")
    require(all(r['status'].startswith('BLOCKED') for r in results.values()),
            'Source unexpectedly passes all gates; independently audit before using.')


if __name__ == '__main__':
    main()
