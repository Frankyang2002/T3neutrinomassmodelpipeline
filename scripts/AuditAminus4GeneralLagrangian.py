"""Read-only A(-4) UV->EFT1 GENERAL operator-inventory audit.

Run from the repository root:
  python -X utf8 -m scripts.AuditAminus4GeneralLagrangian

This is an *inventory* check, not a Matchete/RGBeta Clebsch-Gordan map,
GammaCC-to-Weyl certificate, or dimension-five operator-closure proof.
The audit fails closed if the independently validated two-component Higgs-S2
portal basis is not retained in the matched scalar EFT1 expression.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (ROOT / 'output/general_observable/'
                  'fermion_first_all16_v1/A_m4')
DEFAULT_BETA = (ROOT / 'output/studies/general_scalar_extended/'
                'complete_scalar_symbolic/full16_general_v1/A_m4/eft1_rgbeta.json')
DEFAULT_OUTPUT = ROOT / 'output/eft1_general_lagrangian_audit/A_m4.json'

# Determined from the source's independently computed A(-4) SU(2) basis:
# lambdaH2 has two independent invariants; S2 self has two.
EXPECTED_QUARTIC_SYMBOLS = frozenset((
    'λ', 'lambdaS1Inv1', 'lambdaS2Inv1', 'lambdaS2Inv2',
    'lambdaH1Inv1', 'lambdaH2Inv1', 'lambdaH2Inv2',
    'lambda12Inv1', 'lambdaT3',
))
EXPECTED_EFT1_KEYS = frozenset((
    'gY','g2','g3','yu','yd','ye','mS1Sq','mS2Sq',
    'lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2',
    'lambda12','lambdaT3','lambdaH2Adj','lambdaS2Adj',
))


def _read_json(path: Path) -> dict[str, Any]:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return obj


def audit_source(source_dir: Path, beta_path: Path) -> dict[str, Any]:
    """Validate actual saved Matchete tree inventory against EFT1 RGBeta.

    A successful build request or a successful matching STATUS is not enough;
    independently checked invariant families must survive in the EFT1 tree.
    """
    source_dir = Path(source_dir)
    beta_path = Path(beta_path)
    status = _read_json(source_dir/'native_status.json')
    wilson = _read_json(source_dir/'data/eft1_after_F_wilson_seed.json')
    quartic = _read_json(source_dir/'data/eft1_after_F_scalar_quartic_seed.json')
    beta = _read_json(beta_path)
    native_log = (source_dir/'native.log').read_text(encoding='utf-8-sig')
    matched = (source_dir/'data/eft1_after_F_tree_inputform.txt').read_text(encoding='utf-8-sig')

    # Match exact Coupling[NAME, ...] rather than incidental substring hits.
    present = set(re.findall(r'Coupling\[([^,\]]+)', matched))
    missing = sorted(EXPECTED_QUARTIC_SYMBOLS - present)
    metadata = beta.get('metadata', {})
    beta_keys = set(beta.get('report_betas', {}))
    expected_beta_metadata = {
        'GeneralScalarModel': 'A_-4', 'GeneralScalarStage': 'EFT1',
        'SharedScalar': False, 'dS1': 1, 'dS2': 3, 'dF': 2,
        'alpha': -4, 'GeneralScalarAdditionalQuartics': 0,
        'GeneralScalarNewBetaKeys': [], 'QuadraticKind': 'none',
    }
    beta_meta_ok = all(metadata.get(k) == v for k,v in expected_beta_metadata.items())
    beta_keys_ok = beta_keys == EXPECTED_EFT1_KEYS
    matched_fields_ok = (wilson.get('IntegratedFields') == ['F']
                         and wilson.get('ActiveHeavyFields') == ['S1','S2']
                         and quartic.get('IntegratedFields') == ['F']
                         and quartic.get('ActiveHeavyFields') == ['S1','S2'])
    declared_ok = (status.get('model') == 'A_-4'
                   and status.get('operator_basis_is_full_general_requested') is True
                   and status.get('full_general_extra_yukawa_count') == 0
                   and status.get('full_general_quartic_count') == 0
                   and status.get('full_general_quadratic_count') == 0)
    portal_validated = all(
        re.search(r'Checking candidate: HiggsPortal2Inv%s\s+result: InputForm\[True\]' % k,
                  native_log) is not None for k in ('1','2'))
    portal_basis_two = re.search(
        r'lambdaH2: raw SU\(2\) invariants=2, physical invariants=2', native_log
    ) is not None
    beta_without_F = not ({'MF','y1','y2','y1Conj','y2Conj'} & beta_keys)
    native_success = (status.get('status') == 'FThresholdTreeWilsonExported'
                      and beta.get('status') == 'SuccessCandidateMassBeta'
                      and matched_fields_ok and declared_ok)
    source_complete = (native_success and portal_validated and portal_basis_two
                       and beta_meta_ok and beta_keys_ok and beta_without_F
                       and not missing)
    result = {
        'model': 'A_-4',
        'status': ('GeneralRenormalizableEFT1InventoryPresent'
                   if source_complete else 'Blocked_GeneralEFT1InventoryNotEstablished'),
        'general_lagrangian_input_requested': bool(declared_ok),
        'matchete_matching_completed': bool(native_success),
        'independent_HS2_portal_invariants': 2 if portal_basis_two else None,
        'both_portal_candidates_validated': bool(portal_validated),
        'all_expected_quartics_in_matched_tree': not missing,
        'missing_matched_tree_quartic_symbols': missing,
        'matched_tree_quartic_symbols': sorted(EXPECTED_QUARTIC_SYMBOLS & present),
        'eft1_rgbeta_all_17_keys': bool(beta_keys_ok),
        'eft1_rgbeta_expected_metadata': bool(beta_meta_ok),
        'eft1_rgbeta_no_heavy_fermions': bool(beta_without_F),
        'full_general_renormalizable_eft1_source_ready': bool(source_complete),
        'physical_matchete_to_rgbeta_basis_map_certified': False,
        'physical_gamma_cc_to_weyl_factor_certified': False,
        'full_dimension_five_operator_closure_verified': False,
        'source_paths': {
            'matchete': str(source_dir),
            'rgbeta': str(beta_path),
        },
        'note': 'Full-general request is not proof of full retained matched operator inventory. '
                'This audit does not certify physical C5 even when its inventory passes.',
    }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--beta', type=Path, default=DEFAULT_BETA)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    result = audit_source(args.source_dir, args.beta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n', encoding='utf-8')
    print('General EFT1 matched inventory:', result['status'])
    print('Missing matched-tree quartics:', result['missing_matched_tree_quartic_symbols'])
    print('Audit:', args.output)
    return 0 if result['full_general_renormalizable_eft1_source_ready'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
