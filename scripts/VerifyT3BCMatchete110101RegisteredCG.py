"""Compare *registered* live Matchete B/C 110101 CG components to Stage22.

This completes the missing live-CG component gate after the separately
verified Matchete action and RGBeta Lam reports.  It verifies exact component
identity and the signed L_M=-V_RGBeta polynomial in the same exported SU(2)
component basis, not physical mass-basis EFT matching or a live action
subtraction in one combined Mathematica session.

Run from the repository root:
 python -u -X utf8 -m scripts.VerifyT3BCMatchete110101RegisteredCG \
  --input output/audits/t3_bc_matchete_registered_cg/live_registered_cg_v1.json \
  --run-id bc_live_registered_cg_v1
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

import sympy as sp
import numpy as np

from scripts.AuditT3BC110101MassBasisMap import (
    MODELS, MAT, NATIVE, RGB_NATIVE, RGB_LIVE, ROOT, check_live,
    exact_map, exact_number, load, original_matchete_records,
)
from scripts.AuditT3BCLive110101ActionBridge import exact_action_identity
from scripts.VerifyT3BCNativeLam110101Coefficients import verify as verify_lam
from scripts.VerifyT3BCMatchete110101Action import (
    verify as verify_matchete, verify_source as verify_matchete_source,
)

FORMAT = 't3-bc-live-matchete-110101-registered-cg-v1'
STATUS = 'LiveRegisteredCGExtractionPendingIndependentComparison'
SOURCE = 'installed Matchete GetCGTensor called after T3ExtendedBuildExtraQuarticCandidates'
CG_NAMES = ('T3ExtraCGQ110101I1', 'T3ExtraCGQ110101I2')
MAT_LIVE = ROOT / 'output/audits/t3_bc_matchete_action/matchete_110101_v2.json'
RGB_LIVE_LAM = ROOT / 'output/audits/t3_bc_native_lam/native_lam_coeff_v4.json'
BRIDGE = ROOT / 'output/audits/t3_bc_110101_live_action_bridge/bc_live_action_bridge_v1/summary.json'
MAT_SOURCE = ROOT / 'Lagrangian/interactions/GeneralScalarOperators.wl'
OUT = ROOT / 'output/audits/t3_bc_matchete_registered_cg'


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('B/C live Matchete registered CG audit: ' + message)


def hash_bytes(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def validate_probe(doc: dict) -> dict:
    require(isinstance(doc, dict) and doc.get('format') == FORMAT and
            doc.get('status') == STATUS and doc.get('source') == SOURCE,
            'unexpected native Matchete CG probe format/status/source')
    for key in ('source_basis_tensors_substituted',
                'direct_live_cross_tool_action_identity_verified',
                'physical_coupling_mapping_certified',
                'mass_basis_matching_certified', 'sequential_eft_matching_certified',
                'production_modified'):
        require(doc.get(key) is False, 'unwarranted probe claim: ' + key)
    models = doc.get('models')
    require(isinstance(models, dict) and set(models) == set(MODELS),
            'expected B_-1 and C_-1')
    for model in MODELS:
        row = models[model]
        require(isinstance(row, dict) and row.get('model') == model and
                row.get('status') == 'LiveRegisteredMatcheteCGComponentsCaptured' and
                row.get('extra_quartic_inventory') == 7 and
                row.get('registered_cg_components_extracted') == 32,
                model + ': blocked or incomplete registered CG extraction')
        for key in ('live_registered_cg_components_verified',
                    'direct_live_cross_tool_action_identity_verified',
                    'physical_matching_certified'):
            require(row.get(key) is False, model + ': unwarranted CG claim ' + key)
        rows = row.get('registered_tensors')
        require(isinstance(rows, list) and len(rows) == 2,
                model + ': expected two registered CG tensors')
        for i, record in enumerate(rows):
            require(isinstance(record, dict) and
                    record.get('symbol') == CG_NAMES[i] and
                    record.get('matchete_coupling') == MAT[i] and
                    record.get('term_name') == f'T3ExtraQ_1_1_0_1_0_1_{i+1}' and
                    record.get('cg_present_in_native_plus_hc_interaction') is True and
                    record.get('tensor_accessor') == 'Matchete`GetCGTensor' and
                    record.get('dimensions') == [2, 2, 2, 2] and
                    record.get('component_count') == 16 and
                    isinstance(record.get('components'), list) and
                    len(record['components']) == 16,
                    model + f': invalid registered CG {i+1} metadata')
            require(record.get('returned_head') in ('List', 'SparseArray'),
                    model + ': unexpected tensor return head')
            for value in record['components']:
                exact_number(value)
    return models


def verify_bridge(doc: dict, source_paths: dict) -> None:
    require(doc.get('format') == 't3-bc-live-evidence-gated-110101-action-identity-v1' and
            doc.get('status') == 'ExactExportedComponentActionIdentityWithSeparateLiveGates' and
            doc.get('total_direct_and_hc_coefficients_checked') == 128 and
            doc.get('live_matchete_interactions_checked') == 4 and
            doc.get('live_rgbeta_lam_coefficients_checked') == 128 and
            doc.get('direct_live_matchete_cg_component_equality_verified') is False and
            doc.get('physical_coupling_mapping_certified') is False,
            'previous exported action bridge is missing or not approved')
    rows = doc.get('models')
    require(isinstance(rows, list) and len(rows) == 2 and
            [r.get('model') for r in rows] == list(MODELS) and
            all(r.get('exact_action_identity_in_exported_su2_component_basis') is True and
                r.get('matchete_live_cg_numerical_components_independently_extracted') is False
                for r in rows),
            'previous bridge did not establish exact exported component identity')
    source_hashes = doc.get('source_sha256')
    require(isinstance(source_hashes, dict), 'previous bridge input hashes missing')
    for label, path in source_paths.items():
        require(source_hashes.get(str(path.resolve())) == hash_bytes(path),
                'stale previous bridge source ' + label)


def compare_one(model: str, probe: dict, native: dict, invariants: dict,
                covariance: dict) -> dict:
    reference = original_matchete_records(native, model)
    actual = probe[model]['registered_tensors']
    m_columns = []
    for i, record in enumerate(actual):
        observed = [exact_number(z) for z in record['components']]
        expected = [exact_number(z) for z in reference[i]]
        residuals = [sp.simplify(a-b) for a,b in zip(observed,expected)]
        require(all(x == 0 for x in residuals),
                model + '/' + CG_NAMES[i] + ': live registered CG differs from Stage22 native tensor')
        m_columns.append(sp.Matrix(observed))
    M = sp.Matrix.hstack(*m_columns)
    rgb_tensors = check_live(covariance, invariants, model)
    A = exact_map(native, model, rgb_tensors)
    R = sp.Matrix.hstack(*(
        sp.Matrix([sp.Integer(int(round(float(np.real(z))))) for z in t.flat])
        for t in rgb_tensors
    ))
    action = exact_action_identity(M, R, A)
    require(action['exact_nonzero_component_residuals'] == 0 and
            action['incorrect_sign_rejected'] is True and
            action['spurious_extra_quarter_rejected'] is True,
            model + ': action coefficients/control check failed')
    return {
        'model': model,
        'status': 'ExactLiveMatcheteRegisteredCGAndRGBetaExportedInvariantIdentity',
        'live_matchete_registered_tensors_checked': 2,
        'live_registered_components_equal_native_stage22_exact': 32,
        'live_matchete_registered_cg_direct_and_hc_coefficients_checked': 64,
        'exact_registered_component_disagreements': 0,
        'signed_cross_tool_component_identity_verified': True,
        'live_matchete_registered_cg_extraction_passed': True,
        'native_rgbeta_lam_previously_verified': True,
        'matchete_live_interaction_validation_previously_verified': True,
        'direct_live_matchete_action_subtraction_performed': False,
        'mass_basis_matching_certified': False,
        'sequential_eft_matching_certified': False,
        'physical_coupling_mapping_certified': False,
        **action,
    }


def run(probe_path: Path, mat_live_path: Path, rgb_lam_path: Path,
        bridge_path: Path, native_path: Path, invariant_path: Path,
        covariance_path: Path, mat_source_path: Path) -> dict:
    inputs = {'probe': probe_path,'mat_live': mat_live_path,
              'rgb_lam': rgb_lam_path,'bridge': bridge_path,
              'native':native_path,'invariants': invariant_path,
              'covariance':covariance_path,'mat_source':mat_source_path}
    for name, path in inputs.items():
        require(path.is_file(), 'missing input ' + name + ': ' + str(path))
    verify_matchete_source(mat_source_path)
    vm = verify_matchete(load(mat_live_path))
    vr = verify_lam(load(rgb_lam_path))
    require(vm['matchete_native_terms_validated'] == 4 and
            vm['matchete_combined_checks_passed'] == 2 and
            vr['total_components_checked'] == 128,
            'live Matchete or RGBeta prerequisite gate failed')
    bridge_sources = {k: inputs[k] for k in
                      ('mat_live','rgb_lam','native','invariants','covariance','mat_source')}
    verify_bridge(load(bridge_path), bridge_sources)
    probe = validate_probe(load(probe_path))
    native, invariants, covariance = (load(inputs[x]) for x in
                                    ('native','invariants','covariance'))
    rows = [compare_one(model, probe, native, invariants, covariance)
            for model in MODELS]
    return {
        'format': 't3-bc-live-registered-cg-110101-cross-tool-identity-v1',
        'status': 'LiveMatcheteRegisteredCGExactWithGatedRGBetaIdentityNoEFTMatching',
        'input_sha256': {name: hash_bytes(path) for name,path in inputs.items()},
        'models': rows,
        'live_registered_cg_tensors_checked': 4,
        'live_registered_matchete_components_compared': 64,
        'signed_direct_and_conjugate_coefficients_compared': 128,
        'live_matchete_registered_cg_component_identity_verified': True,
        'cross_tool_signed_component_identity_verified_in_exported_su2_basis': True,
        'direct_live_matchete_action_subtraction_performed': False,
        'mass_basis_matching_certified': False,
        'sequential_eft_matching_certified': False,
        'physical_coupling_mapping_certified': False,
        'production_modified': False,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--matchete-live', type=Path, default=MAT_LIVE)
    p.add_argument('--rgbeta-live', type=Path, default=RGB_LIVE_LAM)
    p.add_argument('--bridge', type=Path, default=BRIDGE)
    p.add_argument('--matchete-native', type=Path, default=NATIVE)
    p.add_argument('--rgbeta-invariants', type=Path, default=RGB_NATIVE)
    p.add_argument('--rgbeta-covariance', type=Path, default=RGB_LIVE)
    p.add_argument('--matchete-source', type=Path, default=MAT_SOURCE)
    p.add_argument('--run-id', required=True)
    p.add_argument('--output-root', type=Path, default=OUT)
    args = p.parse_args()
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id) is not None,
            'invalid run ID')
    output = args.output_root / args.run_id / 'summary.json'
    require(not output.exists(), 'refusing to overwrite ' + str(output))
    result = run(args.input, args.matchete_live, args.rgbeta_live,
                 args.bridge, args.matchete_native, args.rgbeta_invariants,
                 args.rgbeta_covariance, args.matchete_source)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for row in result['models']:
        print(f"{row['model']}: registered Matchete CG components 32/32 exact; "
              'signed direct+h.c. coefficients 64/64')
    print('Direct live Matchete action subtraction and physical EFT matching: NOT CERTIFIED')
    print('Summary:',output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
