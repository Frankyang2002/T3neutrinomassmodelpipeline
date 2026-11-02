"""Preflight/postflight for the native Matchete A_-4 Feynman-rule probe.

Project path: output/general_observable/lo_fermion_boundary/native_feynman_rules_check.py
Read-only verification of existing 16-model source resolved boundary.
Never assigns a physical GammaCC-to-Weyl normalization factor.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
RELEASE = BASE / 'releases' / 'source_021dcc1e8fe8ba0f'


def verified_native(root=ROOT, release=RELEASE):
    import build_boundary as bound
    root, release = Path(root).resolve(), Path(release).resolve()
    base = root / 'output' / 'general_observable'
    output, _ = bound.compile_model(
        'A_-4', root,
        base / 'ordered_wilson_tensor_v4',
        base / 'general_wilson_polynomial_batch_v2',
        base / 'general_wilson_normalization_audit_v1',
    )
    prior = bound.read(release / 'A_m4' / 'boundary.json')
    assert output['checksums_sha256'] == prior['checksums_sha256'], (
        'A_-4 release does not match the currently verified native sources'
    )
    summary = bound.read(release / 'summary.json')
    assert summary['n_models_verified'] == 16
    ordered = bound.read(base / 'ordered_wilson_tensor_v4' / 'A_m4' / 'ordered_tensor.json')
    source_path = bound.source_file('A_-4', root,
        base / 'ordered_wilson_tensor_v4', ordered)
    source = bound.read(source_path)
    assert bound.digest(source_path) == prior['checksums_sha256']['native_source']
    selected = [term for term in source['operator_terms'] if term['source_index'] == 2]
    assert len(selected) == 1
    term = selected[0]
    assert bound.chirality(term['source_term_inputform']) == 'PL'
    assert [x['name'] for x in term['factors']] == ['y1', 'y2']
    assert term['source_term_inputform'].endswith(
        '/Coupling[MF, {Index[d$$1, NFlavor]}, 0]')
    return term['source_term_inputform'], bound.digest(source_path)


def validate_report(rec, native_term):
    assert rec['format'] == 't3-native-matchete-fr-evidence-v1'
    assert rec['model'] == 'A_-4'
    for name in ('matchete_gamma_cc_to_weyl_phase_certified',
                 'independent_connected_uv_four_point_amplitude_computed',
                 'physical_master_wilson_normalization_certified',
                 'original_16_model_release_changed'):
        assert rec[name] is False, f'{name} incorrectly certified'
    if rec['status'] == 'Blocked':
        assert rec['blocked_phase'] and rec['reason']
        return 'Blocked'
    assert rec['status'] == 'NativeUV3pointAndEFT1ContactVerticesExtracted_AmplitudeOpen'
    assert rec['native_source_index'] == 2
    assert rec['native_source_inputform'] == native_term
    assert rec['uv_y1_with_hc'] == 2 and rec['uv_y2_with_hc'] == 2
    assert rec['uv_bsm_yukawa_terms'] == 4
    assert rec['uv_three_point_rules_extracted'] is True
    assert rec['eft1_four_point_rules_extracted'] is True
    assert rec['heavy_propagator_and_both_exchange_orderings_contracted'] is False
    for field in ('uv_yukawa_lagrangian', 'uv_three_point_rules',
                  'eft1_four_point_rules', 'uv_mass_lagrangian_source'):
        content = rec[field]
        assert content['inputform_characters'] > 0
        assert isinstance(content['inputform_sha256'], str) and len(content['inputform_sha256']) == 64
        if content['inputform_complete']:
            assert isinstance(content['inputform'], str)
            assert len(content['inputform']) == content['inputform_characters']
        else:
            assert content['inputform'] is None and content['inputform_characters'] > 100000
    return 'NativeFeynmanRulesExtracted_NoAmplitudeCertification'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=ROOT)
    parser.add_argument('--release', type=Path, default=RELEASE)
    parser.add_argument('--report', type=Path,
                        help='Validate a native Feynman-rule JSON from wolframscript')
    args = parser.parse_args(argv)
    try:
        term, digest = verified_native(args.repo_root, args.release)
        print('PASS: A_-4 source and 16-model release provenance; source sha256:', digest)
        if args.report is None:
            print('PRE-FLIGHT ONLY: no native Feynman-rule amplitude computed')
            return 0
        import build_boundary as bound
        evidence = bound.read(args.report)
        result = validate_report(evidence, term)
        print('REPORT STATUS:', result)
        if result == 'Blocked':
            print(evidence['blocked_phase'] + ': ' + evidence['reason'])
            return 2
        print('NATIVE FR extracted; connected heavy exchange not yet contracted.')
        print('ABSOLUTE Matchete/Weyl normalization remains UNCERTIFIED.')
        return 0
    except (AssertionError, ValueError, KeyError, TypeError, OSError, IndexError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
