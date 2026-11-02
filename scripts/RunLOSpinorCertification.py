#!/usr/bin/env python3
"""All-16 T3 first-fermion-threshold audit and non-physical EFT1 handoff.

Repository-relative destination: scripts/RunLOSpinorCertification.py
Never writes to pre-existing physics/reference output directories.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
import sys


def evidence_report_filename(evidence: dict) -> str:
    """Use a new evidence-v2 namespace, preserving immutable old reports.

    The path is also raw-source-content-addressed: re-exporting the same
    source is idempotent, whereas newly captured Matchete source evidence does
    not overwrite a parallel worker's prior report.
    """
    raw_hash = evidence.get('source_sha256')
    if isinstance(raw_hash, str) and re.fullmatch(r'[0-9a-f]{64}', raw_hash):
        stamp = raw_hash[:16]
    else:
        stamp = 'not_available'
    return f'spinor_evidence_v2_{stamp}.json'


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--source-release', type=Path,
                        default=Path('output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'))
    parser.add_argument('--output', type=Path,
                        default=Path('output/general_observable/lo_certification_spinor'))
    parser.add_argument('--native-evidence', type=Path,
                        default=Path('output/general_observable/lo_fermion_boundary/matchete_gamma_cc_evidence_v2.txt'))
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--require-physical', action='store_true',
                        help='Require independently certified physical normalization and hard matching; fails closed.')
    args = parser.parse_args(argv)
    root = args.repo_root.resolve()
    sys.path.insert(0, str(root))
    from Numerical.lo_certification.spinor import certify_source_release, convention_probe, exchange_probe
    from Numerical.lo_certification.spinor.boundary import BoundaryError, native_evidence_probe
    def within_root(path):
        return path.resolve() if path.is_absolute() else (root / path).resolve()
    try:
        release = within_root(args.source_release)
        destination = within_root(args.output)
        require_output = root / 'output/general_observable/lo_certification_spinor'
        if not args.check_only and not destination.is_relative_to(require_output):
            raise BoundaryError('Output must stay inside exclusive lo_certification_spinor/ directory')
        # Validate all models BEFORE any output is written.
        manifest = certify_source_release(root, release, None, strict_references=True, strict_counts=True)
        manifest['conditional_spinor_probes'] = convention_probe()
        manifest['conditional_zero_momentum_exchange_probes'] = exchange_probe()
        manifest['installed_matchete_source_evidence'] = native_evidence_probe(within_root(args.native_evidence))
        if args.require_physical:
            raise BoundaryError('BLOCKED: no native-Matchete -> master-Weyl phase certificate, '
                                'no all16 native connected-amplitude certificate, and no fermion hard-one-loop matching')
        if not args.check_only:
            from Numerical.lo_certification.spinor.boundary import checked_write
            # Code checks all inputs again while exporting identical versioned data.
            exported = certify_source_release(root, release, destination, strict_references=True,
                                             strict_counts=True)
            assert exported['source_release_summary_sha256'] == manifest['source_release_summary_sha256']
            release_dest = destination / ('source_' + manifest['source_release_summary_sha256'][:16])
            # The previous spinor_evidence_and_status.json is immutable. The
            # v2 decoder can correctly identify UTF-16 source signatures,
            # changing diagnostic metadata without changing the F-tree tensors.
            # Write a NEW content-addressed v2 evidence report beside it.
            evidence_name = evidence_report_filename(manifest['installed_matchete_source_evidence'])
            checked_write(release_dest / evidence_name, {
                'format': 't3-spinor-native-source-evidence-v2',
                'convention_probe': manifest['conditional_spinor_probes'],
                'exchange_probe': manifest['conditional_zero_momentum_exchange_probes'],
                'installed_matchete_evidence': manifest['installed_matchete_source_evidence'],
                'physical_normalization_certified': False,
            })
            manifest['output_location'] = str(release_dest.relative_to(root))
            manifest['evidence_output_file'] = str((release_dest / evidence_name).relative_to(root))
        print(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False))
        print('VERIFIED: 16/16 provisional source numerical reconstructions and common mu_F; '
              'BLOCKED: physical GammaCC/Weyl phase and complete hard-one-loop matching', file=sys.stderr)
        return 0
    except (BoundaryError, OSError, ValueError, KeyError, AssertionError) as exc:
        print(f'BLOCKED: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
