"""Read-only preflight: verify v19/v20 hashes and B/C final mass-basis sources.

This is source metadata verification, not a new mass eigenbasis and not physical
matching. It supports both path separators used by saved v19 input_sha256 on
Windows and POSIX platforms, while rehashing every input byte-for-byte.

    python -m scripts.RunT3ParallelABCMassProvenanceV22 --project-root .
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.bc_mass_eigen_native_bridge_v21 import (
    BC, BCBridgeBlocked, read_json, read_mass_basis, digest, require,
    verify_v19_input_sha256,
    V19_FORMAT, V19_STATUS, V20_FORMAT, V20_STATUS,
)

BASE = Path('output/parallel_t3/wilson')


def preflight(root: Path) -> dict:
    root = Path(root).resolve()
    topology = root / BASE / 'matchete_topology_native_v18.json'
    bc = root / BASE / 'matchete_bc_extra_native_v18.json'
    v19_path = root / BASE / 'complex_scalar_bridge_conditional_v19.json'
    v20_path = root / BASE / 'full_scalar_census_conditional_v20.json'
    v19, v20 = read_json(v19_path), read_json(v20_path)
    require(v19.get('format') == V19_FORMAT and v19.get('status') == V19_STATUS and
            v19.get('baseline_EFT1_complex_slots') == 28 and
            v19.get('certification', {}).get('physical_Wilson_evolution_authorized') is False,
            'preflight: incompatible v19 complex bridge')
    source = verify_v19_input_sha256(root, v19, (topology, bc))
    require(v20.get('format') == V20_FORMAT and v20.get('status') == V20_STATUS and
            v20.get('total_EFT1_scalar_quartic_registration_slots') == 178 and
            v20.get('input_sha256', {}).get('complex_bridge_v19') == digest(v19_path) and
            v20.get('certification', {}).get('physical_Wilson_evolution_authorized') is False,
            'preflight: v20 178-slot census or its v19 SHA256 is stale')
    masses = {}
    for tag in BC:
        path = root / 'output/general_t3_mass_basis' / f'{tag[0]}_m1_final_mass_basis.json'
        basis, saved = read_mass_basis(root, tag, path)
        require(not basis.tachyonic and len(saved.get('source_sha256', {})) >= 1,
                'preflight: incompatible mass-basis source for ' + tag)
        masses[tag] = {
            'squared_masses_gev2': [float(x) for x in basis.eigenvalues_sq],
            'saved_source_files_sha256_rechecked': len(saved['source_sha256']),
        }
    return {'v19': source, 'v20_v19_sha256_chain_verified': True,
            'mass_sources': masses, 'physical_matching_authorized': False}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=Path('.'))
    args = parser.parse_args(argv)
    try:
        result = preflight(args.project_root)
        print('v21/v22 provenance: v19 input SHA256 {}/7; legacy Windows path keys {}; '
              'v20-to-v19 hash verified'.format(
                  result['v19']['source_files_sha256_rechecked'],
                  result['v19']['legacy_windows_separator_keys_accepted']))
        for model, row in result['mass_sources'].items():
            print('  {}: saved B/C final mass eigenvalues {} GeV^2; {} source hashes verified'.format(
                model, row['squared_masses_gev2'],
                row['saved_source_files_sha256_rechecked']))
        print('Read-only provenance PASS; does not certify Matchete matching or Wilson evolution')
        return 0
    except (BCBridgeBlocked, ValueError, KeyError, TypeError, OSError) as exc:
        print('v21/v22 provenance BLOCKED:', exc, file=sys.stderr)
        return 7


if __name__ == '__main__':
    raise SystemExit(main())
