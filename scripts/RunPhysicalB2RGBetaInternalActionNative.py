#!/usr/bin/env python3
r"""B2 native RGBeta launcher: source-verified Git provenance, no Wolfram Git process.

Run from the project root:
  .\.venv\Scripts\python.exe scripts\RunPhysicalB2RGBetaInternalActionNative.py \
    --output output\general_observable\physical_normalisation\b2_rgbeta_action_native_04\internal_action_registry.json

The subprocess receives only fixed environment variable keys populated from
Python's source-verified local Git snapshot. The underlying Wolfram script
still probes RGBeta's internal tensors, not the physical source-action map.
No production code or tracked file is modified.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'Numerical/lo_certification/physical_normalisation'))

from b2_gauge_coupling import AuditBlocked, SOURCE_BLOBS, require
from b2_rgbeta_ward_witness import source_provenance
from b2_rgbeta_internal_action import (
    ALLOWED_B2_PATHS, native_snapshot_guard, read_report, validate_native_data,
)

WL_SCRIPT = Path('scripts/RunPhysicalB2RGBetaInternalActionProbe.wl')
PROVENANCE_ENV = {
    'RGE/running/rgbeta/T3RGBetaModel.wl': 'T3_B2_PROVENANCE_RGBETA',
    'Lagrangian/interactions/T3Topology.wl': 'T3_B2_PROVENANCE_TOPOLOGY',
    'Lagrangian/model/SU2Invariants.wl': 'T3_B2_PROVENANCE_SU2',
}
ENV_KEYS = ('T3_B2_PROVENANCE_HEAD', *PROVENANCE_ENV.values())
RESTRICTED_FLAGS = (
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES',
)


def provenance_environment(provenance: dict, parent_env: dict[str, str]) -> dict[str, str]:
    """Use Python's exact checked source snapshot; never trust inherited values."""
    require(provenance.get('all_six_source_blobs_unchanged') is True,
            'all-six source provenance has not been verified')
    head = provenance.get('repository_head_measured')
    require(isinstance(head, str) and re.fullmatch(r'[0-9a-f]{40}', head) is not None,
            'Git HEAD missing or invalid')
    rows = provenance.get('validated_git_blobs')
    require(isinstance(rows, list) and len(rows) == len(SOURCE_BLOBS),
            'source guard did not cover exactly six tracked files')
    by_path = {}
    for row in rows:
        require(isinstance(row, dict), 'invalid tracked source row')
        path = row.get('path')
        require(path in SOURCE_BLOBS and path not in by_path
                and row.get('git_blob_sha') == SOURCE_BLOBS[path]
                and row.get('validated_with_git_clean_filter') is True,
                'untrusted or duplicate source provenance entry')
        by_path[path] = row
    require(set(by_path) == set(SOURCE_BLOBS), 'incomplete native source provenance')
    require(set(PROVENANCE_ENV) == set(ALLOWED_B2_PATHS),
            'launcher/native validator source list mismatch')
    env = dict(parent_env)
    for key in ENV_KEYS:
        env.pop(key, None)
    env['T3_B2_PROVENANCE_HEAD'] = head
    for path, name in PROVENANCE_ENV.items():
        env[name] = SOURCE_BLOBS[path]
    return env


def run_native(root: Path, output: Path, *, wolfram: str = 'wolframscript',
               env: dict[str, str] | None = None) -> dict:
    root = root.resolve()
    output = output.resolve()
    require(not output.exists(), f'refusing existing native output: {output}')
    script = root / WL_SCRIPT
    require(script.is_file(), f'missing native probe: {script}')
    parent = dict(os.environ if env is None else env)
    require(all(parent.get(flag, '0') != '1' for flag in RESTRICTED_FLAGS),
            'nonrestricted T3_EXPERIMENTAL_GENERAL_* environment')
    # Includes: exact git top-level, HEAD, tree, index, clean-filter worktree
    # hashes, clean tracked paths and unchanged expected Lagrangian contracts.
    proven = source_provenance(root)
    child_env = provenance_environment(proven, parent)
    print('B2 source provenance verified in Python; running native RGBeta.', flush=True)
    try:
        result = subprocess.run(
            [wolfram, '-file', str(script), str(output)],
            cwd=root, env=child_env, check=False,
        )
    except OSError as exc:
        raise AuditBlocked(f'cannot launch wolframscript: {exc}') from exc
    require(result.returncode == 0,
            f'native RGBeta probe exited {result.returncode}; '
            f'inspect diagnostic output if present: {output}')
    data = read_report(output)
    native_snapshot_guard(data, root, proven)
    valid = validate_native_data(data)
    require(valid['native_registry_ratios_measured'] == 30,
            'native coupling inventory did not verify all 30 ratios')
    return {'native_output': str(output), 'git_head': proven['repository_head_measured'],
            'families': valid['families_checked'], 'internal_ratios': 30,
            'physical_source_certified': False}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--wolfram-exe', default='wolframscript',
                        help='full executable path if Wolfram is not on PATH')
    args = parser.parse_args(argv)
    output = args.output if args.output.is_absolute() else ROOT / args.output
    try:
        report = run_native(ROOT, output, wolfram=args.wolfram_exe)
    except (AuditBlocked, ValueError, KeyError, IndexError, OSError) as exc:
        print('BLOCKED_B2_NATIVE_INTERNAL_PROBE:', exc, file=sys.stderr)
        return 3
    print(f"NATIVE_RGBETA_INTERNAL_RATIOS_MEASURED; families="
          f"{report['families']}/5; internal ratios={report['internal_ratios']}/30; "
          f"physical=0/16; output={report['native_output']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
