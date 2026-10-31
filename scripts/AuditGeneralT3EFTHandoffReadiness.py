"""Stage 3: read-only audit of per-particle EFT handoff dependencies.

No RG/matching callback is registered by this script. It validates candidate
Stage 2 events, reports missing branches and constructs a numerical handoff
requirements inventory; NEVER claims a physical EFT step was performed.

  python -m scripts.AuditGeneralT3EFTHandoffReadiness --all
  python -m scripts.AuditGeneralT3EFTHandoffReadiness --model B_-1
  python -m scripts.AuditGeneralT3EFTHandoffReadiness --all --input PATH --output PATH
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from Numerical.experimental.GeneralT3EFTHandoffController import (
    candidate_inventory_from_stage2,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / 'output/studies/general_scalar_extended/continuous_full_scalar/'
          'continuous_v1/workstream_d_fermion_audit_v2/audit_summary.json')
DEFAULT_OUT = (ROOT / 'output/studies/general_scalar_extended/continuous_full_scalar/'
               'continuous_v1/workstream_d_handoff_readiness_v1/summary.json')
FORMAT = 'general-t3-eft-handoff-readiness-v1'


def build_summary(source: dict, *, models: tuple[str, ...] | None = None,
                  source_hash: str = '') -> dict:
    if not isinstance(source, dict) or source.get('format') != \
            'general-t3-running-heavy-mass-crossing-audit-v1' or \
            source.get('status') != 'CandidateAuditComplete' or \
            source.get('physical_matching_certified') is not False:
        raise ValueError('Stage 3: expected complete Stage 2 candidate audit')
    records = source.get('models')
    if not isinstance(records, list) or not records or any(not isinstance(x, dict) for x in records):
        raise ValueError('Stage 3: invalid Stage 2 model records')
    if len({x.get('model') for x in records}) != len(records):
        raise ValueError('Stage 3: duplicated Stage 2 model')
    selected = set(models) if models is not None else {x['model'] for x in records}
    available = {x['model'] for x in records}
    if not selected or not selected.issubset(available):
        raise ValueError('Stage 3: requested models not present in audit')
    rows = []
    for record in records:
        if record['model'] not in selected:
            continue
        inv = candidate_inventory_from_stage2(record)
        rows.append({
            'model': inv.model,
            'status': 'MissingMatchingAndIntermediateRGEs',
            'candidate_mass_inventory_complete': inv.complete,
            'candidate_crossing_count': len(inv.events),
            'unresolved_mass_branches': list(inv.unresolved),
            'first_candidate_event': ({'particle': inv.events[0].particle,
                                       'scale_gev': inv.events[0].scale_gev}
                                      if inv.events else None),
            'ordered_candidate_events': [
                {'particle': x.particle, 'kind': x.kind,
                 'scale_gev': x.scale_gev, 'status': x.source_status}
                for x in inv.events],
            'runtime_handoffs_executed': 0,
            'intermediate_particle_removed_rge_callbacks_registered': 0,
            'threshold_matching_callbacks_registered': 0,
            'physical_matching_certified': False,
            'note': ('These are mass roots on the undecoupled full-UV plus saved '
                     'grouped-F scalar-EFT1 backgrounds, not EFT-matched thresholds.'),
        })
    return {
        'format': FORMAT, 'status': 'CandidateHandoffDependenciesAuditedNotExecuted',
        'source_sha256': source_hash, 'source_format': source['format'],
        'physical_matching_certified': False,
        'models_checked': len(rows),
        'candidate_mass_roots_found': sum(x['candidate_crossing_count'] for x in rows),
        'models_with_all_candidate_branches': sum(x['candidate_mass_inventory_complete'] for x in rows),
        'unresolved_candidate_branches': sum(len(x['unresolved_mass_branches']) for x in rows),
        'matching_callbacks_registered': 0, 'evolution_callbacks_registered': 0,
        'models': rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument('--all', action='store_true')
    selection.add_argument('--model')
    parser.add_argument('--input', type=Path, default=SOURCE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error('Stage 2 input missing: ' + str(args.input))
    if args.output.exists():
        parser.error('Output already exists; use a new --output path: '+str(args.output))
    data = args.input.read_bytes()
    record = json.loads(data.decode('utf-8-sig'))
    summary = build_summary(record, models=((args.model,) if args.model else None),
                            source_hash=sha256(data).hexdigest())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive mode avoids overwriting an existing previous audit, including races.
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')
    for row in summary['models']:
        print(row['model'], 'events=', row['candidate_crossing_count'],
              'missing=', row['unresolved_mass_branches'])
    print('Stage 3 readiness report:', args.output)
    print('No EFT running/matching performed; all callbacks still to be supplied.')


if __name__ == '__main__':
    main()
