"""Opt-in exploratory *executed* B_-1/C_-1 full non-HHSS quartic runner.

Uses the existing mixed production pipeline for UV, F threshold, RGBeta
extended state + bilinear, seven-term Wilson, common scalar threshold,
and SMEFT. Only general_quartic is replaced while integrate_mixed runs.

The Matchete I_k -> RGBeta block-order identification is PROVISIONAL.
Run with a FRESH --run-id. Existing solver files and GitHub are untouched.
"""
from __future__ import annotations

import argparse
import contextlib
import json
from pathlib import Path
import sys

from Numerical.running.T3ExploratoryFullQuartic import build_overlay
from scripts import RunT3SimultaneousEFT1All16 as general
from scripts import RunT3SimultaneousEFT1MixedBC as mixed

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / 'output/general_t3_stage22/matchete_cartesian_candidates.json'
REGISTRY = ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
LABEL = 'EXPLORATORY_Matchete_I_k_equals_RGBeta_registry_block_order'


def make_quartic_override(model, original_quartic, candidate, registry):
    """Return a drop-in general_quartic(model,state,HHSS) function.

    The HHSS tensor comes from the original fixed-Cartesian mixed solver at
    the *same* RHS scale. Model-specific extra quartics come from its native
    extended RGBeta state, not frozen boundary values.
    """
    active = {'hhss': None}
    reg = json.loads(Path(registry).read_text(encoding='utf-8-sig'))
    entry = reg['models'][model]
    required = {
        row['name'] for row in entry['new_rgbeta_symbols']
    } | {
        row['existing_name'] for row in entry['preexisting_rgbeta_symbols']
    }
    overlay = build_overlay(
        model, candidate, registry, original_quartic,
        lambda _state: active['hhss'], required_quartics=required,
    )
    if len(overlay.directions) != 6:
        raise RuntimeError(f'{model}: expected three additional non-HHSS complex quartics')

    def quartic(current_model, running_state, hhss):
        if current_model != model:
            raise RuntimeError(f'Wrong model in exploratory mixed RHS: {current_model}')
        active['hhss'] = hhss
        result = overlay.evaluate(running_state)
        overlay.validate_hhss(running_state, hhss)
        return result

    return quartic, overlay


@contextlib.contextmanager
def full_mixed_quartic(model, candidate=CANDIDATE, registry=REGISTRY):
    """Temporarily inject additional quartics into existing mixed RHS."""
    old = general.general_quartic
    replacement, overlay = make_quartic_override(model, old, candidate, registry)
    general.general_quartic = replacement
    try:
        yield overlay
    finally:
        general.general_quartic = old


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', required=True)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--models', default='B_-1,C_-1')
    parser.add_argument('--candidate', type=Path, default=CANDIDATE)
    parser.add_argument('--registry', type=Path, default=REGISTRY)
    parser.add_argument('--check-only', action='store_true')
    args = parser.parse_args(argv)

    models = tuple(s.strip() for s in args.models.split(',') if s.strip())
    if not models or len(models) != len(set(models)) or not set(models) <= {'B_-1','C_-1'}:
        parser.error('Only B_-1,C_-1 are supported')
    if not args.candidate.is_file() or not args.registry.is_file():
        parser.error('Missing native Matchete candidate export or registry')

    # Fail early on a reused run ID, before the expensive model pipeline.
    output = (ROOT / 'output/studies/general_scalar_extended'
              / 'experimental_simultaneous_eft1_mixed_BC' / args.run_id)
    if output.exists() and not args.check_only:
        parser.error(
            f'Run output already exists: {output}. '
            'Supply a NEW --run-id; existing results will not be overwritten.')

    original = mixed._mixed_simultaneous

    def run_one(model, *values, **kwargs):
        # The existing run invokes this once per model. It retains its own
        # UV -> F matching -> extended EFT1 -> joint scalar -> SMEFT sequence.
        with full_mixed_quartic(model, args.candidate, args.registry) as overlay:
            cases, details = original(model, *values, **kwargs)
        for case in cases.values():
            case['exploratory_full_non_HHSS_quartics_in_Wilson_RHS'] = True
            case['quartic_mapping'] = LABEL
            case['non_HHSS_scalar_ansatz_provisional'] = True
            case['quartic_source'] = overlay.source
            case['physical_rgbeta_mapping_certified'] = False
        details['exploratory_full_non_HHSS_quartics_in_Wilson_RHS'] = True
        details['quartic_mapping'] = LABEL
        details['physical_rgbeta_mapping_certified'] = False
        details['quartic_source'] = overlay.source
        return cases, details

    mixed._mixed_simultaneous = run_one
    try:
        forwarded = ['--source-run', args.source_run, '--run-id', args.run_id,
                     '--models', ','.join(models)]
        if args.check_only:
            forwarded.append('--check-only')
        return mixed.main(forwarded)
    finally:
        mixed._mixed_simultaneous = original


if __name__ == '__main__':
    raise SystemExit(main())
