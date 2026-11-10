#!/usr/bin/env python3
"""Read-only-to-production diagnostic of B_-1/C_-1 UV-frozen F matching.

Runs four new, fixed-Yukawa UV-frozen trials and captures the Workstream B
fermion_boundary_manifest.json before FreshAll16Evaluator removes its scratch.
No production results, fitted Yukawas or physics implementation are modified.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import shutil
import sys
import traceback
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
    BENCHMARK, FreshAll16Evaluator, read_json, save_json
)
from Numerical.parallel_t3.background import runner as b_runner
from Numerical.parallel_t3.rough_inputs import builder as builder_module
from Numerical.parallel_t3.threshold import preflight as preflight_module
from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
from Numerical.parallel_t3.wilson import workflow, rough_stream2
from Numerical.running import UVRunner
from scripts.RunT3All16ExploratoryFit import (
    corrected_quartic_unit_tensors, consistent_project_c5,
    exploratory_preflight, exploratory_rough_case
)
from scripts.BuildT3GeneralisedThesisTables import load_reference
from scripts.RunT3GeneralisedRunningComparisons import _frozen_uv


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tables', type=Path,
                   default=Path('output/parallel_t3/thesis_generalised_tables_v1'))
    p.add_argument('--reference-results', nargs='+', type=Path, required=True,
                   help='Same fitted-reference folders used by the original runner')
    p.add_argument('--output', type=Path, required=True,
                   help='NEW output directory below output/parallel_t3')
    args = p.parse_args()
    out = (ROOT / args.output).resolve()
    base = (ROOT / 'output/parallel_t3').resolve()
    if out.exists() or base not in out.parents:
        p.error('--output must be a NEW directory below output/parallel_t3')
    tables = (ROOT / args.tables).resolve()
    provenance = read_json(tables / 'provenance.json')
    references = {
        ref: load_reference(ref, [(ROOT / x).resolve() for x in args.reference_results])
        for ref in ('D_-2', 'A_-4')
    }
    for ref, obj in references.items():
        if str(obj['source']) != provenance['reference_benchmarks'][ref]['result_json']:
            p.error(f'{ref}: Yukawa reference differs from baseline provenance')
    cfg = read_json(ROOT / BENCHMARK)
    out.mkdir(parents=True)
    original_run_one = b_runner.run_one
    summary = []
    for ref in ('D_-2', 'A_-4'):
        for model in ('B_-1', 'C_-1'):
            case = out / ref / model
            case.mkdir(parents=True)
            def capture_run_one(opts, selected_model):
                result = original_run_one(opts, selected_model)
                source = (ROOT / 'output/parallel_t3/background' /
                          opts.run_id / result['tag'] /
                          'fermion_boundary_manifest.json')
                if source.is_file():
                    shutil.copy2(source, case / 'fermion_boundary_manifest.json')
                return result
            evaluator = FreshAll16Evaluator(
                ROOT, model, cfg,
                batch1='output/batch1_general_contract/batch1_native_v1',
                gamma_factor=float(provenance['gamma_factor_assumed']),
                mu_low=float(provenance['mu_low_gev']),
                vev=float(cfg['vev_GeV']), reduced_potential=True,
                mass_bilinears=provenance['mass_bilinears'],
                wilson_solver=provenance['wilson_solver'],
                wilson_steps=int(provenance['wilson_steps']),
                keep_trial_files=False
            )
            row = {'reference': ref, 'model': model, 'stage_frozen': 'UV'}
            try:
                with ExitStack() as stack:
                    stack.enter_context(patch.object(b_runner, 'run_one', capture_run_one))
                    stack.enter_context(patch.object(
                        builder_module, 'quartic_unit_tensors', corrected_quartic_unit_tensors))
                    stack.enter_context(patch.object(workflow, 'project_c5', consistent_project_c5))
                    stack.enter_context(patch.object(rough_stream2, 'project_c5', consistent_project_c5))
                    stack.enter_context(patch.object(
                        preflight_module, 'preflight_pair', exploratory_preflight))
                    stack.enter_context(patch.object(
                        numerics_module, 'run_rough_case', exploratory_rough_case))
                    stack.enter_context(patch.object(UVRunner, 'run_uv_segment', _frozen_uv))
                    trial = evaluator(references[ref]['x'].copy())
                row.update(status='VALID_EXPLORATORY', observables=trial.observables)
            except Exception as exc:
                row.update(status='BLOCKED', error=f'{type(exc).__name__}: {exc}',
                           traceback=traceback.format_exc())
            manifest_path = case / 'fermion_boundary_manifest.json'
            if manifest_path.is_file():
                fm = json.loads(manifest_path.read_text(encoding='utf-8'))
                row['fermion_matching_status'] = fm.get('status')
                row['underlying_blocker'] = fm.get('blocker')
                row['underlying_traceback_tail'] = fm.get('traceback_tail')
            save_json(case / 'diagnostic.json', row)
            summary.append({k: v for k, v in row.items()
                            if k not in ('observables', 'traceback', 'underlying_traceback_tail')})
            print(ref, model, row['status'], row.get('underlying_blocker',
                                                row.get('error', '')), flush=True)
    save_json(out / 'summary.json', {'cases': summary})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
