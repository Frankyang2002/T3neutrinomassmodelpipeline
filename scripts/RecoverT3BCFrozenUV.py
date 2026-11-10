#!/usr/bin/env python3
"""Targeted B(-1)/C(-1) zero-source recovery for the exploratory UV-off runs.

This is a *scoped* compatibility correction: an individual Matchete LLSS
source term with a provably ZERO numerical scalar prefactor is permitted to
produce an empty component list (the zero tensor). The complete F-boundary
Wilson tensor MUST still be nonzero and finite.

Existing physics files, full-running baselines, successful results, fitted
Yukawas, and non-UV-off model calculations are not modified.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import json
import shutil
from pathlib import Path
import sys
import traceback
from unittest.mock import patch

import numpy as np

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
from Numerical.running import GeneralT3OrderedTensor as ordered_module
from Numerical.running import GeneralT3MixedExactFermionWilson as mixed_module
from RGE.running.intermediate.ScalarOnlyTensorAdapters import _term_scalar_expression
from RGE.running.intermediate.IntermediateMatcheteParsing import parse_mathematica_scalar
from scripts.RunT3All16ExploratoryFit import (
    corrected_quartic_unit_tensors, consistent_project_c5,
    exploratory_preflight, exploratory_rough_case
)
from scripts.BuildT3GeneralisedThesisTables import load_reference
from scripts.RunT3GeneralisedRunningComparisons import _frozen_uv, KEYS, _relative

REFERENCES = ("D_-2", "A_-4")
MODELS = ("B_-1", "C_-1")


def _validated_zero_operator(term: dict, registry: list, model: str,
                             original, zero_sources: list[dict], **kwargs):
    """Only recover an EXACTLY zero scalar-prefactor source contribution.

    A small but nonzero coefficient, cancellation within CG contractions, or
    a malformed source remains blocked until independently diagnosed.
    """
    try:
        return original(term, registry, model, **kwargs)
    except ValueError as exc:
        if str(exc) != "No nonzero real-component Wilson terms":
            raise
        if model not in MODELS:
            raise
        template = term["operator_template_unit_flavor_inputform"]
        raw = _term_scalar_expression({
            "Index": term.get("source_index", 0), "TermInputForm": template
        })
        prefactor = ordered_module._number(parse_mathematica_scalar(raw))
        if prefactor != 0j:
            raise ValueError(
                f"{model} source {term['source_index']}: no real-component "
                "entries but scalar prefactor is NONZERO "
                f"({prefactor!r}); refusing to discard a possible "
                "small or cancelling contribution"
            ) from exc
        basis = ordered_module.scalar_basis(model)
        zero_sources.append({
            "source_index": int(term["source_index"]),
            "reason": "exactly_zero_scalar_prefactor",
            "prefactor_real": float(prefactor.real),
            "prefactor_imag": float(prefactor.imag),
        })
        return {
            "source_index": term["source_index"],
            "basis_scalar_fields": dict(basis.dimensions),
            "real_scalar_dimension_including_higgs": basis.total_dimension,
            "flavour_matrix_axes_external_indices":
                term.get("flavour_matrix_axes_external_indices"),
            "ordered_not_canonical_symmetrized": True,
            "cg_registry_used": [],
            "component_count": 0,
            "components": [],
        }


def _assert_nonzero_full_wilson(model, source, pl, cg, masses, yukawas, U,
                                original, audit):
    result = original(model, source, pl, cg, masses, yukawas, U)
    by_source, indices, details, polynomial_audit = result
    total = np.asarray(by_source).sum(axis=0)
    if not np.isfinite(total).all() or np.max(np.abs(total)) == 0:
        raise ValueError(
            f"{model}: complete F-boundary Wilson tensor is zero/nonfinite; "
            "physical case NOT recovered"
        )
    audit["F_wilson_max_abs"] = float(np.max(np.abs(total)))
    audit["F_wilson_nonzero_sources"] = int(
        sum(np.any(np.asarray(c) != 0) for c in by_source)
    )
    audit["F_wilson_source_count"] = int(len(indices))
    if audit["F_wilson_nonzero_sources"] == 0:
        raise ValueError(f"{model}: no nonzero F-boundary source contributions")
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tables", type=Path,
                   default=Path("output/parallel_t3/thesis_generalised_tables_v1"))
    p.add_argument("--reference-results", type=Path, nargs="+", required=True)
    p.add_argument("--output", type=Path, required=True,
                   help="NEW output directory below output/parallel_t3")
    args = p.parse_args()

    out = (ROOT / args.output).resolve()
    base = (ROOT / "output/parallel_t3").resolve()
    if out.exists() or base not in out.parents:
        p.error("--output must be a NEW directory below output/parallel_t3")
    tables = (ROOT / args.tables).resolve()
    baseline = read_json(tables / "results.json")
    provenance = read_json(tables / "provenance.json")
    reference_dirs = [(ROOT / x).resolve() for x in args.reference_results]
    references = {r: load_reference(r, reference_dirs) for r in REFERENCES}
    for ref in REFERENCES:
        if str(references[ref]["source"]) != (
            provenance["reference_benchmarks"][ref]["result_json"]
        ):
            p.error(f"{ref}: reference source differs from full-running baseline")
    cfg = read_json(ROOT / BENCHMARK)
    out.mkdir(parents=True)
    rows = []

    original_run_one = b_runner.run_one
    original_project = ordered_module.project_ordered_term
    original_mixed_project = mixed_module.project_source

    for ref in REFERENCES:
        for model in MODELS:
            zero_sources = []
            audit = {}
            original = baseline[ref].get(model, {})
            full = original.get("observables") or {}
            row = {
                "reference": ref, "model": model, "stage_frozen": "UV",
                "status": "BLOCKED", "zero_scalar_prefactor_source_terms": zero_sources,
                "full_baseline_source": str(tables / "results.json"),
            }
            try:
                if original.get("status") != "VALID_EXPLORATORY":
                    raise RuntimeError("Missing valid full-running baseline")

                evaluator = FreshAll16Evaluator(
                    ROOT, model, cfg,
                    batch1="output/batch1_general_contract/batch1_native_v1",
                    gamma_factor=float(provenance["gamma_factor_assumed"]),
                    mu_low=float(provenance["mu_low_gev"]),
                    vev=float(cfg["vev_GeV"]), reduced_potential=True,
                    mass_bilinears=provenance["mass_bilinears"],
                    wilson_solver=provenance["wilson_solver"],
                    wilson_steps=int(provenance["wilson_steps"]),
                    keep_trial_files=False,
                )
                def local_project(term, registry, selected_model, **kwargs):
                    return _validated_zero_operator(
                        term, registry, selected_model, original_project,
                        zero_sources, **kwargs
                    )
                def local_mixed_project(m, source, pl, cg, masses, ys, U):
                    return _assert_nonzero_full_wilson(
                        m, source, pl, cg, masses, ys, U,
                        original_mixed_project, audit
                    )
                def capture_run_one(opts, selected_model):
                    result = original_run_one(opts, selected_model)
                    manifest = (ROOT / "output/parallel_t3/background" /
                                opts.run_id / result["tag"] /
                                "fermion_boundary_manifest.json")
                    if manifest.is_file():
                        folder = out / ref / model
                        folder.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(manifest, folder / "fermion_boundary_manifest.json")
                    return result
                with ExitStack() as stack:
                    stack.enter_context(patch.object(
                        b_runner, "run_one", capture_run_one))
                    stack.enter_context(patch.object(
                        ordered_module, "project_ordered_term", local_project))
                    stack.enter_context(patch.object(
                        mixed_module, "project_source", local_mixed_project))
                    stack.enter_context(patch.object(
                        builder_module, "quartic_unit_tensors",
                        corrected_quartic_unit_tensors))
                    stack.enter_context(patch.object(
                        workflow, "project_c5", consistent_project_c5))
                    stack.enter_context(patch.object(
                        rough_stream2, "project_c5", consistent_project_c5))
                    stack.enter_context(patch.object(
                        preflight_module, "preflight_pair", exploratory_preflight))
                    stack.enter_context(patch.object(
                        numerics_module, "run_rough_case", exploratory_rough_case))
                    stack.enter_context(patch.object(
                        UVRunner, "run_uv_segment", _frozen_uv))
                    trial = evaluator(references[ref]["x"].copy())

                obs = trial.observables
                row.update(
                    status="VALID_EXPLORATORY",
                    observables=obs, chi2=trial.chi2, wilson_audit=audit,
                )
                for key in KEYS:
                    row[key + "_full"] = full.get(key)
                    row[key + "_frozen"] = obs.get(key)
                    row[key + "_pct"] = _relative(full.get(key), obs.get(key))
            except Exception as exc:
                row.update(error=f"{type(exc).__name__}: {exc}",
                           traceback=traceback.format_exc(),
                           wilson_audit=audit)
            case_dir = out / ref / model
            manifest_path = case_dir / "fermion_boundary_manifest.json"
            if manifest_path.is_file():
                fm = read_json(manifest_path)
                row["fermion_boundary_status"] = fm.get("status")
                if fm.get("blocker"):
                    row["underlying_blocker"] = fm["blocker"]
                    row["underlying_traceback_tail"] = fm.get("traceback_tail")
            save_json(case_dir / "UV.json", row)
            rows.append(row)
            print(ref, model, row["status"],
                  row.get("error", ""),
                  f"zero sources={len(zero_sources)}", flush=True)

    save_json(out / "results.json", {
        "format": "t3-bc-uv-off-zero-prefactor-recovery-v1",
        "definition": "100*(full-frozen)/frozen",
        "exploratory": True,
        "zero_operator_policy": "exactly zero scalar prefactor only",
        "rows": rows,
    })
    save_json(out / "summary.json", {
        "requested": len(rows),
        "completed": sum(r["status"] == "VALID_EXPLORATORY" for r in rows),
        "failed": sum(r["status"] != "VALID_EXPLORATORY" for r in rows),
        "cases": [
            {"reference": r["reference"], "model": r["model"],
             "status": r["status"], "error": r.get("error"),
             "zero_sources": len(r["zero_scalar_prefactor_source_terms"])}
            for r in rows
        ],
    })
    return 0 if all(r["status"] == "VALID_EXPLORATORY" for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
