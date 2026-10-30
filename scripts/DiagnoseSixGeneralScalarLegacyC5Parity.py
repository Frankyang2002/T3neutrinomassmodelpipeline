"""Read-only B(-3)/C(-3) legacy C5 parity diagnostic.

Reuse the existing general-scalar shared-A(-2) point and exact production
candidate evaluators. Capture the original reference/fast C5 matrices before
calling the original parity guard. NEVER bypass matching verification.

Run from repository root:
 python -X utf8 -m scripts.DiagnoseSixGeneralScalarLegacyC5Parity
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts import CompareT3GeneralScalarSharedA2 as shared
from scripts import RunT3GeneralScalar48Campaign as bulk
from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
from Numerical.running.SMWeinbergStage import evaluate_symbolic_c5

MODELS = ("B_-3", "C_-3")
CASE = "large_uv_frozen_0p5"
LIMIT = 3e-8


def matrix_payload(a: np.ndarray) -> dict:
    a = np.asarray(a, dtype=complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError("Expected a finite 3x3 complex C5 matrix")
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def analyse(reference: np.ndarray, fast: np.ndarray) -> dict:
    reference = np.asarray(reference, dtype=complex)
    fast = np.asarray(fast, dtype=complex)
    if reference.shape != (3, 3) or fast.shape != (3, 3):
        raise ValueError("Unexpected C5 matrix shape")
    scale = max(float(np.max(np.abs(reference))), 1e-30)
    diff = fast - reference
    ref_fro2 = float(np.vdot(reference, reference).real)
    multiplier = (np.vdot(reference, fast) / ref_fro2) if ref_fro2 > 0 else 0j
    resid = fast - multiplier * reference
    ratio2_error = float(np.max(np.abs(fast - 2 * reference)) / scale)
    conj_error = float(np.max(np.abs(fast - np.conjugate(reference))) / scale)
    return {
        "reference_c5": matrix_payload(reference),
        "fast_c5": matrix_payload(fast),
        "difference_fast_minus_reference": matrix_payload(diff),
        "relative_max_error": float(np.max(np.abs(diff)) / scale),
        "relative_frobenius_error": float(np.linalg.norm(diff) / max(float(np.linalg.norm(reference)), 1e-30)),
        "best_complex_multiplier_fast_over_reference": {"real": float(multiplier.real), "imag": float(multiplier.imag)},
        "relative_proportionality_residual": float(np.linalg.norm(resid) / max(float(np.linalg.norm(fast)), 1e-30)),
        "relative_error_if_fast_equals_twice_reference": ratio2_error,
        "relative_error_if_fast_equals_conjugate_reference": conj_error,
        "original_threshold": LIMIT,
        "parity_pass": bool(np.max(np.abs(diff)) / scale <= LIMIT),
        "factor_of_two_consistent_to_original_tolerance": bool(ratio2_error <= LIMIT),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", default="bulk48_v1")
    parser.add_argument("--output-dir", type=Path, default=(
        shared.ROOT / "output/studies/general_scalar_extended/six_model_shared_a2_legacy_c5_diagnostic"))
    args = parser.parse_args()
    source = shared.BASE / args.campaign
    shared.campaign_rows(source)
    donor = shared.donor(source, CASE, allow_poor=False)
    shared.require(donor["status"] == "AcceptedDonor", "accepted large-coupling A(-2) donor required")
    result = {
        "kind": "read-only-general-scalar-legacy-c5-parity-diagnostic-v1",
        "scenario": CASE,
        "source_campaign_sha256": shared.sha(source / "campaign.json"),
        "donor_fit_sha256": donor["source_fit_sha256"],
        "note": "Exact existing shared-A2 benchmark; no fitting, physics source changes, or parity bypass",
        "models": {},
    }
    original = bulk.LegacyKernel.parity_check
    for model in MODELS:
        record: dict = {}

        def instrument(self, uv, eft):
            if not isinstance(self, bulk.LegacyKernel):
                raise TypeError("Unexpected matching kernel")
            view = SimpleNamespace(uv_threshold_state=uv, eft1_threshold_state=eft,
                                   mu_scalar_threshold_gev=1e3)
            config = extract_final_c5_inputs(view).as_config()
            ref = evaluate_symbolic_c5(self.path, config)
            quick = self.evaluate(uv, eft)
            record.update(analyse(ref, quick))
            record["matching_source"] = self.path.relative_to(shared.ROOT).as_posix()
            record["matching_source_sha256"] = shared.sha(self.path)
            # Preserve the guard exactly. An expected failure remains a failure.
            return original(self, uv, eft)

        try:
            bulk.LegacyKernel.parity_check = instrument
            evaluated = shared.evaluate_bulk(CASE, model, donor["parameters"], source,
                                             extra_real=-0.01)
        finally:
            bulk.LegacyKernel.parity_check = original
        record["original_shared_evaluation_status"] = evaluated["status"]
        record["original_shared_evaluation_reason"] = evaluated.get("reason")
        record["mass_screen"] = evaluated.get("screen")
        if not record.get("reference_c5"):
            record["diagnosis_unavailable"] = "Parity check was not reached; inspect original reason"
        result["models"][model] = record
        print(f"{model}: screen={record['mass_screen']}, original_status={evaluated['status']}", flush=True)
        if "relative_max_error" in record:
            print(f"  parity error={record['relative_max_error']:.12g}", flush=True)
            print(f"  best fast/reference={record['best_complex_multiplier_fast_over_reference']}", flush=True)
            print(f"  error for fast=2*reference: {record['relative_error_if_fast_equals_twice_reference']:.12g}", flush=True)
            print("  reference C5:", np.array(record["reference_c5"]["real"]) +
                  1j*np.array(record["reference_c5"]["imag"]), flush=True)
            print("  fast C5:", np.array(record["fast_c5"]["real"]) +
                  1j*np.array(record["fast_c5"]["imag"]), flush=True)
        else:
            print("  ", evaluated.get("reason"), flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dest = args.output_dir / "diagnostic.json"
    dest.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("Saved:", dest)
    return 0 if all("reference_c5" in result["models"][m] for m in MODELS) else 1


if __name__ == "__main__":
    raise SystemExit(main())
