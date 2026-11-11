"""Source-resolved LLSS EFT1-freeze ablation, using the original 16-model inputs.

Intended path: scripts/RunT3CorrectedEFT1Freeze.py

This runner reuses RunT3GeneralSequentialFitAblation's exact fit sources,
renormalisable UV/EFT1 evolution, fermion-tree sources, full scalar pair
response, and final SM+Weinberg evolution. It REPLACES only the original
one-scale-anchored subtraction in `_final_cases`.

For each UV/EFT1 background endpoint, compute independent LLSS components:
    C5(W=0) = C5_scalar_finite_bubble
    C5(W=1) = C5_scalar_finite_bubble + C5_EFT1_leading_log
No difference of separately normalised one-scale amplitudes is used.

The result is an exploratory LLSS-sector fixed-order computation, NOT complete
one-loop matching. A full physical prediction requires independent hard
fermion matching, other EFT1 operator channels and verified normalisations.

PowerShell from project root:
  python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss16_v1 --check-only
  python -X utf8 -m scripts.RunT3CorrectedEFT1Freeze --source-run thesis16_merged_v1 --run-id corrected_llss16_v1

The original pipeline and its saved outputs remain untouched.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

from Numerical.running.GeneralT3SequentialLeadingLog import decode, encode, response, require
from scripts import RunT3GeneralSequentialFitAblation as original


CORRECTED_SCHEMA = "t3-independent-llss-U_E_W_S-ablation-v1"
OUTPUT_BASE = original.ROOT / "output/studies/general_scalar_extended/corrected_llss_eft1_freeze"


def independent_llss_cases(model, c5_one, source, scalar, target, *,
                           uv_on, eft_on, uv_nfev, eft_nfev,
                           spinor_factor, vev):
    """Re-evaluate both W states with one common independently computed LLSS response."""
    llss = response(
        source, scalar, mu_f=original.SCALES[1], mu_s=original.SCALES[2],
        spinor_factor=complex(spinor_factor), check_scale_cancellation=True,
    )
    finite = decode(llss["c5_scalar_finite_bubble"])
    log = decode(llss["c5_eft1_leading_log"])
    full = original.c5_checked(finite + log)
    frozen_wilson = original.c5_checked(finite)
    # IMPORTANT: the original one-scale Matchete coefficient is used ONLY as
    # a diagnostic. It must never set either W=0 or W=1 prediction.
    one_scale = original.c5_checked(c5_one)
    one_scale_discrepancy = float(
        np.linalg.norm(one_scale-full) / max(np.linalg.norm(one_scale), 1e-40)
    )
    if not math.isfinite(one_scale_discrepancy):
        raise ValueError("non-finite one-scale/independent-LLSS discrepancy")
    output = {}
    for w_on, c5 in ((False, frozen_wilson), (True, full)):
        for sm_on in (False, True):
            key = f"{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}"
            result = original.low_energy(
                original.scalar_snapshot_state(scalar), c5, on=sm_on,
                mu_s=original.SCALES[2], mu_low=original.SCALES[3],
                vev=vev, target=target,
            )
            result.update({
                "case": key, "uv_on": uv_on, "eft1_on": eft_on,
                "wilson_log_on": w_on, "smeft_on": sm_on,
                "uv_nfev": uv_nfev, "eft1_nfev": eft_nfev,
                "c5_induced_wilson_log": encode(log),
                "c5_joint_scalar_bubble": encode(finite),
                "c5_at_scalar_threshold": encode(c5),
                "c5_matching_anchored_reference": encode(one_scale),
                "c5_bookkeeping_residual_matching": encode(np.zeros((3, 3), complex)),
                "c5_log_frobenius": float(np.linalg.norm(log)),
                "c5_one_scale_frobenius": float(np.linalg.norm(one_scale)),
                # Retained for compatibility with original write_csv;
                # documented as log / INDEPENDENT LLSS full norm, NOT old anchor.
                "c5_log_over_one_scale_norm": float(
                    np.linalg.norm(log) / max(np.linalg.norm(full), 1e-40)
                ),
                "c5_log_over_independent_llss_norm": float(
                    np.linalg.norm(log) / max(np.linalg.norm(full), 1e-40)
                ),
                "c5_original_one_scale_relative_discrepancy": one_scale_discrepancy,
                "per_general_interaction": llss["per_channel"],
                "quartic_provenance": llss["quartic_provenance"],
                "matching_scale_cancellation": llss["matching_scale_cancellation_residual"],
                "weinberg_projection_residual": max(
                    llss["weinberg_beta_projection_residual"],
                    llss["weinberg_bubble_projection_residual"],
                ),
                "matching_approximation": "independent_F_tree_LLSS_plus_joint_scalar_bubble_only",
                "hard_fermion_matching_included": False,
                "other_EFT1_operator_classes_included": False,
                "full_sequential_one_loop_C5_computed": False,
            })
            output[key] = result
    return output


def main(argv=None):
    # Separate output path protects old results. The upstream controller reads
    # all original model source manifests and performs its own strict preflight.
    # Resume is prohibited, avoiding reuse of the old anchored scheme.
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("--resume", action="store_true")
    known, _ = parser.parse_known_args(argv)
    if known.resume:
        raise ValueError("--resume disabled for corrected run; use a fresh --run-id")
    if "--run-id" not in (sys.argv[1:] if argv is None else argv):
        raise ValueError("A NEW --run-id is required")

    original._final_cases = independent_llss_cases
    original.OUT = OUTPUT_BASE
    original_save = original.save

    def write_with_corrected_provenance(path: Path, data: dict):
        if isinstance(data, dict) and data.get("format") == (
                "t3-full-general-16_U_E_W_S_ablation-v1"):
            data = dict(data)
            data["format"] = CORRECTED_SCHEMA
            data["one_scale_C5_used_as_W_on_reference"] = False
            data["hard_matching_residual_independently_calculated"] = False
            data["matching_scope"] = "source_resolved_LLSS_scalar_pair_only"
            data["independent_W_zero"] = True
            data["not_complete_physical_matching"] = True
        elif isinstance(data, dict) and data.get("status") == "Computed" and "cases" in data:
            data = dict(data)
            data["matching_approximation"] = "independent_LLSS_plus_bubble"
            data["full_sequential_general_C5_certified"] = False
        original_save(path, data)

    original.save = write_with_corrected_provenance
    return original.main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
