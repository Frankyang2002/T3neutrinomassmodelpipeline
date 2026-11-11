"""Read-only A(-4) fermion-threshold LLSS diagnostic.

Project path: scripts/DiagnoseT3FermionThresholdA4.py

Recreates the SAME fitted UV trajectory used in corrected_llss16_v4, only
from mu_UV to mu_F, and invokes its existing dynamic F-tree matching.
No fits, EFT1/scalar/SMEFT evolution, or external repository writes.

The returned LLSS tensor is a source-exact *formal* commuting-spinor tensor;
its absolute Matchete GammaCC -> master-Weyl normalization is uncertified.
No full one-loop hard-fermion C5 matching is calculated by this tool.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from Numerical.orchestration.GeneralRunningBenchmark import initial_state
from Numerical.running.UVRunner import run_uv_segment
from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary
from scripts import RunT3GeneralSequentialFitAblation as sequential
from scripts import RunT3GeneralAll16FromReference as general

ROOT = Path(__file__).resolve().parents[1]
MODEL = "A_-4"
F_SCALE = 1e5


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frob(x) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=np.complex128)))


def encode_complex(z):
    z = complex(z)
    return {"real": float(z.real), "imag": float(z.imag)}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-run", default="thesis16_merged_v1")
    p.add_argument("--reference-run", default="corrected_llss16_v4")
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--wolframscript", default="wolframscript")
    args = p.parse_args(argv)

    source_root = general.OUTPUT / args.source_run
    saved = source_root / sequential.tag(MODEL)
    overrides_path = saved / "general_uv_overrides.json"
    fit_settings = source_root / "settings.json"
    reference_path = (ROOT / "output/studies/general_scalar_extended/"
                      "corrected_llss_eft1_freeze" / args.reference_run /
                      sequential.tag(MODEL) / "cases.json")
    for path in (overrides_path, fit_settings, reference_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required existing campaign input not found: {path}")
    settings = json.loads(fit_settings.read_text(encoding="utf-8-sig"))
    canonical_run = settings.get("canonical_run_id", "full16_general_v1")
    reference = json.loads(reference_path.read_text(encoding="utf-8-sig"))
    if reference.get("model") != MODEL or reference.get("status") != "Computed":
        raise ValueError("A(-4) reference campaign case is not computed")
    reference_case = reference["cases"]["1111"]

    # Exactly the same nonmixed UV initial-state construction and UV ODE settings
    # as scripts/RunT3GeneralSequentialFitAblation.py (the corrected runner's base).
    beta_uv, _beta_eft, _unused, _kernel, _c5path, _kind = general.nonmixed_inputs(
        MODEL, canonical_run, generate_general_c5=False,
        wolframscript=args.wolframscript)
    overrides = json.loads(overrides_path.read_text(encoding="utf-8-sig"))
    _, initial = initial_state(beta_uv, mu_uv=sequential.SCALES[0],
                               mu_f=sequential.SCALES[1],
                               mu_s=sequential.SCALES[2], overrides=overrides)
    sol = run_uv_segment(initial, beta_uv, F_SCALE, rtol=1e-7, atol=1e-10)
    uv_at_f = sol.final_state
    source, diagnostics = nonmixed_boundary(MODEL, uv_at_f,
                                            repo_root=ROOT, mu_f=F_SCALE)
    tensor = np.asarray(source.raw, np.complex128)
    by_source = np.asarray(source.by_source, np.complex128)
    if tensor.shape != (3, 3, 2, 2, 8, 8):
        raise ValueError(f"Unexpected A(-4) LLSS tensor shape {tensor.shape}")
    residual = frob(tensor - by_source.sum(axis=0)) / max(frob(tensor), 1e-40)
    if residual > 1e-11:
        raise ValueError(f"F-tree source additivity failed: {residual}")

    # Check the recomputed heavy masses against the ACTUAL saved corrected
    # run, to protect against accidental use of an incompatible updated fit.
    expected = (reference.get("provenance", {}).get("fermion_tree", {})
                .get("1", {}).get("heavy_mass_eigenvalues_gev"))
    actual = np.asarray(diagnostics["heavy_mass_eigenvalues_gev"], float)
    if expected is None:
        raise ValueError("Saved corrected campaign has no mu_F heavy-mass diagnostic")
    mass_relative = float(np.linalg.norm(actual - np.asarray(expected, float)) /
                          max(np.linalg.norm(expected), 1e-40))
    if mass_relative > 1e-7:
        raise ValueError(f"Fitted mu_F boundary does not reproduce the saved "
                         f"heavy masses (relative discrepancy {mass_relative:.3g})")
    expected_indices = [int(r["source_index"]) for r in
                        reference["provenance"]["fermion_tree"]["1"]["sources"]]
    if list(source.source_indices) != expected_indices:
        raise ValueError("F-tree source indices differ from saved campaign")

    out = args.output or (ROOT / "Reports/output/thesis_qualitative/"
                           "eft1_llss_contributions" / args.reference_run /
                           "fermion_threshold_A_m4")
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "A_m4_LLSS_muF.npz",
                        C_LLSS_formal=tensor,
                        C_LLSS_formal_per_source=by_source,
                        source_indices=np.asarray(source.source_indices, dtype=int),
                        mu_F_gev=np.asarray(F_SCALE))
    sources = []
    for i, idx in enumerate(source.source_indices):
        norm = frob(by_source[i])
        projection = float(np.vdot(tensor, by_source[i]).real /
                           max(frob(tensor)**2, 1e-80))
        sources.append({"source_index": int(idx),
                        "formal_tensor_frobenius_gev_inverse": norm,
                        "signed_projection_fraction": projection,
                        "fraction_of_total_norm": norm/max(frob(tensor),1e-40)})
    with (out / "A_m4_LLSS_sources.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(sources[0]))
        writer.writeheader()
        writer.writerows(sources)
    # Sparse coordinate listing in the *exact* source tensor axes.
    with (out / "A_m4_LLSS_nonzero_components.csv").open(
            "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["flavour_alpha_0based", "flavour_beta_0based",
                         "su2_i_0based", "su2_j_0based", "scalar_a_0based",
                         "scalar_b_0based", "coefficient_real_gev_inverse",
                         "coefficient_imag_gev_inverse"])
        for indices in np.argwhere(np.abs(tensor) > 0):
            coord = tuple(int(k) for k in indices)
            z = complex(tensor[coord])
            writer.writerow([*coord, z.real, z.imag])
    top = sorted(((*[int(k) for k in indices], encode_complex(tensor[tuple(indices)]))
                  for indices in np.argwhere(np.abs(tensor)>0)),
                 key=lambda item: abs(complex(item[-1]["real"], item[-1]["imag"])),
                 reverse=True)[:12]
    report = {
        "model": MODEL, "matching_scale_gev": F_SCALE,
        "source_campaign": args.source_run, "reference_campaign": args.reference_run,
        "reference_case": "1111", "mu_F_UV_running_nfev": int(sol.nfev),
        "tensor_axis_order": ["lepton_flavour_alpha", "lepton_flavour_beta",
                              "lepton_su2_i", "lepton_su2_j",
                              "heavy_real_scalar_a", "heavy_real_scalar_b"],
        "tensor_shape": list(tensor.shape), "scalar_labels": list(source.labels),
        "tensor_frobenius_gev_inverse": frob(tensor),
        "source_sum_relative_residual": residual,
        "source_contributions": sources,
        "largest_components_0based": [{"indices": list(item[:6]),
                                       "coefficient_gev_inverse": item[-1]}
                                      for item in top],
        "saved_muF_heavy_mass_relative_discrepancy": mass_relative,
        "heavy_mass_eigenvalues_gev": actual.tolist(),
        "UV_override_sha256": sha256(overrides_path),
        "reference_cases_sha256": sha256(reference_path),
        "source_seed_sha256": sha256(source.path),
        "source_status": source.source_status,
        "normalization": {
            "value_is_formal_matchete_source": True,
            "spinor_factor_used_by_corrected_campaign": 1.0,
            "GammaCC_to_master_Weyl_absolute_normalization_certified": False,
        },
        "weinberg_muF": {
            "C5_from_tree_fermion_matching": "absent in retained LLSS tree sector",
            "hard_one_loop_C5_at_muF": "not calculated / unknown",
            "full_C5_at_muF": "not determined",
            "do_not_infer_full_C5_equals_zero": True,
        },
        "existing_saved_scalar_threshold": {
            "C5_LLSS_log_frobenius_gev_inverse": float(reference_case["c5_log_frobenius"]),
            "C5_LLSS_log_over_full_norm": float(reference_case["c5_log_over_independent_llss_norm"]),
        },
        "notes": [
            "Replays fitted UV 1e7->1e5 GeV RGE only; no benchmark fit.",
            "Formal LLSS tensor is not a certified physical Weyl Wilson coefficient.",
            "The LLSS index convention is the established source polynomial convention.",
            "Full hard fermion C5 matching has not been independently evaluated.",
        ],
    }
    (out / "A_m4_fermion_matching_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"A(-4) muF LLSS formal tensor: {tensor.shape}; "
          f"Frobenius={frob(tensor):.9e} GeV^-1")
    print(f"Source additivity residual: {residual:.3e}; "
          f"saved muF mass discrepancy: {mass_relative:.3e}")
    for row in sources:
        print(f"  F-tree source {row['source_index']}: "
              f"norm={row['formal_tensor_frobenius_gev_inverse']:.9e} GeV^-1; "
              f"signed projection={100*row['signed_projection_fraction']:.3f}%")
    print("C5 at muF: full one-loop hard matching NOT calculated")
    print("Output:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
