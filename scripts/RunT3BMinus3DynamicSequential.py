"""B(-3) dynamic renormalisable UV/EFT1 -> sequential-C5 -> SMEFT benchmark.

This is Stage 2 of the B(-3) sequential numerical pilot. It RE-RUNS the
canonical, full-general UV and scalar-only EFT1 RGBeta ODEs for every UV input
point. It then numerically substitutes the resulting muF heavy-fermion Yukawas
and masses, and the muS scalar parameters, into the SAVED Matchete F-first /
scalar-second fixed-order one-loop C5 expression, before evolving SMEFT.

What this does NOT do:
  - rerun Matchete symbolically for every UV parameter point (unnecessary for
    a parameter-independent matching kernel);
  - rederive all independent EFT1 dimension-five anomalous dimensions;
  - numerically transport the full LLSS tensor through EFT1 into final C5.
The saved symbolic kernel contains its direct LLSS -> Weinberg leading log.
Thus this is an RG-improved, fixed-order *conditional* sequential treatment,
not a general full-seven-term tensor or all-16 calculation.

Requires scripts/RunT3BMinus3SequentialSMEFT.py from Stage 1.

The canonical RGBeta cache rejects source fingerprints when symbolic source
files have changed. By default, this driver first attempts that strict loader,
and falls back ONLY on its explicit fingerprint-mismatch error to the saved
B(-3) UV/EFT1 beta JSON pair from the original fit campaign. That fallback
checks model, stage, beta inventory, and the absence of fermions in EFT1; it
DOES NOT claim the historic beta functions match the current symbolic source.
Use --beta-source canonical to forbid fallback, or --beta-source archived to
load the historic beta pair directly.

From the repository root (PowerShell):
  python -X utf8 -m scripts.RunT3BMinus3DynamicSequential --check-only
  python -X utf8 -m scripts.RunT3BMinus3DynamicSequential --run-id baseline_v1
  python -X utf8 -m scripts.RunT3BMinus3DynamicSequential --run-id y2_v1 --yukawa-scale 2
  python -X utf8 -m scripts.RunT3BMinus3DynamicSequential --run-id tree_v1 --export-fermion-tree

Additional partial overrides can be supplied as --param-overrides PATH. Values
use the SAME general UV coupling keys and JSON encoding as the saved
"general_uv_overrides.json". Existing saved results are never overwritten.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

from scripts import RunT3BMinus3SequentialSMEFT as stage1


ROOT = Path(__file__).resolve().parents[1]
MODEL = "B_-3"
DEFAULT_C5 = stage1.DEFAULT_C5
DEFAULT_FIT = stage1.DEFAULT_FIT
DEFAULT_OUTPUT = ROOT / "output/general_sequential_pilot/B_m3_dynamic_sequential_smeft"
SCHEMA = "t3-b-minus-3-dynamic-rge-conditional-sequential-smeft-v2"
ARCHIVED_BETA_ROOT = ROOT / "output/studies/general_scalar_extended/complete_scalar_symbolic"
HEAVY_UV_KEYS = frozenset(("y1", "y2", "MF"))
REQUIRED_COMMON_BETAS = frozenset((
    "gY", "g2", "g3", "yu", "yd", "ye", "lambdaH", "mS1Sq", "mS2Sq",
    "lambdaS1", "lambdaS2", "lambdaH1", "lambdaH2", "lambda12", "lambdaT3",
    "lambdaH1Adj", "lambdaH2Adj", "lambda12Adj",
))
REPRO_KEYS_UV = ("gY", "g2", "g3", "lambdaH", "yu", "yd", "ye",
                 "y1", "y2", "MF", "mS1Sq", "mS2Sq", "lambdaT3")
REPRO_KEYS_EFT1 = ("gY", "g2", "g3", "lambdaH", "yu", "yd", "ye",
                   "mS1Sq", "mS2Sq", "lambdaT3")


def local(path: Path) -> Path:
    """Interpret relative user paths from the project root, not cwd."""
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def _sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_archived_beta_pair(run_id: str) -> tuple[dict, dict, dict]:
    """Use the *saved* B(-3) RGBeta pair without asserting current-source identity.

    Only called explicitly (--beta-source archived) or after the strict loader
    reports a source fingerprint mismatch. No global validation is disabled.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_id):
        raise ValueError("Invalid canonical run-id for archived beta lookup")
    folder = ARCHIVED_BETA_ROOT / run_id / "B_m3"
    paths = {
        "UV": folder / "uv_rgbeta.json",
        "EFT1": folder / "eft1_rgbeta.json",
    }
    payloads = {stage: stage1.read_object(path)
                for stage, path in paths.items()}
    required_rep = {"dS1": 2, "dS2": 2, "dF": 1, "alpha": -3}
    keys = {}
    for stage, doc in payloads.items():
        metadata = doc.get("metadata", {})
        if doc.get("status") != "SuccessCandidateMassBeta":
            raise ValueError(f"Archived {stage} RGBeta file has unexpected status")
        if metadata.get("GeneralScalarModel") != MODEL or (
            metadata.get("GeneralScalarStage") != stage
        ) or any(metadata.get(k) != v for k, v in required_rep.items()):
            raise ValueError(f"Archived {stage} RGBeta model/stage/representation mismatch")
        if metadata.get("SharedScalar") is not False:
            raise ValueError(f"Archived {stage} is not the split-scalar model")
        if metadata.get("ReportBetaConvention") != (
            "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
        ):
            raise ValueError(f"Archived {stage} beta convention differs from numerical runner")
        beta = doc.get("report_betas")
        if not isinstance(beta, dict) or not beta or any(
            not isinstance(name, str) or not isinstance(expr, str) or not expr.strip()
            for name, expr in beta.items()
        ):
            raise ValueError(f"Archived {stage} beta expressions are absent/malformed")
        keys[stage] = set(beta)
        expected_extra = metadata.get("GeneralScalarNewBetaKeys", [])
        if not isinstance(expected_extra, list) or not set(expected_extra) <= keys[stage]:
            raise ValueError(f"Archived {stage} is missing a general scalar beta")
    if not REQUIRED_COMMON_BETAS <= keys["UV"] & keys["EFT1"]:
        raise ValueError("Archived UV/EFT1 beta pair lacks common SM/scalar couplings")
    if not HEAVY_UV_KEYS <= keys["UV"] or keys["EFT1"] & HEAVY_UV_KEYS:
        raise ValueError("Archived fermion-first EFT1 inventory is incorrect")
    if keys["UV"] - HEAVY_UV_KEYS != keys["EFT1"]:
        raise ValueError("Archived B(-3) UV/EFT1 coupling inventories disagree")
    return payloads["UV"], payloads["EFT1"], {
        "source": "archived_b_minus_3_rgbeta_pair",
        "model": MODEL,
        "run_id": run_id,
        "uv_json": str(paths["UV"].resolve()),
        "eft1_json": str(paths["EFT1"].resolve()),
        "uv_sha256": _sha256(paths["UV"]),
        "eft1_sha256": _sha256(paths["EFT1"]),
        "UV_beta_count": len(keys["UV"]),
        "EFT1_beta_count": len(keys["EFT1"]),
        "archived_beta_inventories_checked": True,
        "current_symbolic_source_fingerprint_validated": False,
        "historic_beta_functions_may_differ_from_current_symbolic_source": True,
        "physical_matchete_rgbeta_mapping_certified": False,
    }


def load_beta_pair(run_id: str, beta_source: str) -> tuple[dict, dict, dict]:
    """Strict canonical loading or transparent B(-3) historical-cache fallback."""
    if beta_source not in ("auto", "canonical", "archived"):
        raise ValueError("--beta-source must be auto, canonical, or archived")
    if beta_source == "archived":
        uv, eft, info = _load_archived_beta_pair(run_id)
        print("NOTICE: using archived RGBeta beta files; current-source fingerprints "
              "are not asserted.", flush=True)
        info["requested_beta_source"] = beta_source
        return uv, eft, info
    from Numerical.orchestration.GeneralRunningBenchmark import (
        load_canonical_general_exports,
    )
    try:
        uv, eft, info = load_canonical_general_exports(MODEL, run_id)
    except ValueError as exc:
        if beta_source == "canonical" or (
            "source fingerprint differs (or missing provenance)" not in str(exc)
        ):
            raise
        print("NOTICE: canonical RGBeta cache rejected a source fingerprint; "
              "falling back to the archived B(-3) beta pair used by the saved fit.",
              flush=True)
        uv, eft, info = _load_archived_beta_pair(run_id)
        info["canonical_rejection"] = str(exc)
    else:
        info = dict(info)
        info["current_symbolic_source_fingerprint_validated"] = True
    info["requested_beta_source"] = beta_source
    return uv, eft, info


def _matrix_to_entries(matrix: np.ndarray) -> list[list[dict[str, float]]]:
    values = np.asarray(matrix, dtype=complex)
    if values.shape != (3, 3) or not np.isfinite(values).all():
        raise ValueError("Scaled Yukawa must be a finite complex 3x3 matrix")
    return [[{"real": float(z.real), "imag": float(z.imag)} for z in row]
            for row in values]


def build_uv_overrides(base: Mapping[str, Any],
                       changes: Mapping[str, Any] | None,
                       yukawa_scale: float) -> dict[str, Any]:
    """Preserve original full-general UV point and apply optional trial changes."""
    scale = float(yukawa_scale)
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("--yukawa-scale must be positive and finite")
    override = dict(base)
    if changes:
        unknown = set(changes).difference(override)
        if unknown:
            raise ValueError(
                "Unknown UV override key(s), not in fitted B(-3) input: "
                + ", ".join(sorted(unknown))
            )
        override.update(changes)
    if scale != 1.0:
        for key in ("y1", "y2"):
            override[key] = _matrix_to_entries(
                np.asarray(stage1.decode(override[key]), dtype=complex) * scale
            )
    return override


def _relative_difference(current: Any, reference: Any) -> float:
    a = np.asarray(stage1.decode(current), dtype=complex)
    b = np.asarray(stage1.decode(reference), dtype=complex)
    if a.shape != b.shape:
        raise ValueError("Reproduction comparison encountered inconsistent shapes")
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30))


def reproduction_diagnostic(current: Mapping[str, Any],
                            saved: Mapping[str, Any],
                            names: tuple[str, ...]) -> dict[str, float]:
    return {key: _relative_difference(current[key], saved[key])
            for key in names if key in current and key in saved}


def run_dynamic_rges(uv_beta: dict, eft_beta: dict, config: dict,
                     overrides: dict) -> dict:
    """Use EXACT existing RGBeta full-general UV/EFT1 ODE implementations."""
    from Numerical.orchestration.GeneralRunningBenchmark import run

    scales = config["scales"]
    mu_uv = float(config["base_state"]["mu_gev"])
    mu_f = float(scales["mu_fermion_threshold_gev"])
    mu_s = float(scales["mu_scalar_threshold_gev"])
    if not (mu_uv > mu_f > mu_s > 0):
        raise ValueError("Must have UV > all-fermion threshold > all-scalar threshold")
    result = run(
        uv_beta, eft_beta,
        mu_uv=mu_uv, mu_f=mu_f, mu_s=mu_s,
        overrides=overrides, rtol=1e-7, atol=1e-10,
    )
    if result.get("model") != MODEL or result.get("status") != (
        "NumericalRenormalisableRunningOnly"
    ):
        raise RuntimeError("Canonical UV/EFT1 RGE returned unexpected model/status")
    for label in ("uv_solver", "eft1_solver"):
        if result.get(label, {}).get("success") is not True:
            raise RuntimeError(f"Fresh RGBeta {label} was not successful")
        if int(result[label].get("nfev", 0)) <= 0:
            raise RuntimeError(f"Fresh RGBeta {label} performed no RHS evaluations")
    return result


def optional_tree_wilson(uv_endpoint: dict, config: dict, mu_f: float):
    """Separately reconstruct Matchete tree-level LLSS source from THIS muF state.

    The saved hierarchical analytic C5 kernel already embodies the B(-3)
    F-first matching. This tensor is an independent *output* for the next
    workstream; it is not falsely claimed to drive that kernel in Stage 2.
    """
    from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary
    heavy_state = stage1.fitted_fermion_state(uv_endpoint, config, mu_f)
    source, details = nonmixed_boundary(
        MODEL, heavy_state, repo_root=ROOT, mu_f=mu_f
    )
    tensor = np.asarray(source.raw, dtype=complex)
    if tensor.shape[:4] != (3, 3, 2, 2) or not np.isfinite(tensor).all():
        raise ValueError("Fermion-first Matchete LLSS tensor is invalid")
    summary = {
        "calculated_for_this_uv_point": True,
        "source_is_tree_matchete_F_threshold": True,
        "uses_actual_uv_at_muF": True,
        "shape": list(tensor.shape),
        "frobenius_norm_gev_inverse": float(np.linalg.norm(tensor)),
        "source_hash": source.hash,
        "source_diagnostics": details,
        "fed_to_final_C5_by_this_stage": False,
        "reason_not_fed_to_C5": (
            "The one-loop final C5 uses the previously derived symbolic "
            "F-first hard + LLSS-direct-running + joint-scalar matching kernel. "
            "Direct full-Wilson tensor transport remains Stage 3."
        ),
    }
    return tensor, summary


def make_result(c5_path: Path, fit: dict, config: dict,
                target_path: Path, run_data: dict,
                *, baseline: bool) -> dict:
    """Apply the Stage-1 symbolic C5 adapter to freshly evolved endpoints."""
    computed_fit = {
        "model": MODEL,
        "general_rges_used": True,
        "numerical_details": {
            "uv_at_fermion_threshold": run_data["uv_at_fermion_threshold"],
            "eft1_at_scalar_threshold": run_data["eft1_at_scalar_threshold"],
            "uv_nfev": int(run_data["uv_solver"]["nfev"]),
            "eft1_nfev": int(run_data["eft1_solver"]["nfev"]),
        },
    }
    result = stage1.run_numerics(c5_path, computed_fit, config, target_path)
    # The Stage-1 helper's original label referred to SAVED snapshots; here
    # it was given snapshots generated by a new integration above.
    result.pop("saved_renormalisable_running", None)
    result["status"] = "DynamicRGERunSequentialSymbolicC5AndSMEFTComputed"
    result["renormalisable_running"] = {
        "recomputed_from_UV_for_this_parameter_point": True,
        "UV": run_data["uv_solver"],
        "EFT1": run_data["eft1_solver"],
        "fermion_threshold_projection_performed": True,
        "fermion_removed_from_EFT1_beta_inventory": True,
        "mu_initial_to_F_to_S_gev": run_data["scales_gev"],
    }
    result["sequential_matching"].update({
        "symbolic_Matchete_kernel_reused": True,
        "matching_coefficient_evaluated_at_this_point": True,
        "fermion_Yukawas_and_masses_at_fresh_muF": True,
        "scalar_masses_and_lambdaT3_at_fresh_muS": True,
        "F_tree_LLSS_tensor_evolved_numerically_to_muS": False,
        "full_eft1_Wilson_gamma_matrix_solved_in_this_stage": False,
    })
    result["numerical_uv_and_eft1_endpoints"] = {
        "UV_initial": run_data["uv_initial"],
        "UV_at_F": run_data["uv_at_fermion_threshold"],
        "EFT1_at_F": run_data["eft1_at_fermion_threshold"],
        "EFT1_at_S": run_data["eft1_at_scalar_threshold"],
    }
    if baseline:
        previous = fit["numerical_details"]
        uv_diff = reproduction_diagnostic(
            run_data["uv_at_fermion_threshold"],
            previous["uv_at_fermion_threshold"], REPRO_KEYS_UV,
        )
        eft_diff = reproduction_diagnostic(
            run_data["eft1_at_scalar_threshold"],
            previous["eft1_at_scalar_threshold"], REPRO_KEYS_EFT1,
        )
        maximum = max([0.0, *uv_diff.values(), *eft_diff.values()])
        result["baseline_reproduction"] = {
            "compared_with_saved_one_scale_fit_running_only": True,
            "UV_muF_relative_differences": uv_diff,
            "EFT1_muS_relative_differences": eft_diff,
            "maximum_relative_difference": maximum,
            "new_matching_is_not_the_saved_one_scale_fit_matching": True,
        }
    return result


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id", default="baseline_v1")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--canonical-run-id", default="full16_general_v1")
    p.add_argument("--beta-source", choices=("auto", "canonical", "archived"),
                   default="auto", help=(
                       "auto: try current-source fingerprints, then use the saved "
                       "B(-3) beta pair only on fingerprint mismatch; canonical: "
                       "strict; archived: original-fit beta pair"
                   ))
    p.add_argument("--matching-json", type=Path, default=DEFAULT_C5)
    p.add_argument("--fit-root", type=Path, default=DEFAULT_FIT)
    p.add_argument("--param-overrides", type=Path, default=None)
    p.add_argument("--yukawa-scale", type=float, default=1.0)
    p.add_argument("--target", type=Path, default=None)
    p.add_argument("--export-fermion-tree", action="store_true")
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id):
        raise ValueError("--run-id must be one simple, safe folder name")

    c5_path = local(args.matching_json)
    fit_root = local(args.fit_root)
    fit_result_path = fit_root / "result.json"
    fit_config_path = fit_root / "fitted_config.json"
    base_overrides_path = fit_root / "general_uv_overrides.json"
    c5, fit, config, pilot = stage1.source_preflight(
        c5_path, fit_result_path, fit_config_path
    )
    base = stage1.read_object(base_overrides_path)
    change_path = local(args.param_overrides) if args.param_overrides else None
    changes = stage1.read_object(change_path) if change_path else None
    overrides = build_uv_overrides(base, changes, args.yukawa_scale)
    uv_beta, eft_beta, provenance = load_beta_pair(
        args.canonical_run_id, args.beta_source
    )
    if set(overrides).difference(uv_beta.get("report_betas", {})):
        unknown = sorted(set(overrides).difference(uv_beta.get("report_betas", {})))
        raise ValueError(f"Saved initial overrides absent from current RGBeta UV export: {unknown}")
    if args.check_only:
        print("B(-3): saved Matchete F-first and scalar-second symbolic kernel: OK")
        print("B(-3): UV and EFT1 RGBeta beta inventories: OK")
        print("B(-3): beta source = " + str(provenance.get("source")))
        if not provenance.get("current_symbolic_source_fingerprint_validated", False):
            print("NOTICE: archived beta payloads are not fingerprint-matched "
                  "to the current symbolic source")
        print("B(-3): current-UV parameter override keys: OK")
        print(f"Scalar-threshold pilot: {pilot['resumed_scalar_threshold']['status']}")
        print("No fresh RGE integration or output performed by --check-only")
        return 0

    if args.target is None:
        from Numerical.orchestration.RunSmallLambdaRefit import TARGET
        target = Path(TARGET)
    else:
        target = local(args.target)
    if not target.is_file():
        raise FileNotFoundError("Missing oscillation target: " + str(target))

    output = local(args.output_root) / args.run_id
    if output.exists():
        raise FileExistsError(f"Output already exists; choose another --run-id: {output}")

    running = run_dynamic_rges(uv_beta, eft_beta, config, overrides)
    baseline = changes is None and float(args.yukawa_scale) == 1.0
    result = make_result(c5_path, fit, config, target, running, baseline=baseline)
    tensor = None
    if args.export_fermion_tree:
        tensor, tree_summary = optional_tree_wilson(
            running["uv_at_fermion_threshold"], config,
            float(config["scales"]["mu_fermion_threshold_gev"]),
        )
        result["optional_fermion_tree_wilson"] = tree_summary
    else:
        result["optional_fermion_tree_wilson"] = {
            "calculated_for_this_uv_point": False,
            "fed_to_final_C5_by_this_stage": False,
        }
    result["format"] = SCHEMA
    result["model"] = MODEL
    result["symbolic_sequential_pilot"] = {
        "status": pilot["status"],
        "scalar_threshold_status": pilot["resumed_scalar_threshold"]["status"],
        "direct_weinberg_running_loaded": True,
        "hard_pole_direct_log_relation_checked": pilot.get(
            "pole_rge_consistency", {}
        ).get("direct_log_equals_twice_pole_residue"),
    }
    result["sources"] = {
        "canonical_run_id": args.canonical_run_id,
        "beta_source_mode": args.beta_source,
        "beta_source_provenance": provenance,
        "uv_beta_json": provenance.get("uv_json"),
        "eft1_beta_json": provenance.get("eft1_json"),
        "symbolic_sequential_matching_json": str(c5_path.resolve()),
        "original_fitted_overrides": str(base_overrides_path.resolve()),
        "optional_parameter_overrides": str(change_path.resolve()) if change_path else None,
        "oscillation_target": str(target.resolve()),
    }
    paths = (c5_path, c5_path.with_name("pilot_summary.json"),
             fit_result_path, fit_config_path, base_overrides_path, target)
    result["source_sha256"] = {
        str(p.resolve()): _sha256(p) for p in paths
    }
    if change_path:
        result["source_sha256"][str(change_path.resolve())] = _sha256(change_path)
    result["input_controls"] = {
        "historic_RGE_beta_pair_used": provenance.get("source") == "archived_b_minus_3_rgbeta_pair",
        "current_symbolic_source_fingerprint_validated": provenance.get(
            "current_symbolic_source_fingerprint_validated", False),
        "yukawa_scale": float(args.yukawa_scale),
        "additional_overrides_supplied": bool(changes),
        "new_UV_point_used_in_both_RGE_segments": True,
    }
    result["scope"] = (
        "Fresh full-general B(-3) UV/EFT1 RGBeta integrations; F-threshold "
        "parameter evaluation of saved Matchete sequential MSbar matching; "
        "EFT1 LLSS-to-Weinberg fixed one-loop leading log from saved source; "
        "joint scalar matched hard coefficient at S; SMEFT running and "
        "neutrino observables. No new Matchete symbolic run or full LLSS "
        "tensor transport in the final C5."
    )
    stage1.write_outputs(output, result)
    if tensor is not None:
        np.savez_compressed(
            output / "fermion_tree_LLSS_at_MF.npz",
            C_LLSS_gev_inverse=tensor,
            mu_F_gev=np.asarray(float(config["scales"]["mu_fermion_threshold_gev"])),
            model=np.asarray(MODEL),
        )

    print(f"Wrote: {output / 'result.json'}")
    print("Fresh UV/EFT1 RGE evaluations: "
          f"{running['uv_solver']['nfev']} / {running['eft1_solver']['nfev']}")
    if baseline:
        err = result["baseline_reproduction"]["maximum_relative_difference"]
        print(f"Fitted-trajectory endpoint reproduction max relative difference: {err:.3g}")
    print(f"chi2 = {result['neutrinos']['oscillation_chi2']:.8g}")
    print("masses (eV) = " + ", ".join(
        f"{m:.9g}" for m in result["neutrinos"]["masses_ev_m1_m2_m3"]
    ))
    if tensor is not None:
        print(f"Wrote: {output / 'fermion_tree_LLSS_at_MF.npz'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
