"""Run 16 ordinary T3 benchmarks through the conditional fermion-first EFT chain.

This is an additive entrypoint: the existing 14-model runner is unchanged.  For
B_-1 and C_-1 the Matchete-native mixed-doublet F-tree LLSS source is evaluated
at the *current* UV muF endpoint, including all source Yukawa/CG terms.  The
current EFT1 muS scalar state supplies U(muS) to that exploratory F boundary.
The shared runner then evaluates LLSS->C5 fixed-order log, joint scalar finite
bubble, SMEFT running and oscillation observables.  The independently matched
hard one-loop C5 is not generally available.  Neither GammaCC normalization
nor the muF/muS mixed scalar basis identification is certified.

PowerShell, from the repository root:
    python -X utf8 -m scripts.RunT3All16DynamicSequential16 --all --check-only
    python -X utf8 -m scripts.RunT3All16DynamicSequential16 --model B_-1 --run-id mixed_bm1_v1
    python -X utf8 -m scripts.RunT3All16DynamicSequential16 --model C_-1 --run-id mixed_cm1_v1
    python -X utf8 -m scripts.RunT3All16DynamicSequential16 --all --run-id batch16_v1

Requires the unchanged RunT3All16DynamicSequential.py, mixed native Matchete
source/CG files and existing RGBeta and quartic registry data in the project.
No writes to GitHub.  Output folders are write-once, as in the original runner.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3All16DynamicSequential as base

MIXED = frozenset(("B_-1", "C_-1"))
FORMAT = "t3-16-dynamic-native-mixed-llss-scalar-matching-smeft-candidate-v1"
_ORIGINAL_FIT_INPUT = base.fit_input
_ORIGINAL_DYNAMIC_SOFT_SECTOR = base.dynamic_soft_sector


def _mixed_fit_input(model: str, campaign_root: Path) -> tuple[Path, dict, dict, dict]:
    """Restore the B/C historical fit boundary, which lacks overrides.json.

    Their thesis16 fits record the *fitted* y1/y2 in fitted_config, with additional
    scalar quartics and mu2 specified by saved fractions, not by an overrides file.
    Reconstruct the exact declared boundary prescription rather than using
    GeneralRunningBenchmark's nonzero defaults for missing extra couplings.
    """
    if model not in MIXED:
        return _ORIGINAL_FIT_INPUT(model, campaign_root)

    from Numerical.running.GeneralT3MixedExactFermionWilson import read_native
    # Fail at preflight if the source-exact mixed PL/CG inventory is unavailable.
    read_native(model, base.ROOT)

    folder = base.saved_fit_folder(model, campaign_root)
    config = base.stage1.read_object(folder / "fitted_config.json")
    fit = base.stage1.read_object(folder / "result.json")
    expected = {
        "d_s1": base.DIMENSIONS[model[0]][0],
        "d_s2": base.DIMENSIONS[model[0]][1],
        "d_f": base.DIMENSIONS[model[0]][2],
        "alpha": int(model.split("_")[1]),
        "shared_scalar": False,
    }
    if config.get("representation") != expected or fit.get("model") != model:
        raise ValueError(f"{model}: historical fitted representation/model mismatch")
    if fit.get("status") != "ExploratoryGeneralScalarOneScaleConditionalFit":
        raise ValueError(f"{model}: unsupported historical fit status")
    if not isinstance(fit.get("details", {}).get("trajectory"), dict):
        raise ValueError(f"{model}: missing fitted trajectory")

    scales = config.get("scales", {})
    mu_uv = float(config["base_state"]["mu_gev"])
    if not (mu_uv > float(scales["mu_fermion_threshold_gev"]) >
            float(scales["mu_scalar_threshold_gev"]) >
            float(scales["mu_low_gev"]) > 0):
        raise ValueError(f"{model}: invalid UV > F > S > low scales")

    state = config["base_state"]
    ordinary = state["ordinary"]
    sm = state["sm"]
    if not isinstance(ordinary, dict) or not isinstance(sm, dict):
        raise ValueError(f"{model}: missing fitted UV ordinary/SM state")
    overrides = {**sm, **ordinary}
    if len(overrides) != len(sm) + len(ordinary):
        raise ValueError(f"{model}: duplicate SM and ordinary coupling keys")

    # The canonical RGBeta initial-state parser takes 3x3 arrays of {real,imag}
    # entries, not a pair of full {real: matrix, imag: matrix} arrays.
    for name in ("yu", "yd", "ye", "y1", "y2", "MF"):
        if name not in overrides:
            raise ValueError(f"{model}: missing fitted matrix {name}")
        matrix = base.stage1.finite_matrix(overrides[name], name)
        overrides[name] = base.stage2._matrix_to_entries(matrix)

    from scripts import CoupledUVRunning as inventory
    registry = base.stage1.read_object(
        base.ROOT / "RGE/general/generated/general_rgbeta_quartic_registry.json"
    )
    _, new, reused = inventory.inventories(registry, model)
    if len(new) != 1 or len(reused) != 6 or len(set(new + reused)) != 7:
        raise ValueError(f"{model}: unexpected native extra-quartic inventory")
    q_frac = float(fit["extra_quartic_fraction"])
    mu_frac = float(fit["extra_mixed_scalar_bilinear_fraction"])
    if (not math.isfinite(q_frac) or not 0 <= q_frac <= 1 or
            not math.isfinite(mu_frac) or not 0 <= mu_frac <= 0.2):
        raise ValueError(f"{model}: invalid saved extra-quartic or mixing fractions")
    lam = complex(base.stage1.decode(ordinary["lambdaT3"]))
    m1 = base.stage1.finite_real(ordinary["mS1Sq"], "mS1Sq")
    m2 = base.stage1.finite_real(ordinary["mS2Sq"], "mS2Sq")
    if min(m1, m2) <= 0 or abs(lam) <= 0:
        raise ValueError(f"{model}: nonpositive scalar mass or zero topology quartic")
    q = float(q_frac * abs(lam))
    mu2 = float(mu_frac * min(m1, m2))
    for name in (*new, *reused):
        overrides[name] = {"real": q, "imag": 0.0}
    if "muT3GeneralS1S2" in overrides:
        raise ValueError(f"{model}: duplicate mixed quadratic boundary")
    overrides["muT3GeneralS1S2"] = {"real": mu2, "imag": 0.0}
    return folder, config, fit, overrides


def _check_model(model: str) -> None:
    if model not in base.ALL_MODELS:
        raise ValueError(f"Unknown model {model!r}; choices: {base.ALL_MODELS}")
    # Unlike the original 14-model entrypoint, the native mixed boundary is
    # evaluated in _dynamic_soft_sector below.


def _dynamic_soft_sector(prepared: Any, running: dict, *,
                         factor: float, check_scale: bool) -> tuple[Any, Any, dict, dict]:
    if prepared.model not in MIXED:
        return _ORIGINAL_DYNAMIC_SOFT_SECTOR(
            prepared, running, factor=factor, check_scale=check_scale
        )

    from Numerical.running.GeneralT3MixedExactFermionWilson import mixed_exact_boundary
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors, response
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot

    model = prepared.model
    config = prepared.config
    mu_f = float(config["scales"]["mu_fermion_threshold_gev"])
    mu_s = float(config["scales"]["mu_scalar_threshold_gev"])
    at_f = running["uv_at_fermion_threshold"]
    at_s = running["eft1_at_scalar_threshold"]
    heavy = base.stage1.fitted_fermion_state(at_f, config, mu_f)
    source, info = mixed_exact_boundary(
        model, heavy, at_s, repo_root=base.ROOT, mu_f=mu_f
    )
    if tuple(source.labels) != tuple(prepared.reference.labels):
        raise ValueError(f"{model}: native mixed F-tree/reference scalar labels disagree")
    if source.raw.shape != (3, 3, 2, 2, 8, 8):
        raise ValueError(f"{model}: native mixed LLSS tensor has wrong shape")
    if not np.isfinite(source.raw).all():
        raise ValueError(f"{model}: nonfinite native mixed LLSS tensor")

    snapshot = RunningSnapshot(
        model=model,
        source=Path("fresh_full_general_EFT1_endpoint"),
        sha256=base.sha256(json.dumps(at_s, sort_keys=True).encode()).hexdigest(),
        mu_f=mu_f,
        mu_s=mu_s,
        scalar_mass1_sq=base.stage1.finite_real(at_s["mS1Sq"], "mS1Sq"),
        scalar_mass2_sq=base.stage1.finite_real(at_s["mS2Sq"], "mS2Sq"),
        couplings=at_s,
        evidence={"fresh_full_general_EFT1": True,
                  "mixed_native_CG_at_current_fitted_Yukawas": True,
                  "muF_to_muS_mixed_basis_identification_assumed": True},
    )
    bundle = scalar_tensors(
        model, source, prepared.quartics, snapshot, prepared.registry,
        include_general_extra=True, include_portal_completion=True
    )
    result = response(
        source, bundle, mu_f=mu_f, mu_s=mu_s,
        spinor_factor=complex(factor), check_scale_cancellation=check_scale
    )
    return source, bundle, result, info


def _preflight(model: str, args: Any) -> Any:
    prepared = base.preflight(model, args)
    if model in MIXED:
        required = set(prepared.uv_beta["report_betas"])
        provided = set(prepared.overrides)
        if required != provided:
            raise ValueError(
                f"{model}: mixed fitted UV boundary is incomplete/has extra keys: "
                f"missing={sorted(required - provided)}, "
                f"extra={sorted(provided - required)}"
            )
    return prepared


def main(argv: list[str] | None = None) -> int:
    # Use the existing parser, preflight, matching, output and no-overwrite
    # behaviour; only substitute the two source-specific integration hooks.
    base.fit_input = _mixed_fit_input
    base.check_model = _check_model
    base.dynamic_soft_sector = _dynamic_soft_sector
    base.FORMAT = FORMAT
    args = base.parser().parse_args(argv)
    import re
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id) is None:
        raise ValueError("Unsafe --run-id")
    if not math.isfinite(args.yukawa_scale) or args.yukawa_scale <= 0:
        raise ValueError("--yukawa-scale must be finite and positive")
    if not math.isfinite(args.wilson_factor) or args.wilson_factor == 0:
        raise ValueError("--wilson-factor must be finite and nonzero")
    if args.all and args.parameter_overrides:
        raise ValueError("--parameter-overrides cannot be used with --all")
    targets = base.ALL_MODELS if args.all else (args.model,)
    if not args.check_only and not base.absolute(args.target).is_file():
        raise FileNotFoundError(f"Oscillation target not found: {base.absolute(args.target)}")
    base_out = base.absolute(args.output_root) / args.run_id
    if not args.check_only and base_out.exists():
        raise FileExistsError(f"Run directory already exists: {base_out}")

    summary = {
        "format": FORMAT,
        "run_id": args.run_id,
        "beta_source_requested": args.beta_source,
        "models_requested": list(targets),
        "flags": {
            "no_one_scale_C5_substituted": True,
            "B_minus_1_and_C_minus_1_native_mixed_tree_integrated": True,
            "mixed_muF_muS_scalar_basis_identification_assumed": True,
            "independent_one_loop_hard_matching_complete": False,
            "absolute_spinor_normalization_certified": False,
        },
        "models": {},
    }
    failures = 0
    succeeded = 0
    blocked = 0
    for model in targets:
        try:
            entry = _preflight(model, args)
            if args.check_only:
                row = {
                    "status": "ready_to_execute",
                    "UV_beta_count": len(entry.uv_beta["report_betas"]),
                    "EFT1_beta_count": len(entry.eft_beta["report_betas"]),
                    "native_HHSS_channel_count": len(entry.quartics.by_name),
                    "saved_fit": str(entry.fit_root),
                    "model": model,
                }
                print(f"{model}: READY, UV/EFT1 betas "
                      f"{row['UV_beta_count']}/{row['EFT1_beta_count']}, "
                      f"native HHSS={row['native_HHSS_channel_count']}", flush=True)
            else:
                folder = base_out / base.model_slug(model)
                result = base.calculate_model(entry, args, folder)
                row = {
                    "status": "completed",
                    "path": str(folder / "result.json"),
                    "chi2": result["neutrinos"]["oscillation_chi2"],
                    "masses_eV": result["neutrinos"]["masses_ev_m1_m2_m3"],
                    "UV_nfev": result["rge"]["UV"]["nfev"],
                    "EFT1_nfev": result["rge"]["EFT1"]["nfev"],
                    "LLSS_shape": result["fermion_tree"]["shape"],
                    "HHSS_channels": result["scalar_matching"]["HHSS_channel_count"],
                }
                print(f"{model}: COMPLETED, chi2={row['chi2']:.9g}, "
                      f"masses(eV)={', '.join(f'{v:.7g}' for v in row['masses_eV'])}",
                      flush=True)
                succeeded += 1
        except NotImplementedError as exc:
            row = {"status": "unsupported_native_mixed", "reason": str(exc)}
            blocked += 1
            print(f"{model}: UNSUPPORTED — {exc}", flush=True)
        except Exception as exc:
            row = {"status": "failed_preflight_or_run", "reason": str(exc),
                   "type": type(exc).__name__}
            failures += 1
            print(f"{model}: FAILED — {type(exc).__name__}: {exc}", flush=True)
            if not args.check_only:
                error_dir = base_out / base.model_slug(model)
                error_dir.mkdir(parents=True, exist_ok=True)
                (error_dir / "error.json").write_text(
                    json.dumps(row, indent=2) + "\n", encoding="utf-8"
                )
        summary["models"][model] = row

    summary["counts"] = {
        "requested": len(targets), "completed": succeeded,
        "ready": sum(r["status"] == "ready_to_execute"
                     for r in summary["models"].values()),
        "unsupported": blocked, "failed": failures,
    }
    if not args.check_only:
        base_out.mkdir(parents=True, exist_ok=True)
        (base_out / "summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        print("Summary:", base_out / "summary.json")
    else:
        print("Preflight only; no RGEs or files generated.")
    print("Counts:", summary["counts"])
    return 1 if (failures or blocked) else 0


if __name__ == "__main__":
    raise SystemExit(main())
