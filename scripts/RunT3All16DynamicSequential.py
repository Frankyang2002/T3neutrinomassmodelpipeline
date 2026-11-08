"""Conditional all-model T3 fermion-first EFT -> both-scalars -> SMEFT runner.

Supported implementation: 14 *nonmixed* ordinary T3 models.  B_-1 and C_-1
are explicitly blocked: their mixed-doublet fermion boundary needs a separate
native dynamic source, which this script does not approximate.

For EACH point and model:
  UV full-general RGBeta running -> muF;
  integrate out ALL fermions at muF to current-point Matchete tree LLSS;
  remove fermions from the renormalisable state and run EFT1 RGBeta -> muS;
  fixed-one-loop LLSS -> C5 direct leading log using muS HHSS couplings;
  integrate out BOTH scalar multiplets in ONE finite scalar bubble at muS;
  SM + Weinberg evolution -> muLow -> neutrino observables and chi2.

This reuses existing code. It NEVER adds the saved one-scale C5.  For B_-3,
only, the separately available sequential Matchete hard remainder can be
included as an explicit alternative. No equivalent independent hard remainder
is available for the other 13 models; their reported result is the LLSS soft
sector, NOT a complete one-loop UV matching prediction. Hard/soft overlap and
the Matchete GammaCC/Weyl normalization remain unverified for ALL models.

The archived beta option retains exact saved RGBeta payloads WITHOUT claiming
current symbolic source fingerprints.  The original strict canonical loader is
never modified. All model failures are recorded, not substituted.

PowerShell from repository root (after extracting the ZIP):
  python -X utf8 -m scripts.RunT3All16DynamicSequential --model B_-3 --check-only
  python -X utf8 -m scripts.RunT3All16DynamicSequential --model B_-3 --run-id regression_v1
  python -X utf8 -m scripts.RunT3All16DynamicSequential --all --check-only
  python -X utf8 -m scripts.RunT3All16DynamicSequential --all --run-id batch_v1

Requires the existing Stage 1/2/3/4 B(-3) scripts and source-generated
Matchete/RGBeta data already in the repository. No Wolfram invocation and no
GitHub writes. Never overwrites an existing run folder.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from scripts import RunT3BMinus3SequentialSMEFT as stage1
from scripts import RunT3BMinus3DynamicSequential as stage2
from scripts import RunT3BMinus3DynamicLLSSSequential as stage3
from scripts import RunT3BMinus3DynamicScalarThreshold as stage4

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = ROOT / "output/general_sequential_pilot/all16_dynamic_scalar_threshold"
FIT_CAMPAIGN = ROOT / ("output/studies/general_scalar_extended/"
                       "thesis_general_all16/thesis16_general_c5_full_v1")
D0_FIT_CAMPAIGN = ROOT / ("output/studies/general_scalar_extended/"
                          "thesis_general_all16/continue_general16_v1")
BETA_ROOT = ROOT / ("output/studies/general_scalar_extended/"
                    "complete_scalar_symbolic")
FORMAT = "t3-14-nonmixed-dynamic-llss-scalar-matching-smeft-v1"

ALL_MODELS = (
    "A_-4", "A_-2", "A_+0", "B_-3", "B_-1", "B_+1",
    "C_-3", "C_-1", "C_+1", "D_-2", "D_+0", "D_+2",
    "E_-4", "E_-2", "E_+0", "E_+2",
)
MIXED = frozenset(("B_-1", "C_-1"))
DIMENSIONS = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
              "D": (3, 1, 2), "E": (3, 3, 2)}
HEAVY_NAMES = frozenset(("y1", "y2", "MF", "y1Conj", "y2Conj"))
BASE_NAMES = frozenset(("gY", "g2", "g3", "yu", "yd", "ye", "lambdaH",
                        "mS1Sq", "mS2Sq", "lambdaS1", "lambdaS2",
                        "lambdaH1", "lambdaH2", "lambda12", "lambdaT3"))
B3_HARD = stage2.DEFAULT_C5


def model_slug(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def absolute(path: Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def check_model(model: str) -> None:
    if model not in ALL_MODELS:
        raise ValueError(f"Unknown model {model!r}; choices: {ALL_MODELS}")
    if model in MIXED:
        raise NotImplementedError(
            f"{model}: mixed scalar S1/S2 fermion boundary requires its own "
            "dynamic Matchete factorisation; no nonmixed substitution is permitted"
        )


def saved_fit_folder(model: str, campaign_root: Path) -> Path:
    """Use D(+0)'s actual continuation fit, not the failed original fit."""
    folder = campaign_root / model_slug(model)
    if not (folder / "fitted_config.json").is_file() and model == "D_+0":
        if campaign_root.resolve() == FIT_CAMPAIGN.resolve():
            folder = D0_FIT_CAMPAIGN / model_slug(model)
    return folder


def fit_input(model: str, campaign_root: Path) -> tuple[Path, dict, dict, dict]:
    path = saved_fit_folder(model, campaign_root)
    config = stage1.read_object(path / "fitted_config.json")
    fit = stage1.read_object(path / "result.json")
    overrides = stage1.read_object(path / "general_uv_overrides.json")
    expected = {
        "d_s1": DIMENSIONS[model[0]][0],
        "d_s2": DIMENSIONS[model[0]][1],
        "d_f": DIMENSIONS[model[0]][2],
        "alpha": int(model.split("_")[1]),
        "shared_scalar": False,
    }
    if config.get("representation") != expected:
        raise ValueError(f"{model}: saved fitted representation differs from {expected}")
    if fit.get("model") != model:
        raise ValueError(f"{model}: saved fit model name mismatch")
    if fit.get("status") not in ("ExploratoryAcceptedFit", "ExploratoryPoorFit", "Success"):
        raise ValueError(f"{model}: saved result is not a numerical fit output: {fit.get('status')!r}")
    if not isinstance(fit.get("numerical_details"), dict):
        raise ValueError(f"{model}: saved fit has no numerical running snapshots")
    if not isinstance(overrides, dict) or not overrides:
        raise ValueError(f"{model}: saved fitted UV overrides missing")
    scales = config.get("scales", {})
    if not (float(config["base_state"]["mu_gev"]) >
            float(scales["mu_fermion_threshold_gev"]) >
            float(scales["mu_scalar_threshold_gev"]) >
            float(scales["mu_low_gev"]) > 0):
        raise ValueError(f"{model}: expected UV > F > S > low scale")
    return path, config, fit, overrides


def archived_beta_pair(model: str, run_id: str) -> tuple[dict, dict, dict]:
    """Read and validate each model's *archived* UV/EFT1 beta inventory.

    This is not a source-fingerprint certification. It must not be used to
    make claims of current physical Matchete/RGBeta basis equivalence.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_id):
        raise ValueError("Unsafe canonical-run-id")
    folder = BETA_ROOT / run_id / model_slug(model)
    files = {stage: folder / f"{stage.lower()}_rgbeta.json"
             for stage in ("UV", "EFT1")}
    uv = stage1.read_object(files["UV"])
    eft = stage1.read_object(files["EFT1"])
    keys = {}
    dimensions = DIMENSIONS[model[0]]
    for stage, payload in (("UV", uv), ("EFT1", eft)):
        meta = payload.get("metadata", {})
        if payload.get("status") != "SuccessCandidateMassBeta":
            raise ValueError(f"{model}/{stage}: unexpected RGBeta export status")
        expected_meta = {"dS1": dimensions[0], "dS2": dimensions[1],
                         "dF": dimensions[2], "alpha": int(model.split("_")[1]),
                         "GeneralScalarModel": model, "GeneralScalarStage": stage,
                         "SharedScalar": False,
                         "ReportBetaConvention":
                         "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"}
        if any(meta.get(k) != v for k, v in expected_meta.items()):
            raise ValueError(f"{model}/{stage}: archived beta metadata mismatch")
        beta = payload.get("report_betas")
        if not isinstance(beta, dict) or not beta or any(
            not isinstance(k, str) or not isinstance(expr, str) or not expr.strip()
            for k, expr in beta.items()
        ):
            raise ValueError(f"{model}/{stage}: missing/malformed beta expressions")
        keys[stage] = set(beta)
        additional = meta.get("GeneralScalarNewBetaKeys", [])
        if not isinstance(additional, list) or not set(additional) <= keys[stage]:
            raise ValueError(f"{model}/{stage}: missing extra scalar beta entries")
    if not BASE_NAMES <= keys["UV"] & keys["EFT1"]:
        raise ValueError(f"{model}: missing mandatory SM/scalar beta entries")
    if not {"MF", "y1", "y2"} <= keys["UV"]:
        raise ValueError(f"{model}: UV beta file missing heavy-fermion Yukawas/mass")
    if keys["EFT1"] & HEAVY_NAMES:
        raise ValueError(f"{model}: heavy fermion beta functions survive in EFT1")
    if keys["UV"] - HEAVY_NAMES != keys["EFT1"]:
        raise ValueError(f"{model}: unexpected UV-to-EFT1 beta inventory difference")
    return uv, eft, {
        "source": "archived_model_specific_RGBeta_pair",
        "model": model,
        "run_id": run_id,
        "uv_json": str(files["UV"].resolve()),
        "eft1_json": str(files["EFT1"].resolve()),
        "UV_beta_count": len(keys["UV"]),
        "EFT1_beta_count": len(keys["EFT1"]),
        "UV_removed_heavy_beta_keys": sorted(keys["UV"] - keys["EFT1"]),
        "uv_sha256": file_hash(files["UV"]),
        "eft1_sha256": file_hash(files["EFT1"]),
        "current_source_fingerprints_verified": False,
        "physical_matchete_rgbeta_basis_certified": False,
    }


def beta_pair(model: str, run_id: str, source: str) -> tuple[dict, dict, dict]:
    if source not in ("auto", "canonical", "archived"):
        raise ValueError("Invalid beta source")
    if model == "B_-3":
        return stage2.load_beta_pair(run_id, source)
    if source == "archived":
        return archived_beta_pair(model, run_id)
    from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
    try:
        uv, eft, provenance = load_canonical_general_exports(model, run_id)
    except ValueError as exc:
        if source == "canonical" or (
            "source fingerprint differs (or missing provenance)" not in str(exc)
        ):
            raise
        print(f"NOTICE: {model} canonical source fingerprint changed; using archived beta pair")
        uv, eft, provenance = archived_beta_pair(model, run_id)
        provenance["canonical_rejection"] = str(exc)
    else:
        provenance = dict(provenance)
        provenance["current_source_fingerprints_verified"] = True
    provenance["requested_beta_source"] = source
    return uv, eft, provenance


@dataclass(frozen=True)
class SourceFiles:
    boundary: Path
    channels: Path
    quartic_seed: Path
    registry: Path


def prepare_sources(model: str) -> tuple[SourceFiles, Any, Any, dict]:
    from scripts.RunT3GeneralSequentialLeadingLog import (
        DEFAULT_RELEASE, DEFAULT_V4, DEFAULT_SEEDS, DEFAULT_REGISTRY, source_paths,
    )
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_source
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import read_v4_cached_channels

    filemap = source_paths(
        {"release": DEFAULT_RELEASE, "v4": DEFAULT_V4,
         "seeds": DEFAULT_SEEDS, "running": ROOT / "output/general_observable/all16_v1",
         "registry": DEFAULT_REGISTRY},
        model,
    )
    paths = SourceFiles(**{name: filemap[name]
                            for name in ("boundary", "channels", "quartic_seed", "registry")})
    for name in ("boundary", "channels", "quartic_seed", "registry"):
        path = getattr(paths, name)
        if not path.is_file():
            raise FileNotFoundError(f"{model}: missing {name}: {path}")
    reference = load_source(paths.boundary, expected_model=model)
    quartics = read_v4_cached_channels(paths.channels, paths.quartic_seed, reference)
    registry = stage1.read_object(paths.registry)
    if model not in registry.get("models", {}):
        raise ValueError(f"{model}: scalar quartic registry entry absent")
    return paths, reference, quartics, registry


@dataclass
class Prepared:
    model: str
    fit_root: Path
    config: dict
    fit: dict
    overrides: dict
    uv_beta: dict
    eft_beta: dict
    beta_info: dict
    files: SourceFiles
    reference: Any
    quartics: Any
    registry: dict


def current_uv_overrides(saved: dict, changes: dict | None, scale: float) -> dict:
    """Scale ALL allowed UV Yukawa matrices, including extra conjugate channels."""
    values = stage2.build_uv_overrides(saved, changes, 1.0)
    if scale != 1.0:
        for name in ("y1", "y2", "y1Conj", "y2Conj"):
            if name in values:
                scaled = stage1.finite_matrix(values[name], name) * scale
                values[name] = stage2._matrix_to_entries(scaled)
    return values


def preflight(model: str, args: argparse.Namespace) -> Prepared:
    check_model(model)
    campaign = absolute(args.campaign_root)
    fit_root, config, fit, saved_overrides = fit_input(model, campaign)
    changes = None
    if args.parameter_overrides:
        if args.all:
            raise ValueError("--parameter-overrides is available only with --model")
        changes = stage1.read_object(absolute(args.parameter_overrides))
    overrides = current_uv_overrides(saved_overrides, changes, args.yukawa_scale)
    uv_beta, eft_beta, info = beta_pair(model, args.canonical_run_id, args.beta_source)
    missing = set(overrides) - set(uv_beta.get("report_betas", {}))
    if missing:
        raise ValueError(f"{model}: UV override keys absent from beta export: {sorted(missing)}")
    paths, reference, quartics, registry = prepare_sources(model)
    if model != "B_-3" and args.b3_hard_policy != "omit" and not args.all:
        raise ValueError(f"{model}: no model-specific sequential Matchete hard C5 exists")
    if model == "B_-3" and not B3_HARD.is_file():
        raise FileNotFoundError(f"B_-3 saved sequential Matchete hard term missing: {B3_HARD}")
    return Prepared(model, fit_root, config, fit, overrides, uv_beta, eft_beta,
                    info, paths, reference, quartics, registry)


def fresh_running(prepared: Prepared) -> dict:
    """Recalculate the full-general renormalisable trajectory from the UV."""
    from Numerical.orchestration.GeneralRunningBenchmark import run
    cfg = prepared.config
    scales = cfg["scales"]
    result = run(
        prepared.uv_beta, prepared.eft_beta,
        mu_uv=float(cfg["base_state"]["mu_gev"]),
        mu_f=float(scales["mu_fermion_threshold_gev"]),
        mu_s=float(scales["mu_scalar_threshold_gev"]),
        overrides=prepared.overrides,
        rtol=1e-7, atol=1e-10,
    )
    if result.get("status") != "NumericalRenormalisableRunningOnly" or (
        result.get("model") != prepared.model
    ):
        raise ValueError(f"{prepared.model}: renormalisable RGE returned wrong model/status")
    for k in ("uv_solver", "eft1_solver"):
        if result.get(k, {}).get("success") is not True or (
            int(result[k].get("nfev", 0)) <= 0
        ):
            raise RuntimeError(f"{prepared.model}: fresh {k} did not complete")
    return result


def dynamic_soft_sector(prepared: Prepared, running: dict,
                        *, factor: float, check_scale: bool) -> tuple[Any, Any, dict, dict]:
    """Dynamically assemble C_LLSS and its one-loop scalar-contracted C5."""
    from Numerical.running.GeneralT3DynamicFermionWilson import nonmixed_boundary
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot
    from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors, response

    model = prepared.model
    cfg = prepared.config
    mu_f = float(cfg["scales"]["mu_fermion_threshold_gev"])
    mu_s = float(cfg["scales"]["mu_scalar_threshold_gev"])
    at_f = running["uv_at_fermion_threshold"]
    at_s = running["eft1_at_scalar_threshold"]
    heavy = stage1.fitted_fermion_state(at_f, cfg, mu_f)
    source, source_info = nonmixed_boundary(model, heavy, repo_root=ROOT, mu_f=mu_f)
    if tuple(source.labels) != tuple(prepared.reference.labels):
        raise ValueError(f"{model}: F-tree dynamic and reference scalar bases differ")
    n = 2 * sum(DIMENSIONS[model[0]][:2])
    if source.raw.shape != (3, 3, 2, 2, n, n):
        raise ValueError(f"{model}: incorrect C_LLSS shape {source.raw.shape}")
    if not np.isfinite(source.raw).all():
        raise ValueError(f"{model}: nonfinite C_LLSS tensor")
    snapshot = RunningSnapshot(
        model=model,
        source=Path("fresh_full_general_EFT1_endpoint"),
        sha256=sha256(json.dumps(at_s, sort_keys=True).encode()).hexdigest(),
        mu_f=mu_f, mu_s=mu_s,
        scalar_mass1_sq=stage1.finite_real(at_s["mS1Sq"], "mS1Sq"),
        scalar_mass2_sq=stage1.finite_real(at_s["mS2Sq"], "mS2Sq"),
        couplings=at_s,
        evidence={"fresh_full_general_EFT1": True},
    )
    bundle = scalar_tensors(
        model, source, prepared.quartics, snapshot, prepared.registry,
        include_general_extra=True, include_portal_completion=True,
    )
    result = response(
        source, bundle, mu_f=mu_f, mu_s=mu_s,
        spinor_factor=complex(factor), check_scale_cancellation=check_scale,
    )
    return source, bundle, result, source_info


def symmetric_c5(raw: Any, description: str) -> np.ndarray:
    c = np.asarray(stage1.decode(raw), dtype=complex)
    if c.shape != (3, 3) or not np.isfinite(c).all() or not np.allclose(
        c, c.T, rtol=1e-8, atol=1e-28
    ):
        raise ValueError(f"{description}: C5 not finite/symmetric")
    return 0.5 * (c + c.T)


def c5_pieces(response: dict, model: str, running: dict,
              config: dict) -> tuple[dict[str, Any], np.ndarray, np.ndarray | None]:
    direct = symmetric_c5(response["c5_eft1_leading_log"], "direct LLSS log")
    bubble = symmetric_c5(response["c5_scalar_finite_bubble"], "joint scalar bubble")
    soft = symmetric_c5(response["c5_LLSS_total"], "total LLSS matching")
    diff = float(np.linalg.norm(soft - direct - bubble))
    denominator = max(float(np.linalg.norm(soft)),
                      float(np.linalg.norm(direct)),
                      float(np.linalg.norm(bubble)), 1e-30)
    if diff > 1e-9 * denominator + 1e-28:
        raise ValueError(f"LLSS sum does not equal RG log + scalar bubble: {diff}")
    pieces: dict[str, Any] = {
        "direct_log_gev_inverse": stage1.encode_matrix(direct),
        "both_scalar_finite_bubble_gev_inverse": stage1.encode_matrix(bubble),
        "LLSS_soft_total_gev_inverse": stage1.encode_matrix(soft),
        "direct_log_norm": float(np.linalg.norm(direct)),
        "scalar_bubble_norm": float(np.linalg.norm(bubble)),
        "LLSS_soft_total_norm": float(np.linalg.norm(soft)),
        "sum_identity_residual": diff,
        "direct_log_inserted_once": True,
        "scalar_bubble_inserted_once": True,
        "one_scale_Matchete_C5_added": False,
        "hard_matching_independent_for_this_model": False,
    }
    hard: np.ndarray | None = None
    if model == "B_-3":
        hard, saved_log, detail = stage3.saved_threshold_components(
            B3_HARD, running, config
        )
        # The old B_-3 direct log is diagnostic only: never added.
        pieces["saved_B3_analytic_log_not_added"] = stage1.encode_matrix(saved_log)
        pieces["saved_B3_hard_candidate"] = stage1.encode_matrix(hard)
        pieces["saved_B3_hard_source"] = detail
    else:
        pieces["saved_independent_hard"] = None
        pieces["hard_status"] = "not_calculated_for_this_model"
    return pieces, soft, hard


def calculate_model(prepared: Prepared, args: argparse.Namespace,
                    output_dir: Path) -> dict:
    running = fresh_running(prepared)
    source, bundle, response, source_info = dynamic_soft_sector(
        prepared, running, factor=args.wilson_factor,
        check_scale=not args.skip_scale_check,
    )
    pieces, soft, hard = c5_pieces(
        response, prepared.model, running, prepared.config
    )
    chosen = soft
    if args.b3_hard_policy == "retain":
        if hard is None:
            raise RuntimeError("No independent hard candidate for this model")
        chosen = soft + hard
    physics = stage3.smeft_observables(
        chosen, running, prepared.config, absolute(args.target)
    )
    alternative_physics = None
    if prepared.model == "B_-3":
        alternative_c5 = soft + hard if args.b3_hard_policy == "omit" else soft
        alternative_physics = stage3.smeft_observables(
            alternative_c5, running, prepared.config, absolute(args.target)
        )
    result = {
        "format": FORMAT,
        "model": prepared.model,
        "status": "ConditionalDynamicFToBothScalarsToSMEFTCompleted",
        "scales_gev": {
            "UV": float(prepared.config["base_state"]["mu_gev"]),
            **running["scales_gev"],
            "low": float(prepared.config["scales"]["mu_low_gev"]),
        },
        "rge": {
            "UV": running["uv_solver"],
            "EFT1": running["eft1_solver"],
            "renormalisable_parameters_at_F": running["uv_at_fermion_threshold"],
            "renormalisable_parameters_at_S": running["eft1_at_scalar_threshold"],
        },
        "fermion_tree": {
            "model": prepared.model,
            "computed_from_current_UV_at_F": True,
            "source_geometry": source_info,
            "shape": list(source.raw.shape),
            "norm_gev_inverse": float(np.linalg.norm(source.raw)),
            "source_count": len(source.source_indices),
            "source_hash": source.hash,
        },
        "scalar_matching": {
            "HHSS_channel_names": list(bundle.names),
            "HHSS_channel_count": len(bundle.names),
            "HHSS_tensor_provenance": bundle.provenance,
            "mass_eigenvalues_sq_gev2": np.linalg.eigvalsh(bundle.mass2).tolist(),
            "weinberg_beta_projection_residual":
                response["weinberg_beta_projection_residual"],
            "weinberg_bubble_projection_residual":
                response["weinberg_bubble_projection_residual"],
            "scale_cancellation_residual":
                response["matching_scale_cancellation_residual"],
            **pieces,
            "C5_used_at_S_gev_inverse": stage1.encode_matrix(chosen),
            "primary_hard_policy": args.b3_hard_policy if prepared.model == "B_-3"
                                   else "unavailable_soft_only",
            "B3_alternative": None if alternative_physics is None else {
                "hard_policy": "retain" if args.b3_hard_policy == "omit" else "omit",
                "neutrinos": alternative_physics["neutrinos"],
                "smeft": alternative_physics["smeft"],
            },
        },
        **physics,
        "reference": {
            "fit_root": str(prepared.fit_root.resolve()),
            "original_one_scale_fit_chi2": prepared.fit.get("chi2_final"),
            "archived_beta_input": prepared.beta_info,
            "dynamic_source_geometry": str(source.path),
            "quartic_source_cache": str(prepared.files.channels),
            "quartic_seed": str(prepared.files.quartic_seed),
            "scalar_registry": str(prepared.files.registry),
            "assumed_wilson_factor": float(args.wilson_factor),
        },
        "flags": {
            "UV_evolved_fresh": True,
            "all_heavy_fermions_integrated_out_at_F": True,
            "EFT1_RGE_contains_no_heavy_fermions": True,
            "LLSS_F_boundary_evaluated_from_current_Yukawas": True,
            "LLSS_to_C5_direct_fixed_one_loop_log": True,
            "all_scalar_components_integrated_out_at_S": True,
            "SMEFT_evolved_to_low_scale": True,
            "full_LLSS_seven_term_transport_evolved": False,
            "other_EFT1_operator_classes_evolved": False,
            "independent_hard_match_calculated_all_models": False,
            "spinor_normalization_verified": False,
            "scalar_basis_map_verified": False,
            "B3_saved_hard_independence_verified": False,
            "oscillation_fit_performed": True,
        },
        "scope": (
            "Exploratory/conditional full-general renormalisable RGEs, dynamic Matchete "
            "fermion-tree LLSS, one-loop direct LLSS->C5 mixing, simultaneous "
            "two-scalar finite matching, SM+Weinberg RGEs. Independent hard "
            "matching unavailable outside B_-3; physical normalization unverified."
        ),
    }
    # Save only after all numerical work and JSON serialization has succeeded.
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    obs = result["neutrinos"]
    lines = ["quantity,value,unit"]
    for key, value in obs["oscillation_observables"].items():
        unit = "dimensionless" if key.startswith("sin2") else "eV^2"
        lines.append(f"{key},{float(value):.17g},{unit}")
    for i, m in enumerate(obs["masses_ev_m1_m2_m3"], 1):
        lines.append(f"m{i},{float(m):.17g},eV")
    lines.append(f"chi2,{float(obs['oscillation_chi2']):.17g},dimensionless")
    (output_dir / "observables.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    np.savez_compressed(
        output_dir / "dynamic_matching.npz",
        model=np.asarray(prepared.model),
        C_LLSS_muF=source.raw,
        C_LLSS_by_source_muF=source.by_source,
        HHSS_at_muS=bundle.quartic,
        scalar_mass2_at_muS=bundle.mass2,
        C5_direct_log=symmetric_c5(response["c5_eft1_leading_log"], "direct"),
        C5_scalar_bubble=symmetric_c5(response["c5_scalar_finite_bubble"], "bubble"),
        C5_LLSS_soft_total=soft,
        C5_primary=chosen,
        muF_gev=np.asarray(float(prepared.config["scales"]["mu_fermion_threshold_gev"])),
        muS_gev=np.asarray(float(prepared.config["scales"]["mu_scalar_threshold_gev"])),
    )
    if alternative_physics is not None:
        stage4._alternative_csv(output_dir / "observables_alternative.csv", alternative_physics)
    return result


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    choice = p.add_mutually_exclusive_group()
    choice.add_argument("--model", choices=ALL_MODELS, default="B_-3")
    choice.add_argument("--all", action="store_true")
    p.add_argument("--run-id", default="multi_v1")
    p.add_argument("--check-only", action="store_true")
    p.add_argument("--campaign-root", type=Path, default=FIT_CAMPAIGN)
    p.add_argument("--canonical-run-id", default="full16_general_v1")
    p.add_argument("--beta-source", choices=("archived", "canonical", "auto"),
                   default="archived")
    p.add_argument("--yukawa-scale", type=float, default=1.0)
    p.add_argument("--parameter-overrides", type=Path,
                   help="Partial overrides file; only with --model")
    p.add_argument("--target", type=Path,
                   default=ROOT / "configs/nufit_6_1_ic24_no.json")
    p.add_argument("--wilson-factor", type=float, default=1.0)
    p.add_argument("--b3-hard-policy", choices=("omit", "retain"), default="omit")
    p.add_argument("--skip-scale-check", action="store_true")
    p.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", args.run_id) is None:
        raise ValueError("Unsafe --run-id: use 1-64 alphanumeric, underscore or hyphen")
    if not math.isfinite(args.yukawa_scale) or args.yukawa_scale <= 0:
        raise ValueError("--yukawa-scale must be finite and positive")
    if not math.isfinite(args.wilson_factor) or args.wilson_factor == 0:
        raise ValueError("--wilson-factor must be finite and nonzero")
    if args.all and args.parameter_overrides:
        raise ValueError("--parameter-overrides cannot be used with --all")
    targets = ALL_MODELS if args.all else (args.model,)
    if not args.check_only and not absolute(args.target).is_file():
        raise FileNotFoundError(f"Oscillation target not found: {absolute(args.target)}")
    base_out = absolute(args.output_root) / args.run_id
    if not args.check_only and base_out.exists():
        raise FileExistsError(f"Run directory already exists: {base_out}")

    summary = {
        "format": FORMAT,
        "run_id": args.run_id,
        "beta_source_requested": args.beta_source,
        "models_requested": list(targets),
        "flags": {"no_one_scale_C5_substituted": True,
                  "B_minus_1_and_C_minus_1_native_mixed_tree_not_implemented": True},
        "models": {},
    }
    failures = 0
    succeeded = 0
    blocked = 0
    for model in targets:
        try:
            entry = preflight(model, args)
            if args.check_only:
                row = {"status": "ready_to_execute",
                       "UV_beta_count": len(entry.uv_beta["report_betas"]),
                       "EFT1_beta_count": len(entry.eft_beta["report_betas"]),
                       "native_HHSS_channel_count": len(entry.quartics.by_name),
                       "saved_fit": str(entry.fit_root),
                       "model": model}
                print(f"{model}: READY, UV/EFT1 betas "
                      f"{row['UV_beta_count']}/{row['EFT1_beta_count']}, "
                      f"native HHSS={row['native_HHSS_channel_count']}", flush=True)
            else:
                folder = base_out / model_slug(model)
                result = calculate_model(entry, args, folder)
                row = {"status": "completed", "path": str(folder / "result.json"),
                       "chi2": result["neutrinos"]["oscillation_chi2"],
                       "masses_eV": result["neutrinos"]["masses_ev_m1_m2_m3"],
                       "UV_nfev": result["rge"]["UV"]["nfev"],
                       "EFT1_nfev": result["rge"]["EFT1"]["nfev"],
                       "LLSS_shape": result["fermion_tree"]["shape"],
                       "HHSS_channels": result["scalar_matching"]["HHSS_channel_count"]}
                print(f"{model}: COMPLETED, chi2={row['chi2']:.9g}, "
                      f"masses(eV)={', '.join(f'{v:.7g}' for v in row['masses_eV'])}",
                      flush=True)
                succeeded += 1
        except NotImplementedError as exc:
            row = {"status": "unsupported_mixed_boundary", "reason": str(exc)}
            blocked += 1
            print(f"{model}: UNSUPPORTED — {exc}", flush=True)
        except (Exception) as exc:
            row = {"status": "failed_preflight_or_run", "reason": str(exc),
                   "type": type(exc).__name__}
            failures += 1
            print(f"{model}: FAILED — {type(exc).__name__}: {exc}", flush=True)
            if not args.check_only:
                error_dir = base_out / model_slug(model)
                if not error_dir.exists():
                    error_dir.mkdir(parents=True, exist_ok=False)
                (error_dir / "error.json").write_text(
                    json.dumps(row, indent=2) + "\n", encoding="utf-8"
                )
        summary["models"][model] = row

    summary["counts"] = {
        "requested": len(targets), "completed": succeeded,
        "ready": sum(r["status"] == "ready_to_execute"
                     for r in summary["models"].values()),
        "unsupported": blocked,
        "failed": failures,
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
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
