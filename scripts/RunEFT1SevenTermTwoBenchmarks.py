"""Rough two-model, seven-term EFT1 ODE followed by ONE joint scalar threshold.

Run from the T3neutrinomassmodelpipeline repository root (PowerShell)::

    python -X utf8 scripts/RunEFT1SevenTermTwoBenchmarks.py --steps 2

UV -> integrate out ALL heavy fermions at M_F -> EFT1(H,S1,S2) ->
integrate out BOTH heavy scalar multiplets together at M_S -> SMEFT.

Uses the repository's existing source-backed all-16 boundary and 45-Weyl,
12-real-component reduced-quartic inputs, but changes its Wilson evolution
from the Fonseca six-term *candidate* to the seven-term Eq. (4.85)-based
GeneralEFT1SevenTermFlow.master_beta candidate (with independently symmetrised
fermion/scalar pairs in its crossed channel).  In particular this beta
contains the antilinear C.conj() channel.  This does NOT independently
certify the native Matchete GammaCC signed normalization or the complete
model-specific scalar quartic.  Reports are explicitly provisional.

The matching step is C5(MS) = C5_direct_from_ODE(MS) +
C5_joint_heavy_scalar_bubble(MS).  No legacy c5_eft1_leading_log is added.
All inputs are read only, and existing outputs are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running import All16ConditionalFullWilsonTrial as legacy
from Numerical.joint_scalar_matching.core import scalar_bubble_matching
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, master_beta, total_beta, LOOP,
)

MODELS = ("D_-2", "A_-4")
FORMAT = "t3-two-model-seven-term-eft1-joint-scalar-rough-v1"


def require(ok: bool, why: str) -> None:
    if not ok:
        raise ValueError(why)


def absolute(root: Path, relative: Path) -> Path:
    return relative if relative.is_absolute() else root / relative


def encode(z: np.ndarray) -> dict[str, list]:
    z = np.asarray(z, dtype=complex)
    require(np.isfinite(z).all(), "Non-finite coefficient")
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def seven_term_midpoint(
    c0: np.ndarray,
    y_f: np.ndarray,
    y_s: np.ndarray,
    state_f: Mapping[str, Any],
    state_s: Mapping[str, Any],
    hhss: np.ndarray,
    fermion_generators: Mapping[str, tuple],
    scalar_generators: Mapping[str, tuple],
    mu_f: float,
    mu_s: float,
    steps: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """RK2 over ln(mu), full pair-symmetric C_ijab; no projected-only ODE.

    g_Y, g_2, g_3, Higgs Yukawas and lambdaH use the existing two running
    snapshots.  HHSS is supplied at muS and frozen; omitted quartics are
    explicitly zero in this *reduced-potential* numerical trial.
    """
    require(1 <= steps <= 64, "--steps must be between 1 and 64")
    require(0 < mu_s < mu_f, "Require fermion threshold M_F > scalar threshold M_S > 0")
    n = c0.shape[-1]
    require(c0.shape == (45, 45, n, n) and n in (12, 16), "Wrong Wilson shape")
    qhhss = legacy.hhss_full_quartic(hhss)
    qh = legacy.higgs_quartic(n)
    dt = math.log(mu_s / mu_f) / steps
    c = np.asarray(c0, dtype=complex).copy()
    evaluations = 0
    term_norms_start: dict[str, float] = {}

    def rhs(coefficient: np.ndarray, fraction: float) -> np.ndarray:
        nonlocal evaluations, term_norms_start
        y = (1 - fraction) * y_f + fraction * y_s
        h = ((1 - fraction) * float(state_f["lambdaH"])
             + fraction * float(state_s["lambdaH"]))
        q = qhhss + h * qh
        gauge = []
        for group, param in (("SU3", "g3"), ("SU2", "g2"), ("U1", "gY")):
            g = ((1 - fraction) * float(state_f[param])
                 + fraction * float(state_s[param]))
            gauge.append((g, tuple(zip(fermion_generators[group], scalar_generators[group]))))
        terms = master_beta(coefficient, Couplings(y, q, tuple(gauge)))
        if evaluations == 0:
            term_norms_start = {name: float(np.linalg.norm(val)) for name, val in terms.items()}
        evaluations += 1
        return LOOP * total_beta(terms)

    for i in range(steps):
        fraction = i / steps
        k1 = rhs(c, fraction)
        k2 = rhs(c + 0.5 * dt * k1, (i + 0.5) / steps)
        c += dt * k2
        require(np.isfinite(c).all(), f"Non-finite Wilson tensor at step {i+1}")

    denom = max(float(np.max(np.abs(c))), 1e-100)
    pair_residual = max(
        float(np.max(np.abs(c - c.swapaxes(0, 1)))),
        float(np.max(np.abs(c - c.swapaxes(2, 3)))),
    ) / denom
    require(pair_residual < 1e-8, f"Pair symmetry lost: {pair_residual:.4g}")
    return c, {
        "beta": "Eq. (4.85)-based seven-term candidate; crossed channel pair-symmetrised",
        "implementation": "GeneralEFT1SevenTermFlow.master_beta(corrected_crossed=True)",
        "has_conjugate_C_channel": True,
        "full_beta_independently_loop_certified": False,
        "method": "full-tensor explicit midpoint RK2, ln(mu) interpolation",
        "steps": steps,
        "rhs_evaluations": evaluations,
        "mu_F_gev": mu_f,
        "mu_S_gev": mu_s,
        "C_relative_change": float(np.linalg.norm(c - c0) / max(np.linalg.norm(c0), 1e-100)),
        "pair_symmetry_max_relative_residual": pair_residual,
        "first_point_seven_term_beta_norms": term_norms_start,
        "scalar_quartic_coverage": "lambdaH(HdagH)^2/2 and complete supplied HHSS only",
        "additional_quartics_zero_by_rough_approximation": True,
    }


def match_once(c_ms: np.ndarray, background: Any, vev: float) -> dict[str, Any]:
    """Project ODE C_LLHH and one joint eight-real-scalar B0bar bubble.

    Report rather than suppress possible Weinberg projection residuals: the
    seven-term paper template's crossed channel has not independently passed
    gauge covariance on these specific T3 tensors.
    """
    hh = legacy.generic_leptons(c_ms, heavy=False)
    heavy = legacy.generic_leptons(c_ms, heavy=True)
    direct, direct_res = legacy.project_master(hh, strict=False)
    bubble = scalar_bubble_matching(
        heavy, background.quartic, background.mass2, background.mu_s
    )
    bubble_c5, bubble_res = legacy.project_master(bubble["amplitude"], strict=False)
    require(len(bubble["pair_data"]) == 36, "Joint eight-scalar matching lacks 36 unordered pairs")
    total = direct + bubble_c5
    require(np.allclose(total, total.T, rtol=1e-9, atol=1e-27), "Non-symmetric C5")
    masses = legacy.neutrino_mass_proxy(total, vev)
    return {
        "c5_from_eft1_ode_gev_inv": encode(direct),
        "c5_from_joint_scalar_bubble_gev_inv": encode(bubble_c5),
        "c5_smeft_conditional_gev_inv": encode(total),
        "mnu_matrix_eV": encode(-0.5 * vev * vev * 1e9 * total),
        "takagi_masses_eV": masses["takagi_singular_values_eV"],
        "scalar_mass_eigenvalues_gev": np.sqrt(np.linalg.eigvalsh(background.mass2)).tolist(),
        "joint_scalar_eigenpairs": len(bubble["pair_data"]),
        "direct_weinberg_projection_residual": float(direct_res),
        "bubble_weinberg_projection_residual": float(bubble_res),
        "weinberg_projection_is_accurate_1e-6": bool(max(direct_res, bubble_res) <= 1e-6),
        "no_separate_eft1_fixed_order_log_added": True,
        "heavy_scalar_multiplets_integrated_jointly": True,
    }


def run_one(root: Path, model: str, *, release_dir: Path, v6_dir: Path,
            running_dir: Path, steps: int, vev: float,
            full_tensors: Path | None = None) -> dict[str, Any]:
    require(model in MODELS, "This runner is limited to D_-2 and A_-4")
    paths = legacy.load_release_manifest(absolute(root, release_dir))
    source = legacy.load_source(paths[model], expected_model=model)
    bgpath = absolute(root, v6_dir) / "numeric_backgrounds" / (legacy.slug(model) + ".json")
    background = legacy.load_background(bgpath, source)
    runpath = absolute(root, running_dir) / legacy.slug(model) / "running.json"
    running = legacy.read_json(runpath)
    require(running.get("model") == model, "Wrong EFT1 running model")
    state_f, state_s = running["eft1_at_fermion_threshold"], running["eft1_at_scalar_threshold"]
    mu_f, mu_s = float(source.mu_f_gev), float(background.mu_s)
    require(abs(float(running["scales_gev"]["F"]) - mu_f) < 1e-6 and
            abs(float(running["scales_gev"]["S"]) - mu_s) < 1e-6,
            "Coupling snapshots do not match fermion/scalar thresholds")
    n = 4 + len(source.labels)
    require(n == 12 and background.mass2.shape == (8, 8),
            "Expected Higgs + eight real heavy scalars for these models")

    if model in legacy.SPECIAL:
        snapshot = legacy.load_running_snapshot(runpath, model)
        native_to_source = legacy.build_special_basis(model, snapshot, source)["R_original_from_boundary"]
    else:
        native_to_source = np.eye(8)
    scalar_generators = legacy.rotated_scalar_generators(model, native_to_source @ background.U)
    generator_data = legacy.eft.gauge_generators()
    fermion_generators = {g: generator_data[g]["fermion"] for g in ("SU3", "SU2", "U1")}
    c0 = legacy.embedded_boundary(source.raw, background.U)

    # Cheap checks of already implemented source/basis geometry.
    source_matching = np.einsum(
        "pqijab,ak,bl->pqijkl", source.raw, background.U, background.U, optimize=True
    )
    source_ward = legacy.source_ward_residual(source_matching, scalar_generators)
    quartic_ward = legacy.quartic_ward_residual(
        legacy.hhss_full_quartic(background.quartic)
        + float(state_s["lambdaH"]) * legacy.higgs_quartic(n), scalar_generators
    )
    require(max(source_ward.values()) < 2e-9 and max(quartic_ward.values()) < 2e-9,
            "Mismatched source/scalar generator basis")

    y_f = legacy.yukawa_at_endpoint(state_f, n)
    y_s = legacy.yukawa_at_endpoint(state_s, n)
    end, ode = seven_term_midpoint(
        c0, y_f, y_s, state_f, state_s, background.quartic,
        fermion_generators, scalar_generators, mu_f, mu_s, steps
    )
    matching = match_once(end, background, vev)
    report: dict[str, Any] = {
        "format": FORMAT, "model": model,
        "status": "ROUGH_NUMERICAL_SEVEN_TERM_ODE_JOINT_SCALAR_MATCHING",
        "threshold_sequence": "UV -> remove ALL fermions -> EFT1(H,S1,S2) -> jointly remove ALL scalars -> SMEFT",
        "n_SM_Weyl": 45, "n_real_scalar": n,
        "independent_complex_C_coordinates": (45 * 46 // 2) * (n * (n + 1) // 2),
        "master_boundary_factor": legacy.SCALAR_CONVENTION_FACTOR,
        "boundary_factor_status": "declared legacy-to-master projector convention, NOT native GammaCC amplitude certification",
        "source_sha256": source.hash,
        "running_json_sha256": legacy.digest(runpath),
        "scalar_background_sha256": legacy.digest(bgpath),
        "source_ward_max": max(source_ward.values()),
        "quartic_ward_max": max(quartic_ward.values()),
        "evolution": ode, "matching": matching,
        "signed_GammaCC_to_Weyl_factor_independently_known": False,
        "full_scalar_quartic_included": False,
        "hard_one_loop_fermion_matching_included": False,
        "physical_neutrino_mass_prediction_certified": False,
    }
    if full_tensors is not None:
        full_tensors.parent.mkdir(parents=True, exist_ok=True)
        require(not full_tensors.exists(), f"Refusing to replace {full_tensors}")
        np.savez_compressed(
            full_tensors, full_C_at_MF_gev_inv=c0, full_C_at_MS_gev_inv=end,
            C5_SMEFT_gev_inv=(np.asarray(matching["c5_smeft_conditional_gev_inv"]["real"])
                + 1j * np.asarray(matching["c5_smeft_conditional_gev_inv"]["imag"])),
        )
        report["full_tensors_relative_path"] = str(full_tensors.name)
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=2, help="2 midpoint steps by default")
    parser.add_argument("--output-root", type=Path,
                        default=Path("output/general_observable/eft1_seven_term_Dm2_Am4_rough_v1"))
    parser.add_argument("--release-dir", type=Path, default=legacy.DEFAULT_RELEASE)
    parser.add_argument("--v6-dir", type=Path, default=legacy.DEFAULT_V6)
    parser.add_argument("--running-dir", type=Path, default=legacy.DEFAULT_RUNNING)
    parser.add_argument("--vev", type=float, default=246.22)
    parser.add_argument("--save-tensors", action="store_true")
    args = parser.parse_args(argv)
    require(1 <= args.steps <= 64 and args.vev > 0, "Invalid --steps or --vev")
    out = absolute(ROOT, args.output_root)
    require(not out.exists(), f"Will not overwrite existing output directory: {out}")
    out.mkdir(parents=True)
    rows = []
    for model in MODELS:
        print(f"[seven-term EFT1] {model}", flush=True)
        try:
            tensor_path = out / "tensors" / (legacy.slug(model) + ".npz") if args.save_tensors else None
            result = run_one(ROOT, model, release_dir=args.release_dir,
                             v6_dir=args.v6_dir, running_dir=args.running_dir,
                             steps=args.steps, vev=args.vev, full_tensors=tensor_path)
            model_path = out / "per_model" / (legacy.slug(model) + ".json")
            model_path.parent.mkdir(parents=True, exist_ok=True)
            model_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
            masses = result["matching"]["takagi_masses_eV"]
            rows.append({"model": model, "status": "COMPLETE", "m1_eV": masses[0],
                         "m2_eV": masses[1], "m3_eV": masses[2],
                         "C_relative_change": result["evolution"]["C_relative_change"],
                         "max_weinberg_projection_residual": max(
                             result["matching"]["direct_weinberg_projection_residual"],
                             result["matching"]["bubble_weinberg_projection_residual"]),
                         "error": ""})
            print(f"  approximate masses (eV): {masses}", flush=True)
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
            rows.append({"model": model, "status": "FAILED", "m1_eV": None,
                         "m2_eV": None, "m3_eV": None, "C_relative_change": None,
                         "max_weinberg_projection_residual": None,
                         "error": f"{type(exc).__name__}: {exc}"})
            print(f"  failed: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
    with (out / "comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    summary = {
        "format": FORMAT, "status": "COMPLETE" if all(r["status"] == "COMPLETE" for r in rows) else "PARTIAL",
        "models": rows, "steps": args.steps,
        "source": "existing all16 fermion-first boundary and reduced HHSS/Higgs quartic",
        "running": "seven-term full complex Wilson tensor, including C*",
        "matching": "single joint threshold at M_S; no duplicate fixed-order EFT1 log",
        "scalar_quartic_complete": False,
        "Matchete_to_Weyl_sign_certified": False,
        "physical_neutrino_mass_predictions": False,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print("Summary:", out / "summary.json", flush=True)
    return 0 if summary["status"] == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
