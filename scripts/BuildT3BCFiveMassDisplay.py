"""Workstream B: PRESENTATION ONLY for existing B(-1)/C(-1) UV -> EFT1 runs.

Reads the completed native-RGBeta result JSONs. Does not import the RGE
pipeline, evaluate beta functions, integrate ODEs, or modify source reports.

From the repository root:
  python -X utf8 -m scripts.BuildT3BCFiveMassDisplay
  python -X utf8 -m scripts.BuildT3BCFiveMassDisplay --no-plots

Writes a Markdown comparison, compact JSON summary, endpoint and trajectory
CSVs, and seven separate figures. No physical threshold matching, Matchete /
RGBeta basis certification, or neutrino-mass prediction is implied.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "output/workstream_b_bc_coupled_uv_eft1"
OUTPUT_DIR = ROOT / "Reports/output/workstream_b_bc_coupled_uv_eft1"
SOURCE_FORMAT = "workstream-b-bc-coupled-uv-through-eft1-five-mass-v1"
SOURCE_STATUS = "CandidateCoupledUVAndEFT1FiveRealMassesIntegrated"
OUTPUT_FORMAT = "workstream-b-bc-five-mass-display-only-v1"
CONVENTION = "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
MODELS = ("B_-1", "C_-1")
MIXED = "muT3GeneralS1S2"
MASS_SYMBOLS = ("mHSq", "mS1Sq", "mS2Sq", MIXED)
MASS_COORDS = ("mHSq", "mS1Sq", "mS2Sq", f"Re({MIXED})", f"Im({MIXED})")
# Output keys, labels, and filenames. Every quantity has units GeV^2.
QUANTITIES = (
    ("mH2", r"$m_H^2$", "mH2"),
    ("mS12", r"$m_1^2$", "mS1Sq"),
    ("mS22", r"$m_2^2$", "mS2Sq"),
    ("re_mu", r"$\mathrm{Re}\,\mu_R$", "ReMu"),
    ("im_mu", r"$\mathrm{Im}\,\mu_R$", "ImMu"),
    ("m_minus2", r"$m_-^2$", "mMinusSq"),
    ("m_plus2", r"$m_+^2$", "mPlusSq"),
)
KEYS = tuple(q[0] for q in QUANTITIES)


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError("B/C five-mass display: " + reason)


def finite(value: Any, field: str) -> float:
    require(type(value) in (int, float), f"{field} must be a JSON number")
    result = float(value)
    require(math.isfinite(result), f"{field} is not finite")
    return result


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=2e-10, abs_tol=2e-8)


def snapshot(raw: Any, *, model: str, stage: str, index: int) -> dict[str, Any]:
    """Extract only already-saved mass information; perform no physics calculation."""
    tag = f"{model}/{stage}[{index}]"
    require(isinstance(raw, dict), f"{tag} must be an object")
    mixed = raw.get("muT3GeneralS1S2_gev2")
    spectrum = raw.get("scalar_spectrum")
    require(isinstance(mixed, dict) and set(mixed) == {"real", "imag"},
            f"{tag} mixed mass must be a complex real/imag object")
    require(isinstance(spectrum, dict), f"{tag} scalar spectrum missing")
    result = {
        "model": model,
        "stage": stage,
        "mu_gev": finite(raw.get("mu_gev"), f"{tag}.mu_gev"),
        "mH2": finite(raw.get("mHSq_gev2"), f"{tag}.mHSq_gev2"),
        "mS12": finite(raw.get("mS1Sq_gev2"), f"{tag}.mS1Sq_gev2"),
        "mS22": finite(raw.get("mS2Sq_gev2"), f"{tag}.mS2Sq_gev2"),
        "re_mu": finite(mixed.get("real"), f"{tag}.mu.real"),
        "im_mu": finite(mixed.get("imag"), f"{tag}.mu.imag"),
        "m_minus2": finite(spectrum.get("m_minus_squared_gev2"),
                           f"{tag}.m_minus_squared_gev2"),
        "m_plus2": finite(spectrum.get("m_plus_squared_gev2"),
                          f"{tag}.m_plus_squared_gev2"),
    }
    require(result["mu_gev"] > 0 and result["m_minus2"] <= result["m_plus2"],
            f"{tag} invalid scale/eigenvalue order")
    require(spectrum.get("real_degeneracy_per_eigenvalue") == 4,
            f"{tag} unexpected real eigenvalue degeneracy")
    require(type(spectrum.get("tachyonic")) is bool,
            f"{tag} tachyonic flag missing")
    result["tachyonic"] = spectrum["tachyonic"]
    require(result["tachyonic"] == (result["m_minus2"] < 0),
            f"{tag} tachyonic flag inconsistent with saved eigenvalue")
    eigen = spectrum.get("mass_squared_eigenvalues_gev2")
    require(isinstance(eigen, list) and len(eigen) == 8,
            f"{tag} expected eight real scalar Hessian eigenvalues")
    saved = [finite(v, f"{tag}.eigen[{i}]") for i, v in enumerate(eigen)]
    require(all(close(v, result["m_minus2"]) for v in saved[:4]) and
            all(close(v, result["m_plus2"]) for v in saved[4:]),
            f"{tag} saved eight-fold spectrum differs from stored m_+/- values")
    return result


def _same_mass_state(a: Mapping[str, Any], b: Mapping[str, Any], reason: str) -> None:
    require(close(a["mu_gev"], b["mu_gev"]) and
            all(close(a[k], b[k]) for k in KEYS), reason)


def read_model(path: Path, model: str) -> dict[str, Any]:
    """Validate native report provenance, trajectory structure and F-boundary continuity."""
    require(path.is_file(), f"missing source result: {path}")
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(payload, dict), f"invalid result object: {path}")
    expected = {
        "format": SOURCE_FORMAT,
        "status": SOURCE_STATUS,
        "model": model,
        "solver_success": True,
        "one_loop_convention": CONVENTION,
        "loop_factor_applied_exactly_once": True,
        "previous_mass_betas_replaced_in_both_stages": True,
        "physical_matchete_rgbeta_mapping_certified": False,
        "scalar_mass_eigenstate_threshold_matching_performed": False,
        "physical_neutrino_mass_prediction": False,
        "five_real_mass_coordinates": list(MASS_COORDS),
    }
    for k, value in expected.items():
        require(payload.get(k) == value, f"{model} source metadata mismatch: {k}")
    uv_range = payload.get("uv_range_gev")
    eft_range = payload.get("eft1_range_gev")
    require(uv_range == [1e7, 1e5] and eft_range == [1e5, 1e3],
            f"{model} scales differ from verified UV/F/EFT1 benchmark")
    threshold = payload.get("threshold")
    require(isinstance(threshold, dict) and
            threshold.get("scale_gev") == 1e5 and
            threshold.get("tree_level_renormalisable_only") is True and
            threshold.get("finite_mass_matching_performed") is False and
            threshold.get("wilson_coefficient_matching_performed") is False,
            f"{model} unsupported threshold matching assumptions")
    for stage, keys in (("UV", payload.get("uv_beta_keys")),
                        ("EFT1", payload.get("eft1_beta_keys"))):
        require(isinstance(keys, list) and len(set(keys)) == len(keys) and
                set(MASS_SYMBOLS).issubset(set(keys)),
                f"{model}/{stage} scalar beta inventory is incomplete")
        if stage == "UV":
            require({"MF", "y1", "y2"}.issubset(set(keys)),
                    f"{model} UV fermion beta inventory missing")
        else:
            require(not {"MF", "y1", "y2"}.intersection(keys),
                    f"{model} EFT1 retains removed fermions")
    result: dict[str, Any] = {"payload": payload, "source": str(path)}
    for stage, field, limits in (("UV", "uv_trajectory", uv_range),
                                 ("EFT1", "eft1_trajectory", eft_range)):
        rows = payload.get(field)
        require(isinstance(rows, list) and 2 <= len(rows) <= 200,
                f"{model}/{stage} trajectory missing/incomplete")
        decoded = [snapshot(row, model=model, stage=stage, index=i)
                   for i, row in enumerate(rows)]
        require(close(decoded[0]["mu_gev"], limits[0]) and
                close(decoded[-1]["mu_gev"], limits[-1]),
                f"{model}/{stage} trajectory scale endpoints differ")
        require(all(decoded[i]["mu_gev"] > decoded[i+1]["mu_gev"]
                    for i in range(len(decoded)-1)),
                f"{model}/{stage} scales are not strictly decreasing")
        result[stage.lower()] = decoded
    endpoints = (("uv_initial", "uv", 0), ("uv_threshold", "uv", -1),
                 ("eft1_threshold", "eft1", 0), ("eft1_final", "eft1", -1))
    for label, segment, index in endpoints:
        saved = snapshot(payload.get(label), model=model, stage=segment,
                         index=index)
        _same_mass_state(saved, result[segment][index],
                         f"{model} stored {label} differs from trajectory")
    _same_mass_state(result["uv"][-1], result["eft1"][0],
                     f"{model} mass state discontinuity at fermion threshold")
    for name in ("uv_solver_nfev", "eft1_solver_nfev"):
        val = payload.get(name)
        require(type(val) is int and val > 0, f"{model} invalid {name}")
    return result


def _display_num(x: float) -> str:
    return f"{x:,.6g}"


def _write_csv(path: Path, fieldnames: list[str], records: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def _summary_markdown(reports: Mapping[str, dict[str, Any]], comparison: dict,
                      *, plots: bool) -> str:
    b, c = reports["B_-1"], reports["C_-1"]
    rows = [
        "# Workstream B: coupled UV → EFT1 mass benchmark",
        "",
        "Native-RGBeta candidate benchmark; read-only display of previously ",
        "integrated JSON trajectories. All listed mass quantities are in GeV².",
        "No RGEs were recomputed and no physical matching is certified.",
        "",
        "| Quantity (GeV²) | B: 10⁷ GeV | B: 10⁵ GeV | B: 10³ GeV | C: 10⁷ GeV | C: 10⁵ GeV | C: 10³ GeV |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for key, label, _ in QUANTITIES:
        values = [b["uv"][0][key], b["uv"][-1][key], b["eft1"][-1][key],
                  c["uv"][0][key], c["uv"][-1][key], c["eft1"][-1][key]]
        rows.append("| " + label + " | " + " | ".join(map(_display_num, values)) + " |")
    rows += [
        "",
        "## Endpoint comparison at Q = 10³ GeV",
        "",
        "| Quantity (GeV²) | C(-1) − B(-1) |",
        "|---|---:|",
    ]
    for key, label, _ in QUANTITIES:
        rows.append(f"| {label} | {_display_num(comparison[key])} |")
    rows += ["", "## Numerical diagnostics", ""]
    for model in MODELS:
        report = reports[model]
        p = report["payload"]
        all_rows = report["uv"] + report["eft1"][1:]
        tach = sum(row["tachyonic"] for row in all_rows)
        heavy = p["threshold"].get("heavy_fermion_singular_masses_gev")
        require(isinstance(heavy, list) and len(heavy) == 3,
                f"{model} heavy-fermion threshold spectrum missing")
        heavy_str = ", ".join(_display_num(finite(x, "heavy fermion mass")) for x in heavy)
        rows.append(f"- **{model}**: UV/EFT1 nfev = {p['uv_solver_nfev']}/{p['eft1_solver_nfev']}; "
                    f"{len(all_rows)} distinct saved scales; tachyonic sampled points = "
                    f"{tach}; heavy-fermion singular masses at 10⁵ GeV "
                    f"= [{heavy_str}] GeV.")
    rows += [
        "",
        "## Interpretation and limitations",
        "",
        "The extra scalar squared-mass eigenvalues are unbroken-electroweak-symmetry "
        "quadratic-potential eigenvalues, not physical pole masses. Positivity at "
        "sampled scales does not prove the entire trajectory is tachyon-free.",
        "The threshold is an illustrative tree-level projection at 10⁵ GeV; "
        "it does not include finite fermion-threshold corrections or Wilson matching. "
        "The Matchete/RGBeta physical field and coupling map is uncertified, and "
        "the results are not neutrino-mass predictions.",
        "",
    ]
    if plots:
        rows += [
            "Seven per-quantity PNG figures display B(-1) and C(-1) versus the "
            "renormalization scale, with the fermion threshold marked at 10⁵ GeV.",
        ]
    rows.append("")
    return "\n".join(rows)


def make_plots(reports: Mapping[str, dict[str, Any]], folder: Path) -> list[str]:
    """Optional plots from saved arrays only; pyplot is imported only here."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = []
    for key, label, slug in QUANTITIES:
        fig, ax = plt.subplots(figsize=(7.3, 4.1))
        for model in MODELS:
            report = reports[model]
            combined = report["uv"] + report["eft1"][1:]
            ax.plot([r["mu_gev"] for r in combined],
                    [r[key] for r in combined],
                    label="B(-1)" if model == "B_-1" else "C(-1)",
                    linewidth=1.8)
        ax.axvline(1e5, linestyle="--", alpha=0.6, linewidth=1,
                   label="F projection (10⁵ GeV)")
        ax.set(xscale="log", xlim=(1e7, 1e3),
               xlabel="Renormalization scale Q (GeV)",
               ylabel=f"{label} (GeV²)")
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=9)
        ax.set_title(f"B(-1) / C(-1): {label}")
        fig.text(0.5, 0.005, "Native RGBeta benchmark · physical matching uncertified",
                 ha="center", fontsize=8)
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        filename = f"bc_{slug}_uv_eft1.png"
        fig.savefig(folder / filename, dpi=155, bbox_inches="tight")
        plt.close(fig)
        files.append(filename)
    return files


def build(*, input_dir: Path = INPUT_DIR, output_dir: Path = OUTPUT_DIR,
          plots: bool = True) -> dict[str, Any]:
    reports: dict[str, dict[str, Any]] = {}
    hashes = {}
    for model in MODELS:
        tag = model.replace("-", "m").replace("+", "p")
        path = input_dir / f"{tag}_uv_to_eft1_coupled_five_mass_v1.json"
        reports[model] = read_model(path, model)
        hashes[model] = sha256(path.read_bytes()).hexdigest()
    for key in ("initial_higgs_mass_squared_gev2_at_uv", "initial_mu2_fraction_at_uv",
                "initial_mu2_phase_at_uv", "rtol", "atol"):
        b = finite(reports["B_-1"]["payload"].get(key), f"B_-1.{key}")
        c = finite(reports["C_-1"]["payload"].get(key), f"C_-1.{key}")
        require(close(b, c), f"B/C benchmarks differ in input setting {key}")
    output_dir.mkdir(parents=True, exist_ok=True)
    records = [dict(row) for model in MODELS for segment in
               (reports[model]["uv"], reports[model]["eft1"]) for row in segment]
    trajectory_csv = "bc_uv_eft1_mass_trajectories.csv"
    _write_csv(output_dir / trajectory_csv,
               ["model", "stage", "mu_gev", *KEYS, "tachyonic"], records)
    snapshots = []
    for model in MODELS:
        for name, row in (("uv_initial", reports[model]["uv"][0]),
                          ("uv_threshold", reports[model]["uv"][-1]),
                          ("eft1_final", reports[model]["eft1"][-1])):
            snapshots.append({"model": model, "point": name, **row})
    endpoint_csv = "bc_uv_eft1_mass_endpoints.csv"
    _write_csv(output_dir / endpoint_csv,
               ["model", "point", "stage", "mu_gev", *KEYS, "tachyonic"], snapshots)
    delta = {key: reports["C_-1"]["eft1"][-1][key] -
                  reports["B_-1"]["eft1"][-1][key] for key in KEYS}
    figure_files = make_plots(reports, output_dir) if plots else []
    markdown = "bc_uv_eft1_benchmark_summary.md"
    (output_dir / markdown).write_text(
        _summary_markdown(reports, delta, plots=plots), encoding="utf-8")
    summary: dict[str, Any] = {
        "format": OUTPUT_FORMAT,
        "status": "PresentationOnlyNoPhysicsRecomputed",
        "source_sha256": hashes,
        "scale_points_gev": [1e7, 1e5, 1e3],
        "mass_units": "GeV^2",
        "initial_benchmark": {
            key: reports["B_-1"]["payload"][key] for key in
            ("initial_higgs_mass_squared_gev2_at_uv",
             "initial_mu2_fraction_at_uv", "initial_mu2_phase_at_uv")},
        "model_summary": {
            model: {
                "uv_initial": {k: reports[model]["uv"][0][k] for k in KEYS},
                "uv_threshold": {k: reports[model]["uv"][-1][k] for k in KEYS},
                "eft1_final": {k: reports[model]["eft1"][-1][k] for k in KEYS},
                "uv_solver_nfev": reports[model]["payload"]["uv_solver_nfev"],
                "eft1_solver_nfev": reports[model]["payload"]["eft1_solver_nfev"],
                "sampled_tachyonic_count": sum(r["tachyonic"] for r in
                    reports[model]["uv"] + reports[model]["eft1"][1:]),
            } for model in MODELS},
        "eft1_final_C_minus_B": delta,
        "outputs": [markdown, trajectory_csv, endpoint_csv, *figure_files,
                    "bc_uv_eft1_benchmark_summary.json"],
        "physics_reevaluated": False,
        "physical_matching_certified": False,
        "neutrino_mass_predicted": False,
    }
    summary_json = "bc_uv_eft1_benchmark_summary.json"
    (output_dir / summary_json).write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--no-plots", action="store_true",
                        help="Write CSV/JSON/Markdown without requiring matplotlib")
    args = parser.parse_args()
    result = build(input_dir=args.input_dir, output_dir=args.output_dir,
                   plots=not args.no_plots)
    folder = args.output_dir
    print(f"B/C five-mass display: {len(result['outputs'])} files -> {folder}")
    for model in MODELS:
        end = result["model_summary"][model]["eft1_final"]
        print(f"{model}: mH^2={end['mH2']:.6g} GeV^2, "
              f"m_-^2={end['m_minus2']:.6g} GeV^2")
    print("Presentation only; no RGE integration or physical matching performed.")


if __name__ == "__main__":
    main()
