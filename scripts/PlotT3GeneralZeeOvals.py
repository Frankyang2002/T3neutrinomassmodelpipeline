"""Plot the SAVED merged-fit general-T3 Zee sensitivity scans against data.

Read-only post-processing: no RGEs, fitting, Wilson matching, or GitHub writes.
Source: output/studies/general_scalar_extended/zee_style_general_sensitivity/
        zee_general_all16_merged_v1/{points.csv,summary.json}

Plot both NuFIT 6.1 (NO, IC24) Gaussianised marginal inputs and JUNO's
207.2-live-day solar measurement (2026-09-29). Experimental 2D ovals use a
DIAGONAL Gaussian covariance approximation, not official joint contours.
Each model's sensitivity oval is 2 times the empirical centred sample-
standard-deviation covariance of all off-centre T1/T2/T3 gamma points. It is
NOT a statistical confidence region.

PowerShell:
  python -X utf8 -m scripts.PlotT3GeneralZeeOvals
  python -X utf8 -m scripts.PlotT3GeneralZeeOvals --x sin2_theta13 --y sin2_theta23
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (ROOT / "output/studies/general_scalar_extended/"
                  "zee_style_general_sensitivity/zee_general_all16_merged_v1")
DEFAULT_TARGET = ROOT / "configs/nufit_6_1_ic24_no.json"
DEFAULT_OUT = ROOT / "Reports/output/thesis_general_5_3/zee_general_merged_ovals_v1"
MODEL_ORDER = tuple(f"{family}_{a:+d}" for family, aa in {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2)
}.items() for a in aa)
DIRECTIONS = ("T1", "T2", "T3")
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12",
       "sin2_theta13", "sin2_theta23")
SOLAR = frozenset(("delta_m21_sq_ev2", "sin2_theta12"))
LABELS = {
    "delta_m21_sq_ev2": r"$\Delta m^2_{21}\ [10^{-5}\ \mathrm{eV}^2]$",
    "delta_m3l_sq_ev2": r"$\Delta m^2_{3\ell}\ [10^{-3}\ \mathrm{eV}^2]$",
    "sin2_theta12": r"$\sin^2\theta_{12}$",
    "sin2_theta13": r"$\sin^2\theta_{13}$",
    "sin2_theta23": r"$\sin^2\theta_{23}$",
}
UNIT_SCALE = {"delta_m21_sq_ev2": 1e5, "delta_m3l_sq_ev2": 1e3,
              "sin2_theta12": 1., "sin2_theta13": 1., "sin2_theta23": 1.}
# JUNO Collaboration, arXiv:2609.37895, first posted 2026-09-29,
# 207.2 live days (8294 inverse beta-decay candidates); normal ordering.
JUNO_2026 = {
    "central_values": {"sin2_theta12": 0.3036, "delta_m21_sq_ev2": 7.388e-5},
    "one_sigma_errors": {"sin2_theta12": 0.0064, "delta_m21_sq_ev2": 0.078e-5},
    "source": "JUNO Collaboration, 207.2 live days, arXiv:2609.37895 (2026-09-29)",
    "url": "https://arxiv.org/abs/2609.37895",
}
# 68.27% and 99.73% coverage for a *two-dimensional* Gaussian.
JOINT_RADII = {"68% joint": math.sqrt(-2 * math.log(1 - 0.6827)),
               "99.73% joint": math.sqrt(-2 * math.log(1 - 0.9973))}


def _read_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_nufit(path: Path) -> dict:
    data = _read_json(path)
    if data.get("ordering") != "NO" or "NuFIT 6.1" not in data.get("dataset", ""):
        raise ValueError(f"Expected NuFIT 6.1 NO input, not {path}")
    for name in OBS:
        mu = float(data["central_values"][name])
        sigma = float(data["one_sigma_errors"][name])
        if not (math.isfinite(mu) and math.isfinite(sigma) and sigma > 0):
            raise ValueError(f"Invalid NuFIT {name}: {mu}, {sigma}")
    return data


def load_points(path: Path, summary_path: Path, *, models: list[str] | None = None
                ) -> tuple[dict[str, dict[str, list[dict]]], dict[str, dict]]:
    """Validate every gamma direction and reject partial/foreign input files."""
    if not path.is_file() or not summary_path.is_file():
        raise FileNotFoundError("Need both saved points.csv and summary.json in scan folder")
    summary = _read_json(summary_path)
    if not summary.get("schema", "").startswith("t3-general-zee-compensated-sensitivity-"):
        raise ValueError("Unrecognised general-Zee source summary schema")
    meta = {r["model"]: r for r in summary.get("models", [])}
    selected = models or list(MODEL_ORDER)
    if len(selected) != len(set(selected)) or not selected or set(selected) - set(MODEL_ORDER):
        raise ValueError("Invalid/duplicate model selection")
    raw = defaultdict(lambda: defaultdict(list))
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            model = row.get("model")
            if model not in selected:
                continue
            if row.get("status") != "Success" or row.get("direction") not in DIRECTIONS:
                raise ValueError(f"Non-success/unrecognised scan point: {model} {row}")
            parsed = {"model": model, "direction": row["direction"],
                      "gamma": float(row["gamma"]), "chi2": float(row["chi2"])}
            parsed.update({name: float(row[name]) for name in OBS})
            if not all(math.isfinite(v) for v in
                       [parsed["gamma"], parsed["chi2"], *(parsed[n] for n in OBS)]):
                raise ValueError(f"Nonfinite saved values: {model}")
            if parsed["gamma"] <= 0:
                raise ValueError("Nonpositive gamma")
            raw[model][row["direction"]].append(parsed)
    result = {}
    for model in selected:
        if model not in meta or meta[model].get("status") != "Complete" or model not in raw:
            raise ValueError(f"Saved scan incomplete/missing: {model}")
        if meta[model].get("failed_points") != 0:
            raise ValueError(f"Saved scan reports failed points: {model}")
        sets, centres = [], []
        for direction in DIRECTIONS:
            points = sorted(raw[model][direction], key=lambda r: r["gamma"])
            g = [r["gamma"] for r in points]
            if len(g) < 3 or len(g) != len(set(g)) or sum(abs(x-1) < 1e-12 for x in g) != 1:
                raise ValueError(f"Bad gamma grid: {model} {direction}")
            sets.append(g)
            centres.append(next(p for p in points if abs(p["gamma"]-1) < 1e-12))
            raw[model][direction] = points
        if not all(np.allclose(sets[0], g, atol=1e-12, rtol=0) for g in sets[1:]):
            raise ValueError(f"Inconsistent gamma grids: {model}")
        for name in OBS:
            a = centres[0][name]
            if any(abs(a-c[name]) > 1e-6 * max(abs(a), 1e-12) for c in centres[1:]):
                raise ValueError(f"Gamma=1 centres differ: {model} {name}")
        result[model] = dict(raw[model])
    return result, meta


def centre(points: dict[str, list[dict]]) -> dict:
    return next(p for p in points["T1"] if abs(p["gamma"] - 1) < 1e-12)


def centered_covariance(points: dict[str, list[dict]], x: str, y: str) -> np.ndarray:
    """Descriptive covariance around the fitted centre, NOT statistical uncertainty."""
    c = centre(points)
    sample = np.array([[p[x] - c[x], p[y] - c[y]]
                       for direction in DIRECTIONS for p in points[direction]
                       if abs(p["gamma"] - 1) > 1e-12], dtype=float)
    if len(sample) < 3:
        raise ValueError("Need at least three off-centre perturbations")
    return (sample.T @ sample) / (len(sample) - 1)


def central_slope(points: list[dict], observable: str) -> tuple[float, float]:
    """dO/dln gamma at 1, using the nearest two points on opposite sides."""
    down = sorted((p for p in points if p["gamma"] < 1),
                  key=lambda p: 1 - p["gamma"])
    up = sorted((p for p in points if p["gamma"] > 1),
                key=lambda p: p["gamma"] - 1)
    if not down or not up:
        raise ValueError("Gamma range must straddle 1")
    left, right = down[0], up[0]
    h = math.log(right["gamma"]) - math.log(left["gamma"])
    return (right[observable] - left[observable]) / h, right["gamma"] - left["gamma"]


def diag_chi2(point: dict, experiment: dict, names: tuple[str, ...]) -> float:
    return sum(((point[name] - experiment["central_values"][name]) /
                 experiment["one_sigma_errors"][name]) ** 2 for name in names)


def _cov_ellipse(ax, center_xy, cov_xy, *, nstd, **kwargs):
    from matplotlib.patches import Ellipse
    eig, vec = np.linalg.eigh(np.asarray(cov_xy, dtype=float))
    eig = np.maximum(eig, 0)
    if eig[-1] <= 0:
        return False
    major = vec[:, -1]
    width, height = 2*nstd*np.sqrt(eig[-1]), 2*nstd*np.sqrt(eig[0])
    angle = math.degrees(math.atan2(major[1], major[0]))
    ax.add_patch(Ellipse(center_xy, width=width, height=height, angle=angle, **kwargs))
    return True


def _exp_ellipse(ax, exp: dict, x: str, y: str, nstd: float, **kw) -> bool:
    if x not in exp["central_values"] or y not in exp["central_values"]:
        return False
    cx = exp["central_values"][x] * UNIT_SCALE[x]
    cy = exp["central_values"][y] * UNIT_SCALE[y]
    sx = exp["one_sigma_errors"][x] * UNIT_SCALE[x]
    sy = exp["one_sigma_errors"][y] * UNIT_SCALE[y]
    return _cov_ellipse(ax, (cx, cy), np.diag([sx*sx, sy*sy]), nstd=nstd, **kw)


def _plot_panel(ax, model: str, points: dict, nufit: dict, x: str, y: str,
                *, joint_3sigma=True, show_legend=False, title=True):
    colors = {"T1": "#2274bc", "T2": "#259867", "T3": "#9a55ac"}
    mu = centre(points)
    # Only the two solar variables are actually measured by this JUNO result.
    solar = {x, y} == SOLAR
    for tag, exp, color in (("NuFIT 6.1", nufit, "#657899"),
                            ("JUNO 2026", JUNO_2026, "#d47b22")):
        if tag.startswith("JUNO") and not solar:
            continue
        if joint_3sigma:
            _exp_ellipse(ax, exp, x, y, JOINT_RADII["99.73% joint"],
                         edgecolor=color, linewidth=1.0, linestyle=":",
                         facecolor="none", zorder=1)
        _exp_ellipse(ax, exp, x, y, JOINT_RADII["68% joint"],
                     edgecolor=color, linewidth=1.35, facecolor="none",
                     label=tag + " 68% joint (diagonal Gaussian)", zorder=2)
        ax.plot(exp["central_values"][x]*UNIT_SCALE[x],
                exp["central_values"][y]*UNIT_SCALE[y], marker="+",
                markersize=6, color=color, zorder=4)
    cxy = (mu[x]*UNIT_SCALE[x], mu[y]*UNIT_SCALE[y])
    cov = centered_covariance(points, x, y)
    scales = np.diag([UNIT_SCALE[x], UNIT_SCALE[y]])
    cov = scales @ cov @ scales
    _cov_ellipse(ax, cxy, cov, nstd=2., edgecolor="#333333",
                 facecolor="none", linestyle="--", linewidth=1.1,
                 label="Model 2× empirical sensitivity oval", zorder=2)
    for d in DIRECTIONS:
        ps = points[d]
        ax.plot([p[x]*UNIT_SCALE[x] for p in ps],
                [p[y]*UNIT_SCALE[y] for p in ps], '-o', markersize=2.6,
                linewidth=1.3, label=d, color=colors[d], zorder=3)
    ax.scatter([cxy[0]], [cxy[1]], marker="*", s=65, color="black", zorder=5)
    if title:
        label = model.replace("_", "(") + ")"
        ax.set_title(label, fontsize=10)
    ax.set_xlabel(LABELS[x], fontsize=9)
    ax.set_ylabel(LABELS[y], fontsize=9)
    ax.grid(alpha=.18)
    if show_legend:
        ax.legend(fontsize=7, loc="best")


def render(plots: dict, nufit: dict, out: Path, x: str, y: str,
           *, per_model=True) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 4, figsize=(17.6, 17.4), constrained_layout=True)
    for ax, model in zip(axes.ravel(), plots):
        _plot_panel(ax, model, plots[model], nufit, x, y)
        ax.tick_params(labelsize=8)
    for ax in list(axes.ravel())[len(plots):]:
        ax.axis("off")
    fig.suptitle("General T3 merged-fit Zee T1/T2/T3 sensitivity vs NuFIT 6.1 and JUNO 2026\n"
                 "Experimental contours: diagonal-Gaussian approximation; model oval: descriptive, not a confidence region",
                 fontsize=13)
    fig.savefig(out / f"model_grid_{x}_vs_{y}.png", dpi=170)
    plt.close(fig)
    if not per_model:
        return
    individual = out / "models"
    individual.mkdir(exist_ok=True)
    for model, data in plots.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        _plot_panel(ax, model, data, nufit, x, y, show_legend=True)
        fig.suptitle("General-Lagrangian conditional UV/EFT1 RGE + one-scale C5", fontsize=10)
        fig.tight_layout()
        fig.savefig(individual / f"{model.replace('+','p').replace('-','m')}_{x}_vs_{y}.png",
                    dpi=190)
        plt.close(fig)


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{k: row.get(k, "") for k in fields} for row in rows])


def analyse(points: dict[str, dict], meta: dict[str, dict],
            nufit: dict, out: Path, x: str, y: str) -> None:
    slopes = []
    baseline = []
    for model, dirs in points.items():
        c = centre(dirs)
        baseline.append({
            "model": model,
            "saved_fit_chi2_5obs": meta[model].get("baseline_chi2_saved"),
            "nufit_solar_diag_chi2": diag_chi2(c, nufit,
                        ("sin2_theta12", "delta_m21_sq_ev2")),
            "juno_solar_diag_chi2": diag_chi2(c, JUNO_2026,
                        ("sin2_theta12", "delta_m21_sq_ev2")),
            **{name: c[name] for name in OBS},
        })
        for direction in DIRECTIONS:
            for name in OBS:
                deriv, step = central_slope(dirs[direction], name)
                effect = 0.01 * deriv
                row = {"model": model, "direction": direction, "observable": name,
                       "gamma_pair_width": step, "baseline": c[name],
                       "dO_dln_gamma": deriv,
                       "estimated_deltaO_for_1pct_gamma": effect,
                       "abs_effect_over_nufit_sigma": abs(effect) /
                        nufit["one_sigma_errors"][name],
                       "abs_effect_over_juno_sigma":
                           abs(effect) / JUNO_2026["one_sigma_errors"][name]
                           if name in JUNO_2026["one_sigma_errors"] else ""}
                slopes.append(row)
    slopes.sort(key=lambda row: (MODEL_ORDER.index(row["model"]),
                 DIRECTIONS.index(row["direction"]), OBS.index(row["observable"])))
    write_csv(out / "local_sensitivity_slopes.csv", slopes, (
        "model", "direction", "observable", "gamma_pair_width", "baseline",
        "dO_dln_gamma", "estimated_deltaO_for_1pct_gamma",
        "abs_effect_over_nufit_sigma", "abs_effect_over_juno_sigma"))
    write_csv(out / "baseline_nufit_juno_comparison.csv", baseline, (
        "model", "saved_fit_chi2_5obs", "nufit_solar_diag_chi2",
        "juno_solar_diag_chi2", *OBS))
    ranking = sorted((r for r in slopes if r["observable"] in SOLAR),
                     key=lambda row: float(row["abs_effect_over_juno_sigma"]),
                     reverse=True)
    write_csv(out / "solar_sensitivity_ranking.csv", ranking, (
        "model", "direction", "observable", "estimated_deltaO_for_1pct_gamma",
        "abs_effect_over_juno_sigma", "abs_effect_over_nufit_sigma"))


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=DEFAULT_SOURCE,
                   help="Saved general Zee campaign folder (points.csv, summary.json)")
    p.add_argument("--nufit", type=Path, default=DEFAULT_TARGET)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    p.add_argument("--models", default="all", help="Comma-separated model IDs or all")
    p.add_argument("--x", choices=OBS, default="sin2_theta12")
    p.add_argument("--y", choices=OBS, default="delta_m21_sq_ev2")
    p.add_argument("--no-per-model", action="store_true",
                   help="Only combined 4x4 figure (plus CSV diagnostics)")
    args = p.parse_args(argv)
    if args.x == args.y:
        p.error("--x and --y must differ")
    selected = (list(MODEL_ORDER) if args.models == "all" else
                [name.strip() for name in args.models.split(",") if name.strip()])
    rows, meta = load_points(args.source / "points.csv", args.source / "summary.json",
                             models=selected)
    nu = load_nufit(args.nufit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    render(rows, nu, args.output_dir, args.x, args.y,
           per_model=not args.no_per_model)
    analyse(rows, meta, nu, args.output_dir, args.x, args.y)
    manifest = {
        "format": "t3-general-zee-merged-fit-ovals-v1",
        "models": list(rows), "directions": list(DIRECTIONS),
        "gamma_grid": [pt["gamma"] for pt in next(iter(rows.values()))["T1"]],
        "input_files": {str(args.source / "points.csv"):_sha256(args.source / "points.csv"),
                        str(args.source / "summary.json"):_sha256(args.source / "summary.json"),
                        str(args.nufit):_sha256(args.nufit)},
        "nufit": {"source": nu.get("dataset"),
                  "central_values": nu["central_values"],
                  "one_sigma_errors": nu["one_sigma_errors"]},
        "juno_measured_2026": JUNO_2026,
        "experimental_oval_definition":
            "diagonal-covariance Gaussian 68.27% and 99.73% joint 2D; "
            "no published experimental correlation assumed; NOT official confidence contours",
        "model_oval_definition":
            "twice the sample standard-deviation covariance about gamma=1 "
            "of all 30 off-centre gamma points/model; descriptive sensitivity, NOT confidence",
        "physics_scope": "saved conditional general-Lagrangian candidate UV/EFT1 + "
                         "one-scale C5, NOT certified all-fermions-first "
                         "then all-scalars sequential EFT matching",
        "plots_use_saved_points_only": True,
        "juno_solar_chi2_definition": "two marginal Gaussian residuals, diagonal; "
                                       "NOT official JUNO likelihood or combined global chi2",
        "local_slope_definition": "central difference dO/dln(gamma) using nearest gamma pair "
                                  "straddling 1; 1% change approximated linearly",
    }
    (args.output_dir / "provenance.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Plotted {len(rows)}/{len(selected)} complete models: {args.output_dir}")
    print("JUNO measured 2026, NOT the older six-year projection")
    print("No RGEs or matching rerun; no saved scan files modified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
