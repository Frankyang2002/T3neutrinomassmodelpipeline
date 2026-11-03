"""Plot SAVED general-T3 merged-fit sensitivity in Zee Figure-9 style.

Read-only post-processing: no RGEs, fitting, Wilson matching, or GitHub writes.
Source: output/studies/general_scalar_extended/zee_style_general_sensitivity/
        zee_general_all16_merged_v1/{points.csv,summary.json}

Plot both NuFIT 6.1 (NO, IC24) Gaussianised marginal inputs and JUNO's
207.2-live-day solar measurement (2026-09-29). Experimental 2D ovals use a
DIAGONAL Gaussian covariance approximation, not official joint contours.
Zee-style figures use filled 1D marginal experimental bands (not joint
confidence regions), a black fitted point, coloured T1/T2/T3 trajectories,
and optional descriptive model sensitivity oval (--model-oval).
The experimental bands are based on diagonal Gaussian marginal errors.

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
DEFAULT_OUT = ROOT / "Reports/output/thesis_general_5_3/zee_style_all16_v2"
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


# Same visual grammar as Fig. 9 of the Zee-model paper: a blue current-data
# band, a distinct JUNO band and three coloured gamma trajectories.
PALETTE = {
    "nufit_3": "#D6EAF5", "nufit_1": "#59AFE2",
    "juno_1": "#BAD78C", "T1": "#2584BC", "T2": "#D97832",
    "T3": "#478A32", "model_oval": "#6D626A",
}


def _plot_panel(ax, model: str, points: dict, nufit: dict, x: str, y: str,
                *, show_legend=False, title=True, show_model_oval=False):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.set_facecolor("white")
    c = centre(points)
    cx, cy = c[x]*UNIT_SCALE[x], c[y]*UNIT_SCALE[y]
    solar = {x, y} == SOLAR
    # Experimental bands are independent 1D Gaussian marginal intervals.
    # Do NOT label these as 2D joint confidence contours.
    regions = [
        (nufit, 3., PALETTE["nufit_3"], 0.65),
        (nufit, 1., PALETTE["nufit_1"], 0.38),
    ]
    if solar:
        regions.append((JUNO_2026, 1., PALETTE["juno_1"], 0.60))
    for exp, n_sigma, color, alpha in regions:
        for observable, unit, painter in (
            (x, UNIT_SCALE[x], ax.axvspan),
            (y, UNIT_SCALE[y], ax.axhspan),
        ):
            if observable not in exp["central_values"]:
                continue
            mu = exp["central_values"][observable] * unit
            sigma = exp["one_sigma_errors"][observable] * unit
            painter(mu - n_sigma*sigma, mu + n_sigma*sigma,
                    color=color, alpha=alpha, ec="none", zorder=0)

    handles = [
        Patch(facecolor=PALETTE["nufit_3"], edgecolor="none", label="NuFIT 6.1 ±3σ (marginal)"),
        Patch(facecolor=PALETTE["nufit_1"], edgecolor="none", label="NuFIT 6.1 ±1σ (marginal)"),
    ]
    if solar:
        handles.append(Patch(facecolor=PALETTE["juno_1"], edgecolor="none",
                             label="JUNO 2026 ±1σ (marginal)"))
    if show_model_oval:
        cov = centered_covariance(points, x, y)
        scales = np.diag([UNIT_SCALE[x], UNIT_SCALE[y]])
        _cov_ellipse(ax, (cx, cy), scales@cov@scales, nstd=2.,
                     edgecolor=PALETTE["model_oval"], facecolor="none",
                     linestyle=":", linewidth=1.25, zorder=2)
        handles.append(Line2D([], [], color=PALETTE["model_oval"], ls=":",
                              label="Descriptive model response oval"))
    for direction in DIRECTIONS:
        ps = points[direction]
        vx = np.asarray([p[x]*UNIT_SCALE[x] for p in ps], dtype=float)
        vy = np.asarray([p[y]*UNIT_SCALE[y] for p in ps], dtype=float)
        color = PALETTE[direction]
        ax.plot(vx, vy, lw=1.7, color=color, zorder=3,
                marker="o", markersize=3.4, markerfacecolor=color,
                markeredgecolor="white", markeredgewidth=.28)
        ax.scatter([vx[0],vx[-1]], [vy[0],vy[-1]], color=color,
                   s=19, edgecolors="white", linewidths=.35, zorder=4)
        handles.append(Line2D([], [], color=color, lw=1.7, marker="o",
                              markersize=3.5, label=direction))
    ax.plot([cx], [cy], marker="*", markersize=9.2, color="#151515",
            markeredgewidth=.5, markeredgecolor="white", zorder=6)

    # Place the bounds around the experimental 1σ intervals AND the full
    # gamma trajectories so the two most sensitive B/C models are not clipped.
    xvalues = [cx]; yvalues = [cy]
    for direction in DIRECTIONS:
        for point in points[direction]:
            xvalues.append(point[x]*UNIT_SCALE[x]); yvalues.append(point[y]*UNIT_SCALE[y])
    for exp in ([nufit, JUNO_2026] if solar else [nufit]):
        for name, dest in ((x,xvalues),(y,yvalues)):
            if name in exp["central_values"]:
                val = exp["central_values"][name] * UNIT_SCALE[name]
                err = exp["one_sigma_errors"][name] * UNIT_SCALE[name]
                dest.extend([val-err,val+err])
    def limits(values, sigma):
        low, high = min(values), max(values)
        span = max(high-low, 1.6*sigma, 1e-13)
        mid = (high+low)/2
        return (mid-0.60*span, mid+0.60*span)
    ax.set_xlim(*limits(xvalues,nufit["one_sigma_errors"][x]*UNIT_SCALE[x]))
    ax.set_ylim(*limits(yvalues,nufit["one_sigma_errors"][y]*UNIT_SCALE[y]))
    ax.set_xlabel(LABELS[x], fontsize=9)
    ax.set_ylabel(LABELS[y], fontsize=9)
    ax.set_title(model.replace("_", "(")+ ")" if title else "", fontsize=11)
    ax.tick_params(labelsize=8, length=3, direction="in")
    ax.grid(color="#B5BDC4", alpha=.20, linewidth=.5, zorder=.1)
    for edge in ax.spines.values():
        edge.set_color("#686C70"); edge.set_linewidth(.65)
    if show_legend:
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5,1.23),
                  fontsize=8.1, ncol=3, frameon=False, handlelength=1.45,
                  columnspacing=.9)
    return handles


def _save_figure(fig, path: Path, *, dpi: int) -> None:
    """Render in memory: Pillow never opens the Windows destination filename.

    Image.save(path) in earlier versions raised OSError(EINVAL) on some user
    Windows configurations even though plotting was successful. A fresh output
    location and a buffer-to-Path write separate rendering from filesystem I/O.
    """
    from io import BytesIO
    path.parent.mkdir(parents=True,exist_ok=True)
    for fmt in ("png", "pdf"):
        target = path.with_suffix("."+fmt)
        buffer = BytesIO()
        fig.savefig(buffer, format=fmt, dpi=dpi, facecolor="white")
        try:
            target.write_bytes(buffer.getvalue())
        except OSError as exc:
            raise OSError(f"Cannot write {target}: {exc}. Try --output-dir with "
                          "a new writable folder, or --no-grid if only the "
                          "combined figure is affected.") from exc


def render(plots: dict, nufit: dict, out: Path, x: str, y: str,
           *, per_model=True, make_grid=True, model_oval=False) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    # Do individual models FIRST: one large grid file can never prevent the
    # requested 16 thesis figures from being produced.
    if per_model:
        individual = out / "models"
        individual.mkdir(exist_ok=True)
        for model, data in plots.items():
            fig, ax = plt.subplots(figsize=(7.5,5.5))
            _plot_panel(ax, model, data, nufit, x, y,
                        show_legend=True, show_model_oval=model_oval)
            fig.subplots_adjust(top=.77,bottom=.14,left=.14,right=.97)
            fig.text(.52,.98,"General T3: Zee-style compensated UV sensitivity",
                     ha="center",va="top",fontsize=12)
            _save_figure(fig,individual / f"{model.replace('+','p').replace('-','m')}_zee_solar.png",
                         dpi=180)
            plt.close(fig)
            print(f"Saved {model} Zee-style PNG/PDF", flush=True)
    if not make_grid:
        return
    fig, axes = plt.subplots(4,4,figsize=(16.2,15.6))
    handles=None
    for ax, model in zip(axes.flat,plots):
        handles=_plot_panel(ax,model,plots[model],nufit,x,y,
                            show_legend=False,show_model_oval=model_oval)
    for ax in list(axes.flat)[len(plots):]:
        ax.set_axis_off()
    fig.subplots_adjust(left=.065,right=.995,top=.865,bottom=.065,
                        wspace=.37,hspace=.44)
    fig.text(.5,.985,"General T3: Zee-style compensated-coupling sensitivity",
             ha="center",va="top",fontsize=17)
    fig.text(.5,.96,"16 individual best-fit benchmarks • T1/T2/T3 • γ = 0.99–1.01",
             ha="center",va="top",fontsize=11)
    if handles:
        fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.925),
                   ncol=min(6,len(handles)),fontsize=9,frameon=False)
    fig.text(.5,.026,"Shaded regions: one-dimensional Gaussian marginal intervals, not joint contours. "
             "All model predictions use conditional general-Lagrangian one-scale C₅ matching.",
             ha="center",va="bottom",fontsize=9)
    _save_figure(fig,out / "model_grid_zee_solar.png",dpi=170)
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
    p.add_argument("--no-grid", action="store_true",
                   help="Skip the combined grid; write individual Zee-style figures")
    p.add_argument("--model-oval", action="store_true",
                   help="Also show descriptive model-response oval (not confidence region)")
    args = p.parse_args(argv)
    if args.x == args.y:
        p.error("--x and --y must differ")
    selected = (list(MODEL_ORDER) if args.models == "all" else
                [name.strip() for name in args.models.split(",") if name.strip()])
    rows, meta = load_points(args.source / "points.csv", args.source / "summary.json",
                             models=selected)
    nu = load_nufit(args.nufit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.no_per_model and args.no_grid:
        p.error("Cannot disable both individual figures and the combined grid")
    # Write diagnostics first, independent of the image output stage.
    analyse(rows, meta, nu, args.output_dir, args.x, args.y)
    render(rows, nu, args.output_dir, args.x, args.y,
           per_model=not args.no_per_model,
           make_grid=not args.no_grid, model_oval=args.model_oval)
    manifest = {
        "format": "t3-general-zee-style-marginal-bands-v2",
        "models": list(rows), "directions": list(DIRECTIONS),
        "gamma_grid": [pt["gamma"] for pt in next(iter(rows.values()))["T1"]],
        "input_files": {str(args.source / "points.csv"):_sha256(args.source / "points.csv"),
                        str(args.source / "summary.json"):_sha256(args.source / "summary.json"),
                        str(args.nufit):_sha256(args.nufit)},
        "nufit": {"source": nu.get("dataset"),
                  "central_values": nu["central_values"],
                  "one_sigma_errors": nu["one_sigma_errors"]},
        "juno_measured_2026": JUNO_2026,
        "experimental_band_definition":
            "filled 1D Gaussian marginal ±1σ/±3σ for NuFIT 6.1; "
            "filled 1D ±1σ for JUNO 2026 in the solar plane; "
            "NOT joint two-dimensional confidence contours",
        "model_oval_definition":
            "optional: twice the sample standard-deviation covariance about gamma=1 "
            "of all 30 off-centre points/model; descriptive response, NOT confidence",
        "physics_scope": "saved conditional general-Lagrangian candidate UV/EFT1 + "
                         "one-scale C5, NOT certified all-fermions-first "
                         "then all-scalars sequential EFT matching",
        "plots_use_saved_points_only": True,
        "zee_fig9_inspired_bands": True,
        "model_response_oval_drawn": bool(args.model_oval),
        "grid_drawn": not args.no_grid,
        "per_model_figures_drawn": not args.no_per_model,
        "juno_solar_chi2_definition": "two marginal Gaussian residuals, diagonal; "
                                       "NOT official JUNO likelihood or combined global chi2",
        "local_slope_definition": "central difference dO/dln(gamma) using nearest gamma pair "
                                  "straddling 1; 1% change approximated linearly",
    }
    (args.output_dir / "provenance.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved Zee-style figures for {len(rows)}/{len(selected)} complete models: {args.output_dir}")
    print("JUNO measured 2026, NOT the older six-year projection")
    print("No RGEs or matching rerun; no saved scan files modified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
