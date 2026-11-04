"""Plot overview and detailed OVAL figures from saved general-T3 scans.

Sixteen individually fitted models, three compensated T1/T2/T3 variations,
11 gamma points per direction. The overview plots independent-Gaussian
APPROXIMATE joint experimental confidence ellipses (NuFIT 6.1 NO and JUNO
2026 measured solar quantities), plus each model's descriptive response oval.
The separate close-up grid frames the actual gamma trajectories and their
response ovals as deviations O(gamma)-O(1) in labelled powers-of-ten units,
without forcing full-scale experimental regions into view. No arrows or
footnotes are drawn; all gamma points are shown with visible markers.
There is no special star/central marker: gamma=1 is included among the plotted
points. Straight or nearly straight model responses can produce extremely
thin ovals; they are not artificially thickened.

No RGEs, fitting, matching, or other physics calculations are rerun. Existing
C5 results are conditional general one-scale matching, not certified complete
fermion-first sequential EFT matching.

From the repository root:
  python -X utf8 -m scripts.PlotT3GeneralZeeOvals
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
DEFAULT_OUT = ROOT / "Reports/output/thesis_general_5_3/zee_ovals_all16_dualview_v5"
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
# From the repository's 2022 JUNO six-year forecast used in previous plots.
# This is NOT the actual September 2026 measurement. Its illustrative
# central values are explicitly taken from JUNO_2026, not predicted.
JUNO_FUTURE_6Y = {
    "central_values": dict(JUNO_2026["central_values"]),
    "one_sigma_errors": {"sin2_theta12": 0.0016, "delta_m21_sq_ev2": 0.024e-5},
    "source": "JUNO six-year forecast, CPC 46 (2022) 123001, Table 6; centered on JUNO 2026 measured values (illustrative)",
    "url": "https://doi.org/10.1088/1674-1137/ac8bc9",
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


# NuFIT light/dark blue, JUNO orange, and three distinct response colours.
PALETTE = {
    "nufit_3": "#C8E2F4", "nufit_1": "#4D93D0",
    "juno_1": "#F3A36D", "juno_future": "#D57638",
    "T1": "#0077B6", "T2": "#DC653A", "T3": "#289344",
    "model_oval": "#695088",
}


def juno_regions(mode: str) -> list[tuple[dict, str, float, float]]:
    """Return only explicitly labelled measured and/or projected JUNO strips.

    Forecast error bars use the 2022 six-year sensitivity study. Their
    *central value* is the 2026 measured JUNO centre, so these are not a
    measured future likelihood or official JUNO combined fit.
    """
    if mode == "measured":
        return [(JUNO_2026, "JUNO 2026 1σ (measured)", .65, 1.)]
    if mode == "future":
        return [(JUNO_FUTURE_6Y, "JUNO 6-y 1σ (projection)", .72, 1.)]
    if mode == "both":
        return [(JUNO_2026, "JUNO 2026 1σ (measured)", .55, 1.),
                (JUNO_FUTURE_6Y, "JUNO 6-y 1σ (projection)", .66, 1.)]
    raise ValueError(f"Unsupported JUNO mode: {mode}")


def _model_title(model: str) -> str:
    family, sign_alpha = model.split("_")
    return f"{family}({int(sign_alpha):+d})".replace("+0", "0")


# Detail views plot differences from the gamma=1 reference rather than long,
# crowded absolute tick labels (e.g. 0.308799999... or 1e-7 + 3.088e-1).
# The powers of ten are selected independently per model/observable and
# explicitly printed on each axis. The overview retains ABSOLUTE observables.
DELTA_LABELS = {
    "delta_m21_sq_ev2": r"\delta(\Delta m^2_{21})",
    "delta_m3l_sq_ev2": r"\delta(\Delta m^2_{3\ell})",
    "sin2_theta12": r"\delta(\sin^2\theta_{12})",
    "sin2_theta13": r"\delta(\sin^2\theta_{13})",
    "sin2_theta23": r"\delta(\sin^2\theta_{23})",
}


def _detail_frame(points: dict, nufit: dict, x: str, y: str,
                  *, model_oval: bool = True) -> dict:
    """Readable, honest close-up coordinates: (O(gamma)-O(1))/10**power.

    The frame is chosen from the ACTUAL saved trajectory points and, if
    enabled, the two-standard-deviation *descriptive* response oval. No
    covariance or trajectory coordinates are altered: this is an affine
    relabelling of the plotted axes only. Zero-variation directions use a
    numerical-resolution floor so they do not acquire misleading precision.
    """
    c = centre(points)
    cov = centered_covariance(points, x, y) if model_oval else None
    units = {}
    powers = {}
    limits = {}
    for axis, observable in enumerate((x,y)):
        offset = float(c[observable])
        differences = [float(p[observable]) - offset
                       for direction in DIRECTIONS for p in points[direction]]
        radius = (2.0 * math.sqrt(max(0., float(cov[axis, axis])))
                  if cov is not None else 0.0)
        extent = max([radius] + [abs(z) for z in differences])
        # Do not fabricate sub-machine-precision structure for constant data.
        roundoff = max(16 * np.finfo(float).eps * abs(offset), 1e-24)
        meaningful_extent = max(extent, roundoff)
        power = int(math.floor(math.log10(meaningful_extent)))
        unit = 10.0 ** power
        if not (math.isfinite(unit) and unit > 0):
            raise ValueError(f"Invalid close-up power for {observable}")
        plot_values = [v/unit for v in differences]
        if model_oval:
            plot_values.extend((-radius/unit, radius/unit))
        left, right = min(plot_values), max(plot_values)
        width = max(right-left, 1.0)
        midpoint = .5*(left+right)
        margin = .12 * width
        # Keep the gamma=1 origin visible for every model.
        limits[observable] = (min(midpoint-.5*width-margin, -0.08*width),
                              max(midpoint+.5*width+margin, 0.08*width))
        units[observable] = unit
        powers[observable] = power
    return {"origin": {x: float(c[x]), y: float(c[y])},
            "units": units, "powers": powers, "limits": limits}


def _detail_label(observable: str, power: int) -> str:
    suffix = r"\,\mathrm{eV}^2" if observable.startswith("delta_m") else ""
    return rf"${DELTA_LABELS[observable]}\ [10^{{{power}}}{suffix}]$"


def _axis_limits(points: dict, nufit: dict, juno_mode: str, x: str, y: str,
                 *, zoom: str, model_oval: bool = True
                 ) -> tuple[tuple[float, float], tuple[float, float]]:
    """Overview includes experimental ovals; detail frames model points/oval.

    For an ellipse with covariance C and radius r, the projection onto each
    displayed coordinate is r*sqrt(C_ii), regardless of its orientation.
    """
    c = centre(points)
    xs = [float(p[x])*UNIT_SCALE[x] for d in DIRECTIONS for p in points[d]]
    ys = [float(p[y])*UNIT_SCALE[y] for d in DIRECTIONS for p in points[d]]
    if zoom == "experiment":
        experimental = [(nufit, JOINT_RADII["99.73% joint"])]
        if {x, y} == SOLAR:
            experimental += [(e, JOINT_RADII["68% joint"])
                             for e, _, _, _ in juno_regions(juno_mode)]
        for record, radius in experimental:
            for obs, dest in ((x, xs), (y, ys)):
                mu = float(record["central_values"][obs])*UNIT_SCALE[obs]
                sig = float(record["one_sigma_errors"][obs])*UNIT_SCALE[obs]
                dest.extend((mu-radius*sig, mu+radius*sig))
    elif zoom != "trajectory":
        raise ValueError(f"Unknown zoom {zoom!r}")

    if model_oval:
        cov = centered_covariance(points, x, y)
        for obs, dest, idx in ((x, xs, 0), (y, ys, 1)):
            spread = 2.0*math.sqrt(max(float(cov[idx,idx]), 0.0))*UNIT_SCALE[obs]
            origin = c[obs]*UNIT_SCALE[obs]
            dest.extend((origin-spread,origin+spread))

    def lim(vals: list[float], obs: str):
        lo, hi = float(min(vals)), float(max(vals))
        sigma = float(nufit["one_sigma_errors"][obs])*UNIT_SCALE[obs]
        # A small but nonzero floor handles effectively insensitive benchmarks.
        minimum = max(2.5e-5*sigma, abs(float(c[obs])*UNIT_SCALE[obs])*1e-12,
                      1e-15)
        width = max(hi-lo, minimum)
        midpoint = .5*(hi+lo)
        padding = .08 if zoom == "experiment" else .13
        return midpoint-(.5+padding)*width, midpoint+(.5+padding)*width
    return lim(xs,x),lim(ys,y)


def _plot_panel(ax, model: str, points: dict, nufit: dict, x: str, y: str,
                *, show_legend=False, title=True, show_model_oval=True,
                juno_mode="measured", zoom="experiment"):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.set_facecolor("white")
    solar = {x,y} == SOLAR
    detail = (_detail_frame(points,nufit,x,y,model_oval=show_model_oval)
              if zoom == "trajectory" else None)
    if detail is None:
        xlim, ylim = _axis_limits(points,nufit,juno_mode,x,y,zoom=zoom,
                                 model_oval=show_model_oval)
    else:
        xlim, ylim = detail["limits"][x], detail["limits"][y]
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    def coord(obs: str, raw: float) -> float:
        if detail is None:
            return float(raw) * UNIT_SCALE[obs]
        return (float(raw) - detail["origin"][obs]) / detail["units"][obs]
    handles = []
    if zoom == "experiment":
        # Diagonal-covariance Gaussian approximations: we have marginal
        # standard deviations, NOT the official joint NuFIT/JUNO surfaces.
        # The 2D radii give 68.27% and 99.73% probability coverage.
        r1 = JOINT_RADII["68% joint"]
        r3 = JOINT_RADII["99.73% joint"]
        _exp_ellipse(ax,nufit,x,y,r3,
                     facecolor=PALETTE["nufit_3"],edgecolor="#94B9D8",
                     alpha=.38,lw=1.0,zorder=.20)
        _exp_ellipse(ax,nufit,x,y,r1,
                     facecolor=PALETTE["nufit_1"],edgecolor="#28639B",
                     alpha=.32,lw=1.3,zorder=.30)
        handles.extend([
            Patch(facecolor=PALETTE["nufit_3"],edgecolor="#94B9D8",
                  label="NuFIT 6.1 ~3σ (joint Gaussian)"),
            Patch(facecolor=PALETTE["nufit_1"],edgecolor="#28639B",
                  label="NuFIT 6.1 ~1σ (joint Gaussian)"),
        ])
        if solar:
            for i,(record,label,opacity,_) in enumerate(juno_regions(juno_mode)):
                color = (PALETTE["juno_future"] if record is JUNO_FUTURE_6Y
                         else PALETTE["juno_1"])
                _exp_ellipse(ax,record,x,y,r1,
                             facecolor=color,edgecolor="#A14B22",
                             alpha=min(opacity,.28),lw=1.35,zorder=.40+.03*i)
                handles.append(Patch(facecolor=color,edgecolor="#A14B22",
                                     label=label+" (joint Gaussian)"))

    # Model response ellipse is NOT a statistical confidence interval.
    # A nearly rank-one response is genuinely a very thin oval.
    if show_model_oval:
        cc = centre(points)
        centre_xy = (coord(x, cc[x]), coord(y, cc[y]))
        sc = (np.diag([UNIT_SCALE[x], UNIT_SCALE[y]]) if detail is None else
              np.diag([1.0/detail["units"][x], 1.0/detail["units"][y]]))
        model_cov = sc@centered_covariance(points,x,y)@sc
        _cov_ellipse(ax,centre_xy,model_cov,nstd=2.,
                     edgecolor=PALETTE["model_oval"],facecolor="none",
                     linestyle=(0,(4,2)),linewidth=1.4,zorder=4.3)
        handles.append(Line2D([],[],color=PALETTE["model_oval"],
                              ls=(0,(4,2)),lw=1.5,
                              label="Model response oval (descriptive)"))

    for direction in DIRECTIONS:
        ps=points[direction]
        xx=np.asarray([coord(x,r[x]) for r in ps],dtype=float)
        yy=np.asarray([coord(y,r[y]) for r in ps],dtype=float)
        color=PALETTE[direction]
        markersize = 4.9 if detail is not None else 4.2
        ax.plot(xx, yy, color=color,lw=1.6,marker="o",markersize=markersize,
                markerfacecolor=color,markeredgecolor="white",
                markeredgewidth=.4,zorder=3.5)
        handles.append(Line2D([],[],color=color,lw=1.55,marker="o",
                              markersize=4.5,label=direction))
    # No black star, no separate centre marker; gamma=1 is already in the path.
    ax.set_xlabel(_detail_label(x,detail["powers"][x]) if detail else LABELS[x],fontsize=9.3)
    ax.set_ylabel(_detail_label(y,detail["powers"][y]) if detail else LABELS[y],fontsize=9.3)
    if detail is not None:
        # Coordinates are already in order-one, baseline-subtracted units.
        # Disable matplotlib's scientific-offset notation entirely: it was
        # responsible for 1e-7+3.088e-1 overlapping the tick labels.
        from matplotlib.ticker import MaxNLocator, FuncFormatter
        formatter = FuncFormatter(lambda val,_: f"{val:.3g}")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(MaxNLocator(nbins=4,steps=[1,2,2.5,5,10]))
            axis.set_major_formatter(formatter)
    if title:
        ax.set_title(_model_title(model), fontsize=11.1,pad=7)
    ax.grid(alpha=.16,color="#5C6471",linewidth=.55,zorder=.1)
    ax.tick_params(labelsize=8.1,direction="in",length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(.75);spine.set_color("#666666")
    if show_legend:
        ax.legend(handles=handles,ncol=2 if zoom=="experiment" else 4,
                  loc="lower center",bbox_to_anchor=(.5,1.09),
                  fontsize=7.5,frameon=True,framealpha=.92,
                  edgecolor="#D5D5D5",borderpad=.32,labelspacing=.36,
                  handlelength=1.5,columnspacing=.8)
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
           *, per_model=True, make_grid=True, model_oval=True,
           juno_mode="measured", make_detail=True) -> None:
    """Produce both 16-panel views; no numerical re-evaluation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True,exist_ok=True)
    views = ("experiment", "trajectory")
    view_names = {"experiment":"zoom_out","trajectory":"zoom_in"}

    # Individually readable pages, saved before the larger two 4x4 grids.
    if per_model:
        directory = out/"models"
        directory.mkdir(exist_ok=True)
        for model,data in plots.items():
            for zoom in views:
                if zoom == "trajectory" and not make_detail:
                    continue
                fig,ax=plt.subplots(figsize=(7.0,5.4))
                _plot_panel(ax,model,data,nufit,x,y,show_legend=True,
                            show_model_oval=model_oval,juno_mode=juno_mode,
                            zoom=zoom)
                fig.subplots_adjust(top=.77,bottom=.16,left=.15,right=.98)
                tag=model.replace("+","p").replace("-","m")
                _save_figure(fig,directory/f"{tag}_zee_{view_names[zoom]}.png",dpi=175)
                plt.close(fig)
            print(f"Saved {model} (overview + zoom)",flush=True)

    if not make_grid:
        return
    for zoom in views:
        # Both grids have the SAME 16-model ordering and gamma data.
        fig,axes=plt.subplots(4,4,figsize=(16.8,15.4))
        handles=None
        for ax,model in zip(axes.flat,plots):
            handles=_plot_panel(ax,model,plots[model],nufit,x,y,
                                show_model_oval=model_oval,juno_mode=juno_mode,
                                zoom=zoom)
        for ax in list(axes.flat)[len(plots):]:
            ax.set_axis_off()
        fig.subplots_adjust(left=.075,right=.975,top=.855,bottom=.055,
                            wspace=.42,hspace=.45)
        title=("Overview: experimental Gaussian ovals + T1/T2/T3 paths"
               if zoom=="experiment" else
               "Detail: compensated-coupling sensitivity for all 16 models")
        fig.text(.5,.984,title,ha="center",va="top",fontsize=16)
        example=next(iter(plots.values()))["T1"]
        gamma=[p["gamma"] for p in example]
        subtitle = f"16 individually optimised general-T3 fits | UV γ = {min(gamma):.3f}–{max(gamma):.3f}"
        if zoom == "trajectory":
            subtitle += " | δO = O(γ) − O(1); axis powers vary by model"
        fig.text(.5,.955,subtitle,ha="center",va="top",fontsize=10)
        if handles:
            fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.93),
                       ncol=min(len(handles),7),fontsize=8.2,frameon=True,
                       edgecolor="#CCCCCC",fancybox=False)
        _save_figure(fig,out/f"model_grid_zee_{view_names[zoom]}.png",dpi=170)
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
    p.add_argument("--no-model-oval", action="store_true",
                   help="Hide descriptive covariance response ovals (shown by default)")
    p.add_argument("--model-oval", action="store_true",
                   help="Compatibility option: ovals are enabled by default")
    p.add_argument("--juno-mode", choices=("measured", "future", "both"),
                   default="measured",
                   help="JUNO 2026 observed interval, optional re-centred six-year forecast, or both")
    p.add_argument("--no-detail", action="store_true",
                   help="Omit zoomed individual model pages (both 16-model grids still produced)")
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
           make_grid=not args.no_grid, model_oval=not args.no_model_oval,
           juno_mode=args.juno_mode,make_detail=not args.no_detail)
    manifest = {
        "format": "t3-general-zee-dual-view-gaussian-ovals-v5",
        "models": list(rows), "directions": list(DIRECTIONS),
        "gamma_grid": [pt["gamma"] for pt in next(iter(rows.values()))["T1"]],
        "input_files": {str(args.source / "points.csv"):_sha256(args.source / "points.csv"),
                        str(args.source / "summary.json"):_sha256(args.source / "summary.json"),
                        str(args.nufit):_sha256(args.nufit)},
        "nufit": {"source": nu.get("dataset"),
                  "central_values": nu["central_values"],
                  "one_sigma_errors": nu["one_sigma_errors"]},
        "juno_measured_2026": JUNO_2026,
        "juno_forecast_6y": JUNO_FUTURE_6Y,
        "juno_plot_mode": args.juno_mode,
        "experimental_oval_definition":
            "2D independent-error Gaussian approximations: NuFIT 6.1 (NO) "
            "68.27% / 99.73% (Delta chi2=2.30 / 11.83), JUNO 2026 measured "
            "solar 68.27% ellipse or optional six-year projection. "
            "NO covariance measurements supplied; NOT official joint contours.",
        "model_oval_definition":
            "default: 2x sample covariance standard-deviation radius about "
            "gamma=1 for all off-centre T1/T2/T3 points; descriptive, NOT "
            "a statistical confidence region; near-rank-one ovals not enlarged",
        "physics_scope": "saved conditional general-Lagrangian candidate UV/EFT1 + "
                         "one-scale C5, NOT certified all-fermions-first "
                         "then all-scalars sequential EFT matching",
        "plots_use_saved_points_only": True,
        "zee_fig9_inspired_dual_view_gaussian_ovals": True,
        "baseline_star_drawn": False,
        "arrows_drawn": False,
        "grid_footer_text_drawn": False,
        "trajectory_marker_size_pt": 4.9,
        "detail_coordinates": "(observable(gamma)-observable(gamma=1))/10**power; "
                              "power chosen separately for each panel and observable; "
                              "experimental overview remains in absolute coordinates",
        "grid_views": ["zoom_out", "zoom_in"] if not args.no_grid else [],
        "model_response_oval_drawn": not args.no_model_oval,
        "grid_drawn": not args.no_grid,
        "individual_close_up_drawn": not args.no_per_model and not args.no_detail,
        "per_model_figures_drawn": not args.no_per_model,
        "juno_solar_chi2_definition": "two marginal Gaussian residuals, diagonal; "
                                       "NOT official JUNO likelihood or combined global chi2",
        "local_slope_definition": "central difference dO/dln(gamma) using nearest gamma pair "
                                  "straddling 1; 1% change approximated linearly",
    }
    (args.output_dir / "provenance.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved Zee-style dual-view ovals for {len(rows)}/{len(selected)} complete models: {args.output_dir}")
    print(f"JUNO mode: {args.juno_mode}. Measured and projected intervals are never conflated.")
    print("No RGEs or matching rerun; no saved scan files modified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
