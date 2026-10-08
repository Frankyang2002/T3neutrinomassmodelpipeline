"""Plot fixed-order scalar-only LLSS Wilson transport from existing T3 outputs.

This is an *analytic, one-generation/group-component* diagnostic.  It reuses the
existing symbolic transport; it does not solve another RGE, insert the LLSS
self-running into final C5, or reduce full-flavour Yukawa matrices to scalars.

Without explicit numerical assignments only the universal logarithmic transport
kernel is graphable.  Exact model-dependent C0 and beta expressions are saved.
With complete one-generation assignments, actual component curves are plotted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from RGE.running.intermediate.ScalarOnlyWilsonFlow import run_component_wilson_transport

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "output" / "full" / "comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full" / "within_model_comparison"
SCENARIO = "smallY_smallL"  # Symbolic coefficients are not benchmark-evaluated.
REQUIRED = (
    "eft1_after_F_wilson_seed.json",
    "eft1_wilson_rge.json",
    "eft1_rgbeta_rge.json",
)
HBAR = 1.0 / (16.0 * math.pi**2)


def _signature(paths: dict[str, Path], mu_f: float, mu_s: float,
               values: dict[str, Any] | None) -> str:
    hasher = hashlib.sha256()
    for name in REQUIRED:
        hasher.update(name.encode("utf-8"))
        hasher.update(paths[name].read_bytes())
    hasher.update(json.dumps([mu_f, mu_s, values], sort_keys=True).encode("utf-8"))
    return hasher.hexdigest()


def _save_plot(output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output, dpi=200, bbox_inches="tight")
    plt.close()


def plot_transport_kernel(mu_f: float, mu_s: float, output: Path,
                          component_count: int) -> None:
    """Plot the exact common scale factor, not an invented Wilson norm."""
    mu = np.geomspace(mu_f, mu_s, 160)
    eta = np.log(mu / mu_f) * HBAR
    plt.figure(figsize=(8, 4.7))
    plt.plot(mu, eta, linewidth=2)
    plt.axhline(0, linewidth=0.7, linestyle="--")
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$[C_{ijab}(\mu)-C_{ijab}(M_F)]/\beta_{ijab}^{(1)}$")
    plt.title("Scalar-only LLSS: universal fixed-order transport kernel")
    plt.grid(True, alpha=0.25)
    plt.figtext(0.5, -0.015,
                f"{component_count} component(s) with nonzero running; "
                "symbolic one-generation diagnostic, not a Wilson norm",
                ha="center", fontsize=8)
    _save_plot(output)


def plot_component_support(transport: dict[str, Any], output: Path) -> None:
    """Expose the model-dependent LLSS component/mixing structure exactly."""
    counts = [
        int(transport["tree_boundary_component_count"]),
        int(transport["beta_component_count"]),
        int(transport["generated_by_running_component_count"]),
    ]
    labels = ["Tree boundary", "Nonzero one-loop beta", "Generated from zero"]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    bars = ax.bar(labels, counts, width=0.60)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                str(count), ha="center", va="bottom")
    ax.set_ylim(0, max(1, max(counts)) * 1.16)
    ax.set_ylabel("Number of independent LLSS tensor components")
    ax.set_title("Model-dependent LLSS operator support and mixing")
    ax.grid(axis="y", alpha=0.23)
    fig.text(0.5, 0.015,
             "Counted from the existing symbolic Wilson seed and tensor RGE",
             ha="center", fontsize=8)
    _save_plot(output)


def _complex_value(raw: str, substitutions: dict[sp.Symbol, sp.Expr]) -> complex:
    expr = sp.sympify(raw, locals={"conjugate": sp.conjugate})
    val = expr.subs(substitutions).evalf(16)
    if val.free_symbols:
        raise ValueError("Unassigned symbols: " + ", ".join(sorted(map(str, val.free_symbols))))
    number = complex(val)
    if not np.isfinite(number.real) or not np.isfinite(number.imag):
        raise ValueError("A symbolic component evaluated to a non-finite number.")
    return number


def numerical_components(transport: dict[str, Any], mu_f: float, mu_s: float,
                         values: dict[str, Any]) -> dict[str, tuple[complex, complex]]:
    """Numerically evaluate *only* when every symbolic input was specified."""
    substitutions = {sp.Symbol(name): sp.sympify(str(value)) for name, value in values.items()}
    tree = transport["tree_boundary_components"]
    corrections = transport["one_loop_running_components"]
    log_interval = math.log(mu_s / mu_f)
    result: dict[str, tuple[complex, complex]] = {}
    missing: set[str] = set()
    for key in sorted(set(tree) | set(corrections)):
        try:
            c0 = _complex_value(tree.get(key, "0"), substitutions)
            running_at_s = _complex_value(corrections.get(key, "0"), substitutions)
        except ValueError as exc:
            if str(exc).startswith("Unassigned symbols: "):
                missing.update(str(exc).removeprefix("Unassigned symbols: ").split(", "))
                continue
            raise ValueError(f"Component {key}: {exc}") from exc
        result[key] = (c0, running_at_s / log_interval)
    if missing:
        raise ValueError("Missing one-generation assignments: " + ", ".join(sorted(missing)))
    if not result:
        raise ValueError("No Wilson components available to plot.")
    return result


def plot_components(components: dict[str, tuple[complex, complex]],
                    mu_f: float, mu_s: float, output: Path,
                    max_components: int = 8) -> list[str]:
    mu = np.geomspace(mu_f, mu_s, 120)
    log_ratio = np.log(mu / mu_f)
    ranked = sorted(components, key=lambda name: max(
        abs(components[name][0]),
        abs(components[name][0] + HBAR * math.log(mu_s / mu_f) * components[name][1]),
    ), reverse=True)
    chosen = ranked[:max_components]
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    for key in chosen:
        c0, beta = components[key]
        c_mu = c0 + HBAR * log_ratio * beta
        line, = ax.plot(mu, np.abs(c_mu), linewidth=1.6, label=f"{key}")
        ax.plot(mu, np.full(mu.shape, abs(c0)), linestyle="--", linewidth=0.7,
                alpha=0.5, color=line.get_color())
    ax.set_xscale("log")
    if any(max(abs(components[key][0]), abs(components[key][1])) > 0 for key in chosen):
        ax.set_yscale("symlog", linthresh=1e-25)
    ax.set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
    ax.set_ylabel(r"$|C_{ijab}(\mu)|$ [matching-seed units]")
    ax.set_title("Fixed-order LLSS Wilson components (explicit one-generation input)")
    ax.grid(True, alpha=0.22)
    ax.legend(title="Component (i,j,a,b)", fontsize=8, ncol=2)
    fig.text(0.5, 0.01, "Solid: fixed-order; dashed: tree boundary. Not full-flavour C5.",
             ha="center", fontsize=8)
    _save_plot(output)
    return chosen


def _inputs(data_dir: Path) -> dict[str, Path]:
    paths = {name: data_dir / name for name in REQUIRED}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing saved symbolic inputs: " + "; ".join(missing))
    return paths


def _report_model_key(data_dir: Path) -> str:
    """Use the same diagnostic model key as WithinModelComparison."""
    diagnostics = data_dir / "running_diagnostics.json"
    if diagnostics.is_file():
        payload = json.loads(diagnostics.read_text(encoding="utf-8"))
        benchmark = payload.get("benchmark_search", {})
        rep = benchmark.get("model") if isinstance(benchmark, dict) else None
        if isinstance(rep, dict) and rep.get("model_key"):
            return str(rep["model_key"])
    return data_dir.parent.name


def report_model(data_dir: Path, output_dir: Path, *, mu_f: float,
                 mu_s: float, values: dict[str, Any] | None = None,
                 force: bool = False) -> Path:
    if not (math.isfinite(mu_f) and math.isfinite(mu_s) and mu_f > mu_s > 0):
        raise ValueError("Require finite scales M_F > M_S > 0.")
    paths = _inputs(data_dir)
    signature = _signature(paths, mu_f, mu_s, values)
    report_path = output_dir / "llss_wilson_report.json"
    kernel_path = output_dir / "llss_transport_kernel.png"
    support_path = output_dir / "llss_component_support.png"
    component_path = output_dir / "llss_component_running.png"
    readme_path = output_dir / "LLSS_README.md"
    if report_path.is_file() and not force:
        existing = json.loads(report_path.read_text(encoding="utf-8"))
        expected = [report_path, kernel_path, support_path, readme_path]
        if values is not None:
            expected.append(component_path)
        if existing.get("input_signature") == signature and all(p.is_file() for p in expected):
            return report_path

    transport = run_component_wilson_transport(
        wilson_seed_path=paths["eft1_after_F_wilson_seed.json"],
        wilson_rge_path=paths["eft1_wilson_rge.json"],
        rgbeta_path=paths["eft1_rgbeta_rge.json"],
        mu_high=mu_f, mu_low=mu_s, output_path=None,
    )
    if transport.get("status") != "Success":
        raise RuntimeError("The stored symbolic LLSS transport was not successful.")
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_transport_kernel(mu_f, mu_s, kernel_path,
                          int(transport["running_correction_component_count"]))
    plot_component_support(transport, support_path)
    plotted: list[str] = []
    if values is not None:
        components = numerical_components(transport, mu_f, mu_s, values)
        plotted = plot_components(components, mu_f, mu_s, component_path)
    else:
        component_path.unlink(missing_ok=True)  # Avoid stale numerical displays.

    report = {
        "status": "Success", "diagnostic_only": True,
        "representation": "symbolic one-generation LLSS component tensor",
        "rg_improved": False, "used_in_authoritative_final_c5": False,
        "input_signature": signature,
        "source_paths": {key: str(path) for key, path in paths.items()},
        "mu_fermion_gev": mu_f, "mu_scalar_gev": mu_s,
        "numerical_assignments": values,
        "numeric_component_curves": plotted,
        "transport": transport,
    }
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    readme = [
        "# LLSS Wilson running — fixed-order diagnostic", "",
        f"Source: `{data_dir}`", "",
        "The main figure, `llss_transport_kernel.png`, plots only the exact",
        "universal transport factor, **not** a numerical Wilson coefficient.",
        "The companion `llss_component_support.png` gives the model-dependent",
        "component counts at tree level, in the one-loop RGE, and those",
        "generated from a vanishing boundary coefficient:", "",
        r"\[C_{ijab}(\mu)=C_{ijab}^{(0)}+\frac{\ln(\mu/M_F)}{16\pi^2}\beta_{ijab}^{(1)}.\]", "",
        f"Thresholds: M_F={mu_f:g} GeV; M_S={mu_s:g} GeV.",
        f"Tree components: {transport['tree_boundary_component_count']}.",
        f"Nonzero beta components: {transport['beta_component_count']}.",
        f"Running components: {transport['running_correction_component_count']}.",
        f"Generated by mixing: {transport['generated_by_running_component_count']}.", "",
        "The complete *model-dependent* symbolic tree and beta expressions are in",
        "`llss_wilson_report.json` under `transport`. The serialized",
        "`one_loop_running_components` are **beta times ln(M_S/M_F)**,",
        "i.e. the **1/(16pi^2)** factor has not yet been applied.", "",
    ]
    if values is None:
        readme += [
            "No numerical component curves were generated, because the benchmark",
            "uses full-flavour Yukawa matrices whereas this tensor diagnostic",
            "uses one-generation symbolic couplings. No unspecified values were",
            "silently chosen. Pass `--values-json` containing complete symbolic",
            "one-generation parameter assignments to plot actual component curves.", "",
        ]
    else:
        readme += [
            "`llss_component_running.png` contains |C_ijab(mu)| for explicit",
            "one-generation assignments. Solid = fixed-order; dashed = tree.",
            "These assignments are not a three-flavour benchmark projection.",
            f"Plotted components: {', '.join(plotted)}.", "",
        ]
    readme += [
        "LLSS self-running affects the final one-loop C5 only at the next",
        "order if inserted through the scalar loop; this diagnostic does not",
        "change the authoritative low-energy neutrino calculation.", "",
    ]
    readme_path.write_text("\n".join(readme), encoding="utf-8")
    return report_path


def build_reports(input_root: Path = DEFAULT_INPUT, output_root: Path = DEFAULT_OUTPUT,
                  *, scenario: str = SCENARIO, model: str | None = None,
                  mu_f: float = 1e5, mu_s: float = 1e3,
                  values: dict[str, Any] | None = None, force: bool = False) -> list[Path]:
    scenario_root = Path(input_root) / scenario
    if not scenario_root.is_dir():
        raise FileNotFoundError(f"Missing scenario directory: {scenario_root}")
    input_paths = sorted(scenario_root.rglob("eft1_wilson_rge.json"))
    inputs = [(_report_model_key(p.parent), p) for p in input_paths]
    if model is not None:
        inputs = [(key, p) for key, p in inputs
                  if key == model or p.parent.parent.name == model]
    if not inputs:
        raise FileNotFoundError(
            f"No saved LLSS RGE inputs under {scenario_root} for {model or 'any model'}."
        )
    names = [key for key, _ in inputs]
    if len(names) != len(set(names)):
        raise RuntimeError("Multiple saved Wilson RGEs map to the same model report key.")
    reports: list[Path] = []
    failed: list[str] = []
    for key, path in inputs:
        try:
            report = report_model(path.parent, Path(output_root) / key,
                                  mu_f=mu_f, mu_s=mu_s, values=values, force=force)
        except Exception as exc:
            failed.append(f"{key}: {type(exc).__name__}: {exc}")
            print(f"FAILED {failed[-1]}")
            continue
        reports.append(report)
        print(f"LLSS report: {report}")
    print(f"LLSS reports: {len(reports)}/{len(inputs)} succeeded")
    if failed:
        raise RuntimeError("LLSS report failures:\n" + "\n".join(failed))
    return reports


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--scenario", default=SCENARIO)
    parser.add_argument("--model", default=None, help="Optional exact saved model folder name")
    parser.add_argument("--mu-f", type=float, default=1e5, help="F threshold [GeV]")
    parser.add_argument("--mu-s", type=float, default=1e3, help="Scalar threshold [GeV]")
    parser.add_argument("--values-json", type=Path, default=None,
                        help="Explicit one-generation {symbol: numerical value} mapping")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    values = None
    if args.values_json is not None:
        values = json.loads(args.values_json.read_text(encoding="utf-8"))
        if not isinstance(values, dict):
            parser.error("--values-json must contain a JSON object")
    build_reports(args.input_root, args.output_root, scenario=args.scenario,
                  model=args.model, mu_f=args.mu_f, mu_s=args.mu_s,
                  values=values, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
