from __future__ import annotations

"""Generate reproducible tables and figures for the Ma running benchmark.

Run from the repository root with

    python -m RGE.phenomenology.MaRunningReport CONFIG --output-dir OUTPUT_DIR

The report uses the same one-loop UV and EFT integrations as
``MaFullRunningComparison`` and samples each solution only once.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from ..running.MaFullRunningComparison import (
        _threshold_state,
        compare_ma_running,
        state_from_config,
    )
    from ..running.MaUVRGE import sample_ma_uv_trajectory, takagi_majorana
    from ..running.NumericalWeinbergRGE import (
        neutrino_mass_matrix,
        sample_weinberg_trajectory,
    )
    from .NeutrinoDataComparison import mixing_angles_from_pmns_abs
    from .NeutrinoObservables import calculate_neutrino_observables
except ImportError:  # Allow direct execution from a flat development folder.
    from MaFullRunningComparison import _threshold_state, compare_ma_running, state_from_config
    from MaUVRGE import sample_ma_uv_trajectory, takagi_majorana
    from NumericalWeinbergRGE import neutrino_mass_matrix, sample_weinberg_trajectory
    from NeutrinoDataComparison import mixing_angles_from_pmns_abs
    from NeutrinoObservables import calculate_neutrino_observables


CASE_ORDER = ("frozen", "uv_only", "eft_only", "full")
CASE_LABELS = {
    "frozen": "Frozen",
    "uv_only": "UV only",
    "eft_only": "EFT only",
    "full": "Full",
}
COLORS = {
    "frozen": "#7A7A7A",
    "uv_only": "#2878B5",
    "eft_only": "#E07B39",
    "full": "#2A9D6F",
}


def _complex_matrix(payload: dict[str, Any]) -> np.ndarray:
    real = np.asarray(payload["real"], dtype=float)
    imag = np.asarray(payload.get("imag", np.zeros_like(real)), dtype=float)
    return real + 1j * imag


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _observable_rows(comparison: dict[str, Any], ordering: str) -> list[dict[str, Any]]:
    rows = []
    for case in CASE_ORDER:
        payload = comparison["cases"][case]
        mass_ev = _complex_matrix(payload["neutrino_mass_matrix_ev"])
        obs = calculate_neutrino_observables(mass_ev * 1.0e-9, ordering=ordering)
        angles = mixing_angles_from_pmns_abs(obs.pmns_abs)
        rows.append(
            {
                "case": case,
                "label": CASE_LABELS[case],
                "uv_running": case in ("uv_only", "full"),
                "eft_running": case in ("eft_only", "full"),
                "C5_norm_gev_inverse": payload["C5_frobenius_norm_gev_inverse"],
                "mass_norm_ev": payload["mass_frobenius_norm_ev"],
                "m1_ev": obs.masses_ev[0],
                "m2_ev": obs.masses_ev[1],
                "m3_ev": obs.masses_ev[2],
                "delta_m21_sq_ev2": obs.delta_m21_sq_ev2,
                "delta_m31_sq_ev2": obs.delta_m31_sq_ev2,
                "delta_m32_sq_ev2": obs.delta_m32_sq_ev2,
                "sin2_theta12": angles["sin2_theta12"],
                "sin2_theta23": angles["sin2_theta23"],
                "sin2_theta13": angles["sin2_theta13"],
                "takagi_residual": obs.takagi_residual,
            }
        )
    return rows


def _uv_rows(trajectory) -> list[dict[str, Any]]:
    rows = []
    for scale, state in zip(trajectory.scales_gev, trajectory.states):
        _, masses = takagi_majorana(state.M)
        root = np.sqrt(max(state.lambda1 * state.lambda2, 0.0))
        rows.append(
            {
                "scale_gev": scale,
                "log10_scale_gev": np.log10(scale),
                "gY": state.gY,
                "g2": state.g2,
                "g3": state.g3,
                "lambda1": state.lambda1,
                "lambda2": state.lambda2,
                "lambda3": state.lambda3,
                "lambda4": state.lambda4,
                "lambda5": state.lambda5,
                "mH2_gev2": state.mH2,
                "mEta2_gev2": state.mEta2,
                "sqrt_mEta2_gev": np.sqrt(max(state.mEta2, 0.0)),
                "max_abs_h": np.max(np.abs(state.h)),
                "M1_gev": masses[0],
                "M2_gev": masses[1],
                "M3_gev": masses[2],
                "bfb_lambda3_margin": state.lambda3 + root,
                "bfb_lambda345_margin": (
                    state.lambda3 + state.lambda4 - abs(state.lambda5) + root
                ),
            }
        )
    return rows


def _eft_rows(trajectory, scenario: str, vev_gev: float) -> list[dict[str, Any]]:
    rows = []
    for index, scale in enumerate(trajectory.scales_gev):
        mass_ev = neutrino_mass_matrix(trajectory.K[index], vev_gev=vev_gev) * 1.0e9
        masses = np.sort(np.linalg.svd(mass_ev, compute_uv=False))
        rows.append(
            {
                "scenario": scenario,
                "scale_gev": scale,
                "log10_scale_gev": np.log10(scale),
                "gY": trajectory.gY[index],
                "g2": trajectory.g2[index],
                "g3": trajectory.g3[index],
                "lambdaH": trajectory.lambdaH[index],
                "C5_norm_gev_inverse": np.linalg.norm(trajectory.K[index]),
                "mass_norm_ev": np.linalg.norm(mass_ev),
                "m_light_ev": masses[0],
                "m_middle_ev": masses[1],
                "m_heavy_ev": masses[2],
            }
        )
    return rows


def _save_figure(fig: plt.Figure, base_path: Path) -> list[str]:
    outputs = []
    for suffix in (".png", ".pdf"):
        path = base_path.with_suffix(suffix)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        outputs.append(path.name)
    plt.close(fig)
    return outputs


def _plot_effects(rows: list[dict[str, Any]], effects: dict[str, float], output: Path) -> list[str]:
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.6), constrained_layout=True)
    cases = CASE_ORDER
    labels = [CASE_LABELS[case] for case in cases]
    colors = [COLORS[case] for case in cases]

    shift_values = [
        100.0 * effects["uv_only_relative_to_frozen"],
        100.0 * effects["eft_only_relative_to_frozen"],
        100.0 * effects["full_relative_to_frozen"],
    ]
    bars = axes[0, 0].bar(labels[1:], shift_values, color=colors[1:])
    axes[0, 0].bar_label(bars, fmt="%.2f%%", padding=3, fontsize=9)
    axes[0, 0].set_ylabel(r"$\Vert\Delta C_5\Vert/\Vert C_5^{\rm frozen}\Vert$ (%)")
    axes[0, 0].set_title("Running effect on the Weinberg coefficient")

    x = np.arange(len(cases))
    width = 0.23
    for index, mass in enumerate(("m1_ev", "m2_ev", "m3_ev")):
        axes[0, 1].bar(
            x + (index - 1) * width,
            [row[mass] for row in rows],
            width,
            label=rf"$m_{index + 1}$",
        )
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_ylabel("Mass (eV)")
    axes[0, 1].set_title("Neutrino masses")
    axes[0, 1].legend(frameon=False, ncols=3)

    width = 0.34
    axes[1, 0].bar(x - width / 2, [row["delta_m21_sq_ev2"] for row in rows], width, label=r"$\Delta m^2_{21}$")
    axes[1, 0].bar(x + width / 2, [abs(row["delta_m31_sq_ev2"]) for row in rows], width, label=r"$|\Delta m^2_{31}|$")
    axes[1, 0].set_xticks(x, labels)
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylabel(r"Mass splitting (eV$^2$)")
    axes[1, 0].set_title("Mass-squared splittings")
    axes[1, 0].legend(frameon=False)

    for key, label, marker in (
        ("sin2_theta12", r"$\sin^2\theta_{12}$", "o"),
        ("sin2_theta23", r"$\sin^2\theta_{23}$", "s"),
        ("sin2_theta13", r"$\sin^2\theta_{13}$", "^"),
    ):
        axes[1, 1].plot(x, [row[key] for row in rows], marker=marker, label=label)
    axes[1, 1].set_xticks(x, labels)
    axes[1, 1].set_ylabel("Mixing observable")
    axes[1, 1].set_title("Leptonic mixing")
    axes[1, 1].legend(frameon=False)

    for ax in axes.flat:
        ax.grid(axis="y", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Ma benchmark: isolated UV and EFT running effects", fontsize=14)
    return _save_figure(fig, output / "ma_running_effects")


def _plot_uv(rows: list[dict[str, Any]], output: Path) -> list[str]:
    x = np.asarray([row["log10_scale_gev"] for row in rows])
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.6), constrained_layout=True)
    for key, label in (("gY", r"$g_Y$"), ("g2", r"$g_2$"), ("g3", r"$g_3$")):
        axes[0, 0].plot(x, [row[key] for row in rows], label=label)
    axes[0, 0].set_ylabel("Gauge coupling")
    axes[0, 0].legend(frameon=False)

    for index in range(1, 6):
        axes[0, 1].plot(x, [row[f"lambda{index}"] for row in rows], label=rf"$\lambda_{index}$")
    axes[0, 1].set_ylabel("Quartic coupling")
    axes[0, 1].legend(frameon=False, ncols=2)

    for index in range(1, 4):
        axes[1, 0].plot(x, [row[f"M{index}_gev"] for row in rows], label=rf"$M_{index}$")
    axes[1, 0].plot(x, [row["sqrt_mEta2_gev"] for row in rows], "--", label=r"$\sqrt{m_\eta^2}$")
    axes[1, 0].set_ylabel("Mass scale (GeV)")
    axes[1, 0].legend(frameon=False, ncols=2)

    axes[1, 1].plot(x, [row["max_abs_h"] for row in rows], color=COLORS["full"])
    axes[1, 1].set_ylabel(r"$\max |h_{ij}|$")

    titles = ("Gauge couplings", "Scalar quartics", "Heavy and inert scales", "Largest Ma Yukawa")
    for ax, title in zip(axes.flat, titles):
        ax.set_title(title)
        ax.set_xlabel(r"$\log_{10}(\mu/{\rm GeV})$")
        ax.grid(alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
        ax.invert_xaxis()
    fig.suptitle("Full-theory Ma running above the matching scale", fontsize=14)
    return _save_figure(fig, output / "ma_uv_trajectory")


def _plot_eft(rows: list[dict[str, Any]], output: Path) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    for scenario in ("eft_only", "full"):
        subset = [row for row in rows if row["scenario"] == scenario]
        x = [row["log10_scale_gev"] for row in subset]
        axes[0].plot(x, [row["C5_norm_gev_inverse"] for row in subset], label=CASE_LABELS[scenario], color=COLORS[scenario])
        axes[1].plot(x, [row["mass_norm_ev"] for row in subset], label=CASE_LABELS[scenario], color=COLORS[scenario])
    axes[0].set_ylabel(r"$\Vert C_5\Vert_F$ (GeV$^{-1}$)")
    axes[0].set_title("Weinberg coefficient")
    axes[1].set_ylabel(r"$\Vert m_\nu\Vert_F$ (eV)")
    axes[1].set_title("Neutrino mass norm")
    for ax in axes:
        ax.set_xlabel(r"$\log_{10}(\mu/{\rm GeV})$")
        ax.invert_xaxis()
        ax.grid(alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False)
    fig.suptitle("SMEFT running below the common threshold", fontsize=14)
    return _save_figure(fig, output / "ma_eft_trajectory")


def _fmt(value: float, digits: int = 6) -> str:
    return f"{float(value):.{digits}g}"


def _write_markdown(
    path: Path,
    config_name: str,
    comparison: dict[str, Any],
    rows: list[dict[str, Any]],
    generated: list[str],
) -> None:
    effects = comparison["effects"]
    lines = [
        "# Ma full-running comparison",
        "",
        f"Input configuration: `{config_name}`",
        "",
        "This benchmark separates running above and below one common matching threshold. "
        "The fitted point is a benchmark reconstruction of the selected neutrino target; it is not an independent prediction of those target values.",
        "",
        "## Main result",
        "",
        f"- UV-only relative shift in $C_5$: **{100 * effects['uv_only_relative_to_frozen']:.3f}%**",
        f"- EFT-only relative shift in $C_5$: **{100 * effects['eft_only_relative_to_frozen']:.3f}%**",
        f"- Full relative shift in $C_5$: **{100 * effects['full_relative_to_frozen']:.3f}%**",
        f"- UV–EFT interaction relative to frozen $C_5$: **{100 * effects['uv_eft_interaction_relative_to_frozen']:.3f}%**",
        f"- Full/frozen neutrino-mass norm ratio: **{effects['full_mass_norm_ratio_to_frozen']:.6f}**",
        "",
        "## Low-energy comparison",
        "",
        "| Case | UV | EFT | m1 (eV) | m2 (eV) | m3 (eV) | Δm²21 (eV²) | Δm²31 (eV²) | sin²θ12 | sin²θ23 | sin²θ13 |",
        "|---|:---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {label} | {uv} | {eft} | {m1} | {m2} | {m3} | {dm21} | {dm31} | {s12} | {s23} | {s13} |".format(
                label=row["label"],
                uv="yes" if row["uv_running"] else "no",
                eft="yes" if row["eft_running"] else "no",
                m1=_fmt(row["m1_ev"]),
                m2=_fmt(row["m2_ev"]),
                m3=_fmt(row["m3_ev"]),
                dm21=_fmt(row["delta_m21_sq_ev2"]),
                dm31=_fmt(row["delta_m31_sq_ev2"]),
                s12=_fmt(row["sin2_theta12"]),
                s23=_fmt(row["sin2_theta23"]),
                s13=_fmt(row["sin2_theta13"]),
            )
        )
    lines += [
        "",
        "## Reproducibility and limitations",
        "",
        "The calculation uses one-loop Ma RGEs above the threshold, small-$\\lambda_5$ matching, and one-loop SM+Weinberg running below it. "
        "All heavy fields are integrated out at one common scale; split thresholds and finite threshold corrections are not included.",
        "",
        "Generated files:",
        "",
    ]
    lines.extend(f"- `{name}`" for name in generated)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_ma_running_report(
    config: dict[str, Any],
    output_dir: Path,
    *,
    samples: int = 64,
    config_name: str = "config.json",
) -> dict[str, Any]:
    """Generate JSON, CSV, Markdown, PNG and PDF report artifacts."""

    if samples < 2:
        raise ValueError("At least two trajectory samples are required.")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    comparison = compare_ma_running(config)
    comparison_path = output_dir / "ma_full_running_comparison.json"
    comparison_path.write_text(json.dumps(comparison, indent=2) + "\n", encoding="utf-8")

    target = config.get("fit_target", {})
    ordering = str(target.get("ordering", "AUTO"))
    observable_rows = _observable_rows(comparison, ordering)
    _write_csv(output_dir / "ma_case_comparison.csv", observable_rows)

    mu_uv = float(config["mu_uv_gev"])
    mu_match = float(config["mu_matching_gev"])
    mu_low = float(config.get("mu_low_gev", 91.1876))
    vev_gev = float(config.get("vev_gev", 246.22))
    high_state = state_from_config(config["uv"])
    uv_trajectory = sample_ma_uv_trajectory(high_state, mu_uv, mu_match, samples=samples)
    uv_rows = _uv_rows(uv_trajectory)
    _write_csv(output_dir / "ma_uv_trajectory.csv", uv_rows)

    frozen_sm, _ = _threshold_state(high_state)
    evolved_sm, _ = _threshold_state(uv_trajectory.states[-1])
    eft_only = sample_weinberg_trajectory(frozen_sm, mu_match, mu_low, samples=samples)
    full = sample_weinberg_trajectory(evolved_sm, mu_match, mu_low, samples=samples)
    eft_rows = _eft_rows(eft_only, "eft_only", vev_gev) + _eft_rows(full, "full", vev_gev)
    _write_csv(output_dir / "ma_eft_trajectory.csv", eft_rows)

    generated = [
        comparison_path.name,
        "ma_case_comparison.csv",
        "ma_uv_trajectory.csv",
        "ma_eft_trajectory.csv",
    ]
    generated += _plot_effects(observable_rows, comparison["effects"], output_dir)
    generated += _plot_uv(uv_rows, output_dir)
    generated += _plot_eft(eft_rows, output_dir)
    generated.append("ma_running_report.md")
    _write_markdown(
        output_dir / "ma_running_report.md",
        config_name,
        comparison,
        observable_rows,
        generated,
    )

    manifest = {
        "status": "Success",
        "input_config": config_name,
        "trajectory_samples_per_regime": samples,
        "fit_target": target or None,
        "effects": comparison["effects"],
        "files": generated + ["ma_running_report_manifest.json"],
    }
    (output_dir / "ma_running_report_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="fitted Ma benchmark JSON")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=64)
    args = parser.parse_args()

    with args.config.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    result = generate_ma_running_report(
        config,
        args.output_dir,
        samples=args.samples,
        config_name=args.config.name,
    )
    effects = result["effects"]
    print(f"Wrote report to {args.output_dir}")
    print(
        "Relative C5 shifts: "
        f"UV={effects['uv_only_relative_to_frozen']:.6g}, "
        f"EFT={effects['eft_only_relative_to_frozen']:.6g}, "
        f"full={effects['full_relative_to_frozen']:.6g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
