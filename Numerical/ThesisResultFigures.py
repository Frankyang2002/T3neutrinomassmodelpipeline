"""
Generate compact thesis-result figures and tables from BestFitDiagnostics JSON.

Usage
-----
python -m Numerical.ThesisResultFigures BEST_FIT_PHYSICS_JSON

Optional:
python -m Numerical.ThesisResultFigures BEST_FIT_PHYSICS_JSON --output-dir PATH

Outputs
-------
- neutrino_mass_running.png
- c5_running.png
- mnu_matrix_abs_ev.png
- pmns_abs.png
- neutrino_mass_spectrum.png
- best_fit_results.csv
- best_fit_results_table.tex
- best_fit_results_table.md
- figure_captions.txt
- summary.txt

This module does not recompute T3 physics. It formats the outputs already
written by Numerical.BestFitDiagnostics.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    if payload.get("status") != "Success":
        raise ValueError("Best-fit physics JSON must have status='Success'.")
    return payload


def _save(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def _complex_matrix(block: dict[str, Any]) -> np.ndarray:
    real = np.asarray(block["real"], dtype=float)
    imag = np.asarray(block["imag"], dtype=float)
    matrix = real + 1j * imag
    if matrix.shape != (3, 3):
        raise ValueError("Expected a 3x3 complex matrix.")
    return matrix


def plot_neutrino_mass_running(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    running = payload["running"]
    mu = np.asarray(running["mu_gev"], dtype=float)
    masses = np.asarray(running["masses_ev"], dtype=float)

    if masses.shape != (mu.size, 3):
        raise ValueError("running.masses_ev has unexpected shape.")

    plt.figure(figsize=(7.2, 4.8))
    labels = [r"$m_1$", r"$m_2$", r"$m_3$"]
    for index, label in enumerate(labels):
        plt.plot(mu, masses[:, index], label=label, linewidth=1.4)

    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel("Neutrino mass [eV]")
    plt.title("Neutrino-mass running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_c5_running(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    running = payload["running"]
    mu = np.asarray(running["mu_gev"], dtype=float)
    c5_abs = np.asarray(running["c5_abs"], dtype=float)

    if c5_abs.shape != (mu.size, 3, 3):
        raise ValueError("running.c5_abs has unexpected shape.")

    plt.figure(figsize=(7.6, 5.2))
    for i in range(3):
        for j in range(i, 3):
            values = c5_abs[:, i, j]
            positive = values > 0.0
            if not np.any(positive):
                continue
            plt.plot(
                mu[positive],
                values[positive],
                label=rf"$|C_5^{{{i+1}{j+1}}}|$",
                linewidth=1.0,
            )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$|C_5^{ij}|$ [GeV$^{-1}$]")
    plt.title(r"Weinberg-coefficient running")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save(output_path)


def _plot_heatmap(
    matrix: np.ndarray,
    *,
    title: str,
    labels_x: list[str],
    labels_y: list[str],
    output_path: Path,
    value_format: str,
    colorbar_label: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    values = np.asarray(matrix, dtype=float)

    plt.figure(figsize=(5.4, 4.8))
    image = plt.imshow(values, vmin=vmin, vmax=vmax)
    plt.colorbar(image, label=colorbar_label)
    plt.xticks(range(3), labels_x)
    plt.yticks(range(3), labels_y)
    plt.title(title)

    for i in range(3):
        for j in range(3):
            plt.text(
                j,
                i,
                value_format.format(values[i, j]),
                ha="center",
                va="center",
            )

    _save(output_path)


def plot_mnu_matrix(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    block = payload["low_energy"]["mnu_charged_lepton_basis_gev"]
    matrix_ev = np.abs(_complex_matrix(block)) * 1.0e9

    _plot_heatmap(
        matrix_ev,
        title=r"$|m_\nu|$ in the charged-lepton basis",
        labels_x=[r"$e$", r"$\mu$", r"$\tau$"],
        labels_y=[r"$e$", r"$\mu$", r"$\tau$"],
        output_path=output_path,
        value_format="{:.4f}",
        colorbar_label="Magnitude [eV]",
    )


def plot_pmns(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    pmns = np.asarray(payload["low_energy"]["pmns_abs"], dtype=float)

    _plot_heatmap(
        pmns,
        title=r"$|U_{\rm PMNS}|$",
        labels_x=["1", "2", "3"],
        labels_y=[r"$e$", r"$\mu$", r"$\tau$"],
        output_path=output_path,
        value_format="{:.4f}",
        colorbar_label=r"$|U_{\rm PMNS}|$",
        vmin=0.0,
        vmax=1.0,
    )


def plot_mass_spectrum(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    masses = np.asarray(payload["low_energy"]["masses_ev"], dtype=float)

    plt.figure(figsize=(6.2, 4.6))
    x = np.arange(3)
    plt.bar(x, masses)
    plt.xticks(x, [r"$m_1$", r"$m_2$", r"$m_3$"])
    plt.ylabel("Mass [eV]")
    plt.title("Low-energy neutrino mass spectrum")
    plt.grid(True, axis="y", alpha=0.25)
    _save(output_path)


def _table_rows(payload: dict[str, Any]) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []

    rows.append(("fit", r"$\chi^2_{\rm min}$", f"{payload['stored_best_chi2']:.6g}"))
    rows.append(("fit", "ordering", str(payload["ordering"])))

    parameters = payload["best_parameters"]
    rows.append(("parameter", r"$\lambda_{T3}$", f"{float(parameters['lambdaT3_real']):.9g}"))

    for prefix in ("y1", "y2"):
        for i in range(1, 4):
            for j in range(1, 4):
                key = f"{prefix}_{i}{j}"
                label = rf"$({prefix.upper()})_{{{i}{j}}}$"
                rows.append(("parameter", label, f"{float(parameters[key]):.9g}"))

    low = payload["low_energy"]
    masses = low["masses_ev"]
    for index, value in enumerate(masses, start=1):
        rows.append(("neutrino", rf"$m_{index}$ [eV]", f"{float(value):.9g}"))
    rows.append(("neutrino", r"$\sum_i m_i$ [eV]", f"{float(low['sum_masses_ev']):.9g}"))
    rows.append(("neutrino", r"$\Delta m_{21}^2$ [eV$^2$]", f"{float(low['delta_m21_sq_ev2']):.9g}"))
    rows.append(("neutrino", r"$\Delta m_{31}^2$ [eV$^2$]", f"{float(low['delta_m31_sq_ev2']):.9g}"))
    rows.append(("neutrino", r"$\Delta m_{32}^2$ [eV$^2$]", f"{float(low['delta_m32_sq_ev2']):.9g}"))

    thresholds = payload["threshold_diagnostics"]
    for index, value in enumerate(thresholds["fermion_singular_masses_gev"], start=1):
        rows.append(("threshold", rf"$M_{{F,{index}}}$ [GeV]", f"{float(value):.9g}"))
    for index, value in enumerate(thresholds["scalar_running_masses_gev"], start=1):
        rows.append(("threshold", rf"$M_{{S{index}}}(\mu_S)$ [GeV]", f"{float(value):.9g}"))

    return rows


def write_csv(
    rows: list[tuple[str, str, str]],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "quantity", "value"])
        writer.writerows(rows)


def write_markdown(
    rows: list[tuple[str, str, str]],
    output_path: Path,
) -> None:
    lines = [
        "| Section | Quantity | Value |",
        "|---|---|---:|",
    ]
    for section, quantity, value in rows:
        lines.append(f"| {section} | {quantity} | {value} |")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_latex(
    rows: list[tuple[str, str, str]],
    output_path: Path,
) -> None:
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Best-fit parameters and derived neutrino observables for the current T3-B, $\alpha=-1$ benchmark. Scalar masses are running masses evaluated at the scalar matching scale.}",
        r"\label{tab:t3-b-bestfit}",
        r"\begin{tabular}{lll}",
        r"\hline",
        r"Section & Quantity & Value \\",
        r"\hline",
    ]

    previous = None
    for section, quantity, value in rows:
        if previous is not None and section != previous:
            lines.append(r"\hline")
        safe_section = section.replace("_", r"\_")
        lines.append(f"{safe_section} & {quantity} & {value} \\\\")
        previous = section

    lines.extend(
        [
            r"\hline",
            r"\end{tabular}",
            r"\end{table}",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_captions(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    mu_s = float(payload["scales_gev"]["mu_scalar_threshold"])
    mu_low = float(payload["scales_gev"]["mu_low"])

    masses_high = np.asarray(payload["running"]["masses_ev"][0], dtype=float)
    masses_low = np.asarray(payload["running"]["masses_ev"][-1], dtype=float)
    frac = 100.0 * (1.0 - masses_low / masses_high)

    text = f"""Suggested figure captions

neutrino_mass_running.png
Scale dependence of the three neutrino mass eigenvalues in the final SM+Weinberg EFT, evolved from the scalar matching scale mu_S={mu_s:.3e} GeV to mu_low={mu_low:.3e} GeV. For this benchmark each mass decreases by approximately {np.mean(frac):.1f}% over the interval.

c5_running.png
Renormalisation-group evolution of the independent magnitudes |C5^ij| in the final SM+Weinberg EFT between the scalar matching scale and the low scale.

mnu_matrix_abs_ev.png
Absolute value of the low-energy Majorana neutrino mass matrix in the charged-lepton mass basis, in eV.

pmns_abs.png
Absolute values of the PMNS matrix obtained from the Takagi diagonalisation of the low-energy neutrino mass matrix.

neutrino_mass_spectrum.png
Low-energy normal-ordered neutrino mass spectrum at mu_low.
"""
    output_path.write_text(text, encoding="utf-8")


def write_summary(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    low = payload["low_energy"]
    masses = np.asarray(low["masses_ev"], dtype=float)
    running = np.asarray(payload["running"]["masses_ev"], dtype=float)
    changes = 100.0 * (1.0 - running[-1] / running[0])

    lines = [
        "T3-B alpha=-1 best-fit thesis-results summary",
        "",
        f"chi2_min = {float(payload['stored_best_chi2']):.12g}",
        f"ordering = {payload['ordering']}",
        f"m1 [eV] = {masses[0]:.12g}",
        f"m2 [eV] = {masses[1]:.12g}",
        f"m3 [eV] = {masses[2]:.12g}",
        f"sum m_i [eV] = {float(low['sum_masses_ev']):.12g}",
        f"Delta m21^2 [eV^2] = {float(low['delta_m21_sq_ev2']):.12g}",
        f"Delta m31^2 [eV^2] = {float(low['delta_m31_sq_ev2']):.12g}",
        f"Takagi residual = {float(low['takagi_residual']):.12g}",
        "",
        "Final-EFT fractional mass decreases:",
        f"m1: {changes[0]:.6g} %",
        f"m2: {changes[1]:.6g} %",
        f"m3: {changes[2]:.6g} %",
        "",
        "Caveat:",
        "These results characterize the current benchmark configuration and fitted oscillation target.",
        "They should not be presented as a complete phenomenological validation of the model.",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_thesis_results(
    best_fit_json: Path,
    output_dir: Path,
) -> Path:
    payload = _load_json(best_fit_json)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_neutrino_mass_running(
        payload,
        output_dir / "neutrino_mass_running.png",
    )
    plot_c5_running(
        payload,
        output_dir / "c5_running.png",
    )
    plot_mnu_matrix(
        payload,
        output_dir / "mnu_matrix_abs_ev.png",
    )
    plot_pmns(
        payload,
        output_dir / "pmns_abs.png",
    )
    plot_mass_spectrum(
        payload,
        output_dir / "neutrino_mass_spectrum.png",
    )

    rows = _table_rows(payload)
    write_csv(rows, output_dir / "best_fit_results.csv")
    write_markdown(rows, output_dir / "best_fit_results_table.md")
    write_latex(rows, output_dir / "best_fit_results_table.tex")
    write_captions(payload, output_dir / "figure_captions.txt")
    write_summary(payload, output_dir / "summary.txt")

    print(f"Thesis results written to: {output_dir}")
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate compact thesis figures and tables from "
            "BestFitDiagnostics JSON."
        )
    )
    parser.add_argument("best_fit_json", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Default: "
            "output/thesis_results/<input stem>/"
        ),
    )
    args = parser.parse_args()

    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else Path("output") / "thesis_results" / args.best_fit_json.stem
    )

    generate_thesis_results(
        args.best_fit_json,
        output_dir,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
