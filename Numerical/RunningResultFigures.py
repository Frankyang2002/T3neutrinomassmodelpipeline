"""Publication-oriented figures from ``RunningDiagnostics`` JSON.

The plots are chosen to show the quantities most relevant to the multi-threshold T3 running calculation:

- what runs significantly above/between thresholds;
- how the final Weinberg coefficient and m_nu evolve;
- whether physical neutrino mass splittings and mixing angles change.

No T3 physics is recomputed here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    if payload.get("status") != "Success":
        raise ValueError("Diagnostics JSON must have status='Success'.")
    return payload


def _save(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def plot_uv_coupling_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["uv_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(mu, data["y1_frobenius_norm"], label=r"$\|y_1\|_F$", linewidth=1.4)
    plt.plot(mu, data["y2_frobenius_norm"], label=r"$\|y_2\|_F$", linewidth=1.4)
    plt.plot(mu, data["lambdaT3_abs"], label=r"$|\lambda_{T3}|$", linewidth=1.4)
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel("Dimensionless coupling measure")
    plt.title("UV T3 coupling running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_intermediate_scalar_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["intermediate_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    for key, label in (
        ("lambdaH1", r"$\lambda_{H1}$"),
        ("lambdaH2", r"$\lambda_{H2}$"),
        ("lambda12", r"$\lambda_{12}$"),
        ("lambdaT3_abs", r"$|\lambda_{T3}|$"),
    ):
        plt.plot(mu, data[key], label=label, linewidth=1.3)
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel("Quartic coupling")
    plt.title("Scalar-only intermediate-EFT running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_neutrino_mass_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["final_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)
    masses = np.asarray(data["masses_ev"], dtype=float)
    if masses.shape != (mu.size, 3):
        raise ValueError("final_running.masses_ev has unexpected shape.")

    plt.figure(figsize=(7.2, 4.8))
    for index, label in enumerate((r"$m_1$", r"$m_2$", r"$m_3$")):
        plt.plot(mu, masses[:, index], label=label, linewidth=1.4)
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel("Neutrino mass [eV]")
    plt.title("Neutrino-mass running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_mass_splitting_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["final_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)
    dm21 = np.asarray(data["delta_m21_sq_ev2"], dtype=float)
    dm3l = np.asarray(data["delta_m3l_sq_ev2"], dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(mu, dm21, label=r"$\Delta m_{21}^2$", linewidth=1.4)
    plt.plot(mu, np.abs(dm3l), label=r"$|\Delta m_{3\ell}^2|$", linewidth=1.4)
    if np.all(dm21 > 0.0) and np.all(np.abs(dm3l) > 0.0):
        plt.yscale("log")
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"Mass-squared splitting [eV$^2$]")
    plt.title("Neutrino mass-splitting running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_mixing_angle_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["final_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)

    plt.figure(figsize=(7.2, 4.8))
    for key, label in (
        ("sin2_theta12", r"$\sin^2\theta_{12}$"),
        ("sin2_theta13", r"$\sin^2\theta_{13}$"),
        ("sin2_theta23", r"$\sin^2\theta_{23}$"),
    ):
        plt.plot(mu, data[key], label=label, linewidth=1.4)
    plt.xscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$\sin^2\theta_{ij}$")
    plt.title("Neutrino mixing-angle running")
    plt.grid(True, alpha=0.25)
    plt.legend()
    _save(output_path)


def plot_c5_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["final_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)
    c5 = np.asarray(data["c5_abs"], dtype=float)
    if c5.shape != (mu.size, 3, 3):
        raise ValueError("final_running.c5_abs has unexpected shape.")

    plt.figure(figsize=(7.6, 5.2))
    for i in range(3):
        for j in range(i, 3):
            values = c5[:, i, j]
            positive = values > 0.0
            if np.any(positive):
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
    plt.title("Weinberg-coefficient running")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save(output_path)


def plot_mnu_matrix_element_running(payload: dict[str, Any], output_path: Path) -> None:
    data = payload["final_running"]
    mu = np.asarray(data["mu_gev"], dtype=float)
    mnu = np.asarray(data["mnu_charged_lepton_abs_ev"], dtype=float)
    if mnu.shape != (mu.size, 3, 3):
        raise ValueError("final_running.mnu_charged_lepton_abs_ev has unexpected shape.")

    plt.figure(figsize=(7.6, 5.2))
    for i in range(3):
        for j in range(i, 3):
            values = mnu[:, i, j]
            positive = values > 0.0
            if np.any(positive):
                plt.plot(
                    mu[positive],
                    values[positive],
                    label=rf"$|(m_\nu)_{{{i+1}{j+1}}}|$",
                    linewidth=1.0,
                )
    plt.xscale("log")
    if np.any(mnu > 0.0):
        plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$|(m_\nu)_{ij}|$ [eV]")
    plt.title("Running neutrino-mass matrix elements")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save(output_path)


def _complex_matrix(block: dict[str, Any]) -> np.ndarray:
    return np.asarray(block["real"], dtype=float) + 1j * np.asarray(
        block["imag"], dtype=float
    )


def plot_low_energy_mnu_heatmap(payload: dict[str, Any], output_path: Path) -> None:
    matrix = np.abs(
        _complex_matrix(payload["low_energy"]["mnu_charged_lepton_basis_gev"])
    ) * 1.0e9

    plt.figure(figsize=(5.4, 4.8))
    image = plt.imshow(matrix)
    plt.colorbar(image, label="Magnitude [eV]")
    plt.xticks(range(3), [r"$e$", r"$\mu$", r"$\tau$"])
    plt.yticks(range(3), [r"$e$", r"$\mu$", r"$\tau$"])
    plt.title(r"$|m_\nu|$ in the charged-lepton basis")
    for i in range(3):
        for j in range(3):
            plt.text(j, i, f"{matrix[i, j]:.4f}", ha="center", va="center")
    _save(output_path)


def plot_pmns_heatmap(payload: dict[str, Any], output_path: Path) -> None:
    values = np.asarray(payload["low_energy"]["pmns_abs"], dtype=float)

    plt.figure(figsize=(5.4, 4.8))
    image = plt.imshow(values, vmin=0.0, vmax=1.0)
    plt.colorbar(image, label=r"$|U_{\rm PMNS}|$")
    plt.xticks(range(3), ["1", "2", "3"])
    plt.yticks(range(3), [r"$e$", r"$\mu$", r"$\tau$"])
    plt.title(r"$|U_{\rm PMNS}|$")
    for i in range(3):
        for j in range(3):
            plt.text(j, i, f"{values[i, j]:.4f}", ha="center", va="center")
    _save(output_path)


def write_readme(payload: dict[str, Any], output_path: Path) -> None:
    scales = payload["scales_gev"]
    text = f"""T3 numerical result figures

Default figures
---------------
uv_coupling_running.png
    Running of the UV quantities directly entering radiative neutrino-mass
    generation: ||y1||_F, ||y2||_F and |lambdaT3|, from
    mu_UV={scales['mu_uv']:.6g} GeV to
    mu_F={scales['mu_fermion_threshold']:.6g} GeV.

intermediate_direct_weinberg_running.png
    Accumulated direct one-loop LLSS -> Weinberg contribution between the
    fermion and scalar thresholds. This is the intermediate-EFT contribution
    that enters the authoritative one-loop final C5.

c5_threshold_contributions.png
    Component-by-component comparison of the physical hard threshold,
    direct intermediate-running contribution, and their final combined C5
    at the grouped scalar threshold.

c5_running.png
    Independent |C5^ij| entries in the final SM+Weinberg EFT below
    mu_S={scales['mu_scalar_threshold']:.6g} GeV.

neutrino_mass_splitting_running.png
    Solar and atmospheric neutrino mass-squared splittings below mu_S.

neutrino_mixing_running.png
    sin^2(theta12), sin^2(theta13) and sin^2(theta23) below mu_S.

oscillation_observable_sensitivity.png
    Controlled T1/T2/T3 UV-scaling trajectories in three oscillation-observable
    planes. Shaded 1sigma/3sigma regions are Gaussianized regions derived from
    the configured target covariance, not official Delta-chi2 contours.

Not generated by default
------------------------
intermediate_scalar_coupling_running.png
    Generic scalar quartics are useful diagnostics, but this plot does not yet
    isolate the intermediate dimension-five operator that drives the neutrino
    calculation.

neutrino_mass_running.png
    Absolute mass eigenvalues are secondary to the measured mass-squared
    splittings for the present oscillation-focused analysis.

mnu_matrix_element_running.png
    Largely duplicates the information in C5 running through
    m_nu = -(v^2/2) C5, apart from the charged-lepton basis rotation.

mnu_matrix_abs_ev.png
pmns_abs.png
    Static heatmaps are available as diagnostic helpers but are redundant with
    the running mixing observables and sensitivity plots for the main results.

Power-counting note
-------------------
The default intermediate plot shows the direct O(hbar) LLSS -> Weinberg
contribution that actually enters the authoritative one-loop C5. LLSS
self-running itself is O(hbar); inserting that correction into the scalar loop
would first contribute to C5 at O(hbar^2), so it is not promoted to a main
one-loop result figure.
"""
    output_path.write_text(text, encoding="utf-8")

def _target_index(target: dict[str, Any], name: str) -> int:
    names = list(target["observable_names"])
    if name not in names:
        raise KeyError(f"Target does not contain observable {name!r}.")
    return names.index(name)


def _add_gaussian_ellipse(
    ax,
    target: dict[str, Any],
    x_name: str,
    y_name: str,
    *,
    nsigma: float,
    alpha: float,
    label: str | None = None,
) -> None:
    central = np.asarray(target["central_values"], dtype=float)
    covariance = np.asarray(target["covariance"], dtype=float)

    ix = _target_index(target, x_name)
    iy = _target_index(target, y_name)

    cov2 = covariance[np.ix_([ix, iy], [ix, iy])]
    values, vectors = np.linalg.eigh(cov2)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]

    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    width, height = 2.0 * float(nsigma) * np.sqrt(values)

    ax.add_patch(
        Ellipse(
            (central[ix], central[iy]),
            width=width,
            height=height,
            angle=angle,
            fill=True,
            alpha=alpha,
            label=label,
        )
    )


def plot_observable_sensitivity(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    """Plot T3 UV-scaling sensitivity in three oscillation-observable planes."""

    sensitivity = payload.get("sensitivity")
    if not isinstance(sensitivity, dict) or sensitivity.get("status") != "Success":
        return

    target = sensitivity["target"]
    panels = (
        (
            "sin2_theta12",
            "delta_m21_sq_ev2",
            r"$\sin^2\theta_{12}$",
            r"$\Delta m_{21}^2\ [{\rm eV}^2]$",
            "Solar sector",
        ),
        (
            "sin2_theta13",
            "delta_m3l_sq_ev2",
            r"$\sin^2\theta_{13}$",
            r"$\Delta m_{3\ell}^2\ [{\rm eV}^2]$",
            "Atmospheric sector: reactor angle",
        ),
        (
            "sin2_theta23",
            "delta_m3l_sq_ev2",
            r"$\sin^2\theta_{23}$",
            r"$\Delta m_{3\ell}^2\ [{\rm eV}^2]$",
            "Atmospheric sector: atmospheric angle",
        ),
    )

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.5))

    for ax, (x_name, y_name, xlabel, ylabel, title) in zip(axes, panels):
        _add_gaussian_ellipse(
            ax,
            target,
            x_name,
            y_name,
            nsigma=3.0,
            alpha=0.13,
            label=r"Gaussian $3\sigma$",
        )
        _add_gaussian_ellipse(
            ax,
            target,
            x_name,
            y_name,
            nsigma=1.0,
            alpha=0.28,
            label=r"Gaussian $1\sigma$",
        )

        for direction in ("T1", "T2", "T3"):
            points = sensitivity["directions"][direction]["points"]
            x = np.asarray(
                [point["prediction"][x_name] for point in points],
                dtype=float,
            )
            y = np.asarray(
                [point["prediction"][y_name] for point in points],
                dtype=float,
            )
            gamma = np.asarray(
                [point["gamma"] for point in points],
                dtype=float,
            )

            ax.plot(
                x,
                y,
                marker="o",
                markersize=3.0,
                linewidth=1.2,
                label=direction,
            )

            index = int(np.argmin(np.abs(gamma - 1.0)))
            ax.scatter(
                [x[index]],
                [y[index]],
                marker="s",
                s=30,
                zorder=5,
            )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.22)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=5,
        frameon=True,
    )
    fig.suptitle(
        "Sensitivity of low-energy neutrino observables to T3 UV scaling directions",
        y=1.04,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_intermediate_direct_weinberg_running(
    payload: dict[str, Any],
    output_path: Path,
) -> None:
    """Plot the direct LLSS -> Weinberg contribution across the intermediate EFT."""

    running = payload.get("intermediate_direct_weinberg")
    if not isinstance(running, dict):
        return

    mu = np.asarray(running["mu_gev"], dtype=float)
    c5_abs = np.asarray(running["delta_c5_abs"], dtype=float)

    if c5_abs.shape != (mu.size, 3, 3):
        raise ValueError(
            "intermediate_direct_weinberg.delta_c5_abs has unexpected shape."
        )

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
                label=rf"$|\Delta C_5^{{{i+1}{j+1}}}|$",
                linewidth=1.15,
            )

    plt.xscale("log")

    positive_values = c5_abs[c5_abs > 0.0]
    if positive_values.size:
        plt.yscale("log")

    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$|\Delta C_5^{ij}|$ [GeV$^{-1}$]")
    plt.title("Direct Weinberg generation in the intermediate EFT")
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=8, ncol=2)
    _save(output_path)


def generate_running_result_figures(diagnostics_json: Path, output_dir: Path) -> Path:
    payload = _load_json(diagnostics_json)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Focused default result set. Additional plotting helpers remain available
    # for ad-hoc diagnostics but are not emitted automatically by --numerical.
    plot_uv_coupling_running(
        payload,
        output_dir / "uv_coupling_running.png",
    )
    plot_intermediate_direct_weinberg_running(
        payload,
        output_dir / "intermediate_direct_weinberg_running.png",
    )
    plot_c5_threshold_contributions(
        payload,
        output_dir / "c5_threshold_contributions.png",
    )
    plot_c5_running(
        payload,
        output_dir / "c5_running.png",
    )
    plot_mass_splitting_running(
        payload,
        output_dir / "neutrino_mass_splitting_running.png",
    )
    plot_mixing_angle_running(
        payload,
        output_dir / "neutrino_mixing_running.png",
    )
    plot_observable_sensitivity(
        payload,
        output_dir / "oscillation_observable_sensitivity.png",
    )
    write_readme(payload, output_dir / "README.txt")

    print(f"T3 running figures written to: {output_dir}")
    return output_dir


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate T3 running plots from T3 running diagnostics."
    )
    parser.add_argument("diagnostics_json", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    generate_running_result_figures(args.diagnostics_json, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
