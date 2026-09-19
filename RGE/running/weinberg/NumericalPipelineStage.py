from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import sympy as sp

from RGE.matching.FlavorMatchedC5 import (
    extract_t3_loop_kernel,
    is_final_weinberg_json,
    load_final_weinberg_flavor_matrix,
)
from RGE.matching.MatchedEFTRGE import parse_matchete_c5
from RGE.running.weinberg.WeinbergRunning import (
    SMInitialConditions,
    evolve_weinberg,
    neutrino_mass_matrix,
)


def _complex_matrix_from_json(real_part, imag_part) -> np.ndarray:
    real = np.asarray(real_part, dtype=float)
    imag = np.asarray(imag_part, dtype=float)

    if real.shape != imag.shape:
        raise ValueError("Real and imaginary matrix parts must have the same shape.")

    return real + 1j * imag


def _complex_matrix_to_lists(matrix: np.ndarray) -> dict:
    matrix = np.asarray(matrix, dtype=complex)

    return {
        "real": matrix.real.tolist(),
        "imag": matrix.imag.tolist(),
    }


def _sympy_substitutions_from_config(config: dict) -> dict:
    substitutions: dict[sp.Expr, complex | float] = {}

    model = config["t3"]

    scalar_values = {
        "MS1": model["MS1"],
        "MS2": model["MS2"],
        "lambdaT3": model["lambdaT3"],
        "hbar": model.get("hbar", 1.0 / (16.0 * np.pi**2)),
    }

    for name, value in scalar_values.items():
        substitutions[sp.Symbol(name)] = value

    heavy_masses = model["MF"]

    for index, value in enumerate(heavy_masses, start=1):
        substitutions[sp.Symbol(f"MF{index}")] = value

    y1_real = np.asarray(model["y1_real"], dtype=float)
    y1_imag = np.asarray(model["y1_imag"], dtype=float)
    y2_real = np.asarray(model["y2_real"], dtype=float)
    y2_imag = np.asarray(model["y2_imag"], dtype=float)

    if y1_real.shape != y1_imag.shape:
        raise ValueError("y1_real and y1_imag must have the same shape.")

    if y2_real.shape != y2_imag.shape:
        raise ValueError("y2_real and y2_imag must have the same shape.")

    if y1_real.shape != y2_real.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    if y1_real.shape[0] != 3:
        raise ValueError("T3 Yukawa matrices must have three lepton-flavor rows.")

    if y1_real.shape[1] != len(heavy_masses):
        raise ValueError(
            "The number of Yukawa columns must match the number of heavy masses."
        )

    y1 = y1_real + 1j * y1_imag
    y2 = y2_real + 1j * y2_imag

    for p in range(y1.shape[0]):
        for r in range(y1.shape[1]):
            symbol1 = sp.Symbol(f"y1_{p + 1}{r + 1}")
            symbol2 = sp.Symbol(f"y2_{p + 1}{r + 1}")

            substitutions[symbol1] = y1[p, r]
            substitutions[sp.conjugate(symbol1)] = np.conjugate(y1[p, r])

            substitutions[symbol2] = y2[p, r]
            substitutions[sp.conjugate(symbol2)] = np.conjugate(y2[p, r])

    return substitutions



def _hierarchical_scale_substitutions(config: dict) -> dict[sp.Expr, float]:
    """Numerical values for the grouped F -> (S1,S2) threshold symbols."""
    model = config["t3"]

    if "MS" in model:
        ms_value = float(model["MS"])
    else:
        ms1 = float(model["MS1"])
        ms2 = float(model["MS2"])

        if not np.isclose(ms1, ms2):
            raise ValueError(
                "Hierarchical final C5 uses the grouped scalar threshold MS. "
                "Provide t3.MS explicitly, or use equal MS1 and MS2."
            )

        ms_value = ms1

    return {
        sp.Symbol("MS"): ms_value,
    }


def evaluate_final_weinberg_json(
    c5_path: Path,
    config: dict,
) -> np.ndarray:
    """Evaluate the physical Majorana C5 matrix stored in the final JSON."""
    model = config["t3"]
    heavy_masses = model["MF"]

    symbolic = load_final_weinberg_flavor_matrix(
        c5_path,
        n_lepton=3,
        n_heavy=len(heavy_masses),
        split_heavy_masses=True,
    )
    K_symbolic = symbolic["K"]

    substitutions = _sympy_substitutions_from_config(config)
    substitutions.update(_hierarchical_scale_substitutions(config))

    K_numeric = np.empty((3, 3), dtype=complex)

    for p in range(3):
        for q in range(3):
            value = sp.N(K_symbolic[p, q].subs(substitutions), 18)

            if value.free_symbols:
                raise ValueError(
                    "Unresolved symbols remain in hierarchical C5 entry "
                    f"({p},{q}): {sorted(map(str, value.free_symbols))}"
                )

            K_numeric[p, q] = complex(value)

    if not np.allclose(K_numeric, K_numeric.T):
        raise RuntimeError(
            "Numerically evaluated final Weinberg matrix is not symmetric."
        )

    return K_numeric


def evaluate_symbolic_c5(
    c5_path: Path,
    config: dict,
) -> np.ndarray:
    """Evaluate matched flavor C5 without building a symbolic 3x3 matrix."""

    c5_path = Path(c5_path)

    kappa_1g = parse_matchete_c5(
        c5_path.read_text(encoding="utf-8")
    )
    kernel = extract_t3_loop_kernel(kappa_1g)

    model = config["t3"]
    heavy_masses = np.asarray(model["MF"], dtype=float)

    y1 = _complex_matrix_from_json(
        model["y1_real"],
        model["y1_imag"],
    )
    y2 = _complex_matrix_from_json(
        model["y2_real"],
        model["y2_imag"],
    )

    if y1.shape != y2.shape:
        raise ValueError("y1 and y2 must have the same shape.")

    if y1.shape[0] != 3:
        raise ValueError(
            "T3 Yukawa matrices must have three lepton-flavor rows."
        )

    if y1.shape[1] != len(heavy_masses):
        raise ValueError(
            "The number of Yukawa columns must match the number of heavy masses."
        )

    base_substitutions = {
        sp.Symbol("MS1"): model["MS1"],
        sp.Symbol("MS2"): model["MS2"],
        sp.Symbol("lambdaT3"): model["lambdaT3"],
        sp.Symbol("hbar"): model.get(
            "hbar",
            1.0 / (16.0 * np.pi**2),
        ),
    }

    MF = sp.Symbol("MF")
    kernel_values = np.empty(len(heavy_masses), dtype=complex)

    for r, mass in enumerate(heavy_masses):
        kernel_r = kernel.subs(
            {
                **base_substitutions,
                MF: float(mass),
            }
        )
        kernel_values[r] = complex(sp.N(kernel_r, 18))

    K = np.zeros((3, 3), dtype=complex)

    for p in range(3):
        for q in range(p, 3):
            value = 0.0j

            for r in range(len(heavy_masses)):
                value += 0.5 * kernel_values[r] * (
                    np.conjugate(y1[p, r])
                    * np.conjugate(y2[q, r])
                    + np.conjugate(y2[p, r])
                    * np.conjugate(y1[q, r])
                )

            K[p, q] = value
            K[q, p] = value

    return K


def run_numerical_pipeline_stage(
    c5_path: Path,
    output_dir: Path,
    config_path: Path,
) -> dict:
    """Evaluate matched C5 numerically and run it to the requested low scale."""

    c5_path = Path(c5_path)
    output_dir = Path(output_dir)
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    config_path = Path(config_path)

    config = json.loads(
        config_path.read_text(encoding="utf-8")
    )

    if is_final_weinberg_json(c5_path):
        K_initial = evaluate_final_weinberg_json(
            c5_path,
            config,
        )
        c5_input_kind = "final_weinberg_json"
    else:
        K_initial = evaluate_symbolic_c5(
            c5_path,
            config,
        )
        c5_input_kind = "legacy_scalar_c5"

    sm = config["sm"]

    initial = SMInitialConditions(
        gY=sm["gY"],
        g2=sm["g2"],
        g3=sm["g3"],
        lambdaH=sm["lambdaH"],
        ye=np.asarray(sm["ye"], dtype=float),
        yu=np.asarray(sm["yu"], dtype=float),
        yd=np.asarray(sm["yd"], dtype=float),
        K=K_initial,
    )

    result = evolve_weinberg(
        initial,
        config["mu_initial_gev"],
        config["mu_final_gev"],
    )

    mass_low = neutrino_mass_matrix(
        result.K,
        vev_gev=config.get("vev_gev", 246.22),
    )

    initial_path = data_dir / "c5_flavor_matrix_numeric.txt"
    low_path = data_dir / "c5_flavor_matrix_low_scale.txt"
    mass_path = data_dir / "neutrino_mass_matrix_low_scale.txt"

    initial_path.write_text(
        json.dumps(
            _complex_matrix_to_lists(K_initial),
            indent=2,
        ),
        encoding="utf-8",
    )

    low_path.write_text(
        json.dumps(
            _complex_matrix_to_lists(result.K),
            indent=2,
        ),
        encoding="utf-8",
    )

    mass_path.write_text(
        json.dumps(
            _complex_matrix_to_lists(mass_low),
            indent=2,
        ),
        encoding="utf-8",
    )

    summary = {
        "NumericalRGEStatus": "Success",
        "NumericalRGEConfig": str(config_path),
        "C5InputKind": c5_input_kind,
        "C5InputFile": c5_path.relative_to(output_dir).as_posix()
        if c5_path.is_relative_to(output_dir)
        else str(c5_path),
        "NumericalMatchingScaleGeV": result.mu_initial,
        "NumericalLowScaleGeV": result.mu_final,
        "NumericalSolverEvaluations": result.nfev,
        "C5FlavorMatrixNumericFile": initial_path.relative_to(output_dir).as_posix(),
        "C5FlavorMatrixLowScaleFile": low_path.relative_to(output_dir).as_posix(),
        "NeutrinoMassMatrixLowScaleFile": mass_path.relative_to(output_dir).as_posix(),
    }

    summary_path = data_dir / "numerical_rge_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )

    return summary
