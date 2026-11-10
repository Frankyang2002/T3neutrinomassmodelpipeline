"""Numerical general-scalar EFT1 endpoint -> SMEFT C5 -> neutrino observables.

Install at scripts/RunT3GeneralSMEFTBridge.py and run from the project root.
This consumes an ALREADY-MATCHED C5(MS); it does not calculate the F or S
threshold Wilson coefficients. No unrestricted/general matching is claimed.

Example:
 python -m scripts.RunT3GeneralSMEFTBridge --trajectory \
   output/studies/general_scalar_extended/continuous_full_scalar/NAME/A_m4/continuous_trajectory.json \
   --c5-ms path/to/general_matched_c5_ms.npz --c5-key C5_SMEFT_MS \
   --output output/studies/general_scalar_extended/bridged/A_m4
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np


def _load_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _complex(obj: Any, key: str) -> np.ndarray:
    if not isinstance(obj, dict) or "real" not in obj or "imag" not in obj:
        raise ValueError(f"{key}: expected {{real,imag}} numerical data")
    value = np.asarray(obj["real"], dtype=float) + 1j * np.asarray(obj["imag"], dtype=float)
    if not np.isfinite(value).all():
        raise ValueError(f"{key}: nonfinite value")
    return value


def _real(obj: Any, key: str) -> float:
    if isinstance(obj, dict):
        z = _complex(obj, key)
        if z.ndim or abs(z.imag) > 1e-10 * max(1.0, abs(float(z.real))):
            raise ValueError(f"{key}: expected real scalar")
        return float(z.real)
    value = float(obj)
    if not np.isfinite(value):
        raise ValueError(f"{key}: nonfinite scalar")
    return value


def _matrix(obj: Any, key: str) -> np.ndarray:
    a = _complex(obj, key)
    if a.shape != (3, 3):
        raise ValueError(f"{key}: expected 3x3 matrix, got {a.shape}")
    return a


def _load_c5(path: Path, key: str) -> np.ndarray:
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as data:
            if key not in data.files:
                raise ValueError(f"{path}: missing coefficient key {key!r}; available {data.files}")
            c5 = np.asarray(data[key], dtype=complex)
    elif path.suffix.lower() == ".npy":
        c5 = np.asarray(np.load(path, allow_pickle=False), dtype=complex)
    else:
        raise ValueError("C5 input must be .npz or .npy in GeV^-1")
    if c5.shape != (3, 3) or not np.isfinite(c5).all():
        raise ValueError("C5(MS) must be finite complex 3x3")
    scale = max(float(np.max(np.abs(c5))), 1e-50)
    if float(np.max(np.abs(c5 - c5.T))) > 1e-9 * scale:
        raise ValueError("C5(MS) must be complex symmetric, not Hermitian")
    return (c5 + c5.T) / 2


def _array_json(a: np.ndarray) -> dict:
    a = np.asarray(a)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def run(*, trajectory: Path, c5_ms: Path, output: Path, c5_key: str,
        mu_low: float, vev: float) -> dict:
    # Import only after parsing files, allowing this module to be syntax-checked
    # independently of the rest of the project.
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions, evolve_sm_weinberg,
    )
    from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables

    trajectory, c5_ms, output = Path(trajectory), Path(c5_ms), Path(output)
    if not trajectory.is_file() or not c5_ms.is_file():
        raise FileNotFoundError("Both the general trajectory and matched-C5 files are required")
    t = _load_json(trajectory)
    if t.get("format") != "t3-continuous-general-scalar-candidate-v1":
        raise ValueError("Input is not an existing general scalar continuous trajectory")
    if t.get("status") != "CandidateContinuousUVFermionThresholdEFT1Integrated":
        raise ValueError("General UV/EFT1 scalar evolution did not finish")
    model = t.get("model")
    endpoint = t.get("eft1_scalar_endpoint_state")
    sm = endpoint.get("sm") if isinstance(endpoint, dict) else None
    if not isinstance(sm, dict):
        raise ValueError("No saved SM state in the general EFT1 endpoint")
    settings = t.get("settings", {})
    ms = float(settings.get("scalar_endpoint_gev", float("nan")))
    if not (np.isfinite(ms) and 0 < mu_low <= ms and 0 < vev < 1000):
        raise ValueError("Require 0 < mu_low <= MS and physical positive vev")
    if abs(_real(endpoint.get("mu_gev"), "endpoint.mu_gev") - ms) > 1e-8 * ms:
        raise ValueError("General endpoint scale differs from configured MS")
    if not t.get("scalar_endpoint_physically_admissible", False):
        raise ValueError("General scalar endpoint has an inadmissible mass spectrum")
    c5 = _load_c5(c5_ms, c5_key)
    initial = SMWeinbergInitialConditions(
        gY=_real(sm["gY"], "gY"), g2=_real(sm["g2"], "g2"),
        g3=_real(sm["g3"], "g3"), lambdaH=_real(sm["lambdaH"], "lambdaH"),
        ye=_matrix(sm["ye"], "ye"), yu=_matrix(sm["yu"], "yu"),
        yd=_matrix(sm["yd"], "yd"), K=c5,
    ).validated()
    result = evolve_sm_weinberg(initial, ms, float(mu_low))
    if not result.solver_success:
        raise RuntimeError(f"SMEFT RGE integration failed: {result.solver_message}")
    c5_low = np.asarray(result.K, dtype=complex)
    mass = -0.5 * vev * vev * c5_low
    observables, matrices = neutrino_observables(mass, result.ye)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite: {output}")
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "smeft_low.npz", C5_SMEFT_MS=c5,
                        C5_SMEFT_low=c5_low, m_nu_GeV=mass,
                        **matrices)
    summary = {
        "model": model, "status": "GeneralScalarEndpointToSMEFTCompleted",
        "mu_s_gev": ms, "mu_low_gev": float(mu_low), "vev_gev": float(vev),
        "smeft_solver_nfev": int(result.nfev), "observables": observables,
        "c5_at_ms": _array_json(c5), "c5_at_low": _array_json(c5_low),
        "m_nu_gev": _array_json(mass),
        "input_sha256": {"general_scalar_trajectory": _sha(trajectory),
                         "matched_c5_ms": _sha(c5_ms)},
        "included": ["saved general scalar UV/EFT1 endpoint", "SMEFT one-loop running",
                     "Majorana masses and PMNS extraction"],
        "not_computed_here": ["fermion threshold Wilson matching",
                              "EFT1 Wilson RGE", "scalar threshold Wilson matching"],
        "c5_origin": "external_precomputed_matched_coefficient_required",
        "c5_input_key": c5_key,
        "warning": "The C5 source must independently match the SAME general trajectory and scalar threshold; this bridge cannot establish that relation.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--c5-ms", type=Path, required=True,
                        help="Previously matched 3x3 C5(MS) in GeV^-1, NPZ or NPY")
    parser.add_argument("--c5-key", default="C5_SMEFT_MS")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mu-low", type=float, default=100.0)
    parser.add_argument("--vev", type=float, default=246.22)
    args = parser.parse_args()
    result = run(trajectory=args.trajectory, c5_ms=args.c5_ms,
                 output=args.output, c5_key=args.c5_key,
                 mu_low=args.mu_low, vev=args.vev)
    print(json.dumps({"model": result["model"], "status": result["status"],
                      "observables": result["observables"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
