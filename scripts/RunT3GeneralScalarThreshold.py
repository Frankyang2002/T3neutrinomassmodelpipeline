"""General scalar endpoint + canonical EFT1 Wilson tensor -> scalar threshold -> SMEFT.

Uses the repository's existing general scalar mass Hessian, scalar bubble
threshold and RunT3GeneralSMEFTBridge. Does NOT reconstruct a general
Matchete/RGBeta quartic tensor or the fermion-threshold Wilson coefficient.

Required inputs:
 - trajectory: RunT3FullScalarContinuous.py JSON
 - Wilson: full C_psi2phi2_ms tensor in the declared canonical Weyl45,
   interleaved real-scalar coordinates, already evolved from MF to MS
 - quartic: complete V4=lambda_abcd phi^a phi^b phi^c phi^d/4! at the SAME MS
   in those SAME coordinates. The scalar quartic is NOT inferred from scalar
   beta-coordinate names.
The user must explicitly assert the shared canonical frame, after checking
its origin. No phase/sign, field rotation or Wilson normalisation is guessed.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.threshold.matching import (
    match_scalar_threshold, scalar_gauge_ward_residuals,
)


def _json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _complex(value: object, name: str) -> complex:
    if not isinstance(value, dict) or not {"real", "imag"}.issubset(value):
        raise ValueError(f"{name}: expected explicit real/imag object")
    v = complex(float(value["real"]), float(value["imag"]))
    if not np.isfinite([v.real, v.imag]).all():
        raise ValueError(f"{name}: nonfinite")
    return v


def _load_array(path: Path, key: str) -> np.ndarray:
    if path.suffix.lower() == ".npy":
        value = np.load(path, allow_pickle=False)
    elif path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as z:
            if key not in z.files:
                raise ValueError(f"{path}: missing {key}; present={z.files}")
            value = z[key].copy()
    else:
        raise ValueError(f"Input requires .npy/.npz: {path}")
    if not np.isfinite(value).all():
        raise ValueError(f"{path}/{key}: nonfinite")
    return value


def grouped_to_interleaved(mass_grouped: np.ndarray) -> np.ndarray:
    """The general quadratic code orders [R1,...,Rn,I1,...,In]."""
    m = np.asarray(mass_grouped, dtype=float)
    if m.ndim != 2 or m.shape[0] != m.shape[1] or m.shape[0] % 2:
        raise ValueError("Grouped mass Hessian must be real 2n x 2n")
    n = m.shape[0] // 2
    indices = np.array([x for k in range(n) for x in (k, n+k)])
    return m[np.ix_(indices, indices)]


def _mass_from_general(trajectory: dict, model: str) -> tuple[np.ndarray, dict]:
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
        quadratic_mass_hessian, load_manifest,
    )
    from scripts.RunT3FullScalarContinuous import MANIFEST

    state = trajectory.get("eft1_scalar_endpoint_state")
    if not isinstance(state, dict):
        raise ValueError("General scalar endpoint state absent")
    for name in ("mS1Sq", "mS2Sq"):
        if name not in state:
            raise ValueError(f"General endpoint missing {name}")
    mu = _complex(
        trajectory.get("quadratic_mu2_at_scalar_endpoint_gev2"),
        "quadratic_mu2_at_scalar_endpoint_gev2",
    )
    m1, m2 = float(state["mS1Sq"]), float(state["mS2Sq"])
    record = quadratic_mass_hessian(load_manifest(MANIFEST), model, m1, m2, mu)
    m = grouped_to_interleaved(np.asarray(record["real_scalar_hessian_gev2"], float))
    if not np.allclose(m, m.T, rtol=0, atol=1e-9):
        raise ValueError("Real-component Hessian not symmetric")
    if np.linalg.eigvalsh(m)[0] <= 0:
        raise ValueError("Tachyonic or zero scalar eigenmass")
    return m, record


def _ward(mass: np.ndarray, quartic: np.ndarray, model: str) -> dict:
    basis = build_basis(model)
    g = basis.generators()
    gauge = np.stack([g["U1"][0][1], *(pair[1] for pair in g["SU2"])])
    return scalar_gauge_ward_residuals(quartic, mass, gauge)


def run(*, trajectory: Path, wilson: Path, quartic: Path,
        output: Path, wilson_key: str, quartic_key: str,
        assume_canonical: bool, allow_ward_residual: bool,
        allow_projection_residual: bool, mu_low: float, vev: float) -> dict:
    if not assume_canonical:
        raise ValueError(
            "Pass --assume-same-canonical-frame ONLY after confirming both "
            "input tensors are in the same Weyl45/interleaved-real basis"
        )
    trajectory, wilson, quartic, output = map(
        Path, (trajectory, wilson, quartic, output)
    )
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {output}")
    info = _json(trajectory)
    if (info.get("format") != "t3-continuous-general-scalar-candidate-v1"
            or info.get("status") != "CandidateContinuousUVFermionThresholdEFT1Integrated"):
        raise ValueError("Input is not a completed general UV -> EFT1 scalar trajectory")
    model = str(info["model"])
    basis = build_basis(model)
    ms = float(info["settings"]["scalar_endpoint_gev"])
    mf = float(info["settings"]["fermion_threshold_gev"])
    if not (mf > ms >= mu_low > 0):
        raise ValueError("Required hierarchy MF > MS >= mu_low > 0")
    if not info.get("scalar_endpoint_physically_admissible", False):
        raise ValueError("General scalar endpoint is tachyonic")
    m, mass_record = _mass_from_general(info, model)
    c = np.asarray(_load_array(wilson, wilson_key), complex)
    q = np.asarray(_load_array(quartic, quartic_key), float)
    if c.shape != (basis.nf, basis.nf, basis.ns, basis.ns):
        raise ValueError(f"Full Wilson shape should be {(basis.nf,basis.nf,basis.ns,basis.ns)}")
    if q.shape != (basis.ns,)*4:
        raise ValueError(f"Full quartic shape should be {(basis.ns,)*4}")
    check = _ward(m, q, model)
    if not check["passed"] and not allow_ward_residual:
        raise ValueError(
            "U1/SU2 scalar mass-quartic gauge Ward check failed: "
            + json.dumps(check)
            + ". Check real-component order and conjugation; "
            "do not silently force identity maps."
        )

    matched = match_scalar_threshold(
        full_wilson_ms=c,
        lepton_indices=np.array(basis.leptons, dtype=int),
        higgs_indices=np.array(basis.higgs, dtype=int),
        heavy_indices=np.arange(4, basis.ns, dtype=int),
        full_quartic_ms=q,
        scalar_mass2_background=m,
        R_wilson_from_background=np.eye(basis.ns-4),
        mu_s=ms, mode="rge_improved",
        enforce_pure_weinberg=not allow_projection_residual,
    )
    output.mkdir(parents=True, exist_ok=True)
    c5 = np.asarray(matched["c5_smeft_at_ms"], complex)
    c5_file = output / "general_scalar_threshold_c5.npz"
    np.savez_compressed(
        c5_file,
        C5_SMEFT_MS=c5,
        C5_EFT1_MS=matched["c5_eft1"],
        C5_scalar_bubble_MS=matched["c5_scalar_bubble"],
        mass2_real_interleaved_gev2=m,
    )
    # Reuse the previously supplied downstream bridge rather than duplicating
    # SM+Weinberg running and PMNS conventions.
    from scripts.RunT3GeneralSMEFTBridge import run as run_smeft_bridge
    downstream = run_smeft_bridge(
        trajectory=trajectory,
        c5_ms=c5_file,
        c5_key="C5_SMEFT_MS",
        output=output / "smeft",
        mu_low=mu_low,
        vev=vev,
    )
    result = {
        "model": model,
        "status": "GeneralScalarEndpointConditionalScalarThresholdAndSMEFT",
        "mu_f_gev": mf, "mu_s_gev": ms, "mu_low_gev": mu_low,
        "mass_convention": "grouped real/imag -> interleaved; S=(R+iI)/sqrt(2)",
        "wilson_and_quartic_frame_assumption":
            "explicit same canonical Weyl45/real-interleaved declaration; no inferred Matchete mapping",
        "gauge_ward": check,
        "mass_eigenvalues2_gev2": np.linalg.eigvalsh(m).tolist(),
        "general_quadratic_kind": mass_record["quadratic_kind"],
        "scalar_bubble_pairs": matched["bubble"]["pair_count"],
        "weinberg_residual_eft1": matched["rge_input_C5_residual"],
        "weinberg_residual_scalar_bubble": matched["scalar_bubble_projection_residual"],
        "inputs_sha256": {p.name: _sha(p) for p in (trajectory, wilson, quartic)},
        "c5_smeft_at_ms_file": c5_file.name,
        "downstream_summary": "smeft/summary.json",
        "observables": downstream["observables"],
        "not_implemented_here": [
            "generation of full general scalar quartic tensor from RGBeta coordinate names",
            "generation of first-threshold general Wilson coefficient",
            "running the full general EFT1 Wilson tensor",
            "additional finite scalar matching beyond the existing bubble",
        ],
        "warning": "Conditional numerical joining only: same-input provenance and native-to-canonical coefficient normalization must be established by the producer.",
    }
    (output / "manifest.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--trajectory", type=Path, required=True)
    p.add_argument("--wilson", type=Path, required=True)
    p.add_argument("--quartic", type=Path, required=True)
    p.add_argument("--wilson-key", default="C_psi2phi2_ms")
    p.add_argument("--quartic-key", default="quartic_real_full_ms")
    p.add_argument("--assume-same-canonical-frame", action="store_true")
    p.add_argument("--allow-ward-residual", action="store_true",
                   help="Diagnostic only: do not label noncovariant C5 physical")
    p.add_argument("--allow-projection-residual", action="store_true",
                   help="Diagnostic only: permits non-Weinberg projection")
    p.add_argument("--mu-low", type=float, default=100.0)
    p.add_argument("--vev", type=float, default=246.22)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    result = run(
        trajectory=a.trajectory, wilson=a.wilson, quartic=a.quartic,
        output=a.output, wilson_key=a.wilson_key, quartic_key=a.quartic_key,
        assume_canonical=a.assume_same_canonical_frame,
        allow_ward_residual=a.allow_ward_residual,
        allow_projection_residual=a.allow_projection_residual,
        mu_low=a.mu_low, vev=a.vev,
    )
    print(json.dumps({
        "model": result["model"], "status": result["status"],
        "observables": result["observables"]
    }, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
