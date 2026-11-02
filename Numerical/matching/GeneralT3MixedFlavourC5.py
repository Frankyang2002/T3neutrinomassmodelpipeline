"""Provisional one-scale mass-eigenstate C5 for ordinary T3 B(-1)/C(-1).

Uses the EXISTING saved Matchete ordered coefficient and its existing checked
InputForm compiler, RGBeta running parameters and mixed-doublet diagonaliser.
The Stage23/26 conditional map Matchete quartic = -RGBeta quartic is a
DECLARED ASSUMPTION, not an established physical field-basis dictionary.

Nothing in this module implements fermion-first EFT1 Wilson matching, EFT1
Wilson anomalous dimensions, scalar-threshold subtraction/counterterms or
full sequential EFT matching. Outputs must not be labelled physical predictions.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / "output/general_t3_mass_eigen_matchete/one_loop_probe"
STEMS = {"B_-1": "B_m1_final", "C_-1": "C_m1_final"}
DIRECT_RGBETA_TO_MAT = {
    "lambdaHHdagS2S2": "lambdaT3ExtraQ020020I1",
    "lambdaHHdagS1barS1bar": "lambdaT3ExtraQ020200I1",
}
SIGNED_MAP = -1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError("General B/C mixed C5: " + message)


def _complex(value: Any, name: str) -> complex:
    z = complex(value)
    require(math.isfinite(z.real) and math.isfinite(z.imag), name + " is nonfinite")
    return z


def _mp(z: Any) -> mp.mpc:
    z = _complex(z, "numeric coefficient")
    return mp.mpc(str(z.real), str(z.imag))


def _source_files(model: str, directory: Path) -> tuple[Path, Path]:
    require(model in STEMS, "only B_-1 and C_-1 are supported")
    stem = STEMS[model]
    return (directory / (stem + "_one_loop_probe_weinberg_ordered_coefficient.inputform.txt"),
            directory / (stem + "_one_loop_probe.json"))


def validate_source(model: str, coefficient: Path, probe: Path) -> dict[str, Any]:
    """Check archived source metadata; do not promote its certificate to a physical map."""
    require(coefficient.is_file() and coefficient.stat().st_size > 0,
            "missing archived ordered Matchete C5: " + str(coefficient))
    require(probe.is_file(), "missing matching provenance: " + str(probe))
    report = json.loads(probe.read_text(encoding="utf-8-sig"))
    expected = {
        "model": model,
        "format": "t3-matchete-native-heavy-mass-eigenbasis-v1",
        "status": "OneLoopDirectSymbolicWeinbergProbePassedNumericalC5Pending",
        "one_loop_match_stage_completed": True,
        "weinberg_extraction_status": "Success",
        "weinberg_holomorphic_term_count": 118,
        "ordered_coefficient_is_canonical_flavour_symmetric_c5": False,
        "physical_matching_performed": False,
        "production_sequential_matching_performed": False,
        "component_basis_map_certified": False,
        "matchete_rgbeta_quartic_normalization_certified": False,
    }
    for key, value in expected.items():
        require(report.get(key) == value and type(report.get(key)) is type(value),
                f"archived probe mismatch for {key}: {report.get(key)!r}")
    require(report.get("weinberg_ordered_coefficient_inputform_file") == coefficient.name,
            "probe does not identify the selected coefficient file")
    require(report.get("mass_basis_ordinary_check_passed") is True and
            report.get("heavy_mass_basis_detailed_check_passed") is True,
            "archived mass-basis check not passed")
    return {"coefficient_sha256": sha256(coefficient.read_bytes()).hexdigest(),
            "probe_sha256": sha256(probe.read_bytes()).hexdigest()}


def map_direct_quartics(eft_state: Any) -> dict[str, complex]:
    """Stage23/26 CONJECTURAL signed map. No identity/renormalization claim."""
    values = {}
    for rgbeta, matchete in DIRECT_RGBETA_TO_MAT.items():
        raw = getattr(eft_state, rgbeta, None)
        require(raw is not None, "EFT1 lacks running quartic " + rgbeta)
        values[matchete] = SIGNED_MAP * _complex(raw, rgbeta)
    return values


def evaluate_ordered_coefficient(code: Any, masses: np.ndarray,
                                 y1: np.ndarray, y2: np.ndarray,
                                 scalar_basis: Any, lambda_t3: complex,
                                 quartics: dict[str, complex]) -> np.ndarray:
    """Sum archived Matchete A5_ij over heavy indices, then C5=A5+A5^T."""
    require(masses.shape == (3,) and np.all(np.isfinite(masses)) and np.all(masses > 0),
            "invalid three-generation heavy masses")
    require(y1.shape == (3, 3) and y2.shape == (3, 3),
            "Yukawas must both be complex 3x3")
    require(np.isfinite(y1).all() and np.isfinite(y2).all(), "nonfinite Yukawas")
    require(not scalar_basis.tachyonic and not scalar_basis.degenerate and
            np.all(np.asarray(scalar_basis.eigenvalues_sq) > 0),
            "nonpositive/degenerate scalar mass eigenvalues")
    require(set(quartics) == set(DIRECT_RGBETA_TO_MAT.values()),
            "missing/unexpected direct quartic keys")
    mp.mp.dps = 45
    U = np.asarray(scalar_basis.U, dtype=complex)
    require(U.shape == (2, 2), "mixed scalar rotation must be 2x2")
    env = {f"T3MEU{i+1}{j+1}": _mp(U[i, j])
           for i in range(2) for j in range(2)}
    env.update({"MSX": mp.mpf(str(math.sqrt(float(scalar_basis.eigenvalues_sq[0])))),
                "MSY": mp.mpf(str(math.sqrt(float(scalar_basis.eigenvalues_sq[1])))),
                "lambdaT3": _mp(lambda_t3),
                "hbar": mp.mpf(1) / (16 * mp.pi**2),
                "conj": mp.conj, "sqrt": mp.sqrt, "log": mp.log})
    env.update({name: _mp(value) for name, value in quartics.items()})
    ordered = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            total = mp.mpc(0)
            for k in range(3):
                local = dict(env, MF=mp.mpf(str(float(masses[k]))),
                             y1a=_mp(y1[i, k]), y1b=_mp(y1[j, k]),
                             y2a=_mp(y2[i, k]), y2b=_mp(y2[j, k]))
                val = eval(code, {"__builtins__": {}}, local)
                require(bool(mp.isfinite(mp.re(val)) and mp.isfinite(mp.im(val))),
                        "nonfinite ordered Matchete kernel")
                total += val
            ordered[i, j] = complex(total)
    result = ordered + ordered.T  # NOT ordered + ordered.conj().T
    require(np.isfinite(result.real).all() and np.isfinite(result.imag).all(),
            "nonfinite final complex C5")
    return result


@dataclass
class GeneralT3MixedFlavourC5:
    model: str
    source_dir: Path = DEFAULT_SOURCE
    hbar_value: float = 1 / (16 * math.pi**2)
    coefficient_path: Path = field(init=False)
    probe_path: Path = field(init=False)
    provenance: dict[str, Any] = field(init=False)
    last_diagnostics: dict[str, Any] = field(init=False, default_factory=dict)

    def __post_init__(self) -> None:
        require(math.isfinite(self.hbar_value) and self.hbar_value > 0,
                "hbar_value must be finite positive")
        require(math.isclose(self.hbar_value, 1 / (16 * math.pi**2),
                             rel_tol=1e-14, abs_tol=0),
                "B/C archived source uses the fixed one-loop factor 1/(16*pi^2)")
        directory = Path(self.source_dir)
        self.coefficient_path, self.probe_path = _source_files(self.model, directory)
        self.provenance = validate_source(self.model, self.coefficient_path, self.probe_path)

    def evaluate(self, uv_state: Any, eft_state: Any) -> np.ndarray:
        from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data
        from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
        from scripts.RunGeneralT3NumericalC5 import ALLOWED_COUPLINGS, compile_coefficient

        require(uv_state.representation.majorana_fermion,
                "B/C(-1) requires the neutral Majorana heavy mass basis")
        require(eft_state.representation.alpha == -1,
                "EFT1 representation does not match B/C(-1)")
        mu_sq = eft_state.general_couplings.get("muT3GeneralS1S2")
        require(mu_sq is not None, "missing running general mixed bilinear")
        basis = MixedDoubletMassBasis.build(
            self.model, float(eft_state.mS1Sq), float(eft_state.mS2Sq),
            _complex(mu_sq, "muT3GeneralS1S2"),
        )
        basis.verify()
        require(not basis.tachyonic and not basis.degenerate and
                np.min(basis.eigenvalues_sq) > 0, "invalid scalar eigenmass spectrum")
        heavy = _heavy_mass_basis_data(
            uv_state, offdiagonal_atol=1e-7, imaginary_atol=1e-9,
            takagi_residual_rtol=1e-8, takagi_residual_atol=1e-7,
        )
        source = self.coefficient_path.read_text(encoding="utf-8-sig")
        symbols = set(re.findall(r"Coupling\[([A-Za-z][A-Za-z0-9]*)", source))
        require(symbols == ALLOWED_COUPLINGS,
                "archived Matchete coupling inventory does not match checked evaluator: "
                + repr(sorted(symbols ^ ALLOWED_COUPLINGS)))
        compiled = compile_coefficient(source)
        quartics = map_direct_quartics(eft_state)
        result = evaluate_ordered_coefficient(
            compiled, np.asarray(heavy.masses_gev, dtype=float),
            np.asarray(heavy.y1, dtype=complex), np.asarray(heavy.y2, dtype=complex),
            basis, _complex(eft_state.lambdaT3, "lambdaT3"), quartics,
        )
        self.last_diagnostics = {
            "matching_kind": "B/C mixed-doublet mass-eigenstate-Matchete-one-scale-conditional-signed-map",
            "model": self.model,
            "archived_source": str(self.coefficient_path),
            "source_provenance": self.provenance,
            "heavy_basis_residual": float(max(heavy.takagi_residual, heavy.basis_residual)),
            "heavy_masses_at_F_gev": np.asarray(heavy.masses_gev, dtype=float).tolist(),
            "mixed_scalar_basis": basis.info(),
            "quartic_basis_mapping_assumptions": [
                f"UNCERTIFIED {mat} = -{rgb} (Stage23/26 conditional sign)"
                for rgb, mat in DIRECT_RGBETA_TO_MAT.items()],
            "matchete_rgbeta_quartic_map_certified": False,
            "matchete_rgbeta_mixed_bilinear_map_certified": False,
            "complex_mixed_scalar_bilinear_supported_in_candidate_basis": True,
            "sequential_eft_wilson_matching_performed": False,
            "eft1_wilson_running_performed": False,
            "physical_general_c5_certified": False,
        }
        return result
