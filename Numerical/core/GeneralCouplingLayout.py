"""General T3 coupling layout for the existing numerical state/vector pipeline.

This is a registry/validation layer, NOT a second integrator or matching path.
It reads the native UV/EFT1 RGBeta `report_betas` inventories and declares
how the extra parameters must be represented in the existing numerical stack.

No coupling is silently removed.  Couplings without a known numerical type
are rejected; the downstream numerical runner must explicitly support them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

SM_MATRIX = frozenset(("yu", "yd", "ye"))
BSM_MATRIX = frozenset(("y1", "y2", "MF", "y1Conj", "y2Conj"))
REAL_BASE = frozenset((
    "gY", "g2", "g3", "lambdaH", "mS1Sq", "mS2Sq",
    "lambdaS1", "lambdaS2", "lambdaH1", "lambdaH2", "lambda12",
    "lambdaH1Adj", "lambdaH2Adj", "lambdaS1Adj", "lambdaS2Adj",
    "lambda12Adj", "lambda12Cross",
))
COMPLEX_BASE = frozenset((
    "lambdaT3", "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar",
    "lambdaS1bar2S2bar2", "lambdaS1barS2S2bar2",
    "lambdaS1S1bar2S2bar", "lambdaHHdagS1barS2barCross",
    "muT3GeneralS1Sq", "muT3GeneralS2Sq", "muT3GeneralS1S2",
))


@dataclass(frozen=True)
class GeneralCouplingSpec:
    name: str
    kind: str  # real, complex, matrix
    shape: tuple[int, ...] = ()

    def validate(self, value: Any, *, majorana: bool = False) -> float | complex | np.ndarray:
        """Return a validated independent numerical parameter; never discard Im parts."""
        if self.kind == "real":
            arr = np.asarray(value)
            if arr.shape != () or not np.isfinite(arr).all() or abs(complex(arr).imag) > 0:
                raise ValueError(f"{self.name}: expected finite real scalar")
            return float(arr.real)
        if self.kind == "complex":
            arr = np.asarray(value, dtype=complex)
            if arr.shape != () or not np.isfinite(arr).all():
                raise ValueError(f"{self.name}: expected finite complex scalar")
            return complex(arr)
        if self.kind == "matrix":
            arr = np.asarray(value, dtype=complex)
            if arr.shape != self.shape or not np.isfinite(arr).all():
                raise ValueError(f"{self.name}: expected finite matrix {self.shape}")
            if majorana and self.name == "MF" and not np.allclose(
                arr, arr.T, rtol=1e-10, atol=1e-14
            ):
                raise ValueError("MF must be symmetric for Majorana fermions")
            return arr.copy()
        raise ValueError(f"Unknown type for {self.name}: {self.kind}")


def spec_for(name: str) -> GeneralCouplingSpec:
    """Map actual RGBeta names to numerical storage types.

    Extra `lambdaT3GeneralQ...` tensors are registered by the source generator
    as complex non-self-conjugate families. The native export's real
    self-conjugate structural quartics keep their existing real types.
    """
    if name in SM_MATRIX or name in BSM_MATRIX:
        return GeneralCouplingSpec(name, "matrix", (3, 3))
    if name in REAL_BASE:
        return GeneralCouplingSpec(name, "real")
    if name in COMPLEX_BASE or name.startswith("lambdaT3GeneralQ"):
        return GeneralCouplingSpec(name, "complex")
    raise ValueError(f"Unrecognised general T3 RGBeta parameter: {name}")


def layout_from_rgbeta(payload: Mapping[str, Any]) -> tuple[GeneralCouplingSpec, ...]:
    """Use the native export as the *only* inventory of running couplings."""
    betas = payload.get("report_betas")
    if not isinstance(betas, Mapping) or not betas:
        raise ValueError("RGBeta payload missing non-empty report_betas")
    names = list(betas)
    if not all(isinstance(n, str) and n for n in names):
        raise ValueError("RGBeta beta names must be non-empty strings")
    if not all(isinstance(b, str) and b.strip() for b in betas.values()):
        raise ValueError("RGBeta beta expressions must be non-empty strings")
    return tuple(spec_for(n) for n in names)


def validate_general_values(
    payload: Mapping[str, Any],
    values: Mapping[str, Any],
    *,
    majorana: bool = False,
) -> dict[str, float | complex | np.ndarray]:
    """Require exact UV/EFT1 state coverage before an existing ODE run."""
    specs = layout_from_rgbeta(payload)
    expected = {s.name for s in specs}
    actual = set(values)
    if actual != expected:
        raise ValueError(
            f"General T3 state mismatch; missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )
    return {s.name: s.validate(values[s.name], majorana=majorana) for s in specs}


def fermion_first_projection(
    eft1_payload: Mapping[str, Any],
    uv_values: Mapping[str, Any],
    *,
    overrides: Mapping[str, Any] | None = None,
) -> dict[str, float | complex | np.ndarray]:
    """Select exact EFT1 renormalizable parameters; no matching implied.

    Matching corrections are not constructed here. Reject an unknown EFT1
    coupling rather than setting it to zero without a physical boundary.
    """
    source = dict(uv_values)
    if overrides:
        source.update(overrides)
    names = [s.name for s in layout_from_rgbeta(eft1_payload)]
    missing = set(names) - set(source)
    if missing:
        raise ValueError(f"EFT1 threshold values missing: {sorted(missing)}")
    return validate_general_values(eft1_payload, {n: source[n] for n in names})
