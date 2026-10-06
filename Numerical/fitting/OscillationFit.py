"""Neutrino-oscillation fit layer for T3 numerical scans.

This module deliberately separates

    model prediction
        -> oscillation observable vector
        -> experimental target/covariance
        -> chi^2

from any particular NuFIT release.

No oscillation-fit numbers are hard-coded here.  A concrete NuFIT release can
be supplied as JSON once its central values, uncertainties/covariance and mass
ordering convention have been fixed and documented.

The canonical fitted observable vector is

    delta_m21_sq_ev2
    delta_m3l_sq_ev2
    sin2_theta12
    sin2_theta13
    sin2_theta23

with

    delta_m3l^2 = delta_m31^2  for normal ordering,
    delta_m3l^2 = delta_m32^2  for inverted ordering.

The angle extraction uses |U_PMNS|:

    sin^2(theta13) = |U_e3|^2
    sin^2(theta12) = |U_e2|^2 / (1 - |U_e3|^2)
    sin^2(theta23) = |U_mu3|^2 / (1 - |U_e3|^2).
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from physics.NeutrinoObservables import (
    calculate_neutrino_observables,
)


OSCILLATION_OBSERVABLES: tuple[str, ...] = (
    "delta_m21_sq_ev2",
    "delta_m3l_sq_ev2",
    "sin2_theta12",
    "sin2_theta13",
    "sin2_theta23",
)


@dataclass(frozen=True)
class OscillationPrediction:
    """Canonical low-energy oscillation prediction."""

    ordering: str
    delta_m21_sq_ev2: float
    delta_m3l_sq_ev2: float
    sin2_theta12: float
    sin2_theta13: float
    sin2_theta23: float
    takagi_residual: float

    def as_dict(self) -> dict[str, float]:
        return {
            "delta_m21_sq_ev2": self.delta_m21_sq_ev2,
            "delta_m3l_sq_ev2": self.delta_m3l_sq_ev2,
            "sin2_theta12": self.sin2_theta12,
            "sin2_theta13": self.sin2_theta13,
            "sin2_theta23": self.sin2_theta23,
        }

    def vector(
        self,
        names: Sequence[str] = OSCILLATION_OBSERVABLES,
    ) -> np.ndarray:
        values = self.as_dict()

        unknown = [name for name in names if name not in values]
        if unknown:
            raise ValueError(
                "Unknown oscillation observable(s): "
                + ", ".join(unknown)
            )

        return np.asarray(
            [values[name] for name in names],
            dtype=float,
        )


@dataclass(frozen=True)
class OscillationFitTarget:
    """Experimental central values and covariance for one ordering."""

    ordering: str
    observable_names: tuple[str, ...]
    central_values: np.ndarray
    covariance: np.ndarray
    source: str = ""

    def validated(self) -> "OscillationFitTarget":
        ordering = self.ordering.upper()
        if ordering not in {"NO", "IO"}:
            raise ValueError("ordering must be 'NO' or 'IO'.")

        names = tuple(self.observable_names)

        if not names:
            raise ValueError(
                "At least one oscillation observable is required."
            )

        if len(set(names)) != len(names):
            raise ValueError(
                "observable_names cannot contain duplicates."
            )

        unknown = [
            name
            for name in names
            if name not in OSCILLATION_OBSERVABLES
        ]
        if unknown:
            raise ValueError(
                "Unsupported oscillation observable(s): "
                + ", ".join(unknown)
            )

        central = np.asarray(
            self.central_values,
            dtype=float,
        )
        covariance = np.asarray(
            self.covariance,
            dtype=float,
        )

        n = len(names)

        if central.shape != (n,):
            raise ValueError(
                f"central_values must have shape ({n},)."
            )

        if covariance.shape != (n, n):
            raise ValueError(
                f"covariance must have shape ({n}, {n})."
            )

        if not np.all(np.isfinite(central)):
            raise ValueError(
                "central_values must contain finite entries."
            )

        if not np.all(np.isfinite(covariance)):
            raise ValueError(
                "covariance must contain finite entries."
            )

        if not np.allclose(
            covariance,
            covariance.T,
            rtol=1.0e-12,
            atol=1.0e-18,
        ):
            raise ValueError("covariance must be symmetric.")

        eigenvalues = np.linalg.eigvalsh(covariance)

        if np.any(eigenvalues <= 0.0):
            raise ValueError(
                "covariance must be positive definite."
            )

        return OscillationFitTarget(
            ordering=ordering,
            observable_names=names,
            central_values=central.copy(),
            covariance=covariance.copy(),
            source=str(self.source),
        )

    @classmethod
    def from_independent_errors(
        cls,
        *,
        ordering: str,
        central_values: Mapping[str, float],
        one_sigma_errors: Mapping[str, float],
        source: str = "",
    ) -> "OscillationFitTarget":
        """Build a diagonal-covariance target from independent 1-sigma errors."""

        names = tuple(central_values)

        if set(names) != set(one_sigma_errors):
            raise ValueError(
                "central_values and one_sigma_errors must use the same keys."
            )

        central = np.asarray(
            [central_values[name] for name in names],
            dtype=float,
        )
        sigma = np.asarray(
            [one_sigma_errors[name] for name in names],
            dtype=float,
        )

        if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0.0):
            raise ValueError(
                "All one-sigma errors must be finite and positive."
            )

        return cls(
            ordering=ordering,
            observable_names=names,
            central_values=central,
            covariance=np.diag(sigma**2),
            source=source,
        ).validated()

    @classmethod
    def from_json(
        cls,
        path: Path,
    ) -> "OscillationFitTarget":
        """Load a documented oscillation target from JSON."""

        path = Path(path)
        payload = json.loads(
            path.read_text(encoding="utf-8")
        )

        if not isinstance(payload, dict):
            raise ValueError(
                "Oscillation target JSON must contain an object."
            )

        ordering = str(payload["ordering"])
        source = str(payload.get("source", ""))

        if "covariance" in payload:
            names = tuple(payload["observable_names"])
            central = np.asarray(
                payload["central_values"],
                dtype=float,
            )
            covariance = np.asarray(
                payload["covariance"],
                dtype=float,
            )

            return cls(
                ordering=ordering,
                observable_names=names,
                central_values=central,
                covariance=covariance,
                source=source,
            ).validated()

        central_mapping = payload.get("central_values")
        error_mapping = payload.get("one_sigma_errors")

        if not isinstance(central_mapping, dict):
            raise ValueError(
                "JSON central_values must be an object when covariance is absent."
            )

        if not isinstance(error_mapping, dict):
            raise ValueError(
                "JSON one_sigma_errors must be an object when covariance is absent."
            )

        return cls.from_independent_errors(
            ordering=ordering,
            central_values={
                str(name): float(value)
                for name, value in central_mapping.items()
            },
            one_sigma_errors={
                str(name): float(value)
                for name, value in error_mapping.items()
            },
            source=source,
        )


@dataclass(frozen=True)
class OscillationFitResult:
    """Gaussian chi-square result."""

    ordering: str
    chi2: float
    observable_names: tuple[str, ...]
    prediction: np.ndarray
    target: np.ndarray
    residual: np.ndarray
    whitened_residual: np.ndarray
    source: str


def prediction_from_mass_matrix(
    mass_matrix_gev: np.ndarray,
    *,
    ordering: str = "AUTO",
) -> OscillationPrediction:
    """Convert a Majorana mass matrix into the canonical fit observables."""

    observables = calculate_neutrino_observables(
        mass_matrix_gev,
        ordering=ordering,
    )

    pmns_abs = np.asarray(
        observables.pmns_abs,
        dtype=float,
    )

    if pmns_abs.shape != (3, 3):
        raise RuntimeError(
            "PMNS magnitude matrix must be 3x3."
        )

    sin2_theta13 = float(pmns_abs[0, 2] ** 2)
    denominator = 1.0 - sin2_theta13

    if denominator <= 0.0:
        raise RuntimeError(
            "Cannot extract theta12/theta23 because 1-|Ue3|^2 <= 0."
        )

    sin2_theta12 = float(
        pmns_abs[0, 1] ** 2 / denominator
    )
    sin2_theta23 = float(
        pmns_abs[1, 2] ** 2 / denominator
    )

    ordering_label = observables.ordering.upper()

    delta_m3l = (
        observables.delta_m31_sq_ev2
        if ordering_label == "NO"
        else observables.delta_m32_sq_ev2
    )

    result = OscillationPrediction(
        ordering=ordering_label,
        delta_m21_sq_ev2=float(
            observables.delta_m21_sq_ev2
        ),
        delta_m3l_sq_ev2=float(delta_m3l),
        sin2_theta12=sin2_theta12,
        sin2_theta13=sin2_theta13,
        sin2_theta23=sin2_theta23,
        takagi_residual=float(
            observables.takagi_residual
        ),
    )

    values = result.vector()

    if not np.all(np.isfinite(values)):
        raise RuntimeError(
            "Oscillation prediction contains non-finite values."
        )

    return result


def evaluate_oscillation_fit(
    prediction: OscillationPrediction,
    target: OscillationFitTarget,
) -> OscillationFitResult:
    """Evaluate the multivariate Gaussian chi-square."""

    target = target.validated()

    if prediction.ordering != target.ordering:
        raise ValueError(
            "Prediction and target mass orderings differ: "
            f"{prediction.ordering} vs {target.ordering}."
        )

    predicted = prediction.vector(
        target.observable_names
    )
    residual = predicted - target.central_values

    # Cholesky avoids explicitly forming covariance^{-1}.
    chol = np.linalg.cholesky(target.covariance)
    whitened = np.linalg.solve(chol, residual)

    chi2 = float(
        np.dot(whitened, whitened)
    )

    return OscillationFitResult(
        ordering=target.ordering,
        chi2=chi2,
        observable_names=target.observable_names,
        prediction=predicted,
        target=target.central_values.copy(),
        residual=residual,
        whitened_residual=whitened,
        source=target.source,
    )


def evaluate_mass_matrix_fit(
    mass_matrix_gev: np.ndarray,
    target: OscillationFitTarget,
) -> OscillationFitResult:
    """Convenience wrapper from m_nu directly to chi-square."""

    target = target.validated()

    prediction = prediction_from_mass_matrix(
        mass_matrix_gev,
        ordering=target.ordering,
    )

    return evaluate_oscillation_fit(
        prediction,
        target,
    )
