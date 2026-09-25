"""Select the T3 model points requested by one pipeline study.

This module answers only *which* UV models should be calculated.  It keeps the
T3-A...E scan definitions, neutral-state filtering, and direct ``--dims``
selection separate from the Lagrangian/matching implementation.

The selected models are represented by :class:`T3ModelRequest`.  The actual UV
construction and EFT matching are performed by ``Lagrangian.T3ModelMatching``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Literal

from common.T3Model import SMOKE, T3_CLASSES, t3_has_neutral_bsm_component


# Comparison-study definitions kept exactly as in the historical pipeline.
DIMENSION_COMPARISON = tuple(
    (model_class, 0)
    for model_class in ("A", "B", "C", "D", "E")
)

HYPERCHARGE_ALPHAS = tuple(range(-4, 3))
HYPERCHARGE_COMPARISON = tuple(
    (model_class, alpha)
    for model_class in ("A", "B", "C", "D", "E")
    for alpha in HYPERCHARGE_ALPHAS
)


RequestKind = Literal["class", "dimensions", "shared_dimensions"]


@dataclass(frozen=True, slots=True)
class T3ModelRequest:
    """One requested T3 UV model before any Lagrangian calculation is run.

    ``kind='class'`` uses the published T3-A...E representation catalogue.
    ``kind='dimensions'`` stores explicit ``(dS1,dS2,dF)`` dimensions.
    ``kind='shared_dimensions'`` stores the one-physical-scalar ``(dS,dF)``
    branch used for the scotogenic limit.

    The request deliberately contains no output paths, matching artifacts, or
    RGE state.  It is the boundary between study selection and model building.
    """

    kind: RequestKind
    alpha: int
    dimensions: tuple[int, ...] = ()
    model_class: str | None = None

    def __post_init__(self) -> None:
        if self.kind == "class":
            if self.model_class not in T3_CLASSES:
                raise ValueError(f"Unknown T3 model class: {self.model_class}")
            if self.dimensions:
                raise ValueError("Class requests must not provide explicit dimensions.")
            return

        if self.model_class is not None:
            raise ValueError("Explicit-dimension requests must not provide model_class.")

        expected = 3 if self.kind == "dimensions" else 2
        if len(self.dimensions) != expected:
            raise ValueError(
                f"{self.kind} requests require exactly {expected} dimensions."
            )

        if self.kind == "shared_dimensions" and self.alpha != -1:
            raise ValueError("Shared-scalar model requests require alpha=-1.")

    @classmethod
    def class_point(cls, model_class: str, alpha: int) -> "T3ModelRequest":
        """Construct one T3-A...E scan point."""
        return cls(kind="class", model_class=model_class, alpha=alpha)

    @classmethod
    def explicit_dimensions(
        cls,
        dimensions: tuple[int, int, int],
        alpha: int,
    ) -> "T3ModelRequest":
        """Construct one ordinary explicit-dimension T3 point."""
        return cls(kind="dimensions", dimensions=dimensions, alpha=alpha)

    @classmethod
    def shared_dimensions(
        cls,
        dimensions: tuple[int, int],
    ) -> "T3ModelRequest":
        """Construct one physical shared-scalar T3 point."""
        return cls(kind="shared_dimensions", dimensions=dimensions, alpha=-1)


def _neutral_scan_points(
    points: tuple[tuple[str, int], ...],
) -> tuple[tuple[str, int], ...]:
    """Keep scan points containing at least one electrically neutral BSM state."""
    selected: list[tuple[str, int]] = []
    for model_class, alpha in points:
        d_s1, d_s2, d_f = T3_CLASSES[model_class]
        if t3_has_neutral_bsm_component(d_s1, d_s2, d_f, alpha):
            selected.append((model_class, alpha))
    return tuple(selected)


def _scan_definition(
    args: argparse.Namespace,
) -> tuple[str, tuple[tuple[str, int], ...]]:
    """Return the label and benchmark points for a non-``--dims`` scan."""

    if args.smoke:
        label, candidates = "smoke", tuple(SMOKE)
    elif args.hypercharge_comparison:
        label, candidates = "hypercharge comparison", HYPERCHARGE_COMPARISON
    elif args.dimension_comparison:
        label, candidates = "dimension comparison", DIMENSION_COMPARISON
    else:
        return

    if args.force:
        return label, candidates

    return label, _neutral_scan_points(candidates)


def select_study_models(
    args: argparse.Namespace,
    *,
    shared_scalar_mode: bool,
) -> list[T3ModelRequest]:
    """Return the UV model requests selected by one CLI study.

    This function performs no Wolfram calculation and creates no output files.
    It preserves the historical scan definitions and terminal scan summary.
    """

    if args.dims:
        if shared_scalar_mode:
            d_s, d_f = args.dims
            return [T3ModelRequest.shared_dimensions((d_s, d_f))]

        d_s1, d_s2, d_f = args.dims
        return [
            T3ModelRequest.explicit_dimensions(
                (d_s1, d_s2, d_f),
                args.alpha,
            )
        ]

    mode_name, points = _scan_definition(args)
    candidate_count = (
        len(HYPERCHARGE_COMPARISON)
        if args.hypercharge_comparison
        else len(DIMENSION_COMPARISON)
        if args.dimension_comparison
        else len(SMOKE)
    )
    if args.force:
        print(f"T3 scan mode: {mode_name}; force enabled; {len(points)} model(s).")
    else:
        print(
            f"T3 scan mode: {mode_name}; {len(points)}/{candidate_count} "
            "neutral-compatible model(s)."
        )

    return [
        T3ModelRequest.class_point(model_class, alpha)
        for model_class, alpha in points
    ]
