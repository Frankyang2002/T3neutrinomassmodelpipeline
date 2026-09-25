"""Generic effective-field-theory configuration objects.

This module contains small value objects used by the pipeline to describe
what an EFT contains and how one threshold changes that content.  It must not
contain T3-specific matching formulae, beta functions, or assumptions about
which heavy field is integrated out first.

The central distinction is:

- :class:`EFTContent` describes *which heavy fields are active* in one EFT.
- :class:`ThresholdStep` describes *which fields are removed* and at what scale.
- :class:`EFTTransition` joins those two descriptions for one physical matching
  step.
- :class:`EFTRunningInterval` describes the EFT that runs between two adjacent
  matching thresholds.

Detailed couplings and Wilson coefficients remain in the numerical/matching
layers.  These objects are deliberately lightweight so ``pipeline.py`` can use
them as a readable description of the calculation order.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, TypeAlias


ThresholdScale: TypeAlias = str | float


@dataclass(frozen=True, slots=True)
class EFTTruncation:
    """Operator-dimension truncation used when constructing an EFT.

    The production T3 pipeline currently works through mass dimension five.
    Keeping the truncation as explicit run configuration avoids encoding that
    approximation indirectly through a preferred threshold ordering.
    """

    max_operator_dimension: int = 5

    def __post_init__(self) -> None:
        if isinstance(self.max_operator_dimension, bool) or not isinstance(
            self.max_operator_dimension,
            int,
        ):
            raise TypeError("max_operator_dimension must be an integer.")
        if self.max_operator_dimension < 1:
            raise ValueError("max_operator_dimension must be positive.")

    def keeps(self, operator_dimension: int) -> bool:
        """Return whether an operator is retained by this EFT truncation."""
        if isinstance(operator_dimension, bool) or not isinstance(
            operator_dimension,
            int,
        ):
            raise TypeError("operator_dimension must be an integer.")
        if operator_dimension < 1:
            raise ValueError("operator_dimension must be positive.")
        return operator_dimension <= self.max_operator_dimension

    @property
    def label(self) -> str:
        """Return a compact human-readable truncation label."""
        return f"d<={self.max_operator_dimension}"


DEFAULT_EFT_TRUNCATION = EFTTruncation(max_operator_dimension=5)


def _canonical_field_names(fields: Iterable[str]) -> frozenset[str]:
    """Validate physical field names and return an order-independent set."""
    values = tuple(fields)

    if any(not isinstance(field, str) or not field.strip() for field in values):
        raise ValueError("EFT field names must be non-empty strings.")

    stripped = tuple(field.strip() for field in values)
    if len(set(stripped)) != len(stripped):
        raise ValueError("An EFT field cannot be listed more than once.")

    return frozenset(stripped)


@dataclass(frozen=True, slots=True)
class EFTContent:
    """Field content and truncation defining one EFT between thresholds.

    ``active_heavy_fields`` uses *physical* field names and is intentionally
    order-independent.  For ordinary T3 these are drawn from ``F``, ``S1`` and
    ``S2``; shared-scalar mode uses ``F`` and ``S``.  This class itself remains
    model-independent and does not hard-code those names.

    The object does not store numerical couplings.  Its purpose is to let the
    pipeline and dispatch code identify a theory by its actual physical content
    rather than by labels such as ``EFT1`` or ``EFT2``.
    """

    active_heavy_fields: frozenset[str]
    truncation: EFTTruncation = DEFAULT_EFT_TRUNCATION

    def __post_init__(self) -> None:
        canonical = _canonical_field_names(self.active_heavy_fields)
        object.__setattr__(self, "active_heavy_fields", canonical)

        if not isinstance(self.truncation, EFTTruncation):
            raise TypeError("truncation must be an EFTTruncation instance.")

    @classmethod
    def from_fields(
        cls,
        fields: Iterable[str],
        *,
        truncation: EFTTruncation = DEFAULT_EFT_TRUNCATION,
    ) -> "EFTContent":
        """Construct EFT content from any iterable of active field names."""
        return cls(
            active_heavy_fields=_canonical_field_names(fields),
            truncation=truncation,
        )

    @property
    def is_fully_decoupled(self) -> bool:
        """Return whether no heavy fields remain active."""
        return not self.active_heavy_fields

    def has_exactly(self, *fields: str) -> bool:
        """Return whether exactly the requested physical fields are active."""
        return self.active_heavy_fields == _canonical_field_names(fields)

    def after_integrating(self, fields: Iterable[str]) -> "EFTContent":
        """Return the EFT content after removing ``fields``.

        Removing a field that is not active is an invalid threshold transition
        and is rejected here before any matching calculation is attempted.
        """
        to_remove = _canonical_field_names(fields)
        inactive = to_remove - self.active_heavy_fields

        if inactive:
            raise ValueError(
                "Cannot integrate out field(s) that are not active: "
                + ", ".join(sorted(inactive))
            )

        return EFTContent(
            active_heavy_fields=self.active_heavy_fields - to_remove,
            truncation=self.truncation,
        )

    @property
    def heavy_field_label(self) -> str:
        """Return a compact human-readable active-field label."""
        if self.is_fully_decoupled:
            return "none"
        return ", ".join(sorted(self.active_heavy_fields))


@dataclass(frozen=True, slots=True)
class ThresholdStep:
    """One physical threshold in an ordered EFT sequence.

    ``fields_to_integrate`` contains physical field names.  ``scale`` may be a
    numerical value in GeV or a symbolic scale name such as ``MF`` or ``MS1``.
    """

    fields_to_integrate: tuple[str, ...]
    scale: ThresholdScale

    def __post_init__(self) -> None:
        if not self.fields_to_integrate:
            raise ValueError("A threshold step must integrate at least one field.")

        canonical = tuple(field.strip() for field in self.fields_to_integrate)
        _canonical_field_names(canonical)
        object.__setattr__(self, "fields_to_integrate", canonical)

        if isinstance(self.scale, str):
            if not self.scale.strip():
                raise ValueError("A symbolic threshold scale cannot be empty.")
        elif isinstance(self.scale, bool) or not isinstance(self.scale, (int, float)):
            raise TypeError("Threshold scale must be a symbolic string or a number.")
        elif self.scale <= 0:
            raise ValueError("A numerical threshold scale must be positive.")

    @property
    def label(self) -> str:
        """Return the field-group label used in threshold-plan displays."""
        if len(self.fields_to_integrate) == 1:
            return self.fields_to_integrate[0]
        return "(" + ",".join(self.fields_to_integrate) + ")"


@dataclass(frozen=True, slots=True)
class EFTTransition:
    """One run/match boundary connecting two EFT field contents.

    ``before`` and ``after`` are derived from the same physical threshold step.
    The class validates that relationship so later dispatch code can trust the
    transition without re-encoding threshold-order assumptions.
    """

    index: int
    step: ThresholdStep
    before: EFTContent
    after: EFTContent

    def __post_init__(self) -> None:
        if isinstance(self.index, bool) or not isinstance(self.index, int):
            raise TypeError("Transition index must be an integer.")
        if self.index < 1:
            raise ValueError("Transition index must start at 1.")

        expected_after = self.before.after_integrating(
            self.step.fields_to_integrate
        )
        if expected_after != self.after:
            raise ValueError(
                "EFTTransition.after does not match the threshold fields removed "
                "from EFTTransition.before."
            )

        if self.before.truncation != self.after.truncation:
            raise ValueError(
                "One threshold transition must use a consistent EFT truncation."
            )

    @property
    def fields_to_integrate(self) -> tuple[str, ...]:
        """Return the physical fields removed at this threshold."""
        return self.step.fields_to_integrate

    @property
    def scale(self) -> ThresholdScale:
        """Return the matching scale for this transition."""
        return self.step.scale


@dataclass(frozen=True, slots=True)
class EFTRunningInterval:
    """One EFT running region between two adjacent heavy thresholds.

    The interval is identified by its physical field content, not by an ordinal
    label such as ``EFT1``.  ``entered_by`` is the threshold that creates the
    EFT and ``exited_by`` is the next threshold that removes additional heavy
    fields.

    The fully decoupled low-energy theory is intentionally not represented by
    this class because its final running is handled separately from the heavy
    threshold sequence.
    """

    index: int
    content: EFTContent
    high_scale: ThresholdScale
    low_scale: ThresholdScale
    entered_by: EFTTransition
    exited_by: EFTTransition

    def __post_init__(self) -> None:
        if isinstance(self.index, bool) or not isinstance(self.index, int):
            raise TypeError("Running-interval index must be an integer.")
        if self.index < 1:
            raise ValueError("Running-interval index must start at 1.")

        if self.exited_by.index != self.entered_by.index + 1:
            raise ValueError(
                "An EFT running interval must be bounded by adjacent threshold "
                "transitions."
            )

        if self.index != self.entered_by.index:
            raise ValueError(
                "Running-interval index must match the threshold that creates it."
            )

        if self.content != self.entered_by.after:
            raise ValueError(
                "Running-interval content must equal the EFT produced by its "
                "entry threshold."
            )

        if self.content != self.exited_by.before:
            raise ValueError(
                "Running-interval content must equal the EFT entering its next "
                "threshold."
            )

        if self.high_scale != self.entered_by.scale:
            raise ValueError(
                "Running-interval high scale must equal its entry matching scale."
            )

        if self.low_scale != self.exited_by.scale:
            raise ValueError(
                "Running-interval low scale must equal its exit matching scale."
            )

        if self.content.is_fully_decoupled:
            raise ValueError(
                "A heavy-field running interval must contain at least one active "
                "heavy field."
            )

    @property
    def active_heavy_fields(self) -> frozenset[str]:
        """Return the physical heavy fields active throughout the interval."""
        return self.content.active_heavy_fields

    def has_exactly(self, *fields: str) -> bool:
        """Return whether exactly ``fields`` are active in this running region."""
        return self.content.has_exactly(*fields)

