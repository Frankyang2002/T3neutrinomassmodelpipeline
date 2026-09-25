"""Physical heavy-field identities used by the T3 EFT pipeline.

The T3 topology is written formally with scalar roles ``S1`` and ``S2``.  In
ordinary models those roles are distinct physical fields.  In the shared-scalar
(scotogenic) branch they denote one physical scalar ``S`` related by conjugation.

This module is the single place that translates between those two descriptions.
Threshold planning, stage metadata and running backends should work with the
physical field names returned here; only the Wolfram matching boundary should
expand a physical shared scalar back to the formal roles ``S1`` and ``S2``.
"""

from __future__ import annotations

from dataclasses import dataclass


FERMION_FIELD = "F"
FORMAL_SCALAR_FIELDS: tuple[str, str] = ("S1", "S2")
ORDINARY_SCALAR_FIELDS: tuple[str, str] = FORMAL_SCALAR_FIELDS
SHARED_SCALAR_FIELDS: tuple[str, ...] = ("S",)


@dataclass(frozen=True, slots=True)
class T3HeavyFieldScheme:
    """Physical heavy-field naming for one T3 model.

    ``shared_scalar=False`` gives the ordinary physical fields ``F,S1,S2``.
    ``shared_scalar=True`` gives ``F,S`` while preserving the formal Wolfram
    roles ``S1,S2`` behind :meth:`formal_roles_for`.
    """

    shared_scalar: bool = False

    @property
    def scalar_fields(self) -> tuple[str, ...]:
        """Return the physical T3 scalar fields."""
        return SHARED_SCALAR_FIELDS if self.shared_scalar else ORDINARY_SCALAR_FIELDS

    @property
    def heavy_fields(self) -> tuple[str, ...]:
        """Return all physical heavy fields in the established display order."""
        return (FERMION_FIELD, *self.scalar_fields)

    @property
    def scalar_field_set(self) -> frozenset[str]:
        """Return the physical scalar sector as an order-independent set."""
        return frozenset(self.scalar_fields)

    @property
    def heavy_field_set(self) -> frozenset[str]:
        """Return all physical heavy fields as an order-independent set."""
        return frozenset(self.heavy_fields)

    @property
    def threshold_aliases(self) -> dict[str, str]:
        """Return accepted canonicalised CLI aliases for threshold fields."""
        if self.shared_scalar:
            return {
                "F": "F",
                "FERMION": "F",
                "S": "S",
                "SCALAR": "S",
                "S1": "S",
                "S2": "S",
                "SCALAR1": "S",
                "SCALAR2": "S",
            }

        return {
            "F": "F",
            "FERMION": "F",
            "S1": "S1",
            "SCALAR1": "S1",
            "S2": "S2",
            "SCALAR2": "S2",
        }

    @property
    def threshold_field_help(self) -> str:
        """Return the established user-facing list of allowed field names."""
        return "F or S" if self.shared_scalar else "F, S1, or S2"

    def normalise_threshold_field(self, value: str) -> str:
        """Convert one user-facing threshold token to a physical field name."""
        token = value.strip().upper().replace("_", "")
        try:
            return self.threshold_aliases[token]
        except KeyError as exc:
            raise ValueError(
                f"Unknown heavy field {value!r}. Use {self.threshold_field_help}."
            ) from exc

    def formal_roles_for(self, physical_field: str) -> tuple[str, ...]:
        """Map one physical field to the formal T3 roles used by matching.

        Ordinary T3 fields map one-to-one.  The shared physical scalar ``S``
        expands to both formal topology roles because the matching code still
        writes the generic T3 Lagrangian in terms of ``S1`` and ``S2``.
        """
        if physical_field not in self.heavy_field_set:
            raise ValueError(
                f"{physical_field!r} is not a physical heavy field in this T3 scheme."
            )
        if self.shared_scalar and physical_field == "S":
            return FORMAL_SCALAR_FIELDS
        return (physical_field,)


ORDINARY_T3_FIELDS = T3HeavyFieldScheme(shared_scalar=False)
SHARED_SCALAR_T3_FIELDS = T3HeavyFieldScheme(shared_scalar=True)


def t3_heavy_field_scheme(*, shared_scalar: bool = False) -> T3HeavyFieldScheme:
    """Return the immutable physical heavy-field scheme for one model branch."""
    return SHARED_SCALAR_T3_FIELDS if shared_scalar else ORDINARY_T3_FIELDS
