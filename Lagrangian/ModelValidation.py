"""Pure T3 model validation and run-specification construction.

This module owns Python-side checks which decide whether a requested T3 model
is a valid/supported production point and how that point is named on disk.

It does not launch Wolfram, read matching outputs, or construct runtime
result records.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from common.T3Model import (
    T3_CLASSES,
    encode_alpha,
    identify_t3_class,
    shared_scalar_formal_dimensions,
    t3_has_neutral_bsm_component,
    valid_shared_scalar_dimensions,
    valid_shared_scalar_topology_dimensions,
    valid_t3_dimensions,
    valid_t3_topology_dimensions,
)


@dataclass(frozen=True, slots=True)
class ModelRunSpecification:
    """Validated Python-side description of one Wolfram model run."""

    name: str
    alpha: int
    d_s1: int
    d_s2: int
    d_f: int
    output_dir: Path
    model_args: tuple[str, ...]
    shared_scalar: bool = False


def ordinary_run_specification(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    output_root: Path,
    force: bool = False,
) -> ModelRunSpecification:
    """Validate and name one ordinary T3 model point."""

    # ``--force`` bypasses only production-scope and neutrality restrictions.
    # It never bypasses the representation-theory conditions that define T3.
    topology_ok = valid_t3_topology_dimensions(d_s1, d_s2, d_f)
    if not topology_ok:
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) does not form the required T3 topology. "
            "Dimensions must be positive, each scalar must satisfy dS=dF±1, "
            "and S1⊗S2 must contain the triplet."
        )

    if not force and not valid_t3_dimensions(d_s1, d_s2, d_f):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}) is outside the current production support. "
            "Normal mode supports only SU(2) dimensions 1, 2, and 3. "
            "Use --force to attempt a larger representation that still satisfies "
            "the T3 topology conditions."
        )

    if not force and not t3_has_neutral_bsm_component(d_s1, d_s2, d_f, alpha):
        raise ValueError(
            f"({d_s1}, {d_s2}, {d_f}), alpha={alpha} has no electrically neutral "
            "BSM component. Normal mode requires at least one neutral state. "
            "Use --force to run this charged-only point explicitly."
        )

    model_class = identify_t3_class(d_s1, d_s2, d_f)

    if model_class is not None:
        name = f"T3-{model_class}"
        output_dir = output_root / (
            f"T3_{model_class}_alpha_{encode_alpha(alpha)}"
        )
    else:
        name = f"T3-d{d_s1}-d{d_s2}-F{d_f}"
        output_dir = output_root / (
            f"T3_d{d_s1}_d{d_s2}_F{d_f}_alpha_{encode_alpha(alpha)}"
        )

    model_args = (
        "DIMS",
        str(d_s1),
        str(d_s2),
        str(d_f),
        encode_alpha(alpha),
    )

    return ModelRunSpecification(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        output_dir=output_dir,
        model_args=model_args,
        shared_scalar=False,
    )


def shared_run_specification(
    d_s: int,
    d_f: int,
    *,
    output_root: Path,
    force: bool = False,
) -> ModelRunSpecification:
    """Validate and name one physical shared-scalar T3 model point."""

    if not valid_shared_scalar_topology_dimensions(d_s, d_f):
        raise ValueError(
            f"({d_s}, {d_f}) does not form the shared-scalar T3 topology. "
            "Dimensions must be positive, dS=dF±1, and S⊗S must contain the triplet."
        )

    if not force and not valid_shared_scalar_dimensions(d_s, d_f):
        raise ValueError(
            f"({d_s}, {d_f}) is outside supported shared-scalar mode. "
            "Current production support is dS=2 with dF=1 or 3. "
            "Use --force to attempt a larger topology-compatible representation."
        )

    d_s1, d_s2, d_f_formal = shared_scalar_formal_dimensions(
        d_s,
        d_f,
        force=force,
    )
    alpha = -1
    fermion_label = "N" if d_f == 1 else f"F{d_f}"
    name = f"Scotogenic-dS{d_s}-{fermion_label}"
    output_dir = output_root / f"Scotogenic_dS{d_s}_F{d_f}"
    model_args = (
        "DIMS",
        str(d_s1),
        str(d_s2),
        str(d_f_formal),
        encode_alpha(alpha),
    )

    return ModelRunSpecification(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f_formal,
        output_dir=output_dir,
        model_args=model_args,
        shared_scalar=True,
    )


def dimensions_for_class(model_class: str) -> tuple[int, int, int]:
    """Return the formal dimensions of one known T3-A...E model."""

    if model_class not in T3_CLASSES:
        raise ValueError(f"Unknown T3 model class: {model_class}")

    return T3_CLASSES[model_class]
