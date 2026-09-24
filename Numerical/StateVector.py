"""Pack and unpack canonical T3 numerical states for ODE integration.

The ODE independent variable is t = ln(mu), so ``mu_gev`` is metadata rather
than an entry of the state vector.  Every complex matrix/scalar is represented
by real degrees of freedom:

    complex matrix A -> vec(Re A), vec(Im A)
    complex scalar z -> (Re z, Im z)

The block order is deterministic and follows the RGBeta coupling names so a
later beta-function evaluator can map each exported beta directly onto the
corresponding numerical block.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from Numerical.State import (
    SMNumericalState,
    SharedT3UVState,
    T3Representation,
    T3UVState,
)


BlockKind = Literal["real_scalar", "complex_scalar", "complex_matrix"]


@dataclass(frozen=True)
class VectorBlock:
    """One contiguous block in the real ODE vector."""

    name: str
    kind: BlockKind
    start: int
    stop: int
    shape: tuple[int, ...] = ()

    @property
    def size(self) -> int:
        return self.stop - self.start


@dataclass(frozen=True)
class StateVectorLayout:
    """Metadata required to reconstruct one UV state from a real vector."""

    representation: T3Representation
    reference_mu_gev: float
    shared_scalar: bool
    blocks: tuple[VectorBlock, ...]

    @property
    def size(self) -> int:
        if not self.blocks:
            return 0
        return self.blocks[-1].stop

    def block(self, name: str) -> VectorBlock:
        matches = [block for block in self.blocks if block.name == name]
        if len(matches) != 1:
            raise KeyError(
                f"Expected exactly one vector block named {name!r}, "
                f"found {len(matches)}."
            )
        return matches[0]

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(block.name for block in self.blocks)


class _VectorBuilder:
    def __init__(self) -> None:
        self.values: list[np.ndarray] = []
        self.blocks: list[VectorBlock] = []
        self.offset = 0

    def add_real_scalar(self, name: str, value: float) -> None:
        data = np.asarray([float(value)], dtype=float)
        self._append(name, "real_scalar", data, ())

    def add_complex_scalar(self, name: str, value: complex) -> None:
        z = complex(value)
        data = np.asarray([z.real, z.imag], dtype=float)
        self._append(name, "complex_scalar", data, ())

    def add_complex_matrix(self, name: str, value: np.ndarray) -> None:
        matrix = np.asarray(value, dtype=complex)
        data = np.concatenate(
            [
                matrix.real.reshape(-1),
                matrix.imag.reshape(-1),
            ]
        )
        self._append(name, "complex_matrix", data, matrix.shape)

    def _append(
        self,
        name: str,
        kind: BlockKind,
        data: np.ndarray,
        shape: tuple[int, ...],
    ) -> None:
        start = self.offset
        stop = start + int(data.size)
        self.values.append(data)
        self.blocks.append(
            VectorBlock(
                name=name,
                kind=kind,
                start=start,
                stop=stop,
                shape=shape,
            )
        )
        self.offset = stop

    def finish(
        self,
        *,
        representation: T3Representation,
        reference_mu_gev: float,
        shared_scalar: bool,
    ) -> tuple[np.ndarray, StateVectorLayout]:
        vector = (
            np.concatenate(self.values)
            if self.values
            else np.empty(0, dtype=float)
        )
        layout = StateVectorLayout(
            representation=representation,
            reference_mu_gev=float(reference_mu_gev),
            shared_scalar=shared_scalar,
            blocks=tuple(self.blocks),
        )
        return vector, layout


def _pack_sm(builder: _VectorBuilder, sm: SMNumericalState) -> None:
    """Pack the SM variables using RGBeta's coupling names."""

    builder.add_real_scalar("gY", sm.gY)
    builder.add_real_scalar("g2", sm.g2)
    builder.add_real_scalar("g3", sm.g3)

    builder.add_complex_matrix("yu", sm.yu)
    builder.add_complex_matrix("yd", sm.yd)
    builder.add_complex_matrix("ye", sm.ye)

    builder.add_real_scalar("lambdaH", sm.lambdaH)


def pack_uv_state(
    state: T3UVState | SharedT3UVState,
) -> tuple[np.ndarray, StateVectorLayout]:
    """Validate and pack one UV state into a real ODE vector."""

    validated = state.validated()
    builder = _VectorBuilder()
    _pack_sm(builder, validated.sm)

    if isinstance(validated, T3UVState):
        builder.add_complex_matrix("y1", validated.y1)
        builder.add_complex_matrix("y2", validated.y2)
        builder.add_complex_matrix("MF", validated.MF)

        builder.add_real_scalar("mS1Sq", validated.mS1Sq)
        builder.add_real_scalar("mS2Sq", validated.mS2Sq)

        builder.add_real_scalar("lambdaS1", validated.lambdaS1)
        builder.add_real_scalar("lambdaS2", validated.lambdaS2)
        builder.add_real_scalar("lambdaH1", validated.lambdaH1)
        builder.add_real_scalar("lambdaH2", validated.lambdaH2)
        builder.add_real_scalar("lambda12", validated.lambda12)
        builder.add_complex_scalar("lambdaT3", validated.lambdaT3)

        for name in (
            "lambdaH1Adj",
            "lambdaH2Adj",
            "lambdaS1Adj",
            "lambdaS2Adj",
            "lambda12Adj",
            "lambda12Cross",
        ):
            value = getattr(validated, name)
            if value is not None:
                builder.add_real_scalar(name, value)

        for name in (
            "lambdaHHdagS2S2",
            "lambdaHHdagS1barS1bar",
            "lambdaS1bar2S2bar2",
            "lambdaS1barS2S2bar2",
            "lambdaS1S1bar2S2bar",
            "lambdaHHdagS1barS2barCross",
        ):
            value = getattr(validated, name)
            if value is not None:
                builder.add_complex_scalar(name, value)

        return builder.finish(
            representation=validated.representation,
            reference_mu_gev=validated.mu_gev,
            shared_scalar=False,
        )

    builder.add_complex_matrix("h", validated.h)
    builder.add_complex_matrix("MF", validated.MF)
    builder.add_real_scalar("mSSq", validated.mSSq)
    builder.add_real_scalar("lambdaS", validated.lambdaS)
    builder.add_real_scalar("lambda3", validated.lambda3)
    builder.add_real_scalar("lambda4", validated.lambda4)
    builder.add_real_scalar("lambda5", validated.lambda5)

    return builder.finish(
        representation=validated.representation,
        reference_mu_gev=validated.mu_gev,
        shared_scalar=True,
    )


def _read_real_scalar(
    vector: np.ndarray,
    block: VectorBlock,
) -> float:
    data = vector[block.start:block.stop]
    if block.kind != "real_scalar" or data.size != 1:
        raise ValueError(f"Invalid real-scalar block for {block.name}.")
    return float(data[0])


def _read_complex_scalar(
    vector: np.ndarray,
    block: VectorBlock,
) -> complex:
    data = vector[block.start:block.stop]
    if block.kind != "complex_scalar" or data.size != 2:
        raise ValueError(f"Invalid complex-scalar block for {block.name}.")
    return complex(float(data[0]), float(data[1]))


def _read_complex_matrix(
    vector: np.ndarray,
    block: VectorBlock,
) -> np.ndarray:
    if block.kind != "complex_matrix" or len(block.shape) != 2:
        raise ValueError(f"Invalid complex-matrix block for {block.name}.")

    n = int(np.prod(block.shape))
    data = vector[block.start:block.stop]
    if data.size != 2 * n:
        raise ValueError(
            f"Block {block.name} has {data.size} entries; expected {2 * n}."
        )

    real = data[:n].reshape(block.shape)
    imag = data[n:].reshape(block.shape)
    return real + 1j * imag


def _check_vector(vector: np.ndarray, layout: StateVectorLayout) -> np.ndarray:
    result = np.asarray(vector, dtype=float)

    if result.ndim != 1:
        raise ValueError("State vector must be one-dimensional.")

    if result.size != layout.size:
        raise ValueError(
            f"State vector has length {result.size}; "
            f"layout expects {layout.size}."
        )

    if not np.all(np.isfinite(result)):
        raise ValueError("State vector contains non-finite entries.")

    return result


def unpack_uv_state(
    vector: np.ndarray,
    layout: StateVectorLayout,
    *,
    mu_gev: float | None = None,
) -> T3UVState | SharedT3UVState:
    """Reconstruct one validated UV state from a real ODE vector."""

    vector = _check_vector(vector, layout)
    blocks = {block.name: block for block in layout.blocks}

    def real(name: str) -> float:
        return _read_real_scalar(vector, blocks[name])

    def cscalar(name: str) -> complex:
        return _read_complex_scalar(vector, blocks[name])

    def cmatrix(name: str) -> np.ndarray:
        return _read_complex_matrix(vector, blocks[name])

    sm = SMNumericalState(
        gY=real("gY"),
        g2=real("g2"),
        g3=real("g3"),
        lambdaH=real("lambdaH"),
        yu=cmatrix("yu"),
        yd=cmatrix("yd"),
        ye=cmatrix("ye"),
    )

    scale = (
        layout.reference_mu_gev
        if mu_gev is None
        else float(mu_gev)
    )

    if layout.shared_scalar:
        return SharedT3UVState(
            mu_gev=scale,
            representation=layout.representation,
            sm=sm,
            h=cmatrix("h"),
            MF=cmatrix("MF"),
            mSSq=real("mSSq"),
            lambdaS=real("lambdaS"),
            lambda3=real("lambda3"),
            lambda4=real("lambda4"),
            lambda5=real("lambda5"),
        ).validated()

    optional_real_names = (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
    )
    optional_complex_names = (
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    )

    optional_real = {
        name: real(name) if name in blocks else None
        for name in optional_real_names
    }
    optional_complex = {
        name: cscalar(name) if name in blocks else None
        for name in optional_complex_names
    }

    return T3UVState(
        mu_gev=scale,
        representation=layout.representation,
        sm=sm,
        y1=cmatrix("y1"),
        y2=cmatrix("y2"),
        MF=cmatrix("MF"),
        mS1Sq=real("mS1Sq"),
        mS2Sq=real("mS2Sq"),
        lambdaS1=real("lambdaS1"),
        lambdaS2=real("lambdaS2"),
        lambdaH1=real("lambdaH1"),
        lambdaH2=real("lambdaH2"),
        lambda12=real("lambda12"),
        lambdaT3=cscalar("lambdaT3"),
        **optional_real,
        **optional_complex,
    ).validated()
