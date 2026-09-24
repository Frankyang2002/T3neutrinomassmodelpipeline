"""Pack and unpack the renormalisable EFT1 numerical state.

EFT1 is the theory after the heavy fermion F has been integrated out:

ordinary T3:
    SM + S1 + S2

shared-scalar T3:
    SM + S

The real ODE-vector convention matches ``Numerical.StateVector``:
complex matrices are stored as flattened real entries followed by flattened
imaginary entries, and complex scalars as (Re, Im).
"""

from __future__ import annotations

import numpy as np

from Numerical.EFT1State import (
    SharedT3EFT1State,
    T3EFT1State,
)
from Numerical.StateVector import StateVectorLayout, VectorBlock


EFT1State = T3EFT1State | SharedT3EFT1State


class _Builder:
    def __init__(self) -> None:
        self.values: list[np.ndarray] = []
        self.blocks: list[VectorBlock] = []
        self.offset = 0

    def _append(
        self,
        name: str,
        kind: str,
        data: np.ndarray,
        shape: tuple[int, ...] = (),
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

    def real(self, name: str, value: float) -> None:
        self._append(
            name,
            "real_scalar",
            np.asarray([float(value)], dtype=float),
        )

    def complex_scalar(self, name: str, value: complex) -> None:
        z = complex(value)
        self._append(
            name,
            "complex_scalar",
            np.asarray([z.real, z.imag], dtype=float),
        )

    def complex_matrix(self, name: str, value: np.ndarray) -> None:
        matrix = np.asarray(value, dtype=complex)
        self._append(
            name,
            "complex_matrix",
            np.concatenate(
                [
                    matrix.real.reshape(-1),
                    matrix.imag.reshape(-1),
                ]
            ),
            matrix.shape,
        )

    def finish(
        self,
        *,
        state: EFT1State,
    ) -> tuple[np.ndarray, StateVectorLayout]:
        vector = np.concatenate(self.values)
        layout = StateVectorLayout(
            representation=state.representation,
            reference_mu_gev=float(state.mu_gev),
            shared_scalar=state.representation.shared_scalar,
            blocks=tuple(self.blocks),
        )
        return vector, layout


def _pack_sm(builder: _Builder, state: EFT1State) -> None:
    builder.real("gY", state.sm.gY)
    builder.real("g2", state.sm.g2)
    builder.real("g3", state.sm.g3)
    builder.complex_matrix("yu", state.sm.yu)
    builder.complex_matrix("yd", state.sm.yd)
    builder.complex_matrix("ye", state.sm.ye)
    builder.real("lambdaH", state.sm.lambdaH)


def pack_eft1_state(
    state: EFT1State,
) -> tuple[np.ndarray, StateVectorLayout]:
    """Validate and pack one EFT1 state into a real ODE vector."""

    state = state.validated()
    builder = _Builder()
    _pack_sm(builder, state)

    if isinstance(state, SharedT3EFT1State):
        builder.real("mSSq", state.mSSq)
        builder.real("lambdaS", state.lambdaS)
        builder.real("lambda3", state.lambda3)
        builder.real("lambda4", state.lambda4)
        builder.real("lambda5", state.lambda5)
        return builder.finish(state=state)

    builder.real("mS1Sq", state.mS1Sq)
    builder.real("mS2Sq", state.mS2Sq)
    builder.real("lambdaS1", state.lambdaS1)
    builder.real("lambdaS2", state.lambdaS2)
    builder.real("lambdaH1", state.lambdaH1)
    builder.real("lambdaH2", state.lambdaH2)
    builder.real("lambda12", state.lambda12)
    builder.complex_scalar("lambdaT3", state.lambdaT3)

    for name in (
        "lambdaH1Adj",
        "lambdaH2Adj",
        "lambdaS1Adj",
        "lambdaS2Adj",
        "lambda12Adj",
        "lambda12Cross",
    ):
        value = getattr(state, name)
        if value is not None:
            builder.real(name, value)

    for name in (
        "lambdaHHdagS2S2",
        "lambdaHHdagS1barS1bar",
        "lambdaS1bar2S2bar2",
        "lambdaS1barS2S2bar2",
        "lambdaS1S1bar2S2bar",
        "lambdaHHdagS1barS2barCross",
    ):
        value = getattr(state, name)
        if value is not None:
            builder.complex_scalar(name, value)

    return builder.finish(state=state)


def _check_vector(
    vector: np.ndarray,
    layout: StateVectorLayout,
) -> np.ndarray:
    result = np.asarray(vector, dtype=float)

    if result.ndim != 1:
        raise ValueError("EFT1 state vector must be one-dimensional.")

    if result.size != layout.size:
        raise ValueError(
            f"EFT1 state vector has length {result.size}; "
            f"layout expects {layout.size}."
        )

    if not np.all(np.isfinite(result)):
        raise ValueError("EFT1 state vector contains non-finite entries.")

    return result


def unpack_eft1_state(
    vector: np.ndarray,
    layout: StateVectorLayout,
    *,
    mu_gev: float | None = None,
) -> EFT1State:
    """Reconstruct one validated EFT1 state."""

    from Numerical.State import SMNumericalState

    vector = _check_vector(vector, layout)
    blocks = {block.name: block for block in layout.blocks}

    def real(name: str) -> float:
        block = blocks[name]
        data = vector[block.start:block.stop]
        if block.kind != "real_scalar" or data.size != 1:
            raise ValueError(f"Invalid real block for {name}.")
        return float(data[0])

    def cscalar(name: str) -> complex:
        block = blocks[name]
        data = vector[block.start:block.stop]
        if block.kind != "complex_scalar" or data.size != 2:
            raise ValueError(f"Invalid complex-scalar block for {name}.")
        return complex(float(data[0]), float(data[1]))

    def cmatrix(name: str) -> np.ndarray:
        block = blocks[name]
        if block.kind != "complex_matrix":
            raise ValueError(f"Invalid complex-matrix block for {name}.")
        n = int(np.prod(block.shape))
        data = vector[block.start:block.stop]
        return (
            data[:n].reshape(block.shape)
            + 1j * data[n:].reshape(block.shape)
        )

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
        return SharedT3EFT1State(
            mu_gev=scale,
            representation=layout.representation,
            sm=sm,
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

    return T3EFT1State(
        mu_gev=scale,
        representation=layout.representation,
        sm=sm,
        mS1Sq=real("mS1Sq"),
        mS2Sq=real("mS2Sq"),
        lambdaS1=real("lambdaS1"),
        lambdaS2=real("lambdaS2"),
        lambdaH1=real("lambdaH1"),
        lambdaH2=real("lambdaH2"),
        lambda12=real("lambda12"),
        lambdaT3=cscalar("lambdaT3"),
        **{
            name: real(name) if name in blocks else None
            for name in optional_real_names
        },
        **{
            name: cscalar(name) if name in blocks else None
            for name in optional_complex_names
        },
    ).validated()
