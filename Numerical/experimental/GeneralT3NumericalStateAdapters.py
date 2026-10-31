"""Workstream D Stage 4: native T3 state <-> Stage 3 EFT envelope adapters.

Only coordinates already held by the original native states are transported.
No RG callback, matching formula, or physical certification is supplied.

Complex block storage is delegated to the EXISTING production pack/unpack
functions. The optional quadratic mu2 coordinate is explicit because the
continuous runner appends it outside the ordinary RGBeta state layout.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from Numerical.core.State import T3UVState
from Numerical.core.StateVector import pack_uv_state, unpack_uv_state, StateVectorLayout
from Numerical.running.IntermediateScalarState import T3IntermediateScalarState
from Numerical.running.IntermediateScalarStateVector import (
    pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
)
from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, pack_general_scalar_state,
    unpack_general_scalar_state, _model_name,
)
from Numerical.experimental.GeneralT3EFTHandoffController import (
    BlockSpec, StageContract, EFTState, HandoffContractError,
)

QUADRATIC_KINDS = ('none', 'self_s1', 'self_s2', 'mixed')
MASS_SQUARED_NAMES = frozenset(('mS1Sq', 'mS2Sq', 'mSSq'))


def _spec_for_block(block: Any) -> BlockSpec:
    if block.kind not in ('real_scalar', 'complex_scalar', 'complex_matrix'):
        raise HandoffContractError(f'Unsupported native coordinate {block.name}: {block.kind}')
    shape = tuple(block.shape) if block.kind == 'complex_matrix' else ()
    dtype = 'real' if block.kind == 'real_scalar' else 'complex'
    units = ('GeV' if block.name == 'MF' else 'GeV^2' if block.name in MASS_SQUARED_NAMES
             else 'dimensionless')
    return BlockSpec(block.name, shape, dtype, units=units)


def _unpack_blocks(vector: Any, layout: StateVectorLayout) -> dict[str, Any]:
    """Read exactly the existing real/imag production packing convention."""
    x = np.asarray(vector, dtype=float)
    if x.shape != (layout.size,) or not np.isfinite(x).all():
        raise HandoffContractError('Invalid native packed vector')
    values: dict[str, Any] = {}
    previous = 0
    for block in layout.blocks:
        if block.start != previous or block.stop > layout.size or block.name in values:
            raise HandoffContractError('Non-contiguous or duplicate native layout blocks')
        raw = x[block.start:block.stop]
        if block.kind == 'real_scalar' and raw.size == 1:
            value: Any = float(raw[0])
        elif block.kind == 'complex_scalar' and raw.size == 2:
            value = complex(raw[0], raw[1])
        elif block.kind == 'complex_matrix':
            shape = tuple(block.shape)
            n = int(np.prod(shape))
            if len(shape) != 2 or raw.size != 2*n:
                raise HandoffContractError(f'{block.name}: malformed matrix block')
            value = raw[:n].reshape(shape) + 1j*raw[n:].reshape(shape)
        else:
            raise HandoffContractError(f'{block.name}: malformed coordinate block')
        values[block.name] = value
        previous = block.stop
    if previous != layout.size:
        raise HandoffContractError('Native layout does not cover entire packed vector')
    return values


def _repack_blocks(parameters: Mapping[str, Any], layout: StateVectorLayout) -> np.ndarray:
    pieces: list[np.ndarray] = []
    for block in layout.blocks:
        arr = np.asarray(parameters[block.name])
        if block.kind == 'real_scalar':
            pieces.append(np.asarray([float(arr)], dtype=float))
        elif block.kind == 'complex_scalar':
            z = complex(arr)
            pieces.append(np.asarray([z.real, z.imag], dtype=float))
        else:
            a = np.asarray(arr, dtype=complex)
            pieces.append(np.concatenate((a.real.ravel(), a.imag.ravel())))
    vector = np.concatenate(pieces)
    if vector.shape != (layout.size,) or not np.isfinite(vector).all():
        raise HandoffContractError('Native coordinates cannot be repacked')
    return vector


@dataclass(frozen=True)
class RestoredNativeState:
    native: T3UVState | T3IntermediateScalarState | GeneralT3ScalarRGState | SMWeinbergInitialConditions
    mu_gev: float
    quadratic_mu2_gev2: complex | None


@dataclass(frozen=True)
class NativeBridge:
    """Explicitly describes one numerical envelope, not a verified EFT map."""
    contract: StageContract
    state: EFTState
    source_kind: str
    layout: StateVectorLayout | None
    registry: Mapping[str, Any] | None
    quadratic_kind: str

    def restore(self, evolved: EFTState) -> RestoredNativeState:
        """Rebuild a SAME-stage native state; no threshold projection occurs."""
        validated = self.contract.validate(evolved)
        if validated.model != self.state.model:
            raise HandoffContractError('Cannot restore a different T3 model')
        mu = float(validated.mu_gev)
        mu2 = (complex(validated.parameters['mu2'])
               if self.quadratic_kind != 'none' else None)
        if self.source_kind == 'sm_weinberg':
            p = validated.parameters
            native = SMWeinbergInitialConditions(
                gY=float(p['gY']), g2=float(p['g2']), g3=float(p['g3']),
                lambdaH=float(p['lambdaH']),
                ye=np.array(p['ye'], copy=True), yu=np.array(p['yu'], copy=True),
                yd=np.array(p['yd'], copy=True),
                K=np.array(validated.wilson['C5'], copy=True),
            ).validated()
            return RestoredNativeState(native, mu, None)
        if self.layout is None:
            raise HandoffContractError('Native layout absent')
        vector = _repack_blocks(validated.parameters, self.layout)
        if self.source_kind == 'uv':
            native = unpack_uv_state(vector, self.layout, mu_gev=mu)
        elif self.source_kind == 'eft1':
            native = unpack_intermediate_scalar_state(vector, self.layout, mu_gev=mu)
        elif self.source_kind in ('general_uv', 'general_eft1'):
            if self.registry is None:
                raise HandoffContractError('General RGBeta registry missing')
            native = unpack_general_scalar_state(
                vector, self.layout, self.registry,
                uv=self.source_kind == 'general_uv', mu_gev=mu,
            )
        else:
            raise HandoffContractError('Unknown source state kind')
        if _model_name(native.baseline.representation if isinstance(native, GeneralT3ScalarRGState)
                       else native.representation) != validated.model:
            raise HandoffContractError('Restored model representation changed')
        return RestoredNativeState(native, mu, mu2)


def bridge_native_state(
    model: str, native: Any, *, stage: str,
    active_particles: tuple[str, ...],
    registry: Mapping[str, Any] | None = None,
    quadratic_kind: str = 'none', quadratic_mu2_gev2: complex | None = None,
    mu_gev: float | None = None,
    c5_operator_basis: str = '', c5_normalization: str = '',
) -> NativeBridge:
    """Construct a strictly validated Stage 3 envelope from a native state.

    `SMWeinbergInitialConditions` has no intrinsic scale, so `mu_gev` is
    mandatory for that input and its K coefficient is mapped to envelope C5.
    No matching or physical provenance is inferred from the passed strings.
    """
    if quadratic_kind not in QUADRATIC_KINDS:
        raise HandoffContractError('Unsupported scalar quadratic structure')
    if quadratic_kind == 'none' and quadratic_mu2_gev2 is not None:
        raise HandoffContractError('No mu2 coordinate in this model')
    if quadratic_kind != 'none' and quadratic_mu2_gev2 is None:
        raise HandoffContractError('The bilinear mu2 coordinate is required')
    if not isinstance(model, str) or not model.strip() or not isinstance(stage, str) or not stage.strip():
        raise HandoffContractError('Model and stage labels are required')
    layout: StateVectorLayout | None = None
    if isinstance(native, SMWeinbergInitialConditions):
        if quadratic_kind != 'none' or registry is not None:
            raise HandoffContractError('SM+Weinberg stage cannot contain general scalar coordinates')
        if mu_gev is None:
            raise HandoffContractError('SM+Weinberg scale must be explicit')
        src = 'sm_weinberg'
        validated = native.validated()
        p = {name: getattr(validated, name) for name in
             ('gY', 'g2', 'g3', 'lambdaH', 'yu', 'yd', 'ye')}
        specs = tuple(BlockSpec(name, (3,3) if name in ('yu','yd','ye') else (),
                                'complex' if name in ('yu','yd','ye') else 'real')
                      for name in p)
        if not c5_operator_basis.strip() or not c5_normalization.strip():
            raise HandoffContractError('C5 requires explicitly supplied source basis/normalization')
        wspec = (BlockSpec('C5',(3,3),'complex',True,'GeV^-1',
                           operator_basis=c5_operator_basis,
                           normalization=c5_normalization,mass_dimension=-1),)
        wilson = {'C5': validated.K}
        scale = mu_gev
    else:
        if isinstance(native, GeneralT3ScalarRGState):
            if registry is None:
                raise HandoffContractError('General RGBeta state requires the original registry')
            validated = native.validated(registry)
            src = 'general_uv' if isinstance(validated.baseline, T3UVState) else 'general_eft1'
            scale = validated.baseline.mu_gev
            rep = validated.baseline.representation
            vector, layout = pack_general_scalar_state(validated, registry)
        elif isinstance(native, T3UVState):
            validated = native.validated()
            src = 'uv'
            scale = validated.mu_gev
            rep = validated.representation
            vector, layout = pack_uv_state(validated)
        elif isinstance(native, T3IntermediateScalarState):
            validated = native.validated()
            src = 'eft1'
            scale = validated.mu_gev
            rep = validated.representation
            vector, layout = pack_intermediate_scalar_state(validated)
        else:
            raise HandoffContractError('Unsupported native state type')
        if rep.shared_scalar or _model_name(rep) != model:
            raise HandoffContractError('Model label does not match native T3 representation')
        if mu_gev is not None and not np.isclose(float(mu_gev), float(scale), rtol=1e-13, atol=0):
            raise HandoffContractError('Native state already declares a different scale')
        p = _unpack_blocks(vector, layout)
        specs = tuple(_spec_for_block(block) for block in layout.blocks)
        wspec = ()
        wilson = {}
    if quadratic_kind != 'none':
        if 'mu2' in p:
            raise HandoffContractError('Native layout already owns mu2')
        p['mu2'] = complex(quadratic_mu2_gev2)
        specs += (BlockSpec('mu2',(), 'complex',units='GeV^2'),)
    contract = StageContract(stage, active_particles, specs, wspec)
    envelope = contract.validate(EFTState(model, stage, scale, active_particles, p, wilson))
    return NativeBridge(contract, envelope, src, layout, registry, quadratic_kind)
