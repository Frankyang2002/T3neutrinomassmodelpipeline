"""Workstream D Stage 3: explicit, fail-closed numerical EFT handoff controller.

This module does not derive beta functions, project fields, or calculate matching.
It invokes externally registered evolution/matching functions, each with an
explicit source label, and validates numerical states on BOTH sides of each
boundary. No callback is supplied by default. In particular, the candidate
full-UV mass crossings in Stage 2 are not already matched EFT boundaries.

State coordinates and Wilson tensors are keyed arrays with declared shapes and
real/complex symmetry. Scale is metadata (GeV), not an ODE coordinate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Callable, Mapping

import numpy as np


class HandoffContractError(ValueError):
    """Invalid numerical state, event sequence, or callback result."""


class IncompleteCandidateInventory(HandoffContractError):
    """A candidate threshold list omits one or more expected mass branches."""


def _positive_scale(value: float, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise HandoffContractError(f'{label} must be a finite positive scale') from exc
    if not (math.isfinite(number) and number > 0.0):
        raise HandoffContractError(f'{label} must be a finite positive scale')
    return number


def _unique_names(names: tuple[str, ...], description: str) -> None:
    if any(not isinstance(n, str) or not n.strip() for n in names) or len(set(names)) != len(names):
        raise HandoffContractError(f'{description} must contain unique nonempty names')


@dataclass(frozen=True)
class BlockSpec:
    name: str
    shape: tuple[int, ...] = ()
    dtype: str = 'real'                 # real or complex
    symmetric: bool = False             # transpose, NOT Hermitian conjugation
    units: str = ''
    operator_basis: str = ''      # Required when this spec is a Wilson block
    normalization: str = ''       # Required when this spec is a Wilson block
    mass_dimension: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise HandoffContractError('Block name is missing')
        if self.dtype not in ('real', 'complex'):
            raise HandoffContractError(f'{self.name}: block dtype must be real or complex')
        if not isinstance(self.shape, tuple) or any(type(n) is not int or n < 1 for n in self.shape):
            raise HandoffContractError(f'{self.name}: invalid shape')
        if self.symmetric and (len(self.shape) != 2 or self.shape[0] != self.shape[1]):
            raise HandoffContractError(f'{self.name}: symmetric block must be square')

    def validate(self, raw: Any) -> np.ndarray:
        try:
            z = np.asarray(raw, dtype=complex)
        except (TypeError, ValueError, OverflowError) as exc:
            raise HandoffContractError(f'{self.name}: non-numerical block') from exc
        if z.shape != self.shape or not np.isfinite(z).all():
            raise HandoffContractError(f'{self.name}: expected finite shape {self.shape}, got {z.shape}')
        if self.dtype == 'real':
            if np.any(z.imag != 0.0):
                raise HandoffContractError(f'{self.name}: complex value in real coordinate')
            array = z.real.copy()
        else:
            array = z.copy()
        if self.symmetric and not np.allclose(array, array.T, rtol=0, atol=1e-11):
            raise HandoffContractError(f'{self.name}: transpose symmetry violated')
        return array


@dataclass(frozen=True)
class StageContract:
    name: str
    active_particles: tuple[str, ...]
    parameters: tuple[BlockSpec, ...]
    wilson: tuple[BlockSpec, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise HandoffContractError('Stage name is missing')
        _unique_names(self.active_particles, f'{self.name}: active particles')
        _unique_names(tuple(x.name for x in self.parameters), f'{self.name}: parameter blocks')
        _unique_names(tuple(x.name for x in self.wilson), f'{self.name}: Wilson blocks')
        if set(x.name for x in self.parameters) & set(x.name for x in self.wilson):
            raise HandoffContractError(f'{self.name}: parameter/Wilson namespaces collide')
        for x in self.wilson:
            if not x.operator_basis.strip() or not x.normalization.strip() or type(x.mass_dimension) is not int:
                raise HandoffContractError(f'{x.name}: Wilson basis, normalization and mass dimension required')
            if x.name == 'C5' and (x.shape != (3, 3) or x.dtype != 'complex' or not x.symmetric):
                raise HandoffContractError('C5 must be a complex symmetric 3x3 flavour matrix')

    def validate(self, state: 'EFTState') -> 'EFTState':
        if state.stage != self.name or set(state.active_particles) != set(self.active_particles):
            raise HandoffContractError(f'{self.name}: stage or particle inventory disagrees')
        if state.active_particles != self.active_particles:
            raise HandoffContractError(f'{self.name}: particle order differs from declared contract')
        if not isinstance(state.parameters, Mapping) or not isinstance(state.wilson, Mapping):
            raise HandoffContractError(f'{self.name}: state blocks must be mappings')
        ps = {x.name: x for x in self.parameters}
        ws = {x.name: x for x in self.wilson}
        if set(state.parameters) != set(ps) or set(state.wilson) != set(ws):
            raise HandoffContractError(f'{self.name}: state block keys differ from stage contract')
        return EFTState(
            model=state.model, stage=state.stage, mu_gev=_positive_scale(state.mu_gev, 'state scale'),
            active_particles=self.active_particles,
            parameters={name: spec.validate(state.parameters[name]) for name, spec in ps.items()},
            wilson={name: spec.validate(state.wilson[name]) for name, spec in ws.items()},
        )


@dataclass(frozen=True)
class EFTState:
    model: str
    stage: str
    mu_gev: float
    active_particles: tuple[str, ...]
    parameters: Mapping[str, Any]
    wilson: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ThresholdEvent:
    particle: str
    kind: str
    scale_gev: float
    next_stage: str

    def __post_init__(self) -> None:
        if self.kind not in ('fermion', 'scalar') or not self.particle.startswith(self.kind + ':'):
            raise HandoffContractError('Invalid event kind/particle')
        if not isinstance(self.next_stage, str) or not self.next_stage.strip():
            raise HandoffContractError('Event requires explicit next EFT stage')
        _positive_scale(self.scale_gev, 'event scale')


@dataclass(frozen=True)
class Callback:
    name: str
    function: Callable[..., EFTState]
    source: str
    # Explicitly descriptive only: no automatic physical certification.
    status: str = 'candidate'

    def __post_init__(self) -> None:
        if not self.name or not self.source or not callable(self.function):
            raise HandoffContractError('Callback requires name, callable, and source')
        if self.status not in ('candidate', 'externally_verified'):
            raise HandoffContractError('Unrecognized callback provenance status')


@dataclass(frozen=True)
class HandoffLog:
    operation: str
    stage_before: str
    stage_after: str
    scale_start_gev: float
    scale_end_gev: float
    callback: str
    source: str
    declared_callback_status: str
    particle: str = ''


@dataclass(frozen=True)
class HandoffResult:
    status: str                    # Complete or Blocked
    reason: str                    # empty on complete
    state: EFTState                # last validated, never a partially matched state
    next_event: ThresholdEvent | None
    history: tuple[HandoffLog, ...]
    physical_matching_certified: bool = False


class EFTHandoffController:
    """Stage-specific evolution and particle-specific matching callback dispatcher."""

    def __init__(self) -> None:
        self.stages: dict[str, StageContract] = {}
        self.runners: dict[str, Callback] = {}
        self.matchers: dict[tuple[str, str, str], Callback] = {}

    def add_stage(self, contract: StageContract) -> None:
        if contract.name in self.stages:
            raise HandoffContractError(f'Duplicate stage {contract.name}')
        self.stages[contract.name] = contract

    def add_runner(self, stage: str, callback: Callback) -> None:
        if stage not in self.stages or stage in self.runners:
            raise HandoffContractError(f'Unknown or duplicated runner stage {stage}')
        self.runners[stage] = callback

    def add_matcher(self, from_stage: str, particle: str, next_stage: str, callback: Callback) -> None:
        if from_stage not in self.stages or next_stage not in self.stages or from_stage == next_stage:
            raise HandoffContractError('Matching registration needs two different known stages')
        before, after = self.stages[from_stage], self.stages[next_stage]
        if particle not in before.active_particles or set(before.active_particles) - {particle} != set(after.active_particles):
            raise HandoffContractError('Matching transition must remove exactly its declared particle')
        key = (from_stage, particle, next_stage)
        if key in self.matchers:
            raise HandoffContractError('Duplicate matching callback')
        self.matchers[key] = callback

    def run(self, initial: EFTState, events: tuple[ThresholdEvent, ...], *,
            final_scale_gev: float) -> HandoffResult:
        goal = _positive_scale(final_scale_gev, 'final scale')
        if initial.stage not in self.stages or not initial.model:
            raise HandoffContractError('Unregistered initial stage or missing model')
        current = self.stages[initial.stage].validate(initial)
        if not goal < current.mu_gev:
            raise HandoffContractError('Final scale must be below initial scale')
        if not isinstance(events, tuple):
            raise HandoffContractError('Events must be an explicit tuple')
        seen: set[str] = set()
        prev = current.mu_gev
        for event in events:
            scale = _positive_scale(event.scale_gev, 'event scale')
            if event.particle in seen or not goal < scale < prev:
                raise HandoffContractError('Events must be unique, strictly descending, within interval')
            # Independent roots closer than this cannot be safely ordered.
            if (prev - scale) <= 1e-10 * prev:
                raise HandoffContractError('Unresolved simultaneous or near-degenerate events')
            seen.add(event.particle)
            prev = scale
        history: list[HandoffLog] = []

        def blocked(reason: str, upcoming: ThresholdEvent | None) -> HandoffResult:
            return HandoffResult('Blocked', reason, current, upcoming, tuple(history))

        for event in (*events, None):
            next_scale = float(event.scale_gev) if event is not None else goal
            runner = self.runners.get(current.stage)
            if runner is None:
                return blocked(f'Missing RG running callback for stage {current.stage}', event)
            evolved = runner.function(current, next_scale)
            if not isinstance(evolved, EFTState):
                raise HandoffContractError(f'{runner.name}: running callback returned non-EFT state')
            evolved = self.stages[current.stage].validate(evolved)
            if evolved.model != current.model or not math.isclose(evolved.mu_gev, next_scale, rel_tol=1e-11):
                raise HandoffContractError(f'{runner.name}: running result model or endpoint mismatch')
            history.append(HandoffLog('run', current.stage, current.stage, current.mu_gev,
                                      evolved.mu_gev, runner.name, runner.source, runner.status))
            current = evolved
            if event is None:
                return HandoffResult('Complete', '', current, None, tuple(history))
            key = (current.stage, event.particle, event.next_stage)
            before = self.stages[current.stage]
            after = self.stages.get(event.next_stage)
            if event.particle not in before.active_particles or after is None or \
                    set(before.active_particles) - {event.particle} != set(after.active_particles):
                raise HandoffContractError(f'{current.stage}: event is not a legal one-particle transition')
            matcher = self.matchers.get(key)
            if matcher is None:
                return blocked(f'Missing matching callback for {key}', event)
            matched = matcher.function(current, event)
            if not isinstance(matched, EFTState):
                raise HandoffContractError(f'{matcher.name}: matching callback returned non-EFT state')
            matched = after.validate(matched)
            if matched.model != current.model or not math.isclose(matched.mu_gev, current.mu_gev, rel_tol=1e-11):
                raise HandoffContractError(f'{matcher.name}: matching result model or scale mismatch')
            history.append(HandoffLog('match', current.stage, matched.stage,
                                      current.mu_gev, matched.mu_gev, matcher.name, matcher.source,
                                      matcher.status, event.particle))
            current = matched
        raise AssertionError('unreachable')


@dataclass(frozen=True)
class CandidateEvent:
    particle: str
    kind: str
    scale_gev: float
    source_status: str


@dataclass(frozen=True)
class CandidateInventory:
    model: str
    events: tuple[CandidateEvent, ...]
    unresolved: tuple[str, ...]
    complete: bool

    def bind_stages(self, stages_after_particle: Mapping[str, str]) -> tuple[ThresholdEvent, ...]:
        """No implicit EFT field projection: caller must declare EVERY next stage."""
        if not self.complete:
            raise IncompleteCandidateInventory(f'{self.model}: missing {self.unresolved}')
        if set(stages_after_particle) != {e.particle for e in self.events}:
            raise HandoffContractError('Stage binding must cover each event exactly once')
        return tuple(ThresholdEvent(e.particle, e.kind, e.scale_gev,
                                    stages_after_particle[e.particle]) for e in self.events)


def candidate_inventory_from_stage2(row: Mapping[str, Any]) -> CandidateInventory:
    """Validate Stage 2 candidate mass roots; do not infer an EFT or matching.

    Rebuild ordering from source crossing records; reject inconsistencies in the
    cached candidate_event_order list. Model-specific scalar branches are those
    from Stage 1's current four *candidate* quadratic structures, not a claim
    about the certified broken-phase spectrum.
    """
    model = row.get('model')
    kind = row.get('quadratic_kind')
    if not isinstance(model, str) or kind not in ('none', 'mixed', 'self_s1', 'self_s2'):
        raise HandoffContractError('Invalid candidate model/quadratic kind')
    if row.get('status') != 'CandidateRunningFermionCrossingAuditNoMatching' or row.get('uv_scalar_crossings') != []:
        raise HandoffContractError(f'{model}: missing audit data or UV scalar crossings require review')
    expected_scalar = {
        'none': {'scalar:S1', 'scalar:S2'},
        'mixed': {'scalar:Y', 'scalar:X'},
        'self_s1': {'scalar:S1_plus', 'scalar:S1_minus', 'scalar:S2'},
        'self_s2': {'scalar:S1', 'scalar:S2_plus', 'scalar:S2_minus'},
    }[kind]
    expected = {'fermion:0', 'fermion:1', 'fermion:2'} | expected_scalar
    crossings = list(row.get('uv_fermion_crossings', [])) + list(row.get('eft1_scalar_crossings', []))
    unresolved_fermion = list(row.get('unresolved_fermion_crossings', []))
    found: dict[str, CandidateEvent] = {}
    for item in crossings:
        name = item.get('particle')
        mass_kind = item.get('kind')
        if not isinstance(name, str) or name not in expected or name in found or \
                name.split(':', 1)[0] != mass_kind:
            raise HandoffContractError(f'{model}: unexpected/duplicated mass branch {name}')
        found[name] = CandidateEvent(name, mass_kind, _positive_scale(item.get('scale_gev'), name),
                                     str(item.get('status', '')))
    pending: set[str] = set()
    for item in unresolved_fermion:
        name = item.get('particle')
        if name not in expected or name in pending or name in found or not name.startswith('fermion:'):
            raise HandoffContractError(f'{model}: invalid unresolved fermion branch {name}')
        pending.add(name)
    missing = expected - set(found)
    if not pending.issubset(missing):
        raise HandoffContractError(f'{model}: unresolved fermion labels inconsistent')
    if len(found) + len(missing) != len(expected):
        raise HandoffContractError(f'{model}: unexpected branch count')
    events = tuple(sorted(found.values(), key=lambda e: -e.scale_gev))
    if any(events[i].scale_gev <= events[i+1].scale_gev * (1 + 1e-10)
           for i in range(len(events)-1)):
        raise HandoffContractError(f'{model}: nearly coincident candidate mass roots')
    cached = row.get('candidate_event_order', [])
    if not isinstance(cached, list) or len(cached) != len(events):
        raise HandoffContractError(f'{model}: event ordering cache incomplete')
    for entry, expected_entry in zip(cached, events):
        if entry.get('particle') != expected_entry.particle or not math.isclose(
                float(entry.get('scale_gev', -1)), expected_entry.scale_gev, rel_tol=1e-12):
            raise HandoffContractError(f'{model}: event cache differs from recomputed ordering')
    # Stage 2 reports scalar branches only in EFT1 after a grouped F handoff.
    # Both are candidate diagnostics and cannot be concatenated into a valid
    # matched numerical trajectory without additional B/C callbacks.
    return CandidateInventory(model, events, tuple(sorted(missing)), not missing)
