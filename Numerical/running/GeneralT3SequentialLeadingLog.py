"""All-16 exploratory fermion-first LLSS -> Weinberg leading-one-loop bridge.

UV -> integrate out ALL F -> EFT1(S1,S2) -> integrate out ALL scalars -> SMEFT.

Uses the repository's source-resolved Matchete F-tree Wilson tensors and the
joint-scalar master kernel. Constructs SAME-SCALE HHSS tensors from all ordinary
portal directions, the retained native Matchete channels, and *every* extra
HHSS RGBeta field-count sector in the full-general 16-model registry.

Extra invariants use an EXPLICIT, fixed candidate basis identification:
    registry symbol j in a field-count sector -> SU(2) reference invariant j.
This is an algebraic convention for numerical exploration, not a proof that
RGBeta/Matchete native couplings have that basis or overall normalization.
The Matchete GammaCC -> Weyl factor is a user-supplied complex number.
No guessed relative couplings are discarded or silently set to zero.

The scalar-pair contraction yields the leading one-loop log and one simultaneous
finite scalar bubble. It is the correct leading order for tree-level LLSS
sources, but is NOT an all-seven-term physical EFT1 C5 calculation; other
operator sources, hard fermion matching and an independent normalization proof
remain outside this bridge. The seven-term evaluator remains available in
RGE/running/intermediate/GeneralEFT1C5Mixing.py.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

MODELS = tuple(f"{family}_{a:+d}" for family, alphas in {
    'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
    'D': (-2, 0, 2), 'E': (-4, -2, 0, 2),
}.items() for a in alphas)
NATIVE_ALIAS = {
    'lambdaT3': 'lambdaT3',
    'lambdaH1Inv1': 'lambdaH1', 'lambdaH1Inv2': 'lambdaH1Adj',
    'lambdaH2Inv1': 'lambdaH2', 'lambdaH2Inv2': 'lambdaH2Adj',
}
PORTAL_NATIVE = frozenset(NATIVE_ALIAS) - {'lambdaT3'}
DIMS = {'A': (1, 3), 'B': (2, 2), 'C': (2, 2),
        'D': (3, 1), 'E': (3, 3)}
LOOP = 1 / (16 * math.pi**2)


def require(ok: bool, message: str) -> None:
    if not bool(ok):
        raise ValueError('General T3 sequential leading log: ' + message)


def number(raw: Any, name: str) -> complex:
    if isinstance(raw, Mapping) and {'real','imag'} <= set(raw):
        z = complex(float(raw['real']), float(raw['imag']))
    elif isinstance(raw, (int, float, complex, np.number)):
        z = complex(raw)
    else:
        raise ValueError(f'{name}: invalid complex scalar representation')
    require(math.isfinite(z.real) and math.isfinite(z.imag), f'{name}: nonfinite')
    return z


def encode(a: np.ndarray) -> dict:
    v = np.asarray(a, complex)
    require(v.shape == (3, 3) and np.isfinite(v).all(), 'invalid 3x3 matrix')
    return {'real': v.real.tolist(), 'imag': v.imag.tolist()}


def decode(raw: Any) -> np.ndarray:
    if isinstance(raw, Mapping) and {'real','imag'} <= set(raw):
        v = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
    else:
        v = np.asarray(raw, complex)
    require(v.shape == (3, 3) and np.isfinite(v).all(), 'invalid complex matrix')
    return v


def native_aliases(model: str, names: Any, state: Mapping[str, Any],
                   registry: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, str]]:
    """Candidate complete native channel -> running coupling crosswalk.

    Extra names are resolved by exact candidate Matchete symbol from the
    registry, never just by alphabetical order across different monomials.
    """
    from Numerical.lo_certification.couplings.extra_quartic_source_v6 import (
        model_registry_quartics,
    )
    rows = model_registry_quartics(model, registry)
    candidate = {row['Matchete_generator_candidate_symbol']: row['rgbeta_symbol']
                 for row in rows}
    aliases = {}
    values = {}
    for native in sorted(set(names)):
        symbol = NATIVE_ALIAS.get(native, candidate.get(native))
        require(symbol is not None,
                f'{model}: no source-grounded candidate crosswalk for {native}')
        require(symbol in state,
                f'{model}: required scalar coupling {symbol} ({native}) absent')
        aliases[native] = symbol
        values[native] = state[symbol]
    return values, aliases


def geometry_contribution(*, n: int, real_rows: list[dict],
                          imag_rows: list[dict], z: complex,
                          label: str) -> np.ndarray:
    """Real HHSS fourth derivative of V4=-z P-conj(z) conj(P).

    Repository reference geometry provides unit Re(z) and Im(z) directions.
    Each row lists one independent H/H and S/S symmetric component.
    """
    import sympy as sp
    tensor = np.zeros((4, 4, n, n), float)
    assigned: dict[tuple[int,int,int,int], float] = {}
    for rows, factor in ((real_rows, z.real), (imag_rows, z.imag)):
        for item in rows:
            a,b,i,j = item['H1_H2_S1_S2_0based']
            require(0 <= a < 4 and 0 <= b < 4 and 0 <= i < n and 0 <= j < n,
                    f'{label}: bad derivative component')
            expr = sp.sympify(item['value_exact'], locals={'sqrt': sp.sqrt, 'I': sp.I})
            value = complex(sp.N(expr, 17))
            require(abs(value.imag) < 1e-12 and math.isfinite(value.real),
                    f'{label}: nonreal derivative')
            for A,B in {(a,b),(b,a)}:
                for I,J in {(i,j),(j,i)}:
                    key = (A,B,I,J)
                    assigned[key] = assigned.get(key, 0.) + factor * value.real
    for k, value in assigned.items():
        tensor[k] = value
    require(np.isfinite(tensor).all(), f'{label}: nonfinite HHSS tensor')
    return tensor


def general_extra_channels(model: str, registry: Mapping[str, Any],
                           state: Mapping[str, Any], *,
                           excluded: set[str] | frozenset[str] = frozenset(),
                           geometry: dict | None = None
                           ) -> tuple[list[str], list[np.ndarray], list[dict]]:
    """Evaluate ALL extra 2-Higgs, 2-heavy-scalar invariant directions.

    The same generated registry includes quartics without two Higgs legs; those
    have no direct LLSS->LLHH scalar-pair insertion at this order and must stay
    in renormalisable running, not be inserted in this contraction.
    """
    if geometry is None:
        from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
        geometry = model_extra_geometry(model, registry['models'][model],
                                        include_components=True)
    d1,d2 = DIMS[model[0]]
    n = 2*(d1+d2)
    names, tensors, details = [], [], []
    all_hhss = 0
    for sector in geometry['monomial_sectors']:
        symbols = sector['source_registry_symbols']
        basis = sector['exact_reference_real_component_tensor_basis']
        require(len(symbols) == len(basis),
                f'{model}: incomplete reference basis for HHSS sector')
        for sym, ref in zip(symbols, basis):
            label = sym['rgbeta_symbol']
            all_hhss += 1
            if label in excluded:
                continue
            require(label in state, f'{model}: missing GENERAL scalar coupling {label}')
            z = number(state[label], label)
            q = geometry_contribution(
                n=n, real_rows=ref['Re_c_unit_tensor_components'],
                imag_rows=ref['Im_c_unit_tensor_components'], z=z, label=label)
            names.append(label)
            tensors.append(q)
            details.append({'name': label, 'real': z.real, 'imag': z.imag,
                            'reference_invariant': ref['id'],
                            'field_counts': sector['counts_H_BarH_S1_BarS1_S2_BarS2'],
                            'assumption': 'sector_order_candidate_RGBeta_to_SU2_reference_basis'})
    require(all_hhss == geometry['registry_HHSS_extra_symbols'],
            f'{model}: general HHSS registry inventory mismatch')
    return names, tensors, details


@dataclass(frozen=True)
class TensorBundle:
    model: str
    mass2: np.ndarray
    quartic: np.ndarray
    channels: np.ndarray
    names: tuple[str, ...]
    provenance: dict


def scalar_tensors(model: str, boundary: Any, ch: Any,
                   snapshot: Any, registry: Mapping[str, Any], *,
                   include_general_extra: bool = True,
                   include_portal_completion: bool = True,
                   geometry: dict | None = None) -> TensorBundle:
    """Build scalar tensor on the source's fixed real basis, including mixing.

    snapshot is the existing numerical_benchmark_v5.RunningSnapshot; a caller
    may construct one from another same-scale *general* EFT1 endpoint.
    """
    from Numerical.joint_scalar_matching.source_quartic_v4 import evaluate_channels
    from Numerical.joint_scalar_matching.special_mass_v6 import (
        build_special_basis, SPECIAL,
    )
    from Numerical.lo_certification.couplings.portal_completion import (
        completed_portal_numeric,
    )
    require(model in MODELS and boundary.model == model and ch.model == model,
            'scalar model/boundary/channel mismatch')
    state = snapshot.couplings
    d1,d2 = DIMS[model[0]]
    n = 2*(d1+d2)
    require(len(boundary.labels) == n, f'{model}: source real scalar shape mismatch')
    values, aliases = native_aliases(model, ch.by_name, state, registry)
    native_tensors, native_names = evaluate_channels(ch, values)
    names: list[str] = []
    parts: list[np.ndarray] = []
    descriptions: list[dict] = []
    covered_extra = set()
    for k, raw in zip(native_names, native_tensors):
        if k in PORTAL_NATIVE and include_portal_completion:
            continue  # reconstructed in the full RGBeta portal basis below
        mapped = aliases[k]
        if mapped.startswith('lambdaT3GeneralQ') or mapped.startswith('lambdaHHdag'):
            covered_extra.add(mapped)
        names.append(f'native:{k}')
        parts.append(raw)
        descriptions.append({'name': k, 'rgbeta': mapped,
                             'origin': 'native Matchete source CG'})

    if include_portal_completion:
        portal_report, portal_arrays = completed_portal_numeric(
            ch, {'muS': dict(state)})
        status = portal_report['scales']['muS']
        require(not status['missing_rgbeta_couplings'],
                f'{model}: incomplete baseline portal values: '
                f"{status['missing_rgbeta_couplings']}")
        portal_parts = np.asarray(portal_arrays['muS_channel_tensors'], float)
        portal_names = list(portal_arrays['muS_channel_labels'])
        require(portal_parts.shape[1:] == (4,4,n,n),
                f'{model}: baseline portal tensor dimension mismatch')
        names.extend(f'portal:{name}' for name in portal_names)
        parts.extend(portal_parts)
        descriptions.extend({'name': name, 'origin': 'exact RGBeta portal completion'}
                            for name in portal_names)

    if include_general_extra:
        extra_names, extra_parts, extra_descriptions = general_extra_channels(
            model, registry, state, excluded=covered_extra, geometry=geometry)
        names.extend(f'general:{name}' for name in extra_names)
        parts.extend(extra_parts)
        descriptions.extend(extra_descriptions)

    require(bool(parts) and len(names) == len(set(names)),
            f'{model}: no or duplicate scalar HHSS channels')
    original = np.stack(parts, axis=0).astype(float)
    if model in SPECIAL:
        special = build_special_basis(model, snapshot, boundary)
        R = special['R_original_from_boundary']
        mass2 = special['mass2_boundary']
        channels = np.einsum('cABij,ia,jb->cABab', original, R, R,
                             optimize=True)
        basis_status = special['mass_basis_assumption']
    else:
        mass2 = np.diag([snapshot.scalar_mass1_sq]*(2*d1) +
                        [snapshot.scalar_mass2_sq]*(2*d2))
        channels = original
        basis_status = 'original field components equal F-boundary components'
    require(channels.shape == (len(names),4,4,n,n) and
            np.isfinite(channels).all(), f'{model}: nonfinite scalar channels')
    require(np.linalg.eigvalsh(mass2).min() > 0,
            f'{model}: nonpositive scalar mass eigenvalue')
    q = channels.sum(axis=0)
    require(np.allclose(q,q.transpose(1,0,2,3),rtol=1e-9,atol=1e-10) and
            np.allclose(q,q.transpose(0,1,3,2),rtol=1e-9,atol=1e-10),
            f'{model}: quartic not Bose symmetric')
    return TensorBundle(model, mass2, q, channels, tuple(names), {
        'model': model, 'assumed_native_coupling_aliases': aliases,
        'quartic_channels': descriptions,
        'source_basis': basis_status,
        'general_extra_registry_count': registry['models'][model]['total_additional_quartics'],
        'full_general_HHSS_reference_invariant_basis_used': include_general_extra,
        'all_renormalisable_scalar_extra_quartics_are_in_intermediate_RGEs': None,
        'matchete_rgbeta_map_certified': False,
        'spinor_normalization_certified': False,
    })


def response(source: Any, bundle: TensorBundle, *, mu_f: float, mu_s: float,
             spinor_factor: complex = 1+0j,
             check_scale_cancellation: bool = True) -> dict:
    """Compute LLSS direct log + simultaneous finite 2-scalar bubble ONCE.

    Unlike adding an EFT1 log to the one-scale Matchete result, this returns a
    separate LLSS sector. For complete physical C5, a separately matched hard
    contribution is required; this routine does not invent one.
    """
    from Numerical.joint_scalar_matching.large_batch import _response
    require(source.model == bundle.model, 'source/scalar model mismatch')
    require(mu_f > mu_s > 0 and abs(mu_f/source.mu_f_gev-1)<1e-9,
            'invalid matching scales')
    factor = complex(spinor_factor)
    require(math.isfinite(factor.real) and math.isfinite(factor.imag),
            'nonfinite spinor convention factor')
    raw = source.raw * factor
    q = bundle.quartic
    whole = _response(raw, q, bundle.mass2, mu_f, mu_s)
    if check_scale_cancellation:
        shifted = _response(raw, q, bundle.mass2, mu_f,
                            mu_s*math.exp(0.03125))
        diff = float(np.linalg.norm(shifted['sum']-whole['sum']))
        scale = max(float(np.linalg.norm(whole['sum'])),
                    float(np.linalg.norm(whole['running'])),
                    float(np.linalg.norm(whole['bubble'])), 1e-29)
        require(diff / scale < 2e-8,
                'fixed-order running/threshold scale dependence does not cancel')
        cancellation = diff / scale
    else:
        cancellation = None
    entries = []
    for name, quartic in zip(bundle.names, bundle.channels):
        each = _response(raw, quartic, bundle.mass2, mu_f, mu_s)
        entries.append({'name': name,
                        'category': 'quartic_decomposition',
                        'c5_leading_log': encode(each['running']),
                        'c5_joint_scalar_bubble': encode(each['bubble']),
                        'c5_total': encode(each['sum']),
                        'log_frobenius': float(np.linalg.norm(each['running']))})
    # The quartic and fermion-tree entries are two alternative decompositions
    # of the same total, not additive contributions across categories.
    for source_index, tensor in zip(source.source_indices, source.by_source):
        each = _response(tensor*factor, q, bundle.mass2, mu_f, mu_s)
        # Source contributions can fail individual Ward projection: that is
        # acceptable; only the sum is required to be Weinberg shaped.
        entries.append({'name': f'F_tree_source:{source_index}',
                        'category': 'fermion_tree_decomposition',
                        'c5_leading_log': encode(each['running']),
                        'c5_joint_scalar_bubble': encode(each['bubble']),
                        'c5_total': encode(each['sum']),
                        'log_frobenius': float(np.linalg.norm(each['running']))})
    return {
        'model': source.model,
        'format': 't3-general-16-reference-basis-leading-one-loop-v1',
        'mu_F_gev': mu_f, 'mu_S_gev': mu_s,
        'spinor_normalization_assumption': {'real': factor.real, 'imag': factor.imag},
        'c5_beta_16pi2': encode(whole['beta']),
        'c5_eft1_leading_log': encode(whole['running']),
        'c5_scalar_finite_bubble': encode(whole['bubble']),
        'c5_LLSS_total': encode(whole['sum']),
        'c5_wilson_frozen_at_F': encode(whole['bubble']),
        'c5_wilson_running_on': encode(whole['sum']),
        'log_to_total_norm_ratio': float(np.linalg.norm(whole['running']) /
                   max(float(np.linalg.norm(whole['sum'])), 1e-40)),
        'weinberg_beta_projection_residual': whole['beta_residual'],
        'weinberg_bubble_projection_residual': whole['bubble_residual'],
        'matching_scale_cancellation_residual': cancellation,
        'per_channel': entries,
        'quartic_provenance': bundle.provenance,
        'computation_scope': 'LLSS_scalar_pair_leading_one_loop_plus_both_scalar_bubble',
        'full_sequential_one_loop_C5_computed': False,
        'hard_fermion_matching_included': False,
        'other_intermediate_operator_classes_included': False,
    }


def anchor_to_one_scale(one_scale_c5: np.ndarray, leading_log_c5: np.ndarray,
                        *, wilson_running_on: bool) -> dict:
    """Ablation-only decomposition, not independent physical hard matching.

    Defines C_residual = C_one_scale - C_log, so that selecting Wilson running
    ON reproduces the supplied one-scale C5 EXACTLY and OFF removes just the
    identified log. This NEVER adds the log twice. The one-scale and the log
    must refer to compatible source inputs for quantitative interpretation.
    """
    full = decode(one_scale_c5) if isinstance(one_scale_c5, Mapping) else np.asarray(one_scale_c5, complex)
    log = decode(leading_log_c5) if isinstance(leading_log_c5, Mapping) else np.asarray(leading_log_c5, complex)
    require(full.shape == log.shape == (3,3) and np.isfinite(full).all()
            and np.isfinite(log).all(), 'invalid C5 reference/log')
    residual = full - log
    result = residual + (log if wilson_running_on else 0)
    return {'c5_residual_matching': encode(residual),
            'c5_at_scalar_threshold': encode(result),
            'source': 'one_scale_minus_log_anchored_ablation_not_new_hard_matching'}
