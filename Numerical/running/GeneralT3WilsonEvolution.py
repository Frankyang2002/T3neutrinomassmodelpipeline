"""Read-only all-16 Wilson boundary handoff; conditional explicit-Gamma transport.

The ordered Matchete source contributions are NOT a certified physical,
operator-symmetrized C_ijab basis.  This module deliberately does not infer
an anomalous-dimension matrix or a Matchete/RGBeta coupling dictionary.

If an independently validated complete gamma is available in the exact ordered
source-coordinate basis (including zero-boundary closure coordinates), numerical
transport is possible, but is reported as CONDITIONAL rather than physical EFT1.
"""
from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.running.GeneralT3LinearWilsonFlow import LOOP, validate_gamma
from Numerical.running.GeneralT3OrderedTensor import scalar_basis

FORMAT = "t3-ordered-source-flavour-and-real-scalar-tensor-v1"
STATUS = "OrderedNumericalRealComponentTensorFromMatchete"
GAMMA_FORMAT = "t3-explicit-ordered-wilson-gamma-v1"
MODELS = tuple(f"{family}_{a:+d}" for family, alphas in {
    'A':(-4,-2,0), 'B':(-3,-1,1), 'C':(-3,-1,1),
    'D':(-2,0,2), 'E':(-4,-2,0,2)
}.items() for a in alphas)


def tag(model: str) -> str:
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}")
    return model.replace('+', 'p').replace('-', 'm')


def read_json(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return obj


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def save_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _finite_matrix(raw: Any) -> np.ndarray:
    if not isinstance(raw, Mapping):
        raise ValueError('Missing complex flavour matrix')
    try:
        real = np.asarray(raw['real'], dtype=float)
        imag = np.asarray(raw['imag'], dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('Malformed complex flavour matrix') from exc
    if real.shape != (3,3) or imag.shape != (3,3):
        raise ValueError('Flavour matrix must be 3x3')
    result = real + 1j*imag
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite flavour matrix')
    return result


def boundary_coordinates(payload: Mapping[str, Any], model: str) -> tuple[list[str], np.ndarray, dict]:
    """Source-order-preserving coordinates. No physical symmetrization."""
    if payload.get('model') != model or payload.get('format') not in (
        FORMAT, 't3-mixed-native-matchete-ordered-tensor-v1'
    ) or payload.get('status') != STATUS:
        raise ValueError(f'{model}: invalid v3 ordered tensor format/status')
    if (payload.get('scalar_basis_explicit') is not True or
        payload.get('matchete_cg_contractions_retained') is not True or
        payload.get('external_three_family_flavour_order_not_symmetrized') is not True):
        raise ValueError(f'{model}: scalar, CG or flavour provenance not verified')
    if payload.get('canonical_C_ijab_normalization_certified') is not False:
        raise ValueError('Unexpected canonical normalization claim or missing flag')
    mu_f = float(payload.get('mu_F_gev', float('nan')))
    if not math.isfinite(mu_f) or mu_f <= 0:
        raise ValueError('Invalid fermion threshold scale')
    basis = scalar_basis(model)
    operators = payload.get('operators')
    if not isinstance(operators, list) or not operators or len(operators) != payload.get('tree_operator_count_PL'):
        raise ValueError('Incorrect PL operator count')
    labels: list[str] = []
    values: list[complex] = []
    seen_terms = set()
    component_count = 0
    for operator in operators:
        source = operator.get('source_index')
        if not isinstance(source, int) or source <= 0 or source in seen_terms:
            raise ValueError('Invalid or duplicate source operator index')
        seen_terms.add(source)
        if operator.get('ordered_not_canonical_symmetrized') is not True:
            raise ValueError('Missing ordered operator metadata')
        if operator.get('real_scalar_dimension_including_higgs') != basis.total_dimension:
            raise ValueError('Scalar real-coordinate basis mismatch')
        if operator.get('basis_scalar_fields') != basis.dimensions:
            raise ValueError('Scalar species/dimension mismatch')
        components = operator.get('components')
        if not isinstance(components, list) or not components or len(components) != operator.get('component_count'):
            raise ValueError('Invalid real-component count')
        component_count += len(components)
        seen_component = set()
        for component in components:
            i,j = component['lepton_su2']
            a,b = component['scalar_real']
            key = (i,j,a,b)
            if key in seen_component or i not in (1,2) or j not in (1,2) or not (
                isinstance(a,int) and isinstance(b,int) and 1<=a<=basis.total_dimension and 1<=b<=basis.total_dimension
            ):
                raise ValueError('Invalid or duplicate ordered SU2/scalar component')
            seen_component.add(key)
            m = _finite_matrix(component['ordered_flavour_coefficient_gev_inv'])
            for p in range(3):
                for q in range(3):
                    labels.append(f'term{source}:L{i}{j}:S{a},{b}:F{p+1}{q+1}')
                    values.append(m[p,q])
    if component_count != payload.get('total_ordered_real_component_count') or len(labels) != len(set(labels)):
        raise ValueError('Inconsistent total tensor component count')
    return labels, np.asarray(values, dtype=np.complex128), {
        'model':model, 'mu_F_gev':mu_f,
        'tree_operator_count_PL':len(operators),
        'ordered_real_component_count':component_count,
        'nonzero_ordered_flavour_coordinates':len(labels),
        'scalar_real_dimension':basis.total_dimension,
        'boundary_l2_norm_gev_inv':float(np.linalg.norm(values)),
        'coordinate_basis':'ordered_Matchete_real_components_with_flavour',
        'canonical_C_ijab_normalization_certified':False,
        'physical_eft1_wilson_rge_certified':False,
    }



def aggregate_source_coordinates(
    source_labels: list[str], source_values: np.ndarray
) -> tuple[list[str], np.ndarray]:
    """Sum distinct Matchete source terms onto identical ordered LLSS coordinates.

    This is linear addition only. It does not impose i<->j or a<->b
    symmetry, infer physical operator normalization, or discard conjugates.
    A term index is not a physical operator index; retaining it in the gamma
    basis would incorrectly keep separate copies of one operator component.
    """
    if source_values.shape != (len(source_labels),) or len(source_labels) != len(set(source_labels)):
        raise ValueError('Malformed independent source-term coordinates')
    merged: dict[str, complex] = {}
    for label, coefficient in zip(source_labels, source_values):
        if not label.startswith('term') or ':' not in label:
            raise ValueError('Malformed term-labelled coordinate')
        _, physical = label.split(':', 1)
        if not physical.startswith('L') or ':S' not in physical or ':F' not in physical:
            raise ValueError('Malformed ordered physical coordinate')
        key = physical
        merged[key] = merged.get(key, 0j) + coefficient
    coords = sorted(merged)
    return coords, np.asarray([merged[k] for k in coords], dtype=np.complex128)

def verify_batch(root: Path) -> dict[str, dict]:
    """Require exact 16-model v3 inventory; never accept stale/partial summaries."""
    root = Path(root)
    summary = read_json(root/'summary.json')
    rows = summary.get('results')
    if (summary.get('format') != FORMAT or
        summary.get('ordered_real_tensor_count') != 16 or
        summary.get('blocked_count') != 0 or not isinstance(rows, list) or len(rows) != 16):
        raise ValueError('Expected successful complete 16-model v3 batch')
    if len({row.get('model') for row in rows}) != 16 or {row.get('model') for row in rows} != set(MODELS):
        raise ValueError('Incomplete or duplicate v3 model inventory')
    result = {}
    for row in rows:
        model = row['model']
        path = root/tag(model)/'ordered_tensor.json'
        if row.get('status') != STATUS or not path.is_file():
            raise ValueError(f'{model}: ordered_tensor.json absent or blocked')
        if row.get('tensor_sha256') and row['tensor_sha256'] != digest(path):
            raise ValueError(f'{model}: boundary digest changed')
        raw = read_json(path)
        coordinates, initial, report = boundary_coordinates(raw,model)
        if (row.get('tree_operator_count_PL') != report['tree_operator_count_PL'] or
            row.get('real_component_count') != report['ordered_real_component_count']):
            raise ValueError(f'{model}: summary count disagrees with tensor')
        result[model] = {**report, 'tensor_sha256':digest(path),
                         'tensor_file':str(path), 'coordinates':coordinates,
                         'initial':initial}
    return result


def extend_to_gamma_basis(
    initial_labels: list[str], initial_values: np.ndarray, gamma: Mapping[str, Any]
) -> tuple[list[str], np.ndarray, int]:
    """Embed tree support in a gamma-complete basis, leaving new modes zero.

    Crucial: the gamma source must explicitly enumerate *all* outputs which it
    can mix into, even when those coefficients vanish at the matching boundary.
    Merely evolving the sparse support is NOT closed under an RGE.
    """
    full = gamma.get('coordinates')
    if not isinstance(full, list) or len(full) != len(set(full)) or not all(
        isinstance(k, str) and k for k in full
    ):
        raise ValueError('Gamma coordinate inventory invalid')
    if len(initial_labels) != len(set(initial_labels)) or initial_values.shape != (len(initial_labels),):
        raise ValueError('Invalid ordered matching boundary vector')
    if not set(initial_labels).issubset(full):
        raise ValueError('Gamma basis omits source boundary components')
    if gamma.get('closure_includes_zero_boundary_components') is not True:
        raise ValueError('Missing explicit closure inventory declaration')
    positions = {coord:i for i,coord in enumerate(full)}
    init = np.zeros(len(full),dtype=complex)
    for label,z in zip(initial_labels,initial_values):
        init[positions[label]] = z
    return full,init,len(full)-len(initial_labels)


def evolve_given_gamma(
    raw: Mapping[str,Any], model: str, tensor_sha256: str,
    gamma: Mapping[str, Any], *, npoints:int=41, max_dim:int=2500
) -> dict:
    """Numerical evolution only with an explicitly provided gamma source.

    This is *conditional* transport: no physical C_ijab interpretation is made.
    """
    if not isinstance(npoints,int) or not 2<=npoints<=501:
        raise ValueError('Require 2..501 trajectory points')
    source_labels,source_values,report=boundary_coordinates(raw, model)
    labels,base=aggregate_source_coordinates(source_labels,source_values)
    if gamma.get('format') != GAMMA_FORMAT or gamma.get('model') != model:
        raise ValueError('Gamma model/format mismatch')
    if gamma.get('ordered_tensor_sha256') != tensor_sha256:
        raise ValueError('Gamma source not bound to this ordered tensor SHA-256')
    coordinates,c0,generated = extend_to_gamma_basis(labels,base,gamma)
    if len(coordinates)>max_dim or max_dim<1:
        raise ValueError('Gamma dimension exceeds explicit --max-gamma-dim')
    mu_f=report['mu_F_gev']
    mu_s=float(gamma.get('mu_S_gev',float('nan')))
    if not (math.isfinite(mu_s) and 0<mu_s<mu_f):
        raise ValueError('Invalid scalar threshold')
    if 'mu_F_gev' in gamma and float(gamma['mu_F_gev']) != mu_f:
        raise ValueError('Gamma mu_F inconsistent with tensor')
    logs,mats=validate_gamma(gamma,coordinates,mu_f,mu_s)
    start=math.log(mu_f); end=math.log(mu_s)
    def rhs(t, state):
        hi=min(max(int(np.searchsorted(logs,t,side='right')),1),len(logs)-1)
        alpha=(t-logs[hi-1])/(logs[hi]-logs[hi-1])
        return LOOP*((mats[hi-1]*(1-alpha)+mats[hi]*alpha)@state)
    eval_times=np.linspace(start,end,npoints)
    sol=solve_ivp(rhs,(start,end),c0,t_eval=eval_times,rtol=1e-8,atol=1e-11,
                  method='DOP853')
    if not sol.success or sol.y.shape != (len(coordinates),npoints) or not np.isfinite(sol.y).all():
        raise RuntimeError('Conditional gamma solver failed: '+sol.message)
    # Retain all trajectory points, including source-zero modes generated by RG mixing.
    return {
        'format':'t3-conditional-explicit-ordered-source-gamma-trajectory-v1',
        'model':model,'status':'ConditionalExplicitGammaNumericalTransport',
        'ordered_tensor_sha256':tensor_sha256,
        'coordinates':coordinates,'tree_support_coordinates':len(labels),
        'matched_source_flavour_coordinates':len(source_labels),
        'closure_zero_initial_coordinates':generated,
        'equation':'16pi2*dC/dlnmu=Gamma(mu)@C',
        'mu_F_gev':mu_f,'mu_S_gev':mu_s,
        'scales_gev':np.exp(sol.t).tolist(),
        'trajectory_real':sol.y.T.real.tolist(),
        'trajectory_imag':sol.y.T.imag.tolist(),
        'nfev':int(sol.nfev),
        'relative_change':float(np.linalg.norm(sol.y[:,-1]-c0)/max(np.linalg.norm(c0),1e-300)),
        'physical_eft1_wilson_rge_certified':False,
        'canonical_C_ijab_normalization_certified':False,
        'gamma_derived_by_this_module':False,
        'scalar_threshold_matching_performed':False,
    }
