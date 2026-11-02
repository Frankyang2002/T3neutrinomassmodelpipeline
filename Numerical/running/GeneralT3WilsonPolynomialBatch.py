"""All-16 Matchete ordered LLSS -> symmetric *polynomial* Wilson coordinates.

This is a representation/axis bridge, NOT a Matchete-to-master spinor
normalisation certificate. The formal polynomial is

    P(T) = sum T[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b],

where the psi symbols represent symmetric Weyl bilinears and phi are real
commuting scalar coordinates. Every PL source term is aligned independently
against its original Matchete lepton/flavour dummies *before* terms are summed.

The independent orbit coefficient is (sum of source entries in that orbit) /
number of distinct ordered orbit members (1, 2, or 4). Expanding the result
as a separately symmetric tensor reconstructs P exactly, without discarding
antisymmetric physical operators or introducing a guessed EFT normalisation.

No Lagrangian, boundary or RGE files are modified by this module.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import scalar_basis
from Numerical.running.GeneralT3WilsonEvolution import (
    FORMAT as ORDERED_FORMAT, STATUS as ORDERED_STATUS,
    boundary_coordinates, read_json,
)

FORMAT = 't3-general-16-polynomial-wilson-bridge-v1'
STATUS = 'PolynomialRepresentativeVerified_NotSpinorNormalisation'
Key = tuple[int, int, int, int, int, int]  # SU2, flavour, SU2, flavour, real scalar, real scalar

# Flavour indices are attached to the lepton fields, not inferred from the
# order of Yukawa matrices. This pattern also matches nested NCM/Transp/Bar.
_LEPTON = re.compile(
    r'Field\[l,\s*Fermion,\s*\{\s*'
    r'Index\[(?P<su2>[^,\]]+),\s*SU2L\[fund\]\],\s*'
    r'Index\[(?P<flavour>[^,\]]+),\s*Flavor\]'
)


def _require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def _finite_matrix(raw: Mapping[str, Any]) -> np.ndarray:
    _require(isinstance(raw, Mapping) and 'real' in raw and 'imag' in raw,
             'Missing finite complex 3x3 flavour matrix')
    a = np.asarray(raw['real'], dtype=float) + 1j*np.asarray(raw['imag'], dtype=float)
    _require(a.shape == (3, 3) and bool(np.isfinite(a).all()),
             'Expected finite complex 3x3 flavour matrix')
    return a


def _digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def _parse_leptons(template: str) -> tuple[tuple[str, str], tuple[str, str]]:
    _require(isinstance(template, str) and 'GammaCC' in template,
             'Missing exact Matchete LLSS template')
    matched = tuple((m.group('su2').strip(), m.group('flavour').strip())
                    for m in _LEPTON.finditer(template))
    _require(len(matched) == 2, 'Expected exactly two external Matchete lepton fields')
    _require(len({v[0] for v in matched}) == 2 and
             len({v[1] for v in matched}) == 2,
             'External lepton SU2/flavour dummy labels are not distinguishable')
    return matched[0], matched[1]


def source_term_alignment(operator: Mapping[str, Any],
                          source_term: Mapping[str, Any]) -> dict:
    """Verify permutation against exact lepton-field order and Yukawa axes."""
    _require(operator.get('source_index') == source_term.get('source_index'),
             'Source term index does not match ordered tensor')
    axes = source_term.get('flavour_matrix_axes_external_indices')
    _require(isinstance(axes, list) and len(axes) == 2 and
             all(isinstance(x, str) for x in axes) and axes[0] != axes[1],
             'Missing distinct numerical flavour matrix axes')
    _require(operator.get('flavour_matrix_axes_external_indices') == axes,
             'Ordered tensor and original source flavour axes disagree')
    factors = source_term.get('factors')
    _require(isinstance(factors, list) and len(factors) == 2 and
             [f.get('external_index') for f in factors] == axes,
             'Source Yukawa factors do not reproduce matrix axis metadata')
    legs = _parse_leptons(source_term.get('operator_template_unit_flavor_inputform'))
    try:
        permutation = tuple(axes.index(flav) for _, flav in legs)
    except ValueError as exc:
        raise ValueError('Lepton flavour label absent from Yukawa axes') from exc
    _require(sorted(permutation) == [0, 1], 'Lepton-to-flavour axis mapping ambiguous')
    # The projector assigns its L_i,L_j in *first-occurrence* lepton order.
    # Reject changed templates rather than silently permute the wrong legs.
    _require(len({su2 for su2, _ in legs}) == 2,
             'Projector SU2 dummy order ambiguous')
    return {
        'source_index': int(operator['source_index']),
        'lepton_field_dummies': [
            {'su2': s, 'flavour': f} for s, f in legs
        ],
        'matrix_flavour_axes': list(axes),
        'flavour_axis_permutation_for_lepton_order': list(permutation),
    }


def load_source_terms(model: str, ordered: Mapping[str, Any],
                      source: Mapping[str, Any]) -> dict[int, Mapping[str, Any]]:
    """Source metadata only, no stale or mismatched operator inventory."""
    _require(source.get('model') == model, 'Original source model mismatch')
    _require(isinstance(source.get('mu_F_gev'), (int,float)) and
             float(source['mu_F_gev']) == float(ordered['mu_F_gev']),
             'Source and projected tensor fermion threshold scales disagree')
    if model in ('B_-1', 'C_-1'):
        _require(source.get('format') == 't3-mixed-factorised-LLSS-source-v1' and
                 source.get('status') == 'NumericalOrderedSourceFlavourCGExported',
                 'Unverified mixed native-Matchete source inventory')
        # Both JSON artifacts originate in the same successful matching run.
        # A changed U/LLSS export must not be combined with stale CG tensors.
        for old_key, source_key in (
            ('source_status_sha256', 'source_status_sha256'),
            ('source_mass_basis_sha256', 'source_mass_basis_sha256'),
            ('source_native_fresh_llss_sha256', 'llss_sha256'),
        ):
            _require(ordered.get(old_key) == source.get(source_key) and
                     isinstance(ordered.get(old_key), str) and
                     len(ordered[old_key]) == 64,
                     'Mixed source and projected tensor provenance mismatch: ' + old_key)
        terms = [t for t in source.get('operator_terms', [])
                 if t.get('chirality') == 'PL']
        _require(bool(terms) and len(source['operator_terms']) == 2 * len(terms),
                 'Incomplete mixed PL/PR source inventory')
    else:
        _require(source.get('format') == 't3-source-exact-F-tree-flavour-operator-handoff-v1'
                 and source.get('status') ==
                 'SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor',
                 'Unverified ordinary Matchete flavour boundary')
        terms = source.get('operator_terms')
    _require(isinstance(terms, list) and len(terms) == ordered.get('tree_operator_count_PL'),
             'PL source count does not match ordered tensors')
    out: dict[int, Mapping[str, Any]] = {}
    for t in terms:
        ix = t.get('source_index')
        _require(isinstance(ix, int) and ix > 0 and ix not in out,
                 'Duplicate or invalid original Matchete source index')
        out[ix] = t
    _require({op.get('source_index') for op in ordered['operators']} == set(out),
             'Original and projected PL source indices disagree')
    return out


def ordered_aligned_components(model: str, ordered: Mapping[str, Any],
                               source: Mapping[str, Any]) -> tuple[dict[Key, complex], list[dict], dict]:
    """Reassemble aligned three-flavour tensor; rekey BEFORE source summation."""
    _, _, info = boundary_coordinates(ordered, model)  # full v3 validation
    terms = load_source_terms(model, ordered, source)
    raw: dict[Key, complex] = defaultdict(complex)
    alignments = []
    component_count = 0
    for operator in ordered['operators']:
        alignment = source_term_alignment(operator, terms[operator['source_index']])
        perm = alignment['flavour_axis_permutation_for_lepton_order']
        alignments.append(alignment)
        component_count += len(operator['components'])
        for entry in operator['components']:
            i, j = entry['lepton_su2']
            a, b = entry['scalar_real']
            matrix = _finite_matrix(entry['ordered_flavour_coefficient_gev_inv'])
            for p in range(3):
                for q in range(3):
                    z = complex(matrix[p, q])
                    if z == 0j:
                        continue
                    flav = (p + 1, q + 1)
                    raw[(i, flav[perm[0]], j, flav[perm[1]], a, b)] += z
    _require(bool(raw), 'All source Wilson components are zero')
    # Keep actual initial support; exact zero cancellation is allowed.
    aligned = {k: v for k, v in raw.items() if v != 0j}
    _require(bool(aligned), 'All source Wilson components cancel')
    return aligned, alignments, {**info,
        'source_component_count': component_count,
        'source_nonzero_aligned_coordinates': len(aligned)}


def canonical_key(key: Key) -> Key:
    i, p, j, q, a, b = key
    _require(i in (1,2) and j in (1,2) and p in (1,2,3) and q in (1,2,3),
             'Invalid lepton SU2/flavour coordinate')
    _require(a > 0 and b > 0, 'Invalid scalar coordinate')
    if (i, p) > (j, q):
        i, p, j, q = j, q, i, p
    if a > b:
        a, b = b, a
    return (i,p,j,q,a,b)


def distinct_orbit(key: Key) -> tuple[Key, ...]:
    k = canonical_key(key)
    i,p,j,q,a,b = k
    members = {(i,p,j,q,a,b), (j,q,i,p,a,b),
               (i,p,j,q,b,a), (j,q,i,p,b,a)}
    return tuple(sorted(members))


def canonical_polynomial(tensor: Mapping[Key, complex]) -> tuple[dict[Key, complex], dict[Key, complex]]:
    """Return independent canonical tensor and actual polynomial monomial sums.

    If a monomial has orbit size m, C_sym=value/m on each distinct ordered
    coordinate. The output carries only one representative plus orbit size.
    """
    grouped: dict[Key, complex] = defaultdict(complex)
    for k, v in tensor.items():
        _require(len(k) == 6 and math.isfinite(complex(v).real) and
                 math.isfinite(complex(v).imag), 'Invalid finite aligned tensor')
        grouped[canonical_key(k)] += complex(v)
    monomials = {k:z for k,z in grouped.items() if z != 0j}
    independent = {k:z/len(distinct_orbit(k)) for k,z in monomials.items()}
    return independent, monomials


def polynomial_identity(tensor: Mapping[Key, complex],
                        independent: Mapping[Key, complex]) -> dict:
    """Compare all source monomials with expansion of the symmetric tensor."""
    reconstructed: dict[Key, complex] = defaultdict(complex)
    for key, z in independent.items():
        _require(key == canonical_key(key), 'Noncanonical independent coordinate')
        for member in distinct_orbit(key):
            reconstructed[canonical_key(member)] += complex(z)
    original_grouped: dict[Key, complex] = defaultdict(complex)
    for key, value in tensor.items():
        original_grouped[canonical_key(key)] += complex(value)
    keys = set(original_grouped) | set(reconstructed)
    delta = np.asarray([original_grouped.get(k,0j) - reconstructed.get(k,0j)
                        for k in sorted(keys)], dtype=complex)
    norm = np.linalg.norm(list(original_grouped.values()))
    err = np.linalg.norm(delta)
    return {
        'source_monomial_count': sum(v != 0 for v in original_grouped.values()),
        'symmetric_monomial_count': sum(v != 0 for v in reconstructed.values()),
        'max_absolute_monomial_difference_gev_inv':float(max((abs(v) for v in delta), default=0.0)),
        'relative_monomial_difference':float(err / max(float(norm), 1e-300)),
        'preserved_to_numerical_precision': bool(err <= 3e-12 * max(float(norm), 1e-300)),
    }


def encode_independent_tensor(independent: Mapping[Key, complex],
                              dimension: int) -> list[dict]:
    rows = []
    for k, z in sorted(independent.items()):
        _require(k[4] <= dimension and k[5] <= dimension,
                 'Scalar index exceeds real model dimension')
        rows.append({'coordinate':list(k),
                     'coefficient_real':float(z.real), 'coefficient_imag':float(z.imag),
                     'orbit_size':len(distinct_orbit(k))})
    return rows


def build_polynomial_bridge(model: str, ordered: Mapping[str, Any],
                            source: Mapping[str, Any]) -> dict:
    """Build model-independent symmetric representative with fail-closed audit."""
    aligned, alignments, info = ordered_aligned_components(model, ordered, source)
    independent, monomials = canonical_polynomial(aligned)
    proof = polynomial_identity(aligned, independent)
    _require(proof['preserved_to_numerical_precision'],
             'Symmetric representative changed source polynomial')
    dim = scalar_basis(model).total_dimension
    basis = scalar_basis(model)
    return {
        'format':FORMAT, 'model':model, 'status':STATUS,
        'mu_F_gev':info['mu_F_gev'],
        'real_scalar_dimension':dim,
        'scalar_field_offsets_one_based':dict(basis.offsets),
        'scalar_field_complex_dimensions':dict(basis.dimensions),
        'scalar_real_fields':sorted(basis.real_fields),
        'source_PL_term_count':len(alignments),
        'source_ordered_component_count':info['ordered_real_component_count'],
        'nonzero_aligned_flavour_coordinate_count':len(aligned),
        'independent_symmetric_coordinate_count':len(independent),
        'lepton_flavour_alignments':alignments,
        'flavour_permutation_histogram':{
            '[0,1]':sum(a['flavour_axis_permutation_for_lepton_order']==[0,1] for a in alignments),
            '[1,0]':sum(a['flavour_axis_permutation_for_lepton_order']==[1,0] for a in alignments),
        },
        'polynomial_identity':proof,
        'formal_polynomial_convention':'sum C[i,p,j,q,a,b]*psi[i,p]*psi[j,q]*phi[a]*phi[b]',
        'independent_coefficient_convention':'Each entry is one fully symmetric tensor coordinate; reconstruct all orbit members with equal coefficient',
        'independent_symmetric_coordinates':encode_independent_tensor(independent,dim),
        'overall_matchete_spinor_normalization_certified':False,
        'matchete_fierz_eom_mapping_certified':False,
        'physical_master_rge_basis_identification_certified':False,
        'physical_three_family_eft1_gamma_constructed':False,
        'eft1_wilson_running_performed':False,
        'physical_sequential_c5_certified':False,
    }


def master_weinberg_polynomial_check() -> dict:
    """An independent internal master convention check, not a Matchete check."""
    import sympy as sp
    from RGE.matching.WeinbergTensorAdapter import (
        build_sm_eft, weinberg_polynomial, build_weinberg_wilson_tensor,
    )
    scalar, fermion = build_sm_eft()
    original, psi, phi = weinberg_polynomial(scalar, fermion, sp.S.One)
    tensor = build_weinberg_wilson_tensor(scalar, fermion, sp.S.One)
    restored = sum((z*psi[i]*psi[j]*phi[a]*phi[b]
                    for (i,j,a,b),z in tensor.items()), sp.S.Zero)
    residual = sp.expand(original - restored)
    _require(sp.simplify(residual) == 0, 'Master Weinberg tensor normalization changed')
    return {
        'master_weinberg_polynomial_residual':'0',
        'reconstruction_without_additional_polynomial_factor':True,
        'matchete_to_master_operator_normalization_certified':False,
    }
