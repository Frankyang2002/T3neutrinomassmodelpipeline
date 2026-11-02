"""A_-4 ordered Matchete LLSS coefficient -> symmetric polynomial representative.

This is an algebraic change of coefficient *representation*, not a matching
normalization certificate.  In particular, the repository's source matrix is
multiplied by its Matchete tensor/fields as exported; the relative spinor,
Majorana, EOM and Lagrangian sign conventions have not been independently
matched to the master Wilson RGE.

For formal COMMUTING psi placeholders and REAL scalar coordinates phi,

  P(T) = sum_{i,p,j,q,a,b} T[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b].

The simultaneous independent exchanges (i,p)<->(j,q) and a<->b define
an exact Reynolds projection S(T), with P(S(T)) = P(T).  This remains true
for complex tensor coefficients. This module never modifies the ordered input.
"""
from __future__ import annotations

from collections import defaultdict
import re
import math
from typing import Mapping

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import scalar_basis
from Numerical.running.GeneralT3WilsonEvolution import (
    boundary_coordinates, aggregate_source_coordinates,
)
from Numerical.running.GeneralT3WilsonRGEPilot import read_flavour_components

MODEL = "A_-4"
FORMAT = "t3-a-minus4-symmetric-polynomial-bridge-v1"
Key = tuple[int, int, int, int, int, int]


def swap_fermion(key: Key) -> Key:
    i, p, j, q, a, b = key
    return (j, q, i, p, a, b)


def swap_scalar(key: Key) -> Key:
    i, p, j, q, a, b = key
    return (i, p, j, q, b, a)


def orbit(key: Key) -> tuple[Key, ...]:
    flipped = swap_fermion(key)
    return (key, flipped, swap_scalar(key), swap_scalar(flipped))


def validate_a4_components(tensor: Mapping[Key, complex]) -> None:
    if not tensor:
        raise ValueError("A_-4 source tensor empty")
    basis = scalar_basis(MODEL)
    if (basis.dimensions.get("NewScalar1") != 1 or
        basis.dimensions.get("NewScalar2") != 3):
        raise RuntimeError("A_-4 scalar representation changed")
    s1 = set(range(basis.offsets["NewScalar1"],
                   basis.offsets["NewScalar1"] + 2))
    s2 = set(range(basis.offsets["NewScalar2"],
                   basis.offsets["NewScalar2"] + 6))
    for k, z in tensor.items():
        if not isinstance(k, tuple) or len(k) != 6:
            raise ValueError("Expected (SU2,flavour,SU2,flavour,scalar,scalar) key")
        i,p,j,q,a,b = k
        if (i not in (1,2) or j not in (1,2) or
            p not in (1,2,3) or q not in (1,2,3) or
            not ((a in s1 and b in s2) or (a in s2 and b in s1))):
            raise ValueError(f"A_-4 non-LLSS or wrong scalar-field pair: {k}")
        if not (math.isfinite(complex(z).real) and math.isfinite(complex(z).imag)):
            raise ValueError(f"Nonfinite source coefficient: {k}")


def symmetrise_polynomial(tensor: Mapping[Key,complex]) -> dict[Key,complex]:
    """Project onto independent fermion and scalar exchange fixed points.

    This is the *unique orthogonal symmetric representative* of the input
    polynomial under the coordinate-sum convention in the module docstring.
    A missing ordered entry is treated as exactly zero.
    """
    validate_a4_components(tensor)
    keys = {member for k in tensor for member in orbit(k)}
    result = {}
    for k in keys:
        val = sum((complex(tensor.get(other,0j)) for other in orbit(k)),0j)/4
        if val != 0j:
            result[k] = val
    for k, z in result.items():
        if (result.get(swap_scalar(k),0j) != z or
            result.get(swap_fermion(k),0j) != z):
            raise RuntimeError("Polynomial representative not separately symmetric")
    return dict(sorted(result.items()))


def monomial_coefficients(tensor: Mapping[Key,complex]) -> dict[tuple,complex]:
    """Exactly collect coefficients of formal psi_i psi_j phi_a phi_b monomials."""
    result=defaultdict(complex)
    for (i,p,j,q,a,b),z in tensor.items():
        key=(tuple(sorted(((i,p),(j,q)))),tuple(sorted((a,b))))
        result[key]+=complex(z)
    return {k:v for k,v in result.items() if v!=0j}


def polynomial_residual(source: Mapping[Key,complex],
                        representative: Mapping[Key,complex]) -> dict:
    """Identity check in the commuting-placeholder polynomial quotient.

    Test *all* occupied monomials, not only a few random field vectors.
    """
    left=monomial_coefficients(source)
    right=monomial_coefficients(representative)
    all_keys=left.keys()|right.keys()
    error=np.linalg.norm([left.get(k,0j)-right.get(k,0j) for k in all_keys])
    scale=np.linalg.norm(list(left.values()))
    return {"source_polynomial_monomial_count":len(left),
            "projected_polynomial_monomial_count":len(right),
            "max_absolute_monomial_difference_gev_inv":float(max(
                (abs(left.get(k,0j)-right.get(k,0j)) for k in all_keys),default=0.)),
            "relative_monomial_difference":float(error/max(scale,1e-300)),
            "preserved_to_numerical_precision":bool(error<=2e-13*max(scale,1e-300))}


def master_weinberg_polynomial_check() -> dict:
    """Inspect a *live* master-RGE Weinberg example, not a retyped formula.

    Its tensor is generated by four derivatives divided by four.  Verify
    reconstructing P = sum C_ijab psi_i psi_j phi_a phi_b exactly yields
    the example's original polynomial; this fixes the internal polynomial
    coordinate convention, NOT Matchete's overall matching normalization.
    """
    import sympy as sp
    from RGE.matching.WeinbergTensorAdapter import (
        build_sm_eft, weinberg_polynomial, build_weinberg_wilson_tensor,
    )
    scalar,fermion=build_sm_eft()
    poly,psi,phi=weinberg_polynomial(scalar,fermion,sp.Integer(1))
    entries=build_weinberg_wilson_tensor(scalar,fermion,sp.Integer(1))
    reconstructed=sum((v*psi[i]*psi[j]*phi[a]*phi[b]
                       for (i,j,a,b),v in entries.items()),sp.S.Zero)
    residual=sp.simplify(sp.expand(poly-reconstructed))
    if residual!=0:
        raise RuntimeError("Master Weinberg tensor does not reconstruct its source polynomial")
    return {"checked_source":"RGE.matching.WeinbergTensorAdapter.build_weinberg_wilson_tensor",
            "master_reconstruction_residual":"0",
            "master_tensor_reconstructs_polynomial_without_extra_factor":True,
            "weinberg_example_component_count":len(entries),
            "matchete_spinor_normalization_certified":False}


_LEPTON_FIELD = re.compile(
    r"Field\[l,\s*Fermion,\s*\{Index\[(?P<su2>[^,\]]+),\s*SU2L\[fund\]\],"
    r"\s*Index\[(?P<flavour>[^,\]]+),\s*Flavor\]"
)


def a4_lepton_dummy_alignment(payload: Mapping, boundary: Mapping) -> dict:
    """Derive the actual (SU2 component, family index) association.

    _find_lepton_dummies() follows FIRST Field[l] occurrence; a numerical
    Matchete flavour matrix follows the external Yukawa factor axes. Those
    orders are independent. For the checked A_-4 export they are reversed.
    """
    if (boundary.get('model') != MODEL or
        boundary.get('format') != 't3-source-exact-F-tree-flavour-operator-handoff-v1' or
        boundary.get('status') != 'SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor'):
        raise ValueError('Not an A_-4 source-factorised Matchete boundary')
    terms = boundary.get('operator_terms')
    ops=payload.get('operators')
    if not isinstance(terms,list) or len(terms)!=1 or not isinstance(ops,list) or len(ops)!=1:
        raise ValueError('A_-4 expected one Matchete PL source term')
    operator=ops[0]
    term=terms[0]
    if operator['source_index'] != term.get('source_index'):
        raise ValueError('Matchete source index does not match ordered tensor')
    axes=term.get('flavour_matrix_axes_external_indices')
    if (not isinstance(axes,list) or len(axes)!=2 or len(set(axes))!=2 or
        operator.get('flavour_matrix_axes_external_indices')!=axes):
        raise ValueError('Unmatched Matchete numerical flavour axes')
    template=term.get('operator_template_unit_flavor_inputform')
    if not isinstance(template,str):
        raise ValueError('No original Matchete lepton-field template')
    legs=[(m.group('su2').strip(),m.group('flavour').strip())
          for m in _LEPTON_FIELD.finditer(template)]
    if len(legs)!=2 or len(set(legs))!=2:
        raise ValueError('Expected two distinguishable Matchete lepton fields')
    if len({su2 for su2,_ in legs})!=2:
        raise ValueError('Indistinguishable SU(2) lepton dummy ordering')
    if any(flav not in axes for _,flav in legs):
        raise ValueError('Lepton flavour dummy not in numerical matrix axes')
    permutation=[axes.index(flav) for _,flav in legs]
    if sorted(permutation)!=[0,1]:
        raise ValueError('Flavour axes do not uniquely attach to the two leptons')
    factors=term.get('factors')
    if (not isinstance(factors,list) or len(factors)!=2 or
        [factor.get('external_index') for factor in factors]!=axes):
        raise ValueError('Numerical flavour axes do not track source Yukawa factors')
    return {'lepton_fields_in_tensor_su2_order':[
                {'su2_dummy':su2,'flavour_dummy':flav} for su2,flav in legs],
            'matrix_flavour_axis_dummies':axes,
            'flavour_axis_for_lepton_field_order_zero_based':permutation,
            'tensor_su2_and_flavour_orders_differ':permutation!=[0,1],
            'original_matchete_flavour_dummies_checked':True}


def _encode_entry(k:Key,z:complex) -> dict:
    return {"index":[int(x) for x in k],"real":float(z.real),"imag":float(z.imag)}


def bridge_a4(payload: Mapping, source_boundary: Mapping, *,check_master:bool=True) -> dict:
    alignment=a4_lepton_dummy_alignment(payload,source_boundary)
    labels,values,info=boundary_coordinates(payload,MODEL)
    merged_labels,merged_values=aggregate_source_coordinates(labels,values)
    unpaired=read_flavour_components(merged_labels,merged_values,MODEL)
    perm=alignment['flavour_axis_for_lepton_field_order_zero_based']
    original={}
    for (i,p,j,q,a,b),z in unpaired.items():
        flavours=(p,q)
        key=(i,flavours[perm[0]],j,flavours[perm[1]],a,b)
        if key in original:
            raise RuntimeError('Non-unique corrected lepton/flavour coordinate')
        original[key]=z
    validate_a4_components(original)
    representative=symmetrise_polynomial(original)
    proof=polynomial_residual(original,representative)
    if not proof["preserved_to_numerical_precision"]:
        raise RuntimeError("Fermion/scalar symmetric representative changes LLSS polynomial")
    norm_raw=float(np.linalg.norm(list(original.values())))
    norm_proj=float(np.linalg.norm(list(representative.values())))
    diff_norm=float(np.linalg.norm([
        complex(original.get(k,0j))-representative.get(k,0j)
        for k in set(original)|set(representative)]))
    master=master_weinberg_polynomial_check() if check_master else None
    return {"format":FORMAT,"status":"Aminus4PolynomialRepresentativeVerified_NotPhysicalMatchingNormalization",
            "model":MODEL,"mu_F_gev":info["mu_F_gev"],
            "scalar_real_dimension":info["scalar_real_dimension"],
            "source_ordered_coordinate_count":len(merged_labels),
            "lepton_dummy_alignment":alignment,
            "symmetric_nonzero_coordinate_count":len(representative),
            "source_l2_norm_gev_inv":norm_raw,
            "symmetric_representative_l2_norm_gev_inv":norm_proj,
            "relative_representation_norm_difference":diff_norm/max(norm_raw,1e-300),
            "projection_removes_only_polynomial_null_components":proof["preserved_to_numerical_precision"],
            "norm_change_is_not_a_physical_amplitude_loss":True,
            "formal_fields":"commuting psi placeholders (Weyl-bilinear symmetric) and real commuting scalar phi",
            "operator_polynomial":"P=sum C[i,p,j,q,a,b]*psi[i,p]*psi[j,q]*phi[a]*phi[b]",
            "projection_formula":"Cpoly=(Cordered+Cfermion_swap+Cscalar_swap+Cboth_swaps)/4",
            "polynomial_identity":proof,"master_example_convention":master,
            "matchete_to_master_overall_matching_factor_certified":False,
            "spinor_Fierz_EOM_convention_certified":False,
            "full_three_generation_gamma_constructed":False,
            "physical_eft1_wilson_running_certified":False,
            "symmetric_polynomial_coordinates":[_encode_entry(k,z) for k,z in representative.items()]}
