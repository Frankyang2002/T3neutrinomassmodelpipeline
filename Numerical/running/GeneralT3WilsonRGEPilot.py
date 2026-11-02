"""A_-4 diagnostic bridge between matched ordered tensors and master Wilson RGE.

READ-ONLY audit. The raw source tensor is deliberately not identified with
canonical separately symmetric C_ijab: normalizations, Fierz/EOM projection,
and the three-family EFT gauge/Yukawa dictionary have not been certified.
"""
from __future__ import annotations

import math
import re
from typing import Mapping

import numpy as np

from Numerical.running.GeneralT3OrderedTensor import scalar_basis
from Numerical.running.GeneralT3WilsonEvolution import (
    aggregate_source_coordinates, boundary_coordinates,
)

LABEL = re.compile(r'^L([12])([12]):S([1-9]\d*),([1-9]\d*):F([123])([123])$')


def _swap_f(k: tuple[int,...]) -> tuple[int,...]:
    i,p,j,q,a,b=k
    return (j,q,i,p,a,b)


def _swap_s(k: tuple[int,...]) -> tuple[int,...]:
    i,p,j,q,a,b=k
    return (i,p,j,q,b,a)


def read_flavour_components(labels: list[str], values: np.ndarray, model: str) -> dict[tuple[int,...],complex]:
    ns = scalar_basis(model).total_dimension
    if len(labels) != len(set(labels)) or values.shape != (len(labels),):
        raise ValueError('Invalid summed ordered boundary basis')
    result={}
    for label,value in zip(labels,values):
        match=LABEL.fullmatch(label)
        if not match:
            raise ValueError('Non-LLSS or malformed summed boundary coordinate '+label)
        i,j,a,b,p,q=map(int,match.groups())
        if a>ns or b>ns:
            raise ValueError('Scalar coordinate exceeds model basis')
        key=(i,p,j,q,a,b)
        if key in result or not np.isfinite(value):
            raise ValueError('Duplicate/nonfinite Wilson component')
        result[key]=complex(value)
    return result


def symmetry_audit(components: Mapping[tuple[int,...],complex]) -> dict:
    """Evaluate fourfold permutation projection WITHOUT modifying boundary.

    The factors 1/4 are for a mathematical projection only, NOT a physical
    Lagrangian normalization or permission to discard the antisymmetric data.
    """
    all_keys=set(components)
    for key in list(all_keys):
        swapped=_swap_f(key)
        all_keys.update((swapped,_swap_s(key),_swap_s(swapped)))
    if not all_keys:
        raise ValueError('Empty boundary tensor')
    raw_sq=projected_sq=fermion_sq=scalar_sq=lost_sq=0.0
    enlarged=0
    for key in all_keys:
        v=components.get(key,0j)
        vs=components.get(_swap_s(key),0j)
        vf=components.get(_swap_f(key),0j)
        vfs=components.get(_swap_s(_swap_f(key)),0j)
        proj=(v+vs+vf+vfs)/4
        raw_sq+=abs(v)**2
        projected_sq+=abs(proj)**2
        fermion_sq+=abs(v-vf)**2
        scalar_sq+=abs(v-vs)**2
        lost_sq+=abs(v-proj)**2
        if key not in components and abs(proj)>1e-30:
            enlarged+=1
    norm=math.sqrt(raw_sq)
    return {'raw_orbit_component_count':len(all_keys),
            'raw_present_coordinate_count':len(components),
            'new_coordinates_in_separately_symmetric_projection':enlarged,
            'raw_l2_norm_gev_inv':norm,
            'projected_l2_norm_gev_inv':math.sqrt(projected_sq),
            'fermion_exchange_residual_relative':math.sqrt(fermion_sq)/max(norm,1e-300),
            'scalar_exchange_residual_relative':math.sqrt(scalar_sq)/max(norm,1e-300),
            'removed_by_separate_symmetry_projection_relative':math.sqrt(lost_sq)/max(norm,1e-300),
            'symmetry_projection_is_diagnostic_only':True,
            'canonical_operator_normalization_certified':False}


def master_antilinearity_check() -> dict:
    """Evaluate the real master-RGE C* term in an exact one-index toy case.

    It checks the installed source, not a retyped approximation of the term.
    y=1, one Weyl flavour, one real scalar: beta_conjugate(C)=3 C*.
    Thus beta(i)=-3i, while i beta(1)=+3i.
    """
    import sympy as sp
    from types import SimpleNamespace
    from RGE.general.WilsonTensorRGE import (
        WilsonRGEInputs, conjugate_coefficient_yukawa_term,
    )
    model=SimpleNamespace(total_real_scalar_dimension=1)
    inputs=WilsonRGEInputs(fermion_dimension=1,
                            yukawa=lambda *_:sp.Integer(1),
                            quartic=lambda *_:sp.S.Zero)
    class FixedCoefficient:
        def __init__(self, z): self.z=sp.sympify(z)
        def __getitem__(self, key): return self.z
    f=lambda z:sp.simplify(conjugate_coefficient_yukawa_term(
        model,inputs,(1,1,1,1),FixedCoefficient(z)))
    b1=f(1)
    bi=f(sp.I)
    expected=sp.simplify(3*sp.conjugate(sp.I))
    if b1 != 3 or sp.simplify(bi-expected)!=0 or sp.simplify(bi-sp.I*b1)==0:
        raise RuntimeError('Actual master-RGE conjugate source changed')
    return {'master_source_checked':'RGE.general.WilsonTensorRGE.conjugate_coefficient_yukawa_term',
            'unit_yukawa_toy_beta_of_C_1':str(b1),
            'unit_yukawa_toy_beta_of_C_i':str(bi),
            'i_times_beta_of_C_1':str(sp.I*b1),
            'complex_linearity_violated':True,
            'necessary_evolution_structure':'16pi2*dC/dlnmu=A@C+B@conjugate(C)',
            'complete_A_B_constructed':False}


def audit_a_minus_4(payload: Mapping, *, include_master_check: bool=True) -> dict:
    model='A_-4'
    source_labels,source_values,info=boundary_coordinates(payload,model)
    labels,values=aggregate_source_coordinates(source_labels,source_values)
    comp=read_flavour_components(labels,values,model)
    ns=scalar_basis(model).total_dimension
    n_left=6 # 3 generations x 2 SU(2) components
    n_lepton_plus_ec=9 # plus 3 charged-lepton singlet Weyl fields
    def c2(n):return n*(n+1)//2
    closure={'L_only_fermion_pair_by_scalar_pair':c2(n_left)*c2(ns),
             'L_plus_eC_fermion_pair_by_scalar_pair':c2(n_lepton_plus_ec)*c2(ns),
             'note':'illustrative separately symmetric L/eC closure sizes only; not a certified full SM operator closure'}
    report={'format':'t3-a-minus4-master-wilson-rge-basis-audit-v1',
            'model':model,'status':'BoundaryAuditedPhysicalGammaNotConstructed',
            'mu_F_gev':info['mu_F_gev'],
            'source_term_flavour_coordinate_count':len(source_labels),
            'summed_ordered_flavour_coordinate_count':len(labels),
            'real_scalar_dimension':ns,
            'source_symmetry_diagnostic':symmetry_audit(comp),
            'hypothetical_symmetric_closure_sizes':closure,
            'gamma_A_source_verified':False,
            'gamma_B_source_verified':False,
            'physical_eft1_wilson_rge_certified':False,
            'existing_complex_linear_gamma_insufficient_in_general':True}
    if include_master_check:
        report['actual_master_rge_antilinearity']=master_antilinearity_check()
    return report
