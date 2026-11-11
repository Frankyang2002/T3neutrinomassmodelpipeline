"""Exact Matchete-component -> real-Cartesian quartic candidate conversion.

This is strictly the Matchete-side potential.  It does NOT equate Matchete
couplings to RGBeta running parameters, and MUST NOT be loaded as a certified
`t3-running-real-quartic-response-v1` input.

Source convention (subject to independent source-level verification):
  Phi_i = (x_(2i) + I*x_(2i+1))/sqrt(2)
  L_extra = p*c*CG[species...] + h.c., V_extra = -L_extra
  V4 = lambda_abcd*x_a*x_b*x_c*x_d/4!

Every monomial has exactly four real scalar fields.  Values remain exact SymPy
expressions until final serialization and are checked to be real.
"""
from __future__ import annotations

from collections import defaultdict
from itertools import product
from math import factorial
from typing import Any, Mapping

import sympy as sp

from scripts.RGBetaTensorAlgebra import exact, normalized_record, same_species_symmetry

SPECIES = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
MODEL_DIMS = {'A': (1, 3), 'B': (2, 2), 'C': (2, 2),
              'D': (3, 1), 'E': (3, 3)}

class CartesianConversionError(ValueError):
    pass


def _clean(x: sp.Expr) -> sp.Expr:
    y = sp.simplify(sp.expand(x))
    if y.free_symbols or y.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
        raise CartesianConversionError(f'Unresolved or nonfinite coefficient: {y}')
    return y


def _real(x: sp.Expr) -> sp.Expr:
    y = _clean(sp.expand_complex(x))
    if sp.simplify(sp.im(y)) != 0:
        raise CartesianConversionError(f'Coefficient not real: {x}')
    return sp.simplify(sp.re(y))


def _species_slot(model: str, species: int, component: int) -> tuple[int, int, bool]:
    if model[0] not in MODEL_DIMS:
        raise CartesianConversionError(f'Unknown representation: {model}')
    d1, d2 = MODEL_DIMS[model[0]]
    group, conj = species // 2, species % 2 == 1
    dim = (2, d1, d2)[group]
    if component < 0 or component >= dim:
        raise CartesianConversionError('Component outside SU(2) representation')
    offset = (0, 4, 4 + 2*d1)[group] + 2*component
    return offset, offset + 1, conj


def _expand_field_product(model: str, species: list[int], components: tuple[int,...]):
    """Return P(x) as a sparse, 4th-degree, sorted-monomial coefficient map."""
    partial = {(): sp.Integer(1)}
    invroot = 1/sp.sqrt(2)
    for s, component in zip(species, components):
        r, im, conj = _species_slot(model, s, component)
        terms = ((r, invroot), (im, -sp.I*invroot if conj else sp.I*invroot))
        following = defaultdict(lambda: sp.Integer(0))
        for indices, coefficient in partial.items():
            for idx, fac in terms:
                key = tuple(sorted((*indices, idx)))
                following[key] += coefficient * fac
        partial = following
    return {k: _clean(v) for k,v in partial.items() if _clean(v) != 0}


def _offset(indices: tuple[int,...], dims: tuple[int,...]) -> int:
    out = 0
    for i,d in zip(indices,dims):
        out = out*d+i
    return out


def _canonical(value: sp.Expr) -> str:
    return str(sp.simplify(value))


def convert_record(model: str, record: Mapping[str, Any]) -> dict[str,Any]:
    """Convert one native Matchete quartic record to sparse Cartesian responses.

    The resulting two directions are for Re(c), Im(c), where c is the MAtchete
    coefficient in L = +p*c*O + h.c.  They are NOT RGBeta directions.
    """
    if model not in _supported_models():
        raise CartesianConversionError(f'Unsupported model: {model}')
    same_species_symmetry(record)
    normalized = normalized_record(record, model)
    species = [i for i,n in enumerate(normalized['counts']) for _ in range(n)]
    if record['species'] != [SPECIES[i] for i in species]:
        raise CartesianConversionError('Native species ordering mismatch')
    dims = tuple(normalized['dimensions'])
    expected_dims = tuple((2,*MODEL_DIMS[model[0]])[s//2] for s in species)
    if dims != expected_dims:
        raise CartesianConversionError(f'Native SU2 dimensions mismatch {dims} != {expected_dims}')
    prefactor = exact(record['potential_prefactor'])
    if sp.simplify(prefactor) == 0:
        raise CartesianConversionError('Zero Lagrangian prefactor')
    polynomial = defaultdict(lambda: sp.Integer(0))
    for components in product(*(range(d) for d in dims)):
        cg = normalized['components'][_offset(components,dims)]
        if cg == 0:
            continue
        for monomial, val in _expand_field_product(model,species,components).items():
            polynomial[monomial] += prefactor*cg*val
    responses = []
    n = 4+2*sum(MODEL_DIMS[model[0]])
    for component in ('real','imag'):
        sparse = []
        for monomial, coeff in sorted(polynomial.items()):
            # V = -(c P + conjugate(c P)); when c=a+i b,
            # dV/da = -2*Re(P), dV/db = +2*Im(P).
            v = (-2*sp.re(coeff) if component=='real' else 2*sp.im(coeff))
            v = _real(v)
            if v == 0:
                continue
            multiplicity = sp.Integer(1)
            for i in set(monomial):
                multiplicity *= factorial(monomial.count(i))
            tensor_component = _real(v*multiplicity)
            sparse.append({'indices':list(monomial),
                           'lambda_abcd':_canonical(tensor_component),
                           'monomial_coefficient':_canonical(v)})
        responses.append({'component':component,'entries':sparse})
    return {'symbol':normalized['symbol'],
            'native_species':record['species'],
            'cartesian_dimensions':n,
            'tensor_normalization':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
            'potential_source':'V=-[prefactor*c*O+h.c.]',
            'responses':responses}


def _supported_models() -> frozenset[str]:
    return frozenset(('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
                      'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
                      'E_+0','E_+2','E_-2','E_-4'))


def convert_matchete_export(payload: Mapping, registry: Mapping) -> dict:
    if (payload.get('format') != 'stage22-exact-tensors-v1' or
        payload.get('basis') != 'Matchete' or
        payload.get('source') != 'MatcheteInvariantTensors'):
        raise CartesianConversionError('Not a native Matchete Stage 22 tensor export')
    expected_models = _supported_models()
    if set(registry['models']) != set(expected_models) or set(payload['models']) != set(expected_models):
        raise CartesianConversionError('Incomplete 16-model source inventory')
    out = {}
    total = 0
    for model in sorted(expected_models):
        records = payload['models'][model].get('tensors')
        expected = registry['models'][model]['total_additional_quartics']
        if not isinstance(records,list) or len(records) != expected:
            raise CartesianConversionError(f'{model}: native quartic count mismatch')
        produced = [convert_record(model, record) for record in records]
        if len({r['symbol'] for r in produced}) != len(produced):
            raise CartesianConversionError(f'{model}: duplicate Matchete coupling symbols')
        out[model] = {'records':produced,'quartic_slots':expected}
        total += expected
    if total != 50:
        raise CartesianConversionError(f'Expected 50 slots, got {total}')
    return {
        'format':'t3-matchete-candidate-cartesian-quartics-v1',
        'certification_status':'UNCERTIFIED: Matchete-only',
        'rgbeta_mapping_certified':False,
        'cartesian_field_assumption':'complex component=(real+i*imag)/sqrt(2)',
        'potential_convention':'V=-L with L=c*O+h.c.',
        'tensor_convention':'V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
        'models':out,
    }
