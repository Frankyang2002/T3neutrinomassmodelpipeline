"""Source-inventory-checked, SU(2)-invariant reference spans for extra HHSS terms.

The registry gives the allowed field monomials and their number of independent
invariants, *not* a certified mapping from its RGBeta symbols to Matchete CGs.
This module constructs a canonical algebraic basis and its fourth derivatives
without assigning those basis vectors to the RGBeta coupling names.  A basis
change and signed complex phases must be separately certified before numerical
use.  All reference vectors use exact SymPy arithmetic.
"""
from __future__ import annotations

from functools import lru_cache
from math import factorial
import sympy as sp

from .portal_normalization import su2_generators

# Counts: H, BarH, S1, BarS1, S2, BarS2
DIMS={'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}
FORMAT='t3-ws2-extra-hhss-exact-su2-reference-space-v1'


def _check_counts(counts):
    counts=tuple(counts)
    if len(counts)!=6 or any(type(v) is not int or v<0 for v in counts):
        raise ValueError('Expected six nonnegative integer field multiplicities')
    if sum(counts[:2])!=2 or sum(counts[2:])!=2:
        raise ValueError('Not an HHSS monomial')
    return counts


def _kronecker_rep(rep,slot,dims):
    mats=[sp.eye(x) for x in dims]
    mats[slot]=rep
    return sp.kronecker_product(*mats)


def _permute_tensor(v, dims, pos1, pos2):
    """Exchange two identical external boson slots in a tensor vector."""
    if dims[pos1]!=dims[pos2]:
        raise ValueError('Swapping unequal field dimensions')
    shape=tuple(dims)
    a=sp.MutableDenseNDimArray(v,shape)
    # use vector to avoid a SymPy scalar when d=1
    out=[]
    for ix in __import__('itertools').product(*(range(d) for d in dims)):
        swap=list(ix)
        swap[pos1],swap[pos2]=swap[pos2],swap[pos1]
        out.append(a[tuple(swap)])
    return sp.Matrix(out)


def _selected_independent(vecs):
    selected=[]
    rank=0
    for v in vecs:
        if not any(v):
            continue
        k=sp.Matrix.hstack(*selected,v) if selected else v
        rank_new=k.rank()
        if rank_new>rank:
            selected.append(v)
            rank=rank_new
    return selected


@lru_cache(maxsize=None)
def reference_invariants(d1:int,d2:int, counts:tuple[int,...]):
    """Exact SU(2) singlet subspace with identical-boson projection.

    T acts on unbar fields; -T* acts on conjugate fields. Columns carry
    coefficient tensors, hence the transposed representation acts on the
    invariant coefficient vector.  No normalization of a named RGBeta
    `GroupInvariant` or Matchete CG is inferred.
    """
    counts=_check_counts(counts)
    dimensions=(2,2,d1,d1,d2,d2)
    types=tuple(k for k,c in enumerate(counts) for _ in range(c))
    dims=tuple(dimensions[k] for k in types)
    if len(types)!=4:
        raise AssertionError('Quartic must have exactly four fields')
    dimall=__import__('math').prod(dims)
    reps=[]
    for k in types:
        d=dimensions[k]
        t=su2_generators(d)
        reps.append(tuple((-x.conjugate() if k in (1,3,5) else x)
                          for x in t))
    stacked=[]
    for axis in range(3):
        total=sp.zeros(dimall)
        for slot in range(4):
            total+=_kronecker_rep(reps[slot][axis].T,slot,dims)
        stacked.append(total)
    kernel=sp.Matrix.vstack(*stacked).nullspace()
    projected=[]
    for v in kernel:
        w=v
        for k in sorted(set(types)):
            positions=[i for i,x in enumerate(types) if x==k]
            if len(positions)==2:
                w=sp.simplify((w+_permute_tensor(w,dims,*positions))/2)
            elif len(positions)>2:
                raise ValueError('HHSS should not contain >2 identical types')
        projected.append(w)
    basis=_selected_independent(projected)
    # Verify actual Ward identities and exchange symmetry for every basis.
    for v in basis:
        for total in stacked:
            if any(sp.simplify(x)!=0 for x in total*v):
                raise ValueError('Candidate reference basis is not SU2 invariant')
        for k in set(types):
            pos=[i for i,x in enumerate(types) if x==k]
            if len(pos)==2 and any(sp.simplify(x)!=0 for x in v-_permute_tensor(v,dims,*pos)):
                raise ValueError('Nonbosonic HHSS reference basis')
    return types,dims,tuple(basis)


def _realfield_poly(types,dims,v,d1,d2):
    """Convert one invariant monomial to canonical interleaved real fields.

    For each complex component z=(R+iI)/sqrt2. Returns the *complex*
    monomial P; physical V4 for a complex source coefficient c requires
    -c P-c* P* and hence independent Re/Im coefficient directions.
    """
    h=sp.symbols('h0:'+str(4))
    s=sp.symbols('s0:'+str(2*(d1+d2)))
    allfields=h+s
    offsets={0:0,1:0,2:4,3:4,4:4+2*d1,5:4+2*d1}
    data=[]
    for typ,dim in zip(types,dims):
        off=offsets[typ]
        sign=-1 if typ in (1,3,5) else 1
        data.append(tuple((allfields[off+2*i]+sign*sp.I*allfields[off+2*i+1])/sp.sqrt(2)
                          for i in range(dim)))
    poly=sp.S.Zero
    for ix in __import__('itertools').product(*(range(x) for x in dims)):
        row=__import__('functools').reduce(lambda a,b:a*b,
                                         (data[j][idx] for j,idx in enumerate(ix)),sp.S.One)
        poly+=v[__import__('numpy').ravel_multi_index(ix,dims)]*row
    return sp.Poly(sp.expand(poly),*allfields)


def real_component_derivatives(d1:int,d2:int,types,dims,v):
    """Sparse HHSS real/imag coupling fourth derivatives, exact factorials."""
    poly=_realfield_poly(types,dims,v,d1,d2)
    real={};imag={}
    for powers,coefficient in poly.terms():
        if sum(powers)!=4 or sum(powers[:4])!=2 or sum(powers[4:])!=2:
            raise ValueError('Extra invariant produced non HHSS quartic structure')
        hs=tuple(i for i,e in enumerate(powers[:4]) for _ in range(e))
        ss=tuple(i for i,e in enumerate(powers[4:]) for _ in range(e))
        factor=sp.prod(factorial(e) for e in powers)
        d4=sp.simplify(factor*coefficient)
        key=hs+ss
        vr=sp.simplify(-2*sp.re(d4))
        vi=sp.simplify(+2*sp.im(d4))
        if vr!=0:real[key]=sp.sstr(vr)
        if vi!=0:imag[key]=sp.sstr(vi)
    return real,imag


def model_extra_geometry(model:str,registry:dict,*,include_components:bool=True)->dict:
    """Audit all extra monomials, old and new symbols without assigning phases.

    A field-count sector with k independent RGBeta entries must have exactly
    k independent SU2-invariant bosonic tensor directions.  If this is false,
    block instead of dropping or merging symbols.
    """
    d1,d2=DIMS[model[0]]
    sectors={}
    # A missing hypercharge-zero condition invalidates the registry entry.
    alpha=int(model.split('_')[1])
    ys=(sp.Rational(1,2),-sp.Rational(1,2),
        sp.Rational(alpha,2),-sp.Rational(alpha,2),
        sp.Rational(alpha+2,2),-sp.Rational(alpha+2,2))
    for source,items in [('new',registry.get('new_rgbeta_symbols',[])),
                         ('preexisting',registry.get('preexisting_rgbeta_symbols',[]))]:
        for entry in items:
            counts=tuple(entry['counts'])
            if sum(counts[:2])!=2 or sum(counts[2:])!=2:
                continue
            if sum((ys[i]*counts[i] for i in range(6)),sp.S.Zero)!=0:
                raise ValueError(f'{model}/{counts}: nonzero hypercharge extra HHSS structure')
            name=entry['name'] if source=='new' else entry['existing_name']
            sectors.setdefault(counts,[]).append({'rgbeta_symbol':name,'source':source})
    rows=[];candidates=0
    for counts,symbols in sorted(sectors.items()):
        types,dims,basis=reference_invariants(d1,d2,counts)
        if len(symbols)!=len(basis):
            raise ValueError(f'{model}/{counts}: registry lists {len(symbols)} symbols but exact SU2 boson invariant rank is {len(basis)}')
        candidates+=len(symbols)
        reference=[]
        if include_components:
            for index,v in enumerate(basis):
                re,im=real_component_derivatives(d1,d2,types,dims,v)
                if not re and not im:
                    raise ValueError('Zero real-component reference invariant')
                reference.append({
                    'id':f'reference_invariant_{index+1}',
                    'Re_c_unit_tensor_components':[
                        {'H1_H2_S1_S2_0based':list(k),'value_exact':x} for k,x in sorted(re.items())],
                    'Im_c_unit_tensor_components':[
                        {'H1_H2_S1_S2_0based':list(k),'value_exact':x} for k,x in sorted(im.items())],
                })
        rows.append({
            'counts_H_BarH_S1_BarS1_S2_BarS2':list(counts),
            'source_registry_symbols':symbols,
            'independent_complex_invariant_count':len(basis),
            'exact_reference_real_component_tensor_basis':reference,
            'physical_mapping_to_named_RGBeta_coefficients':None,
            'matchete_CG_phase_and_orientation':None,
            'status':'exact_SU2_reference_space_only_not_source_normalized'})
    return {
        'format':FORMAT,'model':model,
        'registry_HHSS_extra_symbols':candidates,
        'monomial_sectors':rows,
        'reference_coordinate_convention':'canonical SU2 highest-weight indices, interleaved real fields (R,I), V4=-cP-c*conjugate(P)',
        'named_source_coupling_map_status':'blocked_no_matchete_RGBeta_CG_basis_change',
        'numerical_extra_couplings_available':False,
        'physical_C5_certified':False,
    }
