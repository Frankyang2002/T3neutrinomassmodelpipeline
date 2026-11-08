"""Independent A(-4) F-only LLSS reconstruction and *blocked* Wilson contract.

READS Matchete raw seed, full tree EFT, source Yukawa/mass boundary, and the
ordered/polynomial exports. Does not identify Matchete GammaCC with the master
Weyl epsilon. Never produces a physically certified Wilson boundary.

Run from repository root:
 python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_boundary_verify

The extended 45-Weyl-fermion ordering is an explicitly PROPOSED three-family
extension of the existing one-family SM ordering, not a currently verified
RGE production-basis constructor.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
import re
import subprocess
from typing import Mapping

import numpy as np

MODEL = 'A_-4'
FORMAT = 't3-a-minus4-fermion-first-conditional-wilson-contract-v1'
RAW = Path('output/general_observable/fermion_first_all16_v1/A_m4/data')
SOURCE = Path('output/general_observable/general_wilson_boundary_v1/A_m4/wilson_boundary.json')
ORDERED = Path('output/general_observable/ordered_wilson_tensor_v4/A_m4/ordered_tensor.json')
POLY = Path('output/general_observable/general_wilson_polynomial_batch_v2/A_m4/polynomial_bridge.json')
RUNNING = Path('output/general_observable/all16_v1/A_m4/running.json')
OUTPUT = Path('output/general_observable/a_minus4_boundary_batch1/conditional_boundary.json')


class BoundaryBlocked(ValueError):
    pass


def demand(ok: bool, reason: str) -> None:
    if not ok:
        raise BoundaryBlocked(reason)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    demand(path.is_file(), f'Missing actual source file: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    demand(isinstance(obj, dict), f'Malformed JSON object: {path}')
    return obj


def compact(txt: str) -> str:
    return re.sub(r'\s+', '', txt)


def signed_real_expr(source: str) -> complex:
    """Parse only arithmetic from exact Matchete native sparse CG values."""
    demand(bool(re.fullmatch(r'[\d\sI\+\-\*/\^().]+', source)),
           'Unexpected native CG arithmetic syntax')
    node = ast.parse(source.replace('^', '**'), mode='eval')
    def walk(x):
        if isinstance(x, ast.Expression):
            return walk(x.body)
        if isinstance(x, ast.Constant) and type(x.value) in (int, float):
            return complex(x.value)
        if isinstance(x, ast.Name) and x.id == 'I':
            return 1j
        if isinstance(x, ast.UnaryOp) and isinstance(x.op, (ast.UAdd, ast.USub)):
            return walk(x.operand) * (-1 if isinstance(x.op, ast.USub) else 1)
        if isinstance(x, ast.BinOp):
            a, b = walk(x.left), walk(x.right)
            if isinstance(x.op, ast.Add): return a+b
            if isinstance(x.op, ast.Sub): return a-b
            if isinstance(x.op, ast.Mult): return a*b
            if isinstance(x.op, ast.Div): return a/b
            if isinstance(x.op, ast.Pow): return a**b
        raise BoundaryBlocked(f'Unsupported native CG arithmetic: {source}')
    value = walk(node)
    demand(math.isfinite(value.real) and math.isfinite(value.imag),
           'Nonfinite native CG value')
    return value


def parse_native_cg(registry: list[dict]) -> tuple[np.ndarray, dict]:
    demand(len(registry) == 1 and registry[0].get('Name') == 'T3MixCG',
           'A(-4) CG registry changed or incomplete')
    item = registry[0]
    demand(compact(item.get('RepsInputForm', '')) ==
           compact('{Bar[SU2L[fund]],Bar[SU2L[fund]],T3BSMd3}'),
           'Native CG representation ordering changed')
    raw = item.get('TensorInputForm', '')
    # This is a deliberately narrow parser for the observed Matchete SparseArray
    # storage format (outer row pointers, inner (j,c) positions, values).
    match = re.fullmatch(
        r'\s*SparseArray\[Automatic,\s*\{2,\s*2,\s*3\},\s*0,\s*'
        r'\{1,\s*\{\{0,\s*3,\s*6\},\s*\{(?P<pairs>(?:\{\s*[12]\s*,\s*[123]\s*\}\s*,?\s*)+)\}\},\s*'
        r'\{(?P<values>[^{}]+)\}\}\]\s*', raw)
    demand(match is not None, 'Native CG SparseArray encoding changed')
    pairs = [(int(a), int(b)) for a, b in
             re.findall(r'\{\s*([12])\s*,\s*([123])\s*\}', match.group('pairs'))]
    vals = [signed_real_expr(s.strip()) for s in match.group('values').split(',')]
    demand(len(pairs) == len(vals) == 6 and len(set(pairs[:3])) == 3 and
           len(set(pairs[3:])) == 3, 'Malformed native six-entry CG')
    out = np.zeros((2, 2, 3), dtype=complex)
    for n, ((b, c), value) in enumerate(zip(pairs, vals)):
        a = n // 3
        demand(out[a, b-1, c-1] == 0, 'Duplicate native CG coordinate')
        out[a, b-1, c-1] = value
    # Independent SU(2) representation diagnostics: Sym^2(2) contains 3,
    # and these three symmetric pair-tensors must be orthogonal and nonzero.
    residual = float(np.linalg.norm(out - np.swapaxes(out, 0, 1)))
    gram = np.einsum('ijc,ijd->cd', out.conj(), out)
    expected = 2/math.sqrt(3)*np.eye(3)
    gram_error = float(np.linalg.norm(gram - expected))
    demand(residual < 1e-13 and gram_error < 1e-13 and
           np.linalg.matrix_rank(out.reshape(4,3), tol=1e-12) == 3,
           'Native CG fails 2x2 symmetric-triplet orthogonality')
    return out, {'cg_symmetric_pair_residual': residual,
                 'cg_gram_residual': gram_error,
                 'native_triplet_rank': 3,
                 'native_triplet_axes_equated_to_RGE_generators': False}


def inspect_native_seed(seed: dict, tree_text: str) -> dict:
    demand(seed.get('IntegratedFields') == ['F'] and
           seed.get('ActiveHeavyFields') == ['S1', 'S2'] and
           seed.get('Level') == 1, 'Not an F-only source')
    terms = seed.get('TreeWilsonTerms')
    demand(isinstance(terms, list) and len(terms) == 2 and
           seed.get('TreeWilsonTermCount') == 2,
           'Unexpected tree-level LLSS operator inventory')
    demand(tree_text.count('GammaCC') == 2,
           'Full matching-generated tree EFT differs from two-term LLSS inventory')
    by_chirality = {t.get('Chirality'): t.get('TermInputForm', '') for t in terms}
    demand(set(by_chirality) == {'PL', 'PR'}, 'Missing conjugate chiral sector')
    pl = compact(by_chirality['PL'])
    pr = compact(by_chirality['PR'])
    for token in ('Bar[Coupling[y1,', 'Bar[Coupling[y2,',
                  'Bar[Field[NewScalar1,Scalar,{},{}]]',
                  'Field[NewScalar2,Scalar,{Index[d$$1,T3BSMd3]},{}]',
                  'CG[T3MixCG,', 'DiracProduct[GammaCC,Proj[-1]]',
                  'NCM[Transp[Field[l,Fermion,{Index[d$$2,SU2L[fund]],Index[d$$2,Flavor]},{}]]',
                  'Field[l,Fermion,{Index[d$$1,SU2L[fund]],Index[d$$1,Flavor]},{}]',
                  'Coupling[MF,{Index[d$$1,NFlavor]},0]'):
        demand(token in pl, f'PL source conjugation/ordering changed: {token}')
    for token in ('Coupling[y1,', 'Coupling[y2,',
                  'Field[NewScalar1,Scalar,{},{}]',
                  'Bar[Field[NewScalar2,Scalar,{Index[d$$1,T3BSMd3]},{}]]',
                  'CG[Bar[T3MixCG],', 'DiracProduct[GammaCC,Proj[1]]',
                  'NCM[Bar[Field[l,Fermion,{Index[d$$1,SU2L[fund]],Index[d$$1,Flavor]},{}]]',
                  'Transp[Bar[Field[l,Fermion,{Index[d$$2,SU2L[fund]],Index[d$$2,Flavor]},{}]]]',
                  'Coupling[MF,{Index[d$$1,NFlavor]},0]'):
        demand(token in pr, f'PR source conjugation/ordering changed: {token}')
    demand('Bar[Coupling[y1,' not in pr and
           'Bar[Coupling[y2,' not in pr and
           'Bar[Field[NewScalar1,Scalar,{},{}]]' not in pr,
           'PR source is not conjugate to PL')
    for t in terms:
        demand(t.get('Scalar1Count') == 1 and t.get('Scalar2Count') == 1 and
               t.get('CGNames') == ['T3MixCG'], 'Unexpected matching field inventory')
        # Full tree is authoritative, not seed counts alone. Matchete sometimes
        # prints parentheses differently when embedding an expanded term.
        fragment = compact(t['TermInputForm']).strip('()')
        demand(fragment in compact(tree_text),
               'Wilson seed term not found verbatim in matching-generated full tree EFT')
    return {'native_full_tree_GammaCC_term_count': 2,
            'PL_term_index': next(t['Index'] for t in terms if t['Chirality']=='PL'),
            'PR_term_index': next(t['Index'] for t in terms if t['Chirality']=='PR'),
            'holomorphic_source': 'PL: conj(y1[p,r])*conj(y2[q,r])/M[r]',
            'conjugate_source': 'PR: y1[p,r]*y2[q,r]/M[r]',
            'source_lepton_order': 'PL left l[q] transposed; right l[p]',
            'source_scalar_order': 'PL S1_bar * S2; PR S2_bar * S1',
            'full_tree_inventory_checked': True}


def dense_flavour(y1: np.ndarray, y2: np.ndarray, masses: np.ndarray) -> np.ndarray:
    a, b = np.asarray(y1, complex), np.asarray(y2, complex)
    m = np.asarray(masses, float)
    demand(a.shape == b.shape == (3,3) and m.shape == (3,) and
           bool(np.all(np.isfinite(a))) and bool(np.all(np.isfinite(b))) and
           bool(np.all(np.isfinite(m))) and bool(np.all(m>0)),
           'Expected 3x3 complex Yukawas and three positive heavy eigenmasses')
    return (a.conj() / m[None,:]) @ b.conj().T


def decode_matrix(item: dict) -> np.ndarray:
    demand(isinstance(item, dict) and set(('real','imag')) <= set(item),
           'Malformed complex 3x3 source matrix')
    mat = np.asarray(item['real'], dtype=float)+1j*np.asarray(item['imag'], dtype=float)
    demand(mat.shape == (3,3) and bool(np.all(np.isfinite(mat))),
           'Expected finite 3x3 matrix')
    return mat


def raw_components(w: np.ndarray, cg: np.ndarray) -> dict[tuple[int, ...], complex]:
    """Independent contraction of PL native SU2/complex scalar source.

    Coordinates: (SU2-of-first-lepton, flavour-first, SU2-second,
                  flavour-second, REAL S1 component, REAL S2 component).
    PL lepton ordering is q,p, but W matrix axes are p,q.
    """
    demand(w.shape == (3,3) and cg.shape == (2,2,3), 'Wrong source matrix dimensions')
    out: dict[tuple[int, ...], complex] = defaultdict(complex)
    for a in range(2):    # Matchete lepton p (right in bilinear)
        for b in range(2): # Matchete lepton q (transposed on left)
            for c in range(3):
                if abs(cg[a,b,c]) < 1e-14:
                    continue
                # S1* = (phi5 - i phi6)/sqrt(2)
                # S2_c = (phi[7+2c] + i phi[8+2c])/sqrt(2)
                for sc1, fac1 in ((5, 1/math.sqrt(2)), (6, -1j/math.sqrt(2))):
                    for sc2, fac2 in ((7+2*c, 1/math.sqrt(2)),
                                      (8+2*c, 1j/math.sqrt(2))):
                        factor = cg[a,b,c]*fac1*fac2
                        for p in range(3):
                            for q in range(3):
                                z = w[p,q]*factor
                                if z != 0:
                                    out[b+1,q+1,a+1,p+1,sc1,sc2] += z
    return dict(out)



def check_pr_formal_hc(w: np.ndarray, cg: np.ndarray) -> dict:
    """PR complex scalar/CG/flavour coefficients vs PL conjugates only.

    This does not test a four-component GammaCC bilinear HC identity.
    """
    maxerr=0.0
    count=0
    for a in range(2):
        for b in range(2):
            for c in range(3):
                if cg[a,b,c] == 0: continue
                for s1,pl1,pr1 in ((5,1/math.sqrt(2),1/math.sqrt(2)),
                                    (6,-1j/math.sqrt(2),1j/math.sqrt(2))):
                    for s2,pl2,pr2 in ((7+2*c,1/math.sqrt(2),1/math.sqrt(2)),
                                        (8+2*c,1j/math.sqrt(2),-1j/math.sqrt(2))):
                        for p in range(3):
                            for q in range(3):
                                pl=w[p,q]*cg[a,b,c]*pl1*pl2
                                # PR source has CG[Bar[T3MixCG]], S2bar*S1,
                                # and unbarred Yukawa products. Scalars commute.
                                pr=w[p,q].conjugate()*cg[a,b,c].conjugate()*pr2*pr1
                                maxerr=max(maxerr,abs(pr-pl.conjugate()))
                                count+=1
    demand(maxerr<=1e-13*max(float(np.max(np.abs(w))),1e-100)+1e-25,
           'PR formal scalar/flavour conjugate differs from PL')
    return {'PL_PR_scalar_CG_flavour_conjugate_components_checked':count,
            'PR_formal_hc_max_abs_residual_gev_inv':maxerr,
            'PR_native_spinor_HC_phase_certified':False}


def global_lepton(su2: int, flavour: int) -> int:
    return 2*(flavour-1) + su2


def orbit(i: int, j: int, a: int, b: int) -> tuple[tuple[int,int,int,int],...]:
    return tuple(sorted({(i,j,a,b),(j,i,a,b),(i,j,b,a),(j,i,b,a)}))


def symmetrise(native: Mapping[tuple[int, ...], complex]) -> tuple[dict, dict]:
    total = defaultdict(complex)
    for (i,p,j,q,a,b), v in native.items():
        fi, fj = global_lepton(i,p), global_lepton(j,q)
        key = (min(fi,fj),max(fi,fj),min(a,b),max(a,b))
        total[key] += v
    independent = {k:v/len(orbit(*k)) for k,v in total.items() if v != 0}
    return independent, dict(total)


def reconstruct(independent: Mapping[tuple[int,int,int,int],complex]) -> dict:
    out = defaultdict(complex)
    for k, v in independent.items():
        for i,j,a,b in orbit(*k):
            out[min(i,j),max(i,j),min(a,b),max(a,b)] += v
    return dict(out)


def max_residual(a: Mapping, b: Mapping) -> tuple[float,float]:
    keys=set(a)|set(b)
    maxerr=max((abs(complex(a.get(k,0))-complex(b.get(k,0))) for k in keys),default=0.)
    scale=max((abs(complex(v)) for v in a.values()),default=0.)
    return maxerr, maxerr/max(scale,1e-300)


def check_existing_ordered(ordered: dict, w: np.ndarray, source: dict) -> dict:
    demand(ordered.get('model') == MODEL and len(ordered.get('operators',[])) == 1,
           'Wrong ordered tensor A(-4) payload')
    op = ordered['operators'][0]
    demand(op['source_index'] == source['operator_terms'][0]['source_index'] and
           op.get('flavour_matrix_axes_external_indices') == ['d$$1','d$$2'],
           'Matchete ordered source index/flavour axes changed')
    x = {}
    for row in op.get('components',[]):
        ii,jj = row['lepton_su2']; a,b=row['scalar_real']
        key= (ii,jj,a,b)
        demand(key not in x,'Duplicate source ordered component')
        x[key]=decode_matrix(row['ordered_flavour_coefficient_gev_inv'])
    demand(len(x)==24, 'Expected 24 real components from the six native CG entries')
    cg, _ = parse_native_cg(source['source_CG_registry'])
    native=raw_components(w,cg)
    projected = defaultdict(lambda:np.zeros((3,3),complex))
    for (i,q,j,p,a,b),z in native.items():
        projected[(i,j,a,b)][p-1,q-1] += z
    keys=set(x)|set(projected)
    residual = max(float(np.max(np.abs(x.get(k,np.zeros((3,3)))-
                                        projected.get(k,np.zeros((3,3)))))) for k in keys)
    scale=max(float(np.max(np.abs(v))) for v in projected.values())
    demand(residual <= 1e-10*scale+1e-25,
           f'Independent native CG projection != existing ordered tensor: {residual}')
    return {'native_ordered_component_count':len(x),
            'native_ordered_matrix_max_abs_residual_gev_inv':residual}


def check_existing_poly(poly: dict, independent: dict) -> dict:
    demand(poly.get('model') == MODEL, 'Polynomial reference model mismatch')
    x = {}
    for row in poly['independent_symmetric_coordinates']:
        a,b,c,d,e,f=row['coordinate']
        fi,fj=global_lepton(a,b),global_lepton(c,d)
        key=(min(fi,fj),max(fi,fj),min(e,f),max(e,f))
        demand(key not in x, 'Duplicate canonical polynomial coordinate')
        x[key]=complex(row['coefficient_real'],row['coefficient_imag'])
        demand(row['orbit_size']==len(orbit(*key)), 'Wrong saved polynomial orbit multiplicity')
    maxerr, rel = max_residual(independent,x)
    demand(maxerr <= max(1e-25,1e-10*max((abs(v) for v in independent.values()),default=0)),
           f'Independent native polynomial != existing polynomial tensor: {maxerr}')
    return {'native_independent_coordinate_count':len(independent),
            'polynomial_max_abs_residual_gev_inv':maxerr,
            'polynomial_max_relative_residual':rel}


def proposed_fermions() -> list[dict]:
    """Proposed three-family 45-Weyl extension of build_sm_eft() block order."""
    blocks = [('L',2,-.5),('eC',1,1.)]
    for colour in range(1,4):
        blocks.extend([(f'Q{colour}',2,1/6),
                       (f'uC{colour}',1,-2/3),(f'dC{colour}',1,1/3)])
    out=[]
    for field, count, hyper in blocks:
        for f in range(1,4):
            for s in range(1,count+1):
                out.append({'index':len(out)+1,'field':field,
                            'generation':f,'component':s,'hypercharge':hyper})
    demand(len(out)==45 and [x['index'] for x in out[:6]]==list(range(1,7)),
           'Three-family extended SM EFT basis bad')
    return out


def scalars() -> list[dict]:
    result=[]
    for f, components, y in [('H',2,.5),('S1',1,-2.),('S2',3,-1.)]:
        for c in range(1,components+1):
            for parity in ('Re','Im'):
                result.append({'index':len(result)+1, 'field':f,
                               'component':c, 'part':parity, 'hypercharge':y})
    return result


def model_basis_report() -> dict:
    return {'fermions':proposed_fermions(), 'real_scalars':scalars(),
            'fermion_dimension':45,'scalar_dimension':12,
            'full_tensor_shape':[45,45,12,12],
            'fermion_order_status':'PROPOSED three-generation extension of source build_sm_eft one-generation species order; not production RGE basis verified',
            'scalar_order_status':'VERIFIED from GeneralT3OrderedTensor.scalar_basis and ScalarBasis.t3',
            'zero_boundary_outside':'i,j in L components 1..6, (a in 5..6 and b in 7..12) or scalar-swapped; all other entries zero',
            'nonzero_entries_format':'independent i<=j,a<=b; reconstruct each distinct fermion/scalar orbit with equal coefficient'}


def verify_from_sources(root: Path, *, check_running: bool=True) -> dict:
    paths = {'native_seed':root/RAW/'eft1_after_F_wilson_seed.json',
             'native_full_tree':root/RAW/'eft1_after_F_tree_inputform.txt',
             'source_flavour':root/SOURCE, 'ordered_tensor':root/ORDERED,
             'polynomial':root/POLY}
    if check_running:
        paths['running']=root/RUNNING
    hashes={key:{'relative_path':str(p.relative_to(root)), 'sha256':sha(p)}
            for key,p in paths.items() if p.is_file()}
    demand(len(hashes)==len(paths), 'Missing required original Matchete or numerical source')
    seed=read(paths['native_seed'])
    source=read(paths['source_flavour'])
    ordered=read(paths['ordered_tensor'])
    poly=read(paths['polynomial'])
    inventory=inspect_native_seed(seed,paths['native_full_tree'].read_text(encoding='utf-8-sig'))
    cg,cg_proof=parse_native_cg(seed['CGRegistry'])
    demand(source.get('model')==MODEL and len(source.get('operator_terms',[]))==1,
           'Expected A(-4) one holomorphic PL source')
    row=source['operator_terms'][0]
    demand(row.get('source_index')==inventory['PL_term_index'] and
           row.get('flavour_matrix_axes_external_indices')==['d$$1','d$$2'] and
           [(f['name'],f['conjugated'],f['external_index'])
            for f in row.get('factors',[])]==[('y1',True,'d$$1'),('y2',True,'d$$2')],
           'Wrong original Yukawa source/conjugations/flavour axes')
    demand(compact(row['source_term_inputform']) ==
           compact(next(t['TermInputForm'] for t in seed['TreeWilsonTerms']
                        if t['Chirality']=='PL')),
           'Source flavour is not raw matching-generated PL operator')
    demand(source.get('source_CG_registry') == seed['CGRegistry'],
           'Source CG differs from native Matchete registry')
    original_w=decode_matrix(row['ordered_flavour_matrix_gev_inv'])
    running_checks={'UV_evolved_3x3_flavour_recomputed':False}
    if check_running:
        from scripts.RunT3GeneralWilsonBoundary import _running_basis
        running,masses,yukawas,basis = _running_basis(MODEL,paths['running'],{'y1','y2'})
        demand(abs(float(source['mu_F_gev'])-float(running['scales_gev']['F']))<=1e-9,
               'F threshold scale mismatch')
        demand(np.allclose(masses,source['heavy_masses_at_F_gev'],rtol=1e-11,atol=1e-7),
               'Three heavy eigenmasses changed')
        expected_w=dense_flavour(yukawas['y1'],yukawas['y2'],masses)
        absdiff=float(np.max(np.abs(expected_w-original_w)))
        demand(absdiff <= 1e-11*max(float(np.max(np.abs(expected_w))),1e-100),
               f'Source Wilson flavour matrix not equal to recomputed UV state: {absdiff}')
        running_checks={'UV_evolved_3x3_flavour_recomputed':True,
                        'flavour_matrix_max_abs_residual_gev_inv':absdiff,
                        'heavy_flavour_mass_basis_residual':basis['heavy_basis_residual'],
                        'three_mass_eigenstates_gev':[float(x) for x in masses],
                        'mu_F_gev':float(source['mu_F_gev'])}
    native=raw_components(original_w,cg)
    conjugate_check=check_pr_formal_hc(original_w,cg)
    independent,total=symmetrise(native)
    reconstructed=reconstruct(independent)
    err, rel=max_residual(total,reconstructed)
    demand(err<=1e-10*max((abs(v) for v in total.values()),default=0)+1e-25,
           'Symmetric polynomial reconstruction changed Matchete source')
    ordered_checks=check_existing_ordered(ordered,original_w,source)
    polynomial_checks=check_existing_poly(poly,independent)
    git_rev=None
    try:
        p=subprocess.run(['git','rev-parse','HEAD'],cwd=root,check=True,
                         capture_output=True,text=True,timeout=5)
        git_rev=p.stdout.strip()
    except (OSError,subprocess.SubprocessError):
        pass
    return {'format':FORMAT,'model':MODEL,
            'required_eft_order':'UV -> integrate out ALL F -> EFT1(SM+S1+S2) -> integrate out S1 and S2 simultaneously -> SMEFT',
            'field_representations':{'F':{'SU2_dimension':2,'Y':-1.5,'Dirac':True,'mass_eigenstates':3},
                                     'S1':{'SU2_dimension':1,'Y':-2.0,'complex_scalar':True},
                                     'S2':{'SU2_dimension':3,'Y':-1.0,'complex_scalar':True}},
            'status':'CONDITIONALLY_RECONSTRUCTED_PHYSICAL_SPINOR_NORMALIZATION_BLOCKED',
            'physical_certification':False, 'git_revision':git_rev,
            'sources':hashes,'native_inventory_checks':inventory,
            'native_CG_checks':cg_proof,'native_PR_formal_conjugate_check':conjugate_check,
            'numeric_flavour_checks':running_checks,
            'independent_component_checks':ordered_checks,
            'independent_polynomial_checks':polynomial_checks,
            'polynomial_reconstruction_max_abs_residual':err,
            'polynomial_reconstruction_relative_residual':rel,
            'basis_contract':model_basis_report(),
            'convention':{
                'lagrangian':'L5 = sum_ijab C_ijab psi_i^alpha psi_j_alpha phi_a phi_b + h.c.; no 1/2 or 1/4 in displayed sum',
                'Weyl_epsilon_12':'+1; psi_i psi_j = psi_j psi_i after Grassmann contraction',
                'fermion_symmetry':'C_ijab=C_jiab (NOT Hermitian conjugation)',
                'scalar_symmetry':'C_ijab=C_ijba',
                'coefficient_reconstruction':'sum over all i,j,a,b; each orbit member assigned same independent coefficient; orbit 1,2,4',
                'native_selected_sector':'PL from NCM[Transp[L_q],GammaCC PL,L_p], PR kept only for conjugate consistency',
                'native_scalar_complex_split':'S1*=(phi5-i phi6)/sqrt(2); S2[c]=(phi[7+2(c-1)]+i phi[8+2(c-1)])/sqrt(2)',
                'native_to_Weyl_factor':None,
                'native_to_Weyl_factor_symbol':'rho_GammaCC (complex sign/phase/prefactor unmeasured)',
                'physical_tensor_relation':'C_physical = rho_GammaCC * C_formal conditional on matching the RGE field/generator axes',
                'unphysical_factor_one_assumed':False},
            'independent_conditional_tensor_entries':[
                {'i':i,'j':j,'a':a,'b':b,'coefficient_real_if_rho_1':float(z.real),
                 'coefficient_imag_if_rho_1':float(z.imag),
                 'distinct_orbit_size':len(orbit(i,j,a,b))}
                for (i,j,a,b),z in sorted(independent.items())],
            'physical_C_ijab_at_MF':None,
            'checks_passed':[
                'native F-only tree exactly two GammaCC sectors',
                'native PL/PR Yukawa and scalar conjugation structure',
                'native CG support/complex phases/orthogonality',
                'PR is the scalar/CG/flavour complex conjugate of PL (spinor phase not established)',
                'three-flavour numeric factorization from evolved UV state' if check_running else 'source numerical matrix accepted without running input',
                'independent native-to-real scalar ordered reconstruction',
                'formal symmetric polynomial source equality'],
            'checks_not_passed':[
                'installed native Matchete GammaCC PL -> declared Weyl epsilon absolute sign/phase',
                'connected UV 4-fermion two-scalar amplitude sign and absolute coefficient',
                'native Matchete SU2 real triplet generator axes -> master RGE generator axes',
                'actual three-generation master fermion constructor',
                'physical normalization of C_ijab in master RGE convention'],
            'native_probe_required':{
                'minimal_action':'In a fresh Matchete kernel, declare the exact A(-4) F Dirac mass and both Yukawa vertices; integrate only F at tree level; evaluate one nonzero PL component between two independent anticommuting external Weyl spinors with explicitly declared C=i gamma^2 gamma^0 and epsilon12=+1. Compare UV exchange and matched contact amplitudes, including the phase, overall sign, normalization, p/q assignment and all LSZ/identical-leg factors.',
                'independent_cross_check':'Repeat with C inverse and with both spinor labels exchanged; verify PR is Hermitian conjugate of PL and measure native SU2 triplet-generator basis matrices.',
                'required_numeric_output':['UV four-point coefficient (Re,Im) in fixed convention',
                                          'Matchete matched coefficient (Re,Im) in same convention',
                                          'GammaCC concrete 4x4 matrix',
                                          'master Weyl epsilon and explicit spinor components',
                                          'native and master SU2 generators in same component basis']}}


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--out',type=Path,default=OUTPUT)
    p.add_argument('--skip-running-check',action='store_true',
                   help='For isolated fixtures only; report will mark UV flavour unchecked')
    args=p.parse_args(argv)
    root=args.root.resolve()
    dest=args.out if args.out.is_absolute() else root/args.out
    try:
        obj=verify_from_sources(root,check_running=not args.skip_running_check)
    except (BoundaryBlocked,FileNotFoundError,ValueError,KeyError,ImportError) as exc:
        print('A(-4) boundary BLOCKED:',exc)
        return 2
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('A(-4) conditional tensor built:',dest)
    print('Physical certification: BLOCKED (native GammaCC phase / field-basis dictionary)')
    print('Independent conditional entries:',len(obj['independent_conditional_tensor_entries']))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
