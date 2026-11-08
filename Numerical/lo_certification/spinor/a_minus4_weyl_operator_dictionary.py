"""A(-4) Batch 1: Grassmann-level Matchete PL -> master Weyl dictionary.

This is a source-anchored, CLOSED-CHOICE dictionary, not an independent
measurement of the installed Matchete GammaCC matrix. It compares the actual
native PL coefficient, with source flavour index ordering and the UV->master
triplet transformation, against the symmetric full-sum 45x45x12x12 master
candidate for 12 independent anticommuting lepton spinor components.

A declared Weyl basis has epsilon_12=+1, C0=i gamma^2 gamma^0, metric +---:
    NCM[Transp[L_q],GammaCC PL,L_p]  ->  rho * (ell_q^T epsilon ell_p).
Both rho=+1 and -1 are tested, and neither is selected as native.

No RGE, scalar matching or physical Wilson normalization is implemented.
Run from repository root:
 python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_weyl_operator_dictionary --root .
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

import numpy as np

from . import a_minus4_boundary_verify as boundary
from . import a_minus4_connected_amplitude_gate as connected
from . import a_minus4_master_basis_gate as master

FORMAT = 't3-a-minus4-grassmann-weyl-operator-dictionary-v1'
DEFAULT_MASTER = master.DEFAULT_OUT
DEFAULT_OUT = Path('output/general_observable/a_minus4_boundary_batch1/weyl_operator_dictionary.json')


class DictionaryBlocked(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise DictionaryBlocked(message)


def load(path: Path) -> dict:
    require(path.is_file(), 'Missing prerequisite '+str(path))
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict),'Invalid JSON prerequisite '+str(path))
    return obj


def grassmann_product(i: int, j: int) -> dict[tuple[int,int], int]:
    """e_i e_j = -e_j e_i; 0 for repeated generator."""
    if i==j:
        return {}
    return {(min(i,j),max(i,j)):1 if i<j else -1}


def weyl_epsilon(i: int, j: int) -> dict[tuple[int,int],int]:
    """ell_i^0 ell_j^1 - ell_i^1 ell_j^0, with anticommuting comps.

    Global Weyl-component generators are indexed 0..11, two per L doublet
    SU2/flavour field i=1..6. W(ell_i,ell_j)=W(ell_j,ell_i).
    For i=j, W(ell_i,ell_i)=2 ell_i^0 ell_i^1 (NOT zero).
    """
    require(1<=i<=6 and 1<=j<=6,'Expected left-lepton master labels 1..6')
    coeff:dict[tuple[int,int],int]=defaultdict(int)
    for sign,a,b in ((1,2*(i-1),2*(j-1)+1),
                     (-1,2*(i-1)+1,2*(j-1))):
        for key,z in grassmann_product(a,b).items():
            coeff[key]+=sign*z
    return {k:v for k,v in coeff.items() if v}


def exterior_summed(rows: Mapping[tuple[int,int],complex]) -> dict[tuple[int,int],complex]:
    """Rows label Weyl *fields*; output is exact degree-two Grassmann basis."""
    out:dict[tuple[int,int],complex]=defaultdict(complex)
    for (i,j),coef in rows.items():
        for key,v in weyl_epsilon(i,j).items():
            out[key]+=complex(coef)*v
    return {k:v for k,v in out.items() if v != 0}


def native_pl_operator(w: np.ndarray, g: np.ndarray, native_phi: np.ndarray) -> dict:
    """DIRECT native PL source; does NOT use the symmetric Wilson tensor.

    PL has (q,b) on the TRANSPOSED left and (p,a) on the right. Its
    Yukawa factor is W[p,q] = Sum_r conj(y1[p,r])*conj(y2[q,r])/Mr.
    The *native* singlet/triplet fields are S1* and S2, respectively.
    """
    require(w.shape==(3,3) and g.shape==(2,2,3) and native_phi.shape==(12,),
            'Wrong native PL component shapes')
    s1star=(native_phi[4]-1j*native_phi[5])/np.sqrt(2.)
    s2=np.array([(native_phi[6+2*c]+1j*native_phi[7+2*c])/np.sqrt(2.)
                 for c in range(3)],complex)
    terms:dict[tuple[int,int],complex]=defaultdict(complex)
    for p in range(3):
        for q in range(3):
            for a in range(2):
                for b in range(2):
                    for c in range(3):
                        z=w[p,q]*g[a,b,c]*s1star*s2[c]
                        if z != 0:
                            # The source's NCM order is q then p.
                            first=boundary.global_lepton(b+1,q+1)
                            second=boundary.global_lepton(a+1,p+1)
                            terms[first,second]+=z
    return exterior_summed(terms)


def master_fullsum_operator(monomials: Mapping[tuple[int,int,int,int],complex],
                            master_phi: np.ndarray, *, coefficient_factor: complex=1) -> dict:
    """Each monomial is ORBIT-SUM, not per-orbit coefficient.

    C_sym for each distinct orbit member is monomial/orbit_size. Thus
    sum_{i,j,a,b} C_sym[i,j,a,b] * (ell_i epsilon ell_j) * phi_a phi_b
    equals exactly this independent-monomial evaluation. No extra /2 or /4.
    """
    require(master_phi.shape==(12,), 'Expected twelve interleaved real scalars')
    terms:dict[tuple[int,int],complex]=defaultdict(complex)
    for (i,j,a,b),z in monomials.items():
        require(1<=i<=j<=6 and 5<=a<=6 and 7<=b<=12,
                'Not an A(-4) independent Wilson coordinate')
        terms[i,j]+=complex(coefficient_factor)*z*master_phi[a-1]*master_phi[b-1]
    return exterior_summed(terms)


def diff(a: Mapping, b: Mapping) -> tuple[float,float]:
    keys=set(a)|set(b)
    absmax=max((abs(a.get(k,0)-b.get(k,0)) for k in keys),default=0.)
    scale=max((abs(v) for v in a.values()),default=0.)
    return float(absmax),float(absmax/max(scale,1e-300))


def full_grassmann_probe(w: np.ndarray, cg_native: np.ndarray,
                          triplet_n: np.ndarray, monomials: Mapping,
                          *, seed: int=20261108) -> dict:
    """Test every degree-2 component of the *full* three-family PL operator.

    Spinors are exterior-algebra generators (NOT commuting trial values).
    Scalar fields are real commuting test values. Test complex W without
    assuming any lepton-flavour symmetry and verify two independent CC signs.
    """
    require(triplet_n.shape==(3,3),'Expected 3x3 complex triplet map')
    r=master.realification(triplet_n)
    rng=np.random.default_rng(seed)
    residual=0.0
    native_norm=0.0
    sign_branch=[]
    flavour_antisym_residual=0.0
    for trial in range(8):
        master_phi=rng.normal(size=12)
        native_phi=master_phi.copy()
        native_phi[6:12]=r.T@master_phi[6:12]
        native=native_pl_operator(w,cg_native,native_phi)
        projected=native_pl_operator((w+w.T)/2,cg_native,native_phi)
        _,flavour_relative=diff(native,projected)
        flavour_antisym_residual=max(flavour_antisym_residual,flavour_relative)
        require(flavour_relative<3e-12,
                'Symmetric lepton CG did not cancel antisymmetric flavour coefficient')
        fitted=master_fullsum_operator(monomials,master_phi)
        require(len(native)>0 and len(fitted)>0, 'Empty Grassmann source/operator')
        native_norm=max(native_norm,max(abs(z) for z in native.values()))
        _,rel=diff(native,fitted)
        residual=max(residual,rel)
        require(rel < 3e-12,
                f'Native PL and complete master Grassmann operator differ, trial {trial}, rel {rel}')
        # GammaCC = rho C0 is a dictionary choice: SAME rho multiplies
        # original Matchete PL and the mapped master Wilson operator.
        for rho in (1,-1):
            left={k:rho*z for k,z in native.items()}
            right={k:rho*z for k,z in fitted.items()}
            _,signed=diff(left,right)
            require(signed < 3e-12, 'Spinor-phase transport mismatch')
            if trial==0:
                sign_branch.append({'rho_C0':rho,'relative_residual':signed})
        if trial==0:
            _,half=diff(native,master_fullsum_operator(monomials,master_phi,
                                                       coefficient_factor=0.5))
            _,cross=diff(native,master_fullsum_operator(monomials,master_phi,
                                                        coefficient_factor=-1))
            require(half > 0.25 and cross > 1.5,
                    'Wrong 1/2 or opposite CC sign unexpectedly accepted')
    return {'samples':8,
            'full_degree_two_Grassmann_component_count':len(native),
            'max_relative_Grassmann_operator_residual':residual,
            'native_source_nonzero_max_abs_gev_inv':native_norm,
            'antisymmetric_flavour_W_cancellation_relative_residual':flavour_antisym_residual,
            'rho_branches':[{'name':('C0' if r['rho_C0']==1 else '-C0'),**r}
                            for r in sign_branch],
            'rho_sign_discriminated':False,
            'extra_half_factor_rejected':True,
            'uncorrelated_opposite_CC_sign_rejected':True}


def verify_evidence(root: Path, native_path: Path, audit_path: Path,
                    connected_path: Path, master_path: Path) -> dict:
    root=root.resolve()
    def actual(p:Path)->Path: return p if p.is_absolute() else root/p
    # Recompute from raw Matchete tree, Yukawas, CGs, mass basis and source
    # polynomial; do not trust the saved master PASS flag.
    rebuilt=master.verify_evidence(root,native_path,audit_path,connected_path)
    saved=load(actual(master_path))
    require(saved.get('format')==master.FORMAT and saved.get('model')=='A_-4' and
            saved.get('source_hashes_sha256')==rebuilt['source_hashes_sha256'] and
            saved.get('physical_C_ijab_at_MF') is None and
            saved.get('physical_certification') is False,
            'Master candidate provenance/status invalid')
    rows_saved=saved.get('conditional_master_weight_basis_independent_entries')
    rows_fresh=rebuilt['conditional_master_weight_basis_independent_entries']
    require(isinstance(rows_saved,list) and len(rows_saved)==len(rows_fresh)==36,
            'Expected 36 independent A(-4) master candidate coordinates')
    require([tuple(x[k] for k in ('i','j','a','b','distinct_orbit_size')) for x in rows_saved]
            ==[tuple(x[k] for k in ('i','j','a','b','distinct_orbit_size')) for x in rows_fresh],
            'Master independent tensor keys changed')
    for a,b in zip(rows_saved,rows_fresh):
        for key in ('coefficient_real_if_rho_1','coefficient_imag_if_rho_1'):
            require(abs(float(a[key])-float(b[key]))<=1e-15*max(abs(float(b[key])),1e-50),
                    'Saved master coefficient differs from source-first reconstruction')

    source=boundary.read(root/boundary.SOURCE)
    w=boundary.decode_matrix(source['operator_terms'][0]['ordered_flavour_matrix_gev_inv'])
    seed=boundary.read(root/connected.MATCHED_SEED)
    cg,_=boundary.parse_native_cg(seed['CGRegistry'])
    n=master.triplet_intertwiner(cg)['triplet_N']
    monomials=master.independent_to_monomials(rows_fresh)
    check=full_grassmann_probe(w,cg,n,monomials)
    # The source benchmark has a diagonal flavour matrix: also check a fully
    # general non-Hermitian 3x3 Yukawa-derived mass contraction with same CG.
    rng=np.random.default_rng(20261109)
    y1=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/5
    y2=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/6
    masses=np.array([111.,173.,239.])
    w_general=boundary.dense_flavour(y1,y2,masses)
    q_master=(cg.reshape(4,3)@n.conj().T).reshape(2,2,3)
    indep,mono_general=boundary.symmetrise(boundary.raw_components(w_general,q_master))
    synthetic=full_grassmann_probe(w_general,cg,n,mono_general,seed=20261109)
    # Numeric check of 4x4 Weyl blocks, separate from Matchete abstract GammaCC.
    z=np.zeros((2,2),complex); ident=np.eye(2,dtype=complex)
    sigma2=np.array([[0,-1j],[1j,0]],complex)
    gam0=np.block([[z,ident],[ident,z]])
    gam2=np.block([[z,sigma2],[-sigma2,z]])
    pl=np.diag([1,1,0,0])
    eps=np.array([[0,1],[-1,0]],complex)
    c0=1j*gam2@gam0
    require(np.max(np.abs((c0@pl)[:2,:2]-eps))<1e-15 and
            np.max(np.abs((-c0@pl)[:2,:2]+eps))<1e-15,
            'Declared Weyl CC sign branches changed')
    require(saved.get('overall_matchete_GammaCC_to_master_Weyl_factor') is None,
            'Saved master source falsely selected native charge-conjugation phase')
    return {
        'format':FORMAT,'model':'A_-4',
        'status':'SOURCE_MATCHED_FULL_GRASSMANN_OPERATOR_DICTIONARY_VERIFIED_NATIVE_CC_CHOICE_OPEN',
        'source_revision':rebuilt['source_revision'],
        'mu_F_gev':rebuilt['mu_F_gev'],
        'input_master_candidate_sha256':sha256(actual(master_path).read_bytes()).hexdigest(),
        'source_native_files_sha256':rebuilt['source_hashes_sha256'],
        'operator_identity':{
            'native_PL':'Sum[p,q,a,b,c] W[p,q] CG[a,b,c] conj(S1) S2_native[c] (ell_[b,q]^T epsilon ell_[a,p])',
            'ordered_W':'W[p,q] = Sum_r conj(y1[p,r])*conj(y2[q,r])/M_r',
            'spinor_order':'q on transposed left, p on right; exterior Weyl contraction is symmetric under exchanging whole lepton fields',
            'flavour_antisymmetric_part':'(W-W^T)/2 cancels because matched SU2 CG and Grassmann Weyl bilinear are symmetric under simultaneous exchange of lepton legs; native W[p,q] ordering still verified independently',
            'weyl_epsilon':'epsilon_12=+1; ell_i^T epsilon ell_j = ell_j^T epsilon ell_i for anticommuting fields',
            'triplet_basis':'S2_master=N*S2_native; S2_native=N^dagger*S2_master',
            'complex_scalar_fields':'conj(S1)=(phi5-i phi6)/sqrt(2); S2_master[c]=(phi_(7+2c)+i phi_(8+2c))/sqrt(2) for zero-based c',
            'master_fullsum':'Sum_{i,j=1..45;a,b=1..12} C_sym^formal[i,j,a,b] (ell_i^T epsilon ell_j) phi_a phi_b',
            'master_independent':'Only i<=j, a<=b entries stored; C_sym^formal=monomial_sum/distinct_orbit_size; no additional prefactor in the explicitly declared full sum',
            'C0':'C0=i gamma^2 gamma^0; (C0 PL)_LH=epsilon; (-C0 PL)_LH=-epsilon',
            'required_relation':'C_sym^declared(rho) = rho * C_sym^formal for rho=+1 or -1; no physical branch selected',
        },
        'checks':{
            'actual_evolved_three_flavour_source_full_Grassmann':check,
            'independent_general_complex_3x3_source_full_Grassmann':synthetic,
            'matched_source_term_and_master_rebuilt':True,
            'master_independent_coordinate_count':len(rows_fresh),
            'no_missing_fermion_swap_half_factor':True,
            'explicit_weyl_CC_candidates_both_valid':True,
        },
        'spinor_dictionary_status':{
            'GammaCC_symbol_in_installed_Matchete_algebra':'abstract, no concrete 4x4 matrix provided',
            'reference_4x4_charge_conjugation':'C0=i gamma^2 gamma^0 chosen explicitly by this dictionary, not native-selected',
            'rho_GammaCC_to_reference_Weyl':None,
            'native_to_master_Yukawa_and_mass_phase_identification_complete':False,
            'source_operator_identity_with_declared_rho_complete':True,
            'full_physical_C_ijab_at_MF':None,
            'physical_certification':False,
        },
        'next_decision':'Absolute native GammaCC sign is a convention choice, not fixed by abstract algebra; do not run further source-only sign tests. A physical coupling dictionary requires a declared consistent convention plus independent native-to-master action alignment.',
        'EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars simultaneously -> SMEFT'
    }


def main(argv:list[str]|None=None) -> int:
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--root',type=Path,default=Path.cwd())
    cli.add_argument('--native',type=Path,default=connected.DEFAULT_NATIVE)
    cli.add_argument('--audit',type=Path,default=connected.DEFAULT_AUDIT)
    cli.add_argument('--connected',type=Path,default=connected.DEFAULT_OUT)
    cli.add_argument('--master',type=Path,default=DEFAULT_MASTER)
    cli.add_argument('--out',type=Path,default=DEFAULT_OUT)
    args=cli.parse_args(argv)
    try:
        result=verify_evidence(args.root,args.native,args.audit,args.connected,args.master)
    except (DictionaryBlocked,master.BasisBlocked,connected.GateBlocked,
            boundary.BoundaryBlocked,ValueError,KeyError,OSError) as e:
        print('A(-4) Weyl dictionary BLOCKED:',e)
        return 2
    out=args.out if args.out.is_absolute() else args.root.resolve()/args.out
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Matchete PL -> full master Grassmann operator: PASS')
    print('Flavour p/q, scalar sqrt(2), fermion exchange and orbit factors: PASS')
    print('Native GammaCC sign: UNMEASURED; both declared choices pass')
    print('Physical C_ijab(M_F): NULL; not certified')
    print('Report:',out)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
