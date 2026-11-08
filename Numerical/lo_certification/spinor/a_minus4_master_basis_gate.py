"""A(-4) Batch-1: Matchete triplet -> master SU(2) *candidate* basis.

Uses the F-only LLSS seed, BOTH registered UV Yukawa CGs and the actual master
RGE SU(2) generator functions. Constructs and verifies the unitary triplet
intertwiner, real 6x6 scalar map and the full 3-flavour symmetric *formal*
Wilson coordinate transport. Verifies the 45-Weyl index extension using the
existing FermionBasis API; does not modify the production build_sm_eft().

The Matchete-to-master fundamental-doublet phase and GammaCC-to-Weyl spinor
factor are NOT measured by these tests; physical C_ijab remains null.
No RGBeta coupling dictionary or Wilson RGE is implemented.

From repository root:
 python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_master_basis_gate \
   --root . \
   --native output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json \
   --audit output/general_observable/a_minus4_boundary_batch1/native_spinor_audit.json \
   --connected output/general_observable/a_minus4_boundary_batch1/connected_amplitude_gate.json \
   --out output/general_observable/a_minus4_boundary_batch1/master_basis_candidate.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import sympy as sp

from . import a_minus4_boundary_verify as boundary
from . import a_minus4_connected_amplitude_gate as connected

FORMAT = 't3-a-minus4-source-anchored-master-triplet-candidate-v1'
DEFAULT_NATIVE = Path('output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json')
DEFAULT_AUDIT = Path('output/general_observable/a_minus4_boundary_batch1/native_spinor_audit.json')
DEFAULT_CONNECTED = Path('output/general_observable/a_minus4_boundary_batch1/connected_amplitude_gate.json')
DEFAULT_OUT = Path('output/general_observable/a_minus4_boundary_batch1/master_basis_candidate.json')
SOURCE_GENERATORS = Path('output/general_observable/lo_direct_one_loop/source_audits/matchete_generator_analysis.json')


class BasisBlocked(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise BasisBlocked(message)


def read(path: Path) -> dict:
    require(path.is_file(), 'Missing input '+str(path))
    val = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(val, dict), 'Expected JSON object '+str(path))
    return val


def digest(path: Path) -> str:
    require(path.is_file(), 'Missing source '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def encode_complex(m: np.ndarray) -> dict:
    x = np.asarray(m, complex)
    require(np.all(np.isfinite(x)), 'Nonfinite complex array')
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def normalized_pair_embedding() -> tuple[np.ndarray,np.ndarray]:
    rt2 = np.sqrt(2.)
    e = np.array([[1,0,0],[0,1/rt2,0],[0,1/rt2,0],[0,0,1]], complex)
    j = np.array([[0,0,1],[0,-1,0],[1,0,0]], complex)
    return e,j


def realification(u: np.ndarray) -> np.ndarray:
    """Complex x'=U x -> interleaved real coordinates x'_R=R(U)x_R."""
    require(u.ndim == 2 and u.shape[0] == u.shape[1], 'Expected square complex basis change')
    n = len(u)
    r = np.zeros((2*n,2*n),float)
    for a in range(n):
        for b in range(n):
            z = complex(u[a,b])
            r[2*a:2*a+2, 2*b:2*b+2] = [[z.real,-z.imag],[z.imag,z.real]]
    return r


def triplet_intertwiner(gmatch: np.ndarray, *, master_generators=None,
                        native_cartan: np.ndarray|None=None) -> dict:
    """Compute S_master = N S_native from *actual* matched CG.

    CG is for L_i L_j S2_native[c], i,j in 2. The S2 field transforms in
    the conjugate of the symmetric doublet-pair representation. E spans the
    symmetric pair; J fixes the standard spin-1 weight basis. This is
    conditional on identifying Matchete's fundamental axes with the master's.
    """
    from RGE.general.GaugeGenerators import su2_complex_generators

    q = np.asarray(gmatch,complex)
    require(q.shape == (2,2,3), 'Expected registered A(-4) 2x2x3 CG')
    require(np.max(np.abs(q - np.swapaxes(q,0,1))) < 3e-13,
            'Matched native CG not symmetric in lepton axes')
    q4 = q.reshape(4,3)
    e,j = normalized_pair_embedding()
    gram = q4.conj().T @ q4
    norm2 = float(np.trace(gram).real/3)
    require(norm2 > 0 and np.max(np.abs(gram-norm2*np.eye(3))) < 3e-13,
            'Native CG triplet Gram matrix changed')
    f = q4/np.sqrt(norm2)
    v = e.conj().T@f
    require(np.max(np.abs(e@v-f)) < 3e-13 and
            np.max(np.abs(v.conj().T@v-np.eye(3))) < 3e-13,
            'Native CG does not span the symmetric pair')
    n = j@v
    require(np.max(np.abs(n.conj().T@n-np.eye(3))) < 3e-13,
            'Native-to-master triplet map not unitary')

    if master_generators is None:
        master_generators = su2_complex_generators(3)
    tmaster = [np.array(x,complex) for x in master_generators]
    t2 = [np.array(x,complex) for x in su2_complex_generators(2)]
    d = [np.kron(t.T,np.eye(2))+np.kron(np.eye(2),t.T) for t in t2]
    native_gens = [-f.conj().T@dk@f for dk in d]
    cgcov = max(float(np.max(np.abs(dk@f+f@gn))) for dk,gn in zip(d,native_gens))
    generr = max(float(np.max(np.abs(n@gn@n.conj().T-tm)))
                 for gn,tm in zip(native_gens,tmaster))
    require(cgcov < 3e-13 and generr < 3e-13,
            'Source CG and actual master generators not covariant')
    # Native-to-master coefficients transform as Q_master = Q_native N^dagger.
    qm = q4@n.conj().T
    canonical_q = np.sqrt(norm2)*e@j
    cgerr = float(np.max(np.abs(qm-canonical_q)))
    require(cgerr < 3e-13, 'Mapped Matchete CG not the standard weight-basis CG')
    r = realification(n)
    require(np.max(np.abs(r.T@r-np.eye(6))) < 3e-13,
            'Interleaved scalar real map is not orthogonal')
    if native_cartan is not None:
        require(np.asarray(native_cartan).shape == (3,3) and
                np.max(np.abs(native_gens[2]-native_cartan)) < 3e-13,
                'Independent Matchete triplet Cartan differs from CG-derived generator')
    return {'triplet_N':n, 'real_R':r, 'q_master':qm, 'native_generators':native_gens,
            'generator_residual':generr, 'cg_covariance_residual':cgcov,
            'cg_standard_residual':cgerr, 'native_CG_norm_sq':norm2}


def independent_to_monomials(records: list[dict]) -> dict[tuple[int,int,int,int],complex]:
    """Extract formal polynomial monomial sums (not physical Weyl amplitudes)."""
    result: dict[tuple[int,int,int,int],complex] = {}
    for rec in records:
        i,j,a,b = (int(rec[x]) for x in ('i','j','a','b'))
        require(1 <= i <= j <= 45 and 1 <= a <= b <= 12 and
                i <= 6 and j <= 6 and 5 <= a <= 6 and 7 <= b <= 12,
                'Not an A(-4) LL(S1bar)(S2) independent coordinate')
        k=(i,j,a,b)
        require(k not in result,'Duplicated independent coordinate')
        orb=len(boundary.orbit(*k))
        require(rec['distinct_orbit_size'] == orb,
                'Source symmetric orbit combinatorics changed')
        z=complex(rec['coefficient_real_if_rho_1'],rec['coefficient_imag_if_rho_1'])
        result[k]=orb*z
    require(result and all(np.isfinite(v) for v in result.values()),
            'Source-independent coefficients empty/nonfinite')
    return result


def transport_monomials(monomials: Mapping[tuple[int,int,int,int],complex],
                        real_map: np.ndarray) -> dict[tuple[int,int,int,int],complex]:
    """Substitute S_native=N^dagger S_master into the formal LLSS polynomial."""
    require(real_map.shape == (6,6), 'Unexpected S2 real-coordinate map')
    n_to_m = real_map.T  # old x_native=R(N^dagger) x_master
    out: dict[tuple[int,int,int,int],complex] = defaultdict(complex)
    for (i,j,a,b), z in monomials.items():
        require(1 <= i <= j <= 45 and a in (5,6) and b in range(7,13),
                'Transport input outside A(-4) light scalar support')
        for new_s2 in range(6):
            val=complex(z)*n_to_m[b-7,new_s2]
            if abs(val)>0:
                out[i,j,a,7+new_s2] += val
    return {k:v for k,v in out.items() if abs(v) > 1e-25}


def monomials_to_independent(monomials: Mapping[tuple[int,int,int,int],complex]) -> list[dict]:
    result=[]
    for k,z in sorted(monomials.items()):
        m=len(boundary.orbit(*k))
        c=z/m
        result.append({'i':k[0],'j':k[1],'a':k[2],'b':k[3],
                       'coefficient_real_if_rho_1':float(c.real),
                       'coefficient_imag_if_rho_1':float(c.imag),
                       'distinct_orbit_size':m})
    return result


def eval_monomials(monomials: Mapping[tuple[int,int,int,int],complex],
                   fermions: np.ndarray, scalars: np.ndarray) -> complex:
    """Commuting-symbol polynomial witness only; not an external spinor amplitude."""
    total=0j
    for (i,j,a,b),z in monomials.items():
        total += z*fermions[i-1]*fermions[j-1]*scalars[a-1]*scalars[b-1]
    return complex(total)


def verify_transport(old: dict, new: dict, r: np.ndarray) -> dict:
    rng=np.random.default_rng(20261108)
    maxres=0.0
    for _ in range(8):
        fermions=rng.normal(size=45)+1j*rng.normal(size=45)
        master=rng.normal(size=12)
        native=master.copy()
        native[6:12]=r.T@master[6:12]
        original=eval_monomials(old,fermions,native)
        transformed=eval_monomials(new,fermions,master)
        residual=abs(original-transformed)/max(abs(original),abs(transformed),1e-50)
        maxres=max(maxres,residual)
        require(residual < 2e-12, 'Master scalar map changes original source polynomial')
    def norm_squared(s):
        return float(sum(abs(z)**2/len(boundary.orbit(*k)) for k,z in s.items()))
    nold,nnew=norm_squared(old),norm_squared(new)
    normerr=abs(nold-nnew)/max(nold,nnew,1e-300)
    require(normerr<2e-12, 'Orthogonal triplet rotation changed full symmetric tensor norm')
    return {'commuting_polynomial_relative_max_residual':maxres,
            'symmetric_full_coordinate_norm_relative_residual':normerr,
            'native_independent_monomials':len(old),
            'master_independent_monomials':len(new),
            'source_tensor_frobenius_norm_gev_inv':float(np.sqrt(nold)),
            'master_tensor_frobenius_norm_gev_inv':float(np.sqrt(nnew)),
            'physical_spinor_operator_identity_established':False}


def three_family_api_basis(expected: dict) -> dict:
    """Construct (do not alter) the master API's three-copy SM basis."""
    from RGE.general.FermionBasis import FermionBasis, WeylFermion
    from RGE.general.ScalarBasis import ScalarBasis
    from RGE.matching.WeinbergTensorAdapter import build_sm_eft

    specifications=[('L',2,sp.Rational(-1,2)),('eC',1,sp.S.One)]
    for colour in (1,2,3):
        specifications.extend([(f'Q{colour}',2,sp.Rational(1,6)),
                               (f'uC{colour}',1,sp.Rational(-2,3)),
                               (f'dC{colour}',1,sp.Rational(1,3))])
    basis=FermionBasis(tuple(WeylFermion(name,d,y,multiplicity=3)
                             for name,d,y in specifications))
    scalar=ScalarBasis.t3(1,sp.Integer(-2),3,sp.Integer(-1))
    require(basis.dimension==45 and scalar.total_real_scalar_dimension==12,
            'Master public API dimension changed')
    require([(x['index'],x['field'],x['generation'],x['component'])
             for x in expected['fermions']] ==
            [(basis.global_index(name,comp,copy),name,copy,comp)
             for name,d,y in specifications for copy in range(1,4)
             for comp in range(1,d+1)],
            'Generated conditional fermion order not master FermionBasis ordering')
    for x in expected['fermions']:
        mat=next(s[2] for s in specifications if s[0]==x['field'])
        require(abs(float(mat)-float(x['hypercharge']))<1e-14,
                'Generated basis hypercharges differ from master API')
    for record in expected['real_scalars']:
        block=scalar.block(record['field'])
        pos=block.local_to_global(2*(record['component']-1)+(1 if record['part']=='Re' else 2))
        require(record['index']==pos,'Generated scalar coordinates differ from master ScalarBasis')
    one_s,one_f=build_sm_eft()
    require(one_f.dimension==15 and one_s.total_real_scalar_dimension==4,
            'Production one-family build_sm_eft changed; inspect before use')
    return {'candidate_45_Weyl_order_matches_FermionBasis_API':True,
            'scalar_12_real_order_matches_ScalarBasis_API':True,
            'production_build_sm_eft_fermion_dimension':one_f.dimension,
            'production_build_sm_eft_is_three_family':False,
            'full_tensor_shape':[45,45,12,12]}


def phase_covariance_probe() -> dict:
    """Test Dirac heavy rephasings & three independent external lepton phases.

    Work with source currents A(q,r), B(p,r) of the declared Weyl saddle,
    so no unstated Matchete/CConj->Weyl phase is inserted. Each term uses
    A_q B_p / M_r and all flavour indices remain active and complex.
    """
    rng=np.random.default_rng(20261108)
    a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    b=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    masses=np.array([157.,193.,251.])
    original=(b/masses[None,:])@a.T
    alpha=rng.uniform(-np.pi,np.pi,size=3)
    beta=rng.uniform(-np.pi,np.pi,size=3)
    theta=rng.uniform(-np.pi,np.pi,size=3)
    # Heavy eta'=exp(i alpha) eta, chi'=exp(i beta) chi,
    # A'=exp(-i alpha) A, B'=exp(-i beta) B, M'=exp(-i(alpha+beta)) M.
    ap=a*np.exp(-1j*alpha)[None,:]
    bp=b*np.exp(-1j*beta)[None,:]
    mp=masses*np.exp(-1j*(alpha+beta))
    transformed=(bp/mp[None,:])@ap.T
    he=np.max(np.abs(transformed-original))/max(np.max(np.abs(original)),1e-300)
    require(he<5e-14, 'Dirac mass/Yukawa heavy rephasing not invariant')
    # Positive mass-basis phases beta=-alpha preserve M>0.
    pos_ap=a*np.exp(-1j*alpha)[None,:]
    pos_bp=b*np.exp(1j*alpha)[None,:]
    positive=(pos_bp/masses[None,:])@pos_ap.T
    poserr=np.max(np.abs(positive-original))/max(np.max(np.abs(original)),1e-300)
    require(poserr<5e-14,'Positive-real heavy Dirac mass basis rephasing failed')
    # Field rephasing L_p'=exp(i theta_p) L_p.
    lep_ap=a*np.exp(-1j*theta)[:,None]
    lep_bp=b*np.exp(-1j*theta)[:,None]
    wl=(lep_bp/masses[None,:])@lep_ap.T
    expected=original*np.exp(-1j*(theta[:,None]+theta[None,:]))
    le=np.max(np.abs(wl-expected))/max(np.max(np.abs(original)),1e-300)
    require(le<5e-14,'Lepton flavour rephasing transport failed')
    return {'arbitrary_complex_dirac_mass_rephasing_relative_residual':float(he),
            'positive_mass_eigenbasis_heavy_rephasing_relative_residual':float(poserr),
            'three_independent_lepton_phase_covariance_relative_residual':float(le),
            'no_convention_dependent_yukawa_sign_selected':True,
            'native_Matchete_to_master_rephasing_dictionary_certified':False}


def verify_evidence(root:Path, native:Path, audit:Path, connected_path:Path|None=None) -> dict:
    root=root.resolve()
    def local(path:Path)->Path: return path if path.is_absolute() else root/path
    native,audit=local(native),local(audit)
    previous=connected.verify_evidence(root,native,audit)  # source-first recomputation
    if connected_path is not None:
        saved=read(local(connected_path))
        require(saved.get('format')==connected.FORMAT and
                saved.get('model')=='A_-4' and
                saved.get('source_files_sha256')==previous['source_files_sha256'] and
                saved.get('native_Cconj_to_master_Weyl_factor') is None and
                saved.get('physical_certification') is False,
                'Connected report stale or falsely certified')
    result=boundary.verify_from_sources(root,check_running=True)
    matched=read(root/connected.MATCHED_SEED)
    uv=read(root/connected.UV_CG_PATH)
    cgs,proof=connected.native_cg_identity(uv,matched)
    # Require separate stored numeric Matchete Cartan as an independent check.
    generator_path=root/SOURCE_GENERATORS
    generator=read(generator_path)
    from Lagrangian.RunT3UVNativeBridge import exact_number
    trip=next((rec for rec in generator.get('representations',[])
               if rec.get('dimension')==3), None)
    require(isinstance(trip,dict) and len(trip.get('T3_in_original_matchete_basis',[]))==3,
            'Missing independent native Matchete triplet Cartan matrix')
    c3=np.array([[complex(exact_number(s)) for s in row]
                 for row in trip['T3_in_original_matchete_basis']],complex)
    transform=triplet_intertwiner(cgs['Gmatch'],native_cartan=c3)
    api=three_family_api_basis(result['basis_contract'])
    old=independent_to_monomials(result['independent_conditional_tensor_entries'])
    # Independent numerical construction directly from the *mapped* CG,
    # rather than thresholding small flavour coefficients from a rotation.
    e,j=normalized_pair_embedding()
    q_exact_zero_pattern=np.sqrt(transform['native_CG_norm_sq'])*(e@j)
    require(np.max(np.abs(q_exact_zero_pattern-transform['q_master']))<3e-13,
            'CG-mapped exact support does not agree with source')
    source_flavour=boundary.read(root/boundary.SOURCE)
    w=boundary.decode_matrix(source_flavour['operator_terms'][0]['ordered_flavour_matrix_gev_inv'])
    master_native=boundary.raw_components(w,q_exact_zero_pattern.reshape(2,2,3))
    master_independent, new=boundary.symmetrise(master_native)
    rotated=transport_monomials(old,transform['real_R'])
    all_keys=set(rotated)|set(new)
    scalar_trans_max=max((abs(rotated.get(k,0)-new.get(k,0)) for k in all_keys),default=0.)
    scale=max((abs(v) for v in new.values()),default=0.)
    require(scale>0 and scalar_trans_max/max(scale,1e-300)<3e-12,
            'Rotated source coefficients differ from independent master CG projection')
    transported=verify_transport(old,new,transform['real_R'])
    transported['direct_mapped_CG_versus_rotated_tensor_relative_max_residual']=float(scalar_trans_max/scale)
    phases=phase_covariance_probe()
    require(previous['source_native_SU2_CG_contraction']['UV_Y1_Y2_CG_to_Matchete_F_only_max_abs_residual']<1e-13,
            'Prior source CG-to-matching identity did not pass')
    return {'format':FORMAT,'model':'A_-4',
            'status':'SOURCE_TRIPLET_TO_MASTER_CANDIDATE_AND_FLAVOUR_PHASE_COVARIANCE_VERIFIED_PHYSICAL_CC_OPEN',
            'source_revision':result['git_revision'],
            'mu_F_gev':result['numeric_flavour_checks']['mu_F_gev'],
            'model_F_only':True,
            'source_hashes_sha256':{'native':digest(native),'audit':digest(audit),
               'UV_registered_CG':digest(root/connected.UV_CG_PATH),
               'F_only_matched_seed':digest(root/connected.MATCHED_SEED),
               'native_matchete_generator_analysis':digest(generator_path)},
            'basis':api,
            'native_to_master_triplet_complex_matrix':encode_complex(transform['triplet_N']),
            'native_to_master_triplet_real_matrix':transform['real_R'].tolist(),
            'transform_equations':{
               'complex_fields':'S2_master = N S2_native',
               'real_fields':'phi_S2_master = R(N) phi_S2_native, indices 7..12',
               'wilson_coefficients':'C_master^formal = C_native^formal with phi_native = R(N)^T phi_master',
               'lepton_SU2':'identified with same sigma_a/2 master fundamental axes (phase alignment still conditional)',
               'master_su2_triplet':'T_a=su2_complex_generators(3)[a]',
               'native_su2_triplet':'derived by covariance from matching-generated symmetric LLSS CG'},
            'checks':{
               'native_triplet_cartan_from_independent_generator_source_equal':True,
               'native_to_master_generators_max_abs_residual':transform['generator_residual'],
               'native_CG_covariance_max_abs_residual':transform['cg_covariance_residual'],
               'mapped_CG_to_master_weight_basis_max_abs_residual':transform['cg_standard_residual'],
               'native_CG_norm_squared':transform['native_CG_norm_sq'],
               'complex_three_flavour_phase_covariance':phases,
               'formal_polynomial_scalar_transport':transported},
            'conditional_master_weight_basis_independent_entries':monomials_to_independent(new),
            'conditional_master_weight_basis_zero_outside_independent_entries':True,
            'conditional_fermion_shape':[45,45], 'conditional_real_scalar_shape':[12,12],
            'tensor_symmetry':'C_ijab=C_jiab=C_ijba; independent coordinates i<=j, a<=b; no extra 1/4 in full sum',
            'overall_matchete_GammaCC_to_master_Weyl_factor':None,
            'physical_C_ijab_at_MF':None,
            'physical_certification':False,
            'remaining':[
               'Installed Matchete GammaCC is abstract; its absolute phase relative to a fixed master Weyl spinor convention is not measured.',
               'Native fundamental doublet and external lepton phase orientation relative to master is conditional on sigma_a/2 conventions.',
               'Master build_sm_eft() is still one-family: the verified 45-Weyl API embedding is not installed in the production RGE.',
               'No physical full three-family EFT1 Wilson RGE or coupling dictionary is asserted.'],
            'EFT_order':'UV -> integrate out ALL F -> EFT1(SM+S1+S2) -> integrate out BOTH scalars -> SMEFT'}


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--native',type=Path,default=DEFAULT_NATIVE)
    p.add_argument('--audit',type=Path,default=DEFAULT_AUDIT)
    p.add_argument('--connected',type=Path,default=DEFAULT_CONNECTED)
    p.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=p.parse_args(argv)
    try: result=verify_evidence(a.root,a.native,a.audit,a.connected)
    except (BasisBlocked, connected.GateBlocked, boundary.BoundaryBlocked,
            FileNotFoundError,ValueError,KeyError,ImportError,TypeError) as exc:
        print('A(-4) master basis gate BLOCKED:',exc)
        return 2
    dest=a.out if a.out.is_absolute() else a.root.resolve()/a.out
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Native->master SU2 triplet CG generator and independent native Cartan: PASS')
    print('Master 45-Weyl API index candidate + 12 real scalars: PASS')
    print('Formal Wilson polynomial, complex flavour/Dirac phase transport: PASS')
    print('Physical GammaCC/master Weyl factor: UNMEASURED (no certification)')
    print('Master-basis conditional coordinates:',len(result['conditional_master_weight_basis_independent_entries']))
    print('Report:',dest)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
