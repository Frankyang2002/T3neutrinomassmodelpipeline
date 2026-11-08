"""A(-4), Batch 1: EXPLICIT declared Weyl/UV coupling dictionary.

This is a reparametrisation of the three-flavour Matchete F-only EFT1
boundary.  It does not determine the installed Matchete GammaCC 4x4 matrix.

Declare C0 = i gamma^2 gamma^0, epsilon_12 = +1.  For a comparison branch
rho in {+1,-1}, GammaCC is *assumed* to represent rho*C0.  Define fields
    L_native[p] = exp(i alpha[p]) L_master[p]
    F_native[r] = exp(i beta[r]) F_master[r]  (both Dirac chiralities)
    S1_native   = exp(i sigma1) S1_master
    S2_native   = exp(i sigma2) N^dagger S2_master.
The source UV vertices are y1 * Bar[L] * CConj[F] * S1 and
                         y2 * Bar[L] * F * Bar[S2].
Consequently, at positive diagonal heavy masses,
    y1_master[p,r] = rho*y1_native[p,r] * exp(-i alpha[p]-i beta[r]+i sigma1)
    y2_master[q,r] =     y2_native[q,r] * exp(-i alpha[q]+i beta[r]-i sigma2)
    MF_master[r] = MF_native[r].
The rho in y1 maps the source's unique GammaCC vertex into declared C0.
The PL Wilson matrix is W[p,q] = sum_r conj(y1[p,r])*conj(y2[q,r])/MF[r].

Only rho=+1 with alpha=beta=sigma1=sigma2=0 is EXPLICITLY DECLARED as
our working numerical convention, not claimed to be Matchete's native C.
No physical Wilson boundary or matched scalar SMEFT is certified here.

Run:
 python -X utf8 -m Numerical.lo_certification.spinor.a_minus4_declared_convention_boundary --root .
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from . import a_minus4_boundary_verify as boundary
from . import a_minus4_connected_amplitude_gate as connected
from . import a_minus4_master_basis_gate as master
from . import a_minus4_weyl_operator_dictionary as dictionary

FORMAT = 't3-a-minus4-declared-weyl-uv-coupling-dictionary-v1'
DEFAULT_WITNESS = dictionary.DEFAULT_OUT
DEFAULT_OUT = Path('output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json')


class ConventionBlocked(ValueError):
    pass


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ConventionBlocked(message)


def valid_phase_vector(value: Any, n: int, label: str) -> np.ndarray:
    arr = np.asarray(value, float)
    require(arr.shape == (n,) and np.all(np.isfinite(arr)), f'{label}: expected {n} finite phases')
    return arr


def complex_symmetric_w(y1: np.ndarray, y2: np.ndarray, mass: np.ndarray) -> np.ndarray:
    """Source-order W, *not* symmetrised until lepton Grassmann contraction."""
    return boundary.dense_flavour(y1, y2, mass)


def uv_couplings_in_declared_basis(
    y1_native: np.ndarray, y2_native: np.ndarray, mass: np.ndarray,
    alpha: np.ndarray, beta: np.ndarray, sigma1: float, sigma2: float,
    *, rho: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full UV vertex phase transport; vectorlike F phases preserve real MF."""
    require(rho in (-1, 1) and type(rho) is int, 'rho must be exactly +1 or -1')
    y1, y2 = np.asarray(y1_native, complex), np.asarray(y2_native, complex)
    m = np.asarray(mass, float)
    a = valid_phase_vector(alpha, 3, 'alpha')
    b = valid_phase_vector(beta, 3, 'beta')
    require(y1.shape == y2.shape == (3,3) and m.shape == (3,) and
            np.all(np.isfinite(y1)) and np.all(np.isfinite(y2)) and
            np.all(np.isfinite(m)) and np.all(m > 0) and
            np.isfinite(sigma1) and np.isfinite(sigma2),
            'A(-4) requires complex 3x3 Yukawas, three positive masses and finite phases')
    phase1 = np.exp(1j*(-a[:,None]-b[None,:]+sigma1))
    phase2 = np.exp(1j*(-a[:,None]+b[None,:]-sigma2))
    return rho*y1*phase1, y2*phase2, m.copy()


def expected_pl_phase(w_native: np.ndarray, alpha: np.ndarray,
                      sigma1: float, sigma2: float, *, rho: int) -> np.ndarray:
    """Directly transform L_native L_native conj(S1_native) S2_native."""
    a=valid_phase_vector(alpha,3,'alpha')
    w=np.asarray(w_native,complex)
    require(w.shape==(3,3) and rho in (-1,1) and type(rho) is int,
            'Invalid native W or GammaCC choice')
    return rho*w*np.exp(1j*(a[:,None]+a[None,:]-sigma1+sigma2))


def scalar_native_from_master(master_phi: np.ndarray, triplet_n: np.ndarray,
                              sigma1: float, sigma2: float) -> np.ndarray:
    """S_native=phase*(N^dagger)S_master; real scalars are interleaved Re/Im."""
    phi=np.asarray(master_phi,float)
    n=np.asarray(triplet_n,complex)
    require(phi.shape==(12,) and n.shape==(3,3) and np.all(np.isfinite(phi)) and
            np.all(np.isfinite(n)) and np.isfinite(sigma1) and np.isfinite(sigma2),
            'Malformed 12-real-scalar/master triplet basis')
    require(np.max(np.abs(n.conj().T@n-np.eye(3)))<5e-13,
            'Nonunitary Matchete/master triplet map')
    out=phi.copy()
    out[4:6]=master.realification(np.array([[np.exp(1j*sigma1)]]))@phi[4:6]
    out[6:12]=master.realification(np.exp(1j*sigma2)*n.conj().T)@phi[6:12]
    return out


def grassmann_rephasing_witness(y1: np.ndarray, y2: np.ndarray, masses: np.ndarray,
                                cg: np.ndarray, triplet_n: np.ndarray,
                                *, seed: int = 20261108, samples: int = 5) -> dict:
    """Compare ALL degree-two exterior components, not commuting Weyl symbols."""
    rng=np.random.default_rng(seed)
    w=complex_symmetric_w(y1,y2,masses)
    qmaster=(cg.reshape(4,3)@triplet_n.conj().T).reshape(2,2,3)
    max_res=0.0
    max_w_res=0.0
    max_beta_res=0.0
    max_branch_res=0.0
    for trial in range(samples):
        alpha=rng.uniform(-2.9,2.9,3)
        beta=rng.uniform(-2.9,2.9,3)
        s1,s2=rng.uniform(-2.9,2.9,2)
        master_phi=rng.normal(size=12)
        native_phi=scalar_native_from_master(master_phi,triplet_n,s1,s2)
        for rho in (1,-1):
            ym1,ym2,mm=uv_couplings_in_declared_basis(y1,y2,masses,alpha,beta,s1,s2,rho=rho)
            wmaster=complex_symmetric_w(ym1,ym2,mm)
            wexp=expected_pl_phase(w,alpha,s1,s2,rho=rho)
            wres=np.max(np.abs(wmaster-wexp))/max(np.max(np.abs(wmaster)),1e-300)
            max_w_res=max(max_w_res,float(wres))
            require(wres<3e-12,'UV Yukawa/MF phases disagree with PL field transport')
            # Source has L_native = e^(i alpha)L_master, hence its external
            # spinors supply alpha[p]+alpha[q]. GammaCC conversion adds rho.
            w_ext=rho*w*np.exp(1j*(alpha[:,None]+alpha[None,:]))
            uv=dictionary.native_pl_operator(w_ext,cg,native_phi)
            _,mono=boundary.symmetrise(boundary.raw_components(wmaster,qmaster))
            eft=dictionary.master_fullsum_operator(mono,master_phi)
            err=dictionary.diff(uv,eft)[1]
            max_res=max(max_res,err)
            require(err<3e-12,'UV-to-EFT full Grassmann phase transport failed')
            # Independent vector-like heavy rephasing cannot alter W in
            # declared basis; it leaves the real positive MF unchanged.
            yb1,yb2,mb=uv_couplings_in_declared_basis(y1,y2,masses,alpha,
                                   beta+0.83,s1,s2,rho=rho)
            bres=np.max(np.abs(complex_symmetric_w(yb1,yb2,mb)-wmaster))/max(
                         np.max(np.abs(wmaster)),1e-300)
            max_beta_res=max(max_beta_res,float(bres))
            require(bres<3e-12,'Vectorlike Dirac phase failed to cancel')
        # Compare equivalent representations of SAME master UV couplings:
        # flipping GammaCC and flipping native y1 simultaneously preserves
        # mapped UV y1_master and the declared EFT contact polynomial.
        yp1,yp2,_=uv_couplings_in_declared_basis(y1,y2,masses,alpha,beta,s1,s2,rho=1)
        yn1,yn2,_=uv_couplings_in_declared_basis(-y1,y2,masses,alpha,beta,s1,s2,rho=-1)
        berr=max(np.max(np.abs(yp1-yn1)),np.max(np.abs(yp2-yn2)))/max(
            np.max(np.abs(yp1)),np.max(np.abs(yp2)),1e-300)
        max_branch_res=max(max_branch_res,float(berr))
        require(berr<3e-12,'Equivalent GammaCC/y1 redefinition not preserved')
    return {'complex_three_flavour_samples':samples,
            'max_full_grassmann_relative_residual':float(max_res),
            'max_w_matrix_relative_residual':float(max_w_res),
            'max_vectorlike_heavy_phase_relative_residual':float(max_beta_res),
            'max_equivalent_rho_y1_branch_relative_residual':float(max_branch_res),
            'independent_Dirac_heavy_phases_tested':3,
            'gamma_cc_branch_physically_discriminated':False}


def encode_rephased_coordinate(rows: list[dict], *, rho: int,
                              expected_count: int|None=None) -> list[dict]:
    """No 1/2 or 1/4 inserted: these are independent full-sum coordinates."""
    require(rho in (-1,1), 'Invalid declared branch')
    result=[]
    for rec in rows:
        i,j,a,b=(int(rec[k]) for k in ('i','j','a','b'))
        require(1<=i<=j<=6 and a in (5,6) and 7<=b<=12,
                'Unsupported non-A(-4) Wilson coordinate')
        require(rec['distinct_orbit_size']==len(boundary.orbit(i,j,a,b)),
                'Orbit normalization inconsistent')
        result.append({'i':i,'j':j,'a':a,'b':b,
            'value_re_gev_inv':float(rho*rec['coefficient_real_if_rho_1']),
            'value_im_gev_inv':float(rho*rec['coefficient_imag_if_rho_1']),
            'distinct_orbit_size':rec['distinct_orbit_size']})
    require(len(result)>0 and len({(q['i'],q['j'],q['a'],q['b']) for q in result})==len(result),
            'Duplicated or empty A(-4) independent coordinates')
    if expected_count is not None:
        require(len(result)==expected_count,
                'Unexpected number of independent A(-4) master coordinates')
    return result


def verify_sources(root:Path, native:Path, audit:Path, connected_file:Path,
                   master_file:Path, weyl_file:Path) -> dict:
    root=root.resolve()
    def abspath(p:Path)->Path: return p if p.is_absolute() else root/p
    # This invokes the prior independent source-first source/CG/Grassmann
    # gates rather than accepting an unproven status string.
    fresh=dictionary.verify_evidence(root,native,audit,connected_file,master_file)
    old=dictionary.load(abspath(weyl_file))
    require(old.get('format')==dictionary.FORMAT and
            old.get('model')=='A_-4' and
            old.get('source_native_files_sha256')==fresh['source_native_files_sha256'] and
            old.get('input_master_candidate_sha256')==fresh['input_master_candidate_sha256'] and
            old.get('spinor_dictionary_status',{}).get('physical_certification') is False,
            'Grassmann witness changed or overclaims physical normalization')
    native_report=dictionary.load(abspath(native))
    require(native_report.get('native_concrete_4x4_cc_matrix') is None and
            native_report.get('native_4x4_CC_to_Weyl_epsilon_ratio') is None,
            'New native GammaCC measurement requires source-first review')
    for key in ('native_yukawa_1_terms_inputform','native_yukawa_2_terms_inputform'):
        require(len(native_report.get(key,[]))==2,'Native UV Yukawa inventory incomplete')
    text_y1=''.join(native_report['native_yukawa_1_terms_inputform'])
    text_y2=''.join(native_report['native_yukawa_2_terms_inputform'])
    mass_text=native_report.get('native_heavy_free_inputform','')
    require('DiracProduct[GammaCC, Proj[1]]' in text_y1 and
            'DiracProduct[Proj[1]]' in text_y2 and
            'Coupling[MF,' in mass_text and
            '-(Coupling[MF,' in mass_text,
            'Native charge-conjugated y1 / standard y2 / Dirac mass changed')
    candidate=dictionary.load(abspath(master_file))
    rows=candidate.get('conditional_master_weight_basis_independent_entries')
    require(isinstance(rows,list) and len(rows)==36 and
            candidate.get('physical_certification') is False and
            candidate.get('overall_matchete_GammaCC_to_master_Weyl_factor') is None,
            'Master-basis boundary certificate changed')
    source=boundary.read(root/boundary.SOURCE)
    w=boundary.decode_matrix(source['operator_terms'][0]['ordered_flavour_matrix_gev_inv'])
    cg,_=boundary.parse_native_cg(boundary.read(root/connected.MATCHED_SEED)['CGRegistry'])
    n=master.triplet_intertwiner(cg)['triplet_N']
    uv_registry=boundary.read(root/connected.UV_CG_PATH)
    tensors,_=connected.native_cg_identity(uv_registry,
                               boundary.read(root/connected.MATCHED_SEED))
    # In the UV y2 vertex the field is Bar[S2_native], hence its CG
    # carries N[k,c] (not conj(N[k,c])).  The PL EFT has S2_native
    # unbarred and therefore carries conj(N[k,c]).
    g2_master=np.einsum('ahc,kc->ahk',tensors['G2'],n)
    gmatch_master=np.einsum('abc,kc->abk',tensors['Gmatch'],np.conjugate(n))
    uv_master=np.einsum('ah,bhk->abk',
                         np.conjugate(tensors['G1']), np.conjugate(g2_master))
    cg_res=float(np.max(np.abs(uv_master-gmatch_master)))
    require(cg_res<3e-13,'Bar[S2] UV triplet CG transport differs from PL contact')
    rng=np.random.default_rng(20261110)
    y1=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/5
    y2=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/7
    m=np.array([113.,173.,229.])
    synthetic=grassmann_rephasing_witness(y1,y2,m,cg,n)
    # Source benchmark W is an independently verified 3x3 matrix; it need
    # not be inverted into arbitrarily nonunique Yukawas.  Its phase law is
    # directly checked, and the full native->master Grassmann witness is
    # already recomputed by dictionary.verify_evidence above.
    a=np.array([0.9,-0.7,1.9]);s1,s2=0.65,-1.12
    actual_phase=expected_pl_phase(w,a,s1,s2,rho=1)
    require(np.all(np.isfinite(actual_phase)) and
            np.linalg.norm(w)>0,'Evolved A(-4) benchmark W invalid')
    status={
      'format':FORMAT,'model':'A_-4',
      'status':'DECLARED_WEYL_UV_COUPLING_DICTIONARY_VERIFIED_PHYSICAL_NATIVE_SIGN_OPEN',
      'source_revision':fresh['source_revision'],
      'mu_F_gev':fresh['mu_F_gev'],
      'declared_gamma_cc_choice':{
          'reference':'C0 = +i gamma^2 gamma^0 in Weyl representation; epsilon_12=+1',
          'working_rho':1,
          'working_rho_is_definition_not_installed_Matchete_measurement':True,
          'rho_minus_one_equivalent_coupling_redefinition':'GammaCC -> -GammaCC together with y1_native -> -y1_native preserves the mapped UV action and its EFT contact',
          'native_GammaCC_to_C0_measured':None,
      },
      'uv_fields_native_in_terms_of_master':{
          'L':'L_native[p] = exp(i alpha[p]) L_master[p], same phase for both SU2 components',
          'F':'F_native[r] = exp(i beta[r]) F_master[r], same phase for left/right, positive diagonal MF unchanged',
          'S1':'S1_native = exp(i sigma1) S1_master',
          'S2':'S2_native = exp(i sigma2) N^dagger S2_master',
          'N_native_to_master_triplet':master.encode_complex(n),
          'native_to_master_CG1':'G1_master[a,h]=G1_native[a,h], with rho absorbed in y1_master',
          'native_to_master_CG2':'G2_master[b,h,k]=sum_c G2_native[b,h,c] N[k,c] for Bar[S2]',
          'native_to_master_PL_CG':'Gmatch_master[a,b,k]=sum_c Gmatch_native[a,b,c] conj(N[k,c]) for S2',
          'chiral_axial_F_rephasings':'excluded because they change positive diagonal MF phase',
          'fundamental_Su2_axis_identification':'declared master reference; independent native phase not measured',
      },
      'uv_couplings_master_in_terms_of_native':{
          'y1':'y1_master[p,r] = rho*y1_native[p,r]*exp(-i alpha[p]-i beta[r]+i sigma1)',
          'y2':'y2_master[q,r] = y2_native[q,r]*exp(-i alpha[q]+i beta[r]-i sigma2)',
          'MF':'MF_master[r]=MF_native[r]>0',
          'PL_W':'W_master[p,q]=rho*W_native[p,q]*exp(i alpha[p]+i alpha[q]-i sigma1+i sigma2)',
          'source_native_PL_W':'sum_r conj(y1_native[p,r])*conj(y2_native[q,r])/MF[r]',
          'source_native_PL_external_fields':'Transp[L[q]] GammaCC PL L[p] * conj(S1) * S2',
      },
      'operator_normalization':{
          'definition':'L_EFT1^(5) = sum_{i,j=1..45; a,b=1..12} C_declared[i,j,a,b] (ell_i^T epsilon ell_j) phi_a phi_b + h.c.',
          'independent_entries':'Only i<=j and a<=b stored; full-sum coefficient = monomial/orbit_size',
          'important':'Do not infer generic-MV 1/4, Weinberg 1/2 or physical coupling relation from this declared full sum',
          'dimension':'GeV^-1',
      },
      'master_full_shape':[45,45,12,12],
      'independent_coordinate_count':36,
      'declared_rho_plus_one_independent_coordinates':encode_rephased_coordinate(rows,rho=1,expected_count=36),
      'rho_minus_one_fixed_native_y1_alternative_coordinates':encode_rephased_coordinate(rows,rho=-1,expected_count=36),
      'numerical_checks':{
          'synthetic_full_three_family_UV_to_EFT_Grassmann_rephasings':synthetic,
          'source_evolved_W_nonzero':True,
          'source_evolved_W_rephased_frobenius_norm_gev_inv':float(np.linalg.norm(actual_phase)),
          'native_action_inventory_checked':True,
          'UV_conjugate_scalar_CG_vs_PL_master_CG_max_abs_residual':cg_res,
          'prior_actual_benchmark_full_Grassmann_revalidated':True,
      },
      'source_files_sha256':{str(p):sha256(abspath(p).read_bytes()).hexdigest()
                             for p in (native,audit,connected_file,master_file,weyl_file)},
      'native_matchete_to_declared_Weyl_factor':None,
      'physical_C_ijab_at_MF':None,
      'physical_certification':False,
      'EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars simultaneously -> SMEFT',
      'no_ODE_or_scalar_threshold_implemented':True,
      'remaining_external_input':'A verified native Matchete -> chosen C0 matrix/vertex convention, or an explicit physics-level adoption of the full UV coupling dictionary; abstract native identities cannot discriminate rho.'
    }
    return status


def main(argv:list[str]|None=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path.cwd())
    parser.add_argument('--native',type=Path,default=master.DEFAULT_NATIVE)
    parser.add_argument('--audit',type=Path,default=master.DEFAULT_AUDIT)
    parser.add_argument('--connected',type=Path,default=master.DEFAULT_CONNECTED)
    parser.add_argument('--master',type=Path,default=master.DEFAULT_OUT)
    parser.add_argument('--weyl',type=Path,default=DEFAULT_WITNESS)
    parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    a=parser.parse_args(argv)
    try:
        result=verify_sources(a.root,a.native,a.audit,a.connected,a.master,a.weyl)
    except (ConventionBlocked, dictionary.DictionaryBlocked, master.BasisBlocked,
            connected.GateBlocked, boundary.BoundaryBlocked,ValueError,KeyError,
            OSError,TypeError) as exc:
        print('A(-4) convention dictionary BLOCKED:',exc)
        return 2
    out=a.out if a.out.is_absolute() else a.root.resolve()/a.out
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Declared C0 Weyl / native UV y1,y2,MF coupling map: PASS')
    print('UV-to-master SU2, complex flavour/Dirac phases and Grassmann operator: PASS')
    print('Declared working branch: rho=+1 by definition, not measured by native Matchete')
    print('Independent declared coordinates:',result['independent_coordinate_count'])
    print('Physical C_ijab(M_F): NULL; not certified')
    print('Report:',out)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
