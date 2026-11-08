"""A(-4) *conditional* EFT1 -> SMEFT matching at one common heavy-scalar scale.

Uses the source-verified 45-Weyl/12-real corrected-RGBeta EFT1 background,
its full evolved six-term-candidate C(M_S), and the existing simultaneous
8-real-scalar MSbar bubble matcher. This is the declared rho=+1 branch only;
the opposite global sign is obtained by real-linearity.

Computed: C5 from already-induced LLHH plus finite O(C_LLSS*lambda_HHSS)
heavy scalar bubble, and Takagi mass singular values with a stated vev.
NOT computed/certified: signed native Matchete GammaCC map, one-loop hard F
threshold, complete one-loop EFT1-to-SMEFT matching, SMEFT running or physical
neutrino masses. No local production source files are changed.
"""
from __future__ import annotations

from hashlib import sha256
from math import exp, log, pi, sqrt
from pathlib import Path
from typing import Any, Mapping
import json
import numpy as np

DEFAULT_WILSON = Path('output/a_minus4_next_batch/conditional_full_wilson_v1')
DEFAULT_OUTPUT = Path('output/a_minus4_next_batch/conditional_joint_scalar_v1')
FORMAT = 'a-minus4-conditional-joint-scalar-smeft-v1'

class ConditionalScalarBlocked(ValueError):
    """Missing physical/numerical contract or an invalid scalar threshold."""


def require(condition: bool, why: str) -> None:
    if not condition:
        raise ConditionalScalarBlocked(why)


def sha(path: Path) -> str:
    require(path.is_file(), 'Missing input: '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    require(path.is_file(), 'Missing JSON input: '+str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data,dict), 'Expected JSON object: '+str(path))
    return data


def enc(value: np.ndarray) -> dict:
    a=np.asarray(value,complex)
    require(bool(np.isfinite(a).all()),'Nonfinite complex tensor')
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def dec(value: Mapping[str,Any], shape: tuple[int,...]) -> np.ndarray:
    require(isinstance(value,Mapping) and set(value)=={'real','imag'},
            'Complex array must have real/imag parts')
    out=np.asarray(value['real'],float)+1j*np.asarray(value['imag'],float)
    require(out.shape==shape and np.isfinite(out).all(),
            'Complex array has the wrong shape or nonfinite entries')
    return out


def maxabs(a:np.ndarray) -> float:
    return float(np.max(np.abs(a))) if np.size(a) else 0.


def relative_error(a:np.ndarray,b:np.ndarray,*,floor:float=1e-100)->float:
    require(np.shape(a)==np.shape(b) and bool(np.isfinite(a).all()) and bool(np.isfinite(b).all()),
            'Cannot compare missing/nonfinite tensor components')
    return maxabs(a-b)/max(maxabs(a),maxabs(b),floor)


def lepton_tensor(c: np.ndarray, positions: np.ndarray, *,heavy:bool) -> np.ndarray:
    """Extract full LLhh or LLSS without truncating the parent 45-Weyl C.

    positions is (3 generations, 2 SU2 components) in 0-based 45-Weyl slots.
    The result is indexed [generation_p,generation_q,weak_i,weak_j,scalar_a,scalar_b].
    """
    c=np.asarray(c,complex)
    require(c.shape==(45,45,12,12) and np.isfinite(c).all(),
            'Full complex 45x45x12x12 Wilson endpoint required')
    positions=np.asarray(positions)
    require(positions.shape==(3,2) and np.issubdtype(positions.dtype,np.integer)
            and np.all((0<=positions)&(positions<45)) and len(set(positions.flat))==6,
            'Explicit distinct 3-family lepton/weak component positions required')
    part=slice(4,12) if heavy else slice(0,4)
    count=8 if heavy else 4
    out=np.empty((3,3,2,2,count,count),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    out[p,q,i,j]=c[int(positions[p,i]),int(positions[q,j]),part,part]
    return out


def master_weinberg_template()->np.ndarray:
    """Master L5 full-sum: C5*T for L5=(1/2)C5(L.H)(L.H)+h.c.

    Exactly the template used by FullSU2SourceResponse, hence
    T[nu,nu,H0R,H0R] = 1/4, and C5 = 4*C[nu,nu,H0R,H0R].
    """
    z=np.zeros((2,4),complex)
    z[0,2],z[0,3]=1/sqrt(2),1j/sqrt(2)
    z[1,0],z[1,1]=-1/sqrt(2),-1j/sqrt(2)
    return (np.einsum('ia,jb->ijab',z,z)+
            np.einsum('ib,ja->ijab',z,z))/4


def project_master(amplitude:np.ndarray, *, strict:bool=True, tol:float=1e-7)->tuple[np.ndarray,float]:
    """Use source implementation's full SU2 tensor projection, not only pivot."""
    from RGE.running.intermediate.FullSU2SourceResponse import project_weinberg
    a=np.asarray(amplitude,complex)
    require(a.shape==(3,3,2,2,4,4) and np.isfinite(a).all(),
            'Full complex LLHH [3,3,2,2,4,4] amplitude is required')
    item=project_weinberg(a)
    c5=dec(item['projected_C5_beta_16pi2_per_NGamma'],(3,3))
    residual=float(item['weinberg_residual_relative'])
    require(item['C5_flavour_symmetric'] and
            relative_error(c5,c5.T)<max(1e-10,tol),
            'Projected C5 violates complex flavour symmetry')
    if strict:
        require(item['SU2_Ward_tensor_algebra_pass'] and residual<tol,
                f'Full LLHH tensor is not Weinberg: residual {residual:.5g}')
    return (c5+c5.T)/2,residual


def gauge_invariant_mass2(m_s1_sq: float,m_s2_sq: float,
                          generators: Mapping[str,Any]|None=None) -> np.ndarray:
    """Unbroken A(-4) S1 (2 real) + S2 (6 real), with no invented mixing."""
    vals=np.asarray([m_s1_sq,m_s2_sq],float)
    require(vals.shape==(2,) and np.isfinite(vals).all() and np.all(vals>0),
            'A(-4) requires positive finite unbroken running mS1Sq,mS2Sq at M_S')
    m=np.diag(np.r_[np.repeat(vals[0],2),np.repeat(vals[1],6)])
    if generators is not None:
        for group in ('SU2','U1'):
            for t in generators[group]['scalar']:
                h=np.asarray(t)[4:,4:]
                residual=relative_error(h@m,m@h,floor=1e-50)
                require(residual<2e-12,
                        'Unbroken scalar mass tensor violates '+group+' Ward identity')
    return m


def neutrino_mass_proxy(c5:np.ndarray,vev_gev:float=246.22)->dict:
    """At scale M_S only: m_nu=-(v^2/2)C5, Takagi singular values.

    Not SMEFT evolved to the electroweak scale; physical masses unverified.
    """
    c5=np.asarray(c5,complex)
    require(c5.shape==(3,3) and np.isfinite(c5).all(),
            'Neutrino-mass proxy requires finite complex symmetric 3x3 C5')
    require(relative_error(c5,c5.T)<2e-11,'C5 not symmetric')
    v=float(vev_gev)
    require(np.isfinite(v) and 0<v<1000,'Choose a positive finite benchmark vev in GeV')
    m=-v*v*c5/2
    e=np.sort(np.linalg.svd(m,compute_uv=False))*1e9
    return {'vev_gev':v,'mass_matrix_GeV':enc(m),
            'takagi_singular_values_eV':e.tolist(),
            'sum_mass_proxies_eV':float(np.sum(e)),
            'delta_m21_sq_eV2':float(e[1]**2-e[0]**2),
            'delta_m31_sq_eV2':float(e[2]**2-e[0]**2),
            'electroweak_scale_SMEFT_running_applied':False,
            'physical_neutrino_masses_certified':False}


def _positions(fi:Mapping[str,int])->np.ndarray:
    return np.asarray([[fi[f'L[{p},{w}]'] for w in (1,2)] for p in (1,2,3)],int)


def _species_summary(bubble:dict, *, master_bubble:np.ndarray) -> dict:
    """Group 36 *eigenstate* pairs; do not mislabel mixed eigenvectors."""
    out={k:np.zeros((3,3),complex) for k in ('S1_S1','S1_S2','S2_S2')}
    ambiguous=0
    for row in bubble['pair_data']:
        key=row['pair']
        a,b=row['eigenstate_indices']
        u=bubble['mass_rotation_U']
        weights=[]
        for k in (a,b):
            frac=float(np.sum(u[:2,k]**2))
            require(abs(frac+(float(np.sum(u[2:,k]**2)))-1)<2e-10,
                    'Heavy scalar eigenvector not normalized')
            if abs(frac-1.)<2e-10: weights.append('S1')
            elif frac<2e-10:weights.append('S2')
            else:weights.append('?')
        if '?' in weights:
            ambiguous+=1
            continue
        group='_'.join(sorted(weights))
        cm,_=project_master(bubble['pair_amplitudes'][key],strict=False)
        out[group]+=cm
    if ambiguous==0:
        require(relative_error(sum(out.values()),master_bubble,floor=1e-55)<2e-9,
                'Eigenstate pair projections do not reconstruct joint bubble')
    return {'pair_count':len(bubble['pair_data']),
            'mixed_species_eigenvector_pairs':ambiguous,
            'grouped_C5_per_eigenpair_gev_inv':{k:enc(v) for k,v in out.items()},
            'grouped_C5_max_abs_gev_inv':{k:maxabs(v) for k,v in out.items()},
            'grouped_species_complete':ambiguous==0}


def compute_from_arrays(c_ms:np.ndarray,quartic_ms:np.ndarray,
                        m_s1_sq:float,m_s2_sq:float,positions:np.ndarray,
                        mu_s:float,*,vev_gev:float=246.22,
                        generators:Mapping[str,Any]|None=None) -> dict:
    """Core pure-numerical matching; full projectors and scalar matcher reused."""
    from Numerical.joint_scalar_matching.core import (
        scalar_bubble_matching, project_weinberg as project_joint)
    mu_s=float(mu_s)
    require(np.isfinite(mu_s) and mu_s>0,'Matching scale must be positive')
    c=np.asarray(c_ms,complex)
    q=np.asarray(quartic_ms)
    require(q.shape==(12,12,12,12) and np.isfinite(q).all() and
            maxabs(np.asarray(q).imag)<1e-11,
            'Corrected real 12^4 quartic tensor required')
    q=q.real
    mass2=gauge_invariant_mass2(m_s1_sq,m_s2_sq,generators)
    d=lepton_tensor(c,positions,heavy=False)
    hh,c5_direct_res=project_master(d,strict=True)
    h=lepton_tensor(c,positions,heavy=True)
    lam=q[:4,:4,4:,4:]
    beta_amp=np.einsum('pqijab,ABab->pqijAB',h,lam,optimize=True)
    beta_c5,beta_res=project_master(beta_amp,strict=True)
    bubble=scalar_bubble_matching(h,lam,mass2,mu_s,
            group_names=['S1']*2+['S2']*6)
    require(len(bubble['pair_data'])==36,
            'All eight heavy real scalar eigenstates must be matched together')
    amp=bubble['amplitude']
    c5_bubble,bubble_res=project_master(amp,strict=True)
    c5_joint,joint_res=project_joint(amp,strict=False)
    factor2=relative_error(c5_bubble,2*c5_joint,floor=1e-55)
    require(factor2<2e-9 and abs(bubble_res-joint_res)<2e-9,
            'Joint-scalar/master Weinberg normalization factor 2 is inconsistent')
    # The static one-loop MSbar logarithm identity must hold at fixed C,lambda,mass:
    # d(delta C5)/d ln mu = - beta_C5|lambda / (16*pi^2).
    step=2e-4
    up=scalar_bubble_matching(h,lam,mass2,mu_s*exp(step))
    down=scalar_bubble_matching(h,lam,mass2,mu_s*exp(-step))
    cup,_=project_master(up['amplitude'],strict=True)
    cdn,_=project_master(down['amplitude'],strict=True)
    derivative=(cup-cdn)/(2*step)
    predicted=-beta_c5/(16*pi*pi)
    log_error=relative_error(derivative,predicted,floor=1e-55)
    require(log_error<2e-8,
            f'Scalar threshold and EFT1 quartic-beta logarithms disagree: {log_error:g}')
    c5=hh+c5_bubble
    require(relative_error(c5,c5.T,floor=1e-55)<1e-10,
            'Combined conditional Weinberg coefficient not complex symmetric')
    grouped=_species_summary(bubble,master_bubble=c5_bubble)
    mass_plus=neutrino_mass_proxy(c5,vev_gev)
    mass_minus=neutrino_mass_proxy(-c5,vev_gev)
    require(np.allclose(mass_plus['takagi_singular_values_eV'],
                        mass_minus['takagi_singular_values_eV'],
                        rtol=2e-14,atol=1e-22),
            'Global rho=-1 branch changed Takagi singular values')
    result={
        'mu_S_gev':mu_s,
        'mS1Sq_at_MS_gev2':float(m_s1_sq),
        'mS2Sq_at_MS_gev2':float(m_s2_sq),
        'heavy_scalar_masses_gev':{'S1':sqrt(m_s1_sq),'S2':sqrt(m_s2_sq)},
        'heavy_real_scalar_eigenvalues_sq_gev2':bubble['eigenvalues_sq_gev2'].tolist(),
        'joint_eigenstate_count':8,'one_common_scalar_threshold':True,
        'mass_basis': 'unbroken gauge invariant S1(2 real) plus S2(6 real); no EWSB mixing',
        'quartic_basis':'corrected native RGBeta MSbar V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
        'direct_C5_from_EFT1_LLHH_gev_inv':enc(hh),
        'finite_scalar_bubble_C5_gev_inv':enc(c5_bubble),
        'conditional_SMEFT_C5_gev_inv':enc(c5),
        'max_abs_direct_C5_gev_inv':maxabs(hh),
        'max_abs_scalar_bubble_C5_gev_inv':maxabs(c5_bubble),
        'max_abs_conditional_C5_gev_inv':maxabs(c5),
        'projected_16pi2_beta_C5_from_quartic_gev_inv':enc(beta_c5),
        'direct_weinberg_tensor_residual':c5_direct_res,
        'quartic_beta_weinberg_tensor_residual':beta_res,
        'scalar_bubble_weinberg_tensor_residual':bubble_res,
        'joint_vs_master_projector_factor2_residual':factor2,
        'matching_scale_log_derivative_relative_residual':log_error,
        'scalar_eigenpair_decomposition':grouped,
        'neutrino_mass_proxy':mass_plus,
        'rho_minus_conditional_C5_gev_inv':enc(-c5),
        'rho_minus_mass_singular_values_equal_rho_plus':True,
        'physical_C5_certified':False,
        'physical_neutrino_masses_certified':False,
        'finite_matching_subset':'O(C_LLSS lambda_HHSS) MSbar bubble at zero external momentum',
        'UV_hard_fermion_one_loop_threshold_included':False,
        'full_EFT1_to_SMEFT_one_loop_matching_certified':False,
        'SMEFT_RG_running_to_electroweak_scale_included':False,
    }
    return result


def run(root:Path, output:Path,*, wilson_root:Path=DEFAULT_WILSON,
        vev_gev:float=246.22)->dict:
    """Read-only source checks then one joint heavy-scalar threshold."""
    from Numerical.running import Aminus4EFT1Couplings as eft
    from Numerical.running.Aminus4ConditionalFullWilson import (
        load_contract, FORMAT as WILSON_FORMAT)
    root=Path(root).resolve()
    output=Path(output)
    if not output.is_absolute():output=root/output
    wilson_root=Path(wilson_root)
    if not wilson_root.is_absolute():wilson_root=root/wilson_root
    require(not output.exists(),'Refusing to overwrite output: '+str(output))
    summary_path=wilson_root/'summary.json'
    endpoint_path=wilson_root/'full_conditional_C_trajectory.npz'
    w=read_json(summary_path)
    require(w.get('format')==WILSON_FORMAT and w.get('model')=='A_-4' and
            w.get('beta_mode')=='fonseca_eq22_six_term' and
            w.get('full_beta_loop_certified') is False and
            w.get('scalar_thresholds_integrated_simultaneously') is False and
            w.get('physical_scalar_mass_eigenbasis_and_CG_map_certified') is False and
            w.get('matched_SMEFT_C5_available') is False,
            'Input must be an unmodified declared conditional full-Wilson EFT1 result')
    ctx=load_contract(root)
    mu_s=float(ctx['scales'][0]);mu_f=float(ctx['scales'][-1])
    require(abs(mu_s-float(w['mu_S_gev']))<1e-8 and
            abs(mu_f-float(w['mu_F_gev']))<1e-8,
            'Wilson EFT1 and corrected RGBeta thresholds do not agree')
    require(ctx['production_summary_sha256']==w['input_provenance']['production_summary_sha256']
            and ctx['tensor_npz_sha256']==w['input_provenance']['tensor_npz_sha256']
            and ctx['boundary_sha256']==w['input_provenance']['boundary_sha256']
            and ctx['fixed_beta_sha256']==w['input_provenance']['fixed_beta_sha256'],
            'Wilson input provenance differs from current corrected RGBeta/background')
    with np.load(endpoint_path,allow_pickle=False) as data:
        require(set(('C_full_rho_plus_MF_gev_inv','C_full_rho_plus_MS_gev_inv',
                     'C_full_rho_minus_MS_gev_inv'))<=set(data.files),
                'Saved full Wilson NPZ lacks both sign branches')
        c0=np.asarray(data['C_full_rho_plus_MF_gev_inv'],complex)
        cm=np.asarray(data['C_full_rho_plus_MS_gev_inv'],complex)
        cn=np.asarray(data['C_full_rho_minus_MS_gev_inv'],complex)
    require(c0.shape==cm.shape==cn.shape==(45,45,12,12) and
            np.isfinite(cm).all() and
            relative_error(c0,ctx['c0'])<1e-12 and
            relative_error(cm,-cn)<1e-12,
            'Saved endpoint is not the verified full declared Wilson trajectory')
    require(relative_error(cm,cm.swapaxes(0,1))<2e-10 and
            relative_error(cm,cm.swapaxes(2,3))<2e-10,
            'Full Wilson endpoint lacks independent pair symmetries')
    expected=dec(w['induced_LLHH_projection_at_MS_not_SMEFT_C5_gev_inv'],(3,3))
    positions=_positions(eft.FI)
    direct,_=project_master(lepton_tensor(cm,positions,heavy=False),strict=True)
    require(relative_error(direct,expected,floor=1e-55)<3e-9,
            'Unchanged Wilson endpoint disagrees with stored full RGE LLHH projection')
    basis_path=root/'output/a_minus4_next_batch/fresh_index_production/running/A_m4_eft1_tensor_basis.json'
    bg=read_json(basis_path)
    require(bg.get('rge_payload_sha256')==ctx['fixed_beta_sha256'] and
            len(bg.get('running_couplings',[]))==7 and
            abs(float(bg['sampled_scales_gev'][-1])-mu_s)<1e-8,
            'Scalar running mass background does not match corrected EFT1')
    ms=eft._point(bg['running_couplings'][-1])
    quartic=ctx['samples'][0].quartic
    # This contract already checks the saved NPZ against independent tensor_point()
    # at M_S, preventing use of unrelated gauge/quartic conventions.
    result=compute_from_arrays(cm,quartic,ms['mS1Sq'],ms['mS2Sq'],positions,
                                mu_s,vev_gev=vev_gev,
                                generators=eft.gauge_generators())
    require(relative_error(dec(result['direct_C5_from_EFT1_LLHH_gev_inv'],(3,3)),expected,
                           floor=1e-55)<3e-9,
            'Matching changed previously evolved EFT1 direct LLHH')
    report={
        'format':FORMAT,'model':'A_-4',
        'status':'CONDITIONAL_EFT1_TO_SMEFT_JOINT_2_SCALAR_BUBBLE_NOT_PHYSICAL',
        'authority':'declared_rho_plus_one_GammaCC_unmeasured',
        'EFT_sequence':'UV -> integrate out ALL fermions at MF -> EFT1(SM+S1+S2) '
                       '-> integrate out ALL scalars at a COMMON MS -> SMEFT',
        'fermion_threshold_MF_gev':mu_f,
        'eight_heavy_real_scalar_components_integrated_at_MS':True,
        'complete_45_Weyl_12_real_Wilson_tensor_used':True,
        'independent_complex_Wilson_coordinates':80730,
        'source_Wilson_ODE_mode':'fonseca_eq22_six_term (candidate gauge extension)',
        'native_Matchete_GammaCC_to_Weyl_sign_and_magnitude_measured':False,
        'matchete_to_RGBeta_scalar_invariant_map_certified':False,
        'physical_neutrino_masses_certified':False,
        'input_sha256':{
            'conditional_wilson_summary':sha(summary_path),
            'conditional_wilson_full_NPZ':sha(endpoint_path),
            'corrected_quartic_basis':sha(basis_path),
            **w['input_provenance']},
        'matching':result,
    }
    output.mkdir(parents=True)
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf8')
    np.savez_compressed(output/'conditional_smeft_c5_and_mnu.npz',
        c5_direct_gev_inv=dec(result['direct_C5_from_EFT1_LLHH_gev_inv'],(3,3)),
        c5_scalar_bubble_gev_inv=dec(result['finite_scalar_bubble_C5_gev_inv'],(3,3)),
        c5_total_rho_plus_gev_inv=dec(result['conditional_SMEFT_C5_gev_inv'],(3,3)),
        c5_total_rho_minus_gev_inv=dec(result['rho_minus_conditional_C5_gev_inv'],(3,3)),
        neutrino_mass_matrix_rho_plus_gev=dec(result['neutrino_mass_proxy']['mass_matrix_GeV'],(3,3)),
        neutrino_mass_singular_values_rho_plus_ev=np.asarray(
            result['neutrino_mass_proxy']['takagi_singular_values_eV'],float))
    return report
