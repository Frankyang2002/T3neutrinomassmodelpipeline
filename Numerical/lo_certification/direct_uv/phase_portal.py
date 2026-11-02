"""Conditional all-16 direct UV phase/portal transport; unchanged exact loop core.

No guessed 4-spinor/Matchete/RGBeta physical dictionary. 'Equal up to a
field redefinition' is proved numerically; alternative coupling identifications
are explicitly *different benchmark hypotheses*, not certified results.
"""
from __future__ import annotations

from copy import deepcopy
from typing import Any
import math
import numpy as np

from Numerical.direct_one_loop.core import _complex, encoded, evaluate
from .certification import CLASSES
from .channels import partition_exact


class ConventionError(ValueError):
    pass


def _require(test: bool, explanation: str) -> None:
    if not test:
        raise ConventionError(explanation)


def _scalar_complex(data: dict[str, Any], field: str) -> np.ndarray:
    return _complex(data[field],field)


def rephase(data: dict[str, Any], fermion_phases: np.ndarray,
            scalar_phases: np.ndarray) -> dict[str, Any]:
    """Exact change of integration coordinates psi=Df psi', phi=Ds phi'.

    For L=-nu psi(Y phi + Z phi*) + h.c.;
    V=phi^dag A phi + 1/2 (phi^T B phi + h.c.).
    The transformed mass/Y/Z tensors are source-basis tensors, not a change
    of physical couplings. No charge indices are permuted or removed.
    """
    out=deepcopy(data)
    M=_scalar_complex(data,'fermion_mass_gev')
    A=_scalar_complex(data,'scalar_A_gev2')
    B=_scalar_complex(data,'scalar_B_gev2')
    pf=np.asarray(fermion_phases,complex)
    ps=np.asarray(scalar_phases,complex)
    _require(pf.shape==(len(M),) and ps.shape==(len(A),), 'Wrong field phase dimensions')
    _require(np.allclose(np.abs(pf),1.,rtol=0,atol=1e-13) and
             np.allclose(np.abs(ps),1.,rtol=0,atol=1e-13),'Field redefinitions must be unit modulus')
    out['fermion_mass_gev']=encoded((pf[:,None]*M)*pf[None,:])
    out['scalar_A_gev2']=encoded((ps.conj()[:,None]*A)*ps[None,:])
    out['scalar_B_gev2']=encoded((ps[:,None]*B)*ps[None,:])
    for original, row in zip(data['vertices'],out['vertices']):
        for name in ('field','conjugate'):
            if name in original:
                v=_complex(original[name],name)
                _require(v.shape==(3,len(M),len(A)),name+' vertex dimensions inconsistent')
                outphase=ps if name=='field' else ps.conj()
                row[name]=encoded(v*pf[None,:,None]*outphase[None,None,:])
    # Reset misleadingly signed input provenance: these are phase-gauge
    # transformations with explicit unitary witnesses, not source file edits.
    out.setdefault('basis_transform',{})['fermion_phase']=encoded(pf)
    out['basis_transform']['scalar_phase']=encoded(ps)
    return out


def portal_geometry(portal_model: dict[str, Any], kind: str,
                    *, tol:float=2e-9) -> dict[str, Any]:
    """Reconstruct identity H0 projectors from the native all-index span.

    The full 2x2xdxd tensor identity is a *source-reported* condition. We
    additionally recompute its neutral-VEV projection from native projectors
    and decomposition coefficients, without identifying coupling signs.
    """
    _require(isinstance(portal_model,dict) and len(portal_model.get('scalars',[]))==2,
             'Missing native portal identity record for both scalar multiplets')
    d1,d2,_,_=CLASSES[kind]
    sectors=[]
    for i,d in enumerate((d1,d2)):
        sector=portal_model['scalars'][i]
        _require(sector.get('scalar')==f'S{i+1}' and sector.get('dimension')==d,
                 'Wrong native scalar portal sector/order')
        _require(sector.get('full_tensor_identity_in_real_native_span') is True,
                 'Full SU(2) native portal identity not established')
        _require(sector.get('does_not_certify_coupling_sign_or_normalization') is True,
                 'Unexpected source certificate state: source audit contract changed')
        coeff=np.asarray(sector['native_coupling_coefficients_for_unit_identity'],float)
        proj=np.asarray(sector['native_vev_projector_diagonals'],float)
        _require(proj.ndim==2 and proj.shape==(len(coeff),d),
                 'Native scalar portal projectors malformed')
        _require(len(sector.get('native_coupling_names',[]))==len(coeff),
                 'Native portal coupling names and projectors mismatch')
        projection=coeff@proj
        residual=float(np.max(np.abs(projection-1.)))
        _require(np.isfinite(projection).all() and residual<tol,
                 'Native full tensor identity does not reconstruct unit H0 portal: '+str(residual))
        sectors.append({'scalar':f'S{i+1}',
                        'native_coupling_names':sector['native_coupling_names'],
                        'unit_identity_native_coefficients':coeff.tolist(),
                        'native_projected_unit_identity':projection.tolist(),
                        'maximum_neutral_projector_error':residual,
                        'full_tensor_residual_as_reported':sector.get('max_abs_four_index_residual'),
                        'physical_Matchete_RGBeta_portal_sign_certified':False})
    return {'status':'native_portal_geometry_verified_physical_sign_pending',
            'sector_witnesses':sectors, 'portal_dictionary_certified':False}


def portal_candidate(data:dict[str,Any], bench:dict[str,Any],
                     portal_model:dict[str,Any], kind:str,
                     *, sign:int) -> tuple[dict[str,Any], dict[str,Any]]:
    """Construct explicit signed RGBeta-to-Matchete *hypothesis*.

    A_s = mS_s^2 - sign * lambdaH_RGBeta * v^2/2 * I_s,
    with sign=+1 meaning numerical lambdaNative=+lambdaRGBeta; sign=-1
    meaning lambdaNative=-lambdaRGBeta. The previous provisional input uses
    sign=+1. Original inter-multiplet A and general B are retained unchanged.
    """
    _require(sign in (-1,1),'Portal identification sign must be +/-1')
    geo=portal_geometry(portal_model,kind)
    ordinary=bench['base_state']['ordinary']; v=float(bench['numerical']['vev_gev'])
    _require(v>0 and math.isfinite(v),'Invalid input Higgs vev')
    out=deepcopy(data)
    A=_scalar_complex(data,'scalar_A_gev2').copy()
    d1,d2,_,_=CLASSES[kind]
    _require(A.shape==(d1+d2,d1+d2),'Wrong scalar dimensions')
    residuals=[]
    for i,(start,d) in enumerate(((0,d1),(d1,d2))):
        l=float(ordinary[f'lambdaH{i+1}'])
        m2=float(ordinary[f'mS{i+1}Sq'])
        _require(np.isfinite(l) and np.isfinite(m2),'Nonfinite portal benchmark')
        projection=np.asarray(geo['sector_witnesses'][i]['native_projected_unit_identity'],float)
        baseline=m2-l*v*v/2*projection
        actual=A.diagonal()[start:start+d]
        residual=float(np.max(np.abs(actual-baseline))/max(1.,abs(m2)))
        _require(residual<2e-9,
                 f'S{i+1}: provisional input diagonal differs from restricted native portal candidate ({residual:.3g})')
        desired=m2-sign*l*v*v/2*projection
        for j in range(d):A[start+j,start+j]+=desired[j]-actual[j]
        residuals.append(residual)
    out['scalar_A_gev2']=encoded(A)
    out['conditional_portal_hypothesis']={
        'native_over_RGBeta_lambdaH_sign':sign,
        'Matchete_source_portal_Lagrangian':'(1/2)*PlusHc[lambdaHInv*Hdag H*Sdag S*CG]',
        'H0_vev':'v/sqrt(2)',
        'physical_cross_package_dictionary_certified':False,
        'original_topology_mixing_unchanged':True,
        'original_B_bilinear_unchanged':True,
        'native_projector_witness':geo,
        'baseline_diagonal_relative_errors':residuals}
    return out, out['conditional_portal_hypothesis']


def _relative(a:np.ndarray,b:np.ndarray,reference:float|None=None)->float:
    denom=max(float(np.linalg.norm(b)),float(np.linalg.norm(a)),
              0. if reference is None else reference,1e-75)
    return float(np.linalg.norm(a-b)/denom)


def mass_record(result:dict[str,Any], baseline:np.ndarray) -> dict[str,Any]:
    m=_complex(result['total']['complex_symmetric_matrix_gev'],'mass',(3,3))
    d=result['diagnostics']
    return {'mass_matrix_gev':encoded(m),
            'takagi_masses_ev':result['total']['neutrino_masses_ev'],
            'mass_matrix_difference_frobenius_relative_to_baseline':_relative(m,baseline),
            'uv_pole_relative':d['relative_uv_pole_coefficient'],
            'cancellation_ratio':d['cancellation_ratio_total_over_sum'],
            'cancellation_factor':d['cancellation_factor_sum_over_total'],
            'vertex_pair_matrices_gev':result['vertex_pair_matrices_gev'],
            'physical_source_certified':False}


def _eval(data:dict[str,Any]) -> dict[str,Any]:
    r=evaluate(data,retain_individual_terms=False)
    _require(r['diagnostics']['uv_finite'] and r['diagnostics']['neutrino_mass_matrix_symmetric'],
             'Conditional loop is not UV finite or symmetric')
    return r


def _complex_action_witness(original:dict[str,Any], transformed:dict[str,Any],
                            pf:np.ndarray,ps:np.ndarray) -> dict[str,Any]:
    """Independent complex-tensor Lagrangian identity, no loop diagonalisation.

    Arbitrarily complex phases may amplify floating cancellation in Takagi
    eigenvectors. This algebraic test does not silently call such cancellation
    failures physical phase violations.
    """
    M=_scalar_complex(original,'fermion_mass_gev')
    A=_scalar_complex(original,'scalar_A_gev2')
    B=_scalar_complex(original,'scalar_B_gev2')
    Mp=_scalar_complex(transformed,'fermion_mass_gev')
    Ap=_scalar_complex(transformed,'scalar_A_gev2')
    Bp=_scalar_complex(transformed,'scalar_B_gev2')
    nf,ns=len(M),len(A)
    psi=(np.arange(nf)+1.)/(nf+1.)+1j*np.cos(np.arange(nf)*.37)
    phi=np.sin(np.arange(ns)+.13)+1j*np.cos(np.arange(ns)*.53)
    nu=np.array([.23+.14j,-.37+.06j,.19-.31j])
    oldpsi=pf*psi;oldphi=ps*phi
    oldmass=oldpsi.T@M@oldpsi;newmass=psi.T@Mp@psi
    oldpotential=np.vdot(oldphi,A@oldphi)+np.real(oldphi.T@B@oldphi)
    newpotential=np.vdot(phi,Ap@phi)+np.real(phi.T@Bp@phi)
    def yukawa(obj,psiv,phiv):
        value=0j
        for row in obj['vertices']:
            y=_complex(row['field'],'Y') if 'field' in row else np.zeros((3,nf,ns),complex)
            z=_complex(row['conjugate'],'Z') if 'conjugate' in row else np.zeros((3,nf,ns),complex)
            value+=np.einsum('a,i,ais,s->',nu,psiv,y,phiv,optimize=True)
            value+=np.einsum('a,i,ais,s->',nu,psiv,z,phiv.conj(),optimize=True)
        return value
    oldy=yukawa(original,oldpsi,oldphi)
    newy=yukawa(transformed,psi,phi)
    residuals={
        'fermion_mass_bilinear':abs(oldmass-newmass)/max(abs(oldmass),abs(newmass),1.),
        'scalar_quadratic_potential':abs(oldpotential-newpotential)/max(abs(oldpotential),abs(newpotential),1.),
        'full_Y_and_Z_complex_vertex_action':abs(oldy-newy)/max(abs(oldy),abs(newy),1.),
    }
    residuals={k:float(v) for k,v in residuals.items()}
    _require(max(residuals.values())<1e-12,'Complex phase transport violates action identity')
    return {'status':'complex_action_identity_verified_not_loop_certificate',
            'residuals':residuals,'physical_source_phase_certified':False}


def inspect_phase_orbit(data:dict[str,Any], baseline:dict[str,Any],
                        *, tol:float=4e-5) -> dict[str,Any]:
    """Exact complex action transport, plus full-loop real sign orbit.

    General complex rephasings are exact algebraically. For M_F >> m_S,
    large cancelling B0 sums combined with arbitrary complex Takagi phase
    choices may hit the unmodified engine's strict symmetry tolerance.
    Therefore the complex witness is an action identity, NOT claimed to be a
    numerical one-loop check; the latter is separately verified for signs.
    """
    M=_scalar_complex(data,'fermion_mass_gev')
    A=_scalar_complex(data,'scalar_A_gev2')
    mm=_complex(baseline['total']['complex_symmetric_matrix_gev'],'baseline',(3,3))
    vals=np.asarray(baseline['total']['neutrino_masses_ev'],float)
    pf=np.exp(1j*(.31+np.arange(len(M))*.173))
    ps=np.exp(1j*(-.67+np.arange(len(A))*.293))
    complex_action=_complex_action_witness(data,rephase(data,pf,ps),pf,ps)
    cases=[]
    for label,pf,ps in [
        ('alternating_real_sign',(-1.)**np.arange(len(M)),(-1.)**np.arange(len(A))),
        ('half_scalar_real_signs',np.ones(len(M)),np.r_[np.ones(len(A)//2),-np.ones(len(A)-len(A)//2)])
    ]:
        # The second sign pattern is just an independent exact signed basis
        # change. No assumption about either scalar multiplet is used.
        transformed=rephase(data,pf,ps)
        test=_eval(transformed)
        test_m=_complex(test['total']['complex_symmetric_matrix_gev'],'test',(3,3))
        abs_error=float(np.linalg.norm(mm-test_m))
        rel=_relative(test_m,mm)
        singular=np.asarray(test['total']['neutrino_masses_ev'],float)
        singular_rel=float(np.linalg.norm(vals-singular)/max(np.linalg.norm(vals),1e-80))
        _require(rel<tol and singular_rel<tol,
                 f'Exact signed field-coordinate change altered mnu: {label}: {rel:.3g}/{singular_rel:.3g}')
        cases.append({'probe':label,'fermion_phases':encoded(pf),
                      'scalar_phases':encoded(ps),
                      'mass_matrix_absolute_error_gev':abs_error,
                      'mass_matrix_relative_error':rel,
                      'takagi_mass_relative_error':singular_rel,
                      'passed':True})
    return {'status':'complex_action_and_two_exact_signed_loop_orbits_verified',
            'complex_phase_action_witness':complex_action,
            'full_loop_sign_rephasing_cases':cases, 'tolerance':tol,
            'general_complex_phase_full_loop_numeric_test':'not_asserted_due_to_cancellation_conditioning',
            'physical_source_vertex_normalisation_certified':False}

def convention_scan(data:dict[str,Any],kind:str,alpha:int,
                    bench:dict[str,Any],portal_model:dict[str,Any],
                    baseline:dict[str,Any]|None=None,
                    *,suppressed_fraction:float=.01) -> dict[str,Any]:
    """Full exact one-loop mass response to unresolved physical dictionaries.

    We do not set *any* absolute phase value or claim a physical result.
    Positive/negative native portal identifications are separate hypotheses.
    The simultaneous S2 -> -S2 mixing/Y2 change is a field redefinition,
    whereas flipping only the mixing or only y2 is a different Lagrangian.
    """
    base=baseline if baseline is not None else _eval(data)
    m0=_complex(base['total']['complex_symmetric_matrix_gev'],'baseline',(3,3))
    orbit=inspect_phase_orbit(data,base)
    plus,proof=portal_candidate(data,bench,portal_model,kind,sign=+1)
    minus,proof_minus=portal_candidate(data,bench,portal_model,kind,sign=-1)
    _require(_relative(_scalar_complex(plus,'scalar_A_gev2'),_scalar_complex(data,'scalar_A_gev2'))<1e-13,
             'Native positive portal hypothesis should reproduce baseline input')
    _require(proof['native_projector_witness']==proof_minus['native_projector_witness'],
             'Changing portal sign altered SU2 geometry')
    # S2 -> -S2 rotates exactly the complete kinetic/M/A/B and *all*
    # independent vertices. This is strictly an unphysical coordinate change.
    d1=CLASSES[kind][0]; ns=_scalar_complex(data,'scalar_A_gev2').shape[0]
    nf=_scalar_complex(data,'fermion_mass_gev').shape[0]
    s2ph=np.ones(ns,complex);s2ph[d1:]=-1.
    redefined=rephase(data,np.ones(nf,complex),s2ph)
    alternative=_eval(redefined)
    altmat=_complex(alternative['total']['complex_symmetric_matrix_gev'],'rephased',(3,3))
    _require(_relative(altmat,m0)<4e-5,'S2 sign-rephasing changed neutrino mass')

    # Counterfactual: flip scalar-mixing sign WITHOUT compensating phase
    # transport of any other term; valid only when no other S1/S2 entries
    # mix and B=0 in the restricted frozen point.
    mixed=deepcopy(data)
    amat=_scalar_complex(data,'scalar_A_gev2').copy()
    amat[:d1,d1:]*=-1;amat[d1:,:d1]*=-1
    mixed['scalar_A_gev2']=encoded(amat)
    topology_result=_eval(mixed)
    # Independent vertex-side hypothesis: multiply only the source-y2
    # scalar vertex by -1. In this restricted benchmark the two hypotheses
    # are exactly related by S2 -> -S2. Their equality is a *basis theorem*,
    # not evidence for choosing the physical source sign.
    y2_flip=deepcopy(data)
    y2_rows=[r for r in y2_flip['vertices'] if r['name']=='y2_HC']
    _require(len(y2_rows)==1,'Expected exactly one restricted y2 source vertex')
    for field in ('field','conjugate'):
        if field in y2_rows[0]:
            y2_rows[0][field]=encoded(-_complex(y2_rows[0][field],field))
    y2_flip_result=_eval(y2_flip)
    topology_m=_complex(topology_result['total']['complex_symmetric_matrix_gev'],'mixflip',(3,3))
    y2_m=_complex(y2_flip_result['total']['complex_symmetric_matrix_gev'],'y2flip',(3,3))
    _require(_relative(topology_m,y2_m)<4e-5,
             'T3 mixing-sign and y2-sign hypotheses not related by S2 rephasing')
    alt_portal_result=_eval(minus)
    cases={
        'baseline_native_plus':mass_record(base,m0),
        'conditional_native_minus_RGBeta_portals':mass_record(alt_portal_result,m0),
        'conditional_opposite_T3_mixing_sign':mass_record(topology_result,m0),
        'conditional_opposite_y2_vertex_sign':mass_record(y2_flip_result,m0),
        'S2_rephasing_with_all_vertices_transported':mass_record(alternative,m0),
    }
    ch=partition_exact(minus,kind,alpha,alt_portal_result,weak_fraction=suppressed_fraction)
    return {'status':'conditional_exact_mass_variants_not_physical_certification',
            'model':data.get('model'),
            'source_portal_geometry':proof['native_projector_witness'],
            'phase_orbit':orbit,
            'hypotheses':{
                'portal_plus_is_previous_provisional_adapter':True,
                'portal_minus_is_RGBeta_to_Matchete_sign_candidate':True,
                'topology_mixing_sign_physically_determined':False,
                'yukawa_source_to_Weyl_phase_physically_determined':False,
                'extra_source_allowed_couplings_fixed_by_input':False,
            },
            'variants':cases,
            'conditional_minus_portal_charge_and_suppression':ch,
            'relative_shift_from_portal_hypothesis':cases['conditional_native_minus_RGBeta_portals']['mass_matrix_difference_frobenius_relative_to_baseline'],
            'relative_shift_from_topology_mixing_hypothesis':cases['conditional_opposite_T3_mixing_sign']['mass_matrix_difference_frobenius_relative_to_baseline'],
            'relative_matrix_discrepancy_mix_flip_vs_y2_flip':_relative(topology_m,y2_m),
            'physical_certified':False}
