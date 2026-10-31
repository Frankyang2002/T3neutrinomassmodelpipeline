"""Source-native B/C complex-heavy-flavour RGE covariance stress.

Numerically compares RGBeta beta functions under a complex fixed unitary
heavy-flavour field redefinition, using existing B/C UV beta exports. A
failure is actionable evidence of incompatible matrix conventions, not a
proof of a wrong physical UV theory. No source exports are modified.
"""
from __future__ import annotations
import math
import numpy as np
from .T3BCMegaPhysicalBridge import require,relative,pack_matrix

MATRIX_NAMES=('MF','y1','y2')


def unitary(seed:int=174)->np.ndarray:
    rng=np.random.default_rng(seed)
    A=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    Q,R=np.linalg.qr(A)
    Q=Q@np.diag(np.conj(np.diag(R))/np.maximum(abs(np.diag(R)),1e-30))
    require(relative(Q.conj().T@Q,np.eye(3))<1e-12,'unitary test matrix failed')
    return Q


def block_complex(v:np.ndarray, layout, name:str)->np.ndarray:
    b=layout.block(name)
    require(b.kind=='complex_matrix' and b.shape==(3,3) and b.stop-b.start==18,
            'unexpected source block '+name)
    return (v[b.start:b.start+9].reshape(3,3)+
           1j*v[b.start+9:b.stop].reshape(3,3))


def replace_complex(v:np.ndarray,layout,name:str,val:np.ndarray)->None:
    b=layout.block(name)
    z=np.asarray(val,complex)
    require(z.shape==(3,3) and np.isfinite(z).all(),'bad replacement matrix '+name)
    v[b.start:b.start+9]=z.real.ravel()
    v[b.start+9:b.stop]=z.imag.ravel()


def change_basis(vector:np.ndarray,layout,U:np.ndarray,*,y1_conjugated:bool=True)->np.ndarray:
    v=np.asarray(vector,float).copy()
    old={name:block_complex(v,layout,name) for name in MATRIX_NAMES}
    replace_complex(v,layout,'MF',U.T@old['MF']@U)
    replace_complex(v,layout,'y1',old['y1']@(U.conj() if y1_conjugated else U))
    replace_complex(v,layout,'y2',old['y2']@U)
    require(np.all(np.isfinite(v)),'nonfinite transformed RGBeta state')
    return v


def residuals_between(reference:np.ndarray,transformed:np.ndarray,layout)->dict:
    rows={}
    for name in MATRIX_NAMES:
        rows[name]=relative(block_complex(transformed,layout,name),
                            block_complex(reference,layout,name))
    names=('gY','g2','g3','lambdaH','mS1Sq','mS2Sq')
    for name in names:
        if name in layout.names:
            b=layout.block(name)
            denom=max(float(np.linalg.norm(reference[b.start:b.stop])),1e-25)
            rows[name]=float(np.linalg.norm(transformed[b.start:b.stop]-reference[b.start:b.stop])/denom)
    return rows


def run_complex_rge_covariance(model:str,*,seed:int=174,integrate:bool=True,
                               points:int=5,rtol:float=5e-8)->dict:
    from scripts import RunT3BCCoupledUVToEFT1 as bc
    from scripts import RunT3BCFullScalarMassEFT1 as ef
    from scripts import RunGeneralT3QuadraticEFT1 as legacy
    from scripts import VerifyT3BCFullScalarMassClosure as closure
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    require(model in ('B_-1','C_-1'),'only B/C(-1)')
    reg=load_registry(legacy.REGISTRY)
    state,y0,layout,_=bc.load_uv_start(model,stage29_dir=bc.STAGE29_DIR,registry=reg)
    mass_scale=min(state.baseline.mS1Sq,state.baseline.mS2Sq)
    mixed0=0.01*mass_scale*complex(math.cos(0.6),math.sin(0.6))
    v=np.concatenate((y0,[mixed0.real,mixed0.imag,-7800.]))
    tag=model.replace('-','m')
    dp=bc.DIMENSIONLESS_DIR/f'{tag}_uv_quadratic_live_rgbeta.json'
    mp=bc.MASS_DIR/f'{tag}_uv_full_mass_v1.json'
    dim=bc.checked_uv_dimensionless(model,registry=reg,layout=layout,path=dp)
    mass=closure.check_payload(closure.read_payload(mp),model,'UV')
    betas=ef.merge_beta_reports(dim,mass,layout.names)
    U=unitary(seed)
    initial_mf=block_complex(v,layout,'MF')
    rotated_mf=block_complex(change_basis(v,layout,U),layout,'MF')
    initial_heavy=np.sort(np.linalg.svd(initial_mf,compute_uv=False))
    rotated_heavy=np.sort(np.linalg.svd(rotated_mf,compute_uv=False))
    spectrum_error=relative(initial_heavy,rotated_heavy)
    require(spectrum_error<1e-12,'heavy spectrum changed under basis transformation')
    t0=math.log(1e7)
    def rhs(t,vec):
        return bc.uv_derivative(vec,mu_gev=float(math.exp(t)),layout=layout,
                                registry=reg,betas=betas)
    tested={}
    original=rhs(t0,v)
    for label,conjugated in [('source_CConj_orientation',True),('same_U_negative_control',False)]:
        rotated=change_basis(v,layout,U,y1_conjugated=conjugated)
        predicted=change_basis(original,layout,U,y1_conjugated=conjugated)
        evaluated=rhs(t0,rotated)
        residual=residuals_between(predicted,evaluated,layout)
        tested[label]={'beta_relative_residual_by_parameter':residual,
                       'max_selected_relative_beta_residual':max(residual.values()),
                       'phase_convention':('y1@conj(U),y2@U' if conjugated else 'y1@U,y2@U')}
    result={'model':model,'status':'NativeRGBetaComplexHeavyBasisCovarianceComputed',
            'source_UV_initial_scale_gev':1e7,
            'unitary_heavy_basis':pack_matrix(U),
            'initial_heavy_singular_masses_gev':initial_heavy.tolist(),
            'heavy_spectrum_rotation_relative_residual':spectrum_error,
            'beta_tests':tested,
            'source_redefinition_is_covariant_numerically':tested['source_CConj_orientation']['max_selected_relative_beta_residual']<1e-6,
            'beta_coupling_basis_independently_verified_by_matchete':False,
            'entire_renormalized_EFT_matching_performed':False}
    if integrate:
        # Use same native full 5-mass UV RHS, not a toy RG approximation.
        base=bc._integrate(v,rhs,1e7,1e5,points=points,rtol=rtol,atol=1e-10,max_step=0.35)
        terminal={}
        for label,conjugated in [('source_CConj_orientation',True),('same_U_negative_control',False)]:
            rotated=change_basis(v,layout,U,y1_conjugated=conjugated)
            sol=bc._integrate(rotated,rhs,1e7,1e5,points=points,rtol=rtol,atol=1e-10,max_step=0.35)
            predicted=change_basis(base.y[:,-1],layout,U,y1_conjugated=conjugated)
            terminal[label]={'endpoint_matrix_relative_residuals':residuals_between(predicted,sol.y[:,-1],layout),
                             'nfev':int(sol.nfev),'solver_success':True}
        result['UV_1e7_to_1e5_real_coupled_ODE_comparison']=terminal
        result['reference_nfev']=int(base.nfev)
    return result
