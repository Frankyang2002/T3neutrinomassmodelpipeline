"""Fast 18-real-Yukawa fit of the already-defined A(-4) rough EFT benchmark.

This deliberately holds the Batch-1 renormalisable background and heavy masses
fixed. The 18 UV real Yukawas are sent to muF via elementwise *baseline*
UV->muF transfer ratios. A full UV RGBeta integration per fit point is NOT
performed; call this a fixed-background approximate fit, not a full RG fit.

At each candidate parameter point:
    UV y1,y2 -> approximate F y1,y2 -> source-derived finite ordered direct
    F-hard + S-soft remainder -> symmetric C5(muS), with factor=1 assumption
    -> EXACT linear SM Weinberg one-loop transport on fixed SM background
    -> Takagi/PMNS -> five Gaussian NuFIT 6.1 observables.

The SMEFT transport matrix is constructed once with six calls to the existing
SMWeinbergEvolution solver on the basis of symmetric complex 3x3 matrices.
Because beta_K is complex-linear in K and SM running does not depend on K,
this is equivalent to rerunning that same solver for each candidate at fixed
SM background, without repeated ODE cost.

The unavailable full counterterms and other dim-5 hard/soft contributions are
not invented or falsely included. No separately computed leading log is added.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import copy
import math

import numpy as np
from scipy.optimize import least_squares

from Numerical.parallel_t3.background.general_batch1_api import decode_native, encode_native
from Numerical.parallel_t3.complete.a4_rough_benchmark import (
    estimate_boundary, run_smeft_observables,
)

SCALE = 1.e-4
PAIRS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2',
       'sin2_theta12','sin2_theta13','sin2_theta23')


class RoughFitError(ValueError):
    """Invalid or unphysical approximate fixed-background fit input."""


def matrix(z: Any) -> np.ndarray:
    x = np.asarray(decode_native(z), dtype=np.complex128)
    if x.shape != (3, 3) or not np.isfinite(x).all():
        raise RoughFitError('Expected a finite complex 3x3 matrix')
    return x


def _real_seed(values: Any, *, scale: float = SCALE) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.shape != (3, 3) or not np.isfinite(x).all() or np.any(x == 0):
        raise RoughFitError('Nonzero finite real 3x3 UV Yukawa seed required')
    return x*scale


def _scalar(z: Any, name: str) -> float:
    c = complex(np.asarray(decode_native(z)))
    if not np.isfinite(c) or abs(c.imag) > 1e-10:
        raise RoughFitError(f'{name} must be real')
    return float(c.real)


def _validate_target(target: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    if target.get('ordering') != 'NO':
        raise RoughFitError('Only the NuFIT NO five-observable fit is implemented')
    try:
        center = np.array([target['central_values'][k] for k in OBS],float)
        sigma = np.array([target['one_sigma_errors'][k] for k in OBS],float)
    except (KeyError, TypeError, ValueError) as exc:
        raise RoughFitError('Missing numerical NuFIT target') from exc
    if not np.isfinite(center).all() or not np.isfinite(sigma).all() or np.any(sigma <= 0):
        raise RoughFitError('Invalid NuFIT target or uncertainties')
    return center,sigma


@dataclass(frozen=True)
class UVToFLinearized:
    """F Yukawas at the saved benchmark and its elementwise UV transfer.

    This is *not* a general RG flow map; it is an explicitly local, simple
    rough-fit approximation to the original UV -> F running.
    """
    original_uv_scaled: np.ndarray  # (18,), values in units of 1e-4
    saved_f_y1: np.ndarray
    saved_f_y2: np.ndarray
    transfer_y1: np.ndarray
    transfer_y2: np.ndarray

    @classmethod
    def from_sources(cls, f: Mapping[str,Any], cfg: Mapping[str,Any]):
        if f.get('model') != 'A_-4' or float(f.get('scale_GeV',0)) != 1e5:
            raise RoughFitError('Wrong F source/threshold')
        c = cfg.get('models', {}).get('A_-4', {})
        uv1 = _real_seed(c.get('y1_UV_times_1e4'))
        uv2 = _real_seed(c.get('y2_UV_times_1e4'))
        values = f['native_uv_rgbeta_couplings_at_muF']
        fy1,fy2 = matrix(values['y1']),matrix(values['y2'])
        t1,t2 = fy1/uv1,fy2/uv2
        if not np.isfinite(t1).all() or not np.isfinite(t2).all():
            raise RoughFitError('Invalid UV-to-F reference transfer')
        return cls(np.concatenate((uv1.reshape(-1),uv2.reshape(-1)))/SCALE,
                   fy1,fy2,t1,t2)

    def f_yukawas(self, x: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
        """18 real fitted UV entries x, reported in units of 1e-4."""
        x=np.asarray(x,dtype=float)
        if x.shape != (18,) or not np.isfinite(x).all():
            raise RoughFitError('Require exactly 18 finite real UV parameters')
        return self.transfer_y1*(x[:9].reshape(3,3)*SCALE),\
               self.transfer_y2*(x[9:].reshape(3,3)*SCALE)

    def modified_f_input(self, f: Mapping[str,Any], x: np.ndarray) -> dict:
        out=copy.deepcopy(f)
        y1,y2=self.f_yukawas(x)
        out['native_uv_rgbeta_couplings_at_muF']['y1']=encode_native(y1)
        out['native_uv_rgbeta_couplings_at_muF']['y2']=encode_native(y2)
        return out


@dataclass(frozen=True)
class LinearSMEFT:
    images: np.ndarray     # shape (6,3,3): image of unit symmetric K at low
    ye_low: np.ndarray
    mu_low: float
    solver_runs: int

    @classmethod
    def from_scalar_input(cls, s: Mapping[str,Any], *, mu_low:float = 100.):
        if s.get('model')!='A_-4' or float(s.get('scale_GeV',0))!=1.e3:
            raise RoughFitError('Require A(-4) scalar threshold at 1 TeV')
        if s.get('both_heavy_scalar_multiplets_to_integrate_out_together') is not True:
            raise RoughFitError('Must integrate both scalars together')
        if not (0<mu_low<1.e3):
            raise RoughFitError('Require mu_low below the 1 TeV scalar threshold')
        values=s['native_eft1_rgbeta_couplings_at_muS']
        from Numerical.running.SMWeinbergEvolution import (
            SMWeinbergInitialConditions,evolve_sm_weinberg,
        )
        kwargs=dict(gY=_scalar(values['gY'],'gY'),g2=_scalar(values['g2'],'g2'),
                    g3=_scalar(values['g3'],'g3'),lambdaH=_scalar(values['lambdaH'],'lambdaH'),
                    ye=matrix(values['ye']),yu=matrix(values['yu']),yd=matrix(values['yd']))
        # K at the model benchmark is O(1e-15) GeV^-1. Integrate unit *physical*
        # basis of 1e-15 GeV^-1 so the ODE absolute tolerance resolves it.
        basis_mag=1e-15
        images=[]
        low_ye=None
        for p,q in PAIRS:
            k=np.zeros((3,3),complex)
            k[p,q]=basis_mag
            k[q,p]=basis_mag
            state=SMWeinbergInitialConditions(**kwargs,K=k)
            low=evolve_sm_weinberg(state,1e3,float(mu_low),rtol=1e-8,atol=1e-24)
            if not low.solver_success or not np.isfinite(low.K).all():
                raise RoughFitError('SMEFT basis transport failed')
            images.append(np.asarray(low.K,complex)/basis_mag)
            if low_ye is None:
                low_ye=np.asarray(low.ye,complex)
            elif not np.allclose(low.ye,low_ye,rtol=1e-11,atol=1e-16):
                raise RoughFitError('SM Yukawa evolution unexpectedly depends on K')
        return cls(np.asarray(images,complex),low_ye,float(mu_low),len(PAIRS))

    def low(self, c5_s: np.ndarray) -> np.ndarray:
        k=np.asarray(c5_s,complex)
        if k.shape!=(3,3) or not np.isfinite(k).all() or not np.allclose(k,k.T,rtol=1e-12,atol=1e-24):
            raise RoughFitError('Source K must be finite symmetric 3x3')
        result=np.einsum('k,kij->ij',np.array([k[p,q] for p,q in PAIRS]),self.images,
                         optimize=True)
        return .5*(result+result.T)


@dataclass
class FitWorkspace:
    f: dict
    s: dict
    cfg: dict
    target: dict
    uvmap: UVToFLinearized
    smeft: LinearSMEFT
    center: np.ndarray
    sigma: np.ndarray
    gamma_factor: float
    evaluations: int=0

    @classmethod
    def from_sources(cls,f:dict,s:dict,cfg:dict,target:dict,*,
                     gamma_factor:float=1.,mu_low:float=100.):
        if cfg.get('schema')!='t3-thesis-benchmark-two-model-pilot-v1' or\
           cfg.get('scales_GeV')!={'UV':1e7,'F':1e5,'S':1e3,'low':100}:
            raise RoughFitError('Unexpected thesis scales or configuration')
        if (isinstance(gamma_factor,bool) or not math.isfinite(gamma_factor)
                or gamma_factor==0):
            raise RoughFitError('GammaCC/Weyl factor must be nonzero finite real')
        center,sigma=_validate_target(target)
        uv=UVToFLinearized.from_sources(f,cfg)
        sm=LinearSMEFT.from_scalar_input(s,mu_low=mu_low)
        return cls(f,s,cfg,target,uv,sm,center,sigma,float(gamma_factor))

    def boundary(self,x:np.ndarray):
        return estimate_boundary(self.uvmap.modified_f_input(self.f,x),self.s,
                                 operator_factor=self.gamma_factor)

    def predict(self,x:np.ndarray) -> tuple[dict,dict]:
        from Numerical.parallel_t3.threshold.rough_numerics import neutrino_observables
        b=self.boundary(x)
        k_low=self.smeft.low(b.c5_at_ms)
        m_low=-.5*float(self.cfg['vev_GeV'])**2*k_low
        obs,_=neutrino_observables(m_low,self.smeft.ye_low)
        self.evaluations+=1
        return obs,{'C5_MS':b.c5_at_ms,'C5_low':k_low,'mnu_low':m_low}

    def residual5(self,obs:Mapping[str,Any]) -> np.ndarray:
        if obs.get('mixing_defined') is not True:
            raise RoughFitError('Takagi/PMNS undefined for fit candidate')
        try:
            pred=np.array([obs['dm21_sq_eV2'],obs['dm31_sq_eV2'],
                           obs['sin2_theta12'],obs['sin2_theta13'],obs['sin2_theta23']],float)
        except (KeyError,ValueError,TypeError) as exc:
            raise RoughFitError('Invalid NuFIT observables') from exc
        if not np.isfinite(pred).all():
            raise RoughFitError('Nonfinite neutrino observables')
        return (pred-self.center)/self.sigma


def fit_18_real(ws:FitWorkspace, *, regularization:float=.02,
                max_nfev:int=180, bounds:float=10.) -> tuple[dict,dict]:
    """Levenberg/trust-region fitting with a mild 18-coordinate seed prior.

    Exact Gaussian five-observable chi2 is reported separately from seed
    regularisation. No expectation of unique 18-dimensional solution.
    """
    if not (0 <= regularization < 1 and 10 <= max_nfev <= 20000 and bounds > 0):
        raise RoughFitError('Invalid fit options')
    base=ws.uvmap.original_uv_scaled.copy()
    obs0,arr0=ws.predict(base)
    dm31=float(obs0['dm31_sq_eV2'])
    if not math.isfinite(dm31) or dm31<=0 or ws.center[1]<=0:
        raise RoughFitError('Nonpositive initial atmospheric mass splitting')
    alpha=float((ws.center[1]/dm31)**.25)
    # Both Yukawas multiplied by alpha => K multiplied by alpha^2,
    # and the mass-squared differences by alpha^4.
    initial=base*alpha
    obs_init,arr_init=ws.predict(initial)
    start_chi2=float(np.dot(ws.residual5(obs_init),ws.residual5(obs_init)))
    anchor=initial.copy()
    if np.max(np.abs(anchor))>=bounds:
        raise RoughFitError('Seed exceeds requested absolute Yukawa bound')

    def residual(x):
        try:
            obs,_=ws.predict(x)
            five=ws.residual5(obs)
        except (ValueError,RuntimeError,np.linalg.LinAlgError,FloatingPointError):
            five=np.full(5,1.e5)
        penalty=regularization*(x-anchor)
        return np.concatenate((five,penalty))

    result=least_squares(residual,anchor,method='trf',
                         bounds=(-bounds*np.ones(18),bounds*np.ones(18)),
                         xtol=1e-11,ftol=1e-11,gtol=1e-11,
                         x_scale='jac',max_nfev=max_nfev)
    fitted=np.asarray(result.x,float)
    obsf,arrf=ws.predict(fitted)
    chi2=float(np.dot(ws.residual5(obsf),ws.residual5(obsf)))
    y1f,y2f=ws.uvmap.f_yukawas(fitted)
    # Report exact SMEFT output for final fitted point, not only accelerated
    # surrogate. This also confirms that the fixed-background transfer agrees.
    b=ws.boundary(fitted)
    true_summary,true_arrays=run_smeft_observables(
        b,ws.s,ws.target,vev_gev=float(ws.cfg['vev_GeV']),
        mu_low_gev=ws.smeft.mu_low)
    direct_obs=true_summary['observables']
    if not np.allclose(true_arrays['C5_low_approx_GeVinv'],arrf['C5_low'],
                       rtol=1e-6,atol=1e-25):
        raise RoughFitError('Linear SMEFT shortcut disagrees with full downstream solver')
    real_chi2=float(np.dot(ws.residual5(direct_obs),ws.residual5(direct_obs)))
    fit={
        'format':'t3-a4-rough-18-real-yukawa-fit-v1',
        'status':'Approximate18RealYukawaFit_NotFullPhysicalMatching',
        'model':'A_-4','parameters':18,
        'integrated_heavy_fields_in_order':[['F'],['S1','S2']],
        'mu_UV_GeV':1.e7,'mu_F_GeV':1.e5,'mu_S_GeV':1.e3,
        'mu_low_GeV':ws.smeft.mu_low,
        'fit_parameter_definition':'18 real UV values, reported as y_UV times 1e4',
        'UV_to_F_map':'elementwise ratios from saved Batch1 benchmark; fixed-background approximation, NOT re-integrated UV RG',
        'renormalizable_EFT1_trajectory_recomputed':False,
        'fermion_hard_soft':'original direct-channel finite pole-paired estimate only',
        'operator_factor_assumed':ws.gamma_factor,
        'SMEFT_transport':'six source-basis one-loop SMEFT solutions; exact complex-linear transport for fixed SM running',
        'SMEFT_transport_basis_solver_runs':ws.smeft.solver_runs,
        'SMEFT_final_direct_solver_runs':1,
        'optimizer':'scipy.optimize.least_squares/trf',
        'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
        'optimizer_nfev':int(result.nfev),'total_prediction_evaluations':int(ws.evaluations),
        'regularisation_coefficient':float(regularization),
        'regularisation_anchor':'uniform initial scale matched to atmospheric Delta m31 squared',
        'initial_uniform_scale_each_Yukawa':alpha,
        'baseline_observables':obs0,
        'uniform_scale_start_observables':obs_init,
        'uniform_scale_start_chi2':start_chi2,
        'final_fit_observables':direct_obs,
        'final_fit_chi2':real_chi2,
        'fast_surrogate_final_chi2':chi2,
        'fitted_y1_UV_times_1e4':fitted[:9].reshape(3,3).tolist(),
        'fitted_y2_UV_times_1e4':fitted[9:].reshape(3,3).tolist(),
        'fitted_y1_muF':encode_native(y1f),'fitted_y2_muF':encode_native(y2f),
        'original_y1_UV_times_1e4':base[:9].reshape(3,3).tolist(),
        'original_y2_UV_times_1e4':base[9:].reshape(3,3).tolist(),
        'five_Gaussian_target_values':dict(zip(OBS,map(float,ws.center))),
        'five_Gaussian_target_1sigma':dict(zip(OBS,map(float,ws.sigma))),
        'full_physical_C5_verified':False,
        'full_Matchete_finite_counterterms_included':False,
        'full_EFT1_Wilson_RGE_included':False,
        'joint_scalar_matching_other_operators_included':False,
        'all_18_UV_Yukawa_couplings_actually_RG_reintegrated':False,
        'not_an_independent_full_physics_fit':True,
    }
    arrays={'C5_MS':arrf['C5_MS'],'C5_low':true_arrays['C5_low_approx_GeVinv'],
            'mnu_100GeV_GeV':true_arrays['mnu_low_approx_GeV'],
            'y1_muF':y1f,'y2_muF':y2f,
            'y1_UV':fitted[:9].reshape(3,3)*SCALE,
            'y2_UV':fitted[9:].reshape(3,3)*SCALE}
    return fit,arrays
