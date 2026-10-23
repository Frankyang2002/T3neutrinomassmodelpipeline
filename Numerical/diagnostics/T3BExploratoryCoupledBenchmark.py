"""Exploratory full-flavour T3-B C12/C5 numerical benchmark.

The off-diagonal coefficients and scalar diagonal coefficient come from the
audited one-generation tensor. The C12 lepton-flavour diagonal action follows
from the generic master RGE collinear fermion term and is checked by a
separate two-generation complex-Yukawa regression. The tested operator sector does not include every generated component.
Heavy-mass logarithms are implemented as a common-scale Wilson matching
boundary, not as sequential partially-decoupled renormalisable EFTs.
Production C5 is never overwritten.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from Numerical.diagnostics.T3BHeavyLogAlignment import flavor_kernel_by_heavy, matrix_from_json, matrix_to_json
from RGE.running.intermediate.T3BC12FlavorRGE import c12_beta_numerator

LOOP = 16.0 * np.pi**2
FWD = 4.0 * np.sqrt(3.0) / 3.0
REV = 8.0 * np.sqrt(3.0) / 3.0


def _complex_scalar(value):
    if isinstance(value, (list, tuple)):
        return complex(*value)
    return complex(value)


def _background(samples):
    mus = np.array([float(x['mu_gev']) for x in samples])
    if len(mus) < 2 or np.any(mus <= 0) or not np.all(np.diff(mus) < 0):
        raise ValueError('Samples must descend from MF to MS')
    t = np.log(mus)[::-1]
    def arr(name, complex_matrix=False):
        if complex_matrix:
            return np.asarray([matrix_from_json(s[name]) for s in samples][::-1])
        return np.asarray([_complex_scalar(s[name]) for s in samples][::-1])
    scalars={k:arr(k) for k in ('g2','lambdaH','lambda12','lambda12Adj','lambdaT3')}
    ye=arr('ye',True); yu=arr('yu',True); yd=arr('yd',True)
    def interp(x, values):
        return np.interp(x,t,values.real)+1j*np.interp(x,t,values.imag)
    def evaluate(x):
        d={k:interp(x,v) for k,v in scalars.items()}
        d.update({k:np.array([[interp(x,v[:,i,j]) for j in range(3)] for i in range(3)])
                  for k,v in [('ye',ye),('yu',yu),('yd',yd)]})
        return d
    return float(t[0]),float(t[-1]),evaluate


def _gamma(c):
    """Numerical coefficients in the lepton-flavour basis.

    The Weinberg diagonal term is lifted *directly* from
    RGE.running.weinberg.WeinbergRGE.beta_weinberg_matrix:
       (2 lambdaH - 3 g2^2 + 2 T) K
         - 3/2 (Ye Ye^dag K + K (Ye Ye^dag)^T).
    It is applied to the intermediate scalar-only EFT because the extra
    scalars have no renormalisable lepton Yukawa interactions; this extension
    still requires a dedicated intermediate-EFT full-flavour derivation.

    The C12 lepton action is derived from the generic master RGE collinear
    fermion term. Two-operator closure remains unverified.
    """
    ye = np.asarray(c['ye'], complex)
    yd = np.asarray(c['yd'], complex)
    yu = np.asarray(c['yu'], complex)
    lepton_left = ye @ ye.conj().T
    lepton_trace = float(np.trace(lepton_left).real)
    full_trace = lepton_trace + 3.0 * float(np.trace(yd @ yd.conj().T).real) + 3.0 * float(np.trace(yu @ yu.conj().T).real)
    scalar12 = -3*c['g2']**2+2*c['lambda12']-c['lambda12Adj']/2
    scalar5 = -3*c['g2']**2+2*c['lambdaH'] + 2*full_trace
    return scalar12, scalar5, lepton_left


def _sym(m):
    return (m+m.T)/2


def _rhs(t, y, bg, *, diagonal=True, reverse=True):
    a=y[:9].reshape(3,3); b=y[9:].reshape(3,3)
    c=bg(t); l=c['lambdaT3']
    da=REV*np.conj(l)*b if reverse else np.zeros_like(a)
    db=FWD*l*a
    if diagonal:
        g12,g5,e=_gamma(c)
        da+=c12_beta_numerator(a, b, ye=c['ye'], g2=c['g2'],
             lambda12=c['lambda12'], lambda12Adj=c['lambda12Adj'],
             lambdaT3=l, reverse=False)
        db+=g5*b-1.5*(e@b+b@e.T)  # source-derived Weinberg diagonal term
    return np.concatenate((_sym(da).ravel(),_sym(db).ravel()))/LOOP


def coupled_heavy_resolved(samples, y1, y2, masses, *, diagonal=True, reverse=True):
    """Each heavy state starts at Mr. For Mr>MF numerical background is held at MF.

    This is explicitly a boundary approximation, not a derived sequential EFT.
    """
    t_low,t_high,bg=_background(samples)
    m=np.asarray(masses,float)
    if m.shape!=(3,) or np.any(m<=0) or np.any(~np.isfinite(m)) or np.log(min(m))<=t_low:
        raise ValueError('Three positive heavy masses must exceed MS')
    source=flavor_kernel_by_heavy(y1,y2,m)
    all_a=[]; all_b=[]
    for mr, a0 in zip(m,source):
        start=np.log(mr)
        y0=np.concatenate((a0.ravel(),np.zeros(9,complex)))
        sol=solve_ivp(lambda t,y:_rhs(t,y,bg,diagonal=diagonal,reverse=reverse),
                      (start,t_low),y0,method='DOP853',rtol=2e-10,atol=1e-14)
        if not sol.success:raise RuntimeError(sol.message)
        all_a.append(_sym(sol.y[:9,-1].reshape(3,3)))
        all_b.append(_sym(sol.y[9:,-1].reshape(3,3)))
    return np.sum(all_a,axis=0),np.sum(all_b,axis=0)



def common_scale_log_boundary(samples, y1, y2, masses, *,
                              diagonal=True, reverse=True,
                              lambda_reference=None):
    """Run solely on the known scalar-only EFT interval [MF, MS].

    Encode the mass-resolved one-loop threshold logs ln(MF/Mr) as an
    effective C5(MF) boundary while seeding C12(MF)=sum_r C12_r. This
    reproduces the production hierarchical direct logarithm to O(1/16pi²)
    when all beta-function coefficients are frozen at the reference value.

    This is a common-scale matching-log approximation, NOT an EFT with
    1 or 2 active heavy fermions between distinct Mr.

    C5 here is the *ordered-tensor* convention. The caller applies the
    established factor of two for the physical Majorana coefficient.
    """
    t_low, t_high, bg = _background(samples)
    mf = float(samples[0]['mu_gev'])
    ms = float(samples[-1]['mu_gev'])
    masses = np.asarray(masses, dtype=float)
    if (masses.shape != (3,) or not np.all(np.isfinite(masses)) or
            np.any(masses <= ms)):
        raise ValueError('Three finite heavy masses must all exceed MS')
    if lambda_reference is None:
        # Production fixed-order direct matching evaluates lambdaT3 at MS.
        lambda_reference = _complex_scalar(samples[-1]['lambdaT3'])
    lambda_reference = complex(lambda_reference)
    if not (np.isfinite(lambda_reference.real) and
            np.isfinite(lambda_reference.imag)):
        raise ValueError('Reference lambdaT3 must be finite')
    sources = flavor_kernel_by_heavy(y1, y2, masses)
    gap_logs = np.log(mf/masses)
    initial_c12 = _sym(sources.sum(axis=0))
    initial_c5 = _sym(FWD*lambda_reference/LOOP *
                      np.einsum('r,rij->ij', gap_logs, sources))
    y0 = np.concatenate((initial_c12.ravel(), initial_c5.ravel()))
    sol = solve_ivp(
        lambda t, y: _rhs(t, y, bg, diagonal=diagonal, reverse=reverse),
        (t_high, t_low), y0, method='DOP853',
        rtol=2e-10, atol=1e-14,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    c12 = _sym(sol.y[:9, -1].reshape(3, 3))
    c5 = _sym(sol.y[9:, -1].reshape(3, 3))
    details = {
        'mass_resolved_log_mf_over_mr': gap_logs.tolist(),
        'heavy_boundary_c5_ordered': matrix_to_json(initial_c5),
        'heavy_boundary_c12': matrix_to_json(initial_c12),
        'lambda_t3_reference_at_ms': [lambda_reference.real,
                                      lambda_reference.imag],
        'evolution_interval_gev': [mf, ms],
        'sequential_partial_decoupling_implemented': False,
    }
    return c12, c5, details


def fixed_order_log_from_common_boundary(samples, y1, y2, masses):
    """Analytic O(hbar) result with lambda frozen at the MS reference.

    This must equal the independently validated production direct term,
    including the physical-Majorana factor two and all three heavy logs.
    """
    mf = float(samples[0]['mu_gev'])
    ms = float(samples[-1]['mu_gev'])
    masses = np.asarray(masses, dtype=float)
    sources = flavor_kernel_by_heavy(y1, y2, masses)
    lam = _complex_scalar(samples[-1]['lambdaT3'])
    boundary_logs = np.log(mf/masses)
    eft_running_log = np.log(ms/mf)
    ordered = FWD*lam/LOOP * (
        np.einsum('r,rij->ij', boundary_logs, sources)
        + sources.sum(axis=0)*eft_running_log
    )
    return _sym(2.0*ordered)

def calculate(control, normalization, y1, y2, masses, *,
              diagonal=True, threshold_scheme="common_scale_log_boundary",
              compare_legacy_threshold=False):
    if control.get('status')!='FixedOrderControlComplete':
        raise ValueError('Require successful fixed-order control')
    if normalization.get('status')!='UnfittedNormalizationPassed':
        raise ValueError('Absolute direct leading-log comparison has not passed')
    samples=control['intermediate_running']['samples']
    mf=float(control['mu_f_gev']); ms=float(control['mu_s_gev'])
    masses=np.asarray(masses,float)
    if np.any(masses<=ms):raise ValueError('Incorrect heavy-mass hierarchy')
    norm_direct=matrix_from_json(normalization['production_direct_c5'])
    predicted=matrix_from_json(normalization['predicted_direct_c5'])
    err=np.linalg.norm(norm_direct-predicted)/np.linalg.norm(norm_direct)
    if err>1e-9:raise ValueError('Production direct leading-log absolute-normalisation failed')
    if not np.isclose(float(samples[0]['mu_gev']),mf) or not np.isclose(float(samples[-1]['mu_gev']),ms):
        raise ValueError('Background scale mismatch')
    if threshold_scheme not in ('common_scale_log_boundary', 'endpoint_extrapolation'):
        raise ValueError('Unknown threshold scheme')
    analytic_direct = fixed_order_log_from_common_boundary(samples,y1,y2,masses)
    ll_residual = float(np.linalg.norm(analytic_direct-norm_direct) /
                        np.linalg.norm(norm_direct))
    if ll_residual > 1e-8:
        raise ValueError('Heavy-boundary leading-order result does not reproduce production direct C5')
    if threshold_scheme == 'common_scale_log_boundary':
        c12,c5,threshold_details=common_scale_log_boundary(
            samples,y1,y2,masses,diagonal=diagonal,reverse=True)
        _,forward,_=common_scale_log_boundary(
            samples,y1,y2,masses,diagonal=diagonal,reverse=False)
    else:
        c12,c5=coupled_heavy_resolved(
            samples,y1,y2,masses,diagonal=diagonal,reverse=True)
        _,forward=coupled_heavy_resolved(
            samples,y1,y2,masses,diagonal=diagonal,reverse=False)
        threshold_details = {
            'sequential_partial_decoupling_implemented': False,
            'approximation': 'Scalar-only background clamped at MF outside its domain',
        }
    # Physical Majorana normalization: factor two relative to ordered tensor.
    evolved=2*c5; forward=2*forward
    baseline=matrix_from_json(control['fixed_order_c5_at_ms'])
    hard=baseline-norm_direct
    candidate=hard+evolved
    denominator=np.linalg.norm(baseline)
    corrected_delta=candidate-baseline
    legacy_comparison=None
    if compare_legacy_threshold and threshold_scheme == 'common_scale_log_boundary':
        _,legacy_c5=coupled_heavy_resolved(
            samples,y1,y2,masses,diagonal=diagonal,reverse=True)
        legacy_evolved=2*legacy_c5
        legacy_candidate=hard+legacy_evolved
        legacy_comparison={
            'legacy_relative_difference_from_fixed_order':float(
                np.linalg.norm(legacy_candidate-baseline)/denominator),
            'common_domain_minus_legacy_relative_to_fixed_order':float(
                np.linalg.norm(candidate-legacy_candidate)/denominator),
            'common_domain_minus_legacy_matrix':matrix_to_json(
                candidate-legacy_candidate),
            'meaning':'Scheme sensitivity, not a calculated matching uncertainty',
        }
    return {
        'status':'ExploratoryCoupledMatrixComplete',
        'interpretation':'NOT a validated physical neutrino mass prediction',
        'c12_diagonal_flavour_formula':('master fermion-collinear: 0.5*(YeYeDagger*C12 + C12*(YeYeDagger).T)' if diagonal else 'disabled'),
        'c12_flavour_source':'RGE/general/WilsonTensorRGE.py + RGE/general/AnomalousDimensions.py',
        'two_operator_closure_verified':False,
        'threshold_scheme':threshold_scheme,
        'threshold_scheme_details':threshold_details,
        'common_boundary_fixed_order_direct_residual':ll_residual,
        'legacy_threshold_comparison':legacy_comparison,
        'c5_diagonal_formula_source':'RGE/running/weinberg/WeinbergRGE.py: beta_weinberg_matrix',
        'heavy_background_above_mf':('not used' if threshold_scheme == 'common_scale_log_boundary'
                                    else 'endpoint held fixed (approximation)'),
        'production_direct_relative_residual':float(err),
        'c12_at_ms':matrix_to_json(c12),
        'coupled_eft1_c5_at_ms':matrix_to_json(evolved),
        'forward_only_eft1_c5_at_ms':matrix_to_json(forward),
        'production_hard_c5_at_ms':matrix_to_json(hard),
        'production_direct_c5_at_ms':matrix_to_json(norm_direct),
        'exploratory_total_c5_at_ms':matrix_to_json(candidate),
        'exploratory_minus_fixed_c5_at_ms':matrix_to_json(corrected_delta),
        'relative_difference_from_fixed_order':float(np.linalg.norm(corrected_delta)/denominator),
        'relative_reverse_feedback':float(np.linalg.norm(evolved-forward)/np.linalg.norm(forward)),
        'c5_symmetry_residual':float(np.linalg.norm(candidate-candidate.T)),
        'validated_physical_prediction_available':False,
        'missing':['Complete three-generation operator-space closure',
                   'Proper sequential heavy-threshold renormalisable matching',
                   'Consistent higher-order scalar threshold matching']
    }



def exploratory_mass_comparison(fixed_c5, candidate_c5, ye_at_ms, *, vev_gev=246.22):
    """Compare Takagi neutrino masses at MS in the charged-lepton mass basis.

    This reuses production C5-to-mass conventions, but does NOT validate
    the exploratory Wilson RGE or propagate the result down to low energy.
    """
    from physics.NeutrinoTrajectory import (
        neutrino_mass_from_c5, charged_lepton_mass_basis_matrix,
    )
    from physics.NeutrinoObservables import takagi_factorization

    result = {}
    for name, c5 in (("fixed_order", fixed_c5), ("exploratory", candidate_c5)):
        c5 = np.asarray(c5, dtype=complex)
        if c5.shape != (3, 3) or not np.all(np.isfinite(c5)):
            raise ValueError(f"{name} C5 must be a finite complex 3x3 matrix")
        mass_gev = neutrino_mass_from_c5(c5, vev_gev=vev_gev)
        mass_flavor_gev = charged_lepton_mass_basis_matrix(mass_gev, ye_at_ms)
        _, masses_gev, residual = takagi_factorization(mass_flavor_gev)
        masses_ev = np.asarray(masses_gev, dtype=float) * 1.0e9
        result[name] = {
            "masses_ev": masses_ev.tolist(),
            "takagi_residual": float(residual),
            "delta_m21_sq_ev2": float(masses_ev[1]**2 - masses_ev[0]**2),
            "delta_m31_sq_ev2": float(masses_ev[2]**2 - masses_ev[0]**2),
        }
    a = np.asarray(result["fixed_order"]["masses_ev"])
    b = np.asarray(result["exploratory"]["masses_ev"])
    result["mass_differences_ev"] = (b-a).tolist()
    result["relative_mass_differences"] = [
        None if aa == 0 else float((bb-aa)/aa)
        for aa, bb in zip(a,b)
    ]
    result["comparison_scale"] = "MS (not low energy)"
    result["validated_physical_prediction_available"] = False
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('config','uv-rgbeta','control','normalization','output'):
        p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--offdiagonal-only',action='store_true',help='Disable diagonal gamma terms')
    p.add_argument('--threshold-scheme',choices=('common_scale_log_boundary',
        'endpoint_extrapolation'),default='common_scale_log_boundary')
    p.add_argument('--compare-legacy-threshold',action='store_true',
        help='Also evaluate the old endpoint-extrapolation approximation')
    a=p.parse_args()
    control=json.loads(a.control.read_text(encoding='utf-8'))
    norm=json.loads(a.normalization.read_text(encoding='utf-8'))
    # Extract physical heavy mass basis exactly from the current canonical adapter.
    # Do not infer Yukawas from the C12 sum (information is lost).
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FinalC5TrajectoryAdapter import _majorana_takagi_mass_basis
    from types import SimpleNamespace
    cfg=json.loads(a.config.read_text(encoding='utf-8'))
    raw=dict(cfg);raw.setdefault('scan',{'bindings':{}})
    state=build_uv_state_from_config(SimpleNamespace(raw=raw),{}).validated()
    uv=run_uv_segment(state,load_rgbeta_payload(a.uv_rgbeta),float(control['mu_f_gev']))
    basis=_majorana_takagi_mass_basis(uv.final_state,symmetry_atol=1e-10,residual_rtol=1e-9,residual_atol=1e-10)
    result=calculate(control,norm,basis.y1,basis.y2,basis.masses_gev,
                     diagonal=not a.offdiagonal_only,
                     threshold_scheme=a.threshold_scheme,
                     compare_legacy_threshold=a.compare_legacy_threshold)
    result['heavy_masses_gev']=basis.masses_gev.tolist()
    # Derive the scalar-threshold charged-lepton Yukawa matrix from the
    # saved intermediate trajectory, not the UV boundary.
    samples=control['intermediate_running']['samples']
    from Numerical.diagnostics.T3BHeavyLogAlignment import matrix_from_json
    ye_at_ms=matrix_from_json(samples[-1]['ye'])
    baseline=matrix_from_json(control['fixed_order_c5_at_ms'])
    exploratory=matrix_from_json(result['exploratory_total_c5_at_ms'])
    result['exploratory_neutrino_mass_comparison_at_ms'] = exploratory_mass_comparison(
        baseline, exploratory, ye_at_ms,
        vev_gev=float(cfg.get('numerical',{}).get('vev_gev',246.22)),
    )
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    summary={k:result[k] for k in ('status','threshold_scheme',
            'common_boundary_fixed_order_direct_residual',
            'relative_difference_from_fixed_order','relative_reverse_feedback',
            'c5_symmetry_residual','validated_physical_prediction_available')}
    legacy=result['legacy_threshold_comparison']
    if legacy is not None:
        summary['legacy_threshold_comparison']={
            key:legacy[key] for key in (
                'legacy_relative_difference_from_fixed_order',
                'common_domain_minus_legacy_relative_to_fixed_order')
        }
    print(json.dumps(summary,indent=2))
    print(json.dumps({'exploratory_neutrino_mass_comparison_at_ms':result['exploratory_neutrino_mass_comparison_at_ms']},indent=2))

if __name__=='__main__':main()
