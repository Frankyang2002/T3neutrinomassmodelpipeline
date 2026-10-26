"""Stage 18: time-resolved extra-quartic contribution to EFT1 LLHH Wilson tensors.

Experimental *source insertion*, not the full homogeneous/mixing EFT1 Wilson
system. The one-generation Stage17 PL seed is lifted by the same full-flavour
symmetric heavy-index contraction as its tree-level Matchete source. The
matched Wilson seed is held fixed at the F boundary; the TWO complex quartics
run with the already-implemented full renormalisable EFT1 RGBeta equations.

No conversion from the real (i,j,a,b) Wilson tensor to a canonical Weinberg
coefficient is made. No low-energy neutrino mass is predicted here.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import sympy as sp
from scipy.integrate import simpson, trapezoid

from Numerical.experimental.GeneralT3ExtraEFT1Mixing import ACTIVE_RGBETA


SOURCE_STATUS = 'ExtraEFT1ScalarPairMixingDerived'
SOURCES = {'y1': sp.Symbol('y1'), 'y2': sp.Symbol('y2')}
MF = sp.Symbol('MF')


def _finite_scalar(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError) as err:
        raise ValueError(f'{label}: finite real scalar required') from err
    if not np.isfinite(result):
        raise ValueError(f'{label}: nonfinite value')
    return result


def _encode(matrix: np.ndarray) -> dict[str, Any]:
    a = np.asarray(matrix, dtype=complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError('Wilson tensor flavour component must be a finite complex 3x3 matrix')
    return {'real': a.real.tolist(), 'imag': a.imag.tolist(),
            'frobenius_norm_gev_inverse': float(np.linalg.norm(a, 'fro'))}


@dataclass(frozen=True)
class SourceComponent:
    higgs_pair: tuple[int, int]
    # The coefficient of conjugate(y1)^2/MF and conjugate(y2)^2/MF.
    # Evaluators take the potential-level RGBeta lambdas (not Matchete signs).
    y1_factor: Any
    y2_factor: Any
    analytic: tuple[str, str]

    def evaluate(self, lambda_s2: complex, lambda_s1: complex,
                 flavor_y1: np.ndarray, flavor_y2: np.ndarray) -> np.ndarray:
        params = (complex(lambda_s2), complex(lambda_s1))
        f1 = complex(self.y1_factor(*params))
        f2 = complex(self.y2_factor(*params))
        result = f1 * flavor_y1 + f2 * flavor_y2
        if not np.isfinite(result).all():
            raise ValueError('Extra-EFT1 Wilson source has nonfinite components')
        return result


def compile_stage17_sources(model: str, payload: Mapping[str, Any]) -> tuple[SourceComponent, ...]:
    """Validate and compile the certified Stage17 formal beta expressions.

    Reject any term that does not factor as conj(y_i)^2 / MF times a function
    linear in the two complex quartics, including their conjugates.
    """
    if model not in ('B_m1', 'C_m1'):
        raise ValueError('Only B/C alpha=-1 supported')
    if payload.get('status') != SOURCE_STATUS or payload.get('model') != model:
        raise ValueError('Missing matching successful Stage17 report')
    if payload.get('source_convention') != 'lambda_Matchete=-lambda_RGBeta, same named scalar field basis':
        raise ValueError('Stage17 source-convention certificate mismatch')
    if payload.get('full_eft1_wilson_running_implemented') is not False:
        raise ValueError('Stage17 input has an unexpected completed-Wilson claim')
    if payload.get('extra_low_energy_neutrino_mass') is not None:
        raise ValueError('Stage17 input should not contain a low-energy mass')

    rows = payload.get('components')
    if not isinstance(rows, list) or len(rows) != payload.get('nonzero_weinberg_real_components'):
        raise ValueError('Inconsistent Stage17 Higgs component count')
    if len(rows) != 4:
        raise ValueError('Stage17 expected four independent nonzero Higgs-real components')
    if payload.get('master_rge_term') != 'delta_beta_ijab^(1) = (1/2) sum_swap(a,b),c,d delta_lambda_abcd C_ijcd':
        raise ValueError('Unrecognized Stage17 master RGE convention')

    y1, y2 = SOURCES['y1'], SOURCES['y2']
    b1, b2 = sp.symbols('yb1 yb2')
    l2, l1 = (sp.Symbol(name) for name in ACTIVE_RGBETA)
    lambda_bars = sp.symbols('l2bar l1bar')
    variables = (l2, l1, *lambda_bars)
    quartic_functions: list[SourceComponent] = []
    seen: set[tuple[int,int]] = set()
    for row in rows:
        indices = row.get('higgs_real_indices')
        if (not isinstance(indices, list) or len(indices) != 2 or
            not all(type(i) is int and 1 <= i <= 4 for i in indices) or
            indices[0] > indices[1]):
            raise ValueError('Invalid Stage17 Higgs real indices')
        pair = (indices[0], indices[1])
        if pair in seen:
            raise ValueError('Repeated Higgs real component')
        seen.add(pair)

        expression = row.get('beta_16pi2_expression')
        if not isinstance(expression, str) or not expression.strip():
            raise ValueError('Missing Stage17 beta expression')
        try:
            symbolic = sp.sympify(expression, locals={
                'conjugate': sp.conjugate, 'sqrt': sp.sqrt,
                'lambdaHHdagS2S2': l2, 'lambdaHHdagS1barS1bar': l1,
                'y1': y1, 'y2': y2, 'MF': MF,
            })
        except (TypeError, ValueError, SyntaxError) as err:
            raise ValueError('Stage17 beta expression could not be parsed') from err
        expression_b = sp.expand(symbolic.xreplace({sp.conjugate(y1): b1, sp.conjugate(y2): b2}))
        polynomial = sp.Poly(expression_b, y1, y2, b1, b2)
        parts = {'y1': sp.S.Zero, 'y2': sp.S.Zero}
        for monomial, coeff in polynomial.terms():
            if monomial == (0, 0, 2, 0):
                family = 'y1'
            elif monomial == (0, 0, 0, 2):
                family = 'y2'
            else:
                raise ValueError(f'Unsupported Stage17 Yukawa monomial {monomial}')
            parts[family] += coeff
        if not any(parts.values()):
            raise ValueError('Zero Stage17 beta component is not allowed in the nonzero list')
        factors = []
        analytical = []
        for family in ('y1','y2'):
            # Stage17's seed contains precisely one inverse MF. Divide it out
            # before lifting a one-generation diagonal term to full flavour.
            coefficient = sp.cancel(parts[family] * MF)
            coefficient = coefficient.xreplace({sp.conjugate(l2):lambda_bars[0],
                                                sp.conjugate(l1):lambda_bars[1]})
            if coefficient.free_symbols - set(variables):
                raise ValueError(f'Unsupported symbols in {family} component: {coefficient.free_symbols-set(variables)}')
            p = sp.Poly(sp.expand(coefficient), *variables)
            if p.total_degree() > 1 or p.coeff_monomial(1) != 0:
                raise ValueError('Extra scalar-pair source must be linear and zero at zero extra couplings')
            # Use conjugation as in Stage17, not a phase guessed at the bridge.
            evaluator = sp.lambdify((l2,l1), coefficient.xreplace({
                lambda_bars[0]:sp.conjugate(l2), lambda_bars[1]:sp.conjugate(l1)}), 'numpy')
            factors.append(evaluator)
            analytical.append(sp.sstr(coefficient))
        quartic_functions.append(SourceComponent(pair, factors[0], factors[1], tuple(analytical)))
    return tuple(sorted(quartic_functions, key=lambda x:x.higgs_pair))


def frozen_seed_flavor_contractions(masses: np.ndarray, y1: np.ndarray,
                                    y2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Full flavour lift of the Stage17 one-generation PL tree Wilson seed.

    PL seed has conjugate(y_i)^2/MF. With diagonal positive mass basis, its
    symmetric flavour extension is Sum_r conj(y_i[p,r]) conj(y_i[q,r])/M_r.
    """
    m = np.asarray(masses, dtype=float)
    if m.shape != (3,) or not np.isfinite(m).all() or np.any(m <= 0):
        raise ValueError('Three finite positive heavy masses required')
    result = []
    for label, yy in (('y1', y1), ('y2', y2)):
        yy = np.asarray(yy, dtype=complex)
        if yy.shape != (3,3) or not np.isfinite(yy).all():
            raise ValueError(f'{label} must be a finite 3x3 matrix')
        c = np.conjugate(yy)
        form = (c / m[None,:]) @ c.T
        if not np.allclose(form,form.T,rtol=1e-12,atol=1e-28):
            raise ValueError('Wilson PL flavour contraction is not symmetric')
        result.append(form)
    return tuple(result)


def integrate_source_over_eft1(components: tuple[SourceComponent, ...],
                                trajectory, masses: np.ndarray,
                                y1: np.ndarray, y2: np.ndarray) -> dict[str, Any]:
    """Integrate *only* the extra-quartic source from F down to scalar threshold.

    The EFT1 state and its quartics run at one loop, while the Stage17 PL
    tree Wilson coefficients are fixed at the F matching boundary. This is
    not the complete Wilson-system evolution or a final C5 matching.
    """
    if len(components) != 4:
        raise ValueError('Expected the four Stage17 nonzero Higgs-real components')
    mu = np.asarray(trajectory.mu_gev, dtype=float)
    if (mu.ndim != 1 or mu.size < 5 or mu.size % 2 != 1 or
        not np.isfinite(mu).all() or np.any(mu <= 0) or
        not np.all(np.diff(mu) < 0)):
        raise ValueError('EFT1 trajectory needs odd >=5 strictly descending positive save scales')
    t = np.log(mu)
    f1, f2 = frozen_seed_flavor_contractions(masses,y1,y2)
    source = np.zeros((len(components),mu.size,3,3),dtype=complex)
    sample_quartics=[]
    for k in range(mu.size):
        state = trajectory.state_at_index(k)
        if state.representation.shared_scalar or (state.representation.d_s1,state.representation.d_s2,state.representation.alpha)!=(2,2,-1):
            raise ValueError('Stage18 EFT1 trajectory is not an ordinary B/C(-1) doublet theory')
        if not np.isclose(float(state.mu_gev),mu[k],rtol=1e-10,atol=1e-8):
            raise ValueError('EFT1 trajectory state scale mismatch')
        l2 = complex(state.lambdaHHdagS2S2)
        l1 = complex(state.lambdaHHdagS1barS1bar)
        if not np.isfinite([l2.real,l2.imag,l1.real,l1.imag]).all():
            raise ValueError('Nonfinite extra quartic along EFT1 trajectory')
        sample_quartics.append({'mu_gev':float(mu[k]), 'lambdaHHdagS2S2':[l2.real,l2.imag],
                                'lambdaHHdagS1barS1bar':[l1.real,l1.imag]})
        for j,component in enumerate(components):
            source[j,k] = component.evaluate(l2,l1,f1,f2)
    loop_factor = 16*np.pi**2
    integrals = simpson(source,x=t,axis=1)/loop_factor
    trapezoids = trapezoid(source,x=t,axis=1)/loop_factor
    errors = integrals-trapezoids
    if not np.isfinite(integrals).all():
        raise ValueError('Integrated Wilson tensor is nonfinite')
    if not np.allclose(integrals,integrals.transpose(0,2,1),rtol=1e-10,atol=1e-25):
        raise ValueError('Integrated Wilson tensor lost Majorana flavour symmetry')
    rows = []
    for j,component in enumerate(components):
        rows.append({'higgs_real_indices':list(component.higgs_pair),
                     'integrated_delta_C_ijab_gev_inverse':_encode(integrals[j]),
                     'quadrature_simpson_minus_trapezoid_gev_inverse':_encode(errors[j]),
                     'coefficient_y1':component.analytic[0],
                     'coefficient_y2':component.analytic[1]})
    return {'status':'ExtraEFT1LLHHTensorSourceIntegrated',
            'source':'Stage17 PL tree-seed scalar-pair beta, Stage16 running EFT1 quartics',
            'scales_gev':{'mu_F':float(mu[0]),'mu_S':float(mu[-1])},
            'eft1_samples':int(mu.size),
            'eft1_source_quartics':sample_quartics,
            'loop_order_of_source_insertion':1,
            'heavy_mass_basis':'positive diagonal masses; PL Yukawa products conjugated',
            'eft1_wilson_seed_frozen_at_mu_F':True,
            'extra_source_only':True,
            'full_eft1_wilson_anomalous_dimension_implemented':False,
            'weinberg_operator_projection_verified':False,
            'scalar_threshold_matching_verified':False,
            'extra_physical_c5_at_mu_s':None,
            'extra_neutrino_mass_low_ev':None,
            'higgs_real_components':rows,
            'complete_two_loop_prediction':False,
            'note':'Numerical leading insertion of certified extra-scalar quartic source into tree PL LLSS Wilson coefficients, with running EFT1 quartics. Not full coupled Wilson running, and real-component coefficients must not be identified with Weinberg C5 without a calibrated projection and threshold matching.'}


def run_stage18_model(model: str, config: Mapping[str,Any],
                      uv_rgbeta: Mapping[str,Any], eft_rgbeta: Mapping[str,Any],
                      stage17: Mapping[str,Any], *, n_save: int = 17) -> dict[str,Any]:
    """Run the actual existing UV + EFT1 ODE layers, then source integration."""
    if n_save < 5 or n_save > 129 or n_save % 2 != 1:
        raise ValueError('n_save must be odd and within [5,129] for Simpson quadrature')
    from Numerical.experimental.GeneralT3CoupledUV import construct_uv_state
    from Numerical.experimental.GeneralT3HierarchicalTrajectory import (
        _threshold_scales,_guard_running_hierarchy,_validate_payload)
    from Numerical.running.T3Trajectory import run_renormalisable_t3_trajectory
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data
    parsed = compile_stage17_sources(model,stage17)
    state=construct_uv_state(model,config)
    scales=_threshold_scales(config,state)
    _validate_payload(uv_rgbeta,'UV',model)
    _validate_payload(eft_rgbeta,'EFT1',model)
    path=run_renormalisable_t3_trajectory(
        state,uv_rgbeta,eft_rgbeta,
        mu_fermion_threshold_gev=scales.mu_f,
        mu_scalar_threshold_gev=scales.mu_s,
        uv_save_scales_gev=np.geomspace(state.mu_gev,scales.mu_f,9),
        eft1_save_scales_gev=np.geomspace(scales.mu_f,scales.mu_s,n_save),
        rtol=3e-8,atol=3e-12,
    )
    masses_report=_guard_running_hierarchy(path,scales)
    heavy=_heavy_mass_basis_data(path.uv_threshold_state,
            offdiagonal_atol=1e-10,imaginary_atol=1e-10,
            takagi_residual_rtol=1e-10,takagi_residual_atol=1e-10)
    result=integrate_source_over_eft1(parsed,path.eft1,heavy.masses_gev,heavy.y1,heavy.y2)
    result.update({'model':model,'uv_nfev':int(path.uv.nfev),
                   'eft1_nfev':int(path.eft1.nfev),
                   'running_masses_gev':masses_report,
                   'three_lepton_flavors':True,
                   'heavy_flavors':3,'source_matching_validation':'Stage17 source report',
                   'mass_basis_residual':float(heavy.takagi_residual)})
    return result
