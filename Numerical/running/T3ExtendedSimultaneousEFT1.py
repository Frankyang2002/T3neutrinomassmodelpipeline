"""Opt-in joint EFT1 evolution with complete running quartic tensor.

The exact existing seven-term Wilson evaluator is reused. No RGBeta scalar
normalisation is inferred. This module is INACTIVE unless explicitly called
with a supplied, validated RunningQuarticBasis.

The joint solution evolves [renormalisable running vector, complex master
Wilson tensor] on ONE DOP853 IVP. In particular, it does NOT evolve the
fermion-decoupled Wilson tensor with frozen quartics.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Mapping, Any
import math
import numpy as np
from scipy.integrate import solve_ivp
from Numerical.running.T3RunningQuarticTensor import (
    RunningQuarticBasis, QuarticEvidenceError)
from Numerical.running.T3QuarticWilsonBridge import running_seven_term_beta


@dataclass(frozen=True)
class JointEFT1Adapter:
    """Existing production state boundary, supplied by the respective model.

    unpack: (ren_vector, mu_gev) -> (running dict, HHSS tensor)
    derivative: (ln_mu, ren_vector) -> d(ren_vector)/dln_mu
    context: (running_dict, HHSS, mu_gev) -> (Wilson Context, extra gamma)
    """
    model: str
    mu_f: float
    mu_s: float
    unpack: Callable
    derivative: Callable
    context: Callable
    required_quartics: tuple[str, ...] = ()


def _pack_complex(c: np.ndarray) -> np.ndarray:
    return np.r_[np.asarray(c.real).ravel(),np.asarray(c.imag).ravel()]


def _unpack_complex(z: np.ndarray, shape: tuple[int,...]) -> np.ndarray:
    n = math.prod(shape)
    if np.asarray(z).shape != (2*n,):
        raise QuarticEvidenceError('Invalid master Wilson packed shape')
    return (z[:n] + 1j*z[n:]).reshape(shape)


def integrate_extended_eft1(
    adapter: JointEFT1Adapter, start_ren: np.ndarray,
    initial_wilson: np.ndarray, basis: RunningQuarticBasis,
    *, loop_factor: float = 1/(16*math.pi**2),
    rtol: float = 2e-7, atol: float = 1e-10,
    max_step: float = 0.35,
    check_hhss: bool = True,
    eft_on: bool = True,
) -> dict[str,Any]:
    """Run certified-response-capable joint DOP853. No extra EFT threshold.

    A caller must supply *all* independent baseline and model-specific scalar
    coupling response tensors. The HHSS check is mandatory for production;
    check_hhss=False is provided ONLY for controlled synthetic tests.
    """
    if adapter.required_quartics and set(basis.symbols) != set(adapter.required_quartics):
        raise QuarticEvidenceError(
            'Full quartic response inventory mismatch: missing=' +
            str(sorted(set(adapter.required_quartics)-set(basis.symbols))) +
            ' excess=' + str(sorted(set(basis.symbols)-set(adapter.required_quartics))))
    if adapter.model != basis.model:
        raise QuarticEvidenceError('Joint adapter and quartic model mismatch')
    if not (0 < adapter.mu_s < adapter.mu_f):
        raise QuarticEvidenceError('Expected EFT1 fermion->scalar flow muF>muS>0')
    if not (math.isfinite(loop_factor) and 0 < loop_factor < 1):
        raise QuarticEvidenceError('Invalid Wilson loop prefactor')
    ren0 = np.asarray(start_ren,dtype=float)
    if ren0.ndim != 1 or not ren0.size or not np.isfinite(ren0).all():
        raise QuarticEvidenceError('Invalid initial renormalisable vector')
    c0=np.asarray(initial_wilson,dtype=complex)
    if c0.ndim!=4 or c0.shape[0]!=c0.shape[1] or c0.shape[2]!=c0.shape[3]:
        raise QuarticEvidenceError('Expected rank-four master Wilson tensor')
    norm=float(np.linalg.norm(c0))
    if not (0<norm<math.inf):
        raise QuarticEvidenceError('Initial Wilson coefficient norm zero/nonfinite')
    nren=ren0.size
    shape=c0.shape
    if nren+2*c0.size>1_000_000:
        raise QuarticEvidenceError('Unexpectedly large combined EFT1 system')
    z0=np.r_[ren0,_pack_complex(c0/norm)]
    calls=0

    def evaluate(t, ren, c):
        state, hhss=adapter.unpack(ren,float(np.exp(t)))
        if not isinstance(state,Mapping):
            raise QuarticEvidenceError('State adapter must expose named coupling map')
        context,gamma=adapter.context(state,hhss,float(np.exp(t)))
        if context.model!=adapter.model:
            raise QuarticEvidenceError('Wilson model mismatch during joint running')
        expected=hhss if check_hhss else None
        beta,terms=running_seven_term_beta(
            context,c,state,basis,expected_hhss=expected,
            extra_scalar_anomalous_dimension=gamma)
        return beta,terms

    _, initial_terms=evaluate(math.log(adapter.mu_f),ren0,c0)
    initial_source=np.array(initial_terms['scalar_pair'],copy=True)
    def rhs(t,z):
        nonlocal calls
        ren=z[:nren]
        c=_unpack_complex(z[nren:],shape)
        beta,_=evaluate(t,ren,c)
        dr=(np.asarray(adapter.derivative(t,ren),float) if eft_on
            else np.zeros_like(ren))
        if dr.shape!=ren.shape or not np.isfinite(dr).all():
            raise QuarticEvidenceError('Invalid renormalisable RGBeta derivative')
        calls+=1
        return np.r_[dr,_pack_complex(loop_factor*beta)]

    sol=solve_ivp(rhs,(math.log(adapter.mu_f),math.log(adapter.mu_s)),z0,
                  method='DOP853',rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success:
        raise RuntimeError('Extended simultaneous EFT1 failed: '+sol.message)
    final_ren=np.array(sol.y[:nren,-1],copy=True)
    final_c=norm*_unpack_complex(sol.y[nren:,-1],shape)
    if not np.isfinite(final_ren).all() or not np.isfinite(final_c).all():
        raise QuarticEvidenceError('Nonfinite joint EFT1 endpoint')
    sf=float(np.linalg.norm(final_c-final_c.swapaxes(0,1))/norm)
    ss=float(np.linalg.norm(final_c-final_c.swapaxes(2,3))/norm)
    if max(sf,ss)>1e-6:
        raise QuarticEvidenceError(f'Master pair symmetry failed: {sf} {ss}')
    return {'end_ren_vector':final_ren,'end_c':final_c,
            'nfev':int(sol.nfev),'steps':int(sol.t.size),'rhs_calls':calls,
            'initial_scalar_pair':initial_source,
            'fermion_symmetry':sf,'scalar_symmetry':ss,
            'tensor_drift':float(np.linalg.norm(final_c-c0)/norm),
            'quartic_tensor_source':basis.source,
            'complete_quartic_response_used':True,
            'ren_beta_active':bool(eft_on)}
