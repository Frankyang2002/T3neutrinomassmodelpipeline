"""Opt-in production-state adapters for the extended EFT1 DOP853 IVP.

The existing entry points remain untouched. To opt in, pass the adapter,
initial renormalisable vector, initial master C tensor, and source-backed
RunningQuarticBasis to integrate_extended_eft1.
"""
from __future__ import annotations
import numpy as np
from Numerical.running.T3ExtendedSimultaneousEFT1 import JointEFT1Adapter


def ordinary_adapter(model, initial_state, eft1_beta, registry, source_c,
                     q_basis, mu_s):
    """Pack *all* new RGBeta scalar quartics, including the ordinary models.

    `initial_state` is a GeneralT3ScalarRGState at muF. RGBeta input must be
    the *general* EFT1 payload, not the baseline-only betas.
    """
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState,pack_general_scalar_state,
        unpack_general_scalar_state,evaluate_general_scalar_derivative)
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    from scripts import RunT3SimultaneousEFT1All16 as general
    if model in ('B_-1','C_-1'):
        raise ValueError('Mixed model must use mixed_adapter')
    if not isinstance(initial_state,GeneralT3ScalarRGState):
        raise TypeError('Pass extended GeneralT3ScalarRGState, not baseline state')
    state=initial_state.validated(registry)
    ren0,layout=pack_general_scalar_state(state,registry)
    mu_f=float(state.baseline.mu_gev)
    def unpack(ren,mu):
        current=unpack_general_scalar_state(
            ren,layout,registry,uv=False,mu_gev=mu)
        running=_state_values(current.baseline)
        running.update(current.extra)
        return running,q_basis.tensor(running)
    def derivative(t,ren):
        current=unpack_general_scalar_state(
            ren,layout,registry,uv=False,mu_gev=float(np.exp(t)))
        result,layout_current=evaluate_general_scalar_derivative(
            eft1_beta,current,registry)
        if tuple(layout_current.names)!=tuple(layout.names):
            raise ValueError('General RGBeta running layout changed')
        return result
    def context(running,hhss,mu):
        return general.general_initial_context(
            model,running,hhss,source_c,mu_f,mu_s)
    required=tuple(sorted(name for name in layout.names if name.startswith('lambda')))
    return JointEFT1Adapter(model,mu_f,float(mu_s),unpack,derivative,context,required),ren0


def mixed_adapter(model, start_ren, eft_layout, legacy_registry,
                  eft1_beta, new, reused, source_original, q_basis,mu_f,mu_s):
    """Reuse the source-verified B/C bilinear, mixed-native RGBeta state.

    Does not repack the bilinear entries; retains their original fixed-basis
    evolution in the existing mixed EFT1 solver's `core.rhs_at`.
    """
    from scripts import RunT3SimultaneousEFT1MixedBC as mixed
    from scripts import RunT3SimultaneousEFT1All16 as general
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    if model not in ('B_-1','C_-1'):
        raise ValueError('Nonmixed model must use ordinary_adapter')
    start_ren=np.asarray(start_ren,float)
    if start_ren.shape!=(eft_layout.size+2,):
        raise ValueError('Mixed bilinear+RGBeta vector shape mismatch')
    def unpack(ren,mu):
        running,_,_=mixed._state_at(
            ren,eft_layout,legacy_registry,mu=mu,new=new,reused=reused)
        return running,q_basis.tensor(running)
    def derivative(t,ren):
        return core.rhs_at(t,ren,layout=eft_layout,
                           registry=legacy_registry,betas=eft1_beta,stage='EFT1')
    def context(running,hhss,mu):
        return general.general_initial_context(
            model,running,hhss,source_original,mu_f,mu_s)
    required=tuple(sorted(name for name in eft_layout.names if name.startswith('lambda')))
    return JointEFT1Adapter(model,float(mu_f),float(mu_s),
                            unpack,derivative,context,required),start_ren
