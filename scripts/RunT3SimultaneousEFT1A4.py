"""Experimental simultaneous RGBeta + seven-term Wilson RGE for T3 A(-4).

Intended project path: scripts/RunT3SimultaneousEFT1A4.py

Run from the project root with the earlier seven-term runner installed:

  python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --self-test
  python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --source-run thesis16_merged_v1 --run-id simultaneous_A_m4_v1 --models A_-4 --check-only
  python -X utf8 -m scripts.RunT3SimultaneousEFT1A4 --source-run thesis16_merged_v1 --run-id simultaneous_A_m4_v1 --models A_-4

The joint DOP853 state is [RGBeta renormalisable variables; Re C5_master;
Im C5_master]. At every RHS call, evaluate RGBeta and all seven terms of
RGE.running.intermediate.GeneralEFT1C5Mixing.calculate_master_beta using
THAT SAME instantaneous renormalisable state.

The A(-4) HHSS tensor is assembled from a once-evaluated linear coupling
basis, checked against scalar_tensors at BOTH matching scales. Non-HHSS
scalar quartics retain the previous seven-term runner's explicit potential
ansatz, evaluated at the running scale. We retain that runner's candidate
left-Weyl SM Yukawa embedding and quark-Higgs anomalous dimension.

This is a genuine simultaneous ODE test for the RETAINED 9-Weyl-field
psi^2 phi^2 basis, not complete verified physical matching: hard one-loop
fermion matching and other operator classes are absent; 7-term tensor and
scalar normalisations are inherited/provisional. Both scalars are matched
jointly at muS using the existing finite bubble, and SMEFT evolves separately.
No original files/runs are modified and no GitHub writes occur.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp

from scripts import RunT3CoupledEFT1SevenTerm as seven

LOOP = 1.0 / (16.0 * math.pi**2)
OUT_SCHEMA = 't3-a-m4-simultaneous-rgbeta-seven-term-wilson-v1'
RTOL = 1.0e-7
ATOL = 1.0e-10


def _complex_number(x):
    if isinstance(x, Mapping) and {'real', 'imag'} <= set(x):
        return complex(float(x['real']), float(x['imag']))
    return complex(x)


def _source_scalar_snapshot(values):
    """Encode the complex T3 coupling in the source-quartic input format.

    The HHSS linear-basis probes use 0j/1j, whereas source_quartic_v4
    accepts real numbers or {'real': ..., 'imag': ...}, not Python complex.
    Keep other EFT1 couplings and the caller's dictionary unchanged.
    """
    out = dict(values)
    if 'lambdaT3' in out and out['lambdaT3'] is not None:
        z = _complex_number(out['lambdaT3'])
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError('lambdaT3: nonfinite source coupling')
        out['lambdaT3'] = {'real': float(z.real), 'imag': float(z.imag)}
    return out


class HHSSLinearBasis:
    """Compute HHSS dependence on running A(-4) portal and T3 couplings once."""

    def __init__(self, model, tree, channels, registry, geometry, source_path,
                 starting_state: Mapping):
        if model != 'A_-4':
            raise ValueError('HHSSLinearBasis supports A_-4 only')
        if registry['models'][model]['total_additional_quartics'] != 0:
            raise ValueError('A(-4) unexpectedly has additional candidate quartics')
        from scripts import RunT3GeneralSequentialFitAblation as original
        self.model = model
        self.tree = tree
        self.channels = channels
        self.registry = registry
        self.geometry = geometry
        self.source_path = source_path
        self.q_shape = (4, 4, 8, 8)
        self.dependencies = tuple(k for k in (
            'lambdaH1', 'lambdaH2', 'lambdaH1Adj', 'lambdaH2Adj', 'lambdaT3'
        ) if k in starting_state and starting_state[k] is not None)
        if not all(k in self.dependencies for k in ('lambdaH1','lambdaH2','lambdaT3')):
            raise ValueError('Expected A(-4) portal/T3 couplings are absent')
        self.directions = tuple(
            (key, component)
            for key in self.dependencies
            for component in (('real','imag') if key == 'lambdaT3' else ('real',))
        )
        base = dict(starting_state)
        for key in self.dependencies:
            base[key] = 0j if key == 'lambdaT3' else 0.
        self.offset = self._direct(base)
        self.tensors = []
        for key, component in self.directions:
            probe = dict(base)
            probe[key] = 1j if component == 'imag' else 1.
            self.tensors.append(self._direct(probe) - self.offset)
        self.tensors = np.asarray(self.tensors, dtype=float)
        if self.tensors.shape != (len(self.directions), *self.q_shape):
            raise ValueError('Unexpected HHSS basis shape')
        self.assert_reproduces(starting_state, 'muF')

    def _direct(self, values):
        from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
        from scripts import RunT3GeneralSequentialFitAblation as original
        snap = original.snapshot(self.model, _source_scalar_snapshot(values),
                                 self.source_path)
        bundle = scalar_tensors(self.model, self.tree, self.channels, snap,
                                self.registry, geometry=self.geometry)
        q = np.asarray(bundle.quartic, dtype=float)
        if q.shape != self.q_shape or not np.isfinite(q).all():
            raise ValueError('A(-4) direct HHSS source gave invalid tensor')
        return q

    def tensor(self, state):
        values = np.array([
            (_complex_number(state[key]).imag if component == 'imag'
             else _complex_number(state[key]).real)
            for key, component in self.directions
        ], dtype=float)
        return self.offset + np.einsum('k,kABab->ABab', values,
                                       self.tensors, optimize=True)

    def assert_reproduces(self, state, label):
        actual = self._direct(state)
        reconstructed = self.tensor(state)
        deviation = float(np.linalg.norm(reconstructed - actual) / max(
            np.linalg.norm(actual), 1e-300))
        if deviation > 1e-9:
            raise ValueError(f'{label}: linear HHSS vs source tensor mismatch '
                             f'{deviation:.3e}')
        return deviation


def _flatten_master(c):
    z = np.asarray(c, complex)
    return np.concatenate((z.real.ravel(), z.imag.ravel()))


def _unflatten_master(real_flat, shape):
    n = int(np.prod(shape))
    v = np.asarray(real_flat, float)
    if v.shape != (2*n,):
        raise ValueError('Invalid complex master tensor vector size')
    return v[:n].reshape(shape) + 1j*v[n:].reshape(shape)


def integrate_simultaneously(eft_start, beta_eft, source_c,
                             q_basis, mu_s, *, eft_on=True,
                             rtol=RTOL, atol=ATOL):
    """Integrate [renormalisable EFT1, full lepton-Weyl master tensor] in ONE IVP.

    Uses the exact existing RGBeta derivative evaluator on renormalisable state,
    while feeding running renormalisable couplings into the full seven-term
    Wilson RHS at every internal DOP853 point.
    """
    from Numerical.running.IntermediateScalarRunner import _rhs_factory
    from Numerical.running.IntermediateScalarStateVector import (
        pack_intermediate_scalar_state, unpack_intermediate_scalar_state,
    )
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        Context, TERMS, calculate_master_beta,
    )
    from scripts import RunT3GeneralSequentialFitAblation as original

    ren0, layout = pack_intermediate_scalar_state(eft_start)
    if set(layout.names) != set(beta_eft.get('report_betas',{})):
        raise ValueError('EFT1 RGBeta inventory differs from packed state')
    c0 = seven.embed_llss(source_c)
    shape = c0.shape
    cscale = float(np.linalg.norm(c0))
    if not math.isfinite(cscale) or cscale <= 0:
        raise ValueError('Zero or nonfinite initial LLSS tensor')
    nren, nc = ren0.size, c0.size
    z0 = np.concatenate((ren0, _flatten_master(c0/cscale)))
    beta_ren = _rhs_factory(beta_eft, layout)
    # Mandatory A(-4) calibration: the initial master scalar-pair term must
    # reproduce the prior LLSS -> Weinberg scalar-quartic source exactly.
    data0 = _state_values(eft_start)
    q0 = q_basis.tensor(data0)
    master0, gamma0 = seven.initial_context(
        data0, q0, source_c, eft_start.mu_gev, mu_s,
    )
    b0 = calculate_master_beta(master0, c0)['scalar_pair']
    c5_master0, _ = seven.project_leptonic_weinberg(b0)
    from Numerical.joint_scalar_matching.core import project_weinberg
    direct0 = np.einsum('fgijab,ABab->fgijAB', source_c, q0,
                        optimize=True)
    c5_direct0, _ = project_weinberg(direct0, strict=False)
    source_parity = float(np.linalg.norm(c5_master0-c5_direct0) /
                          max(np.linalg.norm(c5_direct0),1e-300))
    if source_parity > 1e-8:
        raise ValueError('muF master scalar-pair beta does not match LLSS '
                         f'HHSS source: relative {source_parity:.3e}')
    stats = {'rhs_calls':0}

    def rhs(t, z):
        ren = z[:nren]
        master = _unflatten_master(z[nren:], shape)
        at_t = unpack_intermediate_scalar_state(ren, layout,
                                                mu_gev=float(np.exp(t)))
        data = _state_values(at_t)
        hhss = q_basis.tensor(data)
        # Reuse the current seven-term RGE conventions, now evaluated at t.
        # A(-4) is a no-mixing original scalar basis, hence HHSS directions
        # remain well-defined at all intermediate renormalisation scales.
        context, quark_gamma = seven.initial_context(
            data, hhss, source_c, original.SCALES[1], mu_s,
        )
        terms = calculate_master_beta(context, master)
        gam = np.einsum('bd,ijad->ijab', quark_gamma, master,
                        optimize=True)
        terms['scalar_anomalous_dimension'] += gam + gam.swapaxes(2,3)
        dm = LOOP * sum((terms[key] for key in TERMS),
                         np.zeros(shape, complex))
        dr = beta_ren(t, ren) if eft_on else np.zeros_like(ren)
        stats['rhs_calls'] += 1
        return np.concatenate((dr, _flatten_master(dm)))

    t0 = math.log(eft_start.mu_gev)
    t1 = math.log(float(mu_s))
    sol = solve_ivp(rhs, (t0,t1), z0, method='DOP853',
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError('Combined EFT1 DOP853 failed: '+sol.message)
    end_ren = sol.y[:nren,-1].copy()
    end_c = cscale * _unflatten_master(sol.y[nren:,-1], shape)
    if not (np.isfinite(end_ren).all() and np.isfinite(end_c).all()):
        raise ValueError('Combined EFT1 solver produced nonfinite endpoint')
    end_state = unpack_intermediate_scalar_state(end_ren,layout,
                                                 mu_gev=float(mu_s))
    # Diagnostic tensor symmetry. No silent symmetrisation of the result.
    symf = float(np.linalg.norm(end_c - end_c.swapaxes(0,1)) / cscale)
    syms = float(np.linalg.norm(end_c - end_c.swapaxes(2,3)) / cscale)
    if max(symf,syms) > 1e-6:
        raise ValueError('Combined EFT1 end tensor violates pair symmetry')
    return {
        'end_state':end_state,'end_c':end_c,'end_ren_vector':end_ren,
        'ren_layout':layout,'nfev':int(sol.nfev),
        'steps':int(sol.t.size),'cscale':cscale,
        'tensor_drift':float(np.linalg.norm(end_c-c0)/cscale),
        'fermion_symmetry':symf,'scalar_symmetry':syms,
        'rhs_calls':stats['rhs_calls'],
        'ren_beta_active':bool(eft_on),
        'initial_scalar_quartic_source_parity_rel':source_parity,
    }


def _scalar_endpoint_prediction(model, c5_one, tree, scalar, target, *,
                                uv_on, eft_on, uv_nfev, eft_nfev,
                                spinor_factor, vev, simultaneous):
    """Use the SAME scalar finite matching and SMEFT as prior experiments."""
    from Numerical.running.GeneralT3SequentialLeadingLog import (
        encode, decode, response,
    )
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        TERMS, calculate_master_beta,
    )
    from Numerical.joint_scalar_matching.core import project_weinberg
    from scripts.RunT3CoupledEFT1Quartic import integrate_coupled_quartic
    from scripts import RunT3GeneralSequentialFitAblation as original

    mu_f, mu_s = original.SCALES[1:3]
    state = original.scalar_snapshot_state(scalar)
    fixed = response(tree,scalar,mu_f=mu_f,mu_s=mu_s,
                     spinor_factor=complex(spinor_factor),
                     check_scale_cancellation=True)
    fixed_log=decode(fixed['c5_eft1_leading_log'])
    finite=decode(fixed['c5_scalar_finite_bubble'])
    c5_super,ward=seven.project_leptonic_weinberg(simultaneous['end_c'])
    # Benchmark against the previously successful frozen-coupling seven-term
    # implementation, using the SAME renormalisable endpoint.
    ctx,quark_gamma=seven.initial_context(
        state,np.asarray(scalar.quartic,float),
        np.asarray(tree.raw,complex)*complex(spinor_factor),mu_f,mu_s)
    frozen=seven.solve_master(ctx,quark_gamma=quark_gamma)
    c5_frozen,ward_frozen=seven.project_leptonic_weinberg(frozen['end'])
    quartic=integrate_coupled_quartic(
        np.asarray(tree.raw,complex)*complex(spinor_factor),
        np.asarray(scalar.quartic,float),mu_f,mu_s)
    c5_quartic,quartic_ward=project_weinberg(quartic['c_llhh_end'], strict=False)
    norm=lambda x: float(np.linalg.norm(x))
    rel=lambda a,b: norm(a-b)/max(norm(b),1e-300)
    metrics={
        'simultaneous_vs_frozen_7term_rel':rel(c5_super,c5_frozen),
        'simultaneous_vs_fixed_log_rel':rel(c5_super,fixed_log),
        'simultaneous_vs_quartic_only_rel':rel(c5_super,c5_quartic),
        'frozen_7term_vs_fixed_log_rel':rel(c5_frozen,fixed_log),
        'coupled_llss_master_tensor_drift_rel':simultaneous['tensor_drift'],
        'coupled_rhs_calls':simultaneous['rhs_calls'],
        'combined_ode_nfev':simultaneous['nfev'],
        'initial_scalar_quartic_source_parity_rel':
            simultaneous['initial_scalar_quartic_source_parity_rel'],
        'combined_ode_steps':simultaneous['steps'],
        'combined_ode_fermion_symmetry_residual':simultaneous['fermion_symmetry'],
        'combined_ode_scalar_symmetry_residual':simultaneous['scalar_symmetry'],
        'combined_ode_weinberg_projection_residual':float(ward),
        'frozen_ode_weinberg_projection_residual':float(ward_frozen),
        'quartic_ode_weinberg_projection_residual':float(quartic_ward),
    }
    if not eft_on and metrics['simultaneous_vs_frozen_7term_rel'] > 2e-6:
        raise ValueError('EFT1 renormalisable OFF must reproduce frozen-muS '
                         'seven-term solution; discrepancy=' +
                         str(metrics['simultaneous_vs_frozen_7term_rel']))
    cases={}
    for w_on, c5 in ((False,finite),(True,finite+c5_super)):
        c5=original.c5_checked(c5)
        for sm_on in (False,True):
            code=f'{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}'
            prediction=original.low_energy(state,c5,on=sm_on,mu_s=mu_s,
                                           mu_low=original.SCALES[3],
                                           vev=vev,target=target)
            prediction.update({
                'case':code,'uv_on':uv_on,'eft1_on':eft_on,
                'wilson_log_on':w_on,'smeft_on':sm_on,
                'uv_nfev':uv_nfev,'eft1_nfev':eft_nfev,
                'c5_at_scalar_threshold':encode(c5),
                'c5_induced_wilson_log':encode(c5_super),
                'c5_simultaneous_running':encode(c5_super),
                'c5_frozen_7term':encode(c5_frozen),
                'c5_fixed_order_llss_log':encode(fixed_log),
                'c5_coupled_quartic_only':encode(c5_quartic),
                'c5_joint_scalar_bubble':encode(finite),
                'c5_matching_anchored_reference':encode(original.c5_checked(c5_one)),
                'c5_bookkeeping_residual_matching':encode(np.zeros((3,3),complex)),
                'c5_log_frobenius':norm(c5_super),
                'c5_one_scale_frobenius':norm(c5_one),
                'c5_log_over_one_scale_norm':norm(c5_super)/max(norm(c5_super+finite),1e-300),
                'matching_approximation': 'simultaneous_RGBeta_seven_term_9Weyl_plus_finite_bubble',
                'full_physical_one_loop_matching':False,
                'renormalisable_couplings_frozen_to_muS':False,
                'wilson_fermion_subspace': 'nine lepton Weyl fields; quark Higgs gamma included',
                'per_general_interaction':fixed['per_channel'],
                'quartic_provenance':fixed['quartic_provenance'],
                'matching_scale_cancellation':fixed['matching_scale_cancellation_residual'],
                **metrics,
            })
            cases[code]=prediction
    print(f'{model} U={int(uv_on)} E={int(eft_on)} : '
          f'sim/frozen7={metrics["simultaneous_vs_frozen_7term_rel"]:.6e} '
          f'sim/fixedLL={metrics["simultaneous_vs_fixed_log_rel"]:.6e} '
          f'sim/quarticODE={metrics["simultaneous_vs_quartic_only_rel"]:.6e} '
          f'joint nfev={simultaneous["nfev"]}',flush=True)
    return cases


def _nonmixed_A4(model, saved, canonical_run_id, quartic_source, registry,
                 *, spinor_factor, wolframscript, vev):
    """Retain original UV/F and scalar/SMEFT steps, replace ONLY EFT1 stage."""
    from scripts import RunT3GeneralAll16FromReference as g
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state,_state_values
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry

    if model!='A_-4':raise ValueError('Only A(-4) is implemented')
    muF,muS=original.SCALES[1:3]
    beta_uv,beta_eft,_,kernel,matchete_path,kind=g.nonmixed_inputs(
        model,canonical_run_id,generate_general_c5=False,wolframscript=wolframscript)
    overrides=original.load(saved/'general_uv_overrides.json')
    _,initial=initial_state(beta_uv,mu_uv=original.SCALES[0],
                             mu_f=muF,mu_s=muS,overrides=overrides)
    static_source,channels=quartic_source
    target=OscillationFitTarget.from_json(TARGET)
    geometry=model_extra_geometry(model,registry['models'][model],include_components=True)
    cases={};sources={};endpoint_comparisons={}
    for uv_on in (False,True):
        if uv_on:
            uv_sol=run_uv_segment(initial,beta_uv,muF,rtol=1e-7,atol=1e-10)
            uv_end,uv_nfev=uv_sol.final_state,int(uv_sol.nfev)
        else:
            uv_end,uv_nfev=original.frozen_state(initial,muF),0
        tree,tree_diagnostics=original.nonmixed_boundary(model,uv_end,
                                                         repo_root=original.ROOT,mu_f=muF)
        sources[str(int(uv_on))]=tree_diagnostics
        eft_start=build_intermediate_scalar_boundary(uv_end,matching_scale_gev=muF)
        start_data=_state_values(eft_start)
        q_basis=HHSSLinearBasis(model,tree,channels,registry,geometry,
                                saved/'fitted_config.json',start_data)
        for eft_on in (False,True):
            # Actual joint RGBeta+seven-term solve on the complete retained state.
            combined=integrate_simultaneously(
                eft_start,beta_eft,
                np.asarray(tree.raw,complex)*complex(spinor_factor),
                q_basis,muS,eft_on=eft_on)
            eft_end=combined['end_state']
            end_data=_state_values(eft_end)
            endpoint_reconstruction=q_basis.assert_reproduces(end_data,'muS')
            if eft_on:
                # Compare with the original RGBeta-only trajectory to ensure
                # the enlarged state has not changed renormalisable running.
                separate=run_intermediate_scalar_segment(
                    eft_start,beta_eft,muS,rtol=1e-7,atol=1e-10)
                ref_end=separate.final_vector
                ref_nfev=int(separate.nfev)
                breakdown={
                    block.name:float(np.linalg.norm(
                        (combined['end_ren_vector']-ref_end)[block.start:block.stop]
                    )/max(np.linalg.norm(ref_end[block.start:block.stop]),1e-12))
                    for block in combined['ren_layout'].blocks
                }
                max_endpoint=max(breakdown.values())
                if max_endpoint>5e-5:
                    raise ValueError('Combined and original RGBeta EFT1 endpoints '
                                     f'disagree: worst relative difference={max_endpoint:.3e} '
                                     f'(block {max(breakdown,key=breakdown.get)})')
            else:
                ref_nfev=0
                max_endpoint=0.
            snap=original.snapshot(model,end_data,saved/'fitted_config.json')
            tensor=original.scalar_tensors(model,tree,channels,snap,registry,geometry=geometry)
            scalar=original.ScalarEndpoint(tensor,end_data)
            one_scale=kernel.evaluate(uv_end,eft_end)
            block=_scalar_endpoint_prediction(
                model,one_scale,tree,scalar,target,
                uv_on=uv_on,eft_on=eft_on,
                uv_nfev=uv_nfev,eft_nfev=combined['nfev'],
                spinor_factor=spinor_factor,vev=vev,
                simultaneous=combined)
            for record in block.values():
                record['hhss_linear_basis_parity_rel']=endpoint_reconstruction
                record['separate_rgbeta_endpoint_max_block_rel']=max_endpoint
                record['separate_rgbeta_nfev']=ref_nfev
            cases.update(block)
            endpoint_comparisons[f'{int(uv_on)}{int(eft_on)}']={
                'hhss_basis_parity_rel':endpoint_reconstruction,
                'renormalisable_endpoint_max_block_rel':max_endpoint,
                'combined_nfev':combined['nfev'],
                'separate_nfev':ref_nfev,
            }
    return cases,{
        'fermion_tree':sources,
        'matching_kind':kind,
        'full_general_C5_matchete_path':str(matchete_path),
        'full_general_C5_sha256':original.hash_file(matchete_path),
        'fitted_uv_sha256':original.hash_file(saved/'general_uv_overrides.json'),
        'hhss_linear_basis_directions':list(q_basis.directions),
        'endpoint_comparisons':endpoint_comparisons,
        'simultaneous_rhs':True,
    }


def self_test():
    # Check the source-quartic format used by HHSSLinearBasis offset/probes.
    for value, expected in ((0j, (0., 0.)), (1j, (0., 1.)),
                            (2.-3.j, (2., -3.)),
                            ({'real': 2., 'imag': -3.}, (2., -3.))):
        original = {'lambdaT3': value, 'lambdaH1': 0.125}
        encoded = _source_scalar_snapshot(original)
        if (encoded['lambdaT3'] != {'real': expected[0], 'imag': expected[1]}
                or encoded['lambdaH1'] != 0.125
                or original['lambdaT3'] != value):
            raise AssertionError('HHSS source scalar encoding roundtrip failed')
    print('PASS: lambdaT3 source snapshot scalar encoding')
    # Prior 7-term tensor engine's mathematical and complex-data checks.
    seven.self_test()
    rng=np.random.default_rng(9)
    a=rng.standard_normal((2,3,4))+1j*rng.standard_normal((2,3,4))
    if not np.array_equal(_unflatten_master(_flatten_master(a),a.shape),a):
        raise AssertionError('Wilson complex packing roundtrip failed')
    # Analytic simultaneous running toy: dg/dt=a*g and
    # d(CS,CH)/dt=LOOP*q(t)*(CH,CS), q(t)=q0*g(t).
    g0=0.3; slope=0.04; q0=0.5; dlog=math.log(1e3/1e5)
    def rhs(_t,v):
        g,cs,ch=v
        return [slope*g,LOOP*q0*g*ch,LOOP*q0*g*cs]
    sol=solve_ivp(rhs,(0,dlog),[g0,1.,0.],method='DOP853',
                  rtol=1e-11,atol=1e-13)
    x=LOOP*q0*g0*math.expm1(slope*dlog)/slope
    if not sol.success or not np.allclose(
            sol.y[:,-1],[g0*math.exp(slope*dlog),math.cosh(x),math.sinh(x)],
            rtol=1e-9,atol=1e-12):
        raise AssertionError('Joint toy renormalisable+Wilson ODE failed')
    print('PASS: simultaneous renormalisable + Wilson DOP853 analytic toy')
    print(f'  toy evolving-coupling LLHH={sol.y[2,-1]:.12e}')


def main(argv=None):
    args=list(sys.argv[1:] if argv is None else argv)
    if args==['--self-test']:
        self_test();return 0
    if '--run-id' not in args:
        raise ValueError('Provide a fresh --run-id')
    if '--resume' in args:
        raise ValueError('Resume disabled; use a fresh --run-id')
    if '--models' not in args or args[args.index('--models')+1]!='A_-4':
        raise ValueError('This simultaneous trial requires --models A_-4')
    from scripts import RunT3GeneralSequentialFitAblation as original
    original._nonmixed=_nonmixed_A4
    original.OUT=(original.ROOT/'output/studies/general_scalar_extended/'
                  'experimental_simultaneous_eft1_seven_term')
    old_save=original.save
    def save_with_metadata(path,data):
        if isinstance(data,dict):
            data=dict(data)
            if data.get('format')=='t3-full-general-16_U_E_W_S_ablation-v1':
                data['format']=OUT_SCHEMA
                data['running_background']='simultaneous_RGBeta_and_7term_master'
                data['one_scale_C5_used_as_W_on_reference']=False
                data['hard_matching_residual_independently_calculated']=False
                data['complete_physical_C5']=False
                data['model_scope']='A_-4'
            elif data.get('status')=='Computed' and 'cases' in data:
                data['matching_approximation']='experimental_simultaneous_EFT1_A4'
                data['full_sequential_general_C5_certified']=False
        old_save(path,data)
    original.save=save_with_metadata
    return original.main(args)


if __name__=='__main__':
    raise SystemExit(main())
