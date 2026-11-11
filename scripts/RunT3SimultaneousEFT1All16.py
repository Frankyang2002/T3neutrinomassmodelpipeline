"""Experimental simultaneous RGBeta + seven-term Wilson RGE for 14 T3 nonmixed models.

B_-1 and C_-1 have separate mixed scalar dynamics and are deliberately
reported as Blocked, NOT computed using an invalid nonmixed approximation.

Intended project path: scripts/RunT3SimultaneousEFT1All16.py

Run from the project root with the earlier seven-term runner installed:

  python -X utf8 -m scripts.RunT3SimultaneousEFT1All16 --self-test
  python -X utf8 -m scripts.RunT3SimultaneousEFT1All16 --source-run thesis16_merged_v1 --run-id simultaneous_nonmixed_v1 --models A_-4 --check-only
  python -X utf8 -m scripts.RunT3SimultaneousEFT1All16 --source-run thesis16_merged_v1 --run-id simultaneous_nonmixed_v1 --models A_-4

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
from itertools import permutations

LOOP = 1.0 / (16.0 * math.pi**2)
OUT_SCHEMA = 't3-general-nonmixed-simultaneous-rgbeta-seven-term-wilson-v1'
MIXED = frozenset(('B_-1','C_-1'))
DIMS = {'A':(1,3),'B':(2,2),'C':(2,2),'D':(3,1),'E':(3,3)}
RTOL = 1.0e-7
ATOL = 1.0e-10


def _complex_number(x):
    if isinstance(x, Mapping) and {'real', 'imag'} <= set(x):
        return complex(float(x['real']), float(x['imag']))
    return complex(x)


def _source_scalar_state(data):
    """Source-quartic evaluator requires dict(real,imag) for complex scalars."""
    out = dict(data)
    for name, raw in data.items():
        if not name.startswith('lambda') or isinstance(raw, (list, tuple)):
            continue
        try:
            z = _complex_number(raw)
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(z.real) and math.isfinite(z.imag)):
            raise ValueError(f'Nonfinite scalar coupling {name}')
        if abs(z.imag) > 0 or isinstance(raw, (complex, np.complexfloating)):
            out[name] = {'real': z.real, 'imag': z.imag}
        else:
            out[name] = float(z.real)
    return out


class HHSSLinearBasis:
    """Same-scale linear HHSS response in the original nonmixed scalar basis."""
    def __init__(self, model, tree, channels, registry, geometry, source_path,
                 starting_state: Mapping):
        if model in MIXED:
            raise ValueError(f'{model}: mixed scalar source needs its own moving-basis construction')
        self.model, self.tree, self.channels = model, tree, channels
        self.registry, self.geometry, self.source_path = registry, geometry, source_path
        n = 2*sum(DIMS[model[0]])
        self.q_shape = (4,4,n,n)
        # The source is linear in scalar quartic couplings. All present lambda
        # fields are probed, including source-registry extras. Missing sources
        # cannot be replaced by zero values.
        self.dependencies = tuple(sorted(k for k,v in starting_state.items()
                                          if k.startswith('lambda') and v is not None))
        if not {'lambdaH1','lambdaH2','lambdaT3'} <= set(self.dependencies):
            raise ValueError(f'{model}: required quartic dependencies missing')
        self.directions=[]
        for key in self.dependencies:
            self.directions.append((key,'real'))
            z = _complex_number(starting_state[key])
            if (key == 'lambdaT3' or abs(z.imag)>0 or isinstance(starting_state[key], (complex, np.complexfloating))
                    or isinstance(starting_state[key], Mapping) and {'real','imag'} <= set(starting_state[key])):
                self.directions.append((key,'imag'))
        self.directions = tuple(self.directions)
        base = dict(starting_state)
        for key in self.dependencies: base[key] = 0.
        self.offset = self._direct(base)
        self.tensors=[]
        for key, component in self.directions:
            probe = dict(base)
            probe[key] = (1j if component=='imag' else 1.)
            self.tensors.append(self._direct(probe)-self.offset)
        self.tensors=np.asarray(self.tensors,float)
        if self.tensors.shape != (len(self.directions),*self.q_shape):
            raise ValueError('Unexpected generalized HHSS linear-basis tensor shape')
        self.assert_reproduces(starting_state,'muF')

    def _direct(self, values):
        from Numerical.running.GeneralT3SequentialLeadingLog import scalar_tensors
        from scripts import RunT3GeneralSequentialFitAblation as original
        snap=original.snapshot(self.model,_source_scalar_state(values),self.source_path)
        bundle=scalar_tensors(self.model,self.tree,self.channels,snap,
                              self.registry,geometry=self.geometry)
        q=np.asarray(bundle.quartic,float)
        if q.shape != self.q_shape or not np.isfinite(q).all():
            raise ValueError(f'{self.model}: invalid direct HHSS tensor')
        return q

    def tensor(self,state):
        coords=np.array([(_complex_number(state[key]).imag if component=='imag'
                          else _complex_number(state[key]).real)
                         for key,component in self.directions],float)
        return self.offset + np.einsum('k,kABab->ABab',coords,self.tensors,optimize=True)

    def assert_reproduces(self,state,label):
        actual=self._direct(state)
        computed=self.tensor(state)
        error=float(np.linalg.norm(computed-actual)/max(np.linalg.norm(actual),1e-300))
        if error>1e-9:
            raise ValueError(f'{self.model} {label}: HHSS linear reconstruction residual={error:.3e}')
        return error


def general_embed(model, source):
    """Map LLSS[flavour,flavour,SU2,SU2,realS,realS] to nine-Weyl tensor."""
    n=2*sum(DIMS[model[0]])
    a=np.asarray(source,complex)
    if a.shape != (3,3,2,2,n,n):
        raise ValueError(f'{model}: unexpected LLSS boundary shape {a.shape}')
    c=np.zeros((9,9,4+n,4+n),complex)
    for p in range(3):
        for q in range(3):
            for i in range(2):
                for j in range(2):
                    c[2*p+i,2*q+j,4:,4:]=a[p,q,i,j]
    if not np.allclose(c,c.swapaxes(0,1),rtol=1e-11,atol=1e-20):
        raise ValueError(f'{model}: fermion pair source symmetry failure')
    if not np.allclose(c,c.swapaxes(2,3),rtol=1e-11,atol=1e-20):
        raise ValueError(f'{model}: scalar pair source symmetry failure')
    return c


def general_gauge(model,state):
    d1,d2=DIMS[model[0]]
    n=4+2*(d1+d2)
    alpha=int(model.split('_')[1]); scalar_blocks=((0,2,.5),(4,d1,alpha/2),
                                                     (4+2*d1,d2,(alpha+2)/2))
    g2=seven._get_real(state,'g2'); gy=seven._get_real(state,'gY')
    scalar_su2=[]
    for axis in range(3):
        m=np.zeros((n,n),complex)
        for begin,d,y in scalar_blocks:
            m[begin:begin+2*d,begin:begin+2*d]=seven.scalar_real_generator(
                seven.su2_complex(d)[axis])
        scalar_su2.append(m)
    sy=np.zeros((n,n),complex)
    for begin,d,y in scalar_blocks:
        sy[begin:begin+2*d,begin:begin+2*d]=seven.scalar_real_generator(y*np.eye(d))
    fs=[]
    for axis in range(3):
        f=np.zeros((9,9),complex)
        for gen in range(3):
            f[2*gen:2*gen+2,2*gen:2*gen+2]=seven.su2_complex(2)[axis]
        fs.append(f)
    fy=np.diag([-.5]*6+[1.]*3).astype(complex)
    return ((g2,tuple(zip(fs,scalar_su2))),(gy,((fy,sy),)))


def general_quartic(model,state,hhss):
    """Source HHSS plus explicit provisional baseline non-HHSS potential."""
    d1,d2=DIMS[model[0]]; nS=2*(d1+d2); n=4+nS
    q=np.asarray(hhss,float)
    if q.shape != (4,4,nS,nS):
        raise ValueError(f'{model}: HHSS shape mismatch {q.shape}')
    lam=np.zeros((n,n,n,n),float)
    for A in range(4):
        for B in range(4):
            for a in range(nS):
                for b in range(nS):
                    value=q[A,B,a,b]
                    if value:
                        for ix in set(permutations((A,B,4+a,4+b))):lam[ix]=value
    h=np.diag([1.]*4+[0.]*nS)
    s1=np.diag([0.]*(4)+[1.]*(2*d1)+[0.]*(2*d2))
    s2=np.diag([0.]*(4+2*d1)+[1.]*(2*d2))
    for name,p in (('lambdaH',h),('lambdaS1',s1),('lambdaS2',s2)):
        seven.add_pair_contraction(lam,p,p,seven._get_real(state,name))
    seven.add_pair_contraction(lam,s1,s2,seven._get_real(state,'lambda12'))
    for name,d,start in (('lambdaS1Adj',d1,4),('lambdaS2Adj',d2,4+2*d1)):
        if name not in state or state[name] is None:
            continue
        coeff=seven._get_real(state,name)
        for t in seven.su2_complex(d):
            adj=np.zeros((n,n))
            adj[start:start+2*d,start:start+2*d]=seven.scalar_bilinear_real(t)
            seven.add_pair_contraction(lam,adj,adj,coeff)
    if not np.isfinite(lam).all() or not np.allclose(lam[:4,:4,4:,4:],q):
        raise ValueError(f'{model}: invalid non-HHSS ansatz/HHSS overwrite')
    if not all(np.allclose(lam,lam.swapaxes(0,i),rtol=1e-10,atol=1e-12)
               for i in (1,2,3)):
        raise ValueError(f'{model}: quartic Bose symmetry failure')
    return lam


def general_initial_context(model,state,hhss,source,mu_f,mu_s):
    from RGE.running.intermediate.GeneralEFT1C5Mixing import Context
    n=4+2*sum(DIMS[model[0]])
    yfull=seven.yukawa_context(state)
    # The uploaded A(-4) Yukawa function uses globally fixed N_S=12.
    # Its only nonzero entries are in the four Higgs components: zero-pad.
    y=np.zeros((45,45,n),complex)
    y[:,:,:4]=yfull[:,:,:4]
    qgamma=seven.quark_higgs_anomalous_dimension(y)
    ctx=Context(model=model,nf=9,ns=n,leptons=((0,1),(2,3),(4,5)),
                higgs=(0,1,2,3),mu_f=mu_f,mu_s=mu_s,
                yukawa=y[:9,:9,:],quartic=general_quartic(model,state,hhss),
                gauge=general_gauge(model,state),
                operators=(('source_resolved_fermion_tree_LLSS',general_embed(model,source)),),
                authority='synthetic',certificates={})
    return ctx,qgamma


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
    c0 = general_embed(q_basis.model,source_c)
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
    master0, gamma0 = general_initial_context(
        q_basis.model,data0, q0, source_c, eft_start.mu_gev, mu_s,
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
        context, quark_gamma = general_initial_context(
            q_basis.model,data, hhss, source_c, original.SCALES[1], mu_s,
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
    # Benchmark against the frozen-coupling seven-term implementation
    # with the same generalised candidate context and endpoint.
    ctx,quark_gamma=general_initial_context(
        model,state,np.asarray(scalar.quartic,float),
        np.asarray(tree.raw,complex)*complex(spinor_factor),mu_f,mu_s)
    frozen=seven.solve_master(ctx,quark_gamma=quark_gamma)
    c5_frozen,ward_frozen=seven.project_leptonic_weinberg(frozen['end'])
    if model == 'A_-4':
        quartic=integrate_coupled_quartic(
            np.asarray(tree.raw,complex)*complex(spinor_factor),
            np.asarray(scalar.quartic,float),mu_f,mu_s)
        c5_quartic,quartic_ward=project_weinberg(quartic['c_llhh_end'], strict=False)
    else:
        c5_quartic=np.zeros_like(c5_super)
        quartic_ward=0.
    norm=lambda x: float(np.linalg.norm(x))
    rel=lambda a,b: norm(a-b)/max(norm(b),1e-300)
    metrics={
        'simultaneous_vs_frozen_7term_rel':rel(c5_super,c5_frozen),
        'simultaneous_vs_fixed_log_rel':rel(c5_super,fixed_log),
        'simultaneous_vs_quartic_only_rel':rel(c5_super,c5_quartic) if model=='A_-4' else None,
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
                'c5_coupled_quartic_only':encode(c5_quartic) if model=='A_-4' else None,
                'c5_joint_scalar_bubble':encode(finite),
                'c5_matching_anchored_reference':encode(original.c5_checked(c5_one)),
                'c5_bookkeeping_residual_matching':encode(np.zeros((3,3),complex)),
                'c5_log_frobenius':norm(c5_super),
                'c5_one_scale_frobenius':norm(c5_one),
                'c5_log_over_one_scale_norm':norm(c5_super)/max(norm(c5_super+finite),1e-300),
                'matching_approximation': 'simultaneous_RGBeta_seven_term_9Weyl_plus_finite_bubble',
                'full_physical_one_loop_matching':False,
                'renormalisable_couplings_frozen_to_muS':not eft_on,
                'wilson_fermion_subspace': 'nine lepton Weyl fields; quark Higgs gamma included',
                 'non_HHSS_quartic_ansatz_provisional':True,
                'per_general_interaction':fixed['per_channel'],
                'quartic_provenance':fixed['quartic_provenance'],
                'matching_scale_cancellation':fixed['matching_scale_cancellation_residual'],
                **metrics,
            })
            cases[code]=prediction
    print(f'{model} U={int(uv_on)} E={int(eft_on)} : '
          f'sim/frozen7={metrics["simultaneous_vs_frozen_7term_rel"]:.6e} '
          f'sim/fixedLL={metrics["simultaneous_vs_fixed_log_rel"]:.6e} '
          f'sim/quarticODE={metrics["simultaneous_vs_quartic_only_rel"]} '
          f'joint nfev={simultaneous["nfev"]}',flush=True)
    return cases


def _nonmixed_general(model, saved, canonical_run_id, quartic_source, registry,
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

    if model in MIXED: raise ValueError(f'{model}: mixed scalar RGE requires a distinct coupled state and matching basis')
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
    from scripts import RunT3GeneralSequentialFitAblation as original
    original._nonmixed=_nonmixed_general
    def reject_mixed(model,*args,**kwargs):
        raise NotImplementedError(
            f"{model}: mixed B/C state requires coupled moving scalar basis and "
            "muF-source matching; deliberately not marked Computed")
    original._mixed=reject_mixed
    original.OUT=(original.ROOT/'output/studies/general_scalar_extended/'
                  'experimental_simultaneous_eft1_general_nonmixed')
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
                data['model_scope']='14_nonmixed_models; B_-1_C_-1_blocked'
            elif data.get('status')=='Computed' and 'cases' in data:
                data['matching_approximation']='experimental_simultaneous_EFT1_nonmixed'
                data['full_sequential_general_C5_certified']=False
        old_save(path,data)
    original.save=save_with_metadata
    return original.main(args)


if __name__=='__main__':
    raise SystemExit(main())
