"""Experimental simultaneous EFT1 seven-term Wilson + mixed RGBeta trial: B_-1,C_-1.

Intended path: scripts/RunT3SimultaneousEFT1MixedBC.py

The F-threshold LLSS tensor and all running HHSS tensors are represented in
ONE FIXED original Cartesian S1,S2 basis. In particular the mu-dependent X/Y
mass eigenbasis is NOT treated as a fixed basis when muT3GeneralS1S2 runs.
Both scalars are integrated out together at muS with the existing conditional
B/C scalar matching; SMEFT running is unchanged.

Physical limitations: retained nine lepton Weyl fields only, provisional full
non-HHSS potential and Yukawa dictionary; no independent hard F-loop C5.
This is not a certified complete dimension-five EFT operator calculation.
No GitHub or existing source files are modified.

From repository root:
  python -X utf8 -m scripts.RunT3SimultaneousEFT1MixedBC --self-test
  python -X utf8 -m scripts.RunT3SimultaneousEFT1MixedBC --source-run thesis16_merged_v1 --run-id mixed_BC_v1 --models B_-1,C_-1 --check-only
  python -X utf8 -m scripts.RunT3SimultaneousEFT1MixedBC --source-run thesis16_merged_v1 --run-id mixed_BC_v1 --models B_-1,C_-1
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from scripts import RunT3SimultaneousEFT1All16 as general

LOOP = 1.0 / (16 * math.pi**2)
MIXED = ('B_-1', 'C_-1')
SCHEMA = 't3-mixed-bc-simultaneous-rgbeta-seven-term-wilson-fixed-cartesian-v1'
RTOL, ATOL = 1e-8, 1e-11  # match the existing native mixed RGBeta solver


def _cartesian_source(source, R):
    """C_original[i,j] = R[i,a] R[j,b] C_XY[a,b]."""
    a = np.asarray(source, complex)
    r = np.asarray(R, float)
    if a.shape != (3, 3, 2, 2, 8, 8) or r.shape != (8, 8):
        raise ValueError('Mixed LLSS Cartesian conversion shape mismatch')
    if not np.allclose(r.T @ r, np.eye(8), rtol=1e-10, atol=1e-10):
        raise ValueError('Nonorthogonal mixed scalar rotation')
    return np.einsum('fgpqab,ia,jb->fgpqij', a, r, r, optimize=True)


def _cartesian_hhss(quartic, R):
    q = np.asarray(quartic, float)
    r = np.asarray(R, float)
    if q.shape != (4, 4, 8, 8) or r.shape != (8, 8):
        raise ValueError('Mixed HHSS Cartesian conversion shape mismatch')
    return np.einsum('ABab,ia,jb->ABij', q, r, r, optimize=True)


def _relative(a, b):
    return float(np.linalg.norm(a - b) / max(float(np.linalg.norm(b)), 1e-300))


def _state_at(vector, layout, registry, *, mu, new, reused):
    """Decode full extended mixed RGBeta state, including mu2 and extras."""
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from Numerical.experimental.GeneralT3ScalarRGState import unpack_general_scalar_state
    from Numerical.orchestration.GeneralRunningBenchmark import _state_values

    v = np.asarray(vector, float)
    if v.shape != (layout.size + 2,) or not np.isfinite(v).all():
        raise ValueError('Mixed EFT1 vector must include two real mu2 components')
    extended = unpack_general_scalar_state(v[:layout.size], layout, registry,
                                           uv=False, mu_gev=float(mu))
    values = _state_values(extended.baseline)
    report = core.snapshot(v, layout, registry, stage='EFT1', mu_gev=float(mu),
                           new=new, reused=reused)
    for key, value in report['new_quartics'].items():
        values[key] = value
    for key, value in report['reused_quartics'].items():
        values[key] = value
    values['muT3GeneralS1S2'] = report['quadratic_mu2_gev2']
    return values, extended.baseline, report


def _basis_and_tensor(model, tree_f, channels, registry, geometry,
                      source_path, data):
    """Use source exact XY->original rotation and return ORIGINAL-basis HHSS."""
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
    snap = original.snapshot(model, general._source_scalar_state(data), source_path)
    basis = build_special_basis(model, snap, tree_f)
    R = np.asarray(basis['R_original_from_boundary'], float)
    bundle = original.scalar_tensors(model, tree_f, channels, snap, registry,
                                     geometry=geometry)
    # Preserve EXACTLY the established conditional B/C channel conventions in
    # the Wilson anomalous dimension as well as in muS finite matching.
    from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
        apply_conditional_bc_direct_quartic_map,
    )
    from Numerical.running.GeneralT3ConditionalBCQuarticNormalization import (
        apply_conditional_bc_channel_normalization,
    )
    bundle = apply_conditional_bc_direct_quartic_map(model, bundle)
    bundle = apply_conditional_bc_channel_normalization(model, bundle)
    q = _cartesian_hhss(bundle.quartic, R)
    return R, q, bundle


class FixedCartesianHHSS:
    """Linear quartic response in a fixed original basis, not moving X/Y."""
    def __init__(self, model, tree_f, channels, registry, geometry,
                 source_path, initial_data):
        self.model, self.tree_f, self.channels = model, tree_f, channels
        self.registry, self.geometry, self.source_path = registry, geometry, source_path
        self.dependencies = tuple(sorted(k for k, v in initial_data.items()
                                         if k.startswith('lambda') and v is not None))
        if not {'lambdaH1', 'lambdaH2', 'lambdaT3'} <= set(self.dependencies):
            raise ValueError('Mixed HHSS required portal directions absent')
        directions = []
        for name in self.dependencies:
            directions.append((name, 'real'))
            raw = initial_data[name]
            value = general._complex_number(raw)
            if (name == 'lambdaT3' or abs(value.imag) > 0
                    or isinstance(raw, (complex, np.complexfloating))
                    or isinstance(raw, dict) and {'real', 'imag'} <= set(raw)):
                directions.append((name, 'imag'))
        self.directions = tuple(directions)
        base = dict(initial_data)
        for name in self.dependencies:
            base[name] = 0.0
        self.offset = self._direct(base)
        responses = []
        for name, part in self.directions:
            probe = dict(base)
            probe[name] = 1j if part == 'imag' else 1.0
            responses.append(self._direct(probe) - self.offset)
        self.responses = np.asarray(responses, float)
        if self.responses.shape != (len(self.directions), 4, 4, 8, 8):
            raise ValueError('Mixed fixed-basis HHSS linearisation shape mismatch')
        self.initial_parity = self.assert_reproduces(initial_data, 'muF')

    def _direct(self, data):
        _, q, _ = _basis_and_tensor(self.model, self.tree_f, self.channels,
                                    self.registry, self.geometry, self.source_path,
                                    data)
        if not np.isfinite(q).all():
            raise ValueError('Nonfinite Cartesian HHSS source')
        return q

    def tensor(self, data):
        v = np.array([(general._complex_number(data[k]).imag if part == 'imag'
                       else general._complex_number(data[k]).real)
                      for k, part in self.directions], float)
        q = self.offset + np.einsum('k,kABij->ABij', v, self.responses,
                                     optimize=True)
        if not np.isfinite(q).all():
            raise ValueError('Nonfinite reconstructed HHSS')
        return q

    def assert_reproduces(self, data, label):
        actual = self._direct(data)
        residual = _relative(self.tensor(data), actual)
        if residual > 1e-9:
            raise ValueError(f'{self.model}: {label} HHSS fixed-basis mismatch '
                             f'{residual:.3e}')
        return residual


def integrate_mixed(model, first_vec, eft_layout, registry_legacy, eft_beta,
                    new, reused, source_original, q_basis, mu_f, mu_s, *,
                    eft_on=True):
    """One DOP853 IVP on [mixed native RGBeta, Re(C/scale), Im(C/scale)]."""
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from RGE.running.intermediate.GeneralEFT1C5Mixing import (
        TERMS, calculate_master_beta,
    )
    zren = np.asarray(first_vec, float)
    if zren.shape != (eft_layout.size + 2,):
        raise ValueError('Incorrect mixed RGBeta + bilinear initial vector size')
    c0 = general.general_embed(model, source_original)
    cnorm = float(np.linalg.norm(c0))
    if not 0 < cnorm < math.inf:
        raise ValueError('Mixed initial Wilson norm zero/nonfinite')
    shape = c0.shape
    nren = zren.size
    data0, _, _ = _state_at(zren, eft_layout, registry_legacy, mu=mu_f,
                            new=new, reused=reused)
    q0 = q_basis.tensor(data0)
    initial, _ = general.general_initial_context(model, data0, q0,
                                                   source_original, mu_f, mu_s)
    from Numerical.joint_scalar_matching.core import project_weinberg
    first = calculate_master_beta(initial, c0)['scalar_pair']
    projected, _ = general.seven.project_leptonic_weinberg(first)
    direct = np.einsum('fgpqij,ABij->fgpqAB', source_original, q0,
                       optimize=True)
    direct5, _ = project_weinberg(direct, strict=False)
    parity = _relative(projected, direct5)
    if parity > 1e-8:
        raise ValueError(f'{model}: initial scalar-pair master parity={parity:.3e}')

    count = {'rhs': 0}

    def rhs(t, z):
        ren = z[:nren]
        c = general._unflatten_master(z[nren:], shape)
        data, _, _ = _state_at(ren, eft_layout, registry_legacy,
                               mu=float(math.exp(t)), new=new, reused=reused)
        q = q_basis.tensor(data)
        ctx, qgamma = general.general_initial_context(
            model, data, q, source_original, mu_f, mu_s)
        terms = calculate_master_beta(ctx, c)
        gam = np.einsum('bd,ijad->ijab', qgamma, c, optimize=True)
        terms['scalar_anomalous_dimension'] += gam + gam.swapaxes(2, 3)
        dc = LOOP * sum((terms[k] for k in TERMS), np.zeros(shape, complex))
        dr = (core.rhs_at(t, ren, layout=eft_layout,
                          registry=registry_legacy, betas=eft_beta,
                          stage='EFT1') if eft_on else np.zeros_like(ren))
        count['rhs'] += 1
        return np.r_[dr, general._flatten_master(dc)]

    y0 = np.r_[zren, general._flatten_master(c0 / cnorm)]
    sol = solve_ivp(rhs, (math.log(mu_f), math.log(mu_s)), y0,
                    method='DOP853', rtol=RTOL, atol=ATOL, max_step=0.35)
    if not sol.success:
        raise RuntimeError('Mixed simultaneous DOP853 failed: ' + sol.message)
    end_ren = sol.y[:nren, -1].copy()
    end_c = cnorm * general._unflatten_master(sol.y[nren:, -1], shape)
    if not (np.isfinite(end_ren).all() and np.isfinite(end_c).all()):
        raise ValueError('Mixed simultaneous endpoint nonfinite')
    sf = float(np.linalg.norm(end_c - end_c.swapaxes(0, 1)) / cnorm)
    ss = float(np.linalg.norm(end_c - end_c.swapaxes(2, 3)) / cnorm)
    if max(sf, ss) > 1e-6:
        raise ValueError(f'{model}: C tensor symmetry residual {sf:.3e}/{ss:.3e}')
    return dict(end_ren=end_ren, end_c=end_c, nfev=int(sol.nfev),
                steps=int(sol.t.size), rhs_calls=count['rhs'],
                initial_scalar_parity=parity, fermion_symmetry=sf,
                scalar_symmetry=ss, tensor_drift=float(np.linalg.norm(end_c-c0)/cnorm))


def _mixed_simultaneous(model, saved, quartic_source, reference_running,
                        registry, *, spinor_factor, vev):
    """Original B/C UV trajectory, coupled EFT1 RGE, original muS scalar match."""
    from scripts import RunT3GeneralSequentialFitAblation as original
    from scripts import RunT3BCGeneralFromReference as runner
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunT3ZeroExtraBoundaryBenchmark as core
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state, project_general_after_fermion,
        pack_general_scalar_state,
    )
    from Numerical.running.GeneralT3MixedExactFermionWilson import mixed_exact_boundary
    from Numerical.running.GeneralT3ConditionalBCQuarticMap import (
        apply_conditional_bc_direct_quartic_map,
    )
    from Numerical.running.GeneralT3ConditionalBCQuarticNormalization import (
        apply_conditional_bc_channel_normalization,
    )
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration.RunSmallLambdaRefit import TARGET
    from Numerical.lo_certification.couplings.extra_hhss_geometry import model_extra_geometry
    from Numerical.running.GeneralT3SequentialLeadingLog import encode, decode, response
    from scripts.RunT3SimultaneousEFT1All16 import general_initial_context

    if model not in MIXED:
        raise ValueError('This runner handles only mixed B_-1,C_-1')
    mu_f, mu_s = original.SCALES[1:3]
    fit = original.load(saved/'result.json')
    source = runner.build_reference_source(model)
    source['cfg'] = original.load(saved/'fitted_config.json')
    old_registry, new, reused = (source['registry'], source['new'], source['reused'])
    fraction, mixing_fraction = original.mixed_fit_fractions(fit)
    bundle = bc.beta_bundle(model, source, old_registry)
    coeff = Path(DEFAULT_SOURCE) / ({'B_-1':'B_m1_final',
                                      'C_-1':'C_m1_final'}[model] +
              '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled = bc.compile_source(coeff)
    target = OscillationFitTarget.from_json(TARGET)
    initial, uv_layout, _, _ = pilot.prepare_initial(
        source, old_registry,
        'general_quartics_and_mixing' if mixing_fraction > 0 else 'general_quartics',
        fraction, mixing_fraction)
    _, channels = quartic_source
    geometry = model_extra_geometry(model, registry['models'][model],
                                    include_components=True)
    cases = {}
    source_info = {}
    comparisons = {}
    for uv_on in (False, True):
        if uv_on:
            uv_sol = core.integrate_segment(
                initial,
                lambda t,x: core.rhs_at(t,x,layout=uv_layout,
                                        registry=old_registry,
                                        betas=bundle['UV'],stage='UV'),
                original.SCALES[0],mu_f,n_points=3)
            uv_vec, uv_nfev = uv_sol.y[:,-1].copy(), int(uv_sol.nfev)
        else:
            uv_vec, uv_nfev = initial.copy(), 0
        heavy = unpack_general_scalar_state(
            uv_vec[:uv_layout.size],uv_layout,old_registry,uv=True,
            mu_gev=mu_f).baseline
        eft = project_general_after_fermion(
            unpack_general_scalar_state(uv_vec[:uv_layout.size],uv_layout,
                                       old_registry,uv=True,mu_gev=mu_f),
            old_registry)
        start_vec, eft_layout = pack_general_scalar_state(eft,old_registry)
        start_vec = np.r_[start_vec,uv_vec[uv_layout.size:]]
        eft_beta = core.check_payload(
            bundle['payloads']['EFT1'],model=model,stage='EFT1',
            rep=eft_layout.representation,
            base_names=tuple(eft_layout.names[:-len(new)]),
            new_names=new,reused_names=reused)
        data_f, _, _ = _state_at(start_vec,eft_layout,old_registry,
                                  mu=mu_f,new=new,reused=reused)
        # Crucial difference from old conditional fitting: do NOT rebuild F
        # matching with the final muS scalar basis and call it a muF boundary.
        tree_f, tree_f_info = mixed_exact_boundary(
            model,heavy,data_f,repo_root=original.ROOT,mu_f=mu_f)
        R_f, q_f, _ = _basis_and_tensor(
            model,tree_f,channels,registry,geometry,
            saved/'fitted_config.json',data_f)
        source_c = _cartesian_source(tree_f.raw,R_f)*complex(spinor_factor)
        q_basis = FixedCartesianHHSS(
            model,tree_f,channels,registry,geometry,
            saved/'fitted_config.json',data_f)
        if _relative(q_basis.tensor(data_f),q_f)>1e-9:
            raise ValueError('Initial Cartesian HHSS basis disagrees')
        source_info[str(int(uv_on))] = tree_f_info

        for eft_on in (False,True):
            combined = integrate_mixed(
                model,start_vec,eft_layout,old_registry,eft_beta,new,reused,
                source_c,q_basis,mu_f,mu_s,eft_on=eft_on)
            ren_s = combined['end_ren']
            data_s, _, endpoint = _state_at(ren_s,eft_layout,old_registry,
                                            mu=mu_s,new=new,reused=reused)
            parity_s = q_basis.assert_reproduces(data_s,'muS')
            if eft_on:
                separate = core.integrate_segment(
                    start_vec,
                    lambda t,x: core.rhs_at(t,x,layout=eft_layout,
                                            registry=old_registry,
                                            betas=eft_beta,stage='EFT1'),
                    mu_f,mu_s,n_points=3)
                ref = separate.y[:,-1]
                comparison = max(_relative(ren_s[b.start:b.stop],ref[b.start:b.stop])
                                 for b in eft_layout.blocks)
                comparison = max(comparison,
                                 _relative(ren_s[eft_layout.size:],ref[eft_layout.size:]))
                separate_nfev = int(separate.nfev)
                if comparison > 5e-5:
                    raise ValueError(f'{model}: mixed RGBeta endpoint mismatch '
                                     f'{comparison:.3e}')
            else:
                comparison, separate_nfev = 0.,0
            # Original matching convention is applied ONLY at muS; its scalar
            # mass eigenbasis and conditional channel rescalings are unchanged.
            tree_s, tree_s_info = mixed_exact_boundary(
                model,heavy,data_s,repo_root=original.ROOT,mu_f=mu_f)
            snap = original.snapshot(
                model,general._source_scalar_state(data_s),saved/'fitted_config.json')
            tensor = original.scalar_tensors(
                model,tree_s,channels,snap,registry,geometry=geometry)
            tensor = apply_conditional_bc_direct_quartic_map(model,tensor)
            tensor = apply_conditional_bc_channel_normalization(model,tensor)
            scalar = original.ScalarEndpoint(tensor,data_s)
            fixed = response(tree_s,scalar,mu_f=mu_f,mu_s=mu_s,
                             spinor_factor=complex(spinor_factor),
                             check_scale_cancellation=True)
            finite = decode(fixed['c5_scalar_finite_bubble'])
            fixed_log = decode(fixed['c5_eft1_leading_log'])
            induced, ward = general.seven.project_leptonic_weinberg(combined['end_c'])
            c5_one = bc.c5_at_endpoint(model,endpoint,
                                       bc.heavy_flavour_at_muF(heavy),compiled)
            norm=lambda a: float(np.linalg.norm(a))
            rel=lambda a,b: norm(a-b)/max(norm(b),1e-300)
            metrics = {
                'simultaneous_vs_original_fixed_log_rel':rel(induced,fixed_log),
                'mixed_boundary_muF_vs_muS_relative':rel(source_c,
                    _cartesian_source(tree_s.raw,
                        _basis_and_tensor(model,tree_s,channels,registry,geometry,
                            saved/'fitted_config.json',data_s)[0]) * complex(spinor_factor)),
                'combined_ode_nfev':combined['nfev'],
                'combined_ode_steps':combined['steps'],
                'combined_rhs_calls':combined['rhs_calls'],
                'combined_ode_fermion_symmetry_residual':combined['fermion_symmetry'],
                'combined_ode_scalar_symmetry_residual':combined['scalar_symmetry'],
                'combined_ode_weinberg_projection_residual':float(ward),
                'initial_scalar_quartic_source_parity_rel':combined['initial_scalar_parity'],
                'hhss_original_basis_parity_muF_rel':q_basis.initial_parity,
                'hhss_original_basis_parity_muS_rel':parity_s,
                'separate_rgbeta_endpoint_max_block_rel':comparison,
                'separate_rgbeta_nfev':separate_nfev,
            }
            for w_on, c5 in ((False,finite),(True,finite+induced)):
                c5 = original.c5_checked(c5)
                for sm_on in (False,True):
                    code = f'{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}'
                    result = original.low_energy(
                        original.scalar_snapshot_state(scalar),c5,
                        on=sm_on,mu_s=mu_s,mu_low=original.SCALES[3],
                        vev=vev,target=target)
                    result.update({
                        'case':code,'uv_on':uv_on,'eft1_on':eft_on,
                        'wilson_log_on':w_on,'smeft_on':sm_on,
                        'uv_nfev':uv_nfev,'eft1_nfev':combined['nfev'],
                        'c5_at_scalar_threshold':encode(c5),
                        'c5_induced_wilson_log':encode(induced),
                        'c5_simultaneous_running':encode(induced),
                        'c5_original_fixed_order_llss_log':encode(fixed_log),
                        'c5_joint_scalar_bubble':encode(finite),
                        'c5_matching_anchored_reference':encode(original.c5_checked(c5_one)),
                        'c5_log_over_one_scale_norm':norm(induced)/max(norm(finite+induced),1e-300),
                        'matching_approximation':'mixed_fixed_cartesian_muF_source_simultaneous_7term_plus_original_finite_bubble',
                        'hard_matching_independently_calculated':False,
                        'full_physical_one_loop_matching':False,
                        'fixed_cartesian_wilson_basis':True,
                        'moving_XY_basis_in_wilson_RGE':False,
                        'non_HHSS_scalar_ansatz_provisional':True,
                        'per_general_interaction':fixed['per_channel'],
                        'quartic_provenance':fixed['quartic_provenance'],
                        'matching_scale_cancellation':fixed['matching_scale_cancellation_residual'],
                        **metrics,
                    })
                    cases[code] = result
            comparisons[f'{int(uv_on)}{int(eft_on)}'] = metrics
            source_info[f'{int(uv_on)}{int(eft_on)}_muS_reference'] = tree_s_info
            print(f'{model} U={int(uv_on)} E={int(eft_on)} : '
                  f'mixed joint nfev={combined["nfev"]} '
                  f'sim/originalLL={metrics["simultaneous_vs_original_fixed_log_rel"]:.3e} '
                  f'RGBeta parity={comparison:.3e}',flush=True)
    return cases, {
        'muF_exact_mixed_source':source_info,
        'C5_source_path':str(coeff),
        'C5_source_sha256':original.hash_file(coeff),
        'fitted_config_sha256':original.hash_file(saved/'fitted_config.json'),
        'hhss_original_basis_directions':list(q_basis.directions),
        'endpoint_comparisons':comparisons,
        'simultaneous_rhs':True,
        'fixed_original_cartesian_basis':True,
        'muF_source_rebuilt_at_muS_for_running':False,
        'existing_muS_conditional_scalar_finite_matching_used':True,
        'complete_physical_C5':False,
    }


def self_test():
    """Check field-basis covariance and combined ODE packing without project data."""
    general.self_test()
    rng=np.random.default_rng(314159)
    a=rng.standard_normal((8,8))
    Q,_=np.linalg.qr(a)
    source=(rng.standard_normal((3,3,2,2,8,8))+
            1j*rng.standard_normal((3,3,2,2,8,8)))
    source=.5*(source+source.swapaxes(-1,-2))
    q=rng.standard_normal((4,4,8,8));q=.5*(q+q.swapaxes(-1,-2))
    c_cart=_cartesian_source(source,Q)
    q_cart=_cartesian_hhss(q,Q)
    lhs=np.einsum('fgpqij,ABij->fgpqAB',c_cart,q_cart,optimize=True)
    rhs=np.einsum('fgpqab,ABab->fgpqAB',source,q,optimize=True)
    if not np.allclose(lhs,rhs,rtol=1e-11,atol=1e-11):
        raise AssertionError('Mixed fixed-basis scalar pair contraction not invariant')
    if not np.allclose(
            np.einsum('fgpqij,ia,jb->fgpqab',c_cart,Q,Q,optimize=True),
            source,rtol=1e-11,atol=1e-11):
        raise AssertionError('Mixed LLSS original->XY roundtrip failed')
    print('PASS: mixed scalar fixed-basis tensor rotation and LLSS contraction')


def repair_csv(argv):
    """Repair an already-computed B/C campaign without integrating again.

    Restores the historical CSV schema from existing saved C5 matrices. The
    joint ODE outputs and all other case fields remain unchanged.
    """
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.running.GeneralT3SequentialLeadingLog import decode

    if '--run-id' not in argv:
        raise ValueError('--repair-csv requires --run-id')
    run_id = argv[argv.index('--run-id') + 1]
    if (not run_id or len(run_id)>64 or
            not all(c.isalnum() or c in '_-' for c in run_id)):
        raise ValueError('Invalid run id')
    root = (original.ROOT/'output/studies/general_scalar_extended/'
            'experimental_simultaneous_eft1_mixed_BC'/run_id)
    summary_path = root/'summary.json'
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    summary = original.load(summary_path)
    expected = tuple(summary.get('requested_models', ()))
    if not expected or not set(expected)<=set(MIXED):
        raise ValueError('Stored campaign is not a mixed B/C campaign')
    if '--models' in argv:
        models = tuple(x.strip() for x in
                       argv[argv.index('--models')+1].split(',') if x.strip())
        if models != expected:
            raise ValueError('Requested models differ from saved campaign')
    rows = []
    repaired = 0
    for model in expected:
        tag = original.tag(model)
        file = root/tag/'cases.json'
        if not file.is_file():
            raise FileNotFoundError(file)
        row = original.load(file)
        if row.get('status')!='Computed' or row.get('model')!=model:
            raise ValueError(f'{model}: not a computed case record')
        cases = row.get('cases',{})
        if len(cases)!=16:
            raise ValueError(f'{model}: expected 16 saved cases')
        model_repaired = 0
        for case in cases.values():
            if 'c5_log_over_one_scale_norm' not in case:
                induced = decode(case['c5_induced_wilson_log'])
                bubble = decode(case['c5_joint_scalar_bubble'])
                num = float(np.linalg.norm(induced))
                den = float(np.linalg.norm(induced+bubble))
                if not np.isfinite([num,den]).all():
                    raise ValueError(f'{model}: nonfinite C5 recovery')
                case['c5_log_over_one_scale_norm'] = num/max(den,1e-300)
                repaired += 1
                model_repaired += 1
        # Persist only missing derived metadata; preserve all ODE results.
        if model_repaired:
            original.save(file,row)
        rows.append(row)
    original.write_csv(rows,root)
    print(f'Repaired {repaired} derived C5 ratios; CSV: {root/"all16_cases.csv"}')
    return 0


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if argv==['--self-test']:
        self_test();return 0
    if '--repair-csv' in argv:
        return repair_csv(argv)
    if '--run-id' not in argv or '--resume' in argv:
        raise ValueError('A fresh --run-id is required; --resume disabled')
    if '--models' not in argv:
        raise ValueError('Explicit --models B_-1,C_-1 required')
    names=tuple(x.strip() for x in argv[argv.index('--models')+1].split(',') if x.strip())
    if not names or len(set(names)) != len(names) or not set(names)<=set(MIXED):
        raise ValueError('This mixed runner supports only B_-1,C_-1')
    from scripts import RunT3GeneralSequentialFitAblation as original
    original._mixed = _mixed_simultaneous
    original.OUT = (original.ROOT/'output/studies/general_scalar_extended/'
                    'experimental_simultaneous_eft1_mixed_BC')
    base_save=original.save
    def save(path,payload):
        if isinstance(payload,dict):
            payload=dict(payload)
            if payload.get('format')=='t3-full-general-16_U_E_W_S_ablation-v1':
                payload.update(format=SCHEMA,
                               running_background='simultaneous_RGBeta_and_7term_master',
                               model_scope='B_-1,C_-1',
                               complete_physical_C5=False,
                               one_scale_C5_used_as_W_on_reference=False)
            elif payload.get('status')=='Computed' and 'cases' in payload:
                payload.update(full_sequential_general_C5_certified=False,
                               matching_approximation='experimental_mixed_fixed_basis_simultaneous_EFT1')
        base_save(path,payload)
    original.save=save
    return original.main(argv)


if __name__=='__main__':
    raise SystemExit(main())
