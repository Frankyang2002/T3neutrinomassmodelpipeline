"""Numerical all-16 UV -> all-F-at-MF -> scalar-only EFT1 batch.

Uses existing repository RGBeta expression evaluator, state/vector classes and
DOP853 ODE solvers. B does not perform Wilson evolution or scalar matching.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import traceback
from typing import Any, Mapping

import numpy as np

from .contracts import (MIXED, SELF, MODELS, HEAVY, ContractError,
                        DEFAULT_BATCH1, check_pair, check_finite, encode_complex,
                        hash_file, read_json, safe_run_id, save_json, tag,
                        verify_batch1_provenance)
from .fermion import (evaluate_formal_llss, evaluate_mixed_native_candidate,
                      check_external_hard_matching)


@dataclass(frozen=True)
class Options:
    root: Path
    run_id: str = 'B_all16_qualitative_v1'
    batch1: str = DEFAULT_BATCH1
    mu_uv: float = 1e6
    mu_f: float = 1e5
    mu_s: float = 1e3
    uv_samples: int = 17
    eft1_samples: int = 33
    rtol: float = 1e-8
    atol: float = 1e-11
    hard_matching_root: str | None = None
    force: bool = False

    def validated(self) -> 'Options':
        safe_run_id(self.run_id)
        if not (math.isfinite(self.mu_uv) and math.isfinite(self.mu_f)
                and math.isfinite(self.mu_s)
                and self.mu_uv > self.mu_f > self.mu_s > 0):
            raise ContractError('Required scale hierarchy mu_UV > MF > MS > 0')
        if min(self.uv_samples, self.eft1_samples) < 3:
            raise ContractError('Need at least three samples in each segment')
        if not (0 < self.rtol <= 1e-3 and 0 < self.atol <= 1e-5):
            raise ContractError('Invalid ODE tolerances')
        return self


def _git(root: Path, *args: str) -> str | None:
    try:
        p = subprocess.run(['git', *args], cwd=root, capture_output=True,
                           check=True, timeout=8)
        return p.stdout.decode('utf-8', 'replace').strip()
    except (OSError, subprocess.SubprocessError):
        return None


def _matrix_pattern(name: str, model: str) -> np.ndarray:
    """Three-generation, complex and non-diagonal, deterministic benchmark.

    No fit, no PMNS/CKM claim. Nonzero allowed channels are distinguished;
    MF remains diagonal initially (a valid heavy mass basis at mu_UV).
    """
    if name in ('yu', 'yd', 'ye'):
        diagonal = {'yu': [1.3e-5, 4.8e-3, .79],
                    'yd': [2.4e-5, 4.5e-4, .021],
                    'ye': [2.9e-6, 5.9e-4, .010]}[name]
        y = np.diag(np.asarray(diagonal, dtype=complex))
        y[0,1] = 0.005 * diagonal[1] * (1+0.32j)
        y[1,2] = 0.003 * diagonal[2] * (1-0.24j)
        y[2,0] = 0.002 * diagonal[2] * (1+0.17j)
        return y
    if name in ('y1','y2','y1Conj','y2Conj'):
        fi = MODELS.index(model)
        factor = {'y1': .028,'y2': .023,'y1Conj': .014,'y2Conj': .016}[name]
        phase = 0.15 + 0.015 * fi
        return factor * np.asarray([
            [1+.11j, .20-.14j, .08+.06j],
            [.14+.12j, .86-.06j, .19+.13j],
            [.07-.17j, .12+.05j, 1.04+.09j]], dtype=complex) * np.exp(1j*phase)
    raise ContractError(f'Not a Yukawa matrix: {name}')


def _bench_overrides(model: str, beta_names: list[str], opts: Options) -> dict:
    """Named benchmark values; every unknown required beta is still rejected
    by GeneralRunningBenchmark.initial_state/spec_for (never zero-filled)."""
    result: dict[str, Any] = {}
    for name in beta_names:
        if name in ('yu','yd','ye','y1','y2','y1Conj','y2Conj'):
            result[name] = encode_complex(_matrix_pattern(name, model))
        elif name == 'MF':
            result[name] = encode_complex(np.diag(
                [0.91 * opts.mu_f, 1.00 * opts.mu_f, 1.09 * opts.mu_f]).astype(complex))
        elif name == 'mS1Sq':
            result[name] = (1.32*opts.mu_s)**2
        elif name == 'mS2Sq':
            result[name] = (1.47*opts.mu_s)**2
        elif name.startswith('muT3General'):
            # Do not invent phase-removal arguments for neutral self-bilinears.
            result[name] = {'real': 0.013 * opts.mu_s**2,
                            'imag': 0.004 * opts.mu_s**2}
        elif name == 'lambdaT3':
            result[name] = {'real': 0.015, 'imag': 0.005}
        elif name.startswith('lambdaT3GeneralQ'):
            # Independent quartics, with genuinely complex phases.
            idx = sum(ord(ch) for ch in name) % 13
            result[name] = {'real': .006 + idx * .00035,
                            'imag': .001 + (idx % 5) * .0002}
        elif name in ('gY','g2','g3','lambdaH'):
            result[name] = {'gY': .36,'g2': .64,'g3': 1.04,'lambdaH': .13}[name]
        # Baseline in the source-verified benchmark initializer handles the
        # remaining real and complex scalar potential coordinates. There are
        # NO unrepresented beta keys, because initial_state enforces closure.
    return result


def _snapshot(trajectory: Any, environment: Any, names: list[str]
              ) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for i, _ in enumerate(trajectory.mu_gev):
        state = trajectory.state_at_index(i)
        env = environment(state)
        if set(env) != set(names):
            raise ContractError('Saved trajectory state keys differ from native beta inventory: '
                                f'missing={sorted(set(names)-set(env))}, '
                                f'extra={sorted(set(env)-set(names))}')
        for name in names:
            result.setdefault(name, []).append(np.asarray(env[name], dtype=np.complex128))
    stacked = {name: np.stack(values).astype(np.complex128)
               for name, values in result.items()}
    check_finite(*stacked.values())
    return stacked


def _endpoint_fingerprint(uv: Mapping[str, np.ndarray], names: list[str]) -> str:
    h = hashlib.sha256()
    for name in sorted(names):
        h.update(name.encode('utf-8') + b'\x00')
        a = np.ascontiguousarray(uv[name][-1], dtype=np.complex128)
        h.update(str(a.shape).encode('ascii') + b'\x00')
        h.update(a.tobytes())
    return h.hexdigest()


def _save_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + '.partial.npz')
    np.savez_compressed(tmp, **arrays)
    tmp.replace(path)


def _input_paths(opts: Options, model: str) -> tuple[Path, Path, Path, Path]:
    batch = opts.root / opts.batch1
    t = tag(model)
    return (batch / 'models' / f'{t}.json', batch / 'rgbeta' / t / 'uv.json',
            batch / 'rgbeta' / t / 'eft1.json', batch)


def _manifest_base(model: str, opts: Options) -> dict:
    return {'schema': 't3-parallel-B-background-v1', 'model': model,
            'input_scales_gev': {'UV':opts.mu_uv,'F':opts.mu_f,'S':opts.mu_s},
            'eft_order': ['UV', 'integrate_ALL_F_together_at_MF',
                          'EFT1_with_S1_and_S2_active', 'C_joint_scalar_matching_not_performed'],
            'numerical_scheme': 'MSbar one-loop candidate RGBeta report_betas; t=ln(mu)',
            'renormalisable_matching_at_F': 'tree-level identity projection of common scalar/SM couplings; one-loop finite renormalisable hard threshold not computed',
            'scalar_basis': 'native RGBeta invariant coordinates; no certified Matchete dictionary',
            'mass_basis': 'initial diagonal positive heavy MF; UV one-loop matrix running; heavy singular/Takagi values at MF',
            'SM_yukawa_basis': 'complex 3x3 numerical input; not a fitted physical CKM basis',
            'full_general_symmetry_certified': False,
            'matchete_rgbeta_coupling_dictionary_certified': False,
            'dimension_five_running_performed': False,
            'scalar_matching_performed': False,
            'complete_one_loop_matching': False,
            'physical_neutrino_mass_claim': False}


def run_one(opts: Options, model: str) -> dict:
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    from Numerical.core.RGBetaEvaluator import state_environment
    from Numerical.running.IntermediateScalarRGBetaEvaluator import (
        intermediate_scalar_state_environment)
    from Numerical.running.FermionThresholdBoundary import (
        build_intermediate_scalar_boundary)
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment

    opts.validated()
    if model not in MODELS:
        raise ContractError(f'Invalid model: {model}')
    out = opts.root/'output/parallel_t3/background'/opts.run_id/tag(model)
    if out.exists() and not opts.force:
        # Never silently overwrite an earlier completed/partially blocked run.
        raise ContractError(f'Output already exists: {out}; choose new run-id '
                            'or use --force to rebuild only Workstream B output')
    out.mkdir(parents=True, exist_ok=True)
    if opts.force:
        # Only output files owned by B; clear stale success when retry blocks.
        for filename in ('background.npz','background_manifest.json',
                         'fermion_boundary_formal.npz','fermion_boundary_manifest.json',
                         'blocker.json'):
            (out/filename).unlink(missing_ok=True)
    record_path, uv_path, eft_path, batch = _input_paths(opts, model)
    record, uv_json, eft_json = [read_json(p) for p in (record_path, uv_path, eft_path)]
    names_uv, names_eft = check_pair(model, uv_json, eft_json)
    provenance = verify_batch1_provenance(opts.root, batch, model, record,
                                          uv_path, eft_path)
    overrides = _bench_overrides(model, names_uv, opts)
    declared, state = initial_state(uv_json, mu_uv=opts.mu_uv, mu_f=opts.mu_f,
                                    mu_s=opts.mu_s, overrides=overrides)
    if declared != model:
        raise ContractError('Benchmark initializer returned the wrong model')
    uv_scales = np.geomspace(opts.mu_uv, opts.mu_f, opts.uv_samples)
    eft_scales = np.geomspace(opts.mu_f, opts.mu_s, opts.eft1_samples)
    uv_run = run_uv_segment(state, uv_json, opts.mu_f,
                            save_scales_gev=uv_scales, rtol=opts.rtol, atol=opts.atol)
    state_f = uv_run.final_state
    at_f = build_intermediate_scalar_boundary(state_f, matching_scale_gev=opts.mu_f)
    eft_run = run_intermediate_scalar_segment(
        at_f, eft_json, opts.mu_s, save_scales_gev=eft_scales,
        rtol=opts.rtol, atol=opts.atol)
    uv_snap = _snapshot(uv_run, state_environment, names_uv)
    eft_snap = _snapshot(eft_run, intermediate_scalar_state_environment, names_eft)
    for name in names_eft:
        if not np.allclose(uv_snap[name][-1], eft_snap[name][0], atol=1e-9, rtol=1e-11):
            raise ContractError(f'{model}: discontinuous common coupling at MF: {name}')
    if HEAVY & set(eft_snap):
        raise ContractError(f'{model}: heavy F still active in saved EFT1')
    if not all(np.isfinite(v).all() for v in [uv_run.y, eft_run.y]):
        raise ContractError(f'{model}: nonfinite RGE trajectory')
    if (np.max(np.abs(uv_snap['gY'])) >= 4*math.pi or
        np.max(np.abs(eft_snap['gY'])) >= 4*math.pi):
        raise ContractError(f'{model}: gauge coupling nonperturbative')
    # A rough numerical benchmark veto, not a proof of vacuum stability.
    for stage in (uv_snap, eft_snap):
        for name, value in stage.items():
            if name.startswith('lambda') and np.max(np.abs(value)) >= 4*math.pi:
                raise ContractError(f'{model}: quartic {name} exceeds 4pi')
    uv_hash = _endpoint_fingerprint(uv_snap, names_uv)
    arrays = {'mu_uv_gev': np.asarray(uv_run.mu_gev, dtype=np.float64),
              'mu_eft1_gev': np.asarray(eft_run.mu_gev, dtype=np.float64),
              't_uv': np.asarray(uv_run.t, dtype=np.float64),
              't_eft1': np.asarray(eft_run.t, dtype=np.float64)}
    arrays.update({f'uv__{n}': x for n,x in uv_snap.items()})
    arrays.update({f'eft1__{n}': x for n,x in eft_snap.items()})
    _save_npz(out/'background.npz', arrays)
    manifest = _manifest_base(model, opts)
    manifest.update({
        'status': 'CandidateRenormalisableBackgroundSuccess',
        'input_provenance': provenance,
        'local_git_head': _git(opts.root,'rev-parse','HEAD'),
        'local_git_working_tree': _git(opts.root,'status','--short'),
        'beta_symbol_order_uv': names_uv, 'beta_symbol_order_eft1': names_eft,
        'beta_symbols_uv': {n: {'kind':'complex128', 'shape_per_sample':list(uv_snap[n].shape[1:])}
                            for n in names_uv},
        'beta_symbols_eft1': {n: {'kind':'complex128', 'shape_per_sample':list(eft_snap[n].shape[1:])}
                              for n in names_eft},
        'UV_endpoint_sha256': uv_hash,
        'benchmark_input_overrides': overrides,
        'complete_initial_uv_values': encode_complex(state_environment(state)),
        'SM_at_MS_gauge': {n:float(eft_snap[n][-1].real) for n in ('gY','g2','g3')},
        'scalar_coordinates_at_MS': {
            n: encode_complex(eft_snap[n][-1]) for n in names_eft
            if n.startswith('lambda') or n.startswith('mS') or n.startswith('muT3')},
        'scalar_mass_matrix_at_MS_status':
            'coordinate_values_only; require Workstream C exact CG/basis map',
        'neutral_self_bilinear_complex_phase':
            'not_certified_by_native_real_splitting' if model in SELF else 'not_applicable',
        'mixed_scalar_mass_basis':
            'requires_Workstream_C_native_XY_adapter' if model in MIXED else 'not_applicable',
        'UV_ODE': {'success':uv_run.solver_success, 'nfev':uv_run.nfev,
                   'message':uv_run.solver_message},
        'EFT1_ODE': {'success':eft_run.solver_success, 'nfev':eft_run.nfev,
                     'message':eft_run.solver_message},
        'background_npz_sha256':hash_file(out/'background.npz'),
        'numerical_engine_source_sha256': {
            relative: hash_file(opts.root/relative) for relative in (
                'Numerical/core/RGBetaEvaluator.py',
                'Numerical/core/GeneralCouplingLayout.py',
                'Numerical/orchestration/GeneralRunningBenchmark.py',
                'Numerical/running/UVRunner.py',
                'Numerical/running/IntermediateScalarRunner.py',
                'Numerical/running/FermionThresholdBoundary.py',
            )},
        'last_eft1_sample_is_MS': bool(np.isclose(eft_run.mu_gev[-1], opts.mu_s)),
        'all_fermion_betas_removed_at_MF': sorted(set(names_uv)-set(names_eft)),
        'common_state_continuity_at_MF': True,
    })
    save_json(out/'background_manifest.json', manifest)
    f_manifest = {'schema':'t3-parallel-B-fermion-boundary-v1', 'model':model,
                  'mu_F_gev':opts.mu_f, 'same_run_uv_state_sha256':uv_hash,
                  'formal_LLsS_tree_available':False,
                  'hard_C5_status':'hard_C5_unavailable',
                  'hard_C5_is_physically_zero':None,
                  'complete_one_loop_matching':False,
                  'gammaCC_to_weyl_signed_conversion_certified':False,
                  'active_BSM_fields_below_MF':['S1','S2'],
                  'all_three_heavy_fermions_decoupled_at_single_MF':True,
                  'fermion_matching_source_basis':'Matchete GammaCC source, not canonical Weyl'}
    boundary_status = 'Blocked_MissingFThreshold'
    b_arrays: dict[str, np.ndarray] = {}
    try:
        if model in MIXED:
            b_arrays, detail = evaluate_mixed_native_candidate(
                model, state_f, at_f, root=opts.root, mf=opts.mu_f)
            boundary_status = 'FormalNativeMixedXYFermionTreeCandidate'
        else:
            b_arrays, detail = evaluate_formal_llss(
                model, state_f, root=opts.root, batch=batch,
                record=record, mf=opts.mu_f)
            boundary_status = 'FormalSourceExactFermionTreeLLSS'
        f_manifest.update(detail)
        f_manifest['formal_LLsS_tree_available'] = True
    except Exception as exc:
        # Boundary can fail independently of the completed renormalisable run.
        # Preserve exact reason and all RGE output. NEVER replace by a proxy.
        boundary_status = 'Blocked_'+type(exc).__name__
        f_manifest.update({'blocker':str(exc),
                           'traceback_tail':traceback.format_exc(limit=4)})
    if opts.hard_matching_root is not None:
        try:
            hroot = opts.root / opts.hard_matching_root / tag(model)
            hard, hmeta = check_external_hard_matching(
                hroot/'hard_C5.npz', hroot/'hard_C5_manifest.json', model=model,
                mf=opts.mu_f, uv_state_sha=uv_hash)
            b_arrays['C5_hard_external'] = hard
            f_manifest['hard_C5_status'] = 'external_independent_one_loop_candidate'
            f_manifest['hard_C5_basis'] = hmeta['basis']
            f_manifest['hard_C5_manifest_sha256'] = hash_file(hroot/'hard_C5_manifest.json')
            # Hard LLHH alone is not proof of other hard operators or completeness.
        except Exception as exc:
            f_manifest['hard_C5_status'] = 'hard_C5_requested_but_blocked'
            f_manifest['hard_C5_blocker'] = str(exc)
    if b_arrays:
        _save_npz(out/'fermion_boundary_formal.npz', b_arrays)
        f_manifest['fermion_boundary_npz_sha256'] = hash_file(out/'fermion_boundary_formal.npz')
    f_manifest['status'] = boundary_status
    save_json(out/'fermion_boundary_manifest.json', f_manifest)
    return {'model':model,'tag':tag(model), 'background':manifest['status'],
            'fermion_boundary':boundary_status,
            'hard_C5_status':f_manifest['hard_C5_status'],
            'model_output':str(out.relative_to(opts.root))}


def run_all(opts: Options, models: tuple[str,...] = MODELS) -> dict:
    opts.validated()
    root_out = opts.root/'output/parallel_t3/background'/opts.run_id
    root_out.mkdir(parents=True, exist_ok=True)
    summary = {'schema':'t3-parallel-B-batch-v1','run_id':opts.run_id,
               'created_utc':datetime.now(timezone.utc).isoformat(),
               'local_repo':str(opts.root), 'models_requested':list(models),
               'input_batch1_root':opts.batch1,
               'physics_scope': 'one-loop candidate renormalisable UV and EFT1 + formal F-tree LLSS only',
               'model_results':{}}
    for model in models:
        print(f'B {model}: start', flush=True)
        try:
            row = run_one(opts, model)
        except Exception as exc:
            row = {'model':model, 'tag':tag(model),
                   'background':'Blocked_'+type(exc).__name__,
                   'fermion_boundary':'NotAttempted_BackgroundFailed',
                   'hard_C5_status':'hard_C5_unavailable',
                   'blocker':str(exc),'traceback_tail':traceback.format_exc(limit=5)}
            out = root_out / tag(model)
            out.mkdir(parents=True,exist_ok=True)
            save_json(out/'blocker.json',row)
        summary['model_results'][model] = row
        save_json(root_out/'summary.json', summary)
        print(f"B {model}: {row['background']} / {row['fermion_boundary']}",flush=True)
    summary['completed_utc'] = datetime.now(timezone.utc).isoformat()
    summary['background_success_count'] = sum(
        r['background']=='CandidateRenormalisableBackgroundSuccess'
        for r in summary['model_results'].values())
    summary['formal_f_tree_success_count'] = sum(
        r['fermion_boundary'] in ('FormalSourceExactFermionTreeLLSS',
                                    'FormalNativeMixedXYFermionTreeCandidate')
        for r in summary['model_results'].values())
    save_json(root_out/'summary.json', summary)
    return summary
