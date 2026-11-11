"""Resumable exploratory Yukawa optimisation against the existing coupled EFT1.

Project path: scripts/RunT3CoupledYukawaOvernight.py
No production files are modified. Uses EXACT existing one-model coupled physics
runners; trims only the stage-freeze control loops so the objective is case 1111.
The 18 real y1/y2 coordinates alone vary. Does not change conjugate Yukawas.

This is an exploratory local search, NOT a global optimum or full one-loop match.
"""
from __future__ import annotations

import argparse
import contextlib
import ast
import hashlib
import inspect
import json
import math
import shutil
import sys
import tempfile
import textwrap
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output/studies/general_scalar_extended'
OUTPUT = BASE / 'coupled_yukawa_overnight'
SOURCE = BASE / 'thesis_general_all16'
MIXED = {'B_-1', 'C_-1'}
SPECIAL = {'A_-4', 'D_-2'}


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    tmp.replace(path)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def append_jsonl(path, obj):
    with Path(path).open('a', encoding='utf-8') as f:
        f.write(json.dumps(obj, allow_nan=False) + '\n')
        f.flush()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def narrow_to_1111(func, targets, expected):
    """Compile a guarded, source-identical function with only ablation loops cut.

    Every modified loop has an exact target-name and its pre-existing iterator
    is sliced [-1:], preserving the all-on branch, the body and call graph.
    Abort on upstream source-shape changes instead of silently changing physics.
    """
    tree = ast.parse(textwrap.dedent(inspect.getsource(func)))
    if len(tree.body) != 1 or not isinstance(tree.body[0], ast.FunctionDef):
        raise RuntimeError('Unexpected coupled runner source shape: ' + func.__name__)
    counts = {key: 0 for key in targets}

    class Trim(ast.NodeTransformer):
        def visit_For(self, node):
            self.generic_visit(node)
            if isinstance(node.target, ast.Name):
                variable = node.target.id
            elif isinstance(node.target, ast.Tuple) and node.target.elts and isinstance(node.target.elts[0], ast.Name):
                variable = node.target.elts[0].id
            else:
                return node
            if variable in counts:
                counts[variable] += 1
                # Use the final original all-True alternative, never synthesize a coefficient.
                node.iter = ast.Subscript(value=node.iter,
                    slice=ast.Slice(lower=ast.UnaryOp(op=ast.USub(), operand=ast.Constant(value=1)),
                                    upper=None, step=None), ctx=ast.Load())
            return node

    tree = Trim().visit(tree)
    if counts != expected:
        raise RuntimeError(f'{func.__name__}: ablation source changed: {counts} != {expected}')
    ast.fix_missing_locations(tree)
    env = dict(func.__globals__)
    exec(compile(tree, inspect.getsourcefile(func) or '<runner>', 'exec'), env)
    return env[func.__name__], env


def setup_physics(model, response_root):
    from scripts import RunT3GeneralSequentialFitAblation as original
    from scripts import RunT3SimultaneousEFT1All16 as simultaneous
    from scripts import RunT3SimultaneousEFT1MixedBC as mixed

    scalar_prediction, env_pred = narrow_to_1111(
        simultaneous._scalar_endpoint_prediction, {'w_on', 'sm_on'},
        {'w_on': 1, 'sm_on': 1})
    nonmixed, env_nonmixed = narrow_to_1111(
        simultaneous._nonmixed_general, {'uv_on', 'eft_on'},
        {'uv_on': 1, 'eft_on': 1})
    env_nonmixed['_scalar_endpoint_prediction'] = scalar_prediction
    mixed_runner, _env_mixed = narrow_to_1111(
        mixed._mixed_simultaneous, {'uv_on', 'eft_on', 'w_on', 'sm_on'},
        {'uv_on': 1, 'eft_on': 1, 'w_on': 1, 'sm_on': 1})

    if model in SPECIAL:
        if not response_root:
            raise ValueError(f'{model}: provide --response-root for verified cached non-HHSS context')
        from Numerical.running.T3CoupledQuarticContext import CandidateQuarticContext
        context = CandidateQuarticContext.load(ROOT, Path(response_root).resolve(), model)
        env_nonmixed['general_quartic'] = lambda m, state, hhss: context.tensor(state, hhss)
        # Global accessed by integrate_simultaneously() is on the source module.
        old_general_quartic = simultaneous.general_quartic
        simultaneous.general_quartic = lambda m, state, hhss: context.tensor(state, hhss)
    else:
        old_general_quartic = None

    def restore():
        if old_general_quartic is not None:
            simultaneous.general_quartic = old_general_quartic

    return (mixed_runner if model in MIXED else nonmixed), restore, original


@contextlib.contextmanager
def symbolic_export_compatibility(model, canonical_run_id, allow_stale):
    """Keep validated cached RGEs while tolerating unrelated source edits.

    In stale mode, only the source-fingerprint comparison is satisfied with
    the immutable cached provenance. validate_cached still checks the actual
    cached Lagrangian, operator inventories, beta keys, and source paths.
    No file is modified. This DOES NOT regenerate RGEs from changed sources.
    """
    if model in MIXED:
        yield {'mode': 'mixed_separate_beta_exports'}
        return
    from scripts import RunT3GeneralSymbolicStudy as canonical
    folder = (canonical.DEFAULT_DATA_ROOT / canonical_run_id /
              canonical.model_tag(model))
    provenance_path = folder / 'provenance.json'
    provenance = read_json(provenance_path)
    historical = provenance.get('source_sha256')
    if not isinstance(historical, dict) or not historical:
        raise ValueError(f'{model}: missing historical source hashes in {provenance_path}')
    original_fingerprint = canonical.fingerprint
    current = original_fingerprint()
    all_paths = sorted(set(historical) | set(current))
    differences = [name for name in all_paths
                   if historical.get(name) != current.get(name)]
    report = {
        'model': model,
        'mode': 'cached_export_with_stale_source_override' if differences else 'canonical_source_matching',
        'canonical_run_id': canonical_run_id,
        'cached_provenance_path': str(provenance_path),
        'changed_source_files': differences,
        'cached_provenance_sha256': sha(provenance_path),
        'current_source_fingerprints': current,
        'cached_source_fingerprints': historical,
        'current_sources_recompiled': False,
    }
    if differences and not allow_stale:
        raise ValueError(
            f'{model}: symbolic cache source fingerprint changed in {len(differences)} '
            'file(s). Pass --allow-stale-symbolic-source to reuse the old '
            'cached exports as an explicit approximate calculation. '
            'Changed: ' + ', '.join(differences[:5]))
    if differences:
        print(f'{model}: WARNING reusing existing symbolic exports although '
              f'{len(differences)} source file(s) changed: '
              + ', '.join(differences[:5]), flush=True)
        canonical.fingerprint = lambda: dict(historical)
    try:
        yield report
    finally:
        canonical.fingerprint = original_fingerprint


def initial_values(config):
    o = config['base_state']['ordinary']
    arrays = []
    for name in ('y1', 'y2'):
        item = o[name]
        a = np.asarray(item['real'], float)
        b = np.asarray(item['imag'], float)
        if a.shape != (3, 3) or b.shape != (3, 3) or np.max(abs(b)) > 1e-12:
            raise ValueError('Expected real 3x3 Yukawa matrices for y1 and y2')
        arrays.append(a.ravel())
    out = np.concatenate(arrays)
    if not np.isfinite(out).all():
        raise ValueError('Nonfinite fitted Yukawas')
    return out


def params_from_vector(vector):
    return {f'{name}_{i}{j}': float(vector[k*9+i*3+j])
            for k, name in enumerate(('y1', 'y2'))
            for i in range(3) for j in range(3)}


def candidate(vector0, current, trial, args, model):
    """Seeded, reproducible mixture of local random directions and scale moves."""
    h = int(hashlib.sha256(model.encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(np.random.SeedSequence([args.seed, h, trial]))
    base = np.asarray(current, float).copy()
    amplitudes = np.maximum(abs(vector0)*args.fraction, args.floor)
    mode = trial % 5
    if mode == 0:
        # Bilinear mass-scale compensation: simultaneous common Yukawa scaling.
        base *= math.exp(rng.normal(0, args.fraction))
    elif mode == 1:
        # Relative y1/y2 rescaling keeps leading cross-Yukawa product unchanged.
        scale = math.exp(rng.normal(0, args.fraction))
        base[:9] *= scale
        base[9:] /= scale
    else:
        k = min(args.active, 18)
        indices = rng.choice(18, k, replace=False)
        base[indices] += rng.normal(size=k)*amplitudes[indices]
    return np.clip(base, -args.bound, args.bound)


def run_model(model, args, source_root, result_root, config, overrides, registry,
              release, channels, canonical):
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralSequentialFitAblation as original
    fn, restore, _ = setup_physics(model, args.response_root)
    tag = general.tag(model)
    dest = result_root/tag
    dest.mkdir(parents=True, exist_ok=True)
    best_path = dest/'best.json'
    log_path = dest/'evaluations.jsonl'
    config_sha = sha(source_root/tag/'fitted_config.json')
    override_sha = sha(source_root/tag/'general_uv_overrides.json') if model not in MIXED else None
    x0 = initial_values(config)
    source_hashes = {'config': config_sha, 'override': override_sha}
    if best_path.exists():
        best = read_json(best_path)
        if best['input_sha256'] != source_hashes:
            raise ValueError(f'{model}: fitted input files changed since optimisation started')
        if best['seed'] != args.seed or best['fraction'] != args.fraction or best['floor'] != args.floor or best['active'] != args.active or best['bound'] != args.bound:
            raise ValueError(f'{model}: resume requires unchanged search parameters')
        best_x = np.asarray(best['yukawas_flat'], float)
        start = int(best['last_trial'])+1
    else:
        best = None
        best_x = x0.copy()
        start = 0

    with symbolic_export_compatibility(model, canonical, args.allow_stale_symbolic_source) as symbolic_report, \
         tempfile.TemporaryDirectory(prefix='t3_fit_'+tag+'_', dir=str(dest)) as tmp:
        symbolic_report_path = dest / 'symbolic_export_compatibility.json'
        if symbolic_report_path.is_file():
            previous = read_json(symbolic_report_path)
            if previous != symbolic_report:
                raise ValueError(f'{model}: symbolic export/current source state changed '
                                 'since the earlier optimisation; use a fresh run ID')
        else:
            write_json(symbolic_report_path, symbolic_report)
        saved = Path(tmp)
        if model in MIXED:
            # The production mixed runner reads scalar-fraction metadata here.
            shutil.copy2(source_root/tag/'result.json', saved/'result.json')
        def evaluate(x, label):
            p = params_from_vector(x)
            write_json(saved/'fitted_config.json', general.config_with_yukawas(config,p))
            if model not in MIXED:
                write_json(saved/'general_uv_overrides.json', general.override_yukawas(overrides,p))
            if model in MIXED:
                cases, provenance = fn(model,saved,(release,channels),
                                       args.reference_running,registry,
                                       spinor_factor=args.spinor_factor,vev=args.vev)
            else:
                cases, provenance = fn(model,saved,canonical,(release,channels),registry,
                                       spinor_factor=args.spinor_factor,
                                       wolframscript=args.wolframscript,vev=args.vev)
            if set(cases) != {'1111'}:
                raise RuntimeError(f'{model}: objective must use only 1111, got {list(cases)}')
            record = cases['1111']
            score = float(record['chi2'])
            if not math.isfinite(score) or score < 0:
                raise ValueError('Invalid coupled-pipeline chi2')
            return score, record

        try:
            if best is None:
                score, record = evaluate(x0, 'seed')
                best = {'model':model, 'status':'Computed', 'chi2':score,
                        'seed_chi2':score, 'yukawas_flat':x0.tolist(),
                        'prediction_1111':record['prediction'],
                        'masses_ev':record.get('masses_ev_ascending'),
                        'last_trial':-1,'evaluations':1,'failures':0,
                        'input_sha256':source_hashes,
                        'symbolic_export_compatibility':symbolic_report['mode'],
                        'seed':args.seed,
                        'fraction':args.fraction,'floor':args.floor,
                        'active':args.active,'bound':args.bound,
                        'approximation':'coupled_RGBeta_7term_plus_joint_scalar_bubble',
                        'scalar_context':'cached_Batch3B_nonHHSS' if model in SPECIAL else ('mixed_cartesian' if model in MIXED else 'historical_provisional_nonHHSS')}
                write_json(best_path,best)
                write_json(dest/'best_fitted_config.json', general.config_with_yukawas(config, params_from_vector(x0)))
                if model not in MIXED:
                    write_json(dest/'best_general_uv_overrides.json', general.override_yukawas(overrides, params_from_vector(x0)))
                append_jsonl(log_path,{'trial':-1,'chi2':score,'improved':True,'status':'Computed'})
                print(f'{model} seed coupled chi2={score:.8g}',flush=True)
            consecutive_failures=0
            for trial in range(start,args.trials):
                if args.deadline and time.monotonic() >= args.deadline:
                    print(f'{model}: duration limit reached; saved progress',flush=True)
                    break
                x = candidate(x0,best_x,trial,args,model)
                try:
                    score, record = evaluate(x,f'trial_{trial}')
                    improved = score < float(best['chi2'])
                    consecutive_failures=0
                    append_jsonl(log_path,{'trial':trial,'chi2':score,'improved':improved,'status':'Computed'})
                    if improved:
                        best_x = x
                        best.update({'chi2':score,'yukawas_flat':x.tolist(),
                                     'prediction_1111':record['prediction'],
                                     'masses_ev':record.get('masses_ev_ascending')})
                        write_json(dest/'best_fitted_config.json', general.config_with_yukawas(config, params_from_vector(x)))
                        if model not in MIXED:
                            write_json(dest/'best_general_uv_overrides.json', general.override_yukawas(overrides, params_from_vector(x)))
                        print(f'{model} trial {trial+1}/{args.trials}: improved chi2={score:.8g}',flush=True)
                    else:
                        print(f'{model} trial {trial+1}/{args.trials}: chi2={score:.8g}, best={best["chi2"]:.8g}',flush=True)
                except (KeyboardInterrupt,SystemExit):
                    raise
                except Exception as e:
                    consecutive_failures += 1
                    best['failures'] += 1
                    append_jsonl(log_path,{'trial':trial,'status':'Failed','error':f'{type(e).__name__}: {e}', 'traceback':traceback.format_exc(limit=5)})
                    print(f'{model} trial {trial+1}: ERROR: {type(e).__name__}: {e}',flush=True)
                    if consecutive_failures >= 3:
                        best['last_trial']=trial
                        write_json(best_path,best)
                        raise RuntimeError(f'{model}: three consecutive evaluator failures; stopping') from e
                best['last_trial']=trial
                best['evaluations'] += 1
                write_json(best_path,best)
        finally:
            restore()
    return best


def main(argv=None):
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralSequentialFitAblation as original
    from Numerical.joint_scalar_matching.fermion_boundary_bridge import load_release_manifest,load_source

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id',required=True)
    parser.add_argument('--source-run',default='thesis16_merged_v1')
    parser.add_argument('--models',default=','.join(general.MODELS))
    parser.add_argument('--trials',type=int,default=24)
    parser.add_argument('--active',type=int,default=5)
    parser.add_argument('--fraction',type=float,default=0.20)
    parser.add_argument('--floor',type=float,default=1e-5)
    parser.add_argument('--bound',type=float,default=3.)
    parser.add_argument('--seed',type=int,default=20261112)
    parser.add_argument('--max-hours',type=float,default=0.,help='Stop between evaluations after this many hours (0 = no limit)')
    parser.add_argument('--response-root',type=Path,default=ROOT/'output/audits/batch3b_quartic_fixed_20261111_234059')
    parser.add_argument('--release-dir',type=Path,default=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f')
    parser.add_argument('--v4-channel-dir',type=Path,default=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels')
    parser.add_argument('--quartic-seed-dir',type=Path,default=ROOT/'output/full16_symbolic')
    parser.add_argument('--registry',type=Path,default=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    parser.add_argument('--reference-running',type=Path,default=ROOT/'output/general_observable/all16_v1')
    parser.add_argument('--spinor-factor',type=float,default=1.)
    parser.add_argument('--vev',type=float,default=246.22)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--allow-stale-symbolic-source',action='store_true',
                        help='Reuse saved full-general RGEs despite edits to source files; report differences without modifying canonical validation or exports')
    args=parser.parse_args(argv)

    if not args.run_id.replace('-','').replace('_','').isalnum():
        parser.error('Use an alphanumeric, underscore or hyphen run ID')
    if not 0 <= args.trials <= 100000 or not 1 <= args.active <= 18 or args.floor <= 0 or not 0 < args.fraction < 2 or args.bound <= 0 or args.max_hours < 0:
        parser.error('Invalid search parameters')
    if not math.isfinite(args.spinor_factor) or not math.isfinite(args.vev) or args.vev <= 0:
        parser.error('Invalid spinor factor/vev')
    models=tuple(s.strip() for s in args.models.split(',') if s.strip())
    if not models or len(models)!=len(set(models)) or set(models)-set(general.MODELS):
        parser.error('Invalid models')
    source_root=SOURCE/args.source_run
    if not (source_root/'settings.json').is_file():
        raise FileNotFoundError('Fitted source settings missing: '+str(source_root))
    canonical=read_json(source_root/'settings.json').get('canonical_run_id','full16_general_v1')
    dest=OUTPUT/args.run_id
    if dest.exists() and not (args.resume or args.check_only):
        raise FileExistsError(f'{dest} already exists; specify --resume or fresh run ID')
    if args.resume and not dest.is_dir():
        raise FileNotFoundError(f'Cannot resume missing directory: {dest}')
    if args.check_only:
        ready=True
        for model in models:
            tag=general.tag(model)
            save=source_root/tag
            required=[save/'fitted_config.json',save/'result.json',args.registry,
                      args.release_dir/tag/'boundary.json', args.v4_channel_dir/(tag+'.json')]
            if model not in MIXED: required.append(save/'general_uv_overrides.json')
            if model in SPECIAL: required.append(args.response_root)
            missing=[str(p) for p in required if not p.exists()]
            if not missing and model not in MIXED:
                try:
                    with symbolic_export_compatibility(
                        model, canonical, args.allow_stale_symbolic_source) as info:
                        print(model, 'symbolic mode:', info['mode'], flush=True)
                except (ValueError, FileNotFoundError) as exc:
                    missing.append(f'symbolic export: {exc}')
            print(model,'READY' if not missing else 'MISSING',missing,flush=True)
            ready=ready and not missing
        return 0 if ready else 2
    dest.mkdir(parents=True,exist_ok=True)
    args.deadline = time.monotonic()+args.max_hours*3600 if args.max_hours else None
    registry=read_json(args.registry)
    release_map=load_release_manifest(args.release_dir)
    for model in models:
        tag=general.tag(model)
        saved=source_root/tag
        config=read_json(saved/'fitted_config.json')
        overrides=read_json(saved/'general_uv_overrides.json') if model not in MIXED else None
        ref=load_source(release_map[model], expected_model=model)
        channels=original._load_channels(model,ref,v4=args.v4_channel_dir,seeds=args.quartic_seed_dir)
        best=run_model(model,args,source_root,dest,config,overrides,registry,ref,channels,canonical)
        print(f'{model} saved best chi2={best["chi2"]:.8g} at {dest/tag}/best.json',flush=True)
        if args.deadline and time.monotonic() >= args.deadline:
            break
    return 0


if __name__=='__main__':
    raise SystemExit(main())
