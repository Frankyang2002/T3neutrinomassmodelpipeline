#!/usr/bin/env python3
"""One-command full-general two-threshold native completion attempt.

Runs the real full-general Matchete scalar threshold, not the restricted
analytic triangle or scalar bubble.  Reuses immutable Batch1 and first F
threshold sources.  A numerical C5 is never invented from unsubtracted WXF.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import numpy as np
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.full_general_completion import (
    MODELS, CompletionBlocked, existing_evidence, put_json,
    finish_smeft, from_complex_json,
)
from Numerical.parallel_t3.background.general_batch1_api import decode_native, load_batch1
from Numerical.parallel_t3.complete.final_one_batch_fit import fit_18_real
from Numerical.parallel_t3.background.general_batch1_api import read_json

NATIVE = ROOT / 'Lagrangian/RunT3FullGeneralTwoThresholdNative.wl'
EVALUATOR = ROOT / 'Lagrangian/RunT3FullGeneralNativeC5Numerics.wl'
BENCHMARK = ROOT/'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json'
B1 = ROOT / 'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1'
B2 = ROOT / 'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1'
LAURENT = ROOT / 'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2'
OUT = ROOT / 'output/parallel_t3/full_general_final_one_batch'
FORMAT = 't3-full-general-one-batch-execution-v1'
ENV = ('T3_EXPERIMENTAL_GENERAL_FERMIONS',
       'T3_EXPERIMENTAL_GENERAL_QUARTICS',
       'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
       'T3_EXPERIMENTAL_FLAVOR_MASSES')


def ensure_model_contract(evidence, model: str) -> None:
    if evidence.model != model:
        raise CompletionBlocked('Wrong model at Batch1/Batch2 handoff')
    if abs(evidence.mf_gev-1e5)>1e-8 or abs(evidence.ms_gev-1e3)>1e-10:
        raise CompletionBlocked('Scales must be 1e5/1e3 GeV')
    if set(evidence.stage_statuses) != {'UV_general_RG',
        'F_all_fermion_raw_one_loop', 'F_Laurent'}:
        raise CompletionBlocked('Missing source-derived prior stages')


def run_wolfram(model: str, out: Path, executable: str,
                timeout: int) -> dict:
    env = dict(os.environ)
    env.update({v: '1' for v in ENV})
    cmd = [executable, '-file', str(NATIVE.resolve()), str(out.resolve()), model]
    # A(-4) first threshold completed in the initial run, and exported
    # the exact canonical inputs to stage 2. Reuse those files if present.
    # Do not reuse the failed scalar result or overwrite the old run.
    if model == 'A_-4':
        previous = (OUT / 'thesis_final_one_batch_v1' / model /
                    'debug' / 'fresh_kernel_stage_2')
        if all((previous / name).is_file() for name in
               ('input_tree.wxf', 'input_loop.wxf')):
            cmd.append(str(previous.resolve()))
    try:
        r = subprocess.run(cmd, cwd=ROOT, env=env, text=True,
                           capture_output=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise CompletionBlocked('Native two-threshold Matchete timeout') from exc
    (out/'native_two_threshold.log').write_text(
        'command: '+' '.join(cmd)+'\nexit: '+str(r.returncode)+
        '\nstdout:\n'+r.stdout+'\nstderr:\n'+r.stderr,
        encoding='utf-8')
    path = out/'native_two_threshold_status.json'
    report = read_json(path) if path.is_file() else {}
    if r.returncode != 0 or not report.get('status', '').startswith('NativeTwoThresholdC5SourceExported'):
        raise CompletionBlocked('Native two-threshold matching blocked: '+
            report.get('reason', 'see native_two_threshold.log'))
    for file in ('full_matched_smeft.wxf', 'weinberg_native_ordered_coefficient.wxf',
                 'weinberg_native_ordered_coefficient.txt'):
        if not (out/file).is_file():
            raise CompletionBlocked('Native Matchete did not save '+file)
    return report


def _complex_json_matrix(arr):
    z=np.asarray(arr,complex)
    if z.shape!=(3,3) or not np.isfinite(z).all():
        raise CompletionBlocked('Expected finite complex 3x3 Yukawa matrix')
    return [[{'real':float(v.real),'imag':float(v.imag)} for v in row]
            for row in z]


def numeric_input(evidence, b1: Path) -> dict:
    """Explicit, source-conservative candidate Matchete/RGBeta name map."""
    import numpy as np
    b = read_json(b1/'fermion_matching_input.json')
    f = {k: decode_native(v) for k, v in
         b['native_uv_rgbeta_couplings_at_muF'].items()}
    mf = np.asarray(f['MF'], complex)
    if mf.shape != (3,3):
        raise CompletionBlocked('Require complete three-generation fermion mass matrix')
    diagonal = np.diag(mf)
    ratio = np.linalg.norm(mf-np.diag(diagonal))/np.linalg.norm(mf)
    if ratio >= 1e-8:
        raise CompletionBlocked('Non-diagonal F mass requires explicit bi-unitary Yukawa rotations: '+str(ratio))
    if np.max(abs(diagonal.imag)) >= 1e-8*np.max(abs(diagonal)) or np.min(diagonal.real)<=0:
        raise CompletionBlocked('Nonpositive/complex heavy masses require physical chiral rotation')
    s = evidence.rg_at_ms
    coupling = {}
    for name, value in s.items():
        a = np.asarray(value, complex)
        if a.shape == () and abs(a.imag) <= 1e-10*max(1.,abs(a.real)):
            coupling[name] = float(a.real)
    # This mapping is stated explicitly; it does NOT resolve Matchete CG or
    # independent extra scalar quartic normalisations.
    coupling['λ'] = float(np.real(s['lambdaH']))
    for matchete, native in (('MS1','mS1Sq'),('MS2','mS2Sq')):
        if native not in s:
            raise CompletionBlocked('Missing '+native+' in full general running')
        value = float(np.real(s[native]))
        if value <= 0:
            raise CompletionBlocked('Nonpositive heavy scalar mass squared '+native)
        coupling[matchete] = float(np.sqrt(value))
    # Fields MS1/MS2 at the second threshold use their saved EFT1 values.
    coupling['MuScalar'] = evidence.ms_gev
    return {
        'format':'t3-complete-native-c5-evaluation-input-v1',
        'model':evidence.model,
        'NF':3,
        'mu_S_GeV':evidence.ms_gev,
        'fermion_masses_GeV':[float(v.real) for v in diagonal],
        'heavy_mass_off_diagonal_frobenius_relative':float(ratio),
        'yukawas':{name:_complex_json_matrix(np.asarray(f[name],complex))
                   for name in ('y1','y2')},
        'scalar_couplings':coupling,
        'Matchete_RGBeta_extra_quartic_normalisations_verified':False,
        'symbolic_two_threshold_eft1_wilson_running_included':False,
    }


def numeric_native(evidence, b1: Path, output: Path, args: argparse.Namespace) -> dict:
    """Evaluate native C5 only; an unsubtracted pole is always a hard blocker."""
    payload = numeric_input(evidence,b1)
    input_file=output/'native_numeric_input.json'
    put_json(input_file,payload)
    cmd=[args.wolframscript,'-file',str(EVALUATOR.resolve()),
         str((output/'weinberg_native_ordered_coefficient.wxf').resolve()),
         str(input_file.resolve()),str(output.resolve())]
    env=os.environ.copy()
    env.update({v:'1' for v in ENV})
    try:
        r=subprocess.run(cmd,cwd=ROOT,env=env,text=True,capture_output=True,
                         timeout=args.timeout,check=False)
    except subprocess.TimeoutExpired as exc:
        raise CompletionBlocked('Native numerical C5 conversion timed out') from exc
    (output/'native_numeric.log').write_text(
        'command: '+' '.join(cmd)+'\nexit: '+str(r.returncode)+
        '\nstdout:\n'+r.stdout+'\nstderr:\n'+r.stderr, encoding='utf-8')
    report_file=output/'native_c5_numeric_status.json'
    report=read_json(report_file) if report_file.is_file() else {}
    if r.returncode or report.get('status')!='ThreeFlavourNativeC5FiniteCandidate_NotEFT1RGVerified':
        raise CompletionBlocked('Numerical C5 from full native matching blocked: '+
             str(report.get('reason', report.get('status','see native_numeric.log'))))
    return report


def provisional_fit(evidence, model: str, output: Path, args: argparse.Namespace) -> dict:
    """Re-run FULL UV/EFT1 renormalisable RG per trial; reuse Matchete symbolic C5.

    The candidate remains incomplete without the EFT1 Wilson operator RGE,
    verified finite subtraction and model-specific scalar mass-basis map.
    """
    from scripts.RunT3ThesisGeneralBatch1 import run_one as run_general_rg
    import uuid
    source_config=read_json(args.config)
    target=read_json(ROOT/source_config['fitting']['target'])
    initial=np.concatenate([np.asarray(source_config['models'][model][k],float).ravel()
             for k in ('y1_UV_times_1e4','y2_UV_times_1e4')])*1e-4
    baseline=read_json(args.batch1/model/'manifest.json')
    symbolic_id=baseline['symbolic_run_id']
    scratch=output/'fit_trial_scratch'
    scratch.mkdir(exist_ok=True)
    def evaluate(vec):
        vec=np.asarray(vec,float)
        if vec.shape!=(18,) or not np.isfinite(vec).all():
            raise ValueError('Must supply 18 finite UV Yukawas')
        config=json.loads(json.dumps(source_config))
        config['models'][model]['y1_UV_times_1e4']=(vec[:9]/1e-4).reshape(3,3).tolist()
        config['models'][model]['y2_UV_times_1e4']=(vec[9:]/1e-4).reshape(3,3).tolist()
        nonce=uuid.uuid4().hex[:16]
        config_file=scratch/('config_'+nonce+'.json')
        put_json(config_file,config)
        trial_id='thesis_fit_work_'+nonce
        trial_b1=ROOT/'output/parallel_t3/full_general_batch1'/trial_id
        trial_dir=scratch/('eval_'+nonce)
        trial_dir.mkdir()
        try:
            row=run_general_rg(model,config_file,symbolic_id,trial_id,
                    generate_symbolic=False,samples_uv=5,samples_eft1=5,
                    rtol=1e-7,atol=1e-10,shifts=None)
            if row.get('status')!='CandidateGeneralRenormalisableRGIntegrated':
                raise CompletionBlocked('Full general RG trial did not integrate')
            b1trial=trial_b1/model
            endpoint=load_batch1(b1trial)
            # All other evidence is provenance for the same FORMAL native
            # matching source, not a claim that the trials are verified.
            ev=replace(evidence,rg_at_ms=endpoint.at_scalar_threshold())
            shutil.copyfile(output/'weinberg_native_ordered_coefficient.wxf',
                            trial_dir/'weinberg_native_ordered_coefficient.wxf')
            numeric=numeric_native(ev,b1trial,trial_dir,args)
            c5=from_complex_json(numeric['symmetric_C5_candidate_GeV_inverse'],
                                 label='trial native C5')
            result=finish_smeft(ev,c5,vev=float(config['vev_GeV']),target=target)
            return result['observables']
        finally:
            shutil.rmtree(trial_dir,ignore_errors=True)
            config_file.unlink(missing_ok=True)
            shutil.rmtree(trial_b1,ignore_errors=True)
    fitted=fit_18_real(initial,evaluate,target,max_nfev=args.fit_maxeval,
                       max_yukawa=args.fit_bound)
    put_json(output/'provisional_fit_result.json',fitted)
    return fitted


def execute(model: str, base: Path, args: argparse.Namespace) -> dict:
    out = base/model
    out.mkdir(parents=True)
    record = {
        'format': FORMAT, 'model': model,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'NotExecuted', 'heavy_fermions_integrated_simultaneously': 3,
        'scalar_multiplets_integrated_simultaneously': ['S1', 'S2'],
        'threshold_scales_GeV': {'F': 1e5, 'S': 1e3},
        'analytic_triangle_fallback_used': False,
        'scalar_bubble_fallback_used': False,
        'complete_numerical_C5_available': False,
        'fitted_yukawas_available': False,
    }
    try:
        evidence = existing_evidence(model, args.batch1/model,
                                     args.raw/model, args.laurent/model)
        ensure_model_contract(evidence, model)
        record['source_sha256'] = evidence.sources_sha256
        record['validated_prior_stages'] = evidence.stage_statuses
        if args.preflight:
            record['status'] = 'SourceInputsValidated'
        else:
            native = run_wolfram(model, out, args.wolframscript, args.timeout)
            record['native_matching'] = native
            coeff=numeric_native(evidence,args.batch1/model,out,args)
            record['numeric_native_c5'] = coeff
            record['status'] = 'NativeNumericalC5Candidate_NotCompleteFullGeneralEFT'
            # A candidate low-scale observable is supplied for numerical
            # debugging only.  It cannot be used as the requested final fit.
            config=read_json(args.config)
            target=read_json(ROOT/config['fitting']['target'])
            c5=from_complex_json(coeff['symmetric_C5_candidate_GeV_inverse'],
                                 label='native C5')
            predictions=finish_smeft(evidence,c5,vev=float(config['vev_GeV']),target=target)
            put_json(out/'provisional_neutrino_prediction.json',predictions)
            record['candidate_smeft_neutrino_prediction_file']='provisional_neutrino_prediction.json'
            if args.fit:
                fitted=provisional_fit(evidence,model,out,args)
                record['provisional_fit_result_file']='provisional_fit_result.json'
                record['fit_status']=fitted['status']
            record['reason']=(
                'A numerical candidate was produced using the full native scalar threshold, '
                'but EFT1 Wilson RG, full finite soft subtraction, and native CG normalization '
                'are not verified. Not the requested completed full-general fit.')
    except (CompletionBlocked, ValueError, OSError, KeyError,
            TypeError, RuntimeError) as exc:
        record['status'] = 'Blocked_'+type(exc).__name__
        record['reason'] = str(exc)
        if args.debug:
            record['traceback_tail'] = traceback.format_exc(limit=8)
    put_json(out/'status.json', record)
    print(f'{model}: {record["status"]}', flush=True)
    if record.get('reason'):
        print('  '+record['reason'], flush=True)
    return record


def main(argv: list[str]|None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--all', action='store_true')
    g.add_argument('--model', choices=MODELS)
    p.add_argument('--batch1', type=Path, default=B1)
    p.add_argument('--raw', type=Path, default=B2)
    p.add_argument('--laurent', type=Path, default=LAURENT)
    p.add_argument('--output-root', type=Path, default=OUT)
    p.add_argument('--run-id', default='thesis_final_one_batch_v1')
    p.add_argument('--wolframscript', default='wolframscript')
    p.add_argument('--timeout', type=int, default=7200)
    p.add_argument('--config',type=Path,default=BENCHMARK)
    p.add_argument('--preflight', action='store_true')
    p.add_argument('--debug', action='store_true')
    p.add_argument('--fit',action='store_true',
                   help='Attempt provisional 18-real-Yukawa fit ONLY after a finite native C5 candidate')
    p.add_argument('--fit-maxeval',type=int,default=40)
    p.add_argument('--fit-bound',type=float,default=0.2)
    args = p.parse_args(argv)
    if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}', args.run_id):
        p.error('Invalid run ID')
    if not args.preflight and not (shutil.which(args.wolframscript) or
                                    Path(args.wolframscript).is_file()):
        p.error('wolframscript unavailable')
    if not EVALUATOR.is_file():
        p.error('Missing C5 evaluator '+str(EVALUATOR))
    if not NATIVE.is_file():
        p.error('Missing native matching source '+str(NATIVE))
    out = args.output_root / args.run_id
    if out.exists():
        p.error('Output run-id exists. Use a new run-id to preserve results')
    rows = [execute(model, out, args) for model in (MODELS if args.all else (args.model,))]
    put_json(out/'summary.json', {
        'format': FORMAT+'-summary',
        'models': [{k: x.get(k) for k in ('model','status','reason')} for x in rows],
        'stage': 'source-preflight' if args.preflight else 'native-two-threshold-attempt',
        'all_three_fermions_one_threshold': True,
        'both_scalars_one_threshold': True,
        'physical_full_general_neutrino_predictions_complete': False,
        'fitted_yukawas_complete': False,
    })
    if args.preflight:
        return 0 if all(x['status']=='SourceInputsValidated' for x in rows) else 2
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
