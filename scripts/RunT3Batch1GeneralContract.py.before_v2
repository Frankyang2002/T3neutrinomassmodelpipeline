"""Batch 1 all-16 general T3 source, native-run, and basis contract.

Additive: never replaces an existing production physics module. Uses the
CURRENT LOCAL repository, not handover text, for input. Existing source
folders are never overwritten. Writes output/batch1_general_contract/RUN_ID.

Operational goal: obtain the all-fermions-first symbolic LLSS inventory and
UV/EFT1 native one-loop renormalisable beta evidence for every supported
model, then report unambiguous remaining basis/matching blockers.

This is NOT itself physical GammaCC->Weyl normalization, an RGBeta coupling
map or a complete dimension-five ADM; these require explicit native tensors.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
from typing import Any

from Numerical.basis.T3Batch1TensorMap import (
    solve_dictionary, scan_hermitian_phases, TensorMapError)

ROOT = Path(__file__).resolve().parents[1]
CATALOG = {
    'A': ((1, 3, 2), (-4, -2, 0)),
    'B': ((2, 2, 1), (-3, -1, 1)),
    'C': ((2, 2, 3), (-3, -1, 1)),
    'D': ((3, 1, 2), (-2, 0, 2)),
    'E': ((3, 3, 2), (-4, -2, 0, 2)),
}
MODELS = tuple(f'{k}_{n:+d}' for k, (_, alphas) in CATALOG.items() for n in alphas)
MIXED = frozenset(('B_-1', 'C_-1'))
SELF = frozenset(('A_-2', 'A_+0', 'D_-2', 'D_+0', 'E_-2', 'E_+0'))
SOURCE_FILES = (
    'Lagrangian/model/T3ModelCatalog.wl',
    'Lagrangian/model/T3Fields.wl',
    'Lagrangian/model/LagrangianBuilder.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'Lagrangian/interactions/ScalarPotential.wl',
    'Lagrangian/interactions/GeneralScalarOperators.wl',
    'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl',
    'Lagrangian/RunT3GeneralFThreshold.wl',
    'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl',
    'RGE/running/rgbeta/T3RGBetaModel.wl',
    'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
    'RGE/general/generated/general_rgbeta_quartic_registry.json',
    'RGE/general/WilsonTensorRGE.py',
    'RGE/running/intermediate/GeneralEFT1SevenTermFlow.py',
    'Numerical/running/T3EFT1MasterSourceBridge.py',
)
# The current GitHub main as actually read, NOT a claim about the local tree.
PUBLIC_BASELINE = 'ed8f6ba13c9a803e698ac155d24d8c025c53bc79'
PUBLIC_BLOB_REFERENCE = {
    'Lagrangian/interactions/ScalarPotential.wl': '3dcaf8f49942aa4ce7c3b15687143e43eab38985',
    'Lagrangian/model/LagrangianBuilder.wl': '98d549cda503b751380dc10c4a56212e5312ccd8',
    'Lagrangian/model/T3ModelCatalog.wl': '9a6c47e10fc83d9480f08353682151180c9fbb4a',
    'RGE/running/intermediate/GeneralEFT1SevenTermFlow.py': '571308d2608aba0ce233d87efca0fda4694a6210',
    'RGE/general/generated/general_rgbeta_quartic_registry.json': 'bf674314d60b6b4f8fdb7b9a5d200038611cbdf7',
    'RGE/running/rgbeta/T3RGBetaModel.wl': 'bc9500803f2324ad6504e71d769cb75fbbde9b07',
    'Numerical/running/T3EFT1MasterSourceBridge.py': '5377f82163e60f199b09e598b4e5e2b50790159b',
}
EXTRA_YUKAWAS = SELF
EXPECTED_A4_QUARTICS = frozenset((
    'λ', 'lambdaS1Inv1', 'lambdaS2Inv1', 'lambdaS2Inv2',
    'lambdaH1Inv1', 'lambdaH2Inv1', 'lambdaH2Inv2',
    'lambda12Inv1', 'lambdaT3',
))


def die(message: str) -> None:
    raise RuntimeError(message)


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def utf8_json(path: Path) -> dict[str, Any]:
    if not path.is_file(): die(f'Missing JSON {path}')
    # Accommodate local legacy PowerShell UTF-16 files without overwriting them.
    data = path.read_bytes()
    if data.startswith((b'\xff\xfe', b'\xfe\xff')):
        content = data.decode('utf-16')
    else:
        try: content = data.decode('utf-8-sig')
        except UnicodeDecodeError: content = data.decode('utf-16')
    obj = json.loads(content)
    if not isinstance(obj, dict): die(f'Expected JSON object {path}')
    return obj


def read_log(path: Path) -> str:
    data = path.read_bytes()
    if data.startswith((b'\xff\xfe', b'\xfe\xff')): return data.decode('utf-16')
    try: return data.decode('utf-8-sig')
    except UnicodeDecodeError: return data.decode('utf-16')


def git(args: list[str]) -> str | None:
    try:
        p = subprocess.run(['git', *args], cwd=ROOT, stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, timeout=20, check=True)
        return p.stdout.decode('utf-8', 'replace').strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        return None


def source_manifest() -> dict[str, Any]:
    local = {}
    missing = []
    for relative in SOURCE_FILES:
        path = ROOT/relative
        if path.is_file():
            local[relative] = {'sha256': sha(path), 'size_bytes': path.stat().st_size,
                               'git_blob_sha': git(['hash-object', relative]),
                               'same_as_public_baseline_blob': (
                                   git(['hash-object', relative]) == PUBLIC_BLOB_REFERENCE[relative]
                                   if relative in PUBLIC_BLOB_REFERENCE else None)}
        else: missing.append(relative)
    potential_text = (ROOT/'Lagrangian/interactions/ScalarPotential.wl').read_text(
        encoding='utf-8-sig') if (ROOT/'Lagrangian/interactions/ScalarPotential.wl').is_file() else ''
    portal_corrected = ('T3HiggsPortalHermitianPhase' in potential_text and
                        'phase coupling[] fieldProduct' in potential_text)
    return {'baseline_public_main_sha': PUBLIC_BASELINE,
            'a4_antihermitian_portal_source_correction_present': portal_corrected,
            'local_head': git(['rev-parse','HEAD']),
            'local_branch': git(['branch','--show-current']),
            'git_status_porcelain': git(['status','--short']),
            'files': local, 'missing_source_files': missing}


def validate_registry(reg: dict[str, Any]) -> dict[str, Any]:
    models = reg.get('models', {})
    if set(models) != set(MODELS):
        die(f'Registry model set mismatch: missing {sorted(set(MODELS)-set(models))}, '
            f'extra {sorted(set(models)-set(MODELS))}')
    counts = {'models': 0, 'additional_quartics': 0, 'new_symbols': 0,
              'reused_symbols': 0}
    rows = {}
    for model in MODELS:
        family, aa = model.split('_')
        alpha = int(aa)
        entry = models[model]
        dims = CATALOG[family][0]
        if entry.get('dimensions') != list(dims) or entry.get('alpha') != alpha:
            die(f'{model}: registry dimensions/hypercharge disagree with live catalog')
        new = entry.get('new_rgbeta_symbols', [])
        reused = entry.get('preexisting_rgbeta_symbols', [])
        extra = entry.get('total_additional_quartics')
        if type(extra) is not int or extra != len(new)+len(reused):
            die(f'{model}: inconsistent extra quartic counts')
        for item in new+reused:
            count = item.get('counts', [])
            if len(count) != 6 or sum(count) != 4 or any(type(x) is not int for x in count):
                die(f'{model}: malformed quartic structure')
            y2 = (count[0]-count[1] + alpha*(count[2]-count[3]) +
                  (alpha+2)*(count[4]-count[5]))
            if y2 != 0 or (sum(count[2:]) % 2):
                die(f'{model}: invalid U(1)Y or Z2 symmetry in extra quartic {count}')
        for item in new:
            extension = ROOT/f'RGE/general/generated/{model}_extra_quartics.wl'
            if not extension.is_file(): die(f'Missing native RGBeta extension {extension}')
            if item['name'] not in extension.read_text(encoding='utf-8-sig'):
                die(f'{model}: missing registered quartic {item["name"]}')
        counts['models'] += 1
        counts['additional_quartics'] += extra
        counts['new_symbols'] += len(new)
        counts['reused_symbols'] += len(reused)
        rows[model] = {'dimensions': list(dims), 'alpha': alpha,
                       'hypercharges_2y': {'S1': alpha, 'S2': alpha+2, 'F': alpha+1},
                       'extra_quartics': extra, 'new_quartics': len(new),
                       'reused_quartics': len(reused),
                       'extra_yukawa_expected': model in EXTRA_YUKAWAS,
                       'mixed_mass_basis': model in MIXED,
                       'neutral_self_bilinear': model in SELF}
    if counts != {'models':16,'additional_quartics':50,
                  'new_symbols':38,'reused_symbols':12}:
        die(f'Live registry changed relative to audited expected totals: {counts}')
    return {'totals': counts, 'models': rows}


def wolfram_integer_token(value: int) -> str:
    """WolframScript-safe signed integer token.

    A raw '-4' following -file is consumed as a WolframScript option on the
    user's Windows installation; the native RGBeta reader already supports
    the 'm4' representation for -4 (T3RGBetaParseIntegerToken).
    """
    return f'm{-value}' if value < 0 else str(value)


def native_commands(model: str, root_out: Path, loop_order: int) -> list[tuple[str, list[str], Path]]:
    fam, alpha_s = model.split('_')
    dims = CATALOG[fam][0]
    args = [*(str(x) for x in dims), wolfram_integer_token(int(alpha_s)), model]
    tag = model.replace('+','p').replace('-','m')
    base = ['wolframscript','-file']
    commands = [
        ('rgbeta_uv', base+['RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
                             'UV',*args,str(root_out/'rgbeta'/tag/'uv.json')],
            root_out/'logs'/tag/'rgbeta_uv.log'),
        ('rgbeta_eft1', base+['RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl',
                               'EFT1',*args,str(root_out/'rgbeta'/tag/'eft1.json')],
            root_out/'logs'/tag/'rgbeta_eft1.log'),
    ]
    if model not in MIXED:
        commands.append(('f_threshold', base+['Lagrangian/RunT3GeneralFThreshold.wl',
                           str(root_out/'matchete'/tag),model,str(loop_order)],
                         root_out/'logs'/tag/'f_threshold.log'))
    return commands


def run_local(label: str, command: list[str], log: Path,
              *, timeout: int, resume: bool) -> dict[str, Any]:
    log.parent.mkdir(parents=True,exist_ok=True)
    if log.exists() and not resume:
        return {'stage':label,'status':'Blocked_ExistingLogRequiresResume','command':command}
    exe = command[0]
    try:
        process = subprocess.run(command,cwd=ROOT,env={**os.environ,
                'T3_EXPERIMENTAL_GENERAL_FERMIONS':'1',
                'T3_EXPERIMENTAL_GENERAL_QUARTICS':'1',
                'T3_EXPERIMENTAL_GENERAL_QUADRATICS':'1',
                'T3_EXPERIMENTAL_FLAVOR_MASSES':'1'},
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
            timeout=timeout, check=False)
        # Force UTF-8 output on disk; stdout may contain Windows ACP bytes.
        try: output=process.stdout.decode('utf-8')
        except UnicodeDecodeError: output=process.stdout.decode('cp1252','replace')
        log.write_text(output,encoding='utf-8')
        return {'stage':label,'status':'Completed' if process.returncode == 0 else 'NativeFailed',
                'returncode':process.returncode,'command':command,
                'log':str(log.relative_to(ROOT)) if log.is_relative_to(ROOT) else str(log)}
    except (FileNotFoundError,subprocess.TimeoutExpired,OSError) as e:
        log.write_text(f'{type(e).__name__}: {e}\n',encoding='utf-8')
        return {'stage':label,'status':'NativeFailed','error':str(e),'command':command}


def inspect_f(model: str, root_out: Path, optional_source: Path | None) -> dict[str, Any]:
    tag = model.replace('+','p').replace('-','m')
    folder = optional_source if optional_source and model == 'A_-4' else root_out/'matchete'/tag
    result: dict[str, Any] = {'source_directory':str(folder),'source_exists':folder.is_dir()}
    if not folder.is_dir():
        result['status']='Blocked_MissingFThreshold'; return result
    statusfile = folder/'native_status.json'
    if not statusfile.exists():
        result['status']='Blocked_MissingNativeStatus'; return result
    doc = utf8_json(statusfile)
    result['native_status']=doc.get('status')
    result['requested_general']=doc.get('operator_basis_is_full_general_requested') is True
    result['integrated_fields']=None
    seeds = ('data/eft1_after_F_wilson_seed.json',
             'data/eft1_after_F_scalar_quartic_seed.json')
    for seed in seeds:
        p=folder/seed
        if p.exists():
            q=utf8_json(p)
            result[seed]={'sha256':sha(p),'IntegratedFields':q.get('IntegratedFields'),
                          'ActiveHeavyFields':q.get('ActiveHeavyFields')}
        else: result[seed]={'missing':True}
    valid_seed = all(
        result[seed].get('IntegratedFields')==['F'] and
        result[seed].get('ActiveHeavyFields')==['S1','S2']
        for seed in seeds)
    result['fermions_out_scalars_kept']=valid_seed
    if model=='A_-4':
        path=folder/'data/eft1_after_F_tree_inputform.txt'
        logfile=folder/'native.log'
        if path.exists() and logfile.exists():
            source=read_log(path); trace=read_log(logfile)
            present=set(re.findall(r'Coupling\[([^,\]]+)',source))
            missing=sorted(EXPECTED_A4_QUARTICS-present)
            checks=all(re.search(r'Checking candidate: HiggsPortal2Inv'+i+
                 r'\s+result: InputForm\[True\]',trace) for i in ('1','2'))
            result['a4_portal']={'both_validated':checks,
                'inv1_survives': 'lambdaH2Inv1' in present,
                'inv2_survives': 'lambdaH2Inv2' in present,
                'matched_quartics_missing':missing}
        else: result['a4_portal']={'status':'Blocked_MissingPortalRegressionEvidence'}
    result['status']=('FThresholdInventoryAvailable' if valid_seed and
               result['requested_general'] and doc.get('status') == 'FThresholdTreeWilsonExported'
               else 'Blocked_FThresholdInventoryNotEstablished')
    return result


def inspect_beta(model: str, root_out: Path) -> dict[str, Any]:
    """Inspect native outputs; completed subprocess alone is insufficient."""
    tag = model.replace('+','p').replace('-','m')
    ans = {}
    for name in ('uv','eft1'):
        p=root_out/'rgbeta'/tag/f'{name}.json'
        if not p.is_file():
            ans[name]={'status':'Missing','candidate_ready':False}
            continue
        try:
            doc=utf8_json(p)
        except (OSError,ValueError,UnicodeError) as exc:
            ans[name]={'status':'InvalidJSON','error':str(exc),'candidate_ready':False}
            continue
        betas=doc.get('report_betas',{})
        keys=set(betas) if isinstance(betas,dict) else set()
        remove=sorted(keys & {'MF','y1','y2','y1Conj','y2Conj'})
        metadata=doc.get('metadata',{})
        metadata=metadata if isinstance(metadata,dict) else {}
        correct_identity=(metadata.get('GeneralScalarModel') == model and
                          metadata.get('GeneralScalarStage') == name.upper())
        good_status=doc.get('status') == 'SuccessCandidateMassBeta'
        ready=bool(keys) and good_status and correct_identity and (name!='eft1' or not remove)
        ans[name]={'status':doc.get('status','Unknown'),
                   'sha256':sha(p),'beta_key_count':len(keys),
                   'fermion_keys_unexpected_in_eft1':remove if name=='eft1' else [],
                   'model_metadata':metadata,
                   'model_and_stage_identity_ok':correct_identity,
                   'candidate_ready':ready}
    return ans


def archive_prior(path: Path) -> Path | None:
    """Preserve a failed log/invalid output before an explicit repair."""
    if not path.is_file():
        return None
    k=1
    while True:
        candidate=path.with_name(f'{path.name}.before_batch1_repair_{k}')
        if not candidate.exists():
            shutil.copy2(path,candidate)
            return candidate
        k+=1


def repair_rgbeta(model: str, root_out: Path, timeout: int) -> list[dict[str,Any]]:
    """Only rerun missing/invalid native beta outputs; never rerun F threshold.

    Existing successful UV/EFT1 outputs and both Matchete scalar/fermion seed
    directories remain byte-for-byte untouched.
    """
    before=inspect_beta(model,root_out)
    results=[]
    for stage,command,log in native_commands(model,root_out,0):
        if stage not in ('rgbeta_uv','rgbeta_eft1'):
            continue
        name='uv' if stage=='rgbeta_uv' else 'eft1'
        if before[name]['candidate_ready']:
            results.append({'stage':stage,'status':'ReusedExistingCandidate',
                            'sha256':before[name]['sha256']})
            continue
        dest=Path(command[-1])
        prior_log=archive_prior(log)
        prior_json=archive_prior(dest)
        result=dict(run_local(stage,command,log,timeout=timeout,resume=True))
        result['repair_reason']=before[name]['status']
        result['archived_old_log']=str(prior_log) if prior_log else None
        result['archived_old_json']=str(prior_json) if prior_json else None
        results.append(result)
    return results


def check_repair_provenance(out: Path, current: dict[str, Any]) -> None:
    """Refuse to mix old native results with changed physics source files."""
    old=out/'provenance.json'
    if not old.is_file() or not (out/'summary.json').is_file():
        die('Repair requires an existing Batch 1 run with summary.json and provenance.json')
    prev=utf8_json(old)
    previous_files=prev.get('files',{})
    current_files=current.get('files',{})
    changes=[p for p in SOURCE_FILES if
             previous_files.get(p,{}).get('sha256') != current_files.get(p,{}).get('sha256')]
    if changes:
        die('Source changed since original Batch 1 native run; cannot safely reuse results. '
            'Use a new run-id. Changed files: '+', '.join(changes))


def restore_a4_source_from_report(out: Path) -> Path | None:
    """Preserve the corrected user-provided Matchete A(-4) seed on repair."""
    path=out/'models'/'A_m4.json'
    if not path.is_file():
        return None
    previous=utf8_json(path)
    src=previous.get('f_threshold',{}).get('source_directory')
    if not isinstance(src,str) or not src:
        return None
    value=Path(src)
    return value if value.is_dir() else None


def _extract_tensors(doc: dict[str, Any], model: str) -> list[dict[str, Any]]:
    node=doc.get('models',{}).get(model,{})
    result=node.get('tensors', [])
    if not isinstance(result,list): raise TensorMapError('Tensor export must have a list')
    formatted=[]
    for q in result:
        rec=dict(q)
        if doc.get('format')=='stage22-exact-tensors-v1':
            # Exporter explicitly says PlusHc[contribution], with factor 1.
            rec.setdefault('add_hc',True)
        formatted.append(rec)
    return formatted


def optional_tensor_map(model: str, mat_path: Path | None,
                        rgb_path: Path | None) -> dict[str, Any]:
    if not mat_path or not rgb_path:
        return {'status':'Blocked_ExplicitNativeComponentTensorPairRequired',
                'reason':'RGBeta native declaration metadata is not a component tensor'}
    if not mat_path.is_file() or not rgb_path.is_file():
        return {'status':'Blocked_TensorExportMissing'}
    ma, rb = utf8_json(mat_path), utf8_json(rgb_path)
    if rb.get('status')=='NativeDeclarationMetadataOnly':
        return {'status':'Blocked_RGBetaDeclarationsAreNotComponentTensors'}
    try:
        return solve_dictionary(_extract_tensors(ma,model),_extract_tensors(rb,model))
    except (TensorMapError, ValueError, np.linalg.LinAlgError) as e:  # noqa: F821
        return {'status':'Blocked_TensorMapFailed','error':str(e)}


def one_model(model: str, root_out: Path, source: dict[str, Any],
              registry: dict[str, Any], args: argparse.Namespace,
              stage_runs: list[dict[str, Any]]) -> dict[str, Any]:
    # Neither an RGBeta beta export nor Matchete CheckLagrangian alone proves
    # physical equivalence of their native tensor coordinates.
    fm = ({'status':'Blocked_MixedDoubletMassBasisRequiresNativeEigenDriver'}
          if model in MIXED else inspect_f(model,root_out,args.a4_source))
    betas = inspect_beta(model,root_out)
    dictionary = optional_tensor_map(model,args.matchete_tensors,args.rgbeta_components)
    hermiticity = {'status':'UnavailableWithoutNativeTensorComponents'}
    if args.matchete_tensors and args.matchete_tensors.is_file():
        try:
            rows = _extract_tensors(utf8_json(args.matchete_tensors),model)
            hermiticity = {'status':'NativeTensorNumericalPhaseAudit',
                           'rows':scan_hermitian_phases(rows) if rows else []}
        except (TensorMapError,ValueError) as exc:
            hermiticity = {'status':'Blocked_HermiticityAuditFailed','error':str(exc)}
    blocking = []
    if fm.get('status') != 'FThresholdInventoryAvailable':
        blocking.append(fm['status'])
    if not source['a4_antihermitian_portal_source_correction_present']:
        blocking.append('Blocked_LocalScalarPotentialLacksHermitianPhaseCorrection')
    if model in SELF: blocking.append('Conditional_NeutralSelfBilinearComplexPhaseUnresolved')
    if model in MIXED: blocking.append('Conditional_MixedMassBasisRequiresDedicatedEigenDriver')
    if fm.get('a4_portal') is not None:
        a4=fm['a4_portal']
        if not (a4.get('both_validated') and a4.get('inv1_survives') and
                a4.get('inv2_survives') and not a4.get('matched_quartics_missing', ['missing'])):
            blocking.append('Blocked_A4PortalRegression')
    for stage in ('uv','eft1'):
        if not betas[stage].get('candidate_ready'):
            blocking.append(f'Blocked_{stage.upper()}RGBetaCandidateMissing')
    if dictionary['status'] != 'NumericalPolynomialMapComplete':
        blocking.append(dictionary['status'])
    blocking += ['Blocked_SignedGammaCCWeylConversionUnresolved',
                 'Blocked_DimensionFiveADMDisputed',
                 'Blocked_FOneLoopHardMatchingUncertified']
    return {'model':model,'representation':registry['models'][model],
            'local_source_sha256':{k:v['sha256'] for k,v in source['files'].items()},
            'native_stage_runs':stage_runs, 'f_threshold':fm,'rgbeta':betas,
            'extra_quartic_dictionary':dictionary,
            'native_quartic_hermiticity':hermiticity,
            'spinor_normalization':{'status':'NotDerived','factor':None},
            'dimension_five_adm':{'status':'Disputed','seven_term_physical':False},
            'contract_status':'OperationalCandidate' if not blocking else 'ConditionalOrBlocked',
            'physical_C5_claim_permitted':False,
            'blockers':blocking}


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id',default='batch1_v1')
    ap.add_argument('--models',nargs='+',choices=MODELS,default=list(MODELS))
    mode=ap.add_mutually_exclusive_group()
    mode.add_argument('--run-native',action='store_true',
                    help='Actually run RGBeta UV/EFT1 and 14 F-only Matchete jobs')
    mode.add_argument('--repair-rgbeta',action='store_true',
                    help='Repair only missing/invalid UV/EFT1 beta outputs from an existing run-id')
    ap.add_argument('--resume',action='store_true',
                    help='Explicitly allow replacing output logs and result files in run-id')
    ap.add_argument('--native-timeout',type=int,default=7200)
    ap.add_argument('--f-loop-order',type=int,choices=(0,1),default=0)
    ap.add_argument('--a4-source',type=Path,
                    help='Optional corrected existing A(-4) F-threshold directory')
    ap.add_argument('--matchete-tensors',type=Path,
                    help='Native scalar component tensor JSON (extras only)')
    ap.add_argument('--rgbeta-components',type=Path,
                    help='Same-basis full native RGBeta components (NOT declarations)')
    ap.add_argument('--export-stage22',action='store_true',
                    help='Run existing Matchete tensor and RGBeta declaration exporters')
    return ap.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args=parse_args(argv)
    if any(c in args.run_id for c in '\\/:' ) or not re.fullmatch(r'[A-Za-z0-9_.-]+',args.run_id):
        die('run-id must be one simple filename')
    reg=validate_registry(utf8_json(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'))
    sources=source_manifest()
    if sources['missing_source_files']: die('Missing source files: '+str(sources['missing_source_files']))
    out=ROOT/'output/batch1_general_contract'/args.run_id
    out.mkdir(parents=True,exist_ok=True)
    if args.run_native and not args.resume and (out/'summary.json').exists():
        die('Existing run-id report: choose a NEW run-id or explicitly pass --resume')
    if args.repair_rgbeta:
        if set(args.models)!=set(MODELS):
            die('--repair-rgbeta regenerates the complete all-16 report; omit --models')
        if args.export_stage22:
            die('--export-stage22 is not part of the RGBeta repair')
        check_repair_provenance(out,sources)
        if args.a4_source is None:
            args.a4_source=restore_a4_source_from_report(out)
        if args.a4_source is None:
            die('Cannot recover corrected A(-4) source; pass --a4-source')
        # Preserve previously emitted evidence before overwriting reports.
        archive_prior(out/'summary.json')
        archive_prior(out/'provenance.json')
        archive_prior(out/'registry_audit.json')
    if (args.run_native or args.repair_rgbeta) and not sources['a4_antihermitian_portal_source_correction_present']:
        die('Refusing native run: local ScalarPotential.wl does not include the corrected portal phase')
    (out/'provenance.json').write_text(json.dumps(sources,indent=2)+'\n',encoding='utf-8')
    (out/'registry_audit.json').write_text(json.dumps(reg,indent=2)+'\n',encoding='utf-8')
    if args.export_stage22:
        tasks=[('matchete_extra_tensors', ['wolframscript','-file',
                 'scripts/ExportGeneralT3ScalarTensorMatchete.wl',
                 str(out/'native/matchete_tensors.json')]),
               ('rgbeta_native_declarations',['wolframscript','-file',
                 'scripts/ExportGeneralT3ScalarTensorRGBetaDeclarations.wl',
                 str(out/'native/rgbeta_declarations.json')])]
        for label, command in tasks:
            Path(command[-1]).parent.mkdir(parents=True,exist_ok=True)
            run_local(label,command,out/'logs'/f'{label}.log',
                      timeout=args.native_timeout,resume=args.resume)
        if args.matchete_tensors is None:
            args.matchete_tensors=out/'native/matchete_tensors.json'
    rows={}
    for model in args.models:
        stage_runs=[]
        if args.repair_rgbeta:
            stage_runs=repair_rgbeta(model,out,args.native_timeout)
        elif args.run_native:
            for label,command,log in native_commands(model,out,args.f_loop_order):
                Path(command[-1]).parent.mkdir(parents=True,exist_ok=True)
                stage_runs.append(run_local(label,command,log,
                    timeout=args.native_timeout,resume=args.resume))
        rec=one_model(model,out,sources,reg,args,stage_runs)
        rows[model]=rec
        tag=model.replace('+','p').replace('-','m')
        p=out/'models'/f'{tag}.json';p.parent.mkdir(parents=True,exist_ok=True)
        if args.repair_rgbeta:
            archive_prior(p)
        p.write_text(json.dumps(rec,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(model,rec['contract_status'],
              f"{len(rec['blockers'])} unresolved conditions",flush=True)
    summary={'format':'t3-batch1-general-contract-v1',
         'created_utc':datetime.now(timezone.utc).isoformat(),
         'source_baseline_main':PUBLIC_BASELINE,
         'local_head':sources['local_head'],
         'models_requested':list(args.models),
         'model_count':len(rows),
         'registry_totals':reg['totals'],
         'all_16_requested':set(args.models)==set(MODELS),
         'model_statuses':{m:{'status':v['contract_status'],
                 'physical_C5_claim_permitted':False,'blockers':v['blockers']}
                 for m,v in rows.items()},
         'physics_workflow_order':['UV','all_fermions_at_MF',
              'EFT1_S1_and_S2','all_scalars_simultaneously_at_MS','SMEFT'],
         'complete_two_threshold_pipeline_claimed':False}
    path=out/'summary.json'
    path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Batch 1 report:',path)
    # Structural report succeeds even when physical mapping remains blocked.
    return 0


if __name__=='__main__':
    raise SystemExit(main())
