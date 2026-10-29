"""Stage 20: inspect/generate all 16 experimental general-quartic RGBeta bases.

--export-betas additionally calls the locally installed Wolfram/RGBeta backend.
No changes are made to the production T3RGBetaModel.wl or GitHub.
"""
from __future__ import annotations

import argparse
import os
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scripts.GenerateGeneralT3RGBetaBasis import (MANIFEST, OUTPUT,generate,audit_candidate_gauge_structures)


def verify_existing_rgbeta(registry:dict, source:Path):
    if not source.is_file():
        return {'status':'BlockedMissingProductionRGBetaSource','path':str(source)}
    text=source.read_text(encoding='utf-8')
    names=sorted({item['existing_name'] for body in registry['models'].values()
                  for item in body['preexisting_rgbeta_symbols']})
    missing=[n for n in names if not re.search(r'\bAddQuartic\s*\[\s*'+re.escape(n)+r'\s*,',text)]
    if missing:
        raise ValueError(f'Production RGBeta symbols not present: {missing}')
    return {'status':'SourceStructuralReuseVerified','distinct_reused_symbols':len(names)}


def check_export(path:Path,model:str,stage:str,expected:list[str]) -> dict:
    payload=json.loads(path.read_text(encoding='utf-8-sig'))
    if payload.get('status')!='Success':
        raise ValueError(f'{path}: RGBeta status is not Success')
    meta=payload.get('metadata',{})
    if meta.get('GeneralScalarModel')!=model or meta.get('GeneralScalarStage')!=stage:
        raise ValueError('General scalar RGBeta model/stage mismatch.')
    if sorted(meta.get('GeneralScalarNewBetaKeys',[]))!=sorted(expected):
        raise ValueError('RGBeta beta keys missing one or more new couplings.')
    betas=payload.get('report_betas')
    if not isinstance(betas,dict):
        raise ValueError('No report_betas exported.')
    missing=sorted(set(expected)-set(betas))
    if missing:
        raise ValueError(f'New beta keys missing from serialized report_betas: {missing}')
    if any(not isinstance(betas[k],str) or not betas[k].strip() for k in expected):
        raise ValueError('Empty or non-string extra quartic beta expression.')
    return {'status':'RGBetaBetasExportedCandidateBasis','model':model,'stage':stage,
            'new_complex_quartics':len(expected),'total_beta_keys':len(betas),
            'matchete_coupling_mapping_verified':False}


def normalise_model_selection(models: list[str] | None) -> list[str] | None:
    """Accept comma-separated PowerShell -Models values and repeated --model flags.

    For example, both ['B_-1,C_-1'] and ['B_-1', 'C_-1'] select the
    same two models. Preserve the requested order and avoid duplicate exports.
    None/[] means the caller's default (all models).
    """
    if not models:
        return None
    selected: list[str] = []
    for value in models:
        for part in str(value).split(','):
            name = part.strip()
            if not name:
                raise ValueError('Empty model name in --model/-Models selection.')
            if name not in selected:
                selected.append(name)
    return selected


def wolfram_integer_token(value: int) -> str:
    """Encode negative integers as mN for WolframScript's option parser.

    This matches T3RGBetaParseIntegerToken in T3RGBetaRunnerCommon.wl.
    Literal '-1' can be parsed as a WolframScript switch rather than a
    positional argument (even after -file), producing process exit code 2.
    """
    value = int(value)
    return f'm{-value}' if value < 0 else str(value)


def general_rgbeta_command(
    executable: str, script: Path, stage: str,
    d_s1: int, d_s2: int, d_f: int, alpha: int,
    model: str, destination: Path,
) -> list[str]:
    """Build the script + seven-argument CLI expected by RunGeneralT3RGBeta.wl."""
    if stage not in {'UV', 'EFT1'}:
        raise ValueError(f'Unsupported general RGBeta stage {stage!r}')
    return [
        str(executable), '-file', str(script), stage,
        wolfram_integer_token(d_s1),
        wolfram_integer_token(d_s2),
        wolfram_integer_token(d_f),
        wolfram_integer_token(alpha),
        str(model), str(destination),
    ]


def _process_log(*, command: list[str], returncode: int, stdout: str, stderr: str) -> str:
    """Keep the exact Wolfram invocation and its complete captured output."""
    import subprocess as _subprocess
    return (
        'Command: ' + _subprocess.list2cmdline(command) + '\n'
        + f'ExitCode: {returncode}\n'
        + 'STDOUT:\n' + stdout + '\nSTDERR:\n' + stderr + '\n'
    )


def _error_tail(stdout: str, stderr: str, max_lines: int = 30) -> str:
    lines = [line for line in (stdout + '\n' + stderr).splitlines() if line.strip()]
    return '\n'.join(lines[-max_lines:]) if lines else '(WolframScript returned no output)'


def run_stage20(*, export_betas:bool=False, models:list[str]|None=None,
                stage:str='both',output_dir:Path|None=None,check_gauge:bool=False) -> dict:
    out=OUTPUT if output_dir is None else Path(output_dir)
    manifest=json.loads(MANIFEST.read_text(encoding='utf-8'))
    registry=generate(manifest,out)
    source=ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
    reuse=verify_existing_rgbeta(registry,source)
    requested=normalise_model_selection(models) or sorted(registry['models'])
    unknown=sorted(set(requested)-set(registry['models']))
    if unknown:raise ValueError(f'Unknown model(s): {unknown}')
    report={'status':'CandidateBasisGenerated', 'coverage':registry['totals'],
            'source_reuse':reuse, 'selected_models':requested,
            'quartic_running_validated':False,
            'quadratic_mass_mixing_implemented':False,
            'physical_c5_matching_generalised':False,
            'reports':[],'gauge_invariance_checked':False,
            'gauge_check_method':'NotRequested',
            'native_rgbeta_checkinvariance_used':False}
    if check_gauge and not export_betas:
        raise ValueError('--check-gauge requires --export-betas')
    if check_gauge:
        audit=audit_candidate_gauge_structures(registry)
        report['candidate_gauge_structure_audit']=audit
        report['gauge_check_method']=audit['status']
        # Do not claim that RGBeta performed its rejected native check.
        report['gauge_invariance_checked']=False
    if export_betas:
        if reuse['status']!='SourceStructuralReuseVerified':
            raise RuntimeError('Cannot export without verifying current production RGBeta source.')
        exe=shutil.which('wolframscript')
        if exe is None:
            raise RuntimeError('wolframscript is not on PATH; cannot generate general quartic beta functions.')
        target=ROOT/'output/general_t3_stage20'
        target.mkdir(parents=True,exist_ok=True)
        for model in requested:
            body=registry['models'][model]
            d1,d2,df=body['dimensions']
            alpha=body['alpha']
            expected=[v['name'] for v in body['new_rgbeta_symbols']]
            for name in (['UV','EFT1'] if stage=='both' else [stage]):
                dest=target/f'{model}_{name.lower()}_rgbeta.json'
                command=general_rgbeta_command(
                    exe, ROOT/'RGE/running/rgbeta/RunGeneralT3RGBeta.wl',
                    name,d1,d2,df,alpha,model,dest,
                )
                print(f'Stage20 {model} {name}: exporting {len(expected)} new quartic betas',flush=True)
                proc=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,
                    env={k:v for k,v in os.environ.items()
                         if k != 'T3_GENERAL_RGBETA_VERIFY_GAUGE'})
                log_path=target/f'{model}_{name.lower()}_rgbeta.log'
                log_path.write_text(
                    _process_log(command=command, returncode=proc.returncode,
                                 stdout=proc.stdout, stderr=proc.stderr),
                    encoding='utf-8',
                )
                if proc.returncode:
                    raise RuntimeError(
                        f'{model} {name} RGBeta failed (code {proc.returncode}).\n'
                        f'Log: {log_path}\n'
                        f'Wolfram output (last lines):\n'
                        + _error_tail(proc.stdout,proc.stderr)
                    )
                report['reports'].append(check_export(dest,model,name,expected))
        report['status']='RGBetaCandidateBetasExportedForSelectedModels'
        report['quartic_running_validated']=False  # Beta export alone is not ODE or matched-basis validation.
    summary=ROOT/'output/general_t3_stage20/stage20_summary.json'
    summary.parent.mkdir(parents=True,exist_ok=True)
    summary.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--export-betas',action='store_true')
    ap.add_argument('--check-gauge',action='store_true',
        help='Check U(1) charge, Z2 parity, field orientation and SU(2) candidate rank in Python; RGBeta native CheckInvariance is unsupported here.')
    ap.add_argument('--model',action='append',help='e.g. A_-2; repeat to select multiple models')
    ap.add_argument('--stage',choices=['UV','EFT1','both'],default='both')
    a=ap.parse_args()
    summary=run_stage20(export_betas=a.export_betas,models=a.model,stage=a.stage,
                        check_gauge=a.check_gauge)
    print('Stage20:',summary['status'])
    print('Manifest couplings:',summary['coverage']['matchete_extra_quartics'])
    print('Reused:',summary['coverage']['rgbeta_preexisting_structural'])
    print('New RGBeta vertices:',summary['coverage']['rgbeta_new_vertices'])
    print('Extra quadratic masses still excluded:',summary['coverage']['quadratic_couplings_blocked'])
    print('Beta exports:',len(summary['reports']))
    print('Report: output/general_t3_stage20/stage20_summary.json')
    return 0

if __name__=='__main__':raise SystemExit(main())
