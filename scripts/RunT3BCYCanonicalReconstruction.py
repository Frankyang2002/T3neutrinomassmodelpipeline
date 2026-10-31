"""Experimental independent Y-only EFT1 canonical reconstruction.

This probes a source-defined full renormalisable EFT1 plus the EXACT saved
fermion-threshold LLSS, because the serialized EOMEFT fails Hermiticity.
It deliberately does NOT claim equality to the canonical Matchete EOMEFT,
renormalised scalar matching, or a physical neutrino-mass prediction.

python -u -X utf8 -m scripts.RunT3BCYCanonicalReconstruction --all --check-only
python -u -X utf8 -m scripts.RunT3BCYCanonicalReconstruction --all --match
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl'
PARAM = ROOT / 'output/general_t3_mass_eigen_matchete'
TREE = PARAM / 'eft1_fermion_first_tree'
DEFAULT_OUT = PARAM / 'eft1_scalar_Y_only_canonical_probe'
MODELS = ('B_-1','C_-1')
CHECK_STATUS = 'IndependentCanonicalFullEFTCheckPassedEOMEquivalenceUnproved'
MATCH_STATUS = 'YOnlyIndependentCanonicalRawEFT2NotCertified'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_paths(model: str) -> dict[str,Path]:
    if model not in MODELS:
        raise ValueError('Unsupported B/C model')
    tag=model.replace('+','p').replace('-','m')
    return {
        'wolfram':SCRIPT,
        'parameters':PARAM/(tag+'_final_parameters.json'),
        'tree':TREE/(tag+'_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt'),
    }


def preflight(model: str) -> dict:
    paths=source_paths(model)
    for key,path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f'{model}: required {key} missing: {path}')
    tree=paths['tree'].read_text(encoding='utf-8')
    if len(tree)<1000 or any(t not in tree for t in
        ('Field[NewScalarX','Field[NewScalarY','GammaCC')) or 'Field[NewFermion' in tree:
        raise RuntimeError(f'{model}: canonical EOMEFT inventory changed')
    raw=json.loads(paths['parameters'].read_text(encoding='utf-8'))
    if raw.get('model')!=model:
        raise RuntimeError(f'{model}: mass-basis parameter file model mismatch')
    return {name:sha(path) for name,path in paths.items()}


def run_one(model:str, *, mode:str, outdir:Path, wolframscript:str,
            timeout_sec:int, overwrite:bool=False) -> dict:
    upstream=preflight(model)
    tag=model.replace('+','p').replace('-','m')
    stem=f'{tag}_final_y_only_canonical_{mode}'
    outdir.mkdir(parents=True,exist_ok=True)
    report_path=outdir/(stem+'.json')
    provenance=outdir/(stem+'.provenance.json')
    log=outdir/(stem+'.log')
    if any(p.exists() for p in (report_path,provenance,log)) and not overwrite:
        raise RuntimeError(f'{model}: output exists: {report_path}; use --overwrite for a rerun')
    if mode=='match':
        cert_path=outdir/(f'{tag}_final_y_only_canonical_check.provenance.json')
        if not cert_path.is_file():
            raise RuntimeError(f'{model}: run --check-only first; no certificate: {cert_path}')
        prior=json.loads(cert_path.read_text(encoding='utf-8'))
        if prior.get('source_sha256')!=upstream or prior.get('status')!=CHECK_STATUS:
            raise RuntimeError(f'{model}: source files changed since successful check-only run')
    cmd=[wolframscript,'-file',str(SCRIPT),
         str(source_paths(model)['parameters']),
         str(source_paths(model)['tree']),str(report_path)]
    env=os.environ.copy()
    env.pop('T3_EFT1_LEADING_C5_ONLY',None)
    env.pop('T3_EXPERIMENTAL_GENERAL_QUARTICS',None)
    env['T3_BC_Y_INDEPENDENT_CANONICAL']='1'
    env['T3_BC_Y_INDEPENDENT_CHECK_ONLY']='1' if mode=='check' else '0'
    try:
        proc=subprocess.run(cmd,cwd=ROOT,env=env,check=False,
                            capture_output=True,text=True,encoding='utf-8',
                            errors='replace',timeout=timeout_sec)
        log.write_text('command: '+' '.join(cmd)+'\nexit code: '+str(proc.returncode)+
            '\nstdout:\n'+proc.stdout+'\nstderr:\n'+proc.stderr,encoding='utf-8')
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'{model}: Matchete timed out; {mode} incomplete') from exc
    result=json.loads(report_path.read_text(encoding='utf-8')) if report_path.is_file() else {}
    expected=CHECK_STATUS if mode=='check' else MATCH_STATUS
    reconstruction=outdir/(stem+'_reconstruction_audit.json')
    details=json.loads(reconstruction.read_text(encoding='utf-8')) if reconstruction.is_file() else {}
    status=result.get('status','NO_REPORT')
    certificate={
        'format':'t3-bc-y-only-canonical-reconstruction-experimental-v1',
        'model':model,'mode':mode,'status':status,'returncode':int(proc.returncode),
        'source_sha256':upstream,
        'saved_eomeft_algebraic_equality_proven':False,
        'physical_matching_certified':False,
        'independent_matchete_check_passed':details.get('independent_matchete_check_passed',False),
        'independent_hermiticity':details.get('independent_hermiticity','Unavailable'),
        'raw_matchete_matching_completed':mode=='match' and status==MATCH_STATUS and proc.returncode==0,
        'report':report_path.name,'log':log.name,
        'reconstruction_audit':reconstruction.name if reconstruction.exists() else None,
    }
    provenance.write_text(json.dumps(certificate,indent=2)+'\n',encoding='utf-8')
    print(f'{model}: {mode} status={status}; independent Lagrangian check='
          f'{certificate["independent_matchete_check_passed"]}; EOMEFT equivalence=UNPROVEN',
          flush=True)
    if proc.returncode or status!=expected or not certificate['independent_matchete_check_passed']:
        raise RuntimeError(f'{model}: check/match did not pass. See {log} and {reconstruction}')
    if result.get('physical_matching_performed') is not False or result.get('physical_sequential_C5_computed') is not False:
        raise RuntimeError(f'{model}: unexpected physical-matching flag')
    return certificate


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--all',action='store_true')
    choice.add_argument('--model',choices=MODELS)
    step=parser.add_mutually_exclusive_group(required=True)
    step.add_argument('--preflight-only',action='store_true')
    step.add_argument('--check-only',action='store_true')
    step.add_argument('--match',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--timeout-sec',type=int,default=7200)
    parser.add_argument('--overwrite',action='store_true')
    args=parser.parse_args()
    targets=MODELS if args.all else (args.model,)
    for model in targets:
        if args.preflight_only:
            result=preflight(model)
            print(f'{model}: independent canonical preflight PASS; source SHA256={result["tree"]}',flush=True)
        else:
            run_one(model,mode='check' if args.check_only else 'match',
                    outdir=args.output_dir.resolve(),wolframscript=args.wolframscript,
                    timeout_sec=args.timeout_sec,overwrite=args.overwrite)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
