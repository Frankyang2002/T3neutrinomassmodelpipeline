#!/usr/bin/env python3
"""Live Matchete CG inventory and cross-model consistency, with physical bridge gate.

Usage (repository root, in existing venv):
    python Lagrangian/RunT3UVNativeBridge.py --self-test
    python Lagrangian/RunT3UVNativeBridge.py --all

This is NOT a Matchete/RGBeta physical coupling certification. RGBeta native
triplet matrices and field-component identification remain missing in the
sourced Stage27/28 evidence; the final report must disclose that limitation.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import re


def exact_number(raw: str) -> sp.Expr:
    """Read only exact finite Wolfram-style numeric constants from CG exports."""
    if not isinstance(raw, str) or len(raw)>512:
        raise ValueError('Malformed native CG scalar')
    s=raw.strip().replace('Sqrt[','sqrt(').replace(']',')').replace('^','**')
    if not re.fullmatch(r'[\d\s+\-*/().sqrtI]+',s) or re.search(r'\d+\.\d+',s):
        raise ValueError('Non-exact or unsupported native scalar: '+repr(raw))
    try: v=sp.sympify(s,locals={'sqrt':sp.sqrt,'I':sp.I},evaluate=True)
    except (sp.SympifyError,ValueError,TypeError,SyntaxError) as exc:
        raise ValueError('Could not parse native CG scalar: '+repr(raw)) from exc
    if v.has(sp.Float) or v.free_symbols or v.is_finite is not True:
        raise ValueError('Non-exact or non-finite native CG scalar')
    return v


def exact_zero(v: sp.Expr) -> bool:
    return sp.simplify(v)==0

MODELS = {
    'A': {-4: (1,3,2), -2: (1,3,2), 0: (1,3,2)},
    'B': {-3: (2,2,1), -1: (2,2,1), 1: (2,2,1)},
    'C': {-3: (2,2,3), -1: (2,2,3), 1: (2,2,3)},
    'D': {-2: (3,1,2), 0: (3,1,2), 2: (3,1,2)},
    'E': {-4: (3,3,2), -2: (3,3,2), 0: (3,3,2), 2: (3,3,2)},
}
QUARTICS = {'A': {-4:0,-2:5,0:5}, 'B':{-3:1,-1:7,1:1},
            'C':{-3:1,-1:7,1:1}, 'D':{-2:5,0:5,2:0},
            'E':{-4:0,-2:6,0:6,2:0}}
LABELS = [f'{f}_{"m" if a < 0 else "p"}{abs(a)}' for f in MODELS for a in MODELS[f]]
SCRIPT = ROOT/'Lagrangian/RunT3UVNativeBridge.wl'
OUT = ROOT/'Reports/output/workstream_a/native_cg_bridge'
EVIDENCE = ROOT/'Reports/output/workstream_a/exact_uv_equivalence/exact_uv_equivalence.json'
FINAL = ROOT/'Reports/output/workstream_a/final_uv_evidence/final_uv_evidence.json'
CLOSURE = ROOT/'Reports/output/workstream_a/uv_closure_stage/uv_closure_stage.json'
ST23 = ROOT/'output/general_t3_stage22/stage23_quartic_prefactor_candidates.json'
ST27 = ROOT/'output/general_t3_stage22/stage27_triplet_generator_witness.json'
ST28 = ROOT/'output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json'
RGBETA_SOURCE = ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
CG_KEYS = ('Yukawa1', 'Yukawa2', 'TopologyHHSS')
STATUS = 'LiveMatcheteBaselineCGsPassedRGBetaPhysicalBridgeUncertified'


def require(ok: bool, msg: str) -> None:
    if not ok:
        raise ValueError(msg)


def read(path: Path) -> dict[str, Any]:
    require(path.is_file(), 'Missing prerequisite: '+str(path))
    with path.open(encoding='utf-8') as handle:
        obj=json.load(handle)
    require(isinstance(obj,dict), f'Malformed JSON: {path}')
    return obj


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def model_specs(label: str) -> tuple[str,int,tuple[int,int,int]]:
    require(label in LABELS,'Unknown model '+label)
    family, tag=label.split('_')
    alpha=(-1 if tag.startswith('m') else 1)*int(tag[1:])
    return family,alpha,MODELS[family][alpha]


def expected_shape(d1: int, d2: int, df: int, key: str) -> list[int]:
    if key=='Yukawa1': return [x for x in (2,df,d1) if x>1]
    if key=='Yukawa2': return [x for x in (2,df,d2) if x>1]
    if key=='TopologyHHSS': return [x for x in (2,2,d1,d2) if x>1]
    raise ValueError('Unknown CG type')


def flatten_cg(record: dict[str,Any], dims: list[int], label: str) -> tuple[sp.Expr,...]:
    require(record.get('dimensions')==dims, label+': incorrect component dimensions')
    n=sp.prod(dims)
    raw=record.get('components_inputform')
    require(isinstance(raw,list) and len(raw)==n and record.get('component_count')==n,
            label+': incorrect component count')
    require(record.get('tensor_accessor')=='Matchete`GetCGTensor'
            and record.get('source')=='live registered Matchete CG'
            and record.get('registered_in_assembled_interaction') is True,
            label+': no proof of registered native tensor')
    vals=tuple(exact_number(str(x)) for x in raw)
    require(any(not exact_zero(x) for x in vals),label+': zero CG tensor')
    return vals


def norm_sq(xs: tuple[sp.Expr,...]) -> sp.Expr:
    val=sp.simplify(sum(sp.conjugate(x)*x for x in xs))
    require(val.is_positive is True, 'Non-positive CG Hermitian norm')
    return val


def projectively_equal(xs:tuple[sp.Expr,...],ys:tuple[sp.Expr,...])->bool:
    if len(xs)!=len(ys):return False
    k=next((i for i,v in enumerate(xs) if not exact_zero(v)),None)
    if k is None or exact_zero(ys[k]):return False
    return all(exact_zero(sp.simplify(xs[i]*ys[k]-ys[i]*xs[k])) for i in range(len(xs)))


def hh_symmetric(flat: tuple[sp.Expr,...], dims: list[int]) -> bool:
    # Native first two slots are the identical Higgs doublets.
    assert dims[:2]==[2,2]
    strides=[]
    for i in range(len(dims)):
        p=1
        for d in dims[i+1:]:p*=d
        strides.append(p)
    for idx in product(*(range(d) for d in dims)):
        flipped=(idx[1],idx[0],*idx[2:])
        a=sum(i*s for i,s in zip(idx,strides))
        b=sum(i*s for i,s in zip(flipped,strides))
        if not exact_zero(flat[a]-flat[b]): return False
    return True


def check_model(label:str, report:dict[str,Any]) -> tuple[dict[str,Any],dict[str,tuple[sp.Expr,...]]]:
    f,a,(d1,d2,df)=model_specs(label)
    require(report.get('format')=='t3-workstream-a-live-baseline-cg-v1'
            and report.get('status')=='LiveRegisteredBaselineCGsCaptured'
            and report.get('model')==label
            and report.get('dimensions')==[d1,d2,df]
            and report.get('alpha')==a
            and report.get('full_matchete_check_passed') is True
            and report.get('indexed_flavour_mf') is True
            and report.get('quartics_included')==QUARTICS[f][a]
            and report.get('extra_yukawas_included')==int(a in (-2,0) and f in ('A','D','E'))
            and report.get('rgbeta_component_basis_aligned') is False
            and report.get('signed_physical_coupling_map_certified') is False
            and report.get('matching_performed') is False,
            f'{label}: malformed/unverified native CG report')
    cg=report.get('native_tensors')
    require(isinstance(cg,dict) and set(cg)==set(CG_KEYS),label+': expected 3 native CGs')
    tensors={}
    gram={}
    for key in CG_KEYS:
        expected={'Yukawa1':'T3Y1CG','Yukawa2':'T3Y2CG','TopologyHHSS':'T3MixCG'}[key]
        require(cg[key].get('symbol')==expected, f'{label}: unexpected native symbol for {key}')
        dims=expected_shape(d1,d2,df,key)
        tensors[key]=flatten_cg(cg[key],dims,f'{label}/{key}')
        gram[key]=str(norm_sq(tensors[key]))
    require(hh_symmetric(tensors['TopologyHHSS'],expected_shape(d1,d2,df,'TopologyHHSS')),
            label+': HH exchange symmetry failed for live topology tensor')
    return ({'model':label,'dimensions':[d1,d2,df],
             'registered_baseline_cg_count':3,'component_count':sum(map(len,tensors.values())),
             'native_gram_norms':gram,'higgs_exchange_symmetric':True,
             'matchete_builder_passed':True},tensors)


def check_source_evidence(root:Path=ROOT)->dict[str,Any]:
    prev=read(root/'Reports/output/workstream_a/exact_uv_equivalence/exact_uv_equivalence.json')
    require(prev.get('status')=='ExactScalarIdentityAndConditionalNormalizationPassedNotFullPhysicalCertificate',
            'Prior exact UV equivalence report not passing')
    fin=read(root/'Reports/output/workstream_a/final_uv_evidence/final_uv_evidence.json')
    require(fin.get('status')=='ConditionalUVOperatorEvidencePassedNotFullCompleteness',
            'Prior final UV operator evidence not passing')
    closure=read(root/'Reports/output/workstream_a/uv_closure_stage/uv_closure_stage.json')
    require(closure.get('status')=='StructuralAndNumericalStagePassed',
            'Prior UV closure not passing')
    stage23=read(root/'output/general_t3_stage22/stage23_quartic_prefactor_candidates.json')
    require(stage23.get('status')=='SourceQuarticNormalisationVerifiedConditionalSignedMapsNotPhysicalCertification'
            and stage23.get('totals',{}).get('slots')==50
            and stage23.get('totals',{}).get('blocks')==46
            and stage23.get('signed_physical_coupling_map_certified') is False,
            'Stage23 is not the expected conditional normalization baseline')
    stage27=read(root/'output/general_t3_stage22/stage27_triplet_generator_witness.json')
    stage28=read(root/'output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json')
    require(stage27.get('source_native_matchete_tensor_count')==6
            and stage27.get('global_matchete_rgbeta_cartesian_basis_certified') is False
            and stage27.get('physical_coupling_mapping_certified') is False,
            'Stage27 common-basis status changed: reassess certification')
    require(stage28.get('models_checked')==6
            and stage28.get('native_numerical_triplet_generators_available') is False
            and stage28.get('physical_cross_package_component_basis_certified') is False
            and stage28.get('physical_coupling_mapping_certified') is False,
            'Stage28 triplet native evidence changed: reassess certification')
    source=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    code=source.read_text(encoding='utf-8')
    snippets=('T3RGBetaYukawaInvariant[2, 1]', 'T3RGBetaYukawaInvariant[1, 2]',
              'T3RGBetaYukawaInvariant[2, 3]', 'T3RGBetaYukawaInvariant[3, 2]',
              'T3RGBetaMixInvariant[1, 3]', 'T3RGBetaMixInvariant[3, 1]',
              'T3RGBetaMixInvariant[2, 2]', 'T3RGBetaMixInvariant[3, 3]',
              '{Bar @ S1, l, F}', '{S2, l, secondHeavy}',
              '{H, H, S1, Bar @ S2}')
    require(all(code.count(s)>=1 for s in snippets),
            'RGBeta native baseline Yukawa/topology definitions changed')
    files=[root/p for p in ('Reports/output/workstream_a/exact_uv_equivalence/exact_uv_equivalence.json',
        'Reports/output/workstream_a/final_uv_evidence/final_uv_evidence.json',
        'Reports/output/workstream_a/uv_closure_stage/uv_closure_stage.json',
        'output/general_t3_stage22/stage23_quartic_prefactor_candidates.json',
        'output/general_t3_stage22/stage27_triplet_generator_witness.json',
        'output/general_t3_stage22/stage28_live_rgbeta_triplet_structure.json',
        'RGE/running/rgbeta/T3RGBetaModel.wl')]
    return {'source_sha256':{str(p.relative_to(root)):digest(p) for p in files},
            'prior_exact_scalar_identity_passed':True,
            'native_extra_quartic_slots_conditional_map':50,
            'conditional_quartic_conversion_blocks':46,
            'rgbeta_baseline_operator_definitions_present':11,
            'rgbeta_triplet_numeric_generators_available':False,
            'matchete_rgbeta_common_triplet_basis_certified':False,
            'signed_physical_quartic_dictionary_certified':False}


def make_report(raw:dict[str,dict[str,Any]], *, source:dict[str,Any]|None=None)->dict[str,Any]:
    require(set(raw)==set(LABELS),'Need exactly 16 model-specific Matchete records')
    byfam=defaultdict(list)
    byshape=defaultdict(list)
    out=[]
    allcg={}
    for label in LABELS:
        row,tensors=check_model(label,raw[label])
        out.append(row)
        allcg[label]=tensors
        f,a,(d1,d2,df)=model_specs(label)
        byfam[f].append(label)
        byshape[(d1,d2)].append(label)
    comparisons=[]
    for f,labels in byfam.items():
        ref=labels[0]
        for label in labels[1:]:
            for name in CG_KEYS:
                require(projectively_equal(allcg[ref][name],allcg[label][name]),
                        f'{ref} vs {label} {name}: native CG changed with alpha')
                comparisons.append({'reference':ref,'model':label,'tensor':name,
                                    'projective_CG_identity':True})
    # B/C have the same scalar dimensions but different fermion dimensions;
    # compare the topology tensor only, not Yukawa tensors.
    for i in range(3):
        b=byfam['B'][i];c=byfam['C'][i]
        require(projectively_equal(allcg[b]['TopologyHHSS'],allcg[c]['TopologyHHSS']),
                f'{b}/{c}: scalar topology native CG depends on fermion rep')
        comparisons.append({'reference':b,'model':c,'tensor':'TopologyHHSS',
                            'projective_CG_identity':True})
    return {'format':'t3-workstream-a-native-cg-bridge-v1',
            'status':STATUS, 'models_checked':16,'native_registered_CG_count':48,
            'higgs_exchange_symmetric_models':16,
            'cross_model_projective_CG_comparisons_passed':len(comparisons),
            'model_rows':out,'cross_model_comparisons':comparisons,
            'source_evidence':source or {'synthetic_test_only':True},
            'all_matchete_registered_CG_sources_verified':True,
            'matchete_rgbeta_component_basis_certified':False,
            'physical_signed_coupling_dictionary_certified':False,
            'rgbeta_native_triplet_generators_available':False,
            'explicit_z2_registration_certified':False,
            'physical_matching_performed':False}


def write_report(out:Path,report:dict[str,Any])->None:
    out.mkdir(parents=True,exist_ok=True)
    (out/'native_cg_bridge.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    lines=['# Workstream A — Live baseline Matchete CG bridge', '',
        f'**Status:** `{report["status"]}`', '',
        f'- Live registered baseline CGs: **{report["native_registered_CG_count"]}/48**.',
        f'- Symmetric Higgs-Higgs topology CGs: **{report["higgs_exchange_symmetric_models"]}/16**.',
        f'- Exact projective cross-model comparisons: **{report["cross_model_projective_CG_comparisons_passed"]}**.',
        '- Independent source: installed Matchete `GetCGTensor` (not generated tensor candidates).',
        '- Native RGBeta baseline vertex definitions: checked against `T3RGBetaModel.wl`.', '',
        '| Model | CGs | Components | HH symmetric |',
        '|---|---:|---:|---|']
    for row in report['model_rows']:
        lines.append(f'| {row["model"]} | {row["registered_baseline_cg_count"]} | {row["component_count"]} | PASS |')
    lines+=['','## What is not certified','',
        'Matchete registered CGs now have a direct live extraction and exact internal',
        'consistency checks. This is **not** an independent RGBeta component-level',
        'comparison or a physical signed coupling dictionary. The prior RGBeta native',
        '`S2` generator evidence remains symbolic (no numerical triplet matrices),',
        'and the required conjugate-representation/Bar field intertwiners have not',
        'been identified in a common basis. Therefore **cross-package normalization',
        'and full physical UV equivalence remain unresolved**.', '',
        'The 46 Stage23 / 50-slot quartic map is conditional. This stage neither',
        'calls `Match[]` nor registers an exact discrete `Z2` symmetry in Matchete.',
        'It leaves Workstreams B, C and D and existing output files unchanged.', '']
    (out/'native_cg_bridge.md').write_text('\n'.join(lines),encoding='utf-8')


def self_test()->None:
    # Small exact tensors with genuine HH symmetry. Projective/scaling tests.
    from sympy import sqrt
    v=(sp.Integer(1),sp.Integer(0),sqrt(2))
    require(projectively_equal(v,tuple(3*x for x in v)), 'Projective CG test failed')
    require(not projectively_equal(v,(sp.Integer(1),1,sqrt(2))), 'Projective mismatch was accepted')
    require(hh_symmetric((sp.Integer(1),0,0,1),[2,2]),'HH identity tensor should be symmetric')
    require(not hh_symmetric((0,1,0,0),[2,2]),'HH antisymmetric test was accepted')
    require(norm_sq(v)==3, 'CG norm incorrect')
    for label in LABELS:
        f,a,(d1,d2,df)=model_specs(label)
        for key in CG_KEYS:
            dims=expected_shape(d1,d2,df,key)
            require(len(dims)>=2,'Expected at least 2 live CG indices')
    print('PASS: 16-model shapes; exact projective CG comparison and HH symmetry tests')
    print('NOTE: --all executes the actual installed Matchete and inspects RGBeta source evidence')


def run(args:argparse.Namespace)->int:
    if args.self_test:
        self_test()
        return 0
    if not args.all and not args.from_existing:
        raise ValueError('Use --self-test, --all or --from-existing')
    source=check_source_evidence()
    out=Path(args.output_dir).resolve()
    raw={}
    for label in LABELS:
        family,alpha,_=model_specs(label)
        tag=f'{"m" if alpha<0 else "p"}{abs(alpha)}'
        js=out/'native'/f'{label}.json';log=out/'native'/f'{label}.log'
        if args.from_existing:
            require(js.is_file(),f'{label}: missing native JSON report')
        else:
            require(args.overwrite or (not js.exists() and not log.exists()),
                    f'{label}: native output exists. Use --from-existing or --overwrite')
            js.parent.mkdir(parents=True,exist_ok=True)
            env=os.environ.copy()
            env.update({'T3_EXPERIMENTAL_FLAVOR_MASSES':'1',
                        'T3_EXPERIMENTAL_GENERAL_FERMIONS':'1',
                        'T3_EXPERIMENTAL_GENERAL_QUARTICS':'1',
                        'T3_EXPERIMENTAL_GENERAL_QUADRATICS':'0'})
            cmd=[args.wolframscript,'-file',str(SCRIPT),str(js),family,tag]
            print(f'Matchete registered CG {label}',flush=True)
            try:
                proc=subprocess.run(cmd,env=env,cwd=ROOT,text=True,stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,timeout=args.timeout_seconds,check=False)
                log.write_text(proc.stdout,encoding='utf-8')
                require(proc.returncode==0,f'{label}: Wolfram exit {proc.returncode}; inspect {log}')
            except subprocess.TimeoutExpired as exc:
                log.write_text(f'TIMEOUT after {args.timeout_seconds} seconds\n'+str(exc),encoding='utf-8')
                raise ValueError(f'{label}: Wolfram timed out; inspect {log}') from exc
        raw[label]=read(js)
        row,_=check_model(label,raw[label])
        print(f'  {label}: 3/3 registered native CGs, {row["component_count"]} exact components, HH symmetry PASS',flush=True)
    final=make_report(raw,source=source)
    write_report(out,final)
    print(f'RESULT: {final["status"]} | 48/48 live CGs | '
          f'{final["cross_model_projective_CG_comparisons_passed"]} projective comparisons')
    print('JSON:',out/'native_cg_bridge.json')
    print('REPORT:',out/'native_cg_bridge.md')
    return 0


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    act=parser.add_mutually_exclusive_group(required=True)
    act.add_argument('--self-test',action='store_true')
    act.add_argument('--all',action='store_true')
    act.add_argument('--from-existing',action='store_true')
    parser.add_argument('--overwrite',action='store_true')
    parser.add_argument('--output-dir',default=str(OUT))
    parser.add_argument('--wolframscript',default='wolframscript')
    parser.add_argument('--timeout-seconds',type=int,default=900)
    args=parser.parse_args()
    try:
        return run(args)
    except (ValueError,FileNotFoundError,KeyError,TypeError,subprocess.SubprocessError) as exc:
        print('FAILED:',exc,file=sys.stderr)
        return 1


if __name__=='__main__':
    sys.exit(main())
