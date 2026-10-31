#!/usr/bin/env python3
"""Prepare read-only, fresh-kernel Matchete CG component exports for B(-1), C(-1).

This script does not alter the supplied Wolfram source or call Match. It writes
an instrumented copy to an audit directory, adding a diagnostic exit BEFORE
matching. The copy uses the *original* source directory to load dependencies.

Its result is native CG tensor data, NOT a gamma_CC-to-Weyl convention map,
not a beta function, not a counterterm or a renormalized physical C5.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

SOURCE = Path('Lagrangian/RunT3BCYOnlyCanonicalReconstruction.experimental.wl')
ROOT = Path('output/general_t3_mass_eigen_matchete')
KNOWN = '57afcd888598f86cccc02c069916fee694f0aed3'
ANCHOR = ' T3MECheckpoint["scalar-threshold-match-loader-started"];'
DIR_LINE = 'scriptDir = DirectoryName @ ExpandFileName[$InputFileName];'
SYMBOL = re.compile(r'^[A-Za-z][A-Za-z0-9]*$')
CG_RE = re.compile(r'CG\[(?:Bar\[)?([A-Za-z][A-Za-z0-9]*)')


def git_blob(blob: bytes) -> str:
    return hashlib.sha1(b'blob ' + str(len(blob)).encode() + b'\x00' + blob).hexdigest()


def paths_for(stem: str) -> list[Path]:
    return [
        ROOT / 'eft1_fermion_first_tree' / f'{stem}_final_eft1_fermion_first_tree_llss_tree.inputform.txt',
        ROOT / 'eft1_fermion_first_tree' / f'{stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt',
    ]


def extract_native_names(stem: str) -> list[str]:
    ll, eft = [p.read_text(encoding='utf-8') for p in paths_for(stem)]
    # No assumed CG labels: take every label in each actual source. Only
    # registered tensors are exported by Wolfram; extra SM tensors excluded.
    relevant = sorted(n for n in set(CG_RE.findall(ll + '\n' + eft)) if n.startswith('T3'))
    if not relevant or not all(SYMBOL.fullmatch(n) for n in relevant):
        raise ValueError(f'{stem}: missing or invalid source CG labels')
    return relevant


def wolfram_probe(names: list[str]) -> str:
    # All names come from source tokens with a strict symbol-only grammar.
    literals = '{' + ','.join('"' + n + '"' for n in names) + '}'
    return f'''\n (* BEGIN T3BC READ-ONLY CG COMPONENT EXPORT; NO Match CALLED *)
 Module[{{names={literals}, data, n, symbol, tensor, properties, result}},
   data=Table[
     symbol=ToExpression[n,InputForm];
     tensor=Quiet@Check[
       Matchete`CGManipulations`PackagePrivate`$CGtensors[symbol],$Failed];
     properties=Quiet@Check[
       Matchete`PackageScope`$CGproperties[symbol,Indices],$Failed];
     If[tensor===$Failed || properties===$Failed ||
        !FreeQ[tensor,Matchete`CGManipulations`PackagePrivate`$CGtensors],
       T3MEFail["IndependentCGRegistryUnavailable",n]];
     <|"name"->n,
       "tensor_inputform"->ToString[InputForm[Normal[tensor]]],
       "indices_inputform"->ToString[InputForm[properties]]|>,
     {{n,names}}];
   result=<|
     "format"->"t3bc-native-y-cg-component-export-v1",
     "model"->modelName,
     "source"->"original-canonical-Y-Matchete-kernel-before-Match",
     "canonical_source_checked"->True,
     "CGRegistry"->data,
     "native_Y_matching_performed_in_this_invocation"->False,
     "component_tensor_contract_performed"->False,
     "GammaCC_Weyl_normalization_certified"->False,
     "physical_C5"->Null|>;
   Export[outPath,result,"RawJSON"];
   Print["T3BC CG registry exported: ",Length[data]," -> ",outPath];
   Exit[0];
 ];
 (* END T3BC READ-ONLY CG COMPONENT EXPORT *)\n'''


def prepare(stem: str, outdir: Path, source: Path, allow_changed: bool) -> dict:
    raw=source.read_bytes()
    sha=git_blob(raw)
    if sha != KNOWN and not allow_changed:
        raise ValueError(f'{source} is a changed version (Git blob {sha}); '
                         'review before using --allow-changed-source')
    content=raw.decode('utf-8')
    if content.count(ANCHOR)!=1 or content.count(DIR_LINE)!=1:
        raise ValueError('Required injection anchors not unique in current source; fail closed')
    names=extract_native_names(stem)
    # The script needs the unchanged original Lagrangian/model/* sources.
    # Mathematica's quoted string handles forward-slash absolute paths.
    source_abs=source.resolve().as_posix().replace('\\','\\\\').replace('"','\\"')
    content=content.replace(DIR_LINE, 'scriptDir = DirectoryName @ ExpandFileName["'+source_abs+'"];',1)
    content=content.replace(ANCHOR,wolfram_probe(names)+ANCHOR,1)
    outdir.mkdir(parents=True,exist_ok=True)
    fork=outdir / f'{stem}_ExportNativeCG_NO_MATCH.experimental.wl'
    fork.write_text(content,encoding='utf-8')
    inp=ROOT / f'{stem}_final_parameters.json'
    tree=ROOT / 'eft1_fermion_first_tree' / f'{stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt'
    for file in (inp,tree):
        if not file.is_file():raise FileNotFoundError(file)
    result=outdir / f'{stem}_native_CG.json'
    return {
        'model':{'B_m1':'B_-1','C_m1':'C_-1'}[stem],
        'source_git_blob':sha,
        'input_cg_labels':names,
        'fork_path':str(fork),
        'output_path':str(result),
        'command_powershell':f'$env:T3_BC_Y_INDEPENDENT_CANONICAL="1"; '
          f'Remove-Item Env:T3_BC_Y_INDEPENDENT_CHECK_ONLY -ErrorAction SilentlyContinue; '
          f'wolframscript -file "{fork}" "{inp}" "{tree}" "{result}"',
        'provenance_not_a_counterterm':True,
    }


def self_test() -> None:
    sample=DIR_LINE+'\n'+ANCHOR+'\n'
    assert sample.count(ANCHOR)==1
    assert 'native_Y_matching_performed_in_this_invocation' in wolfram_probe(['GoodCG'])
    assert 'CGRegistry' in wolfram_probe(['A','B'])
    try:
        if not SYMBOL.fullmatch('BadName[]'):raise ValueError('Rejected non-symbol')
    except ValueError:pass
    else:raise AssertionError('Name guard not working')
    print('PASS: strict source anchor, symbol name guard, early exit without Match')


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--prepare',action='store_true')
    ap.add_argument('--source',type=Path,default=SOURCE)
    ap.add_argument('--output-dir',type=Path,
        default=Path('output/audits/t3_bc_y_native_cg_components_v1'))
    ap.add_argument('--allow-changed-source',action='store_true')
    opts=ap.parse_args()
    try:
        if opts.self_test:self_test()
        if opts.prepare:
            manifest={'format':'t3bc-native-y-cg-probe-manifest-v1','models':{}}
            for stem in ('B_m1','C_m1'):
                result=prepare(stem,opts.output_dir,opts.source,opts.allow_changed_source)
                manifest['models'][result['model']]=result
                print(result['model']+': fork created; CG labels='+str(len(result['input_cg_labels'])))
                print('  '+result['command_powershell'])
            (opts.output_dir/'prepare_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        if not (opts.self_test or opts.prepare):ap.error('Use --self-test or --prepare')
        return 0
    except (OSError,ValueError,UnicodeError) as e:
        print('FAIL (no CG export or physical matching certificate): '+str(e),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
