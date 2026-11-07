"""Batch 8: component-certified, fail-closed native X/Y scalar matching.

This launcher uses the unchanged Batch 6 native mixed-mass builder and the
pre-existing Batch 6 component proof. It does NOT integrate out F again.
Only a reconstructed EFT1 verified against that proof, and independently
passing all five Matchete checks, can reach the scalar-only Match[] call.

The saved canonical EOMEFT and the successful mixed_hard_06 proof are read-only.
A successful hard candidate does NOT establish pole cancellation, the complete
three-Majorana-flavour matching, or the physical SMEFT C5 coefficient.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from typing import Any

from scripts import RunT3MixedScalarHardBatch6 as base
from scripts import RunT3MixedComponentAuditBatch6 as audit

MODELS = ('B_-1', 'C_-1')
PROOF_FORMAT = audit.FORMAT
PROOF_STATUS = 'NativeCGComponentConditionsVerifiedNoMatch'


class GateError(ValueError):
    pass


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _run_name(value: str) -> str:
    if not re.fullmatch(r'[A-Za-z0-9_-]+', value):
        raise GateError('Invalid run ID')
    return value


def verify_certificate(model: str, evidence_root: Path, report: dict[str, Any],
                       settings: dict[str, Any], current: dict[str, Any],
                       *, recheck: bool = True) -> dict[str, Any]:
    """Certify the exact saved v6 witness; do not trust a JSON PASS flag alone."""
    if model not in MODELS:
        raise GateError('Only mixed B_-1/C_-1 models are supported')
    if (report.get('format') != PROOF_FORMAT or
        report.get('status') != PROOF_STATUS or
        report.get('source_hard_run') != evidence_root.name or
        report.get('scalar_Match_performed') is not False or
        report.get('fermion_Match_recomputed') is not False or
        report.get('matchete_Hermiticity_check_passed') is not False or
        report.get('physical_C5_certified') is not False):
        raise GateError(model+': incompatible or incomplete component certificate')
    if settings.get('source_fingerprints', {}).get(model) != current.get('source_fingerprints'):
        raise GateError(model+': native first-threshold source inputs changed since component audit')
    rows = [r for r in report.get('models', []) if isinstance(r,dict) and r.get('model') == model]
    if len(rows) != 1:
        raise GateError(model+': exactly one original component-audit row required')
    row = rows[0]
    inputs = row.get('inputs', {})
    native = row.get('native_array_evidence', {})
    if (row.get('status') != PROOF_STATUS or
        row.get('component_equivalence_under_native_CG') is not True or
        row.get('native_arrays_independently_checked') is not True or
        row.get('arrays_from_current_native_run') is not True or
        row.get('scalar_Match_performed') is not False or
        row.get('matchete_Hermiticity_check_passed') is not False or
        row.get('physical_C5_certified') is not False or
        inputs.get('terms') != 137 or inputs.get('sectors') != 25 or
        native.get('native_entries_independently_verified') != 32 or
        row.get('component_conditions', {}).get('exact_component_polynomial_conditions_solved') is not True):
        raise GateError(model+': component identity/native tensor conditions not established')
    stem = audit.MODELS[model]
    folder = evidence_root / stem
    prefix = stem+'_full_eft1_residual'
    witness = {
        'saved_eomeft': folder/(prefix+'_saved_eomeft_full.inputform.txt'),
        'reconstructed_candidate': folder/(prefix+'_reconstructed_candidate_full.inputform.txt'),
        'symbolic_residual': folder/(prefix+'_full_symbolic_residual.inputform.txt'),
    }
    sha = inputs.get('files_sha256')
    if not isinstance(sha,dict) or set(sha)!=set(witness):
        raise GateError(model+': incomplete witness digest inventory')
    for key,path in witness.items():
        if not path.is_file() or file_sha256(path)!=sha[key]:
            raise GateError(model+': witness modified or unavailable: '+key)
    arr = folder/(stem+'_nativeCG_native_cross_cg_arrays.inputform.txt')
    if not arr.is_file() or file_sha256(arr) != native.get('native_tensor_sha256'):
        raise GateError(model+': native tensor array SHA256 differs from audited run')
    # Independently reexecute the repository's exact 25-sector solver and
    # all 32 native tensor-entry comparisons, without any Wolfram execution.
    if recheck:
        verified = audit.audit_model_run(folder,model)
        if (verified.get('status')!=PROOF_STATUS or
            verified.get('component_equivalence_under_native_CG') is not True or
            verified.get('inputs', {}).get('files_sha256') != sha or
            verified.get('native_array_evidence', {}).get('native_tensor_sha256') !=
                native.get('native_tensor_sha256')):
            raise GateError(model+': independent current component recheck failed')
    return {
        'model':model,
        'source_component_run': evidence_root.name,
        'saved_witness':str(witness['saved_eomeft'].resolve()),
        'reconstructed_witness':str(witness['reconstructed_candidate'].resolve()),
        'native_array_witness':str(arr.resolve()),
        'witness_sha256':sha,
        'native_array_sha256':native['native_tensor_sha256'],
        'component_equivalence_under_native_CG':True,
        'matchete_Hermiticity_check_passed':False,
        'physical_C5_certified':False,
    }


def gated_helper() -> str:
    """Alter the *temporary* Wolfram runner, never the saved native driver.

    The previous native helper already constructs the full reconstruction and
    exports 25-sector evidence. Replace only its old literal-equality gate;
    require exact identity to the audited v6 witnesses before validating the
    new reconstruction with Matchete and proceeding to one-loop scalar Match.
    """
    src = base.NATIVE_HELPER
    # Keep all newly introduced Wolfram working symbols local to T3B6Hard.
    local_anchor = '    rebased, sourceCheck, proof, stageName},'
    if src.count(local_anchor) != 1:
        raise GateError('Native Wolfram Module local-variable anchor changed')
    src = src.replace(local_anchor,
        '    rebased, sourceCheck, proof, stageName, proofSavedText, proofNativeText,\n'
        '    proofTensorText, tensorText, proofSavedExpr, proofNativeExpr,\n'
        '    sameSaved, sameNative, sameTensor, componentMatch},', 1)
    old = '    sourceCheck=If[TrueQ[equality],T3B6Validate[reconstructed],'

    end = '  reconstructionPath=FileNameJoin[{out,"hard_scalar_reconstruction_audit.json"}];'
    if src.count(old)!=1 or src.count(end)!=1:
        raise GateError('Existing native reconstruction source changed; refusing unreviewed injection')
    i=src.index(old)
    j=src.index(end,i)
    replacement = r'''
    (* Batch 8: complete exact tensor/component certificate was recomputed
       on the Python side, with source hashes and native arrays verified.
       Tie the newly built two actions to BOTH exact certified witnesses.
       A Matchete-syntactic difference between saved and reconstructed is
       expected; do not claim symbolic equality of those two actions. *)
    T3MECheckpoint["batch8-component-certificate-replay-started"];
    proofSavedText=Quiet@Check[Import[Environment["T3_BATCH8_PROOF_SAVED"],"Text"],$Failed];
    proofNativeText=Quiet@Check[Import[Environment["T3_BATCH8_PROOF_RECONSTRUCTED"],"Text"],$Failed];
    proofTensorText=Quiet@Check[Import[Environment["T3_BATCH8_PROOF_ARRAYS"],"Text"],$Failed];
    tensorText=Quiet@Check[Import[FileNameJoin[{out,
       If[modelName==="B_-1","B_m1","C_m1"]<>
       "_nativeCG_native_cross_cg_arrays.inputform.txt"}],"Text"],$Failed];
    If[!And@@(StringQ/@{proofSavedText,proofNativeText,proofTensorText,tensorText}) ||
       StringLength[proofSavedText]<5000 || StringLength[proofNativeText]<5000,
      Return[T3B6Fail["Batch8 certified witness files missing or empty"]]];
    proofSavedExpr=Quiet@Check[ToExpression[proofSavedText,InputForm],$Failed];
    proofNativeExpr=Quiet@Check[ToExpression[proofNativeText,InputForm],$Failed];
    sameSaved=proofSavedExpr=!=$Failed &&
       (StringTrim[ToString[rebased,InputForm]]===StringTrim[proofSavedText] ||
        TrueQ[T3B6EquivalentQ[rebased,proofSavedExpr]]);
    sameNative=proofNativeExpr=!=$Failed &&
       (StringTrim[ToString[reconstructed,InputForm]]===StringTrim[proofNativeText] ||
        TrueQ[T3B6EquivalentQ[reconstructed,proofNativeExpr]]);
    sameTensor=StringTrim[tensorText]===StringTrim[proofTensorText];
    componentMatch=TrueQ[sameSaved] && TrueQ[sameNative] && TrueQ[sameTensor] &&
       Environment["T3_BATCH8_GATE_VERIFIED"]==="1";
    T3MECheckpoint["batch8-component-certificate-replay-returned"];
    If[!TrueQ[componentMatch],
      Export[FileNameJoin[{out,"hard_scalar_component_gate_native.json"}],<|
        "status"->"WitnessMismatch", "saved_witness_matches"->TrueQ[sameSaved],
        "reconstruction_witness_matches"->TrueQ[sameNative],
        "cross_CG_arrays_match"->TrueQ[sameTensor],
        "fermion_matching_recomputed"->False,"scalar_Match_performed"->False|>,"RawJSON"];
      Return[T3B6Fail["Batch8 component-certificate source/tensor witness mismatch"]]];
    T3MECheckpoint["batch8-reconstruction-independent-CheckLagrangian-started"];
    sourceCheck=T3B6Validate[reconstructed];
    T3MECheckpoint["batch8-reconstruction-independent-CheckLagrangian-returned"];
    reconAudit=<|
      "status"->If[TrueQ[sourceCheck["passed"]],
        "CertifiedComponentReconstructionAndMatcheteChecksPassed",
        "CertifiedComponentReconstructionFailedMatcheteChecks"],
      "algebraic_equality_proven"->TrueQ[equality],
      "component_equivalence_proven"->True,
      "independent_matchete_check_passed"->TrueQ[sourceCheck["passed"]],
      "saved_LLSS_expanded_term_count"->Length[savedMajorana],
      "native_F_free_term_count"->Length[Select[nativeParts,
         FreeQ[#,Verbatim[NewFermion],{0,Infinity}]&]],
      "residual_manifest_file"->If[modelName==="B_-1","B_m1","C_m1"]<>
         "_full_eft1_residual_full_residual_manifest.json",
      "residual_structural_sectors"->Length[proof["coupling_groups"]],
      "residual_additive_terms"->proof["residual_additive_terms"],
      "component_source_run"->Environment["T3_BATCH8_COMPONENT_RUN"],
      "checks"->sourceCheck["checks"]|>;
    If[TrueQ[sourceCheck["passed"]],
      candidate=reconstructed;
      chosen="VerifiedComponentEquivalentNativeReconstruction"];
  ];
'''
    new=src[:i]+replacement+src[j:]
    # Ensure even a hypothetical passing *original* EFT1 cannot silently
    # bypass the component certificate when using this gated entrypoint.
    guard='  tree=candidate;\n'
    if new.count(guard)!=1:
        raise GateError('Native matching entrypoint changed')
    new=new.replace(guard,r'''
  If[chosen=!="VerifiedComponentEquivalentNativeReconstruction",
     Return[T3B6Fail["Batch8 full component/EFT1 native validation gate not passed"]]];
  tree=candidate;
''',1)
    old_report='"algebraic_reconstruction_used"->(chosen==="VerifiedEquivalentNativeReconstruction"),'
    if new.count(old_report)!=1:
        raise GateError('Native hard-report output layout changed')
    new=new.replace(old_report,r'''"algebraic_reconstruction_used"->False,
     "component_equivalent_native_reconstruction_used"->True,
     "component_proof_source_hard_run"->Environment["T3_BATCH8_COMPONENT_RUN"],''',1)
    if (new.index('batch8-reconstruction-independent-CheckLagrangian-started') >=
        new.index('batch6-native-scalar-one-loop-match-started') or
        new.count('Match[tree,EFTOrder->5,LoopOrder->1]') != 1):
        raise GateError('Batch8 scalar matching isolation not established')
    return new


def main(argv: list[str] | None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--component-run',default='mixed_hard_06')
    p.add_argument('--run-id',default='mixed_hard_08')
    p.add_argument('--models',default='B_-1,C_-1')
    p.add_argument('--preflight-only',action='store_true')
    p.add_argument('--timeout-seconds',type=int,default=3600)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--sync-root',type=Path,default=base.SYNC)
    p.add_argument('--running-root',type=Path,default=base.RUNNING_ROOT)
    a=p.parse_args(argv)
    component_id=_run_name(a.component_run)
    output_id=_run_name(a.run_id)
    if component_id==output_id:
        raise GateError('Component evidence and new matching run IDs must differ')
    models=tuple(m.strip() for m in a.models.split(','))
    if not models or len(set(models))!=len(models) or any(m not in MODELS for m in models):
        raise GateError('Only unique B_-1,C_-1 can be selected')
    evidence_root=base.OUTPUT/component_id
    audit_file=evidence_root/'component_equivalence_audit.json'
    report=audit.load_json(audit_file)
    settings=audit.load_json(evidence_root/'settings.json')
    current:dict[str,dict[str,Any]]={}
    gates:dict[str,dict[str,Any]]={}
    errors=0
    for model in models:
        try:
            current[model]=base.preflight(model,a.sync_root,a.running_root)
            gates[model]=verify_certificate(model,evidence_root,report,settings,current[model])
            print(model+': READY: original F-free EOMEFT, independent 25-sector/32-component proof',
                  flush=True)
        except Exception as exc:
            errors+=1
            print(model+': BLOCKED: '+str(exc),flush=True)
    if errors:
        return 2
    if a.preflight_only:
        return 0
    # No post-hoc rewriting of a previously generated Wolfram script or
    # modifications to the native driver. All binding is in this one process.
    proof_by_path={str(Path(v['canonical_eom']).resolve()):gates[m]
                   for m,v in current.items()}
    original_env=base.native_subprocess_environment
    old_helper=base.NATIVE_HELPER
    base.NATIVE_HELPER=gated_helper()
    def env_with_proof(parent:dict[str,str],canonical_eom:str)->dict[str,str]:
        env=original_env(parent,canonical_eom)
        gate=proof_by_path.get(str(Path(canonical_eom).resolve()))
        if gate is None:
            raise GateError('No frozen matching proof corresponds to selected canonical EOMEFT')
        env.update({
            'T3_BATCH8_PROOF_SAVED':gate['saved_witness'],
            'T3_BATCH8_PROOF_RECONSTRUCTED':gate['reconstructed_witness'],
            'T3_BATCH8_PROOF_ARRAYS':gate['native_array_witness'],
            'T3_BATCH8_COMPONENT_RUN':component_id,
            'T3_BATCH8_GATE_VERIFIED':'1',
        })
        return env
    base.native_subprocess_environment=env_with_proof
    args=['--run-id',output_id,'--models',','.join(models),
          '--timeout-seconds',str(a.timeout_seconds),'--wolframscript',a.wolframscript,
          '--sync-root',str(a.sync_root),'--running-root',str(a.running_root)]
    try:
        rc=base.run(args)
    finally:
        base.NATIVE_HELPER=old_helper
        base.native_subprocess_environment=original_env
    # Keep the accepted external evidence and hashes next to the new run;
    # do not pretend these are a completed scalar hard-matching certificate.
    target=base.OUTPUT/output_id
    if target.exists():
        obj={'format':'t3-batch8-gated-reconstruction-source-v1',
             'source_component_run':component_id,
             'source_component_audit_sha256':file_sha256(audit_file),
             'source_settings_sha256':file_sha256(evidence_root/'settings.json'),
             'models':gates,
             'fermion_Match_recomputed':False,
             'component_equivalence_verified':True,
             'full_pole_cancellation_certified':False,
             'physical_C5_certified':False}
        (target/'component_gate.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    return rc


if __name__=='__main__':
    try:
        raise SystemExit(main())
    except Exception as exc:
        print('BATCH8 ERROR: '+str(exc),file=sys.stderr)
        raise SystemExit(2)
