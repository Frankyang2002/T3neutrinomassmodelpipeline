"""Batch 6: guarded native X/Y scalar threshold after SAVED B/C fermion EFT1.

This instruments a temporary copy of the current LOCAL native mixed-doublet
builder at the point immediately BEFORE its fermion Match. No fermion Match is
called. The saved canonical EOMEFT is then used with X and Y declared Heavy.
Do not interpret success of scalar Match as pole/renormalisation certification.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import traceback
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl'
SYNC = ROOT / 'output/general_observable/mixed_f_synchronised_v1'
RUNNING_ROOT = ROOT / 'output/general_observable/mixed_adapter_v3'
OUTPUT = ROOT / 'output/general_sequential_pilot/mixed_scalar_hard_batch6'
MODELS = ('B_-1', 'C_-1')
STAGE = 'MixedScalarHardCandidateExported'
HOOK = 'T3MECheckpoint["fermion-first-production-Match-loader-started"];'
FUNC = 'T3MEPerform[] := Module['
BOTTOM = 'run = CheckAbort['

# Inserted before T3MEPerform; operates inside the native, source-verified
# eigen-doublet CG/field registry. Uses the source's existing U, quartics and
# Matchete tensor conventions. This is intentionally not a newly guessed CG map.
NATIVE_HELPER = r'''
(* Batch 6 deliberately resumes at the saved F-free EOMEFT. Every operation
   below is scalar-only. The source B/C mass-eigenfield builder runs first. *)
T3B6Fail[why_String] := Module[{r},
  r=<|"status"->"Blocked", "model"->modelName,
    "reason"->why, "fermion_matching_recomputed"->False,
    "full_pole_cancellation_certified"->False,
    "physical_C5_certified"->False|>;
  Quiet@Check[Export[FileNameJoin[{DirectoryName[outPath],
    "hard_scalar_report.json"}],r,"RawJSON"],$Failed];
  Print["BATCH6 FAILURE: ",why]; r
];

(* A saved Matchete expression can be mathematically equivalent but fail a
   syntactic Hermiticity check. Never accept a rewrite unless it is exactly
   equal to the saved expression (modulo Matchete's dummy relabelling) AND
   all five detailed native checks pass. No (L+HC[L])/2 substitution. *)
T3B6EquivalentQ[old_,candidate_] := Module[{delta, relabeled, proof},
  delta=Expand[old-candidate];
  If[TrueQ[delta===0],Return[True]];
  relabeled=CheckAbort[TimeConstrained[RelabelIndices[delta],30,$TimedOut],$Aborted];
  If[MemberQ[{$Aborted,$Failed,$TimedOut},relabeled],Return[False]];
  If[TrueQ[relabeled===0],Return[True]];
  proof=CheckAbort[TimeConstrained[Simplify[relabeled==0],30,False],$Aborted];
  TrueQ[proof]
];

T3B6Validate[expr_] := Module[{ord,detail,keys},
  keys={"Hermiticity","CanonicallyNormalized","HeavyMassBasis",
        "GaugeInvariant","AllObjectsDefined"};
  ord=CheckAbort[TimeConstrained[UsingFrontEnd[CheckLagrangian[expr]],900,$TimedOut],$Aborted];
  detail=CheckAbort[TimeConstrained[UsingFrontEnd[
       CheckLagrangian[expr,DetailedOutput->True]],900,$TimedOut],$Aborted];
  <|"ordinary"->ToString[ord,InputForm],
    "checks"->AssociationMap[If[AssociationQ[detail],
       ToString[Lookup[detail,#,"Unavailable"],InputForm],"Unavailable"]&,keys],
    "passed"->(TrueQ[ord] && AssociationQ[detail] &&
       And@@(TrueQ[Lookup[detail,#,False]]& /@ keys))|>
];

T3B6Hard[nativeUV_] := Module[{source, saved, sourceText, nativeParts,
    nMF, mfPattern, nIndexed, nUnindexed, nOther, mode,
    lightFree, heavyFree, originalCounts, rebasedCounts,
    tree, result, raw, current, nxt, steps, hard, ext, out,
    fields, ok, first, trial, attempts={}, chosen="None", candidate,
    savedMajorana, nativeRenorm, reconstructed, equality, reconAudit,
    a, allParts, expectedLLSS, reconstructionPath, checksPath,
    rebased, sourceCheck, proof, stageName},
  out=DirectoryName[outPath];
  source=Environment["T3_BATCH6_EFT1_EOM"];
  If[!StringQ[source] || !FileExistsQ[source],
    Return[T3B6Fail["Saved canonical F-only EOMEFT not found"]]];
  sourceText=Import[source,"Text"];
  If[!StringQ[sourceText] || StringLength[sourceText]<5000 ||
     StringContainsQ[sourceText,"Field[NewFermion,"] ||
     !StringContainsQ[sourceText,"Field[NewScalarX,"] ||
     !StringContainsQ[sourceText,"Field[NewScalarY,"] ||
     !StringContainsQ[sourceText,"GammaCC"],
    Return[T3B6Fail["Saved F-free EOMEFT field/Wilson preflight failed"]]];
  saved=Quiet@Check[ToExpression[sourceText,InputForm],$Failed];
  If[saved===$Failed || AtomQ[saved],
    Return[T3B6Fail["Saved canonical EOMEFT parse failed"]]];
  mfPattern=HoldPattern[Matchete`Coupling[MF,{_},0]];
  nIndexed=Count[saved,mfPattern,Infinity];
  nUnindexed=Count[saved,Matchete`Coupling[MF,{},0],Infinity];
  nOther=Count[saved,Matchete`Coupling[MF,_,_],Infinity]-nIndexed-nUnindexed;
  If[nOther=!=0 || (nIndexed>0 && nUnindexed>0) ||
     (nIndexed===0 && nUnindexed===0),
    Return[T3B6Fail["Ambiguous native MF source type"]]];
  mode=If[nIndexed>0,"IndexedMFConvertedToSingleMassTemplate",
                      "SavedUnindexedSingleMassTemplate"];
  nMF=nIndexed+nUnindexed;
  If[nMF<2 || nMF>100,
    Return[T3B6Fail["Unexpected heavy mass factor inventory"]]];
  If[nIndexed>0,
    saved=saved /. mfPattern :> Matchete`Coupling[MF,{},0]];

  (* In the validated first threshold the scalar fields were LIGHT. Capture
     the entire free action in that registry, before changing their status.
     Preserve the integrated Majorana mass, but do NOT re-integrate F. *)
  T3MECheckpoint["batch6-F-deregister-started"];
  ok=Quiet@Check[RemoveField[NewFermion];True,$Failed];
  If[ok=!=True,Return[T3B6Fail["Cannot deregister integrated fermion"]]];
  ok=Quiet@Check[PreserveIntegratedT3MassCouplings[{"S1","S2"}],$Failed];
  If[ok=!=True,Return[T3B6Fail["Cannot retain integrated real MF coupling"]]];
  T3MECheckpoint["batch6-light-free-capture-started"];
  lightFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]//RelabelIndices,$Aborted];
  If[MemberQ[{$Aborted,$Failed},lightFree],
    Return[T3B6Fail["Light-scalar FreeLag capture failed"]]];

  (* The source is canonical F-first EOMEFT, not MatchedEFT. Its free mass
     is order 1 because X/Y were LIGHT. A mass-token-only substitution is
     insufficient to rebase the whole free action. *)
  originalCounts=(Count[saved,Matchete`Coupling[#,{},1],Infinity]& /@ {MSX,MSY});
  If[!And@@(Positive /@ originalCounts),
    Return[T3B6Fail["Saved light-scalar mass terms missing"]]];
  T3MECheckpoint["batch6-heavy-scalar-registry-started"];
  ok=Quiet@Check[RemoveField[NewScalarX];RemoveField[NewScalarY];True,$Failed];
  If[ok=!=True,Return[T3B6Fail["Cannot reclassify both scalar doublets"]]];
  ok=Quiet@Check[
    DefineField[NewScalarX,Scalar,SelfConjugate->False,
      Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSX}];
    DefineField[NewScalarY,Scalar,SelfConjugate->False,
      Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSY}];True,$Failed];
  If[ok=!=True,Return[T3B6Fail["Cannot define Heavy X/Y"]]];
  fields=GetFields[];
  If[KeyExistsQ[fields,NewFermion] ||
     !And@@(TrueQ[Lookup[Lookup[fields,#,<||>],Heavy,False]]& /@
       {NewScalarX,NewScalarY}),
    Return[T3B6Fail["Scalar-only heavy registry invalid"]]];
  heavyFree=CheckAbort[FreeLag[NewScalarX,NewScalarY]//RelabelIndices,$Aborted];
  If[MemberQ[{$Failed,$Aborted},heavyFree],
    Return[T3B6Fail["Heavy-scalar FreeLag capture failed"]]];
  T3MECheckpoint["batch6-canonical-free-rebase-started"];
  rebased=CheckAbort[RelabelIndices[Expand[saved-lightFree+heavyFree]],$Aborted];
  If[MemberQ[{$Failed,$Aborted},rebased] ||
      !FreeQ[rebased,Verbatim[NewFermion],{0,Infinity}] ||
       FreeQ[rebased,GammaCC],
    Return[T3B6Fail["Canonical scalar free-action rebase failed"]]];
  rebasedCounts=(Count[rebased,Matchete`Coupling[#,{},0],Infinity]& /@ {MSX,MSY});
  If[!And@@(Positive /@ rebasedCounts) ||
      !FreeQ[rebased,Matchete`Coupling[MSX|MSY,{},1]],
    Return[T3B6Fail["Rebased free masses are not pure Heavy EFTOrder 0"]]];
  Export[FileNameJoin[{out,"hard_scalar_input_audit.json"}],<|
      "Status"->"PassedToCheckLagrangian", "Model"->modelName,
      "CanonicalEOMEFTUsed"->True, "FermionMatchingRecomputed"->False,
      "RemainingHeavyFields"->{"NewScalarX","NewScalarY"},
      "CommonMFCount"->nMF, "OriginalIndexedMFCount"->nIndexed,
      "OriginalUnindexedMFCount"->nUnindexed,
      "MFSourceMode"->mode,
      "SourceLightMassCounts"->originalCounts,
      "RebasedHeavyMassCounts"->rebasedCounts,
      "ExactLightHeavyFreeLagRebase"->True,
      "NativeScalarMassEigenbasis"->True,
      "All17NativeQuarticsFromOriginalBuilder"->True,
      "IndependentThreeMassNativeMatching"->False,
      "OriginalFermionSourceUnmodified"->True|>,"RawJSON"];

  (* The existing repository scalar-threshold probe documents syntactic
     Hermiticity failure of the imported EOMEFT. Record detailed checks,
     try ONLY equality-preserving rewrites, and then exact reconstruction
     from the already checked UV eigenbasis with all F terms removed. *)
  T3MECheckpoint["batch6-first-detailed-heavy-check-started"];
  first=T3B6Validate[rebased];
  candidate=rebased;
  If[!TrueQ[first["passed"]],
    Do[
      T3MECheckpoint["batch6-equality-rewrite-"<>stageName<>"-started"];
      trial=CheckAbort[TimeConstrained[Switch[stageName,
         "ExpandRelabel",RelabelIndices[Expand[rebased]],
         "HcSimplify",RelabelIndices[HcSimplify[rebased]],
         "GreensSimplify",RelabelIndices[GreensSimplify[rebased]],
         _, $Failed],120,$TimedOut],$Aborted];
      If[MemberQ[{$Aborted,$Failed,$TimedOut},trial] ||
         !TrueQ[T3B6EquivalentQ[rebased,trial]],
        AppendTo[attempts,<|"method"->stageName,"status"->"RejectedNotEquivalent"|>];
        Continue[]];
      sourceCheck=T3B6Validate[trial];
      AppendTo[attempts,<|"method"->stageName,
         "status"->If[sourceCheck["passed"],"Passed","NotPassed"],
         "checks"->sourceCheck["checks"]|>];
      If[TrueQ[sourceCheck["passed"]],
        candidate=trial;chosen=stageName;Break[]],
      {stageName,{"ExpandRelabel","HcSimplify","GreensSimplify"}}]
  ];

  reconAudit=<|"status"->"NotNeeded", "algebraic_equality_proven"->False,
               "independent_matchete_check_passed"->False|>;
  If[!TrueQ[first["passed"]] && chosen==="None",
    T3MECheckpoint["batch6-native-renormalizable-reconstruction-started"];
    (* Reuse the complete native SU(2)-rotated UV Lagrangian already checked
       by the first-threshold driver. Drop ALL terms containing F. No new
       guessed potential, CG tensor or Yukawa convention enters. *)
    nativeParts=With[{expanded=Expand[nativeUV]},
       If[Head[expanded]===Plus,List@@expanded,{expanded}]];
    nativeRenorm=Total[Select[nativeParts,
       FreeQ[#,Verbatim[NewFermion],{0,Infinity}]&]];
    allParts=With[{expanded=Expand[saved]},
       If[Head[expanded]===Plus,List@@expanded,{expanded}]];
    savedMajorana=Select[allParts,!FreeQ[#,GammaCC]&];
    expectedLLSS=If[modelName==="B_-1",20,24];
    If[Length[savedMajorana]=!=expectedLLSS ||
       !FreeQ[nativeRenorm,GammaCC] ||
       !FreeQ[nativeRenorm,Verbatim[NewFermion],{0,Infinity}],
      Return[T3B6Fail["Reconstruction operator inventory invalid"]]];
    reconstructed=CheckAbort[RelabelIndices[
       Expand[nativeRenorm-lightFree+heavyFree+Total[savedMajorana]]],$Aborted];
    equality=If[reconstructed===$Aborted,False,
      T3B6EquivalentQ[rebased,reconstructed]];
    sourceCheck=If[TrueQ[equality],T3B6Validate[reconstructed],
       <|"ordinary"->"NotAttempted","checks"-><||>,"passed"->False|>];
    reconAudit=<|"status"->If[equality,
          If[sourceCheck["passed"],"EquivalentReconstructionAccepted",
              "EquivalentReconstructionFailedChecks"],
          "ReconstructionNotAlgebraicallyIdentical"],
       "algebraic_equality_proven"->TrueQ[equality],
       "independent_matchete_check_passed"->TrueQ[sourceCheck["passed"]],
       "saved_LLSS_expanded_term_count"->Length[savedMajorana],
       "native_F_free_term_count"->Length[Select[nativeParts,
           FreeQ[#,Verbatim[NewFermion],{0,Infinity}]&]],
       "checks"->sourceCheck["checks"]|>;
    If[TrueQ[sourceCheck["passed"]],
       candidate=reconstructed;chosen="VerifiedEquivalentNativeReconstruction"];
  ];
  reconstructionPath=FileNameJoin[{out,"hard_scalar_reconstruction_audit.json"}];
  Export[reconstructionPath,reconAudit,"RawJSON"];
  checksPath=FileNameJoin[{out,"hard_scalar_eft1_validation.json"}];
  Export[checksPath,<|"model"->modelName,"initial_check"->first,
        "rewrite_attempts"->attempts,"accepted_method"->chosen,
        "reconstruction_audit_file"->FileNameTake[reconstructionPath],
        "source_saved_EOMEFT_unchanged"->True,
        "recomputed_F_matching"->False,"full_EFT_validated"->
          (TrueQ[first["passed"]] || chosen=!="None")|>,"RawJSON"];
  If[!TrueQ[first["passed"]] && chosen==="None",
    Return[T3B6Fail["Full canonical X/Y EFT1 CheckLagrangian failed; see hard_scalar_eft1_validation.json and hard_scalar_reconstruction_audit.json"]]];
  tree=candidate;
  (* From here on, the original first-threshold F matching is unreachable. *)
  T3MECheckpoint["batch6-native-scalar-one-loop-match-started"];
  Get[FileNameJoin[{scriptDir,"RunMatching.wl"}]];
  raw=TimeConstrained[CheckAbort[
    UsingFrontEnd[Match[tree,EFTOrder->5,LoopOrder->1]],$Aborted],3600,$TimedOut];
  If[MemberQ[{$Failed,$Aborted,$TimedOut},raw] || AtomQ[raw] ||
       FreeQ[raw,hbar],
    Return[T3B6Fail["Native scalar-only one-loop Match unavailable"]]];
  Export[FileNameJoin[{out,"hard_raw.wxf"}],raw,"WXF"];
  steps={"GreensSimplify","EOMSimplify","EvaluateLoopFunctions",
     "ReplaceEffectiveCouplings"};
  current=raw;
  Do[
    nxt=TimeConstrained[CheckAbort[UsingFrontEnd[
      Switch[steps[[k]],
      "GreensSimplify",GreensSimplify[current],
      "EOMSimplify",EOMSimplify[current],
      "EvaluateLoopFunctions",EvaluateLoopFunctions[current],
      "ReplaceEffectiveCouplings",ReplaceEffectiveCouplings[current],
      _, $Failed]],$Aborted],900,$TimedOut];
    If[MemberQ[{$Failed,$Aborted,$TimedOut},nxt] || AtomQ[nxt],
      Return[T3B6Fail["Postprocessing failed at "<>steps[[k]]]]];
    current=nxt;
    Export[FileNameJoin[{out,"hard_step_"<>ToString[k]<>".wxf"}],current,"WXF"],
    {k,Length[steps]}];
  hard=Expand[hbar Coefficient[Expand[current],hbar,1]];
  If[hard===0 || !FreeQ[hard,Verbatim[NewScalarX]|Verbatim[NewScalarY],Infinity],
    Return[T3B6Fail["Postprocessed scalar hard sector empty or retains X/Y"]]];
  Export[FileNameJoin[{out,"hard_hbar_sector.wxf"}],hard,"WXF"];
  ext=TimeConstrained[CheckAbort[
     UsingFrontEnd[ExtractWeinbergCoefficient[hard]],$Aborted],900,$TimedOut];
  If[!AssociationQ[ext] || Lookup[ext,"Status",""]=!="Success" ||
     MissingQ[Lookup[ext,"Coefficient",Missing["Missing"]]],
    Return[T3B6Fail["Native scalar hard Weinberg C5 extraction unavailable"]]];
  Export[FileNameJoin[{out,"hard_c5_candidate.txt"}],
    ToString[ext["Coefficient"],InputForm],"Text"];
  result=<|"status"->"MixedScalarHardCandidateExported","model"->modelName,
     "native_scalar_one_loop_matched"->True,
     "weinberg_candidate_terms"->Lookup[ext,"TermCount",0],
     "common_fermion_mass_factors"->nMF,
     "mf_source_mode"->mode,
     "algebraic_reconstruction_used"->(chosen==="VerifiedEquivalentNativeReconstruction"),
     "full_EFT1_CheckLagrangian_passed"->True,
     "independent_three_mass_native_matching"->False,
     "single_heavy_eigenvalue_template"->True,
     "scalar_mass_eigenstates"->{"NewScalarX","NewScalarY"},
     "fermion_matching_recomputed"->False,
     "canonical_EOMEFT_input"->True,
     "full_pole_cancellation_certified"->False,
     "complete_fixed_order_C5"->False,
     "physical_C5_certified"->False|>;
  Export[FileNameJoin[{out,"hard_scalar_report.json"}],result,"RawJSON"];
  result
];
'''

TAIL = r'''
run=CheckAbort[UsingFrontEnd[Catch[T3MEPerform[],"T3MassEigenStop"]],$Aborted];
If[!AssociationQ[run],
  (* T3MEFail already writes the actual failure stage and diagnostic.
     Do not overwrite it with generic "Uncaught native failure". *)
  If[!FileExistsQ[outPath],
    T3MEOutput[<|"status"->"Blocked","model"->modelName,
      "reason"->"Uncaught native failure", "last_stage"->$T3MELastStage,
      "physical_C5_certified"->False|>]];
  Exit[5]];
T3MEOutput[run];
If[Lookup[run,"status",""]==="MixedScalarHardCandidateExported",Exit[0],Exit[5]];
'''


def tag(model: str) -> str:
    if model not in MODELS:
        raise ValueError('Only B_-1 and C_-1 are supported')
    return model.replace('-', 'm')


def stem(model: str) -> str:
    return tag(model) + '_final'


def paths(model: str, sync: Path = SYNC, running_root: Path = RUNNING_ROOT) -> dict[str, Path]:
    base = sync / tag(model)
    s = stem(model)
    return {
        'status': base/'status.json',
        'tree_report': base/'fermion_tree'/f'{s}_eft1_fermion_first_tree.json',
        'canonical_eom': base/'fermion_tree'/f'{s}_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt',
        'matched_llss': base/'fermion_tree'/f'{s}_eft1_fermion_first_tree_llss_tree.inputform.txt',
        'parameters': base/'matchete'/f'{s}_parameters.json',
        'matchete': base/'matchete'/f'{s}_matchete.json',
        'mass_basis': base/'mass_basis'/f'{s}_mass_basis.json',
        'running': running_root/tag(model)/'running.json',
        'driver': DRIVER,
    }


def get_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        raise ValueError(f'JSON object required: {path}')
    return obj


def hash_file(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def preflight(model: str, sync: Path = SYNC, running_root: Path = RUNNING_ROOT) -> dict:
    import numpy as np
    ps = paths(model, sync, running_root)
    for name, path in ps.items():
        if not path.is_file():
            raise ValueError(f'{model}: missing {name} at {path}')
    st, tr, pa, ma, basis, run = (get_json(ps[x]) for x in
        ('status','tree_report','parameters','matchete','mass_basis','running'))
    if st.get('status') != 'CurrentRunningMassBasisFThresholdTreeLLSSExported':
        raise ValueError('F-first synchronized native source not passed')
    if tr.get('status') != 'FermionFirstTreeLLSSExtractedNotSequentialEFT1' or \
        tr.get('integrated_heavy_fields') != ['NewFermion'] or \
        tr.get('scalar_light_fields') != ['NewScalarX','NewScalarY'] or \
        tr.get('llss_term_count') != (20 if model == 'B_-1' else 24) or \
        tr.get('transition_representation') != 'EOMEFT':
        raise ValueError('F-first EOMEFT provenance/multiplicity mismatch')
    if ma.get('all_17_quartics_present') is not True or \
        ma.get('mass_basis_ordinary_check_passed') is not True or \
        ma.get('heavy_mass_basis_detailed_check_passed') is not True or \
        ma.get('all_old_scalar_fields_removed') is not True or \
        ma.get('both_yukawa_interactions_present') is not True:
        raise ValueError('Native Matchete mass-basis/quartic source is unverified')
    if pa.get('quartic_source_count') != 17 or pa.get('yukawa_source_count') != 2 or \
        pa.get('verified_17_tensor_rotation') is not True or \
        pa.get('dimensions_alpha') != [2,2,1 if model.startswith('B') else 3,-1]:
        raise ValueError('Native scalar parameter source inventory incorrect')
    if run.get('model') != model or run.get('status') != 'NumericalRenormalisableRunningOnly':
        raise ValueError('Running-input provenance mismatch')
    if basis.get('model') != model or basis.get('basis',{}).get('tachyonic') or \
        basis.get('basis',{}).get('degenerate_basis_nonunique'):
        raise ValueError('Invalid mass-eigenfield certificate')
    def mat(z):
        return np.asarray(z['real'],float)+1j*np.asarray(z['imag'],float)
    original = mat(pa['original_mass_squared_A_gev2'])
    U = mat(pa['unitary_U'])
    eig = np.asarray(pa['heavy_mass_squared_gev2'],float)
    if original.shape!=(2,2) or U.shape!=(2,2) or eig.shape!=(2,) or \
        not np.all(np.isfinite(eig)) or np.min(eig)<=0 or \
        not np.allclose(U.conj().T@U,np.eye(2),rtol=1e-11,atol=1e-11) or \
        not np.allclose(U.conj().T@original@U,np.diag(eig),rtol=1e-10,atol=1e-7):
        raise ValueError('Independent mixed scalar eigenbasis failed')
    if not np.allclose(U, mat(basis['basis']['unitary_flavour_rotation_U']),rtol=1e-10,atol=1e-10):
        raise ValueError('Matchete and saved F-first scalar mixing matrices mismatch')
    # The canonical EOMEFT and the *expanded* matched LLSS export are NOT
    # term-by-term identical. C_-1 has 8 GammaCC terms in EOMEFT but 24 in
    # the expanded LLSS source (the earlier Batch 6 guard mistakenly required
    # 24 in both). Validate each representation independently.
    eom_report = validate_eom_sources(
        model,
        ps['canonical_eom'].read_text(encoding='utf-8-sig'),
        ps['matched_llss'].read_text(encoding='utf-8-sig'),
        tr,
    )
    source = ps['driver'].read_text(encoding='utf-8-sig')
    driver_anchors(source)  # Structure/ordering check; tolerates indentation.
    if 'T3MENativeHHMassCG' not in source or 'DefineField[NewScalarX' not in source:
        raise ValueError('Native CG rotation or scalar registration not available')
    return {
        'status':'Ready', 'model':model,
        'source_fingerprints':{key:hash_file(value) for key,value in ps.items()},
        'scalar_masses_gev':np.sqrt(eig).tolist(),
        'canonical_eom':str(ps['canonical_eom'].resolve()),
        'parameters':str(ps['parameters'].resolve()),
        'running':str(ps['running'].resolve()),
        'full_pole_cancellation_certified':False,
        'eft1_source_audit': eom_report,
    }


def validate_eom_sources(model: str, eom: str, llss: str, tree_report: dict) -> dict:
    """Validate actual saved canonical EOMEFT and expanded LLSS separately.

    IMPORTANT: synchronized B/C native sources contain MF[] with EMPTY mass
    indices, not MF[{Index[..., NFlavor]}]. This is the one-MF symbolic
    source; the later three-mass lift is a separately declared assumption.
    """
    if model not in MODELS:
        raise ValueError('Unsupported mixed model')
    expected_llss = 20 if model == 'B_-1' else 24
    expected_eom = 20 if model == 'B_-1' else 8
    if tree_report.get('llss_term_count') != expected_llss:
        raise ValueError(f'{model}: first-stage native LLSS count changed')
    for name, content in [('canonical EOMEFT', eom), ('matched LLSS', llss)]:
        if len(content) < 5000 or 'Field[NewFermion,' in content or \
           not all(f'Field[{f},' in content for f in ('NewScalarX','NewScalarY')):
            raise ValueError(f'{model}: {name} lacks F-free X/Y EFT1 fields')

    def mass_counts(content: str) -> tuple[int, int, int]:
        all_mass = len(re.findall(r'Coupling\[MF,', content))
        indexed = len(re.findall(r'Coupling\[MF,\s*\{\s*Index\[', content))
        unindexed = len(re.findall(r'Coupling\[MF,\s*\{\s*\}\s*,\s*0\s*\]', content))
        return all_mass, indexed, unindexed

    eom_gamma = eom.count('GammaCC')
    llss_gamma = llss.count('GammaCC')
    if eom_gamma != expected_eom:
        raise ValueError(f'{model}: canonical EOMEFT GammaCC count {eom_gamma} != {expected_eom}')
    if llss_gamma != expected_llss:
        raise ValueError(f'{model}: matched LLSS GammaCC count {llss_gamma} != {expected_llss}')
    all_eom, indexed_eom, unindexed_eom = mass_counts(eom)
    all_llss, indexed_llss, unindexed_llss = mass_counts(llss)
    # Neither format may contain an unrecognised third MF representation.
    # Reject mixed indexed/unindexed forms, not merely a count difference.
    for label, gamma, all_mass, indexed, unindexed in (
        ('canonical EOMEFT', eom_gamma, all_eom, indexed_eom, unindexed_eom),
        ('expanded LLSS', llss_gamma, all_llss, indexed_llss, unindexed_llss),
    ):
        if all_mass != gamma or indexed + unindexed != gamma or \
           (indexed and unindexed):
            raise ValueError(f'{model}: {label} LLSS/MF source inconsistent: '
                             f'GammaCC={gamma}, indexed MF={indexed}, '
                             f'unindexed MF={unindexed}, all MF={all_mass}')
    if bool(indexed_eom) != bool(indexed_llss):
        raise ValueError(f'{model}: canonical and expanded LLSS MF representations differ')
    fmt = ('IndexedMFRequiresSingleMassConversion' if indexed_eom
           else 'SavedUnindexedSingleMassTemplate')
    return {'canonical_eom_GammaCC_terms': eom_gamma,
            'canonical_eom_indexed_MF_factors': indexed_eom,
            'canonical_eom_unindexed_MF_factors': unindexed_eom,
            'expanded_llss_GammaCC_terms': llss_gamma,
            'expanded_llss_indexed_MF_factors': indexed_llss,
            'expanded_llss_unindexed_MF_factors': unindexed_llss,
            'expanded_llss_reported_terms': expected_llss,
            'mf_source_mode': fmt,
            'three_non_degenerate_MF_mass_lift_inferred_not_source_matched': True,
            'canonical_EOMEFT_used_for_next_threshold': True}


def driver_anchors(source: str) -> dict[str, re.Match[str]]:
    """Fail closed on altered matching structure, not harmless whitespace.

    All matches are required to be entire lines where appropriate. In
    particular the scalar-stage early return must precede the FIRST-stage
    RunT3Matching call, and the old launch tail must be replaced exactly.
    """
    patterns = {
        'scriptdir': r'(?m)^[ \t]*scriptDir[ \t]*=[ \t]*DirectoryName[ \t]*@[ \t]*ExpandFileName\[\$InputFileName\];[ \t]*$',
        'module': r'(?m)^[ \t]*T3MEPerform\[\][ \t]*:=[ \t]*Module\[',
        'hook': r'(?m)^[ \t]*T3MECheckpoint\["fermion-first-production-Match-loader-started"\];[ \t]*$',
        'tail': r'(?m)^[ \t]*run[ \t]*=[ \t]*CheckAbort\[',
    }
    found = {}
    for name, pat in patterns.items():
        occurrences = list(re.finditer(pat, source))
        if len(occurrences) != 1:
            raise ValueError(f'Native F-first driver structural anchor {name}: '
                             f'expected once, found {len(occurrences)}')
        found[name] = occurrences[0]
    if not (found['module'].start() < found['hook'].start() < found['tail'].start()):
        raise ValueError('Native F-first driver structural stage ordering changed')
    body = source[found['hook'].end():found['tail'].start()]
    if body.count('RunT3Matching[lag,5,0]') != 1 or \
       'T3MECheckpoint["fermion-first-tree-Match-started"]' not in body:
        raise ValueError('Native F-first matching call no longer follows checkpoint')
    if 'T3MECheckpoint["fermion-first-lagrangian-checks-returned"]' not in \
       source[found['module'].start():found['hook'].start()]:
        raise ValueError('Cannot certify prior native CheckLagrangian stage')
    return found


def instrument(source: str, lagrangian_directory: Path) -> str:
    a = driver_anchors(source)
    # Apply edits right-to-left so the recorded source offsets remain valid.
    # Truncate the old execution tail: the first-stage matching must not run.
    result = source[:a['tail'].start()] + TAIL
    hook = a['hook']
    result = (result[:hook.start()] +
              hook.group(0).split('T3MECheckpoint',1)[0] +
              'Return[T3B6Hard[lag]];\n' + result[hook.start():])
    # Helper belongs in the native package context, before the procedure.
    module = a['module']
    result = result[:module.start()] + NATIVE_HELPER + '\n' + result[module.start():]
    # Only the temporary runner's scriptDir changes; original remains intact.
    scriptdir = a['scriptdir']
    directory = str(lagrangian_directory.resolve())
    result = (result[:scriptdir.start()] +
              'scriptDir = ExpandFileName[' + json.dumps(directory) + '];' +
              result[scriptdir.end():])
    # The original F matching call remains textually present in the builder,
    # but is unreachable after this early Return. It is never executed.
    if result.count('Return[T3B6Hard[lag]];') != 1 or \
       result.index('Return[T3B6Hard[lag]];') >= result.index('RunT3Matching[lag,5,0]'):
        raise ValueError('Failed to isolate scalar-only branch from F matching')
    return result


# The native B/C builder first validates the *original* scalar fields.
# It cannot pass CheckLagrangian when the extra off-diagonal S1/S2 mass
# operator is enabled. The validated mass-eigenstate driver intentionally
# adds the seven independent quartics itself AFTER its baseline construction.
# The completed synchronized EFT1 source also uses the default unindexed MF.
# Accordingly, reproduce the baseline driver's input environment rather than
# requesting a different, unrotated full-general UV model.
BASELINE_OPT_INS = (
    'T3_EXPERIMENTAL_GENERAL_QUADRATICS',
    'T3_EXPERIMENTAL_GENERAL_QUARTICS',
    'T3_EXPERIMENTAL_GENERAL_FERMIONS',
    'T3_EXPERIMENTAL_FLAVOR_MASSES',
)


def native_subprocess_environment(base: dict[str, str], canonical_eom: str) -> dict[str, str]:
    """Run the existing B/C mass-basis builder with its validated baseline.

    In particular, do NOT inject original-basis mixed scalar masses before
    constructing X/Y. This function does not alter the parent environment.
    """
    env = dict(base)
    for key in BASELINE_OPT_INS:
        env.pop(key, None)
    env['T3_BATCH6_EFT1_EOM'] = str(canonical_eom)
    return env


def write(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def run(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default='B_-1,C_-1')
    p.add_argument('--run-id',default='mixed_hard_01')
    p.add_argument('--preflight-only',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--timeout-seconds',type=int,default=3600)
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--sync-root',type=Path,default=SYNC)
    p.add_argument('--running-root',type=Path,default=RUNNING_ROOT)
    a=p.parse_args(argv)
    names=tuple(x.strip() for x in a.models.split(','))
    if not names or len(set(names))!=len(names) or any(m not in MODELS for m in names):
        p.error('Only unique B_-1,C_-1 can be selected')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',a.run_id) or not 120<=a.timeout_seconds<=7200:
        p.error('Invalid run ID or timeout')
    checks=[]
    for model in names:
        try:
            check=preflight(model,a.sync_root,a.running_root)
            print(f'{model}: READY, saved canonical EOMEFT; no F rerun',flush=True)
        except Exception as exc:
            check={'model':model,'status':'Blocked','reason':str(exc)}
            print(f'{model}: BLOCKED: {exc}',flush=True)
        checks.append(check)
    if a.preflight_only or any(x['status']!='Ready' for x in checks):
        return 0 if all(x['status']=='Ready' for x in checks) else 2
    base=OUTPUT/a.run_id
    cfg={'models':list(names),'timeout_seconds':a.timeout_seconds,
         'source_fingerprints':{x['model']:x['source_fingerprints'] for x in checks}}
    if base.exists():
        if not a.resume or not (base/'settings.json').is_file() or get_json(base/'settings.json')!=cfg:
            raise ValueError('Existing output; use --resume only for identical frozen source inputs')
    else:
        if a.resume: raise ValueError('Cannot resume nonexisting run')
        base.mkdir(parents=True)
        write(base/'settings.json',cfg)
    count=0
    for ch in checks:
        model=ch['model']; dest=base/tag(model); dest.mkdir(exist_ok=True)
        report=dest/'hard_scalar_report.json'
        if a.resume and report.is_file() and get_json(report).get('status')==STAGE:
            if not (dest/'hard_c5_candidate.txt').is_file():
                raise ValueError(f'{model}: successful report missing C5 candidate; refusing overwrite')
            print(model+': reused native scalar hard result')
            count+=1;continue
        native=instrument(DRIVER.read_text(encoding='utf-8-sig'),ROOT/'Lagrangian')
        temp=dest/'native_scalar_stage.generated.wl'
        temp.write_text(native,encoding='utf-8')
        # Never enable original-basis mixed quadratics here: the native B/C
        # driver expects an ordinary 11-vertex baseline, then constructs
        # the seven extra quartics and rotates S1/S2 -> X/Y itself.
        env=native_subprocess_environment(os.environ,ch['canonical_eom'])
        cmd=[a.wolframscript,'-file',str(temp),ch['parameters'],str(dest/'native_result.json')]
        try:
            with (dest/'stdout.log').open('w',encoding='utf-8') as so, \
                 (dest/'stderr.log').open('w',encoding='utf-8') as se:
                r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=so,stderr=se,
                                 timeout=a.timeout_seconds*2,check=False)
            result=get_json(report) if report.is_file() else {}
            if r.returncode!=0 or result.get('status')!=STAGE or \
                not (dest/'hard_c5_candidate.txt').is_file():
                native_path=dest/'native_result.json'
                progress_path=dest/'native_result.progress.json'
                native_report=get_json(native_path) if native_path.is_file() else {}
                progress=get_json(progress_path) if progress_path.is_file() else {}
                raise RuntimeError(f'{model}: native scalar stage incomplete, '
                                   f'return={r.returncode}, hard_report={result.get("status")}, '
                                   f'native_status={native_report.get("status")}, '
                                   f'native_reason={native_report.get("diagnostic",native_report.get("reason"))}, '
                                   f'last_stage={progress.get("last_stage")}; '
                                   f'see {dest / "stdout.log"}')
            count+=1; print(model+': native mixed X/Y hard C5 template exported',flush=True)
        except Exception as exc:
            write(dest/'error.json',{'status':'Blocked','model':model,'reason':str(exc),
                                     'traceback':traceback.format_exc(limit=4)})
            print(model+': BLOCKED: '+str(exc),flush=True)
    write(base/'summary.json',{'status':'HardCandidateExported' if count==len(names) else 'PartialOrBlocked',
            'model_count':len(names),'native_hard_successes':count,
            'fermion_matching_recomputed':False,'ordinary_models_rerun':False,
            'inherited_pole_cancellation_verified':False,'physical_C5_certified':False,
            'generated_utc':datetime.now(timezone.utc).isoformat()})
    return 0 if count==len(names) else 2


if __name__=='__main__':
    try: raise SystemExit(run())
    except Exception as exc:
        print('BATCH6 ERROR: '+str(exc),file=sys.stderr)
        raise SystemExit(2)
