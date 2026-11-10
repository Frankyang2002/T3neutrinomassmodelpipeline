(* Full-operator, source-preserving loop-function / explicit Laurent audit.
   Input MUST be the trusted local Matchete output already generated through
   RunSequentialT3Matching[... {{"F"}}, 5, 1, ...]. Both scalars remain light.

   This is NOT UV-hard-minus-EFT1-soft matching and NOT an MSbar Wilson
   subtraction. In particular, regular terms before counterterms are NOT
   numerical finite Wilson coefficients. No physics is silently substituted.

   wolframscript -file Lagrangian/RunT3Batch2LaurentAudit.wl \
       RAW_ONE_LOOP_INPUTFORM.txt OUTPUT_DIR SYMBOLIC_SECONDS
*)
ClearAll["Global`*"];
b2args = Rest[$ScriptCommandLine];
If[Length[b2args] =!= 3,
 Print["Usage: RAW_ONE_LOOP_INPUTFORM.txt OUTPUT_DIR SYMBOLIC_SECONDS"];
 Exit[2]];
b2source = ExpandFileName[b2args[[1]]];
b2dir = ExpandFileName[b2args[[2]]];
b2seconds = Quiet@Check[ToExpression[b2args[[3]], InputForm], $Failed];
If[!FileExistsQ[b2source] || !IntegerQ[b2seconds] ||
   b2seconds < 30 || b2seconds > 7200,
 Print["Missing source or invalid symbolic seconds (30..7200)"]; Exit[2]];
If[!DirectoryQ[b2dir], CreateDirectory[b2dir,CreateIntermediateDirectories->True]];
b2reportPath = FileNameJoin[{b2dir,"laurent_audit.json"}];
b2hash = FileHash[b2source,"SHA256","HexString"];
b2outputs = <||>;
ClearAll[B2LAWrite,B2LAStop,B2LAExport,B2LAText,B2LAPreview,
         B2LARegQ,B2LAHeads,B2LAUnresolvedQ,B2LACount];
B2LAWrite[data_Association] := Module[{rec},
 rec=Join[<|
   "format"->"t3-batch2-laurent-operator-audit-v1",
   "source_sha256"->b2hash,
   "source_file"->b2source,
   "complete_numerical_flavour_contractions_performed"->False,
   "source_scale_substitution_performed"->False,
   "uv_hard_minus_eft1_soft_subtraction_performed"->False,
   "local_MSbar_counterterms_added"->False,
   "renormalised_MSbar_Wilson_boundary_computed"->False,
   "full_eft1_wilson_RGE_performed"->False,
   "SMEFT_or_neutrino_prediction_computed"->False,
   "output_sha256"->b2outputs
  |>,data];
 Export[b2reportPath,rec,"RawJSON"];
 rec
];
B2LAStop[stage_String,reason_String] := (
 Print["Batch2 loop audit BLOCKED at ",stage,": ",reason];
 B2LAWrite[<|"status"->"Blocked_"<>stage,"reason"->reason,
   "one_loop_hbar_linearity_exact"->False,
   "all_operators_preserved"->False,"native_fermions_absent"->False|>];
 Exit[5]
);
B2LAText[x_] := ToString[InputForm[x]];
B2LAPreview[x_] := StringTake[B2LAText[x],UpTo[1200]];
B2LAExport[filename_String,expr_] := Module[{path,res},
 path=FileNameJoin[{b2dir,filename}];
 res=Quiet@Check[Export[path,expr,"WXF"],$Failed];
 If[res===$Failed || !FileExistsQ[path],
   B2LAStop["WXFExport","Failed exporting "<>filename]];
 AssociateTo[b2outputs, filename->FileHash[path,"SHA256","HexString"]];
 path
];
B2LARegQ[s_Symbol] := MemberQ[{
    "epsilon","eps","\\[Epsilon]","\[Epsilon]", "ε"},
    SymbolName[Unevaluated[s]]];
B2LARegQ[_] := False;
B2LAHeads[expr_] := DeleteDuplicates@Cases[expr,
    HoldPattern[f_Symbol[___]] :> SymbolName[Unevaluated[f]], Infinity];
B2LAUnresolvedQ[name_String] :=
  StringMatchQ[name, "Loop"~~__ | "Integral"~~__ |
    "LoopIntegral" | "A0" | "B0" | "C0" | "D0" |
    "Passarino"~~__ | "PV"~~__];
B2LACount[expr_] := Module[{e=Expand[expr]},
 If[TrueQ[e===0],0,If[Head[e]===Plus,Length[List@@e],1]]];

b2loaded=CheckAbort[Quiet@Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[b2loaded =!= True, B2LAStop["MatcheteLoad","Could not load Matchete"]];
(* Trusted source was produced by the user's own local Matchete driver. *)
b2raw = Import[b2source,"Text"];
If[!StringQ[b2raw] || StringLength[b2raw]===0,
 B2LAStop["Read","Missing InputForm expression"]];
b2held=CheckAbort[Quiet@Check[ToExpression[b2raw,InputForm,HoldComplete],$Failed],$Aborted];
If[Head[b2held]=!=HoldComplete,
 B2LAStop["Parse","Could not parse native generated InputForm"]];
b2expr=CheckAbort[ReleaseHold[b2held],$Aborted];
If[b2expr===$Failed || b2expr===$Aborted,
 B2LAStop["Parse","Could not construct Matchete native expression"]];
If[FreeQ[b2expr,hbar],
 B2LAStop["LoopBookkeeping","Missing Matchete hbar"]];
If[StringContainsQ[B2LAText[b2expr],"Field[NewFermion"],
 B2LAStop["Threshold","Integrated heavy fermion remains in one-loop EFT"]];

Print["Batch2 loop audit: extracting exact hbar^1 from ALL matched operators"];
b2expanded = TimeConstrained[CheckAbort[Expand[b2expr],$Aborted],b2seconds,$Aborted];
If[b2expanded===$Aborted,
 B2LAStop["Expand","Full matched expression expansion timed out"]];
b2linear = Coefficient[b2expanded,hbar,1];
b2difference = Expand[b2expanded - hbar*b2linear];
If[!TrueQ[b2difference===0] || TrueQ[b2linear===0],
 B2LAStop["LoopBookkeeping","Correction is not exactly hbar-linear and nonzero"]];
B2LAExport["raw_hbar1_coefficient.wxf",b2linear];

(* Run Matchete's existing loop evaluator; unlike the raw InputForm export,
   this can expose explicit dimensional poles and finite loop functions. *)
Print["Batch2 loop audit: EvaluateLoopFunctions on full hbar^1 correction"];
b2evaluated = TimeConstrained[
 CheckAbort[Quiet@Check[UsingFrontEnd[EvaluateLoopFunctions[b2linear]],$Failed],$Aborted],
 b2seconds,$Aborted];
If[MemberQ[{$Failed,$Aborted},b2evaluated],
 B2LAStop["LoopEvaluation","Matchete loop-function evaluation failed or timed out"]];
If[!FreeQ[b2evaluated,hbar],
 B2LAStop["LoopEvaluation","Loop evaluator unexpectedly reintroduced hbar"]];
If[StringContainsQ[B2LAText[b2evaluated],"Field[NewFermion"],
 B2LAStop["LoopEvaluation","Heavy fermion reappeared after evaluation"]];
B2LAExport["evaluated_full_hbar1_coefficient.wxf",b2evaluated];

b2heads = B2LAHeads[b2evaluated];
b2preheads = Select[B2LAHeads[b2linear],B2LAUnresolvedQ];
b2unresolved = Select[b2heads,B2LAUnresolvedQ];
b2regs = DeleteDuplicates@Cases[b2evaluated,
   (s_Symbol /; B2LARegQ[s]) :> s, Infinity];
If[Length[b2unresolved]>0,
 B2LAWrite[<|
  "status"->"Blocked_UnresolvedLoopFunctions",
  "reason"->"Matchete evaluation retains potential unreduced loop-function heads",
  "unresolved_function_heads"->b2unresolved,
  "initial_loop_function_heads"->b2preheads,
  "regulator_symbol_count"->Length[b2regs],
  "one_loop_hbar_linearity_exact"->True,
  "all_operators_preserved"->True,"native_fermions_absent"->True|>];
 Print["Unresolved loop heads: ",b2unresolved]; Exit[5]
];
If[Length[b2regs]>1,
 B2LAWrite[<|
  "status"->"Blocked_AmbiguousRegulator",
  "reason"->"Multiple epsilon symbols must be reconciled before Laurent expansion",
  "regulators_inputform"->(B2LAText /@ b2regs),
  "one_loop_hbar_linearity_exact"->True,
  "all_operators_preserved"->True,"native_fermions_absent"->True|>];
 Exit[5]
];
If[Length[b2regs]===0,
 (* Do not call an expression with no explicit epsilon a finite MSbar result:
    the full UV/EFT1 subtraction and renormalisation are still absent. *)
 B2LAWrite[<|
  "status"->"LoopFunctionsEvaluated_NoExplicitRegulator_AuditOnly",
  "one_loop_hbar_linearity_exact"->True,
  "all_operators_preserved"->True,"native_fermions_absent"->True,
  "regulator_symbol_inputform"->Null,
  "explicit_pole_residue_present"->False,
  "evaluated_expression_term_count"->B2LACount[b2evaluated],
  "initial_loop_function_heads"->b2preheads,
  "raw_evaluated_preview_inputform"->B2LAPreview[b2evaluated],
  "note"->"No explicit epsilon does NOT prove finiteness of UV-hard minus EFT1-soft"|>];
 Print["Batch2 loop audit: Matchete evaluated, no explicit epsilon; not finite matching"];
 Exit[0]
];

b2eps=First[b2regs];
Print["Batch2 loop audit: formal epsilon^(-1) / epsilon^0 separation"];
(* Do not silently drop higher poles or regulator-dependent finite terms. *)
b2series=TimeConstrained[
 CheckAbort[Quiet@Check[With[{e=b2eps},Series[b2evaluated,{e,0,0}]],$Failed],$Aborted],
 b2seconds,$Aborted];
If[MemberQ[{$Failed,$Aborted},b2series],
 B2LAWrite[<|"status"->"Blocked_LaurentExpansion",
  "reason"->"Could not obtain complete explicit epsilon Laurent series",
  "regulator_symbol_inputform"->B2LAText[b2eps],
  "one_loop_hbar_linearity_exact"->True,
  "all_operators_preserved"->True,"native_fermions_absent"->True|>];
 Exit[5]
];
b2normal=Expand[Normal[b2series]];
b2pole=With[{e=b2eps},Coefficient[b2normal,e,-1]];
b2regular=With[{e=b2eps},Coefficient[b2normal,e,0]];
b2remainder=Expand[b2normal-b2pole/b2eps-b2regular];
If[!TrueQ[b2remainder===0],
 B2LAWrite[<|"status"->"Blocked_UnresolvedHigherPoles",
  "reason"->"Simple pole plus constant does not reconstruct Laurent series",
  "regulator_symbol_inputform"->B2LAText[b2eps],
  "series_remainder_preview_inputform"->B2LAPreview[b2remainder],
  "one_loop_hbar_linearity_exact"->True,
  "all_operators_preserved"->True,"native_fermions_absent"->True|>];
 Exit[5]
];
B2LAExport["explicit_pole_residue_unsubtracted.wxf",b2pole];
B2LAExport["epsilon_zero_regular_unsubtracted.wxf",b2regular];
b2report=B2LAWrite[<|
 "status"->"LoopFunctionsEvaluated_ExplicitLaurentSeparated_AuditOnly",
 "one_loop_hbar_linearity_exact"->True,
 "all_operators_preserved"->True,"native_fermions_absent"->True,
 "regulator_symbol_inputform"->B2LAText[b2eps],
 "explicit_pole_residue_present"->(!TrueQ[b2pole===0]),
 "explicit_laurent_reconstruction_exact"->True,
 "evaluated_expression_term_count"->B2LACount[b2evaluated],
  "initial_loop_function_heads"->b2preheads,
 "raw_pole_residue_term_count"->B2LACount[b2pole],
 "raw_regular_term_count"->B2LACount[b2regular],
 "pole_preview_inputform"->B2LAPreview[b2pole],
 "regular_preview_inputform"->B2LAPreview[b2regular],
 "regular_part_is_physically_finite_Wilson"->False|>];
Print["Batch2 loop audit: explicit Laurent extraction complete (NOT MSbar matching)"];
Exit[0];
