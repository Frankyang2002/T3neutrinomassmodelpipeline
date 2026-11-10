(* Source-preserving inventory of the COMPLETE, already Matchete-matched
   one-loop first-threshold EFT.  No hard topology is selected/dropped.
   This is not MSbar subtraction or an EFT1 ADM calculation.

   wolframscript -file Lagrangian/RunT3Batch2OneLoopLedger.wl \
     MATCHING/data/eft1_after_F_one_loop_inputform.txt OUTPUT.json
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2, Print["Usage: INPUT_ONE_LOOP_INPUTFORM OUTPUT_LEDGER.json"]; Exit[2]];
source = ExpandFileName[args[[1]]];
out = ExpandFileName[args[[2]]];
If[!FileExistsQ[source], Print["Missing one-loop Matchete source: ",source]; Exit[3]];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
ClearAll[B2Write,B2Fail,B2Text,B2SumFields,B2Sector];
B2Write[x_Association] := Export[out,x,"RawJSON"];
B2Fail[msg_String] := (
 Print["Batch 2 ledger ERROR: ",msg];
 B2Write[<|"format"->"t3-batch2-complete-loop-operator-ledger-v1",
    "status"->"Blocked","reason"->msg,
    "finite_msbar_wilson_coefficients_computed"->False|>]; Exit[5]
);
B2Text[x_] := ToString[InputForm[x]];
B2SumFields[s_String] := StringCount[s,"Field[NewScalar1"] +
  StringCount[s,"Field[NewScalar2"] +
  StringCount[s,"Field[T3Scalar1Minus"] +
  StringCount[s,"Field[T3Scalar1Plus"] +
  StringCount[s,"Field[T3Scalar2Minus"] +
  StringCount[s,"Field[T3Scalar2Plus"];
B2Sector[s_String] := Which[
 StringContainsQ[s,"GammaCC"] &&
    StringCount[s,"Field[H,"]>=2 && B2SumFields[s]===0,
    "CandidateLeptonNumberViolatingLLHH",
 StringContainsQ[s,"GammaCC"] && B2SumFields[s]>=2,
    "CandidateLeptonNumberViolatingLLScalarScalar",
 StringContainsQ[s,"GammaCC"], "OtherGammaCCOperator",
 True, "OtherEFT1Operator"
];

loaded = CheckAbort[Quiet@Check[UsingFrontEnd[Needs["Matchete`"]; True],$Failed],$Aborted];
If[loaded =!= True,B2Fail["Matchete is required to read native operator symbols"]];
raw = Import[source,"Text"];
If[!StringQ[raw] || StringLength[raw]===0,B2Fail["Empty raw Matchete source"]];
(* InputForm was just exported by RunSequentialT3Matching in the same local
   repository; it is a trusted generated source, not arbitrary user input. *)
held = CheckAbort[Quiet@Check[ToExpression[raw,InputForm,HoldComplete],$Failed],$Aborted];
If[Head[held]=!=HoldComplete,B2Fail["Cannot parse trusted generated InputForm"]];
expr = CheckAbort[Quiet@Check[ReleaseHold[held],$Failed],$Aborted];
If[expr===$Failed || expr===$Aborted, B2Fail["Cannot read Matchete one-loop expression"]];
If[FreeQ[expr,hbar] && !StringContainsQ[raw,"hbar"],
 B2Fail["Missing one-loop Matchete bookkeeping hbar"]];
expanded = CheckAbort[Expand[expr],$Aborted];
If[expanded===$Aborted,B2Fail["Expansion of complete one-loop expression aborted"]];
terms = If[expanded===0,{},If[Head[expanded]===Plus,List@@expanded,{expanded}]];
If[Length[terms]===0,B2Fail["Matchete has no nonzero one-loop terms"]];
reconstructed = CheckAbort[Expand[Total[terms] - expanded],$Aborted];
If[reconstructed===$Aborted || !TrueQ[reconstructed===0],
 B2Fail["Exact term reconstruction of complete EFT failed"]];
texts = B2Text /@ terms;
If[AnyTrue[texts,StringContainsQ[#,"Field[NewFermion"]&],
 B2Fail["Heavy fermion field appears in matched EFT1"]];
sectors = B2Sector /@ texts;
records = MapIndexed[
 Function[{s,pos},<|
  "index"->First[pos],
  "sector"->sectors[[First[pos]]],
  "raw_inputform"->s,
  "contains_dim_reg_pole"->(StringContainsQ[s,"\\[Epsilon]^(-1)"] ||
      StringContainsQ[s,"1/\\[Epsilon]"] ||
      StringContainsQ[s,"epsilon^(-1)"] ||
      StringContainsQ[s,"1/epsilon"]),
  "contains_hbar"->StringContainsQ[s,"hbar"],
  "gamma_cc"->StringContainsQ[s,"GammaCC"],
  "coupling_symbols"->DeleteDuplicates@Cases[terms[[First[pos]]],
    HoldPattern[head_Symbol[name_Symbol,___]] /;
      SymbolName[head]==="Coupling" :> SymbolName[name],Infinity],
  "light_scalar_legs_count"->B2SumFields[s],
  "higgs_legs_count"->StringCount[s,"Field[H,"],
  "conjugate_scalar_fields_present"->StringContainsQ[s,"Bar[Field[NewScalar"]
 |>],texts];
counts = Counts[sectors];
report = <|
 "format"->"t3-batch2-complete-loop-operator-ledger-v1",
 "status"->"UnrenormalisedOneLoopTermsEnumerated",
 "source_file"->source,
 "source_sha256"->FileHash[source,"SHA256","HexString"],
 "term_count"->Length[terms],
 "terms"->records,
 "sector_counts"->counts,
 "exact_reconstruction_by_expansion"->True,
 "only_fermions_integrated_out"->True,
 "unselected_terms_retained"->True,
 "selection_for_physical_C5_performed"->False,
 "complete_numeric_flavour_contractions_performed"->False,
 "uv_minus_eft1_soft_subtraction_performed"->False,
 "finite_msbar_wilson_coefficients_computed"->False,
 "eft1_wilson_one_loop_rge_performed"->False,
 "pole_terms_identified_heuristically_only"->True,
 "matchete_to_rgbeta_basis_converted"->False
|>;
If[Quiet@Check[B2Write[report],$Failed]===$Failed,
 Print["Failed to export complete raw operator ledger"]; Exit[5]];
Print["Batch2 all-term ledger: ",Length[terms]," terms (none dropped) -> ",out];
Exit[0];
