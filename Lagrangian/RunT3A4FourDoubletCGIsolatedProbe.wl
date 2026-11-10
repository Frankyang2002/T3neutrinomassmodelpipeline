(* A(-4) INDEPENDENT Matchete CG validator. No T3 source, matching,
   heavy-field changes or Wolfram package internals are overwritten.

   wolframscript -file Lagrangian/RunT3A4FourDoubletCGIsolatedProbe.wl OUTPUT.json

   The six nonzero entries and the four-bar representation list are copied
   verbatim from the user's 2026-11-10 v9 DefineCG trace. The exact tensor
   is neither a physical Wilson coefficient nor a permitted replacement CG.
*)
ClearAll["Global`*"];
t3Args = Rest[$ScriptCommandLine];
If[Length[t3Args] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
t3Output = ExpandFileName[First[t3Args]];
If[!DirectoryQ[DirectoryName[t3Output]],
  CreateDirectory[DirectoryName[t3Output],
    CreateIntermediateDirectories -> True]];
ClearAll[t3String, t3Status, t3Export];
t3String[e_] := ToString[e, InputForm];
t3Status = <|
  "format" -> "t3-a4-isolated-four-doublet-matchete-cg-probe-v10",
  "physics_source" -> "Exact six-component bad CG from v9 DefineCG trace",
  "bad_tensor_is_su2_invariant" -> False,
  "hypothesis_isolated" -> "Source convention vs Matchete CGManipulations undecided",
  "matchete_code_modified" -> False,
  "full_one_loop_Match_replaced" -> False,
  "physical_C5_certified" -> False|>;
t3Export[] := Export[t3Output, t3Status, "RawJSON"];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"];
      ResetAll[]; LoadModel["SM"]; True], $Aborted] =!= True,
  t3Status["status"] = "Blocked_CouldNotLoadMatcheteSM";
  t3Export[]; Exit[3]];
t3Reps = ConstantArray[Bar[SU2L[fund]], 4];
t3Bad = SparseArray[{
  {1,1,1,1} -> 2/Sqrt[3],
  {1,2,1,2} -> 1/Sqrt[3],
  {1,2,2,1} -> 1/Sqrt[3],
  {2,1,1,2} -> 1/Sqrt[3],
  {2,1,2,1} -> 1/Sqrt[3],
  {2,2,2,2} -> 2/Sqrt[3]
}, {2,2,2,2}];
t3Status["bad_component_1111"] = t3String[t3Bad[[1,1,1,1]]];
t3Status["bad_component_2222"] = t3String[t3Bad[[2,2,2,2]]];
t3Status["bad_tensor_reps"] = t3String[t3Reps];
t3Status["bad_tensor_L2_norm_squared"] = t3String[Total[Flatten[Normal[t3Bad]]^2]];
(* Matchete's API rejects the invalid tensor with DefineCG::notInvariant.
   Capture Abort[] in a controlled stand-alone process. *)
t3BadCall = TimeConstrained[CheckAbort[
  Check[UsingFrontEnd[
    DefineCG[T3A4BadFourBarProbeCG, t3Reps, t3Bad]; True],
    $Failed], $Aborted], 120, $TimedOut];
t3Status["bad_DefineCG_result"] = t3String[t3BadCall];
t3Status["bad_DefineCG_failed"] = MemberQ[{$Failed,$Aborted,$TimedOut},t3BadCall];
(* Independent healthy SU(2) invariant construction, same four slots.
   If this fails too, registry mismatch is more likely than a single bad
   composite tensor. Otherwise a genuine singlet can be registered. *)
t3AlgebraReps = ConstantArray[{1},4];
t3Good = TimeConstrained[CheckAbort[UsingFrontEnd[
  InvariantTensors[SU[2],t3AlgebraReps]], $Aborted],120,$TimedOut];
t3Status["valid_singlet_basis_count"] = If[ListQ[t3Good],Length[t3Good],-1];
t3Status["valid_singlet_basis_result_head"] = t3String[Head[t3Good]];
If[ListQ[t3Good] && Length[t3Good] >= 1,
  t3GoodCall = TimeConstrained[CheckAbort[
    Check[UsingFrontEnd[DefineCG[T3A4GoodFourBarProbeCG,
      t3Reps,First[t3Good]]; True],$Failed],$Aborted],120,$TimedOut];
  t3Status["good_DefineCG_result"] = t3String[t3GoodCall];
  t3Status["good_DefineCG_returned_without_error"] = TrueQ[t3GoodCall],
  t3Status["good_DefineCG_result"] = "NotRun";
  t3Status["good_DefineCG_returned_without_error"] = False
];
t3Status["status"] = If[TrueQ[t3Status["bad_DefineCG_failed"]] &&
    TrueQ[t3Status["good_DefineCG_returned_without_error"]],
  "BadCGRejected_IndependentSU2SingletAccepted_SourceUnresolved",
  "IsolatedProbeRequiresRegistryOrVersionInvestigation"];
t3Export[];
Print["A(-4) isolated CG test: ",t3Status["status"]];
Print["Output: ",t3Output];
If[t3Status["status"] =!=
  "BadCGRejected_IndependentSU2SingletAccepted_SourceUnresolved",Exit[4]];
Exit[0];
