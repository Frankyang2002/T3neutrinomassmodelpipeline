(* Workstream B: live Matchete SU(2) fundamental generators for B(-1)/C(-1).
   Does not impose an arbitrary cross-package field rephasing.
   wolframscript -file scripts/ProbeT3BCSU2GeneratorsMatchete.wl NEW_OUTPUT.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ProbeT3BCSU2GeneratorsMatchete.wl NEW_OUTPUT.json"]; Exit[64]];
out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out], Print["Refusing to overwrite ",out]; Exit[64]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
loaded = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True],$Failed],$Aborted];
If[loaded =!= True, Print["Matchete with FrontEnd unavailable"]; Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];

T3BGFail[label_String,why_String] := Throw[<|
  "model"->label,"status"->"Blocked","error"->why,
  "relative_field_phase_certified"->False,
  "physical_mapping_certified"->False|>,"T3BG"];
T3BGOne[label_String,cls_String] := Module[
  {model,sm,fields,raw,gens,err,hermiticity,traceNorm,comm,flat},
  Print["Live Matchete SU(2) generators for ",label];
  model = T3ModelFromClass[cls,-1];
  If[!AssociationQ[model] || model["Alpha"] =!= -1 ||
     model["Scalar1","SU2"] =!= 2 || model["Scalar2","SU2"] =!= 2,
     T3BGFail[label,"Unexpected scalar representation"]];
  ResetAll[];
  sm = CheckAbort[Check[LoadModel["SM"],$Failed],$Aborted];
  If[sm === $Failed || sm === $Aborted,
     T3BGFail[label,"SM model load failed"]];
  fields = CheckAbort[Check[DefineT3Fields[model],$Failed],$Aborted];
  If[fields === $Failed || fields === $Aborted,
     T3BGFail[label,"Matchete field registration failed"]];
  (* Generators[SU[2],{1}] is Matchete's public GroupMagic accessor.
     No hard-coded Pauli matrices substitute for its returned values. *)
  raw = CheckAbort[Check[Generators[SU[2],{1}],$Failed],$Aborted];
  If[raw === $Failed || raw === $Aborted,
     T3BGFail[label,"Matchete native generator accessor failed"]];
  gens = If[Head[raw] === SparseArray,Normal[raw],raw];
  If[Dimensions[gens] =!= {3,2,2} ||
     !AllTrue[Flatten[gens], NumericQ[FullSimplify[#]]&],
     T3BGFail[label,"Native generators do not form three numeric 2x2 arrays"]];
  hermiticity = Max[Abs[Flatten[N[gens - (ConjugateTranspose /@ gens),16]]]];
  traceNorm = Max[Abs[Flatten[N[Table[
       Tr[gens[[i]].gens[[j]]] - KroneckerDelta[i,j]/2,
       {i,3},{j,3}],16]]]];
  comm = Max[Abs[Flatten[N[Table[
     gens[[a]].gens[[b]]-gens[[b]].gens[[a]] -
       I Sum[Signature[{a,b,c}] gens[[c]],{c,3}],
     {a,3},{b,3}],16]]]];
  If[!And@@((NumberQ[#] && # < 10^-11) & /@ {hermiticity,traceNorm}),
    T3BGFail[label,"Native generator hermiticity/normalisation failed"]];
  <|
    "model"->label,
    "status"->"LiveMatcheteSU2GeneratorsCaptured",
    "accessor"->"Matchete`Generators[SU[2], {1}]",
    "generator_dimensions"->{3,2,2},
    "generator_real"->Map[N[Re[#],17]&,gens,{3}],
    "generator_imag"->Map[N[Im[#],17]&,gens,{3}],
    "native_hermiticity_residual"->N[hermiticity],
    "native_trace_normalization_residual"->N[traceNorm],
    "native_pauli_commutator_residual"->N[comm],
    "relative_field_phase_certified"->False,
    "physical_mapping_certified"->False
  |>
];
T3BGRun[] := <|
  "format"->"t3-bc-live-matchete-su2-generators-v1",
  "status"->"LiveMatcheteSU2GeneratorsOnly",
  "models"-><|
   "B_-1"->Catch[T3BGOne["B_-1","B"],"T3BG"],
   "C_-1"->Catch[T3BGOne["C_-1","C"],"T3BG"]|>,
  "rgbeta_numeric_generators_extracted"->False,
  "relative_field_phase_certified"->False,
  "physical_mapping_certified"->False,
  "production_modified"->False|>;
result=CheckAbort[Check[UsingFrontEnd[T3BGRun[]],$Failed],$Aborted];
If[!AssociationQ[result],Print["Matchete generator probe failed"];Exit[5]];
If[!DirectoryQ[DirectoryName[out]],
   CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
If[Check[Export[out,result,"RawJSON"],$Failed] === $Failed || !FileExistsQ[out],
   Print["Matchete generator report write failed"];Exit[6]];
Print["Matchete SU2 generator diagnostic: ",out];
If[!AllTrue[Values[result["models"]],
      Lookup[#,"status",""] === "LiveMatcheteSU2GeneratorsCaptured"&],
   Print["Native SU(2) generator capture blocked; JSON retained"];Exit[7]];
Print["Both models: three native SU2 generators captured. Relative field phase NOT fixed."];
