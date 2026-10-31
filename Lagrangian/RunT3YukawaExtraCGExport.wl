(* Read-only Workstream A: export the ACTUALLY REGISTERED extra-Yukawa CG.
   Usage: wolframscript -file Lagrangian/RunT3YukawaExtraCGExport.wl OUTPUT.json A m2
   Supported: A/D/E x (m2,p0). Requires local Matchete installation.
   This does not run matching, RGBeta, or any production RGE.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] < 3,
 Print["Usage: OUTPUT.json A|D|E m2|p0"]; Exit[2]];
args = Take[$ScriptCommandLine, -3];
out = ExpandFileName[args[[1]]];
cls = ToUpperCase[args[[2]]]; tag = args[[3]];
If[!MemberQ[{"A", "D", "E"}, cls] || !MemberQ[{"m2", "p0"}, tag],
 Print["Unsupported model"]; Exit[2]];
If[Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_FLAVOR_MASSES"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
 Print["Expected general fermions, flavour masses, quartics=1; quadratics !=1"]; Exit[2]];

alpha = If[tag === "m2", -2, 0]; label = cls <> "_" <> tag;
rootL = DirectoryName @ ExpandFileName[$InputFileName];
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[load =!= True, Print["Matchete load failed"]; Exit[3]];
Get[FileNameJoin[{rootL, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{rootL, "model", "LagrangianBuilder.wl"}]];

ClearAll[WardCGRead, WardCGFail, WardCGRun];
WardCGFail[msg_String] := Throw[<|"status" -> "ExtraCGExtractionBlocked", "reason" -> msg|>, "WardCGStop"];
WardCGRead[name_Symbol, dims_List, interaction_] := Module[{arr, raw, flat, shape},
 If[FreeQ[interaction, HoldPattern[CG[name, _List]]],
   WardCGFail[SymbolName[name] <> " not present in assembled interaction"]];
 raw = CheckAbort[Quiet@Check[GetCGTensor[name], $Failed], $Aborted];
 If[raw === $Failed || raw === $Aborted, WardCGFail["GetCGTensor failed for " <> SymbolName[name]]];
 arr = If[Head[raw] === SparseArray, Normal[raw], raw];
 shape = Dimensions[arr]; flat = Flatten[arr];
 If[shape =!= dims || Length[flat] =!= Times @@ dims ||
    !AllTrue[flat, TrueQ[NumericQ[FullSimplify[#]]] &] ||
    !AnyTrue[flat, !TrueQ[PossibleZeroQ[#]] &],
   WardCGFail["Malformed CG tensor for " <> SymbolName[name]]];
 <|"symbol" -> SymbolName[name], "dimensions" -> dims,
   "components_inputform" -> (ToString[FullSimplify[#], InputForm] & /@ flat),
   "registered_in_assembled_interaction" -> True,
   "tensor_accessor" -> "Matchete`GetCGTensor"|>
];
WardCGRun[] := Module[{model, built, df, ds, which, cBase, cExtra, cdims,
                        baseline, extra, same},
 model = T3ModelFromClass[cls, alpha];
 If[!AssociationQ[model], WardCGFail["Model construction failed"]];
 built = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]], $Failed], $Aborted];
 If[!AssociationQ[built] || Lookup[built, "Status", ""] =!= "Success" ||
    !TrueQ[Lookup[built, "FullValidation", False]] ||
    Lookup[built, "AdditionalYukawasIncluded", 0] =!= 1,
  WardCGFail["Matchete full build did not accept exactly one extra Yukawa"]];
 which = If[alpha === 0, 1, 2];
 df = model["Fermion", "SU2"];
 ds = model[If[which === 1, "Scalar1", "Scalar2"], "SU2"];
 cdims = Select[{2, df, ds}, # > 1 &];
 cBase = If[which === 1, T3Y1CG, T3Y2CG];
 cExtra = If[which === 1, T3Y1ConjCG, T3Y2ConjCG];
 baseline = WardCGRead[cBase, cdims, built["LInt"]];
 extra = WardCGRead[cExtra, cdims, built["LInt"]];
 same = TrueQ[And @@ MapThread[TrueQ[FullSimplify[#1 - #2 == 0]] &,
   {Flatten[If[Head[GetCGTensor[cBase]] === SparseArray,
     Normal[GetCGTensor[cBase]], GetCGTensor[cBase]]],
    Flatten[If[Head[GetCGTensor[cExtra]] === SparseArray,
     Normal[GetCGTensor[cExtra]], GetCGTensor[cExtra]]]}]];
 (* This is a test of native tensor EQUALITY, not a statement about the
    equality of the complex-scalar interactions or their physical couplings. *)
 <|"status" -> If[same, "LiveExtraCGEqualityVerified", "LiveExtraCGNotEqualToBaseline"],
   "model" -> label, "alpha" -> alpha, "dimensions" -> {2, df, ds},
   "extra_yukawa" -> If[which === 1, "y1Conj", "y2Conj"],
   "baseline_yukawa" -> If[which === 1, "Yukawa1", "Yukawa2"],
   "baseline_tensor" -> baseline, "extra_tensor" -> extra,
   "extra_equals_baseline_components_exactly" -> same,
   "full_matchete_check_passed" -> True,
   "physical_matchete_rgbeta_normalization_certified" -> False,
   "production_modified" -> False|>
];
result = CheckAbort[Catch[WardCGRun[], "WardCGStop"], $Aborted];
If[!AssociationQ[result], result = <|"status" -> "ExtraCGExtractionBlocked",
 "reason" -> "Wolfram evaluation failed or aborted"|>];
result = Join[<|"format" -> "t3-live-extra-yukawa-cg-v1", "model" -> label|>, result];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out], CreateIntermediateDirectories -> True]];
If[Quiet@Check[Export[out, result, "RawJSON"], $Failed] === $Failed,
 Print["Failed to write ", out]; Exit[7]];
Print[label, ": ", result["status"]];
If[result["status"] =!= "LiveExtraCGEqualityVerified", Exit[6]];
Exit[0];
