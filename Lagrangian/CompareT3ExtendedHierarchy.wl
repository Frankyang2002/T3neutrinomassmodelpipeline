(* Compare seeded versus UV-zeroed native Matchete matching artifacts.
   This compares actual EFT expressions; it is not an EFT1 running calculation.
   Usage: wolframscript -file Lagrangian/CompareT3ExtendedHierarchy.wl
          SEEDED_DIR ZERO_DIR OUT_DIR stage1|hierarchy *)
ClearAll["Global`*"];
cArgs = Rest[$ScriptCommandLine];
If[Length[cArgs] =!= 4 ||
   !MemberQ[{"stage1","hierarchy"},cArgs[[4]]],
   Print["Usage: SEEDED_DIR ZERO_DIR OUT_DIR stage1|hierarchy"]; Exit[2]];
cSeed = ExpandFileName[cArgs[[1]]];
cZero = ExpandFileName[cArgs[[2]]];
cOut = ExpandFileName[cArgs[[3]]];
cMode = cArgs[[4]];
If[!DirectoryQ[cOut],CreateDirectory[cOut,CreateIntermediateDirectories->True]];
cLoad = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[cLoad =!= True,Print["Native Matchete not loaded"]; Exit[3]];
(* Body parsed only after package load, avoiding premature Matchete symbols. *)
cResult = CheckAbort[Module[
  {a,b,required,s,z,d,files=<||>,notProvenZero=<||>,i,
   r,status,base,full,coefficient,nonzero},
  a = Quiet@Check[Import[FileNameJoin[{cSeed,"native.json"}],"RawJSON"],$Failed];
  b = Quiet@Check[Import[FileNameJoin[{cZero,"native.json"}],"RawJSON"],$Failed];
  If[!AssociationQ[a] || !AssociationQ[b], Return[<|
    "status"->"ComparisonBlocked","reason"->"Missing native JSON"|>]];
  If[Lookup[a,"model",None] =!= Lookup[b,"model",None] ||
     Lookup[a,"variant",None] =!= "seeded" ||
     Lookup[b,"variant",None] =!= "zero" ||
     Lookup[a,"mode",None] =!= cMode ||
     Lookup[b,"mode",None] =!= cMode ||
     Lookup[a,"status",""] =!= If[cMode === "stage1",
       "FermionEFT1NativeMatched","NativeHierarchicalC5ExtractedNotRenormalized"] ||
     Lookup[b,"status",""] =!= Lookup[a,"status",None],
    Return[<|"status"->"ComparisonBlocked",
      "reason"->"Source model/variant/mode/success mismatch"|>]];
  required = If[cMode === "stage1",
    {"stage1_tree","stage1_one_loop","stage1_full"},
    {"stage1_tree","stage1_one_loop","stage1_full",
      "stage2_full","stage2_c5_ordered"}];
  Do[
    s = Quiet@Check[Import[FileNameJoin[{cSeed,i<>".wxf"}],"WXF"],$Failed];
    z = Quiet@Check[Import[FileNameJoin[{cZero,i<>".wxf"}],"WXF"],$Failed];
    If[s === $Failed || z === $Failed,
      Return[<|"status"->"ComparisonBlocked",
        "reason"->"Missing seeded/zero WXF for "<>i|>]];
    d = Expand[s-z];
    If[!FreeQ[d, _Missing] || d === $Failed,
      Return[<|"status"->"ComparisonBlocked",
        "reason"->"Unresolved difference for "<>i|>]];
    files[i] = FileNameJoin[{cOut,i<>"_seeded_minus_zero.wxf"}];
    If[Quiet@Check[Export[files[i],d,"WXF"],$Failed] === $Failed,
      Return[<|"status"->"ComparisonBlocked",
        "reason"->"Cannot export difference for "<>i|>]];
    notProvenZero[i] = !TrueQ[PossibleZeroQ[d]];
    If[i === "stage2_c5_ordered",
      Export[FileNameJoin[{cOut,"stage2_c5_difference.inputform.txt"}],
        ToString[d,InputForm],"Text"]
    ];,
    {i,required}];
  <|"status"->"NativeEFTDifferenceExported",
    "model"->a["model"],"mode"->cMode,
    "seeded_vs_uv_zero_expressions_compared"->required,
    "difference_not_proven_zero_by_PossibleZeroQ"->notProvenZero,
    "cross_kernel_dummy_index_canonicalisation_certified"->False,
    "differences_wxf"->files,
    "after_full_native_sequential_matching"->(cMode === "hierarchy"),
    "difference_zero_by_PossibleZeroQ_is_a_sufficient_but_not_necessary_test"->True,
    "absolute_coupling_dictionary_certified"->False,
    "full_renormalized_hierarchical_C5_certified"->False,
    "EFT1_running_inserted"->False|>
],$Aborted];
If[!AssociationQ[cResult],cResult=<|"status"->"ComparisonBlocked",
  "reason"->"Comparison aborted without an association"|>];
Export[FileNameJoin[{cOut,"difference.json"}],cResult,"RawJSON"];
Print["Native EFT comparison: ",Lookup[cResult,"status","Unknown"]];
If[Lookup[cResult,"status",""] =!= "NativeEFTDifferenceExported", Exit[6]];
Exit[0];
