(* Native, standalone one-loop Matchete extraction for the six extra-Yukawa T3 models.
   No RGBeta->Matchete identification is made here. This is a full one-step
   heavy-particle matching pilot, NOT the hierarchical two-threshold EFT.
   Usage: wolframscript -file Lagrangian/RunT3ExtendedC5Native.wl OUT.json A|D|E m2|p0 2|3
   Must enable GENERAL_FERMIONS, GENERAL_QUARTICS, FLAVOR_MASSES and disable GENERAL_QUADRATICS.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 4,
  Print["Usage: OUT.json A|D|E m2|p0 2|3"]; Exit[2]];
out = ExpandFileName[args[[1]]];
cl = ToUpperCase[args[[2]]];
tag = args[[3]];
nHeavy = If[MemberQ[{"2", "3"}, args[[4]]], ToExpression[args[[4]]], $Failed];
allowed = <|"A_m2" -> -2, "A_p0" -> 0, "D_m2" -> -2,
  "D_p0" -> 0, "E_m2" -> -2, "E_p0" -> 0|>;
key = cl <> "_" <> tag;
If[!KeyExistsQ[allowed,key] || nHeavy === $Failed,
  Print["Unsupported model or heavy-flavour count: ",key," / ",args[[4]]];Exit[2]];
If[!And @@ (Environment[#] === "1" & /@ {
    "T3_EXPERIMENTAL_GENERAL_FERMIONS",
    "T3_EXPERIMENTAL_GENERAL_QUARTICS",
    "T3_EXPERIMENTAL_FLAVOR_MASSES"}) ||
    Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] === "1",
  Print["Expected GENERAL_FERMIONS/GENERAL_QUARTICS/FLAVOR_MASSES=1 and GENERAL_QUADRATICS=0"];
  Exit[2]];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]];
dir = DirectoryName @ ExpandFileName[$InputFileName];
progress = FileNameJoin[{DirectoryName[out],FileBaseName[out]<>".progress.json"}];
step = "initialize";
ClearAll[NativeStep, NativeFail, NativeText, NativeSave];
NativeText[x_] := ToString[x, InputForm];
NativeStep[s_String] := (step=s;
  Print["NATIVE_EXT_C5_STAGE: ",s];
  Quiet@Check[Export[progress, <|"model"->key,"heavy_flavours"->nHeavy,
    "step"->s,"physical_dictionary_certified"->False|>,"RawJSON"],$Failed];s);
NativeSave[data_Association] := Module[{p = Join[<|
    "format"->"t3-extended-c5-native-matchete-one-step-v1", "model"->key,
    "heavy_flavours"->nHeavy,"last_step"->step,
    "absolute_matchete_rgbeta_dictionary_certified"->False,
    "hierarchical_matching_certified"->False,
    "all_allowed_scalar_quadratics_included"->False,
    "physical_neutrino_mass_prediction"->False|>, data]},
  Export[out,p,"RawJSON"];
  Print[key,": ",Lookup[p,"status","Unknown"]," -> ",out];p];
NativeFail[msg_String] := Throw[<|"status"->"NativeExtendedC5Blocked",
  "reason"->msg,"failed_step"->step|>,"T3NativeC5Stop"];
(* CRITICAL: Matchete must be loaded as its own COMPLETED top-level
   expression. Wolfram reads/parses the entire following Module BEFORE
   evaluating its body. If Needs["Matchete`"] is placed inside that Module,
   Bar[Coupling[...]] below is parsed as Global`Bar / Global`Coupling,
   producing Bar::shdw, Coupling::shdw and an aborted Matchete load.
   Do not move this package load back into the Module. *)
NativeStep["load-Matchete"];
nativeMatcheteLoad = CheckAbort[
  Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed],
  $Aborted
];
If[nativeMatcheteLoad =!= True,
  NativeSave[<|
    "status" -> "NativeExtendedC5Blocked",
    "reason" -> "Could not load native Matchete in its own top-level expression",
    "failed_step" -> step
  |>];
  Exit[6]
];
(* New top-level expression: now Matchete`Bar and Matchete`Coupling
   resolve in the installed package context, not a premature Global` context. *)
res = CheckAbort[Catch[Module[
  {model,built,matched,ext,coef,withZero,extraPiece,terms,
   interesting, fullPath,zeroPath,diffPath,zeroRules},
  Get[FileNameJoin[{dir,"model","T3ModelCatalog.wl"}]];
  Get[FileNameJoin[{dir,"model","LagrangianBuilder.wl"}]];
  Get[FileNameJoin[{dir,"RunMatching.wl"}]];
  NativeStep["build-full-general-fermion-quartic-Lagrangian"];
  model=T3ModelFromClass[cl,allowed[key],nHeavy];
  If[!AssociationQ[model],NativeFail["Unknown model"]];
  built=CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
  If[!AssociationQ[built] || Lookup[built,"Status",""]=!="Success" ||
      !TrueQ[Lookup[built,"FullValidation",False]] ||
      Lookup[built,"AdditionalYukawasIncluded",-1]=!=1 ||
      Lookup[built,"AdditionalYukawasIndividuallyValidated",-1]=!=1 ||
      !TrueQ[Lookup[built,"MatchingMassBasisReady",False]],
    NativeFail["Full Matchete build and one extra Yukawa did not validate"]];
  If[FreeQ[built["LInt"], If[allowed[key]===0,y1Conj,y2Conj]],
    NativeFail["Independent extra Yukawa not present in assembled LInt"]];
  NativeStep["native-one-loop-simultaneous-heavy-matching"];
  (* The repository's existing matching stages are used, including
     GreensSimplify/EOMSimplify/loop evaluation/effective coupling replacement.
     Do not replace this with a guessed loop formula. *)
  matched=CheckAbort[UsingFrontEnd[RunT3Matching[built["LUV"],5,1]],$Aborted];
  If[!AssociationQ[matched] || Lookup[matched,"Status",""]=!="Success" ||
      !KeyExistsQ[matched,"MatchedEFT"],
    NativeFail["Native Matchete matching stages did not all succeed"]];
  NativeStep["extract-native-symmetric-Weinberg-coefficient"];
  ext=CheckAbort[ExtractWeinbergCoefficient[matched["MatchedEFT"]],$Aborted];
  If[!AssociationQ[ext] || Lookup[ext,"Status",""]=!="Success" ||
      MissingQ[Lookup[ext,"CanonicalCoefficient",Missing["NoC5"]]],
    NativeFail["No certified extractable canonical Weinberg coefficient in full one-step EFT"]];
  coef=ext["CanonicalCoefficient"];
  If[!FreeQ[coef,_Missing] || coef===$Failed || coef===$Aborted,
    NativeFail["Native canonical coefficient unresolved"]];
  NativeStep["isolate-independent-extra-Yukawa-algebraic-sector"];
  (* Simplification of the full coefficient can be expensive. Match already
     canonicalizes; use Expand for the zero-extra algebraic counterfactual. *)
  zeroRules={HoldPattern[Bar[Coupling[y1Conj,___]]]:>0,
    HoldPattern[Bar[Coupling[y2Conj,___]]]:>0,
    HoldPattern[Coupling[y1Conj,___]]:>0,
    HoldPattern[Coupling[y2Conj,___]]:>0};
  withZero=Expand[coef/.zeroRules];
  If[!FreeQ[withZero,y1Conj] || !FreeQ[withZero,y2Conj],
    NativeFail["Zero-extra projection left a neutral-scalar Yukawa in C5"]];
  extraPiece=Expand[coef-withZero];
  interesting=If[allowed[key]===0,y1Conj,y2Conj];
  (* Do not classify a symbolic unevaluated equality as proven zero. *)
  terms=If[TrueQ[PossibleZeroQ[extraPiece]],"AlgebraicallyZero",
    If[FreeQ[extraPiece,interesting],"NoExtraSymbolAfterProjection",
       "ExplicitExtraYukawaDependence"]];
  fullPath=FileNameJoin[{DirectoryName[out],FileBaseName[out]<>".full_C5.inputform.txt"}];
  zeroPath=FileNameJoin[{DirectoryName[out],FileBaseName[out]<>".zero_extra_C5.inputform.txt"}];
  diffPath=FileNameJoin[{DirectoryName[out],FileBaseName[out]<>".extra_difference_C5.inputform.txt"}];
  If[Or @@ (Function[pair, Quiet@Check[
         Export[pair[[1]],NativeText[pair[[2]]],"Text"] === $Failed, True]] /@
      {{fullPath,coef},{zeroPath,withZero},{diffPath,extraPiece}}),
     NativeFail["Cannot export one of the native coefficient text files"]];
  NativeStep["complete"];
  <|"status"->"NativeOneStepExtendedWeinbergExtracted",
    "matching_kind"->"full-one-step-one-loop-EFTOrder5",
    "additional_yukawa"->SymbolName[interesting],
    "extraction_status"->ext["Status"],
    "native_coefficient_full_path"->fullPath,
    "native_coefficient_zero_extra_path"->zeroPath,
    "native_coefficient_extra_difference_path"->diffPath,
    "extra_dependent_sector_diagnostic"->terms,
    "native_full_coefficient_contains_extra_name"->(!FreeQ[coef,interesting]),
    "full_coefficient_flavor_index_structure_independently_verified"->False,
    "coefficient_not_a_direct_3x3_numeric_matrix"->True,
    "weinberg_operator_convention"->Lookup[ext,"OperatorConvention","Unknown"],
    "native_matchete_lagrangian_check_passed"->True,
    "quartic_extensions_included"->Lookup[built,"AdditionalQuarticsIncluded",0],
    "full_scalar_mass_mixing_included"->False,
    "one_loop_matching_performed"->True,
    "matched_c5_in_matchete_native_coordinates"->True,
    "extra_sector_zeroed_only_algebraically_after_matching"->True,
    "cross_package_numeric_substitution_authorized"->False|>
],"T3NativeC5Stop"],$Aborted];
If[!AssociationQ[res],res=<|"status"->"NativeExtendedC5Blocked",
  "reason"->"Matchete aborted or returned an unexpected result", "failed_step"->step|>];
NativeSave[res];
If[Lookup[res,"status",""]=!="NativeOneStepExtendedWeinbergExtracted",Exit[6]];
Exit[0];
