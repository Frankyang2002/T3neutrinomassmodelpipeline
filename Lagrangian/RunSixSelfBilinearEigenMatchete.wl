(* Workstream A: EXPERIMENTAL mass-eigenfield UV audit for six T3 self-bilinear models.

   wolframscript -file Lagrangian/RunSixSelfBilinearEigenMatchete.wl \
       OUTPUT.json A|D|E m2|p0

   Builds the fully extended original UV Lagrangian first, removes precisely its
   one extra holomorphic scalar bilinear, rewrites EVERY other interaction in
   two real self-conjugate scalar multiplets, and rebuilds the free Lagrangian.
   The complex rotation phase is a symbolic unit-modulus spurion: its physical
   value and the heavy eigenmass relations are NOT imposed by Matchete.

   An ordinary/detailed CheckLagrangian pass is a UV consistency check, NOT a
   proof of component normalization, threshold matching, or completeness.
   No existing repository file is modified and Match[] is never called.
*)
ClearAll["Global`*"];
$T3SBArgs = Rest[$ScriptCommandLine];
If[Length[$T3SBArgs] =!= 3,
  Print["Usage: OUTPUT.json A|D|E m2|p0"]; Exit[2]
];
$T3SBOutputPath = ExpandFileName[$T3SBArgs[[1]]];
$T3SBClass = ToUpperCase[$T3SBArgs[[2]]];
$T3SBAlpha = Switch[$T3SBArgs[[3]], "m2",-2,"p0",0,_, $Failed];
If[!MemberQ[{"A","D","E"},$T3SBClass] || $T3SBAlpha === $Failed,
  Print["Only A/D/E with m2 or p0 are supported"]; Exit[2]
];
$T3SBModelName = $T3SBClass <> "_" <>
  If[$T3SBAlpha === -2, "-2", "+0"];
$T3SBDir = DirectoryName @ ExpandFileName[$InputFileName];
$T3SBRequiredFlags = {
  "T3_EXPERIMENTAL_FLAVOR_MASSES",
  "T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS"
};
If[!AllTrue[$T3SBRequiredFlags, Environment[#] === "1" &],
  Print["All four T3 opt-in flags must be 1. Missing: ",
    Select[$T3SBRequiredFlags, Environment[#] =!= "1" &]];
  Exit[2]
];
(* CRITICAL: Matchete must be loaded BEFORE the parser encounters any
   unqualified Matchete public symbol (Field, Coupling, CheckLagrangian,
   DefineField, ...).  Defining the T3SB functions before Needs creates
   Global`Field etc and shadows the package's own API symbols.  The working
   B/C mass-eigenstate driver likewise loads the package first.  This load
   is a standalone top-level expression; don't combine it with the later
   API-symbol preflight in the same held Module/If expression. *)
Print["T3SB_STAGE: load-Matchete-before-API-symbol-parsing"];
Block[{$MessageList = {}},
  $T3SBNeedsResult = CheckAbort[
    UsingFrontEnd[Needs["Matchete`"] ; True], $Aborted
  ];
  $T3SBNeedsMessages = ToString[$MessageList, InputForm];
];

$T3SBStage = "initializing";
$T3SBStages = {};
T3SBStage[s_String] := (
  $T3SBStage = s; AppendTo[$T3SBStages, s];
  Print["T3SB_STAGE: ",s]; s
);
T3SBExport[a_Association] := Module[{result = a},
  If[!DirectoryQ[DirectoryName[$T3SBOutputPath]],
    CreateDirectory[DirectoryName[$T3SBOutputPath],
      CreateIntermediateDirectories -> True]];
  Export[$T3SBOutputPath, Join[<|
    "format"->"t3-six-self-eigen-matchete-v1",
    "model"->$T3SBModelName,
    "last_stage"->$T3SBStage,
    "stage_history"->$T3SBStages,
    "physical_matching_performed"->False,
    "matchete_rgbeta_normalization_certified"->False,
    "native_component_identification_certified"->False,
    "mass_parameter_relation_imposed_in_matchete"->False
  |>,result], "RawJSON"]
];
T3SBFail[code_String, detail_: ""] := (
  Print["T3SB_FAILURE: ",code," ",detail];
  T3SBExport[<|"status"->code,"diagnostic"->detail|>];
  Throw[$Failed,"T3SBStop"]
);

T3SBHasLabelQ[expr_, label_Symbol] := With[{s=label},
  !FreeQ[expr, Verbatim[s], {0,Infinity}]
];
T3SBFieldLabels[expr_] := DeleteDuplicates @ Cases[expr,
  HoldPattern[Field[label_Symbol, ___]] :> label,{0,Infinity}];
T3SBCouplingLabels[expr_] := DeleteDuplicates @ Cases[expr,
  HoldPattern[Coupling[label_Symbol, ___]] :> label,{0,Infinity}];
T3SBMassLists[expr_] := DeleteDuplicates @ Cases[expr,
  HoldPattern[Coupling[MF, ids_List, ___]] :> ids,{0,Infinity}];
T3SBCheckFlavorMass[expr_] := Module[{lists=T3SBMassLists[expr]},
  lists =!= {} && AllTrue[lists, MatchQ[#, {Index[_,NFlavor]}] &] &&
    T3SBHasLabelQ[expr, NewFermion]
];

(* Do not use a public shortcut for a source typed Index: it would wrap twice.
   The field's SU(2) index is preserved verbatim in its replacement. *)
T3SBRealComponent[indices_List] := Field[T3SBReal, Scalar, indices, {}];
T3SBImagComponent[indices_List] := Field[T3SBImag, Scalar, indices, {}];
T3SBRotate[expr_, source_Symbol] := With[{label = source},
  Module[{barred, result},
    barred = expr /. {
      HoldPattern[Bar[Field[label, Scalar, inds_List, rest___]]] :>
        (Bar[T3SBPhase[]] (T3SBRealComponent[inds]-I T3SBImagComponent[inds])/Sqrt[2])
    };
    result = barred /. {
      HoldPattern[Field[label, Scalar, inds_List, rest___]] :>
        (T3SBPhase[] (T3SBRealComponent[inds]+I T3SBImagComponent[inds])/Sqrt[2])
    };
    If[T3SBHasLabelQ[result,label],Return[$Failed]];
    result // Expand // RelabelIndices
  ]
];

T3SBPerform[] := Module[
  {model,built,source,other,dim,neutral,structures,counts,names,data,tensor,
   factor,residual,key,extraMass,allSource,terms,removed,kept,rotated,
   sourceLabels,rotatedLabels,missingLabels,oldFields,free,lag,
   registerReal,registerImag,registerPhase,check,detailed,heavy,countsOK,
   nQuartics,expectedQuartics,nExtraYukawa,allFields,newFields,mfLists,
   report},
  T3SBStage["build-original-full-UV"];
  model = T3ModelFromClass[$T3SBClass,$T3SBAlpha];
  If[model === $Failed,T3SBFail["InvalidModel"]];
  built = BuildT3Lagrangian[model];
  If[!AssociationQ[built] || built["Status"] =!= "Success" ||
     !TrueQ[built["FullValidation"]],
    T3SBFail["OriginalFullLagrangianFailed",ToString[
      If[AssociationQ[built],Lookup[built,"Status","Unknown"],built],InputForm]]
  ];
  nQuartics = Lookup[built,"AdditionalQuarticsIncluded",Missing["Unreported"]];
  expectedQuartics = If[$T3SBClass === "E",6,5];
  nExtraYukawa = Lookup[built,"AdditionalYukawasIncluded",Missing["Unreported"]];
  countsOK = nQuartics === expectedQuartics && nExtraYukawa === 1 &&
    Lookup[built,"AdditionalQuadraticsIncluded",-1] === 1;
  If[!countsOK,T3SBFail["FullOperatorInventoryMismatch",
     ToString[{nQuartics,expectedQuartics,nExtraYukawa,
       Lookup[built,"AdditionalQuadraticsIncluded",Missing["Unknown"]]},InputForm]]];

  neutral = If[$T3SBAlpha === -2, NewScalar2, NewScalar1];
  other = If[$T3SBAlpha === -2, NewScalar1, NewScalar2];
  dim = model[If[$T3SBAlpha === -2,"Scalar2","Scalar1"],"SU2"];
  If[!MemberQ[{1,3},dim], T3SBFail["NotRealSU2Irrep"]];
  T3SBStage["extract-native-invariant"];
  structures = T3ExtendedExtraQuadraticStructures[model];
  If[Length[structures] =!= 1,T3SBFail["QuadraticStructureCountMismatch"]];
  counts = First[structures];
  names = T3ExtendedSpeciesList[counts];
  If[names =!= ConstantArray[If[$T3SBAlpha === -2,"S2dag","S1dag"],2],
     T3SBFail["WrongQuadraticSpecies",ToString[names,InputForm]]];
  data = T3ExtendedBosonicTensorBasis[model,names];
  If[!AssociationQ[data] || data["PhysicalCount"] =!= 1,
     T3SBFail["NativeInvariantExtractionFailed"]];
  tensor=Normal[First[data["Basis"]]];
  If[dim === 1,
    factor = If[Head[tensor] === List, First[Flatten[tensor]],tensor];
    residual = 0,
    If[Dimensions[tensor] =!= {3,3},T3SBFail["UnexpectedNativeTensorShape"]];
    factor = FullSimplify[Tr[tensor]/3];
    residual = FullSimplify[tensor-factor IdentityMatrix[3]]
  ];
  If[TrueQ[PossibleZeroQ[factor]] ||
     !And@@(TrueQ[PossibleZeroQ[#]] & /@ Flatten[{residual}]),
    T3SBFail["NativeTensorNotCartesianIdentity",
      ToString[tensor,InputForm]]
  ];
  key = StringJoin[ToString /@ counts];
  extraMass = Symbol["lambdaT3ExtraM"<>key<>"I1"];
  allSource = Expand[built["LInt"]] // RelabelIndices;
  If[!T3SBHasLabelQ[allSource,extraMass],
    T3SBFail["ExtraMassCouplingMissingFromSource",SymbolName[extraMass]]];

  (* Remove only the holomorphic bilinear, not other interactions.  Every
     monomial carries a coupling; conjugate terms carry the SAME label atom. *)
  T3SBStage["remove-only-holomorphic-bilinear"];
  terms = If[Head[allSource] === Plus,List@@allSource,{allSource}];
  removed = Select[terms,T3SBHasLabelQ[#,extraMass]&];
  kept = Select[terms,!T3SBHasLabelQ[#,extraMass]&];
  If[removed === {} || kept === {} ||
     !TrueQ[Expand[Total[removed]+Total[kept]-allSource]===0] ||
     T3SBHasLabelQ[Total[kept],extraMass],
    T3SBFail["QuadraticOperatorSeparationFailed"]];
  sourceLabels = T3SBCouplingLabels[Total[kept]];
  If[Length[sourceLabels] < 12,
    T3SBFail["UnexpectedSmallInteractionInventory",
      ToString[SymbolName /@ sourceLabels,InputForm]]];

  T3SBStage["deregister-source-scalar"];
  If[Quiet @ Check[With[{label=neutral},RemoveField[label]];True,$Failed] =!= True,
    T3SBFail["NeutralScalarDeregistrationFailed"]];
  T3SBStage["register-two-real-eigenmultiplets"];
  registerPhase = Quiet @ Check[
    DefineCoupling[T3SBPhase,SelfConjugate->False];True,$Failed];
  If[registerPhase =!= True,T3SBFail["RotationPhaseRegistrationFailed"]];
  If[dim === 1,
    registerReal = Quiet @ Check[
      DefineField[T3SBReal,Scalar,SelfConjugate->True,
        Mass->{Heavy,MT3SBReal}];True,$Failed];
    registerImag = Quiet @ Check[
      DefineField[T3SBImag,Scalar,SelfConjugate->True,
        Mass->{Heavy,MT3SBImag}];True,$Failed],
    registerReal = Quiet @ Check[
      DefineField[T3SBReal,Scalar,SelfConjugate->True,
        Indices->{T3BSMd3},Mass->{Heavy,MT3SBReal}];True,$Failed];
    registerImag = Quiet @ Check[
      DefineField[T3SBImag,Scalar,SelfConjugate->True,
        Indices->{T3BSMd3},Mass->{Heavy,MT3SBImag}];True,$Failed]
  ];
  If[registerReal =!= True || registerImag =!= True,
    T3SBFail["RealEigenfieldRegistrationFailed"]];

  T3SBStage["rotate-complete-interaction-inventory"];
  rotated=T3SBRotate[Total[kept],neutral];
  If[rotated === $Failed || T3SBHasLabelQ[rotated,neutral],
    T3SBFail["ComplexNeutralScalarNotEliminated"]];
  rotatedLabels = T3SBCouplingLabels[rotated];
  missingLabels = Complement[sourceLabels,rotatedLabels];
  If[missingLabels =!= {} ||
     !T3SBHasLabelQ[rotated,T3SBPhase] ||
     T3SBHasLabelQ[rotated,extraMass],
    T3SBFail["InteractionCouplingsLostOrExtraMassRetained",
      ToString[SymbolName /@ missingLabels,InputForm]]];
  newFields = T3SBFieldLabels[rotated];
  If[!MemberQ[newFields,T3SBReal] || !MemberQ[newFields,T3SBImag],
    T3SBFail["BothEigenmultipletsNotInInteractions"]];

  (* LFree is REBUILT, not rotated: native FreeLag supplies canonical kinetic
     terms for the real multiplets and diagonal scalar mass terms. *)
  T3SBStage["rebuild-diagonal-free-Lagrangian"];
  free=Quiet @ Check[FreeLag[NewFermion,other,T3SBReal,T3SBImag],$Failed];
  If[free === $Failed,T3SBFail["DiagonalFreeLagFailed"]];
  mfLists = T3SBMassLists[free];
  If[!T3SBCheckFlavorMass[free],
     T3SBFail["IndexedFermionMassMissing",
       ToString[mfLists,InputForm]]];
  lag=(built["LSM"]+free+rotated) // Expand // RelabelIndices;
  If[T3SBHasLabelQ[lag,neutral] ||
     !T3SBHasLabelQ[lag,other] ||
     !T3SBHasLabelQ[lag,T3SBReal] ||
     !T3SBHasLabelQ[lag,T3SBImag] ||
     T3SBHasLabelQ[lag,extraMass],
    T3SBFail["FinalEigenbasisFieldContentMismatch"]];

  T3SBStage["matchete-ordinary-check"];
  check=CheckAbort[CheckLagrangian[lag],$Aborted];
  T3SBStage["matchete-detailed-check"];
  detailed=CheckAbort[CheckLagrangian[lag,DetailedOutput->True],$Aborted];
  heavy=AssociationQ[detailed] &&
    TrueQ[Lookup[detailed,"HeavyMassBasis",False]];
  report=<|
    "status"->If[TrueQ[check] && heavy,
      "MatcheteSelfEigenHeavyBasisPassedMatchingNotPerformed",
      "SelfEigenMatcheteCheckFailed"],
    "native_dimension"->dim,
    "neutral_original_field"->SymbolName[neutral],
    "other_original_field"->SymbolName[other],
    "native_bilinear_species"->names,
    "native_bilinear_coupling"->SymbolName[extraMass],
    "native_bilinear_metric_factor_inputform"->ToString[factor,InputForm],
    "native_bilinear_metric_proportional_identity"->True,
    "neutral_scalar_removed_from_final_expression"->True,
    "quadratic_source_term_count_removed"->Length[removed],
    "interaction_coupling_labels_preserved"->True,
    "interaction_coupling_count_before"->Length[sourceLabels],
    "interaction_coupling_count_after"->Length[rotatedLabels],
    "additional_quartics_retained"->nQuartics,
    "additional_yukawas_retained"->nExtraYukawa,
    "final_indexed_fermion_mass_present"->True,
    "mf_mass_index_lists_inputform"->(ToString[#,InputForm]& /@ mfLists),
    "mass_eigenfield_labels"->{"T3SBReal","T3SBImag"},
    "eigenfield_mass_couplings"->{"MT3SBReal","MT3SBImag"},
    "complex_rotation_phase_spurion"->"T3SBPhase",
    "phase_modulus_one_assumed_not_encoded"->True,
    "ordinary_matchete_check_passed"->TrueQ[check],
    "detailed_heavy_mass_basis_passed"->heavy,
    "detailed_check_inputform"->If[TrueQ[check] && heavy,
      "Passed", ToString[detailed,InputForm]],
    "original_full_matchete_check_passed"->True,
    "no_physical_eigenmass_relation_imposed"->True,
    "independent_operator_normalizations_certified"->False
  |>;
  report
];

(* The package was loaded in its own top-level expression BEFORE any of
   these Matchete API heads were parsed. Never rerun Needs after defining the
   audit helper symbols; a second call would obscure startup diagnostics.
   Do not require $FrontEnd to remain non-Null OUTSIDE UsingFrontEnd. *)
T3SBLoadPreflight[] := Module[
  {probe, probeMessages, functions, ready, result},
  T3SBStage["Matchete-preflight"];
  Print["T3SB_MATCHECK kernel version: ", $Version];
  Print["T3SB_MATCHECK system ID: ", $SystemID];
  Print["T3SB_MATCHECK package lookup: ",
    ToString[Quiet @ Check[FindFile["Matchete`"], $Failed], InputForm]];
  (* Exact context: Names["*`DefineField"] incorrectly accepts the
     Global` symbol that caused the original namespace collision. *)
  functions = Association @ Table[
    name -> (Names["Matchete`" <> name] =!= {}),
    {name, {"CheckLagrangian", "LoadModel", "FreeLag",
            "DefineField", "GetRepresentations"}}
  ];
  Block[{$MessageList = {}},
    probe = CheckAbort[
      Quiet @ Check[UsingFrontEnd[GetRepresentations[]], $Failed],
      $Aborted
    ];
    probeMessages = ToString[$MessageList, InputForm]
  ];
  ready = $T3SBNeedsResult === True &&
    And @@ Values[functions] && AssociationQ[probe];
  result = <|
    "ready" -> TrueQ[ready],
    "kernel_version" -> ToString[$Version, InputForm],
    "system_id" -> ToString[$SystemID, InputForm],
    "matchete_findfile_inputform" -> ToString[
      Quiet @ Check[FindFile["Matchete`"], $Failed], InputForm],
    "front_end_outside_usingfront_end_inputform" ->
      ToString[$FrontEnd, InputForm],
    "needs_result_inputform" -> ToString[$T3SBNeedsResult, InputForm],
    "needs_message_list_inputform" -> $T3SBNeedsMessages,
    "public_api_functions_found_in_matchete_context" -> functions,
    "representations_probe_head" -> ToString[Head[probe], InputForm],
    "representations_probe_message_list_inputform" -> probeMessages
  |>;
  Print["T3SB_MATCHECK READY: ", ready,
        "; Needs result=", result["needs_result_inputform"],
        "; Matchete-context APIs=", ToString[functions, InputForm],
        "; representations head=", result["representations_probe_head"]];
  If[!ready,
    Print["T3SB_MATCHECK diagnostic: ", ToString[result, InputForm]]];
  result
];
$T3SBMatcheteLoadDiagnostic=T3SBLoadPreflight[];
If[!TrueQ[$T3SBMatcheteLoadDiagnostic["ready"]],
  T3SBStage["Matchete-not-ready"];
  T3SBExport[<|"status"->"MatcheteUnavailable",
    "matchete_load_diagnostic"->$T3SBMatcheteLoadDiagnostic|>];
  Exit[3]
];
If[Environment["T3SB_PREFLIGHT_ONLY"] === "1",
  T3SBStage["preflight-only-completed"];
  T3SBExport[<|"status"->"MatcheteReadyPreflight",
    "matchete_load_diagnostic"->$T3SBMatcheteLoadDiagnostic|>];
  Print["T3SB: Matchete preflight passed; no UV physics was run."];
  Exit[0]
];
Get[FileNameJoin[{$T3SBDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{$T3SBDir,"model","LagrangianBuilder.wl"}]];
T3SBStage["run-inside-front-end"];
$T3SBRun=CheckAbort[UsingFrontEnd[Catch[T3SBPerform[],"T3SBStop"]],$Aborted];
If[!AssociationQ[$T3SBRun],
  If[!FileExistsQ[$T3SBOutputPath],
    T3SBExport[<|"status"->"UncaughtOrAbortedError",
      "kernel_result_inputform"->ToString[$T3SBRun,InputForm],
      "message_list_inputform"->ToString[$MessageList,InputForm]|>]];
  Exit[5]
];
T3SBStage["finished"];
T3SBExport[Join[$T3SBRun,<|"matchete_load_diagnostic"->$T3SBMatcheteLoadDiagnostic|>]];
Print[$T3SBModelName,": ",Lookup[$T3SBRun,"status","Unknown"]];
If[Lookup[$T3SBRun,"status",""] =!=
   "MatcheteSelfEigenHeavyBasisPassedMatchingNotPerformed",Exit[6]];
Exit[0];
