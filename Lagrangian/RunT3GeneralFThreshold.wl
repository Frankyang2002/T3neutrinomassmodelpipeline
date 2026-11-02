(* General non-mixed T3 F-first matching, preserving all symmetry-allowed
   Yukawa/quartic interactions. F is Heavy; scalars are Light. For six neutral
   self-bilinear cases use the existing *real-bilinear* mass-eigenfield
   convention, with both real scalar fields Light. Complex self-bilinear
   phases require a separately verified phase-dependent basis and are NOT
   certified here. For B/C(-1) use the independent mixed-doublet driver.

   This does not implement EFT1 Wilson RGEs, scalar-threshold matching, or
   physically certified C5. The Wilson expressions are exact symbolic
   Matchete tree outputs, with their matched SU(2) CG registry.

   wolframscript -file Lagrangian/RunT3GeneralFThreshold.wl OUT_DIR A_-4 0
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3, Print["Usage: OUT_DIR MODEL 0|1"]; Exit[2]];
outDir = ExpandFileName[args[[1]]];
modelName = args[[2]];
loopOrder = Which[args[[3]] === "0", 0, args[[3]] === "1", 1, True, -1];
validModels = <|"A"->{-4,-2,0}, "B"->{-3,1}, "C"->{-3,1},
  "D"->{-2,0,2}, "E"->{-4,-2,0,2}|>;
parsed = StringCases[modelName,
  RegularExpression["^([A-E])_([+-][0-9]+)$"] :> {"$1", "$2"}];
If[Length[parsed] =!= 1 || loopOrder < 0,
  Print["Invalid model/loop order"]; Exit[2]];
class = parsed[[1,1]];
alpha = ToExpression[parsed[[1,2]]];
If[!MemberQ[Lookup[validModels,class,{}],alpha],
  Print["Mixed models require the B/C eigen-doublet driver; invalid model: ",modelName]; Exit[2]];
If[!And@@(Environment[#] === "1" & /@ {
   "T3_EXPERIMENTAL_GENERAL_FERMIONS", "T3_EXPERIMENTAL_GENERAL_QUARTICS",
   "T3_EXPERIMENTAL_GENERAL_QUADRATICS", "T3_EXPERIMENTAL_FLAVOR_MASSES"}),
  Print["All general interactions and heavy flavour masses must be enabled"]; Exit[2]];
If[!DirectoryQ[outDir], CreateDirectory[outDir,CreateIntermediateDirectories->True]];
statusPath = FileNameJoin[{outDir,"native_status.json"}];
sourceDir = DirectoryName@ExpandFileName[$InputFileName];
ClearAll[writeStatus,stop];
writeStatus[rec_Association] := Export[statusPath,
  Join[<|"model"->modelName,"loop_order"->loopOrder,
   "operator_basis_is_full_general_requested"->True,
   "full_sequential_eft1_computed"->False,
   "physical_c5_certified"->False|>,rec],"RawJSON"];
stop[reason_String] := (writeStatus[<|"status"->"Failed", "error"->reason|>];
  Print["F-FIRST FAILED: ",reason]; Exit[5]);
Needs["Matchete`"];
Get[FileNameJoin[{sourceDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir,"model","LagrangianBuilder.wl"}]];
Get[FileNameJoin[{sourceDir,"RunMatching.wl"}]];
Print["F-FIRST model=",modelName," loopOrder=",loopOrder];
model = T3ModelFromClass[class,alpha,3];
If[!AssociationQ[model], stop["Could not construct three-heavy-flavour model"]];
If[SetT3HeavyFields[{"F"}] === $Failed, stop["Heavy field selection failed"]];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]],$Aborted];
If[!AssociationQ[built],stop["Full-general UV Lagrangian builder failed"]];
If[Lookup[built,"Status",None] =!= "Success" ||
   !TrueQ[Lookup[built,"FullValidation",False]] ||
   !TrueQ[Lookup[built,"FermionYukawaInventoryAssembled",False]] ||
   !TrueQ[Lookup[built,"AllScalarExtensionsRequested",False]] ||
   Lookup[built,"AdditionalQuadraticsIncluded",-1] =!= Lookup[built,"AdditionalQuadraticsExpected",-2] ||
   Lookup[built,"AdditionalQuarticsIncluded",-1] =!= Lookup[built,"AdditionalQuarticsExpected",-2],
  stop["Full-general Lagrangian is not validated and operator-complete in its declared inventory: "<>
    ToString[Lookup[built,"Status","Unknown"]]]];
(* Mass-basis preparation required by the six neutral self-bilinear models.
   Copy the real-component orientation from RunT3GeneralSelfMassBasis.wl:
      S = (S_minus + I S_plus)/Sqrt[2]
      Mminus^2 = MS^2 - 2 Re[kappa], Mplus^2 = MS^2 + 2 Re[kappa].
   These are independent mass symbols in the symbolic matching and MUST be
   substituted consistently at a numerical endpoint. All other interactions
   including every extra Yukawa/quartic remain present after rotation. *)
selfModels = {"A_-2", "A_+0", "D_-2", "D_+0", "E_-2", "E_+0"};
selfQ = MemberQ[selfModels, modelName];
lagForThreshold = built["LUV"];
selfSpecies = None;
selfMetadata = <|"mass_basis_kind"->"OriginalComplexScalars",
  "complex_self_bilinear_phase_supported"->False|>;
If[selfQ,
  If[Lookup[built,"AdditionalQuadraticsIncluded",-1] =!= 1 ||
      Lookup[built,"AdditionalYukawasIncluded",-1] =!= 1,
    stop["Self-bilinear case missing one quadratic or extra Yukawa"]];
  selfSpecies = If[alpha === -2, 2, 1];
  originalField = If[selfSpecies === 2, NewScalar2, NewScalar1];
  plusField = If[selfSpecies === 2, T3Scalar2Plus, T3Scalar1Plus];
  minusField = If[selfSpecies === 2, T3Scalar2Minus, T3Scalar1Minus];
  plusMass = If[selfSpecies === 2, T32PlusMass, T31PlusMass];
  minusMass = If[selfSpecies === 2, T32MinusMass, T31MinusMass];
  repDim = model[If[selfSpecies === 2, "Scalar2", "Scalar1"], "SU2"];
  repY = model[If[selfSpecies === 2, "Scalar2", "Scalar1"], "Y"];
  If[!MemberQ[{1,3},repDim] || repY =!= 0,
    stop["Self-bilinear mass basis requires a neutral real SU2 irrep"]];

  (* Isolate only the holomorphic quadratic, never all terms involving S. *)
  hasExtraMassQ[term_] := !FreeQ[term,
    HoldPattern[Coupling[s_Symbol, ___]] /;
      StringStartsQ[SymbolName[s], "lambdaT3ExtraM"]];
  expandedInt = Expand[built["LInt"]];
  allIntTerms = If[Head[expandedInt] === Plus,
    List @@ expandedInt, {expandedInt}];
  extraMassTerms = Select[allIntTerms, hasExtraMassQ];
  If[Length[extraMassTerms] < 2,
    stop["Could not isolate both self-bilinear and Hermitian conjugate"]];
  remainingInt = Expand[Total[Select[allIntTerms, !hasExtraMassQ[#] &]]];
  If[hasExtraMassQ[remainingInt],
    stop["Quadratic spurion remains in interaction Lagrangian"]];
  If[Length[Cases[remainingInt,
      HoldPattern[Field[x_Symbol, Scalar, _, _]] /; x === originalField,
      Infinity]] == 0,
    stop["No evaluated original scalar Field objects in retained interactions"]];

  indexType = BSMIndexType[repDim];
  indices = If[indexType === None, {}, {indexType}];
  If[selfSpecies === 2, RemoveField[NewScalar2], RemoveField[NewScalar1]];
  (* CRITICAL: Sminus and Splus remain Light at the F threshold. *)
  If[selfSpecies === 1,
    DefineField[T3Scalar1Minus,Scalar,SelfConjugate->True,
      Indices->indices,Mass->{Light,T31MinusMass}];
    DefineField[T3Scalar1Plus,Scalar,SelfConjugate->True,
      Indices->indices,Mass->{Light,T31PlusMass}],
    DefineField[T3Scalar2Minus,Scalar,SelfConjugate->True,
      Indices->indices,Mass->{Light,T32MinusMass}];
    DefineField[T3Scalar2Plus,Scalar,SelfConjugate->True,
      Indices->indices,Mass->{Light,T32PlusMass}]
  ];
  (* Bar rules MUST precede Field rules: the imaginary component changes sign. *)
  rotationRules = If[selfSpecies === 1,
    {
      HoldPattern[Bar[Field[NewScalar1,Scalar,ix_,opts_]]] :>
       (Bar[Field[T3Scalar1Minus,Scalar,ix,opts]] -
        I Bar[Field[T3Scalar1Plus,Scalar,ix,opts]])/Sqrt[2],
      HoldPattern[Field[NewScalar1,Scalar,ix_,opts_]] :>
       (Field[T3Scalar1Minus,Scalar,ix,opts] +
        I Field[T3Scalar1Plus,Scalar,ix,opts])/Sqrt[2]
    },
    {
      HoldPattern[Bar[Field[NewScalar2,Scalar,ix_,opts_]]] :>
       (Bar[Field[T3Scalar2Minus,Scalar,ix,opts]] -
        I Bar[Field[T3Scalar2Plus,Scalar,ix,opts]])/Sqrt[2],
      HoldPattern[Field[NewScalar2,Scalar,ix_,opts_]] :>
       (Field[T3Scalar2Minus,Scalar,ix,opts] +
        I Field[T3Scalar2Plus,Scalar,ix,opts])/Sqrt[2]
    }
  ];
  rotated = RelabelIndices[Expand[remainingInt /. rotationRules]];
  If[!FreeQ[rotated,originalField] || hasExtraMassQ[rotated] ||
      FreeQ[rotated,plusField] || FreeQ[rotated,minusField],
    stop["Old scalar, quadratic spurion, or missing mass eigenfield after rotation"]];
  freeEigen = If[selfSpecies === 1,
    FreeLag[NewFermion, NewScalar2, T3Scalar1Minus, T3Scalar1Plus],
    FreeLag[NewFermion, NewScalar1, T3Scalar2Minus, T3Scalar2Plus]
  ];
  lagForThreshold = RelabelIndices[Expand[built["LSM"] + freeEigen + rotated]];
  validEigen = CheckAbort[UsingFrontEnd[CheckLagrangian[lagForThreshold]],$Aborted];
  If[validEigen =!= True,
    Export[FileNameJoin[{outDir,"self_mass_basis_check_failed.inputform.txt"}],
      ToString[InputForm[lagForThreshold]],"Text"];
    stop["Real self-bilinear mass-eigenfield Lagrangian failed CheckLagrangian"]];
  selfMetadata = <|
    "mass_basis_kind"->"RealSelfBilinearSplitLightFields",
    "split_species"->selfSpecies,
    "original_scalar"->ToString[originalField,InputForm],
    "minus_field"->ToString[minusField,InputForm],
    "plus_field"->ToString[plusField,InputForm],
    "minus_mass_symbol"->ToString[minusMass,InputForm],
    "plus_mass_symbol"->ToString[plusMass,InputForm],
    "mass_relation"->"Mminus^2=MS^2-2 Re(kappa); Mplus^2=MS^2+2 Re(kappa)",
    "scalar_mass_components_remain_light"->True,
    "real_bilinear_mass_basis_candidate"->True,
    "complex_self_bilinear_phase_supported"->False,
    "holomorphic_bilinear_absorbed_into_split_free_masses"->True,
    "field_basis_normalization_physically_certified"->False
  |>;
  Export[FileNameJoin[{outDir,"self_mass_basis.json"}],selfMetadata,"RawJSON"];
  Print["F-FIRST: diagonal neutral self-bilinear; Light fields=",
    minusField,", ",plusField,"; original full interactions rotated"];
];

WriteString[$Output,"F-FIRST: model built; calling actual Matchete F-only threshold\n"];
result = CheckAbort[UsingFrontEnd[RunSequentialT3Matching[
  lagForThreshold,model,{{"F"}},5,loopOrder,outDir,
  StringReplace[modelName,{"+"->"p","-"->"m"}]]],$Aborted];
If[!AssociationQ[result] || Lookup[result,"Status",""] =!= "Success" ||
   Length[Lookup[result,"Stages",{}]] =!= 1,
  stop["Matchete F-only matching failed (inspect wolframscript log)"]];
stage = First[result["Stages"]];
If[Lookup[stage,"IntegratedFields",{}] =!= {"F"},
  stop["Unexpected integrated particle inventory"]];
If[!FreeQ[stage["TreeLagrangian"],Field[NewFermion,___]],
  stop["Integrated heavy F remains in tree EFT"]];
(* The pre-existing RunMatching.wl exporter counts only NewScalar1/2.
   A self-bilinear split renames one of those to TWO real, Light scalars;
   otherwise its old Wilson selector silently loses the corresponding terms.
   Re-export the matching-generated full LLSS operator inventory with the
   split field names and exact CG tensors, without modifying RunMatching.wl. *)
If[selfQ,
  originalSeedFiles = FileNames["*wilson_seed.json",
    FileNameJoin[{outDir,"data"}]];
  If[Length[originalSeedFiles] =!= 1,
    stop["Expected one production Wilson seed before split-field recataloguing"]];
  treeTerms = Module[{t = Expand[stage["TreeLagrangian"]]},
    If[Head[t] === Plus, List @@ t, {t}]];
  scalarLabels = {"Field[NewScalar1", "Field[NewScalar2",
    "Field[T3Scalar1Minus", "Field[T3Scalar1Plus",
    "Field[T3Scalar2Minus", "Field[T3Scalar2Plus"};
  splitTerms = Select[treeTerms, Function[term,
    Module[{txt = ToString[InputForm[term]], cnt},
      cnt = Total[StringCount[txt, #] & /@ scalarLabels];
      StringContainsQ[txt, "GammaCC"] && cnt === 2
    ]]];
  If[splitTerms === {},
    stop["No GammaCC with two light scalar factors found in split-field tree EFT"]];
  (* Require a matching-generated operator containing the rotated scalar. *)
  splitLabel = If[selfSpecies === 1,"Field[T3Scalar1","Field[T3Scalar2"];
  If[!AnyTrue[splitTerms,
    StringContainsQ[ToString[InputForm[#]],splitLabel] &],
    stop["No tree Wilson term contains the split scalar mass eigenfields"]];
  relevantCG = DeleteDuplicates @ Cases[HoldComplete @@ splitTerms,
    CG[name_,___] :> Replace[name,Bar[x_] :> x],Infinity];
  splitCGRegistry = Map[Function[name, Module[{reps,tensor},
    reps = Quiet @ Check[Matchete`PackageScope`$CGproperties[name,Indices],$Failed];
    tensor = Quiet @ Check[Matchete`CGManipulations`PackagePrivate`$CGtensors[name],$Failed];
    If[MemberQ[{$Failed},reps] || MemberQ[{$Failed},tensor] ||
       MissingQ[reps] || MissingQ[tensor],
      Nothing,
      <|"Name"->ToString[Unevaluated[name],InputForm],
        "RepsInputForm"->ToString[InputForm[reps]],
        "TensorInputForm"->ToString[InputForm[tensor]]|>
    ]]],relevantCG];
  If[Length[splitCGRegistry] =!= Length[relevantCG],
    stop["Incomplete SU2 CG registry for split-scalar tree Wilson sector"]];
  originalSeed = Import[First[originalSeedFiles],"RawJSON"];
  If[!AssociationQ[originalSeed],stop["Original Wilson seed malformed"]];
  (* Keep the unmodified production export for provenance. This filename
     does NOT match the validator's *wilson_seed.json glob. *)
  Export[FileNameJoin[{outDir,"data","unsplit_seed_archive.json"}],
    originalSeed,"RawJSON"];
  splitSeed = Join[originalSeed,<|
    "TreeWilsonTermCount"->Length[splitTerms],
    "TreeWilsonTerms"->MapIndexed[Function[{term,pos},
      Module[{txt = ToString[InputForm[term]]},<|
        "Index"->First[pos],
        "Chirality"->Which[StringContainsQ[txt,"Proj[-1]"],"PL",
          StringContainsQ[txt,"Proj[1]"],"PR",True,"Unknown"],
        "Scalar1Count"->Total[StringCount[txt,#] & /@
          {"Field[NewScalar1","Field[T3Scalar1Minus","Field[T3Scalar1Plus"}],
        "Scalar2Count"->Total[StringCount[txt,#] & /@
          {"Field[NewScalar2","Field[T3Scalar2Minus","Field[T3Scalar2Plus"}],
        "CGNames"->(ToString[Unevaluated[#],InputForm]& /@
          DeleteDuplicates@Cases[HoldComplete[term],
          CG[name_,___]:>Replace[name,Bar[x_]:>x],Infinity]),
        "TermInputForm"->txt|>]],splitTerms],
    "CGRegistry"->splitCGRegistry,
    "SelfBilinearMassBasis"->selfMetadata,
    "WilsonSeedSource"->"MatcheteFThresholdTreeReclassifiedForSplitLightScalarFields"
  |>];
  Export[First[originalSeedFiles],splitSeed,"RawJSON"];
  Print["F-FIRST: complete split-scalar Wilson terms=",Length[splitTerms]];
];

seedFiles = FileNames["*wilson_seed.json",FileNameJoin[{outDir,"data"}]];
If[Length[seedFiles] =!= 1,stop["Expected exactly one exported EFT1 Wilson seed JSON"]];
seed = Import[First[seedFiles],"RawJSON"];
If[!AssociationQ[seed] || Lookup[seed,"IntegratedFields",{}] =!= {"F"} ||
   !IntegerQ[Lookup[seed,"TreeWilsonTermCount",None]] ||
   Lookup[seed,"TreeWilsonTermCount",0] <= 0 ||
   Length[Lookup[seed,"TreeWilsonTerms",{}]] =!= Lookup[seed,"TreeWilsonTermCount",0],
   stop["F-only Wilson seed export missing or inconsistent"]];
If[AnyTrue[seed["TreeWilsonTerms"],
    StringContainsQ[Lookup[#,"TermInputForm",""],"Field[NewFermion"] &],
  stop["Wilson seed contains integrated fermion"]];
writeStatus[<|
  "status"->"FThresholdTreeWilsonExported",
  "integrated_fields"->{"F"}, "light_bsm_fields"->{"S1","S2"},
  "heavy_flavour_count"->3,
  "tree_wilson_term_count"->seed["TreeWilsonTermCount"],
  "wilson_seed_file"->FileNameTake[First[seedFiles]],
  "full_general_extra_yukawa_count"->built["AdditionalYukawasIncluded"],
  "full_general_quartic_count"->built["AdditionalQuarticsIncluded"],
  "full_general_quadratic_count"->built["AdditionalQuadraticsIncluded"],
  "self_bilinear_mass_basis"->selfMetadata,
  "raw_one_loop_F_matching_included"->(loopOrder==1),
  "hard_one_loop_MSbar_certified"->False,
  "eft1_wilson_anomalous_dimensions_implemented"->False,
  "scalar_threshold_matching_implemented"->False|>];
Print["F-FIRST OK: ",modelName," tree Wilson terms=",seed["TreeWilsonTermCount"]];
Exit[0];
