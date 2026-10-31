(* Workstream A: ISOLATED RGBeta registration for the six neutral-scalar extra Yukawas.
   Native RGBeta v1.2.0 AddYukawa CheckInvariance=True produced unresolved
   symbolic Ward residuals for ALL gauge groups on A_m2, including SU3c
   for colourless fields. This driver deliberately does NOT call or claim
   the native Ward check. It verifies registration/field/invariant metadata
   and exact gauge representation eligibility, not a cross-package CG map.

   wolframscript -file Lagrangian/RunT3RGBetaExtraYukawaBridge.wl OUTPUT.json A m2
   Production RGBeta code and beta functions are never modified/evaluated.
*)
ClearAll["Global`*"];
$T3RYArgs = Rest[$ScriptCommandLine];
If[Length[$T3RYArgs] =!= 3,
  Print["Usage: OUTPUT.json A|D|E m2|p0"]; Exit[2]];
$T3RYOut = ExpandFileName[$T3RYArgs[[1]]];
$T3RYLabel = ToUpperCase[$T3RYArgs[[2]]] <> "_" <> $T3RYArgs[[3]];
$T3RYSpecs = <|
  "A_m2" -> {1,3,2,-2}, "A_p0" -> {1,3,2,0},
  "D_m2" -> {3,1,2,-2}, "D_p0" -> {3,1,2,0},
  "E_m2" -> {3,3,2,-2}, "E_p0" -> {3,3,2,0}|>;
If[!KeyExistsQ[$T3RYSpecs,$T3RYLabel],
  Print["Unsupported model ",$T3RYLabel]; Exit[2]];
$T3RYDir = DirectoryName @ ExpandFileName[$InputFileName];
$T3RYRoot = DirectoryName[$T3RYDir];
$T3RYLoad = CheckAbort[Check[Needs["RGBeta`"] ; True, $Failed],$Aborted];
If[$T3RYLoad =!= True, Print["RGBeta could not load"]; Exit[3]];
Get[FileNameJoin[{$T3RYRoot,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];

ClearAll[T3RYFail,T3RYKeyName,T3RYRead,T3RYRun];
$T3RYStep = "initialization";
T3RYFail[msg_String] := Throw[<|
  "status"->"NativeRGBetaYukawaRegistrationBlocked",
  "reason"->msg,"failed_step"->$T3RYStep,
  "native_ward_check_invoked"->False,"native_ward_check_passed"->False,
  "physical_matchete_rgbeta_coupling_map_certified"->False,
  "one_loop_beta_evaluated"->False|>,"T3RYStop"];
T3RYKeyName[k_] := Last[StringSplit[ToString[k,InputForm],"`"]];
T3RYRead[r_Association,key_String] := Module[{keys},
  keys=Select[Keys[r],T3RYKeyName[#]===key &];
  If[Length[keys]=!=1,T3RYFail["Missing or duplicate native key "<>key]];
  Lookup[r,First[keys],Missing["MissingKey",key]]
];

T3RYRun[] := Module[
  {dims,d1,d2,df,alpha,baseBuild,extra,base,extraFields,baseFields,
   inv,keysBefore,keysAfter,baseRecord,newRecord,baseInv,newInv,
   baseNativeFields,newNativeFields,reg,yScalar,yHeavy,hyperSum,
   extraScalarDimension,couplingLabel,flavourBefore,flavourAfter},
  $T3RYStep="build-native-base-model";
  dims=$T3RYSpecs[$T3RYLabel]; {d1,d2,df,alpha}=dims;
  baseBuild=CheckAbort[Check[T3RGBetaBuild[d1,d2,df,alpha],$Failed],$Aborted];
  If[!AssociationQ[baseBuild],T3RYFail["T3RGBetaBuild failed"]];

  If[alpha===0,
    extra=y1Conj; base=y1; couplingLabel="y1Conj";
    extraFields={S1,l,F}; baseFields={Bar[S1],l,F};
    extraScalarDimension=d1;
    inv=T3RGBetaYukawaInvariant[d1,df];
    yScalar=alpha/2; yHeavy=(alpha+1)/2,
    extra=y2Conj; base=y2; couplingLabel="y2Conj";
    extraFields={Bar[S2],l,Fc}; baseFields={S2,l,Fc};
    extraScalarDimension=d2;
    inv=T3RGBetaYukawaInvariant[d2,df];
    yScalar=-(alpha+2)/2; yHeavy=-(alpha+1)/2
  ];
  $T3RYStep="check-independent-gauge-eligibility";
  hyperSum=FullSimplify[-1/2+yScalar+yHeavy];
  If[!TrueQ[hyperSum===0] || df=!=2 ||
      !MemberQ[{1,3},extraScalarDimension],
    T3RYFail["Hypercharge or SU2 representation content is not eligible"]];

  $T3RYStep="check-existing-native-yukawas";
  keysBefore=Keys[$yukawas];
  If[Length[keysBefore]=!=5 || !KeyExistsQ[$yukawas,base] ||
       KeyExistsQ[$yukawas,extra],
    T3RYFail["Expected five original RGBeta Yukawas and absent extra symbol"]];
  baseRecord=Lookup[$yukawas,base,Missing["NoBase"]];
  If[!AssociationQ[baseRecord],T3RYFail["Missing baseline native record"]];
  baseNativeFields=T3RYRead[baseRecord,"Fields"];
  baseInv=T3RYRead[baseRecord,"Invariant"];
  If[!SameQ[baseNativeFields,baseFields] || !SameQ[baseInv,inv],
    T3RYFail["Baseline field list or SU2 invariant changed"]];
  flavourBefore=T3RYRead[baseRecord,"Indices"];

  (* RGBeta's CheckInvariance option is intentionally OMITTED. Its default is
     False. The observed installed checker printed residuals containing
     unevaluated internal TfLeft and TStructure for even SU3 singlets.
     We must not call that diagnostic a passed Ward identity, and must not
     patch private package implementations in this registration exercise. *)
  $T3RYStep="register-extra-yukawa-without-native-ward-option";
  reg=CheckAbort[Check[
    AddYukawa[extra,extraFields,
      GroupInvariant->inv,
      CouplingIndices->({gen[#2],heavy[#3]} &),
      Chirality->Right]; True,
    $Failed],$Aborted];
  If[reg=!=True,T3RYFail["Native AddYukawa emitted a message or failed"]];

  $T3RYStep="verify-native-registration";
  keysAfter=Keys[$yukawas];
  If[Length[keysAfter]=!=6 || !KeyExistsQ[$yukawas,extra] ||
      Complement[keysAfter,keysBefore]=!={extra},
    T3RYFail["Expected exactly one new registered Yukawa symbol"]];
  newRecord=Lookup[$yukawas,extra,Missing["NoExtra"]];
  If[!AssociationQ[newRecord],T3RYFail["Native added Yukawa record missing"]];
  newNativeFields=T3RYRead[newRecord,"Fields"];
  newInv=T3RYRead[newRecord,"Invariant"];
  flavourAfter=T3RYRead[newRecord,"Indices"];
  If[!SameQ[newNativeFields,extraFields] || !SameQ[newInv,baseInv] ||
       !SameQ[T3RYRead[newRecord,"Chirality"],Right] ||
       !SameQ[flavourAfter,flavourBefore],
    T3RYFail["Registered field order, invariant, chirality or flavour indices differ"]];
  <|
    "status"->"ExtraYukawaRegisteredNativeWardUncertified",
    "model"->$T3RYLabel,"dimensions_alpha"->dims,
    "extra_yukawa"->couplingLabel,
    "baseline_yukawa"->If[alpha===0,"y1","y2"],
    "registered_yukawas_before"->Length[keysBefore],
    "registered_yukawas_after"->Length[keysAfter],
    "extra_fields_inputform"->ToString[newNativeFields,InputForm],
    "baseline_fields_inputform"->ToString[baseNativeFields,InputForm],
    "extra_invariant_inputform"->ToString[newInv[t3ry1,t3ry2,t3ry3],InputForm],
    "baseline_invariant_inputform"->ToString[baseInv[t3ry1,t3ry2,t3ry3],InputForm],
    "same_native_su2_invariant_as_baseline"->True,
    "independent_flavour_matrix_registered"->True,
    "hypercharge_scalar_extra"->0,
    "hypercharge_fermion"->(alpha+1)/2,
    "hypercharge_sum"->ToString[hyperSum,InputForm],
    "independent_u1_charge_selection_passed"->True,
    "independent_su2_singlet_multiplicity_one"->True,
    "bsm_color_singlet_by_model_source"->True,
    "native_ward_check_invoked"->False,
    "native_ward_check_passed"->False,
    "native_ward_check_status"->"NotInvokedAfterUnresolvedA_m2SymbolicDiagnostic",
    "native_ward_check_option_disabled_explicitly"->False,
    "triplet_cross_package_basis_certified"->False,
    "physical_matchete_rgbeta_coupling_map_certified"->False,
    "physical_component_basis_certified"->False,
    "one_loop_beta_evaluated"->False,
    "production_state_modified"->False|>
];

(* Do not wrap Catch in outer Check: it would erase T3RYFail diagnostics. *)
$T3RYRes=CheckAbort[Catch[T3RYRun[],"T3RYStop"],$Aborted];
If[!AssociationQ[$T3RYRes],$T3RYRes=<|
  "status"->"NativeRGBetaYukawaRegistrationBlocked",
  "reason"->"Aborted or non-association result",
  "failed_step"->$T3RYStep,
  "native_ward_check_invoked"->False,
  "native_ward_check_passed"->False|>];
$T3RYRes=Join[<|
  "format"->"t3-workstream-a-rgbeta-extra-yukawa-registration-v1",
  "model"->$T3RYLabel|>,$T3RYRes];
If[!DirectoryQ[DirectoryName[$T3RYOut]],
  CreateDirectory[DirectoryName[$T3RYOut],CreateIntermediateDirectories->True]];
Export[$T3RYOut,$T3RYRes,"RawJSON"];
Print[$T3RYLabel,": ",Lookup[$T3RYRes,"status","Invalid"],
  "; isolated registration only; native Ward not certified"];
If[Lookup[$T3RYRes,"status",""]=!=
   "ExtraYukawaRegisteredNativeWardUncertified",Exit[6]];
Exit[0];
