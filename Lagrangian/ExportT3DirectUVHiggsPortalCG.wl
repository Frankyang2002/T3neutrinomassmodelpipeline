(* Read-only native Matchete Higgs-portal CG export for the direct UV one-loop audit.

   wolframscript -file Lagrangian/ExportT3DirectUVHiggsPortalCG.wl \
       dS1 dS2 dF "(-2)" OUTPUT.json

   Alpha MUST be parenthesized at the CLI, including negative values, because
   wolframscript may interpret a bare negative integer as an option.

   This exports ACTUALLY REGISTERED HiggsPortal{1,2}Inv{k} tensors from
   BuildT3Lagrangian. It DOES NOT identify the native lambdaH{1,2}Inv{k}
   coefficients with the fitted RGBeta lambdaH{1,2} values.
   No matching, RGE, or physics result is overwritten.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] < 5,
 Print["Usage: ExportT3DirectUVHiggsPortalCG.wl dS1 dS2 dF alpha OUTPUT.json"];
 Exit[2]
];
{d1arg,d2arg,dfarg,aarg,outarg} = Take[$ScriptCommandLine,-5];
If[!StringMatchQ[ToString[aarg], "(" ~~ ("-" | "") ~~ DigitCharacter.. ~~ ")"],
 Print["ERROR: alpha must be parenthesized, e.g. \"(-2)\" or \"(0)\"."];
 Exit[2]
];
{d1,d2,df,alpha} = Quiet@Check[ToExpression /@ {d1arg,d2arg,dfarg,aarg},$Failed];
If[!And@@(IntegerQ/@{d1,d2,df,alpha}),Print["Invalid integer input"];Exit[2]];
out = ExpandFileName[outarg];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[!DirectoryQ[DirectoryName[out]],
 CreateDirectory[DirectoryName[out],CreateIntermediateDirectories->True]
];
load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"];True],$Failed],$Aborted];
If[load =!= True,Print["Matchete could not be loaded"];Exit[3]];
Get[FileNameJoin[{root,"Lagrangian","model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{root,"Lagrangian","model","LagrangianBuilder.wl"}]];
model = T3ModelFromDimensions[d1,d2,df,alpha];
If[model === $Failed,Print["Unrecognized model"];Exit[4]];
build = CheckAbort[Check[UsingFrontEnd[BuildT3Lagrangian[model]],$Failed],$Aborted];
If[!AssociationQ[build] || Lookup[build,"Status",""] =!= "Success" ||
   !TrueQ[Lookup[build,"FullValidation",False]],
 Print["Native UV model builder failed"];Exit[5]
];
(* The possible second invariant is model-dependent: detect ACTUAL accepted
   portal candidates and require indices 1..n without a gap. *)
portalNames = Select[Lookup[build,"AllowedInteractions",{}],
 StringMatchQ[#,"HiggsPortal1Inv"~~DigitCharacter..] ||
 StringMatchQ[#,"HiggsPortal2Inv"~~DigitCharacter..]&];
If[Length[portalNames]<2 || Length[portalNames]>4,
 Print["Unexpected native Higgs-portal inventory: ",portalNames];Exit[6]
];
portalRegistry = Table[
 Module[{which,k,cgName,cgSym,reps,tensor,parts},
  parts = StringCases[name,
    "HiggsPortal"~~w:DigitCharacter~~"Inv"~~n:DigitCharacter.. :>
     {ToExpression[w],ToExpression[n]}];
  If[Length[parts] =!= 1,Return[$Failed]];
  {which,k} = First[parts];
  cgName = "T3HiggsPortal"<>ToString[which]<>"CGInv"<>ToString[k];
  cgSym = Symbol[cgName];
  reps = Quiet@Check[Matchete`PackageScope`$CGproperties[cgSym,Indices],$Failed];
  tensor = Quiet@Check[Matchete`CGManipulations`PackagePrivate`$CGtensors[cgSym],$Failed];
  If[reps === $Failed || tensor === $Failed || MissingQ[reps] || MissingQ[tensor],
   Print["Missing registered native tensor: ",cgName];Return[$Failed]
  ];
  <|"Name"->cgName,"Coupling"->("lambdaH"<>ToString[which]<>"Inv"<>ToString[k]),
    "SourceInteraction"->name,"WhichScalar"->which,"Invariant"->k,
    "Dimensions"->{2,2,If[which==1,d1,d2],If[which==1,d1,d2]},
    "RepsInputForm"->ToString[InputForm[reps]],
    "TensorInputForm"->ToString[InputForm[tensor]]|>
 ], {name,portalNames}];
If[MemberQ[portalRegistry,$Failed],Exit[7]];
Do[
 actual=Sort[Lookup[Select[portalRegistry,Lookup[#,"WhichScalar"]===which&],"Invariant",{}]];
 If[actual=!=Range[Length[actual]],
  Print["Missing native invariant or numbering gap: scalar ",which," ",actual];Exit[8]
 ], {which,1,2}];
payload = <|
 "status"->"Success",
 "format"->"t3-native-direct-uv-portal-cg-v1",
 "metadata"-><|"dS1"->d1,"dS2"->d2,"dF"->df,"alpha"->alpha|>,
 "CGRegistry"->portalRegistry,
 "physical_rgbeta_to_matchete_portal_coupling_dictionary_certified"->False,
 "portal_lagrangian_sign"->"+1/2 PlusHc[lambdaH{which}Inv{k} Hdag H Sdag S CG]",
 "vev_quadratic_in_potential"->"-lambdaH{which}Inv{k} v^2/2 * native_CG_H0H0_sbar_s",
 "match_or_running_performed"->False
 |>;
Export[out,payload,"RawJSON"];
Print["Portal CG registry export: Success; tensor count=",Length[portalRegistry]];
Print["Output: ",out];
