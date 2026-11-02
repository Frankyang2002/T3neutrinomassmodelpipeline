(* B2: isolated, READ-ONLY native RGBeta invariant-definition diagnostic.
   Run in fresh kernel, with T3_EXPERIMENTAL_GENERAL_* unset/0.
   wolframscript -file scripts/RunPhysicalB2RGBetaTensorExport.wl OUTPUT.json
   Never interprets a symbolic eps/del/delS2 name as measured numeric axes.
   Does not execute matching, RG evolution or write any production file. *)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] =!= 1,
  Print["Usage: OUTPUT.json"]; Exit[2]];
$B2Out = ExpandFileName[Last[$ScriptCommandLine]];
If[FileExistsQ[$B2Out], Print["Refusing existing output: ", $B2Out]; Exit[2]];
$B2Root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
$B2Step="load";
If[CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted] =!= True,
  Print["BLOCKED: RGBeta package unavailable"]; Exit[3]];
Get[FileNameJoin[{$B2Root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
$B2Models=<|"A"->{1,3,2,0},"B"->{2,2,1,-1},
             "C"->{2,2,3,-1},"D"->{3,1,2,0},"E"->{3,3,2,0}|>;
ClearAll[B2String,B2RegistryKeys,B2Read,B2Run];
B2String[expr_] := ToString[expr,InputForm];
B2RegistryKeys[x_] := If[AssociationQ[x],(ToString[#,InputForm]& /@ Keys[x]),{}];
B2Read[record_Association,n_String] := Module[{matches},
  matches=Select[Keys[record], Last[StringSplit[ToString[#,InputForm],"`"]]===n&];
  If[Length[matches]===1, Lookup[record,First[matches]], Missing["UnknownKey",n]]
];
B2Run[name_String] := Module[{d1,d2,df,alpha,yf,maj,good,rows,inv,exp,expr,
                             rho,nn,base},
  {d1,d2,df,alpha}=$B2Models[name]; yf=(alpha+1)/2;
  maj=(alpha===-1 && OddQ[df]);
  $B2Step="register "<>name;
  good=CheckAbort[Check[
    ResetModel[];T3RGBetaAddSM[];T3RGBetaAddHeavyFermion[df,yf];
    AddScalar[S1,GaugeRep->T3RGBetaGaugeRep[d1,alpha/2]];
    AddScalar[S2,GaugeRep->T3RGBetaGaugeRep[d2,(alpha+2)/2]];
    T3RGBetaAddT3Yukawas[d1,d2,df,yf];
    T3RGBetaAddMasses[d1,d2,df,yf];
    T3RGBetaAddBasicQuartics[d1,d2];
    T3RGBetaAddMixingQuartic[d1,d2]; True,$Failed],$Aborted];
  If[good=!=True,
    Return[<|"family"->name,"status"->"BlockedRegistration",
             "physical_action_certified"->False|>]];
  rows=Table[
    inv=B2Read[Lookup[$yukawas,y, <||>],"Invariant"];
    expr=If[Head[inv]===Missing,Missing["NoInvariant"],
            Quiet@Check[inv[b2a,b2b,b2c],$Failed]];
    <|"coupling"->SymbolName[y],"formal_invariant_inputform"->B2String[expr],
      "actual_registered_fields"->B2String[B2Read[Lookup[$yukawas,y,<||>],"Fields"]]|>,
    {y,{y1,y2}}];
  inv=B2Read[Lookup[$fermionMasses,MF,<||>],"Invariant"];
  expr=If[Head[inv]===Missing,Missing["NoInvariant"],
          Quiet@Check[If[df==1,inv[],inv[b2a,b2b]],$Failed]];
  <|"family"->name,"status"->"RegistryAndFormalInvariantCaptured",
    "dimensions_alpha"->{d1,d2,df,alpha},"majorana"->maj,
    "registered_yukawas"->rows,
    "registered_MF_invariant_formal"->B2String[expr],
    "registered_MF_fields"->B2String[B2Read[Lookup[$fermionMasses,MF,<||>],"Fields"]],
    "mixing_invariant_formal"->B2String[T3RGBetaMixInvariant[d1,d2][b2h,b2h2,b2s,b2t]],
    "basic_portal_singlet_formal"->B2String[
        del[SU2L[fund],b2h,b2h2] * T3RGBetaPairInvariant[d2,b2s,b2t]],
    "physical_action_certified"->False,
    "native_numeric_RGBeta_component_axes_certified"->False|>
];
$B2Data=CheckAbort[Check[
  <|"format"->"b2-native-rgbeta-formal-invariant-probe-v1",
    "status"->"FormalRegistriesCapturedPhysicalNormalisationPending",
    "RGBeta_eps_downvalues_inputform"->B2String[DownValues[eps]],
    "RGBeta_del_downvalues_inputform"->B2String[DownValues[del]],
    "RGBeta_delS2_downvalues_inputform"->B2String[DownValues[delS2]],
    "eps_example_inputform"->B2String[eps[SU2L[fund],1,2]],
    "del_triplet_example_inputform"->B2String[del[SU2L[S2],1,1]],
    "delS2_example_inputform"->B2String[delS2[SU2L,1,1,1]],
    "five_classes"->(B2Run /@ Keys[$B2Models]),
    "numerical_component_axes_verified"->False,
    "physical_source_certified"->False,
    "full_RGBeta_action_normalisation_verified"->False,
    "comment"->"Symbolic names/DownValues do not establish numerical components or physical action prefactors"|>,
  $Failed],$Aborted];
If[!AssociationQ[$B2Data],Print["BLOCKED native probe at ",$B2Step];Exit[5]];
If[!DirectoryQ[DirectoryName[$B2Out]],
 CreateDirectory[DirectoryName[$B2Out],CreateIntermediateDirectories->True]];
If[Quiet@Check[Export[$B2Out,$B2Data,"RawJSON"],$Failed]===$Failed,
 Print["BLOCKED: native JSON export failed"];Exit[6]];
Print["B2 native formal invariant export written; NO physical certificate: ",$B2Out];
Exit[0];
