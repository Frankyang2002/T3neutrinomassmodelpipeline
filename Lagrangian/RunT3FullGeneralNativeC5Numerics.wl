(* Evaluate a complete Matchete-native two-threshold C5 coefficient in the
   fixed source operator convention.  Only valid if all heavy fields are
   removed, all index sums are explicitly contracted and all poles cancel.
   This evaluator does NOT repair a missing EFT1 Wilson RGE or soft subtraction.
   Usage: wolframscript -file SCRIPT COEFFICIENT.wxf INPUT.json OUT_DIR
*)
ClearAll["Global`*"];
nArgs=Rest[$ScriptCommandLine];
If[Length[nArgs]=!=3, Print["Usage: C5.wxf INPUT.json OUT_DIR"]; Exit[2]];
nInput=ExpandFileName[nArgs[[1]]];
nConfig=ExpandFileName[nArgs[[2]]];
nOut=ExpandFileName[nArgs[[3]]];
If[!DirectoryQ[nOut],CreateDirectory[nOut,CreateIntermediateDirectories->True]];
nReport=FileNameJoin[{nOut,"native_c5_numeric_status.json"}];
nStage="init";
nStop[text_String] := (Export[nReport,<|
 "status"->"Blocked","stage"->nStage,"reason"->text,
 "numeric_C5_validated"->False,"full_EFT1_Wilson_RGE_performed"->False|>,
 "RawJSON"]; Print["C5 BLOCKED: ",nStage,": ",text];Exit[5]);
If[!FileExistsQ[nInput]||!FileExistsQ[nConfig],nStop["Missing native WXF or input"]];
nData=Import[nConfig,"RawJSON"];
If[!AssociationQ[nData]||Lookup[nData,"format",""]=!=
   "t3-complete-native-c5-evaluation-input-v1",nStop["Malformed numeric input"]];
nStage="load-Matchete";
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"];True],$Aborted]=!=True,
 nStop["Matchete not loaded"]];
nStage="load-source";
nExpr=Import[nInput,"WXF"];
If[MemberQ[{$Failed,$Aborted},nExpr]||FreeQ[nExpr,Matchete`Coupling],
 nStop["Missing source coefficient or Matchete coupling objects"]];
(* The original C5 convention is ordered A5, with C5(p,q)=A5(p,q)+A5(q,p). *)
nNames[x_]:=ToString[InputForm[x]];
nName[s_Symbol]:=SymbolName[Unevaluated[s]];
nName[x_]:=nNames[x];
nYuk=Lookup[nData,"yukawas",<||>];
nYuk=Map[Map[Function[e,If[AssociationQ[e] &&
   KeyExistsQ[e,"real"] && KeyExistsQ[e,"imag"],
     N[e["real"]]+I N[e["imag"]],e]],#,{2}]&,nYuk];
nScal=Lookup[nData,"scalar_couplings",<||>];
nMass=Lookup[nData,"fermion_masses_GeV",{}];
nMu=Lookup[nData,"mu_S_GeV",0];
If[Length[nMass]=!=3||!VectorQ[nMass,NumericQ]||Min[nMass]<=0||
  !NumericQ[nMu]||nMu<=0,nStop["Missing 3 positive fermion masses or scalar scale"]];
If[!AllTrue[{"y1","y2"},KeyExistsQ[nYuk,#]&],
 nStop["Missing y1/y2 3x3 Yukawa matrices"]];
If[!AllTrue[Values[nYuk],MatrixQ[#,NumericQ]&&Dimensions[#]==={3,3}&],
 nStop["Yukawa entries must be finite numeric 3x3"]];
nStage="prepare-symbolic-flavour";
nIsType[index_,label_String]:=MatchQ[index,
 Matchete`Index[_,r_Symbol]/;SymbolName[r]===label];
nFieldFlav[term_]:=Flatten@Cases[term,
 Matchete`Field[name_Symbol,type_Symbol,indices_,___]/;
   SymbolName[name]==="l"&&SymbolName[type]==="Fermion":>
   Cases[indices,Matchete`Index[v_,rep_Symbol]/;
     SymbolName[rep]==="Flavor":>v,Infinity],Infinity];
nTerms[expr_]:=If[Head[Expand[expr]]===Plus,List@@Expand[expr],{Expand[expr]}];
(* Resolve indices PER MONOMIAL.  Different operators can reuse d$$1 as
   independent bound names; summing over the whole polynomial is incorrect. *)
nExpandSums[expr_]:=FixedPoint[Function[z,
  z/. HoldPattern[Sum[body_,{r_Symbol,lo_Integer,hi_Integer}]]:>
   Total[Table[body/.r->k,{k,lo,hi}]]],expr,10];
nExternal[expr_,p_Integer,q_Integer]:=Module[{labels,ex,px,qx},
 labels=DeleteDuplicates@Cases[expr,
  Matchete`Index[i_Symbol,rep_Symbol]/;
    SymbolName[rep]==="Flavor":>i,Infinity];
 px=Select[labels,nName[#]==="T3D2ExternalFlavorP"&];
 qx=Select[labels,nName[#]==="T3D2ExternalFlavorQ"&];
 If[Length[px]===1&&Length[qx]===1,
  Return[expr/.{
   Matchete`Index[i_Symbol,r_Symbol]/;
     SymbolName[r]==="Flavor"&&nName[i]==="T3D2ExternalFlavorP":>
      Matchete`Index[p,r],
   Matchete`Index[i_Symbol,r_Symbol]/;
     SymbolName[r]==="Flavor"&&nName[i]==="T3D2ExternalFlavorQ":>
      Matchete`Index[q,r]}]];
 (* Standard native Weinberg extractor leaves d$$1 and d$$2 as the two
    external Yukawa-lepton flavours.  Permit this only if the ENTIRE source
    has exactly those two typed Flavor labels; any extra flavour sums need
    a separate explicit projection. *)
 If[Length[labels]===2 &&
    Sort[nName /@ labels] === {"d$$1","d$$2"},
   Return[expr /. {
     Matchete`Index[i_Symbol,r_Symbol]/;
       SymbolName[r]==="Flavor"&&nName[i]==="d$$1":>Matchete`Index[p,r],
     Matchete`Index[i_Symbol,r_Symbol]/;
       SymbolName[r]==="Flavor"&&nName[i]==="d$$2":>Matchete`Index[q,r]}]
 ];
 nStop["Cannot uniquely identify two ordered external lepton flavours"]
];
nContractTerm[term_]:=Module[{inds,res=term,unknown={},r},
 inds=DeleteDuplicates@Cases[term,
   Matchete`Index[s_Symbol,t_Symbol]/;
      SymbolName[t]==="NFlavor":>s,Infinity];
 If[Length[inds]>6,nStop["Excessive independent NFlavor dummy labels"]];
 Do[
   With[{index=inds[[k]]},
      res=Total[Table[res/.Matchete`Index[index,rep_Symbol]/;
       SymbolName[rep]==="NFlavor":>Matchete`Index[v,rep],
       {v,1,3}]]],
   {k,Length[inds]}];
 res
];
nYLookup[name_String, p_Integer,r_Integer]:=Module[{y=Lookup[nYuk,name,{}]},
 If[!MatrixQ[y,NumericQ]||Dimensions[y]!={3,3}||
  !Between[p,{1,3}]||!Between[r,{1,3}],nStop["Unknown indexed Yukawa"]];
 y[[p,r]]];
nReplaceCouplings[expr_]:=Module[{res=expr},
 res=res/.Matchete`Coupling[s_Symbol,indices_List,order_]:>Module[
  {name=nName[s],flav,heavy},
  If[indices==={},
   If[KeyExistsQ[nScal,name],Return[nScal[name]]];
   nStop["Unmapped scalar/mass coupling "<>name]];
  flav=Cases[indices,Matchete`Index[v_Integer,t_Symbol]/;
   SymbolName[t]==="Flavor":>v,Infinity];
  heavy=Cases[indices,Matchete`Index[v_Integer,t_Symbol]/;
   SymbolName[t]==="NFlavor":>v,Infinity];
  If[name==="MF"&&Length[heavy]===1,Return[nMass[[First[heavy]]]]];
  If[KeyExistsQ[nYuk,name]&&Length[flav]===1&&Length[heavy]===1,
    Return[nYLookup[name,First[flav],First[heavy]]]];
  nStop["Unmapped indexed coupling "<>name]
 ];
 FixedPoint[#/.{Matchete`Bar[z_?NumericQ]:>Conjugate[z],
              Matchete`Bar[Matchete`Bar[z_]]:>z}&,res,4]
];
nScaleSymQ[s_Symbol]:=MemberQ[{"Mubar2",FromCharacterCode[956]<>"bar2"},nName[s]];
nScaleSymQ[_]:=False;
nRegQ[s_Symbol]:=MemberQ[{"ϵ","ε","epsilon","eps"},nName[s]];
nRegQ[_]:=False;
nF=Lookup[nData,"NF",3];
If[nF=!=3,nStop["Three fermion generations required"]];
nRaw=nExpandSums[nExpr/.{
 s_Symbol/;nName[s]==="NF":>3,
 s_Symbol/;nName[s]==="hbar":>1/(16 Pi^2),
 s_Symbol/;nScaleSymQ[s]:>nMu^2}];
If[!FreeQ[nRaw,Sum],nStop["Explicit heavy-flavour Sum did not resolve"]];
nStage="evaluate-nine-native-flavour-entries";
nEval[p_Integer,q_Integer]:=Module[{expr=nExternal[nRaw,p,q],terms,contracted,coupled,
  eps,res,fin},
 terms=nTerms[expr];
 contracted=Total[nContractTerm/@terms];
 coupled=nReplaceCouplings[contracted];
 If[!FreeQ[coupled,Matchete`Coupling|Matchete`Index|Matchete`Field|Matchete`CG],
  nStop["Unresolved native field, index, CG or coupling"]];
 eps=DeleteDuplicates@Cases[coupled,s_Symbol/;nRegQ[s]:>s,Infinity];
 If[Length[eps]>1,nStop["More than one dimensional regulator"]];
 If[Length[eps]===0,Return[{N[coupled,20],0}]];
 res=Normal@Series[coupled,{First[eps],0,0}];
 {N[Coefficient[res,First[eps],0],20],
  N[Coefficient[res,First[eps],-1],20]}
];
nValues=Table[nEval[p,q],{p,1,3},{q,1,3}];
nOrdered=nValues[[All,All,1]];
nPole=nValues[[All,All,2]];
If[!MatrixQ[nOrdered,NumericQ]||!MatrixQ[nPole,NumericQ],
 nStop["Non-numeric output after complete three-flavour contraction"]];
nSym=Chop[nOrdered+Transpose[nOrdered],10^-25];
nPoleSym=Chop[nPole+Transpose[nPole],10^-25];
nPoleMax=Max[Abs[Flatten[nPoleSym]]];
nScale=Max[10^-100,Max[Abs[Flatten[nSym]]]];
nPoleRel=N[nPoleMax/nScale,12];
nStage="pole-subtraction-gate";
nEncode[x_]:=<|"real"->N[Re[x],15],"imag"->N[Im[x],15]|>;
nResult=<|
 "format"->"t3-three-flavour-native-two-threshold-c5-numerical-v1",
 "status"->If[TrueQ[nPoleRel<10^-9],
  "ThreeFlavourNativeC5FiniteCandidate_NotEFT1RGVerified",
  "Blocked_UncancelledDimensionalPole"],
 "mu_S_GeV"->nMu,"all_heavy_fermions_integrated_together"->True,
 "both_scalars_integrated_together"->True,
 "ordered_C5_candidate_GeV_inverse"->nEncode[nOrdered],
 "symmetric_C5_candidate_GeV_inverse"->nEncode[nSym],
 "pole_residue_GeV_inverse"->nEncode[nPoleSym],
 "pole_relative_size"->nPoleRel,
 "physical_gamma_CC_normalization_verified"->False,
 "full_eft1_wilson_RGE_performed"->False,
 "complete_MSbar_subtraction_verified"->False,
 "unknown_couplings_defaulted_to_zero"->False,
 "numerical_neutrino_fits_performed"->False
|>;
Export[nReport,nResult,"RawJSON"];
Print["Native C5 pole relative size: ",nPoleRel];
If[nResult["status"] =!= "ThreeFlavourNativeC5FiniteCandidate_NotEFT1RGVerified",
 Exit[5]];
Exit[0];
