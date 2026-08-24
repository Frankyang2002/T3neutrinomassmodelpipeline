(* LagrangianBuilder_TwoScalar.wl
   First genuine three-multiplet T3 builder: Scalar1 (phi'), Scalar2 (phi),
   and one fermion psi.  It keeps only interactions needed to establish the
   T3 topology plus universal norm quartics/portals.  Matchete remains the
   final judge of explicit SU(2) contractions.
*)
ClearAll[
  SU2IndexType, DefineT3TwoScalarFields, DefineT3TwoScalarCouplings,
  ScalarFieldValue, ScalarNorm, ScalarSelf, HiggsPortal, CrossScalarPortal,
  BuildT3YukawaCandidates, BuildT3MixingCandidates,
  ValidateT3Candidate, SelectFirstValidByGroup, BuildT3Lagrangian
];

SU2IndexType[1] := None;
SU2IndexType[2] := SU2L[fund];
SU2IndexType[3] := SU2L[adj];
SU2IndexType[d_] := Missing["UnsupportedSU2", d];

DefineT3TwoScalarCouplings[] := Module[{},
  DefineCoupling[y1, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[y2, Indices -> {Flavor, NFlavor}, SelfConjugate -> False];
  DefineCoupling[lambdaS1, SelfConjugate -> True];
  DefineCoupling[lambdaS2, SelfConjugate -> True];
  DefineCoupling[lambdaH1, SelfConjugate -> True];
  DefineCoupling[lambdaH2, SelfConjugate -> True];
  DefineCoupling[lambda12, SelfConjugate -> True];
  DefineCoupling[lambdaT3, SelfConjugate -> False];
  True
];

DefineT3TwoScalarFields[model_Association] := Module[
  {f, s1, s2, fidx, s1idx, s2idx, findices, s1indices, s2indices, selfConj},
  f = model["Fermion"]; s1 = model["Scalar1"]; s2 = model["Scalar2"];
  fidx = SU2IndexType[f["SU2"]];
  s1idx = SU2IndexType[s1["SU2"]];
  s2idx = SU2IndexType[s2["SU2"]];
  If[AnyTrue[{fidx,s1idx,s2idx}, MissingQ], Return[$Failed]];
  DefineFlavorIndex[NFlavor, Lookup[f,"Multiplicity",1], IndexAlphabet -> {"r","s","t"}];
  findices = If[fidx === None, {NFlavor}, {fidx,NFlavor}];
  s1indices = If[s1idx === None, {}, {s1idx}];
  s2indices = If[s2idx === None, {}, {s2idx}];
  selfConj = TrueQ[PossibleZeroQ[f["Y"]]] && MemberQ[{1,3}, f["SU2"]];
  If[selfConj,
    DefineField[NewFermion,Fermion,SelfConjugate->True,Indices->findices,Mass->{Heavy,Lookup[f,"MassSymbol",MF]}],
    DefineField[NewFermion,Fermion,SelfConjugate->False,Indices->findices,Charges->{U1Y[f["Y"]]},Mass->{Heavy,Lookup[f,"MassSymbol",MF]}]
  ];
  DefineField[NewScalar1,Scalar,SelfConjugate->False,Indices->s1indices,Charges->{U1Y[s1["Y"]]},Mass->{Heavy,Lookup[s1,"MassSymbol",MS1]}];
  DefineField[NewScalar2,Scalar,SelfConjugate->False,Indices->s2indices,Charges->{U1Y[s2["Y"]]},Mass->{Heavy,Lookup[s2,"MassSymbol",MS2]}];
  True
];

ScalarFieldValue[1, False, ___] := NewScalar1[];
ScalarFieldValue[1, True,  ___] := Bar[NewScalar1[]];
ScalarFieldValue[2, False, ___] := NewScalar2[];
ScalarFieldValue[2, True,  ___] := Bar[NewScalar2[]];
ScalarFieldValue[1, False, i_] := NewScalar1[i];
ScalarFieldValue[1, True,  i_] := Bar[NewScalar1[i]];
ScalarFieldValue[2, False, i_] := NewScalar2[i];
ScalarFieldValue[2, True,  i_] := Bar[NewScalar2[i]];

ScalarNorm[which_Integer, d_Integer, i_] := Switch[d,
  1, ScalarFieldValue[which,True] ScalarFieldValue[which,False],
  2|3, ScalarFieldValue[which,True,i] ScalarFieldValue[which,False,i],
  _, $Failed
];
ScalarSelf[which_,d_,coupling_] := Module[{i,j}, coupling[]/2 ScalarNorm[which,d,i] ScalarNorm[which,d,j]];
HiggsPortal[which_,d_,coupling_] := Module[{i,j}, coupling[] Bar[H[i]]H[i] ScalarNorm[which,d,j]];
CrossScalarPortal[d1_,d2_] := Module[{i,j}, lambda12[] ScalarNorm[1,d1,i] ScalarNorm[2,d2,j]];

(* Build one SU(2) contraction for Bar[L] . fermion . scalar, with a few
   equivalent index realizations grouped as alternatives.  fermion_ is either
   NewFermion[...] or CConj[NewFermion[...]]. *)
BuildT3YukawaCandidates[model_Association, which_Integer] := Module[
  {s, ds, df, y, scalarBar, useCConj, i,j,k,I,p,r, f1,f2, baseName},
  s = model[If[which===1,"Scalar1","Scalar2"]]; ds=s["SU2"]; df=model["Fermion","SU2"];
  (* T3 hypercharges imply: vertex 1 uses F^c S1, vertex 2 uses F S2^dagger. *)
  useCConj = which === 1; scalarBar = which === 2;
  y = If[which===1,y1,y2]; baseName = "Yukawa"<>ToString[which];
  f1[args___] := If[useCConj, CConj[NewFermion[args]], NewFermion[args]];
  Switch[{ds,df},
    {1,2},
      {
        <|"Name"->baseName<>"_Direct", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]], f1[i,r]] ScalarFieldValue[which,scalarBar]]|>,
        <|"Name"->baseName<>"_Eps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]], f1[j,r]] CG[eps[SU2L],{i,j}] ScalarFieldValue[which,scalarBar]]|>
      },
    {2,1},
      {
        <|"Name"->baseName<>"_Eps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[r]] CG[eps[SU2L],{i,j}] ScalarFieldValue[which,scalarBar,j]]|>,
        <|"Name"->baseName<>"_Delta", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[r]] ScalarFieldValue[which,scalarBar,i]]|>
      },
    {2,3},
      {
        <|"Name"->baseName<>"_GenEps", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[I,r]] CG[gen[SU2L[fund]],{I,i,k}] CG[eps[SU2L],{k,j}] ScalarFieldValue[which,scalarBar,j]]|>,
        <|"Name"->baseName<>"_Gen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[I,r]] CG[gen[SU2L[fund]],{I,i,j}] ScalarFieldValue[which,scalarBar,j]]|>
      },
    {3,2},
      {
        <|"Name"->baseName<>"_Gen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[j,r]] CG[gen[SU2L[fund]],{I,i,j}] ScalarFieldValue[which,scalarBar,I]]|>,
        <|"Name"->baseName<>"_EpsGen", "Class"->baseName, "AlternativeGroup"->baseName,
          "Expression"->PlusHc[y[p,r] NCM[Bar[l[i,p]],f1[j,r]] CG[eps[SU2L],{j,k}] CG[gen[SU2L[fund]],{I,i,k}] ScalarFieldValue[which,scalarBar,I]]|>
      },
    _, $Failed
  ]
];

(* Essential T3 scalar vertex: H H S1 S2^dagger + h.c.  The hypercharge is
   identically zero for Y1=alpha/2 and Y2=(alpha+2)/2.
   Matchete represents an unbarred fundamental field by an Index that must
   contract a barred CG index.  Hence H/S1 legs use Bar[eps], while the
   triplet projection uses gen with its fundamental/anti-fundamental slots
   oriented explicitly. *)
BuildT3MixingCandidates[model_Association] := Module[
  {d1=model["Scalar1","SU2"], d2=model["Scalar2","SU2"], i,j,k,A,B,C},
  Switch[{d1,d2},
    {2,2}, {
      <|"Name"->"T3ScalarMix_Doublets", "Class"->"T3ScalarMix", "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] H[i] Bar[NewScalar2[i]] CG[Bar[eps[SU2L]],{Bar[j],Bar[k]}] H[j] NewScalar1[k]]|>
    },
    {1,3}, {
      <|"Name"->"T3ScalarMix_1x3", "Class"->"T3ScalarMix", "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[] Bar[NewScalar2[A]] H[i] H[j] CG[gen[SU2L[fund]],{A,k,Bar[i]}] CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}]]|>
    },
    {3,1}, {
      <|"Name"->"T3ScalarMix_3x1", "Class"->"T3ScalarMix", "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[A] Bar[NewScalar2[]] H[i] H[j] CG[gen[SU2L[fund]],{A,k,Bar[i]}] CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}]]|>
    },
    {3,3}, {
      <|"Name"->"T3ScalarMix_3x3", "Class"->"T3ScalarMix", "AlternativeGroup"->"T3ScalarMix",
        "Expression"->PlusHc[lambdaT3[] NewScalar1[A] Bar[NewScalar2[B]] H[i] H[j] CG[gen[SU2L[fund]],{C,k,Bar[i]}] CG[Bar[eps[SU2L]],{Bar[k],Bar[j]}] CG[fStruct[SU2L],{A,B,C}]]|>
    },
    _, {}
  ]
];

ValidateT3Candidate[c_,LSM_,LFree_] := Module[{res},
  Print["Checking candidate: ",c["Name"]];
  res=CheckAbort[Check[CheckLagrangian[(LSM+LFree+c["Expression"])//RelabelIndices],$Failed],$Aborted];
  Print["  result: ",InputForm[res]];
  Association[c,"Valid"->TrueQ[res],"ValidationResult"->res]
];
SelectFirstValidByGroup[list_List] := Module[{seen=<||>,out={},g},
  Do[g=Lookup[x,"AlternativeGroup",None]; If[g===None||!KeyExistsQ[seen,g],If[g=!=None,seen[g]=True];AppendTo[out,x]],{x,list}]; out
];

BuildT3Lagrangian[model_Association] := Module[
  {LSM,LFree,d1,d2,candidates,validated,valid,rejected,LInt,LBSM,LUV,full,hasY1,hasY2,hasMix},
  ResetAll[];
  LSM=LoadModel["SM"]; If[LSM===$Failed,Return[$Failed]];
  If[DefineT3TwoScalarFields[model]===$Failed,Return[$Failed]];
  DefineT3TwoScalarCouplings[];
  d1=model["Scalar1","SU2"]; d2=model["Scalar2","SU2"];
  LFree=FreeLag[NewFermion,NewScalar1,NewScalar2]//RelabelIndices;
  candidates=Join[
    BuildT3YukawaCandidates[model,1], BuildT3YukawaCandidates[model,2],
    {
      <|"Name"->"Scalar1Self","Class"->"ScalarSelf","Expression"->ScalarSelf[1,d1,lambdaS1]|>,
      <|"Name"->"Scalar2Self","Class"->"ScalarSelf","Expression"->ScalarSelf[2,d2,lambdaS2]|>,
      <|"Name"->"HiggsPortal1","Class"->"Portal","Expression"->HiggsPortal[1,d1,lambdaH1]|>,
      <|"Name"->"HiggsPortal2","Class"->"Portal","Expression"->HiggsPortal[2,d2,lambdaH2]|>,
      <|"Name"->"ScalarCrossPortal","Class"->"Portal","Expression"->CrossScalarPortal[d1,d2]|>
    },
    BuildT3MixingCandidates[model]
  ];
  candidates=(Association[#,"Expression"->RelabelIndices[# ["Expression"]]]& /@ candidates);
  validated=ValidateT3Candidate[#,LSM,LFree]& /@ candidates;
  valid=SelectFirstValidByGroup@Select[validated,TrueQ[# ["Valid"]]&];
  rejected=Select[validated,!TrueQ[# ["Valid"]]&];
  LInt=Total[Lookup[valid,"Expression",{}]]//Expand//RelabelIndices;
  LBSM=(LFree+LInt)//Expand//RelabelIndices; LUV=(LSM+LBSM)//Expand//RelabelIndices;
  full=CheckAbort[Check[CheckLagrangian[LUV],$Failed],$Aborted];
  hasY1=AnyTrue[valid,Lookup[#,"Class",""]==="Yukawa1"&];
  hasY2=AnyTrue[valid,Lookup[#,"Class",""]==="Yukawa2"&];
  hasMix=AnyTrue[valid,Lookup[#,"Class",""]==="T3ScalarMix"&];
  <|
    "Status"->If[TrueQ[full],"Success","FullValidationFailed"],
    "Model"->model,"LSM"->LSM,"LFree"->LFree,"LBSM"->LBSM,"LUV"->LUV,
    "AllowedInteractions"->Lookup[valid,"Name",{}],
    "RejectedInteractions"->Lookup[rejected,"Name",{}],
    "T3IngredientsPresent"->TrueQ[hasY1&&hasY2&&hasMix],
    "WeinbergIngredientsPresent"->TrueQ[hasY1&&hasY2&&hasMix],
    "FullValidation"->full
  |>
];
