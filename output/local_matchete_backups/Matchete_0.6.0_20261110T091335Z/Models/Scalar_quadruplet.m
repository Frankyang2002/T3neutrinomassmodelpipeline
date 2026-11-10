(* ::Package:: *)

(* ::Title:: *)
(*Custodial scalar quadruplet model file*)


(* ::Subtitle:: *)
(*Custodial SM extension with two scalar quadruplets: \[CapitalTheta]1~(1,4,1/2) and \[CapitalTheta]3~(1,4,3/2)*)


ParentModel["SM"]


DefineRepresentation[SU2L[quad], SU2L,{3}]


DefineField[\[CapitalTheta]1, Scalar, Indices-> SU2L[quad], Charges-> U1Y[1/2], Mass->{Heavy, M1}];
DefineField[\[CapitalTheta]3, Scalar, Indices-> SU2L[quad], Charges-> U1Y[3/2], Mass->{Heavy, M3}];


DefineCoupling[\[Lambda]\[CapitalTheta]H];


DefineCoupling[\[Lambda]\[CapitalTheta]1, SelfConjugate-> True];
DefineCoupling[\[Lambda]\[CapitalTheta]3, SelfConjugate-> True];
DefineCoupling[\[Lambda]\[CapitalTheta]13, SelfConjugate-> True];


(*DefineCG[\[Epsilon]SU2L[quad],{SU2L[quad],SU2L[quad]}, First@InvariantTensors[SU@2,{{3},{3}},Normalization->4]]
DefineCG[CSU2L,{SU2L[quad],SU2L[adj],Bar@SU2L[fund]}, First@InvariantTensors[SU@2,{{3},{2},CRep@{1}}, Normalization->4]]*)


DefineCG[C4,{SU2L[fund],Bar@SU2L[fund],Bar@SU2L[fund],Bar@SU2L[quad]}, First@InvariantTensors[SU@2,{{1},CRep@{1},CRep@{1},CRep@{3}}, Normalization->2]]


Module[{i,j,k,a,\[Beta],M,L,\[Tau]SU2L,Htilde,L\[CapitalTheta],L\[CapitalTheta]H},	
	(*\[Tau]SU2L[aadj_,ifund_,jfund_]:= 2 CG[gen[SU2L[fund]],{aadj,ifund,jfund}];*)
	Htilde[ifund_]:= Module[{jfund},CG[eps[SU2L],{ifund,jfund}]Bar@H[jfund]];
	
	L\[CapitalTheta] = FreeLag[\[CapitalTheta]1,\[CapitalTheta]3]-\[Lambda]\[CapitalTheta]1[](Bar@\[CapitalTheta]1[M]\[CapitalTheta]1[M])(Bar@\[CapitalTheta]1[L]\[CapitalTheta]1[L])-\[Lambda]\[CapitalTheta]3[](Bar@\[CapitalTheta]3[M]\[CapitalTheta]3[M])(Bar@\[CapitalTheta]3[L]\[CapitalTheta]3[L])-\[Lambda]\[CapitalTheta]13[](Bar@\[CapitalTheta]1[M]\[CapitalTheta]1[M])(Bar@\[CapitalTheta]3[L]\[CapitalTheta]3[L]);
	(*L\[CapitalTheta]H = - \[Lambda]\[CapitalTheta]H[] (-Sqrt[3](Bar@H[i]\[Tau]SU2L[a,i,j]H[j])CSU2L[M,a,\[Beta]]Htilde[\[Beta]] Bar@\[Epsilon]SU2L[quad][M,L] \[CapitalTheta]1[L]+(Bar@H[i]\[Tau]SU2L[a,i,j]Htilde[j])CSU2L[M,a,\[Beta]]Htilde[\[Beta]] Bar@\[Epsilon]SU2L[quad][M,L] \[CapitalTheta]3[L] ) //PlusHc;*)
	L\[CapitalTheta]H = - \[Lambda]\[CapitalTheta]H[] (-Sqrt[3](Bar@H[i]H[j])C4[i,j,k,M]Htilde[k] \[CapitalTheta]1[M]+(Bar@H[i]Htilde[j])C4[i,j,k,M]Htilde[k]  \[CapitalTheta]3[M] ) //PlusHc;

	L\[CapitalTheta] + L\[CapitalTheta]H //RelabelIndices
]
