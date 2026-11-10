(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`CovariantDerivative`*)


(* ::Subtitle:: *)
(*Covariant derivative rules for Matchete paclet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["CD"]
PackageExport["PartialDev"]
PackageExport["PD"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["OpenCD"]


PackageScope["TakeDev"]
PackageScope["TakePartialDev"]
PackageScope["CovDerivativeIndexQ"]


PackageScope["CommuteCDs"]
PackageScope["GAction"]
PackageScope["SeparateGeneratorsFromFS"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


CD::usage = "CD[ind,expr] returns the covariant derivative(s) of a given expression expr. The argument ind can either be a single Lorentz or Spatial index or a list of Lorentz or Spatial indices.";
PartialDev::usage = "PartialDev is a wrapper placed around Lorentz or Spatial indices in the derivative list of Field/FieldStrength to mark them as partial (rather than covariant) derivatives.";
PD::usage = "PD[ind,expr] returns the partial derivative(s) of a given expression expr. Mirrors CD; for gauge singlets it reduces to CD.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


OpenCD::usage = "OpenCD[indices] denotes open covariant derivatives that act on everything to the right of the given expression. The argument indices is must be a list of Lorentz or Spatial indices corresponding to the indices of the covariant derivatives.";


TakeDev::usage = "TakeDev[{Index[\[Mu],Lorentz|Spatial],Index[\[Nu],Lorentz|Spatial],...}, expr] is an internal function encoding the covariant derivative rules.";
TakePartialDev::usage = "TakePartialDev[{Index[\[Mu],Lorentz|Spatial],...}, expr] is an internal function encoding the partial derivative rules.";
CovDerivativeIndexQ::usage = "CovDerivativeIndexQ[ind] returns True if ind is a bare covariant derivative index (Index[_,Lorentz] or Index[_,Spatial]), False otherwise (e.g. for PartialDev-wrapped indices).";


CommuteCDs::usage               = "CommuteCDs[field,n] commutes derivatives n and n+1 (from the left) of a Field object D_{\[Mu]1, ... \[Mu]n, \[Mu](n+1), ...} f.";
GAction::usage                  = "GAction[{Index[\[Mu],Lorentz|Spatial],Index[\[Nu],Lorentz|Spatial]},field] denotes the action of a field strength tensor on a field or FS object.";
SeparateGeneratorsFromFS::usage = "Matchete uses implicit contraction of FS with generators at certain stages in the evaluation. This functions separates F_{mu nu}^a_b = F_{mu nu}^A t^{Aa}_b, with t^{Aa}_b being the appropriate generator.";


(* ::Chapter:: *)
(*Private:*)


(* ::Subsection:: *)
(*Covariant derivative definition and rules*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


CD::indices = "The index or indices '`1`' are invalid. These must be either a symbol or a list of symbols.";
PD::indices = "The index or indices '`1`' are invalid. These must be either a symbol or a list of symbols.";
TakeDev::noncovarg = "Cannot apply a covariant derivative to `1`: its derivative list already contains a PartialDev index. \!\(\*SubscriptBox[\"D\", \"\[Mu]\"]\)\!\(\*SubscriptBox[\"\[PartialD]\", \"\[Nu]\"]\)\[Phi] is not well-defined for charged fields.";


(* ::Subsubsection::Closed:: *)
(*CD definition*)


(* ::Text:: *)
(*Covariant derivative definition*)


CD[ind_, expr_]:= Module[{},
	(* Check that \[Mu] has the correct form *)
	If[
		!MatchQ[ind, _Symbol | {_Symbol..} | Index[_,Lorentz|Spatial|Temporal] | {Index[_,Lorentz|Spatial|Temporal]..} | {} | _Pattern],
		Message[CD::indices,ind];
		Abort[]
	];

	Switch[ind,
		Index[_,Lorentz|Spatial|Temporal],       
			TakeDev[{ind},expr],
		{Index[_,Lorentz|Spatial|Temporal]..},   
			TakeDev[ind,expr],
		_Symbol,                                 
			TakeDev[Index[#, Lorentz]&/@{ind},expr],
		{_Symbol..},                             
			TakeDev[Index[#, Lorentz]&/@ind,expr],
		{},                                      
			expr
	]
]


(* ::Subsubsection::Closed:: *)
(*Partial derivative definition*)


(* ::Text:: *)
(*CovDerivativeIndexQ \[LongDash] True for bare covariant Lorentz or Spatial indices, False for PartialDev-wrapped or anything else*)


CovDerivativeIndexQ[inds_List]:= And@@ CovDerivativeIndexQ/@ inds;
CovDerivativeIndexQ[Index[_, Lorentz|Spatial]]:= True;
CovDerivativeIndexQ[_]:= False;


(* ::Text:: *)
(*PD \[LongDash] user-facing partial derivative, mirrors CD*)


PD[ind_, expr_]:= Block[{},
	If[!MatchQ[ind, _Symbol | {_Symbol..} | Index[_,Lorentz|Spatial] | {Index[_,Lorentz|Spatial]..} | {} | _Pattern],
		Message[PD::indices, ind];
		Abort[]
	];
	Switch[ind,
	Index[_,Lorentz|Spatial],
		TakePartialDev[{ind}, expr]
	,{Index[_,Lorentz|Spatial]..},
		TakePartialDev[ind, expr]
	,_Symbol,
		TakePartialDev[Index[#,Lorentz]&/@{ind}, expr]
	,{_Symbol..},
		TakePartialDev[Index[#,Lorentz]&/@ind, expr]
	,{},
		expr
	]
]


(* ::Subsubsection::Closed:: *)
(*CD Test function*)


(* ::Text:: *)
(*This tells when the covariant derivative is zero*)


DZeroQ[f_? DZeroQ[x__]]:= AllTrue[{x}, DZeroQ](*And@@ DZeroQ/@ {x}*);
DZeroQ[f_[x__]]:= False;
DZeroQ[_]= True;
DZeroQ[Field]= False;
DZeroQ[FieldStrength]= False;
DZeroQ[OpenCD]= False;
DZeroQ[WilsonTerm]= False;


(* ::Subsubsection::Closed:: *)
(*Applying covariant derivative to an expression*)


(* ::Text:: *)
(*Basic rules for the covariant derivative*)


TakeDev[{}, expr_]:= expr;
TakeDev[{\[Mu]_, \[Nu]__}, expr_]:= TakeDev[{\[Mu]}, TakeDev[{\[Nu]}, expr]];
TakeDev[\[Mu]_, _?DZeroQ]= 0;
TakeDev[\[Mu]_, expr:Alternatives[_Plus, _List, _Bar, _Transp]]:= TakeDev[\[Mu], #]&/@ expr;
TakeDev[\[Mu]_, expr:Alternatives[_Times, _NCM, _NCProduct, _FuncNCM]]:=
	Module[{n}, Sum[MapAt[TakeDev[\[Mu], #]&, expr, n], {n, Length@expr}] ];
TakeDev[{\[Mu]_}, Power[expr_, n_]]:= n TakeDev[{\[Mu]}, expr] Power[expr, n-1];


(* Single combined rule for Field|FieldStrength \[LongDash] checks ordering invariant and acts *)
TakeDev[\[Mu]_, f:(Field|FieldStrength)[label_, rest__, \[Nu]_List]]:=
	If[!FreeQ[\[Nu], _PartialDev],
		Message[TakeDev::noncovarg, f]; Abort[]
	,
		Head[f][label, rest, Join[\[Mu], \[Nu]]]
	];


TakeDev[\[Mu]_,X_EoM]:= CD[\[Mu], OperatorToNormalForm@ X];
TakeDev[\[Mu]_, op_Operator]:= Operator@ CD[\[Mu], OperatorToNormalForm@ op];


TakeDev[\[Mu]_, WilsonTerm[a__, \[Nu]_]]:= WilsonTerm[a, Join[\[Mu], \[Nu]]]; 


(* derivatives acting to the right *)
TakeDev[\[Mu]_List,OpenCD[\[Nu]_List]]:= OpenCD[Join[\[Mu],\[Nu]]];


(* ::Subsubsection::Closed:: *)
(*Applying partial derivative to an expression*)


(* ::Text:: *)
(*Basic rules for the partial derivative \[LongDash] parallel to TakeDev*)


TakePartialDev[{}, expr_]:= expr;
TakePartialDev[{\[Mu]_, \[Nu]__}, expr_]:= TakePartialDev[{\[Mu]}, TakePartialDev[{\[Nu]}, expr]];
TakePartialDev[\[Mu]_, _?DZeroQ]= 0;
TakePartialDev[\[Mu]_, expr:Alternatives[_Plus, _List, _Bar, _Transp]]:= TakePartialDev[\[Mu], #]&/@ expr;
TakePartialDev[\[Mu]_, expr:Alternatives[_Times, _NCM, _NCProduct, _FuncNCM]]:=
	Module[{n}, Sum[MapAt[TakePartialDev[\[Mu], #]&, expr, n], {n, Length@expr}] ];
TakePartialDev[{\[Mu]_}, Power[expr_, n_]]:= n TakePartialDev[{\[Mu]}, expr] Power[expr, n-1];

(* Charged fields: prepend PartialDev-wrapped index. \partial_\mu = D_\mu on singlets *)
TakePartialDev[\[Mu]_, (f:(Field|FieldStrength))[label_, rest__, \[Nu]_List]]:=
	f[label, rest, Join[If[GaugeSingletQ@ label, Identity, PartialDev]/@ \[Mu], \[Nu]]];


(* ::Subsection:: *)
(*Canonize derivatives ?*)


(* ::Subsubsection::Closed:: *)
(*Action of field-strength tensor on field*)


(* ::Text:: *)
(*The action of a field strength tensor on a field/FS object*)


GAction[lInds_List, field:(Field|FieldStrength)[lab_, _, inds_, devs_]]:= Module[{charge, charges, gaugeCharges, gaugeInds, group, 
		indc, out, ind, rep, j, A},	
	
	(*Abelian FS*)
	charges= $FieldAssociation[lab, Charges];
	gaugeCharges= Cases[charges, _? (MemberQ[Keys@$GaugeGroups,Head@ #] &)];
	out= Sum[
			{group, charge}= {Head@ charge, First@ charge};
			charge 
			FieldStrength[$GaugeGroups[group, Field], lInds, {}, {}]
		, {charge, gaugeCharges}]* field;
	
	(*Non-Abelian FS*)
	gaugeInds= Cases[inds, _? (MemberQ[Keys@$GaugeGroups,GroupFromInd@ #] &)];
	out+= Sum[
		group= GroupFromInd@ ind;
		rep= Last[ind/. Bar-> Identity];
		indc= If[Head@ ind === Bar, Index[j, rep], Bar@ Index[j, rep]];
		If[rep === group@ adj, 
			I CG[fStruct@ group, {ind, Index[A, group@ adj], indc}], 
			If[Head@ ind === Bar, (*Hermitian generators*)
				- CG[gen@ rep, {Index[A, group@ adj], indc, ind}],
				CG[gen@ rep, {Index[A, group@ adj], ind, indc}] 
			]]
		FieldStrength[$GaugeGroups[group, Field], lInds, {Index[A, group@ adj]}, {}] 
		(field/. ind-> Bar@ indc)
	, {ind, gaugeInds}];
	
	out 
];


(* ::Subsubsection::Closed:: *)
(*Commute derivatives *)


(* ::Text:: *)
(*Commuting derivatives n and n+1 (from the left) of a Field object  D_{\[Mu]1, ... \[Mu]n, \[Mu](n+1), ...} f.*)


CommuteCDs::noncovdev= "CommuteCDs called on `1` with a non-covariant derivative index at position `2`."


CommuteCDs[(f:Field|FieldStrength)[rest__, devs_], n_]:= Module[{},
	If[!CovDerivativeIndexQ@ devs[[n;; n+1]],
		Message[CommuteCDs::noncovdev, f[rest, devs], n];
		Abort[];
	];
	f[rest, Join[devs[[;;n-1]], devs[[{n+1, n}]], devs[[n+2;;]] ] ] + 
		-I TakeDev[devs[[;;n-1]], GAction[devs[[{n, n+1}]], f[rest, devs[[n+2;;]] ] ] ]
]; 


(* ::Subsection:: *)
(*Separating generators from FS*)


(* ::Subsubsection::Closed:: *)
(*SeparateGeneratorsFromFS*)


(* ::Text:: *)
(*Matchete uses implicit contraction of FS with generators at stages in the evaluation. This functions separates F_{mu nu}^a_b = F_{mu nu}^A t^{Aa}_b, with t^{Aa}_b being the appropriate generator*)


SeparateGeneratorsFromFS@ HoldPattern@ FieldStrength[pre__, {ind_, indc_}, devs_]:= Module[{gr, rep, conj, A},
	gr= GroupFromInd@ ind;
	rep= Last[ind/. Bar-> Identity];
	conj= Head@ ind === Bar;
	FieldStrength[pre, {Index[A, gr@ adj]}, devs] *
		If[rep === gr@ adj, 
			I CG[fStruct@ gr, {ind, Index[A, gr@ adj], indc}], 
			If[conj, 
				- CG[gen@ rep, {Index[A, gr@ adj], indc, ind}],
				CG[gen@ rep, {Index[A, gr@ adj], ind, indc}] 
			]]
];


SeparateGeneratorsFromFS@ fs:HoldPattern@ FieldStrength[pre__, {}|{_}, devs_]:= fs;


SeparateGeneratorsFromFS@ expr_:= expr/. fs_FieldStrength:> SeparateGeneratorsFromFS@ fs;
