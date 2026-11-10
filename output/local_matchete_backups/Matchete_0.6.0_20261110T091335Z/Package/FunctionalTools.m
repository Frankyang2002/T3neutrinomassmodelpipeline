(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`FunctionalTools`*)


(* ::Subtitle:: *)
(*Sub package providing the functional tools used in functional one-loop matching.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["FuncNCM"]


PackageScope["DeriveEOM"]
PackageScope["FluctuationOperator"]


PackageScope["VarD"]
PackageScope["FD"]
PackageScope["FuncD"]
PackageScope["FuncDSimplify"]
PackageScope["ShiftCD"]
PackageScope["ShiftFS"]


PackageScope["TransposeThisSpinChain"]


PackageScope["SubstituteFieldStrength"]
PackageScope["SubstituteHeavyFieldStrength"]


PackageScope["GaugeFieldQ"]


PackageScope["ExpandVectorFluctuations"]


PackageScope["ChangeFieldIndex"]


PackageScope["NCProduct"]


PackageScope["FieldTransformsUnderGaugeGroupQ"]


PackageScope["BkgField"]
PackageScope["BackgroundCD"]
PackageScope["BackgroundFS"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


FD::usage =
"FD[expr,field] calculates the (partial) functional derivative of expr with respect to field, where expr is a Lagrangian or any other expression depending on fields. The argument field must be of type Field[label, type, {indices}, {CDerivs}] and the diferent covariant derivatives of the field are treated as independent.
Example: FD[\[ScriptCapitalL], Field[\[Phi], Scalar, {}, {\[Mu],\[Nu]}]] = \!\(\*FractionBox[\(\[Delta]\[ScriptCapitalL]\), \(\[Delta] \((\*SubscriptBox[\(D\), \(\[Mu]\)] \*SubscriptBox[\(D\), \(\[Nu]\)] \[Phi])\)\)]\).";


FuncD::usage =
"FuncD[expr,field] calculates the (partial) functional derivative of expr with respect to field, where expr is a Lagrangian or any other expression depending on fields. Results are returned in a non simplified form. The argument field must be of type Field[label, type, {indices}, {CDerivs}] and the diferent covariant derivatives of the field are treated as independent.
Example: FuncD[\[ScriptCapitalL], Field[\[Phi], Scalar, {}, {\[Mu],\[Nu]}]] = \!\(\*FractionBox[\(\[Delta]\[ScriptCapitalL]\), \(\[Delta] \((\*SubscriptBox[\(D\), \(\[Mu]\)] \*SubscriptBox[\(D\), \(\[Nu]\)] \[Phi])\)\)]\).";


FuncDSimplify::usage =
"FuncDSimplify[expr] simplifies the expression expr that was returned by FuncD. It contracts all indices in Delta[_,_] and Metric[_,_] and transposes spin-chains if necessary.  ";


VarD::usage=
"VarD[expr,field1,field2(optional)] calculates the variational derivative of expr with respect to field2 and then field1, where expr is a Lagrangian or any other expression depending on fields. The arguments field1 and field2 must be of the form Field[label, type, {indices}, {}].";


Grassmann::usage =
"Grassmann is an option for FuncD and VarDraw that specifies whether minus signs from anti-commuting functional derivatives of fermionic fields should be considered. The default value is Grassmann->True.";


DeriveEOM::usage =
"DeriveEOM[\[ScriptCapitalL],\[Phi]] derives the equation of motion for the field \[Phi] from the Lagrangian density \[ScriptCapitalL], i.e. VarD[\[ScriptCapitalL],Bar[\[Phi]]] is computed.";


FluctuationOperator::usage =
"FluctuationOperator[\[ScriptCapitalL],\[Phi],\[Psi]] derives the fluctuation operator, i.e. the X-term, for the fields \[Phi] and \[Psi] from the Lagrangian density \[ScriptCapitalL], i.e. \!\(\*FractionBox[\(\*SuperscriptBox[\(\[Delta]\), \(2\)] \[ScriptCapitalL]\), \(\[Delta]\[Phi]\\\ \[Delta]\[Psi]\)]\) ~ VarD[\[ScriptCapitalL],\[Psi],\[Phi]] is computed.";


RecursiveExpandCD::usage =
"RecursiveExpandCD[indices,group,field] shifts all covariant derivatives by the gauge field fluctuation and expands the result up to leading order in the gauge field fluctuation.";


SubstituteFieldStrength::usage =
"SubstituteFieldStrength[expr] substitues all instances of FieldStrength[...] in expr with the corresponding Fields using \!\(\*SubscriptBox[SuperscriptBox[\(F\), \(a\)], \(\[Mu]\[Nu]\)]\)=(\!\(\*SubscriptBox[\(D\), \(\[Mu]\)]\)\!\(\*SubscriptBox[\(A\), \(\[Nu]\)]\)\!\(\*SuperscriptBox[\()\), \(a\)]\)-(\!\(\*SubscriptBox[\(D\), \(\[Nu]\)]\)\!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\)\!\(\*SuperscriptBox[\()\), \(a\)]\)."


ChangeFieldIndex::usage="ChangeFieldIndex[field, group, indexLabel] returns the index of field that belongs to the given group and returns the field with that index replaced by the given indexLabel.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Functional Non-commutative product *)


(* ::Text:: *)
(*A non-commutative product used in the evaluation on functional derivatives*)


(* ::Subsubsection::Closed:: *)
(*Commutative  check*)


(* ::Text:: *)
(*Default assumption: everything is commutative*)


FuncCommuteQ[f_?FuncCommuteQ[x___]] := AllTrue[{x}, FuncCommuteQ];
FuncCommuteQ[f_[x___]] := False;
FuncCommuteQ[_] := True;


(* ::Text:: *)
(*Define the non-commutative objects*)


FuncCommuteQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold] = False;


(* ::Text:: *)
(*Commutative is a head that can be used to temporarily treat an object as commutative*)


FuncCommuteQ@ Commutative@ _ ^= True;


(* ::Subsubsection::Closed:: *)
(*Non - commutative  symbols  used  in  Matchete  declared  here*)


(* ::Text:: *)
(*Define non-commutative objects*)


((FuncCommuteQ@ # ^= False) &) /@ {
	Field,
	FieldStrength,
	DiracProduct,
	OpenCD,
	Xop, Xterm,
	WilsonLine, WilsonTerm
	};


(* ::Subsubsection::Closed:: *)
(*Properties  of  FuncNCM*)


(* ::Text:: *)
(*Flatness*)


FuncNCM[a___, FuncNCM[b__], c___]:= FuncNCM[a, b, c];
FuncNCM[] = 1;


(* ::Text:: *)
(*Extracting commuting objects*)


FuncNCM[a___, b_?FuncCommuteQ, c___]:= b * FuncNCM[a, c];
FuncNCM[a___, b_?FuncCommuteQ * x_, c___]:= b * FuncNCM[a, x, c];


(* ::Text:: *)
(*Distributivity*)


FuncNCM[a___, b_Plus, c___]:= FuncNCM[a, #, c] & /@ b


(* ::Section:: *)
(*FD = FuncDSimplify@FuncD (partial functional derivative)*)


(* ::Text:: *)
(*This function calls FuncD (which does the actual partial functional derivative) and applies some simplifications to it.*)


(* ::Subsubsection::Closed:: *)
(*Main entry point for computing partial functional derivatives *)


(* ::Text:: *)
(*Power counting parameter*)


$\[Epsilon]FD::usage="$\[Epsilon]FD is a power counting parameter for the expansion of gauge covariant derivatives.";


(* directly drop order $\[Epsilon]FD^3 terms *)
$\[Epsilon]FD /: Power[$\[Epsilon]FD,n_/;n>2] := 0


(* ::Text:: *)
(*Option test*)


Options[FD] = {Grassmann->True};


OptionTest[_, Grassmann] = BooleanQ;


OptionMessage[Grassmann, func:FD, val_] := Message[General::optexpectsval, Grassmann, func, val, "boolean (True or False)"];


FD[x_, y_, OptionsPattern[]] := Module[
	{
		result,
		arg = Expand[x]
	}
	,

	(* computation of the partial functional derivative and simplifications *)
	result = FuncDSimplify@FuncD[arg, y, Grassmann->OptionValue[Grassmann]];

	Return[result]
]


(* ::Text:: *)
(*Check arguments*)


FD::invalidargument =
"The second argument `1` of FuncD is invalid. It must be either of the following forms: Field[...], Bar[Field[...]].";


FD[_,y_,OptionsPattern[]] := Message[FD::invalidargument,y] /; !MatchQ[y, Field[___] | Bar@Field[___] | Transp@Field[___] | Transp@Bar@Field[___] | HoldComplete[___] | NCM[DiracProduct[GammaCC],_] | NCM[_,DiracProduct[GammaCC]]]


(* ::Subsection:: *)
(*Expansion of Field Strength tensors and Covariant Derivatives*)


(* ::Subsubsection::Closed:: *)
(*Expand fluctuations of vector fields*)


ExpandVectorFluctuations[lag_, f1_, f2__]:= ExpandVectorFluctuations[ExpandVectorFluctuations[lag, f2], f1]


(* This functions expands out the vector field fluctuations of the vector field y in Field Strengh tensors and Covariant Derivatives *)
ExpandVectorFluctuations[x_, y_]:=
ExpandVectorFluctuations[x, y]= Module[
	{
		gaugeGroup,
		yField = First@Cases[y,Field[___],All],
		yFieldLabel,
		arg = LagrangianExpand[x]
	}
	,
	yFieldLabel = First[yField];

	(* Shift covariant derivatives and field-strength tensors in the case of gauge fields *)
	If[GaugeFieldQ@yField,

		(* shift explicit gauge fields *)
		arg = LagrangianExpand[arg /. Field[yFieldLabel,type_Vector,rest___]:>(Field[BkgField[yFieldLabel],type,rest]+$\[Epsilon]FD*Field[yFieldLabel,type,rest])];

		(* remove powers for pattern matching *)
		arg = RemovePower[arg];

		(* determine gauge group *)
		gaugeGroup = First@First@FirstPosition[$GaugeGroups, KeyValuePattern[Field->First[yField]]];

		(* shift covariant derivatives iff field is charged under gauge group *)
		arg = arg/.Field[l_,t_,ind_List,deriv_List] :> ShiftCD[deriv,gaugeGroup,Field[l,t,ind,{}],$\[Epsilon]FD] /; (Length[deriv]>0 && FieldTransformsUnderGaugeGroupQ[Field[l,t,ind,{}], gaugeGroup]);
		arg = RemovePower@ReleaseHold[arg];

		(* shift covariant derivatives acting on field-strength tensors *)  (*crosscheck*)
		arg = arg/.FieldStrength[l_,lorentz_,ind_List,deriv_List] :> ShiftCD[deriv,gaugeGroup,FieldStrength[l,lorentz,ind,{}],$\[Epsilon]FD] /; (Length[deriv]>0 && FieldTransformsUnderGaugeGroupQ[FieldStrength[l,lorentz,ind,{}], gaugeGroup]);

		(* shift the field-strength tensors *)
		arg = arg /. FieldStrength[arguments___] :> ShiftFS[FieldStrength[arguments],$\[Epsilon]FD] /; MatchQ[First[{arguments}],$GaugeGroups[gaugeGroup][Field]];
		arg = RemovePower@ReleaseHold[arg];

	];

	(* Shift only the field-strength tensors in the case of vectors that are not gauge *)
	If[VectorFieldQ@yField && !GaugeFieldQ@yField,

		(* shift explicit vector fields *)
		arg = LagrangianExpand[arg /. Field[yFieldLabel,type_Vector,rest___]:>(Field[BkgField[yFieldLabel],type,rest]+$\[Epsilon]FD*Field[yFieldLabel,type,rest])];

		(* remove powers for pattern matching *)
		arg = RemovePower[arg];

		(* shift the field-strength tensors *)
		arg = arg /. FieldStrength[First@yField,args___] :> ShiftFS[FieldStrength[First@yField,args],$\[Epsilon]FD];
	];

	(*Release BackgroundCD*)
	arg = arg/. BackgroundCD[ind__,exp_]:>CD[ind,exp] /.BackgroundFS[all___]:>FS[all];

	(* expand to quadratic order in the fluctuation *)
	arg = LagrangianExpand[ReleaseHold@arg];
	arg = BetterSeries[#,{$\[Epsilon]FD,0,2}]&/@arg;

	RelabelIndices[arg(*, Unique->True*)] (* Uncommenting this option makes SuperSoftMatch extremely slow (investigate) *)
]


(* ::Subsubsection::Closed:: *)
(*Check if a given field transforms under a given gauge group*)


FieldTransformsUnderGaugeGroupQ[(Field|FieldStrength)[label_, _, inds_, _], gaugeGroup_Symbol]:= 
	If[$GaugeGroups[gaugeGroup, Abelian],
		MemberQ[$FieldAssociation[label, Charges][[;;, 0]], gaugeGroup]
	,
		MemberQ[GroupFromRep/@ $FieldAssociation[label, Indices], gaugeGroup]
	];


(* ::Subsubsection::Closed:: *)
(*Auxiliary function to check for gauge fields and heavy vectors*)


(* ::Text:: *)
(*gauge fields*)


GaugeFieldQ@ label_Symbol:= MemberQ[$GaugeGroups, KeyValuePattern[Field-> label] ];


GaugeFieldQ@ Field[label_, _Vector, __]:= GaugeFieldQ@ label;


GaugeFieldQ@ _Field= False;


(* ::Text:: *)
(*heavy vectors*)


VectorFieldQ@ Field[_, type_, __]:= MatchQ[type, Vector[_]]


(* ::Subsubsection::Closed:: *)
(*Shifting all covariant derivatives containing gauge field fluctuations*)


(* returns true if the given gauge group is abelian *)
AbelianQ[gaugeGroup_] := $GaugeGroups[gaugeGroup, Abelian];


(* Expand the quantum in fluctuations of the gauge field associated to the given group in covariant derivatives *)
ShiftCD[indices_List, group_, field_, globalCounter_:1] := Module[{result,counter},

	(* if there is no CD just return the field*)
	If[Length[indices]==0, Return[field]];

	(* else: shift all the CDs *)
	If[AbelianQ[group],
		result = Expand@RecursiveExpandAbelianCD[indices, group, field, counter];
		,
		result = Expand@RecursiveExpandCD[indices, group, field, counter];
	];

	(* fluctuations should not be truncated to leading order when multiple derivatives are present *)
	(*
	(* truncate at leading order in the fluctuation *)
	result = BetterSeries[result,{counter,0,1}];
	*)

	(* remove counting parameter *)
	result = Expand[result/.counter->globalCounter];

	result
]


(* ::Text:: *)
(*Non-Abelian groups*)


RecursiveExpandCD::partialdev = "RecursiveExpandCD does not support PartialDev indices (not implemented).";
RecursiveExpandAbelianCD::partialdev = "RecursiveExpandAbelianCD does not support PartialDev indices (not implemented).";


RecursiveExpandCD[indices_List, group_, initField_, \[Lambda]_:1] := Module[
	{
		\[Mu]1,
		more\[Mu],
		indexAdj=Unique[],
		newRepIndex=Unique[],
		generatorIndices,
		field
	}
	,
	(* replace the group index in field *)
	{field,generatorIndices}=ChangeFieldIndex[initField, group, newRepIndex];

	(* split the Lorentz (or Spatial) indices *)
	{\[Mu]1,more\[Mu]} = TakeDrop[indices,1];
	If[!CovDerivativeIndexQ[First[\[Mu]1]], Message[RecursiveExpandCD::partialdev]; Abort[]];
	(* recursive formula for building up the fluctuation expansion in covariant derivatives *)
	If[Length[more\[Mu]]>=1,
		Plus[
			Times[-1,
				\[Lambda],I,
				(*$GaugeGroups[group][Coupling][],*)
				$GaugeGroups[group, Field][\[Mu]1[[1,1]],indexAdj],
				FieldGenerators[GetFieldLabel[field], group, {indexAdj, First[generatorIndices], Last[generatorIndices]}],
				RecursiveExpandCD[more\[Mu],group,field,\[Lambda]]
			]
			,
			CD[\[Mu]1,RecursiveExpandCD[more\[Mu],group,initField,\[Lambda]]]
		]
		,
		Plus[
			Times[-1,
				\[Lambda],I,
				(*$GaugeGroups[group][Coupling][],*)
				$GaugeGroups[group, Field][\[Mu]1[[1,1]],indexAdj],
				FieldGenerators[GetFieldLabel[field], group, {indexAdj, First[generatorIndices], Last[generatorIndices]}],
				field
			]
			,
			CD[\[Mu]1,initField]
		]
	]
]


(* ::Text:: *)
(*Abelian groups*)


RecursiveExpandAbelianCD[indices_List, group_, initField_, \[Lambda]_:1] := Module[
	{
		\[Mu]1,
		more\[Mu],
		field= initField
	}
	,
	(* split the Lorentz (or Spatial) indices *)
	{\[Mu]1,more\[Mu]} = TakeDrop[indices,1];
	If[!CovDerivativeIndexQ[First[\[Mu]1]], Message[RecursiveExpandAbelianCD::partialdev]; Abort[]];
	(* recursive formula for building up the covariant derivative expansion *)
	If[Length[more\[Mu]]>=1,
		Plus[
			Times[-1, (*crosschecked*)
				\[Lambda],I,
				(*$GaugeGroups[group, Coupling][],*)
				$GaugeGroups[group, Field][\[Mu]1[[1,1]]],
				FieldGenerators[GetFieldLabel[field], group],
				RecursiveExpandAbelianCD[more\[Mu],group,field,\[Lambda]]
			]
			,
			CD[\[Mu]1, RecursiveExpandAbelianCD[more\[Mu],group,initField,\[Lambda]]]
		]
		,
		Plus[
			Times[-1, (*crosschecked*)
				\[Lambda],I,
				(*$GaugeGroups[group, Coupling][],*)
				$GaugeGroups[group, Field][\[Mu]1[[1,1]]],
				FieldGenerators[GetFieldLabel[field], group],
				field
			]
			,
			CD[\[Mu]1, initField]
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*Auxiliary definitions*)


GetFieldLabel[field:(_Field | _FieldStrength)]:=Module[
	{
		label = First[field]
	}
	,
	(* treat Background vector fields *)
	If[Head[label]===BkgField,
		label = First[label]
	];
	(* treat fields expanded in the EFT power counting *)
	If[Head[label]===List,
		label = First[label]
	];
	Return[label]
]


FD::multiplegroupindices="The field `1` has several indices belonging to the group `2`.";


(* Function finding the index of a field that belongs to a specific group, returning that index and the field with the index replaced by a given label *)
ChangeFieldIndex[field_, group_, label_]:=Module[
	{
		(* find all indices of the field *)
		indices = Cases[field,_Index,All],
		groupIndex={}
	}
	,
	(* Find all indices belongin to the given group *)
	Do[
		If[MatchQ[GroupFromRep@Last[ind],group],
			AppendTo[groupIndex,ind]
		]
		,
		{ind,indices}
	];

	(* check that there is only one index *)
	If[Length[groupIndex]>1,
		Message[FD::multiplegroupindices,field,group];
		Abort[]
		,
		groupIndex = First[groupIndex]
	];

	(* return the field with replaced index and the original and new index labels *)
	Return[{
		field/.groupIndex->Index[label,Last[groupIndex]],
		{First[groupIndex],label}
	}]
]


(* ::Subsubsection::Closed:: *)
(*Shifting the field strength by their background configuration*)


(* ::Text:: *)
(*counter is a power-counting variable, default=1 -> no power counting*)


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, {}, lorentzIndices_List], counter_:1] :=
	RelabelIndices[Plus[
		FieldStrength[label, {\[Mu],\[Nu]}, {}, lorentzIndices],
		counter*Field[label, Vector[\[Nu]], {}, Append[lorentzIndices,\[Mu]]],
		-counter*Field[label, Vector[\[Mu]], {}, Append[lorentzIndices,\[Nu]]]
	],Unique->True]


(* ::Text:: *)
(*Nb. this is not a conclusive way of telling if a field strength tensor is for a gauge field or not*)


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, gaugeIndices:{Index[a_,group_[adj]]}, lorentzIndices_List], counter_:1] :=
	Module[{b,c},
	RelabelIndices[Plus[
		FieldStrength[label, {\[Mu],\[Nu]}, gaugeIndices, lorentzIndices],
		counter*Field[label, Vector[\[Nu]], gaugeIndices, Append[lorentzIndices,\[Mu]]],
		-counter*Field[label, Vector[\[Mu]], gaugeIndices, Append[lorentzIndices,\[Nu]]],
		counter^2(**$GaugeGroups[group, Coupling][]*)* CG[fStruct[group],{Index[a,group[adj]],Index[b,group[adj]],Index[c,group[adj]]}]*
		CD[lorentzIndices,Field[label, Vector[\[Mu]], {Index[b,group[adj]]},{}]*Field[label, Vector[\[Nu]], {Index[c,group[adj]]}, {}]]
	],Unique->True]
	]


ShiftFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, gaugeIndices:{Index[a_,group_[adj]],Index[b_,group_[adj]]}, lorentzIndices_List], counter_:1] :=
	Module[{c},
	 RelabelIndices[I CG[fStruct[group],{Index[a,group[adj]],Index[c,group[adj]],Index[b,group[adj]]}]*
	ShiftFS[FieldStrength[label, {\[Mu],\[Nu]}, {Index[c,group[adj]]}, lorentzIndices], counter],Unique->True]
	]


(* ::Section:: *)
(*FuncD \[LongDash] internal functionality*)


Options[FuncD] = {Grassmann->True};


(* ::Subsection:: *)
(*Differentiation rules*)


(* ::Subsubsection::Closed:: *)
(*Sum rule*)


FuncD[Plus[x1_,x2__],y_,OptionsPattern[]] :=
	Plus[
		FuncD[x1,y,Grassmann->OptionValue[Grassmann]],
		FuncD[Plus[x2],y,Grassmann->OptionValue[Grassmann]]
	]


(* ::Subsubsection::Closed:: *)
(*Product rules*)


(* ::Text:: *)
(*Times[]*)


FuncD[Times[x1_,x2__],y_,OptionsPattern[]] :=
	Plus[
		Times[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		Times[x1,FuncD[Times[x2],y,Grassmann->OptionValue[Grassmann]]]
	]


(* ::Text:: *)
(*FuncNCM[]*)


FuncD[FuncNCM[x1_,x2___],y_,OptionsPattern[]] :=
	(*apply the product rule \[LongDash] no Grassmann sign needed: FuncNCM never contains ghost fields*)
	Plus[
		FuncNCM[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		FuncNCM[x1,FuncD[FuncNCM[x2],y,Grassmann->OptionValue[Grassmann]]]
	] /; FreeQ[x1,Plus]


(* ::Text:: *)
(*NCM[] and GrassmannProduct[] \[LongDash] unified rule with head-specific Grassmann field types*)


(*NCM tracks Grassmann signs for Fermions; GrassmannProduct tracks them for Fermions and Ghosts*)
GrassmannSignTypes[NCM]= Fermion;
GrassmannSignTypes[GrassmannProduct]= Fermion|Ghost|AntiGhost;


FuncD[(h:NCM|GrassmannProduct)[x1_,x2__],y_,OptionsPattern[]] := Module[
	{
		fieldType,
		sign = 1,
		gTypes = GrassmannSignTypes[h]
	}
	,
	(* when Grassmann minus sign should be considered *)
	If[OptionValue[Grassmann],
		(*extract field type*)
		fieldType = Cases[y,Field[_,typ_,___]:>typ,All];
		If[MatchQ[fieldType,{gTypes}],
			(*get the number of anti-commuting fields in x1*)
			sign = (-1)^(Plus@@Cases[x1,Field[_,gTypes,___]->1,All]);
		];
	];

	(*apply the product rule*)
	Plus[
		h[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		sign * h[x1,FuncD[h[x2],y,Grassmann->OptionValue[Grassmann]]]
	]
] /; FreeQ[x1,Plus]


FuncD[(h:NCM|GrassmannProduct)[x_],y_,OptionsPattern[]]:= h@FuncD[x,y,Grassmann->OptionValue[Grassmann]]


(* ::Text:: *)
(*NCProduct[]*)


FuncD[NCProduct[x1_,x2__],y_,OptionsPattern[]] := Module[
	{
		fieldType,
		sign = 1
	}
	,
	(*when Grassmann minus sign should be considered*)
	If[OptionValue[Grassmann],
		(*extract field type*)
		fieldType = Cases[y,Field[_,typ_,___]:>typ,All];
		If[MatchQ[fieldType,{Fermion}],
			(*Get the number of anti-commuting fields in x1*)
			sign = (-1)^(Plus@@Cases[x1,Field[_,Fermion,___]->1,All]);
		];
	];

	(*apply the product rule*)
	Plus[
		NCProduct[FuncD[x1,y,Grassmann->OptionValue[Grassmann]],x2],
		sign * NCProduct[x1,FuncD[NCProduct[x2],y,Grassmann->OptionValue[Grassmann]]]
	]
] /; FreeQ[x1,Plus]


FuncD[NCProduct[x_], y_, OptionsPattern[]]:= NCProduct[FuncD[x, y, Grassmann->OptionValue[Grassmann]]]


(* ::Subsubsection::Closed:: *)
(*Power rule*)


FuncD[Power[x_,pow_Integer],y_,OptionsPattern[]]:=
	pow*Power[x,pow-1]*FuncD[x,y,Grassmann->OptionValue[Grassmann]] /; pow>=2


(* ::Subsection:: *)
(*Evaluation of derivatives*)


(* ::Subsubsection::Closed:: *)
(*Vanishing derivatives*)


FuncD[x_,_,OptionsPattern[]] := 0 /; FreeQ[x, Field] (* FieldStrength does not contain vector field fluctuations after the shift in FD *)


FuncD[Bar@Field[___], f:Field[___],OptionsPattern[]] := 0
FuncD[Field[___], Bar[f:Field[___]],OptionsPattern[]] := 0


FuncD[Field[_,_,_,\[Mu]_List], f:Field[_,_,_,\[Nu]_List], OptionsPattern[]] := 0 /; (Length[\[Mu]]!=Length[\[Nu]])
FuncD[Bar@Field[_,_,_,\[Mu]_List], Bar[f:Field[_,_,_,\[Nu]_List]], OptionsPattern[]] := 0 /; (Length[\[Mu]]!=Length[\[Nu]])


FuncD[Field[label1_,__], f:Field[label2_,__], OptionsPattern[]] := 0 /; (label1=!=label2)
FuncD[Bar@Field[label1_,__], Bar[f:Field[label2_,__]], OptionsPattern[]] := 0 /; (label1=!=label2)


(* ::Subsubsection::Closed:: *)
(*Non-vanishing derivatives*)


FuncD::repeatedindices = "The indices `1` in argument x and the indices `2` in argument y of FuncD[x,y] must not have identical labels.";


FuncD::indexnumber = "Unmatched number of group indices found.";


FuncD::unmatchedindex = "Unmatched types of group indices found.";


FuncD[Bar@Field[label_,rest1__], Bar@Field[label_,rest2__], OptionsPattern[]] := FuncD[Field[label,rest1], Field[label,rest2], Grassmann->OptionValue[Grassmann]]


FuncD[Field[label_,type1_,index1_,derivs1_], Field[label_,type2_,index2_,derivs2_], OptionsPattern[]] :=
	Module[
		{ result = 1 }
		,
		(*check that the group indices match up*)
		If[Length[index1]!=Length[index2],
			Message[FuncD::indexnumber];
			Abort[]
			,
			If[(index1/.Index[_,rep_]:>rep)=!=(index2/.Index[_,rep_]:>rep),
				Message[FuncD::unmatchedindex];
				Abort[]
			]
		];

		(*check that there are no duplicate index labels *)
		If[Or@@(!FreeQ[index1,#]&/@index2),
			Message[FuncD::repeatedindices, index1, index2];
			Abort[]
		];
		If[Or@@(!FreeQ[derivs1,#]&/@derivs2),
			Message[FuncD::repeatedindices, derivs1, derivs2];
			Abort[]
		];

		(* add metric for vectors *)
		If[MatchQ[{type1,type2},{Vector[Index[_,Lorentz|Spatial]],Vector[Index[_,Lorentz|Spatial]]}],
				If[type1===type2,
					Message[FuncD::repeatedindices, First@type1, First@type2]; Abort[],
					result *= Metric[First@type1, First@type2]
				]
		];

		(*include Delta[_,_] for all group indices*)
		Do[
			result *= Delta@@pair
			,
			{pair,Transpose[{index1,index2}]}
		];

		(*include Metric[_,_] for all Lorentz (or Spatial) indices of the covariant derivatives*)
		Do[
			result *= Metric@@pair
			,
			{pair,Transpose[{derivs1,derivs2}]}
		];

		Return[result]
	] /; (Length[derivs1]==Length[derivs2])


(*FuncD[Field[label_,type1_,index1_,derivs1_], Field[label_,type2_,index2_,derivs2_], OptionsPattern[]] :=
	Module[
		{ result = 1 }
		,
		(*check that the group indices match up*)
		If[Length[index1]!=Length[index2],
			Message[FuncD::indexnumber];
			Abort[]
			,
			If[(index1/.Index[_,rep_]:>rep)=!=(index2/.Index[_,rep_]:>rep),
				Message[FuncD::unmatchedindex];
				Abort[]
			]
		];

		(*check that there are no duplicate index labels *)
		If[Or@@(!FreeQ[index1,#]&/@index2),
			Message[FuncD::repeatedindices, index1, index2];
			Abort[]
		];
		If[Or@@(!FreeQ[derivs1,#]&/@derivs2),
			Message[FuncD::repeatedindices, derivs1, derivs2];
			Abort[]
		];

		(*multiply by metric tensor for vectors*)
		If[{type1,type2}==={Vector[Index[_,Lorentz]],Vector[Index[_,Lorentz]]},
				If[type1===type2,
					Message[FuncD::repeatedindices, First@type1, First@type2]; Abort[],
					result *= Metric[First@type1, First@type2]
				]
		];

		(*include Delta[_,_] for all group indices*)
		Do[
			result *= Delta@@pair
			,
			{pair,Transpose[{index1,index2}]}
		];

		(*include Metric[_,_] for all Lorentz indices of the covariant derivatives*)
		Do[
			result *= Metric@@pair
			,
			{pair,Transpose[{derivs1,derivs2}]}
		];

		Return[result]
	] /; (Length[derivs1]==Length[derivs2])*)


(* ::Subsection:: *)
(*FuncD with transposed and charge-conjugated fields *)


(* ::Subsubsection::Closed:: *)
(*FuncD w.r.t. Field[...]*)


FuncD[x_, Transp[y_], OptionsPattern[]] := (TransposeThisSpinChain \[CenterDot] FuncD[x, y, Grassmann->OptionValue[Grassmann]])/; FreeQ[x,Plus]


FuncD[Transp[x:(Field[___] | Bar@Field[___])], y:Except[_NCM] , OptionsPattern[]] := (TransposeThisSpinChain \[CenterDot] FuncD[x, y, Grassmann->OptionValue[Grassmann]])


(* ::Text:: *)
(*If y is a NCM[...] then y is the charge conjugate of some field for which the rules below should be used instead, since transposition generates minus signs due to the presence of GammaCC.*)


(* ::Subsubsection::Closed:: *)
(*FuncD  w.r.t. CC@Field[...]*)


FuncD[
	Bar[x_Field],
	NCM[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := NCM[
	DiracProduct[GammaCC],
	TransposeThisSpinChain,
	FuncD[Bar[x], Bar[y], Grassmann->OptionValue[Grassmann]]
]


FuncD[
	Transp@Bar[x_Field],
	NCM[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := (-1)*NCM[
	DiracProduct[GammaCC],
	FuncD[Bar[x], Bar[y], Grassmann->OptionValue[Grassmann]]
]


FuncD[
	x:(Field[___] | Transp@Field[___]),
	NCM[DiracProduct[GammaCC], Transp@Bar[y_Field]],
	OptionsPattern[]
] := 0


(* ::Subsubsection::Closed:: *)
(*FuncD  w.r.t. Bar@CC@Field[...]*)


FuncD[
	x_Field,
	NCM[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := NCM[
	DiracProduct[GammaCC],
	TransposeThisSpinChain,
	FuncD[x, y, Grassmann->OptionValue[Grassmann]]
]


FuncD[
	Transp[x_Field],
	NCM[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := (-1) * NCM[
	DiracProduct[GammaCC],
	FuncD[x, y, Grassmann->OptionValue[Grassmann]]
]


FuncD[
	x:(Bar@Field[___] | Transp@Bar@Field[___]),
	NCM[Transp[y_Field], DiracProduct[GammaCC]],
	OptionsPattern[]
] := 0


(* ::Section:: *)
(*FuncDSimplify*)


Transp[TransposeThisSpinChain] := TransposeThisSpinChain


(* ::Text:: *)
(*Simplifies the expression returned by FuncD by transposing spin-chains where necessary*)


Options@ FuncDSimplify= {Contract-> True}


FuncDSimplify[expr_, OptionsPattern[]] := Module[{solution=expr},
	(* avoid unnecessary transpositions *)
	solution=solution//.HoldPattern@NCM[A___,TransposeThisSpinChain,B___,TransposeThisSpinChain,C___]:>NCM[A,B,C]/;FreeQ[B,TransposeThisSpinChain];
	solution=solution//.HoldPattern@NCM[A___,TransposeThisSpinChain,TransposeThisSpinChain,C___]:>NCM[A,C];
	(* transpose all necessary spin-chains *)
	solution=solution/.HoldPattern@NCM[A___,TransposeThisSpinChain,B___]:>Transp@NCM[A,B];
	(* contract if required *)
	If[OptionValue[Contract],
		Contract[solution],
		solution
	]
]


(* ::Section:: *)
(*Variational derivative \[LongDash] VarD*)


Options[VarDraw] = {Grassmann->True, EFTOrder -> 6};


(* ::Subsubsection::Closed:: *)
(*Check Arguments*)


OptionMessage[Grassmann, func:VarDraw, val_] := Message[General::optexpectsval, Grassmann, func, val, "boolean (True or False)"];


VarDraw::invalidargument =
"The second argument `1` of VarDraw is invalid. It must be either of the following forms:
Field[_,_,_,{}], Bar[Field[_,_,_,{}]], Transp[Field[_,_,_,{}]], Transp[Bar[Field[_,_,_,{}]]], CConj[Field[_,_,_,{}]], Bar[CConj[Field[_,_,_,{}]]].";


SecondArgVarD = Alternatives[
		Field[_,_,_,{}],
		Bar@Field[_,_,_,{}],
		Transp@Field[_,_,_,{}],
		Transp@Bar@Field[_,_,_,{}],
		NCM[DiracProduct[GammaCC], Transp@Bar@Field[_,_,_,{}]],
		NCM[Transp@Field[___], DiracProduct[GammaCC]]
	];


VarDraw[
	_,
	y : Except@SecondArgVarD,
	OptionsPattern[]
] := Message[VarDraw::invalidargument,y]


(* ::Subsection:: *)
(*The total functional derivative*)


(* ::Subsubsection::Closed:: *)
(*Single application of variational derivative*)


CovariantLoop[lagrangian_, fields_List, opts:OptionsPattern[]]? OptionsCheck:=
CovariantLoop[lagrangian, fields, opts]= Module[{},0];


VarDraw[L_, x : SecondArgVarD, opts:OptionsPattern[]]:=
VarDraw[L, x, opts]= Module[
  	{
   		n,
   		arg,
   		result,
   		indices = {}
   	},

  	If[Head@(Cases[x, Field[___], All][[1,2]])===Vector,arg=Expand[L/$\[Epsilon]FD];, arg=L;];

  	(*Add FD w.r.t. Field[...] to result*)
  	result = FD[arg, x, Grassmann -> OptionValue[Grassmann]];

  	(*Add FD w.r.t. CD[Field[...]] to result, included number of CD[] depends on mass dimension*)
  	Do[
   		(*create unique lorentz indices for covariant derivatives*)
   		indices = Table[Unique[],{i, n}];
   		result += ((-1)^n) * FuncDSimplify@CD[
       			Reverse@indices, (* indices need to be reversed due to IBP *)
       			FD[arg, CD[indices, x], Grassmann -> OptionValue[Grassmann]]
       		];
   		,
   		{n, OptionValue[EFTOrder] - 2}
   	];

  	result
  ]


(* ::Subsubsection::Closed:: *)
(*Variational derivatives*)


Options[VarD] = {EFTOrder -> 6};
Options[VarD1] = {EFTOrder -> 6};


(* ::Text:: *)
(*Master function removing all terms which trivially vanish under differentiation*)


VarD[lag_, fs__, opt:OptionsPattern[]]:= Module[{gr, terms, labels},
	labels= Cases[{fs}, Field[l_, __]-> l, All];
	terms= TermsToList@ BetterExpand@ lag;

	(* Remove terms which trivially vanish under differentiation *)
	Do[
		terms= If[MemberQ[$GaugeGroups, KeyValuePattern[Field-> lab]],
			(* gauge field *)
			gr= First@ GaugeGroupByProperty[Field -> lab];
			Select[terms, NotTrivialWRTGaugeFieldQ[gr]]
		,
			(* matter field *)
			Cases[terms, a_/;! FreeQ[a, Field[lab, __]] ]
		];
	, {lab, labels}];

	terms= ExpandVectorFluctuations[Total@ terms, Sequence@@ DeleteDuplicatesBy[{fs},
		FirstCase[{#}, Field[l_, __] :> l, Nothing, All]& ] ];
	VarD1[terms, fs, opt]
];


(* ::Text:: *)
(*Checks if the term is trivially zero under differentiation w.r.t. the field gauge group gr*)


NotTrivialWRTGaugeFieldQ[gr_][term_] := !FreeQ[term, f : (Field[__, {__}] | _FieldStrength) /; FieldTransformsUnderGaugeGroupQ[f, gr]] ||
	!FreeQ[term, (Field|FieldStrength)[$GaugeGroups[gr, Field], __]];


(* ::Text:: *)
(*Variational derivative w.r.t. 1 field*)


VarD1[lag_, f1 : SecondArgVarD, opt:OptionsPattern[]] :=
	RelabelIndices@VarDraw[lag, f1,opt]/.$\[Epsilon]FD->0 /. BkgField[l_]:>l


(* ::Text:: *)
(*Variational derivative w.r.t. 2 fields*)


VarD1[lag_, f1 : SecondArgVarD, f2 : SecondArgVarD, opt:OptionsPattern[]] := Module[{res},
	res = RelabelIndices@VarDraw[Expand[VarDraw[lag, f2,opt] * FuncNCM[OpenCD[{}]]], f1, opt]/.
		OpenCD[{}]-> 1/. $\[Epsilon]FD-> 0/. BkgField[l_]:>l;
	res = res/.Times->FuncNCM;
	res = res/.{FuncNCM[a___,ocd_OpenCD,b___]:>FuncNCM[a,b,ocd]}
]


(* ::Section:: *)
(*Equations of Motion & fluctuation operator*)


(* ::Subsection:: *)
(*EOM*)


Options[DeriveEOM] = {EFTOrder->6};


(* ::Subsubsection::Closed:: *)
(*EoM for one field*)


DeriveEOM[lagrangian_, field:(Field[label_,___,{}] | Bar@Field[label_,___,{}]), OptionsPattern[]] := Module[
	{
		eom,
		deriveField = Bar[field],
		fieldProps = GetFieldsUpdated@ label
	}
	,

	(* check for massive chiral fermions *)
	If[(fieldProps[Type] === Fermion) && (fieldProps[Chiral] =!= False) && (fieldProps[Mass] =!= 0) && !fieldProps[SelfConjugate],
		deriveField = CConj@deriveField
	];

	(* derive EOM *)
	eom = VarD[lagrangian, deriveField, EFTOrder->OptionValue[EFTOrder]];

	Return[eom]
]


(* ::Subsubsection::Closed:: *)
(*EoM for list of fields*)


DeriveEOM[lagrangian_, fields_List, OptionsPattern[]] := (DeriveEOM[lagrangian,#,EFTOrder->OptionValue[EFTOrder]]&/@fields)


(* ::Subsubsection::Closed:: *)
(*Errors*)


DeriveEOM::invalarg = "The argument `1` is invalid; only Field[...] or Bar[Field[...]] are allowed.";


DeriveEOM[_,arg:Except[List[___]|Field[___]|Bar@Field[___]], OptionsPattern[]] := Message[DeriveEOM::invalarg,arg]


(* ::Subsection:: *)
(*Fluctuation operator*)


Options[FluctuationOperator] = {EFTOrder -> 6};


(* ::Text:: *)
(*Acceptable input for functional derivatives *)


fdPattern= Alternatives[
		Field[_,_,_,{}],
		Bar@Field[_,_,_,{}],
		Transp@Field[_,_,_,{}],
		Transp@Bar@Field[_,_,_,{}],
		NCM[DiracProduct[GammaCC], Transp@Bar@Field[_,_,_,{}]],
		NCM[Transp@Field[___], DiracProduct[GammaCC]]
	];


FluctuationOperator::invldField= "`1` is not a valid field variable.";


FluctuationOperator[lag_, field1_, field2_, OptionsPattern[]] := Module[
		{f1, f2, eom, xTerm, f, sign, lagrangian= NCProduct[lag]},
	(*Remove projection operators from functional derivative field (chiral fields always have explicit projectors in the Lagrangian).*)
	{f1, f2}= {field1, field2}/. {DiracProduct[a___, (Transp@ _Proj| _Proj)]:> DiracProduct@ a}/.
		NCM@ x_-> x;

	Do[
		If[! MatchQ[f, fdPattern],
			Message[FluctuationOperator::invldField, f];
			Abort[];
		];
		,
		{f, {f1, f2}}
	];

	(*
	if field2 is fermionic or ghost add a minus sign:
	since VarD[\[ScriptCapitalL],field2] has an odd number of fermions and the fluctuation of field2 has to be anticommuted through the entire expression.
	 *)
	If[FreeQ[field2, Field[_, Fermion|Ghost|AntiGhost, ___], All],
		sign= 1,
		sign= -1
	];

	(* apply functional derivatives *)
	xTerm= sign * VarD[lagrangian,f1,f2, EFTOrder->OptionValue[EFTOrder]];

	xTerm= SortNCProduct[xTerm]/.NCProduct->NCM;

	Return[xTerm]
]


(* ::Section:: *)
(*Auxiliary NCProduct*)


(* ::Subsection:: *)
(*Properties of NCProduct*)


NCProduct[a___, NCProduct[b__], c___]:= NCProduct[a, b, c]
NCProduct[] = 1;
NCProduct@NCProduct[x___]:=NCProduct[x];


(* ::Subsubsection::Closed:: *)
(*General properties*)


(* ::Text:: *)
(*NCProduct contractions*)


(* Extracting commuting objects *)
NCProduct[a___, b_?CommutativeQ, c___]/; FreeQ[b,_NCM]:= b * NCProduct[a, c]
NCProduct[a___, b_?CommutativeQ * x_, c___]/; FreeQ[b,_NCM]:= b * NCProduct[a, x, c]
(* Distributivity *)
NCProduct[a___, b_Plus, c___]:= NCProduct[a, #, c] & /@ b


NCProduct[a___, Times[b_NCM, c___], d___]:= NCProduct[a,b,Times[c],d]


SortNCProduct[expr_]:= Module[
	{
		rules={
			NCProduct[a___,b_NCM,c___]/;ClosedSpinChainQ[b] :> NCProduct[a,c,b],
			NCProduct[a___,x_NCM,y_NCM,b___] /;(ROpenSpinChainQ[x]&&LOpenSpinChainQ[y]) :> (-1)^(GrassmannCount[x]*GrassmannCount[y]) * NCProduct[a,y,x,b]
		}
	}
	,
	expr/.rules
]


(* counts the grassmann parameters in arg *)
GrassmannCount[arg_]:= Count[arg, Fermion|Ghost|AntiGhost, All]
