(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`EFTCounting`*)


(* ::Subtitle:: *)
(*Various methods for handling the EFT expansion and dimension counting *)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping:*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["SeriesEFT"]
PackageExport["EFTOrder"]


PackageExport["HeavyFieldDimension"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["OperatorDimension"]
PackageScope["TypeDim"]


PackageScope["GetMaxEFTOrder"]
PackageScope["GetMinEFTOrder"]
PackageScope["TruncateOperator"]


PackageScope["DenominatorOperatorDimension"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


SeriesEFT::usage= "SeriesEFT[expr, EFTOrder -> 6] truncates the expression expr in the EFT power counting at the mass-dimension specified by the option EFTOrder (6 by default). The EFTOrder can be given as n, to return the expansion up to mass dimension n, or {n} to return exclusively the terms of order exactly n.";
EFTOrder::usage= "EFTOrder is an option for various routines specifying the targeted order in power-counting, as well as defining the power-counting of couplings.";


HeavyFieldDimension::usage= "HeavyFieldDimension is an option of SeriesEFT that determines whether heavy fields are counted by their canonical mass dimension or an extra suppression factor is included to account for EFT counting when the heavy field is replaced in terms of light.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


OperatorDimension::usage = "OperatorDimension[op] returns the mass-dimension of the operator op.";


(* ::Chapter:: *)
(*Private:*)


(* ::Subsection:: *)
(*Dimensionality of an operator*)


(* ::Subsubsection::Closed:: *)
(*Scoped functions*)


(* ::Text:: *)
(*Gives the (minimal) dimension of a (list/sum of) operator(s). *)
(*	<heavyDim> = True indicates that the heavy field gets an extra dimension (minimal from EOM solution).*)
(*	<heavyDim> = False indicates that the canonical dimension is used for heavy fields.*)


OperatorDimension[expr_, heavyDim_:True]:= Block[{},
	$heavyFieldDim= heavyDim; (*Flag that's only used in FieldDimension*)
	OperatorDimensionInternal@ expr
];


(* ::Text:: *)
(*Returns the maximum/minimum operator dimension among all terms in a Lagrangian. The optional second argument <heavyDim> is passed to OperatorDimension.*)


GetMaxEFTOrder[lag_, heavyDim_:True]:= Max[OperatorDimension[#, heavyDim]& /@ TermsToList@ lag]

GetMinEFTOrder[lag_, heavyDim_:True]:= Min[OperatorDimension[#, heavyDim]& /@ TermsToList@ lag]


(* ::Subsubsection::Closed:: *)
(*Internal helper functions*)


(* ::Text:: *)
(*Internal function threading over Operator, List and Plus*)


OperatorDimensionInternal[x_Operator] := OperatorDimensionInternal[OperatorToNormalForm@ x];
OperatorDimensionInternal[c_ x_Operator] := OperatorDimensionInternal[c] + OperatorDimensionInternal[OperatorToNormalForm@ x];


OperatorDimensionInternal[AtomicOp[{fieldList_,devs_},___]] :=
	Total[(TypeDim@ $FieldAssociation[#][Type] + Boole@ $FieldAssociation[#][Heavy])&/@
	(fieldList /. Conj->Identity)] + devs

OperatorDimensionInternal[CompOp[_,{fieldList_,devs_},___]] :=
	Total[(TypeDim@ $FieldAssociation[#][Type] + Boole@ $FieldAssociation[#][Heavy])&/@
	(fieldList /. Conj->Identity)] + devs


(* ::Text:: *)
(*How OperatorDimensionInternal calls recursively into an expression*)


OperatorDimensionInternal[expr_Plus]:= Min[OperatorDimensionInternal/@ (List@@ expr)];
OperatorDimensionInternal[expr_List]:= Min[OperatorDimensionInternal/@ expr];


OperatorDimensionInternal@ (Bar|Transp)@ expr_:= OperatorDimensionInternal@ expr;
OperatorDimensionInternal@ Power[expr_, n_]:= n* OperatorDimensionInternal@ expr;
(*All kinds of products*)
OperatorDimensionInternal@ expr:_Times|_NCM|_GrassmannProduct|_FuncNCM:= 
	Total[OperatorDimensionInternal/@ List@@ expr];


(* ::Text:: *)
(*Rules for specific objects*)


OperatorDimensionInternal[0]= 100;
OperatorDimensionInternal@ Coupling[__, n_]:= n;
OperatorDimensionInternal@ Matchete`FeynmanRules`PackagePrivate`PDint[\[Mu]_, expr_]:= 
	Length@ \[Mu]+ OperatorDimensionInternal@ expr;
OperatorDimensionInternal[f:_Field|_FieldStrength]:= FieldDimension@ f;
OperatorDimensionInternal@ _= 0


(* ::Text:: *)
(*Counting rule for IR fields*)


FieldDimension[Field[f:Except[List[___]], type_, _, derivs_List]]:=
	Length[derivs] + TypeDim[type] + Boole[$heavyFieldDim && $FieldAssociation[f][Heavy]]; (* heavy fields have at least one suppression factor*)


(* ::Text:: *)
(*Counting rule for expanded UV fields*)


FieldDimension[Field[{_,n_,___}, type_, _, derivs_List]]:= n + Length[derivs] + TypeDim[type];


(* ::Text:: *)
(*Counting rule for FieldStrengthTensors fields*)


FieldDimension[FieldStrength[f_, _, _, derivs_List]]:= Length[derivs] + 2 + If[$heavyFieldDim && GetFieldsUpdated[f, Heavy], 2, 0];


(* ::Subsubsection::Closed:: *)
(*TypeDim*)


(* ::Text:: *)
(*Dimension for field types*)


TypeDim[Scalar|Vector[_]|Ghost|AntiGhost]= 1;
TypeDim[Fermion]= 3/2;


TypeDim[Field[_,type_,___]]:= TypeDim@ type;


(* ::Subsubsection::Closed:: *)
(*DenominatorOperatorDimension*)


(* ::Text:: *)
(*Is this actually necessary, or can we replace it by the updated version of OperatorDimension?*)


(*
DenominatorOperatorDimension[arg_] := Module[{expr,$\[Lambda]},
	expr= arg/. c:Coupling[_,_,pow_]:> c * $\[Lambda]^pow;
	expr= expr/. f_Field:> f * $\[Lambda]^FieldDimension[f];
	expr= expr/. f_FieldStrength:> f * $\[Lambda]^FieldDimension[f];
	expr= expr/. d_Matchete`FeynmanRules`PackagePrivate`PDint:> Power[$\[Lambda], Length@First[d]] * d;
	Exponent[expr, $\[Lambda]]
]
*)


(* ::Subsection:: *)
(*SeriesEFT*)


(* ::Subsubsection::Closed:: *)
(*SeriesEFT*)


(* ::Text:: *)
(*SeriesEFT truncates an expression at the desired EFT-order. Supports exclusive with {<order>} option*)


SeriesEFT::notarationale = "The order for the EFT series expansion specified `1` is neither an integer nor a rationale.";


Options[SeriesEFT]= {
		EFTOrder-> 6,
		HeavyFieldDimension-> True,
		ExpandDenominator-> True
	};


SeriesEFT[arg_, OptionsPattern[]]:= Module[
	{
		order= OptionValue[EFTOrder],
		expr= LagrangianExpand[arg],
		heavyDim= OptionValue@ HeavyFieldDimension,
		result,
		TruncateOp, TruncateOpEx
	},
	(* determine whether a Taylor expansion should be applied or not *)
	If[OptionValue@ ExpandDenominator,
		(* Taylor expand denominators if required *)
		TruncateOp   = TaylorTruncateOperator;
		TruncateOpEx = TaylorTruncateOperatorExact
		,
		(* Use minimum EFT dimension for each denominator, not performing and Taylor expansion *)
		TruncateOp   = TruncateOperator;
		(*TruncateOpEx = TruncateOperatorExact*)
		(* exact truncation can only be achieved with a Taylor expansion *)
		TruncateOpEx = TaylorTruncateOperatorExact
	];
	
	(* determine whether EFTOrder is exclusive or inclusive *)
	Switch[order,
		(* select terms up to the given order *)
		_Integer | _Rational,
			If[Head[expr]===Plus,
				result = TruncateOp[#,order,heavyDim]&/@ expr
				,
				result = TruncateOp[expr,order,heavyDim]
			],
		(* select terms of the exact order *)
		{_Integer | _Rational},
			order = First[order];
			If[Head[expr]===Plus,
				result = TruncateOpEx[#,order,heavyDim]&/@ expr
				,
				result = TruncateOpEx[expr,order,heavyDim]
			],
		_,
			Message[SeriesEFT::notarationale, order];
			Abort[]
	];
	result
]


(* ::Subsubsection::Closed:: *)
(*Operator truncation functions*)


TruncateOperator[op_, dim_, heavyDim_:True]:= If[OperatorDimension[op, heavyDim] > dim, 0, op]


TruncateOperatorExact[op_, dim_, heavyDim_:True]:= If[OperatorDimension[op, heavyDim] === dim, op, 0]


(* ::Subsubsection::Closed:: *)
(*Operator truncation functions with Taylor expansions*)


(* ::Text:: *)
(*This function performs a SeriesExpansion if required and it should be applied only to individual terms (not on sums).*)


TaylorTruncateOperator[arg_, ord_, heavyDim_]:=Module[
	{
		expr,
		$\[Lambda],
		denominator= Denominator[arg]
	}
	,
	If[FreeQ[denominator, _Plus, All],
		(* no Taylor expansion required *)
		TruncateOperator[arg, ord, heavyDim]
		,
		(* Taylor expansion required *)
		$heavyFieldDim= heavyDim;
		expr= arg/. op_Operator-> OperatorToNormalForm[op];
		expr= expr/. c:Coupling[_,_,pow_]:> c * $\[Lambda]^pow;
		expr= expr/. f_Field:> f * $\[Lambda]^FieldDimension[f];
		expr= expr/. f_FieldStrength:> f * $\[Lambda]^FieldDimension[f];
		expr= expr/. d_Matchete`FeynmanRules`PackagePrivate`PDint:> Power[$\[Lambda], Length@First[d]] * d;
		BetterSeries[expr, {$\[Lambda], 0, ord}]/. $\[Lambda]->1
	]
]


TaylorTruncateOperatorExact[arg_, ord_, heavyDim_]:=Module[
	{
		expr,
		$\[Lambda],
		denominator= Denominator[arg]
	}
	,
	If[FreeQ[denominator, _Plus, All],
		TruncateOperatorExact[arg, ord, heavyDim]
		,
		$heavyFieldDim= heavyDim;
		expr= arg/. c:Coupling[_,_,pow_]:> c * $\[Lambda]^pow;
		expr= expr/. f_Field:> f * $\[Lambda]^FieldDimension[f];
		expr= expr/. f_FieldStrength:> f * $\[Lambda]^FieldDimension[f];
		expr= expr/. d_Matchete`FeynmanRules`PackagePrivate`PDint:> Power[$\[Lambda], Length@First[d]] * d;
		SeriesCoefficient[expr, {$\[Lambda], 0, ord}]
	]
]
