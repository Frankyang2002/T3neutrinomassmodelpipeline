(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`TreeLevelMatching`*)


(* ::Subtitle:: *)
(*Sub package providing the tree-level matching.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["PrintEOM"]


PackageScope["DetermineEOMs"]


PackageScope["ReplaceHeavyEOM"]


PackageScope["IndexToPattern"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PrintEOM::usage = 
"PrintEOM[L] computes the equations of motion of the Lagrangian L and prints them using NiceForm. By default the EOM of all heavy fields are printed. 
PrintEOM[L, Field -> {l1, l2, ...}] prints only the EOM of a specific set of fields defined by their labels l1, l2, ...";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Expanding the EoMs in the EFT power counting*)


(* ::Subsubsection::Closed:: *)
(*Expands fields in the EFT power counting*)


(* ::Text:: *)
(*Expand field produces a series expansion of the form Field[f, ...]:> Field[{f, 1}, ...] + Field[{f, 2}, ...]+ ... for heavy fields up to the highest order allowed for on general grounds. *)


ExpandField::notype = "The Field type `1` is unknown.";


Options[ExpandField]={EFTOrder -> 6};


ExpandField[Field[label_,type_,indices_,derivs_],OptionsPattern[]]:=Module[
	{
		expansion,
		cd = Unique[],
		myPattern,
		maxOrder
	}
	,
	(* Determine the max EFT order up to which the fields have to be expanded *)
	(*
		From any tree-level or one-loop vertex there is at least two external fermions (dim 3), 
		two external vectors, or 1 scalar. This follows from conservation of Lorentz symmetry.
	*)
	maxOrder = OptionValue[EFTOrder]+ Switch[type,
		Fermion, 
			- 3,
		Vector[_],
			- 2,
		Scalar ,
			- 1,
		_, 
			Message[ExpandField::notype,type];
			Abort[]	
	];
	
	(* write the field expansion, where n\[GreaterEqual]1 is the number of supression factors (m/M)^n *)
	expansion = Sum[
		Field[{label,n},type,indices,cd]
		,
		{n, 1, maxOrder}
	];
	(* construct the replacement rules *)
	With[{res=expansion},
		Return[
			(Field[label,IndexToPattern[type],IndexToPattern[indices],myPattern[cd,_]] :> res)/.myPattern->Pattern
		]
	]
]


(* ::Subsubsection::Closed:: *)
(*Expand EOMs in the EFT power counting*)


Options[ExpandEOM]={EFTOrder -> 6};


(* ::Text:: *)
(*ExpandEOM expands the UV fields in the equations of motion in the EFT power counting*)


ExpandEOM[eoms_List, fieldsUV_List, OptionsPattern[]]:=Module[
	{
		expandedEOM,
		order = OptionValue[EFTOrder]
	}
	,
	(* plug in expanded fields (heavy FieldStrength tensors have to be substituted before) *)
	expandedEOM = Expand[eoms /. (ExpandField[#, EFTOrder->order]&/@fieldsUV)];
	Return[
		SeriesEFT[#, EFTOrder->order]& /@ expandedEOM
	]
]


(* ::Subsubsection::Closed:: *)
(*Expand heavy fields in operators*)


Options[ExpandOperators]={EFTOrder -> 6};


ExpandOperators[arg_, uvFields_, OptionsPattern[]] := Module[
	{
		expr = LargeSumExpand[arg],
		return
	}
	,
	(* expand operator up to desired dimension *)
	If[Head[expr]===Plus,
		Return[ExpandOp[#, uvFields, EFTOrder->OptionValue[EFTOrder]]& /@ expr]
		,
		Return[ExpandOp[expr, uvFields, EFTOrder->OptionValue[EFTOrder]]]
	];
]


Options[ExpandOp]={EFTOrder -> 6};


ExpandOp[arg_, uvFields_, OptionsPattern[]] := Module[
	{
		expr,
		order = OptionValue[EFTOrder],
		inclusiveOrder
	}
	,
	(* find the inclusive EFT order *)
	If[MatchQ[order, {_}],
		inclusiveOrder = First[order],
		inclusiveOrder = order
	];
	(* expand uv fields in EFT power counting *)
	expr = Expand[arg /. (ExpandField[#, EFTOrder->inclusiveOrder]&/@uvFields)];
	(* expand up to desired order *)
	expr = SeriesEFT[expr, EFTOrder->order];
	(* remove the power *)
	If[Head[expr]===Plus,
		expr = RemovePower[#]& /@ expr,
		expr = RemovePower[expr]
	];
	Return[expr]
]


(* ::Subsubsection::Closed:: *)
(*Substitute EOM*)


Options[ReplaceHeavyEOM]={EFTOrder -> 6};


ReplaceHeavyEOM[arg_, OptionsPattern[]] := Module[
	{
		uvFields = FindUvFields[arg],
		order = OptionValue[EFTOrder],
		expr = SubstituteHeavyFieldStrength[arg]
	}
	,
	expr = ExpandOperators[expr, uvFields, EFTOrder->order];
	expr = LargeSumExpand@ReleaseHold[expr /. Last[$currentEOMs]];
	
	Return[expr]
]


(* ::Subsection:: *)
(*Auxiliary functionalities*)


(* ::Subsubsection::Closed:: *)
(*Substitute heavy field-strength tensors*)


(* ::Text:: *)
(*Rewrite: Subscript[F^a, \[Mu]\[Nu]]=(Subscript[D, \[Mu]] Subscript[A, \[Nu]])^a-(Subscript[D, \[Nu]] Subscript[A, \[Mu]])^a*)


SubstituteHeavyFieldStrength[expr_] := Expand[expr/.(FieldStrength[l_,x___]:>SubstituteFS@FieldStrength[l,x]/;GetFieldsUpdated[l][Heavy])]


SubstituteFS[FieldStrength[label_, {\[Mu]_Index,\[Nu]_Index}, gaugeIndices_List, lorentzIndices_List]] :=
	Field[label, Vector[\[Nu]], gaugeIndices, Append[lorentzIndices,\[Mu]]] - Field[label, Vector[\[Mu]], gaugeIndices, Append[lorentzIndices,\[Nu]]]


(* ::Subsubsection::Closed:: *)
(*Make indices to patterns*)


(* ::Text:: *)
(*Replace every Index[a,rep] with Index[a_,rep]*)


IndexToPattern[expr_]:=Module[
	{rules, myPattern}
	,
	rules = (Index[a_,rep_]:>Index[myPattern[a,_],rep]);
	Return[(expr/.rules)/.myPattern->Pattern]
]


(* ::Subsubsection::Closed:: *)
(*Fast expansion of large sums*)


(* This might be much faster than Mathematica's built in Expand[...] when dealing with large sums *)
LargeSumExpand[sum_Plus] := Expand[#]& /@ sum
LargeSumExpand[expr:Except[_Plus]] := Expand[expr]


(* ::Subsubsection::Closed:: *)
(*Finds the UV fields in the given expression*)


FindUvFields::error = "The field `1` is not part of the Lagrangian.";


FindUvFields[lagrangian_]:=Module[
	{
		fieldAssociation = GetFieldsUpdated[],
		uvFields
	},
	(* find all field labels *)
	uvFields = DeleteDuplicates@Join[
		Cases[lagrangian, Field[label_,___]:>label, All],
		Cases[lagrangian, FieldStrength[label_,___]:>label, All]
	];
	(* select labels of all heavy fields *)
	uvFields = If[fieldAssociation[#][Heavy]===True,#,Nothing[]]&/@uvFields;
	(* pick the apropriate fields *)
	uvFields = FirstCase[lagrangian,Field[#,___],Message[FindUvFields::error,#]; Abort[],All]&/@uvFields;
	(* make the indices unique *)
	uvFields = uvFields /. Index[_,rep_]:>Index[Unique[],rep];
	Return[uvFields]
]


(* ::Section:: *)
(*Solving the EoM *)


(* ::Subsection:: *)
(*Solve EOM at a fixed order in the EFT power counting*)


(* ::Subsubsection::Closed:: *)
(*Solve single EOM at a given order*)


SolveEOM::nosolution="The field `1` appears more than once in the equation of motion `2`. Tyring to find a solution...";
SolveEOM::multiplesolutions="The equation of motion `1` has either no or more than one solutions. The solutions are `2`.";


SolveOneEOMfixedOrder[eom_, f:(Field[l_,_,_,{}] | Bar@Field[l_,_,_,{}]), n_?IntegerQ]:=Module[
	{
		solution,
		field = Cases[eom,Field[{l,n},_,_,{}],All],
		cdPat = Unique[],
		left,
		right,
		relabelInd,
		ruleDelayed,
		cdDelayed,
		equation, const, constRule={}
	}
	,
	(* if the given field order appears more than once throw an error *)
	If[
		Length[field]>1,
		Message[SolveEOM::nosolution,f,eom]
	];
	(* if given field order does not appear in the eom then return zero solution *)
	If[
		Length[field]==0,
		Return[RuleDelayed[Field[{l,n},_,_,_],0]]
	];
	
	(* finding the solution *)
	If[Head[f]===Bar,
		field=Bar[field]
	];
	
	(* remove the NCM head for fermions in the eom *)
	equation= ((eom/.NCM[x_] :> x /; !FreeQ[x, First@field, All]) /. NCM[DiracProduct[_Proj],First@field] /; ($FieldAssociation[First@First@field][Chiral]=!=False) -> First@field );
	
	(* Solve and Reduce are not able of solving the equation m*x+y1+y2+...+yn==0 for x, where all yi are independent of x, if n is in the few hundrets. *)
	(* Thus we first define const=y1+y2+...+yn, then solve the equation m*x+const==0, before replacing back const->y1+y2+...+yn. *)
	(* For complicated EOM systems this gives SEVERAL orders of magnitude improvement for the computation time! *)
	
	(* replace all terms independent of First@field by const and create the rule constRule which allows to replace them back *)
	equation= equation/.Plus[a_,b___]/;(FreeQ[{b},First@field,All]&&!FreeQ[{a},First@field,All]):>(constRule=(const->Plus[b]); a+const);
	
	(* solve the eom for the given field and replace back the constants *)
	solution= Flatten@Solve[equation==0, First@field]/. constRule;
	(* expand the solutions *)
	solution= MapAt[LagrangianExpand, solution, {All,2}];
	
	(* throw error if there is not exactly one solution *)
	If[Length[solution]!=1,
		Message[SolveEOM::multiplesolutions,eom,solution];
		Abort[]
	];
	
	(* change the value of $currentHeavyDims for the heavy fields *)
	If[(!KeyExistsQ[$currentHeavyDims,l]) && (Last@First[solution] =!= 0),
		AppendTo[$currentHeavyDims, l -> n]
	];
	
	(* changing the solutions to patterns *)
	solution = SolutionToPattern[First@solution];
	Return[solution]
]


(* ::Subsubsection::Closed:: *)
(*Transform solution to appropriate delayed rule*)


SolutionToPattern[Rule[lhs_,rhs_]]:=Module[
	{
		left,
		right = RelabelIndices[rhs,Unique->True],
		delayedCD,
		indexRelabel,
		result
	}
	,
	(* simplified rule if solution vanishes *)
	If[rhs==0,
		With[{indexPatternCD = Unique["\[Mu]"]},
			left = lhs /. (Field[a_, b_, c_, {}] :> Field[a, IndexToPattern[b], IndexToPattern[c], Pattern[indexPatternCD,_]]);
			Return@RuleDelayed[left,0]
		];
	];
	
	(* non-vanishing solutions *)
	With[{indexPatternCD = Unique["\[Mu]"]},
		left = lhs /. (Field[a_, b_, c_, {}] :> Field[a, IndexToPattern[b], IndexToPattern[c], Pattern[indexPatternCD,_]]);
		right = indexRelabel[
			delayedCD[indexPatternCD, right],
			Unique->True
		];
		result = RuleDelayed[left, Evaluate@right];
		result = result /. delayedCD->CD /. indexRelabel->RelabelIndices;
	];
	
	Return[result]
]


(* ::Subsubsection::Closed:: *)
(*Solve all EOMs at a given order*)


SolveAllEOMfixedOrder[eoms_List,fields_List,n_?IntegerQ]:=Module[
	{
		solution,
		fullSet = Transpose[{eoms,fields}],
		eomSet
	}
	,	
	(*
	n is the number of suppression factors;
	thus the EFT order is n + mass-dimension (and not the number of suppression factos) of the field.
	*)
	eomSet = SeriesEFT[First[#],EFTOrder->{n+TypeDim@Last[#]}]& /@ fullSet;
	
	eomSet = Transpose[{eomSet,fields}];
	
	(* solve all the eom at the given order *)
	solution = Table[
		SolveOneEOMfixedOrder[First[eq],Last[eq],n]
		,
		{eq,eomSet}
	];
	Return[solution]
]


(* ::Subsection:: *)
(*Full EOM function*)


(* ::Subsubsection::Closed:: *)
(*Full solution to EOM expanded to all required orders*)


Options[FullEomSolution]={EFTOrder -> 6};


OptionMessage[EFTOrder, func:FullEomSolution, val_] := Message[General::optexpectsval, EFTOrder, func, val, "interger value"];


FullEomSolution[uvEOMs_List, uvFields_List, OptionsPattern[]] := Module[
	{
		eom,
		order = OptionValue[EFTOrder],
		solutionOrders = {},
		fullSolution
	}
	,
	
	(* Expand the EOM in the power counting *)
	eom = ExpandEOM[
		(* For heavy FieldStrength substitute the corresponding Vectorfields *)
		SubstituteHeavyFieldStrength[uvEOMs],
		uvFields,
		EFTOrder->order
	];
	
	(* solve the EOM order by order up to order -1, which is the highest EFTOrder of any heavy field *)
	Do[
		(* solve the eom at order i *)
		AppendTo[solutionOrders, SolveAllEOMfixedOrder[eom, uvFields, i]];
		(* plug in solutions of order i into the eom *)
		eom = RemovePower[eom] /. Last[solutionOrders];
		eom = LagrangianExpand@ReleaseHold[eom];
		,
		{i, 1, order-1}
	];
	solutionOrders = Flatten[solutionOrders];
	
	(* construct full solution *)
	fullSolution = (ExpandField[#,EFTOrder->order]&/@uvFields) /. solutionOrders;

	Return[{fullSolution, solutionOrders}]
]


(* ::Subsubsection::Closed:: *)
(*Main EOM function*)


DetermineEOMs::usage="DetermineEOMs[lagrangian] derives and solves the equations of motion for all heavy fields from the given lagrangian.";


Options[DetermineEOMs]={EFTOrder -> 6};


DetermineEOMs[lag_, OptionsPattern[]] := Module[
	{
		uvFields,
		eoms,
		order = OptionValue[EFTOrder]
	}
	,
	(* find all UV fields in the given lag *)
	uvFields = FindUvFields[lag];
	(* derive the eom for the UV fields from lag *)
	eoms = DeriveEOM[lag, uvFields, EFTOrder->order];
	(* solve eoms, and store the solutions in $currentEOMs *)
	(* $currentEOMs[[1]] is the inclusive solution up to EFTOrder and $currentEOMs[[2]] is the specific solution for each order *)
	$currentEOMs = FullEomSolution[eoms, uvFields, EFTOrder->order];
	(* only return inclusive solution *)
	Return@First[$currentEOMs]
]


(* ::Section:: *)
(*Printing EOMs*)


(* ::Subsubsection::Closed:: *)
(*Print EOM*)


Options[PrintEOM]={Field -> All, EFTOrder -> 6};


PrintEOM::FieldOptionValue = "The option value for Field must be All, or specified as a List."


PrintEOM[expr_, OptionsPattern[]] := Module[
	{
		lagrangian = Expand[expr],
		order = OptionValue[EFTOrder],
		inclusiveOrder,
		uvFields
	}
	,
	inclusiveOrder = If[MatchQ[order,{_}],
		First[order],
		order
	];
	
	uvFields = FindUvFields[lagrangian];
	(* derive eom *)
	SetCurrentLagrangian[lagrangian, 0, inclusiveOrder];
	
	(* print EOMs *)
	Switch[OptionValue[Field],
		All  , Nothing,
		{___}, uvFields = DeleteCases[uvFields, Field[Except[Alternatives@@OptionValue[Field]],___]],
		_    , Message[PrintEOM::FieldOptionValue]; Abort[]
	];
	Module[
		{
			fields = uvFields /. Field[a___,b:Except[{}]]:>Field[a,{}],
			aux
		}
		,
		aux = (#->SeriesEFT[RelabelIndices[#/. First[$currentEOMs]], EFTOrder->order]) &/@ fields;
		Return[Format[#,NiceForm]&/@aux]
	];
]
