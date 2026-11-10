(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`FeynmanRules`*)


(* ::Subtitle:: *)
(*Calculation of Feynman Rules for Matchete.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["FeynmanRules"]
PackageExport["Ext"]
PackageExport["ExtMom"]
PackageExport["Legs"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* Public Function Declarations *)
FeynmanRules::usage = "FeynmanRules[lag] computes all Feynman rules for the given Lagrangian lag.
FeynmanRules[lag, fields] computes the Feynman rules for the given ensemble of fields from the given Lagrangian lag.";


Ext::usage= "Ext[field, no] denotes the external spinors (for fermions), polraization vectors (for vectors), or labels the fields (for scalars) for Feynman rules. The argument field is the label of the corresponding fields, where for conjugated fields field is set to Conj[label]. The argument no is an integer used to number and distinguish the different fields in one vertex."


ExtMom::usage= "ExtMom[\[Mu], Ext[field, no]] denotes the momentum of the external field Ext[field, no] with Lorenty index \[Mu]."


Legs::usage= "Is an option used in Matchete."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*FeynmanRules main functions*)


(* ::Subsection:: *)
(*FeynmanRules*)


(* ::Subsubsection::Closed:: *)
(*Options, Messages, and Attributes*)


SetAttributes[FeynmanRules, HoldRest]; (*Don't evaluate Bars on field list*)


Options@ FeynmanRules= {
	Legs    -> All, 
	(*Mode    -> "Standard", *)
	Verbose -> True
	,
	"ExternalMomentum" -> "incoming",
	"FermionFlow"      -> "Standard"
};


OptionTest[FeynmanRules, "ExternalMomentum"] = MatchQ["incoming"|"outgoing"];
OptionTest[FeynmanRules, "FermionFlow"]      = MatchQ["Standard"|"DEHK"];


(* move to Matchete.m ? *)
OptionMessage[Legs, FeynmanRules, val_]               := Message[General::optexpectsval, Legs, FeynmanRules, val, "possitive integer, a list with a single positive integer entry, or 'All'"];
(*OptionMessage[Mode, FeynmanRules_,val_]               := Message[General::optexpectsval, Mode, FeynmanRules, val, "\"Standard\" or \"MadGraph\""];*)
OptionMessage[Verbose, FeynmanRules, val_]            := Message[General::optexpectsval, Verbose, FeynmanRules, val, "Boolean"];

OptionMessage["ExternalMomentum", FeynmanRules, val_] := Message[General::optexpectsval, "ExternalMomentum", FeynmanRules, val, "\"incoming\" or \"outgoing\""];
OptionMessage["FermionFlow", FeynmanRules, val_]      := Message[General::optexpectsval, "FermionFlow", FeynmanRules, val, "\"Standard\" or \"DEHK\""];


FeynmanRules::list= "The second argument must be a list of fields."; 
FeynmanRules::unkwnfield= "The second argument must be a list of fields. `1` is/are not recognized field(s).";
FeynmanRules::legopt= "Cannot specify 'Legs' option when giving a particular set of external fields.";
FeynmanRules::legs= "The option Legs can only take the following values: integer, a list with a single entry which is an integer, or All."
FeynmanRules::dupl= "The same external legs `1` have been detected in multiple Feynman rules."


(* ::Subsubsection::Closed:: *)
(*Single Term*)


FeynmanRules[lagrangian_, inputFields_List, opts:OptionsPattern[]]? OptionsCheck:= Module[
	{
		fields,
		expandedExpr,
		interactionTerms,
		positionPowerList, 
		unbaredFields,
		indices,
		InteractionTermQ, 
		feynmanRule
	},
	
	(* Legs option invalid when looking at specific vertex *)
	If[MemberQ[{opts}, Rule[Legs, _]],
		Message[FeynmanRules::legopt];
		Abort[];
	];
	
	(*Check fields input*)
	If[Head @inputFields =!= List,
		Message[FeynmanRules::list];
		Abort[];
	];
	If[!SubsetQ[Keys@ $FieldAssociation, inputFields/. Conj-> Identity],
		Message[FeynmanRules::unkwnfield, Complement[inputFields, Keys@ $FieldAssociation]];
		Abort[];
	];
	
	(*Check that the lagrangian terms make sense individually*)
	If[!CheckLagrangian[lagrangian, 
			CanonicallyNormalized-> False, 
			GaugeInvariant->        False,
			GaugeAnomalies->        False,
			HeavyMassBasis->        False,
			Hermiticity->           False
	], 
		Abort[]
	];
	
	(*Get the fields with indices*)
	fields= Hold@ inputFields/. Bar-> Conj// ReleaseHold;
	fields= Table[FieldDoFs[fields[[n]]]@ Symbol["d$$" <> ToString@ n], {n, Length@ fields}];
	fields= If[FreeQ[#, Bar@ _Field], Identity, Bar]@ FirstCase[#, _Field| _Field, {}, All] &/@ fields;
	
	(* drop theta terms *)
	(* 
	Note: 
	While the sum of all contributions to the Feynman rules by theta terms adds up to zero, 
	this is not easily seen or identified since individual term show up and 
	GreensSimplify and similar functions cannot be used on Feynman rules level. 
	Thus we drop theta term contributions from the start. 
	*)
	expandedExpr= RelabelIndices[lagrangian - GetThetaTerms[lagrangian]];
	
	(*Expand out all fields in the Lagrangian*)
	expandedExpr= ExpandLagrangianFields@ expandedExpr;
	
    (*Now get all the relevant interaction terms which include the given fields*)
    unbaredFields = fields /. Bar[x_] :> x;
	InteractionTermQ[term_]:= AllTrue[
			Tally[unbaredFields[[All, 1]]], 
			#[[2]] == Count[term, #[[1]], {0, Infinity}]&
		] && Length[unbaredFields] === Count[term, Field[___], {0, Infinity}];
	
	(*Extract terms that include all symbols in the list with the correct multiplicities*)
	interactionTerms = Select[expandedExpr, InteractionTermQ];
	If[OptionValue@ Verbose,
		Print["There is/are ",Length[interactionTerms]," interaction term(s) containing the fields ", unbaredFields[[All, 1]]];
	];
	
	(* remove GammaCC in MadGraph mode *)
	If[OptionValue["FermionFlow"]==="DEHK",
		interactionTerms = interactionTerms/. GammaCC-> GammaCCFR;
		interactionTerms= interactionTerms//. DiracProduct[A___, GammaCCFR, B___]:> DiracProduct[A, B];
	];
	
	(*Calculate the feynman rules*)
	feynmanRule= GetFeynmanRules[Total@ RelabelIndices[interactionTerms, Unique->True], 
		fields, "FermionFlow"-> OptionValue@ "FermionFlow"];
	
		(*Flip sign of external momentum for MadGraph*)
	If[OptionValue@ "ExternalMomentum" === "outgoing",
		feynmanRule= feynmanRule/. mom_ExtMom:> -mom;
	];
	
	feynmanRule
];


(* ::Subsubsection::Closed:: *)
(*All Rules*)


FeynmanRules[lagrangian_, OptionsPattern[]]? OptionsCheck:= Module[
	{
		expandedExpr,
		positionPowerList,
		interactionL,
		rules= <||>,
		startComputation,
		fields,
		fieldCounts,
		selectedTerms,
		interactionTerms,
		maxLegs= OptionValue@ Legs,
		terms,
		fieldKey
	}
	,
	(*Check that the lagrangian terms make sense individually*)
	If[!CheckLagrangian[lagrangian, 
			GaugeInvariant        -> False, 
			GaugeAnomalies        -> False, 
			CanonicallyNormalized -> False
		], 
		Abort[]
	];	
	
	(* drop theta terms *)
	(* 
	Note: 
	While the sum of all contributions to the Feynman rules by theta terms adds up to zero, 
	this is not easily seen or identified since individual term show up and 
	GreensSimplify and similar functions cannot be used on Feynman rules level. 
	Thus we drop theta term contributions from the start. 
	*)
	expandedExpr= RelabelIndices[lagrangian - GetThetaTerms[lagrangian]];
	
	(* TO DO: add check that Lagrangian is in mass basis *)

	(*Expand out all fields in the Lagrangian*)
	expandedExpr= ExpandLagrangianFields@ expandedExpr;
	
	(*Split interaction terms off the lagrangian*)
	(*maxLegs= If[OptionValue@ Legs === All, 100, OptionValue@ Legs];*)
	Switch[maxLegs,
		All, 
			interactionL= Plus@@ Select[expandedExpr, Count[#, _Field, All] >= 3 &],
		{_Integer}, 
			interactionL= Plus@@ Select[expandedExpr, First[maxLegs] == Count[#, _Field, All] >= 3 &],
		_Integer,
			interactionL= Plus@@ Select[expandedExpr, maxLegs >= Count[#, _Field, All] >= 3 &],
		_,
			Message[FeynmanRules::legs]; Abort[]
	];
	
	(*Compute all Feynman rules*)
	
	(* remove GammaCC in MadGraph mode *)
	If[OptionValue["FermionFlow"]==="DEHK",
		interactionL = interactionL/. GammaCC-> GammaCCFR;
		interactionL= interactionL//. DiracProduct[A___, GammaCCFR, B___]:> DiracProduct[A, B];
	];
	
	While[
		(* run this loop until no terms are left *)
		interactionL =!= 0,
		
		(* get fields appearing in first term in specific order *)
		fields = GetFieldsInOrder[TermsToList[interactionL][[1]]];
		
		(* count how often which field appears in this term *)
		fieldCounts = Counts[fields /. Field[label_,___] :> Field[label,___]];
		
		(* select all terms which contain at least the same number of each field *)
		selectedTerms = TermsToList[interactionL];
		Do[
			If[Head@field===Field && !$FieldAssociation[First@field, SelfConjugate],
				(* for complex fields we need to subtract the number of Bar@Field appearing to avoid double counting *)
				selectedTerms = Select[selectedTerms, (Count[#, field, All]-Count[#, Bar@field, All]) >= fieldCounts[field]&]
				,
				selectedTerms = Select[selectedTerms, Count[#, field, All] >= fieldCounts[field]&]
			],
			{field, Keys[fieldCounts]}
		];
		
		(* take only terms with the same total number of fields *)
		interactionTerms = Plus@@ Select[selectedTerms, Count[#, Field[___], All] === Total[Values[fieldCounts]]&];
		
		(* subtract selected terms from interaction Lagrangian *)
		interactionL -= interactionTerms;
		
		(* go to standard form for interations *)
		interactionTerms= RelabelIndices[interactionTerms, Unique -> True];
		
		(* determine field labels as key *)
		fieldKey= Cases[fields, f:(Bar[_Field]|Except[Bar[_Field],_Field]):> If[Head@f===Bar, Conj[f[[1,1]]],f[[1]]], 1];
		
		(* compute the Feynman rules for the selected terms *)
		If[KeyExistsQ[rules, fieldKey],
			(* sum up terms if keys are duplicate \[Rule] should never happen *)
			Message[FeynmanRules::dupl, fieldKey]; (* warning *)
			rules[fieldKey] += GetFeynmanRules[interactionTerms, RelabelIndicesInOrder[fields], "FermionFlow"-> OptionValue["FermionFlow"]]
		,
			rules[fieldKey]= GetFeynmanRules[interactionTerms, RelabelIndicesInOrder[fields], "FermionFlow"-> OptionValue["FermionFlow"]]
		];
		
		(*AppendTo[rules, GetFeynmanRules[interactionTerms, RelabelIndicesInOrder[fields], Mode-> OptionValue[Mode]]];*)		
	];
	
	(*Flip sign of external momentum for MadGraph*)
	If[OptionValue@ "ExternalMomentum" === "outgoing",
		rules= rules/. mom_ExtMom:> -mom;
	];
	
	(*,"Compute all Feynman rules"];*)

	(*Parallel implementation:
	(*Get all field labels occuring in the interaction lagrangian*)
	allFields = EchoTiming[DeleteDuplicates@Cases[interactionL, Field[a_, ___] -> a, Infinity]];

	(*Convert to list*)
	terms = EchoTiming[TermsToList@interactionL];

	(*Gather terms together by  field counts*)
	terms = EchoTiming[GatherBy[terms, Table[{Count[#, Field[field, ___], Infinity], Count[#, Bar@Field[field, ___], Infinity]}, {field, allFields}]&]];

	(*Sum all terms with the same field content*)
	terms = EchoTiming[Total /@ terms];

	(*Relabel indices*)
	interactionTerms = EchoTiming[RelabelIndices[#, Unique->True]& /@ ReleasePseudoTimes@terms];

	(*Calculate the Feynman rules for the terms in parallel*)
	ParallelNeeds["Matchete`"];
	DistributeDefinitions[$FieldAssociation, $FlavorIndices, $DropDiagonalCouplings];
	rules = ParallelMap[GetFeynmanRules, interactionTerms];*)
	
	(* sort field labels in each Key *)
	rules= KeyMap[SortFieldLabels, rules];
	(* sort the keys by field content *)
	rules= KeySortBy[rules, SortFunctionKeysFR];
	
	rules
]


(* ::Subsubsection::Closed:: *)
(*Sorting function for Feynman rules*)


(* ::Text:: *)
(*Takes a list of field labels (a key of the association returned by FeynmanRules) as input and assigns a score to it by which association can be sorted into a canonical ordering. *)


SortFunctionKeysFR[expr_List]:= Module[{weight},
	weight= 100*Length[expr];
	weight= weight+ Total[(WeightFieldTypes/@ expr)[[;;,1]]] (* the 1st returned value from WeightFieldTypes is the weight whereas the 2nd value is not needed here *)
]


(* ::Text:: *)
(*Takes a field label as input and returns an assigned score depending on the field type. The second argument is optional and allows to manually add to the score (if not given it defaults to zero) . *)


WeightFieldTypes[l_, num_:0]:= {Switch[$FieldAssociation[l/. (*OverBar|*)Conj->Identity, Type],
	Vector,  1,
	Scalar,  2+ If[$FieldAssociation[l/. (*OverBar|*)Conj->Identity, GoldstoneBoson], 0.1, 0], (* list Goldstones after physical scalars *)
	(* 
	Add the off-set counter num (defaults to 0) such that fermions are not sorted alphabetically, which potentially changes the global minus sign.
	*)
	Fermion, 3+ num,
	_,       5+ num (* Is this necessary for ghosts? *)
], l/. (*OverBar|*)Conj->Identity}


(* ::Text:: *)
(*Function to sort the labels in one key.*)


SortFieldLabels[expr_List]:= Module[{tmp=0.},
	(* tmp is used such that fermions are not sorted alphabetically *)
	SortBy[expr, WeightFieldTypes[#,tmp=tmp+0.001]&, ComparisonForSorting]
]


(* ::Text:: *)
(*This function ensures that the field labels are sorted according to their score, and if the score is identical, they are sorted canonically, but ignoring the Conj (OverBar) heads...*)


ComparisonForSorting[a_,b_]:= If[a[[1]]<b[[1]], 
	True, 
	If[a[[1]]==b[[1]],
		OrderedQ[{a[[2]],b[[2]]}]
		,
		False
	]
]


(* ::Subsubsection::Closed:: *)
(*ExpandLagrangianFields*)


(* ::Text:: *)
(*Expands out all gauge fields from CDs and FSs and then rescales the gauge fields to move the gauge couplings from the gauge kinetic term to interaction terms. *)
(*Returns a list of Lagrangian terms held by PseudoTimes *)


ExpandLagrangianFields@ lag_:= Module[
	{couplingRules, expandedLag= LagrangianExpand@ lag, maxDim}
	,
	(* determine maximum EFT order *)
	maxDim= GetMaxEFTOrder[expandedLag, False];
	
	(* Expand the covariant derivatives and field strength tensors *)
	expandedLag = ExpandFS@ ExpandCD@ expandedLag;
	
	(* Reintroduce gauge couplings *)
    couplingRules = GetGaugeCouplingConventionReplacement@ lag;
    (* code below uses default gauge coulpings only, which potentially do not match the couplings in front of the gauge kinetic terms *)
    (*couplingRules = {Field[#[Field],rest___]->#[Coupling][]Field[#[Field],rest]}&/@ Values[GetGaugeGroups[]];*)
    
    (* Series necessary if gauge couplings receive NP shifts *)
    expandedLag = SeriesEFT[expandedLag /. couplingRules, EFTOrder->maxDim, ExpandDenominator-> True, HeavyFieldDimension-> False];
    
    (* Expand the fields such that they can be counted properly *)
    TermsToList@ PseudoTimes@ RelabelIndices@ expandedLag
];


(* ::Subsubsection::Closed:: *)
(*GetFeynmanRules*)


Options@ GetFeynmanRules= {"FermionFlow"-> "Standard"}


(*GetFeynmanRules[expr_] := GetFeynmanRules[expr, RelabelIndicesInOrder @ GetFieldsInOrder @ PseudoTimes @ First @ TermsToList[expr]];*)

GetFeynmanRules[interactionTerms_, fields_List, OptionsPattern[]] := Module[
	{
		feynmanRules,
		expr= interactionTerms, 
		listedExpr,
		fieldNumber = 1,
		expandPDint,
		expandBarPDint
	}
	,
	(* take functional derivatives w.r.t. all fields appearing in the given terms *)
	Do[
		listedExpr= TermsToList@ expr;
		expr= Sum[FunctionalD[summand, field, fieldNumber, "FermionFlow"-> OptionValue["FermionFlow"]], {summand, listedExpr}];
		(* increase counter *)
		fieldNumber++
		,
		{field, fields}
	];
	
	(* reactivate Times *)
	feynmanRules = ReleasePseudoTimes[expr];
	
	(*Replace every partial derivative with the momentum of the particle (implementation of fieldnumber index )*)
	(*For incomming momentat: \[PartialD] \[Rule] -i k *)
	expandPDint= PDint[{\[Mu]_}, Field[f_Ext, type_, {indices___, rest___}, {lIndices___}]] :> 
		-I ExtMom[\[Mu], f]Field[f, type, {indices, rest}, {lIndices}];
	
	expandBarPDint= PDint[{\[Mu]_}, Bar[Field[f_Ext, type_, {indices___, rest___}, {lIndices___}]]] :> 
		-I ExtMom[\[Mu], f]Bar@Field[f, type, {indices, rest}, {lIndices}];
	
	(*The i here comes from the expansion of exp{iSint}*)
	feynmanRules= ReleaseHold@ feynmanRules;
	
	(*Simplify@RelabelIndices@ReplaceCGs@GreensSimplify@Contract[I feynmanRules //. expandPDint //. expandBarPDint]*)
	(*EchoTiming@(*Simplify@*)RelabelIndices@Contract[I EchoTiming[feynmanRules //. expandPDint //. expandBarPDint,"rules"]]*)
	Collect[
		RelabelIndices@ SimplifySymmetries@ Contract[I* feynmanRules//. expandPDint//. expandBarPDint], 
		_Field | Bar[_Field] | Transp[_Field] | Transp[Bar[_Field]] | _NCM | _GrassmannProduct
	]
]


(* ::Subsection:: *)
(*FunctionalDerivative*)


(* ::Subsubsection::Closed:: *)
(*Option*)


Options@ FunctionalD= {"FermionFlow"-> "Standard"}


(* ::Subsubsection::Closed:: *)
(*Errors*)


FunctionalD::mode = "'FermionFlow \[Rule] `1`' is invalid. Only 'FermionFlow \[Rule] \"Standard\"' and 'FermionFlow \[Rule] \"DEHK\"' are currently supported.";


FunctionalD::unsupported = "Unsupported field type: `1`.";


FunctionalD[expr_, field_, fieldNumber_, OptionsPattern[]] := Module[
    {},
    Message[FunctionalD::unsupported, field];
    Null
]


(* ::Subsubsection::Closed:: *)
(*Scalars*)


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Field[_, Scalar, _, _]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
		label = First[field], 
		positionField,
		positionBarField
	}
	,
	(*get position of all fields (The only important variable is positionField)*)
	positionField    = Position[expr, Field[label, ___]];
	positionBarField = Position[expr, Bar[Field[label, ___]]];
	positionField    = Select[positionField, Not[MemberQ[positionBarField,Most[#]]]&];
	
	(*replace all occurences of field with an external leg*)
	Total[
		ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
				(* Create the product of Delta terms *)
				deltaProduct = Times @@ Table[Delta[Extract[expr, #][[3, i]], field[[3, i]]], {i, Length[field[[3]]]}];
				
				(* Replace with modified expression *)
				deltaProduct * Field[Ext[label, fieldNumber], Scalar, field[[3]], Extract[expr, #][[4]]]
			] 
		]& /@ positionField
	]
]


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Bar[Field[_, Scalar, _, _]]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
		label = First[field][[1]],
		positionBarField
	}
	,
	(*get position of all bared fields *)
	positionBarField = Position[expr, Bar[Field[label, ___]]];
	
	(*replace all occurences of field with an external leg*)
	Total[
		ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
			
			(* Create the product of Delta terms *)
			deltaProduct = Times @@ Table[Delta[First[Extract[expr, #]][[3, i]], First[field][[3, i]]], {i, Length[First[field][[3]]]}];
			
			(* Replace with modified expression *)
			deltaProduct * Bar[Field[Ext[(*OverBar*)Conj[label], fieldNumber], Scalar, First[field][[3]], First[Extract[expr, #]][[4]]]]
			] 
		] & /@ positionBarField
	]
]


(* ::Subsubsection::Closed:: *)
(*Fermions/Ghosts*)


(* ::Text:: *)
(*TranpMarker ensures that if a derivative of a transposed field is taken w.r.t. the non-transposed field, the remaining spin-chain is transposed:  df^T/df X = X^T*)


Transp[TranspMarker]:= TransposeThisSpinChain
(* if the TranspMarker gets transposed, it indicates that the entire corresponding spin-chain should be transposed *)
TransposeThisSpinChain/: Times[TransposeThisSpinChain, x___]:= NCM[TransposeThisSpinChain, Times[x]]


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Field[_, Fermion|Ghost|AntiGhost, _, _]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
		label = First[field],
		positionGrassmann,
		positionField, 
		positionBarField,
		minusGrassmann= Switch[OptionValue["FermionFlow"],
			"Standard", -1,
			"DEHK",     +1,
			_, Message[FunctionalD::mode]; Abort[]
		],
		res
	}
	,
	(*get position of all grassmann fields to account for anticommutativity*)
	positionGrassmann = Position[
		expr, 
		(*_?(MatchQ[#, Field[_,Fermion|Ghost|AntiGhost,___]] && Not@MatchQ[#, Field[Ext[_],___]] &)*)
		_?(MatchQ[#, Field[Except[_Ext],Fermion|Ghost|AntiGhost,___]] &)
	];
	(* get position of all fields *)
	positionField    = Position[expr, Field[label, ___]];
	positionBarField = Position[expr, Bar[Field[label,___]]];
	positionField    = Select[positionField,Not[MemberQ[positionBarField,Most[#]]]&];
	
	(*replace all occurences of field with an external leg*)
	res= Flatten@@ Flatten@@ Total[
		(minusGrassmann)^(1+Position[positionGrassmann, #]) * ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
				(* Create the product of Delta terms *)
				deltaProduct = Times @@ Table[Delta[Extract[expr, #][[3, i]], field[[3, i]]], {i, Length[field[[3]]]}];
				
				(* Replace with modified expression *)
				(* if inserted into a Trans[...] the TranspMarker ensures that this spin-chain gets transposed *)
				TranspMarker * deltaProduct * Field[Ext[label, fieldNumber], field[[2]], field[[3]], Extract[expr, #][[4]]]
			] 
		]& /@ positionField
	];
	
	(* remove remaining TranpMarker and transpose required spin-chains with FuncDSimplify *)
	res= res/. TranspMarker-> 1;
	FuncDSimplify[res, Contract-> False]
]


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Bar[Field[_, Fermion|Ghost|AntiGhost, _, _]]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
		label = First[field][[1]],
		positionAll, 
		positionGrassmann,
		positionBarGrassmann,
		positionBarField,
		minusGrassmann= Switch[OptionValue["FermionFlow"],
			"Standard", -1,
			"DEHK",     +1,
			_, Message[FunctionalD::mode]; Abort[]
		],
		res
	}
	,
	(*get position of all fermionic fields to account for anticommutativity (The only important variable is positionBarField and positionBarFermions)*)
	positionGrassmann = Position[
		expr, 
		_?(MatchQ[#, Field[Except[_Ext],Fermion|Ghost|AntiGhost,___]] &)
	];
    positionBarGrassmann = Most /@ positionGrassmann; (*position of fermions without bar*)
	positionBarField = Position[expr, Bar[Field[label,___]]];
	
	(*replace all occurences of field with an external leg*)
	res= Flatten@@ Flatten@@ Total[
		(minusGrassmann)^(1+Position[positionBarGrassmann, #]) * ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
				(* Create the product of Delta terms *)
				deltaProduct = Times @@ Table[Delta[First[Extract[expr, #]][[3, i]], First[field][[3, i]]], {i, Length[First[field][[3]]]}];
				
				(* Replace with modified expression *)
				(* if inserted into a Trans[...] the TranspMarker ensures that this spin-chain gets transposed *)
				TranspMarker * deltaProduct * Bar[Field[Ext[(*OverBar*)Conj[label], fieldNumber], First[field][[2]], First[field][[3]], First[Extract[expr, #]][[4]]]]
			] 
		]& /@ positionBarField
	];
	
	(* remove remaining TranpMarker and transpose required spin-chains with FuncDSimplify *)
	res= res/. TranspMarker-> 1;
	FuncDSimplify[res, Contract-> False]
]


(* ::Subsubsection::Closed:: *)
(*Vectors*)


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Field[_, Vector[___], _, _]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
		label = First[field],
		positionField,
		positionBarField
	}
	,
	(* get position of all fields *)
	positionField    = Position[expr, Field[label, ___]];
	positionBarField = Position[expr, Bar[Field[label,___]]];
	positionField    = Select[positionField,Not[MemberQ[positionBarField,Most[#]]]&];
	
	(*replace all occurences of field with an external leg*)
	Total[
		ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
				(* Create the product of Delta terms *)
				deltaProduct = Times @@ Table[Delta[Extract[expr, #][[3, i]], field[[3, i]]], {i, Length[field[[3]]]}];
				
				(* Replace with modified expression *)
				deltaProduct * Metric[Extract[expr, #][[2, 1]], field[[2, 1]]] * Field[Ext[label, fieldNumber], field[[2]], field[[3]], Extract[expr, #][[4]]]
			] 
		]& /@ positionField
	]
]


FunctionalD[expr_, field_, fieldNumber_, opt:OptionsPattern[]] /; MatchQ[field, Bar[Field[_, Vector[___], _, _]]] := FunctionalD[expr, field, fieldNumber, opt] = Module[
	{
	label = First[field][[1]],
	positionBarField
	}
	,
	(*get position of all fields (The only important variable is positionField)*)
	positionBarField = Position[expr, Bar@Field[label, ___]];
	
	(*replace all occurences of field with an external leg*)
	Total[
		ReplacePart[
			expr, 
			# -> Module[{deltaProduct},
				(* Create the product of Delta terms *)
				deltaProduct = Times @@ Table[Delta[First[Extract[expr, #]][[3, i]], First[field][[3, i]]], {i, Length[First[field][[3]]]}];
				
				(* Replace with modified expression *)
				deltaProduct * Metric[First[Extract[expr, #]][[2, 1]], First[field][[2, 1]]] * Bar@Field[Ext[(*OverBar*)Conj[label], fieldNumber], First[field][[2]], First[field][[3]], First[Extract[expr, #]][[4]]]
			] 
		]& /@ positionBarField
	]
]


(* ::Section:: *)
(*Non-covariant expansions*)


(* ::Subsubsection::Closed:: *)
(*Expand covariant derivatives*)


ExpandCD::usage = "ExpandCD[expr] expands each covariant derivative in momentum space.";


ExpandCD[arg_] := Module[{
		expr= PseudoTimes[arg], 
		nonGaugeGroups= Alternatives@@ Join[Keys@ $FlavorIndices, Blank/@ Keys@ $GlobalGroups]},
	(* expand CD of Fields and FieldStrengths *)
	expr= expr//. f:(Field|FieldStrength)[__, {__}]:> ExpandCDSingleField[f, nonGaugeGroups];
	
	ReleasePseudoTimes@ expr	
]


ExpandCDSingleField[(f:Field|FieldStrength)[label_, aux_, indices_List, {\[Mu]_, lIndices___}], nonGaugeGroups_]:= Module[
		{abelianTerms, (*antifundSign,*) charge, entry, fieldIndices, gaugeField, gaugeIndices, generator, group, nonAbelianTerms, partial, x, y},
	
	partial= PDint[{\[Mu]/. PartialDev-> Identity}, f[label, aux, indices, {lIndices}]];
	
	(*Shortcut if derivative was partial*)
	If[Head@ \[Mu] === PartialDev, Return@ partial];
	
	(* abelian terms of the covariant derivative *)
	abelianTerms = Table[
		(*coupling = GetGaugeGroups[Head[group]][Coupling][];*)(* gauge couplings are associated to propagators in our convention *)
		charge = First[group];
		gaugeField = GetGaugeGroups[Head[group]][Field][First@ \[Mu]];
		-I * charge * gaugeField * f[label, aux, indices, {lIndices}]
		,
		{group, $FieldAssociation[label, Charges]}
	];
	
	(* disregard flavor and global indices *)
	gaugeIndices = DeleteCases[
		indices, 
		Index[_, nonGaugeGroups] | _Integer
	];
	
	(* non-abelian terms of the covariant derivative *)
	nonAbelianTerms = Table[
		(* this is now taken care of directly by FieldGenerators[...] *)
		(*
		antifundSign = Which[
			MatchQ[entry, Bar[_]], -1, True, +1
		]; (* Identify transformation behaviour based on the Bar in front of the gauge group*)
		*)
		entry = entry /. Bar[x_] :> x; (* set every index (fund or antifund) equal *)
		(* get some unique index labels *)
		x = Unique[]; (* adj index *)
		y = Unique[]; (* genric rep index *)
		(*coupling = GetGaugeGroups[Head @ entry[[2]]][Coupling][];*) (* gauge couplings are associated to propagators in our convention *)
		generator = FieldGenerators[label, Head @ entry[[2]], {x, entry[[1]], y}];
		gaugeField = GetGaugeGroups[Head @ entry[[2]]][Field][First@ \[Mu], x];
		fieldIndices = indices /. entry :> Index[y, entry[[2]]];
		(*-I coupling generator gaugeField Field[field, type, fieldIndices, {}]*);
		-I* (*antifundSign* *) generator* gaugeField* f[label, aux, fieldIndices, {lIndices}]
		,
		{entry, gaugeIndices} 
	];
	(* recursive expansions from outside to inside as otherwise there is an issue with due to the non-commutativety of PDint and CD *)
	partial + Total[nonAbelianTerms] + Total[abelianTerms]	
]


(* ::Subsubsection::Closed:: *)
(*Expand field-strength tensors*)


ExpandFS::usage = "ExpandFS[expr] expands each field-strength tensor in momentum space.";


ExpandFS[arg_] := Module[
	{
		flavors,
		group,
		x,
		y,
		coupling,
		generator,
		nonGaugeGroups= Alternatives@@Join[Keys@GetFlavorIndices[], Blank/@Keys[GetGlobalGroups[]]],
		expr= PseudoTimes[arg]
	}
	,
	(* first expand CD if derivatives of FS are present *)
	expr= expr //. fs:FieldStrength[_, {_,_}, _List, {__}] :> ExpandCD[fs];
			
	(* expand FS *)
	expr= expr //. FieldStrength[field_, {index1_, index2_}, indices_List, {}] :> Module[{},
		If[FreeQ[Lookup[Values[GetGaugeGroups[]],Field], field],
			ExpandCD[CD[index1, field[index2, Sequence @@ indices]] - CD[index2, field[index1, Sequence @@ indices]]]
			,
			Switch[indices,
				{},
					ExpandCD[CD[index1, field[index2]] - CD[index2, field[index1]]],
				_,
					(* delete non-gauge indices *)
					group = DeleteCases[
						indices, 
						Index[_, nonGaugeGroups] | _Integer
					];
					
					PDint[{index1}, field @@ Join[{index2}, indices]] - PDint[{index2}, field @@ Join[{index1}, indices]] +
					If[Length[group] == 0,
						0
						,
						x = Unique[];
						y = Unique[];
						(*coupling = GetGaugeGroups[Head @ group[[1, 2]]][Coupling][];*) (* gauge couplings are associated to propagators in our convention *)
						generator = FieldGenerators[field, Head @ group[[1, 2]], {x, group[[1, 1]], y}];
						- I generator field[index1, x] field[index2, y]
					]
			]
		]
	];
	
	ReleasePseudoTimes[expr]
];


(* ::Section:: *)
(*PartialDerivative*)


(* ::Subsubsection::Closed:: *)
(*Rules*)


(* this is potentially problematic since partial and covariant derivatives do not commute, so one should never use CD[{\[Mu]}, PDint[{\[Nu]}, Field[...]]] *)
(*TakeDev[{lindices_}, PDint[{\[Mu]_}, Field[label_, type_, indices_, {}]]] := PDint[{\[Mu]}, Field[label, type, indices, {lindices}]];*)


PDintZeroQ[f_? PDintZeroQ[x__]]:=And@@ PDintZeroQ/@ {x};
PDintZeroQ[f_[x__]]:= False;
PDintZeroQ[_]= True;
PDintZeroQ[Field]= False;
PDintZeroQ[FieldStrength]= False;
PDintZeroQ[OpenCD]= False;
PDintZeroQ[WilsonTerm]= False;


PDint[{}, expr_]:= expr;
PDint[{\[Mu]_, \[Nu]__}, expr_]:= PDint[{\[Mu]}, PDint[{\[Nu]}, expr]];
PDint[\[Mu]_, _?PDintZeroQ]= 0;
PDint[\[Mu]_, expr:Alternatives[_Plus, _List]]:= PDint[\[Mu], #]&/@ expr;
PDint[\[Mu]_, expr:Alternatives[_Times, _NCM, _NCProduct, _FuncNCM, _GrassmannProduct]]:=
	Module[{n}, Sum[MapAt[PDint[\[Mu], #]&, expr, n], {n, Length@expr}] ];
PDint[{\[Mu]_}, Power[expr_, n_]]:= n PDint[{\[Mu]}, expr] Power[expr, n-1];


Bar@PDint[\[Mu]_,Field[label_,type_,indices_,\[Nu]_]]        := PDint[\[Mu], Bar@Field[label,type,indices,\[Nu]]];
Transp@PDint[\[Mu]_,Field[label_,type_,indices_,\[Nu]_]]     := PDint[\[Mu], Transp@Field[label,type,indices,\[Nu]]];
Transp@Bar@PDint[\[Mu]_,Field[label_,type_,indices_,\[Nu]_]] := PDint[\[Mu], Transp@Bar@Field[label,type,indices,\[Nu]]];
Bar@PDint[\[Mu]_, FieldStrength[label_,li_,indices_,\[Nu]_]] := PDint[\[Mu], Bar@FieldStrength[label,li,indices,\[Nu]]];


(* ::Section:: *)
(*Charge conjugation matrix that positions itself next to the charge conjugated field*)


(* ::Subsubsection::Closed:: *)
(*GammaCCFR*)


(*Assumes that GammaCCFR stands at the left side of the DiracProduct initially*)
NCM[Bar[f1_Field], DiracProduct[GammaCCFR, a: Except[_Proj|Gamma5] ..], Transp[Bar[f2_Field]]]:= PositionGammaCCFRatRight[Bar[f1], DiracProduct[GammaCCFR, a], Transp[Bar[f2]]];

PositionGammaCCFRatRight[Bar[f1_Field], DiracProduct[a___, GammaCCFR, Transp@g_GammaM, b___], Transp[Bar[f2_Field]]]:= Power[-1, Ceiling[Length@ g/2]] PositionGammaCCFRatRight[Bar[f1], DiracProduct[a, g, GammaCCFR, b], Transp[Bar[f2]]];
PositionGammaCCFRatRight[Bar[f1_Field], DiracProduct[a___, GammaCCFR, g_GammaM, b___], Transp[Bar[f2_Field]]]:= Power[-1, Ceiling[Length@ g/2]] PositionGammaCCFRatRight[Bar[f1], DiracProduct[a, Transp@ g, GammaCCFR, b], Transp[Bar[f2]]];
PositionGammaCCFRatRight[Bar[f1_Field], DiracProduct[a___, GammaCCFR], Transp[Bar[f2_Field]]]:= NCM[Bar[f1], DiracProduct[a, GammaCCFR], Transp[Bar[f2]]];

(*Projectors and Gamma5 always stand at the right side*)
PositionGammaCCFRatRight[Bar[f1_Field], DiracProduct[a___, GammaCCFR, Transp@s: _Proj | Gamma5], Transp[Bar[f2_Field]]]:= NCM[Bar[f1], DiracProduct[a, GammaCCFR, s, b], Transp[Bar[f2]]];
PositionGammaCCFRatRight[Bar[f1_Field], DiracProduct[a___, GammaCCFR, s: _Proj | Gamma5], Transp[Bar[f2_Field]]]:= NCM[Bar[f1], DiracProduct[a, GammaCCFR, s], Transp[Bar[f2]]];


(* ::Section:: *)
(*Help Functions*)


(* ::Subsubsection::Closed:: *)
(*Get fields with ordered canonical index names*)


RelabelIndicesInOrder::usage= "RelabelIndicesInOrder[fields] takes a list of Fields (fields) and returns the same fields with ordered canonical index labels."


RelabelIndicesInOrder[fields_List] := MapIndexed[
	Module[
		{
			field = #1, 
			idx   = First[#2]
		}
		,
		field/. Index[_,group_]:>Index[ToExpression["d$$"<>ToString[idx]],group]
	] &, fields
];


(* ::Subsubsection::Closed:: *)
(*Get ordered fields *)


FeynmanRules::unknowntype= "The field type `1`is unknown. "


typeOrder[type_] := Switch[
	type,
	Fermion, 1,
	Vector[_], 2,
	Scalar, 3,
	Ghost|AntiGhost, 4,
	_, Message[FeynmanRules::unknowntype,type]; 5
]


(* Resort fields with typeOrder to match helas convention *)
sortFieldsByType[fields_List] := Module[
	{keys}
	,
	keys = MapIndexed[
		{ 
			typeOrder[Replace[#1, Bar[x_] :> x][[2]]],
			First[#2]
		} &,
		fields
	];
	fields[[Ordering[keys]]]
];


GetFieldsInOrder::usage= "GetFieldsInOrder[term] returns an ordered List of all fields appearing in the given Lagrangian term.";


GetFieldsInOrder[expr_] := Module[
	{
		raw,
		keep,
		final
	}
	,
	
	raw = Cases[expr, x: (Field[___]| Bar[Field[___]]), All];
	keep = MapIndexed[
		Function[
			{field, idx},
			!(MatchQ[field, Field[___]] && idx[[1]] < Length[raw] && If[$FieldAssociation[First[field], SelfConjugate], False, MatchQ[raw[[idx[[1]] + 1]], Bar[field]]])
		],
		raw
	];
	
	final = Reverse[Pick[raw, keep] /. Field[label_, type_, indices_, ___] :> Field[label, type, indices, {}]];

	sortFieldsByType[final]
];


(* ::Subsubsection::Closed:: *)
(*Get couplings of gauge fields*)


GetGaugeCouplingConventionReplacement::usage= "GetGaugeCouplingConventionReplacement[lag] reads off the couplings (g) in front of gauge field-strength tensors -1/(4g^2)F^2 in the given Lagrangian (lag) and returns a replacement rule F->g*F, which can be used to map the usual convention.";


GetGaugeCouplingConventionReplacement@ lag_:= Module[{assoc, new, requiredGaugeFields, missingGaugeFields},
	(* Get coupling values for all gauge fields from kinetic terms in the Lagrangian. *)
	new= GetGaugeCouplingsFromLagrangian[lag, Power-> -1/2];
	(* get gauge fields appearing in the present CD and FS *)
	requiredGaugeFields= #[Field]&/@ (GetGaugeGroupsFromCDandFS[lag]/. $GaugeGroups);
	(* get the missing gauge fields and throw a warning potentially *)
	missingGaugeFields= Complement[requiredGaugeFields, Keys[new]];
	If[Length[missingGaugeFields]>0,
		Print["No kinetic term could be found for the gauge fields " <> ToString@ missingGaugeFields <> " to determine their gauge coupling normalization. Standard normalization for the corresponding gauge couplings is assumed."];
	];
	
	(*Get default coupling values for all gauge fields that do not have kinetic term in the Lagrangian (term)*)
	assoc= Association@@ (#[Field]-> #[Coupling][]&)/@ List@@ $GaugeGroups;
	
	(*Overwrite the updated values*)
	Do[
		assoc@ field= new@ field;
	, {field, Keys@ new}];
	
	(*Return replacement association for the gauge fields *)
	KeyValueMap[f:Field[#1, __]-> f #2&, assoc]
];


GetGaugeGroupsFromCDandFS::usage= "GetGaugeGroupsFromCDandFS[arg] returns a list of the gauge groups that have indices in the given argument (arg).";


GetGaugeGroupsFromCDandFS[arg_]:= Module[{expr= arg, gaugeGroups},
	(* determine non-Abelian indices *)
	expr= Cases[expr, (Field[_,_,inds_,{__}] | FieldStrength[_,_,inds_,_])-> inds, All];
	(* get representations from indies *)
	expr=  DeleteDuplicates[Cases[expr, Index[_,x_]-> x, All]]/. $Representations;
	(* drop Lorentz, flavor and global indices *)
	expr= expr/. Alternatives@@ Join[{Lorentz}, Keys@ $FlavorIndices, Keys@ $GlobalGroups]-> Nothing;
	(* get non-Abelian gauge groups *)
	gaugeGroups= #[Group]&/@ expr;
	
	(* determine Abelian representations from charges *)
	expr= DeleteDuplicates@ Cases[arg, (Field[f_,_,_,{__}] | FieldStrength[f_,___])-> f, All]/. $FieldAssociation;
	(* find charge labels *)
	expr= DeleteDuplicates[Head/@ Flatten[#[Charges]&/@ expr]];
	
	(* add Abelian gauge groups *)
	gaugeGroups= Join[expr, gaugeGroups]
]


(* ::Subsubsection::Closed:: *)
(*Select terms up to a maximum number of fields*)


SelectFieldNumber[expr_, num_] := Module[{lag, number, terms, number},
	lag = Distribute[expr, Plus, Times];
    terms = TermsToList@PseudoTimes@lag;
	If[Head@num === List,
	number = First@num;
	terms = If[Count[#,Field[__]|FieldStrength[__],Infinity] === number, #, 0]& /@ terms;,
	terms = If[Count[#,Field[__]|FieldStrength[__],Infinity] <= num, #, 0]& /@ terms;
	];
    Total@terms
]


(* ::Subsubsection::Closed:: *)
(*Get theta terms*)


(* ::Text:: *)
(*GetThetaTerms takes a Lagrangian(-like) expression as input and returns all theta terms contained.*)


GetThetaTerms[lag_]:= SelectOperatorClass[
	Total@ DeleteCases[
		TermsToList[lag/. _Field-> 0],
		x_/; FreeQ[x, _LCTensor]
	],
	{},
	4
]


(* ::Subsubsection::Closed:: *)
(*Remove terms that vanish by contraction of Lorentz indices*)


(* Only looks at terms containing at least two times the same momentum. *)
SimplifySymmetries[sum_Plus]:= SimplifySymmetries/@ sum;

SimplifySymmetries[term_]:= Module[
	{out = term},
	If[Max[Counts[Cases[out, ExtMom[_, x_] :> x, All]], 0] < 2, Return@ out;];
	out = NormalToOperatorForm[out];
	out = out - ReplaceAll[out, o_Operator:> First[IdentitiesSymmetry[o], 0]];
	OperatorToNormalForm[out, RelabelIndices->False]
]
