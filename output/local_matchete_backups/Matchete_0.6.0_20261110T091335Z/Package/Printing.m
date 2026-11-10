(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Printing`*)


(* ::Subtitle:: *)
(*Sub package providing the printing routine NiceForm.*)


(* ::Text:: *)
(*Nb. "ToBoxes[expr, NiceForm/StandardForm]" can be used to see how an expression is formatted in the internal Mathematica Box forms.*)
(*Tutorial on Box structures and notation: https://www.wolfram.com/broadcast/video.php?sx=notation&c=104&v=3570 *)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["NiceForm"]
PackageExport["$PrintIndexLabels"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["$LabelsNiceForm"]
PackageScope["AuxForm"]
PackageScope["Overline"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


NiceForm::usage          = "NiceForm[expr] prints the expression expr in a human readable form.";
$PrintIndexLabels::usage = "$PrintIndexLabels=True|False determines whether to print the representation of indices as subscripts. It is set to False by default.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* turn off an error when calling Definition[...] *)
Off[ParentForm::deflt];


(* ::Section:: *)
(*Definition of Format*)


(* ::Subsubsection::Closed:: *)
(*Properties of the NiceForm format*)


(* ::Text:: *)
(*Inherit from StandardForm*)


ParentForm[NiceForm]^= StandardForm;


(* ::Text:: *)
(*Add to list of allowed BoxForms*)


AppendTo[$BoxForms, NiceForm];


(* ::Text:: *)
(*Sanitize output to remove a wrapping FormBox[..., NiceForm] from the output cell (which is anyway not displayed). This might reduce the option for manipulating with the displayed output in other ways.*)


NiceForm/: FormBox[expr_, NiceForm]:= expr;


(* ::Text:: *)
(*Allow printing in NiceForm*)


NiceForm /: Print[A___, NiceForm[arg_], B___] := Print[A, Format[arg, NiceForm], B]


(* ::Text:: *)
(*Allow Echo@NiceForm[...]*)


NiceForm/:Echo[NiceForm[arg_]]:=Echo@Format[arg,NiceForm]


(* ::Text:: *)
(*Allow EchoFunction[NiceForm][...]*)


NiceForm/:EchoFunction[NiceForm]:=(Echo@NiceForm[#];#)&;


(* ::Text:: *)
(*Allow Message[NiceForm[...]]*)


NiceForm/:Message[msg_,NiceForm[arg_]]:=Message[msg,Format[arg,NiceForm]];


(* ::Text:: *)
(*Variables are not  defined with the head NiceForm*)


NiceForm /: Set[a_, NiceForm@b_] := (a = b; NiceForm@b)


(* ::Subsubsection::Closed:: *)
(*Auxiliary formatting form*)


(* ::Text:: *)
(*Passing the output through an Auxiliary form means that we can control the font, w/o having to modify the stylesheet.*)


AppendTo[$BoxForms, AuxForm];
ParentForm[AuxForm]^= StandardForm;


(* ::Text:: *)
(*Fix the output font irrespective of the StyleSheet (note that default seems to be Times New Roman, and users wouldn't have predefined any rules for Niceform).*)
(*Default font is "Source Code Pro" on Mac, "Consolas" on Windows, and "Bitstream Vera Sans Mono" on Linux. *)


(*Picks the default code font based on the defautl at the time of loading the package*)
(*$NiceFormFont= AbsoluteCurrentValue@ {StyleHints, "CodeFont"}; *)


(*NiceForm/: MakeBoxes[expr_, NiceForm]:= StyleBox[MakeBoxes[expr, AuxForm], FontFamily-> $NiceFormFont];*)
NiceForm/: MakeBoxes[expr_, NiceForm]:= StyleBox[MakeBoxes[expr, AuxForm], "StandardForm"];


(* ::Text:: *)
(*AuxForm should never called directly. We can use the FormBox that get's introduced for post-processing: *)


AuxForm/: FormBox[expr_, AuxForm]:= expr;


(* ::Section:: *)
(*General printing rules*)


(* ::Subsubsection::Closed:: *)
(*General style*)


(* ::Text:: *)
(*Formatting functions*)


SubscriptStyle2@ x_:= StyleBox[x, FontSize-> 10];


Overline[label_]:= Switch[
	Head@label,
	String,
		OverscriptBox[label, StyleBox["_", FontWeight-> Bold]],
	Ext,
		label /. Ext[x_, y_] :> OverscriptBox[RowBox[{OverscriptBox[DownIndices[x/.OverBar|Conj->Identity, ""],StyleBox["_", FontWeight-> Bold]]}],""], 
		(* w/o the empty OverscriptBox and DownIndices above the printing does not work\[Ellipsis] *)
	_,
		OverscriptBox[ToString@label, StyleBox["_", FontWeight-> Bold]]
]


(* ::Text:: *)
(*Remove Skeletons from arguments given to Message function*)


Unset[$MessagePrePrint]


(* ::Subsubsection::Closed:: *)
(*DiracProduct*)


AuxForm/: MakeBoxes[DiracProduct[x_, y__], AuxForm]:=
	RowBox@ Map[MakeBoxes[#, AuxForm]&, {x, y}]


AuxForm/: MakeBoxes[DiracProduct@ x_, AuxForm]:= MakeBoxes[x, AuxForm]


(* ::Subsubsection::Closed:: *)
(*Transpose*)


AuxForm/: MakeBoxes[Transp[x:Except[DiracProduct[_, __]]], AuxForm]:= SuperscriptBox[MakeBoxes[x, AuxForm], "T"];


AuxForm/: MakeBoxes[Transp[d:DiracProduct[_, __]], AuxForm]:= 
	SuperscriptBox[FormBox[
		StyleBox[
			RowBox@ {"(", MakeBoxes[d, AuxForm], ")"} /. RowBox@ {x___, RowBox@ {a___}, y___}:> RowBox@ {x, a, y},
			SpanSymmetric-> False
		], StandardForm], "T"];


(* ::Subsubsection::Closed:: *)
(*Indices*)


(* ::Text:: *)
(*Global flag to set whether index representations are printed*)


$PrintIndexLabels = False;


(* ::Text:: *)
(*UpDownIndices[label,indices] prints the label with regular indices as superindices and barred indices as subindices.*)


UpDownIndices[label_, indices_]:= Block[{labelPrint},
	labelPrint= Switch[label
		,_OverscriptBox| _String,
			label
		,_Ext,
			label /. Ext[x_, y_] :> RowBox[{DownIndices[x, ""]}]
		,_,
			ToString@ label
	];
	
	If[!MatchQ[indices, _List], 
		Return[SubscriptBox[labelPrint, MakeBoxes[indices, AuxForm]]//. 
			RowBox[{a___, RowBox@ {x__}, b___}]:> RowBox@ {a, x, b}];
	];
	
	Which[
		indices === {},
			labelPrint,
		Select[indices, Head[#] === Bar &] === {},
			SuperscriptBox[labelPrint, SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, indices], "RowDefault"] ],
		Select[indices, Head[#] =!= Bar &]==={},
			SubscriptBox[labelPrint, SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, indices[[;;, 1]]], "RowDefault"] ],
		True, 
			SubsuperscriptBox[labelPrint, SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, Cases[indices, _Bar][[;;, 1]]], "RowDefault"], 
				SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, Cases[indices, _Index|_Integer]], "RowDefault"] ]
	]/. TemplateBox[{a_}, _]-> a
];


DownIndices[label_, indices_]:= Block[{labelPrint, indexList},
	indexList = If[ListQ[indices], indices, {indices}];
	labelPrint= Switch[label
		,_OverscriptBox| _String,
			label
		,_,
			ToString@ label
		];
	SubscriptBox[labelPrint, SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, indexList], "RowDefault"] ]/. TemplateBox[{a_}, _]-> a
]


AuxForm/: MakeBoxes[Index[label_Pattern, rep_], AuxForm]:= 
		If[$PrintIndexLabels,
			(* print with rep label *)
			SuperscriptBox[label, ToString@ rep]
		,
			(* print without rep label *)
			ToString@ label				
		];


AuxForm/: MakeBoxes[Index[label_, rep_], AuxForm]:= Block[{printLabel},
	(*Determine the appropriate printed label based on the alphabet of the representation*)
	printLabel= Lookup[Lookup[$IndexAlphabets, rep, $UndefinedIndexAlphabet], label, ToString@ label];
	
	If[$PrintIndexLabels,
		(* Print w/ rep label *)
		SuperscriptBox[printLabel, ToString@ rep]
	,
		(* Print w/o rep label *)
		printLabel
	]
];


(* ::Subsubsection::Closed:: *)
(*NCM*)


(* ::Text:: *)
(*Multiple arguments*)


AuxForm/: MakeBoxes[x: NCM[_,__], AuxForm] := Module[
		{
			product=List@@x,
			rowBox={}
		}
	,
	(* put spin-chains in brackets *)
	AppendTo[rowBox, "("];
	(* print the arguments of NCM *)
	Do[
		With[{y=factor},
			If[Head[y]===Plus,
				AppendTo[rowBox, " "];
				AppendTo[
					rowBox, FormBox[
						StyleBox[
							RowBox[{"(",MakeBoxes[y, AuxForm],")"}],
							SpanSymmetric->False
						],
						StandardForm
					]
				]
				,
				AppendTo[rowBox,MakeBoxes[y, AuxForm]]
			];
		];
		(* merge the arguments with a \cdot *)
		AppendTo[rowBox,"\[CenterDot]"];
		,
		{factor, product}
	];
	(* remove the last \cdot and close brackets *)
	rowBox = Drop[rowBox,-1];
	AppendTo[rowBox,")"];
	
	(* return as StandardForm box *)
	FormBox[StyleBox[RowBox[rowBox], SpanSymmetric->False], StandardForm]
]


(* ::Text:: *)
(*Single Argument*)


AuxForm/: MakeBoxes[NCM@ arg_, AuxForm]:= MakeBoxes[arg, AuxForm];


(* ::Text:: *)
(*GrassmannProduct \[LongDash] implicit product notation (thin space, no CenterDot)*)


AuxForm/: MakeBoxes[x: GrassmannProduct[_,__], AuxForm] := Module[
		{
			product= List@@ x,
			rowBox= {"\[NegativeThinSpace]"} (*Negative thin space*)
		}
	,
	AppendTo[rowBox, "("];
	Do[
		With[{f= factor},
			AppendTo[rowBox, MakeBoxes[f, AuxForm]];
		];
		AppendTo[rowBox, "\[ThinSpace]"];
		,
		{factor, product}
	];
	rowBox= Drop[rowBox, -1];
	AppendTo[rowBox, ")"];
	AppendTo[rowBox, "\[NegativeThinSpace]"]; (*Negative thin space*)
	FormBox[StyleBox[RowBox[rowBox], SpanSymmetric-> False], StandardForm]
]


AuxForm/: MakeBoxes[GrassmannProduct@ arg_, AuxForm]:= MakeBoxes[arg, AuxForm];


(* ::Subsection:: *)
(*Plus*)


(* ::Subsubsection::Closed:: *)
(*AuxForm for Plus*)


(* ::Text:: *)
(*Definition of how sums are displayed*)


AuxForm/: MakeBoxes[a_Plus, AuxForm]:= Module[
		{
			terms= List@@ a,
			rowBox= {},
			startFlag= True
		},
		
	(* sort the factors of Times *)
	terms = SortPlus[terms];
	(* Add minus sign where appropriate *)
	Do[
		If[startFlag,
			(* First term *)
			startFlag=False;
			(* First term in the sum does not get a sign if positive *)
			If[MinusSignQ[term],
				AppendTo[rowBox,"-"];
				term=-term
			]
			,
			(* Later terms *)
			If[MinusSignQ[term],
				AppendTo[rowBox,"-"];
				term=-term
				,
				AppendTo[rowBox,"+"]
			];
		];
		(* Apply MakeBoxes to the term *)
		With[{x=term},
			AppendTo[rowBox, MakeBoxes[x, AuxForm]]
		]
		,
		{term,terms}
	];
	(* Display as RowBox *)
	RowBox@ rowBox
]


(* ::Subsubsection::Closed:: *)
(*Sign*)


NiceForm::sign="The sign of a term could not be detected.";


(* Returns true (false) if the argument contains a (no) Minus sign *)
MinusSignQ[arg_Times]:=Module[
	{
		factors=List@@arg,
		sign=1
	}
	,
	Do[
		(* for each negative number multiply by -1 *)
		If[NumberQ[x],
			If[Im[x]==0,
				(* for real numbers *)
				sign*=Sign[x]
				,
				(* for complex numbers *)
				If[Re[x]===0,
					(* for imaginary numbers *)
					sign*=Sign@Im[x]
					,
					(* if not zero, take the sign of the real part for complex numbers *)
					sign*=Sign@Re[x]
				]
			];
		]
		,
		{x, factors}
	];
	Switch[sign,
		1, Return[False],
		-1,Return[True],
		_,Message[NiceForm::sign]
	]
]


MinusSignQ[x:Except[Times]]:=False


(* ::Subsubsection::Closed:: *)
(*SortingFunctionPlus*)


SortPlus[terms_]:=Module[
	{
		constants={},
		tadpoles={},
		propagators={},
		marginal={},
		remainingTerms,
		maxEFTdimension= 100,
		eftOperators={}
	}
	,
	(*For the future: consider exchanging all of these special cases for a heuristic scoring function*)
	
	(* find all terms with no fields *)
	Do[
		If[CountFields[term]===0,AppendTo[constants,term]]
		,
		{term,terms}
	];
	remainingTerms=Complement[terms,constants];
	
	(* find all terms with one fields *)
	Do[
		If[CountFields[term]===1,AppendTo[tadpoles,term]]
		,
		{term,remainingTerms}
	];
	remainingTerms=Complement[remainingTerms,tadpoles];
	
	(* find all terms with two fields *)
	propagators=FindPropagators[remainingTerms];
	remainingTerms=Complement[remainingTerms,propagators];
	
	(* marginal interactions *)
	marginal=FindMarginalInteractions[remainingTerms];
	remainingTerms=Complement[remainingTerms,marginal];
	marginal=SortInteractions[marginal];
	
	(* EFT operators *)
	(* find maximum EFT dimension *)
	maxEFTdimension=If[Length[remainingTerms]==0,
		4,
		Max[MaxDim/@(remainingTerms/. Ext[x_,_]:>x/. (*OverBar*)Conj[y_]:>y)]
	];
	
	(* find EFT operators at all orders *)
	Do[
		Module[{aux=FindEFTOperators[remainingTerms,dim]},
			remainingTerms=Complement[remainingTerms,aux];
			aux=SortInteractions[aux];
			AppendTo[eftOperators,aux];
		]
		,
		{dim,5,maxEFTdimension}
	];
	eftOperators=Flatten[eftOperators];
	
	(* return the sorted list *)
	Flatten@Join[constants,tadpoles,propagators,marginal,eftOperators,remainingTerms]
]


(* ::Subsubsection::Closed:: *)
(*Ancillary functions*)


(* ::Text:: *)
(*Counts the number of covariant derivatives acting on fields*)


CountDerivatives[expr_]:=Module[
	{term=RemovePower[expr]},
	Plus@@Cases[term,Field[_,_,_,cds_List]:>Length[cds],All]+Plus@@Cases[term,FieldStrength[_,_,_,cds_List]:>Length[cds],All]
]


(* ::Text:: *)
(*Counts the number of fields and field-strength tensors*)


CountFields[expr_]:=Module[
	{term=RemovePower[expr]},
	Count[term, _Field|_FieldStrength, All]
]


(* ::Text:: *)
(*Function that returns the field type*)


GetFieldTypes[expr_]:=Module[
	{term=RemovePower[expr]},
	First[
		Cases[term,Field[_,type_,___]:>type,All],
		First[Cases[term,FieldStrength[___]:>FieldStrength,All],None]
	]
]


(* ::Text:: *)
(*Returns the field content of an expression*)


GetFieldContent[arg_]:=Module[
	{expr=RemovePower[arg]}
	,
	Flatten@Join[
		Cases[expr,Field[l_,t_,___]:>Field[l,t,___],All],
		Cases[expr,FieldStrength[l_,___]:>FieldStrength[l,___],All]
	]
]


(* ::Text:: *)
(*Returns the spin type of the fields in expr*)


GetFieldSpin[expr_]:=Flatten@Join[
	Cases[expr,Field[_,t_,___]:>t,All],
	Cases[expr,FieldStrength[l_,___]:>FieldStrength,All]
]/.Vector[_]->Vector


(* ::Text:: *)
(*Returns the maximum dimension of operators in term*)


MaxDim[term_]:=If[FreeQ[term,Plus],
	OperatorDimension[term],
	Max[OperatorDimension/@TermsToList[term]]
]


(* ::Text:: *)
(*Returns the minimum dimension of operators in term*)


MinDim[term_]:=If[FreeQ[term,Plus],
	OperatorDimension[term],
	Min[OperatorDimension/@TermsToList[term]]
]


(* ::Subsubsection::Closed:: *)
(*Find kinetic and mass terms*)


(* ::Text:: *)
(*Extracts an print-ordered list of all kinetic and mass terms from the list of terms*)


FindPropagators[expr_]:=Module[
	{
		terms=expr,
		kineticTerms={},
		remainingTerms,
		propagatorTerms
	}
	,
	
	terms= Select[terms, CountFields@ # === 2 &]; 
	
	(* Find all kinetic terms *)
	Do[
		If[KineticQ[term],
			AppendTo[kineticTerms,term]
		]
		,
		{term,terms}
	];
	(* sort the kinetic terms *)
	kineticTerms=SortBy[kineticTerms,KineticSortingFunction[#]&];
	
	(* find all remaining terms *)
	remainingTerms=Complement[terms,kineticTerms];
	(* for each kinetic term find the corresponding mass terms *)
	propagatorTerms= Table[
			FindMassTerm[term, remainingTerms] (*{kin term, massterm1, ...}*)
		, {term, kineticTerms}];
	(*Gather by identical mass terms (chiral kinetic terms can share masses). 
	Terms w/o masses should be gathered. GatherBy preserves ordering*)
	propagatorTerms= GatherBy[propagatorTerms, If[Length@ # > 1, Rest@ #, First@ #]&]; 
	propagatorTerms= Join@@(Join[#[[;;, 1]], #[[1, 2;;]]]&/@ propagatorTerms);
	
	(* Return all operators contributing to the two point function *)
	(* Add all remaining terms with two fields at the end *)
	remainingTerms= Complement[remainingTerms, propagatorTerms];
	Join[propagatorTerms, remainingTerms]
]


(* ::Text:: *)
(*Returns True if expr is a kinetic term and False otherwise*)


KineticQ[expr_]:=Module[
	{
		nFields=CountFields[expr],
		type,
		nCD
	}
	,
	If[nFields!=2,
		(* If there are not exactly 2 fields return False *)
		False
	,
		(* If there are exatly 2 fields *)
		type=GetFieldTypes[expr];
		nCD=CountDerivatives[expr];
		(* whether expr is kinetic term depends on the number of derivatives *)
		Switch[type,
			Scalar, nCD === 2,
			FieldStrength, nCD === 0,
			Fermion, nCD === 1,
			Ghost|AntiGhost, nCD === 2,
			None, False,
			_, False
		]
	]
]


(* ::Text:: *)
(*The function below might not yet find correctly the mass term for heavy vector fields*)


(* finds the mass term belonging to kinTerm in the list of terms *)
FindMassTerm[kinTerm_,terms_]:=Module[
	{
		fieldContent=GetFieldContent[kinTerm],
		massTerms={kinTerm}
	}
	,
	(* find all terms with two fields *)
	Do[
		If[CountFields[term]==2 && CountDerivatives[term]==0,
			If[MatchQ[GetFieldContent[term],fieldContent/.FieldStrength->Field],
				AppendTo[massTerms,term]
			]
		]
		,
		{term,terms}
	];
	massTerms
]


(* ::Text:: *)
(*This function assigns a weight to the kinetic term given as argument which can be used to sort all kinetic terms*)


KineticSortingFunction[term_]:=Module[
	{
		fieldContent=GetFieldContent[term],
		weight
	}
	,
	(* assign weight according to field type *)
	weight= Switch[fieldContent,
			{FieldStrength[___],FieldStrength[___]}, 1000,
			{Field[_,_Vector,___],Field[_,_Vector,___]}, 2000,
			{Field[_,Scalar,___],Field[_,Scalar,___]}, 3000,
			{Field[_,Fermion,___],Field[_,Fermion,___]}, 4000,
			{Field[_,AntiGhost,___],Field[_,Ghost,___]}, 5000,
			_, 100000
		];
	(* If both fields do not match increase weight *)
	If[!MatchQ@@fieldContent,
		weight+=900
	];

	weight
]


(* ::Text:: *)
(*Sorting interactions*)


FindMarginalInteractions[terms_]:=Module[
	{inertactions={}}
	,
	Do[
		If[MinDim[term]<=4,
			AppendTo[inertactions,term]
		]
		,
		{term,terms}
	];
	inertactions
]


(* ::Text:: *)
(*Finds all EFT operators of dimension dim in terms*)


FindEFTOperators[terms_,dim_]:=Module[
	{operators={}}
	,
	Do[
		If[MinDim[term]==dim,
			AppendTo[operators,term]
		]
		,
		{term,terms}
	];
	operators
]


(* ::Text:: *)
(*Sorts the given interactions by field content*)


SortInteractions[interactionTerms_List]:=Module[
	{
		terms,
		pureGauge={},
		pureVector={},
		pureScalar={},
		mixGaugeVector={},
		mixGaugeScalar={},
		mixVectorScalar={},
		mixGaugeVectorScalar={},
		mixFermionGauge={},
		mixFermionVector={},
		mixFermionGaugeVector={},
		mixFermionScalar={},
		mixFermionGaugeVectorScalar={},
		pureFermion={},
		remainder={},
		aux
	}
	,
	Do[
		aux=DeleteDuplicates@GetFieldSpin[term];
		Switch[aux,
			{FieldStrength},
				AppendTo[pureGauge,term],
			{Vector},
				AppendTo[pureVector,term],
			{Scalar},
				AppendTo[pureScalar,term],
			{OrderlessPatternSequence[FieldStrength,Vector]},
				AppendTo[mixGaugeVector,term],
			{OrderlessPatternSequence[FieldStrength,Scalar]},
				AppendTo[mixGaugeScalar,term],
			{OrderlessPatternSequence[Scalar,Vector]},
				AppendTo[mixVectorScalar,term],
			{OrderlessPatternSequence[FieldStrength,Vector,Scalar]},
				AppendTo[mixGaugeVectorScalar,term],
			{OrderlessPatternSequence[Fermion,FieldStrength]},
				AppendTo[mixFermionGauge,term],
			{OrderlessPatternSequence[Fermion,Vector]},
				AppendTo[mixFermionVector,term],
			{OrderlessPatternSequence[Fermion,FieldStrength,Vector]},
				AppendTo[mixFermionGaugeVector,term],
			{OrderlessPatternSequence[Fermion,Scalar]},
				AppendTo[mixFermionScalar,term],
			{OrderlessPatternSequence[Fermion,FieldStrength,Vector,Scalar]},
				AppendTo[mixFermionGaugeVectorScalar,term],
			{OrderlessPatternSequence[Fermion,FieldStrength,Scalar]},
				AppendTo[mixFermionGaugeVectorScalar,term],
			{OrderlessPatternSequence[Fermion,Vector,Scalar]},
				AppendTo[mixFermionGaugeVectorScalar,term],
			{Fermion},
				AppendTo[pureFermion,term],
			_,
				AppendTo[remainder,term]
		]
		,
		{term,interactionTerms}
	];
	(* sort each category by number of CDs then combine everything *)
	Flatten@Join[SortByCD/@{
		pureGauge,pureVector,pureScalar,
		mixGaugeVector,mixGaugeScalar,mixVectorScalar,mixGaugeVectorScalar,
		mixFermionGauge,mixFermionVector,mixFermionGaugeVector,
		mixFermionScalar,mixFermionGaugeVectorScalar,pureFermion,remainder
	}]
]


(* ::Text:: *)
(*Sorts given terms by number of derivatives*)


SortByCD[terms_List]:=Module[
	{tab,res}
	,
	(* return empty list if emtpy list is given as argument *)
	If[Length[terms]==0,Return[{}]];
	(* count CDs for each term *)
	tab=Table[
		{term,CountDerivatives[term]}
		,
		{term,terms}
	];
	(* Create empty table with length given by the maximum number of derivatives *)
	res=Table[{},{i,0,Max[tab[[;;,2]]]}];
	(* Sort terms according to number of derivatives *)
	Do[
		AppendTo[res[[Last[term]+1]],First@term]
		,
		{term,tab}
	];
	Flatten@ res
]


(* ::Subsection:: *)
(*Times*)


(* ::Subsubsection::Closed:: *)
(*AuxForm for Times*)


(* ::Text:: *)
(*Determines how products are printed *)


AuxForm/: MakeBoxes[a_Times, AuxForm]:= Module[
		{
			product= List@@ a,
			minusFlag= False,
			denomCouplings, numCouplings, otherParts, out
		}
	,
	(* check if this term is negative *)
	If[MinusSignQ[a],
		minusFlag= True;
		If[Head[-a] === Times,
			product= Flatten[{List@@(-a)}],
			product= Flatten[{-a}]
		];
	];

	(* sort the factors of Times *)
	product= SortBy[product, SortingFunctionTimes[#]&];

	(* partition into index-free coupling numerator/denominator and everything else *)
	denomCouplings= Select[product, MatchQ[#,
		Power[Coupling[_, {}, _] | Bar[Coupling[_, {}, _]], _?Negative]]&];
	numCouplings= Select[product, MatchQ[#,
		Coupling[_, {}, _] | Bar[Coupling[_, {}, _]] |
		Power[Coupling[_, {}, _] | Bar[Coupling[_, {}, _]], _?Positive]]&];
	otherParts= Select[product, !MemberQ[Join[numCouplings, denomCouplings], #]&];
	
	(*Printing the product*)
	out= If[denomCouplings === {},
			TimesRowBox@ product
		,
			(*If a denominator is present, insert the fraction in the product*)
			Block[{numBox, denomBox, fracBox, fracScore= 4.5, beforeFrac, afterFrac},
				numBox= If[numCouplings === {}, "1", TimesRowBox@ numCouplings];
				denomBox= TimesRowBox[denomCouplings /. Power[base_, n_] :> Power[base, -n]];
				fracBox= FormBox[FractionBox[numBox, denomBox], StandardForm];
				(*fracScore= Min[SortingFunctionTimes/@ Join[numCouplings, denomCouplings]];*)
				beforeFrac= Select[otherParts, SortingFunctionTimes[#] < fracScore &];
				afterFrac= Select[otherParts, SortingFunctionTimes[#] >= fracScore &];
				If[otherParts === {},
					fracBox
				,
					RowBox@ Riffle[
						Join[
							If[beforeFrac =!= {}, {TimesRowBox@ beforeFrac}, {}],
							{fracBox},
							If[afterFrac =!= {}, {TimesRowBox@ afterFrac}, {}]
						], " "
					]
				]
			]
		];
	
	(*Adding an overall minus*)
	If[minusFlag, RowBox@ {"-", out}, out]
]


(* ::Text:: *)
(*Creates a RowBox out of a list of factors *)


TimesRowBox[factors_List]:= Block[{rowBox= {}, startFlag= True},
	Do[
		If[!startFlag, AppendTo[rowBox, " "]];
		With[{x= factor},
			If[MatchQ[x, _Plus | Complex[Except[0|0.`], _]],
				AppendTo[rowBox, FormBox[
					StyleBox[RowBox@ {"(", MakeBoxes[x, AuxForm], ")"}, SpanSymmetric-> False],
					StandardForm]]
			,
				AppendTo[rowBox, MakeBoxes[x, AuxForm]]
			]
		];
		startFlag= False
	, {factor, factors}];
	RowBox@ rowBox
];


(* ::Subsubsection::Closed:: *)
(*Sorting function for times*)


(* ::Text:: *)
(*This function determines the printing order of NiceForm in each Times[...]*)
(*ToDo: Add ClebschGordan coefficients*)


SortingFunctionTimes[x_] := Switch[x,
	(* numbers *)
	_Complex, 0,
	_?NumberQ, 1,
	Power[_?NumberQ, _], 1.1,
	(* hbar *)
	hbar|\[Epsilon], 1.5,
	(* gauge-fixing *)
	_\[Xi]Gauge, 1.7,
	(* symbols *)
	Bar[_Symbol], 2,
	_Symbol, 2,
	(* FlavorSum *)
	_FlavorSum, 3.5,
	(* couplings w/o indices *)
	Bar[Coupling[_,{},_]], 4 + AlphabeticLabelOrder[First@First@x],
	Coupling[_,{},_], 4 + AlphabeticLabelOrder[First@x],
	(* couplings w/ indices *)
	Bar[Coupling[_,{__},_]], 6 + AlphabeticLabelOrder[First@First@x],
	Coupling[_,{__},_], 6 + AlphabeticLabelOrder[First@x],
	(* Plus *)
	_Plus, Max[SortingFunctionTimes /@ List @@ x] + 1,
	_Times, Max[SortingFunctionTimes /@ List @@ x], (*Can happen in nested expressions*)
	(* loop functions *)
	_LF|_Log, 10,
	(* Delta, Metric, CG *)
	_Delta, 11,
	_Metric, 12,
	_CG, 13,
	(* Dirac product *)
	_DiracProduct, 14,
	(* external momenta *)
	_ExtMom, 15,
	(* - - - - - - fields - - - - - - *)
	(* scalars *)
	Bar@Field[_,Scalar,___], 20 + AlphabeticLabelOrder[First@First@x],
	Field[_,Scalar,___], 20 + AlphabeticLabelOrder[First@x],
	(* Vectors *)
	Bar@Field[_,_Vector,___], 22 + AlphabeticLabelOrder[First@First@x],
	Field[_,_Vector,___], 22 + AlphabeticLabelOrder[First@x],
	(* Field-Strength *)
	Bar[_FieldStrength], 24 + AlphabeticLabelOrder[First@First@x],
	_FieldStrength, 24 + AlphabeticLabelOrder[First@x],
	(* fermion spin-chains *)
	_NCM, 26,
	(* operators *)
	_Operator, 30,
	_AtomicOp, 30,
	_CompOp, 30,
	(* Power: score of base + 1 *)
	Power[_, _], SortingFunctionTimes[First@ x] + 1,
	(* unknown *)
	_, 100000
]


(* ::Subsubsection::Closed:: *)
(*Alphabetic order of a label*)


(* ::Text:: *)
(*This auxiliary function gives a number in [0,1] determining the alphabetic order of the given argument*)


AlphabeticLabelOrder[x_] := Module[
	{str=ToString[x],numbers,max=4}
	,
	(* transform to characters to numbers *)
	numbers = ToString/@ToCharacterCode[str];
	(* fill missing 0's such that each number has 'max=4' digits *)
	numbers = Table[
		If[(max-StringLength[number])>0,
			StringJoin@Flatten@Join[{number},ConstantArray["0",max-StringLength[number]]],
			number
		]
		,
		{number,numbers}
	];
	(* make decimal *)
	PrependTo[numbers,"0."];
	ToExpression@StringJoin[numbers]
]


(* ::Text:: *)
(*Options for StyleBox to change the Parenthesis.*)


(*SpanMaxSize*)(*SpanAdjustments*)(*SpanSymmetric*)


(* ::Section:: *)
(*Specific printing definitions*)


(* ::Text:: *)
(*Global association of NiceForm printing for fields and couplings.*)


$LabelsNiceForm= <|Coupling-> <||>, Field-> <||>|>;


(* ::Subsubsection::Closed:: *)
(*Clebsch Gordan coefficients*)


AuxForm/: MakeBoxes[CG[label_, indices_], AuxForm]:= UpDownIndices[CGLabel@ label, indices];


CGLabel[Bar@ _del]= "\[Delta]";
CGLabel[Bar@ label_]:= Overline@ CGLabel@ label;
CGLabel[_eps]= "\[CurlyEpsilon]";
CGLabel[_del]= "\[Delta]";
CGLabel[_dSym]= "d";
CGLabel[_fStruct]= "f";
CGLabel[_gen]= "T";
CGLabel[label_]:= ToString@ label;


(* ::Subsubsection::Closed:: *)
(*Coupling*)


AuxForm/: MakeBoxes[Bar@ Coupling[label_, indices_, _], AuxForm]:= 
	UpDownIndices[Overline[label/. $LabelsNiceForm[Coupling]], Bar/@indices];


AuxForm/: MakeBoxes[Coupling[label_, indices_, _], AuxForm]:= 
	UpDownIndices[label/. $LabelsNiceForm[Coupling], indices];


(* ::Subsubsection::Closed:: *)
(*Constants*)


AuxForm/: MakeBoxes[\[Mu]bar2, AuxForm]:= SuperscriptBox[Overline@ "\[Mu]", 2];
AuxForm/: MakeBoxes[hbar, AuxForm]:= "\[HBar]";


AuxForm/: MakeBoxes[\[Xi]Gauge[fieldLabel_], AuxForm]:= 
	SubscriptBox["\[Xi]", SubscriptStyle2@ ToString[fieldLabel/. $LabelsNiceForm[Field]]];


(* ::Subsubsection::Closed:: *)
(*DiracAlgebra*)


AuxForm/: MakeBoxes[Gamma5, AuxForm]:= SubscriptBox["\[Gamma]", "5"];
AuxForm/: MakeBoxes[GammaM@ mu_, AuxForm]:= DownIndices["\[Gamma]", {mu}];
AuxForm/: MakeBoxes[GammaM[mu_, nu__], AuxForm]:= DownIndices["\[CapitalGamma]", {mu, nu}];


AuxForm/: MakeBoxes[Proj[1], AuxForm]:= SubscriptBox["P", "R"];
AuxForm/: MakeBoxes[Proj[-1], AuxForm]:= SubscriptBox["P", "L"];
AuxForm/: MakeBoxes[GammaCC, AuxForm]:= "C";


(* ::Subsubsection::Closed:: *)
(*Evanescent operator*)


AuxForm/: MakeBoxes[EvaOp[class_, id_, inds_List], AuxForm]:= 
	SubsuperscriptBox["E" <> ToString@ id, 
		OpClassFormatting@ class,
		SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, inds], "RowDefault"]];


(* ::Subsubsection::Closed:: *)
(*Fields*)


$ExtColor= Darker[Gray,0.1]; (* this is the color in which Ext[] is printed *)


(* ::Text:: *)
(*Fields*)


AuxForm/: MakeBoxes[Field[label_:Except[_Ext], type_, indices_, CDerivs_], AuxForm]:= Switch[type
	,_Vector,
		MakeFieldBox[label/. $LabelsNiceForm[Field], Append[indices, First@ type], CDerivs]
	,_ ,
		MakeFieldBox[label/. $LabelsNiceForm[Field], indices, CDerivs]
	];
	
AuxForm/: MakeBoxes[Field[label_Ext, type_, indices_, CDerivs_], AuxForm]:= Switch[type
	,_Vector,
		StyleBox[SubscriptBox[RowBox@{"\[LeftAngleBracket]",MakeFieldBox[label/. $LabelsNiceForm[Field], Append[indices, First@ type], CDerivs],"\[RightAngleBracket]"},label[[2]]],FontColor->$ExtColor]
	,_ ,
		StyleBox[SubscriptBox[RowBox@{"\[LeftAngleBracket]",MakeFieldBox[label/. $LabelsNiceForm[Field], indices, CDerivs],"\[RightAngleBracket]"},label[[2]]],FontColor->$ExtColor]
	];


AuxForm/: MakeBoxes[Bar@ Field[label_:Except[_Ext], type_, indices_, CDerivs_], AuxForm]:= Switch[type
	,_Vector,
		MakeFieldBox[Overline[label/. $LabelsNiceForm[Field]], Append[Bar/@ indices, First@ type], CDerivs]
	,_ ,
		MakeFieldBox[Overline[label/. $LabelsNiceForm[Field]], Bar/@ indices, CDerivs]
	];

AuxForm/: MakeBoxes[Bar@ Field[label_Ext, type_, indices_, CDerivs_], AuxForm]:= Switch[type
	,_Vector,
		StyleBox[SubscriptBox[RowBox@{"\[LeftAngleBracket]",MakeFieldBox[Overline[label/. $LabelsNiceForm[Field]], Append[Bar/@ indices, First@ type], CDerivs],"\[RightAngleBracket]"},label[[2]]],FontColor->$ExtColor]
	,_ ,
		StyleBox[SubscriptBox[RowBox@{"\[LeftAngleBracket]",MakeFieldBox[Overline[label/. $LabelsNiceForm[Field]], Bar/@ indices, CDerivs],"\[RightAngleBracket]"},label[[2]]],FontColor->$ExtColor]
	];


(* ::Text:: *)
(*Field strengths*)


AuxForm/: MakeBoxes[FieldStrength[label_, lorInds_, indices_, CDerivs_], AuxForm]:= 
	MakeFieldBox[label/. $LabelsNiceForm[Field], Join[lorInds, indices], CDerivs];
AuxForm/: MakeBoxes[HoldPattern@ Bar@ FieldStrength[label_, lorInds_, indices_, CDerivs_], AuxForm]:= 
	MakeFieldBox[Overline@ ToString[label/. $LabelsNiceForm[Field]], Join[lorInds, Bar/@ indices], CDerivs];


(* ::Text:: *)
(*Formatting functions*)


MakeFieldBox[lab_, indices_, CDinds_]:= CDFormatting[UpDownIndices[lab, indices], CDinds];


CDFormatting[lab_, {}]:= lab;
CDFormatting[lab_, CDinds:{__}]:= TemplateBox[Append[FormatDerivList@ CDinds, lab], "RowDefault"];
CDFormatting[lab_, other_]:= TemplateBox[{SubscriptBox["D", MakeBoxes[other, AuxForm]], lab}, "RowDefault"];


FormatDerivList[inds_]:= Module[{tagged, collapsed},
	tagged= inds/. {i_Index:> {"D", i}, PartialDev[i_]:> {"\[PartialD]", i}};
	collapsed= tagged//. {
		{x___, {s_, a_}, {s_, a_}, y___}/; !IntegerQ[a]:> {x, {s, 2}, y},
		{x___, {s_, n_Integer}, {s_, m_Integer}, y___}:> {x, {s, n+m}, y}};
	Replace[collapsed, {
			{s_, n_Integer}:> SuperscriptBox[s, n],
			{s_, i_}:> SubscriptBox[s, MakeBoxes[i, AuxForm]]
		}, {1}]
];


(* ::Text:: *)
(*Powers of fields/field strengths*)


AuxForm/: MakeBoxes[Power[f:(Field|FieldStrength)[__, derivs_], n_], AuxForm]:= If[derivs==={},
	SuperscriptBox[MakeBoxes[f, AuxForm], n],
	(*FormBox[..., StandardForm] ensures large parenthesis*)
	(*StyleBox[..., SpanSymmetric-> False] ensures that the parenthesis do not reach further down that is necessary*)
	SuperscriptBox[FormBox[StyleBox[
			RowBox@ {"(", MakeBoxes[f, AuxForm], ")"},
			SpanSymmetric-> False
		], StandardForm], n]
]


(* ::Subsubsection::Closed:: *)
(*FlavorSums*)


AuxForm/: MakeBoxes[FlavorSum[ind_Index], AuxForm]:= UnderscriptBox[StyleBox["\[Sum]", 18], SubscriptStyle2@ MakeBoxes[ind, NiceForm]]


(* ::Subsubsection::Closed:: *)
(*Logarithms*)


AuxForm/: MakeBoxes[Log[Times[\[Mu]bar2, Power[a_, -2]]], AuxForm]:= FormBox[RowBox@ {
			Log, "[",
				FractionBox[
					MakeBoxes[\[Mu]bar2, AuxForm],
					MakeBoxes[Power[a, 2], AuxForm]
				], "]"
	}, StandardForm];


AuxForm/: MakeBoxes[Log[Times[Power[a_,2],Power[b_,-2]]], AuxForm]:=FormBox[RowBox@ {
			Log,"[",
				FractionBox[
					MakeBoxes[Power[a,2], AuxForm],
					MakeBoxes[Power[b,2], AuxForm]
				], "]"
	}, StandardForm]


(* ::Subsubsection::Closed:: *)
(*Loops*)


AuxForm/: MakeBoxes[LF[masses_, powers_], AuxForm]:= 
	RowBox@{SubscriptBox["LF", RowBox@ Riffle[powers,","]], "[", Sequence@@ Riffle[(MakeBoxes[#, AuxForm]&/@ masses), ","], "]"};


(* ::Subsubsection::Closed:: *)
(*Metric, Delta & LCTensor*)


AuxForm/: MakeBoxes[Delta@ ind__, AuxForm]:= DownIndices["\[Delta]", {ind}];


AuxForm/: MakeBoxes[LCTensor@ ind__, AuxForm]:= DownIndices["\[CurlyEpsilon]", {ind}];


AuxForm/: MakeBoxes[Metric@ ind__, AuxForm]:= DownIndices["g", {ind}];


(* ::Subsubsection::Closed:: *)
(*Numbers*)


AuxForm/: MakeBoxes[z_Complex, AuxForm]:= MakeBoxes[z, StandardForm]


(* ::Subsubsection::Closed:: *)
(*Open covariant derivatives*)


AuxForm/: MakeBoxes[OpenCD[inds:{Index[_, Lorentz]..}], AuxForm]:= FormBox[StyleBox[
		RowBox[Flatten@ {"(", Replace[inds//. {x___, a_, a_, y___}:> {x, 2, y}, 
	{2:> SuperscriptBox[OverscriptBox["D", "\[LongRightArrow]"], 2], i_:> SubscriptBox[OverscriptBox["D", "\[LongRightArrow]"], MakeBoxes[i, AuxForm]]}, {1}],
	")"} ]
	, SpanSymmetric-> False], StandardForm]


(* ::Subsubsection::Closed:: *)
(*External momenta*)


AuxForm/: MakeBoxes[ExtMom[ind_Index, Ext[lab_, no_]], AuxForm]:= FormBox[
	SubsuperscriptBox["p", MakeBoxes[ind, AuxForm], no (*SubscriptBox[ToBoxes[lab/. $LabelsNiceForm[Field]], no]*)] (* use ToBoxes here since MakeBoxes has attribute HoldComplete*)
, StandardForm]


(* ::Subsubsection::Closed:: *)
(*Simplifications*)


(* ::Text:: *)
(*Print + H.c. for HcTerms *)


AuxForm/: MakeBoxes[HcTerms@ expr_, AuxForm]:= 
	FormBox[RowBox[{"[", MakeBoxes[expr, AuxForm], "+", StyleBox["H.c.", FontSize-> 14], "]"}]/. 
		RowBox@ {x___,  RowBox@ {a__}, y___}:> RowBox@ {x, a, y}, StandardForm];


(* ::Subsubsection::Closed:: *)
(*Finite Temperature*)


AuxForm/: MakeBoxes[\[CapitalLambda]\[CapitalTau], AuxForm]:= SubscriptBox["\[CapitalLambda]","\[CapitalTau]"];
AuxForm/: MakeBoxes[\[CapitalLambda]bar2, AuxForm]:= SuperscriptBox[Overline@ "\[CapitalLambda]", 2];


AuxForm/: MakeBoxes[Zeta[x_], AuxForm]:= SubscriptBox["\[Zeta]",ToString@x];


AuxForm/: MakeBoxes[TemporalProjector[x___], AuxForm]:= SuperscriptBox["\!\(\*SubscriptBox[\(\[DoubleStruckCapitalP]\), \(0\)]\)",MakeBoxes[{x},AuxForm]];


AuxForm/: MakeBoxes[Thermalhbar, AuxForm]:= "\!\(\*SubscriptBox[\(\[HBar]\), \(\[CapitalTau]\)]\)";


AuxForm/: MakeBoxes[Log[Times[c___, \[CapitalLambda]bar2, Power[a_, -2]]], AuxForm]:= FormBox[RowBox@ {
			Log, "[",
				FractionBox[
					If[Length@{c}===0, MakeBoxes[\[CapitalLambda]bar2, AuxForm], MakeBoxes[Times@c \[CapitalLambda]bar2, AuxForm]],
					MakeBoxes[Power[a, 2], AuxForm]
					(*If[Length@{c}===0, MakeBoxes[Power[a, 2], AuxForm], MakeBoxes[1/Times@c Power[a, 2], AuxForm]]*)
				], "]"
	}, StandardForm];


(* ::Subsubsection::Closed:: *)
(*Conj*)


AuxForm/: MakeBoxes[Conj[arg_], AuxForm]:= 
	SuperscriptBox[MakeBoxes[arg, AuxForm], "*"];


(* ::Subsection:: *)
(*Debugging  *)


(* ::Subsubsection::Closed:: *)
(*Operator*)


(* ::Text:: *)
(*Ignore "Operator"*)


AuxForm/: MakeBoxes[HoldPattern@ Operator@ op__, AuxForm]:= 
	With[{temp= Times@ op}, MakeBoxes[temp, AuxForm] ];


(* ::Subsubsection::Closed:: *)
(*Atomic operators*)


(* ::Text:: *)
(*AtomicOp*)


AuxForm/: MakeBoxes[AtomicOp[class_, id_, inds_List], AuxForm]:= 
	SubsuperscriptBox["\[ScriptCapitalA]" <> ToString@ id, 
		OpClassFormatting@ class,
		SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, inds], "RowDefault"]];


AuxForm/: MakeBoxes[AtomicOp[class_, id_, inds_Pattern], AuxForm]:= 
	SubsuperscriptBox["\[ScriptCapitalA]" <> ToString@ id, OpClassFormatting@ class, SubscriptStyle2@ MakeBoxes[inds, StandardForm]];


(* ::Text:: *)
(*CompOp*)


AuxForm/: MakeBoxes[CompOp[idenType_, class_, id_, inds_List], AuxForm]:= 
	SubsuperscriptBox["\[ScriptCapitalC]" <> ToString@ id, 
		MapAt[Prepend[$identityFormat@ idenType <> ","], OpClassFormatting@ class, 1],
		SubscriptStyle2@ TemplateBox[Map[MakeBoxes[#, AuxForm]&, inds], "RowDefault"]];


$identityFormat= <|Evanescent-> "ev", FourDimensional-> "4", dDimensional-> "\[ScriptD]"|>;


(* ::Text:: *)
(*Format the operator class*)


OpClassFormatting@ {fields_List, 0}:= FieldListFormatting@ fields;
OpClassFormatting@ {fields_List, devs_Integer}:= Block[{out}, 
	out= FieldListFormatting@ fields;
	MapAt[Prepend[If[devs > 1,
			SuperscriptBox["\[PartialD]", ToString@ devs]
		,
			"\[PartialD]"
		]], out, 1]
];


FieldListFormatting@ fields_List:= Block[{tally, label},
	tally= Tally@ fields;
	TemplateBox[Table[
		label= If[Head@ f[[1]] === Conj, 
			Overline@ f[[1, 1]]
		,
			ToString@ f[[1]]
		];
		If[f[[2]]> 1,
			SuperscriptBox[label, ToString@ f[[2]]]
		,
			label
		]
	, {f, tally}], "RowDefault"]
]
