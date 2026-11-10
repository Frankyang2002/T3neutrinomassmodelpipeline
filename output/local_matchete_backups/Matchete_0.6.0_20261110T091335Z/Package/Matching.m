(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Matching`*)


(* ::Subtitle:: *)
(*Various modules for defining the EFT power counting and matching.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["Match"]
PackageExport["CovariantLoop"]


(* ::Text:: *)
(*Options *)


PackageExport["LoopOrder"]
PackageExport["Matching"]
PackageExport["Divergence"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["$FieldTypes"]
PackageScope["hScalar"]
PackageScope["lScalar"]
PackageScope["hFermion"]
PackageScope["lFermion"]
PackageScope["hVector"]
PackageScope["lVector"]
PackageScope["hGhost"]
PackageScope["lGhost"]
PackageScope["hAntiGhost"]
PackageScope["lAntiGhost"]
PackageScope["fieldFormat"]


PackageScope["ResetAuxiliaryMatchingInformation"]
PackageScope["ResetSubstitutionInformation"]
PackageScope["$currentEOMs"]
PackageScope["$currentHeavyDims"]
PackageScope["$Xsubs"]
PackageScope["$XsubsMom"]
PackageScope["$LoopMomOrders"]
PackageScope["$Msubs"]
PackageScope["$Gsubs"]
PackageScope["$XFieldDofs"]
PackageScope["$XOrders"]
PackageScope["$XOrdMin"]
PackageScope["$EvaOrdMin"]
PackageScope["$EvaOrders"]
PackageScope["$EvaSubs"]
(*PackageScope["$currentXdims"]*)


PackageScope["SetCurrentLagrangian"]
PackageScope["$currentLagrangian"]


PackageScope["ListPowerTypeTraces"]
PackageScope["MatchReduce"]


PackageScope["FieldDoFs"]


PackageScope["LoopMatch"]


PackageScope["$ValidationRun"]


PackageScope["MassiveVectorsInLag"]


(* ::Text:: *)
(*Options *)


PackageScope["Simplifications"]


PackageScope["WhichTraces"]


(* ::Text:: *)
(*Dimensional Reduction*)


PackageScope["DimReduction"]


(* ::Section:: *)
(*Usage  messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


Match::usage=
	"Match[Lag, EFTOrder -> 6, LoopOrder-> 1] returns the EFT Lagrangian resulting from integrating out the heavy fields from the Lagrangian Lag. The optional arguments EFTOrder and LoopOrder determine, respectively, the desired order in the EFT and loop expansion, with LoopOrder admitting the values 0 (tree), {1} (loop) or 1 (tree+loop).";


LoopOrder::usage = "LoopOrder is an option for Match routine specifying the targeted order in the loop expansion. It can take the values 0 (tree), {1} (loop) or 1 (tree+loop).";


CovariantLoop::usage=
	"CovariantLoop[Lag, {fields}] returns the value of the supertraces involving exactly the degrees of freedom specified in the field list. It takes the option EFTOrder.";


Matching::usage= "Option value used in Matchete."
Divergence::usage= "Option value used in Matchete."


(* ::Subsubsection::Closed:: *)
(*Internal*)


IsolateMassTerms::usage="IsolateMassTerms[\[ScriptCapitalL], Heavy -> All] returns all massterms of the Lagrangian \[ScriptCapitalL]. Replacing All by True or False only returns the heavy or light mass terms respectively.";


ListPowerTypeTraces::usage= "ListPowerTypeTrace[order, <lightFieldsOnly>] returns a list of all the propagator set that can contribute in the supertraces up to a given order and given a Lagrangian.";


LoopMatch::usage=
	"LoopMatch[Lag] returns \!\(\*SubsuperscriptBox[\(L\), \(EFT\), \((1)\)]\) given a UV Lagrangian as input.";


$ValidationRun::usage= "Global flag for vor validation mode."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*EFT power counting *)


(* ::Subsubsection::Closed:: *)
(*Global flag for the mass dimension*)


$currentEFTOrder = 6;


(* ::Section:: *)
(*Extract theory information*)


(* ::Subsubsection::Closed:: *)
(*Field type info*)


(* ::Text:: *)
(*Properties of the field types propagating in the supertraces*)


$FieldTypes= <|
	hScalar-> <|Type-> Scalar, Heavy-> True|>,
	lScalar-> <|Type-> Scalar, Heavy-> False|>,
	hFermion-> <|Type-> Fermion, Heavy-> True|>,
	lFermion-> <|Type-> Fermion, Heavy-> False|>,
	hVector-> <|Type-> Vector, Heavy-> True|>,
	lVector-> <|Type-> Vector, Heavy-> False|>,
	hGhost-> <|Type-> Ghost, Heavy-> True|>,
	lGhost-> <|Type-> Ghost, Heavy-> False|>,
	hAntiGhost-> <|Type-> AntiGhost, Heavy-> True|>,
	lAntiGhost-> <|Type-> AntiGhost, Heavy-> False|>
|>;


(* ::Subsection:: *)
(*Substitution rules *)


(* ::Text:: *)
(*To determine substitution rules given Lagrangian input *)


(* ::Subsubsection::Closed:: *)
(*Reset substitution information*)


(* ::Text:: *)
(*Resets the information derived in SetSubstitutions*)


ResetSubstitutionInformation[]:= Block[{},
	{$XOrders, $XOrdMin, $Xsubs, $Msubs, $Gsubs}= {<||>, <||>, {}, {}, {}};
	{$LoopMomOrders, $XsubsMom}= {<||>, {}};
	{$EvaOrders, $EvaOrdMin, $EvaSubs}= {<||>, <||>, {}};
	$XFieldDofs= {};
]


(* ::Subsubsection::Closed:: *)
(*Derive Xterm and Mterm substitutions as well their EFT orders*)


(* ::Text:: *)
(*lag: UV Lagrangian*)
(*	Compute fluctuation operator from given lag.*)
(*	Sets the global flags {$XFieldDofs, $XSubs, $XOrders, $XOrdMin, $MSubs}.*)


Options[SetSubstitutions] = {EFTOrder -> 6, Mode -> Matching};


SetSubstitutions[lag_, OptionsPattern[]] := Module[
	{tmp, xTemp, terms, evaTemp, fieldDofs, Xterms,effLag,allFieldLabels, gfTerms, \[Alpha], \[Beta], i, j, iPattern,jPattern,
		devs,allEFTorders,l1,l2,f1,f2, minXOrderByFieldType, xSubsTMP,xSubsMomTMP,RelabelI,
		minEvaOrderByFieldType, orders, \[Lambda], \[Kappa], nMom, nOpen}
	,
	(*Reset all substitution information*)
	ResetSubstitutionInformation[];

	(*Determine fields in lag sorted by classes, accounting for complex dofs as in [2012.08506, (2.33)]*)
	$XFieldDofs= LagrangianDofs@ lag;

	(*prepare patterns corresponding to indices i and j*)
	iPattern=ToExpression[ToString[i]<>"_"];
	jPattern=ToExpression[ToString[j]<>"_"];

(* Determine field content *)
	(*Association: Keys = field classes, Values =  field variable functions*)
	fieldDofs=(FieldDoFs/@ $XFieldDofs);
	(*List of all field labels in lag*)
	allFieldLabels=DeleteCases[Flatten@ Values@ $XFieldDofs, _Conj];

(*Subtract kinetic terms that go into propagators*)
	effLag= lag-KinOpLagrangian@@ allFieldLabels;

(*Determine X-terms*)
	(*Iterate over field types*)
	Do[
		(*save all EFT orders for this type*)
		{minXOrderByFieldType, minEvaOrderByFieldType}={{}, {}};

		(*Iterate over fields of given type*)
		Table[
			{l1,f1,l2,f2}=Flatten@{Keys[field1],Values[field1],Keys[field2],Values[field2]};
			(*perform functional derivatives*)
			xTemp=-FluctuationOperator[effLag,Bar@f1@i,f2@j,EFTOrder->OptionValue[EFTOrder]];
			xTemp= xTemp/. {Field[_, Ghost, __]-> 0, Field[_, AntiGhost, __]-> 0, Field[_, _Vector, __] -> 0};
			
			If[OptionValue@Mode===DimReduction, xTemp= xTemp/. Field[_,Fermion,__]->0];

			(*transform to momentum space: OpenCD[a]->OpenCD[a]-\[ImaginaryI]*LoopMom[a]*)
			xTemp=Expand[
				xTemp/.OpenCD[inds_List]:>FuncNCM@@Table[OpenCD[{\[Mu]}]-I*LoopMom[\[Mu]],{\[Mu],inds}]
			];
			
			(*X-terms*)
			(*Find all EFTOrder*)
			allEFTorders=GetAllOperatorDimensions[xTemp];
			minXOrderByFieldType= Join[minXOrderByFieldType, allEFTorders];
			If[OptionValue@ Mode === Evanescent,
				minEvaOrderByFieldType= Join[minEvaOrderByFieldType, GetAllOperatorDimensions@ Coefficient[xTemp, ev]];
			];
			
			(*Write replacements rules order-by-order*)
			(*XSubs for matching*)
			If[OptionValue@ Mode === Matching||OptionValue@ Mode === DimReduction,
				{xSubsTMP, orders}= Reap@ Flatten@ Table[
						terms= SeriesEFT[xTemp, EFTOrder-> {ord}]/. 
							{lmom_LoopMom-> \[Lambda]* lmom, Prop[0]->\[Lambda]^(-2) Prop[0], o_OpenCD-> \[Lambda]* \[Kappa]* o};
						terms= CoefficientList[terms, {\[Lambda], \[Kappa]}, (Exponent[terms, \[Lambda]]+ 1){1, 1}];
						Table[
							If[terms[[nMom, nOpen]] === 0, 
								Nothing
							,
								(*Populate order (EFT-base order, # cov. mom. operators)*)
								Sow@ {ord, nMom -1};
								(*X-terms are organized as 
									Xterm[labels, ind. patterns, base EFT order, #Momentum operators, #open CDs]*)
								Xterm[{l1,l2}, {iPattern, jPattern}, ord, nMom- 1, nOpen- 1]:>
									Evaluate@RelabelI[terms[[nMom, nOpen]], Unique-> True]
							]
						, {nMom, Length@ terms}, {nOpen, nMom}]
					, {ord, DeleteCases[allEFTorders, n_/;n>99]} ];
				
				xSubsTMP= xSubsTMP/. {RelabelI[0, _]-> 0, RelabelI-> RelabelIndices};
				If[Length@ xSubsTMP > 0,
					AppendTo[$XOrders, {l1,l2}-> DeleteDuplicates@ First@ orders];
				,
					AppendTo[$XOrders, {l1,l2}-> {{100, 0}}];
				];
			,
				(*Xorders for the divergent calculations are ordered by *)
				{xSubsTMP, orders}= Reap@ Flatten@ Table[
						terms= SeriesEFT[xTemp, EFTOrder-> {ord}]/. 
							{lmom_LoopMom-> \[Lambda]* lmom, Prop[0]->\[Lambda]^(-2) Prop[0], o_OpenCD-> \[Lambda]* \[Kappa]* o};
						terms= CoefficientList[terms, {\[Lambda], \[Kappa], ev}, {Exponent[terms, \[Lambda]]+ 1, Exponent[terms, \[Lambda]]+ 1, 2}];
						Table[{
							If[terms[[nMom, nOpen, 1]] === 0, 
								Nothing
							,
								(*Populate order (EFT-base order, # cov. mom. operators, ev. True/False)*)
								Sow@ {ord, nMom -1, False};
								Xterm[{l1,l2}, {iPattern, jPattern}, ord, nMom- 1, nOpen- 1, False]:>
									Evaluate@RelabelI[terms[[nMom, nOpen, 1]], Unique-> True]
							]
						,
							If[OptionValue@ Mode =!= Evanescent || terms[[nMom, nOpen, 2]] === 0, 
								Nothing
							,
								(*Populate order (EFT-base order, # cov. mom. operators, ev. True/False)*)
								Sow@ {ord, nMom -1, True};
								Xterm[{l1,l2}, {iPattern, jPattern}, ord, nMom- 1, nOpen- 1, True]:>
									Evaluate@RelabelI[terms[[nMom, nOpen, 2]], Unique-> True]
							]
							(*X-terms are organized as Xterm[labels, ind. patterns, 
								base EFT order, #Momentum operators, #open CDs, ev.]*)
						}, {nMom, Length@ terms}, {nOpen, nMom}]
					, {ord, DeleteCases[allEFTorders, n_/;n>99]} ];
				xSubsTMP= xSubsTMP/. {RelabelI[0, _]-> 0, RelabelI-> RelabelIndices};
				
				If[Length@ xSubsTMP > 0,
					AppendTo[$XOrders, {l1,l2}-> DeleteDuplicates@ First@ orders];
				,
					AppendTo[$XOrders, {l1,l2}-> {{100, 0, False}}];
				]
			];
			$Xsubs= Join[$Xsubs, xSubsTMP];
		, {field1,fieldDofs@fieldType1},{field2,fieldDofs@fieldType2} ];

		(*determine min. EFT order for field type*)
		If[minXOrderByFieldType==={},
			minXOrderByFieldType=100,
			minXOrderByFieldType=Min@minXOrderByFieldType
		];
		AppendTo[$XOrdMin, {fieldType1,fieldType2}-> minXOrderByFieldType];

		If[OptionValue@ Mode === Evanescent,
			If[minEvaOrderByFieldType === {},
				minEvaOrderByFieldType= 100
			,
				minEvaOrderByFieldType= Min@ minEvaOrderByFieldType
			];
			AppendTo[$EvaOrdMin, {fieldType1, fieldType2}-> minEvaOrderByFieldType];
		]
	,
	{fieldType1, Keys@ $XFieldDofs}, {fieldType2, Keys@ $XFieldDofs}
	];

	(*Test if the X terms are consistent with the EFT expansion*)
	tmp= MinimalBy[$XOrders, Min@ #[[;;,1]]&];
	If[Min[First/@ First/@ tmp] <= 0,
		Message[SetCurrentLagrangian::Xdims, First@ Keys@ tmp];
		Abort[];
	];

(* Determine Mass substitutions *)
	$Msubs= Flatten@ Table[
			Table[
				tmp= GetFieldsUpdated[f, Mass];
				(*give indices to masses if necessary*)
				If[Length@ $CouplingAssociation[tmp, Indices]=== 1,
					Mterm[f,iPattern]->tmp[i]
				,
					Mterm[f,___]->tmp[]
				]
			,
			{f, DeleteDuplicates[$XFieldDofs[fieldType]/. Conj[x_]->x]}
			]
		,
		{fieldType,{hScalar, hFermion, hVector, hGhost, hAntiGhost}} (*NB: Loop only over heavy DoFs*)
		];

(* Determine gauge coupling substitutions *)
	$Gsubs= Table[
			GaugeCTerm@ f-> BetterSeries[1/GetFieldsUpdated[f, Coupling], {hbar, 0, 2}]
		, {f, Intersection[$XFieldDofs@ lVector, List@@ Query[All, Key@ Field]@ $GaugeGroups]}]; (*NB: Loop over gauge fields*)
];




(* ::Subsubsection::Closed:: *)
(*Lagrangian DoFs*)


(* ::Text:: *)
(*Gives an association with the DoFs of a Lagrangian for each field type .*)


LagrangianDofs[lag_]:= LagrangianDofs[lag]= Module[{fields},
	(*Get all fields from Lagrangian*)
	fields= Complement[OccurringFields@ lag, FieldByProperty[BackgroundField-> True]];
	(*Create association by field type*)
	fields= Association@@
		KeyValueMap[(#1-> Intersection[GetFieldsUpdatedByProperty[#2], fields]&), $FieldTypes];
	(*Add conjugate DoFs*)
	LagrangianDofsAux/@ fields
 ];


(* ::Text:: *)
(*Extends list of fields with conjugate DoFs if they are not selfconjugate *)


LagrangianDofsAux@ fieldList_List:= LagrangianDofsAux/@ fieldList// Flatten;
LagrangianDofsAux@ field_Symbol:= If[$FieldAssociation[field, SelfConjugate], field, {field, Conj@ field}];


(* ::Text:: *)
(*"Conj" is used to describe the conjugated DoF.*)


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


(* ::Text:: *)
(*Auxiliary function for constructing the relevant DoFs of a field with indices to go in functional derivatives*)


FieldDoFs[fields_List]:= (#->FieldDoFs[#])&/@ fields;
FieldDoFs[Conj@ f_Symbol]:= FieldDoFs[f, True];
FieldDoFs[f_Symbol]:= FieldDoFs[f, False];
FieldDoFs[f_Symbol, conj_]:= Block[{props= $FieldAssociation[f], inds, i},
	inds= Sequence@@ ConstantArray[i, Length@ props@ Indices + If[$FieldAssociation[f, Type]=== Vector, 1, 0]];
	ReplacePart[Switch[{props@ Type, conj}
		,{Fermion, False},
			Function[i, Evaluate@ f@ inds]
		,{Fermion, True},
			Function[i, Evaluate@ CConj@ f@ inds]
		,{_, False},
			Function[i, Evaluate@ f@ inds]
		,{_, True},
			Function[i, Evaluate@ Bar@ f@ inds]
	], {1}-> i]
];


(* ::Text:: *)
(*Auxiliary function for determining the maximum EFT order of all terms in a sum/list*)


AllOperatorDimension[0]={0};
AllOperatorDimension[0.]={0};
AllOperatorDimension[{}]={0};
AllOperatorDimension[expr_Plus]:=Flatten[ AllOperatorDimension/@ (List@@ expr)];
AllOperatorDimension[expr_List]:= Flatten[AllOperatorDimension/@ expr];

GetAllOperatorDimensions[term_]:=Module[{tmp=BetterExpand[term]},
	tmp=AllOperatorDimension[tmp];
	tmp=DeleteDuplicates@Flatten@{tmp/.AllOperatorDimension->OperatorDimension};
	tmp/.{0}->{100}
]


(* ::Text:: *)
(*Returns a list of the labels of all fields occurring in the Lagrangian*)


OccurringFields[L_]:= DeleteDuplicates@ Cases[L, (Field|FieldStrength|WilsonTerm)[lab_, __]-> (lab/.Conj->Identity), All];


(* ::Subsection:: *)
(*Update the current Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*SetCurrentLagrangian*)


(* ::Text:: *)
(*Initialize loop and EFT order of the current Lagrangian*)


$currentLoopOrder = -1;


ResetAuxiliaryMatchingInformation[]:= Block[{},
	$currentEOMs = <||>;
	$currentFieldAssociation = <||>;
	$currentHeavyDims = <||>;
]
ResetAuxiliaryMatchingInformation[];


(* ::Text:: *)
(*Sets the global state for the the current Lagrangian and associated objects*)


Options@ SetCurrentLagrangian= {Verbose-> True, Mode-> Matching};


SetCurrentLagrangian::Xdims= "The \"`1`\" X-term was found to have EFT-order 0. This is not gonna turn out well!"


SetCurrentLagrangian[lagrangian_, loopOrder_, eftOrder_, OptionsPattern[]] := Module[
	{eftOrd, newLag, lag = HcExpand@ lagrangian, lagFixed, gfTerms, tmp, gVectors, gAbelian}
	,
	(* expand Lagrangian and write mass terms in canonical format *)
	lag = CanonizeFermionMassTerms@ ContractCGs@ lag;

	If[OptionValue@ Mode === Matching,
		lag = IntroduceEffectiveMasses@ lag;
	];

	lag = RelabelIndices@ Contract@ lag;
	
	(* Add gauge-fixing and ghost terms for all gauge vectors *)
	gVectors= Select[DeleteDuplicates@ Cases[lag, (Field|FieldStrength)[f_, __]:> f, Infinity],GaugeFieldQ];
	gAbelian= Select[gVectors, $GaugeGroups[First@GaugeGroupByProperty[Field->#], Abelian]&];

	(* Assumes ghosts and vectors are quantum (removed in X terms) *)
	(* WARNING: The GF terms need to be written in this very specific way or it will fail to give right X-terms*)
	gfTerms=Sum[
			If[MemberQ[gAbelian,gVec],
					1/2* DetermineCouplingOfGaugeField[gVec, lag][Coupling]* Module[{\[Alpha],\[Beta]}, Bar@gVec[\[Alpha]]BackgroundCD[{\[Beta],\[Alpha]}, gVec[\[Beta]]] ]
			,
					1/2* DetermineCouplingOfGaugeField[gVec, lag][Coupling]* Module[{\[Alpha],\[Beta],a,b,c}, Bar@gVec[\[Alpha],a]BackgroundCD[{\[Beta],\[Alpha]}, gVec[\[Beta],a]] - I FieldGenerators[gVec, GroupFromRep@@$FieldAssociation[gVec][Indices], {a,b,c}]gVec[\[Alpha],a]gVec[\[Beta],b] BackgroundFS[gVec,\[Alpha],\[Beta],c] ]
				+
				Module[{gh, ghbar, indK, \[Mu], a, b, c},
					{gh, ghbar}= $GhostFields@ gVec;
					FuncNCM[CD[\[Mu],ghbar[indK]],CD[\[Mu],gh[indK]]] + I FieldGenerators[gVec, GroupFromRep@@$FieldAssociation[gVec][Indices], {a,b,c}]gVec[\[Mu],b]FuncNCM[ghbar[a],CD[\[Mu],gh[c]]]
				]
			]
			,
			{gVec,gVectors}
		];
	lagFixed= lag + gfTerms;
				
	If[(newLag = lag =!= $currentLagrangian),
		If[!CheckLagrangian@ lag, (*Abort[];*)Null ];
		(* set the new Lagrangian *)
		$currentLagrangian= lag;
		$currentMode= OptionValue@ Mode;
		(* determine its fields with their Masses *)
		DetermineCurrentFieldAssociation[lagFixed, Mode-> OptionValue@ Mode];
	];

	eftOrd= If[Head@ eftOrder === List, First@ eftOrder, eftOrder];
	If[newLag || eftOrder > $currentEFTOrder,
		(* Reset $currentHeavyDims *)
		$currentHeavyDims = <||>;
		(* Derive and solve the EOMs *)
		DetermineEOMs[lag, EFTOrder-> eftOrder]; (* This internally sets $currentEOMs and $currentHeavyDims *)
	];
	$currentEFTOrder= eftOrder;
	
	If[(newLag || loopOrder > $currentLoopOrder) && ($currentLoopOrder= loopOrder) > 0,
		OptionalMonitor[OptionValue@ Verbose,
			SetSubstitutions[
				If[loopOrder<2, lagFixed/.hbar->0, lagFixed], (* drop loop-suppressed contributions to Xterms *)
				EFTOrder-> eftOrder, 
				Mode-> OptionValue@ Mode
			];
		, "Evaluating X-terms"];
	];

	lag
];


(* ::Subsection:: *)
(*Matching routines*)


(* ::Subsubsection::Closed:: *)
(*List the power-type supertraces*)


(* ::Text:: *)
(*To determine all the kinds of power type tracers labelled by their propagator types*)
(*Future improvements possible to evanescent and counterterm traces is possible if including momentum order.*)


ListPowerTypeTraces[{order_Integer}, opt___]:= ListPowerTypeTraces[order, opt];


ListPowerTypeTraces[order_Integer, lightOnly_:False]:= Module[{possibilities, seed, temp, Xords},
	Xords= $XOrdMin (*$currentXdims*);
	(*Recursively determines candidates for traces with order \[LessEqual] order*)
	possibilities= Flatten[Last@ Reap[
		(*Determines whether to only consider light traces or traces with at least one heavy field*)
		If[lightOnly,
			seed= Sow@ {{lScalar}, {lFermion}, {lVector}, {lGhost}};
			While[(Length@ seed> 0),
				temp= Flatten/@ Tuples[{seed, {lScalar, lFermion, lVector}}];
				seed= Sow@ Select[temp, Total@ BlockMap[Xords, #, 2, 1] < order&];
			];
		,
			seed= Sow@ {{hScalar}, {hFermion}, {hVector}, {hGhost}, {hAntiGhost}};
			While[(Length@ seed> 0),
				temp= Flatten/@ Tuples[{seed, {hScalar, lScalar, hFermion, lFermion, hVector, lVector, hGhost, lGhost, hAntiGhost, lAntiGhost}}];
				seed= Sow@ Select[temp, Total@ BlockMap[Xords, #, 2, 1] < order && Length@ # <= order &];
			];
		];
	], 2];

	(*Eliminate based on no. of X-terms*)
	possibilities= Select[possibilities, Length@ # <= order &];

	(*Eliminates traces with order > order*)
	possibilities= DeleteCases[possibilities,
		_? (Total[Min/@ BlockMap[Xords, #, 2, 1]]+ Min@ Xords@ #[[{-1, 1}]] > order&), {1}];

	(*Delete duplicates under cyclic permutations*)
	DeleteDuplicatesBy[possibilities, (First@ Sort@ NestList[RotateLeft, #, Length@# - 1] &)]
];


(* ::Subsubsection::Closed:: *)
(*All one-loop contributions*)


(* ::Text:: *)
(*The full 1-loop matching function evaluating all log- and power-type traces relevant to the Lagrangian.*)


Options[LoopMatch]= {
	EFTOrder -> 6,
	Verbose -> True,
	WhichTraces -> All
	};


LoopMatch[opt:OptionsPattern[]]? OptionsCheck:= Module[
		{field, fields, out, logTraces, powerTraces, i=0, myTraces},

	(* determine log and power traces *)
	logTraces = {hScalar, hFermion, hVector, hGhost, hAntiGhost};
	powerTraces = ListPowerTypeTraces[OptionValue@ EFTOrder];

	(* reduce to manually selected traces if applicable *)
	If[OptionValue[WhichTraces] =!= All,
		myTraces = OptionValue[WhichTraces]/.fieldFormatReverse;
		logTraces = Intersection[logTraces, {#}]&/@myTraces;
		logTraces = logTraces /. {{{a___}}:>{a}} /. {{}->Nothing};
		powerTraces = Intersection[powerTraces, {#}]&/@myTraces;
		powerTraces = powerTraces /. {{{a___}}:>{a}} /. {{}->Nothing};
		Print["Computing the Log-Type traces: ", logTraces/.fieldFormat];
		Print["Computing the Power-Type traces: ", powerTraces/.fieldFormat];
	];
	
	out= OptionalMonitor[OptionValue@ Verbose,

			Sum[
				(*Check if a field with non zero charges exists*)
				If[Or[
					Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field, Charges-> {__}] > 0,
					Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field,
						Indices-> inds_/; IntersectingQ[GroupFromRep/@ inds, Keys@ $GaugeGroups]] > 0
				],
					LogTypeSTr[field, OptionValue@ EFTOrder]
				,
					0
				]
			, {field, logTraces}]
		, StringForm["Evaluating log-type supertrace: `1`", field/. fieldFormat] ];

	out+= OptionalMonitor[OptionValue@ Verbose,
			Sum[i++;
				Sow@ StringReplace[ToString@fields,{", " -> "-", "{" -> "", "}" -> "", "Matchete`PackageScope`" -> ""}];
				Sow@ PowerTypeSTr[fields, OptionValue@ EFTOrder, FilterOptions[opt, PowerTypeSTr]]
			, {fields, powerTraces}]
		, StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)",
			fields/. fieldFormat, i, Length@ powerTraces] ];
	
	out
];


fieldFormat= {hScalar-> "\[CapitalPhi]", lScalar-> "\[Phi]", hFermion-> "\[CapitalPsi]", lFermion-> "\[Psi]", hVector-> "V", lVector-> "A", lGhost-> "cA", hGhost-> "cV", lAntiGhost-> "cbarA", hAntiGhost-> "cbarV"};
fieldFormatReverse= {"\[CapitalPhi]"-> hScalar, "\[Phi]"-> lScalar, "\[CapitalPsi]"-> hFermion, "\[Psi]"-> lFermion, "V"-> hVector, "A"-> lVector, "cA"-> lGhost, "cV"-> hGhost,"cbarA"-> lAntiGhost, "cbarV"-> hAntiGhost};


(* ::Subsubsection::Closed:: *)
(*Routine for specific covariant loops*)


(* ::Text:: *)
(*The user function for computing specific supertrace contributions by specifying the internal fields. *)


CovariantLoop::bkgfld= "The field(s) `1` is a non-propagating background field."
CovariantLoop::ukwnfld= "The field(s) `1` is not part of the Lagrangian."
CovariantLoop::noheavy= "There are no heavy fields in the loop."


Options@ CovariantLoop= {
		EFTOrder-> 6,
		Mode-> Matching
	};


OptionTest[CovariantLoop, Mode]= MatchQ[Matching| Divergence];
OptionMessage[Mode, CovariantLoop, val_]:= Message[General::optexpectsval, Mode, CovariantLoop, val, "value Matching or Divergence"];


CovariantLoop[lag_, field_Symbol, opts:OptionsPattern[]]:= CovariantLoop[lag, {field}, opts];


CovariantLoop[lagrangian_, fields_List, opts:OptionsPattern[]]? OptionsCheck:=
CovariantLoop[lagrangian, fields, opts]= Module[
		{bkgFields, lagFields, loopMode, n, types, dofNumbers, ord, out, lag= lagrangian},
	loopMode= OptionValue@ Mode;
	lagFields= DeleteDuplicates@ Join[
			Cases[lag, Bar[(Field|FieldStrength)[label_, __]]:> Conj@ label, All],
			Cases[lag, (Field|FieldStrength)[label_, __]:> label, All]
		];

	(*Check fields*)
	If[!SubsetQ[lagFields, fields],
		Message[CovariantLoop::ukwnfld, Complement[fields, Join@@ List@@ lagFields]];
		Abort[];
	];
	bkgFields= Intersection[fields, FieldByProperty[BackgroundField-> True]];
	If[Length@ bkgFields> 0,
		Message[CovariantLoop::bkgfld, bkgFields];
		Abort[];
	];

	(*Update Lagrangian*)
	ord= OptionValue@ EFTOrder;
	lag = SetCurrentLagrangian[lag, 1, If[Head@ ord === List, First @ord, ord], Mode-> loopMode];
	(* check that there is at least one heavy field, must be done after setting the Lagrangian *)

	If[loopMode === Matching && Intersection[GetFieldsUpdatedByProperty[Heavy-> True], fields] === {},
		Message[CovariantLoop::noheavy];
		Abort[];
	];

	(*Number the fields *)
	types= FieldType/@ fields;

	(*Add log type STr*)
	out= If[Length@ fields === 1,
			LogTypeSTr[First@ types, ord, Fields-> fields, Mode-> loopMode]
		,
			0
		];

	out+ PowerTypeSTr[types, ord, Fields-> fields, Mode-> loopMode]//ContractCGs//MatchReduce
];


FieldType[f_]:= Switch[Lookup[$FieldAssociation[f], {Type, Heavy}]
	,{Scalar, True}, hScalar
	,{Scalar, False}, lScalar
	,{Fermion, True}, hFermion
	,{Fermion, False}, lFermion
	,{Vector, True}, hVector
	,{Vector, False}, lVector
	,{Ghost, True}, hGhost
	,{Ghost, False}, lGhost
	,{AntiGhost, True}, hAntiGhost
	,{AntiGhost, False}, lAntiGhost
];


(* ::Subsubsection::Closed:: *)
(*Main Matchete routine for integrating out fields*)


(* If True, the Match routine writes its results to an internal database *)
$ValidationRun= False;


Options[Match]={EFTOrder -> 6, LoopOrder-> 1, Simplifications -> All, Verbose -> Monitor, WhichTraces -> All};


Match::heavyvec= "Lagrangian contains massive vector(s) `1`. Massive vectors are not yet supported for loop-level matching. Stay tuned for future updates\[Ellipsis]";


Match[lag_, opts:OptionsPattern[]]? OptionsCheck :=
Match[lag, opts] = Module[{
		(*lagrangian = BetterExpand[ContractCGs@lag],*)
		lagrangian = lag,
		eftOrder=OptionValue@EFTOrder,
		loopOrder=OptionValue@LoopOrder,
		LagrangianEFT,
		LagrangianEFT1,
		heavyVectors,
		ReplaceHeavyEOMOpts = FilterOptions[opts, ReplaceHeavyEOM],
		VerboseOption = (OptionValue@Verbose===Print||OptionValue@Verbose===Monitor),
		traceResults,
		time
	},

	(*Check that no heavy vectors appear in loop-level matching*)
	If[MatchQ[loopOrder, 1|{1}],
		heavyVectors= MassiveVectorsInLag@ lag;
		If[Length@ heavyVectors =!= 0,
			Message[Match::heavyvec, heavyVectors];
			Abort[];
		];
	];

	time= First@ Timing[

	(* Set global variables for the given Lagrangian *)
	lagrangian = SetCurrentLagrangian[lagrangian, If[loopOrder === {1}, 1, loopOrder],
		If[Head@ eftOrder === List, First @eftOrder, eftOrder], Verbose-> VerboseOption];

	MyPrint["Integrating out the fields: ", Sequence@@Riffle[Intersection[Matchete`PackageScope`OccurringFields[lagrangian],GetFieldsUpdatedByProperty[Heavy->True]],", "], Verbose->OptionValue@Verbose===Print];

	LagrangianEFT= If[MatchQ[loopOrder, 0|1],
			(*Tree-level Lagrangian*)
			OptionalMonitor[VerboseOption,
				RelabelIndices@ RefineDiracProducts@ ReplaceHeavyEOM[lagrangian, ReplaceHeavyEOMOpts]
			, "Matching at tree level\[Ellipsis]"]
		, 0]+ If[MatchQ[loopOrder, 1|{1}],
			(*One-loop Lagrangian*)
			OptionalMonitor[VerboseOption,
				traceResults= Reap[
					LagrangianEFT1 = LoopMatch[EFTOrder->eftOrder,Verbose->VerboseOption, WhichTraces->OptionValue[WhichTraces]];
				];
				traceResults=Quiet[traceResults[[2,1]]]; (* this produces an error when the WhichTrace option is used, which however should never be used in validation *)
				LagrangianEFT1
			, "Matching at 1-loop level..."]
		, 0];

	LagrangianEFT= LagrangianEFT//ContractCGs//MatchReduce;
	];
	
	(* This is only performed if the result of the matching should be saved for validation *)
	If[$ValidationRun,
		QuietEcho@ SaveValidationResults[traceResults, LagrangianEFT, time, lag, "Directory"->"MatchingResults"];
	];
	
	(* ensure that the output Lagrangian ia manifestly Hermitian *)
	(* this might not be the case due to mass/index permutations in LF that are not recognized. *)
	LagrangianEFT = (LagrangianEFT + Bar[LagrangianEFT])/2;
	
	LagrangianEFT
];


(* ::Text:: *)
(*Checks if what massive vectors are present in the Lagrangian. Ignores Light vectors unless Light-> True.*)


Options@ MassiveVectorsInLag= {Light-> False};


MassiveVectorsInLag[lag_, OptionsPattern[]]:= Module[{massiveVectors, lagVectors},
	massiveVectors= If[OptionValue@ Light,
			FieldByProperty[Type-> Vector, Mass-> Except[0]]
		,
			FieldByProperty[Type-> Vector, Heavy-> True]
		];
	lagVectors= DeleteDuplicates@ Cases[lag, (Field[lab_, _Vector, __] | FS[lab_, __])-> lab, All];
	Intersection[massiveVectors, lagVectors]
];


(* ::Subsubsection::Closed:: *)
(*Routine to simplify the output*)


(* ::Text:: *)
(*Function to convert CGs to Matchete Deltas and to evaluate single scale loop integrals*)


MatchReduce[expr_]:= Module[{res,summedInd},
	res= Contract[ContractCGs[expr/.
	{
		(*\[Epsilon]^-1-> 0,*)(*Needs to be done after Gamma reduction*)
		FieldStrength[label_,linds_,{Bar@ind1_,ind2_},{CDer___}]:> Module[{A},
				-CG[Bar[gen[ind2[[2]]]],{Index[A,GroupFromRep[ind2[[2]]][adj]],Bar@ind1,ind2}] FieldStrength[label,linds,{Index[A,GroupFromRep[ind2[[2]]][adj]]},{CDer}]
			],
		FieldStrength[label_,linds_,{ind1_,Bar@ind2_},{CDer___}]:> Module[{A},
				CG[gen[ind2[[2]]],{Index[A,GroupFromRep[ind2[[2]]][adj]],ind1,Bar@ind2}] FieldStrength[label,linds,{Index[A,GroupFromRep[ind2[[2]]][adj]]},{CDer}]
			],
		FieldStrength[label_,linds_,{ind1_Index,ind2_Index},{CDer___}]:> Module[{A},
				CG[gen[ind2[[2]]],{Index[A,GroupFromRep[ind2[[2]]][adj]],ind1,ind2}] FieldStrength[label,linds,{Index[A,GroupFromRep[ind2[[2]]][adj]]},{CDer}]
			]
	}]//.
	{
		LF[masses_?(Length[#]==1&),powers_]:> EvaluateLoopFunctions@LF[masses,powers],
		HoldPattern@CG[del[_],{a_,b_}]:> Delta[a,b],
		HoldPattern[CG[eps[x_],a:{_Index..}]CG[Bar@eps[x_],b:{_Bar..}]]:> Det[Outer[Delta,a,b]]
	}];

	(* expand results to reduce auxiliary FlavorSum expressinos *)
	res= TermsToList@ res;
	res= Sum[
		If[FreeQ[term,_FlavorSum,All],
			term
			,
			summedInd= FirstCase[term,FlavorSum[ind_]:>ind,Nothing,All];
			If[FreeQ[term/._FlavorSum->1,summedInd,All],
				term/._FlavorSum->DimRep[Last@summedInd],
				term
			]
		]
		,
		{term,res}
	];
	res

]


(* ::Subsubsection::Closed:: *)
(*Canonize fermion mass terms*)


CanonizeFermionMassTerms[lagrangian_] := Module[
	{
		lag = TermsToList@ lagrangian (*needs to be already expanded*)
	}
	,
	lag = Plus@@Table[
		Which[
		(* canonize Dirac masses of vectorlike fermions *)
		MatchQ[term, aux_ * Bar@Field[l_,Fermion,___]\[CenterDot] DiracProduct[Proj[-1]]\[CenterDot] Field[l_,Fermion,___] /; (FreeQ[aux, _Field, All] && !($FieldAssociation[l][SelfConjugate]) && $FieldAssociation[l][Chiral]===False)]
		,
			(term /. {Proj[-1]->1}) - (term /. {Proj[-1]->Proj[+1]})
		,
		(* canonize Majorana masses *)
		MatchQ[term, aux_ * Transp[Field[l_,Fermion,___]]\[CenterDot] DiracProduct[GammaCC,Proj[-1]]\[CenterDot] Field[l_,Fermion,___] /; (FreeQ[aux, _Field, All] && $FieldAssociation[l][SelfConjugate] && $FieldAssociation[l][Chiral]===False)],
			(term /. {Proj[-1]->1}) - (term /. {Proj[-1]->Proj[+1]})
		,
		True
		,
			term
		]
		,
		{term,lag}
	];

	RelabelIndices[lag]
];
