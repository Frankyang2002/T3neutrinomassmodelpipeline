(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Counterterms`*)


(* ::Subtitle:: *)
(*Tools needed to calculate finite and divergent counterterms*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsection:: *)
(*Exported*)


PackageExport["OnShell"]


PackageExport["DR2MS"]


(* ::Subsection:: *)
(*Internal*)


PackageScope["UVDivergentAction"]


PackageScope["EvanescentContribution"]


PackageScope["TruncateLagrangian"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


OnShell::usage = "Option for UVDivergentAction. Can be set to: Automatic, True, or False";


(* ::Subsection:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Divergence*)


(* ::Subsection:: *)
(*UVDivergentAction*)


Options@ UVDivergentAction = { EFTOrder -> Automatic , Verbose -> True, OnShell -> Automatic };


UVDivergentAction::orderArgument = "The option EFTOrder-> `1`  is neither an integer nor Automatic and cannot be used as an order specification.";


UVDivergentAction[expr_,OptionsPattern[]]:=Module[{powerTraces,logTraces,isOnShell=True,L=expr,order,
								ptr,LDiv,field,i},

	(* determine order argument and cut the Lagrangian to the given order *)
	{order, L}= TruncateLagrangian[OptionValue@ EFTOrder, L];
	
	(* check if this Lagrangian contains redundant operators - if so, this Lagrangian is not in the on-shell basis *)
	If[Length @ FieldsToShift @ InternalSimplify @ L > 0, isOnShell = False];
	If[OptionValue@OnShell===False, isOnShell = False];
	
	(* set the expanded Lagrangian and determine X terms *)
	SetCurrentLagrangian[L, 1, order, Mode -> Divergence];
	 
	$MonitorString="";
	
	(* determine the Log traces *)
	logTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
					Sum[
						If[Or[
							Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field, Charges-> {__}] > 0,
							Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field,
								Indices-> inds_/; IntersectingQ[GroupFromRep/@ inds, Keys@ $GaugeGroups]] > 0
						],
							$MonitorString = StringForm["Evaluating log-type supertrace: `1`", field/. fieldFormat];
							LogTypeSTr[field, order, Mode -> Divergence]
						,
							0
						]
					,{field,{hScalar,hFermion,hVector,hGhost,hAntiGhost,lScalar,lFermion,lVector,lGhost,lAntiGhost}}
					],
					$MonitorString					
				];
	
	(* list the possible power traces, including light fields *)
	ptr = ListPowerTypeTraces[ order, True ];
	
	powerTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
						Sum[
							$MonitorString = StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)",
											ptr[[i]] /. fieldFormat, i, Length@ ptr];
							PowerTypeSTr[ptr[[i]],order,Mode->Divergence],
						{i, Length@ptr}]
					,
						$MonitorString];
	
	(* ExtractePoles can be significantly faster than Series *)
	LDiv= SeparateOutConstants[logTraces + powerTraces][[2]];
	LDiv = ExtractPoles[LDiv];

	If[isOnShell  ,
		ExtractPoles[EOMSimplify[L + LDiv,DummyCoefficients->True]], 
		GreensSimplify@ LDiv
	]
]


TruncateLagrangian[eftOrderOption_, lag_]:= Module[{order}, 
	Switch[eftOrderOption
	,_Integer,
		{eftOrderOption, SeriesEFT[lag, EFTOrder-> eftOrderOption]}
	,Automatic,
		{order= GetMaxEFTOrder@ lag, SeriesEFT[lag, EFTOrder-> order]}
	,_,
		Message[UVDivergentAction::orderArgument, eftOrderOption]; Abort[];
	]
];


(* function to extractt 1/\[Epsilon] poles in an efficient way *)
ExtractPoles[expr_]:=Module[{x},
	Coefficient[expr/.{\[Epsilon]->\[Epsilon]/x}, x]
]


(* ::Section:: *)
(*Evanescent contribution*)


(* ::Subsection:: *)
(*EvanescentContribution*)


Options@ EvanescentContribution = {EFTOrder-> Automatic, Verbose-> True};


EvanescentContribution::orderArgument = "The option EFTOrder-> `1`  is neither an integer nor Automatic and cannot be used as an order specification.";


EvanescentContribution[lagrangian_, OptionsPattern[]]? OptionsCheck:= Module[
		{lag, order, out, traceTypes, fields, i= 0},
		
	(* Determine order argument and cut the Lagrangian to the given order *)
	{order, lag}= TruncateLagrangian[OptionValue@ EFTOrder, lagrangian];
	
	(* Set the expanded Lagrangian and determine X terms *)
	SetCurrentLagrangian[lag, 1, order, Mode-> Evanescent, Verbose-> OptionValue@ Verbose];
	
	(* Loop-level evanescent insertions *)
	traceTypes= ListEvanescentTypeTraces[order];
	out= OptionalMonitor[OptionValue@ Verbose,
			Sum[i++;
				PowerTypeSTr[fields, order, Mode-> Evanescent]
			, {fields, traceTypes}]
		, StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)",
			fields/. fieldFormat, i, Length@ traceTypes] ];
	
	out// ContractCGs// MatchReduce
]


(* ::Subsubsection::Closed:: *)
(*Listing evanescent trace types*)


(* ::Text:: *)
(*ListEvanescentTypeTraces provides the list of all power type traces beginning the reading from an evanescent operator, while being indifferent as to the nature of the remaining insertions*)


ListEvanescentTypeTraces@ {order_Integer}:= ListEvanescentTypeTraces@ order;


ListEvanescentTypeTraces@ order_Integer:= Module[{possibilities, seed, temp},
	(*Recursively determines candidates for traces with order \[LessEqual] order*)
	possibilities= Flatten[Last@ Reap[
		seed= Sow@ {{lScalar}, {lFermion}, {lVector}};
		While[(Length@ seed> 0),
			temp= Flatten/@ Tuples[{seed, {lScalar, lFermion, lVector}}];
			seed= Sow@ Select[temp, $EvaOrdMin@ #[[;;2]]+ Total[BlockMap[$XOrdMin, #[[2;;]], 2, 1]] <= order&];
		];
	], 2];

	(*Eliminates traces with order > order*)
	possibilities= Select[possibilities,
		If[Length@# === 1, $EvaOrdMin@ #[[{-1, 1}]], 
			$EvaOrdMin@ #[[;;2]]+ Total@ BlockMap[$XOrdMin, #[[2;;]], 2, 1]+ $XOrdMin@ #[[{-1, 1}]]] <= order&]
];


(* ::Section:: *)
(*Finite scheme-change counterterms*)


(* ::Subsection:: *)
(*DRbar-MSbar*)


(* ::Subsubsection::Closed:: *)
(*Extracting the Lagrangian terms relevant for tree-level matching*)


GetTreeLevelTerms::usage=" returns all terms of the Lagrangian that can potentially contribute to the matching at tree level.";

GetTreeLevelTerms::sanitycheck= "Could not correctly identify tree-level terms."

GetTreeLevelTerms[lag_]:=Module[{L,Ltree1,Ltree2,Ltree,Lloop},
	(* only keep tree-level contributions *)
	L = lag/. hbar->0;
	
	(* keep only terms with at most one heavy field *)
	L = Expand[L];
	
	(* get all Lagrangian terms with at most one heavy field *)
	Ltree1 = Plus@@Cases[
		L,
		x_/;Length[Cases[RemovePower[x],f_Field/;$FieldAssociation[First@f,Heavy],All]]<=1,
		1
	];
	
	(* get kinetic and mass terms of the heavy fields that contribute to the tree-level matching *)
	Ltree2 = Plus@@Cases[
		L,
		x_/;(
		(Length[Cases[RemovePower[x],_Field,All]]==2)
		&&
		MatchQ[Cases[RemovePower[x],f_Field/;$FieldAssociation[First@f,Heavy]:>First[f],All],{y1_,y2_}/;((!FreeQ[Ltree1,y1,All])&&(!FreeQ[Ltree1,y2,All]))]),
		1
	];
	
	(* get all Lagrangian terms relevant for tree-level matching *)
	Ltree = Ltree1 + Ltree2;
	
	(* loop part *)
	Lloop = L-Ltree;
	
	(* sanity check *)
	If[((Ltree+Lloop)-L)=!=0, 
		Message[GetTreeLevelTerms::sanitycheck]; 
		Abort[]
	];
	
	Return@{Ltree,Lloop}
]


(* ::Subsubsection::Closed:: *)
(*Helper functions*)


(* ::Text:: *)
(*Determine gauge coupling from representations*)


GaugeCouplingsFromRep::usage= "returns the gauge coupling corresponding to the gauge group of a given representation.";


GaugeCouplingsFromRep[rep_]:= $GaugeGroups[Head[rep],Coupling][]


(* ::Text:: *)
(*Determine all representations of a field*)


GetFieldGaugeReps::usage= "returns all representations of non-Abelian groups under which a given field transforms together with all its charges under Abelian groups.";


GetFieldGaugeReps[\[Psi]_]:= $FieldAssociation[\[Psi],Indices]~Join~$FieldAssociation[\[Psi],Charges]/.Alternatives@@(Keys@$GlobalGroups~Join~Keys@$FlavorIndices)->Nothing


(* ::Text:: *)
(*Calculate the quadratic Casimir operator for a representation*)


QuadraticCasimir::usage= "returns the quadratic Casimir operator for a given representation.";
QuadraticCasimir::noCasimir= "The quadratic Casimir could not be de termined for `1`.";


QuadraticCasimir[rep_]:=If[MemberQ[Keys@GetGroups[],Head[rep]],
	(* non-Abelian groups *)
	Casimir2[Head[rep]/.$Groups,$Representations[rep,DynkinCoefficients]]
	,
	(* Abelian groups *)
	If[MatchQ[rep,(_[n_]/;NumberQ[n])],
		First[rep]^2
		,
		Message[QuadraticCasimir::noCasimir,rep];
		Abort[]
	]
]


(* ::Text:: *)
(*Calculate the squared gauge coupling times the quadratic Casimir operator summed over all representations of a field*)


SumC2g2::usage= "for a given field label this functions returns the quadratic Casimir operator times the square of the associated gauge coupling summed over all reresentations of the field.";


SumC2g2[label_]:=Module[{tmp=GetFieldGaugeReps[label],C2,g},
	C2 = QuadraticCasimir/@ tmp;
	g  = GaugeCouplingsFromRep/@ tmp;
	Plus@@(C2*(g)^2)
]


(* ::Text:: *)
(*Introduce flavor deltas such that operators do not contain internally contracted flavor indices*)


UnifyOpIndices::usage= "for every dummy index contracted within a given operator a Kronecker delta is introduced.";


(* unifies the dummy indices in an operator by introducing new indices and factoring out appropriate Delta functions *)
UnifyOpIndices[oper_Operator]:=Module[{op=oper,inds,newInds,deltas,rules},
	inds    = FindDummyIndices[op];
	newInds = Table[ind/.First[ind]->Unique[],{ind,inds}];
	deltas  = Times@@(Delta[First@#,Last@#]&/@Transpose@{inds,newInds});
	rules   = (First@#->Last@#)&/@Transpose@{inds,newInds};
	
	(* ensure every index replacement is only performed once *)
	Do[
		op= ReplaceFirst[op,rule]
		,
		{rule,rules}
	];
	op*deltas
]


(* ::Subsubsection::Closed:: *)
(*Scheme change shift for quartic scalar couplings*)


QuarticScalarShift::usage= "takes four scalar fields as input in the form [Bar@\!\(\*SubscriptBox[\(\[Phi]\), \(1\)]\),Bar@\!\(\*SubscriptBox[\(\[Phi]\), \(2\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(3\)]\),\!\(\*SubscriptBox[\(\[Phi]\), \(4\)]\)] and returns the corresponding finite counterterm originating from the \!\(\*OverscriptBox[\(DR\), \(_\)]\)-\!\(\*OverscriptBox[\(MS\), \(_\)]\) scheme change.";


QuarticScalarShift[\[Phi]1_,\[Phi]2_,\[Phi]3_,\[Phi]4_]:= Module[
	{
		groups,
		A,B,a1,b1,c1,d1,e1,f1,a2,b2,c2,d2,e2,f2,
		ind1,ind2,ind3,ind4,id13$24,id14$23,
		res,
		includeDeltas
	}
	,
	(* includes a Delta[a_i,b_i] for indA={a_1,a_2,...} and indB={b_1,b_2,...} if indA and indB have same length, otherwise 0 is returned *)
	includeDeltas[indA_,indB_]:= If[Length[indA]==Length[indB],Times@@(Delta[First[#],Last[#]]&/@Transpose@{indA,indB}),0];
	
	(* markers for when \[Phi]i=\[Phi]j *)
	id13$24=If[(First@\[Phi]1===First@\[Phi]3)&&(First@\[Phi]2===First@\[Phi]4),1,0];
	id14$23=If[(First@\[Phi]1===First@\[Phi]4)&&(First@\[Phi]2===First@\[Phi]3),1,0];
	
	(* list of gauge groups commen between all fields *)
	groups=Intersection[Head/@GetFieldGaugeReps[First@\[Phi]1],
		Head/@GetFieldGaugeReps[First@\[Phi]2],
		Head/@GetFieldGaugeReps[First@\[Phi]3],
		Head/@GetFieldGaugeReps[First@\[Phi]4]
	];
	
	res=-(hbar/4)*Sum[
		(* field indices that belong to the current gauge rep. *)
		a1=FirstCase[\[Phi]1,Index[l_,gr1[_]]:>l,0,All];
		c1=FirstCase[\[Phi]3,Index[l_,gr1[_]]:>l,0,All];
		d1=FirstCase[\[Phi]2,Index[l_,gr1[_]]:>l,0,All];
		f1=FirstCase[\[Phi]4,Index[l_,gr1[_]]:>l,0,All];
		a2=FirstCase[\[Phi]1,Index[l_,gr2[_]]:>l,0,All];
		c2=FirstCase[\[Phi]3,Index[l_,gr2[_]]:>l,0,All];
		d2=FirstCase[\[Phi]2,Index[l_,gr2[_]]:>l,0,All];
		f2=FirstCase[\[Phi]4,Index[l_,gr2[_]]:>l,0,All];
		(* field indices that do not belong to the current gauge reps and that should be contracted by Deltas *)
		ind1=Cases[\[Phi]1,Index[Except[a1|a2],_],All];
		ind2=Cases[\[Phi]2,Index[Except[d1|d2],_],All];
		ind3=Cases[\[Phi]3,Index[Except[c1|c2],_],All];
		ind4=Cases[\[Phi]4,Index[Except[f1|f2],_],All];
		(* include proper power of gauge couplings *)
		$GaugeGroups[gr1,Coupling][]^2*$GaugeGroups[gr2,Coupling][]^2*
		(* distinguish Abelian and non-Abelian groups *)
		Switch[{$GaugeGroups[gr1,Abelian],$GaugeGroups[gr2,Abelian]},
			(* {Abelian, Abelian} *)
			{True,True},
				(id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]+id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3])(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1]
				),
			(* {Abelian, non-Abelian} *)
			{True,False},
				id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2,{B,a2,c2}]+
					FieldGenerators[First@\[Phi]1,gr2,{B,a2,c2}]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2,{B,d2,f2}]+
					FieldGenerators[First@\[Phi]2,gr2,{B,d2,f2}]FieldGenerators[First@\[Phi]4,gr1]
				)+
				id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
					FieldGenerators[First@\[Phi]1,gr1]FieldGenerators[First@\[Phi]3,gr2,{B,d2,c2}]+
					FieldGenerators[First@\[Phi]1,gr2,{B,d2,c2}]FieldGenerators[First@\[Phi]3,gr1]
				)*(
					FieldGenerators[First@\[Phi]2,gr1]FieldGenerators[First@\[Phi]4,gr2,{B,a2,f2}]+
					FieldGenerators[First@\[Phi]2,gr2,{B,a2,f2}]FieldGenerators[First@\[Phi]4,gr1]
				),
			(* {non-Abelian, Abelian} *)
			{False,True},
				id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
					FieldGenerators[First@\[Phi]1,gr1,{A,a1,c1}]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1,{A,a1,c1}]
				)*(
					FieldGenerators[First@\[Phi]2,gr1,{A,d1,f1}]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1,{A,d1,f1}]
				)+
				id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
					FieldGenerators[First@\[Phi]1,gr1,{A,d1,c1}]FieldGenerators[First@\[Phi]3,gr2]+
					FieldGenerators[First@\[Phi]1,gr2]FieldGenerators[First@\[Phi]3,gr1,{A,d1,c1}]
				)*(
					FieldGenerators[First@\[Phi]2,gr1,{A,a1,f1}]FieldGenerators[First@\[Phi]4,gr2]+
					FieldGenerators[First@\[Phi]2,gr2]FieldGenerators[First@\[Phi]4,gr1,{A,a1,f1}]
				),
			(* {non-Abelian, non-Abelian} *)
			{False,False},
				If[gr1===gr2,
					(* if both groups are identical *)
					id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
						FieldGenerators[First@\[Phi]1,gr1,{A,a1,b1}]FieldGenerators[First@\[Phi]3,gr1,{B,b1,c1}]+
						FieldGenerators[First@\[Phi]1,gr1,{B,a1,b1}]FieldGenerators[First@\[Phi]3,gr1,{A,b1,c1}]
					)*(
						FieldGenerators[First@\[Phi]2,gr1,{A,d1,e1}]FieldGenerators[First@\[Phi]4,gr1,{B,e1,f1}]+
						FieldGenerators[First@\[Phi]2,gr1,{B,d1,e1}]FieldGenerators[First@\[Phi]4,gr1,{A,e1,f1}]
					)+
					id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
						FieldGenerators[First@\[Phi]1,gr1,{A,d1,b1}]FieldGenerators[First@\[Phi]3,gr1,{B,b1,c1}]+
						FieldGenerators[First@\[Phi]1,gr1,{B,d1,b1}]FieldGenerators[First@\[Phi]3,gr1,{A,b1,c1}]
					)*(
						FieldGenerators[First@\[Phi]2,gr1,{A,a1,e1}]FieldGenerators[First@\[Phi]4,gr1,{B,e1,f1}]+
						FieldGenerators[First@\[Phi]2,gr1,{B,a1,e1}]FieldGenerators[First@\[Phi]4,gr1,{A,e1,f1}]
					)
					,
					(* if both groups are different *)
					id13$24*includeDeltas[ind1,ind3]*includeDeltas[ind2,ind4]*(
						2*FieldGenerators[First@\[Phi]1,gr1,{A,a1,c1}]FieldGenerators[First@\[Phi]3,gr2,{B,a2,c2}]
					)*(
						2*FieldGenerators[First@\[Phi]2,gr1,{A,d1,f1}]FieldGenerators[First@\[Phi]4,gr2,{B,d2,f2}]
					)+
					id14$23*includeDeltas[ind1,ind4]*includeDeltas[ind2,ind3]*(
						2*FieldGenerators[First@\[Phi]1,gr1,{A,d1,c1}]FieldGenerators[First@\[Phi]3,gr2,{B,d2,c2}]
					)*(
						2*FieldGenerators[First@\[Phi]2,gr1,{A,a1,f1}]FieldGenerators[First@\[Phi]4,gr2,{B,a2,f2}]
					)
				]
		]
		,
		{gr1,groups},{gr2,groups}
	]//ContractCGs
]


(* ::Subsubsection::Closed:: *)
(*Scheme conversion function: DR-MS*)


DR2MS::usage= "DR2MS[lag] takes the renormalizable Lagrangian \!\(\*
StyleBox[\"lag\", \"TI\"]\) that is renormalized in dimensional reduction using the \!\(\*OverscriptBox[\(DR\), \(_\)]\)  scheme as input and returns the version of that Lagrangian renormalized using dimensional regularization in the \!\(\*OverscriptBox[\(MS\), \(_\)]\) scheme by incorporating the appropriate finite counterterms required for the scheme change.";


DR2MS::nonren= "The given Lagrangian contains non-renormalizable operators, which are currently not supported by DR2MS.";
DR2MS::noncanonicnorm= "The kinetic terms of the gauge bosons are not canonically normalized.";


DR2MS[lagrangian_]:= Module[
	{
		lagAux,
		lagTree,lagLoop,
		lag,
		LagGaugeKin,
		LagFermionMass,
		LagYukawa,
		LagQuartic,
		LagNonShifted
		,
		fieldContent,
		gaugeShifts,
		gaugeShiftsInverse={},
		fermionMassShift,
		yukawaShifts,
		fields,
		baredFields,
		quarticShifts
	}
	,
	(* ensure this function is only applied to renormalizable Lagrangians *)
	
	If[GetMaxEFTOrder[LagrangianExpand@ lagrangian, False]>4,
		Message[DR2MS::nonren];
		Abort[]
	];
	
	(* separate terms relevant for tree level matching *)
	{lagTree,lagLoop}= GetTreeLevelTerms[lagrangian];
	
	(* group by operators *)
	lag= CollectOperators[lagTree];
	
	(* get kinetic terms of all gauge bosons and check they are canonically normalized *)
	LagGaugeKin= KineticTerms[lag]/. _Field->0// OperatorToNormalForm;
	If[!KineticCanonicalQ[LagGaugeKin],
		Message[DR2MS::noncanonicnorm];
		Abort[]
	];
	LagGaugeKin= TermsToList[LagGaugeKin];
	
	(*  TermsToList@lag is problematic here since it expands lag before going to the List *)
	lag= If[Head[lag]===Plus,
		List@@lag,
		{lag}
	];
	
	
	(* get all Fermion mass terms *)
	LagFermionMass= Cases[
		lag/. _FieldStrength->0 /. Field[__,Except[{}]]->0,
		x_/;((Count[x,Field[_,Fermion,___],All]==2)&&(Count[x,Field[_,Except[Fermion],___],All]==0)),
		1
	];
	
	(* get all scalar-fermion-fermion terms *)
	LagYukawa= Cases[
		lag/. _FieldStrength->0,
		x_/;((Count[x,Field[_,Fermion,___],All]==2)&&(Count[RemovePower[x],Field[_,Scalar,___],All]==1)),
		1
	];
	
	(* get all scalar quartics *)
	(* NOTE: instead of strating from the present scalar quartics, it might be better to start from the gauge interactions in the availavle sclar kinetic terms! *)
	LagQuartic= Cases[
		lag/._FieldStrength->0,
		x_/;((Count[RemovePower[x],Field[_,Scalar,___],All]==4)&&(Count[x,Field[_,Except[Scalar],___],All]==0)),
		1
	];
	
	(* determine all terms that will not be shifted *)
	LagNonShifted= LagrangianExpand[(Plus@@lag)-(Plus@@LagGaugeKin+Plus@@LagFermionMass+Plus@@LagYukawa+Plus@@LagQuartic)];

(* compute scheme changes  *)
	(* proper gauge coupling shifts *)
	If[Length[LagGaugeKin]>0&&LagGaugeKin=!={0},
		(*gaugeShifts= (FirstCase[#,_Coupling,0,All]->FirstCase[#,_Coupling,0,All](1+hbar/6 SumC2g2[FirstCase[#,FieldStrength[l_,___]:>l,0,All]]))&/@LagGaugeKin;*)
		(*corrected sign*)
		gaugeShifts= (FirstCase[#,_Coupling,0,All]->FirstCase[#,_Coupling,0,All](1-hbar/6 SumC2g2[FirstCase[#,FieldStrength[l_,___]:>l,0,All]]))&/@LagGaugeKin;
		(* This shift leads to NON-canonically normalized kinetic terms for the gauge bosons *)
		LagGaugeKin= (First[#]/.Last[#])&/@Transpose@{LagGaugeKin,gaugeShifts};
		(* compute the inverse shift of the MSbar gauge coupling to regain a canonical normalization at the end *)
		gaugeShiftsInverse=gaugeShifts/.hbar->-hbar;
	];

	(* Fermion mass shifts *)
	If[Length[LagFermionMass]>0,
		fermionMassShift= hbar*SumC2g2[FirstCase[#,Field[l_,Fermion,___]:>l,0,All]]&/@LagFermionMass;
		LagFermionMass= (First[#]/.Last[#])&/@Transpose@{
			NormalToOperatorForm[LagFermionMass],
			(op_Operator:>op(1+#))&/@fermionMassShift
		};
	];

	(* Yukawa mass shifts *)
	If[Length[LagYukawa]>0,
		fieldContent= Join[Cases[#,Field[l_,Fermion,___]:>l,All],Cases[#,Field[l_,Scalar,___]:>l,All]]&/@LagYukawa;
		yukawaShifts= Simplify[hbar/2 (SumC2g2[#[[1]]]+SumC2g2[#[[2]]]-2SumC2g2[#[[3]]])]&/@fieldContent;
		(*LagYukawa= (First[#]/.Last[#])&/@Transpose@{NormalToOperatorForm[LagYukawa],(op_Operator:>op(1-#))&/@yukawaShifts};*)
		(*corrected sign*)
		LagYukawa= (First[#]/.Last[#])&/@Transpose@{NormalToOperatorForm[LagYukawa],(op_Operator:>op(1+#))&/@yukawaShifts};
		LagYukawa= OperatorToNormalForm[LagYukawa];
	];
	
	(* shifting scalar quartics *)
	If[Length[LagQuartic]>0,
		LagQuartic=    CollectOperators[LagQuartic, NormalForm->False]/.op_Operator:>UnifyOpIndices[op];
		fields=        Cases[#,Field[_,Scalar,___],All]&/@LagQuartic;(*all fields*)
		baredFields=   Cases[#,Bar[f:Field[_,Scalar,___]]:>f,All]&/@LagQuartic;(*all bared fields*)
		fields=        Complement[First@#,Last@#]&/@Transpose@{fields,baredFields};(*all non-bared fields*)
		fields=        Join[First@#,Last@#]&/@Transpose@{baredFields,fields};(*sort: {Bar@\[Phi],Bar@\[Phi],\[Phi],\[Phi]}*)
		quarticShifts= QuarticScalarShift[Sequence@@#]&/@fields;
		(*LagQuartic=    (First@#+(Last@#)*FirstCase[First@#,_Operator,All])&/@Transpose@{LagQuartic,quarticShifts};*)
		(*corrected sign*)
		LagQuartic=    (First@#-(Last@#)*FirstCase[First@#,_Operator,All])&/@Transpose@{LagQuartic,quarticShifts};
	];
	
	(* combine results *)
	Collect[
		ExpandDenominator[
			Contract@RelabelIndices[
				Plus@@LagGaugeKin+Plus@@LagFermionMass+Plus@@LagYukawa+ContractCGs[Plus@@LagQuartic]+LagNonShifted
			]/.gaugeShiftsInverse
		],
		hbar,
		CollectOperators
	]+lagLoop
]
