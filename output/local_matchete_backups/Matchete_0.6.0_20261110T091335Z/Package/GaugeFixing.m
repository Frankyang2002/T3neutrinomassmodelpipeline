(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`GaugeFixing`*)


(* ::Subtitle:: *)
(*Procedures for gauge-fixing, both ordinary, background, and partially fixed*)


(* ::Text:: *)
(*We mainly follow the general notation used in https://arxiv.org/abs/2404.11640*)


(* ::Chapter:: *)
(*Public:*)


(* ::Subsection:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["GaugeFixLagrangian"]
PackageExport["R\[Xi]"]
PackageExport["\[Xi]Gauge"]
PackageExport["Gauge"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["IndividualRXiGaugeFixingTerms"]


(* ::Subsection:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


GaugeFixLagrangian::usage= "GaugeFixLagrangian[<lag>, Gauge-> <choice>] returns the gauge fixed version of the Lagrangian. This is only supported for broken Lagrangians that are already at their vacuum (cf. ImplementVacuumConditions)."
\[Xi]Gauge::usage= "\[Xi]Gauge[V] is the gauge parameter of the R\[Xi] gauge (associated with the vector V)";
R\[Xi]::usage= "Used as an option value in Matchete";
Gauge::usage= "Is an option used in Matchete.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* ::Subsection:: *)
(*Background gauge-fixing (for functional calculations)*)


(* ::Subsection:: *)
(*Ordinary gauges*)


(* ::Subsubsection::Closed:: *)
(*GaugeFixLagrangian *)


(* ::Text:: *)
(*Function to return the gauge-fixed version of a Lagrangian. It only works for broken-phase Lagrangians if they have been expanded about their vacuum (ImplementVacuumConditions). *)


GaugeFixLagrangian::unkwngauge= "`1` is not a valid gauge.";
GaugeFixLagrangian::partGF= "The input Lagrangian looks to have already been (partially) gauge-fixed.";
GaugeFixLagrangian::partGFnonU= "The input Lagrangian appears to have already been partially gauge-fixed not using the unitary gauge This is curently not supported.";
GaugeFixLagrangian::unbrokenunitary= "Input Lagrangian found to be in the symmetric phase. Cannot go to the unitary gauge for proper (massless) gauge fields.";
GaugeFixLagrangian::alreadysubstituted= "The vector masses in the input Lagrangian are already written in terms of the mass symbols. SubstituteMasses-> False is not supported in this case; use Automatic or True.";
GaugeFixLagrangian::partGFmix= "The input Lagrangian appears to have already been partially gauge-fixed, but there is mixing between the fixed and the unfixed sector.";
GaugeFixLagrangian::partial= "The input Lagrangian appears to be already partially gauge-fixed with gauges `1` which are retained.";


Options@ GaugeFixLagrangian= {
		EFTOrder-> Automatic,
		Full-> True,
		Gauge-> R\[Xi],
		SubstituteMasses-> Automatic,
		Verbose-> True
	};


OptionMessage[EFTOrder, GaugeFixLagrangian, val_]:= Message[General::optexpectsval, EFTOrder, GaugeFixLagrangian, val, "positive Integer, 'Automatic', or 'None'"];
OptionMessage[Full, GaugeFixLagrangian, val_]:= Message[General::optexpectsval, Full, GaugeFixLagrangian, val, "'True' or 'False'"];
OptionMessage[Gauge, GaugeFixLagrangian, val_]:= Message[General::optexpectsval, Gauge, GaugeFixLagrangian, val, "'R\[Xi]' or 'Unitary'"];
OptionMessage[SubstituteMasses, GaugeFixLagrangian, val_]:= Message[General::optexpectsval, SubstituteMasses, GaugeFixLagrangian, val, "Automatic, True, or False"];
OptionMessage[Verbose, GaugeFixLagrangian, val_]:= Message[General::optexpectsval, Verbose, GaugeFixLagrangian, val, "'True' or 'False'"];


OptionTest[GaugeFixLagrangian, EFTOrder]= MatchQ[(_Integer? Positive) | Automatic];
OptionTest[GaugeFixLagrangian, Full]= BooleanQ;
OptionTest[GaugeFixLagrangian, Gauge]= MatchQ[Unitary | R\[Xi]];
OptionTest[GaugeFixLagrangian, SubstituteMasses]= MatchQ[Automatic| True| False];
OptionTest[GaugeFixLagrangian, Verbose]= BooleanQ;


GaugeFixLagrangian[lag_, OptionsPattern[]] ? OptionsCheck:= Module[
		{brokenPhase, gfVectors, gfGhosts, gfParam, individualGauges, nonGF, 
			gfLag= 0, maxOrder= OptionValue@ EFTOrder, workingLag= lag, workingLagGF, echo= If[OptionValue[Verbose], Echo, Identity]},
	If[!CheckLagrangian[workingLag, GaugeAnomalies-> False, GaugeInvariant-> False, HeavyMassBasis-> False],
		Abort[];
	];
	
	brokenPhase= BrokenPhaseLagQ@ lag;
	
	(* Depending on whether partial GF is allowed. 
	Partial GF is only allowed to be unitary gauge (and heavy vectors) *)
	If[OptionValue[Full] || !brokenPhase,
		(* abort if there is partial of full GF present *)
		If[!NonGaugeFixedLagQ@ workingLag,
			Message[GaugeFixLagrangian::partGF];
			Abort[];
		];
	,
		(*When partial GF is allowed, substract the terms which already have GF*)
		(*Check which fields are already gauge fixed*)
		{nonGF, individualGauges}= Quiet[NonGaugeFixedLagQ[workingLag, DetailedOutput-> True], 
			{NonGaugeFixedLagQ::nongaugekinfix, NonGaugeFixedLagQ::noghostkin}];
		If[!nonGF,
			{workingLag, gfLag, gfGhosts, gfParam}= 
				SubtractPreviouslyFixedTerms[workingLag, individualGauges, OptionValue@ Verbose];
		];
	];

	(*Determine the representation of the vector masses in the gauge-fixing terms (the gauge-fixing terms inherit the form of the masses in the input Lagrangian).*)
	If[brokenPhase,
		Switch[OptionValue@ SubstituteMasses
		,True, (*Express all masses via the field mass symbols*)
			workingLag= SubstituteMassesInLagrangian@ workingLag;
			gfLag= SubstituteMassesInLagrangian[gfLag, "CheckVectorMassConsistency"-> False] (* this part of the lagrangian is already gauge fixed *);
			
		,False, (*Leave the masses as couplings; unsupported on an already-substituted Lagrangian*)
			If[VectorMassesSubstitutedQ@ workingLag,
				Message[GaugeFixLagrangian::alreadysubstituted];
				Abort[];
			];
		(*Automatic: keep the form of the input Lagrangian*)
		];
	];

	(*Determine working EFT order and what expansion function should be applied to the g.f. terms.*)
	If[maxOrder === Automatic,
		maxOrder= GetMaxEFTOrder[lag, False];
	];

	(*Select the appropriate function for constructing the gauge*)
	workingLag= Switch[{OptionValue@ Gauge, brokenPhase}
	,{Unitary, True},
		UnitaryGauge[workingLag, maxOrder]
	,{Unitary, False},
		Message[GaugeFixLagrangian::unbrokenunitary];
		Abort[]
	,{R\[Xi], True},
		workingLag+ BrokenPhaseRXiTerms[workingLag, maxOrder]
	,{R\[Xi], False},
		workingLag+ UnbrokenRXiGaugeFixingTerms[workingLag]
	,_,
		Message[GaugeFixLagrangian::unkwngauge, OptionValue@ Gauge];
		Abort[];
	];
	
	(*Remove unwanted gauge fixing terms that were incorrectly introduced for the fields that were already gauge fixed*)
	If[!OptionValue[Full] && gfLag =!= 0,
		(* the code below avoids expanding workingLag such that one does not have to call CollectOperators or GreensSimplify in the end *)
		(* convert to list without expanding *)
		workingLag= TermsToList[workingLag, Expand-> False];
		(* separate gauge param dependent and independent terms for performance *)
		workingLagGF= Cases[workingLag, x_/; !FreeQ[x, Alternatives@@ gfParam], 1];
		workingLag= Plus@@ Cases[workingLag, x_/; FreeQ[x, Alternatives@@ gfParam], 1];
		(* drop gauge param dependence *)
		workingLagGF= workingLagGF/. Plus[A_, rest___]/; (FreeQ[A, _Plus] && !FreeQ[A, Alternatives@@ gfParam]):> Plus[rest];
		workingLagGF= workingLagGF/. Times[A_, rest___]/; (FreeQ[A, _Plus] && !FreeQ[A, Alternatives@@ gfParam])-> 0;
		(* recombine *)
		workingLag= (Plus@@ workingLagGF) + workingLag;
		(* both CollectOperators and GreensSimplify have severe performance issues below *)
		workingLag= workingLag+ gfLag/. Field[Alternatives@@ gfGhosts, ___]-> 0(*// CollectOperators*)(*GreensSimplify*);
	];
	 
	workingLag
];


(* ::Text:: *)
(*True if the massive-vector masses in the Lagrangian are already written in terms of the field mass symbols.*)


VectorMassesSubstitutedQ[lag_]:= AnyTrue[HeavyVectorLabels[], !FreeQ[lag, $FieldAssociation[#, Mass]]&];


SubtractPreviouslyFixedTerms[lag_, individualGauges_, verbose_?BooleanQ]:= Module[
		{gfVectors, gfGhosts, gfLag, gfParam, workingLag},
	
	(*Find previously gauge fixed fields*)
	gfVectors= Select[individualGauges, !MatchQ[#, Indeterminate]&];
	
	If[!ContainsOnly[Values[gfVectors], {Unitary}],
		Message[GaugeFixLagrangian::partGFnonU];
		Abort[];
	];
	
	If[verbose, 
		Echo@ StringForm[GaugeFixLagrangian::partial, gfVectors];
	];
	
	(*Get (anti-)ghosts and goldstones associated to the already fixed fields*)
	gfGhosts= Flatten[List@@ Lookup[$GhostFields, Keys@ gfVectors] ];
	(*Get gauge parameter to already fixed fields*)
	gfParam= \[Xi]Gauge/@ Keys[gfVectors];
	
	(*Substract the part of workingLag that is already GF*)
	workingLag= TermsToList@ lag;
	gfLag= Cases[workingLag, x_/; !FreeQ[x, Alternatives@@ Keys[gfVectors], All], 1];
	workingLag= lag - Total[gfLag];
	
	(*Check that there is no mixing between GF and non-GF sectors*)
	If[!FreeQ[
			DeleteCases[PseudoTimes/@ gfLag, x_/; (Count[x, _Field|_FieldStrength, All] > 2), 1], 
			Alternatives@@ Complement[Keys[individualGauges], Keys@ gfVectors], 
			All
		],
		Message[GaugeFixLagrangian::partGFmix];
		Abort[];
	];
	
	(*{other terms, gauge-fixed terms, ghosts, gauge parameters}*)
	{workingLag, Total@ gfLag, gfGhosts, gfParam}
]


(* ::Subsubsection::Closed:: *)
(*Unbroken R\[Xi] gauge*)


(* ::Text:: *)
(*Returns the R\[Xi] gauge-fixing and ghost terms for an unbroken Lagrangian*)
(*Should only be called on a Lagrangian in a symmetric phase *)


UnbrokenRXiGaugeFixingTerms[lag_]:= Module[{gfTerms},
	gfTerms= KeyValueMap[IndividualRXiGaugeFixingTerms, GetGaugeCouplingsFromLagrangian[lag, Power-> 1]];
	RelabelIndices@ Total@ gfTerms
];


(* ::Text:: *)
(*The gauge-fixing + ghost terms for a non-participating gauge group*)


IndividualRXiGaugeFixingTerms[gaugeField_, invCoupling_]:= Module[{rest, a, mu, nu},
	rest= Sequence@@ If[Length@ $FieldAssociation[gaugeField, Indices] === 1, {a}, {}];
	-1/(2 \[Xi]Gauge@ gaugeField)* invCoupling* PD[mu, gaugeField[mu, rest]]* PD[nu, gaugeField[nu, rest]]+ 
	If[Length@ {rest} =!= 0, (*Add ghosts if non-Abelian*)
		-I GrassmannProduct[$GhostFields[gaugeField][[2]][rest], PD[mu, CD[mu, $GhostFields[gaugeField][[1]][rest] ] ] ]
	,
		0
	]
]
