(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`SuperTrace`*)


(* ::Subtitle:: *)
(*Paclet for evaluating supertraces in the functional matching methods.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["PowerTypeSTr"]
PackageScope["LogTypeSTr"]


PackageScope["EvaluateSTr"]


(* ::Text:: *)
(*Options*)


PackageScope["Fields"]


(* ::Text:: *)
(*What are those??*)


PackageScope["Xop"]
PackageScope["Xterm"]
PackageScope["Mop"]
PackageScope["Mterm"]
PackageScope["GaugeCop"]
PackageScope["GaugeCTerm"]
PackageScope["WilsonLine"]
PackageScope["WilsonTerm"]
PackageScope["WilsonExpand"]
PackageScope["ExpandGenFSs"]


PackageScope["ActWithOpenCDs"]
PackageScope["CloseFermionLoop"]
PackageScope["WilsonExpand"]
PackageScope["ExpandGenFSs"]
PackageScope["RemoveSymmetryVanishingWilsonTerms"]


(* ::Text:: *)
(*Specific functions for Finite Temperature*)


PackageScope["PropSp"]
PackageScope["PropTp"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Chapter:: *)
(*SuperTrace*)


(* ::Section:: *)
(*Power-type traces*)


(* ::Text:: *)
(*Evaluation of power-type supertraces*)


(* ::Subsection:: *)
(*Evaluate power-type traces*)


(* ::Subsubsection::Closed:: *)
(*Power-type traces*)


Options@ PowerTypeSTr= {Fields-> All, Mode-> Matching};


PowerTypeSTr[propTypes_List, eftOrder_, OptionsPattern[]]:= Module[
		{maxEFTOrd, maxPropExpansionOrder, propOrder, result, genexp, insertions, preFactor},
	(* propTypes: list of propagator types as obtained by ListPowerTypeTraces. *)
	preFactor= -I hbar/ 2 Switch[First@ propTypes,
			hScalar| hVector| lScalar| lVector, 1,
			hFermion| hGhost| hAntiGhost| lFermion| lGhost| lAntiGhost, - 1
		];
	maxEFTOrd= (eftOrder/. List-> Identity);

	(* Determine maximal expansion order for propagators: MaxPropExpansionOrder. *)
	maxPropExpansionOrder = maxEFTOrd - Total[Partition[propTypes, 2, 1, 1]/. $XOrdMin];

	(* Expand and evaluate STr *)
	(*Sum exclusively over all possible propagator expansion orders*)
	result = Sum[
		(* Expand all propagators in STr to the given propOrder (exclusive) *)
		genexp = GenericPropagatorExpansion[propTypes, propOrder]; (* This could be cached *)
		
		(* Determine all allowed insertions of X-terms, Masses, ... such that the total EFT order of the STr is \[LessEqual] EFTOrder.
		Only return unique insertions and include count how often they are repeated. *)
		insertions = Switch[OptionValue@ Mode
			,Matching|DimReduction,
				(*Inclusive expansion up to max EFT order*)
				DeterminePowerInsertions[propTypes, maxEFTOrd- propOrder, OptionValue@ Fields]
			,Divergence,
				(*Exclusive or inclusive as determined by eftOrder*)
				DetermineDivergentInsertions[propTypes, eftOrder- propOrder,
					PropMomentumPower@ propTypes + propOrder - 4, OptionValue@ Fields]
			,Evanescent,
				(*Exclusive or inclusive as determined by eftOrder*)
				DetermineEvanescentInsertions[propTypes, eftOrder- propOrder,
					PropMomentumPower@ propTypes + propOrder - 4, OptionValue@ Fields]
			];
		(* Sum over all possible insertions and evaluate the traces *)
		Sum[
			EvaluateSTr[genexp, ins, propTypes, OptionValue@ Mode]
		,
			{ins, insertions}
		]
	,
		{propOrder, 0, maxPropExpansionOrder}
	];

	If[OptionValue@ Mode === Matching,
		result= ReplaceHeavyEOM[result, EFTOrder-> eftOrder];
	];


	(* ensure that the loop factor is multiplied into each term ensuring an expanded output.
	This is required in complicated models where result can become extremely long *)
	If[Head[result]===Plus,
		(preFactor*#)&/@result,
		preFactor*result
	]

]


(* ::Subsubsection::Closed:: *)
(*Power of loop momenta*)


(* ::Text:: *)
(*Returns the total inverse power of loop momenta associated with the propagator types (assumes massless fermions)*)


PowerTypeSTr::hFermionMom= "Cannot determine momenta associated with heavy fermion.";


(* ::Text:: *)
(*Counts the total momentum power carried by the *)


PropMomentumPower@ propTypes_List:= Total@ Replace[propTypes, {lFermion-> 1,
	hFermion:> (Message[PowerTypeSTr::hFermionMom]; Abort[]), _-> 2}, {1}];


(* ::Subsection:: *)
(*Generic power-type hard-region expansion of propagators*)


(* ::Text:: *)
(*Function to determine the generic hard-region expanded trace blue-print based on the propagator types*)
(*	propTypes: list of propagator types as obtained by GenerateSuperTraces.*)
(*	Expand each propagator in propTypes such that the overall EFT order of the STr is exactly expOrder.*)


(* ::Subsubsection::Closed:: *)
(*GenericPropagatorExpansion*)


GenericPropagatorExpansion[propTypes_List, expOrder_]:= Module[
		{gaugeCouplings, genInds, n, masses, props, Xops, types, orders, propNo= Length@ propTypes},
	genInds=  Table[Symbol["$i"<> ToString@ n], {n, propNo+ 1}];

	(*Create the list of X-op (+ WilsonLine) templates to go in template*)
	Xops= MapThread[Xop,
		{Partition[propTypes, 2, 1, 1], Partition[genInds, 2, 1], Table[n, {n, propNo}]}];
	AppendTo[Xops, WilsonLine[First@ propTypes, genInds[[{-1, 1}]]]];

	(*Determine the corresponding propagator masses *)
	types= RotateLeft[propTypes];
	masses= MapThread[Mop, {types, genInds[[2;;]], Table[n, {n, propNo}]}];
	(*And the gauge coupling insertions from gauge fields*)
	gaugeCouplings= Times@@ Table[If[types[[n]] === lVector, GaugeCop[n], 1], {n, propNo}];

	(*Sum over all ways to expand the propagators*)
	gaugeCouplings* Sum[
		props= MapThread[PropAtOrder, {types, masses, orders}];
		(*Interweave Xops and expanded propagators *)
		FuncNCM@@ (Riffle[Xops, props]/. PropAtOrder-> PropExpand)
	, {orders, IntegerSets[expOrder, propNo]}]// Expand
]


(*Return[
	GaugeCop[2]*... FuncNCM[Xop[{\[CapitalPhi],\[CapitalPhi]},{i,j},1], (Prop@Mop[\[CapitalPhi],j,1])^k, LoopMom[___], CD[___], .. , WilsonLine[\[CapitalPhi],{i,j}]]
	+GaugeCop[2]*... FuncNCM[Xop[{\[CapitalPhi],\[CapitalPhi]},{i,j},2], (Prop@Mop[\[CapitalPhi],j,2])^k, LoopMom[___], CD[___], .. , WilsonLine[\[CapitalPhi],{i,j}]]
	+ ...
]*)


(* ::Subsubsection::Closed:: *)
(*Expansion of propagators*)


(* ::Text:: *)
(*Determine expansion type*)


PropExpand[fType_, mass_, ord_, OptionsPattern[]]:= Switch[$FieldTypes[fType, Type]
	,Fermion,
		PropFermionExpand[mass, ord]
	,Vector,
		-PropBosonExpand[mass, ord]
	,_,
		PropBosonExpand[mass, ord]
]


(* ::Subsubsection::Closed:: *)
(*Expansion of boson propagator*)


(* ::Text:: *)
(*The expansion is the n'th EFT order term of 1/[(k+ P)^2 - M^2] = 1/[k^2 - M^2] \sum_{n=0} ( -1)^n ( [2 k.P + P^2] / [k^2 - M^2] )^n  *)


PropBosonExpand[mass_, 0]:= Prop@ mass;
PropBosonExpand[mass_, ord_]:= Module[{indices, m, set, singleCDs,  pairCDs},
	indices= Index[#, Lorentz]&/@ Table[Unique@ "mu", {m, ord}];

	(*Sum over the number of D^2 insertions*)
	Power[-I, ord] Sum[
		Power[-1, m]* Power[2, ord- 2m]* Power[Prop@ mass, ord+ 1- m]*
		 Times@@ LoopMom/@ indices[[;;-2 m-1]]*
		(*Sum the ways to position the D^2 insertions*)
		Sum[
			singleCDs= OpenCD/@ List/@ indices[[;;-2 m-1]];
			pairCDs= FuncNCM@@ OpenCD/@ {{#}, {#}} &/@  indices[[-2m;; -m-1]];
			pairCDs= FuncNCM@@@ VariableLengthPartition[pairCDs, set];
			(*interweave the list of D_mu's and D^2's*)
			FuncNCM@@ Riffle[pairCDs, singleCDs]
		, {set, IntegerSets[m, 1+ord- 2m]}]
	, {m, 0, Floor[ord/2]}]
]


(* ::Subsubsection::Closed:: *)
(*Expansion of fermion propagator*)


(* ::Text:: *)
(*The expansion is the n'th EFT order term of 1/[(k.\[Gamma]+ P.\[Gamma]) - M] = [k.\[Gamma] + P.\[Gamma]+ M] [k^2 +2 k.P +P.\[Gamma] P.\[Gamma] - M^2]*)
(*=    [k.\[Gamma] + P.\[Gamma]+ M]/[k^2 - M^2] \sum_{n=0} ( -1)^n ( [2 k.P + P.\[Gamma] P.\[Gamma]] / [k^2 - M^2] )^n *)


PropFermionExpand[mass_, ord_]:= Module[{ind= Index[Unique@ "mu", Lorentz]},
	FuncNCM[(DiracProduct@ GammaM@ LoopMom+ mass), PropFermionExpandHelper[mass, ord]]+
	I* FuncNCM[DiracProduct@ GammaM@ ind, OpenCD@ {ind}, PropFermionExpandHelper[mass, ord -1]]
]


PropFermionExpandHelper[mass_, -1]:= 0;
PropFermionExpandHelper[mass_, 0]:= Prop@ mass;
PropFermionExpandHelper[mass_, ord_]:= Module[{indices, m, set, singleCDs,  pairCDs},
	indices= Index[#, Lorentz]&/@ Table[Unique@ "mu", {m, ord}];

	(*Sum over the number of D^2 insertions*)
	Power[-I, ord] Sum[
		Power[-1, m]* Power[2, ord- 2m]* Power[Prop@ mass, ord+ 1- m]*
		 Times@@ LoopMom/@ indices[[;;-2 m-1]]*
		(*Sum the ways to position the D^2 insertions*)
		Sum[
			singleCDs= OpenCD/@ List/@ indices[[;;-2 m-1]];
			pairCDs= FuncNCM@@@ Partition[FuncNCM[DiracProduct@ GammaM@ #, OpenCD@ {#}]&/@ indices[[-2m;;]], 2];
			pairCDs= FuncNCM@@@ VariableLengthPartition[pairCDs, set];
			(*interweave the list of D_mu's and D^2's*)
			FuncNCM@@ Riffle[pairCDs, singleCDs]
		, {set, IntegerSets[m, 1+ ord- 2m]}]
	, {m, 0, Floor[ord/2]}]
]


(* ::Subsection:: *)
(*Determine insertions of X *)


(* ::Subsubsection::Closed:: *)
(*Determine matching insertions*)


(* ::Text:: *)
(*Determine all insertions (with total EFT order <= order) for Xop, Mop, and Wilson line.	 For equivalent insertions only return one of them and the count.*)
(*propTypes: list of propagator types as obtained by ListPowerTypeTraces. *)
(*order: EFT order to which the insertions are allowed to sum maximally.*)


DeterminePowerInsertions[propTypes_List, order_, propFields_]:= Module[
		{fieldPattern, insertionOrders, newOrds, ordSet, sample, XSamples, Xords, Xtypes, 
		count, aux, props},

	Xtypes= Partition[propTypes, 2, 1, 1];

	props = propTypes/. $XFieldDofs;
	Xords = Min/@ Map[First, $XOrders, {2}];

	(* figure out all possible Xterm combinations that are allowed to be inserted up to the given order in the EFT power counting *)
	XSamples = Table[
			(* one of the possible combination of propagators *)
			aux = Table[props[[j]][[count[j]]], {j, Length[props]}];
			(* partition the list to obtain Xterms *)
			aux = Partition[aux, 2, 1, 1];
	
			(* drop the term if it is too high in the EFT power counting *)
			If[Total[aux/. Xords]> order, Nothing, aux] 

			(* loop over all fields for each propagator type *)
		, Evaluate[Sequence@@ Table[{count[i], 1, Length[props[[i]]]}, {i, Length[props]}]] ];
	
	(* Flatten the nexted list returned by Table*)
	XSamples = Flatten[XSamples, Length@ propTypes- 1];

	(*Remove fields if a particular loop is considered*)
	If[propFields =!= All,
		fieldPattern= Partition[propFields, 2, 1, 1];
		XSamples= Cases[XSamples, x_/; MatchQ[x/. Conj-> Identity, fieldPattern]];
	];

	If[Length@ XSamples === 0, Return@ {};];

	(*List all concrete combinations Xterms+order not greater than "order"*)
	insertionOrders= Select[Tuples[#/. $XOrders], Function[x, (Total@ x[[;;, 1]] <= order)]]&/@ XSamples;
	XSamples= Flatten[MapThread[Function[{x, y}, Transpose/@ Thread[{x, y}, List, {2}]],
		{XSamples, insertionOrders}], 1];

	(*Delete cyclically identical insertions and add sym. factor*)
	XSamples= {CyclicSymFactor@ #, #}&/@ DeleteDuplicatesByCylcicity@ XSamples;
	
	(*Apply possible expansions for the open CDs. Super traces are only cyclic before going to momentum space.*)
	XSamples= Join@@ Table[
			newOrds= PopulateCovMomOps[order, sample[[2, ;;, 2]]];
			Table[
				sample[[2, ;;, 2]]= ordSet;
				sample
			, {ordSet, newOrds}]
		, {sample, XSamples}];
	
	(*Create replacement patterns*)
	(*Xterm format: Xterm[{field labels}, {indices}, {base EFT ord, # cov. mom. ops., # of those being open CDs}]*)
	MapAt[Function[Xs, Flatten@ MapIndexed[
		{Xop[Xtypes[[First@ #2]], {i_, j_}, First@ #2]-> Xterm[First@ #1, {i, j}, Sequence@@ Last@ #1],
		If[$FieldTypes[Xtypes[[First@ #2, 2]], Heavy],
			Mop[Xtypes[[First@ #2, 2]], i_, First@ #2]-> Mterm[#1[[1, 2]], i],
			Nothing],
		If[Xtypes[[First@ #2, 2]] === lVector,
			GaugeCop[First@ #2]-> GaugeCTerm@ #1[[1, 2]],
			Nothing],
		If[First@ #2=== Length@ propTypes,
			WilsonLine[Xtypes[[First@ #2, 2]], {i_, j_}]-> WilsonTerm[#1[[1, 2]], {i, j}, {}]
		, Nothing]} &,
		Xs]], XSamples, {All, 2}]
];


(* ::Text:: *)
(*Adjust for cyclic permutations*)


DeleteDuplicatesByCylcicity@ list_:= DeleteDuplicatesBy[list,
	First@ Sort@ NestList[RotateRight, #, Length@ #- 1]& ];
CyclicSymFactor@ list_:= Length@ DeleteDuplicates@ NestList[RotateRight, list, Length@ list- 1]/
	Length@ list;


(* ::Text:: *)
(*Expansion of the covariant momentum operators (up to max EFT order)*)


PopulateCovMomOps[maxOrd_, orders_]:= Block[{expOrd, openCDs},
	expOrd= maxOrd- Total@ orders[[;;, 1]];
	openCDs= Select[Tuples[Range[0, #]&/@ orders[[;;, 2]]], Total@ # <= expOrd&];
	Transpose@ Join[Transpose@ orders, #]&/@ List/@ openCDs
];


(* ::Subsubsection::Closed:: *)
(*Determine divergent insertions*)


(* ::Text:: *)
(*Determine all insertions (with total eftOrder and momOrder ) for Xop,  and Wilson line.*)
(*For equivalent insertions only return one of them and the count.*)


DetermineDivergentInsertions[propTypes_List, eftOrder_, momOrder_, propFields_]:= Module[
		{evaSubs, fieldPattern, insertionOrders, newOrds, ordSet, sample, XSamples, Xords, 
			Xtypes, maxEFTOrd, props, aux, count, myList},
	maxEFTOrd= (eftOrder/. List-> Identity);

	(*List all possible fields for the X terms*)
	Xtypes= Partition[propTypes, 2, 1, 1];
	props = propTypes/. $XFieldDofs;
	Xords= Min/@ First/@ Transpose/@ $XOrders;

	(* figure out all possible Xterm combinations that are allowed to be inserted up to the given order in the EFT power counting *)
	XSamples = Table[
		(* one of the possible combination of propagators *)
		aux = Table[props[[j]][[count[j]]], {j, Length[props]}];
		(* partition the list to obtain Xterms *)
		aux = Partition[aux, 2, 1, 1];

		(* drop the term if it is too high in the EFT power counting *)
		aux = If[Total[(aux/. Xords)]> eftOrder, Nothing, aux];

		(* drop the term if it is too high in the loop momentum counting *)
		(* ... to do ... *)
		aux
		,
		(* loop over all fields for each propagator type *)
		Evaluate[Sequence@@ Table[{count[i], 1, Length[props[[i]]]}, {i, Length[props]}]]
	];
	(* Flatten the nested list returned by Table*)
	XSamples = Flatten[XSamples, Length@ propTypes -1];

	(*Remove fields if a particular loop is considered*)
	If[propFields =!= All,
		fieldPattern= Partition[propFields, 2, 1, 1];
		XSamples= Cases[XSamples, x_/; MatchQ[x/. Conj-> Identity, fieldPattern]];
	];

	If[Length@ XSamples === 0, Return@ {};];

	(*List all concrete combinations Xterms not greater than eftOrder or momOrder*)
	insertionOrders= If[Head@ eftOrder === List,
			Select[Tuples[#/. $XOrders],
				Function[x, (Total[x[[;;, ;;2]], 2] === maxEFTOrd+ momOrder && Total@ x[[;;, 2]] >= momOrder)]]&/@ XSamples
		,
			Select[Tuples[#/. $XOrders],
				Function[x, (Total[x[[;;, ;;2]], 2] <= maxEFTOrd+ momOrder && Total@ x[[;;, 2]] >= momOrder)]]&/@ XSamples
		];
	XSamples= Flatten[MapThread[Function[{x, y}, Transpose/@ Thread[{x, y}, List, {2}]],
		{XSamples, insertionOrders}], 1];

	(*Delete cyclically identical insertions and add sym. factor*)
	XSamples= {CyclicSymFactor@ #, #}&/@ DeleteDuplicatesByCylcicity@ XSamples;

	(*Apply possible expansions for the open CDs. Super traces are only cyclic before going to momentum space.*)
	XSamples= Join@@ Table[
			newOrds= PopulateCovMomOpsLoopExact[momOrder, sample[[2, ;;, 2]]];
			Table[
				sample[[2, ;;, 2]]= ordSet;
				sample
			, {ordSet, newOrds}]
		, {sample, XSamples}];

	(*Create replacement patterns*)
	(*Xterm format: Xterm[{field labels}, {indices}, {base EFT ord, # cov. mom. ops., # of those being open CDs, ev. T/F}]*)
	MapAt[Function[Xs, Flatten@ MapIndexed[
		{Xop[Xtypes[[First@ #2]], {i_, j_}, First@ #2]-> Xterm[First@ #1, {i, j}, Sequence@@ Last@ #1],
		If[Xtypes[[First@ #2, 2]] === lVector,
			GaugeCop[First@ #2]-> GaugeCTerm@ #1[[1, 2]],
			Nothing],
		If[First@ #2=== Length@ propTypes,
			WilsonLine[Xtypes[[First@ #2, 2]], {i_, j_}]-> WilsonTerm[#1[[1, 2]], {i, j}, {}]
		, Nothing]} &,
		Xs]], XSamples, {All, 2}]
];


(* ::Text:: *)
(*Expansion of the covariant momentum operators (keeping a specific number of loop momenta)*)


PopulateCovMomOpsLoopExact[momOrder_, orders_]:= Block[{expOrd, openCDs},
	expOrd= Total@ orders[[;;, 2]]- momOrder;
	openCDs= Select[Tuples[Range[0, #]&/@ orders[[;;, 2]]], Total@ # === expOrd&];
	Transpose@ Insert[Transpose@ orders, #, 3]&/@ openCDs
];


(* ::Subsubsection::Closed:: *)
(*Determine evanescent insertions*)


(* ::Text:: *)
(*Determine all evanescent insertions (with total eftOrder and momOrder ) for Xop,  and Wilson line.*)
(*The first Xterm is always assumed to be of evanescent type*)
(*For equivalent insertions only return one of them and the count.*)


DetermineEvanescentInsertions[propTypes_List, eftOrder_, momOrder_, propFields_]:= Module[
		{evaSubs, fieldPattern, insertionOrders, XSamples, Xords, evOrds, Xtypes, maxEFTOrd},
	maxEFTOrd= (eftOrder/. List-> Identity);

	(*List all possible fields for the X terms*)
	Xtypes= Partition[propTypes, 2, 1, 1];
	XSamples= Partition[#, 2, 1, 1]&/@ Tuples[propTypes/. $XFieldDofs];

	(*Remove fields if a particular loop is considered*)
	If[propFields =!= All,
		fieldPattern= Partition[propFields, 2, 1, 1];
		XSamples= Cases[XSamples, x_/; MatchQ[x/. Conj-> Identity, fieldPattern]];
	];

	(*Selecting only possible fields*)
	evaSubs= DeleteCases[Select[Last]/@ $XOrders, {}];
	XSamples= Select[XSamples, MatchQ[Alternatives@@ Keys@ evaSubs]@* First];

	{evOrds, Xords}= {Min/@ First/@ Transpose/@ evaSubs, Min/@ First/@ Transpose/@ $XOrders};
	XSamples= Select[XSamples, (First@ #/. evOrds) + Total@ (#[[2;;]]/. Xords) <= maxEFTOrd &] ;
	If[Length@ XSamples === 0, Return@ {};];

	(*List all concrete combinations Xterms not greater than eftOrder or momOrder*)
	insertionOrders= If[Head@ eftOrder === List,
			Select[Tuples@ Prepend[#[[2;;]]/. $XOrders, #[[1]]/. evaSubs],
				Function[x, (Total[x[[;;, ;;2]], 2] === maxEFTOrd+ momOrder && Total@ x[[;;, 2]] >= momOrder)]]&/@ XSamples
		,
			Select[Tuples@ Prepend[#[[2;;]]/. $XOrders, #[[1]]/. evaSubs],
				Function[x, (Total[x[[;;, ;;2]], 2] <= maxEFTOrd+ momOrder && Total@ x[[;;, 2]] >= momOrder)]]&/@ XSamples
		];
	XSamples= Flatten[MapThread[Function[{x, y}, Transpose/@ Thread[{x, y}, List, {2}]],
		{XSamples, insertionOrders}], 1];

	(*Delete cyclically identical insertions and add sym. factor*)
	XSamples= {CyclicSymFactor@ #, #}&/@ DeleteDuplicatesByCylcicity@ XSamples;

	(*Apply possible expansions for the open CDs. Super traces are only cyclic before going to momentum space.*)
	XSamples= Join@@ Table[
			newOrds= PopulateCovMomOpsLoopExact[momOrder, sample[[2, ;;, 2]]];
			Table[
				sample[[2, ;;, 2]]= ordSet;
				sample
			, {ordSet, newOrds}]
		, {sample, XSamples}];

	(*Create replacement patterns*)
	(*Xterm format: Xterm[{field labels}, {indices}, {base EFT ord, # cov. mom. ops., # of those being open CDs, ev. T/F}]*)
	MapAt[Function[Xs, Flatten@ MapIndexed[
		{Xop[Xtypes[[First@ #2]], {i_, j_}, First@ #2]-> Xterm[First@ #1, {i, j}, Sequence@@ Last@ #1],
		If[Xtypes[[First@ #2, 2]] === lVector,
			GaugeCop[First@ #2]-> GaugeCTerm@ #1[[1, 2]],
			Nothing],
		If[First@ #2=== Length@ propTypes,
			WilsonLine[Xtypes[[First@ #2, 2]], {i_, j_}]-> WilsonTerm[#1[[1, 2]], {i, j}, {}]
		, Nothing]} &,
		Xs]], XSamples, {All, 2}]
];


(* ::Section:: *)
(*Log traces *)


(* ::Text:: *)
(*Evaluation of log-type supertraces*)


(* ::Subsection:: *)
(*Evaluate log-type traces*)


(* ::Subsubsection::Closed:: *)
(*Log-type traces*)


(* ::Text:: *)
(*Evaluate the particular log-type supertrace as indicated by the propagator class *)


Options@ LogTypeSTr= {Fields-> All, Mode-> Matching};


(*Inclusive*)
LogTypeSTr[propType_, order_Integer, opt:OptionsPattern[]]:= Module[{m},
	Sum[LogTypeSTr[propType, {m}, opt], {m, order}]
];

(*Exclusive*)
LogTypeSTr[propType_, {eftOrder_}, OptionsPattern[]]:= Module[
		{result, genexp, insertions, preFactor},

	preFactor= I hbar/ 2 Switch[propType,
			hScalar| hVector| lScalar| lVector, 1,
			hFermion| hGhost| hAntiGhost| lFermion| lGhost| lAntiGhost, - 1
		];

	(* Determine all allowed insertions of X-terms, Masses, ... such that the total EFT order of the STr is \[LessEqual] EFTOrder.
	Only return unique insertions and include count how often they are repeated. *)
	insertions = DetermineLogInsertions[propType, OptionValue@ Fields];
	If[Length@ insertions === 0, Return@ 0;];

	(* Expand all propagators in STr to the given propOrder (exclusive) *)
	genexp = GenericLogExpansion[propType, eftOrder]; (* This could be cached *)
	(* Sum over all possible insertions and evaluate the traces *)
	result= Sum[
		EvaluateSTr[genexp, ins, {propType}, OptionValue@Mode]
	, {ins, insertions}];

	preFactor result
]


(* ::Subsubsection::Closed:: *)
(*Determine log insertions *)


(* ::Text:: *)
(*Determine the log insertions corresponding to the propagator type*)


DetermineLogInsertions[propType_, propFields_]:= Module[
		{f, lagFields, chargedFields, fields, factors, replacements},
	lagFields= Lookup[$XFieldDofs, propType, {}];

	(*Filter out uncharged fields*)
	chargedFields= Union[
		FieldByProperty[Sequence@@ Normal@ $FieldTypes@ propType, Charges-> {__}],
		FieldByProperty[Sequence@@ Normal@ $FieldTypes@ propType, Indices-> inds_/;
			IntersectingQ[GroupFromRep/@ inds, Keys@ $GaugeGroups]]
		];
	fields= Cases[lagFields, f_/; MemberQ[chargedFields, f/. Conj-> Identity]];

	(*Remove fields if a particular loop is considered*)
	If[propFields =!= All,
		fields= Cases[fields, f_/; MatchQ[f/. Conj-> Identity, First@ propFields]];
	];

	If[Length@ fields === 0, Return@ {};];

	replacements= {
		If[MemberQ[{hScalar, hFermion, hVector, hGhost, hAntiGhost}, propType], Mop[propType, i_, 1]-> Mterm[#, i], Nothing]
		,
		WilsonLine[propType, {i_, j_}] -> WilsonTerm[#, {i, j}, {}]
	}&/@ fields;

	(*Introduce a factor of 1/2 for chiral fermion fields, to compensate for the projectors in the trace*)
	factors= Table[If[MemberQ[{LeftHanded, RightHanded}, $FieldAssociation[f, Chiral]], 1/2, 1],
		{f, fields/. Conj-> Identity}];

	Transpose@ {factors, replacements}
]


(*Return[{
	{1 (*count*), {Mop[\[CapitalPhi],j,1]->Mterm[h,j], .., WilsonLine[\[CapitalPhi],{i,j}]->WilsonTerm[\[Phi],{i,j},{}]}},
	.. (* other insertions *)
}]*)


(* ::Subsection:: *)
(*Generic log- type hard-region expansion of the log term*)


(* ::Text:: *)
(*Function to determine the generic hard-region expanded log blue-print based on field type*)
(*	Expand the log such that the overall EFT is exactly expOrder.*)


(* ::Subsubsection::Closed:: *)
(*GenericLogExpansion*)


GenericLogExpansion[propType_, expOrder_] := Module[
		{genInd, n, masses, prop, Xops, types, orders, propNo= Length@ propTypes},
	genInd= Symbol["$i"<> ToString@ 1];

	prop= PropLogExpand[propType, Mop[propType, genInd, 1], expOrder];
	FuncNCM[prop, WilsonLine[propType, {genInd, genInd}]]
]


(* ::Subsubsection::Closed:: *)
(*Expansion of log*)


(* ::Text:: *)
(*Determine expansion type*)


PropLogExpand[fType_, mass_, ord_, OptionsPattern[]]:= Switch[$FieldTypes[fType, Type]
	,Fermion,
		LogFermionExpand[mass, ord]
	,_,
		LogBosonExpand[mass, ord]
]


(* ::Subsubsection::Closed:: *)
(*Expansion of boson log*)


(* ::Text:: *)
(*The expansion is the n'th EFT order term of log[(k+ P)^2 - M^2] = cst. - \sum_{n=1} ( -1)^n ( [2 k.P + P^2] / [k^2 - M^2] )^n  *)


LogBosonExpand[mass_, 0]:= 0;
LogBosonExpand[mass_, ord_]:= Module[{indices, m, set, singleCDs,  pairCDs},
	indices= Index[#, Lorentz]&/@ Table[Unique@ "mu", {m, ord}];

	(*Sum over the number of D^2 insertions*)
	Power[-I, ord] Sum[
		Power[-1, m+1]/ (ord- m)* Power[2, ord- 2m]* Power[Prop@ mass, ord- m]*
		 Times@@ LoopMom/@ indices[[;;-2 m-1]]*
		(*Sum the ways to position the D^2 insertions*)
		Sum[
			singleCDs= OpenCD/@ List/@ indices[[;;-2 m-1]];
			pairCDs= FuncNCM@@ OpenCD/@ {{#}, {#}} &/@  indices[[-2m;; -m-1]];
			pairCDs= FuncNCM@@@ VariableLengthPartition[pairCDs, set];
			(*interweave the list of D_mu's and D^2's*)
			FuncNCM@@ Riffle[pairCDs, singleCDs]
		, {set, IntegerSets[m, 1+ ord- 2m]}]
	, {m, 0, Floor[ord/2]}]
]


(* ::Subsubsection::Closed:: *)
(*Expansion of fermion log*)


(* ::Text:: *)
(*The n'th term of the squared fermion log reads tr log[(k.\[Gamma]+ P.\[Gamma]) - M] = cst.  - 1/2 \sum_{n=1} ( -1)^n ( [2 k.P + P^2 - i \sigma_{ab} P_a P_b] / [k^2 - M^2] )^n*)


LogFermionExpand[mass_, 0]:= 0;
LogFermionExpand[mass_, ord_]:= Module[{indices, m, set, singleCDs,  pairCDs},
	indices= Index[#, Lorentz]&/@ Table[Unique@ "mu", {m, ord}];

	(*Sum over the number of D^2 insertions*)
	Power[-I, ord] Sum[
		Power[-1, m+1]/ (ord- m)* Power[2, ord- 2m- 1]* Power[Prop@ mass, ord- m]*
		 Times@@ LoopMom/@ indices[[;;-2 m-1]]*
		(*Sum the ways to position the D^2 insertions*)
		Sum[
			singleCDs= OpenCD/@ List/@ indices[[;;-2 m-1]];
			pairCDs= FuncNCM@@@ Partition[FuncNCM[DiracProduct@ GammaM@ #, OpenCD@ {#}]&/@ indices[[-2m;;]], 2];
			pairCDs= FuncNCM@@@ VariableLengthPartition[pairCDs, set];
			(*interweave the list of D_mu's and D^2's*)
			FuncNCM@@ Riffle[pairCDs, singleCDs]
		, {set, IntegerSets[m, 1+ ord- 2m]}]
	, {m, 0, Floor[ord/2]}]
]


(* ::Section:: *)
(*Shared evaluation functions*)


(* ::Subsection:: *)
(*Main evaluation of supertraces*)


(* ::Subsubsection::Closed:: *)
(*Evaluate Supertrace*)


(* ::Text:: *)
(*Function for performing a supertrace given a generic propagator expansion and an X-term substitution*)


EvaluateSTr[expr_, {factor_, replacement_}, propTypes_, mode_]:= Module[{out,propcount},
	out= factor expr/. replacement/. $Xsubs/. $Msubs/. $Gsubs;
	
	(*Modify propagators if we are working with a spatial Lagrangian*)
	If[!FreeOfZeroModes[out],
		out= out/.{Prop-> PropSp, Lorentz-> Spatial};
	];
	
	out= ActWithOpenCDs@ out/. FuncNCM-> NCM;
	
	Switch[mode,
		DimReduction,
		out= out//ExtractMomenta;
		out= CloseFermionLoop[out, propTypes];	
		out = out// CollectMomentaSpatial//CanonizeSpatial//WilsonToStatic// RemoveSymmetryVanishingWilsonTerms;
	,
		_,
		out= out// GatherLoopMomenta// RemoveSymmetryVanishingWilsonTerms;
		out= CloseFermionLoop[out, propTypes];
		];
	
	out= out// EvaluateSymmetricLorentzInds//ContractMetric// WilsonExpand;
		
	Switch[mode,
		Matching,
		out= LoopIntegrate[out];
	,
		DimReduction,
		out= If[MatchQ[First@ propTypes, hFermion| lFermion], $FSumInt out, out];
		out= LoopSumIntegrate[out];
	,
		_,		 
		out= Expand[out /. Prop[x_] :> propcount Prop[x]];
		out= Expand[I/\[Epsilon] Coefficient[out, propcount^2]/. _Prop ->1]
	];
	
	out= out// RelabelIndices// ExpandGenFSs// ContractDelta;
	
	Switch[mode,
		DimReduction,
		out= out // RefineDiracProducts;
	,
		_,
		out= out// ContractCGs// ContractDelta// RefineDiracProducts;
	];
	
	If[mode === Divergence,
		EpsExpand[out, Order-> -1]
	,
		out// EpsExpand
	]// RelabelIndices
]


(* ::Subsubsection::Closed:: *)
(*Identify vanishing WilsonTerms*)


(* ::Text:: *)
(*A fully symmetric tensor from the CDs contracting into the WilsonTerm vanish identically.*)
(*Assumes one WilsonTerm and at most one SymmetricLorentzInds in each term. It should be called after acting with all CDs.*)


RemoveSymmetryVanishingWilsonTerms@ expr_:= Module[{out= BetterExpand@ expr, symInds, wilsonInds},
	If[Head@ out === Plus, Return[RemoveSymmetryVanishingWilsonTerms/@ out];];
	symInds= FirstCase[out, SymmetricLorentzInds@ inds___-> {inds}, {}, All];
	wilsonInds= FirstCase[out, WilsonTerm[__, inds_]-> inds, {}, All];
	If[Length@ wilsonInds === 0, Return@ out;];

	(*Checks for overlap between wilsonInds and symInds*)
	If[SubsetQ[symInds, wilsonInds] || MatchQ[wilsonInds, {a_, a_}], 0, out]
]


(* ::Subsubsection::Closed:: *)
(*Close fermion lines*)


(* ::Text:: *)
(*Traces fermion line if necessary and performs Dirac product simplifications.*)


CloseFermionLoop[expr_, propagators_]:= Module[{out= expr},
	out= If[MatchQ[First@ propagators, hFermion| lFermion],
			FermionTrace@ out,
			out
		];
	ContractMetric@ CanonizeSpinorLines@ out
]


(* ::Subsubsection::Closed:: *)
(*Act with the open covariant derivatives*)


(* ::Text:: *)
(*Terminate all open CDs to the right in the functional expression. *)


MapThread


(*
ActWithOpenCDs@ expr_:= expr//. {
	FuncNCM[a___, cd:OpenCD[{mu_}], x:Except[_OpenCD], b___]:>
		FuncNCM[a, x, cd, b] + FuncNCM[a, CD[mu, x], b],
	FuncNCM[a___, _OpenCD]:> 0
};
*)


ActWithOpenCDs@ expr_:= expr//. {
	FuncNCM[a___, OpenCD[{mu_}], b___/;FreeQ[{b},_OpenCD]]:> Sum[
		FuncNCM[a, Sequence@@ MapAt[CD[mu,#]&, {b}, i]]
		,
		{i, Length@{b}}
	],
	FuncNCM[a___, _OpenCD]:> 0
};


(* ::Subsection:: *)
(*Wilson lines*)


(* ::Text:: *)
(*WilsonTerms are understood as WilsonTerm[\[Phi], {i, j}, {\[Mu], \[Nu], ...}] = D_x^\[Mu] D_x^\[Nu]... U^i_j(x, y) |_{y=x}, where index i transforms in the representation of field \[Phi], and j in the conjugate representation. *)


(* ::Subsubsection::Closed:: *)
(*Derivatives on the Wilson lines (parallel displacement propagators)*)


(* ::Text:: *)
(*Compute derivatives of all Wilson lines in an expression in the coincidence limit*)


WilsonExpand@ expr_:= expr/. w_WilsonTerm:> WilsonTermExpand@@ w;


(* ::Text:: *)
(*Memoize in case of more than six derivatives*)


WilsonTermExpand[field_, inds:{_,_}, devInds: {___}]:= If[Length@devInds>4,
	 With[{devPats= devInds/. Index[a_, t_]:> Index[Pattern[#, Blank[]]&@ a, t], indPats= Pattern[#, Blank[]]&/@ inds},
		WilsonTermExpand[field, indPats, devPats] = WilsonTermExpandEvaluate[field,inds,devInds]
		],
		
		WilsonTermExpandEvaluate[field,inds,devInds]
	]


(* ::Text:: *)
(*Compute derivatives acting on single line, depending on the gauge representation of the underlying object*)


WilsonTermExpandEvaluate[field_, {ind1_, ind2_}, devInds_List]:= Module[
		{conj, devSet, fieldCharges, fsSum, indices, fieldLabel, flavorIndices,
		gaugeIndSet, gaugeIndices, flavorDeltas, fieldStrengthFactor, lorentz},
	fieldLabel= If[(conj= MatchQ[field, _Conj]), First@ field, field];
	lorentz= If[$FieldAssociation[fieldLabel, Type] === Vector, Metric[ind1, ind2], 1];
	indices= If[conj, Bar, Identity]@ Map[{Index[ind1, #1], Bar@ Index[ind2, #1]} &,
		$FieldAssociation[fieldLabel, Indices]];

	(*Returns the trivial line*)
	If[Length@ devInds === 0,
		Return[lorentz* Times@@ Delta@@@ indices];
	];

	(*Determines gauge and flavor indices*)
	gaugeIndices= Cases[indices, _? (MemberQ[Keys@ $GaugeGroups, GroupFromInd@ First@ #] &)];
	flavorIndices= Complement[indices, gaugeIndices];
	fieldCharges= If[conj, MapAt[Minus, #, {All, 1}]&, Identity]@ $FieldAssociation[fieldLabel, Charges];

	flavorDeltas= Times@@ Delta@@@ flavorIndices;

	(*Determines the sum over various field strength tensors*)
	fieldStrengthFactor= Sum[
		gaugeIndSet= GaugeIndexSet[gaugeIndices, Length@ devSet];
		DevPreFact@ devSet* ContractDelta[Times@@
			MapThread[FSWilsonFactor[#1, #2, fieldCharges]&, {devSet, gaugeIndSet}]]
	, {devSet, DerivativeSubLists@ devInds}];

	lorentz* flavorDeltas* fieldStrengthFactor
]


DevPreFact[indSets_]:= Times@@ (-I (Length@ # -1)/ Length@ #&)/@indSets;


(* ::Text:: *)
(*Determines all index sets relevant to the computation of the CDs on the Wilson line. Returns a List of all ways to partition the original list in sublists.*)


DerivativeSubLists@ l_List:=Block[{sets, len= Length@ l},
	If[len< 2, Return@ {}; ];
	If[len< 4, Return@ {{l}}; ];

	(*Constructs a list of partitions removing a valid subset from the original list*)
	sets= Subsets[Range[len -1], {2, len -2}] ~Join ~{{}};
	sets= PartDrop[l, #]&/@ sets;

	(*If possible make further partitions of the first subset of each set*)
	Join@@ Table[
		If[!MatchQ[set, {{}, _}],
			Join[#, set[[{-1}]]]&/@ DerivativeSubLists@ First@ set
		,
			{{Last@ set}}
		]
	, {set, sets}]
];


(* ::Text:: *)
(*Partitions out a subset of the list *)


PartDrop[l_List, i_]:= {l[[i]], l[[Complement[Range@ Length@ l, i]]]};


(* ::Subsubsection::Closed:: *)
(*Generic field strength tensors*)


(* ::Text:: *)
(*Properties*)


FSWilsonFactor[{a_, a_}, _, _]= 0;
FSWilsonFactor[_, {___, {(Index[a_, rep_]|Bar@ Index[a_, rep_])..}, ___}, _]= 0;


(* ::Text:: *)
(*Expand out all generic field strength tensors in an expression*)


ExpandGenFSs@ expr_:= expr/. FSWilsonFactor:> DevTermOnWilson// ContractDelta// SeparateGeneratorsFromFS


(* ::Text:: *)
(*A particular index set in the formula for the derivative set on a Wilson line *)


DevTermOnWilson[devInds_, gaugeInds_, charges_]:= Block[
		{abelFS, nonAbelFS, charge, group, n, perm, permutations, prefact},
	permutations= Permutations@ devInds[[;; -2]];

	(*Abelian*)
	abelFS= Times@@ Delta@@@ gaugeInds* Sum[
			First@ charge* (*$GaugeGroups[Head@ charge, Coupling][]* *)
			Sum[
				(*FieldStrength[$GaugeGroups[Head@ charge, Field], devInds[[-2;;]], {}, devInds[[;; -3]]]*)
				FieldStrength[$GaugeGroups[Head@ charge, Field], {perm[[-1]], devInds[[-1]]},
					{}, perm[[;; -2]]]
			, {perm, permutations}]
		, {charge, charges}];

	(*Non-Abelian FS*)
	nonAbelFS= Sum[
			group= GroupFromInd@ gaugeInds[[n, 1]];
			(*$GaugeGroups[group, Coupling][]**) Times@@ Delta@@@ Delete[gaugeInds, n]*
			Sum[
				FieldStrength[$GaugeGroups[group, Field], {perm[[-1]], devInds[[-1]]},
					gaugeInds[[n]], perm[[;; -2]]]
			, {perm, permutations}]
		, {n, Length@ gaugeInds}];

	(*Sum over all permutations of the of all but the last index of indSet*)
	(abelFS+ nonAbelFS)/ Length@ permutations
]


PermutationReplacements@ list_:= Thread@ Rule[list, #]&/@ Permutations@ list;


(* ::Text:: *)
(*Creates several sets of gauge indices that contract to the original set, e.g., {{i, j}} -> { {{i, i[1]}}, {{i[1], j}} }*)


GaugeIndexSet[originalSet_, 1]:= {originalSet};
GaugeIndexSet[originalSet_, multiples_Integer]:= Block[{lab, out},
	If[Length@ originalSet === 0, Return@ ConstantArray[{}, multiples];];
	lab= FirstCase[originalSet[[1, 1]], Index[ind_, _]:> ind, {}, All];
	out= ConstantArray[originalSet, multiples];
	out[[;;-2, ;;, 2]]= MapIndexed[#1/. Index[_, rep_]:> Index[lab@ First@ #2, rep] &, out[[;;-2, ;;, 2]] ];
	out[[2;;, ;;, 1]]= MapIndexed[#1/. Index[_, rep_]:> Index[lab@ First@ #2, rep] &, out[[2;;, ;;, 1]] ];
	out
]


(* ::Section:: *)
(*General properties*)


(* ::Subsubsection::Closed:: *)
(*Mass term properties*)


(* ::Text:: *)
(*Masses of light field types vanish*)


Mop[lScalar|lFermion|lVector|lGhost|lAntiGhost, __]:= 0;


(* ::Text:: *)
(*Heavy masses are assumed to be real*)


Mterm[Conj@ lab_, i_]:= Mterm[lab, i];


(* ::Subsubsection::Closed:: *)
(*Ancillary functions*)


(* ::Text:: *)
(*Partition a list into lists of variabel length  *)


VariableLengthPartition[list_, lengths_]:=
	list[[#1 +1;; #2]]&@@@ Partition[Accumulate@Prepend[lengths, 0], 2, 1];
