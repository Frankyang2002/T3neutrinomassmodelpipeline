(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`CouplingSimplifications`*)


(* ::Subtitle:: *)
(*Paclet for simplifying coupling and coefficient expressions*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["SimplifyCouplings"]


PackageExport["DefineCoSinePair"]


(* ::Subsubsection::Closed:: *)
(*Scoped*)


PackageScope["CollectCoefficients"]
PackageScope["CoefficientSimplify"]
PackageScope["ContractOverallDelta"]


PackageScope["CouplingContractionPattern"]
PackageScope["$nonTrivCouplingPattern"]


PackageScope["TrigSimplify"]
PackageScope["ResetTrigonometrics"]
PackageScope["$trigPairs"]
PackageScope["$unitaryMatrixSymbols"]


PackageScope["CSSimplify"]


PackageScope["DefineMixingMatrix"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


SimplifyCouplings::usage= "SimplifyCouplings[expr] simplify coupling contractions with dummy indices.";


DefineCoSinePair::usage= "DefineCoSinePair[<cos>, <sin>] defines a set of couplings to behave as though they are cosine and sine functions of the same angle."


(* ::Subsubsection::Closed:: *)
(*Scoped*)


DefineMixingMatrix::usage= "DefineMixingMatrix[<string>, <size>, <type>] defines a matrix of symbols behaving as the entries of a matrix of the requested type under trigonometric simplifications."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Coupling collection*)


(* ::Text:: *)
(*Provides methods to collect groups of couplings accounting for naming of repeated indices and coupling symmetries*)


(* ::Subsubsection::Closed:: *)
(*Basic coupling patterns*)


(* ::Text:: *)
(*Basic pattern to check if coefficient object  contains indices.*)
(*This must be loaded first to be used in the following functions*)


$nonTrivCouplingPattern=
	Coupling[_, {__}, _]| Bar@ Coupling[_, {__}, _]|
	Power[Coupling[_, {_}, _], _]| Power[Bar@ Coupling[_, {_}, _], _]|
	_Delta| Bar[_Delta]| Log[Times[_, Power[Coupling[_, {_}, _], _]]]|
	LF[{___, Coupling[_, {_}, _], ___}, _];


(* ::Subsection:: *)
(*Collect couplings with operator symmetries*)


(* ::Text:: *)
(*This version of coefficient collecting accounts for the flavor symmetries of the operator*)


(* ::Subsubsection::Closed:: *)
(*Collect coefficient of operators*)


(* ::Text:: *)
(*Main function to apply to a Lagrangian on internal form (involving AtomicOp, and CompOp)*)


Options[CollectCoefficients] = {
		CoefficientSimplify-> True,
		ContractOverallDelta -> True, 
		Simplify-> False,
		TrigSimplify-> False
	}


CollectCoefficients[expr_, OptionsPattern[]]:= Module[{out, simplificationFunc, trigVars= Join@@ $trigPairs},
	out=  Collect[ContractUnitaryCouplings@ CanonizeLoopFunctions@ expr, {_AtomicOp|_CompOp}];
	(*Should be updated when we implement flavor symmetries in the simplification module*)
	out= out/. {
		Times[coef___, op:AtomicOp[class_, id_, inds_]]:> op*
			CanonizeCouplingContractions[Times@ coef, inds, $operators[class, id, Symmetries]],
		Times[coef___, op:CompOp[red_, class_, id_, inds_]]:> op*
			CanonizeCouplingContractions[Times@ coef, inds, $compoundOperators[red, class, id, Symmetries]]
	};

	(*Determine what (if any) simplification function to apply to the coefficient*)
	simplificationFunc= Which[
		OptionValue@ CoefficientSimplify && OptionValue@ TrigSimplify,
			Simplify[Collect[#, c_Coupling/; !MemberQ[trigVars, c], TrigSimplify], 
				TimeConstraint-> 0.1]&
		,OptionValue@ CoefficientSimplify,
			Simplify[#, TimeConstraint-> 0.1]&
		,OptionValue@ TrigSimplify,
			CollectCoefficients[#, c_Coupling/; !MemberQ[trigVars, c], TrigSimplify]&
		,True, 
			Identity
		];
	
	Quiet[
		out= Collect[out, {_AtomicOp|_CompOp, hbar, \[Epsilon]},
			Collect[RelabelIndices@ #, {Power[c_Coupling, _?Negative]/; !MemberQ[trigVars, c], _LF, _Log}, simplificationFunc]& ]; 
	, Simplify::time];

	If[OptionValue[ContractOverallDelta],
		(* convert output to list without expanding it -> do NOT use TermsToList here! *)
		out= If[Head@ out === Plus, List@@ out, List@ out];

		(* contract back overall deltas *)
		out= Plus@@( ContractDeltaSingleTerm[#,"overall-only"->True]&/@ out )
	];
	
	If[OptionValue@ Simplify,
		out= Quiet[Collect[out, {_AtomicOp|_CompOp, hbar, \[Epsilon]}, Simplify[#,TimeConstraint->1]&]
			,Simplify::time];
	];
	
	out
]


(* ::Subsubsection::Closed:: *)
(*Contract unitary couplings*)


Options@ ContractUnitaryCouplings= {Expand-> True} (* Expand-> False might miss some possible simplifications due to unitarity relations, but still gives correct results *)


ContractUnitaryCouplings[expr_, OptionsPattern[]]:= Module[
	{
		out,
		(* find complex and real unitary couplings *)
		unitaryComplex= GetCouplingByProperty[Unitary-> True, SelfConjugate-> False],
		unitaryReal=    GetCouplingByProperty[Unitary-> True, SelfConjugate-> True]
	}
	,
	If[Length[unitaryComplex]+Length[unitaryReal]===0,
		(* no unitary couplings: do nothing *)
		out= expr
		,
		(* unitary couplings present: expand output (by default) *)
		out= TermsToList[expr, Expand-> OptionValue@ Expand];
		(* contract unitary couplings in each term individually *)
		out= ContractUnitaryCouplingsSingleTerm[#, unitaryComplex, unitaryReal]&/@ out;
		(* sum results *)
		out= Plus@@ out
	];
	
	out
]


ContractUnitaryCouplingsSingleTerm[expr_, unitaryComplex_, unitaryReal_]:= Module[{mark, out= expr, subs, inds, diagInds},
	(* find diagonal indices in a term *)
	(* NOTICE: do not use DiagonalInidicesSingleTerm here *)
	out= PseudoTimes@ out;
	(* Check which inidces appear a different number of times before/after applying $DropDiagonalCouplings *)
	diagInds= Plus@@Cases[out, _Index, All] - Plus@@Cases[out/. $DropDiagonalCouplings, _Index, All];
	diagInds= Cases[diagInds, _Index, All];
	
	(* reset out *)
	out= expr;
	
	(* create replacement patterns; these should not trigger, if the contracted index is shared with a diagonal coupling *)
	subs= Switch[{Length@ unitaryComplex, Length@ unitaryReal}
		,{0, 0},
			Return@ out;
		,{_, 0},
			{
			Coupling[lab_, {p_, r_}, 0] Bar@ Coupling[lab_, {p_, s_}, 0]/; (MemberQ[unitaryComplex, lab]&& !MemberQ[diagInds, p]):> mark* Delta[r, s],
			Coupling[lab_, {p_, r_}, 0] Bar@ Coupling[lab_, {s_, r_}, 0]/; (MemberQ[unitaryComplex, lab]&& !MemberQ[diagInds, r]):> mark* Delta[p, s]
			}
		,{0, _},
			{
			Coupling[lab_, {p_, r_}, 0] Coupling[lab_, {p_, s_}, 0]/; (MemberQ[unitaryReal, lab]&& !MemberQ[diagInds, p]):> mark* Delta[r, s],
			Coupling[lab_, {p_, r_}, 0] Coupling[lab_, {s_, r_}, 0]/; (MemberQ[unitaryReal, lab]&& !MemberQ[diagInds, r]):> mark* Delta[p, s]
			}
		,{_, _},
			{
			Coupling[lab_, {p_, r_}, 0] Bar@ Coupling[lab_, {p_, s_}, 0]/; (MemberQ[unitaryComplex, lab]&& !MemberQ[diagInds, p]):> mark* Delta[r, s],
			Coupling[lab_, {p_, r_}, 0] Bar@ Coupling[lab_, {s_, r_}, 0]/; (MemberQ[unitaryComplex, lab]&& !MemberQ[diagInds, r]):> mark* Delta[p, s],
			Coupling[lab_, {p_, r_}, 0] Coupling[lab_, {p_, s_}, 0]/; (MemberQ[unitaryReal, lab]&& !MemberQ[diagInds, p]):> mark* Delta[r, s],
			Coupling[lab_, {p_, r_}, 0] Coupling[lab_, {s_, r_}, 0]/; (MemberQ[unitaryReal, lab]&& !MemberQ[diagInds, r]):> mark* Delta[p, s]
			}
	];
	
	(* apply unitarity relations *)
	out= out/. subs; 
	
	(* contract resulting deltas and re-apply unitarity relations until no simplifications happen *)
	While[!FreeQ[out, mark],
		out= ContractDelta@ out/. mark-> 1/. subs 
	];
	
	out
];


(* ::Subsubsection::Closed:: *)
(*Simplify coefficient given symmetries of open flavor indices*)


(* ::Text:: *)
(*Simplifies a coefficient given symmetries of the operator it is contracted with*)


CanonizeCouplingContractions[coef_, opInds_List, syms_List]:= Module[
		{aSymMatches, contractions, couplingTerms, couplingContractions, newPats, out, rules, pos, symMatches, symmetryPermutations},
	out= PseudoTimes@ BetterExpand@ coef;
	(*Remove everything but valid coupling patterns from PseudoTimes*)
	out= out//.
		PseudoTimes[a___, b: Except[$nonTrivCouplingPattern], c___]:> b PseudoTimes[a, c]/.
		PseudoTimes[]->1;

	(*All non-trivial coupling contractions*)
	couplingContractions= Cases[{out}, _PseudoTimes, All];

	(*Construct matching patterns, while gathering the coupling contractions that matches these patterns*)
	symmetryPermutations= GroupBy[syms, Last, Map[First]];
	{rules, contractions}= {{}, {}};
	While[Length@ couplingContractions> 0,
		AppendTo[rules, CouplingContractionPatternWithSyms[First@ couplingContractions, opInds, symmetryPermutations]];
		(*rules[[-1]] is now {rep. rules w/ + signature, rep. rules w/ - signature (or None)} *)
		pos= Position[couplingContractions,
			If[rules[[-1, 2]] === None, rules[[-1, 1, 1]], Alternatives@@ rules[[-1, ;;, 1]]] ];
		AppendTo[contractions, Extract[couplingContractions, pos]];
		couplingContractions= Delete[couplingContractions, pos];
	];
	
	(*For each pattern w/ associated coupling contractions, match patterns of both signatures in all possible ways.
		This ensures proper cancellation with e.g. symmetric coupling and antisymmetric open indices.*)
	(*The result is replacement rules for each coupling contraction in coef*)
	rules= Flatten@ Table[
		Thread@ Rule[contractions[[n]],
			Table[
				symMatches= ReplaceList[contraction, rules[[n, 1]]];
				aSymMatches= If[rules[[n, 2]] === None, {}, ReplaceList[contraction, rules[[n, 2]]]];
				(Total@ symMatches- Total@ aSymMatches)/ (Length@ symMatches+ Length@ aSymMatches)// Expand
		, {contraction, contractions[[n]]}] ]
	, {n, Length@ rules}];

	out/. rules// ReleasePseudoTimes
]


(* ::Subsubsection::Closed:: *)
(*Construct coupling contraction patterns*)


(* ::Text:: *)
(*Creates matching patterns out of a coupling contraction given symmetries of the open indices. *)
(*Return:*)
(*	{pattern that matches with + signature, pattern that matches with - signature/None (if not applicable)}*)


CouplingContractionPatternWithSyms[coups:PseudoTimes[$nonTrivCouplingPattern..], opInds_List, syms_Association]:= Module[
		{openInds, perm, repInds, out, couplingRules, openInPattern, openPerms, pos, repl},
	(*Replace repeated indicies with patterns*)
	repInds= Cases[coups/. $DropDiagonalCouplings, Index[__], All];
	(*If the coefficient does not fully contract the operator (non-Lagrangian), symmetries are not applicable!*)
	If[!SubsetQ[repInds, opInds],
		Return@ {CouplingContractionPattern@ coups, None};
	];
	repInds= DeleteCases[Tally@ repInds, {i_, 1}][[;;, 1]];

	{repInds, openInds}= Map[With[{temp= Unique@ First@ #}, Rule[#, Pattern[temp, Blank[]] ] ]&,
		{repInds, opInds}, {2}];

	openPerms= Map[opInds[[#]]&, syms, {2}];
	openInPattern= openInds[[;;, 2]]/. Verbatim[Pattern][s_, Blank[]]-> s;

	(*Make coupling pattern*)
	couplingRules= coups/. repInds/. openInds;
	couplingRules= CouplingPattern/@ List@@ couplingRules;
	
	(*Identify delta indices in diagonal couplings, e.g., \[Delta]_ab M_b = \[Delta]_ab M_a*)
	pos= Position[couplingRules, _Delta];
	If[Length@ pos> 0,
		repl= Extract[couplingRules, pos]/. Delta[args__]:> 
			(Verbatim/@ Alternatives[args]-> Alternatives[args]);
		pos= Delete[Range@ Length@ couplingRules, pos[[;;, {1}]]];
		couplingRules[[pos]]= couplingRules[[pos]]/. repl;
	];
	
	(*Build a signed replacement rule*)
	AppendTo[couplingRules, {1, coups}];
	couplingRules= RuleDelayed@@ MapThread[PseudoTimes, couplingRules]// Activate;

	{ (*Pattern matching '+' signature permutation of open inds*)
		With[{oi= openInPattern, op= openPerms@ 1},
			MapAt[Condition[#, MemberQ[op, oi]]&, couplingRules, 1]
		]
	, (*Pattern matching '-' signature permutation of open inds*)
		If[KeyExistsQ[openPerms, -1],
			With[{oi= openInPattern, op= openPerms[-1]},
				MapAt[Condition[#, MemberQ[op, oi]]&, couplingRules, 1]
			]
		,
			None
		]
	}
]


(* ::Text:: *)
(*Makes a replacement pattern corresponding to a given combination of the couplings. *)


CouplingContractionPattern[coups:PseudoTimes[$nonTrivCouplingPattern..]]:= Module[
		{indices, out, couplingRules},
	(*Replace repeated indicies with patterns*)
	indices= DeleteCases[Tally@ Cases[coups, Index[__], All], {i_, 1}][[;;, 1]];
	indices= With[{temp= Unique@ First@ #}, Rule[#, Pattern[temp, Blank[]]] ]&/@ indices;
	out= coups/. indices;

	(*Build a signed replacement rule*)
	couplingRules= CouplingPattern/@ List@@ out;
	AppendTo[couplingRules, {1, coups}];
	RuleDelayed@@ MapThread[PseudoTimes, couplingRules]// Activate
]


(* ::Subsubsection::Closed:: *)
(*Constructs patterns out of coefficient objects*)


(* ::Text:: *)
(*Constructs the pieces for a substitution rule with all symmetry permutations of the indices  of a coupling pattern given the patterns. *)
(*CouplingPattern returns {LHS, RHS}.*)


CouplingPattern@ c:Coupling[lab_, indsPat_, ord_]:= Module[{symsPat, inds, syms, x},
	syms= $CouplingAssociation[lab, Symmetries];
	If[Length@ syms === 1,
		Return[{c, 1}];
	];

	symsPat= GroupBy[KeyMap[indsPat[[#]]&, syms], Identity, Keys];
	If[KeyExistsQ[symsPat, -1],
		{Coupling[lab, (x: Alternatives@@ Lookup[symsPat, 1]) | Alternatives@@ Lookup[symsPat, -1], ord],
			Inactive[DidMatchSwitch][1, -1][x]}
	,
		{Coupling[lab, Alternatives@@ Lookup[symsPat, 1], ord], 1}
	]
]


CouplingPattern[c:Bar@ Coupling[lab_, indsPat_, ord_]]:= Module[{symsPat, inds, syms, x},
	syms= $CouplingAssociation[lab, Symmetries];
	If[Length@ syms === 1,
		Return[{c, 1}];
	];

	symsPat= GroupBy[KeyMap[indsPat[[#]]&, syms], Identity, Keys];
	If[KeyExistsQ[symsPat, -1],
		{Bar@ Coupling[lab, (x: Alternatives@@ Lookup[symsPat, 1]) | Alternatives@@ Lookup[symsPat, -1], ord],
			Inactive[DidMatchSwitch][1, -1][x]}
	,
		{Bar@ Coupling[lab, Alternatives@@ Lookup[symsPat, 1], ord], 1}
	]
]


CouplingPattern@ pow:Power[_Coupling, _]:= {pow, 1};
CouplingPattern@ log_Log:= {log, 1};


(*Account for the symmetries of the Loop functions*)
CouplingPattern@ LF[masses_, subscripts_]:= Module[{perms, n},
	perms= Transpose@ {masses, Most@ subscripts, Table[n, {n, Length@ masses}]};
	perms= GroupBy[perms, {#[[1, 1]], #[[2]]} &];
	perms= KeyValueMap[If[GetCouplings[#1, DiagonalCoupling] === {}, {#2[[;;, -1]]}, Permutations@ #2[[;;, -1]]]&,
		perms];

	With[{lhs= Alternatives@@ Table[LF[masses[[n]], Append[subscripts[[n]], subscripts[[-1]]]], {n, Join@@@ Tuples@ perms}]},
		{lhs, 1}]
];


(*Put deltas with flavor indices on the same form*)
CouplingPattern@ Delta[indsPat__]:= {Delta@ indsPat, 1};
CouplingPattern@ Bar@Delta[indsPat__]:= {Bar@ Delta@ indsPat, 1};


(* ::Subsection:: *)
(*Collect couplings generic*)


(* ::Subsubsection::Closed:: *)
(*SimplifyCouplings*)


(* ::Text:: *)
(*Collect like coupling contractions in a  generic expression *)


SimplifyCouplings@ l_List:= SimplifyCouplings/@ l;


SimplifyCouplings@ expr_:= Module[{out= PseudoTimes@ Expand@ expr, couplingContractions, pats= {}},
	(*If only one term, nothing to collect*)
	If[Head@ out =!= Plus, Return@ ReleasePseudoTimes@ out];
	(*All non-trivial coupling contractions*)
	couplingContractions= Cases[out, _PseudoTimes, All];
	couplingContractions= PseudoTimes@@@ DeleteCases[Cases[#,
		$nonTrivCouplingPattern]&/@ couplingContractions, {}];

	(*Construct patterns*)
	(*Optimization may be required*)
	While[Length@ couplingContractions> 0,
		AppendTo[pats, CouplingContractionPattern@ First@ couplingContractions];
		couplingContractions= DeleteCases[couplingContractions, pats[[-1, 1]]];
	];
	(*Remove everything but the coupling patterns from PseudoTimes*)
	out= out//. PseudoTimes[a___, b: Except[$nonTrivCouplingPattern]]:>
		b PseudoTimes@ a;
	out/. pats// ReleasePseudoTimes
]


(* ::Section:: *)
(*Trigonometric and orthogonal simplifications*)


(* ::Subsection:: *)
(*Simplify trigonometric expressions*)


(* ::Subsubsection::Closed:: *)
(*TrigonometricSimplify*)


ResetTrigonometrics[]:= Block[{}, 
	$trigPairs= {};
	$unitaryMatrixSymbols= <||>;
];
ResetTrigonometrics[];


(* ::Text:: *)
(*Initialize a cosine-sine pair of couplings*)


DefineCoSinePair[cos_Symbol, sin_Symbol]:= Module[{c, s},
	DefineCoupling[cos, EFTOrder-> 0, SelfConjugate-> True];
	DefineCoupling[sin, EFTOrder-> 0, SelfConjugate-> True];
	{c, s}= {cos[], sin[]};
	AppendTo[$trigPairs, {c, s}];
];


(* ::Text:: *)
(*Initialize a matrix of couplings, satisfying SO/O/SU/U(N) constraints*)


DefineMixingMatrix::size= "A mixing matrix should be of size <= 5.";
DefineMixingMatrix::type= "The input mixing matrix should a square matrix of symbols.";


DefineMixingMatrix[name_String, size_Integer, type_Symbol]:= Module[{matrix},
	If[size > 5 || size < 2,
		Message[DefineMixingMatrix::size];
		Abort[];
	];
	matrix= Table[
			Symbol[name <> ToString@ m <> ToString@ n]
		, {m, size}, {n, size}];
	DefineMixingMatrix[matrix, type];
];


DefineMixingMatrix[symbolMatrix_, type_Symbol]:= Module[{basis, relations, matrix, size, symbols},
	If[!MatchQ[symbolMatrix, {{__Symbol}..}] || !SquareMatrixQ@ symbolMatrix,
		Message[DefineMixingMatrix::type];
		Abort[];
	];
	size= Length@ symbolMatrix;
	If[size > 5 || size < 2,
		Message[DefineMixingMatrix::size];
		Abort[];
	];
	
	DefineCoupling[#, SelfConjugate-> True, NiceForm-> CanonicalNiceForm@ #]&/@ Flatten@ symbolMatrix;
	matrix= Map[#[]&, symbolMatrix, {2}];
	
	Switch[type
	,O|SO,
		relations = Flatten@ Table[
				{matrix[[m]] . matrix[[n]] - KroneckerDelta[m, n],
				matrix[[All, m]] . matrix[[All, n]] - KroneckerDelta[m, n]}
			, {n, size}, {m, n, size}];
	,_,
		Abort[];	
	];
	
	If[type === SO,
		AppendTo[relations, Det@ matrix - 1]
	]; 
	
	symbols= Flatten@ matrix;
	basis= GroebnerBasis[relations, symbols, MonomialOrder-> DegreeReverseLexicographic];
	
	$unitaryMatrixSymbols[symbols]= basis;
];


CanonicalNiceForm@ symbol_Symbol:= Module[{str= ToString@ symbol, index, base, split},
	split= StringPosition[str, RegularExpression["\\D\\d*$"]][[1, 1]];
	{base, index}= StringTakeDrop[str, split];
	"\!\(\*SubscriptBox[\("<> base<> "\), \("<> index<> "\)]\)"
];


(* ::Text:: *)
(*Simplification with cosine-sine pair*)


TrigSimplify@ expr_:= Module[{cspattern, varPattern, vars, out= expr},
	(*Simplify with all co-sine pairs*)
	Do[
		cspattern= Alternatives@@ csPair;
		If[!FreeQ[out, cspattern],
			out= Collect[out, Except[cspattern, _Coupling|Bar@ _Coupling], CSSimplify[#, csPair]&];
		];
	, {csPair, $trigPairs}];
	(*Simplify using orthogonality/unitarity relations*)
	Do[
		varPattern= Alternatives@@ vars;
		If[!FreeQ[out, varPattern],
			out= Collect[out, Except[varPattern, _Coupling|Bar@ _Coupling], 
				UnitaryMatrixSimplify[#, vars, $unitaryMatrixSymbols@ vars]&];
		];
	, {vars, Keys@ $unitaryMatrixSymbols}];
	out
]


UnitaryMatrixSimplify[expr_, vars_List, basis_List]:= Module[{},
	Last@ PolynomialReduce[expr, basis, vars, MonomialOrder-> DegreeReverseLexicographic]
]


(* ::Subsubsection::Closed:: *)
(*CSSimplify*)


(* ::Text:: *)
(*Heuristic simplification that's applicable to rational functions of cos and sine (c, s)*)


CSSimplify[expr_, {c_, s_}]:= Module[{num, den, out= Together@ expr},
	{num, den}= {Numerator@ out, Denominator@ out};
	num= CSPolySimplify[num, {c, s}];
	den= CSPolySimplify[den, {c, s}];
	num/ den// Simplify
]


(* ::Text:: *)
(*Simplification function for a polynomial in cos and sine*)


CSPolySimplify[poly_, {c_, s_}]:= Module[
		{case1, case2, candidate, coefs, out, seq, sequences, subSequence, temp, n, x, y},
	out= Expand@ poly/. {
			Power[s, n_]:> If[OddQ@ n, s Power[x, (n- 1)/2], Power[x, n/2]],
			Power[c, n_]:> If[OddQ@ n, c Power[y, (n- 1)/2], Power[y, n/2]] 
		};
	
	(*Sine first simplification *)
	candidate= out/. y-> 1- x// Expand;
	candidate= FixedPoint[IterativeCoSineReduction[x, y], candidate];
	
	(*Try cosine first-simplification*)
	case1= Expand@ candidate;
	If[Head@ case1 === Plus, 
		temp= out/. x-> 1- y// Expand;
		temp= FixedPoint[IterativeCoSineReduction[y, x], temp];
		case2= Expand@ temp;
		If[Head@ case2 =!= Plus || Length@ case2 < Length@ case1, candidate= temp];
	];
		
	candidate/. {x-> s^2, y-> c^2}
]


(* ::Text:: *)
(*Attempts heuristic simplification steps on subsequences of the sin^2 series*)


$maxCoSineSearchOrd= 5; (*Sets the limit for the number of terms to include in brute force simplification*)
IterativeCoSineReduction[x_, y_][poly_]:= Module[
		{coefs, seq, sequences, subSequence, n},
	coefs= CoefficientList[poly, x];
	sequences= SequencePosition[coefs, {Except[0]..}, Overlaps-> False];
	(*Sum over subsequences --- 
	assumption is that one can rarely simplify with between terms separated by s^4 or more*)
	Sum[
		subSequence= coefs[[Span@@ seq]];
		If[seq[[2]]- seq[[1]] >= $maxCoSineSearchOrd, (*Too long to bruteforce*)
			Sum[subSequence[[n]] Power[x, n-1], {n, Length@ subSequence}]
		,
			BruteCoSineReduction[subSequence, x, y]
		] * Power[x, First@ seq- 1]
	, {seq, sequences}]
]


(* ::Text:: *)
(*Attempts substitutions of s^2 -> 1- c^2 at each step in the series (but not e.g. s^4 -> s^2 (1- c^2) to improve performance) *)


BruteCoSineReduction[coefs_List, x_, y_]:= Module[
		{best, minLen, aList, bList, perm, perms, n, score, temp, len= Length@ coefs},
	If[len === 1, First@ coefs];
	
	aList= Table[Power[x, n], {n, 1, len-1}];
	bList= Expand@ Table[Power[1-y, n], {n, 1, len-1}];
	perms= Tuples@ Transpose@ {aList, bList};
	
	minLen= $maxCoSineSearchOrd + 10; 
	Catch[
		Do[
			temp= coefs[[1]] + coefs[[2;;]] . perm// Expand;
			If[Head@ temp =!= Plus, Throw@ temp]; (*Single term*)
			score= Length@ temp;
			If[score < minLen,
				{minLen, best}= {score, temp};
			];
		, {perm, perms}];
		best
	]
]
