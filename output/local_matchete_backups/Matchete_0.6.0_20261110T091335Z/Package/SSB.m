(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`SSB`*)


(* ::Subtitle:: *)
(*Procedures for spontaneous symmetry breaking*)


(* ::Text:: *)
(*We mainly follow the general notation used in https://arxiv.org/abs/2404.11640*)


(* ::Chapter:: *)
(*Public:*)


(* ::Subsection:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["SetSymmetryBreakingPattern"]
PackageExport["RepresentationDecomposition"]
PackageExport["FieldDecomposition"]
PackageExport["CGDecomposition"]
PackageExport["ToBrokenPhase"]
PackageExport["ImplementVacuumConditions"]
PackageExport["SubstituteMasses"]
PackageExport["GetVacuumConditions"]


PackageExport["SetSSBReplacements"]


PackageExport["Singlet"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["$SSBInfo"]
PackageScope["ResetSymmetryBreaking"]


PackageScope["BrokenPhaseLagQ"]
PackageScope["SubstituteMassesInLagrangian"]
PackageScope["HeavyVectorLabels"]


PackageScope["UnitaryGauge"]
PackageScope["BrokenPhaseRXiTerms"]


(* ::Subsection:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


SetSymmetryBreakingPattern::usage= "SetSymmetryBreakingPattern[<group>, <stability group>] sets the current symmetry breaking pattern.";
RepresentationDecomposition::usage= "RepresentationDecomposition[<original rep>, <decomposition>] dictates how each original representation decomposes under the stability group."
FieldDecomposition::usage= "FieldDecomposition[<field>, <decomposition>] specifies how a field decomposes under the stability group.";
CGDecomposition::usage= "CGDecomposition[<CG>, <decomposition>] specifies how a CG of the original group decomposes into CGs of the stability group.";
ToBrokenPhase::usage= "ToBrokenPhase[<Lagrangian>] takes the Lagrangian to the broken phase as dictated by the set symmetry breaking pattern."
ImplementVacuumConditions::usage= "ImplementVacuumConditions[<Lagrangian>] simplifies a broken phase Lagrangian, with assumptions based on the vacuum being true and the broken phase fields being properly identified as eigenstates."
SubstituteMasses::usage= "SubstituteMasses is an option that controls whether the mass terms (and vector-Goldstone kinetic mixing terms) are written in terms of the field mass symbols rather than the underlying couplings."
GetVacuumConditions::usage= "GetVacuumConditions[<Lagrangian>] returns an association containing all vacuum conditions that would be applied by ImplementVacuumConditions on the same Lagrangian.";
SetSSBReplacements::usage= "SetSSBReplacements[rules] lets one fix the coupling replacement rules that are applied to the Lagrangian at the beginning of the routine ToBrokenPhase.";


Singlet::usage= "Used to indicates a singlet representation in the decomposition of broken symmetry representations."


(* ::Subsubsection::Closed:: *)
(*Internal*)


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Setting up the breaking pattern*)


(* ::Subsubsection::Closed:: *)
(*Index projectors*)


(* ::Text:: *)
(*The BrokenInds facilitate contraction between  the right broken components from various fields and CGs *)


ResetBrokenInds[]:= Module[{},
	Clear@ BrokenInds;
	BrokenInds/: Bar@ b_BrokenInds:= b;
	BrokenInds[Bar@ ind_, m_]:= BrokenInds[ind, m];
]


(* ::Subsection:: *)
(*Breaking pattern*)


(* ::Subsubsection::Closed:: *)
(*Groups*)


(* ::Text:: *)
(*Defines the current symmetry-breaking structure pattern.*)
(*	- Consider adding a test if the stability group is a subgroup of the group*)


ResetSymmetryBreaking[]:= Block[{}, 
	$SSBInfo= <|
			Charges-> <||>,
			CGSubstitutions-> <||>,
			FieldSubstitutions-> <||>,
			OriginalGroup-> {},
			Representations-> <||>,
			StabilityGroup-> {}
		|>;
	ResetBrokenInds[];
	$SSBCouplingReplacements= {};
];
ResetSymmetryBreaking[];


SetSymmetryBreakingPattern::notgroup= "The `1` received invalid argument `2`. It must be a known group or list of groups.";
SetSymmetryBreakingPattern::overlap= "The group(s) `1` appear(s) multiple times among the group breaking pattern.";


SetSymmetryBreakingPattern[group_, stabilityGroup_]:= Module[{allGroups, orgGroup, stabGroup},
	orgGroup= If[Head@ group =!= List, {group}, group];
	stabGroup= If[Head@ stabilityGroup =!= List, {stabilityGroup}, stabilityGroup];
	
	If[!MatchQ[orgGroup, {Alternatives@@ Join[Keys@ $GaugeGroups, Keys@ $GlobalGroups]..}],
		Message[SetSymmetryBreakingPattern::notgroup, "Group", group];
		Abort[];	
	];
	
	If[!MatchQ[stabGroup, {Alternatives@@ Join[Keys@ $GaugeGroups, Keys@ $GlobalGroups]...}],
		Message[SetSymmetryBreakingPattern::notgroup, "stability group", stabilityGroup];
		Abort[];	
	];
	
	allGroups= Join[orgGroup, stabGroup];
	If[!DuplicateFreeQ@ allGroups,
		allGroups= DeleteCases[Tally@ allGroups, {_, 1}][[;;, 1]];
		Message[SetSymmetryBreakingPattern::overlap, allGroups];
		Abort[];
	];
	
	(*Check inputs are list of defined gauge groups*)
	ResetBrokenInds[];
	$SSBInfo= <|
			Charges-> <||>,
			CGSubstitutions-> <||>,
			FieldSubstitutions-> <||>,
			OriginalGroup-> orgGroup,
			Representations-> <||>,
			StabilityGroup-> stabGroup
		|>;	
];


(* ::Subsubsection::Closed:: *)
(*Representations*)


(* ::Text:: *)
(*Defines the decomposition of a representation from the original group in terms of representations of the stability groups*)


RepresentationDecomposition::invldorgrep= "`1` is not a recognized representation of any of the groups of the original symmetry, nor is it the unit charge from an abelian component.";


RepresentationDecomposition[rep_, decomp_]:= Module[{abelians},
	If[MemberQ[Keys@ $Representations, rep],
		If[!MemberQ[$SSBInfo@ OriginalGroup, $Representations[rep, Group]],
			Message[RepresentationDecomposition::invldorgrep, rep];
			Abort[];
		];
		NonAbelianRepresentationDecomposition[rep, decomp];
	,
		abelians= Keys@ Select[$GaugeGroups, #@ Group === U1 &];
		If[!MatchQ[rep, (Alternatives@@ abelians)[1]] || !MemberQ[$SSBInfo@ OriginalGroup, Head@ rep],
			Message[RepresentationDecomposition::invldorgrep, rep];
			Abort[];
		];
		AbelianRepresentationDecomposition[rep, decomp];
	];
];


RepresentationDecomposition::invldcharge= "The unit charge `1` of the original group must decompose as a rational charge of an Abelian component of the stability (remnant) group or Nothing.";


(* ::Text:: *)
(*Specify how a unit charge of the original Abelian factor decomposes in terms of the stability group charges.*)


AbelianRepresentationDecomposition[gr_[1], decomp_]:= Module[
	{c, decomposition, q, replacements, validCharges},
	If[decomp === Nothing,
		$SSBInfo[Charges, gr]= gr[_]-> Nothing;
		Return[];
	];
	
	(*Ensure we have a list of charges*)
	decomposition= Flatten@ {decomp};
	
	validCharges= Intersection[Keys@ Select[$GaugeGroups, #@ Group === U1 &], 
		$SSBInfo@ StabilityGroup];
	replacements= Table[
			If[!MatchQ[c, (Alternatives@@ validCharges)[_Integer|_Rational]],
				Message[RepresentationDecomposition::invldcharge, c];
				Abort[];
			];
			With[{lhs= gr[q_], rhs= (q # &)/@ c},
				{lhs-> rhs}]
		, {c, decomposition}];
	
	$SSBInfo[Charges, gr]= replacements;
];


RepresentationDecomposition::invldsumrepord= "The embedding of the complex representation '`1`' must either precede or follow its conjugate in the decomposition of the real representation '`2`'.";
RepresentationDecomposition::incmpdims= "The dimension `1` of the decomposition does not agree with the dimension `2` of the original representation.";
RepresentationDecomposition::complexunfollowed= "The embedding of the complex representation '`1`' cannot be unpaired in the decomposition of the real representation '`2`'.";
RepresentationDecomposition::realfollowcomplex= "The embedding of the real representation '`1`' directly follows the unpaired complex representation '`2`', which is not its conjugate, in the decomposition of the real representation '`3`'.";
RepresentationDecomposition::repord= "Make sure that any complex representation in the embedding is immediately followed by its conjugate.";


NonAbelianRepresentationDecomposition[rep_, decomp_]:= Module[
		{abelians, brkRep, componentDim, componentReal, decomposition, n, nextConj, pairs, originalReal, dim= 0},
	abelians= Keys@ Select[$GaugeGroups, #@ Group === U1 &];
	
	(*Check original representation*)
	If[!MemberQ[Keys@ $Representations, rep] || 
		!MemberQ[$SSBInfo@ OriginalGroup, $Representations[rep, Group]]
	,
		Message[RepresentationDecomposition::invldorgrep, rep];
		Abort[];
	];
	originalReal= $Representations[rep, Reality] === +1;
	
	(*Ensure that the decomposition is {{r1, r2, ...}, ...}, where each sublist signifies a tensor product 
	of representations from the stability group*)
	decomposition= If[Head@ decomp =!= List, List@ decomp, decomp]/. Singlet-> {};
	decomposition= If[Head@ # =!= List, List@ #, #]&/@ decomposition;
	decomposition= Sort/@ decomposition;
	
	(*Check representations of the decomposition*)
	(*pairs contain information for which blocks in the metric of real original reps are non-zero*)
	(*Complex representations embedded into real representations of the original symmetry are 
		assumed to form conjugate pairs with the neighboring rep*)
	nextConj= None; pairs= {};
	Do[
		brkRep= decomposition[[n]];
		{componentReal, componentDim}= ProductRepsCheckAndInfo@ brkRep;
		dim += componentDim;
		
		(*No further checks if original was complex*)
		If[!originalReal, 
			Continue[];
		];
		nextConj;
		If[componentReal,
			If[nextConj =!= None, (*A real component cannot follow and unpaired complex one*)
				Message[RepresentationDecomposition::realfollowcomplex, brkRep, decomposition[[n -1]], rep];
				Message[RepresentationDecomposition::repord];
				Abort[];
			];
			AppendTo[pairs, {n, n}-> 1];
		,
			Which[
			nextConj === None, (*First in a complex pair*)
				nextConj= ConjugateProductRep@ brkRep;
				AppendTo[pairs, {n, n+ 1}-> 1];
			,brkRep =!= nextConj, (*Doesn't match the previous complex component*)
				Message[RepresentationDecomposition::complexunfollowed, decomposition[[n -1]], rep];
				Message[RepresentationDecomposition::repord];
				Abort[];
			,True,
				nextConj= None;
				AppendTo[pairs, {n, n- 1}-> 1];
			];
		];
	, {n, Length@ decomposition}];
	If[originalReal && nextConj =!= None,
		Message[RepresentationDecomposition::complexunfollowed, Last@ decomposition, rep];
		Message[RepresentationDecomposition::repord];
		Abort[];
	];
	
	(*Check that dimension agrees*)
	If[dim =!= $Representations[rep, Dimension],
		Message[RepresentationDecomposition::incmpdims, dim, $Representations[rep, Dimension]];
		Abort[];
	];	
	
	(*Update the symmetry breaking pattern*)
	$SSBInfo[Representations, rep]= decomposition;
	
	(*Update the behavior of BrokenInds associated with the original representation*)
	If[originalReal,
		pairs= Normal@ SparseArray@ pairs;
		With[{array= pairs, diag= Diagonal@ pairs},
			BrokenInds/: BrokenInds[Index[a_, rep], m_]BrokenInds[Index[a_, rep], n_]:= 
				array[[n, m]];
			BrokenInds/: Power[BrokenInds[Index[a_, rep], m_], 2]:= diag[[m]]; 
		];
	,
		BrokenInds/: BrokenInds[Index[a_, rep], m_]BrokenInds[Index[a_, rep], n_]:= 
			KroneckerDelta[n, m];
		BrokenInds/: Power[BrokenInds[Index[a_, rep], m_], 2]:= 1;
	];
	
	(*Add the automatically infered decomposition of the kronecker delta of the original representation*)
	AddDeltaDecomposition@ rep;
];


(* ::Text:: *)
(*Conjugate a product representation *)


ConjugateProductRep@ reps_List:= Sort[ConjugateProductRep/@ reps];
ConjugateProductRep@ u1_[q:_Integer|_Rational]:= u1[-q];
ConjugateProductRep@ rep_:= Bar@ rep;


(* ::Text:: *)
(*Checks if the representations in the representation product are all valid and associated with unique groups. Returns the reality and dimension of the representation product .*)


RepresentationDecomposition::invldcharge= "`1` is not a rational charge of the stability (remnnant) groups.";
RepresentationDecomposition::invldsumrep= "`1` is not a valid representation the stability (remnnant) groups.";
RepresentationDecomposition::repeatedfactorgroup= "The representation product `1` cannot contain multiple factors associated with the same group.";


ProductRepsCheckAndInfo@ reps_List:= Module[{abelians, group, groups= {}, dim= 1, reality= True},
	abelians= Keys@ Select[$GaugeGroups, #@ Group === U1 &];
	Do[
		Switch[r
		,u1_[_]/; MemberQ[abelians, u1],
			If[!MemberQ[$SSBInfo@ StabilityGroup, Head@ r] || !MatchQ[r, _[_Rational|_Integer]],
				Message[RepresentationDecomposition::invldcharge, r];
				Abort[];
			];
			AppendTo[groups, Head@ r];
			reality= False;
		,b_/; MemberQ[Keys@ $Representations, b/. Bar-> Identity],
			group= $Representations[r/. Bar-> Identity, Group];
			If[!MemberQ[$SSBInfo@ StabilityGroup, group],
				Message[RepresentationDecomposition::invldsumrep, r];
				Abort[];
			];
			If[$Representations[r/. Bar-> Identity, Reality] =!= 1, 
				reality= False;
			];
			AppendTo[groups, group];
			dim*= $Representations[r/. Bar-> Identity, Dimension];
		,_,
			Message[RepresentationDecomposition::invldsumrep, r];
			Abort[];
		];
	, {r, reps}];
	
	If[!DuplicateFreeQ@ groups,
		Message[RepresentationDecomposition::repeatedfactorgroup, reps];
		Abort[];
	];
	
	{reality, dim}
];


(* ::Subsection:: *)
(*Decomposition of objects*)


(* ::Subsubsection::Closed:: *)
(*Field decomposition*)


(* ::Text:: *)
(*Defines replacement rules for substitution of the original fields to constituent parts under the stability group.*)


FieldDecomposition::brkrepfield= "The original field (`1`) is charged under the stability/remnant groups `2` of the symmetry breaking indicating that it is already a broken phase field.";
FieldDecomposition::charge= "The field was specified with a charge `1` at element `2`, which doesn't match the inferred value `3`.";
FieldDecomposition::chargesum= "The charges of the different terms in element `1` do not all agree. They take unique values `2`.";
FieldDecomposition::delta= "The field decomposition does not support the use of 'Delta' (found at position `1`). Use the CG objects 'del' of the appropriate representation instead.";
FieldDecomposition::derivs= "The original field (`1`) must be specified w/o any derivatives.";
FieldDecomposition::dimmismatch= "The dimensionality of the decomposition array does not match that infered from the original field (`1`).";
FieldDecomposition::decompfield= "The field `1` used in the decomposition is charged under the broken symmetry group.";
FieldDecomposition::field= "The first argument (`1`) must be a recognized field (or a projector multiplying a field).";
FieldDecomposition::indsum= "The open indices in the different terms in element `1` do not all agree (`2`).";
FieldDecomposition::indmismatch= "The open indices in element `1` are not consistent with those infered by the decomposition of the original field (`2`).";


FieldDecomposition[NCM[DiracProduct@ _Proj, field_Field], decomposition_]:=
	FieldDecomposition[field, decomposition];


FieldDecomposition[field_, decomp_]:= Module[
		{brokenInds, decompCharges, decompInds, decomposition, depth, dims, fields, fieldInds, indReplace, inferredIndices, 
		inferredCharges, infered, lhs, originalGroups, part, rhs, repInds, \[Mu]},	
	(*Check the field to be decomposed*)
	If[Head@ field =!= Field || !KeyExistsQ[$FieldAssociation, First@ field],
		Message[FieldDecomposition::field, field];
		Abort[];
	];
	originalGroups= Join[
		GroupFromRep/@ $FieldAssociation[First@ field, Indices],
		Head/@ $FieldAssociation[First@ field, Charges] ];
	originalGroups= Intersection[originalGroups, $SSBInfo[StabilityGroup]];
	If[Length@ originalGroups > 0,
		Message[FieldDecomposition::brkrepfield, First@ field];
		Abort[];
	];
	If[Last@ field =!= {},
		Message[FieldDecomposition::derivs, field];
		Abort[];
	];
	
	(*Standardize the decomposition*)
	decomposition= CanonizeDecompositionTensor@ decomp;
	
	(*Check the fields in the decomposition*)
	fields= DeleteDuplicates@ Cases[decomposition, Field[f_, __]-> f, All];
	Do[
		If[IntersectingQ[$SSBInfo@ OriginalGroup, 
				Join[GroupFromRep/@ $FieldAssociation[f, Indices], Head/@ $FieldAssociation[f, Charges] ] ],
			Message[FieldDecomposition::decompfield, f];
			Abort[];
		];
	, {f, fields}];
	
	(*Broken indices, inferred new indices, and inferred new charges*)
	fieldInds= Append[field[[3]], Replace[field[[2]], {Vector@ ind_-> ind, _-> Nothing}]];
	{brokenInds, inferredIndices, inferredCharges}= 
		InferIndicesAndCharges[fieldInds, $FieldAssociation[First@ field, Charges]];
	depth= Length@ brokenInds;
	
	(*Check dimensions of decomposition*)
	dims= Dimensions[inferredCharges, depth];
	If[If[Head@ decomposition === List, Dimensions@ decomposition, {}] =!= dims,
		Message[FieldDecomposition::dimmismatch, dims];
		Abort[];
	];
	
	(*Check that no Deltas occur*)
	If[!FreeQ[decomposition, _Delta],
		Message[FieldDecomposition::delta, FirstPosition[decomposition, expr_/; !FreeQ[expr, Delta], {}, {depth}]];
		Abort[];
	];
	
	(*Checks that open indices and charges match in all components*)
	Do[
		decompInds= If[part =!= {}, Extract[decomposition, part], decomposition];
		If[decompInds === 0, Continue[] ];
		
		(*Check that the charges are consistent between the terms*)
		decompCharges= TermCharges/@ TermsToList@ decompInds// DeleteDuplicates;
		If[Length@ decompCharges > 1, 
			Message[FieldDecomposition::chargesum, part, decompCharges];
			Abort[];
		]; 
		
		(*Check that the charges agree with the inffered ones*)
		infered= If[part =!= {}, Extract[inferredCharges, part], inferredCharges];
		If[First@ decompCharges =!= infered, 
			Message[FieldDecomposition::charge, First@ decompCharges, part, infered];
			Abort[];
		]; 
		
		(*Get the indices of all the terms in the decomposition and check that they are identical*)
		decompInds= TermOpenIndices/@ TermsToList@ decompInds// DeleteDuplicates;
		If[Length@ decompInds > 1, 
			Message[FieldDecomposition::indsum, part, decompInds];
			Abort[];
		]; 
		
		(*Check that the infered and found indices for the decomposition agrees*)
		infered= If[part =!= {}, Extract[inferredIndices, part], inferredIndices];
		If[First@ decompInds =!= infered, 
			Message[FieldDecomposition::indmismatch, part, infered];
			Abort[];
		]; 
	, {part, Tuples[Range/@ dims]}];
	
	(*Make LHS of field substitution pattern*)
	fieldInds= Cases[field, _Index, All];
	indReplace= fieldInds/. ind:Index[lab_, rep_]:> 
		With[{temp= lab}, Rule[ind, Index[Pattern[temp, Blank[]], rep] ] ];
	fieldInds= fieldInds/. Index[l_, _]:> l;
	
	lhs= field/. indReplace;
	lhs[[-1]]= \[Mu]___;
	
	(*Make RHS of field substitution pattern*)
	If[depth > 0,
		rhs= MapIndexed[(#1 Times@@ Thread@ BrokenInds[brokenInds, #2]&), decomposition, {depth}];
		rhs= Plus@@ Flatten[rhs, depth];
	,
		rhs = decomposition;
	];
	rhs= rhs/. Field[rest__, {}]-> Field[rest, \[Mu]];
	
	(*Makes repeated indices in the replacement unique in substitution rule*)	
	repInds= Complement[Cases[decomposition, Index[l_, r_]:> l, All], fieldInds];
	$SSBInfo[FieldSubstitutions, First@ field]= With[
		{temp= rhs, sub= (Index[#, rep_]-> Hold@ Index[Unique@ #, rep]&)/@ repInds},
			RuleDelayed[lhs, temp/. ReleaseHold@ sub]
		];
];


DecomposeFields@ expr_:= expr/. List@@ $SSBInfo@ FieldSubstitutions;


(* ::Text:: *)
(*Determines the open indices of a term, counting conjugated indices *)


TermOpenIndices[term_]:= Module[{barInds, inds, temp, tField},
	temp= term/. Bar@ Field[x1__, indices_, x2_]:> tField[x1, Bar/@ indices, x2];
	barInds= Cases[temp, Bar@ i_Index-> i, All];
	inds= Cases[Tally@ Cases[temp, _Index, All], {i_, 1}-> i];
	inds/. {i_/; MemberQ[barInds, i]:> Bar@ i}// Sort
];


(* ::Text:: *)
(*Determines the charges in a term*)


TermCharges[term_]:= Module[{tField, temp},
	temp= term/. Bar@ Field[x1___]:> tField[x1];
	MergeChargeList@ Flatten@ Join[Cases[temp, Field[f_, __]:> $FieldAssociation[f, Charges], All],
		Cases[temp, tField[f_, __]:> Map[Minus, $FieldAssociation[f, Charges], {2}], All] ]
];


(* ::Text:: *)
(*Canonize decompositions tensor: accepts scalars, tensors, SparseArray, and a list in ArrayRules format*)


CanonizeDecompositionTensor@ decomp_:= 
	Switch[decomp
	,{Rule[{_Integer...}, _]..},
		Normal@ SparseArray@ decomp
	,_SparseArray,
		Normal@ decomp
	,_,
		decomp
	];


(* ::Subsubsection::Closed:: *)
(*CG decomposition*)


(* ::Text:: *)
(*Defines replacement rules for CGs of the original symmetry group. *)
(*The decomposition supports both SparseArray and regular Arrays (nested lists).*)


CGDecomposition::charge= "The CG was specified with a non-zero component at `1` with non-vanishing charge(s) `2`.";
CGDecomposition::delta= "The CG decomposition does not support the use of 'Delta' (found at position `1`). Use the CG objects 'del' of the appropriate representation instead.";
CGDecomposition::dimmismatch= "The dimensionality of the decomposition array does not match that infered from the original CG (`1`). If the indices decompose trivially ensure that the depth of the array is correct.";
CGDecomposition::indmismatch= "The open indices in element `1` are not consistent with those infered by the decomposition of the original CG (`2`).";
CGDecomposition::indsum= "The open indices in the different terms in element `1` do not all agree (`2`).";
CGDecomposition::invldcg= "The CG `1` is not associated with an unbroken group.";
CGDecomposition::invldrep= "No decomposition has been specified for the representation(s) `1` of the original CG. First provide a decomposition with RepresentationDecomposition.";
CGDecomposition::normmismatch= "The squared norm of original CG is `1`, but the squared norm of the decomposition is `2`. The norm must agree.";
CGDecomposition::notcg= "The first argument must be an un-bared CG.";


CGDecomposition[cg_, decomp_]:= Module[
		{brokenInds, cmpl, decompInds, decomposition, depth, dims, indReplace, inferredIndices, 
			inferredCharges, lhs, lhsNorm, part, repInds, rhs, rhsNorm},
	(*Check the original CG*)
	If[Head@ cg =!= CG || Head@ cg[[1]] === Bar,
		Message[CGDecomposition::notcg];
		Abort[];
	];
	If[!MemberQ[$SSBInfo@ OriginalGroup, $CGproperties[cg[[1]], Group]],
		Message[CGDecomposition::invldcg, cg[[1]]];
		Abort[];
	];
	cmpl= Complement[$CGproperties[cg[[1]], Indices]/. Bar-> Identity, Keys@ $SSBInfo@ Representations];
	If[Length@ cmpl > 0,
		Message[CGDecomposition::invldrep, cmpl];
		Abort[];
	];
	
	(*Broken indices, inferred new indices, and inferred new charges*)
	{brokenInds, inferredIndices, inferredCharges}= InferIndicesAndCharges[cg[[-1]], {}];
	depth= Length@ brokenInds;
	
	(*Standardize the decomposition*)
	decomposition= CanonizeDecompositionTensor@ decomp;
	
	(*Check dimensions of decomposition*)
	dims= Dimensions[inferredCharges, depth];
	If[If[Head@ decomposition === List, Dimensions@ decomposition, {}] =!= dims,
		Message[CGDecomposition::dimmismatch, dims];
		Abort[];
	];
	(*Check that no Deltas occur*)
	If[!FreeQ[decomposition, _Delta],
		Message[CGDecomposition::delta, FirstPosition[decomposition, expr_/; !FreeQ[expr, Delta], {}, {depth}]];
		Abort[];
	];
	
	(*Checks open indices and charges match in all components*)
	Do[
		decompInds= Extract[decomposition, part];
		If[decompInds === 0, Continue[] ];
		(*Check that the charges vanish for any non-zero part of the decomposed CG*)
		If[Length@ Extract[inferredCharges, part] > 0, 
			Message[CGDecomposition::charge, part, Extract[inferredCharges, part]];
			Abort[];
		]; 
		
		(*Get the indices of all the terms in the decomposition and check that they are identical*)
		decompInds= TermOpenIndices/@ TermsToList@ decompInds// DeleteDuplicates;
		If[Length@ decompInds > 1, 
			Message[CGDecomposition::indsum, part, decompInds];
			Abort[];
		]; 
		
		(*Check that the infered and found indices for the decomposition agrees*)
		If[First@ decompInds =!= Extract[inferredIndices, part], 
			Message[CGDecomposition::indmismatch, part, Extract[inferredIndices, part]];
			Abort[];
		]; 
	, {part, Tuples[Range/@ dims]}];
	
	(*LHS of the substitution pattern*)
	indReplace= brokenInds/. ind:Index[lab_, rep_]:> 
		With[{temp= lab}, Rule[ind, Index[Pattern[temp, Blank[]], rep] ] ];
	lhs= cg/. indReplace;
	
	(*RHS of the substitution pattern*)
	rhs= MapIndexed[(#1 Times@@ Thread@ BrokenInds[brokenInds, #2]&), decomposition, {depth}];
	rhs= Plus@@ Flatten[rhs, depth];
	
	(*Make repeated indices in the replacement unique*)
	repInds= Complement[Cases[decomposition, Index[l_, r_]-> l, All], Cases[cg[[2]], Index[l_, r_]-> l, All]];
	$SSBInfo[CGSubstitutions, First@ cg]= With[
		{temp= rhs, sub= (Index[#, rep_]-> Hold@ Index[Unique@ #, rep]&)/@ repInds},
			(*Make conjugate replacement pattern if CG has unique conjugate:*)
			If[$CGproperties[First@ cg, UniqueConj],
				{RuleDelayed[lhs, temp/. ReleaseHold@ sub], RuleDelayed[Bar@ lhs, Bar@ temp/. ReleaseHold@ sub]}
			,
				{RuleDelayed[lhs, temp/. ReleaseHold@ sub]}
			]
		];
	
	(*Check normalization of CG substitution*)
	lhsNorm= Simplify@ ContractCGs[cg Bar@cg];
	rhsNorm= Simplify@ ContractCGs@ ReleasePseudoTimes@ DecomposeCGs@ PseudoTimes[cg Bar@cg];
	If[lhsNorm =!= rhsNorm, 
		Message[CGDecomposition::normmismatch, lhsNorm, rhsNorm];
		KeyDropFrom[$SSBInfo@ CGSubstitutions, First@ cg];
		Abort[];
	];
	
	(*Automatically infer decomposition of the structure constants, from the decomposition of generators*)
	If[Head@ cg[[1]] === gen && !MemberQ[$SSBInfo, fStruct@ cg[[1, 1, 0]]],
		AddFStructDecomposition[cg[[1, 1, 0]], cg[[1, 1]]];
	];
];


DecomposeCGs@ expr_:= expr/. Flatten[List@@ $SSBInfo@ CGSubstitutions];


(* ::Subsubsection::Closed:: *)
(*InferIndicesAndCharges*)


(* ::Text:: *)
(*Extract tensors of new indices and charges from the decomposition of an original set of indices and charges (of the broken symmetry)*)


InferIndicesAndCharges::unkwnqdecomp= "The charge(s) `1` transform(s) under the original symmetry but no breaking has been assigned with RepresentationDecomposition.";
InferIndicesAndCharges::unkwnrepdecomp= "The representation '`1`' transforms under the original symmetry but no breaking has been assigned with RepresentationDecomposition.";


InferIndicesAndCharges[inds_List, charges_List]:= Module[
		{abelians, baseCharges, brokenReps, brokenInds, indices, inferredCharges, inferredIndices, 
		newReps, part, qs, rep, unspecified},
	abelians= Keys@ Select[$GaugeGroups, #@ Group === U1 &];
	
	(*Determine the representations involved in the decomposition*)	
	brokenInds= Cases[inds, _Index?(MemberQ[$SSBInfo@ OriginalGroup, GroupFromInd@ #] &), All];
	brokenReps= RepFromInd/@ brokenInds;
	newReps= Table[
			Lookup[$SSBInfo[Representations], rep,
				Message[InferIndicesAndCharges::unkwnrepdecomp, rep]; Abort[];]
		, {rep, brokenReps}];
	
	(*Determine base charges*)
	unspecified= Complement[Intersection[Head/@ charges, $SSBInfo@ OriginalGroup], Keys@ $SSBInfo@ Charges];
	If[Length@ unspecified > 0,
		Message[InferIndicesAndCharges::unkwnqdecomp, unspecified]; Abort[];
	];
	baseCharges= #/. Lookup[$SSBInfo@ Charges, Head@ #, {}]&/@ charges;
	baseCharges= MergeChargeList@ Flatten@ baseCharges;
	 
	(*Iterate over all parts of the tensor decomposition*)
	inferredCharges= inferredIndices= ConstantArray[0, Length/@ newReps];
	Do[
		(*Replace broken indices with the new representation given to the current part*)
		indices= inds/. Thread@ Rule[brokenInds, 
			MapThread[#1/. Index[i_, _]:> Index[i, #2]&, {brokenInds, MapThread[Part, {newReps, part}]}]]; 
		(*Treat and separate singlet representations*)
		indices= indices/. 
			ind:Index[_, _List]:> Thread@ ind/. 
			Bar@ Index[a_, u1_[q_]]/; MemberQ[abelians, u1]:> Index[a, u1[-q]]// Flatten;
		{qs, indices}= SelectAndDeleteCases[indices, Index[_, u1_[q_]]/; MemberQ[abelians, u1]:> u1[q]];
		(*Total charges*)
		qs= MergeChargeList@ Join[qs, baseCharges]; 
		
		inferredIndices[[Sequence@@ part]]= Sort@ indices;
		inferredCharges[[Sequence@@ part]]= qs;
	, {part, Tuples[Range/@ Length/@ newReps]}];
	
	{brokenInds, inferredIndices, inferredCharges}
];


MergeChargeList@ charges_List:= Module[{qs},
	qs= Total/@ GroupBy[charges, Head-> First];
	qs= KeyValueMap[#1@ #2&, qs]/. _[0]-> Nothing;
	Sort@ qs
];


(* ::Subsubsection::Closed:: *)
(*Decomposition of the structure constants*)


(* ::Text:: *)
(*Use the known decomposition of the generators of a representation to determine the decomposition of the structure constants*)


AddFStructDecomposition[group_, rep_]:= Module[{decomposition, dynk, a, b, c, X, Y, Z},
	(*Check group, representation, and generator decomposition*)
	
	dynk= DynkinIndex[$GaugeGroups[group, Group], 
		$Representations[rep, DynkinCoefficients] ];
	decomposition= -I/ dynk ContractCGs@ Expand@ DecomposeCGs[
		CG[gen@ rep, {X, a, b}] (CG[gen@ rep, {Y, b, c}] CG[gen@ rep, {Z, c, a}]- 
			CG[gen@ rep, {Z, b, c}] CG[gen@ rep, {Y, c, a}])]// Expand;
		
	$SSBInfo[CGSubstitutions, fStruct@ group]= With[
		{temp= decomposition, inds= Index[#, group@ adj]&/@ {X_, Y_, Z_}},
			CG[fStruct@ group, inds]:> temp
		];
	
	(*Check normalization of (f^ABC)^2 before/after symmetry breaking*)

]


(* ::Subsubsection::Closed:: *)
(*Decomposition of the Kronecker deltas*)


(* ::Text:: *)
(*Use the knowledge of a representation decomposition to determine the decomposition of a delta*)


AddDeltaDecomposition[brokenRep_]:= Module[
		{abelians, bi1, bi2, conj, decomposition, factor, rep, isRealRep, realRep, repDecomposition, n, a, b},
	abelians= Keys@ Select[$GaugeGroups, #@ Group === U1 &];
	
	isRealRep= $Representations[brokenRep, Reality] === 1;
	repDecomposition= $SSBInfo[Representations, brokenRep];
	{bi1, bi2}= {Index[a, brokenRep], Bar@ Index[b, brokenRep]};
	
	(*Decompose the original delta as a sum of projectors (BrokenInds) multiplying remnant deltas*)
	conj= False; 
	decomposition= Sum[
			rep= repDecomposition[[n]];
			
			realRep= First@ ProductRepsCheckAndInfo@ rep;
			If[isRealRep && !realRep, (*off-diagonal if complex in real rep*)
				conj= !conj;
				If[conj,
					BrokenInds[bi1, n] BrokenInds[bi2, n+ 1]
				,
					BrokenInds[bi1, n] BrokenInds[bi2, n- 1]
				]
			,
				BrokenInds[bi1, n] BrokenInds[bi2, n]
			] *
			Product[(*Ensures delta contractions of the factors of remnant decomposition*)
					If[MatchQ[factor, u1_[_]/; MemberQ[abelians, u1]], 
						1
					,
						If[Head@ factor === Bar,
							CG[Bar[del@@ factor], {a, b}]
						,
							CG[del@ factor, {a, b}]
						]
					]
				, {factor, rep}]
		, {n, Length@ repDecomposition}];
	
	(*Add the substitution pattern to the symmetry breaking information*)
	$SSBInfo[CGSubstitutions, del@ brokenRep]= With[
		{dec= decomposition, inds= {Index[a_, brokenRep], Bar@ Index[b_, brokenRep]}},
			CG[del@ brokenRep, inds]:> dec
		];
]


(* ::Subsection:: *)
(*Derivatives*)


(* ::Subsubsection::Closed:: *)
(*Heavy vector action*)


(* ::Text:: *)
(*The action of a heavy vector on an unbroken field, i.e., a field with indices in the original gauge group*)


VAction[\[Mu]_, f:(Field[lab_, _, inds_, _]| FieldStrength[lab_, _, inds_, _])]:= Module[{ind, charge, newInd, newLab},
	Sum[(*Non-abelian indices*)
		ContractDelta@ If[MemberQ[$SSBInfo@ OriginalGroup, GroupFromInd@ ind],
			newInd= If[Head@ ind === Bar, Bar, Identity]@ Index[newLab, Last[ind/. Bar-> Identity]];
			- I VActionOnIndex[\[Mu], ind, newInd] BrokenPhaseCD[f/. ind-> newInd] 
		,
			0
		]
	, {ind, inds}] + Sum[(*Charges*)
		If[MemberQ[$SSBInfo@ OriginalGroup, Head@ charge],
			- I VActionOnCharge[\[Mu], charge]
		,
			0
		]
	,{charge, $FieldAssociation[lab, Charges]}] * BrokenPhaseCD@ f
]


(* ::Text:: *)
(*Gives the action of the massive vectors on a single index of the original group*)


VActionOnIndex[\[Mu]_, brokenInd_, newInd_]:= Module[{generator, group, l, adjInd},
	group= GroupFromInd@ brokenInd;
	adjInd= Index[l, group@ adj];
	(*Determine unbroken generator associated with the index*)
	generator= If[Head@ brokenInd === Bar,
			-Bar@ CG[gen@ brokenInd[[1, 2]], {adjInd, Bar@ brokenInd, newInd}]
		,
			CG[gen@ brokenInd[[2]], {adjInd, brokenInd, Bar@ newInd}]
		];
	(*Get the unbroken gauge field * generator -> break it and select heavy only*)
	ContractDelta@ Expand@ LightGaugeFieldsToZero@ DecomposeCGs@ DecomposeFields[
		$GaugeGroups[group, Field][\[Mu], l] generator
	]
]


(* ::Text:: *)
(*Gives the action of the massive vectors on a single charge of the original group*)


VActionOnCharge[\[Mu]_, brokenGroup_@ charge_]:= LightGaugeFieldsToZero@ DecomposeFields[
		charge $GaugeGroups[brokenGroup, Field][\[Mu]]
	]


(* ::Text:: *)
(*Remove the unbroken vector fields*)


LightGaugeFieldsToZero@ expr_:= expr/.
	(Field[#, __]:> 0& )/@ List@@ Query[Key/@ $SSBInfo@ StabilityGroup, Key@ Field]@ $GaugeGroups;


(* ::Subsubsection::Closed:: *)
(*Decomposing (derivatives of) fields*)


(* ::Text:: *)
(*We identify the action of a heavy vector in the decomposition [D_mu, f] = [d_mu, f] -i [V_mu, f] of the full covariant derivative on a field. This is the main function for decomposing unbroken fields, as they may have derivatives.*)


(*Fields*)
BrokenPhaseCD[Field[lab__, inds_, {\[Mu]_, \[Nu]___}]]:= Module[{f= Field[lab, inds, {\[Nu]}]},
	CD[\[Mu], BrokenPhaseCD@ f] + VAction[\[Mu], f]
];
BrokenPhaseCD[f:Field[___, {}]]:= DecomposeFields@ f;


(*FieldStrengths*)
BrokenPhaseCD[FieldStrength[rest__, inds_, {\[Mu]_, \[Nu]___}]]:= Module[{f= FieldStrength[rest, inds, {\[Nu]}]},
	CD[\[Mu], BrokenPhaseCD@ f] + VAction[\[Mu], f] 
];
BrokenPhaseCD[FieldStrength[lab_, {\[Mu]_, \[Nu]_}, inds_, {}]]:= Module[{FSs, $a, $b},
	(*Sum of all FS over all constituen vectors. FS defined w.r.t. cov. dev. of the remnant group*)
	FSs= DecomposeFields@ Field[lab, Vector@ \[Mu], inds, {}]/. 
		Field[l_, _, gInds_, {}]:> FieldStrength[l, {\[Mu], \[Nu]}, gInds, {}];
	(*Include [V, V] terms*)
	FSs+ LightGaugeFieldsToZero@ VAction[\[Mu], Field[lab, Vector@ \[Nu], inds, {}]]
];


(* ::Section:: *)
(*Prepare Lagrangian for the broken phase*)


(* ::Subsection:: *)
(*Going to the broken phase*)


(* ::Subsubsection::Closed:: *)
(*Coupling replacements*)


(* ::Text:: *)
(*Function to set a global list of coupling replacement rules invoked during, e.g., symmetry breaking.*)


SetSSBReplacements::notrule= "The argument of SetSSBReplacements must be either a rule or a list of rules.";
SetSSBReplacements::notcoupling= "The LHS of the rule `1` is expected to be a coupling (pattern) of a couopling.";
SetSSBReplacements::indexmismatch= "The open indices of the LHS and RHS of the rule `1` are incompatible.";
SetSSBReplacements::invldobj= "The rule `1` contains Fields/FieldStrengths/Dirac or Lorentz objets. This most certainly shouldn't happen!";
SetSSBReplacements::repeatedcoupling= "The coupling '`1`' occurs on both left- and right-hand sides of the replacement rules. Use a different (new) coupling for replacements to clearly distinguish the new ones.";


SetSSBReplacements[arg_]:= Module[{c, couplings, rhsCouplings, lhsOpen, r, rhsOpen, rules},
	rules= If[Head@ arg =!= List, {arg}, arg];
	If[!MatchQ[rules, {(_Rule|_RuleDelayed)...}],
		Message[SetSSBReplacements::notrule];
		Abort[];
	];
	
	rules= Table[
			If[Head@ First@ r =!= Coupling,
				Message[SetSSBReplacements::notcoupling, r];
				Abort[];
			];
			If[!FreeQ[r, Field|FieldStrength|Lorentz|DiracProduct],
				Message[SetSSBReplacements::invldobj, r];
				Abort[];
			];
			
			(*Check that the LHS and RHS patterns match.*)
			lhsOpen= FindOpenIndices@ First @r ;
			(* Take care of sums on the rhs *)
			rhsOpen= Sort/@ FindOpenIndices/@ TermsToList@ Last@ r;
			rhsOpen= DeleteDuplicates[rhsOpen];
			(* Warning if rhs has mismatch of open indices in itself *)
			If[Length@ rhsOpen > 1, 
				Message[SetSSBReplacements::indexmismatch, r];
				Abort[]
			,
				rhsOpen= First@ rhsOpen;
			];
			
			lhsOpen= lhsOpen/. Index[Verbatim[Pattern][ind_, Blank[]], rep_]:> Index[ind, rep];
			If[Sort@ lhsOpen =!= Sort@ rhsOpen,
				Message[SetSSBReplacements::indexmismatch, r];
				Abort[];
			];
			
			(*Make appropriate patterns out of the open lhs indices.*)
			lhsOpen= (#-> Replace[#, Index[ind_, rep_]:> 
				With[{indPat= ind}, (Index[indPat_, rep])] ])&/@ lhsOpen;
			r[[1]]= r[[1]]/. lhsOpen;
			r
		, {r, rules}];
		
	couplings= Cases[rules[[;;, 1]], Coupling[c_, __]-> c];
	rhsCouplings= Cases[rules[[;;, 2]], Coupling[c_, __]-> c, All];
	couplings= Intersection[couplings, rhsCouplings];
	If[Length@ couplings > 0,
		Do[
		Message[SetSSBReplacements::repeatedcoupling, c];
		, {c, couplings}];
		Abort[];
	];
	
	$SSBCouplingReplacements= rules;
];


(* ::Subsubsection::Closed:: *)
(*ToBrokenPhase*)


(* ::Text:: *)
(*Takes a Lagrangian to the broken phase, decomposing all relevant fields and CGs charged under the original symmetry *)


Options@ ToBrokenPhase= {
		EFTOrder-> Automatic,
		GoldstoneBoson-> True
	};


OptionMessage[EFTOrder, ToBrokenPhase, val_]:= Message[General::optexpectsval, EFTOrder, ToBrokenPhase, val, "positive Integer or 'Automatic'"];


OptionTest[ToBrokenPhase, EFTOrder]= MatchQ[(_Integer? Positive) | Automatic];


ToBrokenPhase::nobreak= "No symmetry breaking pattern specified (use SetSymmetryBreakingPattern).";
ToBrokenPhase::nocgdecomp= "Cannot go to the broken phase without a decomposition for the CG `1` (use CGDecomposition).";
ToBrokenPhase::nofielddecomp= "Cannot go to the broken phase without a decomposition for the field `1` (use FieldDecomposition).";


ToBrokenPhase[lagrangian_, OptionsPattern[]] ? OptionsCheck:= Module[
	{
		cgs, fields, fieldSubstitutionslab, order, out, replaceableCGs,
		lag= lagrangian,
		maxOrder= OptionValue@ EFTOrder
	}
,
	If[$SSBInfo@ OriginalGroup === {}, 
		Message[ToBrokenPhase::nobreak];
		Abort[];
	]; 
	
	If[!CheckLagrangian[lag, 
			CanonicallyNormalized-> False, 
			GaugeAnomalies-> False,
			HeavyMassBasis-> False,
			Hermiticity-> False,
			UndefinedObject-> False], 
		Abort[]
	];
	
	(*Check that all fields have decompositions*)
	fields= DeleteDuplicates@ Cases[lag, Field[lab_, __]|FieldStrength[lab_, __]:> lab, Infinity];
	Do[
		(*Only fields that need decomposisitions proceed*)
		If[!IntersectingQ[$SSBInfo@ OriginalGroup,
			Join[GroupFromRep/@ $FieldAssociation[f, Indices], Head/@ $FieldAssociation[f, Charges] ] ],
			Continue[];
		];
		If[!KeyExistsQ[$SSBInfo@ FieldSubstitutions, f],
			Message[ToBrokenPhase::nofielddecomp, f];
			Abort[];
		];
	, {f, fields}];

	(*Remove compound CGs w/o defined replacements*)
	replaceableCGs= Complement[Keys@ $CGreplacements, Keys@ $SSBInfo@ CGSubstitutions];
	lag= ReplaceCGs[lag, CGsToReplace-> replaceableCGs];
	
	(*Check that all remaining CGs have decompositions*)
	cgs= DeleteDuplicates@ Cases[lag, CG[lab_, __]:> (lab/. Bar-> Identity), Infinity];
	Do[
		(*Only CGs that need decomposisitions proceed*)
		If[!MemberQ[$SSBInfo@ OriginalGroup, $CGproperties[c, Group]],
			Continue[];
		];
		If[!KeyExistsQ[$SSBInfo@ CGSubstitutions, c],
			Message[ToBrokenPhase::nocgdecomp, c];
			Abort[];
		];
	, {c, cgs}];

	If[maxOrder === Automatic, 
		maxOrder= GetMaxEFTOrder[lag, False];
	];
	
	(*For high mass dimensions, compute symmetry breaking per mass dimension for better performance, 
	utilizing suitably truncated expansions.*)
	Monitor[
		out= Sum[
			(* truncate the field substitutions at EFT order maxOrder-order to reduce replacement complexity *)
				BreakLagrangianByOrder[SeriesEFT[lag, EFTOrder-> {order}, HeavyFieldDimension-> False, ExpandDenominator-> True], 
					order, maxOrder, OptionValue@ GoldstoneBoson]
			, {order, maxOrder}];
	, StringForm["Breaking terms at order `1`\[Ellipsis]", order]];
	
	(*Apply basic simplifications *)
	Monitor[
		InternalCollectOperators[out, TrigSimplify-> True, NormalForm-> True]
	, "Collecting operators\[Ellipsis]"]
];


(* ::Text:: *)
(*Break the part Lagrangian of a Lagrangian with terms up to <order> by truncating the expansion and all decompositions to <maxOrder>. With <includGBs> one can decide if GBs should be removed already at level of the field substitutions.*)


BreakLagrangianByOrder[lagrangianPart_, order_Integer, maxOrder_Integer, includeGBs_?BooleanQ]:= Module[
		{fieldSubstitutions, out, replacementOrder},
	replacementOrder= maxOrder- order;
	
	(* Perform specified coupling replacements *)
	out= ReplaceInLagrangian[lagrangianPart, TruncateReplacements[$SSBCouplingReplacements, replacementOrder], 
		Simplify-> False];

	(* Sums of couplings in denominators have to be expanded *)
	out= SeriesEFT[out, EFTOrder-> maxOrder, HeavyFieldDimension-> False, ExpandDenominator-> True];
	
	(*Decompose fields and CGs*)
	(*Temporarily changes global field substitutions to improve efficiency*)
	(*Consider finding some non-global way of doing this*)
	fieldSubstitutions= $SSBInfo@ FieldSubstitutions;
	$SSBInfo@ FieldSubstitutions= TruncateFieldSubstitutions[$SSBInfo@ FieldSubstitutions, replacementOrder, includeGBs];
	out= DecomposeExpr@ out;
	$SSBInfo@ FieldSubstitutions= fieldSubstitutions;
	
	(*Truncate EFT series again *)
	SeriesEFT[out, EFTOrder-> maxOrder, HeavyFieldDimension-> False, ExpandDenominator-> True]
];


(* ::Text:: *)
(*Perform the decomposition and expansion of all operators in the expression. *)
(*- New version groups products of symmetric-phase fields to facilitate fast contraction of the resulting tensor-network*)
(*N.b. PseudoTimes instead of Operator to prevent rewriting EOM objects. *)


DecomposeExpr@ expr_:= Module[{out},
	out= PseudoTimes@ RelabelIndices[expr, Unique-> True];
	Monitor[
		out= out/. op_PseudoTimes:> DecomposeOperator@ op;
	, "Decomposing CGs and Fields\[Ellipsis]"];
	Monitor[
		Contract@ ContractCGs@ ContractDelta@ out
	, "Contracting along broken indices\[Ellipsis]"]
]


(*DecomposeExpr@ expr_:= Module[{out},
	out= PseudoTimes@ RelabelIndices[expr, Unique-> True];
	Monitor[
		out= DecomposeCGs@ out/. f:(_Field| _FieldStrength):> BrokenPhaseCD@ f;
	, "Decomposing CGs and Fields\[Ellipsis]"];
	Monitor[
	Contract@ ContractCGs@ ContractDelta@ LagrangianExpand@ ReleasePseudoTimes@ out
	, "Contracting along broken indices\[Ellipsis]"]
];*)


(* ::Text:: *)
(*Groups products and successively unwraps for improved performance in large contractions. *)


DecomposeOperator@ op_PseudoTimes:= Module[{out},
	out= GroupBrokenContractions@ op;
	out= DecomposeCGs@ out/. f:(_Field| _FieldStrength):> BrokenPhaseCD@ f;
	out//. c_$contraction/; FreeQ[c, $contraction, {2, Infinity}]:> Expand[Times@@ c]
];


(* ::Text:: *)
(*Truncates the SSB replacements at EFTOrder order*)


TruncateReplacements[replacements_List, order_Integer]:= Map[
	RuleDelayed[First@ #, SeriesEFT[Last@ #, EFTOrder-> GetMinEFTOrder[Last@ #] +order]]&,
	replacements
];


(* ::Text:: *)
(*Truncates the field substitutions at EFTOrder order*)


(*TruncateFieldSubstitutions[replacements_Association, order_Integer]:= Module[{rhs, keys},
	rhs = ReplaceAll[
		Map[
			SeriesEFT[Last[#], EFTOrder-> GetMinEFTOrder[Last[#]]+order]&,
			(List@@replacements /. HoldPattern@ReplaceAll[expr_,_]:>expr) /. Field[rest__, inds_, \[Mu]_]:>Field[rest, List[inds,\[Mu]], {}]//Echo
	 	],
		Field[rest__, {inds_List,\[Mu]_}, {}]:> Field[rest, inds, \[Mu]]
	]//Echo;
	ReplacePart[replacements, Echo@ Thread[Table[{Key@ key, 2, 1}, {key, Keys@ replacements}]-> rhs]]
];*)


TruncateFieldSubstitutions[replacements_Association, order_Integer, includeGBs_?BooleanQ]:= Module[
		{gbSub, keys, rhs},
	gbSub= If[includeGBs, {}, Thread[GetFieldsByProperty[GoldstoneBoson-> True]-> 0] ];
	
	(*The replacement part for the field is part {2, 1} (wrapped in a Replace all for the repeated indices)*)
	(*The derivative pattern used in the fields, is a symbol on the r.h.s. SeriesEFT does not trigger 
	properly on it and it temporarily moved to the other indices while perfoming the truncation.*)
	rhs= KeyValueMap[Function[{field, replacement}, 
			{Key@ field, 2, 1}-> SeriesEFT[#, EFTOrder-> GetMinEFTOrder[#]+ order, HeavyFieldDimension-> False]&@
				(replacement[[2, 1]]/. gbSub/. Field[rest__, inds_, \[Mu]_]:> Field[rest, {inds, \[Mu]}, {}])
				/. Field[rest__, {inds_List,\[Mu]_}, {}]:> Field[rest, inds, \[Mu]]
		], replacements];
	
	ReplacePart[replacements, rhs]
];


(* ::Subsubsection::Closed:: *)
(*Group broken index contractions*)


(* ::Text:: *)
(*Performs a nested groupings of the factors in an operator that contain broken indices. *)
(*- The nested products are wrapped in a "$contraction" head. *)


GroupBrokenContractions@ op_PseudoTimes:= Module[
		{adjacencies, grouping, ind, indices, nestingRules, objects, pos, representations, trivial},
	representations= Length/@ $SSBInfo@ Representations;
	objects= List@@ op; 
	(*Separate out trivial objects*)
	indices= DeleteDuplicates@ Cases[objects, Index[_, rep_]/; MemberQ[Keys@ representations, rep], All];
	pos= DeleteDuplicates[First/@ Position[objects, ind_Index/; MemberQ[indices, ind]]];
	{trivial, objects}= {objects[[Complement[Range@ Length@ objects, pos]]], objects[[pos]]};
	
	(*Construct adjacency association*)
	(*Each edge carries weight equal to the log of the number of representation in the corresponding index 
	contraction. This let's us determine the "size" of subgraphs by adding up the boundary weights to the rest of the graph.
	Efficient contractions should avoid large boundaries, which results in a proliferation of intermediate terms.*)
	adjacencies= Association@@ Array[#-> <||>&, Length@ objects];
	Do[
		pos= First/@ Position[objects, ind];
		adjacencies[pos[[1]], pos[[2]]]= adjacencies[pos[[2]], pos[[1]]]= N@ Log@ representations@ Last@ ind
	, {ind, indices}];
	
	(*Obtains a nesting of the product. It replaces subgraphs with new effective nodes. 
	The nesting rules can then unwrap the nested products*)
	{grouping, nestingRules}= RecursiveGrouping[adjacencies, {}, Length@ adjacencies+ 1];
	grouping= Keys@ grouping//. nestingRules; (*e.g. {{1,3}, {2, {{4, 6}, 5}}}*)
	grouping= Replace[grouping, List-> $contraction, {2, Infinity}, Heads-> True]/. n_Integer:> objects[[n]];
	
	Times@@ trivial* Times@@ grouping
];


(* ::Text:: *)
(*Groups the vertices for a given adjacency association.*)
(*Greedy algorithm recursively replaces to neighboring nodes with a new effective node (with number "next"), based on how much the product size is reduced. *)


RecursiveGrouping[adjacencyAssoc_Association, rules_List, next_Integer]:= Module[
		{min, newAdj, nodes, pair, newRules, score, len= Length@ adjacencyAssoc},
	nodes= Keys@ adjacencyAssoc;
	
	(*Finds the best pair (greedy algorithm)*)
	min= 10000.;
	Do[(*Repeat accross all possible pairs*)
		score= PairingScore[adjacencyAssoc, nodes[[{m, n}]]];
		If[score < min, 
			min= score; 
			pair= nodes[[{m, n}]];
		];
	, {m, len}, {n, m+ 1, len}];
	
	(*Stops if no neighboring pair exists; disconected components are possible*)
	If[min > 1000., Return@ {adjacencyAssoc, rules}];
	
	(*Replaces the pair with an effective node*)
	newAdj= MergePairs[adjacencyAssoc, pair, next];
	newRules= Append[rules, next-> pair]; 
	
	(*Continues recursion*)
	RecursiveGrouping[newAdj, newRules, next+ 1]
];


(* ::Text:: *)
(*Determines merge score for selected pair. Impossibly large default if not neighboring*)


PairingScore[adjacencyAssoc_Association, {u_Integer, v_Integer}]:= Block[{},
	If[!KeyExistsQ[adjacencyAssoc[u], v], Return[10000.];]; (*Very high if no edge found: we shouldn't pair unconnect vertices*)
	(*The size of the boundary of the new effective vertex*)
	Total@ adjacencyAssoc@ u + Total@ adjacencyAssoc@ v - 2 adjacencyAssoc[u, v]
];


(* ::Text:: *)
(*Merges a pair from an adjacency association to form a single new vertex*)


MergePairs[adjacencyAssoc_Association, {u_Integer, v_Integer}, new_Integer]:= Module[{assoc, newEntry, temp},
	newEntry= KeyDrop[Merge[{adjacencyAssoc[u], adjacencyAssoc[v]}, Total], {u, v}];
	assoc= KeyDrop[adjacencyAssoc, {u, v}];
	Do[
		temp= Total@ Lookup[assoc[n], {u, v}, 0]; 
		assoc[n]= KeyDrop[assoc[n], {u, v}];
		assoc[n][new]= temp
	, {n, Keys@ newEntry}];
	assoc[new]= newEntry;
	assoc
];


(* ::Subsection:: *)
(*Prepare Lagrangian *)


(* ::Text:: *)
(*Functions to prepare the broken phase Lagrangian for further computations*)


(* ::Subsubsection::Closed:: *)
(*ImplementVacuumConditions *)


ImplementVacuumConditions::incKinTerm= "The Lagrangian cannot be canonically normalized. The kinetic term involving the fields `1` is normalized to `2`."; 
ImplementVacuumConditions::invldout= "Could not bring the Lagrangian to a valid vacuum form.";
ImplementVacuumConditions::symmetricPhase= "The input Lagrangian looks to be in the symmetric phase. Vacuum conditions should only be applied to a Lagrangian already expanded out in the broken phase eigenstates. Maybe you fogot to apply ToBrokenPhase?";


(* ::Text:: *)
(*Takes a Lagrangian in the broken phase and determines the conditions on the coefficients required to put it on a sensible form.*)
(*ImplementVacuumConditions removes all forbidden terms from the Lagrangian, implicitly assuming that a consistent solution of for the parameters exists.  *)


Options@ ImplementVacuumConditions= {
		SubstituteMasses-> True
	};


OptionMessage[SubstituteMasses, ImplementVacuumConditions, val_]:= Message[General::optexpectsval, SubstituteMasses, ImplementVacuumConditions, val, "True or False"];


OptionTest[ImplementVacuumConditions, SubstituteMasses]= MatchQ[True| False];


ImplementVacuumConditions[lag_, OptionsPattern[]] ? OptionsCheck:= Module[{out, effConditions, flat},
	{out, effConditions}= Most@ VacuumConditionsCore@ lag;
	
	(*Constrain Lagrangian*)
	Monitor[
		flat= Flatten[List@@ List@@@ effConditions@ "Conditions"];
		out= out/. flat// ReplaceEffectiveCouplings;
		(*Optionally write the mass and vector-GB mixing terms in terms of the field mass symbols*)
		If[OptionValue@ SubstituteMasses,
			out= SubstituteMassesInLagrangian@ out;
		];
		out= out// GreensSimplify;
	, "Simplifying Lagrangian\[Ellipsis]"];
	
	(*Check that the Lagrangian now has canonically normalized kinetic terms*)
	If[!CheckLagrangian[out, 
		GaugeAnomalies-> False, 
		HeavyMassBasis-> False
	],
		Message[ImplementVacuumConditions::invldout];
		(*Output likely preferred with the warning to help debug.*)
	];

	out
];


(* ::Text:: *)
(*Core module memoizes the input Lagrangian*)


VacuumConditionsCore[lag_]:= VacuumConditionsCore[lag]= Module[
		{out, terms, outLeading, conditions, effConditions= <|"Conditions"-> <||>|>},
	(*Remove vacuum*)
	
	out= lag- (lag/. (_Field| _FieldStrength)-> 0);
	If[!CheckLagrangian[out,
			CanonicallyNormalized-> False,
			FreeOfHeavyTadpoles-> False,
			GaugeAnomalies-> False,
			HeavyMassBasis-> False
	],
		Abort[]
	];
	
	(*Check if Lagrangian is in a broken phase*)
	If[!BrokenPhaseLagQ@ lag,
		Message[ImplementVacuumConditions::symmetricPhase];
		Abort[];
	];

	(* introduce effective couplings and simplify again *)
	out= IntroduceEffectiveCouplings@ out;
	terms= InternalSimplify[out, InternalOpRepresentation-> True, VLFBasis-> Automatic, TrigSimplify-> True];

	(*Determine coefficients*)
	Monitor[
		(*Kinetic terms*)
		effConditions["Conditions", "KineticTerms"]= ExtractKineticTermConditions@ AtomicToOperatorForm@ terms;
		(*Tadpole coefficients*)
		effConditions["Conditions", "Tadpoles"]= ExtractTadpoleConditions@ terms;
		(*Off-diagonal masses*)
		{effConditions["Conditions", "Masses"], effConditions@ "MassSubstitutions"}= 
			ExtractMassConditionsAndMasses[out, terms];
		(*Coefficients of non-GB mixing terms*)
		effConditions["Conditions", "ScalarVectorMixings"]= ExtractScalarVectorMixing[out, terms];
		(*Goldstone boson masses*)
		effConditions["Conditions", "GBMasses"]= ExtractGBMassConditions[out, terms];
	, "Determining conditions\[Ellipsis]"];

	(*Put simplification conditions in terms of the original couplings*)
	conditions= effConditions;
	conditions@ "Conditions"= Query[All, All, ReplaceEffectiveCouplings]@ conditions@ "Conditions";
	conditions@ "MassSubstitutions"= Query[All, ReplaceEffectiveCouplings]@ conditions@ "MassSubstitutions";
	
	{out, effConditions, conditions}
];


(* ::Subsubsection::Closed:: *)
(*Conditions on tadpole coefficients*)


(* ::Text:: *)
(*Finds the coefficients of all tadpole terms and set them to zero*)


ExtractTadpoleConditions@ expr_:= Module[{terms},
	(*terms= expr/. (AtomicOp| CompOp)[{Except@{{_}, 0}, _}, _]-> 0;*)
	terms= expr/. {AtomicOp[Except@ {{_}, 0}, __]-> 0, CompOp[_, Except@ {{_}, 0}, __]-> 0};
	If[terms === 0, Return@ <||>;];
	terms= TermsToList@ terms;
	terms= terms/. Times[AtomicOp[{{f_}, 0}, __], rest___]:> Rule[f, Times@ rest];
	Rule[#, 0]&/@ Association@@ terms
]


(* ::Subsubsection::Closed:: *)
(*Conditions on GB masses*)


(* ::Text:: *)
(*Finds all instances of GB mass terms and determines conditions for elimination*)


ExtractGBMassConditions[lag_, lagInt_]:= Module[{terms, gbs},
	gbs= Intersection[GetFieldsByProperty[GoldstoneBoson-> True], Cases[lag, Field[f_, __]:>f, All]];
	
	terms= lagInt/. {AtomicOp[Except@ {fields:{_, _} /; SubsetQ[gbs, fields/. Conj-> Identity], 0}, __]-> 0,
		_CompOp-> 0};
	If[terms === 0, Return@ <||>;];
	terms= TermsToList@ terms;
	
	terms= terms/. Times[AtomicOp[{fields:{_, _}, _}, __], rest___]:> (fields-> (Times@ rest-> 0));
	Association@@ terms
];


(* ::Subsubsection::Closed:: *)
(*Conditions on off-diagonal masses*)


(* ::Text:: *)
(*Finds all instances of off-diagonal masses and determines the coefficients required to eliminate them. *)
(*Also finds all mass substitutions for the radial scalar modes and massive vectors.*)


General::masslessSSBField= "The radial scalar/massive vector field(s) `1` acquire a mass term in the broken phase but were defined as massless (Mass-> 0) in DefineField. Their mass substitutions are omitted; define them with a non-zero mass symbol to use the mass-substitution machinery (SubstituteMasses / the \"MassSubstitutions\" output of GetVacuumConditions).";


ExtractMassConditionsAndMasses[lag_, terms_]:= Module[
		{conditions, expr, fields, fieldLabels, masses, massSubstitutions, i, j},
	(*Group fields by quantum numbers*)
	fields= DeleteDuplicates@ Cases[lag, Field[f_, __]:>f, All];
	fields= $FieldAssociation[[Key/@ fields, Key/@ {Type, Indices, Charges, Chiral}]];
	fields= DeleteCases[fields, KeyValuePattern[Chiral-> LeftHanded|RightHanded]]; (*removes chiral fermions*)
	(*fields= List@@(Keys/@ GroupBy[MapAt[Sort, fields, {All, Key/@{Indices, Charges}}], Identity]);*)
	(*The above does not work in Mathematica v14.0, therefore replaced it with a nested MapAt*)
	fields= List@@(Keys/@ GroupBy[MapAt[Sort, MapAt[Sort, fields, {All, Key@Indices}], {All, Key@Charges}], Identity]);
	
	expr= AtomicToNormalForm[terms/. 
		{AtomicOp[Except@ {{_, _}, 0}, __]-> 0, CompOp[_, Except@ {{_, _}, 0}, __]-> 0}];
	masses= Association@@ (#-> Simplify@ MassMatrix[#, expr, Global`d$$1, Global`d$$2]&)/@ fields;
	
	(*Mass substitutions for heavy vectors and scalar radial modes (no GBs)*)
	(*Doesn't support flavor indices in the current version*)
	fieldLabels= Join[HeavyVectorLabels[], RadialScalarLabels[]];
	massSubstitutions= Flatten@ KeyValueMap[
			Table[
				If[#2[[i, i]]=!= 0 && MemberQ[fieldLabels, #1[[i]]],
					(*The field must carry a mass symbol to be rewritten; warn and skip if it was defined massless.*)
					If[$FieldAssociation[#1[[i]], Mass] === 0,
						Message[General::masslessSSBField, #1[[i]]];
						Nothing
					,
						(*The mass matrix eigenvalue is the squared mass; the field mass symbol denotes the (linear) mass.*)
						Rule[#1[[i]], $FieldAssociation[#1[[i]], Mass][]^2-> #2[[i, i]]*
							If[$FieldAssociation[#1[[i]], Type] === Scalar, -1, 1]]
					]
				,
					Nothing
				]
			, {i, Length@ #1}]&
		, masses];
	
	conditions= Rule[#, 0]&/@ Association@@ Flatten@ KeyValueMap[
			Table[
				If[#2[[i, j]]=!= 0, {#1[[i]], #1[[j]]}-> #2[[i,j]], Nothing]
			, {j, Length@ #1}, {i, j- 1}]&
		, masses];
	
	{conditions, massSubstitutions}
];


(* ::Text:: *)
(*Mass matrix*)


MassMatrix[fields_, lag_, i_, j_]:= Module[{deltas, f1, f2, fs, fBars, gaugeInds, fieldType, indNo},
	(*For projecting over open gauge indices*)
	gaugeInds= DeleteCases[$FieldAssociation[First@ fields, Indices], 
		_? (!MemberQ[Keys@ $GaugeGroups, GroupFromRep@ #]&)];
	deltas= Times@@(Delta[Index[i, #],Index[j, #]]&/@ gaugeInds);
	deltas= deltas/ Power[deltas, 2];
	
	(*Fields used in taking functional derivatives*)
	fieldType= $FieldAssociation[First@ fields, Type];
	{fBars, fs}= Transpose@ Table[
			indNo= Length@ $FieldAssociation[f1, Indices]+ If[fieldType === Vector, 1, 0];
			{Bar[f1@@ ConstantArray[i, indNo]], f1@@ ConstantArray[j, indNo]}
		, {f1, fields}];
	
	Map[Contract, Table[VarD[RelabelIndices[lag, Unique->True], f1, f2]/. (_Field| _FieldStrength)-> 0
		, {f1, fBars}, {f2, fs}] deltas If[fieldType === Vector, Metric[i, j]/\[ScriptD], 1], {2}]
];


(* ::Subsubsection::Closed:: *)
(*Conditions on kinetic mixing of non-GB scalars with vectors*)


(* ::Text:: *)
(*Finds all cases of kinetic mixing between heavy vectors and non-GB scalars and determine what coefficients are put to zero*)


ExtractScalarVectorMixing[lag_, expr_]:= Module[{vectors, terms, nonGBscalars},
	vectors= DeleteDuplicates@ Cases[lag, Field[f_, _Vector, __]:> f, All];
	(*Keep all scalar-vector mixing terms*)
	(*terms= expr/. (AtomicOp| CompOp)[{Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, _}, _]-> 0;*)
	terms= expr/. {AtomicOp[Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, __]-> 0,
		CompOp[_, Except@ {{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, __]-> 0};
	If[terms === 0, Return@ <||>;];
	(*Remove the Goldstone terms*)
	nonGBscalars= Alternatives@@ GetFieldsByProperty[Type-> Scalar, GoldstoneBoson-> False];
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];

	terms= terms/. Times[AtomicOp[{f_, _}, __], rest___]:> 
		If[FreeQ[f, nonGBscalars], Nothing, Rule[f, Times@ rest]];
	Rule[#, 0]&/@ Association@@ terms	
]


(* ::Subsubsection::Closed:: *)
(*Conditions on non-canonical kinetic terms*)


(* ::Text:: *)
(*Finds all kinetic terms and determine the coefficients required for canonical normalization and elimination of all off-diagonal terms*)


ExtractKineticTermConditions[lag_]:= Module[{gaugeFields, inconsistencies, terms, var},
	(*Might have to use simplified Lagrangian*)
	
	(*Remove all gauge-kinetic terms: their normalization is unkown due to the freedom of the gauge coupling*)
	gaugeFields= Values[$GaugeGroups][[All, Key@ Field]];
	terms= Operator@ lag/. Operator[FieldStrength[lab1_, __], FieldStrength[lab2_, __], ___]/; 
		SubsetQ[gaugeFields, {lab1, lab2}]-> 0;
	terms= Collect[terms, _Operator];
	If[terms === 0, Return@ <||>;];
	terms= If[Head@ terms === Plus, List@@ terms, {terms}];
	
	(*Finds all instances of kinetic terms*)
	terms= Cases[terms, op_Operator|Times[op_Operator, rest__]/; KineticOpQ@ op:> 
		{op, Times@ rest}];
	(*Set coefficients equal to their canonical values*)
	terms= Rule[Cases[#1, Field[f_, __]|FieldStrength[f_, __]:> f, All], 
		#2 == KineticTermNormalization@ #1]&@@@ terms;
	
	(*Check for inconsistencies*)
	inconsistencies= DeleteCases[terms, Except@ Rule[_, False]];
	If[Length@ inconsistencies> 0,
		KeyValueMap[(Message[ImplementVacuumConditions::incKinTerm, #2, First@ #1]&), Association@@ inconsistencies];
		Abort[];
	];
	
	(*Solve non-trivial conditions*)
	terms= DeleteCases[terms, HoldPattern@ Rule[_, True]];
	terms[[;;, 2]]= SmartSolve/@ terms[[;;, 2]];
	
	Association@@ terms
];


(* ::Text:: *)
(*Standard normalization for kinetic terms *)


KineticTermNormalization@ op_:= op/.{
		(*Scalars*)
		HoldPattern@ Operator[Bar@Field[f_, Scalar, _, {\[Mu]_}], Field[f_, Scalar, _, {\[Mu]_}]]-> 1,
		HoldPattern@ Operator[Field[f_, Scalar, _, {\[Mu]_}], Field[f_, Scalar, _, {\[Mu]_}]]-> 1/2,
		HoldPattern@ Operator[Bar@Field[f_, Scalar, _, {}], EoM@ Field[f_, Scalar, _, {}]]-> -1,
		HoldPattern@ Operator[Field[f_, Scalar, _, {}], EoM@ Field[f_, Scalar, _, {}]]-> -1/2,
		HoldPattern@ Operator[Field[f_, Scalar, _, {}], EoM@ Bar@ Field[f_, Scalar, _, {}]]-> -1,
		(*Fermions*)
		(* N.b. one cannot use  CenterDot with HoldPattern *)
		HoldPattern@ Operator@ NCM[Bar@ Field[f_, Fermion, _, {}], EoM@ Field[f_, Fermion, _, {}] ]-> I,
		HoldPattern@ Operator@ NCM[Bar@ Field[f_, Fermion, _, {}], DiracProduct@ _Proj,
			EoM@ Field[f_, Fermion, _, {}] ]-> I,
		HoldPattern@ Operator@ NCM[Transp@ Field[_, Fermion, _, {}], DiracProduct[GammaCC], 
			EoM@ Field[_, Fermion, _, {}] ]-> I/2,
		(* self-conjugate fields cannot be chiral *)
		(*HoldPattern@ Operator@NCM[Transp@ Field[_, Fermion, _, {}], DiracProduct[GammaCC, _Proj],
			EoM@ Field[_, Fermion, _, {}] ]-> \[ImaginaryI],*)
		(*Vectors*)
		HoldPattern@ Operator[FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}],
			FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}]] -> -1/4,
		HoldPattern@ Operator[Bar@ FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}],
			FieldStrength[f_, {\[Mu]_, \[Nu]_}, a_, {}]] -> -1/2,
		HoldPattern@ Operator[Field[f_, Vector@ \[Mu]_, a_, {}], EoM@ Field[f_, Vector@ \[Mu]_, a_, {}]] -> 1/2,
		HoldPattern@ Operator[Bar@ Field[f_, Vector@ \[Mu]_, a_, {}], EoM@ Field[f_, Vector@ \[Mu]_, a_, {}]] -> 1,
		HoldPattern@ Operator[Field[f_, Vector@ \[Mu]_, a_, {}], EoM@ Bar@ Field[f_, Vector@ \[Mu]_, a_, {}]] -> 1,
		(*Off-diagonals*)
		_Operator-> 0
		}; 


(* ::Text:: *)
(*A custom solve that attempts to find a good variable to solve in, by looking for an equality of the form *)
(*	0 = a + b x^n *)


SmartSolve@ equality:Equal[lhs_, rhs_]:= Module[
		{candidates, eq, expr, powers, var, vars, temp, tempExpr, terms},
	expr= lhs- rhs;
	terms= TermsToList@ expr;
	(*Consider adjusting for the posibility of coliding symbol and coupling names*)
	vars= DeleteDuplicates@ Cases[equality, _Coupling, All];
	tempExpr= equality/. {_Coupling-> 1, _Delta-> 1};
	
	(*TEMPORARY*)
	If[Length@ vars === 0, 
		Return@ Nothing
	];
	
	(*Check if expression is of the form a+ b var^n = 0 for some var*)
	candidates= Table[Catch[
			If[!PolynomialQ[expr, var], Throw[Nothing]];
			(*What powers of the variable is present in the terms*)
			powers= DeleteDuplicates@ Cases[terms, Power[var, n_]-> n, All];
			Switch[Length@ powers
			,0,
				{var, 1}
			,1,
				If[FreeQ[expr/. Power[var, _]-> 1, var],
					{var, First@ powers}
				,
					Nothing
				]
			,_,
				Nothing
			]
		], {var, vars}];
	
	(*If not candidate found*)
	If[Length@ candidates === 0, 
		(*Default variable: Improve this*)
		Return@ First@ Solve[equality, First@ var]
	];
	
	(*Solve in the smalest power of the symbols*)
	var= First@ MinimalBy[candidates, Last];
	If[Last@ var > 1,
		eq= equality/. Power@@ var-> temp;
		First@ First@ Solve[eq, temp]/. temp:> Power@@ var
	,
		First@ First@ Solve[equality, First@ var]
	]
]


(* ::Subsubsection::Closed:: *)
(*GetVacuumConditions *)


(* ::Text:: *)
(*GetVacuumConditions returns only the conditions.*)


Options@ GetVacuumConditions= {
		TimeConstraint-> 60,
		Variables-> Automatic
	};


GetVacuumConditions::dbg= "Debugging information: `1`";
GetVacuumConditions::eqform= "Invalid form of the orignal equations: assuming equations of the form <expr> == <number>.";
GetVacuumConditions::invldcoupling= "The symbol `1` does not correspond to any coupling label. Invalid in Variables option.";
GetVacuumConditions::inconsistent= "Found a solution but it seems to be inconsistent. Verify for yourself that the result is consistent.";
GetVacuumConditions::time= "Solution routine exceeded time constraint of `1` and aborted. You can adjust option TimeConstraint to allocate additional solution times.";
GetVacuumConditions::tensorcoupling= "`1` is tensor coupling. Can only solve for scalar couplings."; 
GetVacuumConditions::unexpct= "Encountered unexpected equation form: check manually if equation system is consistent.";


(*Converts the substitutions to equalities*)
GetVacuumConditions[lag_, OptionsPattern[]]? OptionsCheck:= Module[
		{out, conditions, eftOrd, priorityCouplings= OptionValue@ Variables},
	conditions= Last@ VacuumConditionsCore@ lag;
	conditions@ "Conditions"= Apply[Equal, conditions@ "Conditions", {2}];
	eftOrd= GetMaxEFTOrder[lag, False];
	
	priorityCouplings= Switch[priorityCouplings
		, Automatic,
			{}
		, _List, 
			priorityCouplings
		, _, 
			{priorityCouplings}
		];
	priorityCouplings= Table[
			If[!KeyExistsQ[$CouplingAssociation, c],
				Message[GetVacuumConditions::invldcoupling, c];
				Abort[];
			];
			If[Length@ $CouplingAssociation[c, Indices] > 0, 
				Message[GetVacuumConditions::tensorcoupling, c];
				Abort[];
			];
			c[]
		, {c, priorityCouplings}];
	
	(*Solve the conditions*)
	Monitor[
		conditions@ "Solution"= TimeConstrained[ReduceVacuumConditions[conditions@ "Conditions", eftOrd, priorityCouplings],
			OptionValue@ TimeConstraint, Message[GetVacuumConditions::time, OptionValue@ TimeConstraint]; $TimedOut];
	, "Solving system\[Ellipsis]"];
	
	(*Present the mass substitutions*)
	conditions@ "MassSubstitutions"= Replace[conditions@ "MassSubstitutions", 
		Rule[_, r_Rule]:> r, {1}];
	
	conditions
];


(* ::Text:: *)
(*Function to find a consistent solution to the vacuum conditions *)


ReduceVacuumConditions[vacCond_Association, eftOrd_Integer, priorityCouplings_List]:= Module[
		{delta, dim, eqs, loopResult, nonZero, ord, originalEqs, partialSol, solution, 
			substitutions, type, variables},
	(*delta is a Head used to denote perturbative shifts of couplings*)
	delta@ delta@ x_:= delta@ x; 
	
	(*Gather equations from all conditions*)
	originalEqs= Flatten[List@@@ List@@ vacCond];
	
	(*Check that all equations are of the form <something> == 0*)
	If[Length@ originalEqs === 0, 
		Return@ True;
	];
	If[!MatchQ[originalEqs, {Equal[_, _?NumberQ]..}],
		Message[GetVacuumConditions::eqform];
		Abort[];
	];
	originalEqs= Replace[originalEqs, Equal[a_, b_?NumberQ]:> Equal[a- b, 0], {1}];
	
	(*Determines the leading EFT order of each equation*)
	eqs= originalEqs= {#, GetMinEFTOrder@ #[[1]]}&/@ originalEqs;
	
	(*Get sorted list of variables to solve for*)
	{variables, nonZero}= GetEquationVariables@ eqs[[;;, 1]];
	variables= Complement[variables, priorityCouplings];
	variables= Join[variables, priorityCouplings];
	
	(*Solve equation system perturbatively in EFT order. Up to EFTOrder - 4 above the leading order. 
	Solved variables get replaced by solution + perturbative off-set for next order*)
	partialSol= <||>; (*association w/ solutions up to current order*)
	loopResult= Do[
		solution= SolveSystemAtOrder[eqs, variables, nonZero, ord];
		(*delta-variables == expr, preffered over other parts of 'Or' conditions *)
		solution= solution//. Or[a___, eq:Equal[Coupling[_delta, __]|Bar@ Coupling[_delta, __], _], b___]:> eq;
		
		(*Remove degenerate solution from squared variables*)
		solution= solution//. Or[a___, x_ == y1_, x_ == y2_, b___]/; y1 == -y2:> 
			Module[{y},
				y= If[Count[y1, _?Negative, All] > Count[y2, _?Negative, All], y2, y1];
				Or[a, x == y, b]
			];
		
		(*Check form of the solutions*)
		Switch[solution
		,False, 
			Return[False];
		,_Or, 
			Message[GetVacuumConditions::unexpct];
				Return[solution];
		,True,
			Continue[];
		];
		solution= If[Head@ solution === And, List@@ solution, List@ solution];
		
		Do[
			If[!MatchQ[s, _Coupling == _], 
				Message[GetVacuumConditions::unexpct];
				Return@ solution;
			];
			(*Update the expression for the solved-for couplings*)
			If[MatchQ[s[[1]], Coupling[_delta, __]|Bar@ Coupling[_delta, __]],
				partialSol[s[[1, 1, 1]]]+= s[[2]];
			,
				partialSol[s[[1, 1]]]= s[[2]];
			];
		, {s, solution}];
		
		(*Substitute couplings for solutions + perturbative version of the same coupling for next order*)
		substitutions= Replace[solution, {
			(c:Coupling[lab_, inds_, order_]) == expr_ :> 
				{Rule[c, RelabelIndices@ expr + Coupling[delta@ lab, inds, order+ 1]], 
					Rule[c, Coupling[delta@ lab, inds, order+ 1]]},
			(c:Bar@ Coupling[lab_, inds_, order_]) == expr_ :> 
				{Rule[c, RelabelIndices@ expr + Bar@ Coupling[delta@ lab, inds, order+ 1]], 
					Rule[c, Bar@ Coupling[delta@ lab, inds, order+ 1]]}
			}, {1}];
		eqs= eqs//. substitutions[[;;, 1]];
		(*Update variables to ensure that the perturbative deltas are the ones solved for in next iteration*)
		variables= variables/. substitutions[[;;, 2]];
	, {ord, 0, eftOrd- 4}];
	
	(*If no solution was found, return False*)
	Switch[loopResult
	, False, 
		Return@ False;
	, Except[Null],
		Message[GetVacuumConditions::dbg, Iconize@ loopResult];
		Return@ False;
	];
	
	partialSol= Normal@ KeyMap[#[]&, partialSol];
	
	(*Check solution*)
	originalEqs= DeleteCases[originalEqs//. partialSol, {True, _}];
	originalEqs= SeriesEFT[#[[1, 1]], EFTOrder-> #[[2]]]&/@ originalEqs;
	If[!MatchQ[originalEqs, {0...}],
		Message[GetVacuumConditions::inconsistent];
	];
	
	partialSol
];


(* ::Text:: *)
(*Get the system up to the relevant order, keeping only the equations that can be valid up to that order*)


SolveSystemAtOrder[equationSystem_List, variables_List, nonZero_List, order_Integer]:= Module[
		{system, solution},
	(*Expand each equation up to the required order above leading*)
	system= Function[{x}, SeriesEFT[x, EFTOrder-> #[[2]] + order]]/@ #[[1]]&/@ equationSystem;
	
	solution= Reduce[Join[system, nonZero], variables];
	(*We pursue the general solution, so we safely disregard the inequalities*)
	solution/. _Unequal-> True (*Thread@ Rule[nonZero, True]*)
];


(* ::Text:: *)
(*Get the coupling variables in the equation system ordered by how easy they will be to solve for.*)


GetEquationVariables@ equations:{__Equal}:= Module[
		{flavored, linCouplings, nonZero, powerCouplings, productCouplings, terms},
	terms= Map[TermsToList, List@@@ equations, {2}]// Flatten; 
	
	(*We do not currently attempt to solve for flavored couplings*)
	flavored= Cases[terms, Coupling[_, {__}, _], All];
	flavored= DeleteDuplicates@ Flatten[{#, Bar@ #}&/@ flavored];
	
	(*Everywhere else we look for couplings w/o indices*)
	(*Get linear couplings*)
	linCouplings= DeleteDuplicates@ Join[Cases[terms, _Coupling], 
		Cases[terms, Times[c:Coupling[_, {}, _]|Bar@ Coupling[_, {}, _], rest__]/; FreeQ[{rest}, _Coupling]:> c] ];
	
	(*Get power couplings*)
	powerCouplings= Complement[Cases[terms, Power[c:Coupling[_, {}, _]|Bar@ Coupling[_, {}, _], Except[-1]]:> c, All], linCouplings];
	
	(*Get product couplings*)
	productCouplings= Cases[terms, prod_Times/; Count[prod, _Coupling, All] > 1:> prod, All];
	productCouplings= DeleteDuplicates@ Cases[productCouplings, Coupling[_, {}, _], All];
	(*And get their conjugates*)
	productCouplings= DeleteDuplicates@ Flatten[{#, Bar@ #}&/@ productCouplings];
	
	(*Assume that couplings in product terms are non-zero*)
	nonZero= Thread[Complement[productCouplings, linCouplings] != 0];
	
	productCouplings= Complement[productCouplings, Join[powerCouplings, linCouplings]];
	
	{Join[powerCouplings, productCouplings, linCouplings], nonZero}
];


(* ::Subsubsection::Closed:: *)
(*Substitute the field mass symbols*)


(* ::Text:: *)
(*Rewrites the vector and radial scalar mass terms and the vector-Goldstone kinetic mixing terms of a broken-phase Lagrangian in terms of the field mass symbols: diagonal vector/radial-scalar mass terms become M_F^2 and the GB-vector kinetic mixing coefficients become (signed) M_V. *)
(*It is idempotent once the masses are already substituted. *)


SubstituteMassesInLagrangian::massconflict= "The diagonal mass term for `1` is in disagreement with the vector-Goldstone kinetic mixing coefficient.";
SubstituteMassesInLagrangian::massless= "The vector/radial scalar `1` was defined as massless. Mass substitution cannot be completed."


Options@ SubstituteMassesInLagrangian= {"CheckVectorMassConsistency"-> True}


SubstituteMassesInLagrangian[0, OptionsPattern[]]:= 0 (* code below throws errors for lag=0 *)


SubstituteMassesInLagrangian[lag_, OptionsPattern[]]:= Module[
		{gbTerms, gbMassSq, massTerms, relevantTerms, out= HcExpand@ lag},
	{relevantTerms, out}= SeparateRelevantTerms@ out; 
	(*Ensure terms are single coefficients * operator*)
	relevantTerms= HcExpand@ IntroduceEffectiveCouplings@ relevantTerms; 
	
	(*Linear vector masses M_V read off the GB-vector kinetic mixing coefficients*)
	gbTerms= IsolateVectorGBMixingTerms@ relevantTerms;
	massTerms= relevantTerms- gbTerms;
	gbTerms= DeleteCases[TermsToList@ gbTerms, 0, 1]; 
	gbMassSq= Association[#[[2]]-> #[[1, 1]]^2&/@ (GBVectorMassRule/@ gbTerms)];
	
	(*Check that each diagonal vector mass M_V^2 agrees with the squared mixing coefficient*)
	If[OptionValue["CheckVectorMassConsistency"],
		CheckVectorMassConsistency[massTerms, gbMassSq];
	];
	
	out+ Total[GBMixingTermSubstituted/@ gbTerms]
		+ Total[SubstituteMassInTerm/@ TermsToList@ massTerms]
];


(* ::Text:: *)
(*Isolates the relevant terms from the rest of the Lagrangian w/o changing the form/grouping of the other terms*)


SeparateRelevantTerms@ lag_:= Module[{terms},
	terms= If[Head@ lag === Plus, List@@ lag, {lag}];
	Plus@@@ SelectAndDeleteCases[terms, x_/; VectorMassTermQ@ x || VectorGBMixingTermQ@ x || RadialMassTermQ@ x]
]


RadialMassTermQ[term_]:= MatchQ[term/. Bar@ f_Field-> f, 
	(c_.)* Field[s_, Scalar, _, {}]^2/; MemberQ[RadialScalarLabels[], s] && FreeQ[c, Field|FieldStrength]];


VectorMassTermQ[term_]:= MatchQ[term/. Bar@ f_Field-> f, 
	(c_.)* Field[v_, _Vector, _, {}]^2 /; FreeQ[c, Field|FieldStrength]];


VectorGBMixingTermQ[term_]:= MatchQ[term/. Bar@ f_Field-> f, 
	(c_.)* Field[_, Vector@ mu_, _, {mu_}] Field[_, Scalar, _, {}] | (c_.)* Field[_, Vector@ mu_, _, {}] Field[_, Scalar, _, {mu_}]/; FreeQ[c, Field|FieldStrength]];


(* ::Text:: *)
(*Builds the rule replacing the kinetic mixing coefficient of a single vector-GB term by the (signed) linear vector mass symbol, together with the associated vector. *)


GBVectorMassRule[term_]:= Replace[Operator@ term/. Bar@ f:_Field-> f, {
	Times[Operator[Field[_, Scalar, __], Field[vec_, Vector@ a_, _, {a_}]], rest__]:>
		{Rule[Times@ rest, -$FieldAssociation[vec, Mass][]], vec},
	Times[Operator[Field[_, Scalar, _, {a_}], Field[vec_, Vector@ a_, __]], rest__]:>
		{Rule[Times@ rest, $FieldAssociation[vec, Mass][]], vec}
}];


(* ::Text:: *)
(*Rewrites the coefficient of a single GB-vector mixing term as the linear field mass symbol. Each term is processed individually.*)


GBMixingTermSubstituted[term_]:= term/. First@ GBVectorMassRule@ term;


(* ::Text:: *)
(*Rewrites the coefficient of a single diagonal mass term (vector or radial scalar) as the squared field mass symbol. Idempotent once the masses are already substituted.*)


SubstituteMassInTerm[term_]:= Module[{label, coeff},
	label= FirstCase[term, Field[l_, ___]:> l, Missing[], All];
	If[$FieldAssociation[label, Mass] === 0, 
		Message[SubstituteMassesInLagrangian::massless, label];
		Abort[];
	];
	coeff= term/. Field[label, rest___]:> 1;
	term/ coeff* DiagonalMassFactor[label]* $FieldAssociation[label, Mass][]^2
];


(* ::Text:: *)
(*The prefactor relating a field's diagonal mass-term coefficient to its squared mass symbol: 1/2 for self-conjugate fields (from the symmetric bilinear) and an extra -1 for scalars (whose mass terms enter with a relative minus sign).*)


DiagonalMassFactor[label_]:= If[$FieldAssociation[label, SelfConjugate], 1/2, 1]*
	If[$FieldAssociation[label, Type] === Scalar, -1, 1];


(* ::Text:: *)
(*Checks that the diagonal mass of a vector matches the squared coefficient of its GB-vector kinetic mixing term (compared at the relevant EFT order).*)


CheckVectorMassConsistency[massTerms_, gbMassSq_Association]:= Module[{vectorTerms, term},
	(*If GB mixing term not present, the key will not be present in gb Mass association*)
	vectorTerms= Cases[TermsToList@ massTerms/. Bar@ f_Field:> f, c_.* Field[v_, _Vector, __]^2:>
		If[KeyExistsQ[gbMassSq, v], ReplaceEffectiveCouplings@ {v, c, DiagonalMassFactor@ v * gbMassSq@ v}, Nothing] ];
	
	Table[
		If[SeriesEFT[term[[3]], EFTOrder-> GetMaxEFTOrder@ term[[2]]]- term[[2]] =!= 0,
			Message[SubstituteMassesInLagrangian::massconflict, term[[1]]];
			Abort[];
		]
	, {term, vectorTerms}];
];


(* ::Text:: *)
(*The radial (physical) scalar modes of the breaking pattern.*)


RadialScalarLabels[]:= Keys@ DeleteCases[
	$FieldAssociation[[Key/@ BrokenPhaseScalarLabels[]]], KeyValuePattern[GoldstoneBoson-> True]];


(* ::Subsection:: *)
(*Get information from Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Check to see if Lagrangian is in a broken phase *)


(* ::Text:: *)
(*Checks if Lagrangian is consistent with the broken phase Lagrangian after current symmetry breaking pattern is triggered.*)


BrokenPhaseLagQ::missmatchgbs= "The Lagrangian contains `1` GBs (`2`) and `3` massive vectors (`4`). Matching numbers are expected.";
BrokenPhaseLagQ::partialSSB= "The Lagrangian contains some but not all of the massive vectors involved in the current SSB pattern. Didn't find `1` in the Lagrangian.";


BrokenPhaseLagQ@ lag_:= Module[{vectors, gbs},
	(*Vectors appearing explicitly in the Lagrangian: massive*)
	vectors= DeleteDuplicates@ Cases[lag, Field[f_, _Vector, __]:> f, All];
	
	(*The gb-vector matching is now done elsewhere*)
	(*gbs= Intersection[DeleteDuplicates@ Cases[lag, Field[f_, Scalar, __]:> f, All],
		GetFieldsByProperty[GoldstoneBoson-> True]];	
	(*If any GBs are present we are not in unitary gauge: there's a 1:1 of massive vectors to GBs*)
	If[Length@ gbs > 0 && Length@ vectors - Length@ gbs =!= 0,
		Message[BrokenPhaseLagQ::missmatchgbs, Length@ gbs, gbs, Length@ vectors, vectors];
		Abort[];
	];*)
	
	Which[
	Length@ vectors === 0,
		Return@ False;
	,Sort@ vectors === Sort@ HeavyVectorLabels[],
		Return@ True;
	,IntersectingQ[vectors, HeavyVectorLabels[]],
		Message[BrokenPhaseLagQ::partialSSB, Complement[HeavyVectorLabels[], vectors]];
		Abort[];
	,True,
		(*All massive vectors are outside SSB context, i.e., manually defined massive vectors*)
		Return@ False;
	];
];


(* ::Subsubsection::Closed:: *)
(*Determine masses *)


(* ::Text:: *)
(*Get the kinetic-mixing terms between GBs and heavy vectors from a Lagrangian*)


IsolateVectorGBMixingTerms@ lag_:=Module[{terms= TermsToList@ lag},
	Plus@@ Select[terms, FreeQ[#, FieldStrength| Fermion] && 
		!FreeQ[Operator@ #/. Bar-> Identity,
		Operator[Field[_, Scalar, _, {a_| PartialDev@ a_}], Field[_, Vector@ a_, _, {}] ]|
		Operator[Field[_, Scalar, _, {}], Field[_, Vector@ a_, _, {a_| PartialDev@ a_}] ] ]& ]
];


(* ::Text:: *)
(*Works only for ungauged, broken-phase Lagrangians*)


DetermineFieldsForMassRedefinitions@ lag_:= Module[{fields, gbVectorParis, mixingTerms},
	fields= DeleteDuplicates@ Cases[lag, Field[label_, __]|FieldStrength[label_, __]:> label, All];
	mixingTerms= TermsToList@ IsolateVectorGBMixingTerms@ lag;
	gbVectorParis= {FirstCase[#, Field[lab_, _Vector,  __]:> lab, {}, All], 
		FirstCase[#, Field[lab_, Scalar, __]:> lab, {}, All]}&/@ mixingTerms/. Bar-> Identity// DeleteDuplicates;
	gbVectorParis
]


(* ::Section:: *)
(*Broken-phase gauge fixing*)


(* ::Subsection:: *)
(*Group algebra of the broken phase*)


(* ::Text:: *)
(*To compute how the structure constant of the product gauge group decomposes in the broken phase*)


(* ::Subsubsection::Closed:: *)
(*Determine heavy vectors *)


(* ::Text:: *)
(*The heavy vectors of the breaking pattern are*)


HeavyVectorLabels[]:= Module[{fields},
	fields= DeleteDuplicates@ Cases[List@@ $SSBInfo[FieldSubstitutions][[;;, 2]], Field[lab_, __]-> lab, All];
	Intersection[GetFieldsByProperty[Type-> Vector, Mass-> Except[0]], fields]
];


HeavyVectorDOFs[]:= Module[{vecs},
	vecs= HeavyVectorLabels[];
	Flatten[If[$FieldAssociation[#, SelfConjugate], {#}, {#, Conj@ #}]&/@ vecs]
]


(* ::Text:: *)
(*The light vectors are the gauge fields themself *)


LightVectorLabels[]:= Module[{fields},
	fields= DeleteDuplicates@ Cases[List@@ $SSBInfo[FieldSubstitutions][[;;,2]], Field[lab_, __]-> lab, All];
	Intersection[GetFieldsByProperty[Type-> Vector, Mass-> 0], fields]
];


(* ::Text:: *)
(*Returns the block structure of \kappa_{ij} *)


HeavyIndexMetric[]:= Module[{conj, count, vec, dofs},
	dofs= HeavyVectorDOFs[]/. Conj-> Identity;
	{conj, count}= {False, 0};
	Normal@ SparseArray@ Table[
			count++; 
			If[$FieldAssociation[vec, SelfConjugate], 
				{count, count}-> 1
			,
				conj= !conj;
				If[conj,
					{count, count+ 1}-> 1
				,
					{count, count- 1}-> 1
				]
			]
		, {vec, dofs}] 	
];


(* ::Subsubsection::Closed:: *)
(*Rotation matrix L^{i(\[Alpha])}_A on the Lie algebra space  *)


(* ::Text:: *)
(*Returns (L^{-1})^A_i as a pure function from the label of a generic gauge unbroken index, A (taking values in the adjoint of all unbroken groups), and returns an association with the blocks from each broken product group to the heavy vector basis*)


GaugeIndexDecompositionHeavy[]:= Module[{unbrokenGaugeFields, f, mu, inds, out, vec, heavies},
	unbrokenGaugeFields= Through@ Lookup[$GaugeGroups, $SSBInfo[OriginalGroup]][Field];
	heavies= HeavyVectorDOFs[];
	
	(*Construct pure function*)
	out= A|-> Evaluate[Table[(*loop over each original gauge field*)
			inds= Sequence@@ ConstantArray[A, Length@ $FieldAssociation[f, Indices]];
			First@ GaugeGroupByProperty[Field-> f]-> 
				Table[(*Make a list with components for each massive vector DOF*)
					f[mu, inds]/. $SSBInfo[FieldSubstitutions, f]/. 
						If[Head@ vec === Conj, {Bar@ Field[First@ vec, __]-> 1}, 
							{Bar@ Field[vec, __]-> 0, Field[vec, __]-> 1}]/. 
						_Field-> 0
				, {vec, heavies}]
		, {f, unbrokenGaugeFields}]];
	
	Apply[Association]@* out
];


(* ::Text:: *)
(*Similar but for (L^{-1})^A_\[Alpha]*)


GaugeIndexDecompositionLight[]:= Module[{unbrokenGaugeFields, f, mu, inds, out, vec, lights},
	unbrokenGaugeFields= Through@ Lookup[$GaugeGroups, $SSBInfo[OriginalGroup]][Field];
	lights= LightVectorLabels[];
	
	out= A|-> Evaluate[Table[
		inds= Sequence@@ ConstantArray[A, Length@ $FieldAssociation[f, Indices]];
		First@ GaugeGroupByProperty[Field-> f]-> 
			Table[(*Make a list with components for each light vector DOF (all real)*)
				f[mu, inds]/. $SSBInfo[FieldSubstitutions, f]/. 
					Field[vec, __]-> 1/. _Field-> 0
			, {vec, lights}]
	, {f, unbrokenGaugeFields}]];
	
	Apply[Association]@* out
];


(* ::Subsubsection::Closed:: *)
(*Gauge coupling matrix of the symmetric-phase gauge fields*)


(* ::Text:: *)
(*Recover the original symmetric phase gauge couplings (squared) given the vector-field rotations and the ahat matrix.  *)
(*a^{AB} = (L^-1)^A_\alpha  (L^-1)^B_\beta \hat{a}^{\alpha\beta} + (L^-1)^A_i (L^-1)^B_j \kappa^{ij}*)


SymmetricPhaseGaugeCouplings[ahat_]:= Module[
		{heavyDeltas, heavyIndex, lightDeltas, lightIndex, kappa, originalGroups, a, gr, A, B},
	heavyIndex= GaugeIndexDecompositionHeavy[];
	heavyDeltas= KroneckerDeltaOfDOF[#][A, B]&/@ HeavyVectorDOFs[]// DiagonalMatrix;
	lightIndex= GaugeIndexDecompositionLight[];
	lightDeltas= KroneckerDeltaOfDOF[#][A, B]&/@ LightVectorLabels[]// DiagonalMatrix;
	kappa= HeavyIndexMetric[];
	originalGroups= $SSBInfo@ OriginalGroup;
	
	{A, B}|-> Evaluate[Table[
			lightIndex[A][gA] . lightDeltas . ahat . lightIndex[B][gB]
			+ heavyIndex[A][gA] . heavyDeltas . kappa . heavyIndex[B][gB]
		, {gA, originalGroups}, {gB, originalGroups}]]
];


(* ::Text:: *)
(*Provides a Kronecker delta in the representation of a field DOF*)


KroneckerDeltaOfDOF[dof_]:= Module[{label, conj},
	{label, conj}= {dof/. Conj-> Identity, Head@ dof === Conj};
	{A, B}|-> Evaluate@ Product[
			If[KeyExistsQ[$FlavorIndices, ind],
				Delta[ind][A, B]
			,
				If[Xor[conj, Head@ ind === Bar], Bar, Identity]@ 
					del[ind/. Bar-> Identity][A, B]
			]
		, {ind, $FieldAssociation[label, Indices]}]
];


(* ::Text:: *)
(*Inverts the coupling matrix perturbatively in EFT order, assuming the LO matrix is diagonal (or with diagonal pairs in the new basis for the adjoint representation). *)


InvertGaugeCouplingMatrix::inconsistentLO= "The LO coupling of `1` is not compatible with a gauge a single gauge coupling: `2`";
InvertGaugeCouplingMatrix::inconsistentOrd= "Couldn't determine perturbation order.";
InvertGaugeCouplingMatrix::nonDiag= "The LO coupling matrix is non-diagonal; `1`. Cannot invert it with available methods.";
InvertGaugeCouplingMatrix::zeroLO= "The leading order gauge coupling associated with `1` was found to be 0";


InvertGaugeCouplingMatrix[couplingMatrixFunc_Function, eftOrd_Integer]:= Module[
		{coupling, delOrder, gr, inds, leading, invCouplingMatrix, invLeading, originalGroups, out, 
			perturbation, symGaugeCouplings, test, A, B, n},
	originalGroups= $SSBInfo@ OriginalGroup;
	
	(*Determine the LO part of the coupling matrix*)
	leading= Map[SeriesEFT[#, EFTOrder-> 0]&, couplingMatrixFunc[A, B], {2}];
	If[!DiagonalMatrixQ@ leading, 
		Message[InvertGaugeCouplingMatrix::nonDiag, leading];
		Abort[];
	];
	(*The leading order matrix should simply be the symmetric phase gauge-coupling matrix, but in the 
	broken-phase representation decomposition*)
	(*Project out the gauge couplings*)
	symGaugeCouplings= Table[
			gr= originalGroups[[n]];
			coupling= leading[[n, n]]* If[$GaugeGroups[gr, Abelian], 1, 
				DecomposeCGs@ del[gr@ adj][B, A]/$Representations[gr@ adj, Dimension]]// 
				Expand// Contract// TrigSimplify;
			If[coupling === 0,
				Message[InvertGaugeCouplingMatrix::zeroLO, gr];
				Abort[];
			];
			(*Check that the LO result is consistent with a single gauge coupling in the new representation basis*)
			test= leading[[n, n]]- coupling* If[$GaugeGroups[gr, Abelian], 1, DecomposeCGs@ del[gr@ adj][B, A]];
			test= Collect[test, _BrokenInds, TrigSimplify];
			If[test =!= 0,
				Message[InvertGaugeCouplingMatrix::inconsistentLO, gr, leading[[n, n]]];
				Abort[];
			];
			coupling
		, {n, Length@ leading}];
	
	(*Determine the inverse matrix perturbatively*)
	perturbation= {a, b}|-> Evaluate[couplingMatrixFunc[a, b]- leading/. {A-> a, B-> b}// Expand];
	delOrder= OperatorDimension@ perturbation[A, B];
	If[delOrder === 0, 
		Message[InvertGaugeCouplingMatrix::inconsistentOrd];
		Abort[];
	];
	
	invLeading= {A, B}|-> Evaluate@ DiagonalMatrix@ Table[
			gr= originalGroups[[n]];
			symGaugeCouplings[[n]]^-1* If[$GaugeGroups[gr, Abelian], 1, DecomposeCGs@ del[gr@ adj][B, A]]
		, {n, Length@ leading}];
	
	invCouplingMatrix= Sum[
			inds= Join[{A}, Table[Unique[], 2n], {B}];
			out= (-1)^n* invLeading@@ inds[[;;2]];
			Do[
				out= out . perturbation@@ inds[[2k;; 2k+1]] . invLeading@@ inds[[2k+1;; 2k+2]];
				out= Map[SeriesEFT[#, EFTOrder-> eftOrd]&, out, {2}];
			, {k, n}];
			out
		, {n, 0, Floor[eftOrd/delOrder]}];
	invCouplingMatrix= Map[Collect[#, _BrokenInds, TrigSimplify]&, invCouplingMatrix, {2}];
	
	(*Return as a function*)
	{a, b}|-> Evaluate[invCouplingMatrix/. {A-> a, B-> b}]
];


(* ::Subsubsection::Closed:: *)
(*Broken-phase structure constants*)


(* ::Text:: *)
(*Gives the structure constants restricted to the space and basis of the heavy vectors, f^{i}_{jk}. It is returned as a pure function of the three index labels, and returns a rank-3 tensor in the space of massive vector DOFs (arranged according to HeavyVectorDOFs). *)
(*- The conversion from original to massive vector indices (f_{ABC} -> f_{ijk}) is done purely with lower indices (using L^-1). The first index is raised with a kappa afterwards. *)
(*- To lower the indices on the original structure constants use f_{ABC} = a^{-1}_{AD} f^D_{BC}*)


StructureConstantHHH[aInvMatrixFunc_Function, eftOrder_]:= Module[{d, gr, hi, groups, n, m, p},
	hi= GaugeIndexDecompositionHeavy[];
	groups= $SSBInfo@ OriginalGroup;
	(*Tensor multiplying the heavy vector metric to raise the first index.*)
	{a, b, c}|-> Evaluate@ Map[TrigSimplify@ SeriesEFT[#, EFTOrder-> eftOrder]&, TensorContract[
		TensorProduct[Expand[HeavyIndexMetric[] . Transpose[List@@ hi[a]] . aInvMatrixFunc[a, d]], 
			Table[
				gr= groups[[m]]; 
				If[m === n === o && !$GaugeGroups[gr, Abelian],
					DecomposeCGs@ CG[fStruct@ gr, {d, b, c}]
				,
					0
				]
			, {m, Length@ groups}, {n, Length@ groups}, {o, Length@ groups}],
		List@@ hi[b], List@@ hi[c]], {{2, 3}, {4, 6}, {5, 8}}], {3}]
]


(* ::Text:: *)
(*Returns the structure constants restricted to the space and basis of the heavy vectors, f^{i}_{j \[Alpha]}.*)
(*Same as above*)


StructureConstantHHL[aInvMatrixFunc_Function, eftOrder_]:= Module[
		{d, gr, groups, heavyIndex, lightIndex, n, m, p},
	heavyIndex= GaugeIndexDecompositionHeavy[];
	groups= $SSBInfo@ OriginalGroup;
	lightIndex= GaugeIndexDecompositionLight[];
	(*Tensor multiplying the heavy vector metric to raise the first index.*)
	{a, b, c}|-> Evaluate@ Map[TrigSimplify@ SeriesEFT[#, EFTOrder-> eftOrder]&, TensorContract[
		TensorProduct[Expand[HeavyIndexMetric[] . Transpose[List@@ heavyIndex[a]] . aInvMatrixFunc[a, d]], 
			Table[
				gr= groups[[m]]; 
				If[m === n === o && !$GaugeGroups[gr, Abelian],
					DecomposeCGs@ CG[fStruct@ gr, {d, b, c}]
				,
					0
				]
			, {m, Length@ groups}, {n, Length@ groups}, {o, Length@ groups}],
		List@@ heavyIndex[b], List@@ lightIndex[c]], {{2, 3}, {4, 6}, {5, 8}}], {3}]
]


(* ::Text:: *)
(*Returns the structure constants restricted to the space and basis of the heavy vectors, f^{\[Alpha]}_{i j}.*)
(*Same as above*)


StructureConstantLHH[aInvMatrixFunc_Function, ahat_, eftOrder_]:= Module[
		{d, gr, groups, heavyIndex, lightIndex, n, m, p},
	heavyIndex= GaugeIndexDecompositionHeavy[];
	groups= $SSBInfo@ OriginalGroup;
	lightIndex= GaugeIndexDecompositionLight[];
	(*Tensor multiplying the heavy vector metric to raise the first index.*)
	{a, b, c}|-> Evaluate@ Map[TrigSimplify@ SeriesEFT[#, EFTOrder-> eftOrder]&, TensorContract[
		TensorProduct[Expand[ahat . Transpose[List@@ lightIndex[a]] . aInvMatrixFunc[a, d]], 
			Table[
				gr= groups[[m]]; 
				If[m === n === o && !$GaugeGroups[gr, Abelian],
					DecomposeCGs@ CG[fStruct@ gr, {d, b, c}]
				,
					0
				]
			, {m, Length@ groups}, {n, Length@ groups}, {o, Length@ groups}],
		List@@ heavyIndex[b], List@@ heavyIndex[c]], {{2, 3}, {4, 6}, {5, 8}}], {3}]
]


(* ::Subsection:: *)
(*Scalar fields and VEVs*)


(* ::Subsubsection::Closed:: *)
(*Determine scalar multiplets receiving a VEV *)


(* ::Text:: *)
(*A list of all the scalar fields receiving a VEV \[CurlyPhi]^{\prime a}. Other scalars are irrelevant for gauge fixing.  *)


SymmetryBreakingScalars[]:= Module[{scalars, f, vev},
	scalars= Intersection[GetFieldsByProperty[Type-> Scalar], Keys@ $SSBInfo[FieldSubstitutions]];
	(*Select scalars with a non-zero VEV*)
	Table[
		vev= $SSBInfo[FieldSubstitutions, f][[2]]/. _Field-> 0;
		If[vev === 0, Nothing, f]
	, {f, scalars}]
];


SymmetryBreakingScalarDOFs[]:= Module[{scals},
	scals= SymmetryBreakingScalars[];
	Flatten[If[$FieldAssociation[#, SelfConjugate], {#}, {#, Conj@ #}]&/@ scals]
];


(* ::Text:: *)
(*VEVs of all the scalars developing a VEV (v^a)*)


ScalarVEVs[]:= Module[{scalars, f, inds, out, vev},
	scalars= SymmetryBreakingScalars[];
	out= a|-> Evaluate[Table[
		inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
		f-> DecomposeFields[f@ inds]/. _Field-> 0
	, {f, scalars}]];
	Apply[Association]@* out
];


(* ::Text:: *)
(*All the scalar degrees of freedom embedded in those fields (VEV-subtracted)*)


ScalarFields[]:= Module[{scalars, f, inds, out, vev, fields},
	scalars= SymmetryBreakingScalars[];
	out= a|-> Evaluate[Table[
		inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
		(*Subtract of the VEV*)
		fields= DecomposeFields[f@ inds];
		f-> Expand[fields- (fields/. _Field-> 0)]
	, {f, scalars}]];
	Apply[Association]@* out
];


(* ::Text:: *)
(*Returns the block structure of h_{ab} for the collected scalar fields*)


SymPhaseCanonicalScalarMetric[]:= Module[{conj, count, vec, dofs},
	dofs= SymmetryBreakingScalarDOFs[]/. Conj-> Identity;
	{conj, count}= {False, 0};
	Normal@ SparseArray@ Table[
			count++; 
			If[$FieldAssociation[vec, SelfConjugate], 
				{count, count}-> 1
			,
				conj= !conj;
				If[conj,
					{count, count+ 1}-> 1
				,
					{count, count- 1}-> 1
				]
			]
		, {vec, dofs}] 	
];


(* ::Subsubsection::Closed:: *)
(*Background(VEV)-dependent scalar metric*)


(* ::Text:: *)
(*The mass-basis broken-phase scalars (radial modes and would-be Goldstones) that the symmetry-breaking scalars decompose into.*)


BrokenPhaseScalarLabels[]:= Module[{fields},
	fields= DeleteDuplicates@ Cases[
			KeyTake[$SSBInfo[FieldSubstitutions], SymmetryBreakingScalars[]][[;;, 2]],
			Field[lab_, __]-> lab, All];
	Intersection[GetFieldsByProperty[Type-> Scalar], fields]
];


BrokenPhaseScalarDOFs[]:= Module[{scals},
	scals= BrokenPhaseScalarLabels[];
	Flatten[If[$FieldAssociation[#, SelfConjugate], {#}, {#, Conj@ #}]&/@ scals]
];


(* ::Text:: *)
(*The canonical metric \[Sigma]_{ab} (1's on- or slightly off-diagonal) over the broken-phase mass basis. Same block structure as ScalarMetric and HeavyIndexMetric but over the broken-phase DOFs.*)


BrokenScalarMetric[]:= Module[{conj, count, scalar, dofs},
	dofs= BrokenPhaseScalarDOFs[]/. Conj-> Identity;
	{conj, count}= {False, 0};
	Normal@ SparseArray@ Table[
			count++;
			If[$FieldAssociation[scalar, SelfConjugate],
				{count, count}-> 1
			,
				conj= !conj;
				If[conj,
					{count, count+ 1}-> 1
				,
					{count, count- 1}-> 1
				]
			]
		, {scalar, dofs}]
];


(* ::Text:: *)
(*The rotation matrix K^a_b relating the original symmetry-breaking scalars \[CurlyPhi]^a= K^a_b \[Phi]^b to the mass-basis broken-phase scalars \[Phi]^b. Returned as a pure function of a generic index, yielding the matrix with rows ordered by SymmetryBreakingScalarDOFs and columns by BrokenPhaseScalarDOFs. K is field-independent (the EFT order resides in its light/heavy mass-ratio couplings).*)


ScalarIndexDecomposition[]:= Module[
		{scalarDOFs, brokenDOFs, label, expr, scalar, scalarFields, brokenScalar},
	scalarDOFs= SymmetryBreakingScalarDOFs[];
	scalarFields= ScalarFields[]; (*VEV-subtracted scalar field functions decomposed in broken phase components*)
	brokenDOFs= BrokenPhaseScalarDOFs[];
	a|-> Evaluate@ Table[
			label= scalar/. Conj-> Identity;
			expr= If[Head@ scalar === Conj, Bar, Identity]@ scalarFields[a][label];
			Table[(*Pick out the coefficient of each broken-phase scalar DOF*)
				expr/. If[Head@ brokenScalar === Conj, {Bar@ Field[First@ brokenScalar, __]-> 1},
						{Bar@ Field[brokenScalar, __]-> 0, Field[brokenScalar, __]-> 1}]/.
					_Field-> 0
			, {brokenScalar, brokenDOFs}]
		, {scalar, scalarDOFs}]
];


(* ::Text:: *)
(*The inverse background metric \overline{h}^{-1}= K \[Sigma] K^T over the original symmetry-breaking scalar DOFs, returned as a pure function of two generic indices. Analogous to SymmetricPhaseGaugeCouplings.*)


InverseSymPhaseScalarMetric[]:= Module[{kDecomp, scalarDeltas, sigma, a, b},
	kDecomp= ScalarIndexDecomposition[];
	scalarDeltas= KroneckerDeltaOfDOF[#][a, b]&/@ BrokenPhaseScalarDOFs[]// DiagonalMatrix;
	sigma= BrokenScalarMetric[];
	{a, b}|-> Evaluate[kDecomp[a] . scalarDeltas . sigma . Transpose@ kDecomp[b]]
];


(* ::Text:: *)
(*Inverts the inverse background metric perturbatively in EFT order to give \overline{h}_{ab}, using that \overline{h} = h(0) at LO. Analogous to InvertGaugeCouplingMatrix, but the LO metric is block, should be the canonical scalar metric (self inverse).*)


InvertScalarMetric::leadingord= "The leading order scalar metric `1` doesn't agree with the canonical form.";
InvertScalarMetric::inconsistentOrd= "Couldn't determine perturbation order.";


InvertScalarMetric[invMetricFunc_Function, eftOrd_Integer]:= Module[
		{delOrder, dims, dofs, h0Metric, inds, invMetric, leading, out,
			perturbation, A, B, n, m, k},

	(*The canonical h0 metric, is its own inverse*)
	h0Metric= {a, b}|-> Evaluate[
			DiagonalMatrix[DecomposeCGs@ KroneckerDeltaOfDOF[#][A, B]&/@ SymmetryBreakingScalarDOFs[]] .
				SymPhaseCanonicalScalarMetric[]/. {A-> a, B-> b}
		];
	
	(*Check that the leading order inverse metric agrees with h0*)
	leading= Map[SeriesEFT[#, EFTOrder-> 0]&, invMetricFunc[A, B], {2}];
	If[leading =!= h0Metric[A, B],
		Message[InvertScalarMetric::leadingord, leading];
		Abort[];
	];

	(*Determine the inverse matrix perturbatively*)
	perturbation= {a, b}|-> Evaluate[invMetricFunc[a, b]- h0Metric[a, b]// Expand];
	delOrder= OperatorDimension@ perturbation[A, B];

	invMetric= Sum[
			inds= Join[{A}, Table[Unique[], 2n], {B}];
			out= (-1)^n* h0Metric@@ inds[[;;2]];
			Do[
				out= out . perturbation@@ inds[[2k;; 2k+1]] . h0Metric@@ inds[[2k+1;; 2k+2]];
				out= Map[SeriesEFT[#, EFTOrder-> eftOrd]&, out, {2}];
			, {k, n}];
			out
		, {n, 0, If[delOrder === 0, 0, Floor[eftOrd/delOrder]]}];
	invMetric= Map[Collect[#, _BrokenInds, TrigSimplify]&, invMetric, {2}];

	(*Return as a function*)
	{a, b}|-> Evaluate[invMetric/. {A-> a, B-> b}]
];


(* ::Subsubsection::Closed:: *)
(*Broken generators *)


(* ::Text:: *)
(*Returns the broken generator x^a_ {ib} as a function of dummy indices into a rank-3 tensor (heavy vector dofs (i), scalar dofs (a), conj scalar dofs (b)).*)


BrokenGenerators[]:= Module[{fields, nullVec, m, n},
	fields= SymmetryBreakingScalarDOFs[];
	nullVec= ConstantArray[0, Length@ HeavyVectorDOFs[]];
	{i, a, b}|-> Evaluate@ Transpose[Table[
			If[n === m, BrokenGeneratorsOnField[fields[[n]]][i, a, b], nullVec]
		, {m, Length@ fields}, {n, Length@ fields}], {2, 3, 1}]
]


(* ::Text:: *)
(*Returns the broken generators acting on the representation of an original scalar field (charged under the original group) in terms of broken index projectors in scalar space and a vector in the space of heavy vector components. *)


BrokenGeneratorsOnField[field_]:= Module[{charges, rep, reps, indNo, indI, indJ, n, k},
	(*Adjust the index representations and the charges depending on conjugation of field*)
	{reps, charges}= If[Head@ field === Conj, 
			{Bar/@ $FieldAssociation[First@ field, Indices], MapAt[Minus, $FieldAssociation[First@ field, Charges], {All, 1}]}
		,
			{$FieldAssociation[field, Indices], $FieldAssociation[field, Charges]}
		];
	indNo= Length@ reps;
	
	{A, i, j}|-> Evaluate[
		{indI, indJ}= {ConstantArray[i, indNo], ConstantArray[j, indNo]};
		(*Non-Abelian*)
		Expand@ Sum[
			Product[
				If[k === n,
					BrokenGeneratorOnRep[reps[[k]]][A, indI[[k]], indJ[[k]]]
				,
					DecomposeCGs@ If[Head@ reps[[k]] === Bar, 
						rep= First@ reps[[k]];
						Bar
					, 
						rep = reps[[k]];
						Identity
					]@ CG[del@ rep, {Index[indI[[k]], rep], Bar@ Index[indJ[[k]], rep]}]
				]
			, {k, indNo}]
		, {n, indNo}] + 
		(*Abelian*)
		Expand[
			Sum[
				BrokenGeneratorOnCharge[ch]
			, {ch, charges}] *
			Product[
				DecomposeCGs@ If[Head@ reps[[k]] === Bar, 
					rep= First@ reps[[k]];
					Bar
				, 
					rep = reps[[k]];
					Identity
				]@ CG[del@ rep, {Index[indI[[k]], rep], Bar@ Index[indJ[[k]], rep]}]
			, {k, indNo}]
		]
	]
]


BrokenGeneratorOnCharge@ group_@ q_:= Module[{A},
	q GaugeIndexDecompositionHeavy[][A][group]
]


BrokenGeneratorOnRep[representation_]:= Module[{rep, group, conj},
	rep= If[(conj= Head@ representation === Bar), First@ representation, representation];
	group= GroupFromRep@ rep;
	{A, i, j}|-> Evaluate@ Table[
		hA DecomposeCGs@ If[conj, Minus@* Bar, Identity]@ CG[gen@ rep, {A, i, j}]//Expand
	, {hA, GaugeIndexDecompositionHeavy[][A][group]}]
]


(* ::Text:: *)
(*The decay constant defined by f_i^a = - i x^a_{ib} v^b*)


DecayConstantMatrix[]:= Module[{scalar, vevs, b},
	vevs= ScalarVEVs[][b];
	{i, a}|-> Evaluate@ Transpose@ Table[
			-I* BrokenGeneratorsOnField[scalar][i, a, b]* 
				If[Head@ scalar === Conj, Bar@ vevs[scalar[[1]]], vevs[scalar]]// Expand
		, {scalar, SymmetryBreakingScalarDOFs[]}]
]


(* ::Subsection:: *)
(*Construct gauge-fixing terms*)


(* ::Subsubsection::Closed:: *)
(*Unitary gauge*)


(* ::Text:: *)
(*Unitary gauge removes all GoldstoneBosons and fixes the remnant gauge-fields as per the R\[Xi]. *)
(*N.b. Assumes broken phase of lag.*)


UnitaryGauge[lag_, eftOrd_Integer]:= Module[{gaugeAssoc, gaugeGFterms, replacement},
	(*Gauge fixes the gauge fields of the stability group as per R\[Xi]*)
	gaugeAssoc= GaugeStatesAssociation[lag, True];
	gaugeGFterms= RXiRemnant[gaugeAssoc];
	gaugeGFterms= ReplaceInLagrangian[gaugeGFterms, $SSBCouplingReplacements, Simplify-> False];
	(*Apply series expansion to expression as perscribed by option in GaugeFixLagrangian*)
	gaugeGFterms= SeriesEFT[gaugeGFterms, EFTOrder-> eftOrd];
	gaugeGFterms= InternalCollectOperators[gaugeGFterms, TrigSimplify-> True, NormalForm-> True];
	
	(*Removes all gbs*)
	replacement= Thread@ Rule[Values[gaugeAssoc@ GoldstoneBosons]/. Conj-> Identity, 0];
	gaugeGFterms+ (lag/. replacement)
];


(* ::Subsubsection::Closed:: *)
(*Ordinary R\[Xi] gauge *)


(* ::Text:: *)
(*The ordinary R\[Xi] gauge (not the BF gauge).*)
(*We have redefined anti-ghost -> I* anti-ghost to ensure that the field does not get a sign under Hermitian conjugation.*)
(*N.b. Assumes broken phase of lag.*)


BrokenPhaseRXiTerms[lag_, eftOrder_]:= Module[
		{aHatMatrix, aMatrix, fdecay, fhhl, fhhh, fields, flhh, gaugeAssoc, gaugeGFterms, gfTerms, ghostTerms, 
		hMet, kappa, massMatrix, out, xgen, xih, xihInv, xilInv, mu, nu, a, b, c, i, j, k},
	gaugeAssoc= GaugeStatesAssociation[lag, True];
	fields= CollectiveFields@ gaugeAssoc@ GoldstoneBosons;
	
	kappa= HeavyIndexMetric[];
	(*hMet= SymPhaseCanonicalScalarMetric[];*)
	hMet= InvertScalarMetric[InverseSymPhaseScalarMetric[], eftOrder- 4];
	massMatrix[]:= RelabelIndices[gaugeAssoc[VectorMasses], Unique-> True]; 
	
	(*Determine group sructure based on the global SSB information*)
	aHatMatrix= InvertMetric@ gaugeAssoc@ InverseAHat;
	aMatrix= InvertGaugeCouplingMatrix[SymmetricPhaseGaugeCouplings@ aHatMatrix, eftOrder- 4]; (*original couplings*)
	fhhh= StructureConstantHHH[aMatrix, eftOrder- 4];
	fhhl= StructureConstantHHL[aMatrix, eftOrder- 4];
	flhh= StructureConstantLHH[aMatrix, aHatMatrix, eftOrder- 4];
	fdecay= DecayConstantMatrix[];
	xgen= BrokenGenerators[];
	
	(*Gauge parameter matrices of the heavy vectors: \xi_i^{1(-1)}*)
	xih= \[Xi]Gauge/@ HeavyVectorDOFs[]/. Conj-> Identity;
	
	xihInv= DiagonalMatrix[Power[#, -1]&/@ xih];
	xih= DiagonalMatrix@ xih;
	
	(*gf for the non SSB-participating gauge fields*)
	gaugeGFterms= RXiRemnant[gaugeAssoc];
	
	(*Massive vectors fields:*)
	(*(\partial_\mu V_\mu^i)^2*)
	gfTerms= -1/2* PD[mu, fields[MassiveVectors][mu, a]] . xihInv . kappa . PD[nu, fields[MassiveVectors][nu, a]];
	(*M_i \partial_\mu V_\mu^i \chi_i*)
	gfTerms+= PD[mu, fields[MassiveVectors][mu, a]] . kappa . massMatrix[] . fields[GoldstoneBosons][a];
	(*M_i^2 \chi^2_i*)
	gfTerms+= -1/2 fields[GoldstoneBosons][a] . xih . kappa . massMatrix[] . massMatrix[] . fields[GoldstoneBosons][a];
	
	(*Massive ghosts, kinetic*)
	ghostTerms= -I Inner[GrassmannProduct, fields[MassiveGhostBar][a], PD[mu, CD[mu, fields[MassiveGhost][a] ] ] ]+
		-I Inner[GrassmannProduct, fields[MassiveGhostBar][a], xih . massMatrix[] . massMatrix[] . fields[MassiveGhost][a]];

	(*Ghost-vector interactions:*)
	(*\parital_\mu \bar{c}_a V_\mu^i f^{a}_{ij} u^j*)
	ghostTerms+= I Inner[GrassmannProduct, PD[mu, fields[GhostBar][a]], TensorContract[TensorProduct[
		fields[MassiveVectors][mu, i], flhh[a, i, j], fields[MassiveGhost][j]], {{1, 3}, {4, 5}}]];
	(*\parital_\mu \bar{u}_i V_\mu^j f^{i}_{ja} c^a*)
	ghostTerms+= I Inner[GrassmannProduct, PD[mu, fields[MassiveGhostBar][i]], TensorContract[TensorProduct[
		fields[MassiveVectors][mu, j], fhhl[i, j, a], fields[Ghost][a]], {{1, 3}, {4, 5}}]];
	(*\parital_\mu \bar{u}_i V_\mu^j f^{i}_{jk} u^k*)
	ghostTerms+= I Inner[GrassmannProduct, PD[mu, fields[MassiveGhostBar][i]], TensorContract[TensorProduct[
		fields[MassiveVectors][mu, j], fhhh[i, j, k], fields[MassiveGhost][k]], {{1, 3}, {4, 5}}]];

	(*Ghost-scalar interactions:*)
	(*M_j \chi^j \bar{u}_i f^i_{j a} c^a*)
	ghostTerms+= I Inner[GrassmannProduct, xih . fields[MassiveGhostBar][i], TensorContract[TensorProduct[
		massMatrix[] . fields[GoldstoneBosons][j], fhhl[i, j, a], fields[Ghost][a]], {{1, 3}, {4, 5}}]];
	(*\bar{u}_i f^i_a x^a_{jb} \varphi^b u^j*)
	ghostTerms+= -Inner[GrassmannProduct, fields[MassiveGhostBar][i] . xih . kappa . fdecay[i, a] . hMet[a, b], 
		TensorContract[TensorProduct[xgen[j, b, c], fields[Scalars][c], fields[MassiveGhost][j]], {{3, 4}, {1, 5}}]];
		
	out= gaugeGFterms+ gfTerms+ ghostTerms;
	out= ContractDelta@ ReplaceInLagrangian[out, $SSBCouplingReplacements, Simplify-> False];
	(*Apply series expansion to expression as perscribed by option in GaugeFixLagrangian*)
	out= SeriesEFT[out, EFTOrder-> eftOrder];
	InternalCollectOperators[out, TrigSimplify-> True, NormalForm-> True]
];


(* ::Subsubsection::Closed:: *)
(*R\[Xi] gauge-fixing for the unbroken gauge fields  *)


(* ::Text:: *)
(*Gives the ordinary R\[Xi] terms for the unbroken (and the non-participating) gauge fields. *)


RXiRemnant[stateAssoc_]:= Module[{fields, xilInv, out, mu, nu, a},	
	(*gf for the non SSB-participating gauge fields*)
	out= Total@ KeyValueMap[Function[{fs, met}, IndividualRXiGaugeFixingTerms[fs[[1]], met]], 
		stateAssoc@ NonParticipatingGaugeFields];
	
	(*Stability group fields of the symmetry breaking*)
	fields= CollectiveFields@ stateAssoc@ GoldstoneBosons;
	xilInv= DiagonalMatrix[Power[#, -1]&/@ \[Xi]Gauge/@ LightVectorLabels[]];
	
	(*Unbroken gauge fields: (\partial_\mu B_\mu^\alpha)^2 *)
	out+= -1/2* PD[mu, fields[Vectors][mu, a]] . xilInv . stateAssoc@ InverseAHat . PD[nu, fields[Vectors][nu, a]];
	
	(*Unbroken ghosts kinetic: \bar{c}_\alpha \partial_\mu d^\mu c^\alpha*)
	out+= -I Inner[GrassmannProduct, fields[GhostBar][a], PD[mu, CD[mu, fields[Ghost][a] ] ] ];
	
	out
];


(* ::Subsubsection::Closed:: *)
(*Find Goldstone bosons and decay constants*)


(* ::Text:: *)
(*We determine the GBs corresponding to the heavy vectors in the broken phase along with their associated decay constants.*)


GaugeStatesAssociation[lag_, brokenPhase_? BooleanQ]:= Module[{gbs, invaHat, vecMasses, nonParticipating, terms, vectors},
	(*Coupling matrices of the (remaining) gauge fields*)
	{invaHat, nonParticipating}= InverseaHatMatrix[lag, brokenPhase];
	
	(*Shortcut if in symmetric phase - no massive vectors or GB*)
	If[!brokenPhase,
		Return@ <|NonParticipatingGaugeFields-> nonParticipating|>;
	];
	
	(*Get the massive vectors and GBs with which they mix*)
	terms= GBMixingTerms@ lag;
	vectors= DeleteDuplicates@ Cases[terms, Field[f_, _Vector, __]:> f, All];
	{gbs, vecMasses}= If[Length@ vectors > 0, 
			Association@@@ Transpose[Thread[#-> GBAndDecayConstant[#, terms]]&/@ vectors]
		,
			{<||>, <||>}
		];
		
	(*Construct diagonal mass matrix: M_i*)
	vecMasses= DiagonalMatrix@ Flatten@ Table[
			If[$FieldAssociation[vec, SelfConjugate], {1}, {1, 1}]* vecMasses@ vec
		, {vec,  HeavyVectorLabels[]}];
	vecMasses= vecMasses/. Missing["KeyAbsent",_]-> 0; (* work-around if HeavyVectorLabels[] contains fields not present in  lag *)
	
	<|
		GoldstoneBosons-> gbs,
		InverseAHat-> invaHat,
		NonParticipatingGaugeFields-> nonParticipating,
		VectorMasses-> vecMasses
	|>
];


(* ::Text:: *)
(*Isolate the Lagrangian terms w/ kinetic mixing between massive vectors and GBs*)


GBMixingTerms@ expr_:= Module[{vectors, terms},
	vectors= DeleteDuplicates@ Cases[expr, Field[f_, _Vector, __]:> f, All];
	terms= InternalCollectOperators[expr, InternalOpRepresentation-> True];
	terms= terms/. AtomicOp[Except@{{OrderlessPatternSequence[Alternatives@@ Join[vectors, Conj/@ vectors], _]}, 1}, __]-> 0;
	AtomicToNormalForm@ terms
];


(* ::Text:: *)
(*Search for each vector*)


GBAndDecayConstant::mult= "Multiple potential GB candidates for `1`";


GBAndDecayConstant[vector_, expr_]:= Module[{a, gb, inds, term, vec, vecMass},
	inds= ConstantArray[a, Length@ $FieldAssociation[vector, Indices] +1];
	term= VarD[RelabelIndices[expr, Unique->True], vector@@ inds];
	term= Collect[term, _Field|Bar@ _Field];
	
	(* Check if multiple potential GB candidates to the Vector *)
	If[Count[term, _Field, All] > 1,
		Message[GBAndDecayConstant::mult, vector];
		Abort[];
	];
	
	(*Determine GBs. N.b. V^i \kappa \chi^i, so the GBs associated w/ V^i is the conjugate of the one appearing in the mixing term*)
	gb= FirstCase[term, Bar@ Field[f_, __]-> f, 
		FirstCase[term, Field[f_, __]-> Conj@ f, All], All];

	(*Determine masses of the heavy vectors*)
	vecMass= term/. _Field-> 1;

	{gb, vecMass}
];


(* ::Text:: *)
(*Function to obtain the diagonal of the unbroken gauge coupling matrix \hat{a}_{\[Alpha] \[Beta]}^{-1}.*)
(*Gauge fields that are not involved in the breaking (nonParticipating) simply get an association of inverse couplings.*)
(*Non-broken Lagrangians only have non-participating gauge fields.*)


InverseaHatMatrix::overlap= "Some of the kinetic terms mix gauge fields from unbroken gauge fields from the SSB with gauge fields that do not participate in SSB in the first place.";
InverseaHatMatrix::nonpresentgaugefields= "The remnant gauge fields `1` of the SSB pattern are not present in the Lagrangian.";


InverseaHatMatrix[lag_, brokenPhase_:True]:= Module[
		{ahatAssoc, fields, nonParticipating, participating, stabilityFields, terms, f1, f2},
	(*Select the gauge kinetic terms and read of their coefficients*)
	terms= Replace[
		TermsToList@ SelectOperatorClass[lag/. _LCTensor->0 (*Remove topological gauge terms*), {}, 4], 
		{
			Times[FieldStrength[lab_, __]^2 , coef__]:> {{lab, lab}, Times@ coef},
			Times[FieldStrength[lab1_, __], FieldStrength[lab2_, __], coef__]:> {Sort@ {lab1, lab2}, Times@ coef}
		}, {1}];
	fields= DeleteDuplicates@ Flatten@ terms[[;;, 1]];
	
	(*Stability fields (remnant gauge fields) of the breaking*)
	stabilityFields= If[brokenPhase, 
			If[!SubsetQ[fields, LightVectorLabels[]],
				Message[InverseaHatMatrix::nonpresentgaugefields, Complement[LightVectorLabels[], fields]];
				Abort[];
			];
			LightVectorLabels[]
		, 
			{}
		];
	
	(*Normalize the couplings and group by fields*)
	ahatAssoc= -4* Plus@@ #&/@ GroupBy[terms, First-> Last];
	
	(*Gauge fields and aHat that where not involved in the SSB*)
	nonParticipating= Complement[fields, stabilityFields];
	nonParticipating= KeySelect[ahatAssoc, SubsetQ[nonParticipating, #]&];
	
	(*Gauge fields and aHat that are involved in SSB*)
	participating= KeySelect[ahatAssoc, SubsetQ[stabilityFields, #]&];
	
	If[Length@ ahatAssoc =!= Length@ nonParticipating+ Length@ participating,
		Message[InverseaHatMatrix::overlap];
		Abort[];
	];
	
	(*Construct ahat matrix for the SSB (in the collective B^\alpha space)*)
	participating= Table[
			If[f1 =!= f2, 1/2, 1]* participating[Sort@ {f1, f2}]
		, {f1, stabilityFields}, {f2, stabilityFields}];
	
	{participating, nonParticipating}
];


(* ::Subsubsection::Closed:: *)
(*Collected fields*)


(* ::Text:: *)
(*Produce reducible multiplets of ghosts, vectors (light/massive), and scalars in the broken phase set by the current SSB patterns *)


CollectiveFields[gbInfo_]:= Module[{B, c, cbar, chi, f, inds, phi, u, ubar, V, vev},
	(*Unbroken/gauge ghosts*)
	c= a|-> Evaluate[Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			$GhostFields[f][[1]][inds]
		, {f, LightVectorLabels[]}] ];
	(*Unbroken/gauge anti-ghosts*)
	cbar= a|-> Evaluate[Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			$GhostFields[f][[2]][inds]
		, {f, LightVectorLabels[]}] ];
	(*Unbroken vectors*)
	B= {mu, a}|-> Evaluate[Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			f[mu, inds]
		, {f, LightVectorLabels[]}] ];
	
	(*Massive ghosts*)
	u= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			If[$FieldAssociation[f, SelfConjugate], 
				{$GhostFields[f][[1]][inds]}, {$GhostFields[f][[1]][inds], Bar@ $GhostFields[f][[1]][inds]}
			]
		, {f, HeavyVectorLabels[]}] ];
	(*Massive anti-ghosts*)
	ubar= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			If[$FieldAssociation[f, SelfConjugate], 
				{$GhostFields[f][[2]][inds]}, {$GhostFields[f][[2]][inds], Bar@ $GhostFields[f][[2]][inds]}
			]
		, {f, HeavyVectorLabels[]}] ];
	(*Massive vectors*)
	V= {mu, a}|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			If[$FieldAssociation[f, SelfConjugate], 
				{f[mu, inds]}, {f[mu, inds], Bar@ f[mu, inds]}
			]
		, {f, HeavyVectorLabels[]}] ];
	
	(*Goldstone bosons*)
	chi= a|-> Evaluate[Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			chi= gbInfo[f];
			chi= If[Head@ chi === Conj, Bar@ chi[[1]][inds], chi[inds]];
			If[$FieldAssociation[f, SelfConjugate], 
				{chi}, {chi, Bar@ chi}
			]
		, {f, HeavyVectorLabels[]}] ];
	(*Scalars w/o VEV (phi, not phi')*)
	phi= a|-> Evaluate[Expand[#- (#/. _Field-> 0)]&/@ DecomposeFields@ Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			If[$FieldAssociation[f, SelfConjugate], 
				{f@ inds}, {f@ inds, Bar@ f@ inds}
			]
		, {f, SymmetryBreakingScalars[]}] ];
	(*VEV*)
	vev= a|-> Evaluate[DecomposeFields@ Flatten@ Table[
			inds= Sequence@@ ConstantArray[a, Length@ $FieldAssociation[f, Indices]];
			If[$FieldAssociation[f, SelfConjugate], 
				{f@ inds}, {f@ inds, Bar@ f@ inds}
			]
			, {f, SymmetryBreakingScalars[]}]/. _Field-> 0
		];

	<|
		Vectors-> B, 
		Ghost-> c,
		GhostBar-> cbar,
		GoldstoneBosons-> chi,
		MassiveGhost-> u,
		MassiveGhostBar-> ubar,
		MassiveVectors-> V,
		Scalars-> phi,
		VEV-> vev
	|>
];
