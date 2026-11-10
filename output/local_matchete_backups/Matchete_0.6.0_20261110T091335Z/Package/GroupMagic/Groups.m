(* ::Package:: *)

Package["Matchete`"] 


(* ::Title:: *)
(*GroupMagic`Groups`*)


(* ::Subtitle:: *)
(*Contains the necessary information for constructing groups  *)


(* ::Text:: *)
(*Largely based on [Cahn]: "Semi-Simple Lie Algebras and Their Representations" by Robert N. Cahn (1984).*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Export*)


PackageExport["CartanMatrix"]
PackageExport["AlgebraRank"]
PackageExport["GroupDimension"]


PackageExport["AdjointRepresentation"]
(*PackageExport["SmallestRepresentation"]*)
PackageExport["FundamentalRepresentation"]
PackageExport["PositiveRoots"]


PackageExport["RepresentationWeights"]
PackageExport["RepresentationDimension"]
PackageExport["RepresentationsUpToDimension"]
PackageExport["FSIndicator"]
PackageExport["DualRepresentation"]
PackageExport["RepresentationProduct"]
PackageExport["DynkinIndex"]
PackageExport["Casimir2"]


PackageExport["Levels"]
PackageExport["SymmetricIndices"]
PackageExport["AntisymmetricIndices"]
PackageExport["Orthonormalize"]
PackageExport["Normalization"]


PackageExport["CRep"]


PackageExport["Alg"]
PackageExport["SU"]
PackageExport["SO"]
PackageExport["Sp"]


PackageExport["Generators"]
PackageExport["InvariantTensors"]
PackageExport["StructureConstants"]


(* ::Subsubsection::Closed:: *)
(*Scope*)


PackageScope["TestAlg"]
PackageScope["RepresentationCheck"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Export*)


CartanMatrix::usage   = "CartanMatrix[alg] returns the Cartan matrix of a Lie algebra.";
AlgebraRank::usage    = "AlgebraRank[alg] returns the rank of a Lie algebra.";
GroupDimension::usage = "GroupDimension[alg] returns the number of generators.";
MatCommutator::usage  = "MatCommutator[A,B] = A.B-B.A";


AdjointRepresentation::usage     = "AdjointRepresentation[alg] provides the Dynkin coefficients of the adjoint representation of a Lie algebra."; 
FundamentalRepresentation::usage = "FundamentalRepresentation[alg] provides the Dynkin coefficients of the fundamental representation of a Lie algebra."; 
PositiveRoots::usage             = "PositiveRoots[alg] provides the weights of the positive roots of the Lie algebra."; 


RepresentationWeights::usage        = "RepresentationWeights[alg, DynkinCoef] returns the weights of a representation with multiplicity.";
RepresentationDimension::usage      = "RepresentationDimension[alg, DynkinCoef] returns the dimension of a representation.";
RepresentationsUpToDimension::usage = "RepresentationsUpToDimension[alg, n] returns all representations with dimensions \[LessEqual] n, as a list with {dim, Dynkin coefficients}.";
FSIndicator::usage                  = "FSIndicator[alg, DynkinCoef] returns the Frobenius-Schur indicator of the representation -- +1 for real, 0 for complex, and -1 for pseudo-real representations." ;
DualRepresentation::usage           = "DualRepresentation[alg, rep] returns the Dynkin coefficients of the dual representation to rep (the same if rep is (pseudo)real)." 
RepresentationProduct::usage        = "RepresentationProducut[alg, listOfReps] returns an association of the representations in the prouct of the representations in the input and their multiplicities.";
DynkinIndex::usage                  = "DynkinIndex[alg,rep] returns the Dynkin index of the representation.";
Casimir2::usage                     = "Casimir2[alg,rep] returns the quadratic Casimir of the representation.";


Levels::usage               = "Levels -> True/False is an option to indicate whether RepresentationWeights should be sorted by their level.";
SymmetricIndices::usage     = "SymmetricIndices[{set of numbers representing indices}] restricts to CGs with those indices completely symmetrized." ;
AntisymmetricIndices::usage = "AntisymmetricIndices[{set of numbers representing indices}] restricts to CGs with those indices completely antisymmetrized."  ;
Orthonormalize::usage       = "Orthonormalize -> True/False determines whether the to orthogonalize a set of CGs.";
Normalization::usage        = "Normalization -> num, normalizes the CGs to have norm num.";


CRep::usage = "CRep[DynkinCoef] denotes the conjugate of the representation with the given heighest weight." ;


Alg::usage = "Alg[class, index] is used to identify the various compact simple Lie algebras, e.g., Alg[\"A\", 2].";
SO::usage  = "SO[n] returns the Lie algebra associated with the ordinary Lie group";
Sp::usage  = "Sp[n] returns the Lie algebra associated with the ordinary Lie group";
SU::usage  = "SU[n] returns the Lie algebra associated with the ordinary Lie group";


Generators::usage         = "Generators[alg, rep] returns Hermitian and normalised generators of the representation in a basis common to all representations of the algebra.";
InvariantTensors::usage   = "InvariantTensors[alg, listOfReps] returns a basis of CG coefficients with indices contracted to representations in the input (order preserved). Possible options are SymmetricIndices/AntisymmetricIndices-> {set of indices by number} and Normalization -> num." ;
StructureConstants::usage = "StructureConstants[alg] returns the structure constants of the algebra in the basis common to the generators of all representations: [T^a, T^b]= i f^{abc} T^c.";


SmallestRepresentation::usage = "SmallestRepresentation[alg] provides the Dynkin coefficients of the representation with the smallest dimension.";
CSMatrices::usage             = "CSmatrices[alg, DynkinCoef] returns a set of matrices satisfying the Chevalley-Serre relations.";


(* ::Chapter:: *)
(*Private:*)


(* ::Subsection:: *)
(*Group/Algebra information*)


(* ::Subsubsection::Closed:: *)
(*Algebras *)


(* ::Text:: *)
(*Algebras of the simple Lie groups. The following algebras are isomorphic:*)
(*	A1 = C1 ~ B1*)
(*	D2 ~ A1 x A1 (not simple)*)
(*	B2 ~ C2*)
(*	D3 ~ A3*)


General::unkwnGroup= "`1`[`2`] is not a known simple group.";
General::notsimpleSO4= "SO(4) is not a simple group: its algebra is isomorphic to \!\(\*SubscriptBox[\(A\), \(1\)]\)\[Times]\!\(\*SubscriptBox[\(A\), \(1\)]\). Most of our methods fail for non-simple groups and so it is not supported. You can consider SU(2\!\(\*SuperscriptBox[\()\), \(2\)]\) instead.";  


SU@ n_Integer/; n>= 2:= Alg["A", n -1];
SU@ x_:= (Message[General::unkwnGroup, "SU", x]; Abort[]);


Sp@ n_Integer/; (EvenQ@ n && n>= 2):= 
	Switch[n
	,2,
		Alg["A", 1] (*SU(2) = Sp(2)*)
	,_, 
		Alg["C", n/2]
	]; 
Sp@ x_:= (Message[General::unkwnGroup, "Sp", x]; Abort[]);


SO@ n_Integer/; (n>= 3):= 
	Switch[n
	,4,
		Message[General::notsimpleSO4];
		Abort[];
	,_? OddQ,
		Alg["B", (n- 1)/ 2]
	,_, 
		Alg["D", n/2]
	]; 
SO@ x_:= (Message[General::unkwnGroup, "SO", x]; Abort[]);


(* ::Text:: *)
(*Algebra test *)


Alg::unkwn= "`1` is not a valid simple Lie algebra.";


LieAlgQ@ Alg[name_String, n_Integer]:= 
	Switch[name, 
		"A", n> 0,
		"B"|"C", n> 0,
		"D", n> 1, 
		"E", MatchQ[n, 6|7|8],
		"F", n === 4, 
		"G", n === 2,
		_, False];
LieAlgQ@ _= False;  


TestAlg@ alg_:= If[!LieAlgQ@ alg, Message[Alg::unkwn, alg]; Abort[];];   


(* ::Subsubsection::Closed:: *)
(*Determining Cartan matrices *)


(* ::Text:: *)
(*Returns the Cartan Matrix associated with an algebra*)


CartanMatrix::unkwn= "`1` is not a known simple Lie algebra.";


CartanMatrix[alg: Alg[Name_String, n_Integer]]:= CartanMatrix[alg]= Block[{temp},
	TestAlg@ alg;
	Switch[alg
	,Alg["A", 1], 
		temp= SparseArray[{1, 1} -> 2, {1, 1}]
	,Alg["A", m_]/; m> 1,  
		temp= SparseMatrixBand[n]
	,Alg["B", 1], 
		temp= CartanMatrix@ Alg["A", 1]
	,Alg["B", m_]/; m> 1, (*Note that this convention is different from the isomorphic C_2*)
		temp= SparseMatrixBand[n];
		temp[[-2, -1]]= -2;
	,Alg["C", 1],
		temp= CartanMatrix@ Alg["A", 1]
	,Alg["C", m_]/; m> 1,   
		temp= SparseMatrixBand[n];
		temp[[-1, -2]]= -2;
	,Alg["D", 2],
		temp= DiagonalMatrix@ {2, 2};
	,Alg["D", 3],
		temp= CartanMatrix@ Alg["A", 3]
	,Alg["D", m_]/; m> 3, 
		temp= SparseMatrixBand[n];
		temp[[-1, -2]]= temp[[-2, -1]]= 0;
		temp[[-1, -3]]= temp[[-3, -1]]= -1;
	,Alg["E", 6| 7| 8],
		temp= SparseMatrixBand[n];
		temp[[-1, -2]]= temp[[-2, -1]]= 0;
		temp[[-4, -1]]= temp[[-1, -4]]= -1;
	,Alg["F", 4],
		temp= SparseMatrixBand[4];
		temp[[2, 3]]= -2;
	,Alg["G", 2],
		temp= SparseMatrixBand[2];
		temp[[1, 2]]= -3;
	];
	Normal@ temp
];
CartanMatrix[x_]:= (Message[CartanMatrix::unkwn, x]; Abort[];)


SparseMatrixBand@ n_:= SparseArray[{Band@ {1, 1}-> 2, Band@ {2, 1}-> -1, Band@ {1, 2}-> -1}, {n, n}];


(* ::Subsubsection::Closed:: *)
(*Rank*)


(* ::Text:: *)
(*The rank of an algebra*)


AlgebraRank@ alg_:= Length@ CartanMatrix@ alg;


(* ::Subsubsection::Closed:: *)
(*Algebra dimension *)


GroupDimension[Matchete`U1]:= 1;
GroupDimension[alg_]:= RepresentationDimension[alg, AdjointRepresentation[alg]];


(* ::Section:: *)
(*Roots*)


(* ::Subsection:: *)
(*Weights of the roots *)


(* ::Subsubsection::Closed:: *)
(*Positive roots*)


(* ::Text:: *)
(*[Cahn, VII]: The positive roots are determined  starting from the zero weight, adding simple roots step by steps. The weights of negative roots are minus that of the positive (if \[Alpha] is a root then so is -\[Alpha]).  *)


PositiveRoots[alg_Alg]:= PositiveRoots[alg]= Block[{rank, weights, prevWeights, i, m, p, pos,
	cartanMat= CartanMatrix@ alg},
	rank= Length@ cartanMat;
	
	weights= Last@ Reap[
		prevWeights= Association@@ Table[cartanMat[[i]]-> -UnitVector[rank, i] (*m_i vector*)
			,{i, rank}]// Sow;
		While[Length@ prevWeights> 0,
			prevWeights= Sow@ Merge[Table[
				m= prevWeights@ w;
				p= m- w;
				pos = Position[p, x_/; x> 0,{1}][[;;, 1]];
				Association@@ Table[w+ cartanMat[[i]]-> (m[[i]]+ 1) UnitVector[rank, i]
					,{i, pos}]
			,{w, Keys@ prevWeights}], Total];
		];
	];
	Join@@ weights[[1]]// Keys 
];


PositiveRoots::arg= "PositiveRoots[alg] takes a Lie algebra as argument."
PositiveRoots@___:= (Message[PositiveRoots::arg]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Adjoint representation*)


(* ::Text:: *)
(*Returns the Dynkin coefficient of the adjoint representation of a given Lie Algebra*)


AdjointRepresentation@ alg_Alg:= AdjointRepresentation@ alg= Block[{invCM, posRoots},
	TestAlg@ alg;
	posRoots= PositiveRoots@ alg;
	invCM= Inverse@ CartanMatrix@ alg;
	First@ MaximalBy[posRoots, (Total[# . invCM]&)]
];


AdjointRepresentation::arg= "AdjointRepresentation[alg] takes a Lie algebra as argument."
AdjointRepresentation@ ___:= (Message[AdjointRepresentation::arg]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Fundamental representation*)


(* ::Text:: *)
(*Returns the Dynkin coefficient of the fundamental representation of a given Lie Algebra*)


FundamentalRepresentation@ alg: Alg[_, n_Integer]:= Block[{},
	TestAlg@ alg;
	Switch[alg
	,Alg["B", 1],
		{2}
	,Alg["B", 2],
		{1, 0}
	,Alg["B", 2],
		{1, 1}
	,Alg["D", 3],
		{0, 1, 0}
	,Alg["A"|"B"|"C"|"D", _],
		Join[{1}, ConstantArray[0, AlgebraRank@ alg -1]]
	,Alg["E"|"F"|"G", _],
		SmallestRepresentation@ alg
	]
];


FundamentalRepresentation::arg= "FundamentalRepresentation[alg] takes a Lie algebra as argument."
FundamentalRepresentation@ ___:= (Message[FundamentalRepresentation::arg]; Abort[];) 


SmallestRepresentation@ alg_Alg:= 
	RepresentationsUpToDimension[alg, RepresentationDimension[alg, AdjointRepresentation@ alg]][[2, 2]];


(* ::Subsection:: *)
(*Lengths of the simple roots and metric  *)


(* ::Subsubsection::Closed:: *)
(*Simple root lengths *)


(* ::Text:: *)
(*The Cartan matrix is given by A_ij = 2 <\[Alpha]_ i, \[Alpha]_j>/ <\[Alpha]_ j, \[Alpha]_ j> [Cahn, VII], so when <\[Alpha]_i, \[Alpha]_ j> != 0,  A_ ij/ A_ji = <\[Alpha]_ i, \[Alpha]_i>/ <\[Alpha]_j, \[Alpha], j>. *)


(* ::Text:: *)
(*Gives the length of the simple roots, assuming that the longest is length Sqrt[2] in the default case.*)
(*For trace normalization purpose, one can choose the convention where the shortest is 1 by using the option RootSqMin->1.*)


Options@ RootLengthSquared= {RootSqMax -> 2, RootSqMin -> 0}
Options@ CartanKillingMetric= {RootSqMax -> 2, RootSqMin -> 0}


(*Special caes for the non-simple algebra*)
RootLengthSquared[Alg["D", 2], opt:OptionsPattern[]]:= RootLengthSquared[Alg["D", 2], opt]= Block[{},
	If[OptionValue[RootSqMin]=== 0,
		OptionValue[RootSqMax],
		OptionValue[RootSqMin]
	]* IdentityMatrix[2]
];


(* ::Text:: *)
(*This algorithm only works for simple algebras*)


RootLengthSquared[alg_,opt:OptionsPattern[]]:= RootLengthSquared[alg,opt]= Block[{cm, lengthUnkown, lengthsSq, rank, currentRoot, 
	pos, i, j, next, n},
	cm= CartanMatrix@ alg;
	rank= Length@ cm; 
	lengthsSq= ConstantArray[0, rank];
	lengthUnkown= ConstantArray[True, rank];
	currentRoot= ConstantArray[0, rank];
	lengthsSq[[1]]= 2; lengthUnkown[[1]]= False; 
	currentRoot[[1]]= 1; next= 2;
	Do[
		i= currentRoot[[n]];
		pos= Flatten@ Position[cm[[i]], _?(# < 0&), {1}];
		Do[
			If[lengthUnkown[[j]],
				lengthsSq[[j]]= lengthsSq[[i]] cm[[j, i]]/ cm[[i, j]];
				lengthUnkown[[j]]= False;
				currentRoot[[next++]]= j;
			];
		,{j, pos}];		
	 ,{n, rank}]; 
	
	If[OptionValue[RootSqMin]===0,
		lengthsSq *OptionValue[RootSqMax]/ Max@ lengthsSq,
		lengthsSq *OptionValue[RootSqMin]/ Min@ lengthsSq
	]
];


(* ::Subsubsection::Closed:: *)
(*Cartan-Killing metric *)


(* ::Text:: *)
(*The Cartan-Killing metric, G_ij, on root space in the basis of Dynkin coefficients, so <M,N> = (M)_i G_ij (N)_j*)


CartanKillingMetric[alg_, opt:OptionsPattern[]]:= CartanKillingMetric[alg,opt]= 
	1/2 Inverse@ CartanMatrix@ alg . DiagonalMatrix@ RootLengthSquared[alg,opt];


(* ::Section:: *)
(*Representations*)


(* ::Subsection:: *)
(*Properties*)


(* ::Text:: *)
(*CRep is used to denote representations whose vectors are obtained by complex conjugating vectors in a representation. *)


CRep@ CRep@rep_:= rep;


(* ::Subsubsection::Closed:: *)
(*Representation check *)


(* ::Text:: *)
(*Conjugate representations*)


General::realRep= "`1` is a real representation of `2`. CRep cannot be used here.";


RepresentationCheck[alg_, CRep@ \[CapitalLambda]_]:= Block[{},
	RepresentationCheck[alg, \[CapitalLambda]]; 
	If[FSIndicator[alg, \[CapitalLambda]] === 1,
		Message[General::realRep, \[CapitalLambda], alg];
		Abort[];
	];
];


(* ::Text:: *)
(*Other representations*)


General::dynkCoef= "`1` is not a valid Dynkin Coefficient.";
General::rankLab= "The length of `1` does not match the rank of `2`.";


RepresentationCheck[alg_, \[CapitalLambda]_]:= Block[{},
	TestAlg@ alg;
	If[!MatchQ[\[CapitalLambda], {(_Integer? NonNegative)..}],
		Message[General::dynkCoef, \[CapitalLambda]];
		Abort[];		
	]; 
	If[Length@ \[CapitalLambda] =!= AlgebraRank@ alg,
		Message[General::rankLab, \[CapitalLambda], alg];
		Abort[];
	];
	DoubleCoverRepCheck[alg, \[CapitalLambda]];
];


(* ::Text:: *)
(*Representations for the double-covers*)


General::doublecover= "The representation `1` exists only in the double cover (`2`) of `3`.";


DoubleCoverRepCheck[alg:Alg["B", 1], rep:{_Integer? OddQ}]:= Block[{},
	Message[General::doublecover, rep, "Alg[A, 1]", alg];
	Abort[];
];
DoubleCoverRepCheck[alg:Alg["B", 2], rep:{_, _Integer? OddQ}]:= Block[{},
	Message[General::doublecover, rep, "Alg[C, 2]", alg];
	Abort[];
];
DoubleCoverRepCheck[alg:Alg["D", 2], rep:{a_, b_}]/; OddQ[a+ b]:= Block[{},
	Message[General::doublecover, rep, "Alg[A, 1] \[Times] Alg[A, 1]", alg];
	Abort[];
];
DoubleCoverRepCheck[alg:Alg["D", 3], rep:{a_, _, b_}]/; OddQ[a+ b]:= Block[{},
	Message[General::doublecover, rep, "Alg[A, 3]", alg];
	Abort[];
];


(* ::Subsubsection::Closed:: *)
(*Dynkin index and quadratic Casimir *)


(* ::Text:: *)
(*Casimir2 gives the quadratic Casimir of the representation.*)
(*DynkinIndex gives the trace normalization of the representation.*)
(*Both use the normalization condition where the trace over fundamental generators is 1/2, which is equivalent to requiring the length squared of the shortest root to be 1.*)


Casimir2[alg_,rep_]:= Block[{},
	RepresentationCheck[alg, rep];
	rep . CartanKillingMetric[alg, RootSqMin->1] . (rep+ ConstantArray[2, Length@ CartanMatrix[alg]]) 
];
Casimir2[alg_, CRep@ rep_]:= Casimir2[alg, rep]


Casimir2::arg= "Casimir2[alg, rep] takes a Lie algebra and a representation as arguments."
Casimir2@ ___:= (Message[Casimir2::arg]; Abort[];) 


DynkinIndex[alg_,rep_]:= Block[{},
	RepresentationCheck[alg, rep];
	RepresentationDimension[alg, rep]/ GroupDimension[alg] Casimir2[alg,rep]
];
DynkinIndex[alg_, CRep@ rep_]:= DynkinIndex[alg, rep]


DynkinIndex::arg= "DynkinIndex[alg, rep] takes a Lie algebra and a representation as arguments."
DynkinIndex@ ___:= (Message[DynkinIndex::arg]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Weights and dimensions of the corresponding weight spaces  *)


(* ::Text:: *)
(*Returns an association with all the weight and dimensionality of the associated weight spaces of representation with highest weight (\[CapitalLambda])_i.*)


Options@ RepresentationWeights= {Levels-> False};
Options@ RepresentationWeightsInternal= {Levels-> False}


(* ::Text:: *)
(*The weights of the conjugate representation (produced by taking \[Phi]^* of weight vector \[Phi]) are minus the weight of the original vector.*)


RepresentationWeights[alg_, \[CapitalLambda]_, opt:OptionsPattern[] ]? OptionsCheck:= Block[{},
	RepresentationCheck[alg, \[CapitalLambda]];
	RepresentationWeightsInternal[alg, \[CapitalLambda], opt]
];


RepresentationWeightsInternal[alg_, CRep@ \[CapitalLambda]_List, opt:OptionsPattern[]]:= 
	KeyMap[-#&, RepresentationWeightsInternal[alg, \[CapitalLambda], opt]];


RepresentationWeightsInternal[alg_, \[CapitalLambda]_, opt:OptionsPattern[] ]:= RepresentationWeights[alg, \[CapitalLambda], opt]=
Block[{rank, weights, prevWeights, m, p, w, pos, i, M,
	weightDimensions, posRoots, dMetric, level, root, dim, num,
	cartanMat= CartanMatrix@ alg},	
	rank= Length@ cartanMat; 
	
	(*Determine the weights of the representation [Cahn, X]*)
	weights= Reap[
		prevWeights= Sow@ <|\[CapitalLambda]-> ConstantArray[0, rank] (*p_i vector*)|>;
		While[Length@ prevWeights> 0,
			prevWeights= Sow@ Merge[Table[
				p= prevWeights@ w;
				m= w+ p;
				pos = Position[m, x_/; x> 0,{1}][[;;, 1]];
				Association@@ Table[w- cartanMat[[i]]-> (p[[i]]+ 1) UnitVector[rank, i]
					,{i, pos}]
			,{w, Keys@ prevWeights}], Total];
		];
	][[-1, 1, ;;-2]]; (*Last element is an empty association*)
	weights= Keys/@ weights;
	
	(*Determine dimensionality of the weights with Freudenthal's formula [Cahn, XI]*)	
	posRoots= PositiveRoots@ alg;
	dMetric= CartanKillingMetric@ alg;	
	
	(*Setup association with dimensions of the weight spaces*)
	weightDimensions= Join@@ weights; 
	weightDimensions= Association@@ (Rule[#, 0]&)/@ weightDimensions;
	weightDimensions@ \[CapitalLambda]= 1; 
	
	Do[(*Determine dimension of weight space layer by layer*)
		num= 0;
		Do[(*Sum of all positive roots in Freudenthtal's formula *)
			w= M+ root;
			dim= Lookup[weightDimensions, Key@ w, 0];
			While[dim> 0,
				num+= dim (w . dMetric . root);
				w+= root;
				dim= Lookup[weightDimensions, Key@ w, 0];
			];
		,{root, posRoots}];
		(*Multiply by demoniator*)
		weightDimensions@ M= 2 num/ ((\[CapitalLambda]+ M+ ConstantArray[2, rank]) . dMetric . (\[CapitalLambda]- M));
	,{level, weights[[2;;]]}, {M, level}];
	
	If[OptionValue@ Levels,
		Association@@@ Map[(#-> weightDimensions@ # &), weights, {2}]
	,
		weightDimensions
	]
]; 


(* ::Subsubsection::Closed:: *)
(*Dimension of the representation*)


(* ::Text:: *)
(*The dimension of a representation can be determined efficiently with Weyl's dimension formula [Cahn, XIII]*)


RepresentationDimension[alg_, \[CapitalLambda]_]:= Block[{},
	RepresentationCheck[alg, \[CapitalLambda]];
	RepresentationDimensionInternal[alg, \[CapitalLambda]]
] 


RepresentationDimensionInternal[alg_, \[CapitalLambda]_]:= Block[{roots, dMetric, r, \[Delta], temp},
	roots= PositiveRoots@ alg;
	dMetric= CartanKillingMetric@ alg;
	\[Delta]= ConstantArray[1, Length@ \[CapitalLambda]]; 
	Product[
		temp= r . dMetric;
		(temp . (\[CapitalLambda]+ \[Delta]))/ (temp . \[Delta])
	, {r, roots}] 
];


RepresentationDimensionInternal[alg_, CRep@ \[CapitalLambda]_List]:= RepresentationDimensionInternal[alg, \[CapitalLambda]];


(* ::Subsubsection::Closed:: *)
(*Representations up to given dimension*)


(* ::Text:: *)
(*Function for returning all representations up to a given dimension*)


RepresentationsUpToDimension[alg_, dimLimit_Integer]:= Block[{carry, dim, i, rank, rep, representations},
	TestAlg@ alg;
	If[dimLimit< 1, Return@{}; ];
	rank= Length@ CartanMatrix@ alg;
	representations= Reap[
		rep= ConstantArray[0, rank];
		Sow@ {1, rep};
		i=1;
		While[i<= rank,
			rep[[i]]++;
			dim= RepresentationDimensionInternal[alg, rep];
			If[dim> dimLimit,
				rep[[i++]]= 0;
				Continue[];
			];
			Sow@ {dim, rep};
			i= 1;
		];
	][[2, 1]];
	
	representations= Select[representations, !InvalidRep[alg, #[[2]]]&];
	
	ReverseSortBy[representations, (-First@# &)]
];


InvalidRep[Alg["B", 1|2], label_]:= OddQ[label[[-1]]];
InvalidRep[Alg["D", 2|3], label_]:= OddQ[label[[1]]+ label[[-1]]];
InvalidRep[__]:= False;


(* ::Subsection:: *)
(*Conjugate (dual) representations*)


(* ::Subsubsection::Closed:: *)
(*Frobenius-Schur indicator*)


(* ::Text:: *)
(*The Frobenius-Schur indicator indicates whether a given representation is +1 for real, 0 for complex, and -1 for quateronic (pseudo-real) representations.  *)
(*https://mathoverflow.net/questions/227869/is-there-a-formula-for-the-frobenius-schur-indicator-of-a-rep-of-a-lie-group*)


FSIndicator[alg_, CRep@ \[CapitalLambda]_]:= FSIndicator[alg, \[CapitalLambda]];


FSIndicator[alg_, \[CapitalLambda]_]:= FSIndicator[alg, \[CapitalLambda]]= Block[{weights, temp},
	RepresentationCheck[alg, \[CapitalLambda]];
	
	weights= RepresentationWeights[alg, \[CapitalLambda]];
	(*If not selfdual, the representation is complex*)
	temp= KeyMap[Minus, weights];
	If[! And@@ Merge[{weights, temp}, (Length@ # === 2 && #[[1]] === #[[2]] &)],
		Return@ 0;
	];
	
	temp= Total@ PositiveCoroots@ alg . CartanKillingMetric@ alg . \[CapitalLambda]; 
	Power[-1, temp]
];


PositiveCoroots@ alg_:= Block[{roots, dMetric= CartanKillingMetric@ alg},
	roots= PositiveRoots@ alg;
	2 #/ (# . dMetric . #)&/@ roots
]  


(* ::Subsubsection::Closed:: *)
(*Determine conjugate representation*)


(* ::Text:: *)
(*Function for finding the Dynkin coefficients of the dual representation. *)
(*https://en.wikipedia.org/wiki/Dual_representation#General_semisimple_Lie_algebras*)


DualRepresentation[alg_, rep_]:= Block[{weights, invCM},
	RepresentationCheck[alg, rep];
	
	(*(Pseudo)real representations are self dual*)
	If[FSIndicator[alg, rep] =!= 0, Return@ rep;];
	
	(*Find the largest weight in the dual representaiton*)
	invCM= Inverse@ CartanMatrix@ alg;
	weights= Keys@ KeyMap[Minus, RepresentationWeights[alg, rep]];
	First@ MaximalBy[weights, (# . invCM &)]
];


(* ::Subsection:: *)
(*Product representations *)


(* ::Subsubsection::Closed:: *)
(*RepresentationProduct*)


(* ::Text:: *)
(*The irreducible representations in a product representation can be determined by successively removing the irreducible representations of the largest remaining weight from the weight space of the product. *)
(*Alternative: Konstant's formula?  *)


RepresentationProduct::invldarg= "`1` is not a valid argument; a list of representations is expected."


(* ::Text:: *)
(*Returns an association of representations in the product and their multiplicities. *)


RepresentationProduct[alg_, reps_]:= Block[{productWeights, invCM, rep, dim},	
	If[!MatchQ[reps, {(_List|CRep@ _List)..}], 
		Message[RepresentationProduct::invldarg, reps];
		Abort[];
	];
	RepresentationCheck[alg, #]&/@ reps; 
	
	(*Determine weights and their multiplicity in the product*)
	productWeights= RepresentationWeights[alg, #]&/@ reps; 
	productWeights= Rule[Total@ #[[;;, 1]], Times@@ #[[;;, 2]]]&/@ Tuples[Normal/@ productWeights];
	productWeights= Normal@ Merge[productWeights, Total];
	
	(*Order weights*)
	invCM= Inverse@ CartanMatrix@ alg;
	productWeights= Association@ Sort[productWeights, (FirstCase[(#1[[1]]- #2[[1]]) . invCM, Except@ 0, 0] >=0 &)];
	
	(*Successively removes hegihest weight representations from the weight space to determine irreducible factors*)
	Association@Reap[
		Do[
			dim= productWeights@ rep;
			If[dim<= 0, Continue[];];
			Sow[rep-> dim];
			productWeights= Merge[{productWeights, -dim RepresentationWeights[alg, rep]}, Total]; 
		,{rep, Keys@ productWeights}];
	][[2, 1]]
];


(* ::Section:: *)
(*Invariants*)


(* ::Text:: *)
(*The various invariant tensors and generators of representations are all dependent on the particular choice of basis states.*)


(* ::Subsection:: *)
(*Representation matrices/generators (Chevalley-Serre)*)


(* ::Subsubsection::Closed:: *)
(*Adjoint representation*)


(* ::Text:: *)
(*Returns the sets of {H_i, E_i, F_i} for the simple  roots, obeying the Chevalley-Serre relations of the adjoint representation. *)
(*The method uses the CS matrices for the fundamental representation to derive the CS matrices for the adjoint representation. *)


CSMatrices[alg_, \[CapitalLambda]_List]/; \[CapitalLambda] === AdjointRepresentation@ alg:= CSMatrices[alg, \[CapitalLambda]]= Block[
	{Hi, Ei, Fi, invMetric, repMatrices, invRot, minRep,
		cartanMat= CartanMatrix@ alg,
		rank= AlgebraRank@ alg
		},
	
	(*Relies on some other representation than the adjoint*)
	minRep= SmallNonAdjRepresentation@ alg;
	
	{Hi, Ei, Fi}= Transpose@ CSMatrices[alg, minRep];
	Ei= Association@@ Thread[Rule[cartanMat, Ei]];
	Fi= Association@@ Thread[Rule[-cartanMat, Fi]]; 
	
	repMatrices = GeneratorBasis[alg, minRep];
	
	(*Some performance might be gained here:*)
	invMetric= SparseArray@ Inverse@ Outer[MatInnerProd, repMatrices, repMatrices, 1];
	
	(*Determine the action of {Hi, Ei, Fi} on the basis of the generators, which are weight vectors*) 
	Table[Transpose/@{
		 SparseArray@ Outer[RepAction@ Hi[[i]], repMatrices, repMatrices, 1] . invMetric,
		 SparseArray@ Outer[RepAction@ Ei@ cartanMat[[i]], repMatrices, repMatrices, 1] . invMetric ,
		 SparseArray@ Outer[RepAction@ Fi[- cartanMat[[i]]], repMatrices, repMatrices, 1] . invMetric
		}
	,{i, rank}] //Simplify 
]


RepAction[L_][A_, B_]:= MatInnerProd[B, MatCommutator[L, A]];


SmallNonAdjRepresentation@ alg_:= Switch[alg
	,Alg["E", 8],
		{0, 0, 0, 0, 0, 0, 1, 0}(*Suitable representation to compute CS matrices*)
	,Alg["B", 1],
		{4}
	,_,
		SmallestRepresentation@ alg
	];


(* ::Subsubsection::Closed:: *)
(*Generic representations*)


(* ::Text:: *)
(*Returns a the set of {H_i, E_i, F_i} for the simple  roots, obeying the Chevalley-Serre relations, in the given representation [Fonseca, B1]. The basis of the representation space is always chosen with weight vectors, i.e. the basis states are eigenvectors of the adjoint generators. *)


CSMatrices::decomp = "Inconsistent decomposition. CSMatrices failed."
CSMatrices::singlet = "Cannot be applied on singlet representations."


CSMatrices[alg_, \[CapitalLambda]_List]:= CSMatrices[alg, \[CapitalLambda]]= Block[{weights, weightDims, level, M, rank, emptyMatrices,
	raisingMatrices, X, i, j, targetWeights, E1, E2, targetDims, \[CapitalOmega], offSets, repDimension, temp,
	cartanMat= CartanMatrix@ alg},
	If[MatchQ[\[CapitalLambda], {0..}],
		(*In this case the representation is a singlet *)
		Message[CSMatrices::singlet];
		Abort[];
	];
	
	rank= Length@ cartanMat;
	weights= RepresentationWeights[alg, \[CapitalLambda], Levels-> True];
	
	(*The goal here is to construct the raising matrices [\[Omega]]^i of [Fonseca, (B.29)], 
	such that E_i |\[Omega]> = ([\[Omega]]^i)^T |\[Omega]+ \[Alpha]_i> *)
	emptyMatrices= ConstantArray[{}, rank];
	weightDims= raisingMatrices= Join@@ weights;
	raisingMatrices[[;;]]= ConstantArray[emptyMatrices, {Length@ raisingMatrices}];
	Do[ (*Perhaps this loop should be implemented with sparse arrays*)
		targetWeights= Table[w + cartanMat[[i]], {i, rank}];
		targetDims= Lookup[weightDims, targetWeights, 0];
		(*Determine the symmetric X of [Fonseca, (B.29)]*)
		X= Table[
			If[targetDims[[i]] === 0, 
				{}
			,
				Table[
					If[targetDims[[j]] === 0, 
						{}
					,		
						E1= Lookup[raisingMatrices, Key@ targetWeights[[i]], emptyMatrices][[j]];
						E2= Lookup[raisingMatrices, Key@ targetWeights[[j]], emptyMatrices][[i]];
						If[E1 === {} || E2 === {},
							ConstantArray[0, targetDims[[{i, j}]]]
						,
							Transpose@ E1 . E2
						]+ 
						If[i === j,
							targetWeights[[i, i]] IdentityMatrix@ targetDims[[i]]
						,
							0
						]
					]
				, {j, rank}]
			]
		,{i, rank}];
		X= ArrayFlatten@ DeleteCases[X, {}, 2];
		
		(*Determine appropriate rasing matrix*)
		\[CapitalOmega]= CholeskyLikeDecomposition@ X;
		\[CapitalOmega]= Transpose@ DeleteCases[Transpose@ \[CapitalOmega], {0..}];
		\[CapitalOmega]= PadRight[\[CapitalOmega], {Length@ \[CapitalOmega], weightDims@ w}];
		If[Simplify@ X =!= Simplify[\[CapitalOmega] . Transpose@ \[CapitalOmega]], (*Test consistency of solution*)
			Message[CSMatrices::decomp];
			Abort[];
		]; 
		raisingMatrices@ w= TakeList[\[CapitalOmega], targetDims];
	,{layer, weights[[2;;]]}, {w, Keys@ layer}];
	
	(*Construct the Sparse Arrays*)
	weights= Keys@ weightDims;
	offSets= WeightPositions[alg, \[CapitalLambda]];
	repDimension= Total@ weightDims;
	Table[{ (*H_i*)
		SparseArray[ 
			Flatten@ Table[({#, #}-> w[[i]] &)/@ (offSets@ w + Range@ weightDims@ w)
				, {w, weights}]
		,{repDimension, repDimension}]
	, (*E_i, [Fonseca, (B.29)]*)
		E1= SparseArray[ Flatten@ Table[
			temp= raisingMatrices[w][[i]];
			If[temp=== {}, {},
				temp= ArrayRules[temp][[;; -2]];
				temp[[;;, 1]] += ConstantArray[{offSets[w+ cartanMat[[i]]], offSets@ w}, Length@ temp];
				temp
			]
		, {w, weights}]
		,{repDimension, repDimension}]
	, (*F_i*)
		Transpose@ E1
	},{i, rank}]
	
];


(* ::Text:: *)
(*The CS matrices of the conjugate representation (produced by taking \[Phi]^* of weight vector \[Phi]) are the negative transposed of the original ones. This parallels the generators -T^* acting on \[Phi]^* to the T of \[Phi] (using T\[ConjugateTranspose]=T).  *)


CSMatrices[alg_, CRep@ \[CapitalLambda]_List]:= -Map[Transpose, CSMatrices[alg, \[CapitalLambda]], {2}];


(* ::Subsection:: *)
(*Generators*)


(* ::Subsubsection::Closed:: *)
(*Generator basis (in the weight basis)*)


(* ::Text:: *)
(*Determines a set of generators for a representation corresponding to a the weight-vector basis in the Lie algebra, by extending the set of CS matrices for the representation. *)


GeneratorBasis[alg_, \[CapitalLambda]_]:=Block[{Hi, Ei, Fi, E\[Alpha], repMatrices, zeroWeight,
		rank= AlgebraRank@ alg,
		weightPos= WeightPositions[alg, AdjointRepresentation@alg],
		cartanMat= CartanMatrix@ alg
		},
	RepresentationCheck[alg, \[CapitalLambda]];
	
	zeroWeight= ConstantArray[0, rank];

	{Hi, Ei, Fi}= Transpose@ CSMatrices[alg, \[CapitalLambda]];
	
	(*Construct a basis of generators*)
	Ei= Association@@ Thread[Rule[cartanMat, Ei]];
	E\[Alpha]= BasisForRoots[Ei, Keys@ weightPos, ReverseCommutator->True];
	Fi= Association@@ Thread[Rule[-cartanMat, Fi]]; 
	E\[Alpha]= Join[E\[Alpha], BasisForRoots[Fi, Keys@ weightPos]];
	
	(*Matches the ladder and Cartan generators to the adjoint index in the weight basis. *)
	repMatrices= ConstantArray[SparseArray[{}, {1, 1} RepresentationDimensionInternal[alg, \[CapitalLambda]]], GroupDimension[alg]];
	Do[
		repMatrices[[weightPos@ root+ 1]]= E\[Alpha]@ root; 
	, {root, Keys@ E\[Alpha]}]; 
	repMatrices[[weightPos@ zeroWeight+ 1;; weightPos@ zeroWeight+ rank]]= Hi; 
	
	SparseArray@repMatrices
]


(* ::Text:: *)
(*Takes the representation matrices of +/- the simple roots (e.g. E_i or F_i from CS representations) in an association <|root-> mat|> and a list of all roots in the algebra. Returns representation matrices for all roots that can be obtained with commutators of the input *)


BasisForRoots[simpBasis_, roots_,OptionsPattern[{ReverseCommutator->False}]]:= Block[{rootCheck, rootBasis, lvl, \[Alpha]},
	rootCheck= Association@@ (Rule[#, False]&/@ roots);
	rootBasis= {};
	lvl= simpBasis;
	While[Length@ lvl> 0,
		AppendTo[rootBasis, lvl];
		lvl= Reap[Do[
			Do[
				\[Alpha]= root+ simpRoot;
				If[Lookup[rootCheck, Key@ \[Alpha], True], Continue[]; ];
				Sow[\[Alpha]-> If[OptionValue@ ReverseCommutator, -1, 1] MatCommutator[simpBasis@ simpRoot, lvl@ root]];
				rootCheck@ \[Alpha]= True; 
			, {simpRoot, Keys@ simpBasis}];
		, {root, Keys@ lvl}] ][[2]];
		If[Length@ lvl> 0, lvl= First@ lvl;];
		lvl= Association@ lvl;
	];
	Join@@ rootBasis 
];


(* ::Subsubsection::Closed:: *)
(*BasisRotation*)


(* ::Text:: *)
(*Given the weight basis u^a of the Lie algebra, the generators of a representation T^a  = \[Rho](u^a) are not Hermitian and orthonormal. A linear transformation L^a_b is found such that T'^a = L^a_b T^b satisfy the usual properties for the basis generator. The corresponding basis of the Lie algebra u'^a = L^a_b u^b does not consist of weight vectors. *)


BasisRotation[alg_]:= BasisRotation[alg]= Block[
	{rotation, generatorBasis, Ei, FEcorrespondance, Hi, metric, norm, minRep,
	nG=GroupDimension@alg, 
	rank = AlgebraRank@alg,
	roots=Keys@RepresentationWeights[alg,AdjointRepresentation@alg]},
	
	minRep= SmallNonAdjRepresentation@ alg;
	
	generatorBasis= GeneratorBasis[alg, minRep];
	
	(*Matches basis vectors in ladder pairs (step up and down).*)
	FEcorrespondance=(Flatten@{Position[roots, #], Position[roots,-#]+rank-1})&/@ (roots[[;;(nG-rank)/2]]);
	(*Ei=generatorBasis[[;;(nG-rank)/2]];
	FEcorrespondance=(Flatten@{Position[generatorBasis, #], Position[generatorBasis, Transpose@ #]})&/@ Ei;*)
	(*Position[generatorBasis, #] does not work anymore because generatorBasis is a full SparseArray*)
		
	norm= Sqrt[DynkinIndex[alg, minRep]]/ Sqrt[2 MatNorm[#]&/@ (generatorBasis[[;;(nG-rank)/2]]) (*Ei*) ];
	norm= Join[norm, ConstantArray[1,rank], Reverse@norm];
	
	Hi=generatorBasis[[ (nG-rank)/2 + 1 ;; (nG+rank)/2 ]] ;
	metric= SparseArray[Table[ MatInnerProd[i, j],{i,Hi},{j,Hi}]/DynkinIndex[alg,minRep]];
	
	rotation= norm SparseArray[Flatten@({ConstantArray[#[[1]],2] -> 1 , # -> 1 , Reverse@# -> -I ,ConstantArray[#[[2]],2] -> I }&/@FEcorrespondance), {nG,nG}];
	rotation[[ (nG-rank)/2 + 1 ;; (nG+rank)/2 , (nG-rank)/2 + 1 ;; (nG+rank)/2 ]]=
		SparseArray@ GramSchmidtMetric@ metric;
	
	rotation
]


(* ::Subsubsection::Closed:: *)
(*Generators*)


(* ::Text:: *)
(*Returns the generator matrices of a representation in a basis satisfying the usual properties of Hermitianity and orthonormality *)


Generators[alg_, rep_]:= Generators[alg, rep]= Block[{generators, L, 
	rot = BasisRotation[alg]},
	RepresentationCheck[alg, rep];
	generators= GeneratorBasis[alg, rep];
	(*Rotate generators of real representaiton according to the standard basis of the representaiton.*)
	If[FSIndicator[alg, rep] === 1,
		L= RealBasisRotation[alg, rep];
		generators= (Transpose@ Inverse@ L . # . Transpose@ L & /@ generators)// Simplify;
	];
	SparseArray@ Simplify[rot . generators]
	
]


Generators::arg= "Generators[alg, rep] takes a Lie algebra and a representation as arguments."
Generators@ ___:= (Message[Generators::arg]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Structure constants*)


(* ::Text:: *)
(*The structure constant of the algebra in the basis used for generators of all representations*)


StructureConstants[alg_]:= Block[{},
	TestAlg@ alg; 
	I Generators[alg, AdjointRepresentation@ alg]
];


StructureConstants::arg= "StructureConstants[alg] takes a Lie algebra as argument."
StructureConstants@ ___:= (Message[StructureConstants::arg]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Auxiliary functions*)


NormalizeGenerator::error1="Normalization could not be completed since the basis choice for the generators yields Tr[\!\(\*SuperscriptBox[\(T\), \(a\)]\)\!\(\*SuperscriptBox[\(T\), \(b\)]\)]=0 in some cases."


NormalizeGenerator[gen_,alg_,rep_]:= If[MatNorm[gen]=!=0, 
										Return[gen  Sqrt[DynkinIndex[alg,rep] / MatNorm[gen]]],
										Message[NormalizeGenerator::error1];
										Return[gen]
										];


Generator2Root[CartanGen_,mat_]:=Module[{},
		Return[MatInnerProd[MatCommutator[#,mat],mat]/ MatInnerProd[mat,mat] &/@CartanGen //Simplify]
]


(* ::Subsubsection::Closed:: *)
(*Checks on the generators*)


(* ::Text:: *)
(*Structure Constants obtained from explicit generator matrices*)


SpecificStructureConstants[generators_,a_,b_,c_]:=Block[{},
		-I MatInnerProd[generators[[c]],MatCommutator[generators[[a]] , generators[[b]]]]/ MatNorm[generators[[c]]] //Simplify
]
SpecificStructureConstants[generators_]:=Block[{groupdim},
		groupdim= Length[generators];
		SparseArray@ Table[SpecificStructureConstants[generators,a,b,c],{a,1,groupdim},{b,1,groupdim},{c,1,groupdim}] 
]


(* ::Text:: *)
(*Linear Independence check (to decide whether to add to the basis or not)*)


CheckLinearIndependence::setnotorthogonal="The input set is not linearly independent."


CheckLinearIndependence[set_List]:=Module[{rank,len},
	rank=MatrixRank@(Flatten/@Normal/@set);
	len=Length[set];
	If[len =!= rank, Message[CheckLinearIndependence::setnotorthogonal]; Return[False];, Return[True]];
]
CheckLinearIndependence[set_List,mat_]:=Module[{rank,newrank},
	rank=MatrixRank@(Flatten/@Normal/@set);
	If[CheckLinearIndependence[set],
		newrank=MatrixRank[Flatten/@Normal/@Append[set,mat]];
		If[newrank =!= rank, Return[True],Return[False]],
		Message[CheckLinearIndependence::setnotorthogonal];
		Abort[];
	]
]


(* ::Text:: *)
(*Check several conditions on the generators : Hermiticity, Dynkin Index, quadratic Casimir, Basis complete and linearly independent + Option to print the generators.*)


CheckGeneratorsConditions[alg_,rep_,OptionsPattern[{DisplayGenerators->False}]]:=Block[{generators=Generators[alg,rep],groupdim,BoolCheck},
		groupdim=Length[generators];
		Print["Representation ", rep, " of ", alg];
		
		BoolCheck=AllTrue[Flatten[Table[generators[[a]]-generators[[a]]\[ConjugateTranspose],{a,1,groupdim}]],# === 0 &];
		Print["Hermiticity condition T^a = (T^a)\[ConjugateTranspose] : ", BoolCheck];
		If[!BoolCheck,Print[{Table[generators[[a]]-generators[[a]]\[ConjugateTranspose],{a,1,groupdim}]}//MatrixForm] ];
		
		BoolCheck= (Table[Tr[generators[[a]] . generators[[b]]],{a,1,groupdim},{b,1,groupdim}] === DynkinIndex[alg,rep]IdentityMatrix[groupdim]);
		Print["Dynkin Index: Tr[\!\(\*SuperscriptBox[\(T\), \(a\)]\)\!\(\*SuperscriptBox[\(T\), \(b\)]\)] = ",DynkinIndex[alg,rep]," \!\(\*SuperscriptBox[\(\[Delta]\), \(ab\)]\) : ", BoolCheck];
		If[!BoolCheck,Print[Table[Tr[generators[[a]] . generators[[b]]],{a,1,groupdim},{b,1,groupdim}]]];
		
		BoolCheck= (Sum[generators[[a]] . generators[[a]],{a,1,groupdim}] === Casimir2[alg,rep]IdentityMatrix[RepresentationDimensionInternal[alg,rep]]);
		Print["Quadratic Casimir: (\!\(\*SuperscriptBox[\(T\), \(a\)]\)\!\(\*SubscriptBox[\()\), \(ik\)]\)(\!\(\*SuperscriptBox[\(T\), \(a\)]\)\!\(\*SubscriptBox[\()\), \(kj\)]\) = ",Casimir2[alg,rep]," \!\(\*SubscriptBox[\(\[Delta]\), \(ij\)]\) :",BoolCheck];
		If[!BoolCheck,Print[Sum[generators[[a]] . generators[[a]],{a,1,groupdim}]//MatrixForm] ];
		
		BoolCheck = CheckLinearIndependence[List@@generators] && GroupDimension[alg]=== groupdim;
		Print["Generator basis complete and linearly independent = ", BoolCheck];
		
		(*BoolCheck= (Count[Normal/@generators, _?DiagonalMatrixQ]==Length@CartanMatrix@alg);
		Print["Number of diagonal matrices is ",Length@CartanMatrix@alg," : " BoolCheck];*)
		
		BoolCheck = (Total[Abs/@Flatten[SpecificStructureConstants[generators]-SpecificStructureConstants[Generators[alg,SmallestRepresentation@alg]]]]===0);
		Print["Structure constants same as for the fundamental representation = ", BoolCheck];
				
		If[OptionValue[DisplayGenerators],
			Print["Generators:"];
			{generators} //MatrixForm]
]


(* ::Text:: *)
(*Also a condition on the complex conjugate representation, i.e. T^a(R* ) = - (T^a(R))**)


(* ::Subsection:: *)
(*Invariants *)


(* ::Subsubsection::Closed:: *)
(*Invariant tensors*)


Options@ InvariantTensors= {
	AntisymmetricIndices-> {},
	Normalization-> Default,
	SymmetricIndices-> {} 
	};


InvariantTensors::invldarg= "`1` is not a valid argument; a list of representations is expected."
InvariantTensors::nosing= "There are no singlets in the product of the representations."
InvariantTensors::syminp= "`1` is not a valid list (of lists) of enumerated indices to (anti)symmetrize."
InvariantTensors::impsym= "Representaitons of indices `1` are not identical and (anti)symmetrization is imposible."


InvariantTensors[alg_, reps_List, OptionsPattern[]]? OptionsCheck:= Block[{canonicalOrd, inds, inputOrd, inv, 
		invariants, n, perm, repDimensions, rot, syms, asyms, transp, ToSingleIndex, ToMultiIndex},
	(*Check arguments*)
	If[!MatchQ[reps, {(_List|CRep@ _List)..}], 
		Message[InvariantTensors::invldarg, reps];
		Abort[];
	];
	RepresentationCheck[alg, #]&/@ reps;
	
	(*Verify that the representations allow for a singlet*)
	If[Lookup[RepresentationProduct[alg, reps], Key@ ConstantArray[0, AlgebraRank@ alg], 0] === 0,
		Message[InvariantTensors::nosing, reps];
		Return@ {};
	];
	
	(*Check symmetrization is between identical representations*)
	syms= Symmetrizations[reps, OptionValue@ SymmetricIndices];
	asyms= Symmetrizations[reps, OptionValue@ AntisymmetricIndices];
	
	(*Cannonically sort the representations*)
	canonicalOrd= OrderingBy[reps, RepresentationDimensionInternal[alg, #]& ];
	inputOrd= Ordering@ canonicalOrd; 
	
	(*Determine invariants*)
	invariants= ComputeInvariants[alg, reps[[canonicalOrd]], Orthonormalize-> False];
	invariants= Transpose[#, canonicalOrd]&/@ invariants;

	(*Symmetrize the invariants. Normalization is fixed after.*)
	invariants= Table[
			Do[
				inv= Signature@ inds SparseArray@ Simplify@ Sum[
						transp= Range[Length@ reps]; transp[[inds]]= transp[[perm]]; 
						Signature@ perm Transpose[inv, transp]
					,{perm, Permutations@ inds}];
			,{inds, asyms}];
			Do[
				inv= SparseArray@ Simplify@ Sum[
						transp= Range[Length@ reps]; transp[[inds]]= transp[[perm]]; 
						Transpose[inv, transp]
					,{perm, Permutations@ inds}];
			,{inds, syms}];
			inv
		,{inv, invariants}];
	invariants= DeleteCases[SortBy[ArrayRules[#][[;; -2]]&/@ invariants, Length], _?(Length@ # === 0 &)];
	If[Length@ invariants=== 0, Return@ {}];
	
	(*Orthogonalize invariants*)
	{ToSingleIndex, ToMultiIndex}= TensorIndexMappings[RepresentationDimensionInternal[alg, #]&/@ reps];
	
	repDimensions= RepresentationDimensionInternal[alg, #]&/@ reps;
	invariants= MapAt[ToSingleIndex, invariants, {All, All, 1}];
	invariants= SparseArray[#, {Times@@ repDimensions}]&/@ invariants;
	
	invariants= ArrayRules[#][[;; -2]]&/@ SparseArray/@ GramSchmidtComplex@ invariants;
	invariants= DeleteCases[invariants, _?(Length@ # === 0 &)];
	invariants= MapAt[ToMultiIndex@ #[[1]]&, invariants, {All, All, 1}];
	invariants= SparseArray[#, repDimensions]&/@ invariants; 
	
	(*Apply rotations*)
	invariants= Table[
			Do[
				If[FSIndicator[alg, reps[[n]] ] < 1, Continue[]; ];
				rot= RealBasisRotation[alg, reps[[n]]]; 
				inv= Simplify@ TensorContract[Transpose[TensorProduct[rot, inv], 1 <-> n+ 2], {1, 2}];
			, {n, Length@ reps}];
			inv= SparseArray@ inv;
			(*Fix first element to have Arg 0, using U(1) ambiguity*)
			inv= Conjugate@ Sign@ ArrayRules[inv][[1, 2]] inv//Simplify;
			(*Normalization*)
			inv/TensorNorm@ inv * If[OptionValue@ Normalization === Default,
				Power[Times@@ repDimensions, 1/4], Sqrt@ OptionValue@ Normalization] 
		, {inv, invariants}];
	
	invariants 
];


InvariantTensors::arg= "InvariantTensors[alg, repList] takes a Lie algebra and a list of representations as arguments."
InvariantTensors@ ___:= (Message[InvariantTensors::arg]; Abort[];) 


TensorNorm@ T_:= Conjugate@ Flatten@ T . Flatten@ T// Sqrt// Simplify; 


(* ::Text:: *)
(*Ancillary function for checking the (anti)symmetrization option*)


Symmetrizations[reps_, optSym_]:= Block[{syms, len= Length@ reps},
	syms= Switch[optSym
		,{_Integer..},
			{optSym}
		,{{_Integer..}..}|{},
			optSym
		,_,
			Message[InvariantTensors::syminp, optSym];
			Abort[];
		];
	Do[
		If[Max@ sym > len || Min@ sym < 1,
			Message[InvariantTensors::syminp, optSym];
			Abort[];
		];
		If[Length@ DeleteDuplicates@ reps[[sym]] > 1,
			Message[InvariantTensors::impsym, sym];
			Abort[];
		];
	, {sym, syms}];	
	
	syms
]; 


(* ::Subsubsection::Closed:: *)
(*Compute tensors *)


(* ::Text:: *)
(*Determines a basis of invariant tensors in a product of representations [Fonseca, B2].*)


Options@ ComputeInvariants= {Orthonormalize-> True};


ComputeInvariants[alg_, reps_List, opt:OptionsPattern[]]:= ComputeInvariants[alg, reps, opt]= Block[
	{csMats, dimOffsets, ET, h, i, identities, inv, invariants, n, nullPositions, nullWeights, 
		repDims, repPosistion, repWeights, temp, ToMultiIndex, ToSingleIndex, ws, weightOffsets, 
		rank= Length@ CartanMatrix@ alg},
	
	repWeights= RepresentationWeights[alg, #]&/@ reps;
	repDims= Total/@ repWeights;
	weightOffsets= WeightPositions[alg, #]&/@ reps;
	{ToSingleIndex, ToMultiIndex}= TensorIndexMappings@ repDims;	
	
	(*Determine all components of the tensor product with weight 0*)
	nullWeights= Cases[Tuples[Keys/@ repWeights], _? (MatchQ[Total@ #, {0..}] &)];
	nullWeights= Association@@((#-> Inner[Function[{x, y}, x@ y], h@@ repWeights, h@@ #, List])&/@ nullWeights);
	(*The position of elements in the tensor product with weight 0*)
	nullPositions= ToSingleIndex/@ Flatten[Table[
			ConstantArray[Inner[(#1@ #2&), h@@ weightOffsets, h@@ ws, List], Times@@ nullWeights@ ws]
			+ Tuples[Range/@ nullWeights@ ws]
		, {ws, Keys@ nullWeights}], 1];
	
	(*Determine the raising matrices in the reducible product representation*)
	csMats= CSMatrices[alg, #]&/@ reps;
	identities= Table[SparseArray[{n_, n_}->1, {i, i}], {i, repDims}];
	ET= Table[Sum[temp= identities; temp[[n]]= Transpose@ csMats[[n, i, 2]];
			TensorProductMatrix@ temp
		,{n, Length@ reps}] ,{i, rank}];
	
	(*Find the null space corresponding to the invariant tensors*)
	ET= (Join@@ ET)[[;;, nullPositions]];
	n= DeleteDuplicates@ ArrayRules[ET][[;; -2, 1, 1]];
	invariants= NullSpace@ ET[[n]];
	If[OptionValue@ Orthonormalize,  
		invariants= (*Orthogonalize@*) GramSchmidtComplex@ invariants * Sqrt@ Max@ repDims;
	];
	
	(*Convert back from nullspace vectors in the 0-weights to tensors in the product*)
	Table[
		SparseArray[MapAt[ToMultiIndex@ nullPositions[[#[[1]] ]]&, ArrayRules[inv][[;; -2]], {All, 1}], 
			repDims] 
	,{inv, invariants}] 
];


(* ::Subsubsection::Closed:: *)
(*Auxilliary *)


(* ::Text:: *)
(*Makes a tensor product and then reshapes it into a matrix *)


TensorProductMatrix[matrices_List]:= Block[{prod, dimOffsets, rules, f},
	prod= TensorProduct@@ matrices;
	dimOffsets= FoldList[Times, 1, Length/@ matrices];
	f= (1+ dimOffsets[[;; -2]] . (# -1) &);
	
	rules= ArrayRules[prod][[;; -2]];
	rules= MapAt[({f@ #[[1;; ;; 2]], f@ #[[2;; ;; 2]]} &), rules, {All, 1}];
	Length@ rules;
	SparseArray[rules, dimOffsets[[-1]] {1, 1}]	
];


(* ::Text:: *)
(*Returns functions for mapping the coordinates of a multi-index tensor to a vector and back *)


TensorIndexMappings[dims_List]:= With[{
		dimOffsets= FoldList[Times, 1, dims][[;;-2]]
	},
	{(1+ dimOffsets . (# -1) &), (*To single index*)
	Function[{ind}, FoldPairList[{Mod[#1, #2], Floor[#1/ #2]}&, ind-1, dims]+ 1] (*To multi-index*)
	}
];


(* ::Subsection:: *)
(*Vector-space rotations *)


(* ::Text:: *)
(*To rotate the basis of a real representation to a basis where the identity matrix is invariant. *)


RealBasisRotation::unrlrep= "The representation `1` of `2` is not real.";


RealBasisRotation[alg_, \[CapitalLambda]_]:= RealBasisRotation[alg, \[CapitalLambda]]= Block[{},
	RepresentationCheck[alg, \[CapitalLambda]];
	If[FSIndicator[alg, \[CapitalLambda]]< 1, 
		Message[RealBasisRotation::unrlrep, \[CapitalLambda], alg];
		Return@ SparseArray[{{i_, i_}-> 1}, RepresentationDimensionInternal[alg, \[CapitalLambda]] {1, 1}];
	];
	
	If[\[CapitalLambda] === AdjointRepresentation@ alg,
		Return@ BasisRotation@ alg
	];
	
	ConjugateTranspose@ TakagiFactorization@ First@ ComputeInvariants[alg, {\[CapitalLambda], \[CapitalLambda]}]
 ]; 


(* ::Subsection:: *)
(*Other functions*)


(* ::Text:: *)
(*Gives the first position of the weights in the weight vectors of the representation*)


WeightPositions[alg_, \[CapitalLambda]_]:= Block[{i= 0},
	((i+= #)- #&)/@ RepresentationWeights[alg, \[CapitalLambda]]
];


(* ::Section:: *)
(*Auxiliary function*)


(* ::Subsection:: *)
(*Matrix functions*)


(* ::Subsubsection::Closed:: *)
(*Matrix commutators and inner product *)


MatCommutator[a_,b_]:=a . b-b . a;
MatNorm[A_]:= MatInnerProd[A,A];
MatInnerProd[A_, B_]:= Tr[ConjugateTranspose@ A . B];  


(* ::Subsubsection::Closed:: *)
(*CholekyDecomposition*)


(* ::Text:: *)
(*Cholesky-like decomposition for finding lower triangular L given symmetric A such that A=L L\[Transpose].*)


CholeskyLikeDecomposition[A_]:= Block[{L, dim= Length@ A, i ,j},
	(*Cholesky\[Dash]Banachiewicz, algorithm *)
	L= ConstantArray[0, {dim, dim}];
	Do[
		If[i=== j,
			L[[j, j]]= Simplify@ Sqrt[A[[j, j]]- L[[j]] . L[[j]] ];
		,
			If[L[[j, j]]=== 0, Continue[];];
			L[[i, j]]= (A[[i, j]]- L[[i]] . L[[j]])/ L[[j, j]]// Simplify;
		];	
	,{i, dim}, {j, i}];
	L
];


(* ::Subsubsection::Closed:: *)
(*Gram-Schmidt reduction*)


(* ::Text:: *)
(*Gram-Schmidt procedure for orthogonalizing vectors*)


GramSchmidt[A_]:= Block[{dim= Dimensions@ A, i, out, proj, temp},
	out= SparseArray@ Simplify@ A;
	out[[1]]= Normalize@ out[[1]];
	Do[
		temp= out[[i]];
		out[[i]]= temp - (temp . Transpose@ out[[;; i-1]]) . out[[;; i-1]]// Simplify;
		temp= out[[i]] . out[[i]]//Simplify;
		If[temp> 0, out[[i]]/= Sqrt@ temp];
	,{i, 2, dim[[1]]}]; 
	out
];


GramSchmidtComplex[A_]:= Block[{dim= Dimensions@ A, i, out, proj, temp},
	out= Simplify@ A;
	out[[1]]= Normalize@ out[[1]];
	Do[
		temp= out[[i]];
		out[[i]]= temp - (temp . ConjugateTranspose@ out[[;; i-1]]) . out[[;; i-1]]// Simplify;
		temp= out[[i]] . Conjugate@ out[[i]]//Simplify;
		If[temp> 0, out[[i]]/= Sqrt@ temp];
	,{i, 2, dim[[1]]}]; 
	out
];


(* ::Text:: *)
(*Does a Gram-Schmidt transformation on a basis to orthonormalize it given a metric.*)


GramSchmidtMetric@ met_:= Block[{temp, proj, dim},
	dim= Length@ met;
	proj= 0 met; 
	Table[
		temp= UnitVector[dim, i];
		temp-= proj . met . temp;
		temp/= Sqrt[Conjugate@ temp . met . temp];
		proj+= TensorProduct[temp, Conjugate@temp];
		temp
	, {i, dim}]
];


(* ::Subsubsection::Closed:: *)
(*Takagi Factorization*)


(* ::Text:: *)
(*Takagi factorization of matrix A=A\[Transpose] is A = U. \[CapitalSigma]. U\[Transpose], where U is unitary and \[CapitalSigma] diagonal with non-negative entries. This method returns U*)


TakagiFactorization::nonSym= "Takagi factorization can only be performed on Symmetric matrices."
TakagiFactorization::method= "An appropriate method for Takagi factorizing the matrix has not been implemented."


TakagiFactorization[A_]:= Block[{Ap, block, blocks, dim, rules, U},
	If[! SymmetricMatrixQ@ A, 
		Message[TakagiFactorization::nonSym];
		Abort[];
	];
	
	(*Takagi factorize blocks individually*)
	dim = Length@ A; U= SparseArray[{}, {dim, dim}];
	blocks= MatrixBlocks@ A;
	Do[
		Ap= A[[block, block]];
		If[UnitaryMatrixQ@ Ap,
			U[[block, block]]= UnitaryTakagi@ Ap;
			Continue[];
		];		
		Message[TakagiFactorization::method];
		Abort[];
	,{block, blocks}];
	
	U
	(*dim = Length@ A;
	rules= ArrayRules[A][[;;-2]];
	(*If off-diagonal*)
	If[And@@ (Total@ First@ # === dim+ 1 &/@ rules),
		Return@ TakagiMethod2x2@ A;
	]; 
	Message[TakagiFactorization::method];
	Abort[];
	*)	
];


(* ::Text:: *)
(*Method for rotating away the 2x2 blocks *)


TakagiMethod2x2@ A_:= Block[{dim= Length@ A, i, j, rules, temp},
	rules= Flatten@ Table[
		temp= A[[i, j= dim+ 1- i]];
		If[temp =!= 0,
			temp= If[temp> 0, 1, I]/ Sqrt@ 2;
			{{i, i}-> I temp, {i, j}-> temp,
				{j ,i}-> -I temp, {j, j}-> temp} 
		,
			{{i, i}-> 1, {j, j}-> 1}
		]  
	, {i, Floor[dim/ 2]}];
	If[OddQ@ dim, 
		i= (dim +1)/ 2;
		AppendTo[rules, {i, i}-> If[A[[i, i]]< 0, I, 1]];
	];
	SparseArray[rules, {dim, dim}]	
];


(* ::Text:: *)
(*Takagi factorization for unitary, symmetric A: http://link.springer.com/article/10.1134/S0965542512010034*)


UnitaryTakagi@ A_:= Block[{dim= Length@ A, An= A, n, U, x, y, temp, mu, V},
	U= IdentityMatrix@ dim;
	Do[
		x= UnitVector[dim+ 1- n, 1];
		temp= An . Conjugate@ x; mu= Simplify[Conjugate@ x . temp]; 
		If[Abs@ mu === 1,
			y= x/ Sqrt@ mu;
		,
			y= temp + x; y/= Sqrt[Conjugate@ y . y];
		];
		y= Simplify@ y;
		V= Transpose@ ExtendToOrthonormalBasis@ y;
		An= ConjugateTranspose[V][[2;;]] . An . Conjugate[V][[;;, 2;;]]// Simplify;
		U[[;;, n;;]]= U[[;;, n;;]] . V; 
	, {n, dim- 1}];
	
	U[[;;, -1]]*= Sqrt@ An[[1, 1]]/ Sqrt@ Abs@ An[[1, 1]]// Simplify;
	U
];


(*ASqrt@ Complex[a_, b_]:= Block[{c},
	c= Sign[b]/Sqrt[2]*Sqrt[a+ Sqrt[a^2+ b^2]];
	c+ b I/ (2 c)
];*)
(*ASqrt@ z_:= Sqrt@ Abs@ z Exp[I/2 Arg@ z]; *)


(* ::Text:: *)
(*Extend vector to Orthonormal basis*)


ExtendToOrthonormalBasis@ v_:= Block[{B},
	B= Join[{v}, NullSpace@ {v}];
	GramSchmidtComplex@ B 
];


(* ::Text:: *)
(*Determine the blocks of a symmetric matrix *)


MatrixBlocks[A_]:= Block[{block, entries, i, j, mat, taken, queue, qTake, qPut, dim= Length@ A},
	entries= Table[
		Flatten[First/@ ArrayRules[A[[i]]][[;; -2]]]
	, {i, dim}];
	
	taken= ConstantArray[False, dim];
	queue= ConstantArray[0, dim+ 1]; qPut= qTake= 1; 
	block= 1;
	Last@ Reap@ Do[
		If[taken[[i]], Continue[];];
		queue[[qPut++]]= i; taken[[i]]= True;
		While[queue[[qTake]] =!= 0,
			Do[
				If[!taken[[j]], 
					queue[[qPut++]]= j;
					taken[[j]]= True;
				];
			, {j, entries[[queue[[qTake]] ]]}];
			Sow[queue[[qTake++]], block]; 
		];
		block++;
	, {i, dim}]
] 
