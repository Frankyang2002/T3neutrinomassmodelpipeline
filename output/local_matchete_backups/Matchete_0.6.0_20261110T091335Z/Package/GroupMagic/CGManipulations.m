(* ::Package:: *)

Package["Matchete`"] 


(* ::Title:: *)
(*GroupMagic`CGManipulations`*)


(* ::Subtitle:: *)
(*Implements the CG coefficients and their contractions*)


(* ::Section:: *)
(*Scoping & usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["CG"]
PackageExport["DefineCG"]
PackageExport["RenameCG"]
PackageExport["RemoveCG"]
PackageExport["RemoveRepresentation"]
PackageExport["RemoveGroup"]
PackageExport["ClearGroups"]
PackageExport["DefineCompositeCG"]
PackageExport["ContractCGs"]
PackageExport["DefineCGGroup"]
PackageExport["DefineRepresentation"]
PackageExport["GetCGTensor"]


PackageExport["GetGroups"]
PackageExport["GetRepresentations"]


PackageExport["ReplaceCGs"]
PackageExport["CGsToReplace"]


PackageExport["del"]


PackageExport["Dimension"]
PackageExport["DynkinCoefficients"]
PackageExport["Group"]
PackageExport["Reality"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["$Groups"]
PackageScope["$Representations"]
PackageScope["$CGproperties"]
PackageScope["$CGreplacements"]


PackageScope["InBasis"]
PackageScope["UniqueConj"]
PackageScope["DeltaDecomposable"]


(* ::Subsubsection::Closed:: *)
(*Definitions*)


DefineCGGroup::usage=
	"DefineCGGroup[grName, alg] defines a new group with which to create and manipulate CG objects. It is referenced by the name grName and its Lie algebra determined by alg.";


DefineRepresentation::usage=
	"DefineRepresentation[repName, groupName, DynkinCoef] defines a representation with reference name repName of the group and the specified Dynkin coefficients.";


CG::usage= 
	"CG[symb, {indices}] denotes the Clebsch-Gordan coefficient symb, with indices given in the list.";


DefineCG::usage= 
	"DefineCG[symb, {index types}, tensor] defines a new Clebsch-Gordan coefficient, symb, with the specified index types (ordered). tensor specifies the explicit tensor assignment of the CG.";


RenameCG::usage=
	"RenameCG[oldName-> newName] assigns a new symbol to an already defined CG."


RemoveCG::usage=
	"RemoveCG[symb] removes the CG of symb.";


RemoveRepresentation::usage=
	"RemoveRepresentation[symb] removes the corresponding representation and all associated CGs.";


RemoveGroup::usage=
	"RemoveGroup[symb] removes the group and all associated representations and CGs.";


DefineCompositeCG::usage=
	"DefineCompositeCG[symb, {cgs}, {{indices}, {indices}, ...}, defines a new CG from the contraction of known cgs. Repeated indices describe determine how the contraction is done (as per Einstein sumation). Open indices give determine the indices of the new CG (in order). The option ReplaceCGs-> True (default) creates an automatic replacement for this contraction used in ReplaceCGs. ";


ContractCGs::usage= 
	"ContractCGs[expr] evaluates all contractions of CGs occuring in the expression.";


ClearGroups::usage=
	"ClearGroups[] clears all definitions of CGs, groups, and representations.";


GetCGTensor::usage=
	"CGtensor[symb] returns the tensor associated with the CG symbol.";


GetGroups::usage= "GetGroups[] returns an association of all the groups that have been defined by DefineCGGroup.";
GetRepresentations::usage= "GetRepresentations[group] returns an association of all the representations of the group that have been defined by DefineRepresentation. If group is nor specified, then it returns all groups and representations defined.";


ReplaceCGs::usage= 
	"ReplaceCGs[expr] replaces all the automatically generated CGs with their constituent CGs. The option CGsToReplace-> All can also take a list of Specific CGs to replace.";
CGsToReplace::usage= "Option used in ReplaceCGs.";


del::usage     = "Standard head for delta Clebsch-Gordan coefficients. E.g., CG[del@ rep, {a, b}] for \!\(\*SubscriptBox[\(\[Delta]\), \(ab\)]\).";


Dimension::usage= 
	"Dimension is a property dictating the size of a representation."
DynkinCoefficients::usage= 
	"DynkinCoefficients is a property of a representation, the highest weight associated with it."
Reality::usage= 
	"Reality is a property of a representation determining if it is Real, Pseudoreal, or Complex."


Group::usage= 
	"Group is a property indicating the name of the associated Lie group.";


(* ::Section:: *)
(*Tensor contractions *)


(* ::Subsection:: *)
(*Clebsch-Gordan definitions *)


(* ::Subsubsection::Closed:: *)
(*Representations *)


(* ::Text:: *)
(*Get the group and representation properties*)


GetGroups[]:= $Groups;


GetRepresentations[args___]:= 
	(ReplaceAt[{1-> "Real", 0-> "Complex", -1-> "PseudoReal"}, Key@Reality]/@ $Representations)[args];


(* ::Text:: *)
(*For assigning a new group to the current instance*)


DefineCGGroup::args= "DefineCGGroup takes arguments group name (Symbol) and a Lie algebra."


DefineCGGroup[groupName_Symbol, alg_]:= Module[{},
	TestAlg@ alg;
	AppendTo[$Groups, groupName-> alg];  
 ]; 
DefineCGGroup@ ___:= (Message[DefineCGGroup::args]; Abort[];) 


(* ::Text:: *)
(*For assigning a new representation to the current instance*)


DefineRepresentation::unkwngroup= "`1` has not been defined as a group."
DefineRepresentation::repinuse= "`1` is already used to denote a representation of `2`."
DefineRepresentation::args= "Invalid arguments."


Options[DefineRepresentation]= {IndexAlphabet-> None};


DefineRepresentation[repName_, groupName_Symbol, dynkCoef_, opts: OptionsPattern[]] ? OptionsCheck:=
Module[{fsIndicator, lieAlg},
	If[KeyFreeQ[$Groups, groupName], 
		Message[DefineRepresentation::unkwngroup, groupName];
		Abort[];
	];
	RepresentationCheck[$Groups@ groupName, dynkCoef];
	
	(*If repName is already used for a representation, abort or do nothing*)
	If[!KeyFreeQ[$Representations, repName],
		If[{groupName, dynkCoef} === Lookup[$Representations@ repName, {Group, DynkinCoefficients}],
			(* If provided, overwrite the index alphabet for the representation *)
			If[OptionValue@ IndexAlphabet =!= None,
				AppendTo[$IndexAlphabets, repName-> BuildIndexAssoc@ OptionValue@ IndexAlphabet]
			];
			Return[];
		,
			Message[DefineRepresentation::repinuse, repName, $Representations[repName, Group]];
			Abort[];		
		];
	];
	
	lieAlg= $Groups@ groupName;
	fsIndicator= FSIndicator[lieAlg, dynkCoef];
	
	(*Set up the representation*)
	AppendTo[$Representations, repName-> <|
		Group-> groupName, 
		DynkinCoefficients-> dynkCoef, 
		Dimension-> RepresentationDimension[$Groups@ groupName, dynkCoef],
		Reality-> fsIndicator
		|>];
	If[$Representations[repName, Reality]=== 1, 
		Bar@ ind:Index[_, repName]:= ind;
	,
		SetBarable@ repName;
	];
	
	(*Add delta CG to the representation*)
	If[fsIndicator=== 1,
		DefineCG[del@ repName, {repName, repName}, 
			First@ InvariantTensors[lieAlg, {dynkCoef, dynkCoef}, SymmetricIndices-> {1, 2}] ];
	,
		DefineCG[del@ repName, {repName, Bar@ repName}, 
			First@InvariantTensors[lieAlg, {dynkCoef, CRep@ dynkCoef}] ];
	];
	
	(*Incorporation with Matchete*)
	If[MemberQ[Join[Keys@ $GlobalGroups, Keys@ $GaugeGroups], groupName],
		(* Define generator for the representation if the group is a gauge or global group -- these always define an adjoint represenetation*)
		If[fsIndicator === 1,
			DefineCG[gen@ repName, {groupName@ adj, repName, repName}, Generators[lieAlg, dynkCoef]];
		,
			DefineCG[gen@ repName, {groupName@ adj, repName, Bar@repName}, Generators[lieAlg, dynkCoef]];
		];
	];
	
	(* If provided, add index alphabet for the representation *)
	If[OptionValue@ IndexAlphabet =!= None,
		AppendTo[$IndexAlphabets, repName-> BuildIndexAssoc@ OptionValue@ IndexAlphabet]
	];
	
	(*Rules of contraction of deltas to improve efficiency*)
	Delta /: Delta[Index[a_,repName],Index[b_,repName]] * Delta[Index[b_,repName],Index[c_,repName]] :=
		Delta[Index[a,repName],Index[c,repName]];
	Delta /: Delta[Index[a_,repName],Index[b_,repName]]^2 :=
		Delta[Index[a,repName],Index[a,repName]];
	Delta[Index[a_,repName],Index[a_,repName]] := DimRep[repName];
];
DefineRepresentation@ ___:= (Message[DefineRepresentation::args]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Properties*)


(* ::Text:: *)
(*Conjugation of CGs *)


(*This is a tentative function for the action of the Bar on the CG, it needs more extensive checks.*)


TensorsOverlap[ten1_,ten2_]:=CGProduct[ten1, ten2]/Sqrt[CGProduct[ten1, ten1] CGProduct[ten2, ten2]];


Bar@ cg_CG:= CGBarInternal@ cg;


CGBarInternal@ CG[Bar@ symb_, inds_]:= CG[symb, Bar/@ inds]; (*cg is necessarily barable*)


CGBarInternal@ CG[symb_/; $CGproperties[symb, UniqueConj], inds_]:= CG[Bar@ symb, Bar/@ inds];


CGBarInternal@ CG[symb_, inds_]:= Block[
		{reps, tens, cReps, cTens, p, perms, coeff, pout, coeffout},
	reps= $CGproperties[symb, Indices];
	tens= $CGtensors[symb];
	cTens= Conjugate@ tens;
	cReps= Bar/@ reps; 
	perms= DegeneratePerms[cReps, reps];
	(*Check individual permutations*)
	{coeffout, pout}= Catch[
		Do[
			coeff= TensorsOverlap[tens, PermuteArray[cTens, p]];
			If[Abs@coeff === 1, Throw@ {coeff, p};]; 
		,{p, perms}];
	];
	coeffout CG[symb, Bar/@inds[[pout]]]
];


(* ::Text:: *)
(*CG shorthand input *)


CG::unkwn= "`1` is not a known CG. Verify that this CG has been loaded into the kernel.";
CG::inds= "The CG `1` has `2` indices.";


CG[symb_, inds: {_Symbol..}]:= Block[{indexTypes},
	If[!MemberQ[Keys@ $CGtensors, symb], 
		Message[CG::unkwn, symb]; Abort[]; ];
	indexTypes= $CGproperties[symb, Indices];
	If[Length@ inds =!= Length@ indexTypes,
		Message[CG::inds, symb, Length@ indexTypes]; Abort[]; ];
	CG[symb, MapThread[Index, {inds, indexTypes}] ]
];


CG[Bar@ symb_, inds: {_Symbol..}]:= Bar@ CG[symb, inds]


CG::args= "CG should have 2 arguments.";
CG[_, _, __]:= (Message[CG::args]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*Utility functions *)


(* ::Text:: *)
(*The inner product from fully contracting two CGs of the same type*)


CGProduct[cg1_, cg2_]:= Conjugate@ Flatten[cg1/. $CGtensors] . Flatten[cg2/. $CGtensors]// Simplify; 


(*Floating point sparse arrays*)
CGNumProduct[cg1_SparseArray, cg2_SparseArray]:= Conjugate@ Flatten[cg1] . Flatten[cg2]; 


(*w/o conjugation*)
SAProduct[cg1_SparseArray, cg2_SparseArray]:= Flatten[cg1] . Flatten[cg2]// Simplify


(* ::Text:: *)
(*Permutations*)


ComposePerm[\[Sigma]1_List, \[Sigma]2_List]:= \[Sigma]2[[\[Sigma]1]]; (*Gives \[Sigma]1\[SmallCircle]\[Sigma]2 *)


InversePerm[\[Sigma]_List]:= Ordering@ \[Sigma];


FindPerm[set_List]:= Ordering@ set;
FindPerm[seti_List, setf_List]:= ComposePerm[InversePerm@ FindPerm@ setf, FindPerm@ seti];


(* ::Text:: *)
(*For taking e.g. T_{a, b, c, d} -> T_{c, b, d, a} when \[Sigma]{a, b, c, d} = {c, b, d, a} *)


PermuteArray[tens_, \[Sigma]_List]:= Transpose[tens, InversePerm@ \[Sigma]];


(* ::Text:: *)
(*All permutations under which the list is unchanged *)


DegeneratePerms@ set_List:= Block[{groupings, perms},
	groupings= List@@ PositionIndex@ set;
	perms= Flatten/@ Tuples[Permutations/@ groupings];
	perms[[;;, InversePerm@ Flatten@ groupings]] 	
 ]; 


DegeneratePerms[set_List, targetSet_List]:= 
	ComposePerm[FindPerm[set, targetSet], #]&/@ DegeneratePerms@ set;


(* ::Text:: *)
(*Factor out symmetry subgroup from symmetry group (assuming full groups)*)


FactorOutSymmetrySubgroup[group_List, subGroup_List]:= 
	DeleteDuplicatesBy[group, Function[{x}, First@Sort[x[[#]]&/@ subGroup]]];


(* ::Text:: *)
(*Determine if two tensors are linearly dependent*)


(*LinDepQ[ten1_, ten2_]:= Abs[CGProduct[ten1, ten2]]^2/(CGProduct[ten1, ten1] CGProduct[ten2, ten2]) === 1;*)
(*Floating poit evaluation for speed up*)
LinDepQ[ten1_SparseArray, ten2_SparseArray]:= Abs[CGNumProduct[ten1, ten2]]^2 === CGNumProduct[ten1, ten1] CGNumProduct[ten2, ten2];


(* ::Text:: *)
(*Determines if the conjugate of the CG is linearly independent from itself (with permutation of indices)*)


UniqueConjQ[reps_, tens_SparseArray]:= Block[{cReps, cTens= Conjugate@ tens, p, perms},
	cReps= Bar/@ reps; 
	If[Sort@cReps =!= Sort@reps, Return@ True;];
	perms= DegeneratePerms[cReps, reps];
	Catch[
		Do[
			If[LinDepQ[tens, PermuteArray[cTens, p]],
				Throw@ False;]; 
		,{p, perms}];
	True
	]
];


(* ::Text:: *)
(*Determine whether a list is a subset of another counting multiplicity  *)


SubsetMultQ[subset_, set_]:= AllTrue[Merge[{Association@@ Rule@@@ Tally@ set, 
	- Association@@ Rule@@@ Tally@ subset}, Total], NonNegative];
	(*And@@ NonNegative/@ Merge[{Association@@ Rule@@@ Tally@ set, 
	- Association@@ Rule@@@ Tally@ subset}, Total];*)


(* ::Text:: *)
(*Test if an explicit tensor, with indices in representation set, is invariant*)


TensorInvariantQ[indexTypes_, tensor_]:= Module[{alg, rep, reps, generators, 
		n= Length@ indexTypes},
	alg= $Groups@ $Representations[First@ indexTypes/. Bar-> Identity, Group];
	(*Find the generators of associated to the representation of all external indicies*)
	reps= Table[
			If[Head@ rep === Bar, 
				CRep@ $Representations[First@ rep, DynkinCoefficients],
				$Representations[rep, DynkinCoefficients]
			]
		,{rep, indexTypes}];
	generators= Generators[alg, #]&/@ reps;
	
	(*Determines infinitesimal action of the algebra on the tensor T^{ai_k}_j CG^{i_1...j...i_n}, 
		and checks if it vanishes*)
	SAZeroQ@ Sum[
		EinsteinSummation[{generators[[i]], tensor},
			{{n+2, i, n+1}, Join[Range[i -1], {n+1}, Range[i+1, n]]}-> Prepend[Range@ n, n+2] ]
	, {i, n}] 
];


(* ::Text:: *)
(*Get all the elements of an association that satisfy certain properties*)


GetByProperty[assoc_Association, props_Association]:= GetCGsByProperty[assoc, Normal@ props];
GetByProperty[assoc_Association, props:_List|_Rule]:=
	Keys@ Select[assoc, MatchQ[#, KeyValuePattern[props]]&];


(* ::Subsection:: *)
(*Defining CG tensors*)


(* ::Text:: *)
(*List of explicit CG tensors and their properties properties *)


$CGtensors= <||>; (*CG tensors, symbolic*)
$CGNumTensors= <||>; (*CG tensors, floats*)
$CGproperties= <||>;


(* ::Subsubsection::Closed:: *)
(*Definition from tensor*)


(* ::Text:: *)
(*For  defining new CGs *)


DefineCG::args= "Invalid Arguments.";
DefineCG::cginuse= "`1` is already used to denote a cg of `2`."
DefineCG::invldIndType= "`1` is not defined as a representation.";
DefineCG::notTensor= "`1` is not an explicit tensor.";
DefineCG::notSameGroup= "Not all indices are representations of the same group.";
DefineCG::notInvariant= "The tensor is not invariant under group transformations.";
DefineCG::dimensions= "The dimensions of the indexTypes do not match those of the tensor representation.";


DefineCG[symb: _Symbol | _Symbol[_], indexTypes_List, tensorIn_]:= Module[
		{decomposable, temp, inBasis, indices, rep, reps, tensor, i},
	(*If symb is already used for a cg, abort or do nothing*)
	If[!KeyFreeQ[$CGproperties, symb],
		If[indexTypes === $CGproperties[symb, Indices] && SAZeroQ[tensorIn- $CGtensors@ symb],
			Return[];
		,
			Message[DefineCG::cginuse, symb, $CGproperties[symb, Group]];
			Abort[];		
		];
	];
	
	reps= indexTypes/. {Bar@ x_:> x};
	Do[
		If[KeyFreeQ[$Representations, rep],
			Message[DefineCG::invldIndType, rep];
			Abort[];
		];
	, {rep, reps}]; 
	
	If[Length@ DeleteDuplicates[$Representations[#, Group]&/@ reps] > 1,
		Message[DefineCG::notSameGroup];
		Abort[];
	];
	
	If[!ArrayQ@ tensorIn,
		Message[DefineCG::notTensor, tensorIn];
		Abort[];
	];
	tensor= SparseArray@ tensorIn;
	
	If[Dimensions@ tensor =!= Lookup[Lookup[$Representations, reps], Dimension],
		Message[DefineCG::dimensions];
		Abort[];
	];
	
	If[!TensorInvariantQ[indexTypes, tensorIn],
		Message[DefineCG::notInvariant];
		Abort[];
	];
	
	(*We determine wether the CG redundant: if it can be decomposed using other CGs,
		it does not need to be included in the basis.*)
	indices= Table[Index[Unique@ u, rep], {rep, indexTypes}];
	If[(inBasis = FindCGMatch[tensor, indices] === None), 
		ResetCGProjectors[];
	];
	
	(*Determine whether cg\[CircleTimes] cgConj should be included in tensor basis*)
	decomposable= If[Length@ indexTypes === 2,
			DeltaDecomposableQ[indexTypes, tensor]
		,
			False
		];
	
	$CGtensors@ symb= tensor;
	$CGNumTensors@ symb= N@ $CGtensors@ symb;
	$CGproperties@ symb= <|
			Group-> $Representations[First@ reps, Group], 
			Indices-> indexTypes, 
			UniqueConj-> UniqueConjQ[indexTypes, tensor],
			Symmetries-> FindTensorSymmetries[symb/. $CGNumTensors, indexTypes], 
			Real-> Element[Normal@ tensor, Reals] === True,
			InBasis-> inBasis,
			DeltaDecomposable-> decomposable
		|>; 
	
	(*Symbol is complex unless tensor and indices are real*)
	If[!($CGproperties[symb, Real] && (Bar/@ indexTypes === indexTypes) ), 
		SetBarable@ symb;
		$CGtensors@ Bar@ symb= Conjugate@ tensor;
	];		
	
	(*Create CG shorthand*)
	symb[inds__]:= CG[symb, List@ inds];
];
DefineCG@ ___:= (Message[DefineCG::args]; Abort[];) 


(* ::Text:: *)
(*Checks whether the tensor product of a CG and its conjugate is decomposable in terms of deltas*)


DeltaDecomposableQ[indexTypes_, tensor_]:= Block[{indProd, match},
	indProd= Table[Index[Unique@ u, rep], {rep, Join[indexTypes, Bar/@ indexTypes]}];
	match= Expand@ FindCGMatch[TensorProduct[tensor, Conjugate@ tensor], indProd];
	
	If[match === None, Return@ False; ];
	If[Head@ match === Plus, match= First@ match; ];
	If[Count[match, _CG, {0, Infinity}] === Length@ indexTypes, True, False]
]


(* ::Subsubsection::Closed:: *)
(*Renaming*)


(* ::Text:: *)
(*Assigns a new name to a previously (potentially automatically) defined CG*)


RenameCG::unkwn= "The CG `1` has not been definined.";
RenameCG::inuse= "The name `1` is already in use.";
RenameCG[oldName_-> newName_]:= Module[{},
	If[!MemberQ[Keys@ $CGtensors, oldName],
		Message[RenameCG::unkwn, oldName];
		Abort[];
	];
	If[MemberQ[Keys@ $CGtensors, newName],
		Message[RenameCG::inuse, newName];
		Abort[];
	];
	
	$CGtensors= KeyRename[$CGtensors, oldName-> newName];
	$CGNumTensors= KeyRename[$CGNumTensors, oldName-> newName];
	$CGproperties= KeyRename[$CGproperties, oldName-> newName];
	(*Symbol is complex unless tensor and indices are real*)
	If[$CGproperties[newName, UniqueConj], 
		SetBarable@ newName;
		(*barred cgs*)
		$CGtensors= KeyRename[$CGtensors, Bar@ oldName-> Bar@ newName];
	];
	
	ResetCGProjectors[];
	(*Change replacement rules if any*)
	If[KeyExistsQ[$CGreplacements, oldName],
		$CGreplacements@ oldName= $CGreplacements@ oldName/. 
			HoldPattern@ RuleDelayed[CG[x:(oldName|Bar@ oldName), inds_], rhs_]:> 
				RuleDelayed[CG[x/. oldName-> newName, inds], rhs];
		$CGreplacements= KeyRename[$CGreplacements, oldName-> newName];
	];
	
	If[$CGproperties[newName, UniqueConj], 
		UnsetBarable@ oldName;
	];
	
	(*Fix CG shorthands*)
	newName[inds__]:= CG[newName, List@ inds];
	Switch[oldName
	, _Symbol,
		Clear@ oldName
	, _Symbol[_],
		SubValues[Evaluate@ Head@ oldName]= DeleteCases[SubValues[Evaluate@ Head@ oldName], 
			RuleDelayed[HoldPattern[Verbatim[HoldPattern][oldName[_] ] ], _] ];
	];
];


KeyRename[association_, old_-> new_]/; KeyExistsQ[association, old]:= 
	KeyDrop[old]@ Insert[association, new-> association@ old, Key@ old];


(* ::Subsubsection::Closed:: *)
(*Composite CG*)


(* ::Text:: *)
(*Define a composite CG made from a contraction of other CGs*)


DefineCompositeCG::args= "Invalid Arguments.";
DefineCompositeCG::unkwnCG= "`1` is not an assigned CG.";
DefineCompositeCG::indCount= "The number of specified indices do not match the CGs"; 
DefineCompositeCG::indMult= "An index occurs more than twice.";
DefineCompositeCG::invldContraction= "The two indices with label `1` are incompatible types.";


Options@ DefineCompositeCG= {ReplaceCGs-> True};


DefineCompositeCG[name_, cgs_List, indicesIn:{_List..}, OptionsPattern[]]:= Block[
	{indexTypes, indices, eleminate, cg, tensor, contractions, ind, inds},
	Do[
		If[KeyFreeQ[$CGproperties, cg],
			Message[DefineCompositeCG::unkwnCG, cg];
			Abort[];
		];
	, {cg, cgs/. Bar-> Identity}];
	
	(*Test index contractions*)
	indexTypes= If[Head@ # === Bar, Bar/@ $CGproperties[Bar@ #, Indices], $CGproperties[#, Indices]]&/@ cgs;
	If[Length/@ indexTypes =!= Length/@ indicesIn, 
		Message[DefineCompositeCG::indCount]; 
		Abort[];
	];
	contractions= DeleteCases[Tally[Flatten@ indicesIn/. Bar@ x_-> x], {_, 1}];
	If[MemberQ[contractions, {_, x_}/; x> 2],
		Message[DefineCompositeCG::indMult];
		Abort[];
	]; 
	contractions= contractions[[;;, 1]];
	Do[ 
		inds= indexTypes[[Sequence@@ #]]&/@ Position[indicesIn, ind];
		If[inds =!= Bar/@ Reverse@ inds,
			Message[DefineCompositeCG::invldContraction, ind];
			Abort[];
		];
	, {ind, contractions}];
	
	(*Determine the index types of the resultant CG*)
	indices= MapThread[Index, {Flatten@ indicesIn, Flatten@ indexTypes}];
	eleminate= DeleteCases[Tally[indices/. HoldPattern@ Bar@ x_-> x], {_, 1}][[;;, 1]];
	eleminate= Alternatives@@ Join[eleminate, Bar/@ eleminate];
	indices= DeleteCases[indices, eleminate];
	indexTypes= indices/. Index[_, rep_]:> rep;
	
	tensor= EinsteinSummation[cgs/. $CGtensors, indicesIn-> Cases[Tally@ Flatten@ indicesIn, {lab_, 1}-> lab]];
	
	DefineCG[name, indexTypes, tensor];
	If[OptionValue@ ReplaceCGs,
		NewCGReplacement[CG[name, indices], 
			Times@@ Thread@ CG[cgs, indicesIn] ];
	];
]; 
DefineCompositeCG@ ___:= (Message[DefineCompositeCG::args]; Abort[];) 


(* ::Subsubsection::Closed:: *)
(*GetCGTensor*)


(* ::Text:: *)
(*Return the tensor representation of a CG *)


GetCGTensor::unkwn= "`1` is not a known CG. Verify that this CG has been loaded into the kernel.";


GetCGTensor[symb_]:= Block[{},
	If[!MemberQ[Keys@ $CGtensors, symb], 
		Message[GetCGTensor::unkwn, symb]; Abort[]; ];
	$CGtensors@ symb
];


(* ::Subsection:: *)
(*Clean up*)


(* ::Subsubsection::Closed:: *)
(*Remove CGs*)


RemoveCG::unkwn= "The CG `1` has not been definined.";
RemoveCG[name_]:= Block[{},
	If[!MemberQ[Keys@ $CGtensors, name],
		Message[RemoveCG::unkwn, name];
		Abort[];
	];
	
	KeyDropFrom[$CGtensors, name];
	KeyDropFrom[$CGtensors, Bar@ name];
	KeyDropFrom[$CGNumTensors, name];
	KeyDropFrom[$CGproperties, name];
	KeyDropFrom[$CGreplacements, name];
	ResetCGProjectors[];
	
	(*Remove Shorthand*)
	Switch[name
	, _Symbol,
		Clear@ name
	, _Symbol[_],
		SubValues[Evaluate@ Head@ name]= DeleteCases[SubValues[Evaluate@ Head@ name], 
			RuleDelayed[HoldPattern[Verbatim[HoldPattern][name[_] ] ], _] ];
	];
	Quiet@ UnsetBarable@ name;
];


(* ::Subsubsection::Closed:: *)
(*Remove representation*)


RemoveRepresentation::unkwn= "The Representation `1` has not been definined.";
RemoveRepresentation[name_]:= Block[{group, cg},
	If[!MemberQ[Keys@ $Representations, name],
		Message[RemoveRepresentation::unkwn, name];
		Abort[];
	];
	group= $Representations[name, Group];
	
	(*Removing representations*)
	Do[
		If[$CGproperties[cg, Group] =!= group, Continue[] ];
		If[MemberQ[$CGproperties[cg, Indices]/. {Bar@ x_:> x}, name],
			RemoveCG@ cg;
		];
	, {cg, Keys@ $CGproperties}];
	KeyDropFrom[$CGproperties, name];
	
	(*Unsetting the associated Bar properties*)
	If[$Representations[name, Reality]=== 1, 
		Bar@ ind:Index[_, name]=.;
		,
		UnsetBarable@ name;
	];
	KeyDropFrom[$Representations, name];
];


(* ::Subsubsection::Closed:: *)
(*Remove group*)


RemoveGroup::unkwn= "The Group `1` has not been definined.";
RemoveGroup[name_]:= Block[{rep},
	If[!MemberQ[Keys@ $Groups, name],
		Message[RemoveGroup::unkwn, name];
		Abort[];
	];
	
	(*If gauge or global symmetry group call the specific removal functions instead*)
	Which[KeyExistsQ[$GaugeGroups, name],
		RemoveGaugeGroup@ name;
		Return[];
	, KeyExistsQ[$GlobalGroups, name],
		RemoveGlobalGroup@ name;
		Return[];
	];
	
	KeyDropFrom[$Groups, name];
	Do[
		If[$Representations[rep, Group] === name, 
			RemoveRepresentation@ rep; 
		];
	, {rep, Keys@ $Representations}];
];


(* ::Subsubsection::Closed:: *)
(*Clear all CGs, groups, and representations*)


ClearGroups[]:= Block[{},
	If[Head@ $Groups === Association,
		RemoveGroup/@ Keys@ $Groups;
	,
		$Groups= <||>;
		$Representations= <||>;
		$CGtensors= <||>;
		$CGproperties= <||>;
		$CGreplacements= <||>;
	];
	ResetCGProjectors[];
	$CGcount= 1;
];
ClearGroups[]


(* ::Subsection:: *)
(*CG matching and contractions*)


(* ::Subsubsection::Closed:: *)
(*Construct CG basis *)


(* ::Text:: *)
(*Find the elementary tensors relevant for a set of indices *)


ElementaryTensors[indexTypes_List]:= Block[{cg, inds, out},
	out= Reap[Do[
		If[!$CGproperties[cg, InBasis], Continue[]]; (*Don't include redundant tensors*)
		If[SubsetMultQ[inds= $CGproperties[cg, Indices], indexTypes],
			Sow[cg-> inds];
		];
		(*Checks what conjugate tensors should be added *)
		If[$CGproperties[cg, UniqueConj] && SubsetMultQ[Bar/@ inds, indexTypes],
			Sow[Bar@ cg-> Bar/@ inds];
		];
	, {cg, Keys@ $CGproperties}]][[2]];
	
	If[Length@ out> 0, 
		Association@@ out[[1]],
		<||>
	] 
];


(* ::Text:: *)
(*Finds the full set of basis tensors for given open indices and their unique permutations.*)
(*This is not guaranteed to be a true basis as there might be some redundancy left over.*)


TensorBasis[indexTypes_List]:= Module[{compTensors, types, eTensInds, len, cur, target, type, 
	countMat, v, temp, proceed= True, tens, inds, t, fromOrigPerm, perms, sa},
	(*Finds elementary CGs and determine what index types are available from among them*)
	eTensInds= ElementaryTensors@ indexTypes; 
	types= DeleteDuplicates@ indexTypes;
	
	(*Sets a matrix counting what indices are available from a given set of elementary CGs*)
	countMat= Table[Count[#, type]&/@ (List@@ eTensInds), {type, types}];
	target= Table[Count[indexTypes, type], {type, types}];
	
	(*Loops through all possible combinations of ECGs consistent with indexTypes*)
	len= Length@ eTensInds;
	If[len < 1, Return@ {};]; 
	v= ConstantArray[0, len];
	cur= 1; 
	(*Counts the tensors as though they were numbers, increasing digits and carrying the remainder*)
	compTensors= Reap[
		While[proceed,
			v[[cur]]++; temp= countMat . v;
			While[proceed && VecCompare[Or, Greater, temp, target],
				If[cur< len,
					v[[cur++]]= 0; v[[cur]]++; temp= countMat . v;
				,
					proceed= False; 
				];
			];
			If[VecCompare[And, SameQ, temp, target], Sow@ v; ];
			cur= 1;
		];
	][[2]];
	
	(*Recast the vectors to groupings of ECGs*)
	If[Length@ compTensors < 1, Return@ {}; ];
	compTensors= Flatten@ MapThread[ConstantArray, {Keys@ eTensInds, #}]&/@ compTensors[[1]];
	
	(*Remove posibility of e.g., \[CurlyEpsilon]* \[CurlyEpsilon] combinations for SU(2) groups*)
	compTensors= RemoveDeltaDecomposables@ compTensors;
	
	(*For each ECG combination determines the linearly independent permutations*)
	compTensors= Table[
		inds= Flatten[tens/. eTensInds];
		(*DEPRECATED: all (non-unique) permutations*)
		(*perms= DegeneratePerms[inds, indexTypes];*)

		(*"Unique" permutation group of of the CG tensor product *)		
		perms= FactorOutSymmetrySubgroup[DegeneratePerms@ inds, SymmetriesOfTensorProduct[tens, eTensInds]];
		
		(*Map unto permutations of the original index set*)
		fromOrigPerm= FindPerm[indexTypes, inds];
		perms= InversePerm@ ComposePerm[#, fromOrigPerm]&/@ perms;
		
		sa= TensorProduct@@ tens/. $CGNumTensors;
		perms= Reap[
			While[Length@ perms>= 1,
				Sow@ perms[[1]];
				(*Check numerically for for tensor overlap*)
				t= PermuteArray[sa, perms[[1]] ];
				perms= DeleteCases[perms[[2;;]], _? (LinDepQ[t, PermuteArray[sa, #]]&)];
			];
		 ][[2, 1]];
		 Thread[{tens, perms}, List, -1] 
	, {tens, compTensors}];
	Flatten[compTensors, 1]	
];


(* ::Text:: *)
(*Determines full symmetry group inherent in tensor product of CGs *)


SymmetriesOfTensorProduct[tens_List, tensIndices_Association]:= Module[{i= 1, inds, p, perms, toOrigPerm},
	(*Enumerate the indices *)
	inds= Map[(i++ &), tens/. tensIndices, {2}];
	(*Find symmetries of the individual CGs*)
	perms= Keys@ $CGproperties[#, Symmetries]&/@ (tens/. Bar-> Identity);
	inds= MapThread[Table[#1[[p]], {p, #2}]&, {inds, perms}];
	(*Symmetries of like CG exchanges*)
	perms= Flatten/@ Tuples[Permutations/@ List@@ PositionIndex@ tens];
	toOrigPerm= FindPerm[First@ perms, Range@ Length@ First@ perms];
	perms= #[[toOrigPerm]]&/@ perms;
	Join@@ Table[Flatten/@ Tuples@ inds[[p]], {p, perms}]
]


(* ::Text:: *)
(*For comparing all numbers in two vectors *)


VecCompare[logOp_, comp_, v1_, v2_]:= logOp@@ MapThread[comp, {v1, v2}];


(* ::Text:: *)
(*Remove delta-decomposable  products of tensors from the basis sets*)


RemoveDeltaDecomposables@ tensorSets_:= Module[{decomposableCGs, out= tensorSets},
	decomposableCGs= Intersection[DeleteDuplicates@ Flatten[tensorSets/. Bar-> Identity],
		GetByProperty[$CGproperties, DeltaDecomposable-> True] ];
	decomposableCGs= {#, Bar@ #}&/@ decomposableCGs;
	Do[
		out= Select[out, (!ContainsAll[#, dcg] &)];
	, {dcg, decomposableCGs}];
	out
]


(* ::Subsubsection::Closed:: *)
(*Projection matrices for CGs*)


(* ::Text:: *)
(*Determines the pieces necessary  for projecting tensor to CG basis in the vector space of CGs*)


CGProjector[indexTypes_List]:= CGProjector[indexTypes]= Block[
		{tensors, m, n, invmetric, sa},
	tensors= TensorBasis@ indexTypes;
	If[Length@ tensors < 1, Return@ {{}, {{}}}; ];
	
	(*Find overlap matrix. Utilizes hermitianity to save half the computation time*)
	invmetric= Table[
			sa= TensProdPermutationToSA[Bar/@ tensors[[m, 1]], tensors[[m, 2]]];
			Table[
				If[n >= m, 
					SAProduct[sa, TensProdPermutationToSA[tensors[[n, 1]], tensors[[n, 2]]] ]
				,
					0
				]
			, {n, Length@ tensors}]
		, {m, Length@ tensors}];
	invmetric+= ConjugateTranspose@ ReplacePart[invmetric, {i_, i_}-> 0];
	
	{tensors, invmetric}
];


(* ::Text:: *)
(*Call to remove memoized down-values from CGProjector*)


ResetCGProjectors[]:= RemoveAssociatedDownValues@ CGProjector[{__}];


(* ::Text:: *)
(*Function to convert a symbolic list of CGs with a permutation into a SparseArray, *)


TensProdPermutationToSA[cgList_List,perm_List]:=
	PermuteArray[TensorProduct@@ Replace[cgList, $CGtensors, {1}], perm];


(* ::Subsubsection::Closed:: *)
(*Match to CG basis *)


(* ::Text:: *)
(*Function for matching a tensor with certain indices to the predefined set of CG coefficients. *)


MatchToCGs::failed= "Failed to add new CG basis tensor from `1`"


(*Default behavior for numbers\[LongDash]rank-zero tensors*)
MatchToCGs[num_, {}, originalTensors_: None]:= num;


MatchToCGs[tensor_SparseArray, inds_, originalTensors_: None]:= Block[{out, indTypes},
	out= FindCGMatch[tensor, inds];
	
	(*If no match found, add define new CG*)
	If[out === None,
		indTypes= inds/. Index[_, rep_]-> rep;
		DefineNewCG[indTypes, tensor];
		out= FindCGMatch[tensor, inds];
		If[out === None,
			Message[MatchToCGs::failed, originalTensors];
			Abort[];
		];
		NewCGReplacement[out, originalTensors];
	];
	out
];


(* ::Text:: *)
(*Attempts to decompose tensor on the basis of known tensors*)


FindCGMatch[tensor_, inds_]:= Block[{indTypes, indTypesOut, indPerm, tens, indsOut, indices,
		cgBasis, cgMetric, overlaps, basis, compTensor, composition, n, temp, order},
	(*If input tensor vanishes*)
	If[SAZeroQ@ tensor, Return@ 0;]; 
	
	(*Find index types and associated permutation to bring it to canonical form*)
	indTypes= inds/. Index[_, rep_]-> rep;
	indPerm= FindPerm[indTypes];
	indsOut= inds[[indPerm]];
	indTypesOut= indTypes[[indPerm]];
	(*Conjugate, permuted tensor*)
	tens= PermuteArray[tensor, indPerm];
	
	(*Find basis tensors*)
	{cgBasis, cgMetric}= CGProjector[indTypesOut]; 
	If[Length@ cgBasis< 1,
		Return@ None;
	]; 
	
	(*Project tensor: The tensors with least overlap in each redundant set are removed 
		to reduce the relevant CG tensors to a basis for the projection.*)
	overlaps= SAProduct[TensProdPermutationToSA@@ MapAt[Bar, #, {1, All}], tens]&/@ cgBasis;
	order= Ordering[- Abs@ overlaps/ Sqrt@ Diagonal@ cgMetric// N];	
	basis= DeleteCases[RowReduce@ cgMetric[[order, order]], {0..}];
	basis= FirstPosition[#, 1]&/@ basis// Flatten;
	basis= order[[basis]];
	
	(*Project tensor*)
	overlaps= overlaps[[basis]]; 
	cgBasis= cgBasis[[basis]];
	(*Solve g_ij x_j = y_i*)
	cgMetric= cgMetric[[basis, basis]];
	composition= RowReduce[Join[cgMetric, List/@ overlaps, 2]][[;;, -1]];
		
	(*Test if projection succesfull: otherwise add new CG*)
	compTensor= Sum[(*Evaluate only the SAs with non-trivial overlap with the CG*)
			If[composition[[n]] =!= 0,
				composition[[n]]* TensProdPermutationToSA@@ cgBasis[[n]]
			,
				0
			]
		, {n, Length@ composition}];
	If[!SAZeroQ[tens- compTensor],
		Return@ None;
	];
	
	(*Produce CG decomposition*)
	Sum[
		indices= TakeList[indsOut[[InversePerm@ cgBasis[[n, 2]] ]],
			(Length@ $CGproperties[#/. Bar-> Identity, Indices]&)/@ cgBasis[[n, 1]]];
		composition[[n]] Times@@ MapThread[CG, {cgBasis[[n ,1]], indices}]
	,{n, Length@ composition}]
];


(* ::Text:: *)
(*Add new CG for the unknown tensor*)


$CGcount= 1;
DefineNewCG[indexTypes_List, tensor_SparseArray]:= Module[{lab},
	lab= Symbol["cg"<> ToString@ $CGcount++];
	DefineCG[lab, indexTypes, tensor];
	Echo@ StringForm["Added new CG `1` with indices `2`", lab, indexTypes];
];


(* ::Subsubsection::Closed:: *)
(*Contracting CGs in an expression*)


(* ::Text:: *)
(*For bringing a single CG onto canonical form*)


IdentifyCG@ CG[tens_, inds_]/; !DuplicateFreeQ[inds/. Bar@ x_-> x]:= 
	MatchToCGs@@ EinsteinSummation[{tens/. $CGtensors}, {inds/. Bar@ x_-> x}, KeepIndices-> True];
IdentifyCG@ CG[tens_SparseArray, inds_]/; DuplicateFreeQ[inds/. Bar@ x_-> x]:= MatchToCGs[tens, inds];
IdentifyCG@ CG[tens_, inds_]/; DuplicateFreeQ[inds/. Bar@ x_-> x]:= CG[tens, inds];
IdentifyCG@ CG[Bar@ tens_, inds_]/; DuplicateFreeQ[inds/. Bar@ x_-> x]:= 
	If[$CGproperties[tens, UniqueConj], 
		CG[Bar@ tens, inds], 
		MatchToCGs[Bar@ tens/. $CGtensors, inds]
	]; 


(* ::Text:: *)
(*Contracts CGs in an expression according to their indices. The resulting tensors are matched to the CG basis. *)


ContractCGs::repInds= "Index label \"`1`\" occurs `2` times in the expression.";


ContractCGs@ expr:(_List|_Plus|_HcTerms):= ContractCGs/@ expr; 


ContractCGs@ expr_:= Block[{terms= LagrangianExpand@ expr},
	If[Head[terms] === Plus,
		ContractCGsSingleTerm/@ terms
	,
		ContractCGsSingleTerm@ terms
	]
];


(* ::Text:: *)
(*Applies contractions to single terms. Both Times and Operator threads over sums that may be generated from contracted CGs. *)


ContractCGsSingleTerm@ cg_CG:= IdentifyCG@ cg;


ContractCGsSingleTerm[term_Times]:= Block[{out= term},
	If[FreeQ[out, _CG], Return@ out; ];
	If[FreeQ[out, _Operator], 
		Times@@ ContractCGsProduct@@ out
	,
		out/. op_Operator:> Operator@@ ContractCGsProduct@@ op
	]// Expand 
];


ContractCGsSingleTerm[term_Power]:= Block[{out= term},
	If[FreeQ[out, _CG], Return@ out; ];
	Times@@ ContractCGsProduct[term]// Expand
];


ContractCGsSingleTerm[term_Operator]:= Block[{out= term},
	If[FreeQ[out, _CG], Return@ out; ];
	Operator@@ ContractCGsProduct@@ out// Expand
];


ContractCGsSingleTerm@ other_:= other;


(* ::Text:: *)
(*Contracts all CGs among factors of a product (or equivalent)*)


ContractCGsProduct[factors__]:= Block[
		{cg, cgs, head, indRules, repeatedInds, contractedCGs, set, uncontractedCGs, out= {factors}},
		
	(*Extract all CGs*)
	out= out/. Power[c_CG, 2]:> CGProduct[c[[1]], c[[1]]];
	cgs= Cases[out, _CG];
	out= DeleteCases[out, _CG];
	
	(*Identify conjugate repreated indices*)
	repeatedInds= Tally@ Cases[cgs, _Index, {3, 4}];
	If[MemberQ[repeatedInds, {_, n_/; n> 2}],
		Message[ContractCGs::repInds, FirstCase[repeatedInds, {ind_, n_}:> Sequence[ind, n]]];
		Abort[];
	]; 
	repeatedInds= DeleteCases[repeatedInds, {_, 1}][[;;, 1]];
	If[Length@ repeatedInds> 0,
		indRules= MapThread[Rule, {Bar/@ repeatedInds, repeatedInds}];
		cgs= cgs/. indRules;
	];
	
	(*Determine the disjoint CG components*)
	contractedCGs= Position[cgs, #][[;;, 1]]&/@ repeatedInds;
	contractedCGs= DeleteDuplicates/@ (contractedCGs//. 
		{OrderlessPatternSequence[{a___, x_, b___}, {c___, x_, d___}, rest___]}-> {{x, a, b, c, d}, rest});
	
	(*All the uncontracted CGs are identified individually*)
	uncontractedCGs= IdentifyCG/@ cgs[[Complement[Range@ Length@ cgs, Flatten@ contractedCGs] ]];
	
	(*Contract the CGs within the disjoint sets and match to basis*)
	contractedCGs= Table[
			set= cgs[[set]];
			cg= EinsteinSummation[set[[;;, 1]]/. $CGtensors, set[[;;, 2]], KeepIndices-> True];
			If[Length@ cg[[2]]=== 0,
				cg[[1]]
			,
				MatchToCGs[cg[[1]], cg[[2]], Times@@ set]
			]
		, {set, contractedCGs}];
	
	(*Returns list of factors*)
	Join[out, uncontractedCGs, contractedCGs]
];


(* ::Subsection:: *)
(*CG replacements*)


(* ::Text:: *)
(*To automatically replace composite CGs, especially the automatically constructed ones, with the  factorized form. *)


(* ::Subsubsection::Closed:: *)
(*Make new replacement rule*)


(* ::Text:: *)
(*Function for adding a new CG replacement rule. `contraction' is a product of CGs.*)


NewCGReplacement[_, None]:= Null; 
NewCGReplacement[cg_CG, contraction_]:= Module[{inds},
	$CGreplacements@ First@ cg= If[$CGproperties[First@ cg/. Bar-> Identity, UniqueConj],
			{MakeCGReplacement[cg, contraction], MakeCGReplacement[Bar@ cg, Bar@ contraction]}
		,
			{MakeCGReplacement[cg, contraction]}
		];
]


(* ::Text:: *)
(*Construct replacement rule*)


MakeCGReplacement[cg_, contraction_]:= Module[{inds},
	inds= Tally@ Cases[contraction, _Index, Infinity];
	With[{
			org= contraction/. c_CG:> ReconstructBars@ c,
			indSing= Cases[inds, {i_, 1}:> Alternatives[Bar@ i, i]-> Unique@"v$"],
			indDup= Cases[inds, {Index[i_, rep_], 2}:> Index[i, rep]-> Hold@ Index[Unique@"u$", rep]]
		},
		(cg/. MapAt[Pattern[#, Blank[]]&, indSing, {All,2}]):> (org/. ReleaseHold@ indDup/. indSing)
	]
]


(* ::Text:: *)
(*Reconstruct the Bars that have previously been removed on contracted indices*)


ReconstructBars@ CG[cg_, inds_]:= Module[{conj, symb, bars, newInds},
	symb= If[(conj= Head@ cg === Bar),
			Bar@ cg,
			cg
		];
	bars= If[conj, Bar, Identity]/@ $CGproperties[symb, Indices];
	newInds= MapThread[(If[Head[#1/. Index[_, rep_]-> rep] =!= Head@ #2, Bar@ #1, #1] &),{inds, bars}];
	CG[cg, newInds]
] 


(* ::Subsubsection::Closed:: *)
(*Substitution of the CGs*)


(* ::Text:: *)
(*Functions that replace the CGs*)


Options@ ReplaceCGs= {CGsToReplace-> All}


ReplaceCGs::uknwnReplacement= "There is no replacement rule assocated with one/some of the keys `1`."


ReplaceCGs[expr_, OptionsPattern[]]? OptionsCheck:= Module[{keys, rules},
	rules= Flatten@ Switch[keys= OptionValue@ CGsToReplace
		,All,
			List@@ $CGreplacements
		,None, 
			{}
		,_List,
			If[!SubsetQ[Keys@ $CGreplacements, keys],
				Message[ReplaceCGs::uknwnReplacement, keys];
				Abort[];
			]; 
			Lookup[$CGreplacements, keys]
	];
	expr/. rules
]


(*ReplaceCGs indices on HcTerms*)
ReplaceCGs[HcTerms[expr_], opt:OptionsPattern[]]:=HcTerms[ReplaceCGs[expr, opt]];


(* ::Section:: *)
(*Auxiliary function*)


(* ::Subsection:: *)
(*Einstein summation*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


(* ::Text:: *)
(*Routine for contracting a list of SparseArray tensors, with indices {{a, b}, {c,b, ...},  ...} to final tensor with indices {c, ...}. *)


Options@ EinsteinSummation= {KeepIndices-> False};


EinsteinSummation[tensorsIn_, in_, opt: OptionsPattern[]]:= EinsteinSummation[tensorsIn, in-> {}, opt];


EinsteinSummation[tensorsIn_List, in_-> out_, OptionsPattern[]]:= Block[{tinds= in, tensors= tensorsIn, 
		contract, ind, n, traced},
	(*Trace any indvidual tensors if necessary*)
	Do[
		{contract, traced}= IndexContractions@ tinds[[n]];
		tensors[[n]]= SARemove0@ TensorContract[tensors[[n]], contract];
		tinds[[n]]= DeleteCases[tinds[[n]], Alternatives@@ traced];
	, {n, Length@ in}];
	
	(*Contract pairwise*)
	{tensors, tinds}= EinsteinSummationByPairs[tensors, tinds];
	
	(*Rearrange indices*)
	If[Length@ out> 0,
		tensors= Transpose[tensors, FindPermutation[tinds, out]];
		tinds= out;
	];
	
	If[OptionValue@ KeepIndices,
		{tensors, tinds}
	,
		tensors
	]	
];


(* ::Subsubsection::Closed:: *)
(*Assets*)


(* ::Text:: *)
(*The strategy for contracting multiple tensors is to contract them pairwise, starting from the ones multiplying to give a tensor of lowest possible rank. *)


EinsteinSummationByPairs[tensors_, tinds_]:= Block[{contract, pairs, tens, inds},
	If[Length@ tinds < 2,
		Return[{First@ tensors, First@ tinds}];
	]; 
	
	(*The next contraction is that which makes for the fewest open indices*)
	pairs= Subsets[MapIndexed[{First@ #2, #1}&, tinds], {2}];
	contract= SortBy[pairs, 
		({Length@ Complement[Union@@ #[[;;, 2]], Intersection@@ #[[;;, 2]] ],
			Total[Length/@ #[[;;, 2]]]}&)][[1, ;;, 1]];
	
	(*Contract relevant pair*)
	{tens, inds}= EinSum2[tensors[[contract]], tinds[[contract]]];
	
	(*Reapply on the remaining tensors*)
	EinsteinSummationByPairs[
		Append[Delete[tensors, List/@ contract], tens],
		Append[Delete[tinds, List/@ contract], inds] ]
];


(* ::Text:: *)
(*Compute the contraction of two tensors *)


EinSum2[tensors_, in_]:= Block[{indices, tensor, contract, traced},
	indices= Flatten@ in;
	{contract, traced}= IndexContractions@ indices;
	(*The Activate@ TensorContract[Inactive[TensorProduct][...], ...] call is not documented but mentioned here: https://mathematica.stackexchange.com/questions/77511/efficient-tensor-product-followed-by-contraction 
	It seems to behave more as a dot product than a TensorProduct followed by traces, as would be the case w/o the inactivation and which is extremely inefficient.*)
	tensor= Inactive[TensorProduct]@@ tensors;
	{SARemove0@ Activate@ TensorContract[tensor, contract], DeleteCases[indices, Alternatives@@ traced]}
];


(* ::Text:: *)
(*Given a list of indices, determine which are to be contracted*)


IndexContractions@ inds_List:= Block[{doubles},
	doubles= Flatten@ DeleteCases[Tally@ inds, {_, 1}][[;;, 1]];
	{Flatten@ Position[inds, #]&/@ doubles, doubles}
];


(* ::Subsection:: *)
(*Sparse arrays*)


(* ::Subsubsection::Closed:: *)
(*Sparse array properties*)


(* ::Text:: *)
(*To remove zeros from sparse arrays*)


SARemove0@ expr_SparseArray:= SparseArray@ expr;
SARemove0@ expr_:= expr;


(* ::Text:: *)
(*Checks if a SparseArray or Tensor vanishes identically*)


SAZeroQ@ expr_:= MatchQ[{0}]@ DeleteDuplicates@ Simplify@ DeleteDuplicates@ ArrayRules[expr][[;;,2]];


(* ::Subsection:: *)
(*Tensor functions*)


(* ::Subsubsection::Closed:: *)
(*Tensor symmetries *)


(* ::Text:: *)
(*Returns an association of symmetries and corresponding signatures*)


FindTensorSymmetries[tens_SparseArray, indTypes_]:= Block[{couplingPermutations, cycles, newPerms},
	(*Determine the generators (with signature) for the tensor symmetries*)
	cycles= TensorSymmetry@ tens/.{
			Symmetric[{a_, b_, c___}]:> {{Cycles[{{a, b}}], 1}, {Cycles[{{a, b, c}}], 1}},
			Antisymmetric[{a_, b_, c___}]:> {{Cycles[{{a, b}}], -1}, 
				{Cycles[{{a, b, c}}], -Power[-1, Length@ {a, b, c}]}}
		}// DeleteDuplicates;
		
	(*Delete symmetries not between identical index types*)
	cycles= DeleteCases[cycles, {Cycles[sets_], _}/; (Nand@@ (MatchQ[indTypes[[#]], {(a_)..}]&/@ sets)) ];
	
	(*Generate all permutations from the symmetry generators*)
	couplingPermutations=<|Range@ Length@ indTypes-> 1|>;	
	newPerms= couplingPermutations;
	While[Length@ newPerms > 0,
		newPerms= Association@@ (KeyValueMap[Function[{key, val}, Permute[key, #[[1]]]-> val #[[2]]], 
			newPerms]&/@ cycles); 
		KeyDropFrom[newPerms, Keys@ couplingPermutations];
		couplingPermutations= Join[couplingPermutations, newPerms];
	];
	couplingPermutations
];


(* ::Text:: *)
(*Gives the symmetries of a tensor as {List of symmetries, List of anti-symmetries}. *)


(*FindTensorSymmetries[tens_SparseArray]:= Block[{rank= TensorRank@ tens, i, j, taken, temp, 
	out, dims= Dimensions@ tens, mtens= -tens},
	(*Determines all indices with symmetries between them*)
	taken= ConstantArray[False, rank];
	out= Reap[Do[
		If[taken[[i]], Continue[];];
		Sow[i, i];
		Do[
			If[dims[[i]] =!= dims[[j]] || taken[[j]], Continue[];];
			temp= SparseArray@ Transpose[tens, i<-> j];
			If[temp === tens || temp === mtens, 
				Sow[j, i];
				taken[[j]]= True;
			]; 
		,{j, i+ 1, rank}]; 
	,{i, rank}]][[2]];
	
	(*Sorts the symmetries according to sign*)
	out= {tens === SparseArray@ Transpose[tens, #[[1]]<-> #[[2]]], #} &/@ DeleteCases[out, _? (Length@ # === 1&)];
	DeleteCases[GatherBy[Join[{{True, {}}, {False, {}}}, out], First][[;;, ;;, 2]], {}, {2}]
];*)
