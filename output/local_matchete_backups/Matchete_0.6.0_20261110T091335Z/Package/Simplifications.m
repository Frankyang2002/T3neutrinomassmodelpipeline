(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Simplifications`*)


(* ::Subtitle:: *)
(*Paclet for simplifying and reducing EFT operators to an off-shell basis*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["GreensSimplify"]
PackageExport["IBPIdentities"]
PackageExport["PreferInputOperators"]


PackageExport["CollectOperators"]
PackageExport["SelectOperatorClass"]
PackageExport["GetOperatorCoefficient"]


PackageExport["HcSimplify"]
PackageExport["HcExpand"]
PackageExport["HcTerms"]


PackageExport["ReductionIdentities"]
PackageExport["dDimensional"]
PackageExport["FourDimensional"]
PackageExport["EvanescenceFree"]
PackageExport["Evanescent"]
PackageExport["VLFBasis"]
PackageExport["VectorLike"]


PackageExport["Conj"]


(* ::Text:: *)
(*Evanescent*)


PackageExport["EvaOp"]
PackageExport["GetEvanescentTerm"]
PackageExport["ResetEvanescentOperators"]
PackageExport["ExpandEvanescentOperators"]
PackageExport["Origin"]
PackageExport["RedundantOperator"]
PackageExport["BasisOperator"]
PackageExport["RemoveLCTensor"]
PackageExport["Schouten"]


(* ::Subsubsection::Closed:: *)
(*Scoped*)


PackageScope["Operator"]
PackageScope["NormalToOperatorForm"]
PackageScope["OperatorToNormalForm"]
PackageScope["ToAtomicForm"]
PackageScope["AtomicToOperatorForm"]
PackageScope["AtomicToNormalForm"]
PackageScope["ReleaseOperators"]


PackageScope["NormalForm"]


PackageScope["EoM"]
PackageScope["SeparateOutConstants"]


PackageScope["CanonicalKineticOpQ"]
PackageScope["KineticOpQ"]


PackageScope["InternalSimplify"]
PackageScope["InternalOpRepresentation"]
PackageScope["InequivalentPermutations"]
PackageScope["ConjugateIndexExchange"]
PackageScope["Score"]
PackageScope["PromoteOperatorsToBasis"]
PackageScope["UpdateChiralBasis"]


PackageScope["AtomicOp"]
PackageScope["CompOp"]
PackageScope["ConstructOperatorIdentities"]
PackageScope["ResetOperatorAssociations"]
PackageScope["MatchOperatorPatterns"]
PackageScope["LookupOperatorFlavorProperties"]
PackageScope["$operators"]
PackageScope["$compoundOperators"]
PackageScope["HcSimplifyInternal"]

PackageScope["IdentitiesSymmetry"]


PackageScope["EOMDevs"]


PackageScope["CanonizeKinetic"]


PackageScope["SeparateInteractionTerm"]


PackageScope["InternalCollectOperators"]


(* ::Text:: *)
(*Evanescent operators*)


PackageScope["DefineEvanescentOperator"]
PackageScope["FindEvanescentOperator"]
PackageScope["$EvanescentTerms"]
PackageScope["Reducible4dOpQ"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


CollectOperators::usage= "CollectOperators[expr] collects all identical operators (under index relabeling and obvious symmetries) in an expression."
SelectOperatorClass::usage= "SelectOperatorClass[expr, fields, devNo] selects all terms in the expression with operators consisting of the specified matter fields (including Bars) and number of derivatives. Field strength tensors count as 2 derivatives."
GetOperatorCoefficient::usage= "GetOperatorCoefficient[expr, operator] extracts the coefficient of the operator matching exactly (except for index names) the form given by operator."


GreensSimplify::usage= "GreensSimplify[expr] simplifies an expression with IBP identities and commutation relations."
GreensSimplify::partialdev= "The expression contains partial derivative indices (PartialDev). GreensSimplify with ReductionIdentities -> EvanescenceFree requires a gauge-invariant expression."
IBPIdentities::usage= "IBPIdentities[fields, devNo] returns a list with all the integration-by-part identities used by GreensSimplify for the class of operators consisting of the exact matter fields (including Bars) and the specified number of derivatives. Field strength tensors count as 2 derivatives."
PreferInputOperators::usage= "Option for simplification functions to decide if the input operators should be prioritized for the output basis."


HcSimplify::usage="HcSimplify[Lagrangian] identifies terms with their hermitian conjugate and retains only one of them under the head HcTerms.";
HcExpand::usage="HcExpand[Lagrangian] expands the head HcTerms to return both the operators and their hermitian conjugate.";
HcTerms::usage="HcTerms[expr] is a placeholder that contain terms whose hermitian conjugate have been removed form the Lagrangian with HcSimplify.";


ReductionIdentities::usage="Option for GreensSimplify to decide whether or not to use the subset of identities valid only in 4 dimensions (fierz, gamma reduction and Levi-Civita relations) in the simplification. The value can be dDimensional (4D identies are not included), FourDimensional (4D identities are included), Evanescent (4D identities are included with the corresponding evanescent operators) or EvanescenceFree (4D identities are included and evanescent contributions are evaluated)."


EvaOp::usage                     = "EvaOp[class, id, (indices)] denotes the evanescent operator stored in the association $EvanescentTerms. To obtain the definition, use GetEvanescentTerm[class, id]."
Evanescent::usage                = "Evanescent is an option for Fierz, GammaReduction and RemoveLCTensor. If set to True (default), it stores the evanescent operator (if generated) in an association and includes it in the result."
GetEvanescentTerm::usage         = "GetEvanescentTerm[class, (id)] outputs the informations (Operators, Origin) about the (class of) evanescent operator(s) given as an argument. If no argument is given, it returns the full list of generated evanescent operators."
ResetEvanescentOperators::usage  = "ResetEvanescentOperators[] resets the association containing the definition of all evanescent operators generated."
ExpandEvanescentOperators::usage = "ExpandEvanescentOperators[expr] replaces EvaOp[class, id, (indices)] by the corresponding definition of the evanescent operator."
FindEvanescentOperator::usage    = "Identifies an evanescent operator already defined and replace it with an EvaOp object."
VLFBasis::usage                  = "Option name for simplificaiton functions to determine the gamma matrix basis for vector-like fermions."
VectorLike::usage                = "Option used in Matchete."


Conj::usage= "Conj is used as a head that abstractly denotes the conjugate of, e.g., field labels for conjugate degrees of freedom."


RedundantOperator::usage= "RedundantOperator is an attribute of evanescent operators. It stores the operator that generated the evanescent structure before the 4D operation."
BasisOperator::usage= "BasisOperator is an attribute of evanescent operators. It stores the operator generated by RedundantOperator after applying the 4D operation."
Origin::usage= "Origin is an attribute of evanescent operators. It is either Fierz, GammaReduction, LCTensor, or Schouten."
Schouten::usage= "Denotes that an evanescent operator originates from applying Schouten identities to Lorentz indices.";


RemoveLCTensor::usage            = "Get rid of Levi-Civita tensors either by combining with a Lorentz sigma to give a \!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) or by expanding the product of two Levi-Civita tensor to give the fully antisymmetric product of metrics."


(* ::Subsubsection::Closed:: *)
(*Scoped*)


EoM::usage="EoM[A_] is a field operator in the form of the equations of motion of a field.";
Operator::usage="Operator[x] separates prefactors from field operators and puts the expressions into forms that can more easily be processed by the reduction algorithms.";


LookupOperatorFlavorProperties::usage= "LookupOperatorFlavorProperties[op] returns the flavor properties of the operator in question."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Operator part*)


(* ::Subsection:: *)
(*Operator object*)


(* ::Subsubsection::Closed:: *)
(*Properties*)


(* ::Text:: *)
(*Operator is separates prefactors from the fields and group structures of an object. It distributes over an expression, so after applying Operator to an expression, all operators are contained in an Operator object.*)


SetAttributes[Operator, {Orderless}];


(*Flat*)
HoldPattern@ Operator[Operator@ x___, y___]:= Operator[x, y];
(*Distributive *)
Operator@ l_List:= Operator/@ l;
op:Operator[_Plus, ___]:= Distribute@ Unevaluated@ op;

(*Operator@ prod_Times:= Operator@@ prod;*)
Operator[HoldPattern@ Times@ a__, b___]:= Operator[a, b];

Operator[a_/; FreeQ[a, Field| FieldStrength| CG| LCTensor| ExtMom], b___]:= a Operator@ b;

Operator[Power[a_, n_Integer? Positive], b___]:= Operator[Sequence@@ ConstantArray[a, n], b];

Operator[]:= 1;


(* ::Subsubsection::Closed:: *)
(*Recognize EoM*)


(* ::Text:: *)
(*Setting fields to 0 sets EoM to zero*)


EoM@ 0= 0;


(* ::Text:: *)
(*Behavior of EoM under conjugation*)


EoM/: Bar@ EoM@ X_:=EoM@ Bar@X


(* ::Text:: *)
(*Structures occurring with Majorana fermions*)


EoM[(f:Transp@ Field[_, Fermion, _, {}])\[CenterDot] CC]:= -EoM@ f \[CenterDot] CC;
EoM[-CC\[CenterDot] (f:Field[_, Fermion, _, {}])]:= CC\[CenterDot] EoM@ f ;


(* ::Text:: *)
(*Recognize EoM-type structures...*)
(*	- Matching on, e.g., a_Index ensures that we exclusively match on covariant derivatives.*)


(* Scalars *)
Operator[eom:(Field[_, Scalar, _, {mu_Index, mu_Index}]| Bar@ Field[_, Scalar, _, {mu_Index, mu_Index}]), rest___]:=
	Operator[eom/. Field[l_, Scalar, i_, {b_, b_}]:> EoM@ Field[l, Scalar, i, {}], rest];

(* Fermions *)
Operator[X___\[CenterDot] DiracProduct[A___, GammaM[mu_]]\[CenterDot] Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	Operator[X\[CenterDot] DiracProduct[A]\[CenterDot] EoM[Field[f, Fermion, i, {}]], rest];
Operator[X___\[CenterDot] DiracProduct[A___, GammaM[mu_], Proj@ p_]\[CenterDot] Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	Operator[X\[CenterDot] DiracProduct[A, Proj[-p]]\[CenterDot] EoM[Field[f, Fermion, i, {}]], rest];
Operator[X___\[CenterDot] DiracProduct[A___, GammaM[mu_], Gamma5]\[CenterDot] Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	-Operator[X\[CenterDot] DiracProduct[A, Gamma5]\[CenterDot] EoM[Field[f, Fermion, i, {}]], rest];
(* Bar@ Fermions *)
Operator[Bar@ Field[f_, Fermion, i_, {mu_Index}]\[CenterDot] DiracProduct[GammaM[mu_], A___]\[CenterDot] X___, rest___]:=
	Operator[EoM[Bar@ Field[f, Fermion, i, {}]]\[CenterDot] DiracProduct[A]\[CenterDot] X, rest];
Operator[Bar@ Field[f_, Fermion, i_, {mu_Index}]\[CenterDot] DiracProduct[GammaCC,Transp@GammaM[mu_], A___]\[CenterDot] X___, rest___]:=
	-Operator[EoM[Bar@ Field[f, Fermion, i, {}]]\[CenterDot] DiracProduct[GammaCC,A]\[CenterDot] X, rest];
(* Transp@ Fermions *)
Operator[Transp@ Field[f_, Fermion, i_, {mu_Index}]\[CenterDot] DiracProduct[GammaCC, GammaM[mu_], A___]\[CenterDot] X___, rest___]:=
	-Operator[EoM[Transp@ Field[f, Fermion, i, {}]]\[CenterDot] DiracProduct[GammaCC, A]\[CenterDot] X, rest];
(*Transp@ Bar@ Fermions*)
Operator[X___\[CenterDot] DiracProduct[A___, Transp@ GammaM[mu_]]\[CenterDot] Transp@ Bar@ Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	Operator[X\[CenterDot] DiracProduct[A]\[CenterDot] EoM[Transp@ Bar@ Field[f, Fermion, i, {}]], rest];
Operator[X___\[CenterDot] DiracProduct[A___, Transp@ GammaM[mu_], Proj@ p_]\[CenterDot] Transp@ Bar@ Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	Operator[X\[CenterDot] DiracProduct[A, Proj[-p]]\[CenterDot] EoM[Transp@ Bar@ Field[f, Fermion, i, {}]], rest];
Operator[X___\[CenterDot] DiracProduct[A___, Transp@ GammaM[mu_], Gamma5]\[CenterDot] Transp@ Bar@ Field[f_, Fermion, i_, {mu_Index}], rest___]:=
	-Operator[X\[CenterDot] DiracProduct[A, Gamma5]\[CenterDot] EoM[Transp@ Bar@ Field[f, Fermion, i, {}]], rest];

(*Vectors*)
Operator[FieldStrength[V_, linds:{OrderlessPatternSequence[a_, b_]}, ind_, {a_Index}], rest___]:=
	Signature@ linds Signature@ {a, b} Operator[EoM@ Field[V, Vector[b], ind, {}], rest];
Operator[Bar@ FieldStrength[V_, linds:{OrderlessPatternSequence[a_, b_]}, ind_, {a_Index}], rest___]:=
	Signature@ linds Signature@ {a, b} Operator[Bar@ EoM@ Field[V, Vector[b], ind, {}], rest];


(* ::Subsubsection::Closed:: *)
(*NormalForm, rewriting Operator and EoM objects*)


Options@ OperatorToNormalForm= {Unique-> False, CanonizeKinetic -> False, RelabelIndices->True};


OperatorToNormalForm[expr_,OptionsPattern[]]:= 
	If[OptionValue@ CanonizeKinetic, CanonizeKineticalOps@ expr, expr]/.
		op_Operator:> Activate@ If[OptionValue@ RelabelIndices, RelabelNonFlavorIndices[#, Unique-> OptionValue@Unique]&, Identity][
			Inactive[Times]@@ op/. EoM@ x_:> EoMStandardForm@ x]


(* ::Text:: *)
(*Function for bringing Operator objects back into normal form while expanding any EoM objects. *)


(*Options@ NormalForm= {Unique-> True, CanonizeKinetic -> True};*)


(* ::Text:: *)
(*The call of OperatorFlavorSeparate is important to ensure a proper treatment of all diagonal flavor indices in the call of RelabelIndices. The final contract removes the deltas.*)


(*NormalForm[expr_, OptionsPattern[]]:= If[OptionValue@CanonizeKinetic,CanonizeKineticalOps[expr],expr] /.
		op_Operator:> Activate@ RelabelNonFlavorIndices[(Inactive[Times]@@ op)/. EoM@ x_:> EoMStandardForm@ x, Unique-> OptionValue@ Unique]*)


(*
NormalForm[expr_, OptionsPattern[]]:= (*ContractDelta[*)
	CanonizeKineticalOps[expr, OptionValue@CanonizeKinetic]/.
		op_Operator:> OperatorFlavorSeparate[op]/.
		op_Operator:> Activate@ RelabelIndices[(Inactive[Times]@@ op)/. EoM@ x_:> EoMStandardForm@ x,
			Unique-> OptionValue@ Unique]/.
		x: Times[_Operator, __]:> ContractDelta[x, Expand-> False]
(*,Expand->False]*) (* The overall ContractDelta can lead to incorrect reults when using the option Expand->False *)
*)


(* ::Text:: *)
(*Lightweight function for releasing operator form in a Lagrangian*)


ReleaseOperators[expr_]:= expr/.
	op_Operator:> Activate[Inactive[Times]@@ op/. EoM@ x_:> EoMStandardForm@ x]


(* ::Text:: *)
(*Rewrite EOMs as standard field expressions*)


EoMStandardForm[f:Field[_, Scalar, __]]:= Module[{a}, CD[{a, a}, f] ];
EoMStandardForm[Bar@ f:Field[_, Scalar, __]]:= Bar@ Module[{a}, CD[{a, a}, f] ];
EoMStandardForm[f: Field[_, Fermion, __]]:= Module[{a}, \[Gamma][a]\[CenterDot] CD[a, f]];
EoMStandardForm[Bar@ f:Field[_, Fermion, __]]:= Bar@ Module[{a}, \[Gamma][a]\[CenterDot] CD[a, f]];
EoMStandardForm[Transp@ f:Field[_, Fermion, __]]:= Transp@ Module[{a}, \[Gamma][a]\[CenterDot] CD[a, f]];
EoMStandardForm[Transp@ Bar@ f:Field[_, Fermion, __]]:= Module[{a}, Transp@ \[Gamma][a]\[CenterDot] Transp@ Bar@ CD[a, f]];
EoMStandardForm[Field[f_, Vector[Index[mu_, type_]], ind_, {}]]:=
	Module[{nu}, FieldStrength[f, {Index[nu, type], Index[mu,type]}, ind, {Index[nu, type]}] ];
EoMStandardForm[Bar@ Field[f_, Vector[Index[mu_, type_]], ind_, {}]]:=
	Module[{nu}, Bar@ FieldStrength[f, {Index[nu, type], Index[mu,type]}, ind, {Index[nu, type]}] ];


(* ::Text:: *)
(*To put scalar kinetic terms on textbook form *)


CanonizeKineticalOps[expr_] :=
	expr/. {
		(*Scalars*)
		Operator[Bar@ Field[f1_, Scalar, inds1_, {}], EoM@ Field[f2_, Scalar, inds2_, {}]]:>
			-Operator[Bar@ CD[Global`d$$1, f1@@inds1], CD[Global`d$$1,f2@@inds2]],
		Operator[Field[f1_, Scalar, inds1_, {}], EoM@ Field[f2_, Scalar, inds2_, {}]]:>
			-Operator[CD[Global`d$$1, f1@@inds1], CD[Global`d$$1,f2@@inds2]],
		(*Vectors*)
		Operator[Bar@ Field[f1_, Vector@ \[Mu]_, inds1_, {}], EoM@ Field[f2_, Vector@ \[Mu]_, inds2_, {}]]:>
			-1/2* Operator[Bar@ FieldStrength[f1, {Index[Global`d$$1, Lorentz], Index[Global`d$$2, Lorentz]}, inds1, {}], FieldStrength[f2, {Index[Global`d$$1, Lorentz], Index[Global`d$$2, Lorentz]}, inds2, {}]],
		Operator[Field[f1_, Vector@ \[Mu]_, inds1_, {}], EoM@ Field[f2_, Vector@ \[Mu]_, inds2_, {}]]:>
			-1/2* Operator[FieldStrength[f1, {Index[Global`d$$1, Lorentz], Index[Global`d$$2, Lorentz]}, inds1, {}], FieldStrength[f2, {Index[Global`d$$1, Lorentz], Index[Global`d$$2, Lorentz]}, inds2, {}]]
	};


(* ::Subsubsection::Closed:: *)
(*Converting to Operator Form*)


NormalToOperatorForm[expr_]:=Operator@expr;


(* ::Subsubsection::Closed:: *)
(*Operator invariants *)


(* ::Text:: *)
(*The operators collected into $operators are assumed to satisfy certain invariants: *)
(*1) There cannot be multiple GammaM in the same DiracProduct (in Operator form, in normal form a lone gamma matrix that's part of a Fermion EOM is allowed)*)
(*2) There cannot be any contractable CGs present in an operator, i.e., two CGs that share a dummy index. *)
(*3) There cannot be any Metrics or Deltas in the operator*)


OperatorSatisfiesInvariantsQ@ op_Operator:= Module[{},
	FreeQ[op, Metric|Delta] && 
		FreeQ[op, CG[_del, _]] &&
		FreeQ[op, DiracProduct[___, _GammaM, _GammaM, ___]] &&
		DuplicateFreeQ[Cases[op, CG[_, inds_]:> inds, Infinity]/. Bar-> Identity]
];


(* ::Text:: *)
(*Enforce operator invariants for all parts of an expression*)


EnforceOperatorInvariants= RefineDiracProducts@* ContractDelta@* ContractCGs;


(* ::Section:: *)
(*Operator construction*)


(* ::Subsubsection::Closed:: *)
(*Operator association*)


(* ::Text:: *)
(*Global associations with operator pattern substitutions (both ways)*)
(*$operators = <| opType -> <| operator ID -> OperatorProperties |> |>*)


$operators= <||>;
ResetOperatorAssociations[]:= Module[{class},
	Do[
		UnsetHermitianIdentities[class, _];
	,{class, Keys@ $operators}];
	$operators= <||>;
	$compoundOperators= <|dDimensional-> <||>, Evanescent-> <||>, FourDimensional-> <||> |>;
];
ResetOperatorAssociations[];


(* ::Subsubsection::Closed:: *)
(*General chiral basis setting *)


(* ::Text:: *)
(*Set the operator basis preference for operators involving VLF fermion lines. *)
(*	Automatic: VLF fermions for the kinetic and mass terms; chiral elsewhere.*)
(*	Chiral: Chiral spinor basis everywhere.*)
(*	VectorLike: All spin chains with VLFs are in the vector basis*)


$VLFBasis= Automatic;


UpdateChiralBasis::unkwnbasis= "`1` is not an accepted option for the VLFBasis.";
UpdateChiralBasis::notimpl= "VectorLike option for VLF basis is not yet implemented. Keep an eye out for future updates.";


UpdateChiralBasis[newVLFBasis_]:= Module[{},
	If[newVLFBasis === $VLFBasis, Return[]; ];
	If[newVLFBasis === VectorLike,
		Message[UpdateChiralBasis::notimpl];
		Abort[];
	];
	If[!MatchQ[newVLFBasis, Chiral|Automatic],
		Message[UpdateChiralBasis::unkwnbasis, newVLFBasis];
		Abort[];
	];
	$VLFBasis= newVLFBasis;
	ResetOperatorAssociations[];	
];


(* ::Subsection:: *)
(*Utility functions*)


(* ::Subsubsection::Closed:: *)
(*Separate a term by operator and coefficient, returns {coupling, operator}*)


SeparateInteractionTerm::IncompatibleForm = "The object `1` is not of the form (Operator[a]*b)"


SeparateInteractionTerm[dL_]:=Module[{L = Collect[Operator@Expand@dL,_Operator],cpl,op},
	If[MatchQ[L, _Operator b_:1],
		cpl= L/._Operator -> 1;
		op = L/cpl;
		Return[{cpl, (*RelabelIndices@*)OperatorToNormalForm[op,Unique->True,CanonizeKinetic->True]}];
	,
		Message[ SeparateInteractionTerm::IncompatibleForm ,Format[L,NiceForm]];
		Return[{1,dL}]
	];
]


(* ::Subsubsection::Closed:: *)
(*Indices*)


(* ::Text:: *)
(*Return all open/contracted indices appearing in an object *)


OpenIndices@ obj_:= Cases[Tally@ Cases[obj, Index[__], All], {ind_, 1}-> ind];
ContractedIndices@ obj_:= Cases[Tally@ Cases[obj, Index[__], All], {ind_, 2}-> ind];


(* ::Text:: *)
(*Constructs canonical dummy indices of the given types *)


ConstructDummyIndices@ types_List:= ConstructDummyIndices@ types=
	MapIndexed[(Index[ToExpression["d$$"<> ToString[First@ #2]], #1]&), types];


(* ::Subsubsection::Closed:: *)
(*FirstElementBy*)


(* ::Text:: *)
(*The performance of MinimalBy for finding a single smallest element is much worse than sort would be (at least for small examples)*)


FirstElementBy[list_, func_]:= First@ SortBy[list, func];
FirstElementBy[func_]@ list_:= FirstElementBy[list, func];


(* ::Subsubsection::Closed:: *)
(*Conjugate operator class*)


(* ::Text:: *)
(*Returns the conjugate of an operator class*)


OpClassConjugate@ {fieldTypes_List, devs_}:= {Sort@ Conj@ fieldTypes, devs};


(* ::Subsubsection::Closed:: *)
(*Self-conjugate operator class *)


(* ::Text:: *)
(*Determines if the operator class is of the self-conjugate type*)


SelfConjugateClassQ@ {fieldTypes_List, _}:=
	Sort@ fieldTypes === Sort@ Conj@ fieldTypes;
SelfConjugateClassQ@ fieldTypes_List:= SelfConjugateClassQ@ {fieldTypes, 0};


(* ::Subsubsection::Closed:: *)
(*Conjugate Atomic/CompOp*)


ConjugateAtomicOperator@ AtomicOp[opClass_, opID_, inds_]:=
	AtomicOp[OpClassConjugate@ opClass, opID, inds[[$operators[opClass, opID, ConjugateFieldExchange]]]];


ConjugateAtomicOperator@ CompOp[red_, opClass_, opID_, inds_]:= Module[{opAssociation},
	opAssociation= $compoundOperators[red, opClass, opID];
	If[opAssociation@ SelfConjugate === False,
		CompOp[red, OpClassConjugate@ opClass, opID, inds[[opAssociation@ ConjugateFieldExchange]]]
	,
		opAssociation@ SelfConjugate* CompOp[red, opClass, opID, inds[[opAssociation@ ConjugateFieldExchange]]]
	]
];


(* ::Subsection:: *)
(*Construct operator matching patterns*)


(* ::Subsubsection::Closed:: *)
(*Construct op ID matching patterns from an operator *)


SetBarable[AtomicOp]


(* ::Text:: *)
(*Set up a pattern to identify an Operator object with its AtomicOp object*)


patternLabels= Table[Unique@ ordPat, 20];


OperatorMatchingPatterns[opClass_, opID_, op_Operator]:= Module[
		{indices, indReplace, pattern, signature},

	{pattern, signature}= MakeOperatorPattern@ op;

	(*Flavor indices *)
	indices= Flatten@ OperatorFieldsAndFlavors[op][[;;, ;;, 2]];
	indReplace= IndexPatternReplace/@ indices;

	(*Produce substitution rules from operator pattern*)
	{
		(*Replace operator with AtomicOp object*)
		With[{temp= signature,
				inds= indReplace[[;;, 2, 1]],
				pat= pattern/. indReplace},
			RuleDelayed[pat, Activate@ temp AtomicOp[opClass, opID, inds]]
		]
	,
		(*Function for matching operator with specific flavor indices*)
		With[{temp= signature,
				pat= pattern/. Thread@ Rule[indices, Slot/@ Range@ Length@ indices]},
			Function[Evaluate@ RuleDelayed[pat, Activate@ temp]]
		]
	}
];


IndexPatternReplace@ ind:Index[lab_, rep_]:= With[{temp= Unique@ lab},
	Rule[ind, temp:Index[_, rep]]
];


(* ::Subsubsection::Closed:: *)
(*Make operator matching pattern*)


(* ::Text:: *)
(*Function for creating pattern and associated signature from an Operator*)


MakeOperatorPattern@ op_Operator:= Module[{antisyms, cgsyms, indices, pattern, lorentzSign, counter= 1},
	(*Make pattern out of contracted indices*)
	indices= IndexPatternReplace/@ ContractedIndices@ op;
	pattern= op/. indices/. $lorentzSubs;
	lorentzSign= LorentzSign@ op;

	(*Name the antisymmetrized pattern indices*)
	pattern= pattern/. {obj: (_LCTensor|_FieldStrength|_GammaM):>
		(obj/. x_OrderlessPatternSequence:> Pattern[Evaluate@ patternLabels[[counter++]], x])};

	(*Provide sign signatures for all the anti-symmetric indices*)
	antisyms= Times@@ Cases[pattern,
			Verbatim[Pattern][name_, pats_OrderlessPatternSequence]:>
				Inactive[Signature]@ {name} Inactive[Signature][List@@ pats/. Verbatim[Pattern][x_, Blank[]]:>x]
		, All];

	(*Account for the symmetries of the CGs and spin-chain transposition (for identical fermions)*)
	cgsyms= Times@@ Reap[Sow@ Nothing;
			pattern= pattern/. cg_CG:> CGPattern@ cg;
			pattern= pattern/. spinChain_NCM:> SpinLinePattern@ spinChain;
		][[2, 1]];

	{pattern, lorentzSign* antisyms* cgsyms}
]


(* ::Text:: *)
(*Create orderless Patterns for the anti-symmetric Lorentz constructions*)


$lorentzSubs= {
	LCTensor[\[Mu]___]:> LCTensor@ OrderlessPatternSequence@ \[Mu],
	FieldStrength[A_, {\[Mu]__}, a__]:> FieldStrength[A, {OrderlessPatternSequence@ \[Mu]}, a],
	GammaM[\[Mu]_, \[Nu]__]:> GammaM@ OrderlessPatternSequence[\[Mu], \[Nu]]
};


(* ::Text:: *)
(*Compensate for the canonical ordering of patterns when creating OrderlessPatternSequence*)


LorentzSign@ op_:=
	Times@@ Cases[op, LCTensor[\[Mu]___]:> Signature@ List@ \[Mu], Infinity] *
	Times@@ Cases[op, FieldStrength[_, {\[Mu]__}, __]:> Signature@ List@ \[Mu], Infinity] *
	Times@@ Cases[op, GammaM[\[Mu]_, \[Nu]__]:> Signature@ {\[Mu], \[Nu]}, Infinity]


(* ::Text:: *)
(*Identify symmetric contractions of the CGPatterns *)


CGPattern@ CG[cg_, indsPat_]:= Module[{symsPat, inds, syms, x},
	syms= $CGproperties[cg/. Bar-> Identity, Symmetries];
	If[Length@ syms === 1,
		Return@ CG[cg, indsPat]
	];

	symsPat= GroupBy[KeyMap[indsPat[[#]]&, syms], Identity, Keys];
	(*Sow signature for the RHS of the AtomicOp if anti-symmetries are present*)
	If[KeyExistsQ[symsPat, -1],
		Sow@ Inactive[DidMatchSwitch][1, -1][x]
	];

	CG[cg, (x: Alternatives@@ Lookup[symsPat, 1, {}]) | Alternatives@@ Lookup[symsPat, -1, {}] ]
]


(* ::Text:: *)
(*Identify (anti-)symmetries from transposition of spin-chains with identical fermions *)


SpinLinePattern[spinLine:NCM[Transp@ Field[lab_, rest1__],
	dirac:RepeatedNull[DiracProduct[GammaCC,___], 1], Field[lab_, rest2__]] ]:= 
	Module[{l},
		(l: spinLine)| NCM[Transp@ Field[lab, rest2], TranspDiracPattern[l, dirac], Field[lab, rest1]]
	];
SpinLinePattern[spinLine:NCM[Bar@ Field[lab_, rest1__],
	dirac:RepeatedNull[DiracProduct[GammaCC,___], 1], Bar@ Transp@ Field[lab_, rest2__]] ]:= 
	Module[{l},
		(l: spinLine)| NCM[Bar@ Field[lab, rest2], TranspDiracPattern[l, dirac], Transp@ Bar@ Field[lab, rest1]]
	];
(*EOMs: both fields must be EOMS for an ambiguious ordering*)
SpinLinePattern[spinLine:NCM[EoM@ Transp@ Field[lab_, rest1__],
	dirac:RepeatedNull[DiracProduct[GammaCC,___], 1], EoM@ Field[lab_, rest2__]] ]:= 
	Module[{l},
		(l: spinLine)| NCM[EoM@ Transp@ Field[lab, rest2], TranspDiracPattern[l, dirac], EoM@ Field[lab, rest1]]
	];
SpinLinePattern[spinLine:NCM[EoM@ Bar@ Field[lab_, rest1__],
	dirac:RepeatedNull[DiracProduct[GammaCC,___], 1], EoM@ Bar@ Transp@ Field[lab_, rest2__]] ]:= 
	Module[{l},
		(l: spinLine)| NCM[EoM@ Bar@ Field[lab, rest2], TranspDiracPattern[l, dirac], EoM@ Transp@ Bar@ Field[lab, rest1]]
	];
SpinLinePattern@ else_:= else;


(* ::Text:: *)
(*Creates patterns for the transposed Dirac algebra, and sowing the relative sign (if any)*)


TranspDiracPattern[symb_Symbol, dp_DiracProduct]:= Module[{gammas, chiral, proj, sign},
	{gammas, chiral}= Replace[dp, {
			DiracProduct[GammaCC, g___, c:(_Proj|Gamma5)]:> {Reverse@ {g}, c},
			DiracProduct[GammaCC, g___]:> {Reverse@ {g}, 1}
		}];
	sign= Times@@ Cases[gammas, (GammaM@ Verbatim[Pattern][_, ind_]|Transp@ GammaM@ Verbatim[Pattern][_, ind_]):> 
		 If[Head@ ind === OrderlessPatternSequence, Power[-1, Ceiling[Length@ ind/2]], -1 ]];
	
	chiral= ReleaseHold@ Switch[chiral
		,1, 
			Hold@ Sequence[]
		,Gamma5,
			sign*= Times@@ Cases[gammas, (GammaM@ Verbatim[Pattern][_, ind_]|Transp@ GammaM@ Verbatim[Pattern][_, ind_]):> 
				If[Head@ ind === OrderlessPatternSequence, Power[-1, Length@ ind], -1]];
			Gamma5
		,_Proj,
			proj= First@ chiral;
			proj*= Times@@ Cases[gammas, (GammaM@ Verbatim[Pattern][_, ind_]|Transp@ GammaM@ Verbatim[Pattern][_, ind_]):> 
				If[Head@ ind === OrderlessPatternSequence, Power[-1, Length@ ind], -1]];
			Proj[proj]
		];
	
	If[sign === -1, Sow@ Inactive[DidMatchSwitch][1, -1][symb]];
	DiracProduct[GammaCC, Sequence@@ gammas, chiral]
]


(* ::Subsubsection::Closed:: *)
(*Return Operator from AtomicOp *)


(* ::Text:: *)
(*Set up a pattern to expand an  AtomicOp object to the corresponding operator *)


OperatorExpansionPattern[class_, id_, op_Operator]:= Module[{indices},
	indices= Flatten@ OperatorFieldsAndFlavors[op][[;;, ;;, 2]];
	With[{inds= indices},
		RuleDelayed[AtomicOp[class, id, flavInds_], op/. Thread@ Rule[inds, flavInds] ]
	]
];


(* ::Subsection:: *)
(*Operator properties*)


(* ::Subsubsection::Closed:: *)
(*All operator properties *)


(* ::Text:: *)
(*Characterizes all properties of an operator object *)


OperatorProperties::nonop= "Passed a non-operator \"`1`\".";
OperatorProperties::nonminop= "The operator `1` is non minimal and cannot be used internally.";


OperatorProperties[_, _, x:Except[_Operator]]:= (Message[OperatorProperties::nonop, x]; Abort[];);


OperatorProperties[opClass_, opID_, op_Operator]:= Module[{count= 1, conjIndExchange, conjFieldExchange, equivClasses,
		fieldTypes, flavorInds, flavorPerms, indexGrouping, indexTypes, iPower, opBar, opIDpattern, permPattern, reexpPattern,
		selfConjugate, selfconjugateType, symmetries},
	If[!OperatorSatisfiesInvariantsQ@ op, 
		Message[OperatorProperties::nonminop, op];
		Abort[];
	];
	
	fieldTypes= OperatorFieldsAndFlavors@ op;
	flavorInds= fieldTypes[[;;, ;;, 2]];
	fieldTypes= fieldTypes[[;;, ;;, 1]];
	indexTypes= Flatten@ flavorInds/. Index[_, type_]-> type;
	(*The (flavor) index grouping contains the fields by type and the enumerated flavor indices for each field.
		{{{1,2},{3,4}},{{5},{6}}} means that there are two identical fields with indices {1,2} and {3,4}, and two identical fields with indices {5} and {6}.*)
	indexGrouping= flavorInds/. _Index:> count++;

	(*Get patterns for matching the operator*)
	{opIDpattern, permPattern}= OperatorMatchingPatterns[opClass, opID, op];
	reexpPattern= OperatorExpansionPattern[opClass, opID, op];

	(*Valid permutation of flavor indices*)
	flavorPerms= Flatten/@ Tuples[Permutations/@ indexGrouping];
	(*Flavor symmetries*)
	flavorInds= Flatten@ flavorInds;
	symmetries= DeleteCases[Rule[#1, op/. {permPattern@@ #2, _-> 0} ]& @@@ ({#, flavorInds[[#]]}&/@ flavorPerms),
		Rule[_, 0]];

	(*Equivalence classes for the flavor indices*)
	equivClasses= Sort/@ Outer[Part, flavorPerms, symmetries[[;;, 1]], 1];
	equivClasses= DeleteDuplicates@ equivClasses[[;;, 1]];

	(*Selfconjugate*)
	selfconjugateType= SelfConjugateClassQ@ Flatten@ fieldTypes;
	opBar= Bar@ op;
	selfConjugate= !FreeQ[opBar, First@ opIDpattern];
	conjFieldExchange= If[selfconjugateType, FindPermutationOrder[fieldTypes, Conj@ fieldTypes], {}];
	conjIndExchange= If[selfConjugate,
			FindPermutationOrder@@ Cases[{op, Bar@ op}/. opIDpattern, AtomicOp[__, inds_]-> inds, All]
		,
			(*For selfconjugate-type operators, the permutation order can only be found once the conjugate operator has been added.
			This information is then added after the fact.*)
			{}
		];
	
	(*Determine if the operator definition should include a factor of i*)
	iPower= If[selfConjugate,
			MatchQ[opBar/. opIDpattern, Times[-1, _AtomicOp]] 
		,
			OperatorIPower@ op
		];
	If[iPower,
		reexpPattern[[2, 1]]*= I;
		(*Prevent evaluation of r.h.s. of the RuleDelayed while multiplying with -I*)
		opIDpattern= opIDpattern/. a_AtomicOp-> -I* a;
	];

	<|
		AtomicOpExpansionPattern-> reexpPattern,
		ConjugateFieldExchange-> conjFieldExchange,
		ConjugateIndexExchange-> conjIndExchange,
		ID-> {opClass, opID},
		IndexGrouping-> indexGrouping,
		IndexType-> indexTypes,
		InequivalentPermutations-> equivClasses,
		OperatorMatchingPattern-> opIDpattern,
		SelfConjugate-> selfConjugate,
		SelfconjugateType-> selfconjugateType,
		Score-> OpScore[op, selfConjugate],
		Subclass-> OperatorSubclass@ op,
		Symmetries-> symmetries
	|>
]


(* ::Text:: *)
(*Determines the appropriate power of i for an operator based on its factors*)


OperatorIPower@ op_Operator:= Module[{p},
	p= Count[op, _FieldStrength, Infinity]; (*Number of field-strength tensors*)
	p+= Total@ Cases[op, (Field|FieldStrength)[__, devs_]:> Length@ devs, All];
	p+= Total@ Count[op, eom_EoM/; FreeQ[eom, Scalar|Vector], All];
	Mod[p, 2] === 1
];


(* ::Text:: *)
(*Gives the conjugate index permutation of an operator (may be inferable from the properties in $operators now)*)


OpConjFlavorPermutation@ op_Operator:= Module[{},
	FindPermutationOrder@@ Cases[MatchOperatorPatterns@ {op, Bar@ op}, AtomicOp[__, inds_]-> inds, All]
]


(* ::Subsubsection::Closed:: *)
(*Specific operator tests *)


(* ::Text:: *)
(*Tests if an Operator is a canonical kinetic term*)


CanonicalKineticOpQ[Times[op_Operator, __]]:= CanonicalKineticOpQ@ op;


CanonicalKineticOpQ[op_Operator]:= MatchQ[op, Alternatives[
		(*Scalars*)
		Operator[Bar@Field[_, Scalar, _, {}], EoM@ Field[_, Scalar, _, {}]],
		Operator[Field[_, Scalar, _, {}], EoM@ Field[_, Scalar, _, {}]],
		(*Fermions*)
		Operator[Bar@ Field[_, Fermion, _, {}]\[CenterDot] EoM@ Field[_, Fermion, _, {}] ],
		Operator[Bar@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct@ (_Proj|Gamma5)\[CenterDot] 
			EoM@ Field[_, Fermion, _, {}] ],
		Operator[Transp@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct[GammaCC]\[CenterDot] EoM@ Field[_, Fermion, _, {}] ],
		Operator[Transp@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct[GammaCC, _Proj]\[CenterDot] 
			EoM@ Field[_, Fermion, _, {}] ],
		(*Vectors*)
		Operator[Bar@ Field[_, Vector@ \[Mu]_, _, {}], EoM@ Field[_, Vector@ \[Mu]_, _, {}]],
		Operator[Field[_, Vector@ \[Mu]_, _, {}], EoM@ Field[_, Vector@ \[Mu]_, _, {}]]
		(*Operator[FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}],
			FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}]],
		Operator[Bar@ FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}],
			FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}]]*)
		] ];


(* ::Text:: *)
(*Tests if an operator is generally a kinetic operator*)


KineticOpQ[Times[op_Operator, __]]:= KineticOpQ@ op;
KineticOpQ[op_Operator]:= CanonicalKineticOpQ@ op || MatchQ[op, Alternatives[
		(*Scalars*)
		HoldPattern@ Operator[Bar@Field[_, Scalar, _, {\[Mu]_}], Field[_, Scalar, _, {\[Mu]_}]],
		HoldPattern@ Operator[Field[_, Scalar, _, {\[Mu]_}], Field[_, Scalar, _, {\[Mu]_}]],
		HoldPattern@ Operator[Field[_, Scalar, _, {}], EoM@ Bar@ Field[_, Scalar, _, {}]],
		(*Fermions*)
		Operator[EoM@ Bar@ Field[_, Fermion, _, {}]\[CenterDot] Field[_, Fermion, _, {}] ],
		Operator[EoM@ Bar@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct@ (_Proj|Gamma5)\[CenterDot] 
			Field[_, Fermion, _, {}] ],
		Operator[EoM@ Transp@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct[GammaCC]\[CenterDot] Field[_, Fermion, _, {}] ],
		Operator[EoM@ Transp@ Field[_, Fermion, _, {}]\[CenterDot] DiracProduct[GammaCC, _Proj]\[CenterDot] 
			Field[_, Fermion, _, {}] ],
		(*Vectors*)
		Operator[Field[_, Vector@ \[Mu]_, _, {}], EoM@ Bar@ Field[_, Vector@ \[Mu]_, _, {}]],
		Operator[FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}], FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}]],
		Operator[Bar@ FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}], FieldStrength[_, {\[Mu]_, \[Nu]_}, _, {}]]
	] ];


(* ::Text:: *)
(*Tests if operator is fermion mass term is canonical *)


CanonicalFermionMassOpQ= MatchQ[
	(*VL fermion*)
	Operator@ NCM[Bar@ Field[l1_, Fermion, _, {}], Field[l1_, Fermion, _, {}]] |
	(*Majorana fermion*)
	Operator@ NCM[Transp@ Field[l1_, Fermion, _, {}], DiracProduct@ GammaCC, Field[l1_, Fermion, _, {}]]
];


(* ::Subsubsection::Closed:: *)
(*Operator Type*)


(* ::Text:: *)
(*Head for denoting Conjugate field type*)


Conj@ Conj@ f_:= f;
Conj@ l_List:= Conj/@ l;


(* ::Text:: *)
(*Operator class discriminator *)


OperatorClass@ op_Operator:= Module[{fields, temp,
		devs= 0, gaugeFields= List@@ Query[All, Key@ Field]@ $GaugeGroups},
	devs+= Plus@@ Cases[op, EoM@ (Bar|Transp)@ Field[_, type_, __]:> EOMDevs@ type, All];
	devs+= Plus@@ Cases[op, EoM@ Transp@ Bar@ Field[_, Fermion, __]:> EOMDevs@ Fermion, All];
	devs+= Plus@@ Cases[op, EoM@ Field[_, type_, __]:> EOMDevs@ type, All];
	devs+= Plus@@ Cases[op, FieldStrength[___, linds_]:>
		1+ Length@ linds, All];
	devs+= Plus@@ Cases[op, Field[___, linds_]:> Length@ linds, All];

	temp= op/. EoM-> Identity/. FieldStrength[lab_, inds_, rest__]:> Field[lab, Vector@ First@ inds, rest];
	devs+= Plus@@ Cases[temp, Field[lab_, _Vector, __]/; MemberQ[gaugeFields, lab]-> 1, All];
	temp= DeleteCases[temp, Field[lab_, _Vector, __]/; MemberQ[gaugeFields, lab], All];

	fields= Cases[temp, Bar@ Field[lab_, (Scalar|Fermion|_Vector),__]:> Conj@ lab, All];
	temp= DeleteCases[temp, Bar@ Field[_, (Scalar|Fermion|_Vector),__], All];
	fields= Join[fields, Cases[temp, Field[lab_, (Scalar|Fermion|_Vector),__]:> lab, All]];
	{Sort@ fields, devs}
];
EOMDevs[type_]:= Switch[type, Scalar, 2, Fermion, 1, _Vector, 2];


(* ::Text:: *)
(*Returns an ordered list of field types and corresponding flavor indices *)


OperatorFieldsAndFlavors@ op_Operator:= Module[
		{fields, flavorReps= Keys@ $FlavorIndices, temp= op},
	fields= Cases[temp, (Bar@ Field[lab_, _, inds_, _] | Bar@ FieldStrength[lab_/;!GaugeFieldQ@ lab, _, inds_, _]):>
		{Conj@ lab, Cases[inds, Index[_, rep_/; MemberQ[flavorReps, rep]] ]}, All];
	temp= DeleteCases[temp, Bar@ (_Field | _FieldStrength), All];
	fields= Join[fields,
		Cases[temp, (Field[lab_, _, inds_, _] | FieldStrength[lab_/;!GaugeFieldQ@ lab, _, inds_, _]):>
			{lab, Cases[inds, Index[_, rep_/; MemberQ[flavorReps, rep]] ]}, All] ];
	GatherBy[SortBy[fields, First], First]
]


(* ::Text:: *)
(*Returns the subclass of the operator. This is all the "local" information of the operator, that is, what the unique details that can be determined w/o knowing about dummy-index contractions *)


OperatorSubclass@ op_Operator:= Module[{cgs, eoms, fields, fss, temp},
	cgs= Cases[op, CG[cg_, _]:> cg, All];
	
	eoms= Join[Cases[op, (EoM@ Bar@ Field[f_, __] | EoM@ Transp@ Bar@ Field[f_, Fermion, __]):> Conj@ f, All],
		Cases[op, (EoM@ Field[f_, __] | EoM@ Transp@ Field[f_, Fermion, __]):> f, All]];
	temp= DeleteCases[op, _EoM, All];
	
	(*Fields and their derivatives*)
	fields= Cases[temp, Bar@ Field[f_, _, _, devs_]:> {Conj@ f, Length@ devs}, All];
	temp= DeleteCases[temp, Bar@ _Field, All];
	fields= Join[fields, Cases[temp, Field[f_, _, _, devs_]:> {f, Length@ devs}, All]];
	
	(*FieldStrengths  and their derivatives*)
	fss= Cases[temp, Bar@ FieldStrength[f_, _, _, devs_]:> {Conj@ f, Length@ devs}, All];
	temp= DeleteCases[temp, Bar@ _FieldStrength, All];
	fss= Join[fss, Cases[temp, FieldStrength[f_, _, _, devs_]:> {f, Length@ devs}, All]];
	
	Sort/@ {cgs, eoms, fields, fss}
];


(* ::Subsubsection::Closed:: *)
(*Operator ranking*)


(* ::Text:: *)
(*Ranks operators by preference (higher is more preferable)*)
(*Canonical terms: +10000 (You want these in the basis!)*)
(*EOMS: +1000 (should always be preferred)*)
(*Evanescent: +100 (keep if possible, they are also removable) *)
(*Other adjustments: [-10, +10] (ordering by preference)*)


OpScore::unexp= "OpScore received unexpected argument `1`"


OpScore[op_Operator, selfConj_]:= Module[{cgs, cgInds, score= 0, composedCGs= Keys@ $CGreplacements},
	(*Score by chiral basis*)
	score= ChiralScoring@ op;
	
	(*Check for canonical terms*)
	If[CanonicalKineticOpQ@ op, Return[score+ 20000]] (*Should rank higher than just EOM count to ensure canonical kinetic terms*);

	(*EoMs can be removed by field redefinitions and any EoM should always be included in the basis*)
	score+= 10000 Count[op, _EoM, Infinity];
	(*Field strength tensors are preferable to more derivatives*)
	score+= 1 Count[op, _FieldStrength, Infinity];
	(*Penalize longer Dirac structures*)
	score+= -.01 Total@ Cases[op, DiracProduct[___, g:HoldPattern@ _GammaM, ___]:> Length@ g, Infinity];
	(*Penalize non-ordered CDs*)
	score+= -.1 Count[op, Field[__, {___, a_, __, a_, ___}], Infinity];
	(*Prefer when derivatives and vectors share indices*)
	score+= .1 Count[op, Field[_, Vector[a_], _, {___, a_}], Infinity];
	
	(*Penalize CGs, to prefer simpler cgs and deltas, i.e., no cgs above all*)
	cgs= Cases[op, _CG]; 
	If[cgs =!= {}, 
		score+= -.01 Length@ cgs;
		{cgs, cgInds}= Transpose[List@@@ cgs/. Bar->Identity];
		
		score+= -.1 Count[cgs, _eps];
		score+= -.1 Count[cgs, name_/; MemberQ[composedCGs, name]];
		
		(*Penalize CGs with derivative fields*)
		cgInds= Flatten@ cgInds;
		score+= -.03 Total@ Cases[op, Field[_, _, grInds_, devInds_]/; IntersectingQ[grInds/. Bar->Identity, cgInds]:>
			Length@ devInds, Infinity];
	];
	
	(*Preferance of self-conjugate operators*)
	If[selfConj, score+= .3];
	(*Score for different combinations of fermions bilinears*)
	score+= FierzScore[op];
	
	score
]


OpScore@ x_:= (Message[OpScore::unexp, x]; Abort[];);


(* ::Text:: *)
(*Give preference to chiral or vector/axial basis depending on current VLFBasis*)


ChiralScoring[op_Operator]:= Switch[$VLFBasis
	,Automatic,
		MassOrKinVLFTerm@ op* Count[op, _Proj, Infinity]
	,Chiral,
		(*Preference to chiral projectors*)
		.5* Count[op, _Proj, Infinity]
	,VectorLike,
		(*Penalty to chiral projectors in VL fermion lines*)
		-.5* Count[op, ncm_NCM? VectorLikeSpinChainQ/; !FreeQ[ncm, _Proj], All]
	];


(*-1 for mass/kin terms of vector-like fermions, +1 for everything else*)
MassOrKinVLFTerm@ op_Operator:= Which[
	Count[op, _Field|_FieldStrength, All] === 2 && (*No more than 2 "fields"*)
	FreeQ[op, Field[_, Fermion, _, {__}], All] && (*No non-EOM derivatives*)
	!FreeQ[op, _NCM? VectorLikeSpinChainQ], (*Contains spin-chain with VLF fermions*)
		-.1,
	Count[op, _NCM, All] >= 2, (*Chiral four-fermion operators must be in basis (weighs higher than self-conjugate)*)
		.5, 
	True, 
		.1
	]


(* ::Subsection:: *)
(*Match all operators in an expression to patterns  *)


(* ::Subsubsection::Closed:: *)
(*Search expression*)


(* ::Text:: *)
(*Substitutes all operators in an expression with AtomicOp objects for further manipulations  *)


Options@ MatchOperatorPatterns= {ResetIdentities-> {dDimensional, FourDimensional, Evanescent}};


MatchOperatorPatterns[expr_, OptionsPattern[]]:= Module[{class, classes, ops, out},
	(*Separate out flavor index contractions within each operator*)
	(*Ensure that all operators are DiracRefined*)
	out= OperatorFlavorSeparate/@ Operator/@ TermsToList@ expr;

	(*Groups all un-matched operators in expression*)
	classes= GroupBy[out, FirstCase[#, op_Operator:> OperatorClass@ op, {}, All]& ];

	(*Considering operator classes individually, the replacements are better targeted, which speeds up performance*)
	Sum[Catch[
		If[class === {}, (*Constants are kept unchanged*)
			Throw[Plus@@ classes@ class];
		];
		ops= Cases[classes@ class, _Operator, All];
		(*Construct new patterns out of any new operators encounterd in the sets*)
		MakeNewOperatorPatterns[class, ops, OptionValue@ ResetIdentities];

		Plus@@ classes@ class/. OpToAtomicReplacementPattern@ class
	], {class, Keys@ classes}]

];


(* ::Text:: *)
(*Separates flavor contracted indices within an operator to extract a delta*)


OperatorFlavorSeparate@ expr_:= expr/. op_Operator:> OperatorFlavorSeparate@ op


OperatorFlavorSeparate@ op_Operator:= Module[{flavInds, deltas, out= op},
	(*Find all internally contracted flavor indices*)
	flavInds= DeleteCases[ContractedIndices@ op,
		Index[_, type_/; !MemberQ[Keys@ $FlavorIndices, type]] ];

	(*Extract a delta for every contraction*)
	deltas= Delta[#, MapAt[Unique, #, 1]]& /@ flavInds;
	Do[
		out= ReplaceFirst[out,  replacement]
	, {replacement, Rule@@@ deltas}];
	out Times@@ deltas
];


(* ::Text:: *)
(*For matching patterns in IBP identities (performance optimized)*)
(*identityTypeResets is a subset of {dDimensional, FourDimensional, Evanescent}*)


MatchOperatorPatternsInIdentities[opClass_, identities_List, identityTypeResets_List]:= Module[
		{len, newOps, ops, out, repl, subclasses, temp},
	out= Operator@ identities;

	(*Speed up replacements by first substituting in dummy symbols for each unique operator*)
	ops= DeleteDuplicates@ Cases[out, _Operator, All];
	ops= MapIndexed[Rule[#1, temp@@ #2]&, ops];
	out= out/. ops;

	(*Make replacements of dummies with atomic operators that are already known*)
	len= Length@ $operators[opClass];
	
	ops= GroupBy[ops[[;;, {2, 1}]], OperatorSubclass@* Last];
	repl= OpToAtomicReplacementPattern@ opClass;
	subclasses= GroupBy[Range@ len, $operators[opClass, #, Subclass]&];
	
	ops= Flatten@ Table[
		ops@ subclass/. repl[[Lookup[subclasses, Key@ subclass, {}]]]
	, {subclass, Keys@ ops}];
	
	(*Construct new patterns out of any new operators encounterd*)
	newOps= Cases[ops, _Operator, All];
	MakeNewOperatorPatterns[opClass, newOps, identityTypeResets];

	(*Make replacements of the final dummies with the new replacement patterns*)
	ops= ops/. OpToAtomicReplacementPattern[opClass][[len+1;;]];

	out/.ops
];


(* ::Subsubsection::Closed:: *)
(*Creates new operator classes*)


(* ::Text:: *)
(*Makes new pattern rules for a list of operators of a given type.*)
(*resetIdentities is a subset of {dDimensional, Evanescent, FourDimensional}*)


MakeNewOperatorPatterns::outofsync= "Operator Class `1` and its conjugate is out of sync."


MakeNewOperatorPatterns[opClass_, opList_List, resetIdentities_List]:= Module[{conjOpClass, identifiers, newOps, nextID, op, perm,
		opBar, classAssociation, classAssociationConj, pat, notSelfConjugateClass, remainingOps= opList, props, propsConj},
	(*Find patterns for the previously encountered operators*)
	classAssociation= Lookup[$operators, Key@ opClass, <||>];

	(*Remove operators that have already been added to $operators*)
	(*This has much better performance in some cases than the equivalent DeleteCases*)
	If[Length@ classAssociation =!= 0,
		remainingOps= Cases[DeleteDuplicates@ remainingOps/. OpToAtomicReplacementPattern[opClass], _Operator];
		If[Length@ remainingOps === 0, Return[];];
	];

	(*Check class is syncronized with its conjugate*)
	If[(notSelfConjugateClass= !SelfConjugateClassQ@ opClass),
		classAssociationConj= Lookup[$operators, Key[conjOpClass= OpClassConjugate@ opClass], <||>];
		If[Length@ classAssociation =!= Length@ classAssociationConj,
			Message[MakeNewOperatorPatterns::outofsync, opClass];
		];
	];

	(*Loops through the unidentitfied operators, creating new patterns for them*)
	nextID= Length@ classAssociation+ 1; pat= Alternatives[];
	newOps= Association@@ Reap[Do[
			If[MatchQ[pat]@ op, Continue[]; ];
			props= OperatorProperties[opClass, nextID, op];
			AppendTo[pat, First@ props@ OperatorMatchingPattern];
			opBar= Bar@ op; (*May be +/- the Bar operator*)
			opBar= Replace[opBar, Times[-1, o_Operator]-> o];
			(*If valid, appends the conjugate operator to the conjugate class*)
			Which[notSelfConjugateClass,
				classAssociationConj[nextID]=
					OperatorProperties[conjOpClass, nextID, opBar];

				Sow[nextID++-> props];
			,!props@ SelfConjugate, (*Otherwise ensure the conjugate operator is added to present class*)
				propsConj= OperatorProperties[opClass, nextID+1, opBar];
				AppendTo[pat, First@ propsConj@ OperatorMatchingPattern];

				(*Rectify the conjugate permutation order between the operator and its conjugate*)
				perm= FindPermutationOrder@@ Cases[{op/. props@ OperatorMatchingPattern, opBar/. propsConj@ OperatorMatchingPattern},
					AtomicOp[__, inds_]-> inds, All];
				props@ ConjugateIndexExchange= perm;
				propsConj@ ConjugateIndexExchange= InversePermutationOrder@ perm;

				Sow[nextID++-> props];
				Sow[nextID++-> propsConj];
			,True,
				Sow[nextID++-> props];
			];
		, {op, remainingOps}]][[2, 1]];

	(*Add the conjugate operators*)
	If[notSelfConjugateClass,
		$operators@ conjOpClass= classAssociationConj;
	];

	$operators@ opClass= classAssociation~ Join~ newOps;

	(*Reset simplification identities*)
	UnsetHermitianIdentities[opClass, Alternatives@@ resetIdentities];
];


(* ::Subsubsection::Closed:: *)
(*Conversion function*)


ToAtomicForm[expr_]:= MatchOperatorPatterns[expr];


(* ::Subsection:: *)
(*Replacement rules*)


(* ::Subsubsection::Closed:: *)
(*Replacement patterns *)


(* ::Text:: *)
(*Collect all patterns for substituting Operator objects for AtomicOp objects*)


OpToAtomicReplacementPattern[]:= Block[{},
	Flatten@ Apply[List, Query[All, All, Key@ OperatorMatchingPattern]@ $operators, {0, 1}]
]
OpToAtomicReplacementPattern[opType_]:= Block[{},
	Flatten@ Apply[List, Query[{Key@ opType}, All, Key@ OperatorMatchingPattern]@ $operators, {0, 1}]
]


(* ::Text:: *)
(*Collect all patterns for substituting AtomicOp objects back to Operator objects*)


AtomicToOpReplacementPattern[]:= Block[{},
	Flatten@ Apply[List, Query[All, All, Key@ AtomicOpExpansionPattern]@ $operators, {0, 1}]
]
AtomicToOpReplacementPattern[opType_]:= Block[{},
	Flatten@ Apply[List, Query[{Key@ opType}, All, Key@ AtomicOpExpansionPattern]@ $operators, {0, 1}]
]


(* ::Subsubsection::Closed:: *)
(*Operator Canonization*)


(* ::Text:: *)
(*Returns the signed canonical form of an AtomicOp with the given open indices, using the obvious symmetries of the operator*)


FlavorCanonizeInternalOp@ expr_:= expr/. op:(_AtomicOp|_CompOp):> FlavorCanonizeInternalOp@ op;


FlavorCanonizeInternalOp@ AtomicOp[opClass_, opID_, inds_]:= Block[{perms},
	perms= $operators[opClass, opID, Symmetries];
	perms= FirstElementBy[inds[[First@ #]]-> Last@ #&/@ perms, First];
	Last@ perms AtomicOp[opClass, opID, First@ perms]
];


FlavorCanonizeInternalOp@ CompOp[id__, inds_]:= Block[{perms},
	perms= $compoundOperators[id, Symmetries];
	perms= FirstElementBy[inds[[First@ #]]-> Last@ #&/@ perms, First];
	Last@ perms CompOp[id, First@ perms]
];


(* ::Subsubsection::Closed:: *)
(*Canonize flavor indices in expression*)


(* ::Text:: *)
(*Ensure that the AtomicOp/CompOp object in each terms of an expression carries canonical indices*)


ExprFlavorCanonize@ expr:HoldPattern@ Times[couplings__?(FreeQ[Operator]), op:((AtomicOp|CompOp)[__, inds_])]:= Block[
		{indReplace= Thread@ Rule[inds, ConstructDummyIndices@ inds[[;;, 2]]]},
	If[OpenIndicesQ@ expr, Return@ expr];
	RelabelIndices[RelabelIndices[Times@ couplings, Unique-> True]/. indReplace] (op/. indReplace)
];
ExprFlavorCanonize@ expr_:= Block[{out=BetterExpand@ expr},
	If[Head@ out === Plus, ExprFlavorCanonize/@ out, out]
];


(* ::Text:: *)
(*Checks if there are open indices in a single term (assumes that there is no sum over terms)*)


OpenIndicesQ@ expression_:= Block[{expr= PseudoTimes@ expression, list},
	expr = expr /. $DropDiagonalCouplings; (* We need to remove open indices *)
	list = Tally[Cases[expr, _Index, All]];
	Length@ Cases[list, {x_, 1}] > 0
]


(* ::Subsubsection::Closed:: *)
(*Expand all operators from AtomicOp/CompOp to Operator Form*)


(* ::Text:: *)
(*Puts AtomicOp and CompOp on Operator form*)


AtomicToOperatorForm[expr_]:= Module[{class, out, reductionType, terms},
	(*The expression may be collected by operators, so we don't want to expand*)
	terms= If[Head@ expr === Plus, List@@ expr, {expr}];

	(*Consider adding a check that no terms has multiple operators, e.g., from nested collection*)

	(*Group terms by operator class*)
	terms= GroupBy[terms,
		FirstCase[#, AtomicOp[class_, __] | CompOp[_, class_, __] :> class, {}, All]& ];

	reductionType= FirstCase[expr, CompOp[redType_, __]:> redType, dDimensional, All];

	(*Perform the replacements by operator class for performance*)
	out= Sum[Catch[
			If[class === {}, (*Constants are kept unchanged*)
				Throw[Plus@@ terms@ class];
			];
			Plus@@ terms@ class/.
				Flatten@ Apply[List, Query[All, Key@ CompOpReplacement]@
					Lookup[$compoundOperators@ reductionType, Key@ class, <||>]]/.
				AtomicToOpReplacementPattern@ class
		], {class, Keys@ terms}];

	out
]


AtomicToNormalForm[expr_]:=OperatorToNormalForm[AtomicToOperatorForm@ expr,CanonizeKinetic->True];


AtomicToNormalForm[assoc_Association]:=Map[AtomicToNormalForm@# &, assoc, Infinity];


(* ::Subsection:: *)
(*Compound operators *)


(* ::Subsubsection::Closed:: *)
(*Construct relevant compound operators for a class *)


(* ::Text:: *)
(*Constructs all compound operators for an operator class and updates the global association for class *)


ConstructCompoundBasis[reductionType_, opClass_, opsToSymmetrize_]:= Module[{atomicOperators, compoundOperators, indexPattern, symmetryReps},
	(*Don't make compounds for kinetic or mass operators; they get special treatment*)
	If[KineticOrMassClassQ@ opClass,
		$compoundOperators[reductionType, opClass]= <||>;
		Return[];
	];
	atomicOperators= List@@ $operators@ opClass;

	compoundOperators= Join@@ MapIndexed[ConstructCompoundsForOp[#1, MemberQ[opsToSymmetrize, First@ #2]]&,
		atomicOperators];

	(*Remove duplicates (from H.c.) *)
	compoundOperators= DeleteDuplicatesBy[compoundOperators, First@ Sort@ {#, -#}& @*First];

	indexPattern= ConstructDummyIndices@ First[atomicOperators]@ IndexType;
	indexPattern= Index[Pattern[Evaluate@ First@ #, Blank[]], Last@ #]&/@ indexPattern;

	(*Save compounds in association*)
	$compoundOperators[reductionType, opClass]= Association@@
		MapIndexed[(First@ #2->
			Append[Last@ #1, CompOpReplacement-> (CompOp[reductionType, opClass, First@ #2, indexPattern]-> First@ #1)]&),
		compoundOperators];
]


(* ::Text:: *)
(*Determine if operator is of the Kinetic term class *)


KineticOrMassClassQ@ {{f1_, Conj@ f2_Symbol| f2_Symbol}, n_}:=
	Switch[List@@ $FieldAssociation[[{Key@ f1, Key@ f2}, Key@ Type]] 
	,{Fermion..},
		MatchQ[n, 0| 1]
	,{Scalar..},
		MatchQ[n, 0| 2]
	,{Vector..},
		MatchQ[n, 0| 2]
	]
KineticOrMassClassQ@ _:= False;


(* ::Subsubsection::Closed:: *)
(*Construct compounds from atomic operator*)


(* ::Text:: *)
(*Determine all relevant symmetrizations and/or (anti-)Hermitian combinations involving the operator, i.e., the compoundOperators*)


ConstructCompoundsForOp[opProperties_, flavorSymmetrization_]:= Block[{dummyInds, conjIndexOrder, conjIndexReplacement,
		conjugated, indexSyms, indexSymsFull, inequivPerms, n, pos, sym, symmetryCombinations, temp},

	If[(!opProperties@ SelfconjugateType || opProperties@ SelfConjugate)&& !flavorSymmetrization,
		Return@ {};
	];
	(*All flavor permutations of the original operator*)
	inequivPerms= Flatten/@ Tuples[Permutations/@ opProperties@ IndexGrouping];

	If[flavorSymmetrization,
		(*Find all relevant symmetrization of the flavor indices in the operator*)
		indexSymsFull= Tuples[FlavorSymmetryReps/@ opProperties[IndexGrouping]];
		(*Expand out the symmetrized indices*)
		indexSymsFull= Map[Map[Join@@ #[[;;, 1]]-> Times@@ #[[;;, 2]]& ]@* Tuples,
			ExpandSymmetrization@ indexSymsFull];

		(*Normalize symmetries*)
		indexSyms= Function[{perm}, MapAt[#/ Length@ perm&, perm, {All, 2}]]/@ indexSymsFull;
		(*Canonically sort indices in each term by the symmetries of the Op*)
		indexSyms= Outer[(Part[First@ #1, First@ #2]-> Last@#1 Last@#2&),
			indexSyms, opProperties@ Symmetries];
		indexSyms= Map[FirstElementBy[First], indexSyms, {2}];
		(*Combine identical elements*)
		indexSyms= DeleteCases[Normal@ Map[Merge[Total], indexSyms, {1}], Rule[_, 0], {2}];
		pos= Flatten@ Position[indexSyms, {__}, {1}];
		indexSymsFull= indexSymsFull[[pos]];
		indexSyms= indexSyms[[pos]];

		(*Expand with the herited symmetries from the parent operator*)
		indexSymsFull= MergeSymmetries[opProperties@ Symmetries, #]&/@ indexSymsFull;
	,
		indexSymsFull= {opProperties@ Symmetries};
		indexSyms= {{indexSymsFull[[1, 1]]}};
	];

	(*Construct operator form*)
	dummyInds= ConstructDummyIndices@ opProperties@ IndexType;
	symmetryCombinations= Plus@@@ Apply[#2* AtomicOp[Sequence@@ opProperties@ ID, dummyInds[[#1]]]&, indexSyms, {2}];

	(*Return: the operator is not self-conjuagte type, we don't need to Hermitianize it*)
	If[!opProperties@ SelfconjugateType,
		(*Determine the unique flavor permutations, i.e., those that are not redundant by symmetry*)
		inequivPerms= Map[Sort, Transpose@ Outer[Part, inequivPerms, indexSymsFull[[;;, ;;, 1]], 1, 2], {2}];
		inequivPerms= DeleteDuplicates/@ inequivPerms[[;;, ;;, 1]];
		{#[[1]], <|
			ConjugateIndexExchange-> False,
			InequivalentPermutations-> #[[3]],
			Score-> opProperties@ Score +.01,
			SelfConjugate-> False,
			Symmetries-> #[[2]]
		|>}&/@ Transpose@ {symmetryCombinations, indexSymsFull, inequivPerms}// Return;
	];

	(*Construct complex-conjugated version of the operators in hcCombinations*)
	conjugated= ConjugateAtomicBySubclass@ AtomicOp[Sequence@@ opProperties@ ID, dummyInds];
	conjugated= Replace[conjugated, {Times[_, atom_AtomicOp]:> atom, atom_AtomicOp:> atom}];
	conjIndexReplacement= Thread@ Rule[Last@ conjugated, dummyInds];
	conjIndexOrder= FindPermutationOrder[Last@ conjugated, dummyInds];

	If[conjIndexOrder === Range@ Length@ dummyInds && conjugated[[2]] === Last@ opProperties@ ID,
		(*Operator is real: No need to Hermitianize operator*)
		symmetryCombinations= List/@ symmetryCombinations;
	,
		(*Construct the conjugated of the operator, with the proper index-permutations*)
		conjugated= ConjugateAtomicBySubclass@ symmetryCombinations/. conjIndexReplacement;
		symmetryCombinations= Table[
				(*The anti-Hermitian combination is made Hermitian through multiplication with I*)
				Expand@ FlavorCanonizeInternalOp@ Replace[{1, I}* symmetryCombinations[[n]] + {1, -I}* conjugated[[n]],
					{___, 0, ___}-> {symmetryCombinations[[n]]}]
			, {n, Length@ symmetryCombinations}];

		(*Keep only the flavor symmetries of both the op and its h.c. part*)
		indexSymsFull= Table[
			temp= sym;
			temp[[;;, 1]]= conjIndexOrder[[#[[conjIndexOrder]] ]]&/@ temp[[;;, 1]];
			Intersection[temp, sym]
		, {sym, indexSymsFull}];
	];

	(*Determine the unique flavor permutations, i.e., those that are not redundant by symmetry*)
	inequivPerms= Map[Sort, Transpose@ Outer[Part, inequivPerms, indexSymsFull[[;;, ;;, 1]], 1, 2], {2}];
	inequivPerms= DeleteDuplicates/@ inequivPerms[[;;, ;;, 1]];
	(*Return the compunds*)
	Join@@ Table[
		{symmetryCombinations[[m, n]], <|
			ConjugateIndexExchange-> conjIndexOrder,
			InequivalentPermutations-> inequivPerms[[m]],
			Score-> opProperties@ Score + (3-n)* .01, (*Prioritze Hermitian over anti-Hermitian combinations*)
			SelfConjugate-> True,
			Symmetries-> indexSymsFull[[m]]
		|>}
	, {m, Length@ symmetryCombinations}, {n, Length@ symmetryCombinations[[m]]}]
]


(* ::Text:: *)
(*Flavor symmetries: *)
(*To be expanded with multi-index representations *)


(*FlavorSymmetryReps@ inds_List:= Symmetrization[inds, #]&/@ IntegerPartitions@ Length@ inds;*)
FlavorSymmetryReps@ {}:= {Symmetrization[{}, {}]}
FlavorSymmetryReps@ {ind_}:= {Symmetrization[{ind}, {1}]};
FlavorSymmetryReps@ {i1_, i2_}:= {Symmetrization[{i1, i2}, {2}], Symmetrization[{i1, i2}, {1, 1}]};
FlavorSymmetryReps@ inds_List:= {Symmetrization[inds, {}]};


(*To be generalized*)
ExpandSymmetrization@ expr_:= expr/. {
		Symmetrization[inds_, {}]:> {Flatten@ inds-> 1},
		Symmetrization[inds_, {1}]:> {Flatten@ inds-> 1},
		Symmetrization[{inds1_, inds2_}, {2}]:> {Flatten@ {inds1, inds2}-> 1, Flatten@ {inds2, inds1}-> 1},
		Symmetrization[{inds1_, inds2_}, {1, 1}]:> {Flatten@ {inds1, inds2}-> 1, Flatten@ {inds2, inds1}-> -1}
	}


MergeSymmetries::incmp= "The symmetries are incompatible.";
MergeSymmetries[sym1_, sym2_]:= Module[{newSym},
	newSym= DeleteDuplicates@ Flatten@ Outer[#1[[1, #2[[1]]]]-> #1[[2]]* #2[[2]]&, sym1, sym2];
	If[!DuplicateFreeQ[First/@ newSym], Message[MergeSymmetries::incmp]; Abort[]; ];
	newSym
];


(* ::Text:: *)
(*For efficient conjugation of AtomicOps*)


ConjugateAtomicBySubclass::notfound= "The conjugate operator not found.";


ConjugateAtomicBySubclass@ expr:Except[_AtomicOp]:= expr/. op_AtomicOp:> ConjugateAtomicBySubclass@ op;
ConjugateAtomicBySubclass@ op:AtomicOp[class_, id_, _]:= Module[{out, repl, subclass},
	out= Bar[op/. $operators[class, id, AtomicOpExpansionPattern]];
	(*The conjugate might include a sign*)
	subclass= OperatorSubclass@ FirstCase[{out}, _Operator, 
		Abort[], All];
	(*All operator matching patterns of the relevant subclass*)
	repl= List@@ Query[Select[#[Subclass] === subclass &], Key@ OperatorMatchingPattern][
		$operators@ OpClassConjugate@ class];
	out= out/. repl;
	If[FreeQ[out, AtomicOp],
		Message[ConjugateAtomicBySubclass::notfound]; 
		Abort[];
	];
	out
];


(* ::Subsubsection::Closed:: *)
(*Construct compound basis*)


ConstructConjugateCompounds[reductionType_, opClass_]:= Block[{},
	$compoundOperators[reductionType, OpClassConjugate@ opClass]=
		Association@@ KeyValueMap[#1->
			 ConjugateCompound[reductionType, opClass, #1, #2]&, $compoundOperators[reductionType, opClass]];
]


(* ::Text:: *)
(*Constructs the conjugate compound operator from a compound operator*)


ConjugateCompound[reductionType_, opClass_, compID_, compProps_]:= Module[
		{atomBar, atomicID, compReplacement, dummies, indPermutation, sign},
	(*Determines the atomic operator ID used in the definition
	(NB. there is only one atomicID for compounds of a complex class)*)
	{atomicID, dummies}= FirstCase[compProps@ CompOpReplacement,
		AtomicOp[_, id_, inds_]:> {id, ConstructDummyIndices@ inds[[;;, 2]]}, 0, Infinity];
	
	(*Determines the conjugation properties of the constituent atomic operators*)
	atomBar= Bar[AtomicOp[opClass, atomicID, dummies]/. AtomicToOpReplacementPattern@ opClass]/.
		OpToAtomicReplacementPattern@ OpClassConjugate@ opClass;
	sign= If[MatchQ[atomBar, Times[-1, _AtomicOp]], -1, 1];
	indPermutation= FindPermutationOrder[dummies, 
		FirstCase[{atomBar}, AtomicOp[__, inds_]:> inds, {}, All] ];
	
	(*Adjust the indices and operatorIDs used in the original compound operator*)
	$compoundOperators[reductionType, opClass, compID, ConjugateIndexPermutation]= indPermutation;
	
	<|
		Symmetries-> MapAt[Ordering[indPermutation][[#[[indPermutation]]]]&, compProps@ Symmetries, {All, 1}],
		Score-> compProps@ Score,
		SelfConjugate-> False,
		CompOpReplacement-> compProps@ CompOpReplacement/.{
				CompOp[__, indPats_]:> CompOp[reductionType, OpClassConjugate@ opClass, compID, indPats[[indPermutation]]],
				AtomicOp[__, inds_]:> sign* AtomicOp[OpClassConjugate@ opClass, atomicID, inds[[indPermutation]]]
			},
		ConjugateIndexPermutation-> Ordering@ indPermutation
	|>
]


(* ::Subsection:: *)
(*Manipulate expression *)


(* ::Subsubsection::Closed:: *)
(*Separate out constants *)


(* ::Text:: *)
(*Splits constants off from a Lagrangian. *)
(*	Return: {constants, remaining expression}*)


SeparateOutConstants@ expr_:= Module[{consts, remainder},
	consts= Last[Last@ Reap[
		remainder= Replace[expr,
			a_/;FreeQ[a, _Field| _FieldStrength, All]:> (Sow[a]; 0),
			If[Head@ expr === Plus, 1, 0] ]
	],{}];
	{Plus@@ consts, remainder}
]


(* ::Subsubsection::Closed:: *)
(*Collect atomic operators in expr*)


(* ::Text:: *)
(*Collects all identical operators in an expression to the same form *)


Options@ CollectOperators= {ContractOverallDelta->True,NormalForm-> True, Simplify-> True(*, CoefficientFullSimplify->False*)};


CollectOperators[arg_, OptionsPattern[]]? OptionsCheck:= Block[
		{out, consts, expr=LagrangianExpand@ HcExpand@ arg},
	(*subtract off constants*)
	If[FreeQ[expr, _Field| _FieldStrength, All], Return[expr]];
	{consts, expr}= SeparateOutConstants@ expr;

	LagrangianLikeCheck@ expr;

	out= MatchOperatorPatterns@ EnforceOperatorInvariants@ expr;

	out= CollectCoefficients[ExprFlavorCanonize@ out, ContractOverallDelta->OptionValue@ ContractOverallDelta,Simplify-> OptionValue@ Simplify(*, CoefficientFullSimplify->OptionValue@CoefficientFullSimplify*)];

	out= out/. AtomicToOpReplacementPattern[];

	If[OptionValue@ NormalForm, 
		OperatorToNormalForm[out, CanonizeKinetic-> True], out] + consts
];


(* ::Text:: *)
(*For internal use keeping Atomic form*)


Options@ InternalCollectOperators= {
		InternalOpRepresentation-> False,
		NormalForm-> False, 
		TrigSimplify-> False
	}


InternalCollectOperators[arg_, OptionsPattern[]]:= Block[{out},
	out= CollectCoefficients[ExprFlavorCanonize@ MatchOperatorPatterns@ EnforceOperatorInvariants@ arg,
		TrigSimplify-> OptionValue@ TrigSimplify];
		
	If[OptionValue@ InternalOpRepresentation,
		Return@ out
	];
	
	out= AtomicToOperatorForm@ out;
	If[OptionValue@ NormalForm,
		OperatorToNormalForm[out, CanonizeKinetic-> True]
	,
		out
	]
];


(* ::Subsubsection::Closed:: *)
(*Extract operator class*)


(* ::Text:: *)
(*Extract all operators of a certain class from an expression*)


SelectOperatorClass::args= "SelectOperatorClass[expr, fields, devNo] takes three arguments.";


SetAttributes[SelectOperatorClass, HoldRest];


Options[SelectOperatorClass]={Evanescent->True};


SelectOperatorClass[expr_ /; MemberQ[expr, _HcTerms],fields_,devs_, opt:OptionsPattern[]] := SelectOperatorClass[ HcExpand @ expr , fields, devs,opt];


SelectOperatorClass[arg_, fields_, devs_:0, OptionsPattern[]]:= Block[{gaugeFields, matterFields, out, opType, conjType, FSCount, FSPower, pw, derivatives, expr=arg, consts, \[Psi]irrelevant, evaExpr },
	(* First, set all fields not appearing in fields_ to zero, for better performance *)
	\[Psi]irrelevant= Alternatives@@ Complement[
		DeleteDuplicates@Cases[expr, (Field[label_,___]|FieldStrength[label_,___]):>label, All],
		DeleteDuplicates[fields/.Bar->Identity]
	];
	
	(* FS also count as derivatives here *)
	If[devs>1,
		expr= expr/.{Field[\[Psi]irrelevant,___]->0}
	,
		expr= expr/.{Field[\[Psi]irrelevant,___]->0, FieldStrength[\[Psi]irrelevant,___]->0}
	];

	(* subtract off constants *)
	If[FreeQ[expr,_Field,All]&&FreeQ[expr,_FieldStrength,All], Return[expr]];
	expr=  TermsToList@ HcExpand@ expr;
	(* TermstoList is required since the Select below returns 1 if the head of expr is Times *)
	(*Extract evanescent operators*)
	evaExpr= Total@ Select[expr, !FreeQ[#,_EvaOp]&];
	expr= Total[expr];
	consts= Last[Last@Reap[
		expr= Replace[
			expr,
			{a_/;(FreeQ[a,_Field,All]&&FreeQ[a,_FieldStrength,All]):>(Sow[a]; 0),_EvaOp->0},
			If[MatchQ[Head@expr, Plus|List], 1, 0]
		]
	],{}];
	
	LagrangianLikeCheck@ expr;
	If[OptionValue@Evanescent,expr=expr+evaExpr];

	opType= Hold@ fields/. Bar-> Conj// ReleaseHold;
	MatterFieldListCheck@ opType;
	DevNoCheck@ devs;

	gaugeFields= List@@ Query[All, Key@ Field]@ $GaugeGroups;
	FSPower=Product[pw[gLabel]^Count[opType, gLabel], {gLabel, gaugeFields}];
	FSCount=2 Sum[Count[opType, gLabel], {gLabel, gaugeFields}];
	opType=DeleteCases[opType, Alternatives@@ gaugeFields];
	derivatives=devs+Total@FSCount;

	opType= {Sort@ opType, derivatives};
	conjType= OpClassConjugate@ opType;
	out= MatchOperatorPatterns@ EnforceOperatorInvariants@ expr;
	out= out/. AtomicOp[class_, _, _]/; !MatchQ[class, opType|conjType]-> 0;
	out= out/. EvaOp[class_,___]/; !MatchQ[class, opType|conjType]->0;
	out= CollectCoefficients@ ExprFlavorCanonize@ out/. AtomicToOpReplacementPattern[];
	out= If[FSPower===1, out, Coefficient[out/.FieldStrength[label_,rest___]:>pw[label] FieldStrength[label,rest], FSPower]/._pw:>0];
	out= (*NormalForm[out, Unique-> False]*)OperatorToNormalForm[out,CanonizeKinetic->True]
]


SelectOperatorClass@ ___:= (Message[SelectOperatorClass::args]; Abort[];);


General::devNo= "`1` is not a valid number of derivatives. A non-negative integer is expected.";
General::fields= "`1` is not a valid list of matter fields. A list of matter field labels (allowing for Bar) is expected.";


MatterFieldListCheck@ fields_:= Block[{matterFields},
	matterFields= Keys@ Query[Select[(MatchQ[#[Type], Scalar|Fermion|Vector]&)]]@ $FieldAssociation;
	If[!MatchQ[fields, f_List /; SubsetQ[matterFields, f/. Conj-> Identity]],
		Message[General::fields, fields];
		Abort[];
	];
]


DevNoCheck@ devs_:=
	If[!MatchQ[devs, _Integer? NonNegative],
		Message[General::devNo, devs];
		Abort[];
	];



(* ::Subsubsection::Closed:: *)
(*GetOperatorCoefficient*)


(* ::Text:: *)
(*Function the extracts the coefficients of an operator matching exactly the form specified*)


GetOperatorCoefficient::unkwnop= "`1` is not regonized as an operator";
GetOperatorCoefficient::nonminop= "The operator `1` is not a minimal operator and its coefficient cannot be extracted.";


GetOperatorCoefficient[lagrangian_, operator_]:= Module[{expr, atomOp, preFact, class, id},
	Operator/@ operator/. op_Operator:> If[!OperatorSatisfiesInvariantsQ@ op, 
			Message[GetOperatorCoefficient::nonminop, op];
			Abort[];
		]; 
	
	(*Cast to atomic operators to compare operators*)
	atomOp= MatchOperatorPatterns@ operator;
	Switch[atomOp
	,_AtomicOp,
		preFact= 1;
	,Times[_? NumberQ, _AtomicOp],
		{preFact, atomOp}= Replace[atomOp, Times[x_, a_AtomicOp]:> {1/x, a}];
	,_,
		Message[GetOperatorCoefficient::unkwnop, operator];
		Abort[];
	];
	{class, id}= Replace[atomOp, AtomicOp[c_, i_, _]:> {c, i}];

	expr= TermsToList@ MatchOperatorPatterns@ EnforceOperatorInvariants@ HcExpand@ lagrangian;

	(*Tests if flavor indices are present:*)
	preFact* Plus@@ If[Length@ Last@ atomOp === 0,
		Cases[expr, Times[couplings__?(FreeQ[Operator]), AtomicOp[class, id, {}]]| AtomicOp[class, id, {}]:>
			Times@ couplings]
	,
		Cases[expr, Times[couplings__?(FreeQ[Operator]), op:(AtomicOp[class, id, currentInds_])]:>
			TransformCouplingFlavors[Times@ couplings, currentInds, Last@ atomOp]]
	]// Simplify
]


(* ::Text:: *)
(*Change the labels of the open flavor indices of a coupling term to some new ones*)


TransformCouplingFlavors[couplings_, currentInds_, newInds_]:= Module[{replacement},
	replacement= Thread@ Rule[currentInds, newInds];
	RelabelIndices[couplings, Unique-> True]/. replacement// RelabelIndices
]


(* ::Subsubsection::Closed:: *)
(*Lookup operator properties*)


(* ::Text:: *)
(*Extracts the operator information *)


LookupOperatorProperties::nonminop= "The operator `1` is not minimal and cannot be mapped.";


LookupOperatorProperties[operator_Operator]:= Module[{op},
	If[!OperatorSatisfiesInvariantsQ@ operator,
		Message[LookupOperatorProperties::nonminop, operator];
		Abort[];
	];
	op= MatchOperatorPatterns@ operator;
	FirstCase[op, AtomicOp[class_, id_, _]:> $operators[class, id], <||>, All]
]


LookupOperatorProperties[operator_AtomicOp]:= Module[{op},
	op= operator;
	FirstCase[{op}, AtomicOp[class_, id_, _]:> $operators[class, id], <||>, All]
]


(* ::Text:: *)
(*Looks op flavor properties of an operator*)


LookupOperatorFlavorProperties[operator:_AtomicOp|_Operator]:= Module[{props, selfConj, selfConjExchange},
	props= LookupOperatorProperties@ operator;

	selfConj= props@ SelfConjugate;
	selfConjExchange= If[selfConj,
			If[Length@ First@ props@ InequivalentPermutations > 0,
				OpConjFlavorPermutation@ FirstCase[Operator@ operator, _Operator|_AtomicOp, Operator[], All]
			,
				True
			]
		,
			False
		];
	<|
		SelfConjugate-> selfConj,
		ConjugateIndexExchange-> selfConjExchange,
		InequivalentPermutations-> props@ InequivalentPermutations,
		Symmetries-> props@ Symmetries
	|>
]


(* ::Text:: *)
(*Special treatment for CompOp*)


LookupOperatorFlavorProperties[operator_CompOp]:= Module[{selfConj,inequivPerms, selfConjExchange},
	selfConj = $compoundOperators[Sequence@@ operator[[;;3]], SelfConjugate];

	(* take inequivalent permutations from one of the AtomicOp that make up this CompOp *)
	inequivPerms = LookupOperatorFlavorProperties[FirstCase[ToAtomicForm@RelabelIndices@AtomicToNormalForm@operator,_AtomicOp,Nothing,All]]@ InequivalentPermutations;

	selfConjExchange= $compoundOperators[Sequence@@ operator[[;;3]], ConjugateIndexExchange];
	(*selfConjExchange= If[TrueQ@ selfConj || selfConj === -1,
			If[Length@ First@ inequivPerms > 0,
				OpConjFlavorPermutation@ operator
			,
				True
			]
		,
			False
		];*)

	<|
		SelfConjugate-> selfConj,
		ConjugateIndexExchange-> selfConjExchange,
		InequivalentPermutations-> inequivPerms,
		Symmetries-> $compoundOperators[Sequence@@ operator[[;;3]], Symmetries]
	|>
]


(*OpConjFlavorPermutation@ op_AtomicOp:= OpConjFlavorPermutation@ FirstCase[Operator@ RelabelIndices@AtomicToNormalForm@ op, _Operator, Operator[], All]*)
OpConjFlavorPermutation@ AtomicOp[class_, id_, _]:= $operators[class, id, ConjugateIndexExchange];


OpConjFlavorPermutation@ op_CompOp:= Module[{bop},
	bop = InternalSimplify[Bar@RelabelIndices@AtomicToNormalForm@op, InternalOpRepresentation->True];
	FindPermutationOrder@@ Cases[{op, bop} , CompOp[__, inds_]-> inds, All]
]


OpConjFlavorPermutation@ op_Operator:= Module[{},
	FindPermutationOrder@@ Cases[MatchOperatorPatterns@ {op, Bar@ op}, AtomicOp[__, inds_]-> inds, All]
]


(* ::Section:: *)
(*IbP Simplify*)


(* ::Text:: *)
(*Exact simplification procedure*)


(* ::Subsection:: *)
(*ConstructOperatorIdentities*)


(* ::Subsubsection::Closed:: *)
(*Main function*)


(* ::Text:: *)
(*Function to make a minimal viable set of identities (in the form of substitution rules) governing an operator class*)
(*reduction= <dDimensional| Evanescent| FourDimensional> determines the type of identities used *)


ConstructOperatorIdentities::basis= "Somethings gone terribly wrong: we could not map the operators into the simplification vector space. ";
ConstructOperatorIdentities::indices= "Somethings gone terribly wrong: the operator contain open gauge or Lorentz indices.";
ConstructOperatorIdentities::flavor= "It looks like different flavor permutations among the identities for operator class `1` are incompatible. The result of the simplification should be considered suspect.";


Options[ConstructOperatorIdentities]={ShowScores->False};


ConstructOperatorIdentities[opClass_, reduction_, OptionsPattern[]]:= Module[
	{opID= 1, diff, dummies, fourDIdentities, identities, indexToPattern, lhs, nonRedIdentities, op, 
		opIdentities, ordering, operators, indPerms, opsToSymmetrize, rhs, evaOperators, evaopt},
	If[!KeyExistsQ[$operators, opClass], Return@ {}; ];

	evaopt= reduction === Evanescent;
	fourDIdentities= (reduction === FourDimensional) || evaopt;

	(*Initial set of operators to generate identities from*)
	dummies= ConstructDummyIndices@ First[$operators@ opClass]@ IndexType;
	indexToPattern= #-> Index[Pattern[Evaluate@ First@ #, Blank[]], Last@ #]&/@ dummies;

	(*Loop to make all identities of the operator class*)

	identities= Flatten@ Reap[While[opID <= Length@ $operators@ opClass,
		op= Replace[AtomicOp[opClass, opID, dummies]/. $operators[opClass, opID, AtomicOpExpansionPattern],
			Times[_, o_Operator]:> o];
		indPerms= $operators[opClass, opID++, InequivalentPermutations];
		(*If[(a=(Select[(FindOpenIndices@op),Index[_,Lorentz], Infinity]))=!={}, Echo@NiceForm@{a,op}; Message[ConstructOperatorIdentities::indices]; (*Abort[];*)];*)
		
		(*Don't produce identities if 4d identities are applied and operator is 4d reducible*)
		If[fourDIdentities && Reducible4dOpQ@ op, Sow@ {}; Continue[];];

		indPerms= Thread[dummies-> dummies[[#]]]&/@ indPerms;

		(*Construct all identities generated from the operator*)
		opIdentities= RelabelIndices/@ Join@@ Through[{
				 IdentitiesIBP
				,IdentitiesCDCommutation
				,IdentitiesJacobi
				,IdentitiesSpinorLineDerivatives (*May cause infinite loop with Schouten*)
				,IdentitiesSpinorLineTransposition
				,IdentitiesDiracCommutation
				,IdentitiesHeavyFieldStrengths
				,IdentitiesChirality
				,IdentitiesSymmetry
				,IdentitiesGroupEpsilons
				,IdentitiesGroupFierz
				,IdentitiesGroupSchouten
				,If[fourDIdentities, IdentitiesFierz[#, Evanescent-> evaopt]&, Nothing]
				,If[fourDIdentities, IdentitiesLorentzSchouten[#, evaopt]&, Nothing]
			}@ op];

		(*Identify operators in the identities with the canonical form ones*)
		opIdentities= MatchOperatorPatternsInIdentities[opClass, opIdentities,
			Complement[{dDimensional, FourDimensional, Evanescent}, {reduction}]];

		(*Use identities with all inequivalent index permutations of the original operator*)
		Sow@ Flatten@ FlavorCanonizeInternalOp[opIdentities/. indPerms];
	]][[2, 1]];
	
	(*Remove trivial identities (not involving any operators)*)
	identities= DeleteDuplicates@ DeleteCases[identities, 0]; 
	(*Determine operators with non-trivial flavor permutations*)
	opsToSymmetrize= FlavorPermutedOperators@ identities;
	
	(*Add compound identities with all inequivalent flavor permutations*)
	ConstructCompoundBasis[reduction, opClass, opsToSymmetrize];
	
	identities= Join[identities, Join@@ KeyValueMap[
		Function[{id, prop}, Table[
			CompOp[reduction, opClass, id, dummies]- (CompOp[reduction, opClass, id, dummies]/. prop@ CompOpReplacement)/.
				Thread@ Rule[dummies, dummies[[perm]]]
			, {perm, prop@ InequivalentPermutations}]
			],
		$compoundOperators[reduction, opClass] ] ];
	
	If[Length@ identities === 0,
		Return@ {};
	];

	(*Determine non-redundant set of operators and indexpermutations*)
	operators= Join[
		Join@@ KeyValueMap[Function[{id, prop},
				Thread@ {AtomicOp[opClass, id, dummies[[#]]]&/@ prop@ InequivalentPermutations, prop@ Score}],
			$operators@ opClass],
		Join@@ KeyValueMap[Function[{id, prop},
				Thread@ {CompOp[reduction, opClass, id, dummies[[#]]]&/@ prop@ InequivalentPermutations, prop@ Score}],
			$compoundOperators[reduction, opClass] ]
		];

	evaOperators= DeleteDuplicates@Cases[identities, _EvaOp,Infinity];
	evaOperators= {#, Replace[#, EvaOp[type_, id_, _]:> $EvanescentTerms[type, id, Score]]}&/@ evaOperators;

	operators= operators~ Join~ evaOperators;
		
	If[OptionValue@ShowScores,
		Echo@NiceForm@(MatrixForm@(operators//AtomicToNormalForm//RelabelIndices));];

	(*Order the operators according to IBPScore*)
	{operators, ordering}= Transpose@ operators;
	ordering= Ordering@ ordering;
	operators= operators[[ordering]];

	(*Canonize flavor indices in operators so that they match the inequivalent permutations*)
	identities= identities/. AtomicOp[opClass, id_, inds_]:>
		FirstCase[List@@@ MapAt[inds[[#]]&, $operators[opClass, id, Symmetries], {All, 1}],
			{i:Alternatives@@ (dummies[[#]]&/@ $operators[opClass, id, InequivalentPermutations]), sign_}:>
			sign* AtomicOp[opClass, id, i] ];

	(*Identify the operators with basis vectors in a vectorspace*)
	identities= (identities/. MapIndexed[(#1-> UnitVector[Length@ operators, First@ #2]&), operators]);
	
	(*Verify that the mapping into vector space was succesfull*)
	If[!MatchQ[identities, {{(_Integer| _Rational| _Complex)..}..}],
		Message[ConstructOperatorIdentities::basis];
		Print[opClass];
		EchoFunction[Iconize]@identities;
	];

	(*Use row reduction to generate a list of substitution rules for redundant operators*)
	identities= DeleteCases[RowReduce@ identities, {0..}];

	(*Recast the identities of the redundant operators in terms of patterns for AtomicOp/CompOp objects*)
	identities= Table[
		lhs= operators[[First@ FirstPosition[identity, 1]]];
		rhs= lhs- identity . operators;
		(lhs/. indexToPattern)-> rhs
	, {identity, identities}];

	(*Keep only substitutions from AtomicOp; CompOp and EvaOp are not present in the expression where the simplifications are applied.*)
	identities= DeleteCases[identities, Rule[_CompOp| _EvaOp, _]];
	(*Return@ identities;*)
	(*Don't include multiple flavor permutations *)
	nonRedIdentities= DeleteDuplicatesBy[identities, #[[1, 2]]&];
	
	(*There has previously been problems where different flavor permutations did not give manifestly identical identities*)
	(*Check that this is actually the case (no abort!)*)
	(*We don't have good control for >2 repeated equivalent flavors (e.g. qqql class), so no check then.*)
	If[!MemberQ[Length/@ $operators[opClass, 1, IndexGrouping], _?(# > 2&)], 
		diff= (identities[[;;, 1]]/. Verbatim[Pattern][x_, Verbatim[Blank[]]]:> x/. nonRedIdentities) - identities[[;;, 2]];
		diff= FlavorCanonizeInternalOp@ diff/. _EvaOp-> 0; (*the symmetries of the evanescent operators are not well-known*)
		If[!VectorQ[diff, # === 0 &],
			Message[ConstructOperatorIdentities::flavor, Format[opClass, NiceForm]];
		];
	];
	
	nonRedIdentities 
];


(* ::Text:: *)
(*Determines the set of all AtomicOps that have different flavor permutations across a "connected" set of operator identities*)


FlavorPermutedOperators@ identities_List:= Module[{ops, tally},
	ops= Cases[#, _AtomicOp, All]&/@ identities;
	ops= JoinOverlappingSets@ ops; 
	tally= Flatten[Tally/@ ops[[;;, ;;, 2]], 1];
	DeleteDuplicates@ Select[tally, #[[2]] > 1&][[;;, 1]]
];


(* ::Text:: *)
(*Determines the connected components of a list of sets, i.e., gather the sets with overlapping components. Works by constructing a graph with the minimal set of edges to connect each set internally *)


JoinOverlappingSets@ sets_List:= Module[{ss, vertices, edges}, 
	ss= DeleteCases[DeleteDuplicates/@ sets, {}];
	vertices= DeleteDuplicates@ Flatten@ ss;
	edges=Flatten[Thread[UndirectedEdge[First@ #, Rest@ #]]&/@ ss];
	ConnectedComponents@ Graph[vertices, edges]
];


(* ::Subsubsection::Closed:: *)
(*Reducible operator in identities*)


(* ::Text:: *)
(*Unless ReductionIdentities-> dDimensional there are operators that should never survive the reduction step with GammaReduction and RemoveLCTensors. As such we don't have to build new identities from these reducible operators. Reducible operators are operators with either *)
(*	i) multiple LC tensors *)
(*	ii) a LC tensor contracted with any gamma matrix *)
(*	iii) four-fermion operators with >2 shared Lorentz indices between the spin chains*)
(*	iv) four-fermion operators with left and right tensor currents*)


Reducible4dOpQ@ op_Operator:= Module[{gammaInds, lcInds, spinChains, spinChainIndices, nproj},
	(*Condition i)*)
	lcInds= Cases[op, _LCTensor, All];
	If[Length@ lcInds > 1, Return@ True; ];

	(*Condition ii)*)
	If[Length@ lcInds === 1,
		lcInds= List@@ First@ lcInds;
		gammaInds= Union@@ Cases[op, g_GammaM:> List@@ g, All];
		If[IntersectingQ[lcInds, gammaInds], Return@ True; ];
	];

	(*Condition iii) and iv)*)
	spinChains= Cases[op, NCM[_, d_DiracProduct, _]-> d, All];
	If[Length@ spinChains < 2, Return@ False; ];
	spinChainIndices= Join@@ Cases[#, g_GammaM|Transp@ g_GammaM :> List@@ g]&/@ spinChains;
	spinChainIndices=Length[Intersection@@ spinChainIndices];
	nproj=Length@DeleteDuplicates[Cases[spinChains, _Proj,All]];
	spinChainIndices>2 || (spinChainIndices===2 && nproj>1)
]


(* ::Subsection:: *)
(*Operator Identities *)


(* ::Text:: *)
(*Functions returning the results of applying various identities to an operator*)


(* ::Text:: *)
(*With every identity that (may introduce) introduces new CGs locally, it  is important to apply a ContractCG to ensure complete contractions of CGs in all operators. Otherwise the assumption of CG reduction as part of identifying identical operators is violated. *)


(* ::Subsubsection::Closed:: *)
(*IdentitiesIBP*)


(* ::Text:: *)
(*Gives identities from IBP relations *)


IdentitiesIBP@ op_Operator:= Module[{positions, \[Mu]ind},
	If[!FreeQ[op, _PartialDev], Return[{}]]; (*Break in case of partial derivatives*)
	(*Positions of all outmost derivatives*)
	positions= Join[#, {-1, 1}]&/@ Position[op, (Field|FieldStrength)[___, {__}], All, Heads-> False];
	Table[
		RefineDiracProducts@ CD[op[[Sequence@@ pos]], Delete[op, pos]]
	, {pos, positions}]
	(*Account for EoMs...*)
	~ Join~
	(RefineDiracProducts@ CD[\[Mu]ind, #]&/@ ReplaceListSubExprs[op, EoM[f_]:> EoMSplitter[\[Mu]ind, f]])
];


(* ::Text:: *)
(*Remove derivative with index \[Mu] from EoM (for IbP)*)


EoMSplitter::unexp= "EoMSplitter applied to unexpected expression `1`.";
EoMSplitter[\[Mu]_, f:Field[_, Scalar, __] ]:= CD[\[Mu], f];
EoMSplitter[\[Mu]_, f:Bar@ Field[_, Scalar, __] ]:= CD[\[Mu], f];
EoMSplitter[\[Mu]_, f:Field[_, Fermion, __] ]:= \[Gamma]@ \[Mu]\[CenterDot] f;
EoMSplitter[\[Mu]_, f:Bar@ Field[_, Fermion, __] ]:= f\[CenterDot] \[Gamma]@ \[Mu];
EoMSplitter[\[Mu]_, f:Transp@ Field[_, Fermion, __] ]:= f\[CenterDot] Transp@ \[Gamma]@ \[Mu];
EoMSplitter[\[Mu]_, f:Transp@ Bar@ Field[_, Fermion, __] ]:= Transp@ \[Gamma]@ \[Mu]\[CenterDot] f;
EoMSplitter[\[Mu]_, Field[f_, Vector[Index[\[Nu]_, type_]], inds_, {}] ]:= FieldStrength[f, {Index[\[Mu], type], Index[\[Nu], type]}, inds, {}];
EoMSplitter[\[Mu]_, Bar@ Field[f_, Vector[Index[\[Nu]_, type_]], inds_, {}] ]:=Bar@ FieldStrength[f, {Index[\[Mu], type], Index[\[Nu], type]}, inds, {}];
EoMSplitter[\[Mu]_, x_]:= (
	Message[EoMSplitter::unexp, x];
	Abort[];
);


(* ::Subsubsection::Closed:: *)
(*IdentitiesCDCommutation*)


(* ::Text:: *)
(*Commutation relations from commuting CD on fields/field strengths*)


IdentitiesCDCommutation@ op_Operator:= Block[{positions},
	(*Position patterns match only on covariant derivatives (e.g., no partials)*)
	positions= Position[op, (Field|FieldStrength)[___, {_Index, __Index}], All, Heads-> False];
	Flatten@ Table[
		pos= Sequence@@ pos;
		Table[
			If[op[[pos, -1, n]] === op[[pos, -1, n+ 1]],
				Nothing
			,
				OperatorToNormalForm[ReplacePart[op, {pos}-> CommuteCDs[op[[pos]], n] ]- op]//
					ContractCGs// ContractDelta// Operator
			]
		, {n, Length@ op[[pos, -1]] -1}]
	, {pos, positions}]
]


(* ::Subsubsection::Closed:: *)
(*IdentitiesJacobi*)


(* ::Text:: *)
(*Jaccobi identities: D_a A_{bc} + D_b A_{ca} + D_c A_{ab} = 0*)


IdentitiesJacobi@ op_Operator:= Block[{positions, indices},
	positions= Position[op, FieldStrength[lab_/; GaugeFieldQ@ lab, ___, {__}], All, Heads-> False];
	Table[
		pos= Sequence@@ pos;
		indices= op[[pos, 2]]~ Join~ op[[pos, -1, {-1}]];
		If[DuplicateFreeQ@ indices,
			indices= Thread@ Rule[indices, #]&/@ NestList[RotateRight, indices, 2];
			ReplacePart[op, {pos}-> CD[op[[pos, -1, ;;-2]],
				Plus@@ (op[[pos, ;;-2]]~ Join~ op[[pos, {-1}, {-1}]]/. indices)] ]
		,
			Nothing
		]
	, {pos, positions}]
]


(* ::Subsubsection::Closed:: *)
(*IdentitiesSpinorLineDerivatives*)


(* ::Text:: *)
(*The application of GAction on the fermion field may create new CGs that may be contractable with other CGs*)


IdentitiesSpinorLineDerivatives@ op_Operator:=
	ReplaceListSubExprs[op, f: Field[_, Fermion, _, {a_Index, a_Index}]:> D2Fermion@ f]//
		OperatorToNormalForm// ContractCGs// ContractDelta// Operator// RefineDiracProducts
(*The RefineDiracProducts ensures that all non-eom gammas are embedded into a single antisymmetric gamma*)


(* ::Text:: *)
(*D2 identity on a fermion (D2 = Dslash Dslash -  1/2 F_ab sigma_ab)*)


D2Fermion[f:Field[lab_, Fermion, inds_, {\[Mu]_Index, \[Mu]_Index}] ]:= Module[{\[Nu]= Index[Unique[lInd], Lorentz]},
	- f +
	\[Gamma]@ \[Nu] \[CenterDot] \[Gamma]@ \[Mu] \[CenterDot] Field[lab, Fermion, inds, {\[Nu], \[Mu]}] +
	I/2 \[Gamma][\[Mu], \[Nu]] \[CenterDot] GAction[{\[Mu], \[Nu]}, Field[lab, Fermion, inds, {}]]
]


(* ::Subsubsection::Closed:: *)
(*IdentitiesSpinorLineTransposition*)


(* ::Text:: *)
(*Identities from all possible transpositions of spinor lines involving identical spinors (but potentially different flavor indices).*)


IdentitiesSpinorLineTransposition@ op_Operator:= Module[{positions},
	positions= Position[op, NCM[Bar@ Field[lab_, __], _, Transp@ Bar@ Field[lab_, __] ] |
		NCM[Transp@ Field[lab_, __], _, Field[lab_, __] ]];
	Table[
		op- MapAt[Transp, op, lines]
	, {lines, Subsets@ positions}]
]


(* ::Subsubsection::Closed:: *)
(*IdentitiesDiracCommutation*)


IdentitiesDiracCommutation@ op_Operator:= Module[{rule},
	rule= ncm: Alternatives[
			NCM[_, dp_DiracProduct, Field[__, {\[Mu]_}]|Transp@ Bar@ Field[__, {\[Mu]_}] ]/; 
				ExpandableDiracProductQ[dp, \[Mu]],
			NCM[Bar@ Field[__, {\[Mu]_}]|Transp@ Field[__, {\[Mu]_}], dp_DiracProduct, _]/; 
				ExpandableDiracProductQ[dp, \[Mu]]
		]:> ncm- ExpandGammas@ ncm;
	RefineDiracProducts@ ReplaceListSubExprs[op, rule]
];


ExpandableDiracProductQ[dp_DiracProduct, index_]:= Module[{},
	If[Count[dp, _GammaM| Transp@ _Gamma] =!= 1, Return@ False; ];
	If[FreeQ[dp, GammaM@ OrderlessPatternSequence[index, __]], Return@ False;];
	True
];


(* ::Text:: *)
(*We might need more identities here *)


(*IdentitiesDiracCommutation@ op_Operator:=
	If[Count[op,_NCM,Infinity]=!=2,
	Join[
		ReplaceListSubExprs[op, d: DiracProduct[___, GammaM[_, __], ___]:>
			ASymGammaExpand@ d] -op,
		ReplaceListSubExprs[op, d: DiracProduct[___, Transp@ GammaM[_, __], ___]:>
			Transp@ ASymGammaExpand@ Transp@ d] -op,
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, GammaM@ \[Mu]_, GammaM@ \[Nu]_, b___]:>
			-DiracProduct[a, GammaM@ \[Nu], GammaM@ \[Mu], b] + 2 Metric[\[Mu], \[Nu]] DiracProduct[a, b]] -op,
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, Transp@ GammaM@ \[Mu]_, Transp@ GammaM@ \[Nu]_, b___]:>
			-DiracProduct[a, Transp@ GammaM@ \[Nu], Transp@ GammaM@ \[Mu], b] + 2 Metric[\[Mu], \[Nu]] DiracProduct[a, b]] -op,
		(*Additional*)
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, GammaM[\[Mu]_, \[Nu]_], b___]:>
			DiracProduct[a, GammaM@ \[Mu], GammaM@ \[Nu], b] - Metric[\[Mu], \[Nu]] DiracProduct[a, b]] -op,
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, GammaM[\[Mu]_, \[Nu]_], b___]:>
			Metric[\[Mu], \[Nu]] DiracProduct[a, b] - DiracProduct[a, GammaM@ \[Nu], GammaM@ \[Mu], b]] -op,
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, Transp@GammaM[\[Mu]_, \[Nu]_], b___]:>
			DiracProduct[a, Transp@GammaM@ \[Nu], Transp@GammaM@ \[Mu], b] - Metric[\[Mu], \[Nu]] DiracProduct[a, b]] -op,
		ContractMetric/@ ReplaceListSubExprs[op,
			DiracProduct[a___, Transp@GammaM[\[Mu]_, \[Nu]_], b___]:>
			Metric[\[Mu], \[Nu]] DiracProduct[a, b] - DiracProduct[a, Transp@GammaM@ \[Mu], Transp@GammaM@ \[Nu], b]] -op
		],
		Nothing
	]*)


(* ::Subsubsection::Closed:: *)
(*IdentitiesChirality*)


(* ::Text:: *)
(*Identities for matching the chirality for VL fermions. They should only be applied when a spinor line involves two vector-like fields. Otherwise the chirality is fixed by whatever chiral fermion is present.  *)


IdentitiesChirality@ op_Operator:= 
	ReplaceListSubExprs[op, {
		NCM[f1: Alternatives[
				Bar@ Field[l1_, Fermion, __], Transp@ Field[l1_, Fermion, __],
				EoM@ Bar@ Field[l1_, Fermion, __], EoM@ Transp@ Field[l1_, Fermion, __] ],
			dp___DiracProduct,
			f2: Alternatives[
				Field[l2_, Fermion, __], Transp@ Bar@ Field[l2_, Fermion, __],
				EoM@ Field[l2_, Fermion, __], EoM@ Transp@ Bar@ Field[l2_, Fermion, __] ] ]/;
			!$FieldAssociation[l1, Chiral] && !$FieldAssociation[l2, Chiral]:>
				f1\[CenterDot] Replace[NCM@ dp, {
					NCM@ DiracProduct[g___, Gamma5]:> DiracProduct@ g \[CenterDot] (\[Gamma]@ 5+ PL- PR),
					NCM@ DiracProduct[g___, Proj@ p_]:> DiracProduct@ g \[CenterDot] (DiracProduct@ Proj@ p - 1/2 -p/2* \[Gamma]@ 5),
					x_:> x \[CenterDot] (1- PL- PR) (*default: no \[Gamma]5 or projector*)
				}]\[CenterDot] f2
	}];


(*IdentitiesChirality@ op_Operator:= 
	ReplaceListSubExprs[op, {
		NCM[f1: Alternatives[
				Bar@ Field[l1_, Fermion, __], Transp@ Field[l1_, Fermion, __],
				EoM@ Bar@ Field[l1_, Fermion, __], EoM@ Transp@ Field[l1_, Fermion, __] ],
			dp___DiracProduct,
			f2: Alternatives[
				Field[l2_, Fermion, __], Transp@ Bar@ Field[l2_, Fermion, __],
				EoM@ Field[l2_, Fermion, __], EoM@ Transp@ Bar@ Field[l2_, Fermion, __] ] ]/;
			!$FieldAssociation[l1, Chiral] && !$FieldAssociation[l2, Chiral]:>
				f1\[CenterDot] (NCM@ dp/. DiracProduct[g___, _Proj|Gamma5]->DiracProduct@ g)\[CenterDot] (1- PL- PR)\[CenterDot] f2,
		NCM[f1: Alternatives[
				Bar@ Field[l1_, Fermion, __], Transp@ Field[l1_, Fermion, __],
				EoM@ Bar@ Field[l1_, Fermion, __], EoM@ Transp@ Field[l1_, Fermion, __] ],
			dp___DiracProduct,
			f2: Alternatives[
				Field[l2_, Fermion, __], Transp@ Bar@ Field[l2_, Fermion, __],
				EoM@ Field[l2_, Fermion, __], EoM@ Transp@ Bar@ Field[l2_, Fermion, __] ] ]/;
			!$FieldAssociation[l1, Chiral] && !$FieldAssociation[l2, Chiral]:>
				f1\[CenterDot] (NCM@ dp/. DiracProduct[g___, _Proj|Gamma5]->DiracProduct@ g)\[CenterDot] (\[Gamma]@ 5+ PL- PR)\[CenterDot] f2
	}];*)


(* ::Subsubsection::Closed:: *)
(*IdentitiesHeavyFieldStrengths*)


(* ::Text:: *)
(*Decomposes heavy field strength tensors in terms of derivatives.*)


IdentitiesHeavyFieldStrengths@ op_Operator:= ReplaceListSubExprs[op,
		fs:FieldStrength[lab_/; !GaugeFieldQ@ lab, {\[Mu]_, \[Nu]_}, inds_, devs_]:> 
			CD[devs, (Field[lab, Vector@ \[Nu], inds, {\[Mu]}]- Field[lab, Vector@ \[Mu], inds, {\[Nu]}])]- fs
	];


(* ::Subsubsection::Closed:: *)
(*IdentitiesSymmetry*)


(* ::Text:: *)
(*Identifies identically vanishing operators from anti-symmetric objects, e.g., the Levi-Civita in H_i eps_ij H_j in the SM.*)


IdentitiesSymmetry@ op_Operator:= Module[{patterns},
	(*Create all patterns from anti-symmetric exchanges of indices for various objects*)
	(*Mark potential places for these anti-symmetric exchanges*)
	patterns= ReplaceListSubExprs[op, {
			cg:CG[lab_, inds_]/; MemberQ[$CGproperties[lab/. Bar-> Identity, Symmetries], -1]:> MapAt[Inactive, cg, 0],
			x:(_LCTensor | _FieldStrength | _GammaM):> MapAt[Inactive, x, 0]
		}];
	If[Length@ patterns === 0, Return@ {};];

	(*Create full operator matching patterns (w/o the symmetry of the marked Inactive head)*)
	patterns= MakeOperatorPattern/@ patterns;

	(*Create full pattern replacement rules *)
	patterns[[;;, 1]]= patterns[[;;, 1]]/. x:Inactive[_][__]:> AntiPattern@ x;
	patterns= MapAt[Activate, RuleDelayed@@@ patterns, {All, 2}];

	Catch[
		Do[
			(*Try all matches of each pattern. Signature 1 indicates incompatible symmetries of the operator*)
			If[MemberQ[ReplaceList[op, pat], 1],
				Throw@ {op}
			];
		, {pat, patterns}];
		{}
	]
];


AntiPattern@ Inactive[CG][lab_, inds_]:= Module[{syms},
	syms= Keys@ Select[# === -1 &]@ $CGproperties[lab/. Bar-> Identity, Symmetries];
	syms= Alternatives@@ (inds[[#]]&)/@ syms;
	CG[lab, syms]
];


AntiPattern@ Inactive[head:LCTensor|GammaM][inds__]:= Module[{sig},
	sig= Signature@ List@ inds;
	With[{syms= Alternatives@@ PatternSequence@@@ Select[Permutations@ List@ inds, sig* Signature@ # === -1 &]},
		HoldPattern@ head[syms]
	]
];


AntiPattern@ Inactive[FieldStrength][lab_, lor_, gr_, devs_]:= With[{lorSiwtch= lor[[{2, 1}]]},
	HoldPattern@ FieldStrength[lab, lorSiwtch, gr, devs] 
];


(* ::Subsubsection::Closed:: *)
(*IdentitiesGroupEpsilons*)


(* ::Text:: *)
(*Identity related to CG coefficients that can simplify even without any contraction *)
(*\[Epsilon][i1, i2,...] \[Epsilon][j1, j2,...] = \[CapitalSigma] \[Sigma](perm) \[Delta][i1, j1] \[Delta][i2, j2] ...*)


IdentitiesGroupEpsilons@ op_Operator:= Module[{cgs, newop= op},
	(*Screen if we should attempt CG replacements*)
	cgs= Select[Cases[newop, CG[name_, __]:> $CGproperties[name/. Bar-> Identity, Group]],
			MatchQ[$Groups@ #, Alg["A", _]]&];
	If[DuplicateFreeQ@ cgs, 
		Return@ {}
	];
	
	(*Perform epsilon replacements*)
	ReplaceList[ReplaceCGs@ newop,
		o: HoldPattern@ Operator[CG[ep_eps, inds1_], CG[Bar@ ep_eps, inds2_], rest__]:>
		o- Signature@ inds2 Plus@@ (Signature@ # Times@@ Thread@ Delta[inds1, #]&/@ Permutations@ inds2) *
			Times@ rest
	]// ContractGroupIndices
];


(* ::Subsubsection::Closed:: *)
(*IdentitiesFierz*)


(* ::Text:: *)
(*Identity related to Fierzing of two fermions bilinears*)


IdentitiesFierz[ op_Operator, opt:OptionsPattern[]]:=
	Join[
		ReplaceList[op, 
			o: HoldPattern@Operator[SP1_NCM,SP2_NCM, rest___/;
				FreeQ[{rest},_NCM]]/;
				FbasisQ[SP1,SP2]/;
				OrderedQ@ {SP1, SP2}:>
			o- Operator[Fierz[OperatorToNormalForm@o, Order->{1,4,3,2},opt]]
		],
		ReplaceList[op, 
			o: HoldPattern@Operator[SP1_NCM,SP2_NCM, rest___/;
				FreeQ[{rest},_NCM]]/;
				FbasisQ[SP1,SP2]/;
				OrderedQ@ {SP1, SP2}:>
			o- Operator[Fierz[OperatorToNormalForm@o, Order->{1,3,4,2},opt]]
		]
	];


(* ::Subsubsection::Closed:: *)
(*IdentitiesLorentzSchouten*)


(* ::Text:: *)
(*The identity constructed by the Schouten identity will be the four-dimensional identity to prevent errant factors of \[ScriptD] in the row reduction step. *)


IdentitiesLorentzSchouten[op_Operator, evanecentOp_]:= Module[
		{ind, indices, positions, rel, relations, subset, dim=If[FreeQ[op,Spatial], $Dimensions, $Dimensions-1]},
	indices= DeleteDuplicates@ Cases[op, Index[_, Lorentz|Spatial], All];
	If[Length@ indices < dim+1, Return@ {};];

	(*Anti-symmetrize accross all subsets dim+1 unique Lorentz indices*)
	relations= Table[
		positions= Table[FirstPosition[op, ind], {ind, subset}];
		(*Ensure that the original operator enters with '+' (for defining evanescent ops) *)
		Signature@ subset* Sum[
				Signature@ perm * ReplacePart[op, Thread@ Rule[positions, perm]]
			, {perm, Permutations@ subset}]
	, {subset, Subsets[indices, {dim+1}] }];
	
	(*Take \[ScriptD]=$Dimensions limit for the relation*)
	relations= RelabelIndices/@ ContractMetric/@ relations/. \[ScriptD]-> $Dimensions;
	
	(*Add evanescent operators*)
	If[evanecentOp,
		Table[
			rel- DefineEvanescentOperator[op, op- rel, Schouten]
		, {rel, relations}]
	,
		relations
	]
];




(* ::Subsubsection::Closed:: *)
(*IdentitiesGroupFierz*)


(* ::Text:: *)
(*Identities shuffling around singular fundamental SU(2) generators*)
(*	0= \[Delta]^i_j t^{ak}_l + \[Delta]^k_l t^{ai}_j - \[Delta]^k_j t^{ai}_l - \[Delta]^i_l t^{ak}_j*)
(*or alternatively*)
(*	0 = f_i g^i t^{ak}_l + f_i g^j t^{ai}_j \[Delta]^k_l - f_i g^k t^{ai}_l - f_l g^i t^{ak}_i*)
(*TODO:*)
(*Expand this identity to SU(N) by relying on *)
(*	\[Delta]^i_j t^{ak}_l = N {t_a, t_b}^i_j t^{bk}_l - N d_{abc} t^{ci}_j t^{bk}_l*)
(*and the ordinary SU(N) Fierz identity to simplify the generator contractions where possible*)


IdentitiesGroupFierz@ operator_Operator:= Module[
		{fundInds, otherInds, genPos, op, out, su2Groups,
			iInd, jInd, kInd, lInd, iBarPos, iPos, kPos, lPos},
	su2Groups= Keys@ Select[$Groups, # === Alg["A", 1]&];
	(*Ensure that indices on fields are Bared or not locally*)
	op= operator/. Bar@ (f:Field|FieldStrength)[lab_, t_, inds_, devs_]:>
		f[Conj@ lab, t, Bar/@ inds, devs]// ReplaceCGs;

	out= Flatten@ Table[Catch[
		(*Determine positions of prospective SU(2) generators*)
		genPos= Position[op, CG[gen@ gr@ fund, _]];
		If[Length@ genPos === 0, Throw@ Nothing;];

		(*Determine all fundamental SU(2) indices of same group*)
		fundInds= DeleteDuplicates@ Cases[op, Index[_, gr@ fund], All];
		If[Length@ fundInds < 3, Throw@ Nothing;];

		jInd= Index[Unique[j], gr@ fund];

		(*Loop over generators and 3rd prospective index*)
		Table[
			{kInd, lInd}= op[[Sequence@@ gPos]][[2, 2;;]];
			{kPos, lPos}= {Join[gPos, {2, 2}], Join[gPos, {2, 3}]};

			otherInds= Complement[fundInds, {kInd, lInd}/.Bar->Identity];
			Table[
				iBarPos= FirstPosition[op, Bar@ iInd];
				iPos= First@ DeleteCases[Position[op, iInd], Append[iBarPos, 1]];
				(*Identity*)
				op- ReplacePart[op, {kPos-> iInd, iPos-> kInd}]- ReplacePart[op, {iBarPos-> lInd, lPos-> Bar@ iInd}] +
					ContractDelta[Delta[kInd, lInd]ReplacePart[op, {kPos-> iInd, lPos-> Bar@ jInd, iPos-> jInd}]]
			, {iInd, otherInds}]
		, {gPos, genPos}]
	], {gr, su2Groups}];

	(*Restore ordinary Bar notation for Field and FieldStrength objects in output identities*)
	out= out/. (f:Field|FieldStrength)[Conj@ lab_, t_, inds_, devs_]:>
		Bar@ f[lab, t, Bar/@ inds, devs];

	ContractGroupIndices@ out 
]


(* ::Subsubsection::Closed:: *)
(*IdentitiesGroupSchouten*)


(* ::Text:: *)
(*Applies Schouten (anti-symmetrization) identity to group indices. For instance \[CurlyEpsilon]_{[ab} \[CurlyEpsilon]_{c]d} = 0 for fundamental SU(2) indices.*)


IdentitiesGroupSchouten@ operator_Operator:= Module[
		{indPos, indSet, indsToSymmetrize, op, opInds, out, rep, reps},
	(*For now only look for indices in small representations*)
	reps= Keys@ Select[$Representations, #[Dimension] < 5&];
	
	(*Ensure that indices on fields are Bared or not locally*)
	op= operator/. Bar@ (f:Field|FieldStrength)[lab_, t_, inds_, devs_]:>
		f[Conj@ lab, t, Bar/@ inds, devs];
	
	(*Group indices by representation*)
	opInds= GroupBy[DeleteDuplicates@ Cases[op, Index[_, Alternatives@@ reps], All], Last];
	
	(*Apply all possible symmetrizations*)
	out= Flatten@ Table[
			indsToSymmetrize= If[$Representations[rep, Reality] =!= 1, Bar/@ opInds@ rep, opInds@ rep];
			Table[
				indPos= Table[FirstPosition[op, ind, All], {ind, indSet}];
				Sum[
					Signature@ perm* ReplacePart[op, Thread@ Rule[indPos, perm]]
				, {perm, Permutations@ indSet}]
			, {indSet, Subsets[indsToSymmetrize, {$Representations[rep, Dimension] + 1}]}]
		, {rep, Keys@ opInds}];	
	
	(*Restore ordinary Bar notation for Field and FieldStrength objects in output identities*)
	out= out/. (f:Field|FieldStrength)[Conj@ lab_, t_, inds_, devs_]:>
		Bar@ f[lab, t, Bar/@ inds, devs];
	ContractGroupIndices@ out
];


(* ::Subsection:: *)
(*Simplify  using  rules*)


(* ::Text:: *)
(*Apply the simplification identities to a Lagrangian expression*)


(* ::Subsubsection::Closed:: *)
(*Construct Hermitian identities*)


(* ::Text:: *)
(*To ensure the same ordering between Hermitian operator types, the corresponding identities are simultaneously constructed*)
(*reduction= <dDimensional| Evanescent| FourDimensional> determines the type of identities used *)


ConstructHermitianSimplificationIdentities[opClass_, reduction_]:=
		ConstructHermitianSimplificationIdentities[opClass, reduction]= Module[{identities},
	identities= ConstructOperatorIdentities[opClass, reduction];

	If[SelfConjugateClassQ@ opClass, Return@ identities];

	(*For complex-type classes the conjugate identities are constructed*)
	ConstructConjugateCompounds[reduction, opClass];

	Join[identities,
		identities/.
			atom_AtomicOp:> ConjugateAtomicOp@ atom/.
			eva_EvaOp:> Bar@ eva/.
			c_Complex:> Conjugate@ c/.
			CompOp[reduction, opClass, id_, inds_]:> CompOp[reduction, OpClassConjugate@ opClass, id,
				inds[[$compoundOperators[reduction, opClass, id, ConjugateIndexPermutation]]]]/.
			HoldPattern[Times[-1, x_]-> y_] :> Rule[x, -y] (*account for possibility of -1 from conjugation*)
	]
]


ConjugateAtomicOp@ atom:AtomicOp[class_, id_, inds_]:=
	Bar@ Replace[atom, $operators[class, id, AtomicOpExpansionPattern]]/.
		op_Operator:> op/. $operators[OpClassConjugate@ class, id, OperatorMatchingPattern]


(* ::Subsubsection::Closed:: *)
(*Unset the identities*)


(* ::Text:: *)
(*Remove all memoized values of ConstructHermitianSimplificationIdentities associated to the operator class and its conjugate*)


UnsetHermitianIdentities[opClass_, reduction_] := Block[{},
	With[{class= opClass},
		RemoveAssociatedDownValues[ConstructHermitianSimplificationIdentities[class, reduction]];
	];
	If[!SelfConjugateClassQ@ opClass,
		With[{class= OpClassConjugate@ opClass},
			RemoveAssociatedDownValues[ConstructHermitianSimplificationIdentities[class, reduction]];
		];
	];
]


(* ::Subsubsection::Closed:: *)
(*IBP Simplification function*)


(* ::Text:: *)
(*Applies the operator identities to a given pattern *)


Options@ IBPSimplify= {
		ReductionIdentities-> dDimensional,  
		Verbose->True
	};


IBPSimplify[expr_, OptionsPattern[]]:= Module[
		{class, evTerms, out, physTerms, subs, terms, treeLag, evaopt, redID, n= 0},
	{redID, evaopt}= Switch[OptionValue@ ReductionIdentities
		,dDimensional,
			{dDimensional, False}
		,Evanescent|EvanescenceFree,
			{Evanescent, True}
		,FourDimensional,
			{FourDimensional, False}
	];
	out= If[redID === dDimensional,
			expr
		,
			(*Reduce Lorentz/Dirac structures that should always be reduced, i.e.,
			structures that are never part of the basis.*)
			(*It seems more reliable to first remove LC tensors, but maybe we need to consider the ordering,
			or even unifying the two.*)
			(*Consider collecting the operators first (Nb. w/o going to AtomicOps!)*)
			GammaReductionInternal[RemoveLCTensorInternal[OperatorFlavorSeparate@ Operator@ expr,
				evaopt], evaopt, redID === FourDimensional]
		];
		
	out= MatchOperatorPatterns@ out;

	(*Determine identities for each group of operator types*)
	(*Collect terms by Hermitian operator class*)
	terms= Normal@ GroupBy[TermsToList@ out, FirstCase[#, AtomicOp[class_, __] :> class, {}, All]& ];
	terms= Merge[MapAt[If[# === {}, {}, First@ Sort@ {#, OpClassConjugate@#}] &,
		terms, {All, 1}], Apply@ Join];
	(*Perform the replacements by operator class for performance*)
	OptionalMonitor[OptionValue@ Verbose,
		out= Sum[Catch[
				n++;
				If[class === {}, (*Constants are kept unchanged*)
					Throw[Plus@@ terms@ class];
				];
				subs= ConstructHermitianSimplificationIdentities[class, redID];
				Plus@@ terms@ class/. subs// ExprFlavorCanonize
				
			], {class, Keys@ terms}];
	, 
		Row[{"Simplifying operators in class ", Format[class, NiceForm], 
			"\t("<> ToString@ n<> "/"<> ToString@ Length@ terms<> ")"}, " "] 
	];
		
	If[redID === Evanescent,
		out= out/. {hbar* _EvaOp-> 0, hbar* Bar@ _EvaOp-> 0};
	];

	(*Calculate the finite shift from the evanescent operators*)
	If[OptionValue@ ReductionIdentities === EvanescenceFree && !FreeQ[out, _EvaOp],
		out= Catch[
			(*Separate physical and evanescent terms*)
			{physTerms, evTerms}= SelectAndDeleteCases[TermsToList@ out, _? (FreeQ[#, EvaOp]&)];
			physTerms= Plus@@ physTerms;
			evTerms= Plus@@ evTerms/. evaop_EvaOp :> CollectOperators[ev*  RefineDiracProducts@ ExpandEvanescentOperators@ evaop, Simplify->False];
			If[evTerms === 0, (*this may occur e.g. for Cee where the coefficient symmetries kill evanescent contribution*)
				Throw[physTerms];
			];

			(*Pass tree-level lag + evanescent terms to evanescent loops*)
			treeLag= AtomicToNormalForm@ CollectCoefficients[physTerms/. hbar->0];
			evTerms= OptionalMonitor[OptionValue@ Verbose,
				EvanescentContribution[evTerms+ treeLag, Verbose-> OptionValue@ Verbose]
			, "Computing evanescent shifts\[Ellipsis]"];

			physTerms+ IBPSimplify[evTerms, ReductionIdentities-> Evanescent]
		];
	];

	(*Expands the dimensional dependence from the Dirac/Lorentz identities*)
	If[OptionValue@ ReductionIdentities =!= FourDimensional && !FreeQ[out, \[ScriptD]],
		out= EpsExpand@ out;
	];

	out
]


(* ::Subsection:: *)
(*Output functions*)


(* ::Subsubsection::Closed:: *)
(*Scoped simplification functions*)


(* ::Text:: *)
(*Simplification function outputting the Lagrangian either in Operator form or in the internal representation (e.g. AtomicOp) *)


Options@ InternalSimplify= {
		CoefficientSimplify-> True,
		ContractOverallDelta-> True,
		InternalOpRepresentation-> False, 
		ReductionIdentities-> dDimensional,
		Verbose-> True,
		VLFBasis-> None,
		TrigSimplify-> False
	};


InternalSimplify[expr_, opts:OptionsPattern[]]:= Module[{out},
	(*Update chiral basis - internally - only by explicit choice*)
	If[OptionValue@ VLFBasis =!= None, 
		UpdateChiralBasis@ OptionValue@ VLFBasis;
	];
	
	out= CollectCoefficients[IBPSimplify[EnforceOperatorInvariants@ expr,
			ReductionIdentities-> OptionValue@ ReductionIdentities, Verbose-> OptionValue@ Verbose],
		FilterOptions[opts, CollectCoefficients] ];

	If[OptionValue@ InternalOpRepresentation,
		out
	,
		AtomicToOperatorForm[out]
	]
];


(* ::Subsubsection::Closed:: *)
(*Main Exported simplification function*)


(* ::Text:: *)
(*Output to NormalForm for the operators *)


Options[GreensSimplify]= {
		PreferInputOperators-> False,
		ReductionIdentities-> dDimensional,
		Verbose-> True,
		VLFBasis-> Automatic
	};


GreensSimplify[arg_,OptionsPattern[]]? OptionsCheck:= Block[
		{consts, expr=LagrangianExpand@HcExpand[arg], redID= OptionValue@ ReductionIdentities},
	UpdateChiralBasis@ OptionValue@ VLFBasis;
	
	If[redID === EvanescenceFree,
		AddToBibliography["EvanescentTreatment", "Simplified expression to evanescent-free scheme (with GreensSimplify)"];
		If[!FreeQ[expr, _PartialDev], Message[GreensSimplify::partialdev]; Abort[]]
	];

	(* subtract off constants *)
	{consts, expr}= SeparateOutConstants@ expr;

	(* simplifications *)
	LagrangianLikeCheck@ expr;
	
	expr= AtomicToNormalForm@ CollectCoefficients@ IBPSimplify[EnforceOperatorInvariants@ expr, 
		ReductionIdentities-> redID, Verbose-> OptionValue@ Verbose];

	(* add back the constants *)
	expr + consts
]


(* ::Subsubsection::Closed:: *)
(*Return simplification identities *)


(* ::Text:: *)
(*Returns simplifications rules determined by ConstructOperatorIdentities to reduce operators from an operator class to the basis*)


Options[IBPIdentities]={ReductionIdentities->dDimensional, ShowScores->False};


IBPIdentities::args= "IBPIdentities[fields, devNo] takes two arguments.";


SetAttributes[IBPIdentities, HoldFirst];
IBPIdentities[fields_List, devs_Integer, OptionsPattern[]](*? OptionsCheck*):= Module[{opClass, identities},
	opClass= Hold@ fields/. Bar-> Conj// ReleaseHold;
	MatterFieldListCheck@ opClass;
	DevNoCheck@ devs;

	opClass= {opClass// Sort, devs};
	If[!MemberQ[Keys@ $operators, opClass],
		Return@ {};
	];

	(*Remove compound definitions:*)
	identities= DeleteCases[ConstructOperatorIdentities[ opClass, OptionValue@ReductionIdentities, ShowScores->OptionValue@ShowScores], Rule[_CompOp, _]];

	identities/. Verbatim[Pattern][x_, Blank[]]-> x// AtomicToNormalForm
]


IBPIdentities@ ___:= (Message[IBPIdentities::args]; Abort[];);


(* ::Subsection:: *)
(*HcSimplification*)


(* ::Subsubsection::Closed:: *)
(*HcTerms and HcExpand*)


(* ::Text:: *)
(*Collect Hermitian conjugate terms in a Lagrangian*)


HcTerms@ 0= 0;
Bar@HcTerms[arg___]:=HcTerms[arg]


HcExpand[Lag_]:= Lag /.{HcTerms[expr_]:> expr + Bar@expr} ;


(* ::Subsubsection::Closed:: *)
(*HcSimplify*)


HcSimplify::nothermitian = "The Lagrangian is not hermitian. The non-Hermitian part \[ScriptCapitalL]-\[ScriptCapitalL]\[ConjugateTranspose] is given by: `1`"


(* ::Text:: *)
(*This implementation of HcSimplify proceeds under the assumption that the Lagrangian is explicitly Hermitian after application of InternalSimplify*)


Options[HcSimplify]= {
		Simplify-> True,
		VLFBasis-> Automatic
	};


HcSimplify[arg_, opts:OptionsPattern[]]:= 
	HcSimplifyInternal[arg, InternalOpRepresentation -> False, Simplify-> OptionValue@ Simplify, VLFBasis-> OptionValue@ VLFBasis];


Options[HcSimplifyInternal]= {
		ContractOverallDelta-> True,
		InternalOpRepresentation-> False(*True*), 
		ReductionIdentities-> dDimensional,
		Simplify-> False,
		VLFBasis-> Automatic
	};


HcSimplifyInternal[arg_,OptionsPattern[]]:= Module[{complex, out, real, expr=arg, hermite, hermiteManifest},
	(* check hermiticity and manifest hermiticity *)
	{hermite, hermiteManifest} = HermitianQ[expr,"manifest"->True];
	If[!hermite,
		Message[HcSimplify::nothermitian, Format[GreensSimplify[expr- Bar@ expr],NiceForm]];
		Abort[];
	];
	(* if the Lagrangian is Hermitian, but not manifestly Hermitian, manifestly Hermitianize it *)
	If[!hermiteManifest,
		expr = (expr+Bar[expr])/2;
	];

	out= InternalSimplify[EnforceOperatorInvariants@ expr, ContractOverallDelta->OptionValue@ ContractOverallDelta, InternalOpRepresentation-> True, VLFBasis-> OptionValue@ VLFBasis, ReductionIdentities -> OptionValue@ ReductionIdentities];
	If[out === 0, Return@ 0;];
	out= If[Head@ out =!= Plus, {out}, List@@ out];

	real= Cases[out, (Times[op:(_AtomicOp|_CompOp), rest__]|op:(_AtomicOp|_CompOp) )/; RealOpQ@ op:>
		{op, Times@ rest}];
	complex= Cases[out, (Times[op:(_AtomicOp|_CompOp), rest__]|op:(_AtomicOp|_CompOp))/; !RealOpQ@ op:>
		{op, Times@ rest}];
	
	complex= DeleteDuplicatesBy[SortBy[complex,Count[#,_Bar,Infinity]&], (Sort@ {#, MapAt[OpClassConjugate, #, If[Head@# === AtomicOp,1,2]]}&@* First)];

	If[TrueQ@ OptionValue[InternalOpRepresentation],
		out= Plus@@ Times@@@ real + HcTerms[Plus@@ Times@@@ complex]
		,
		out= CollectOperators[AtomicToNormalForm[Plus@@ Times@@@ real],ContractOverallDelta->OptionValue@ ContractOverallDelta, Simplify->OptionValue@Simplify] + HcTerms[CollectOperators[AtomicToNormalForm[Plus@@ Times@@@ complex],ContractOverallDelta->OptionValue@ ContractOverallDelta,Simplify->OptionValue@Simplify]] 
	];

	out
]


(* ::Text:: *)
(*Test if an AtomicOp is real or complex (exception if the operator is a kinetic term)*)


RealOpQ@ op:AtomicOp[opClass_, id_, _]:= $operators[opClass, id, SelfConjugate] ||
	CanonicalKineticOpQ@ AtomicToOperatorForm[op];
RealOpQ@ CompOp[red_, opClass_, id_, _]:=
	TrueQ@ $compoundOperators[red, opClass, id, SelfConjugate];


(* ::Subsection:: *)
(*Manipulate operator basis (scoring)*)


(* ::Subsubsection::Closed:: *)
(*PromoteOperatorsToBasis*)


(* ::Text:: *)
(*Set operators as preferred*)


PromoteOperatorsToBasis::nonminop= "The operator `1` is not minimal and cannot be used in the basis."


PromoteOperatorsToBasis[ops_]:= Module[{atomics, class, classes, operators},
	operators= If[Head@ ops === List, ops, TermsToList@ ops];
	Operator/@ operators/. op_Operator:> If[!OperatorSatisfiesInvariantsQ@ op, 
			Message[PromoteOperatorsToBasis::nonminop, OperatorToNormalForm@ op];
			Abort[];
		]; 

	atomics= DeleteDuplicates@ Cases[MatchOperatorPatterns@ ops, _AtomicOp, All];
	classes= DeleteDuplicates@ Replace[atomics, AtomicOp[class_, __]-> class, {1}];
	Do[
		UnsetHermitianIdentities[class, _];
	, {class, classes}];
	PromoteOperatorScore/@ atomics;
];


(* ::Text:: *)
(*Sets preferential scoring to an atomic operator*)


PromoteOperatorScore@ atom:AtomicOp[class_, id_, _]:= Module[{newScore},
	(*If not self-conjugate sets score for the conjugate operator*)
	If[!$operators[class, id, SelfConjugate],
		PromoteIndividualOperatorScore@ ConjugateAtomicOp@ atom;
	];
	PromoteIndividualOperatorScore@ atom;
];
PromoteIndividualOperatorScore@ atom:AtomicOp[class_, id_, _]:= Module[{newScore},
	(*Recalculate the OpScore to prevent accidental repeated preferences*)
	newScore= 30+ OpScore[FirstCase[{atom/. $operators[class, id, AtomicOpExpansionPattern]}, _Operator, 0, All],
		$operators[class, id, SelfConjugate]];
	$operators[class, id, Score]= newScore;
];


(* ::Section:: *)
(*Evanescent terms*)


(* ::Subsection:: *)
(*Evanescent operators*)


(* ::Subsubsection::Closed:: *)
(*Handle functions*)


(* ::Text:: *)
(*$EvanescentTerms stores the evanescent operators by class and the operator ID:*)
(*	<|class1-> <|id1-> properties, ...|> , ...|>*)


ResetEvanescentOperators[]:= Block[{},
	$EvanescentTerms= <||>;
];
ResetEvanescentOperators[];


(* ::Text:: *)
(*Expand the evanescent operators in terms of AtomicOps *)


ExpandEvanescentOperators[expr_]:=
	expr/. EvaOp[class_, id_, inds_]:> AtomicToNormalForm@
		(($EvanescentTerms[class, id, RedundantOperator]- $EvanescentTerms[class, id][BasisOperator])/.
		AssociationThread[FindOpenIndices[$EvanescentTerms[class, id, RedundantOperator]], inds]);


(* ::Text:: *)
(*Looking up an evanescent operator in the global association*)


GetEvanescentTerm::args= "GetEvanescentTerm takes either no argument, one argument: the class of operators ({fields, nodevs}), two argument: the class and the id of the operator, or an EvaOp object direcly.";


SetAttributes[GetEvanescentTerm, HoldAll];


GetEvanescentTerm[]:= $EvanescentTerms //AtomicToNormalForm;
GetEvanescentTerm[op_EvaOp]:= GetEvanescentTerm@@ op[[;;2]] ;
GetEvanescentTerm[class_, id_]:= GetEvanescentTerm[class][id];


GetEvanescentTerm[{fields_,devs_}]:= Block[{gaugeFields, opType, conjType, FSCount, FSPower, pw, derivatives },

	opType= Hold@ fields/. Bar-> Conj// ReleaseHold;
	MatterFieldListCheck@ opType;
	DevNoCheck@ devs;

	gaugeFields= List@@ Query[All, Key@ Field]@ $GaugeGroups;
	FSPower=Product[pw[gLabel]^Count[opType, gLabel], {gLabel, gaugeFields}];
	FSCount=2 Sum[Count[opType, gLabel], {gLabel, gaugeFields}];
	opType=DeleteCases[opType, Alternatives@@ gaugeFields];
	derivatives=devs+Total@FSCount;

	opType= {Sort@ opType, derivatives};
	
	$EvanescentTerms[opType] //AtomicToNormalForm
]


GetEvanescentTerm@ ___:= (Message[GetEvanescentTerm::args]; Abort[];);


(* ::Text:: *)
(*Only consider tree-level evanescent operators (otherwise it's a two loop contribution)*)


(*EvaOp/:(hbar*EvaOp[___]):=0;
EvaOp/:(hbar*Bar @ EvaOp[___]):=0;*)


(* ::Subsubsection::Closed:: *)
(*Lookup an evanescent operator by operator difference *)


FindEvanescentOperator::suscoef= "No scaling of found for evanescent operator recognition (operator class `1`). Error suspected:"


FindEvanescentOperator[opDifference_]:= Module[
		{class, coef, inputIDs, inputInds, inputOp= opDifference, evaIDs, evaInds, evaOp},
	class= FirstCase[inputOp, AtomicOp[c_, __]:> c, {}, All];
	If[MatchQ[$EvanescentTerms@ class, _Missing], (*If no evs of this class have been created*)
		Return@ None;
	];
	
	(*All ops showing up in input difference*)
	inputIDs= SortBy[Cases[inputOp, AtomicOp[class, ID_, inds_]-> {ID, inds}, All], First];
	inputInds= inputIDs[[1, 2]];
	inputIDs= inputIDs[[;;, 1]];
	
	Catch[
		Do[
			(*All ops showing up in current evanescent operator*)
			evaOp= $EvanescentTerms[class, id, RedundantOperator]- $EvanescentTerms[class, id, BasisOperator];
			evaIDs= SortBy[Cases[evaOp, AtomicOp[class, ID_, inds_]-> {ID, inds}, All], First];
			evaInds= evaIDs[[1, 2]];
			evaIDs= evaIDs[[;;, 1]];
			(*Check if all atomics coincide*)
			If[!MatchQ[Cases[evaOp, AtomicOp[class, ID_, _]-> ID, All], {OrderlessPatternSequence@@ inputIDs}],
				Continue[];
			];

			(*Checks if input is number times the evanescent op*)
			coef= inputOp/ (evaOp/. Thread[evaInds-> inputInds])// FlavorCanonizeInternalOp// Simplify;
			If[FreeQ[coef, _AtomicOp],
				(*If so, return the input mapped to that operator*)
				evaInds= inputInds[[FindPermutationOrder[evaInds,
					FindOpenIndices@ $EvanescentTerms[class, id, RedundantOperator]]]];
				Throw[coef* EvaOp[class, id, evaInds]]
(*			,
				Message[FindEvanescentOperator::suscoef, class];
				Echo@ Iconize@ coef;*)
			];
		, {id, Length@ $EvanescentTerms@ class}];
		None
	]
];


(* ::Subsubsection::Closed:: *)
(*Define evanescent operators*)


(* ::Text:: *)
(*Define a global evanescent operator as the difference between initial and final operator, unless the operator already exists (modulo scaling) *)


DefineEvanescentOperator[iniOperator_, finOperator_, origin_]:= Module[
		{basisOp, class, evaOperator, id, redundantOp, result, score},
	(*Map to atomic operator space*)
	{redundantOp, basisOp}= MatchOperatorPatterns/@ Contract/@ {iniOperator, finOperator};
	evaOperator= Collect[redundantOp- basisOp, _AtomicOp];
	(*Return immediately if the purported evanescent operator is trivially zero*)
	If[evaOperator === 0, Return@ 0; ];
	 
	(*Check if evaOperator already exists*)
	result= FindEvanescentOperator[evaOperator];
	If[result =!= None, Return@ result; ];

	class= FirstCase[redundantOp, AtomicOp[c_, __]:> c, {}, All];

	(*Creates a sub-association in $EvanescentTerms if the class is not already pressent.*)
	If[MatchQ[$EvanescentTerms@ class, _Missing],
		$EvanescentTerms@ class= <||>;
	];
	id= Length@ $EvanescentTerms@ class+ 1;

	{basisOp, redundantOp}= InternalCollectOperators[#, InternalOpRepresentation-> True]&/@
		{basisOp, redundantOp};

	(*Give simplification preference to evanescent operators:*)
	score= 100
			(*Take the sum of the minimums of the redundant operator and the basis operator*)
			+ Total[Min/@ (Cases[#, AtomicOp[class, ID_, _]:> $operators[class, ID, Score], All]&/@{redundantOp, basisOp}/.{}->0)]
			(*favor if the basis operator is 0*)
			+ If[basisOp===0, 50,0];

	$EvanescentTerms[class, id]= <|
			RedundantOperator-> redundantOp,
			BasisOperator-> basisOp,
			Origin-> origin,
			Score-> score
		|>;
	EvaOp[class, id, FindOpenIndices[redundantOp]]
];


(* ::Subsection:: *)
(*Levi-Civita relations*)


(* ::Subsubsection::Closed:: *)
(*Internal RemoveLCTensor*)


(* ::Text:: *)
(*This function assumes operator form of input*)


RemoveLCTensorInternal[expr_, produceEva_]:= expr/. a_Operator:> RemoveLCTensorInternal[a, produceEva];


RemoveLCTensorInternal[op_Operator, _]/; FreeQ[op, _LCTensor]:= op;


RemoveLCTensorInternal[op_Operator, produceEva_? BooleanQ]:= Module[
		{opexp= ReleaseOperators@ op, evaOperator, result},
	result= ContractMetric@ LC2Delta@ LC2Gamma5@ opexp;
	evaOperator=If[produceEva,
			DefineEvanescentOperator[opexp, result, LCTensor]
		,
			0
		];
	evaOperator+ result// RelabelIndices// Operator
]


(* ::Subsubsection::Closed:: *)
(*Public RemoveLCTensor*)


(* ::Text:: *)
(*RemoveLCTensor uses four dimensional relations to simplify LC tensors from an expression*)
(*RemoveLCTensor assumes NormalForm input -> Returns NormalForm output*)


Options@ RemoveLCTensor= {Evanescent-> True};


RemoveLCTensor[expr_, OptionsPattern[]]? OptionsCheck:= 
	RemoveLCTensorInternal[OperatorFlavorSeparate@ Operator@ expr, OptionValue@ Evanescent]// 
	OperatorToNormalForm// ContractDelta// RelabelIndices;
