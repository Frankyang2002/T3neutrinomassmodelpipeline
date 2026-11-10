(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Lagrangians`*)


(* ::Subtitle:: *)
(* Checks on validity of Lagrangians and methods to extract information from Lagrangian.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["CheckLagrangian"]


PackageExport["Hermiticity"]
PackageExport["ContractedIndices"]
PackageExport["ClosedSpinChains"]
PackageExport["CanonicallyNormalized"]
PackageExport["HeavyMassBasis"]
PackageExport["FreeOfHeavyTadpoles"]
PackageExport["ChargeNeutral"]
PackageExport["GaugeInvariant"]
PackageExport["UndefinedObject"]
PackageExport["GaugeAnomalies"]
PackageExport["DetailedOutput"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["LagrangianLikeCheck"]


PackageScope["HermitianQ"]
PackageScope["UncontractedIndices"]
PackageScope["OpenSpinChains"]
PackageScope["ComplexSpinChains"]
PackageScope["HeavyMassBasisQ"]
PackageScope["HeavyTadpoleQ"]
PackageScope["HeavyTadpoles"]
PackageScope["ChargeNeutralQ"]
PackageScope["KineticCanonicalQ"]
PackageScope["GaugeVectorOccurenceQ"]
PackageScope["GaugeSingletQ"]
PackageScope["GaugeAnomaliesQ"]


PackageScope["NonGaugeFixedLagQ"]


PackageScope["GetCharge"]
PackageScope["GaugeAnomalyContribution"]


PackageScope["OccurringFields"]
PackageScope["IsolateMassTerms"]
PackageScope["DetermineCouplingOfGaugeField"]


PackageScope["DetermineCurrentFieldAssociation"]
PackageScope["GetFieldsUpdated"]
PackageScope["GetFieldsUpdatedByProperty"]
PackageScope["$currentFieldAssociation"]


PackageScope["GetGaugeCouplingsFromLagrangian"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


CheckLagrangian::usage = "CheckLagrangian[Lagrangian, options] performs a series of checks on the input Lagrangian.
The options {Hermiticity, ContractedIndices, ClosedSpinChains, CanonicallyNormalized, HeavyMassBasis, ChargeNeutral, FreeOfHeavyTadpoles, GaugeInvariant, UndefinedObject, GaugeAnomalies} take values True/False and determine which of the checks are performed.
If the option DetailedOutput->True, an association with the result of each test is returned.";


Hermiticity::usage    = ContractedIndices::usage = ClosedSpinChains::usage = CanonicallyNormalized::usage = HeavyMassBasis::usage = UndefinedObject::usage = ChargeNeutral::usage = GaugeInvariant::usage = FreeOfHeavyTadpoles::usage = GaugeAnomalies::usage = "Option for the CheckLagrangian function. It takes the values True/False to determine whether the check is performed.";
DetailedOutput::usage = "Option for the CheckLagrangian function. It takes the values True/False. An association with the result of each test is returned when set to True.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


HermitianQ::usage               = "HermitianQ[expr] returns True if the expression is hermitian, False otherwise.";
UncontractedIndices::usage      = "UncontractedIndices[expr] returns the list of open indices and mismatch between bared and not bared indices.";
OpenSpinChains::usage           = "OpenSpinChains[expr] returns the list of open spin chains in expr."
ComplexSpinChains::usage        = "ComplexSpinChains[expr] returns the list of spin chains in expr with more than 2 fermions."
HeavyMassBasisQ::usage          = "HeavyMassBasisQ[Lag] returns True if the Lagrangian is in the heavy fields mass basis, False otherwise.";
HeavyTadpoleQ::usage            = "HeavyTadpoleQ[operator] identifies terms with only one heavy field."
HeavyTadpoles::usage            = "HeavyTadpoles[Lag] returns the list of heavy tadpole terms with mass dimension less than 1. "
ChargeNeutralQ::usage           = "ChargeNeutralQ[term1 (+ term2 + ...)] returns True if the term(s) are neutral under all Abelian gauge groups defined, False otherwise.";
KineticCanonicalQ::usage        = "KineticCanonicalQ[Lag] returns True if the Lagrangian's kinetic terms are canonically normalized and diagonal, False otherwise.";
GaugeVectorOccurenceQ::usage    = "GaugeVectorOccurenceQ[expr] returns True if the expression contains gauge vectors outside field strength tensors or covariant derivatives, False otherwise";
GaugeSingletQ::usage            = "GaugeSingletQ[field] returns True if the field is a gauge singlet."
GaugeAnomaliesQ::usage          = "GaugeAnomaliesQ[Lag] returns True if the Lagrangian is gauge anomaly free, False otherwise.";


GetCharge::usage                = "GetCharge[term, Agroup, ChargeNeutral -> True/False] returns the total charge of the term under the abelian group Agroup. It issues a warning if not neutral when ChargeNeutral -> True.";
GaugeAnomalyContribution::usage = "GaugeAnomalyContribution[field/fieldList] returns the contributions to the gauge anomaly triangles from the field 'field' or the fields in the 'fieldList' list.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Validating correctness of Lagrangian*)


(* ::Subsection:: *)
(*Check modules *)


(* ::Subsubsection::Closed:: *)
(*Lagrange-like expression*)


(* ::Text:: *)
(*Throws an error if the overall form of an expression can not conceivably correspond to a Lagrangian or an operator.  *)


LagrangianLikeCheck::nofield= "The following term(s) in the expression does not contain fields or field strength tensors and the expression is not compatible with a Lagrangian or an operator: `1`";


LagrangianLikeCheck@ 0:= Null;
LagrangianLikeCheck@ expr_:= Module[{temp= TermsToList@ expr},
	If[(Or@@ FreeQ[Field|FieldStrength]/@ temp) && temp=!={0} ,
		Message[LagrangianLikeCheck::nofield, Format[Cases[temp,a_/;FreeQ[a,_Field|_FieldStrength,All],1],NiceForm]];
		Abort[];
	];
	Null
];


(* ::Subsubsection::Closed:: *)
(*Hermiticity*)


(* ::Text:: *)
(*Applies a FullSimplify to each term in a sum, while assuming that the arguments of all Log are real positive numbers, as long as the arguments do not contain explicit complex numbers.*)


SimplifyLogs[expr_]:=If[Head[expr]===Plus,
	FullSimplify[#,And@@DeleteDuplicates@Cases[#,Log[arg_/;FreeQ[arg,_Complex,All]]:>(arg>0),All]]&/@ CollectOperators[expr]
	,
	FullSimplify[expr,And@@DeleteDuplicates@Cases[expr,Log[arg_/;FreeQ[arg,_Complex,All]]:>(arg>0),All]]
]
SimplifyLogs[0]=0


(* ::Text:: *)
(*If "manifest" -> True return 2 booleans for the Lagrangian being Hermitian, and it being manifestly Hermitian.*)


HermitianQ[arg_,OptionsPattern[{"manifest"->False}]]:= Module[{expr=AtomicToNormalForm[arg]},
	If[OptionValue["manifest"],
		Module[
			{
				(* 
				Even though only GreensSimplify would be sufficient below, 
				there is a very significant performance gain by including CollectOperators and RelabelIndices before,
				especially for very long expressinos, such as the broken-phase D=8 SMEFT Lagrangian. 
				*)
				tmp= InternalSimplify[CollectOperators@ RelabelIndices[expr - Bar[expr]], TrigSimplify-> True]
			}
		,
			{
				SimplifyLogs[tmp/. lf_LF:> EvaluateLoopFunctions[lf]] === 0,
				tmp === 0
			}
		]
	,
		SimplifyLogs[InternalSimplify[CollectOperators@ RelabelIndices[expr - Bar[expr]], TrigSimplify-> True]/. 
			lf_LF:> EvaluateLoopFunctions[lf]] === 0
	]
]


(* ::Subsubsection::Closed:: *)
(*Contracted indices*)


(*UncontractedIndices@ 0= {};*)


(*UncontractedIndices[Lag_]:= Module[{indOpen=Flatten[FindOpenIndices/@ TermsToList@ Lag], indDummyBar},
	If[indOpen=!={}, Return[indOpen]];
	indDummyBar=Flatten[Complement[Bar/@ FindDummyIndices@ #,
			Flatten[(List@@ #)//. {
				Power[x_,n_]/; Element[n, NonNegativeIntegers]:> ConstantArray[x, n],
				FieldStrength[_, ind_, indG_, indL_]:> Flatten@ Join[ind, indG, indL],
				_?NumericQ-> Nothing,
				_Proj-> Nothing, 
				NCM-> List,
				DiracProduct-> List,
				GrassmannProduct-> List,
				GammaM[ind__]:> List@ ind, 
				Field[_,Vector[indV_], ind_, indL_]:> Join[List@ indV, ind, indL],
				Field[_,_, ind_, indL_]:> Join[ind,indL], 
				Coupling[_, ind_, _]:> ind, 
				CG[_,ind_]:> ind
			}/. PartialDev-> Identity] ]
		& /@TermsToList@ Lag];
	If[indDummyBar=!={}, Return[Bar/@indDummyBar]];
	{}
]*)


(* ::Text:: *)
(*Returns a list w/ all problem indices across all terms in a Lagrangian/expression*)


UncontractedIndices::open= "The indices `1` in are uncontracted in the term `2`.";
UncontractedIndices::nonconj= "The complex index '`1`' is not contracted with a CONJUGATE index in the term `2`.";


UncontractedIndices@ lag_:= Flatten[UncontractedIndicesInTerm/@ TermsToList@ lag]


UncontractedIndicesInTerm@ term_:= Module[{conjugates, dummyInds, ind, problemInds, explicitTerm},
	(*Tripple indices are directly flagged w/ FindIndices*)
	{dummyInds, problemInds}= FindIndices@ term;
	If[Length@ problemInds > 0, 
		Message[UncontractedIndices::open, Format[problemInds, NiceForm], Format[term, NiceForm]];
	];
	
	(*The complex dummy indices*)
	conjugates= Complement[Bar/@ dummyInds, dummyInds];
	(*Flatten out Bars so they appear directly on the conjugate indices (not in front of the entire field object)*)
	explicitTerm= PseudoTimes@ term/. Bar@ f:_Field|_FieldStrength:> Map[Bar, f, 3];
	Do[
		(*Check if there are both index and it's bared*)
		Switch[Count[explicitTerm, ind, All]
		,2,
			Message[UncontractedIndices::nonconj, Format[ind, NiceForm], Format[term, NiceForm]];
			AppendTo[problemInds, Bar@ ind];
		,0,
			Message[UncontractedIndices::nonconj, Format[Bar@ ind, NiceForm], Format[term, NiceForm]];
			AppendTo[problemInds, Bar@ ind];
		,_,
			Null
		];
	, {ind, conjugates}];
	problemInds
]


(* ::Subsubsection::Closed:: *)
(*Open spin chains*)


(* ::Text:: *)
(*Identifies and returns open and complex spin chains in an expression*)


OpenSpinChains[expr_]:= Select[Flatten@Join[Cases[{expr}, _NCM,Infinity],Cases[{{expr}/._NCM->0},Field[_,Fermion, ___] ,All]], OpenSpinChainQ@# &]


ComplexSpinChains[expr_]:= Select[Cases[{expr}, _NCM,Infinity], Count[#,Fermion, Infinity]>2 &]


(* ::Text:: *)
(*Identifies inconsistent spin chains *)


InconsistentSpinChains@ expr_:= Cases[{expr},
		NCM[Transp@ _Field, DiracProduct[Except[GammaCC], ___], ___]|
		NCM[___, DiracProduct[Except[GammaCC], ___], Transp@ Bar@ _Field]
	, Infinity]


(* ::Subsubsection::Closed:: *)
(*Canonically normalized kinetic term*)


(* ::Text:: *)
(*Checks if all kinetic terms (not masses!) are canonically normalized.*)


KineticCanonicalQ[L_]:=Module[{bkgFields, fields, gaugeFields, LShouldBe, LIs},
	bkgFields= FieldByProperty[BackgroundField-> True];
	fields= Complement[OccurringFields@ L, bkgFields];
	gaugeFields= Intersection[fields, List@@( #@ Field &)/@ $GaugeGroups[]];
	
	LShouldBe= KineticTerms@ Sum[FreeLag@ f, {f, fields}];
	LIs= KineticTerms@ L/. Field[lab_, __]/; MemberQ[bkgFields, lab]-> 0;
	(*Do not compare the gauge field normalization, since the coupling coulg be different from expected*)
	GreensSimplify[LShouldBe -LIs/. Thread@ Rule[gaugeFields, 0]] === 0
];


(* ::Subsubsection::Closed:: *)
(*HeavyMassBasisQ*)


HeavyMassBasisQ::nomasscoupling  = "There is no Coupling in the coefficient of the heavy mass terms `1`."
HeavyMassBasisQ::mixedfields     = "There should not be any mixed field mass terms with coupling with EFT order 0."


(*check if heavy mass is diagonal in the input*)
HeavyMassBasisQ[Lag_]:=Module[{terms,coupling,fields,indicestypes},
	terms = Collect[Contract@ Operator@ IsolateMassTerms[Lag, Heavy->True], _Operator] //TermsToList;
	Catch[
	If[terms=={0},Throw[True]];
	(*for each term...*)
	(fields = Cases[#,_Field,Infinity];
	(*...check if it is twice the same field*)
	If[!MatchQ[fields, {Field[f_,__],Field[f_,__]}], Message[HeavyMassBasisQ::mixedfields]; Throw[False]];
	(*...check there is at least one Coupling[...] in the coefficient*)
	If[FreeQ[#, _Coupling, All], Message[HeavyMassBasisQ::nomasscoupling,#]; Throw[False]]
	) & /@ terms;
	Throw[True]]
];


(* ::Subsubsection::Closed:: *)
(*HeavyTadpoleQ*)


HeavyTadpoleQ[(c_:1)op_Operator, cond_]:=
		MatchQ[OperatorType[op,CountEoMDerivatives->True], {{Field[f_,_,__]},0}/; (cond[f])]


HeavyTadpoles[Lag_]:=Module[{tadterms, massterms=IsolateMassTerms[Lag,Heavy->True], heavyQ},
	heavyQ = !FreeQ[massterms,#,All]&; (*heavy masss mixing was checked before*)
	tadterms = (Operator[Lag]//.o_Operator /;(!HeavyTadpoleQ[o, heavyQ]):>0) //TermsToList;
	tadterms = Select[tadterms, OperatorDimension[#/._Operator:>1]<1 &];
	tadterms
];


(* ::Subsubsection::Closed:: *)
(*TotalCharge*)


GetCharge::ChargeNeutral = "The term `1` is not invariant under `2` gauge group."
GetCharge::Symbols = "The term `1` is not invariant under `2` gauge group."


GetCharge[field_Symbol,Agroup_]:=Cases[$FieldAssociation[field][Charges],_Agroup]/.{Agroup[x_]}:>x/.{}->0;
GetCharge[Field[f_,___],Agroup_]:=GetCharge[f,Agroup]
GetCharge[Bar@Field[f_,___],Agroup_]:=-GetCharge[f,Agroup]
GetCharge[FieldStrength[f_,___],Agroup_]:=GetCharge[f,Agroup]
GetCharge[Bar@FieldStrength[f_,___],Agroup_]:=-GetCharge[f,Agroup]
GetCharge[list_List,Agroup_]:=GetCharge[#,Agroup]& /@list


GetCharge[term:Except[_Symbol|_Field|_FieldStrength|_Bar|_Power|_List],Agroup_,OptionsPattern[{ChargeNeutral->False}]]:=Module[{tot},
			tot=Total[GetCharge[#,Agroup]& /@ Cases[RemovePower@term, Bar[_Field|_FieldStrength], Infinity]] //Expand;
			tot=1/2 Total[GetCharge[#,Agroup]& /@ Cases[RemovePower@term, _Field|_FieldStrength, Infinity]]+  tot //Expand;
			If[ tot =!= 0 && OptionValue@ChargeNeutral == True, Message[GetCharge::ChargeNeutral,Format[term,NiceForm], Agroup]];
			Return[tot]
			]


TotalCharge[Lag_]:= Total@(Abs/@Flatten@Outer[GetCharge, TermsToList@ Lag,Keys@Select[$GaugeGroups,#[Abelian] === True&],{ChargeNeutral->True}]);


(* ::Subsubsection::Closed:: *)
(*Gauge non-invariant objects*)


(* ::Text:: *)
(*NonGaugeInvariant objects*)


NonHomogeneousTransformingQ@ expr_:= GaugeVectorOccurenceQ[expr] || !FreeQ[Lag, _PartialDev];


(* ::Text:: *)
(*Checks for isolated gauge fields *)


GaugeVectorOccurenceQ[expr_]:= !FreeQ[expr, Field[Alternatives@@(Query[Transpose]@GetGaugeGroups[])[Field],___], Infinity];


(* ::Subsubsection::Closed:: *)
(*Is the Lagrangian gauge-fixed?*)


(* ::Text:: *)
(*NonGaugeFixedLagQ returns True if and only if we are sure that the Lagrangian has not been gauge fixed. That is partial fixing returns False!*)


NonGaugeFixedLagQ::nongaugekinfix= "Could not find the kinetic term fixing for the gauge field `1`.";
NonGaugeFixedLagQ::noghostkin= "Could not find a ghost kinetic for the ghost associated to `1`.";
NonGaugeFixedLagQ::vectorkinfix= "Could not find the kinetic term fixing for the massive vector `1`.";
NonGaugeFixedLagQ::unitaryghost= "Found a ghost for a massive vector `1` in a presumed unitary gauge.";


Options@ NonGaugeFixedLagQ= {DetailedOutput-> False};


NonGaugeFixedLagQ[lag_, OptionsPattern[]]:= Module[{agh, gauges, gbs, field, gh, massiveVectors, messages,
		nonFixed= True, terms= TermsToList@ Operator@ lag},
	massiveVectors= MassiveVectorsInLag[lag, Light-> True];
	
	(*If the Lagrangian is gauge invariant w/o massive vectors, it isn't gauge-fixed*)
	If[Length@ massiveVectors === 0 && !NonHomogeneousTransformingQ@ lag,
		If[OptionValue@ DetailedOutput,
			Return@ {True, <||>};
		,
			Return@ True;
		];
	];
	
	(*Reaps messages for missing gauge-fixing components*)
	{messages, gauges}= Reap[
			(*Check that each gauge fields has been fixed*)
			Do[
				If[GaugeFixedGaugeBosonQ[terms, field], nonFixed= False]
			, {field, GaugeFieldsInLag@ lag}];
			
			(*Check that each massive vector has been fixed*)
			Do[
				If[GaugeFixedMassiveVectorQ[terms, field], nonFixed= False]
			, {field, massiveVectors}];
		, {$msg, $gauge}, Sequence@@ #2&][[2]];
	
	(*Display messages if the Lagrangian is not found to be non-fixed*)
	If[!nonFixed, ReleaseHold@ messages];
	
	If[OptionValue@ DetailedOutput,
		{nonFixed, Association@@ gauges}
	,
		nonFixed
	]
];


(* ::Text:: *)
(*Returns all the gauge fields of the Lagrangian*)


GaugeFieldsInLag@ lag_:= Module[{vectorlabels},
	vectorlabels= DeleteDuplicates@ Cases[lag, FieldStrength[label_, __]-> label, All];
	Complement[vectorlabels, MassiveVectorsInLag[lag, Light-> True]]
];


(* ::Text:: *)
(*Returns all the Goldstone bosons of the Lagrangian*)


GoldstoneBosonsInLag@ lag_:= Module[{scalars},
	scalars= DeleteDuplicates@ Cases[lag, Field[label_, Scalar, __]-> label, All];
	Intersection[scalars, GetFieldsByProperty[GoldstoneBoson-> True]]
];


(* ::Text:: *)
(*Checks if the kinetic gauge-fixing and ghost kinetic terms are present for a gauge field (label). Otherwise sows errors and returns false.*)


GaugeFixedGaugeBosonQ[terms_, label_]:= Module[{abelian, gauge, gh, agh, pass= True},
	If[!MemberQ[terms, Times[Operator[Field[label, _Vector, __, l1_], Field[label, _Vector, __, l2_]], _.]/; Length@ l1 + Length@ l2 === 2],
		Sow[Hold@ Message[NonGaugeFixedLagQ::nongaugekinfix, label], $msg];
		pass= False;
	];
	
	gauge= If[MemberQ[terms, Times[Operator[Field[label, Vector@ mu_, _, {_@ mu_}], Field[label, Vector@ nu_, _, {_@ nu_}]], 
			Power[\[Xi]Gauge[label], -1], _.]],
		R\[Xi], Indeterminate];
	
	(*Terminate if Abelian *)
	If[$GaugeGroups[First@ GaugeGroupByProperty[Field-> label], Abelian], 
		Sow[label-> gauge, $gauge];
		Return@ pass;
	];
	
	(*Look for ghost fields*)
	{gh, agh}= $GhostFields@ label;
	If[FreeQ[terms, Times[Operator[GrassmannProduct[Field[agh, AntiGhost, __, l1_], Field[gh, Ghost, __, l2_]]], _.]/; 
		Length@ l1 + Length@ l2 === 2]
	,
		Sow[Hold@ Message[NonGaugeFixedLagQ::noghostkin, label], $msg];
		gauge= Indeterminate;
		pass= False
	];
	
	Sow[label-> gauge, $gauge];
	pass
];


(* ::Text:: *)
(*Checks if the kinetic gauge-fixing and ghost kinetic terms are present for a massive vector field (label). Otherwise sows errors and returns false.*)


GaugeFixedMassiveVectorQ[terms_, label_]:= Module[{abelian, gh, agh, gaugeGFterm, gbs, ghostKinTerms, messages={}},
	If[!MemberQ[terms, Times[Operator[Bar@ Field[label, _Vector, __, l1_] | Field[label, _Vector, __, l1_], Field[label, _Vector, __, l2_]],
		_.]/; Length@ l1 + Length@ l2 === 2]
	,
		AppendTo[messages, Hold@ Message[NonGaugeFixedLagQ::vectorkinfix, label]];
		gaugeGFterm= None
	];  
	
	If[gaugeGFterm =!= None,
		gaugeGFterm= If[MemberQ[terms, Times[Operator[Field[label, Vector@ mu_, _, {_@ mu_ | mu_}], Field[label, Vector@ nu_, _, {_@ nu_| nu_}]]| 
			Operator[Bar@ Field[label, Vector@ mu_, _, {_@ mu_ | mu_}], Field[label, Vector@ nu_, _, {_@ nu_| nu_}]], 
			Power[\[Xi]Gauge[label], -1], _.]],
		R\[Xi], Indeterminate];
	];
	
	(*Look for ghost fields*)
	{gh, agh}= $GhostFields@ label;
	ghostKinTerms= True;
	If[FreeQ[terms, Times[Operator[GrassmannProduct[Field[agh, AntiGhost, __, l1_], Field[gh, Ghost, __, l2_]]], 
		_.]/; Length@ l1 + Length@ l2 === 2]
	,
		AppendTo[messages, Hold@ Message[NonGaugeFixedLagQ::noghostkin, label]];
		ghostKinTerms= False
	];  
	(*Check if bared-term is present (for complex ghosts)*)
	If[!$FieldAssociation[label, SelfConjugate] && 
			FreeQ[terms, Times[Operator[GrassmannProduct[Bar@ Field[agh, AntiGhost, __, l1_], Bar@ Field[gh, Ghost, __, l2_]]], 
				_.]/; Length@ l1 + Length@ l2 === 2]
	,
		AppendTo[messages, Hold@ Message[NonGaugeFixedLagQ::noghostkin, Conj@ label]];
		If[ghostKinTerms, Indeterminate];
	];  
	
	(*Consider checking if GBs are present*)
	gbs= KinMixingScalars[terms, label];
	
	Switch[{gaugeGFterm, ghostKinTerms, gbs}
	, {None, False, {}}, 
		Sow[label-> Unitary, $gauge];
		True
	, {R\[Xi], True, {}|{_}}, 
		Sow[label-> R\[Xi], $gauge];
		True
	, {None, False, {_}}, (*Non-fixed*)
		Sow[label-> Indeterminate, $gauge];
		False
	, _, (*Probably (partially) fixed*)
		Do[
			Sow[msg, $msg];
		, {msg, messages}];
		Sow[label-> Indeterminate, $gauge];
		True
	]
];


(* ::Text:: *)
(*Returns the kinetic-mixing scalars associated with the vector*)


KinMixingScalars[terms_List, vecLabel_]:= DeleteDuplicates@ Cases[terms/. Bar@ f:_Field-> f, 
	(Times[Operator[Field[s_, Scalar, __], Field[vecLabel, Vector@ a_, _, {a_}]], _.] |
	Times[Operator[Field[s_, Scalar, _, {a_}], Field[vecLabel, Vector@ a_, __]], _.]):> s
];


(* ::Subsection:: *)
(*Gauge anomalies*)


(* ::Subsubsection::Closed:: *)
(*GaugeSingletQ*)


(* ::Text:: *)
(*Function returning True/false if a field is a gauge singlet *)


GaugeSingletQ::NotAField = "'`1`' is not a previously defined field.";


GaugeSingletQ[field_]:=Module[{},
(* Check that the argument is a previously defined field *)
	If[!MemberQ[Keys@$FieldAssociation,field],
		Message[GaugeSingletQ::NotAField,field];
		Abort[]
	];

	If[Length@ GaugeGroupByProperty[Field-> field] > 0, 
		Return@ False;
	];

	Join[$FieldAssociation[field][Charges], 
		Intersection[Keys@ $GaugeGroups, GroupFromRep/@ $FieldAssociation[field][Indices]]] === {}
]


(* ::Subsubsection::Closed:: *)
(*Function that compute gauge anomalies*)


GaugeAnomalyContribution[fieldList_List]:=Merge[GaugeAnomalyContribution[#]&/@fieldList,Total];


GaugeAnomalyContribution[field_Symbol]:=Module[
	{
	i,j,k,
	A=Global`d$$1,B=Global`d$$2,C=Global`d$$3,
	ind=$FieldAssociation[field][Indices],
	gaugeInd=Select[$FieldAssociation[field][Indices],MemberQ[Keys@$GaugeGroups,GroupFromRep@#]&],
	charges=$FieldAssociation[field][Charges],
	chiralFactor=If[$FieldAssociation[field][Chiral]===LeftHanded, 1, -1],
	nonAbelianAnomalies,
	nonAbelianAnomalyFactors,
	mixedAnomalies,
	mixedAnomalyFactors,
	AbelianAnomalies,
	AbelianAnomalyFactors,
	gravityAnomalies,
	gravityAnomalyFactors
	},

	If[$FieldAssociation[field][Chiral]===False||GaugeSingletQ[field]===True,Return[<||>]];

	(* Contruct all possible non-Abelian triangles and compute the corresponding anomalies, given by Tr[{T[rep],T[rep]}T[rep]]*)
	nonAbelianAnomalies=Table[{i,i,i},{i,gaugeInd}];
	nonAbelianAnomalyFactors=Head/@#->If[Length[ind]>1,Times@@DimRep/@DeleteCases[ind,#[[1]]],1]chiralFactor ContractCGs[(CG[gen[#[[1]]],{A,i,j}]CG[gen[#[[3]]],{B,j,k}]+CG[gen[#[[2]]],{B,i,j}]CG[gen[#[[1]]],{A,j,k}])CG[gen[#[[3]]],{C,k,i}]]&/@nonAbelianAnomalies;

	(* Contruct all possible mixed non-Abelian/Abelian triangles and compute the corresponding anomalies, given by Tr[T[rep] T[rep]] U1[Charge] *)
	mixedAnomalies=Flatten[Table[{i,i,j},{i,gaugeInd},{j,charges}],1];
	mixedAnomalyFactors=Head/@#->If[Length[ind]>1,Times@@DimRep/@DeleteCases[ind,#[[1]]],1]chiralFactor ContractCGs[CG[gen[#[[1]]],{A,i,j}]CG[gen[#[[2]]],{B,j,i}]] #[[3]][[1]]&/@mixedAnomalies;

	(* Contruct all possible Abelian triangles and compute the corresponding anomalies, given by U1[Charge]^3 *)
	(*AbelianAnomalies=If[Length[charges]===1,
		Table[{i,i,i},{i,charges}],
		Flatten[DeleteCases[#/.{},_?NumericQ]]&/@(Expand[(Plus@@charges)^3]/.{Times->List,Plus->List,Power->Apply[List]@*ConstantArray})
	];*)
	AbelianAnomalies=Switch[Length[charges],
		0,
			{}, (* Exception for purely non-Abelian theories *)
		1,
			Table[{i,i,i},{i,charges}],
		_,
			Flatten[DeleteCases[#/.{},_?NumericQ]]&/@(Expand[(Plus@@charges)^3]/.{Times->List,Plus->List,Power->Apply[List]@*ConstantArray})
	];
	AbelianAnomalyFactors=Head/@#->If[Length[ind]>0,Times@@DimRep/@ind,1]chiralFactor #[[1]][[1]] #[[2]][[1]]#[[3]][[1]]&/@AbelianAnomalies;

	(* Contruct all possible mixed Gravity/Abelian triangles and compute the corresponding anomalies, given by U1[Charge] *)
	gravityAnomalies=Table[{"Gravity"[],"Gravity"[],i},{i,charges}];
	gravityAnomalyFactors=Head/@#->If[Length[ind]>0,Times@@DimRep/@ind,1]chiralFactor #[[3]][[1]]&/@gravityAnomalies;

	Association@Join [nonAbelianAnomalyFactors,mixedAnomalyFactors,AbelianAnomalyFactors,gravityAnomalyFactors]
]


(* ::Subsubsection::Closed:: *)
(*Function that tests if the Lagrangian has anomalies*)


GaugeAnomaliesQ[lagrangian_]:= Total@ Values@ GaugeAnomalyContribution@ 
	Complement[OccurringFields@ lagrangian, FieldByProperty[BackgroundField-> True]] === 0;


(* ::Subsection:: *)
(*Check Lagrangian function*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


CheckLagrangian::Hermiticity          = "The input Lagrangian is not hermitian. The non-Hermitian part \[ScriptCapitalL]-\[ScriptCapitalL]\[ConjugateTranspose] is `1`";
CheckLagrangian::ContractedIndices    = "The input Lagrangian contains open indices or a mismatch of bared and not bared indices.";
CheckLagrangian::ClosedSpinChains     = "The input Lagrangian contains open, complex, or inconsistent (using transposed spinors with no corresponding C-matrix) spin chains."
CheckLagrangian::CanonicallyNormalized = "The input Lagrangian kinetic part is not canonically normalized."
CheckLagrangian::HeavyMassBasis       = "The input Lagrangian is not in the mass basis for heavy fields."
CheckLagrangian::FreeOfHeavyTadpoles  = "The input Lagrangian contains heavy tadpoles that should be removed."
CheckLagrangian::ChargeNeutral        = "The input Lagrangian is not fully neutral."
CheckLagrangian::GaugeInvariant = "The Lagrangian is not gauge invariant: it contains either gauge vector bosons outside of covariant derivatives or field strength tensors, or partial derivative indices (PartialDev)."
CheckLagrangian::UndefinedObject      = "The object `1` is undefined."
CheckLagrangian::GaugeAnomalies       = "The input Lagrangian is not gauge anomaly free."


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[CheckLagrangian]={
		CanonicallyNormalized-> True,
		ChargeNeutral-> True,
		ContractedIndices-> True,
		ClosedSpinChains-> True,
		DetailedOutput-> False,
		GaugeInvariant-> True,
		GaugeAnomalies-> True,
		HeavyMassBasis-> True,
		FreeOfHeavyTadpoles-> True,
		Hermiticity-> True,
		UndefinedObject-> True,
		ZeroModeMixing-> True
	};


(* ::Subsubsection::Closed:: *)
(*CheckLagrangian*)


(* ::Text:: *)
(*Function for checking if the Lagrangian is valid for use in the matching routines *)


CheckLagrangian[Lagrangian_,opt:OptionsPattern[]]? OptionsCheck:= Module[
		{Lag= RelabelIndices@ HcExpand@ Lagrangian, DetOutput=<||>, out, OSpinChains, HeavTadpoles, 
		TotCharge, UncIndices, ExtraHeads, mHermiticity= False, mContractedIndices= False, 
		mClosedSpinChains= False, mCanonicallyNormalized= False, mMassBasis= False, mHeavyTadpoles= False,
		mChargeNeutral= False, mGaugeInvariant= False, mUndefinedObject= False, mGaugeAnomalies= False, mZeroModeMixing= False},
	
	(* Remove ghosts from the check for now *)
	Lag= Lag/. Field[_, Ghost, __]-> 0;
	
	(*Checks that Lag is a series of terms with fields or FS tensors*)
	LagrangianLikeCheck@ Lag;

	(*check if all spin chains are closed and contain only 2 fermions*)
	If[OptionValue@ ClosedSpinChains,
		OSpinChains=Join[OpenSpinChains[Lag],ComplexSpinChains[Lag], InconsistentSpinChains[Lag]];
		If[ (mClosedSpinChains = OSpinChains =!= {}) ,
			Message[CheckLagrangian::ClosedSpinChains];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"Open/Complex/InconsistentSpinChains"->(Format[#,NiceForm]&/@OSpinChains)]]
		];

	(*check if L is hermitian *)
	(*HermitianQ cannot evaluate improperly contracted spinchains*)
	If[OptionValue@Hermiticity && !mClosedSpinChains,
		If[(mHermiticity = !HermitianQ[Lag]),
			Message[CheckLagrangian::Hermiticity,Format[(Lag - Bar[Lag] //GreensSimplify),NiceForm]];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"Hermiticity"->!mHermiticity]]
		];

	(*check if all indices are contracted *)
	If[OptionValue@ ContractedIndices,
		UncIndices = UncontractedIndices[Lag];
		If[ (mContractedIndices = UncIndices =!= {}),
			Message[CheckLagrangian::ContractedIndices,UncIndices];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"UncontractedIndices"->UncIndices]]
		];

	(*check if kinetaic part is canonically normalized*)
	If[OptionValue@ CanonicallyNormalized,
		If[ (mCanonicallyNormalized = !KineticCanonicalQ[Lag]) ,
			Message[CheckLagrangian::CanonicallyNormalized];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"CanonicallyNormalized"->!mCanonicallyNormalized]]
		];

	(*check if the heavy mass is diagonal*)
	If[OptionValue@ HeavyMassBasis,
		If[ (mMassBasis= !HeavyMassBasisQ[Lag]),
			Message[CheckLagrangian::HeavyMassBasis];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"HeavyMassBasis"->!mMassBasis]]
		];

	If[OptionValue@ FreeOfHeavyTadpoles,
		HeavTadpoles=HeavyTadpoles[Lag];
		If[ (mHeavyTadpoles= HeavTadpoles =!= {}),
			Message[CheckLagrangian::FreeOfHeavyTadpoles];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"HeavyTadpoles"->(Format[#,NiceForm]&/@HeavTadpoles)]]
		];

	(*check if L is neutral under all charges*)
	If[OptionValue@ ChargeNeutral,
		TotCharge=TotalCharge[Lag];
		If[ (mChargeNeutral= ((TotCharge=!=0) && FreeQ[TotCharge,_Symbol])) ,
			Message[CheckLagrangian::ChargeNeutral];
			];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"Charges"->TotCharge]]
		];

	(*check gauge invariance: no bare gauge vectors, no partial derivative indices*)
	If[OptionValue@ GaugeInvariant,
		mGaugeInvariant= NonHomogeneousTransformingQ@ Lag;
		If[mGaugeInvariant, Message[CheckLagrangian::GaugeInvariant]];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"GaugeInvariant"->!mGaugeInvariant]]
		];

	(*check if all objects are defined*)
	If[OptionValue@ UndefinedObject,
		ExtraHeads= DeleteDuplicates@ DeleteCases[Flatten[{Lag}//.
				{Times| Plus| NCM| GrassmannProduct|LF->List,
				Power[a_,___]:> a, Bar[a_]:> a , Transp[a_]:> a, Log[a_]:>a,
				hbar->1, \[Mu]bar2->1, ev-> 1, \[Epsilon]->1, \[Pi]->1, \[ScriptD]->1,
				Thermalhbar->1, \[CapitalLambda]bar2->1, \[CapitalLambda]\[CapitalTau]->1, \[CapitalTau]->1, Zeta[a_]:>a}],
			_Field | _CG | _Coupling | _?NumberQ | _DiracProduct  | _FieldStrength | _LCTensor |  _Delta | _\[Xi]Gauge | _FlavorSum];
		If[(mUndefinedObject)= (ExtraHeads=!={}),
			Message[CheckLagrangian::UndefinedObject,#]& /@ ExtraHeads;
		];
		If[OptionValue@DetailedOutput, AppendTo[DetOutput,"AllObjectsDefined"->!mUndefinedObject]]
	];

    (*check if the Lagrangian is gauge anomaly free *)
    If[OptionValue@ GaugeAnomalies,
		If[(mGaugeAnomalies = !GaugeAnomaliesQ@Lag),
			Message[CheckLagrangian::GaugeAnomalies];
		];
		If[OptionValue@DetailedOutput,
			AppendTo[DetOutput,"GaugeAnomalies"-> Format[Normal@GaugeAnomalyContribution@OccurringFields@Lag, NiceForm]]
		];
	];
	
	(*check is Lagrangian contaings terms with a mixing of zero mode and regular fields*)
	If[OptionValue@ ZeroModeMixing,
		(mZeroModeMixing=#[[1]];
		If[OptionValue@DetailedOutput, 
			AppendTo[DetOutput,"ZeroModeMixing"->Format[#[[2]],NiceForm]]
		])& @CheckZeroModeMixing[Lag, OptionValue@DetailedOutput]
	];

	(* output the test results combined, we can output more details if wanted *)
	If[OptionValue@DetailedOutput,
		Return[DetOutput]
	];
	out= Nor[mHermiticity, mContractedIndices, mClosedSpinChains, mCanonicallyNormalized,mMassBasis,
		mChargeNeutral, mHeavyTadpoles, mGaugeInvariant, mUndefinedObject, mGaugeAnomalies, mZeroModeMixing];
	(*Cache results if Lagrangian is verified True*)
	If[out,
		CheckLagrangian[Lagrangian]= out;
	];
	out
]


(* ::Section:: *)
(*Get information from Lagrangian*)


(* ::Subsection:: *)
(*Get Lagrangian information*)


(* ::Text:: *)
(*Returns a list of the labels of all fields occurring in the Lagrangian*)


OccurringFields[L_]:= DeleteDuplicates@ Cases[L, (Field|FieldStrength)[lab_, __]-> lab, All];


(* ::Subsubsection::Closed:: *)
(*Get gauge couplings*)


(* ::Text:: *)
(*Returns an association of the gauge couplings for the gauge fields of the Lagrangian.*)
(*Power-> -1/2:     g  *)
(*Power-> 1:     1/g^2  *)


Options@ GetGaugeCouplingsFromLagrangian= {
		Power-> -1/2, (*Power applied to the coefficients of the gauge-kinetic terms*)
		EFTOrder-> Automatic, (*Automatic| All| Integer*)
		LoopOrder-> 1
	};


GetGaugeCouplingsFromLagrangian::mixing= "The Lagrangian has a kinetic mixing term between the gauge fields `1` and `2`. This is not currently supported and is ignored.";


GetGaugeCouplingsFromLagrangian[lag_, OptionsPattern[]]:= Module[
		{coefficients, gaugeFields, mix, mixings, maxOrder= OptionValue@ EFTOrder, expr= lag},
	(*This approach should be much fast than e.g., GetOperatorCoefficient because no mapping to AtomicOps is needed*)
	expr= TermsToList[expr/. {_Field-> 0, _LCTensor-> 0}];
	coefficients= Replace[expr, {
			Times[FieldStrength[lab_, __]^2 , coef__]/; FreeQ[{coef}, FieldStrength]:> {{lab, lab}, Times@ coef},
			Times[FieldStrength[lab1_, __], FieldStrength[lab2_, __], coef___]/; FreeQ[{coef}, FieldStrength]:> {Sort@ {lab1, lab2}, Times@ coef},
			FieldStrength[lab_, __]^2:> {{lab, lab}, 1},
			_-> Nothing (*Non-matching operators are not included, e.g. X^3-class*)
		}, {1}];
	
	(*Remove the coefficients of non-gauge fields*)
	gaugeFields= List@@( #@ Field &)/@ $GaugeGroups[]; 
	coefficients= Replace[coefficients, {fields_, _}/; !SubsetQ[gaugeFields, fields]-> Nothing, {1}];
	If[coefficients === {}, Return@ <||>;];
	
	(*Gather and normalize terms*)
	coefficients= -4* Plus@@ #&/@ GroupBy[coefficients, First-> Last];
	
	(*Check for mixing fields*)
	mixings= Cases[Keys@ coefficients, {lab1_, lab2_}/; lab1 =!= lab2];
	If[Length@ mixings > 0, 
		Do[
			Message[GetGaugeCouplingsFromLagrangian::mixing, mix[[1]], mix[[2]]];
		, {mix, mixings}];
	];
	(*Remove mixing fields*)
	coefficients= KeyMap[First, KeySelect[coefficients, SameQ@@ #&]];
	
	(*Determine EFT expansion order*)
	If[maxOrder === Automatic, maxOrder= GetMaxEFTOrder@ lag- 4];
	ExpandGaugeCouplingPower[#, OptionValue@ Power, maxOrder, OptionValue@ LoopOrder]&/@ coefficients
];


(* ::Text:: *)
(*Function for expanding out a power of a gauge coupling (from the kinetic term) (1/g^2 + \delta g)^x*)


ExpandGaugeCouplingPower[Power[c_Coupling, -2], pow:(_Intger|_Rational), __]:= Power[c, -2 pow];
ExpandGaugeCouplingPower[expr_, 1, __]:= expr;


ExpandGaugeCouplingPower[expr_, pow:(_Integer|_Rational), eftOrder_, loopOrd_Integer]:= Module[{series, \[Lambda]},
	series= BetterSeries[expr^pow, {hbar, 0, loopOrd}];
	If[eftOrder =!= All,
		 series= series/. c:Coupling[__, ord_]:> c \[Lambda]^ord;
		 series= BetterSeries[series, {\[Lambda], 0, eftOrder}]/. \[Lambda]-> 1;
	];
	series/.Power[Power[c_Coupling ,n_Integer],r_Rational]:> Power[c, n r]
];


(* ::Subsubsection::Closed:: *)
(*Extract all mass terms from a Lagrangian*)


Options@ IsolateMassTerms= {Heavy -> All};


IsolateMassTerms[lag_, OptionsPattern[]]:= Module[{terms= TermsToList@ HcExpand@ lag},
	(*A term is a mass term iff it contains no FSs or field derivatives and has exactly two fields*)
	terms= Select[terms, FreeQ[#, _FieldStrength|Field[__, {__}]] && 
		Count[PseudoTimes@ #, _Field, All] === 2 &];
	Plus@@ Switch[OptionValue@ Heavy
		,True,
			Select[terms, OperatorDimension[#, False] < 4 &]
		,False,
			Select[terms, OperatorDimension[#, False] >= 4 &]
		,All,
			terms
		]
];


(* ::Subsection:: *)
(*Update the current Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Determine  Fields  and  Masses*)


(* ::Text:: *)
(*Define  the  function  $currentFieldAssociation[]  which  works  similar  to  $FieldAssociation, but  containing  the  definitions  specific  to  the  current  Lagrangian . Furthermore, light  masses  are  treated  as  interactions  and  thus  for  a  light  but  massive  field  \[Phi]  we  have  Mass -> 0  in  $currentFieldAssociation[\[Phi]], contrary  to  $FieldAssociation[\[Phi]]*)


Options@ DetermineCurrentFieldAssociation= {Mode-> Matching};


DetermineCurrentFieldAssociation::massterm="Could not identify the mass coupling `2` for the field `1`.";


DetermineCurrentFieldAssociation[lagrangian_, OptionsPattern[]]:=Module[
	{
		allLight,
		lag = lagrangian,
		fieldLabels,
		\[ScriptCapitalL]Mass,
		\[ScriptCapitalL]Gauge,
		massTerms,
		fieldAssoc
	}
	,
	allLight= MatchQ[OptionValue@ Mode, Divergence | Evanescent | DimReduction];

	(* determine fields in new Lagrangian *)
	fieldLabels = Complement[OccurringFields@ lag, FieldByProperty[BackgroundField-> True]];

	(* extract all heavy mass terms *)
	\[ScriptCapitalL]Mass = IsolateMassTerms[lag, Heavy -> True];
	\[ScriptCapitalL]Mass = CollectOperators[\[ScriptCapitalL]Mass, NormalForm->False];

	(* Get all gauge kinetic terms*)
	\[ScriptCapitalL]Gauge = SelectOperatorClass[lag, {}, 4];

	(* loop over fields *)
	$currentFieldAssociation = <||>;
	Do[
		AssociateTo[
			$currentFieldAssociation,
			field -> If[MemberQ[Query[All, Key@ Field]@ $GaugeGroups, field],
				DetermineCouplingOfGaugeField[field, \[ScriptCapitalL]Gauge]
			,
				DetermineFieldMassProperties[field, \[ScriptCapitalL]Mass/.Except[Field[field,___], _Field]->0, allLight]
			]
		]
		,
		{field, fieldLabels}
	];

	(* Throw warnings if Heavy attribute changed *)
	If[OptionValue@Mode=!=DimReduction,
		Do[
			If[$FieldAssociation[f,Heavy]=!=$currentFieldAssociation[f,Heavy],
				Echo["The field "<> $LabelsNiceForm[Field]@ f<> " was defined as "<> If[$FieldAssociation[f,Heavy],"heavy","light"]<> 
					", but in the given Lagrangian it is "<> If[$currentFieldAssociation[f,Heavy],"heavy","light"]<> "."
				, "Warning:"];
			];
			If[$currentFieldAssociation[f,Heavy],
				If[$currentFieldAssociation[f,Mass]=!=$FieldAssociation[f,Mass],
					Echo["The heavy field "<> $LabelsNiceForm[Field]@ f<> " was defined with mass "<> $LabelsNiceForm[Coupling]@ $FieldAssociation[f,Mass]<> 
						", but the mass found in the Lagrangian is "<> $LabelsNiceForm[Coupling]@ $currentFieldAssociation[f,Mass]<> "."
					, "Warning:"];
				]
			]
			,
			{f, Keys@$currentFieldAssociation}
		];
	];
]


DetermineFieldMassProperties[l_, massterm_, allLight_:False]:=Module[
	{
		defProps = $FieldAssociation[l],
		heavy, mass,
		massTerm = Contract@ massterm
	}
	,

	(* determine mass-coupling depending on field type *)
	If[massTerm===0,
		(* massless/light fields -> light masses are treated as interactions *)
		mass = 0;
		heavy = False
	,
		(* heavy/massive fields *)
		Switch[defProps[Type],
			Scalar,
				If[defProps[SelfConjugate],
					mass = Sqrt[-2*massTerm/._Operator->1], (*real*)
					mass = Sqrt[-massTerm/._Operator->1] (*complex*)
				],
			Vector,
				If[defProps[SelfConjugate],
					mass = Sqrt[2*massTerm/._Operator->1], (*real*)
					mass = Sqrt[massTerm/._Operator->1] (*complex*)
				],
			Fermion,
				If[defProps[SelfConjugate],
					mass = -2*massTerm/._Operator->1 (*Majorana*)
				,
					Switch[defProps[Chiral],
						LeftHanded | RightHanded, (*Weyl*)
							mass = -massTerm/._Operator->1,
						False, (*Dirac*)
							mass = -massTerm/._Operator->1
					]
				]
		];

		mass = Contract@ mass;

		mass = mass /. {Sqrt[x_Coupling^2]:>x /; $CouplingAssociation[First@x, SelfConjugate]}; (* simplify Sqrt[] for real couplings *)
		If[MatchQ[mass, _Coupling],
			(* extract mass label *)
			If[$CouplingAssociation[First@ mass, EFTOrder]==0,
				mass = First@ mass;
				heavy = True
				,
				mass = 0; (* light field masses are treated as interactions *)
				heavy = False
			]
			,
			Message[DetermineCurrentFieldAssociation::massterm, l, mass];
			Abort[]
		];
	];

	(* return field properties in given Lagrangian *)
	<|
		Type          -> defProps[Type],
		Indices       -> defProps[Indices],
		Charges       -> defProps[Charges],
		SelfConjugate -> defProps[SelfConjugate],
		Chiral        -> defProps[Chiral],
		Mass          -> mass,
		Heavy         -> If[allLight, False, heavy],
		UFO           -> defProps[UFO]
	|>
]


(* ::Text:: *)
(*Function to determine the (square) gauge coupling of a given gauge field as it appears in the Lagrangian.*)


DetermineCouplingOfGaugeField[fieldLab_, lag_]:= Module[{coupling2, kinOperator, mu1, mu2, A,
		defProps= $FieldAssociation@ fieldLab},
	kinOperator= -1/4* If[Length@ $FieldAssociation[fieldLab, Indices] > 0,
			FS[fieldLab, mu1, mu2, A]^2
		,
			FS[fieldLab, mu1, mu2]^2
		];
	coupling2= GetOperatorCoefficient[lag, kinOperator];

	defProps@ Coupling= coupling2;
	defProps
]


(* ::Subsubsection::Closed:: *)
(*Return most recently determined fields association*)


$currentFieldAssociation = 0;


GetFieldsUpdated[] := If[Head[$currentFieldAssociation]===Association, $currentFieldAssociation, $FieldAssociation]


GetFieldsUpdated[f_] := If[Head[$currentFieldAssociation]===Association && KeyExistsQ[$currentFieldAssociation,f],
	$currentFieldAssociation@ f,
	$FieldAssociation@ f
]


GetFieldsUpdated[f_,prop_] := If[Head[$currentFieldAssociation]===Association && KeyExistsQ[$currentFieldAssociation,f],
	$currentFieldAssociation[f,prop],
	$FieldAssociation[f,prop]
]


GetFieldsUpdatedByProperty[]:= GetFieldsUpdated[];
GetFieldsUpdatedByProperty[propsSeq__]:= GetFieldsUpdatedByProperty@ List@ propsSeq;
GetFieldsUpdatedByProperty[props_Association]:= GetFieldsUpdatedByProperty@ Normal@ props;


GetFieldsUpdatedByProperty[props:_List|_Rule]:= If[Head[$currentFieldAssociation]===Association,
	Keys@ Select[$currentFieldAssociation, MatchQ[#, KeyValuePattern[props]]&],
	Keys@ Select[$FieldAssociation, MatchQ[#, KeyValuePattern[props]]&]
]
