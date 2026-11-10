(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`FlavorSymmetries`*)


(* ::Subtitle:: *)
(*Determine flavor invariants.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["ImposeFlavorSymmetry"]


PackageExport["FlavorInvariantCoupling"]


PackageExport["FlavorInvariants"]


(*PackageExport["GetFlavorStructures"]*)


(*PackageExport["RemoveFlavorStructures"]*)


(* ::Subsubsection::Closed:: *)
(*Scoped*)


(*PackageScope["$couplingFlavorStructure"]*)


PackageScope["DetermineConjugationPropertiesOfFlavorStructure"]


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Assigning flavor invariants to couplings*)


(* ::Subsection:: *)
(*Package scoped variable: $couplingFlavorStructure*)


(* ::Subsubsection::Closed:: *)
(*Definition, accessing, resetting*)


(*$couplingFlavorStructure= <||>;*)


(*GetFlavorStructures::usage= "GetFlavorStructures[] returns an association of all flavor invariants currently in use for all couplings.
GetFlavorStructures[c] returns a list of all flavor invariants currently associated to the coupling labeld by c."

GetFlavorStructures[coupl___]:= $couplingFlavorStructure[coupl]*)


(*RemoveFlavorStructures::usage= "RemoveFlavorStructures[] removes all flavor invariants currently in use.
RemoveFlavorStructures[c1, c2, ...] removes the flavor invaraints currently associated to the couplings labeld by c1, c2, \[Ellipsis]."

RemoveFlavorStructures[]:= ($couplingFlavorStructure= <||>);
RemoveFlavorStructures[coupl__]:= KeyDropFrom[$couplingFlavorStructure, #]&/@ {coupl}*)


(* ::Subsection:: *)
(*ImposeFlavorSymmetry*)


(* ::Subsubsection::Closed:: *)
(*Options*)


Options@ ImposeFlavorSymmetry= {
	Except->   {},
	Simplify-> True
}


(* ::Subsubsection::Closed:: *)
(*Messages*)


ImposeFlavorSymmetry::input=          "The 2nd argument of ImposeFlavorSymmetry must be a list or association.";
ImposeFlavorSymmetry::flavoredFields= "Only flavored fields can be assigned flavor symmetries.";
ImposeFlavorSymmetry::symFormat=      "The given symmetries are not of the right format.";
ImposeFlavorSymmetry::noporp=         "Could not identify the flavor transformation properties of `1` in the term: `2`";
ImposeFlavorSymmetry::multprop=       "Could not assign consistent flavor-symmetry transformation properties to the coupling `1`.";
ImposeFlavorSymmetry::multflav=       "Fields with multiple flavor indices are currently not supported.";


(* ::Subsubsection::Closed:: *)
(*Usage*)


ImposeFlavorSymmetry::usage= "ImposeFlavorSymmetry[L, sym] determine the flavor invariants for all couplings present in the expression L, according to the given flavor symmetries (sym) and and associates the invariants to the coupling.";


(* ::Subsubsection::Closed:: *)
(*Main function*)


ImposeFlavorSymmetry[L_, symmetries_, OptionsPattern[]]:= Module[
	{
		syms= symmetries,
		coups= DeleteDuplicates@Cases[L, Coupling[label_,{__},_]-> label, All], (* ignore unflavored couplings *)
		expr= TermsToList[L],
		tmp, fields, fieldsBar, ind, symsBar,
		globalGroups= GetGlobalGroups[],
		globalU1, globalAbelianGroups, globalNonAbelianGroups,
		globalReps,
		reps= GetRepresentations[],
		currentRep,
		flavoredFields= KeySelect[
			$FieldAssociation,
			!FreeQ[$FieldAssociation[#], Alternatives@@(Keys[$FlavorIndices])]&
		],
		rhsPattern, auxTensorProduct,
		indDim, indGiven
	}
	,
	(* drop couplings that should be excluded *)
	coups= DeleteCases[coups, x_/; MatchQ[x, Alternatives@@ OptionValue[Except]], 1];
	
	(* get list of global groups *)
	globalAbelianGroups=    KeySelect[globalGroups,  globalGroups[#,Abelian]&];
	globalNonAbelianGroups= KeySelect[globalGroups, !globalGroups[#,Abelian]&];
	globalReps= Flatten[(#[Representations])&/@ Values[globalNonAbelianGroups]];
	globalU1= Keys@ globalAbelianGroups;
	
	(* check input is of the right form *)
	If[!MatchQ[Head[syms], List|Association],
		(* check format of symmetries *)
		Message[ImposeFlavorSymmetry::input];
		Abort[]
		,
		If[Head[syms]===List, syms= Association@ syms]; (* transform List -> Association *)
		If[!AllTrue[Keys[syms], MatchQ[#, Alternatives@@(Keys@flavoredFields)]&],
			(* check that all fields are flavored *)
			Message[ImposeFlavorSymmetry::flavoredFields];
			Abort[]
			,
			(* generate pattern for allowed symmetry input *)
			rhsPattern= Alternatives@@ Join[
				globalReps,
				CRep/@ globalReps,
				Bar/@ globalReps
			];
			rhsPattern= Flatten@ Alternatives[
				rhsPattern, 
				Sequence@@ (#[_?NumericQ]&/@ globalU1),
				auxTensorProduct[rhsPattern],
				auxTensorProduct[OrderlessPatternSequence[rhsPattern, (Alternatives@@ (#[_?NumericQ]&/@ globalU1))..]],
				auxTensorProduct[(Alternatives@@ (#[_?NumericQ]&/@ globalU1))..],
				Singlet
			];
			With[{rhsPatternAUX= rhsPattern},
				rhsPattern= HoldPattern[rhsPatternAUX]/. auxTensorProduct-> (List|TensorProduct);
			];
			(* check symmetries satisfy the pattern *)
			If[!MatchQ[Values@ syms, {{({(_Integer)..}->rhsPattern)..}..}],
				Message[ImposeFlavorSymmetry::symFormat];
				Abort[]
			];
		]
	];
	
	(* include flavored fields of the Lagrangian for which no flavor symmetry has been specified *)
	AssociateTo[syms, #-> {}]&/@ Complement[
		Cases[Keys@flavoredFields, x_/; !FreeQ[L, x, All], All], (* all flavored fields in Lagrangian *)
		Keys@syms (* fields for which flavor symmetry was defined *)
	];
	
	(* include Singlet for omitted flavors *)
	Do[
		(* flavors for which a transformation was specified *)
		indGiven= Sort[Flatten[syms[f][[;;,1]]]];
		(* determine number of flavors *)
		indDim= #[IndexDimension]&/@ (Cases[$FieldAssociation[f,Indices], x_/; !FreeQ[x, Alternatives@@ Keys[$FlavorIndices], All]]/. $FlavorIndices);
		(* abort if a field with multiple flavor indices is present *)
		If[Length@indDim > 1, Message[ImposeFlavorSymmetry::multflav]; Abort[]];
		(* add Singlet for missing flavors *)
		syms[f]= syms[f]~Join~Table[
			{i}-> Singlet
			,
			{i, Complement[Range[First@indDim], indGiven]}
		];
		,
		{f, Keys@ syms}
	];
	
	(* compute the symmetry replacements of the bared fields *)
	symsBar= syms/. {
		(rep_/; MemberQ[globalReps,rep]&&reps[rep,Reality]=!=Real):> CRep[rep],
		(u1_[x_]/; MemberQ[Keys@globalAbelianGroups,u1]):> u1[-x]
	};
	
	(* reset flavor structures *)
	(*$couplingFlavorStructure=<||>;*)
	
	(* loop over all coupling labels *)
	Association@@ Table[
		(* get all terms containing this coupling, but ignore the terms with the bared couplings for convenience *)
		With[{coupl=coup},
			tmp= Cases[expr, x_/;(!FreeQ[x,coup,All]&&FreeQ[x,HoldPattern@Bar@Coupling[coupl,__],All]), 1];
		];
		
		(* loop over all terms containing that one coupling *)
		currentRep= Table[
			(* get indices of the coupling *)
			ind= FirstCase[term, Coupling[coup,inds_,_]:> inds, All];
			
			(* get rid of powers *)
			term= RemovePower[term];
			
			(* get the bared fields contracted into the coupling *)
			fieldsBar= Cases[
				Hold[Evaluate@term]/. Coupling[coup,___]-> 1,
				(f:(Bar@Field[l_,__]|Bar@FieldStrength[l_,__]))/;!FreeQ[f,#,All]:>l,
				All
			]&/@ ind;
			
			(* get the non-bared fields contracted into the coupling *)
			fields= Cases[
				Hold[Evaluate@term]/. Coupling[coup,___]-> 1/. Bar@_Field-> 1/. Bar@_FieldStrength-> 1,
				(f:(Field[l_,__]|FieldStrength[l_,__]))/; !FreeQ[f,#,All]:> l,
				All
			]&/@ ind;
			
			(* replace field label by given flavor transformation proprties *)
			fields=    fields/.    syms;
			fieldsBar= fieldsBar/. symsBar;
			
			(* combine reps of bared and non-bared fields *)
			Table[
				Switch[{Length[fields[[i]]],Length[fieldsBar[[i]]]},
					{1,0}, First@fields[[i]],
					{0,1}, First@fieldsBar[[i]],
					_, Message[ImposeFlavorSymmetry::noporp, coup, Format[ReleaseHold@term,NiceForm]];
					   Abort[]
				]
				,
				{i,Length@ fields}
			]
			,
			{term, tmp}
		];
		
		(* all terms for one coupling must lead to the same representations, if not throw an error *)
		currentRep= DeleteDuplicates@ currentRep;
		If[Length@currentRep > 1,
			Message[ImposeFlavorSymmetry::multprop, coup];
			Continue[]
			,
			currentRep= First@ currentRep
		];
		
		(* assign flavor structures to the coupling *)
		FlavorInvariantCoupling[
			coup,
			FlavorInvariants[currentRep, Simplify-> OptionValue[Simplify]] (* determines flavor invariants *)
		]
		,
		{coup, coups}
	]
]


(* ::Subsection:: *)
(*Assigning a flavor invariant to a coupling*)


(* ::Subsubsection::Closed:: *)
(*Messages*)


FlavorInvariantCoupling::undefcoupl= "The symbol `1` does not match a defined coupling label.";
FlavorInvariantCoupling::dims=       "The dimensions of all given flavor tensors and/or of the coupling do not match.";
FlavorInvariantCoupling::existing=   "The coupling `1` already has associated flavor structures. Overwriting...";
FlavorInvariantCoupling::conjprop=   "Unknown conjugation property `1` encountered for coupling `2`.";
FlavorInvariantCoupling::noconj=     "No conjugate could be found for coupling `1` with flavor structure `2`.";


(* ::Subsubsection::Closed:: *)
(*Usage*)


FlavorInvariantCoupling::usage= "FlavorInvariantCoupling[c, struct] associates the given flavor invariant (struct) to the coupling (c).";


(* ::Subsubsection::Closed:: *)
(*FlavorInvariantCoupling*)


(* ::Text:: *)
(*Handle the case of no invariants*)


FlavorInvariantCoupling[coup_Symbol, {}]:= coup-> {}
(*If[KeyExistsQ[$couplingFlavorStructure, coup],
	(* print warning if key already exists *)
	Message[FlavorInvariantCoupling::existing, coup]
	,
	(* otherwise associate empty list *)
	AssociateTo[$couplingFlavorStructure, coup-> {}]
]*)


(* ::Text:: *)
(*Handle the case with invariants*)


FlavorInvariantCoupling[coup_Symbol, sa_SparseArray]:= FlavorInvariantCoupling[coup, {sa}]


FlavorInvariantCoupling[coup_Symbol, structuresInput:{(_SparseArray)..}]:= Module[
	{
		coups= $CouplingAssociation,
		structures= structuresInput,
		currentStruct, symStruct, permStruct,
		syms, conj,
		res= {},
		checkmark
	}
	,
	(* check that coup matches a defiend coupling symbol *)
	If[!MemberQ[Keys@coups,coup],
		Message[FlavorInvariantCoupling::undefcoupl, coup];
		Abort[]
	];
	(* check that all SparseArray have same dimensions and these match the indices on the coupling *)
	If[Length@ DeleteDuplicates[{$FlavorIndices[#, IndexDimension]&/@ (coups[coup, Indices])}~Join~(Dimensions/@ structures)]=!=1,
		Message[FlavorInvariantCoupling::dims];
		Abort[]
		,
		(*(* print warning if key already exists *)
		If[KeyExistsQ[$couplingFlavorStructure, coup],
			Message[FlavorInvariantCoupling::existing, coup];
		];*)
		
		(* NOTE: we also need to check the SelfConjugate propertie! *)
		
		(* check for symmetries *)
		syms= coups[coup, Symmetries];
		
		(* combine the invariants according to the symmetries of the coupling *)
		While[Length[structures]>0,
			(* take first structure and remove it from list *)
			currentStruct= structures[[1]];
			structures= structures[[2;;]];
			(* create empty list for symmetry related structures *)
			symStruct= {};
			(* loop over symmetries of the current structure *)
			Do[
				(* if this symmetry permutation is equal to another structure, add it to symStruct and delete it from structures *)
				If[MemberQ[structures, permStruct= Transpose[syms[perm]* currentStruct, perm]],
					AppendTo[symStruct, permStruct];
					structures= DeleteCases[structures, permStruct, 1];
				]
				,
				{perm, Keys@ syms}
			];
			(* if there are symmetry related structures, add them up *)
			If[Length[symStruct]>0,
				currentStruct= currentStruct + Total[symStruct]
			];
			(* add the new symmetrized structure to the result *)
			AppendTo[res, currentStruct]
		];

		(* canonize all SparseArray so one can use MatchQ/SameQ below *)
		res= SparseArray/@ res;
				
		(* identify conjugate structures *)
		conj= coups[coup, SelfConjugate];
		Switch[conj,
			False, res= List/@ res, (* for intrinsically complex couplings only wrap every structure in am individual list *)
			True,  res= List/@ res, (* for intrinsically real couplings only wrap every structure in am individual list *)
			{__},     (* if a non-trivial conjugation symmetry exists loop over all structures *)
				(* set structures to res, then reset res *)
				structures= res;
				res= {};
				(* loop over all structures *)
				While[Length[structures]>0,
					currentStruct= First[structures];
					(* determine conjugate of currentStruct and drop struct from structures *)
					permStruct= SparseArray@ Transpose[currentStruct, conj]; (* canonize again using SparseArray again *)
					structures= DeleteCases[structures, currentStruct, 1];
					(* check if currentStruct is self conjugate *)
					If[MatchQ[permStruct, currentStruct | -currentStruct],
						(* currentStruct is already (anti-)selfconjugate -> add it to res *)
						AppendTo[res, {currentStruct}]
						,
						(* currentStruct is not yet (anti-)selfconjugate -> loop over all structures to find its conjugate partner *)
						checkmark= False;
						Do[
							If[permStruct===str,
								(* if str is the conjugate, append it together with currentStruct to res and delete it from structures *)
								checkmark= True;
								AppendTo[res, {currentStruct, str}];
								structures= DeleteCases[structures, str, 1];
								(* exit do-loop if the conjugate was found [there won't be a 2nd one] *)
								Break[]
							]
							,
							{str, structures}
						];
						(* throw an error if no conjugate could be found *)
						If[!checkmark, 
							Message[FlavorInvariantCoupling::noconj, coup, currentStruct]; Abort[]
						];
					]
				],
			_, Message[FlavorInvariantCoupling::conjprop, conj, coup]; Abort[]
		];
		
		(* less reliable version of the code above *)
		(*
		(* symmetrize all flavor structures *)
		res= DeleteDuplicates@ Table[
			(* DeleteDuplicates ensures that already symmetric flavor structures keep their normalization *)
			SparseArray[Plus@@ DeleteDuplicates[Table[Transpose[syms[perm]*struc, perm], {perm, Keys@ syms}]]]
			,
			{struc, structures}
		];
		*)
		
		(* delete empty entries [this should no longer happen] *)
		res= DeleteCases[#, x_/; Length[x["NonzeroValues"]]==0, 1]&/@ res;
		
		(* save to $couplingFlavorStructure *)
		(*AssociateTo[$couplingFlavorStructure, coup-> res]*)
		coup-> res
	]
]


(* ::Subsection:: *)
(*Determine the conjugation properties of a flavor invariant given certain symmetries*)


(* ::Subsubsection::Closed:: *)
(*Messages*)


DetermineConjugationPropertiesOfFlavorStructure::sym= "The flavor structure `1` is inconsitent with the conjugation property `2` and the symmetries `3`.";
DetermineConjugationPropertiesOfFlavorStructure::re= "The flavor structure `1` contains some real and some complex entries, considering the given conjugation property `2` and symmetries `3`. Its coupling is taken as real for consistentcy.";
DetermineConjugationPropertiesOfFlavorStructure::im= "The flavor structure `1` contains some imaginary and some complex entries, considering the given conjugation property `2` and symmetries `3`. Its coupling is taken as imaginary for consistentcy.";


DetermineConjugationPropertiesOfFlavorStructure::missconj= "The flavor structure `1` for the coupling `2` was determined to be complex but no conjugate was found."


(* ::Subsubsection::Closed:: *)
(*DetermineConjugationPropertiesOfFlavorStructure*)


(* ::Text:: *)
(*This function takes a flavor invariant (given as SparseArray), an index permutation under conjugation (for example {2,1,4,3}), and a list of symmetries (for example {{1,2,3,4}->1, {3,4,1,2}->1}) as input.*)
(*It then determines if the coupling associated to this flavor structure is real (Re), imaginary (Im), or complex (ReIm).*)
(*It returns a list {Re/Im/ReIm, {a,b}}, where a and b are empty lists for Re and Im, while for ReIm a contains all non-conjugated entries and b al conjugated entries.*)


DetermineConjugationPropertiesOfFlavorStructure[sa_SparseArray, conj_, symmetries_]:= 
DetermineConjugationPropertiesOfFlavorStructure[sa, conj, symmetries]= Module[
	{
		conjPerm, tab, conjProp, res, param= {}, paramConj= {}, symRelated, conjRelated,
		(* get the position and values of non-zero-entires *)
		nonZero= DeleteCases[ArrayRules@ sa, HoldPattern[_-> 0]],
		syms= Normal[symmetries]
	}
	,
	(* loop over all non-zero entries *)
	tab= Table[
		(* initialize with general complex case *)
		conjProp= ReIm;
	
		(* index permutation under conjugation *)
		conjPerm= Permute[First@ entry, conj];
		
		If[conjPerm === First@ entry,
			(* the entry is (anti-)selfconjugate *)
			Switch[Conjugate[Last@ entry],
				+Last@ entry, conjProp= Re, (* selfconjugate *)
				-Last@ entry, conjProp= Im, (* anti-selfconjugate *)
				_, Message[DetermineConjugationPropertiesOfFlavorStructure::sym, sa, conj, syms]; Abort[]
			]
			,
			(* loop over all symmetries *)
			Do[
				If[Permute[conjPerm, First@ sym] === First@ entry,
					(* if a fitting symmetry is bound check if it is real or imaginary *)
					Switch[Last[sym]*Conjugate[Last@ entry],
						+Last@ entry, conjProp= Re; Break[], (* conjugate entry equals a symmetry related entry -> it is real *)
						-Last@ entry, conjProp= Im; Break[], (* conjugate entry equals minus a symmetry related entry -> it is imaginary *)
						_, Message[DetermineConjugationPropertiesOfFlavorStructure::sym, sa, conj, syms]; Abort[]
					]
				]
				,
				{sym, syms}
			];
		];
		conjProp
		,
		{entry, nonZero}
	];
	
	(* determine return value depending on the results *)
	res= Switch[DeleteDuplicates@ tab,
		{Re},   Re,   (* fully real *)
		{Im},   Im,   (* fully imaginary *)
		{ReIm}, ReIm, (* fully complex *)
		(* mixed *)
		x_/; (FreeQ[x,Im] && !FreeQ[x,Re]), 
			Message[DetermineConjugationPropertiesOfFlavorStructure::re, sa, conj, syms]; 
			Re, (* Re and ReIm -> Re *)
		x_/; (FreeQ[x,Re] && !FreeQ[x,Im]), 
			Message[DetermineConjugationPropertiesOfFlavorStructure::im, sa, conj, syms]; 
			Im, (* Im and ReIm -> Im *)
		_, 
			Message[DetermineConjugationPropertiesOfFlavorStructure::sym, sa, conj, syms]; 
			Abort[] (* Re and Im -> inconsitent *)
	];
	
	(* if complex (ReIm) param found, find also entries that are conjugated *)
	If[res===ReIm,
		(* 
		WARNING: 
			- if the While does a 2nd run, it might use the opposite sym/conj association compared to the 1st run;
			- for exact flavor symmetries this is no problem, but once spurions are included it might cause inconsistencies.
		*)
		While[Length[nonZero]>0,
			(* find enrtries related by symmetry to 1st entry *)
			symRelated=  Permute[nonZero[[1,1]], #]&/@ syms[[;;,1]];
			(* find conjugate of 1st entry *)
			conjRelated= Permute[nonZero[[1,1]], conj];
			(* find enrtries related by symmetry to the conjugate of the 1st entry *)
			conjRelated= Permute[conjRelated, #]&/@ syms[[;;,1]];
			(* append params and conjugated params *)
			AppendTo[param, #]&/@ symRelated;
			AppendTo[paramConj, #]&/@ conjRelated;
			(* drop these terms *)
			nonZero= DeleteCases[nonZero, y_/; MemberQ[symRelated~Join~conjRelated, First[y]], 1];
		]
	];
	
	{res, {param, paramConj}}
]


(* ::Section:: *)
(*Determination of flavor invariants*)


(* ::Subsection:: *)
(*Outer product of invariants*)


(* ::Subsubsection::Closed:: *)
(*Options*)


Options@ FlavorInvariants= {
	Dimensions-> Automatic, 
	Simplify-> True
};


(* ::Subsubsection::Closed:: *)
(*Messages*)


FlavorInvariants::symFormat=    "The given symmetries are not of the right format.";
FlavorInvariants::dimMissmatch= "The dimensions on the left- and right-hand side of `1` do not match!";


(* ::Subsubsection::Closed:: *)
(*Usage*)


FlavorInvariants::usage= "FlavorInvariants[arg] receives a list of representations as argument (arg) and returns all singlets that are in teh product of all possible combinatinos of these representations.";


(* ::Subsubsection::Closed:: *)
(*FlavorInvariants*)


FlavorInvariants[arg_List, OptionsPattern[]]:= Module[
	{
		reps,
		$FlavorInvariants, (* dummy function such that Flatten does not act on SparseArray *)
		indsFinal, indsInitial, inds,
		res,
		maxDims,
		globalReps= Flatten[#[Representations]&/@ Values[GetGlobalGroups[]]],
		globalU1= Keys@ KeySelect[GetGlobalGroups[], GetGlobalGroups[#,Abelian]&],
		rhsPattern,
		auxTensorProduct
	}
	,
	(* check input is of the right form *)
	rhsPattern= Alternatives@@ Join[
		globalReps,
		CRep/@ globalReps,
		Bar/@ globalReps
	];
	rhsPattern= Flatten@ Alternatives[
		rhsPattern, 
		Sequence@@ (#[_?NumericQ]&/@ globalU1),
		auxTensorProduct[rhsPattern],
		auxTensorProduct[OrderlessPatternSequence[rhsPattern, (Alternatives@@ (#[_?NumericQ]&/@ globalU1))..]],
		auxTensorProduct[(Alternatives@@ (#[_?NumericQ]&/@ globalU1))..],
		Singlet
	];
	With[{rhsPatternAUX=rhsPattern},
		rhsPattern= HoldPattern[rhsPatternAUX]/. auxTensorProduct-> (List|TensorProduct);
	];
	If[!MatchQ[arg,{{({(_Integer)..}->rhsPattern)..}..}],
		Message[FlavorInvariants::symFormat];
		Abort[]
	];
	
	(* check dimenions of reps in input *)
	Do[
		Do[
			(* check that number of indices matches dimension of the representation *)
			(* Below we use `{z_}-> z` to drop lists of single terms and `Nothing-> 1` to simplify TensorProduct[a_,1]=TensorProduct[a] *)
			If[Length@First[y]=!=(Last[y]/.((Blank[#]-> Nothing)&/@ globalU1)/. {z_}-> z/. Nothing-> 1/. (CRep|Bar-> Identity)/. GetRepresentations[])[Dimension],
				(* check singlets have only one index *)
				If[MatchQ[Last@y, Singlet | (Alternatives@@ (#[_]&/@ globalU1))] && (Length@First[y]===1),
					Continue[]
					,
					Message[FlavorInvariants::dimMissmatch, y];
					Abort[]
				]
			]
			,
			{y, x}
		]
		,
		{x, arg}
	];
	
	(* if dimensions are not specified take the larges occuring index *)
	If[OptionValue[Dimensions]===Automatic,
		maxDims= Max/@ Flatten[arg[[;;,;;,1]]//. {a___,l1:{_Integer..},l2:{_Integer..},b___}:> {a,l1~Join~l2,b}, 1],
		maxDims= OptionValue[Dimensions]
	];
	
	(* get all combinatinos of representatinos *)
	res= Sort[Flatten@ Outer[$FlavorInvariants@{##}&, Sequence@@ arg]];
	
	(* keeps only the full U1 singlets *)
	res= DeleteCases[res, x_/;! MatchQ[Sum[Abs@Total@Cases[x,u1[n_]:>n,All],{u1,globalU1}],(0|0.)], 1];
	
	(* drop U1 charges from now on *)
	res= res/. ((Blank[#]-> Singlet)&/@ globalU1);
	
	(* remove singlets from TensorProduct *)
	res= res/. tp_TensorProduct:> If[(tp/. Singlet-> 1)===1, Singlet (* handle fully abelian case *), tp/. Singlet-> 1];
	
	(* the same but if lists are used instead of TensorProduct *)
	res= res/. ls:Except[{_Integer..}, List[__]]/; FreeQ[ls, Rule|RuleDelayed]:> If[(ls/. Singlet-> Nothing)==={}, Singlet (* handle fully abelian case *), ls/. Singlet-> Nothing/. {a_}-> a];
	
	(* compute flavor invariants and delete all with no singlets *)
	res= DeleteCases[res/. $FlavorInvariants[x__]:> {x, FlavorInvariantsInternal[x[[;;,2]], Simplify-> OptionValue[Simplify]]}, {_,{}}, 1];
	
	(* fix dimension of tensors *)
	res= Table[(* loops over all rep-products *)
		(* get all combinations of allowed indices for the embedding, these are the indices of the final full tensor *)
		indsFinal= Flatten[Outer[List, Sequence@@(First[tens][[;;,1]])], Length@First[tens]-1];
		
		(* the tensor one starts with is smaller and the final indices might be out of range *)
		indsInitial= Table[
			(* this replaces the indices into canonical ordering, e.g., if the final indices are {1,3}->{1,2} and {3,2}->{1,2} *)
			ten[[1]]/. (Rule[First@#, Last@#]&/@ Transpose@{ten[[1]], Range[Length@ten[[1]]]})
			,
			{ten, First@tens}
		];
		(* this generates the equivalent of indsFinal but with indices in canonical ordering which can be used as part specification for the initial tensor *)
		indsInitial= Flatten[Outer[List, Sequence@@indsInitial], Length@First[tens]-1];
		(* combine final and initial indices *)
		inds= Transpose@ {indsFinal, indsInitial};
		
		(* for every invariant of this product embedd it in larger structure *)
		Table[
			SparseArray[(First[#]-> y[[Sequence@@(Last[#][[Flatten@Position[First[tens][[;;,2]],Except[Singlet|List],1]]])]])&/@ inds, maxDims]
			,
			{y, Last@ tens}
		]
		,
		{tens, res}
	];
	
	(* unify into single list *)
	Flatten[res, 1]
]


(* ::Subsection:: *)
(*Individual flavor invariants*)


(* ::Subsubsection::Closed:: *)
(*Options*)


Options@ FlavorInvariantsInternal= {
	Dimensions-> Automatic, 
	Simplify-> True
};


(* ::Subsubsection::Closed:: *)
(*Messages*)


FlavorInvariantsInternal::undefRep= "The input `1` is not a list of defined representations of global groups.";


(* ::Subsubsection::Closed:: *)
(*FlavorInvariantsInternal*)


FlavorInvariantsInternal[arg_List, OptionsPattern[]]:= Module[
	{
		expr= arg/. Bar-> CRep,(* for user-friendliness allow for Bar as input but immediately replace by CRep *)
		globalGroups= GetGlobalGroups[],
		allReps= $Representations,
		globalReps,
		presentReps,
		repsToDynkinRule,
		count= 1,
		res, dims
	}
	,
	(* get all representations of global groups *)
	globalReps= Flatten[#[Representations]&/@ Values[globalGroups]];
	
	(* handle full singlet case *)
	If[MatchQ[expr, {Singlet..}], 
		Return@ {1}
	];
	
	(* drop singlets since InvariantTensors cannot handle them *)
	expr= expr/. Singlet-> Nothing;
	
	(* check argument is a list of global reps *)
	If[!MatchQ[expr/. CRep-> Identity, {(Alternatives@@globalReps)..}|{}],
		Message[FlavorInvariantsInternal::undefRep, expr];
		Abort[]
	];
	
	(* save the ordering of the representations *)
	expr= {#, count++}&/@ expr;
	
	(* group all present representations by group *)
	presentReps= GroupBy[expr, allReps[First[#]/. CRep-> Identity, Group]&];
	
	(* replacement rules from representation names to Dynkin indices *)
	repsToDynkinRule= #[DynkinCoefficients]&/@ allReps;
	
	(* compute the invariant for every group algebra present, the keys mark which indices these correspond to in the final tensor *)
	res= Association@ Table[
		(*  position of indices  \[Rule]  invariant tensor structure for the specifed group *)
		(presentReps[gr][[;;,2]])-> Quiet[
			InvariantTensors[globalGroups[gr,Group], (presentReps[gr]/. repsToDynkinRule)[[;;,1]]]
			, {InvariantTensors::nosing}
		]
		,
		{gr, Keys@ presentReps}
	];
	
	(* return empty list if no singlet was found *)
	If[MemberQ[Values[res],{}] || MatchQ[res,<||>], Return@ {}];

	(* take the product representation of the different groups *)
	res= Flatten[
		Outer[
			Transpose[ (* this traspose ensure teh indices are brought back into the right ordering *)
				Outer[Times, ##],(* this computes the tensor product between the different algebras *)
				Flatten@Keys@ res
			]&, 
			Sequence@@ Values[res],(* this considers all possible combinations if multiple invariants have been found for any of the group algebras *)
			1
		],
		Length@ Values[res]-1 (* this is to ensure the SparseArray does not get flattened out as well *)
	];
	
	(* attempt to take linear combinations of the obtained tensors  *)
	If[OptionValue[Simplify],
		res= SimplifyInvariantBasis[res]
	];
	
	(* pad with zeros to match dimenions and add singlets *)
	(* NOTE: only adding zeros at the end -> need option for inserting zero rows and columns in the middle as  well *)
	If[OptionValue[Dimensions]=!=Automatic,
		(* get starting dimensions *)
		dims= Dimensions[First@ res];
		
		(* include zeros to match required dimensions *)
		res= Table[
			SparseArray[arr["NonzeroPositions"]-> arr["NonzeroValues"], OptionValue[Dimensions]]
			,
			{arr, res}
		];
	];

	res
]


(* ::Subsection:: *)
(*Simplifications of flavor invariants*)


(* ::Subsubsection::Closed:: *)
(*Messages*)


SimplifyInvariantBasis::failed= "Failed to simplify the basis of invariant tensors. Returning the unsimplified basis."


(* ::Subsubsection::Closed:: *)
(*SimplifyInvariantBasis*)


SimplifyInvariantBasis[arg:{_SparseArray..}]:= Module[
	{
		vars= Table[ToExpression["x"<>ToString[i]], {i, Length@ arg}], (* linar combination variables *)
		terms,
		sol0, sol1, sol,
		expr,
		ordering
	}
	,
	
	(* list of all entries in the linear combination that are not zero *)
	terms= DeleteCases[Flatten@ Normal@ Sum[vars[[i]]*arg[[i]], {i, Length@ arg}], 0, 1];
	
	(* all solutions that set entries to zero *)
	sol0= DeleteDuplicates@ Flatten@ Table[
		Quiet@ Solve[entry==0, vars]
		,
		{entry, terms}
	];
	
	(* all solutions that set entries to on1 *)
	sol1= DeleteDuplicates@ Flatten@ Table[
		Quiet@ Solve[entry==1, vars]
		,
		{entry, terms}
	];
	
	(* if there si only one arg, we cannot allow for zero solutions *)
	If[Length@arg==1,
		sol0= sol0/. Rule[_,0]-> Nothing;
		sol1= sol1/. Rule[_,0]-> Nothing;
	];
	
	(* combinations of solutions *)
	If[Length@sol0 > 0,
		If[Length@sol1 > 0,
			sol= Flatten@ Outer[And, sol0/. Rule-> Equal, sol1/. Rule-> Equal]
			,
			sol= sol0/. Rule-> Equal
		],
		If[Length@sol1 > 0,
			sol= sol1/. Rule-> Equal
			,
			Message[SimplifyInvariantBasis::failed];
			Return[arg];
			Abort[] (* for emergencies *)
		]
	];
	
	(* find common solutions *)
	sol= DeleteDuplicates@ Flatten[FindInstance[#, vars]&/@ sol, 1];
	
	(* catch the case where only one solution was found *)
	If[!MatchQ[sol, {_List..}],
		sol= {sol}
	];
	
	(* get a list of original tensors and all linear combinations just found as solutions *)
	expr= DeleteDuplicates[arg~Join~ReleaseHold[Sum[vars[[i]]*Hold[Evaluate[arg[[i]]]], {i, Length@ arg}]/. sol]];
	
	(* find the ones containing the least amount of entries not matching 0 or 1 *)
	ordering= OrderingBy[expr, Count[Flatten@ Normal@ #, Except[0|1]]&][[;;Length@arg]];
	
	(* get the best resluts *)
	expr= SparseArray/@ expr[[ordering]];
	
	(* check linear independence *)
	If[MatrixRank[Flatten/@ Normal/@ expr]=!=Length[expr],
		Message[SimplifyInvariantBasis::failed];
		expr= arg
	];
	
	expr
]
