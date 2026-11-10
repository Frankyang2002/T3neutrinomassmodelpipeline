(* ::Package:: *)

(* ::Title:: *)
(*SM model file with electroweak symmetry breaking*)


(* ::Section:: *)
(*Loading the unbroken phase Lagrangian*)


ParentModel["SM"]


(* ::Section:: *)
(*General broken EW phase definition*)


(* ::Subsection:: *)
(*Broken phase gauge groups*)


(* ::Text:: *)
(*Define QED as the gauge group of the unbroken phase:*)


DefineGaugeGroup[U1em, U1, \[ScriptE], \[ScriptCapitalA], 
	UFO -> <|
		"pdg"->      22, 
		"name"->     "A", 
		"coupling"-> "ge"
	|>
];


(* ::Subsection:: *)
(*Broken phase fields*)


(* ::Subsubsection:: *)
(*Massive vectors*)


(* ::Text:: *)
(*Define the massive \[ScriptCapitalZ] and \[ScriptCapitalW] vector bosons:*)


DefineField[\[ScriptCapitalZ], Vector,
	SelfConjugate -> True,
	Mass          -> {Light,M\[ScriptCapitalZ]},
	UFO           -> <|
		"pdg"->   23, 
		"name"->  "Z", 
		"mass"->  "MZ", 
		"width"-> "WZ"
	|>
];


DefineField[\[ScriptCapitalW], Vector,
	Charges       -> U1em[1],
	Mass          -> {Light,M\[ScriptCapitalW]},
	UFO           -> <|
		"pdg"->   24, 
		"name"->  "W", 
		"mass"->  "MW", 
		"width"-> "WW"
	|>
];


(* ::Subsubsection:: *)
(*Scalars*)


(* ::Text:: *)
(*Define the physical Higgs field and the Goldstone bosons:*)


DefineField[\[ScriptH], Scalar,
	SelfConjugate -> True,
	Mass          -> {Light,M\[ScriptH]},
	UFO           -> <|
		"pdg"->   25, 
		"name"->  "h", 
		"mass"->  "Mh", 
		"width"-> "Wh"
	|>
];


DefineField[\[Chi]0, Scalar,
	Mass           -> 0,
	GoldstoneBoson -> True,
	SelfConjugate  -> True,
	NiceForm       -> "\!\(\*SuperscriptBox[\(\[Chi]\), \(0\)]\)",
	UFO            -> <|
		"pdg"->  250, 
		"name"-> "G0"
	|>
];


DefineField[\[Chi], Scalar,
	Charges        -> U1em[1],
	Mass           -> 0,
	GoldstoneBoson -> True,
	UFO            -> <|
		"pdg"->  251, 
		"name"-> "Gp"
	|>
];


(* ::Subsubsection:: *)
(*Fermions*)


(* ::Text:: *)
(*Define the broken phase vector-like fermions:*)


DefineField[\[GothicU], Fermion,
	Indices     -> {SU3c[fund], Flavor},
	Charges     -> {U1em[+2/3]},
	Mass        -> {Light, M\[GothicU], {Flavor}},
	NiceForm    -> {Default, "\!\(\*SubscriptBox[\(M\), \(\[GothicU]\)]\)"},
	UFO         -> <|
		"pdg"->   {2,    4,    6}, 
		"name"->  {"u",  "c",  "t"}, 
		"mass"->  {"Mu", "Mc", "Mt"}, 
		"width"-> {"Wu", "Wc", "Wt"}
	|>
];


DefineField[\[GothicD], Fermion,
	Indices     -> {SU3c[fund], Flavor},
	Charges     -> {U1em[-1/3]},
	Mass        -> {Light, M\[GothicD], {Flavor}},
	NiceForm    -> {Default, "\!\(\*SubscriptBox[\(M\), \(\[GothicD]\)]\)"},
	UFO         -> <|
		"pdg"->   {1,    3,    5}, 
		"name"->  {"d",  "s",  "b"}, 
		"mass"->  {"Md", "Ms", "Mb"}, 
		"width"-> {"Wd", "Ws", "Wb"}
	|>
];


DefineField[\[GothicE], Fermion,
	Indices     -> {Flavor},
	Charges     -> {U1em[-1]},
	Mass        -> {Light, M\[GothicE], {Flavor}},
	NiceForm    -> {Default, "\!\(\*SubscriptBox[\(M\), \(\[GothicE]\)]\)"},
	UFO         -> <|
		"pdg"->   {11,   13,    15}, 
		"name"->  {"e",  "mu",  "ta"}, 
		"mass"->  {"Me", "Mmu", "Mta"}, 
		"width"-> {"We", "Wmu", "Wta"}
	|>
];


DefineField[\[Nu], Fermion,
	Indices     -> {Flavor},
	Mass        -> 0,
	Chiral      -> LeftHanded,
	UFO         -> <|
		"pdg"->  {12,   14,   16}, 
		"name"-> {"ve", "vm", "vt"}
	|>
];


(* ::Subsection:: *)
(*Broken phase couplings*)


(* ::Text:: *)
(*Define the Higgs vacuum expectation value (v) and the CKM matrix as new couplings:*)


DefineCoupling[v, 
	EFTOrder      -> 1,
	SelfConjugate -> True,
	UFO           -> <|"name"->"VEV"|>
];


(* ::Text:: *)
(*Define the couplings sw and cw corresponding to the Sine and Cosine of the Weinberg angle, respectively (this facilitates automatic simplifications based on cw^2+sw^2 = 1):*)


DefineCoSinePair[cw, sw];


(* ::Text:: *)
(*Mass-basis rotation matrices for left - handed quark fields*)


ParameterDefault[QuarkAlignment-> "DownAlignment"]


Switch[QuarkAlignment
, "UnitCKM",
	Nothing
, "DownAlignment"|"UpAlignment",
	DefineCoupling[CKM,
	Indices       -> {Flavor, Flavor},
	NiceForm      -> "\!\(\*SubscriptBox[\(V\), \(CKM\)]\)",
	Unitary       -> True,
	UFO           -> <|
		"name"->{
			{"CKM_11", "CKM_12", "CKM_13"},
			{"CKM_21", "CKM_22", "CKM_23"},
			{"CKM_31", "CKM_32", "CKM_33"}
		}
	|>
	];
, "NoAlignment",
	DefineCoupling[VuL,
		Indices       -> {Flavor, Flavor},
		NiceForm      -> "\!\(\*SubscriptBox[\(V\), SubscriptBox[\(u\), \(L\)]]\)",
		Unitary       -> True,
		UFO           -> <|
			"name"->{
				{"VuL_11", "VuL_12", "VuL_13"},
				{"VuL_21", "VuL_22", "VuL_23"},
				{"VuL_31", "VuL_32", "VuL_33"}
			}
		|>
	];
	DefineCoupling[VdL,
		Indices       -> {Flavor, Flavor},
		NiceForm      -> "\!\(\*SubscriptBox[\(V\), SubscriptBox[\(d\), \(L\)]]\)",
		Unitary       -> True,
		UFO           -> <|
			"name"->{
				{"VdL_11", "VdL_12", "VdL_13"},
				{"VdL_21", "VdL_22", "VdL_23"},
				{"VdL_31", "VdL_32", "VdL_33"}
			}
		|>
	];
, _,
	Echo["Invalid value for QuarkAlignment parameter specified"]
];


(* ::Section:: *)
(*Symmetry breaking pattern*)


(* ::Subsection:: *)
(*General symmetry breaking information*)


(* ::Subsubsection:: *)
(*Group breaking pattern*)


(* ::Text:: *)
(*Define the breaking pattern of the gauge groups:*)


SetSymmetryBreakingPattern[
	{SU2L, U1Y},
	{U1em}
];


(* ::Subsubsection:: *)
(*Decomposition of representation*)


(* ::Text:: *)
(*Define how the representations decompose under this breaking pattern:*)


RepresentationDecomposition[ SU2L@adj,
	{U1em[1], U1em[-1], Singlet}
];
RepresentationDecomposition[ SU2L@fund,
	{U1em[1/2], U1em[-1/2]}
];
RepresentationDecomposition[U1Y[1], U1em[1]];


(* ::Subsubsection:: *)
(*Clebsch\[Dash]Gordan decomposition*)


(* ::Text:: *)
(*Define how the Clebsch\[Dash]Gordan coefficients decompose:*)


Module[{J, i, j}, 
	(* generators of teh fundamental representation *)
	CGDecomposition[ gen[SU2L@fund][J,i,j], {
		(* indices below correspond to {J,i,j} above *)
		{2,1,2} -> +(1/Sqrt[2])(*T+*),
		{1,2,1} -> +(1/Sqrt[2])(*T-*),
		{3,1,1} -> +(1/2)      (*T3*),
		{3,2,2} -> -(1/2)      (*T3*)
	}];
	(* epsilon tensor *)
	CGDecomposition[ eps[SU2L][i,j], {
		{1,2} -> +1,
		{2,1} -> -1
	}];
];


(* ::Text:: *)
(*The decomposition of the adjoint Clebsch-Gordan coefficients/adjoint generators is deduced automatically from the decomposition of the fundamental.*)


(* ::Subsection:: *)
(*Decomposition of fields*)


(* ::Subsubsection:: *)
(*Vector decomposition*)


(* ::Text:: *)
(*Define how the weak gauge bosons decompose:*)


Module[{\[Mu], a},
	FieldDecomposition[ W[\[Mu],a], {
		\[ScriptE][]/sw[]\[ScriptCapitalW][\[Mu]],              (* U1em[1] *)
		\[ScriptE][]/sw[]Bar@\[ScriptCapitalW][\[Mu]],          (* U1em[-1] *)
		(cw[]/sw[]\[ScriptE][]\[ScriptCapitalZ][\[Mu]] + \[ScriptCapitalA][\[Mu]]) (* Singlet *)
	}];	
	FieldDecomposition[ B[\[Mu]],
		(-sw[]/cw[]\[ScriptE][]\[ScriptCapitalZ][\[Mu]] + \[ScriptCapitalA][\[Mu]]) (* Singlet *)
	];
];


(* ::Subsubsection:: *)
(*Scalar decomposition*)


(* ::Text:: *)
(*Decompose the Higgs into spin up and down:*)


Module[{i},
	FieldDecomposition[ H[i], {
		-I*\[Chi][],
		(v[] + \[ScriptH][] + I*\[Chi]0[])/Sqrt[2]
	}];
];


(* ::Subsubsection:: *)
(*Fermion decomposition*)


(* ::Text:: *)
(*Decompose the doublet fermions in spin up and down, and the trivial decomposition for the singlets:*)


Module[{a, i, p},
	FieldDecomposition[ l[i, p], {
		\[Nu][p],
		\[GothicE][p]
	}];
	FieldDecomposition[ u[a, p], \[GothicU][a, p] ];
	FieldDecomposition[ d[a, p], \[GothicD][a, p] ];
	FieldDecomposition[ e[p], \[GothicE][p] ];
	(*The decomposition of the left-handed quarks is determined by the QuarkAlignment*)
	Switch[QuarkAlignment
	, "UnitCKM",
		FieldDecomposition[ q[a, i, p], {
			\[GothicU][a, p],
			\[GothicD][a, p]
		}];
	, "DownAlignment",
		FieldDecomposition[ q[a, i, p], {
			Bar@CKM[r, p]\[GothicU][a, r],
			\[GothicD][a, p]
		}];
	, "UpAlignment",
		FieldDecomposition[ q[a, i, p], {
			\[GothicU][a, p],
			CKM[p, r]\[GothicD][a, r]
		}];
	, "NoAlignment",	
		FieldDecomposition[ q[a, i, p], {
			VuL[p, r]\[GothicU][a, r],
			VdL[p, r]\[GothicD][a, r]
		}];
	];
];


(* ::Subsection:: *)
(*Broken phase coupling replacements*)


(* ::Text:: *)
(*Substitution rules for the SM couplings in terms of convenient broken phase parameters. The QuarkAlignment variable specifies if the Yukawa coupling matrices are aligned to either up or down quarks.*)


Module[{yukawas, s, t},
	(*Replace Yukawa couplings in favor of fermion masses and alignment matrices*)
	yukawas= Switch[QuarkAlignment
		, "UnitCKM", {
			Yu[s_, t_]-> (Sqrt[2]/v[]) Delta[Flavor][s, t]M\[GothicU][t],
			Yd[s_, t_]-> (Sqrt[2]/v[]) Delta[Flavor][s, t]M\[GothicD][t]
		}, "DownAlignment", {
			Yu[s_, t_]-> (Sqrt[2]/v[]) Bar@ CKM[t, s]M\[GothicU][t],
			Yd[s_, t_]-> (Sqrt[2]/v[]) Delta[Flavor][s, t]M\[GothicD][t]
		}, "UpAlignment", {
			Yu[s_, t_]-> (Sqrt[2]/v[]) Delta[Flavor][s, t] M\[GothicU][t],
			Yd[s_, t_]-> (Sqrt[2]/v[]) CKM[s, t] M\[GothicD][t]
		}, "NoAlignment", {
			Yu[s_, t_]-> (Sqrt[2]/v[]) VuL[s, t] M\[GothicU][t],
			Yd[s_, t_]-> (Sqrt[2]/v[]) VdL[s, t] M\[GothicD][t]
		}];
		
	SetSSBReplacements[Join[yukawas, {
		(*Repalce lepton Yukawa*)
		Ye[s_, t_]-> (Sqrt[2]/v[]) Delta[Flavor][s, t] M\[GothicE][t],
		(*Replace electroweak gauge couplings *)
		gL[]-> \[ScriptE][]/sw[],
		gY[]-> \[ScriptE][]/cw[]
	}]];
];


(* ::Section:: *)
(*Return*)


(* ::Text:: *)
(*The "0" in the last line ensures that only the parent Lagrangian in the unbroken phase is returned.*)
(*The broken phase Lagrangian can be obtained using: *)
(*ImplementVacuumConditions@ToBrokenPhase@LoadModel["SM+breaking"]*)


0
