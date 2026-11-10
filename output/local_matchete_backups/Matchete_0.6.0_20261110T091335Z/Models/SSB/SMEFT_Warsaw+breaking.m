(* ::Package:: *)

(* ::Title:: *)
(*SMEFT model file with electroweak symmetry breaking*)


(* ::Section:: *)
(*Loading the unbroken phase Lagrangian*)


ParentModel["SMEFT_Warsaw"]
ParentModel["SM+breaking", "UseLagrangian"-> False]


(* ::Text:: *)
(*Turn off B- and L-violating operators by default:*)


ParameterDefault[BLViolation-> False];


(* ::Section:: *)
(*Broken phase definitions*)


(* ::Subsection:: *)
(*Additional couplings*)


(* ::Text:: *)
(*SMEFT VEV*)


DefineCoupling[vT, 
	EFTOrder      -> 1,
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubscriptBox[\(v\), \(T\)]\)",
	UFO           -> <|"name"->"VEVT"|>
];


(* ::Text:: *)
(*Shifted strong gauge coupling*)


DefineCoupling[\[ScriptG]s, 
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubscriptBox[\(\[ScriptG]\), \(s\)]\)",
	UFO           -> <|"name"-> "gs"|>
];


(* ::Subsection:: *)
(*Overwriting the ordinary SM field decomposition*)


(* ::Subsubsection:: *)
(*Vector decomposition*)


(* ::Text:: *)
(*Define how the weak gauge bosons decompose:*)


Module[{\[Mu], a, \[Epsilon]WB = -(1/2) vT[]^2 cHWB[]},
	FieldDecomposition[ W[\[Mu],a],
		{
			(*This is equivalen to gL[] \[ScriptCapitalW][\[Mu]] and expressing gL in terms of \[ScriptE] *)
			\[ScriptE][]/sw[] (1 - \[Epsilon]WB cw[]/sw[]) \[ScriptCapitalW][\[Mu]],         (* U1em[1] *)
			\[ScriptE][]/sw[] (1 - \[Epsilon]WB cw[]/sw[]) Bar@ \[ScriptCapitalW][\[Mu]],     (* U1em[-1] *)
			(\[ScriptE][]/sw[] (cw[] - \[Epsilon]WB/sw[]) \[ScriptCapitalZ][\[Mu]] + \[ScriptCapitalA][\[Mu]]) (* Singlet *)
		}
	];
	
	FieldDecomposition[ B[\[Mu]],
		(-(\[ScriptE][]/cw[])(sw[] - \[Epsilon]WB/cw[]) \[ScriptCapitalZ][\[Mu]] + \[ScriptCapitalA][\[Mu]]) (* Singlet *)
	];
];


(* ::Subsubsection:: *)
(*Scalar decomposition*)


(* ::Text:: *)
(*Decompose the Higgs into spin up and down:*)


Module[
	{
		i,
		cHkin=  (cHBox[]-cHD[]/4)vT[]^2,
		c\[Chi]0kin= -(cHD[]/4) vT[]^2
	}
,
	FieldDecomposition[ H[i],
		{
			-I*\[Chi][],
			(vT[] + (1+cHkin)\[ScriptH][] + I*(1+c\[Chi]0kin)*\[Chi]0[])/Sqrt[2]
		}
	]
];


(* ::Subsection:: *)
(*New couplings*)


(* ::Text:: *)
(*Substitution rules for the SM couplings in terms of convenient broken phase parameters. The QuarkAlignment variable specifies if the Yukawa coupling matrices are aligned to either up or down quarks. *)
(*N.b. the default value of QuarkAlignment is inherited from "SM+breaking".*)


Module[{yukawas, s, t, \[Epsilon]WB = -(1/2) vT[]^2 cHWB[]},
	(*Replace Yukawa couplings in favor of fermion masses and alignment matrices*)
	yukawas= Switch[QuarkAlignment
		, "UnitCKM", {
			Yu[s_,t_]-> (Sqrt[2]/vT[]) Delta[Flavor][s, t]M\[GothicU][t] + 1/2 vT[]^2 cuH[s, t],
			Yd[s_,t_]-> (Sqrt[2]/vT[]) Delta[Flavor][s, t] M\[GothicD][t] + 1/2 vT[]^2 cdH[s, t]
		}, "DownAlignment", {
			Yu[s_,t_]-> (Sqrt[2]/vT[]) Bar@CKM[t, s]M\[GothicU][t] + 1/2 vT[]^2 cuH[s, t],
			Yd[s_,t_]-> (Sqrt[2]/vT[]) Delta[Flavor][s, t] M\[GothicD][t] + 1/2 vT[]^2 cdH[s, t]
		}, "UpAlignment", {
			Yu[s_,t_]-> (Sqrt[2]/vT[]) Delta[Flavor][s, t] M\[GothicU][t] + 1/2 vT[]^2 cuH[s, t],
			Yd[s_,t_]-> (Sqrt[2]/vT[]) CKM[s, t]M\[GothicD][t] + 1/2 vT[]^2 cdH[s, t]
		}, "NoAlignment", {
			Yu[s_, t_]-> (Sqrt[2]/vT[]) VuL[s, t]M\[GothicU][t] + 1/2 vT[]^2 cuH[s, t],
			Yd[s_, t_]-> (Sqrt[2]/vT[]) VdL[s, t]M\[GothicD][t] + 1/2 vT[]^2 cdH[s, t]
		}];
		
	SetSSBReplacements[Join[yukawas, {
		(*Replace lepton Yukawa*)
		Ye[s_,t_]-> (Sqrt[2]/vT[]) Delta[Flavor][s, t] M\[GothicE][t] + 1/2 vT[]^2 ceH[s, t],
		(*Replace electroweak gauge couplings*)
		gL[]-> \[ScriptE][]/sw[] (1 -vT[]^2 cHW[] -\[Epsilon]WB cw[]/sw[]),
		gY[]-> \[ScriptE][]/cw[] (1 -vT[]^2 cHB[] -\[Epsilon]WB sw[]/cw[]),
		gs[]-> \[ScriptG]s[] (1 -vT[]^2 cHG[])
	}]];
];


(* ::Section:: *)
(*Return*)


(* ::Text:: *)
(*The "0" in the last line ensures that only the parent Lagrangian in the unbroken phase is returned.*)
(*The broken phase Lagrangian can be obtained using: *)
(*ImplementVacuumConditions@ToBrokenPhase@LoadModel["SMEFT_Warsaw+breaking"]*)


0
