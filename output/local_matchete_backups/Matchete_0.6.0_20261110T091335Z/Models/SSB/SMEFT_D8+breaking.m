(* ::Package:: *)

(* ::Title:: *)
(*SMEFT model file with electroweak symmetry breaking*)


(* ::Section:: *)
(*Load parent model*)


ParentModel["SMEFT_D8"]
ParentModel["SMEFT_Warsaw+breaking", "UseLagrangian"-> False]


ParameterDefault[BLViolation-> False];


(* ::Section:: *)
(*Broken phase definitions*)


(* ::Subsection:: *)
(*Broken phase couplings*)


(* ::Subsubsection:: *)
(*Shifted gauge couplings*)


(* ::Text:: *)
(*\[ScriptG]s is defined in SMWEFT_Warsaw*)


DefineCoupling[\[ScriptG]L, 
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubscriptBox[\(\[ScriptG]\), \(L\)]\)",
	UFO           -> <|"name"-> "gL"|>
];

DefineCoupling[\[ScriptG]Y, 
	SelfConjugate -> True,
	NiceForm      -> "\!\(\*SubscriptBox[\(\[ScriptG]\), \(Y\)]\)",
	UFO           -> <|"name"-> "gY"|>
];


(* ::Subsection:: *)
(*Overwriting the ordinary SM field decomposition*)


(* ::Subsubsection:: *)
(*Vector decomposition*)


(* ::Text:: *)
(*Define how the weak gauge bosons decompose:*)


Module[{x, X, alpha, beta, gamma, delta, W3, \[ScriptCapitalB], Wp, Wm},
	
	alpha = 1 - 2 * x^2 * cHW[] * vT[]^2 - x^4 * cW2H4$1[] * vT[]^4;
    beta = 1 - 2 * x^2 * cHW[] * vT[]^2 - x^4 * cW2H4$1[] * vT[]^4 - x^4 * cW2H4$3[] * vT[]^4;
    gamma = 1 - 2 * x^2 * cHB[] * vT[]^2 - x^4 * cB2H4$1[] * vT[]^4;
    delta = x^2 * cHWB[] * vT[]^2 + 1/2 * x^4 * cWBH4$1[] * vT[]^4;

    (* Transformation matrix *)
    X = 1/2 * {{Sqrt[alpha]/Sqrt[beta], 0}, {0, 1}} . {{1, -1}, {1, 1}} .
        {{1/Sqrt[1 + delta/Sqrt[beta*gamma]], 0},{0, 1/Sqrt[1 - delta/Sqrt[beta*gamma]]}} .
        {{1, 1}, {-1, 1}} .
        {{cw[], sw[]}, {-sw[], cw[]}};

    X = Normal[Series[X, {x, 0, 4}]] /. x->1;

    {W3, \[ScriptCapitalB]} = {{\[ScriptG]L[], 0},{0, \[ScriptG]Y[]}} . X . {\[ScriptCapitalZ][\[Mu]], \[ScriptCapitalA][\[Mu]] / \[ScriptE][]};
    {Wp, Wm} = \[ScriptG]L[] * {\[ScriptCapitalW][\[Mu]], Bar@\[ScriptCapitalW][\[Mu]]};


	FieldDecomposition[ W[\[Mu],a],
		{
			Wp, (* U1em[1] *)
			Wm, (* U1em[-1] *)
			W3  (* Singlet *)
		}
	];
	
	FieldDecomposition[ B[\[Mu]],
		\[ScriptCapitalB]       (* Singlet *)
	];
];



(* ::Subsubsection:: *)
(*Scalar decomposition*)


(* ::Text:: *)
(*Decompose the Higgs into spin up and down:*)


Module[
	{
		cHkin=  cHBox[] * vT[]^2 - cHD[] * vT[]^2 / 4 - cH6$1[] * vT[]^4 / 8 - cH6$2[] * vT[]^4 / 8 + 3/2 * cHBox[]^2 * vT[]^4 + 3/32 * cHD[]^2 * vT[]^4 - 3 * cHBox[] * cHD[] * vT[]^4 / 4,
		c\[Chi]0kin= -1/4 cHD[] vT[]^2 - 1/8 cH6$1[] vT[]^4 - 1/8 cH6$2[] vT[]^4 + 3/32 cHD[]^2 vT[]^4,
		c\[Chi]kin= -1/8 cH6$1[] vT[]^4 + 1/8 cH6$2[] vT[]^4
	},
	FieldDecomposition[ H[i],
		{
			-I*(1+c\[Chi]kin)*\[Chi][],
			(vT[] + (1+cHkin)\[ScriptH][] + I*(1+c\[Chi]0kin)*\[Chi]0[])/Sqrt[2]
		}
	];
];


(* ::Subsubsection:: *)
(*Parameter substitutions*)


(* ::Text:: *)
(*Substitution rules for the SM couplings in terms of convenient broken phase parameters*)


Module[
    {	
		yukawas,
        p,r,s,t,
        alpha = 1 - 2 * cHW[] * vT[]^2 - cW2H4$1[] * vT[]^4,
        gamma = 1 - 2 * cHB[] * vT[]^2 - cB2H4$1[] * vT[]^4
    },
	yukawas= Switch[QuarkAlignment
		, "UnitCKM", {
			Yu[s_,t_]-> Sqrt[2] / vT[] * M\[GothicU][t] * Delta[Flavor][s,t] + cuH[s,t] / 2 * vT[]^2 + cquH5[s,t] / 4 * vT[]^4,
			Yd[s_,t_]-> Sqrt[2] / vT[] * M\[GothicD][t] * Delta[Flavor][s,t] + cdH[s,t] / 2 * vT[]^2 + cqdH5[s,t] / 4 * vT[]^4
		}, "DownAlignment", {
			Yu[s_,t_]-> Sqrt[2] / vT[] * M\[GothicU][t] * Bar@CKM[t,s] + cuH[s,t] / 2 * vT[]^2 + cquH5[s,t] / 4 * vT[]^4,
			Yd[s_,t_]-> Sqrt[2] / vT[] * M\[GothicD][t] * Delta[Flavor][s,t] + cdH[s,t] / 2 * vT[]^2 + cqdH5[s,t] / 4 * vT[]^4
		}, "UpAlignment", {
			Yu[s_,t_]-> Sqrt[2] / vT[] * M\[GothicU][t] * Delta[Flavor][s,t] + cuH[s,t] / 2 * vT[]^2 + cquH5[s,t] / 4 * vT[]^4,
			Yd[s_,t_]-> Sqrt[2] / vT[] * M\[GothicD][t] * CKM[s,t] + cdH[s,t] / 2 * vT[]^2 + cqdH5[s,t] / 4 * vT[]^4
		}, "NoAlignment", {
			Yu[s_, t_]-> Sqrt[2] / vT[] * M\[GothicU][t] * VuL[s,t] + cuH[s,t] / 2 * vT[]^2 + cquH5[s,t] / 4 * vT[]^4,
			Yd[s_, t_]-> Sqrt[2] / vT[] * M\[GothicD][t] * VdL[s,t] + cdH[s,t] / 2 * vT[]^2 + cqdH5[s,t] / 4 * vT[]^4
		}];

	SetSSBReplacements[Join[yukawas, {
		(*Replace lepton Yukawa*)
		Ye[s_,t_] :> Sqrt[2] / vT[] * M\[GothicE][t] * Delta[Flavor][s,t] + ceH[s,t] / 2 * vT[]^2 + cleH5[s,t] / 4 * vT[]^4,
		
		(* replace gauge couplings *)
		gs[] -> \[ScriptG]s[] (1 - cHG[] vT[]^2 - 1/2 cG2H4$1[] vT[]^4 - 1/2 cHG[]^2 vT[]^4),
		gL[] -> \[ScriptG]L[] Sqrt[alpha],
		gY[] -> \[ScriptG]Y[] Sqrt[gamma]
		}]
	];
];


(* ::Section:: *)
(*Return*)


(* ::Text:: *)
(*The "0" in the last line ensures that only the parent Lagrangian in the unbroken phase is returned.*)
(*The broken phase Lagrangian can be obtained using: *)
(*ImplementVacuumConditions@ToBrokenPhase@LoadModel["SMEFT_D8+breaking"]*)


0
