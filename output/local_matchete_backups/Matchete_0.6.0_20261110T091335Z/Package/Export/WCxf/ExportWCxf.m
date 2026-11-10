(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Export to WCxf files*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsection:: *)
(*Exported*)


PackageExport["ExportWCxf"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsection:: *)
(*Exported*)


ExportWCxf::usage = "ExportWCxf[arg, input, output] generates WCxf files for the specified BSM parameter points.
- arg: either a SMEFT Lagrangian obtained, e.g., by the Match function, or directly a list of matching conditions obtained by MapEffectiveCouplings[\[ScriptCapitalL]EFT, LoadModel[\"SMEFT_Warsaw\"], ShiftRenCouplings->True]. Notice that only BSM models matching directly onto the default SMEFT Warsaw basis Lagrangian used by Matchete are currently supported and the option ShiftRenCouplings hast to be set to True for the matching;
- input: path to a JSON input file containing the following Keys:
	- \"matching scale\" : matching scale in GeV;
	- \"parameters\" : list of all parameter names for which numerical input is provided;
	- \"dictionary\" : dictionary mapping from Matchete notation to the parameter names given in \"parameters\";
	- \"parameter points\" : dictionary containing all numerical values, Keys are integer IDs for the points, and values are lists of numbers for the parameters specified in \"parameters\" in the same order;
	- \"SM\" : optional dictionary for overwriting default SM values or scheme choices.
- output: path to the directory where all WCxf files should be stored."


(* ::Chapter:: *)
(*Determine numeric values of all SM and BSM parameters*)


(* ::Section:: *)
(*Default SM parameters [taken from the PDG]*)


(* ::Subsection:: *)
(*SM input at matching scale*)


(* ::Text:: *)
(*Scaleful quantities in units of GeV*)


Options@ParametersSM = {
	(* scheme choices *)
	"EW scheme"    -> 1, (* allowed values are: 1, 2, 3, 4, 5 *)
	"flavor basis" -> "down", (* flavor alignement in quark sector; allowed values are: "down", "up" *)
	(* numerical input at Z pole; only a subset of these will be used depedning on specified scheme *)
	"mH"     -> 125.20,             (* physical Higgs boson mass *)
	"GFermi" -> 1.1663788*10^-5,    (* Fermi's constant *)
	"mZ"     -> 91.1880,            (* Z-boson mass *)
	"mW"     -> 80.3692,            (* W-boson mass *)
	"sW2"    -> 0.22305,            (* Sine of weak mixing angle *)
	"\[Alpha]EM"    -> 7.2973525693*10^-3, (* electromagnetic fine-structure constant *)
	"\[Alpha]S"     -> 0.1180,             (* strong fine-structure constant *)
	(* fermion masses *)
	"me"     -> 0.511*10^-3,        (* electron mass *)
	"m\[Mu]"     -> 105.66*10^-3,       (* muon mass *)
	"m\[Tau]"     -> 1.77693,            (* tau mass *)
	"md"     -> 4.70*10^-3,         (* down-quark mass *)
	"ms"     -> 93.5*10^-3,         (* strange-quark mass *)
	"mb"     -> 4.183,              (* bottom-quark mass *)
	"mu"     -> 2.16*10^-3,         (* up-quark mass *)
	"mc"     -> 1.2730,             (* charm-quark mass *)
	"mt"     -> 162.5,              (* top-quark mass (MSbar value) *)
	(* CKM Wolfenstein parameters *)
	"\[Lambda]c"     -> 0.22501,
	"A"      -> 0.826,
	"\[Rho]Bar"   -> 0.1591,
	"\[Eta]Bar"   -> 0.3523
}


ParametersSM::UnknownScheme = "The specified EW input scheme (\"EW scheme\" -> `1`) is unknown. Currently allowed values are 1, 2, 3, 4, 5."
ParametersSM::alignment     = "The specified quark flavor alignment (\"flavor basis\" -> `1`) is unknown. Currently allowed choices are \"down\" or \"up\"."


ParametersSM[scale_, OptionsPattern[]] := Module[
	{
		mH, GFermi, mZ, mW, sW, \[Alpha]EM, vev, g2mZ, g1mZ, ge, \[Alpha]SmZ, g3mZ,
		\[Beta]0g3, \[Beta]0g2, \[Beta]0g1, g3Sol, g3func,
		mass, yTop, yTopLow, yTopSol,
		\[Lambda]c, A, \[Rho], \[Eta], WolfensteinReplacements, VCKM, YukawaReplace
	}
	,
	(* Electroweak & Higgs sector *)
	(* various possible inputs *)
	mH     = OptionValue["mH"];      (* physical Higgs boson mass *)
	GFermi = OptionValue["GFermi"];  (* Fermi's constant *)
	mZ     = OptionValue["mZ"];      (* Z-boson mass *)
	mW     = OptionValue["mW"];      (* W-boson mass *)
	sW     = Sqrt[OptionValue["sW2"]];   (* Sine of weak mixing angle *)
	\[Alpha]EM    = OptionValue["\[Alpha]EM"];     (* electromagnetic fine-structure constant *)
	
	(* required output: \[Mu]EW, \[Lambda]H, g1, g2 *)
	(* 
	NOTE: currently no SMEFT corrections are implemented; for details see 
	 - https://arxiv.org/abs/2111.12515 
	 - https://arxiv.org/abs/2305.03763
	 - https://arxiv.org/abs/2312.08446 
	*)
	Switch[OptionValue["EW scheme"],
		1, (* {Subscript[G, F],Subscript[m, Z],Subscript[m, W]} *)
			vev  = (Sqrt[2]*GFermi)^(-1/2);    (* electroweak vacuum expectation value *)
			g2mZ = 2*mW/vev;                   (* SU(2)L gauge coupling *)
			g1mZ = Sqrt[(2*mZ/vev)^2-g2mZ^2],  (* U(1)Y gauge coupling *)
		2, (* {Subscript[G, F],Subscript[m, Z],Subscript[\[Alpha], em]} *)
			vev  = (Sqrt[2]*GFermi)^(-1/2);                         (* electroweak vacuum expectation value *)
			ge   = Sqrt[4\[Pi]*\[Alpha]EM];                                    (* QED gauge coupling *)
			g1mZ = Sqrt[2*mZ^2/vev^2 (1-Sqrt[1-ge^2*vev^2/mZ^2])]; (* U(1)Y gauge coupling *)
			g2mZ = 2 mZ/vev ge/g1mZ,                                (* SU(2)L gauge coupling *)
		3, (* {Subscript[\[Alpha], em],Subscript[m, Z],Subscript[m, W]} *)
			ge   = Sqrt[4\[Pi]*\[Alpha]EM];              (* QED gauge coupling *)
			g2mZ = ge/Sqrt[1-(mW/mZ)^2];      (* SU(2)L gauge coupling *)
			vev  = 2*mW/g2mZ;                 (* electroweak vacuum expectation value *)
			g1mZ = Sqrt[4*(mZ/vev)^2-g2mZ^2], (* U(1)Y gauge coupling *)
		4, (* {Subscript[G, F],Subscript[s, W],Subscript[m, Z]} *)
			vev = (Sqrt[2]*GFermi)^(-1/2);  (* electroweak vacuum expectation value *)
			g1mZ = 2*sW*mZ/vev;             (* U(1)Y gauge coupling *)
			g2mZ = 2*(mZ/vev)Sqrt[1-sW^2],  (* SU(2)L gauge coupling *)
		5, (* {Subscript[s, W],Subscript[m, Z],Subscript[\[Alpha], em]} *)
			ge   = Sqrt[4\[Pi]*\[Alpha]EM];             (* QED gauge coupling *)
			g2mZ = ge/sW;                    (* SU(2)L gauge coupling *)
			g1mZ = ge/Sqrt[1-sW^2];          (* U(1)Y gauge coupling *)
			vev  = 2*mZ/Sqrt[g1mZ^2+g2mZ^2], (* electroweak vacuum expectation value *)
		_,
			Message[ParametersSM::UnknownScheme, OptionValue["EW input scheme"]];
			Abort[]
	];
	(* scheme independent formulae *)
	\[Mu]EW = Sqrt[1/2 mH^2];  (* Higgs doublet mass term *)
	\[Lambda]H  = mH^2/vev^2;      (* quartic Higgs doublet coupling *)
	
	(* QCD *)
	\[Alpha]SmZ = OptionValue["\[Alpha]S"]; (* strong fine structure constant *)
	g3mZ = Sqrt[4\[Pi]*\[Alpha]SmZ];     (* SU(3)c gauge coupling *)
	
	(* RGE evolution to matching scale *)
	(* only consider pure SM one loop running *)
	(* beta functions \[Beta]0 from [1704.04504] *)
	\[Beta]0g3 = -(7/(16\[Pi]^2));         (* SU(3)c *)
	\[Beta]0g2 = 1/(16\[Pi]^2)*(-(19/6));  (* SU(2)L *)
	\[Beta]0g1 = 1/(16\[Pi]^2)*(41/6);     (* U(1)Y *)
	(* gauge couplings at matching scale *)
	{g3, {g3func,g3Sol}} = RGEvolve[\[Beta]0g3, mZ, scale, g3mZ, "solution" -> True];
	g2 = RGEvolve[\[Beta]0g2, mZ, scale, g2mZ];
	g1 = RGEvolve[\[Beta]0g1, mZ, scale, g1mZ];
	
	(* fermion masses *)
	mass["e"] = OptionValue["me"];
	mass["\[Mu]"] = OptionValue["m\[Mu]"];
	mass["\[Tau]"] = OptionValue["m\[Tau]"];
	mass["d"] = OptionValue["md"];
	mass["s"] = OptionValue["ms"];
	mass["b"] = OptionValue["mb"];
	mass["u"] = OptionValue["mu"];
	mass["c"] = OptionValue["mc"];
	mass["t"] = OptionValue["mt"];
	
	(* RG evolution of the top mass/Yukawa *)
	(* low scale input *)
	yTopLow = Sqrt[2]/vev*mass["t"];
	(* 1-loop SM running to matching scale *)
	yTopSol = First@NDSolve[
		{
			(16\[Pi]^2)*\[Mu]*yTop'[\[Mu]] == (9/2)*yTop[\[Mu]]^3-8(g3func/.g3Sol)[\[Mu]]^2 yTop[\[Mu]],
			yTop[mass["t"]] == yTopLow
		},
		yTop,
		{\[Mu], mZ-5, scale+50}
	];
	(* matching scale output *)
	mass["t"] = vev/Sqrt[2] (yTop[scale]/.yTopSol);
	
	(* Wolfenstein parametrization *)
	WolfensteinReplacements = {
		\[Lambda]c -> OptionValue["\[Lambda]c"],
		A  -> OptionValue["A"],
		\[Rho]  -> OptionValue["\[Rho]Bar"] (1 + OptionValue["\[Lambda]c"]^2/2), (* \[Rho]Bar=\[Rho](1-\[Lambda]c^2/2) *)
		\[Eta]  -> OptionValue["\[Eta]Bar"] (1 + OptionValue["\[Lambda]c"]^2/2)  (* \[Eta]Bar=\[Eta](1-\[Lambda]c^2/2) *)
	};
	(* CKM matrix *)
	VCKM = {
		{ 1-\[Lambda]c^2/2,           \[Lambda]c,         A*\[Lambda]c^3*(\[Rho]-I*\[Eta]) },
		{ -\[Lambda]c,                1-\[Lambda]c^2/2,   A*\[Lambda]c^2         },
		{ A*\[Lambda]c^3*(1-\[Rho]-I*\[Eta]),   -A*\[Lambda]c^2,    1              }
	};
	VCKM = VCKM //. WolfensteinReplacements;
	(* 
	NOTE: noe CKM input schemes are currently implemented; for details see:
	 - https://arxiv.org/abs/1812.08163
	*)
	
	(* Determine SMEFT Yukawas in up or down basis *)
	(* NOTE: the corrections of d=6 Wilson coefficients proportional to cuH, cdH, and ceH are effectively d=8 contributions. Thus they can and must be neglected. *)
	Switch[OptionValue["flavor basis"],
		"down",
			YukawaReplace[] := With[
				{
					$Yu = Sqrt[2]/vev (VCKM\[ConjugateTranspose] . DiagonalMatrix[{mass["u"],mass["c"],mass["t"]}]) (*+ vev^2/2 Table[Coupling[Global`cuH,{i,j},0],{i,1,3},{j,1,3}]*)
				}
				,
				{
					Coupling[Global`celH,{1,1},0]            -> Sqrt[2]/vev mass["e"] (*+ vev^2/2 Coupling[Global`ceH,{1,1},0]*),
					Coupling[Global`celH,{2,2},0]            -> Sqrt[2]/vev mass["\[Mu]"] (*+ vev^2/2 Coupling[Global`ceH,{2,2},0]*),
					Coupling[Global`celH,{3,3},0]            -> Sqrt[2]/vev mass["\[Tau]"] (*+ vev^2/2 Coupling[Global`ceH,{3,3},0]*),
					Coupling[Global`celH,{p_,r_},0]/;(p=!=r) -> 0 (*+ vev^2/2 Coupling[Global`ceH,{p,r},0]*)
					,
					Coupling[Global`cdqH,{1,1},0]            -> Sqrt[2]/vev mass["d"] (*+ vev^2/2 Coupling[Global`cdH,{1,1},0]*),
					Coupling[Global`cdqH,{2,2},0]            -> Sqrt[2]/vev mass["s"] (*+ vev^2/2 Coupling[Global`cdH,{2,2},0]*),
					Coupling[Global`cdqH,{3,3},0]            -> Sqrt[2]/vev mass["b"] (*+ vev^2/2 Coupling[Global`cdH,{3,3},0]*),
					Coupling[Global`cdqH,{p_,r_},0]/;(p=!=r) -> 0 (*+ vev^2/2 Coupling[Global`cdH,{p,r},0]*)
					,
					Coupling[Global`cquH,{p_,r_},0]          :> $Yu[[p,r]]
				}
			]
			,
		"up",
			YukawaReplace[] := With[
				{
					$Yd = Sqrt[2]/vev (DiagonalMatrix[{mass["d"],mass["s"],mass["b"]}] . VCKM (*+ vev^2/2Table[Coupling[Global`cdH,{i,j},0],{i,1,3},{j,1,3}]*))
				}
				,
				{
					Coupling[Global`celH,{1,1},0]            -> Sqrt[2]/vev mass["e"] (*+ vev^2/2Coupling[Global`ceH,{1,1},0]*),
					Coupling[Global`celH,{2,2},0]            -> Sqrt[2]/vev mass["\[Mu]"] (*+ vev^2/2Coupling[Global`ceH,{2,2},0]*),
					Coupling[Global`celH,{3,3},0]            -> Sqrt[2]/vev mass["\[Tau]"] (*+ vev^2/2Coupling[Global`ceH,{3,3},0]*),
					Coupling[Global`celH,{p_,r_},0]/;(p=!=r) -> 0 (*+ (vev^2/2)Coupling[ceH,{p,r},0]*)
					,
					Coupling[Global`cquH,{1,1},0]            -> Sqrt[2]/vev mass["u"] (*+ vev^2/2Coupling[Global`cuH,{1,1},0]*),
					Coupling[Global`cquH,{2,2},0]            -> Sqrt[2]/vev mass["c"] (*+ vev^2/2Coupling[Global`cuH,{2,2},0]*),
					Coupling[Global`cquH,{3,3},0]            -> Sqrt[2]/vev mass["t"] (*+ vev^2/2Coupling[Global`cuH,{3,3},0]*),
					Coupling[Global`cquH,{p_,r_},0]/;(p=!=r) -> 0 (*+ (vev^2/2)Coupling[Global`cuH,{p,r},0]*)
					,
					Coupling[Global`cdqH,{p_,r_},0]          :> $Yd[[p,r]]
				}
			]
			,
		_,
			Message[ParametersSM::alignment, OptionValue["flavor basis"]];
			Abort[]	
	];
	
	Return@{
		\[Mu]EW, \[Lambda]H, g3, g2, g1, YukawaReplace
	}
]


(* ::Subsection:: *)
(*Function for RGE evolution*)


(* ::Text:: *)
(*RG evolution at one-loop with beta function \[Beta]0, from scale0 to scale1 with initial condition value0.*)


Options@RGEvolve = {"solution" -> False}


RGEvolve[\[Beta]0_, scale0_, scale1_, value0_, OptionsPattern[]] := Module[
	{g,gSol}
	,
	gSol = First@NDSolve[
		{
			\[Mu]*g'[\[Mu]]   == \[Beta]0*g[\[Mu]]^3, 
			g[scale0] == value0
		}, 
		g, 
		{\[Mu],scale0-10,scale1+100}
	];
	
	If[OptionValue["solution"],
		{g[scale1]/.gSol, {g,gSol}}
		,
		g[scale1]/.gSol
	]
]


(* ::Section:: *)
(*(B)SM replacement rules*)


(* ::Subsection:: *)
(*SM*)


(* same options as ParametersSM *)
Options@NumericReplacementsSM = Options@ParametersSM;


NumericReplacementsSM::usage = "NumericReplacementsSM[] returns the numerical replacment rule for all SM parameteres."


NumericReplacementsSM[matchingScale_, opt:OptionsPattern[]] := Module[{\[Mu]EW, \[Lambda]H, g3, g2, g1, VCKM, YukawaReplace},
	{\[Mu]EW, \[Lambda]H, g3, g2, g1, YukawaReplace} = ParametersSM[matchingScale, opt];
	Join[
		{
			ToExpression["hbar"]                                  -> 1/(16\[Pi]^2),
			Coupling[ToExpression["cB2"]|ToExpression["gY"],{},0] -> g1,
			Coupling[ToExpression["cW2"]|ToExpression["gL"],{},0] -> g2,
			Coupling[ToExpression["cG2"]|ToExpression["gs"],{},0] -> g3,
			Coupling[ToExpression["cH2"],{},0]                    -> \[Mu]EW^2,
			Coupling[ToExpression["cH4"],{},0]                    -> \[Lambda]H (* ?WRONG?: there is an offset between cH4 and \[Lambda]H due to EOM redundant operators *)
		},
		YukawaReplace[]
	]
]


(* ::Subsection:: *)
(*BSM*)


NumericReplacementsBSM::usage = "NumericReplacementsBSM[params, vals, dict] Returns a list of numerical replacement rules for the BSM parameters (specified in the List params) with one entry per paramerter point consisting of its own table of replacment rules for this point.
The numerical values for all paramerter points are taken from vals.
The argument dictionary is a list of replacement rules mapping from the Matchete convention to the specified parameters."


ExportWCxf::unknownInput = "The numerical input for the BSM paramerter points could not be determined. Allowed values are either an association or the path to a file containing all pramerter points."


NumericReplacementsBSM[params_, numerics_, dict_] := Module[{repl},
	(* create replacement rule for each parameter point *)
	repl = Table[
		Table[
			params[[j]] -> Values[numerics][[i,j]]
			,
			{j,Length@params} (* loop over keys for each paramerter points *)
		]
		,
		{i, Length@numerics} (* loop over number or paramerter points *)
	];

	Normal[dict /. repl]
]


(* ::Section:: *)
(*Main function for WCxf generation*)


ExportWCxf::blviolation = "WARNING: WCxf export is currently not supported for Baryon and/or Lepton number violating operators. Non-vanishing matching condition for the Wilson coefficient `1` detected, but ignored for WCxf export."


ExportWCxf::defaultSMEFT = "ExportWCxf currently only supports the default SMEFT Lagrangian and parameter name changes are allowed. Please load and use the default SMEFT definitions with LoadModel[\"SMEFT_Warsaw\"]."


ExportWCxf::renCouplings = "ExportWCxf currently only supports matching conditions determined with MapEffectiveCouplings with the Option ShiftRenCouplings set to True."


ExportWCxf::incompleteDict = "The following BSM parameter(s) were not specified in the input dictionary: `1`."


ExportWCxf::nonrenormalized = "WARNING: Matching conditions containing \[Epsilon] dependence detected. Dropping all \!\(\*FractionBox[\(1\), \(\[Epsilon]\)]\) poles and setting \[Epsilon]->0 afterwards."


ExportWCxf::missingparams = "Missing input parameters detected. The BSM parameters present in the matching conditions are `1`. The paramertes specifed in the input dictionary are `2`."


Options@ExportWCxf= {Point-> 1}


(*ExportWCxf[matchingCoditions:{Rule[_,_]..}, matchingScale_?Positive, numericInput_, dict:({(_Rule|_RuleDelayed)..}|_Association), outputDir_, opt:OptionsPattern[]] := Module[*)
ExportWCxf[arg_, inputFile_, outputDir_, OptionsPattern[]] := Module[
	(* 
	- arg:        either a BSM Lagrangian that can be matched onto the default SMEFT Lagrangian, or a list of matching conditions
	- inputFile:  path to the input JSON file containing all paramerte points and further information (matching scale, non-default SM parameter values/settings, notation dictionary)
	- outputDir:  path to the directory where the WCxf files will be stored
	*)
	{
		ReplaceBSM,
		ReplaceSM,
		MatchingConditions, matchingScale,
		smParams = {
			Global`gY, Global`gL, Global`gs, Global`Ye, Global`Yu, Global`Yd, Global`\[Mu]2, Global`\[Lambda],
			Global`cB2, Global`cW2, Global`cG2, Global`celH, Global`cquH, Global`cdqH, Global`cH2, Global`cH4, Global`CH2
		},
		dict, smOptions, params, numerics,
		NumericMatching,
		DefaultAssoc,
		tmp$$assoc,
		inputParam, rhsParam, lhsParam, inds, flavors=GetFlavorIndices[],
		smeftWC,
		flavorBasis,
		stratingPoint= OptionValue[Point]
	}
	,
	(* add literature *)
	AddLiterature[
	"WCxf export",
	"Exporting WCxf files from Matchete",
	"@article{Belfatto:2025ids,
    author = {Belfatto, Benedetta and Blanke, Monika and Heisig, Jan and Kr{\\\"a}mer, Michael and Rathmann, Lena and Wilsch, Felix},
    title = \"{Toward a Comprehensive Exploration of Flavored Dark Matter Models}\",
    eprint = \"2511.10490\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    reportNumber = \"P3H-25-091, TTP25-043, TTK-25-37\",
    month = \"11\",
    year = \"2025\"
}"
	];
	
	If[MatchQ[arg,{_Rule..}|_Association],
	(* if matching conditions directly provided as argument *)
		(* set MatchingConditions variable as a list  *)
		MatchingConditions = Normal[arg]; (* transform association to list if required *)
		(* check if default SMEFT Lagangian is loaded *)
		If[FreeQ[DownValues[LoadModel], HoldPattern[LoadModel["SMEFT_Warsaw"]]],
			Message[ExportWCxf::defaultSMEFT];
			Abort[]
		];
		(* check that ShiftRenCouplings was used and the RHS of the matching conditions does not contain default d=4 SM parameteres *)
		If[
			!FreeQ[
				Association@@Table[
					If[FreeQ[First[mc],Alternatives@@smParams,All],
						mc,
						Nothing
					],
					{mc, MatchingConditions}
				],
				Alternatives[Global`gY, Global`gL, Global`gs, Global`Ye, Global`Yu, Global`Yd, Global`\[Mu]2, Global`\[Lambda]],
				All
			],
			Message[ExportWCxf::renCouplings];
			Abort[]
		];
		(* check matching conditions are renormalized and drop poles if not so *)
		If[!FreeQ[MatchingConditions,\[Epsilon],All],		
			Message[ExportWCxf::nonrenormalized];
			MatchingConditions = MatchingConditions /. \[Epsilon]^-1->0 /. \[Epsilon]->0
		];
		,
	(* if BSM Lagrangian is provided as argument *)
		(* determine matchign conditions *)
		MatchingConditions = GetMatchingConditionsForWCxf[arg];
	];
	
	(* exclude SM rules from matching conditions *)
	MatchingConditions = MatchingConditions /. Rule[Coupling[Alternatives@@smParams,__],_]->Nothing;
	
	(* exclude B- & L-violating operators from matching conditions *)
	MatchingConditions = MatchingConditions /. Rule[Coupling[l:Global`cllHH|Global`cduu|Global`cduq|Global`cqqq|Global`cqqu,__],Except[0|0.]]:>(Message[ExportWCxf::blviolation,l]; Nothing);
	MatchingConditions = MatchingConditions /. Rule[Coupling[Global`cllHH|Global`cduu|Global`cduq|Global`cqqq|Global`cqqu,__],_]->Nothing;
	
	(* make matching conditions an association *)
	If[Head[MatchingConditions]===List, MatchingConditions=Association@@MatchingConditions];
	
	(* determine the parameters in matchign conditions *)
	lhsParam = DeleteDuplicatesBy[Cases[Keys@MatchingConditions,  _Coupling,All], First];
	rhsParam = DeleteDuplicatesBy[Cases[Values@MatchingConditions,_Coupling,All], First];
	
	(* determine parameters required as input parameters *)
	inputParam = Complement[rhsParam,lhsParam]/.Coupling[Alternatives@@smParams,___]->Nothing;
	
	(* determine all flavor index combinations *)
	inputParam = Flatten@Table[
		(* find indices and their dimensions *)
		inds = Cases[param,_Index,All];
		inds = Table[
			{ind,flavors[Last[ind]][IndexDimension]}
			,
			{ind,inds}
		];
		Table[param,Evaluate[Sequence@@inds]]
		,
		{param,inputParam}
	];
	Print["The following parameters where determined as input parameters: \n", Format[inputParam,NiceForm]];
	
	(* read in JSON input file *)
	{matchingScale, params, dict, numerics, smOptions} = ReadInputJSON[inputFile];
	If[Length[smOptions]==0, 
		smOptions = Nothing,
		smOptions = Normal[smOptions]
	];
	
	(* check that all required parameters are specifed in input *)
	If[Sort[inputParam]=!=(Sort[Keys[dict]]/.Index[n_?IntegerQ,_]->n),
		Message[ExportWCxf::missingparams, Format[Sort[inputParam],NiceForm], Format[(Sort[Keys[dict]]/.Index[n_?IntegerQ,_]->n),NiceForm]];
		Abort[]
	];
	
	(* derive replacements for SM and BSM parameters *)
	ReplaceBSM = NumericReplacementsBSM[params, numerics, dict]/.Index[n_?IntegerQ,_]->n;
	(* SM replacements must come after BSM due to overloading parameter names *)
	ReplaceSM  = NumericReplacementsSM[matchingScale,Sequence@@smOptions];
	
	(* make flavor indices explicit in matching conditions; after this all flavor indices are pure integers w/o Index Head *)
	MatchingConditions = ExpandFlavorInMatchingConditions[MatchingConditions];
	
	(* substitute SM values *)
	MatchingConditions = MatchingConditions /. ReplaceSM;
	
	(* derive function to numerically evaluate matching conditions *)
	NumericMatching = BuildEvaluationFunction[MatchingConditions, inputParam];
	
	flavorBasis = If[KeyExistsQ[smOptions,"SM"] && KeyExistsQ[smOptions["SM"],"flavor basis"],
		smOptions["SM"]["flavor basis"],
		"flavor basis" /. Options@ParametersSM
	];
	
	DefaultAssoc = <|
		"eft"   -> "SMEFT",
		"basis" -> Switch[flavorBasis,
			"down", "Warsaw",
			"up",   "Warsaw up"
		],
		"scale" -> matchingScale
	|>;
	
	Do[
		(* generate numeric values *)
		tmp$$assoc = NumericMatching[
			Sequence@@(Join[ReplaceBSM[[i]],{\[Mu]bar2->matchingScale^2}](*/.Index[i:1|2|3,_]:>i*))
		];
		
		(* covert to independent basis *)
		(*Echo@NiceForm[KeySelect[tmp$$assoc,!FreeQ[#,Global`cee,All]&]];*)
		tmp$$assoc = ToIndependentBasis@tmp$$assoc;
		(*Echo@NiceForm[KeySelect[tmp$$assoc,!FreeQ[#,Global`cee,All]&]];*)
		
		(* express Wilson coefficients in WCxf convention *)
		KeyDropFrom[tmp$$assoc,{Global`\[Mu]2[],Global`\[Lambda][],Global`cH2[]}];
		tmp$$assoc = KeyMap[MatcheteToWCxf,tmp$$assoc];
		
		(* replace real and imaginary parts *)
		tmp$$assoc = Evaluate/@(tmp$$assoc/.Complex[r_,i_]:><|"Re"->r,"Im"->i|>);
		
		(* create final association *)
		tmp$$assoc = Append[DefaultAssoc,"values"->tmp$$assoc];
				
		(* export the results *)
		(*ExportYAML[
			FileNameJoin[{outputDir,"ParamPoint_"<>ToString[i]}],
			tmp$$assoc
		];*)
		ExportJSON[
			FileNameJoin[{outputDir,"ParamPoint_"<>ToString[stratingPoint+i-1]}],
			tmp$$assoc
		];

		,
		{i, Length@ReplaceBSM (* loop over all BSM parameter points *)}
	];

]


(* ::Section:: *)
(*Reading input file*)


ExportWCxf::missingkey       = "The input JSON file must contain all of the following required keys: `1`"
ExportWCxf::positivescale    = "The specified matching scale (`1`) is not a positive number."
ExportWCxf::paramlist        = "The key \"parameters\" must have a list of all parameter names as value."
ExportWCxf::dictassoc        = "The key \"dictionary\" must have a dictionary as value."
ExportWCxf::intkeys          = "The keys of the \"parameter points\" dictionary must be strings of integer numbers."
ExportWCxf::listvalues       = "The values of the \"parameter points\" dictionary must be lists of numbers."
ExportWCxf::listvalueslength = "The values of the \"parameter points\" dictionary must be lists with length equal to the number of parametres specifed in \"parameters\"."
ExportWCxf::smdict           = "The value of \"SM\" must be a dictionary."
ExportWCxf::smdictkeys       = "The following keys of the \"SM\" dictionary are ignored: `1`. Only the following keys can be used: `2`."


ReadInputJSON[path_] := Module[
	{
		input = Import[path,"RawJSON"], (* RawJSON ensures that Associations are used instead of Lists *)
		requiredKeys = {"matching scale", "parameters", "dictionary", "parameter points"}
	}
	,
	(* check if all required keys exist *)
	If[!(And@@(KeyExistsQ[input,#]&/@requiredKeys)),
		Message[ExportWCxf::missingkey, requiredKeys];
		Abort[]
	];
	
	(* convert Strings to Expressions *)
	input["parameters"]       = ToExpression/@ input["parameters"];
	input["dictionary"]       = ToExpression/@ KeyMap[ToExpression, input["dictionary"]];
	input["parameter points"] = KeyMap[ToExpression, input["parameter points"]];
	
	(* check JSON is sensible *)
	If[!(NumberQ[input["matching scale"]]&&Positive[input["matching scale"]]),
	Message[ExportWCxf::positivescale, input["matching scale"]];
		Abort[]
	];
	If[Head[input["parameters"]]=!=List,
		Message[ExportWCxf::paramlist];
		Abort[]
	];
	If[Head[input["dictionary"]]=!=Association,
		Message[ExportWCxf::paramlist];
		Abort[]
	];
	If[!(And@@(IntegerQ/@Keys[input["parameter points"]])),
		Message[ExportWCxf::intkeys];
		Abort[]
	];
	If[!(And@@(((Head[#]===List) && (#\[Element]Reals))&/@Values[input["parameter points"]])),
		Message[ExportWCxf::listvalues];
		Abort[]
	];
	If[Or@@((#=!=Length[input["parameters"]])& /@ Length/@Values[input["parameter points"]]),
		Message[ExportWCxf::listvalueslength];
		Abort[]
	];
	If[KeyExistsQ[input["SM"]],
		If[Head[input["SM"]]===Association,
			If[!(And@@(MemberQ[Options[ParametersSM][[;;,1]],#]&/@Keys[input["SM"]])),
				Message[ExportWCxf::smdictkeys, 
					Complement[Keys[input["SM"]],Options[ParametersSM][[;;,1]]],
					Options[ParametersSM][[;;,1]]
				]
			]
			,
			Message[ExportWCxf::smdict];
			Abort[]
		];
	];
	
	
	(* return data *)
	{
		input["matching scale"],
		input["parameters"],
		input["dictionary"],
		input["parameter points"],
		If[KeyExistsQ[input,"SM"], input["SM"], <||>]
	}
]


(* ::Section:: *)
(*Determine matching conditions*)


Options@GetMatchingConditionsForWCxf= {ReductionIdentities-> EvanescenceFree}


GetMatchingConditionsForWCxf[arg_, OptionsPattern[]] := Module[
	{
		\[ScriptCapitalL]Warsaw = LoadModel["SMEFT_Warsaw"],
		matchingConditions
	},
	
	matchingConditions = MapEffectiveCouplings[arg, \[ScriptCapitalL]Warsaw
		,EOMSimplify                  -> True
		,ReductionIdentities          -> OptionValue[ReductionIdentities]
		,ShiftRenCouplings            -> True
		,AppendEffectiveCouplingsDefs -> True
	]/.\[Epsilon]^-1->0; (* renormalize *)
	
	matchingConditions
]


(* ::Section:: *)
(*Expand in explicit flavor indices*)


(* ::Subsection:: *)
(*ExpandFlavor*)


ExpandFlavorInMatchingConditions[mc_] := Module[
	{
		expr = mc, res,
		inds, flavors=GetFlavorIndices[],
		pattern, myCoeff
	}
	,
	(* delete vanishing coefficients *)
	expr = DeleteCases[expr,0];
	
	(* expand rhs and canonize flavor indices *)
	expr = RelabelIndices[#, Expand->True, Unique->False]& /@ expr;
	
	(* introduce explicit sums for Einstein convention *)
	expr = ExplicitEinsteinSums/@expr;
	
	(* consider all combinations of flavor indices for the Wilson coefficients *)
	res = Flatten@Table[
		myCoeff = (coeff/.Pattern->pattern)/.pattern[x_,_]:>x;
		inds    = Cases[myCoeff,_Index,All];
		inds    = Table[
			{ind,flavors[Last[ind]][IndexDimension]}
			,
			{ind,inds}
		];
		Flatten@Table[
			myCoeff -> expr[coeff]
			,
			Evaluate[Sequence@@inds]
		]
		,
		{coeff,Keys[expr]}
	];
	(* after this point all Index[n,Flavor] are replaced by n for integer n. *)
	
	(* replace the deltas *)
	res = res/.Delta->KroneckerDelta;
	
	(* make association *)
	res = Association@@res;
	
	(* delete vanishing coefficients *)
	res = DeleteCases[res,0];
	
	Return[res]
]


(* ::Section:: *)
(*Substituting Einstein summation convention by explicit sums*)


(* ::Text:: *)
(*NOTICE: This is now defined in ExportUFO.m*)


(* ::Subsection:: *)
(*ExplicitEinsteinSums*)


(*ExplicitEinsteinSums[sum_Plus]:=ExplicitEinsteinSums/@sum*)


(*ExplicitEinsteinSums[term:Except[_Plus]]:=Module[{repeatedInds,flavors=GetFlavorIndices[]},
	(* determine repeated indices *)
	repeatedInds=Matchete`PackageScope`FindDummyIndices[term];
	(* determine dimensions of these indices *)
	repeatedInds=Table[
		{ind,flavors[Last[ind]][IndexDimension]}
		,
		{ind,repeatedInds}
	];
	(* perform explicit sums over flavor indices *)
	If[Length[repeatedInds]>0, Sum[term,Evaluate[Sequence@@repeatedInds]],term]/.FlavorSum[_?NumericQ]->1
]*)


(* ::Section:: *)
(*Mapping to independent basis used by WCxf*)


(* ::Subsection:: *)
(*Categorize SMEFT Wilson coefficients*)


$Hermitian2Fermion=Alternatives[
	Global`cHl1, Global`cHl3, Global`cHe, Global`cHq1, Global`cHq3 , Global`cHu, Global`cHd
];


$Hermitian4FermionNonIdentical=Alternatives[
	Global`clq1, Global`clq3, Global`ceu, Global`ced, Global`cud1, Global`cud8, Global`cle, Global`clu, Global`cld, Global`cqe, Global`cqu1, Global`cqu8, Global`cqd1, Global`cqd8
];


$Hermitian4FermionIdentical=Alternatives[
	Global`cll, Global`cqq1, Global`cqq3, Global`cuu, Global`cdd, Global`cee
];


$Hermitian4Fermion=Flatten[$Hermitian4FermionNonIdentical|$Hermitian4FermionIdentical];


(* ::Subsection:: *)
(*Canonize SMEFT Wilson coefficients*)


CanonizeWC[arg_]:=CanonizeWC[arg]=Module[
	{expr = arg},
	(* 2-fermion *)
	expr=expr//.(
		Coupling[
			label:$Hermitian2Fermion,
			{p_?IntegerQ,r_?IntegerQ},
			aux_
		]/;p>r
	):>Conjugate@Coupling[label,{r,p},aux];
	
	(* 4-fermion [cee] *)
	expr=expr//.(
		Coupling[
			label:Global`cee,
			{p_?IntegerQ,r_?IntegerQ,s_?IntegerQ,t_?IntegerQ},
			aux_
		]/;(
			(((p+t)<(p+r))&&((p+t)<(s+t)))
			||
			(((s+r)<(p+r))&&((s+r)<(s+t)))
		)
	):>Coupling[label,{p,t,s,r},aux];
	
	(* 4-fermion [identical] *)
	expr=expr//.{
		(
			Coupling[
				label:$Hermitian4FermionIdentical,
				{p_?IntegerQ,r_?IntegerQ,s_?IntegerQ,t_?IntegerQ},
				aux_
			]/;((s<p&&s<r)||(t<p&&t<r)||(s==p&&t<r)||(s==r&&t<p)||(t==p&&s<r)||(t==r&&s<p))
		):>Coupling[label,{s,t,p,r},aux]
	};
	
	(* 4-fermion [all] *)
	expr=expr//.{
	(
		Coupling[
			label:$Hermitian4Fermion,
			{p_?IntegerQ,r_?IntegerQ,s_?IntegerQ,t_?IntegerQ},
			aux_
		]/;r<p
	):>Conjugate@Coupling[label,{r,p,t,s},aux]
	,
	(
		Coupling[
			label:$Hermitian4Fermion,
			{p_?IntegerQ,r_?IntegerQ,s_?IntegerQ,t_?IntegerQ},
			aux_
		]/;(r==p&&t<s)
	):>Conjugate@Coupling[label,{r,p,t,s},aux]
	}
	;
	
	(* 4-fermion [identical] remove conjugates *)
	expr=expr//.{
		Conjugate@Coupling[
			label:$Hermitian4FermionIdentical,
			{p_?IntegerQ,r_?IntegerQ,r_?IntegerQ,p_?IntegerQ},
			aux_
		]:>Coupling[label,{p,r,r,p},aux]
	};
	
	expr
]


(* ::Subsection:: *)
(*Map to independent basis*)


ToIndependentBasis::nonSymmetric="The coefficients `1` and `2` are not in the symmetric basis.";
ToIndependentBasis::nonHermitian="The given coefficients are not Hermitian.";


ToIndependentBasis[assoc_Association]:=Module[
	{
		wcLabel,
		canonicLabels,
		independentAssoc=<||>
	},
	(* create list of WCs and their canonized versions *)
	canonicLabels=Do[
		wcLabel=CanonizeWC[wc];
		
		(* check if conjugate *)
		If[Head[wcLabel]=!=Conjugate,
			(* if label is not conjugated add it to the list of indepednent coefficients *)
			If[!KeyExistsQ[independentAssoc,wcLabel],
				(* add canonic coefficient if it did not exist before *)
				AssociateTo[independentAssoc,wcLabel->assoc[wc]]
				,
				(* sum values of all coefficients that correspond to the same canonic coefficient *)
				AssociateTo[independentAssoc,wcLabel->(independentAssoc[wcLabel]+assoc[wc])]
			]
			,
			(* for conjugated labels check that a term and its conjugate are the hermitian of each other *)
			If[assoc[wc]!=Conjugate[assoc[wcLabel]],
				Message[ToIndependentBasis::nonHermitian];
				Abort[]
			]
		];
		,
		{wc,Keys[assoc]}
	];
	
	
	(* chop small entries *)
	independentAssoc=Chop[#,10^(-20)]&/@independentAssoc;
	(* delete vanishing coefficients *)
	independentAssoc=DeleteCases[independentAssoc,0];
	
	independentAssoc
]


(* ::Section:: *)
(*Efficient numeric evaluation*)


(* ::Subsection:: *)
(*Numeric evaluation*)


BuildEvaluationFunction[arg_,inputParam_] := Module[
	{
		res = arg,
		lf, $EvaluateLoopFunctions,
		internalParams,
		tmpCouplings = $TempCouplings,
		func, $Association
		,lfRules= {} (* for exporting the LF Taylor expansions as rules *)
	}
	,
	(* get all LF *)
	lf = DeleteDuplicates@Cases[res,_LF,All];
	(* associate LF to their evaluated form *)
	internalParams = Association@@Flatten[{#->$EvaluateLoopFunctions[#]}&/@lf];
	
	(* 
	 - collect all LF structures 
	 - determine all possible degeneracy limits
	 - determine all corresponding Taylor expansions once
	 - make $EvaluateLoopFunctions chose the Taylor expansion in the degeneracy limit instead of using the full LF
	*)
	Module[
		{  			
			taylorOrder = 2, (* order of the Taylor expansion *)  			
			acc = 0.02, (* degeneracy threshold when to use Taylor expansion of LF *)
			lfIndexList, lfList, degeneracyLimits, expansionsList, mNonDegenerate, mDegenerate, 
			(* patterns for LF with up to 6 masses *) 
			massesToPattern    = Table[ToExpression["m"<>ToString[i]] -> Pattern[Evaluate[ToExpression["m"<>ToString[i]]], Blank[]], {i,6}],
			massesToNumPattern = Table[ToExpression["m"<>ToString[i]] -> PatternTest[Pattern[Evaluate[ToExpression["m"<>ToString[i]]], Blank[]], NumericQ], {i,6}]
		},
		(* get indices of all LF *)
		lfIndexList = DeleteDuplicates@Cases[res,LF[_,ind_]:>ind,All];
		(* generate minimal set of generic LF *)
		lfList = LF[Table[ToExpression["m"<>ToString[i]],{i,Length[#]-1}],#]&/@lfIndexList;
		(* determine all degeneracy limits *)
		degeneracyLimits = Table[
			Subsets[Table[ToExpression["m"<>ToString[i]],{i,Length[inds]-1}],{2,Length[inds]-1}]
			,
			{inds,lfIndexList}
		];
		(* determine all Taylor expansions *)
		expansionsList = Do[
			Do[
				(* conditions for (non-)degeneracy *)
				mNonDegenerate = And@@Table[Abs[mNonDeg/degeneracyLimits[[i,j,-1]]]<(1-acc)||(1+acc)<Abs[mNonDeg/degeneracyLimits[[i,j,-1]]],{mNonDeg,Complement[First[lfList[[i]]],degeneracyLimits[[i,j]]]}];
				mDegenerate    = And@@Table[(1-acc)<=Abs[degeneracyLimits[[i,j,k]]/degeneracyLimits[[i,j,-1]]]<=(1+acc),{k,Length[degeneracyLimits[[i,j]]]-1}];
				(* set the corresponding DownValues *)
				With[
					{
						lhs = $EvaluateLoopFunctions[lfList[[i]] /. massesToNumPattern] /; Evaluate@And[mDegenerate,mNonDegenerate],
						rhs = BetterSeries[EvaluateLoopFunctions[lfList[[i]]],Sequence@@Table[{degeneracyLimits[[i,j,k]]}~Join~{degeneracyLimits[[i,j,-1]]}~Join~{taylorOrder},{k,Length[degeneracyLimits[[i,j]]]-1}]]
					}
					,
					AppendTo[lfRules, lhs:> rhs]; (* for exporting the Taylor exapnsions *)
					lhs := rhs
				]
				,
				{j,Length[degeneracyLimits[[i]]]}
			]
			,
			{i,Length@lfList}
		];
	];
	
	(* exporting Taylor expansions *)
	(*
	Export[
		FileNameJoin@{$MatchetePath, "Package", "Export", "LFExpansions.m"},
		(lfRules/. $EvaluateLoopFunctions-> Identity)/. pat_PatternTest:> First[pat]
	];
	*)
	
	(* add temporary couplings as internal parameters and remove them from matching conditions *)
	Do[
		If[!FreeQ[Keys[res],cp,All],
			(AssociateTo[internalParams,#->res[#]]; KeyDropFrom[res,#])&[FirstCase[Keys[res],Coupling[cp,___]]]
		]
		,
		{cp,tmpCouplings}
	];
	
	(* speed-up replacements later *)
	res = res/.hbar->N[1/(16\[Pi]^2)];
	res = Collect[#,_Coupling,Simplify]&/@res;
	
	(* generate a function that can evaluate the WC numerically *)
	With[{int=internalParams,ext=res},
		Options[func]=(#->0)&/@Join[inputParam,{\[Mu]bar2}];
		func[opt___:OptionsPattern[]] := Block[
			{
				preEval,
				rule={}
			},
			PrependTo[rule,#]&/@{opt};
			preEval = int/.hbar->N[1/(16\[Pi]^2)]/.rule;
			(* need to temporarly remove Association for DownValues of $EvaluateLoopFunctions to trigger *)
			preEval = preEval /. Association -> $Association;
			preEval = preEval /. $EvaluateLoopFunctions -> EvaluateLoopFunctions;
			preEval = preEval /. $Association -> Association;
			preEval = N/@preEval;
			N/@(Evaluate/@ext/.preEval/.rule)
		]
	];

	(* return this function *)
	func
]


(* ::Section:: *)
(*Translate Matchete to WCxf conventions*)


(* ::Subsection:: *)
(*Dictionary*)


WCxfDictionary[]:={
	ToExpression["cllHH"]->"cllHH",
	
	ToExpression["cG"]->"G",
	ToExpression["cGt"]->"Gtilde",
	ToExpression["cW"]->"W",
	ToExpression["cWt"]->"Wtilde",
	
	ToExpression["cH"]->"phi",
	
	ToExpression["cHBox"]->"phiBox",
	ToExpression["cHD"]->"phiD",
	
	ToExpression["cHG"]->"phiG",
	ToExpression["cHB"]->"phiB",
	ToExpression["cHW"]->"phiW",
	ToExpression["cHWB"]->"phiWB",
	ToExpression["cHGt"]->"phiGtilde",
	ToExpression["cHBt"]->"phiBtilde",
	ToExpression["cHWt"]->"phiWtilde",
	ToExpression["cHWtB"]->"phiWtildeB",
	
	ToExpression["cuH"]->"uphi",
	ToExpression["cdH"]->"dphi",
	ToExpression["ceH"]->"ephi",
	
	ToExpression["ceW"]->"eW",
	ToExpression["ceB"]->"eB",
	ToExpression["cuG"]->"uG",
	ToExpression["cuW"]->"uW",
	ToExpression["cuB"]->"uB",
	ToExpression["cdG"]->"dG",
	ToExpression["cdW"]->"dW",
	ToExpression["cdB"]->"dB",
	
	ToExpression["cHl1"]->"phil1",
	ToExpression["cHl3"]->"phil3",
	ToExpression["cHe"]->"phie",
	ToExpression["cHq1"]->"phiq1",
	ToExpression["cHq3"]->"phiq3",
	ToExpression["cHu"]->"phiu",
	ToExpression["cHd"]->"phid",
	ToExpression["cHud"]->"phiud",
	
	ToExpression["cll"]->"ll",
	ToExpression["cqq1"]->"qq1",
	ToExpression["cqq3"]->"qq3",
	ToExpression["clq1"]->"lq1",
	ToExpression["clq3"]->"lq3",
	
	ToExpression["cee"]->"ee",
	ToExpression["cuu"]->"uu",
	ToExpression["cdd"]->"dd",
	ToExpression["ceu"]->"eu",
	ToExpression["ced"]->"ed",
	ToExpression["cud1"]->"ud1",
	ToExpression["cud8"]->"ud8",
	
	ToExpression["cle"]->"le",
	ToExpression["clu"]->"lu",
	ToExpression["cld"]->"ld",
	ToExpression["cqe"]->"qe",
	ToExpression["cqu1"]->"qu1",
	ToExpression["cqu8"]->"qu8",
	ToExpression["cqd1"]->"qd1",
	ToExpression["cqd8"]->"qd8",
	
	ToExpression["cledq"]->"ledq",
	
	ToExpression["cquqd1"]->"quqd1",
	ToExpression["cquqd8"]->"quqd8",
	ToExpression["clequ1"]->"lequ1",
	ToExpression["clequ3"]->"lequ3"
}


(* ::Subsection:: *)
(*Translator*)


MatcheteToWCxf[arg_]:=MatcheteToWCxf[arg]=arg/.{
	Coupling[label_,{},_]:>(label/.WCxfDictionary[]),
	Coupling[label_,{p_,r_},_]:>(label/.WCxfDictionary[])<>"_"<>ToString[p]<>ToString[r],
	Coupling[label_,{p_,r_,s_,t_},_]:>(label/.WCxfDictionary[])<>"_"<>ToString[p]<>ToString[r]<>ToString[s]<>ToString[t]
}


(* ::Section:: *)
(*Exporting to JSON and YAML*)


(* ::Subsection:: *)
(*JSON export*)


ExportJSON[name_String,assoc_]:=Export[name<>".json",assoc,"JSON"]


(* ::Subsection:: *)
(*YAML export (mapping JSON to YAML)*)


ExportYAML[name_String,assoc_]:=Export[name<>".yaml",
	StringTrim@StringReplace[
		StringReplace[
			StringReplace[
				ExportString[assoc,"JSON"],
				{"\""|"{"|"}"|","->"",":"->": ","\n\t"->"\n"}
			]
			,
			"\n\t\n"->"\n"
		],
		"\t"->"  "
	]
,"Text"]
