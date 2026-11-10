(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`*)


(* ::Subtitle:: *)
(*Paclet for routines options and utility functions*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["$MatchetePath"]
PackageExport["CheckForUpdate"]
PackageExport["SuggestBibliography"]
PackageExport["AddLiterature"]


(* ::Text:: *)
(*Deprecated*)


PackageExport["DefineGroupRepresentation"]
PackageExport["DefineGroup"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["OptionsCheck"]
PackageScope["OptionTest"]
PackageScope["OptionMessage"]
PackageScope["FilterOptions"]


PackageScope["Defined"]
PackageScope["OptionalMonitor"]
PackageScope["PrintMessages"]
PackageScope["RemoveAssociatedDownValues"]
PackageScope["RemoveAssociatedUpValues"]


PackageScope["SubscriptStyle"]


PackageScope["PseudoTimes"]
PackageScope["ReleasePseudoTimes"]


PackageScope["FindPermutationOrder"]
PackageScope["InversePermutationOrder"]


PackageScope["ReplaceListSubExprs"]
PackageScope["ReplaceFirst"]
PackageScope["ReplaceShieldSubexpressions"]


PackageScope["SelectAndDeleteCases"]


PackageScope["TermsToList"]


PackageScope["ComplementWithDuplicates"]


PackageScope["DidMatchSwitch"]


PackageScope["IntegerSets"]
PackageScope["NonOverlappingPairs"]


PackageScope["BetterExpand"]
PackageScope["LagrangianExpand"]
PackageScope["FastExpand"]
PackageScope["LayeredExpand"]


PackageScope["BetterSeries"]
PackageScope["InvertMetric"]


PackageScope["BetterContractCGs"]


PackageScope["MyPrint"]


PackageScope["$PrintMessages"]


PackageScope["AddToBibliography"]


PackageScope["CheckVersionCompatibility"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


$MatchetePath::usage  = "$MatchetePath is the path to the Matchete package."


CheckForUpdate::usage = "CheckForUpdate[] compares the local Matchete version to the one in the git repository."


SuggestBibliography::usage = "SuggestBibliography[] automatically collects and prints the suggested bibliography based on the usage of functions in the current Matchete session. Option \"Explanation\"-> True can be used to get the reason behind the suggestion of the individual papers. Option \"References\"-> All can be used to get a list of all papers associated with the Matchete package.";


(* ::Text:: *)
(*Deprecated*)


DefineGroupRepresentation::usage= "\"DefineGroupRepresentation\" is deprecated as of v0.3.0. Please refer to \"DefineRepresentation\" instead."
DefineGroup::usage= "\"DefineGroup\" is deprecated as of v0.3.0. Please refer to \"DefineCGGroup\" instead."


(* ::Subsubsection::Closed:: *)
(*Internal*)


OptionsCheck::usage  = "OptionsCheck can be applied on any function with optional arguments and checks the validity of optional arguments.";
OptionTest::usage    = "OptionTest[opt,func,val] must be implemented such that it returns True if val is an allowed value for the option opt in the function func, and False otherwise. ";
MyPrint::usage       = "MyPrint[message,Verbose -> True/False] is a printing function that can be deactivated when Verbose is set to False. The Verbose option is set to the global flag $PrintMessages, which can be changed using the PrintMessages routine."


PrintMessages::usage  = "PrintMessages[True/False] sets whether information messages are displayed by some of the routines."


TermsToList::usage=
"Transform a sum of terms into a list, or convert a single term into a list."


(* ::Chapter:: *)
(*Private:*)


$MatchetePath=DirectoryName[$InputFileName,2];


(* ::Section:: *)
(*Version and updates*)


(* ::Subsubsection::Closed:: *)
(*Version information*)


$MatcheteVersionURL = "https://gitlab.com/matchete/matchete/-/raw/master/version";


GetVersionString[]:=StringTrim[First@StringSplit[Import@FileNameJoin[{$MatchetePath,"version"}],"\n"]]


(* ::Subsubsection::Closed:: *)
(*Update package*)


CheckForUpdate[]:=Module[{importString,nrOnly,yn},
	importString = Import[$MatcheteVersionURL];
	If[importString === $Failed || Head[importString]=!=String, Print["Could not fetch version number from repository."]; Return[]];
	nrOnly = First @ StringSplit[importString, "\n"];

	If[StringTrim[nrOnly] === GetVersionString[],
		MessageDialog["Matchete is up-to-date.\nYour version: "<>GetVersionString[]],
		yn = ChoiceDialog["A new Matchete version is available!\nYour version: "<> GetVersionString[]<>"\nNew version: "<>StringTrim@nrOnly<>" ("<>StringSplit[importString, "\n"][[2]]<>")\n\n"<>"Do you want to update?",{"Yes"->True,"No"->False}];
		If[yn,UpdateMatchete[]]
	];
]


UpdateMatchete[]:=Module[{},
	Import["https://gitlab.com/matchete/matchete/-/raw/master/install.m"]
]


(* ::Subsubsection::Closed:: *)
(*Check version compatibility*)


CheckVersionCompatibility::usage="checks that the current MAtchete version is the same or newer than the version given as argument.";


CheckVersionCompatibility::inval="Invalid version number (`1`) given. Version numbers must follow the format \"a.b.c\", where a, b, c can each be only a sequence of integers. Cannot check whether your Matchete version is compatible with this code. I am still running the code, but watch our for problems.";


CheckVersionCompatibility::incomp="In compatible Matchete version detected. The feature you are using requires Mathematica v`1`, whereas you are currently using v`2`. You can update to the latest version using the function: CheckForUpdate[]";


GetVersion[]:=Import@FileNameJoin[{$MatchetePath,"version"}];


CheckVersionCompatibility[v_String]:= Module[{currentVersion=GetVersion[],requestedVersion=v},
	(* check that the version numbers have a valid format *)
	If[!StringMatchQ[currentVersion,(DigitCharacter..)~~"."~~(DigitCharacter..)~~"."~~(DigitCharacter..)],
		Message[CheckVersionCompatibility::inval,currentVersion]; Return[True]
	];
	If[!StringMatchQ[requestedVersion,(DigitCharacter..)~~"."~~(DigitCharacter..)~~"."~~(DigitCharacter..)],
		Message[CheckVersionCompatibility::inval,requestedVersion]; Return[True]
	];
	
	(* split version numbers and convert to integers *)
	currentVersion= ToExpression/@StringSplit[currentVersion,"."];
	requestedVersion= ToExpression/@StringSplit[requestedVersion,"."];
	
	(* compare *)
	If[currentVersion[[1]]>requestedVersion[[1]],
		Return[True]
		,
		If[currentVersion[[1]]==requestedVersion[[1]],
			If[currentVersion[[2]]>requestedVersion[[2]],
				Return[True]
				,
				If[currentVersion[[2]]==requestedVersion[[2]],
					If[currentVersion[[3]]>=requestedVersion[[3]],
						Return[True]
					];
				];
			];
		];
	];
	Message[CheckVersionCompatibility::incomp,v,GetVersion[]];
	Return[False]
]


(* ::Section:: *)
(*Meta level functions*)


(* ::Subsection:: *)
(*Options*)


(* ::Subsubsection::Closed:: *)
(*OptionsChecker*)


(* ::Text:: *)
(*General function for checking Options of functions*)


SetAttributes[OptionsCheck, HoldFirst];
OptionsCheck @ func_[___, opts : OptionsPattern[]] := And[
	And@@ (Message[General::invalidopt, #1, func] &)@@@ FilterRules[List@ opts, Except@ Options@ func],
	And@@ (OptionTest[func, #1][#2] || OptionMessage[#1, func, #2] &)@@@ FilterRules[List@ opts, Options@ func]
];


General::invalidopt = "Invalid option `1` given for function `2`.";
General::invalidarg = "Option `1` for function `2` received invalid value `3`.";
General::optexpectsval = "Option `1` for function `2` received invalid value `3`. A `4` is expected.";
OptionMessage[opt_, func_, val_] := Message[General::invalidarg, opt, func, val];


(* ::Subsubsection::Closed:: *)
(*Option tests*)


(* ::Text:: *)
(*Specific tests for options*)


(*KEEP ALPHABETICAL in option name please*)
OptionTest[_, AdjAlphabet]                      = #===None||(ListQ[#]&&And@@(Head[#1]===String&/@#))&;
OptionTest[_, AppendEffectiveCouplingsDefs]     = BooleanQ;
OptionTest[_, CanonicallyNormalized]            = BooleanQ;
OptionTest[_, ChargeNeutral]                    = BooleanQ;
OptionTest[_, Chiral]                           = MatchQ[False|LeftHanded|RightHanded];
OptionTest[_, ClosedSpinChains]                 = BooleanQ;
OptionTest[CollectOperators, NormalForm]        = BooleanQ;
OptionTest[_, ContractedIndices]                = BooleanQ;
OptionTest[_, ContractOverallDelta]             = BooleanQ;
OptionTest[_, "Database"]                       = BooleanQ;
OptionTest[_, DefaultValues]                    = BooleanQ;
OptionTest[_, DetailedOutput]                   = BooleanQ;
OptionTest[_, DummyCoefficients]                = BooleanQ;
OptionTest[_, EffectiveCouplingSymbol]          = StringQ;
OptionTest[EOMSimplify, EFTOrder]               = MatchQ[Automatic| a_Integer/; a>= 4];
OptionTest[_, EFTOrder]                         = MatchQ[{_Integer?Positive}| (_Integer?Positive)];
OptionTest[ExportUFO, FeynmanRules]             := Module[{pat= Alternatives@@ Keys[$FieldAssociation]}, MatchQ[#, a_Association/; And@@ (MatchQ[#, l_List/; (And@@ (MatchQ[#/. (*OverBar|*)Conj-> Identity, pat] &/@l)) ]&/@ Keys[a])]]&;
OptionTest[_, FlavorInvariants]                 = AssociationQ;
OptionTest[_, FundAlphabet]                     = #===None||(ListQ[#]&&And@@(Head[#1]===String&/@#))&;
OptionTest[_, GaugeInvariant]                   = BooleanQ;
OptionTest[_, FreeOfHeavyTadpoles]              = BooleanQ;
OptionTest[DefaultParamCard | ExportUFO, Gauge] = MatchQ[Automatic| Unitary| Manual| R\[Xi]];
OptionTest[_, GaugeAnomalies]                   = BooleanQ;
OptionTest[_, GoldstoneBoson]                   = BooleanQ;
OptionTest[_, HeavyMassBasis]                   = BooleanQ;
OptionTest[_, Hermiticity]                      = BooleanQ;
OptionTest[_, IndexAlphabet]                    = #===None||(ListQ[#]&&And@@(Head[#1]===String&/@#))&;
OptionTest[ExportUFO, Information]              = AssociationQ;
OptionTest[_, InputFile]                        = FileExistsQ;
OptionTest[_, "Internal"]                       = BooleanQ;
OptionTest[_, KeepTrivalReplacements]           = BooleanQ;
OptionTest[_, Legs]                             = MatchQ[All|_Integer?Positive|{_Integer?Positive}];
OptionTest[EOMSimplify, LoopOrder]              = MatchQ[0| 1| Automatic];
OptionTest[_, LoopOrder]                        = MatchQ[0| 1| {1}];
(*OptionTest[_, Mode]                             = MatchQ["Standard"|"MadGraph"];*)
OptionTest[_, ModelParameters]                  = (ListQ[#] && And@@(Head[#1]===Rule &/@#))&;
OptionTest[_, OutputDirectory]                  = MatchQ[Automatic|_?PathFormQ];
OptionTest[_, Path]                             = BooleanQ;
OptionTest[_, PreferInputOperators]             = BooleanQ;
OptionTest[_, ReductionIdentities]              = MatchQ[dDimensional|Evanescent|EvanescenceFree|FourDimensional];
OptionTest[_, Rules]                            = BooleanQ;
OptionTest[DefaultParamCard | ExportUFO, R\[Xi]]    = MatchQ[1| General];
OptionTest[_, SelfConjugate]                    = Or[BooleanQ[#],VectorQ[#,Positive]]&;
OptionTest[_, Simplifications]                  = MatchQ[All| None];
OptionTest[_, Simplify]                         = BooleanQ;
OptionTest[_, SortByEFTOrder]                   = BooleanQ;
OptionTest[_, Symmetries]                       = ListQ;
OptionTest[_, ThermalPowerCounting]             = MatchQ[Automatic| {_Integer?Positive}| (_Integer?Positive) | ThermalPowerCounting[_]];
OptionTest[_, TimeConstraint]                   = Positive;
OptionTest[_, UFO]                              = AssociationQ;
OptionTest[_, UndefinedObject]                  = BooleanQ;
OptionTest[_, Variables]                        = MatchQ[_Symbol| {_Symbol...}];
OptionTest[Match, Verbose]                      = MatchQ[Print|Monitor|None];
OptionTest[_, Verbose]                          = BooleanQ;
OptionTest[_, VLFBasis]                         = MatchQ[Automatic| Chiral(*| VectorLike*)];
OptionTest[_, WhichTraces]                      = MatchQ[{_List..}|All];
OptionTest[_, Dimensions]                       = IntegerQ;


(* ::Subsubsection::Closed:: *)
(*Error messages*)


(* ::Text:: *)
(*Error messages for the options*)


(*KEEP ALPHABETICAL in option name please*)
OptionMessage[AdjAlphabet, func_, val_]                         := Message[General::optexpectsval, AdjAlphabet, func, val, "list of strings or None"];
OptionMessage[AppendEffectiveCouplingsDefs, func_, val_]        := Message[General::optexpectsval, Verbose, func, val, "Boolean"];
OptionMessage[Chiral, func_, val_]                              := Message[General::optexpectsval, Chiral, func, val, "value False, LeftHanded, or RightHanded"];
OptionMessage[DummyCoefficients, func_, val_]                   := Message[General::optexpectsval, DummyCoefficients, func, val, "Boolean"];
OptionMessage[EffectiveCouplingSymbol, func_, val_]             := Message[General::optexpectsval, EffectiveCouplingSymbol, func, val, "String"];
OptionMessage[EFTOrder, EOMSimplify, val_]                      := Message[General::optexpectsval, EFTOrder, EOMSimplify, val, "integer >=4 or the value Automatic"];
OptionMessage[EFTOrder, func_, val_]                            := Message[General::optexpectsval, EFTOrder, func, val, "positive integer or List with one positive integer"];
OptionMessage[FeynmanRules, ExportUFO, val_]                    := Message[General::optexpectsval, FeynmanRules, ExportUFO, Format[Short[val,20],NiceForm], "association with every key being  a list of field labels"];
OptionMessage[FlavorInvariants, func_, val_]                    := Message[General::optexpectsval, FlavorInvariants, func, val, "association"];
OptionMessage[FundAlphabet, func_, val_]                        := Message[General::optexpectsval, FundAlphabet, func, val, "list of strings or None"];
OptionMessage[Gauge, func:(_DefaultParamCard|_ExportUFO), val_] := Message[General::optexpectsval, Gauge, func, val, "value Automatic, Unitary, R\[Xi], or Manual"];
OptionMessage[GoldstoneBoson, func_, val_]                      := Message[General::optexpectsval, GoldstoneBoson, func, val, "Boolean"];
OptionMessage[IndexAlphabet, func_, val_]                       := Message[General::optexpectsval, IndexAlphabet, func, val, "list of strings or None"];
OptionMessage[Information, ExportUFO, val_]                     := Message[General::optexpectsval, Information, ExportUFO, val, "association"];
OptionMessage[KeepTrivalReplacements, func_, val_]              := Message[General::optexpectsval, Verbose, func, val, "Boolean"];
OptionMessage[Legs, func_, val_]                                := Message[General::optexpectsval, Legs, func, val, "possitive integer, a list with a single positive integer entry, or All"];
OptionMessage[LoopOrder, func_, val_]                           := Message[General::optexpectsval, LoopOrder, func, val, "value 0, 1 or {1}"];
(*OptionMessage[Mode, func_, val_]                                := Message[General::optexpectsval, Mode, func, val, "value \"MadGraph\" or \"Standard\""];*)
OptionMessage[ModelParameters, func_, val_]                     := Message[General::optexpectsval, ModelParameters, func, val, "list of replacement rules"];
OptionMessage[NormalForm, func_, val_]                          := Message[General::optexpectsval, NormalForm, func, val, "Boolean"];
OptionMessage[OutputDirectory, func_, val_]                     := Message[General::optexpectsval, OutputDirectory, func, val, "path to the directory or Automatic"];
OptionMessage[Path, func_, val_]                                := Message[General::optexpectsval, Path, func, val, "String"];
OptionMessage[PreferInputOperators, func_, val_]                := Message[General::optexpectsval, PreferInputOperators, func, val, "Boolean"];
OptionMessage[ReductionIdentities, func_, val_]                 := Message[General::optexpectsval, ReductionIdentities, func, val, "value dDimensional, Evanescent, EvanescenceFree, or FourDimensional"];
OptionMessage[R\[Xi], func:(_DefaultParamCard|_ExportUFO), val_]    := Message[General::optexpectsval, R\[Xi], func, val, "value 1 or General"];
OptionMessage[SelfConjugate, func_, val_]                       := Message[General::optexpectsval, SelfConjugate, func, val, "boolean (True or False) or a list of positive integers indicating the index positions"];
OptionMessage[Simplifications, func_, val_]                     := Message[General::optexpectsval, Simplifications, func, val, "value All or None"];
OptionMessage[Simplify, func_, val_]                            := Message[General::optexpectsval, Simplify, func, val, "Boolean"];
OptionMessage[SortByEFTOrder, func_, val_]                      := Message[General::optexpectsval, Verbose, func, val, "Boolean"];
OptionMessage[Symmetries, func_, val_]                          := Message[General::optexpectsval, Symmetries, func, val, "list of index symmetries"];
OptionMessage[TimeConstraint, func_, val_]                      := Message[General::optexpectsval, TimeConstraint, func, val, "positive number"];
OptionMessage[UFO, func_, val_]                                 := Message[General::optexpectsval, UFO, func, val, "Association (typically with Keys \"pdg\" and \"name\")"];
OptionMessage[Variables, func_, val_]                           := Message[General::optexpectsval, Variables, func, val, " (list of) coupling label(s)"];
OptionMessage[Verbose, Match, val_]                             := Message[General::optexpectsval, Verbose, Match, val, "Print, Monitor or None"];
OptionMessage[Verbose, func_, val_]                             := Message[General::optexpectsval, Verbose, func, val, "Boolean"];
OptionMessage[VLFBasis, func_, val_]                            := Message[General::optexpectsval, VLFBasis, func, val, "value Automatic or Chiral"(*"value Automatic, Chiral, or VectorLike"*)];


(* ::Subsubsection::Closed:: *)
(*Filter options*)


(* ::Text:: *)
(*Function to filter a set option arguments given valid options for a function*)


FilterOptions[opts___Rule, func_Symbol]:= 
	Sequence@@ FilterRules[{opts}, Options@ func];


(* ::Subsubsection::Closed:: *)
(*PathFormQ*)


PathFormQ[s_String] := Module[{},
  Or[
    (* Unix/macOS absolute path: starts with / *)
    StringMatchQ[s, "/" ~~ ___],

    (* Windows absolute path: drive letter like C:\ or C:/ *)
    StringMatchQ[s, RegularExpression["[A-Za-z]:[/\\\\].*"]],

    (* UNC network path: \\server\share *)
    StringMatchQ[s, "\\\\" ~~ __],

    (* Relative path: contains / or \ somewhere (but not just a plain word) *)
    StringContainsQ[s, "/" | "\\"]
  ]
]

PathFormQ[_] := False  (* Non-strings are never paths *)


(* ::Subsection:: *)
(*Check and change functions*)


(* ::Subsubsection::Closed:: *)
(*Defined*)


(* An auxiliary function that determines if a variable has already been defined *)
Defined[label_] := ValueQ[label] ||
               Head@label =!= Symbol ||
               Attributes[label] =!= {} ||
               DownValues[label] =!= {} ||
               SubValues[label] =!= {} ||
               Head[label::usage]=!=MessageName


(* ::Subsubsection::Closed:: *)
(*Remove associated DownValues and UpValues*)


(* ::Text:: *)
(*Selectively unsets all DownValues of a symbol, whose arguments matches the given pattern*)


SetAttributes[RemoveAssociatedDownValues, HoldAllComplete];


RemoveAssociatedDownValues[func_[arg___]]:= Module[{values},
	DownValues@ func= DeleteCases[DownValues@ func, 
		HoldPattern[Verbatim[HoldPattern][func[p:PatternSequence@ arg]]:> _]];
];


(* ::Text:: *)
(*Selectively unsets all UpValues of a symbol, whose arguments matches the given pattern*)


SetAttributes[RemoveAssociatedUpValues, HoldAllComplete];


RemoveAssociatedUpValues[func_[arg___]]:= Block[{},
	UpValues@ func= DeleteCases[UpValues@ func,
		HoldPattern[Verbatim[HoldPattern][_[___, func[PatternSequence@ arg], ___]]:> _]];
];


(* ::Subsection:: *)
(*Print function that can be globally deactivated*)


(* ::Subsubsection::Closed:: *)
(*OptionalMonitor*)


(* ::Text:: *)
(*A version of monitor that can be disabled (while still executing the code), controlled by the 3rd Boolean argument*)


SetAttributes[OptionalMonitor, HoldRest];
OptionalMonitor[True, expr_, mon_]:= Monitor[expr, mon];
OptionalMonitor[False, expr_, mon_]:= expr;


(* ::Subsubsection::Closed:: *)
(*Setter function for the printing flag*)


PrintMessages::flag="The flag '`1`' is not a boolean. Please use either True or False."


$PrintMessages=False;
PrintMessages[flag_]:=Module[{},
	If[!BooleanQ[flag],
		Message[PrintMessages::flag,flag];
		Abort[]
	];

	$PrintMessages=flag;
]


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[MyPrint]={Verbose :> $PrintMessages}


(* ::Subsubsection::Closed:: *)
(*MyPrint*)


MyPrint[string__,OptionsPattern[]]? OptionsCheck:=Module[{},
	If[OptionValue@Verbose,Print[string]];
];


(* ::Subsection:: *)
(*Other*)


(* ::Subsubsection::Closed:: *)
(*Auto-completion function*)


AddAutoCompletion[function_String][args___]:=Module[{processed},
	processed=ReplaceAll[{args},
	{
		None->0,
		"AbsoluteFileName"->2,
		"RelativeFileName"->3,
		"Color"->4,
		"PackageName"->7,
		"DirectoryName"->8,
		"InterpreterType"->9
	}
	];
	Function[FE`Evaluate@FEPrivate`AddSpecialArgCompletion@#][function->processed]
]


(* ::Subsubsection::Closed:: *)
(*Change Matchete Logo*)


ChangeMatcheteLogo[logo_String]:=Switch[logo,
	"Matchotter"
	,
	Global`$logoMatchete= Import[FileNameJoin[{$MatchetePath,"Kernel","LogoMatchotter.svg"}]];
	,
	_| "Matchete"
	,
	Global`$logoMatchete= Import[FileNameJoin[{$MatchetePath,"Kernel","LogoMatchete.eps"}]];
];


(* ::Section:: *)
(*Utility functions*)


(* ::Subsection:: *)
(*Expansion*)


(* ::Text:: *)
(*Smarter versions of the Expand function, which can be incredibly slow at times *)


(* ::Subsubsection::Closed:: *)
(*BetterExpand*)


(* ::Text:: *)
(*BetterExpand distributes over Plus, for MUCH better performance on long expressions. Why is this  not default Mathematica behaviour?*)


BetterExpand@ expr_Plus:= Expand/@ expr;
BetterExpand@ expr_:= Expand@ expr;


(* ::Subsubsection::Closed:: *)
(*LagrangianExpand*)


(* ::Text:: *)
(*Distributes all sums that may appear in Lagrangian terms (including positive powers)*)


LagrangianExpand@ expr_Plus:= LagrangianExpand/@ expr;

LagrangianExpand@ expr_Times:= expr//. {
	pwr:Power[_Plus, _Integer? Positive]:> Expand@ pwr, (* ListConvolve mght be better for large sums, but only then *)
	(*y:Times[rest___,sum1_Plus,sum2_Plus]:> ((rest*#)&/@ (Plus@@ ListConvolve[List@@sum1, List@@sum2, 1]));*) (* unclear if this really improves the performance -> for very large sums it might be better, but not for small ones *)
	x_Times/; !FreeQ[x, _Plus] && FreeQ[x, Power[_Plus, _Integer? Positive]]:> Distribute@ x
}
	
LagrangianExpand[expr:Power[_, _Integer? Positive]]:= Expand@ expr; (* ListConvolve mght be better for large sums, but only then *)

LagrangianExpand@ expr_:= expr;


(*Dist@ expr_Plus:= Dist/@ expr;
Dist@ term_Times:= Distribute[term, Plus, Times];*)


(* ::Subsubsection::Closed:: *)
(*FastExpand (experimental)*)


(* ::Text:: *)
(*Fast expansion on all levels for very large expressions (with several layers of Times & Plus) *)


(* ::Text:: *)
(*Can throw an error if $RecursionLimit is set too low, but might crash the Wolfram kernel if $RecursionLimit is set too high.*)


FastExpand[arg_]:=Module[{res,myPlus,recLimit=$RecursionLimit},
	$RecursionLimit=10000; (* this function is fast but uses many recursions *)
	(* Wrap all sums on all level with a MultiplicationBox *)
	res=arg//.sum_Plus:>MultiplicationBox[myPlus@@sum];
	(* the change Plus <-> myPlus is necessary for ReplaceRepeated to terminate *)
	res=res//.myPlus->Plus;
	(* remove MultiplicationBox after all UpValues triggered *)
	$RecursionLimit=recLimit;
	res//.MultiplicationBox->Identity
]


(* fast multiplication of 2 sums *)
MultiplicationBox/:MultiplicationBox[sum1_Plus]*MultiplicationBox[sum2_Plus]:=MultiplicationBox[
	Plus@@ListConvolve[List@@sum1,List@@sum2,1]
]

(* fast expansions of powers of sums *)
MultiplicationBox/:Power[MultiplicationBox[sum_Plus],pow_/;(IntegerQ[pow]&&pow>=2)]:=Power[MultiplicationBox[sum],pow-2]*MultiplicationBox[Plus@@ListConvolve[List@@sum,List@@sum,1]]

(* fast multiplication with prefactors *)
MultiplicationBox/:(coeff:Except[_MultiplicationBox|_Plus|Power[_MultiplicationBox|_Plus,_]])*MultiplicationBox[sum_Plus]:=MultiplicationBox[
	Plus@@ListConvolve[{coeff},List@@sum,1]
]

(* combine sums - should never be necessary *)
MultiplicationBox/:x_+MultiplicationBox[sum_Plus]:=MultiplicationBox[x+sum]

(* flatten MultiplicationBox *)
MultiplicationBox@MultiplicationBox[x_]:=MultiplicationBox[x]

(* simplify vanishing expressions *)
MultiplicationBox[0]=0;

(* some additional definitions *)
MultiplicationBox/:Bar[MultiplicationBox[arg_]]:=MultiplicationBox[Bar[arg]]
MultiplicationBox/:Transp[MultiplicationBox[arg_]]:=MultiplicationBox[Transp[arg]]
MultiplicationBox/:CConj[MultiplicationBox[arg_]]:=MultiplicationBox[CConj[arg]]


(* ::Subsubsection::Closed:: *)
(*LayeredExpand (very experimental)*)


(* thread over sums *)
LayeredExpand[arg_Plus]:=LayeredExpand/@arg

(* determine depth of expression *)
LayeredExpand[arg_]:=LayeredExpand[arg,Depth[arg]]

(* expand sums starting from deepest level *)
LayeredExpand[arg_,depth_]:=LayeredExpand[
	Replace[arg,{Times[coeff__,sum_Plus]:>Plus@@ListConvolve[{coeff},List@@sum,1]},{depth}],
	depth-1
]

(* stop expansion once level 0 is reached *)
LayeredExpand[arg_,-1]:=arg


(* ::Subsection:: *)
(*Series expansion*)


(* ::Text:: *)
(*Version of the Series function that fixes a bug since Mathematica 14.3*)


(* ::Subsubsection::Closed:: *)
(*BetterSeries*)


BetterSeries[f_, {x_, x0_, nmax_}]:= Normal[ Series[f, {x, x0, nmax}] + O[x, x0]^(nmax+1) ]


BetterSeries[f_, {x_, x0_, nmax_}, multiExp:{_,_,_}..]:= Normal[ Series[f, {x, x0, nmax}, multiExp] + O[x, x0]^(nmax+1) + Plus@@(O[#[[1]],#[[2]]]^(#[[3]]+1)&/@{multiExp})]


(* ::Subsection:: *)
(*Matrices *)


(* ::Subsubsection::Closed:: *)
(*Inverse metric*)


(* ::Text:: *)
(*Invert a metric (e.g. a gauge coupling matrix). Symmetric, positive definite, square matrix*)


InvertMetric::blocksize= "The block `1` of the `2` is too large for symbolic inversion.";
InvertMetric::nonsym= "The matrix `1` is not symmetric.";


InvertMetric@ matrix_List:= Module[{block, blockComponents, out},
	out= ConstantArray[0, Dimensions@ matrix];
	If[SimplifyCouplings[matrix- Transpose@ matrix] =!= out, 
		Message[InvertMetric::nonsym, matrix];
		Abort[];
	];
	blockComponents= MatrixBlockComponents@ matrix;
	Do[
		If[Length@ block > 2, 
			Message[InvertMetric::blocksize, block, matrix];
			Abort[];
		];
	, {block, blockComponents}];
	
	(*Invert each block separately*)
	Do[
		out[[block, block]]= Inverse@ matrix[[block, block]];
	, {block, blockComponents}];
	
	out
];


(* ::Text:: *)
(*Returns a list of list with the connected components of a matrix. E.g. {{1,3}, {2, 4}}.*)


MatrixBlockComponents::nonsquare= "The matrix `1` is not square.";


MatrixBlockComponents@ matrix_List:= Block[{adjacencyMatrix},
	If[!SquareMatrixQ@ matrix,
		Message[MatrixBlockComponents::nonsquare, matrix];
		Abort[];
	];
	(*Make unit connections*)
	adjacencyMatrix= Map[Boole[# =!= 0]&, matrix, {2}];
	ConnectedComponents@ AdjacencyGraph@ adjacencyMatrix
];


(* ::Subsection:: *)
(*Permutations*)


(* ::Subsubsection::Closed:: *)
(*FindPermutationOrder*)


(* ::Text:: *)
(*Returning the ordering list needed to make  permutation[[ordering list]] === target*)


FindPermutationOrder[permutation_List, target_List]:=
	Permute[Range@ Length@ target, FindPermutation[permutation, target]]


(* ::Subsubsection::Closed:: *)
(*InversePermutation*)


(* ::Text:: *)
(*Produces the inverse of the given permutation*)


InversePermutationOrder[permutation_List]:= Ordering@ permutation;


(* ::Subsection:: *)
(*Replacement functions*)


(* ::Subsubsection::Closed:: *)
(*ReplaceListSubExprs*)


(* ::Text:: *)
(*Function mimicking ReplaceList but on all subexpressions*)


(*This implementation does not apply the rules in all posible ways at each subexpression*)
(*ReplaceListSubExprs[expr_, rule_Rule|rule_RuleDelayed]:=
	MapAt[Function[{x}, x/. rule], expr, #]&/@ Position[expr, First@ rule, Infinity];*)


ReplaceListSubExprs[expr_, rule_Rule|rule_RuleDelayed]:= Module[{op, pos, rep},
	Flatten[Table[
			op= expr;
			op[[Sequence@@ pos]]= rep;
			op
		, {pos, Position[expr, First@ rule]}
		, {rep, ReplaceList[expr[[Sequence@@ pos]], rule]}]
	, 1]
]


ReplaceListSubExprs[expr_, rules:{_Rule|_RuleDelayed...}]:= Module[{rule},
	Join@@ Table[ReplaceListSubExprs[expr, rule], {rule, rules}]
]


(* ::Subsubsection::Closed:: *)
(*ReplaceFirst*)


(* ::Text:: *)
(*Function applying a replacement rule once on the first match to the rule*)


ReplaceFirst[expr_, rule_Rule|rule_RuleDelayed]:=
	MapAt[Function[{x}, x/. rule], expr, FirstPosition[expr, First@ rule, {}] ];


(* ::Subsubsection::Closed:: *)
(*ReplaceShieldSubexpressions*)


ReplaceShieldSubexpressions[expr_, rule_Rule|rule_RuleDelayed, shieldPattern_]:=
	ReplaceShieldSubexpressions[expr, {rule}, shieldPattern];


ReplaceShieldSubexpressions[expr_, rules:{(_Rule|_RuleDelayed)..}, shieldPattern_]:= Module[{shieldPos, replacePos},
	shieldPos= Position[expr, shieldPattern];
	shieldPos= Alternatives@@ (Append[#, ___]&)/@ shieldPos;
	replacePos= Position[expr, Alternatives@@ First/@ rules];
	(*Remove all replacement positions that are in the shielded part of the expression*)
	replacePos= Select[replacePos, Not@* MatchQ[shieldPos]];
	ReplaceAt[expr, rules, replacePos]
]


(* ::Subsection:: *)
(*Combinatorics function*)


(* ::Subsubsection::Closed:: *)
(*IntegerSet*)


(* ::Text:: *)
(*IntegerSet[s,n] returns all ordered sets of n integers {Subscript[\[Mu], 1],...,Subscript[\[Mu], n]}, such that Subscript[\[CapitalSigma], k] Subscript[\[Mu], k]=s and Subscript[\[Mu], k]>=0.*)


IntegerSets[sum_, ints_]:= Flatten[Permutations@ PadRight[#, ints]&/@
	DeleteCases[IntegerPartitions@ sum, _?(Length@ # > ints &)], 1];


(* ::Subsubsection::Closed:: *)
(*Non-overlapping pairs*)


(* ::Text:: *)
(*NonOverlappingPairs returns all possible ways of splitting a set into an unordered set non-overlapping unordered pairs. For a set of size n there are (n-1)!! such splittings. There are no ways for set of an odd size. *)


NonOverlappingPairs@ {}:= {{}};
NonOverlappingPairs@ set_List:= Block[{len= Length@ set},
	(*If[OddQ@ len, Abort[]];*)
	Flatten[Table[
		Join[{set[[{1, n}]]}, #]&/@ NonOverlappingPairs@ Join[set[[2;; n-1]], set[[n+1;; len]]]
	,{n, 2, len}], 1]
]


(* ::Subsection:: *)
(*CG algebra simplification*)


(* ::Subsubsection::Closed:: *)
(*BetterContractCGs*)


(* ::Text:: *)
(*A minor improvement of ContractCGs: It removes the contraction of symmetric objects with antisymmetric structure factors before using ContractCGs, and optimizes the contractions of repeated CG structures*)


BetterContractCGs[expr_]:=Module[{out=LagrangianExpand@expr,termNoCgs,cgsList,termsCgs},
	out= out//ContractDelta//RelabelIndices;
	out= out//.CG[fStruct[_],{_,A_,B_}|{A_,B_,_}|{A_,_,B_}]Field[W_,Scalar,{A_},d_]Field[W_,Scalar,{B_},d_]:>0;
	
	termNoCgs= out/._CG->0;
	cgsList= DeleteCases[DeleteDuplicates[Times@@Cases[#,_CG,All]&/@(TermsToList[out])],1];
	termsCgs= Coefficient[out,cgsList]/._CG->0;
	cgsList= ContractCGs/@cgsList;
	
	out= termsCgs . cgsList+termNoCgs//ContractDelta
]


(* ::Subsection:: *)
(*Other*)


(* ::Subsubsection::Closed:: *)
(*PseudoTimes*)


(* ::Text:: *)
(*A Times-like head to expand out powers *)


SetAttributes[PseudoTimes, {Orderless}];
PseudoTimes@ expr_Plus:= PseudoTimes/@ expr;
PseudoTimes@ expr_Times:= PseudoTimes@@ expr;
PseudoTimes[a___, PseudoTimes@ b___]:= PseudoTimes[a, b]
PseudoTimes[a___, n_Integer]:= n PseudoTimes@ a;
PseudoTimes[a___, b_Plus]:= PseudoTimes[a, #]&/@ b;
PseudoTimes[a___, Power[b_, n_Integer/; n > 1]]:= PseudoTimes[a, Sequence@@ ConstantArray[b, n]];


ReleasePseudoTimes@ expr_:= expr/. PseudoTimes-> Times;


(* ::Subsubsection::Closed:: *)
(*Select and delete cases in an expression *)


SelectAndDeleteCases[expr_, rule:(Rule|RuleDelayed)[lhs_, _], args___]:=
	{Cases[expr, rule, args], DeleteCases[expr, lhs, args]};
SelectAndDeleteCases[expr_, args__]:= {Cases[expr, args], DeleteCases[expr, args]};


(* ::Subsubsection::Closed:: *)
(*TermsToList*)


(* ::Text:: *)
(*Transform a sum of terms into a list, or convert a single term into a list*)


Options@ TermsToList= {Expand -> True};


TermsToList[expr_, OptionsPattern[]]:= Module[{temp= If[OptionValue@Expand, LagrangianExpand@expr, expr]},
	If[Head@ temp === Plus, List@@ temp, List@ temp]
];


(* ::Subsubsection::Closed:: *)
(*ComplementWithDuplicates*)


ComplementWithDuplicates[a_List,b_List]:=Fold[DeleteCases[#,#2,1,1]&,a,b]


(* ::Subsubsection::Closed:: *)
(*DidMatchSwitch*)


(* ::Text:: *)
(*Use any pattern variable, e.g., (x: a|...)| (b|...) to discriminate between a pattern match with pattern(s) a or pattern(s) b. x is a NullSequence if the b pattern(s) get matched.*)


DidMatchSwitch[t_, f_][]:= f;
DidMatchSwitch[t_, f_][___]:= t;


(* ::Section:: *)
(*Generate bibliography*)


(* ::Subsubsection::Closed:: *)
(*SuggestBibliography function*)


(* ::Text:: *)
(*Function to generate a suggested bibliography to the user  *)


Options@ SuggestBibliography= {
		"References"-> Default,
		"Explanation"-> False
	};


OptionTest[SuggestBibliography, "Explanation"]= BooleanQ;
OptionTest[SuggestBibliography, "References"]= MatchQ[All| Default];


OptionMessage["References", SuggestBibliography, val_]:=
	Message[General::optexpectsval, "References", SuggestBibliography, val, "value 'All' or 'Default'"];


SuggestBibliography[OptionsPattern[]] ? OptionsCheck:= Module[{abbreviation, reasons, bib, refs},
	(*What references should be included in the bibliography*)
	refs= If[OptionValue@ "References" === All,
			Keys@ LiteratureList
		,
			Keys@ $RelevantCitations
		];

	(*What are the reasons behind the suggestions of the individual references*)
	If[OptionValue@ "Explanation",
		reasons= Table[
			abbreviation= StringReplace[LiteratureList@ ref, 
				RegularExpression["(^.*\\s)*@article\\{(.*),(.|\\s)*"]:> "$2"];
			abbreviation= "Reference '" <> abbreviation <> "' is suggested for the following reasons:";
			StringJoin@@ Riffle[Prepend[$RelevantCitations@ ref, abbreviation], "\n - "]
		, {ref, refs}];

		StringJoin@@ Riffle[reasons, "\n\n"]// Echo
	];

	(*Merge the references and present the result*)
	bib= StringJoin@@ Riffle[Lookup[LiteratureList, refs], "\n\n"];
	CellPrint[Cell[bib, "Output", "PageWidth"-> Infinity]];
	Button["Copy to clipboard", CopyToClipboard@ bib]
];


(* ::Subsubsection::Closed:: *)
(*Function for collecting the bibliography*)


(* ::Text:: *)
(*Function for adding references to the suggested bibliography *)


AddToBibliography[paper_String, reason_String]:= Block[{},
	If[!KeyExistsQ[$RelevantCitations, paper],
		$RelevantCitations@ paper= {};
	];

	If[FreeQ[$RelevantCitations@ paper, reason],
		AppendTo[$RelevantCitations@ paper, reason];
	];
];


$RelevantCitations= <||>;
AddToBibliography["ProofOfConcept", "Use of the core Matchete functionality"];


(* ::Subsubsection::Closed:: *)
(*Collection of all papers (Matchete and others)*)


MatchetePapers= <|
	"EvanescentTreatment"->
"%General treatment of evanescent operators in EFT matching (particularly to the SMEFT)
@article{Fuentes-Martin:2022vvu,
    author = {Fuentes-Mart\\'\\i{}n, Javier and K\\\"onig, Matthias and Pag\\`es, Julie and Thomsen, Anders Eller and Wilsch, Felix},
    title = \"{Evanescent operators in one-loop matching computations}\",
    eprint = \"2211.09144\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    reportNumber = \"MITP-22-091, TUM-HEP-1428/22, ZU-TH-48/22\",
    doi = \"10.1007/JHEP02(2023)031\",
    journal = \"JHEP\",
    volume = \"02\",
    pages = \"031\",
    year = \"2023\"
}",
	"ProofOfConcept"->
"%Introduction of the core Matchete package and v0.1
@article{Fuentes-Martin:2022jrf,
    author = {Fuentes-Mart\\'\\i{}n, Javier and K\\\"onig, Matthias and Pag\\`es, Julie and Thomsen, Anders Eller and Wilsch, Felix},
    title = \"{A proof of concept for matchete: an automated tool for matching effective theories}\",
    eprint = \"2212.04510\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    reportNumber = \"MITP-22-105, TUM-HEP-1443/22, ZU-TH-58/22\",
    doi = \"10.1140/epjc/s10052-023-11726-1\",
    journal = \"Eur. Phys. J. C\",
    volume = \"83\",
    number = \"7\",
    pages = \"662\",
    year = \"2023\"
}",
	"SuperTracer"->
"%Early implementation of the functional tools at the heart of Matchete
@article{Fuentes-Martin:2020udw,
    author = {Fuentes-Martin, Javier and K\"onig, Matthias and Pag\`es, Julie and Thomsen, Anders Eller and Wilsch, Felix},
    title = \"{SuperTracer: A Calculator of Functional Supertraces for One-Loop EFT Matching}\",
    eprint = \"2012.08506\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    reportNumber = \"MITP-20-076, TUM-HEP-1302/20, ZU-TH-54/20\",
    doi = \"10.1007/JHEP04(2021)281\",
    journal = \"JHEP\",
    volume = \"04\",
    pages = \"281\",
    year = \"2021\"
}",
	"Matchotter"->
"%Functional treatment of thermal field theory
@article{Fuentes-Martin:2026???,
    author = {Fuentes-Martin, Javier and Lopez Miras, Javier and Moreno-S\[AAcute]nchez, Adri\[AAcute]n},
    title = \"{Matchotter: An Automated Tool for Dimensional Reduction at Finite Temperature}\",
    eprint = \"2604.21972\",
    archivePrefix = \"arXiv\",
    primaryClass = \"hep-ph\",
    month = \"4\",
    year = \"2026\"
}"
|>;


LiteratureList=MatchetePapers;


(* ::Subsubsection::Closed:: *)
(*Adding papers*)


AddLiterature::usage="AddLiterature[name, info, bibtex] adds a paper labeled by name to the list of references used by SuggestBibliography[]. The argument info provides a description of the paper and bibtex should give the corresponding bibtex code."


AddLiterature[name_String, info_String, bibtex_String]:=Module[
	{entry}
	,
	entry = name->"%"<>info<>"\n"<>bibtex;
	AssociateTo[LiteratureList, entry];
	AddToBibliography[name, info]
]


(* ::Section:: *)
(*DumpSave*)


(* ::Subsubsection::Closed:: *)
(*W.I.P.*)


(*SetAttributes[SaveMatcheteSession, HoldRest]*)


(*SaveMatcheteSession[fileName_String, symbols_List:{}]:=Module[
	{
		(* For some reason averything crashes when you DumpSave GroupMagic`... *)
		builtInDefs={"Matchete`", (*"GroupMagic`",*) Hold[Global`$MatcheteVersion], NCM, Format, NiceForm},
		file,
		mySymbols,
		$DumpSave,
		tmp
	}
,
	(* determine file for saving *)
	Switch[FileExtension[fileName],
		"",   file=fileName<>".mx",
		"mx", file=fileName,
		_,    (Print["Changing file extension from \"",FileExtension[fileName],"\" to \"mx\"."];file=DirectoryName[fileName]<>FileBaseName[fileName]<>".mx")
	];

	(* keep symbols unevaluated *)
	mySymbols= Map[Hold, Hold@symbols, {2}][[1]];

	(* list of everything that must be saved *)
	mySymbols= Echo@Join[builtInDefs, mySymbols];

	(* dummy function that does not evaluate and keeps its arguments unevaluated *)
	SetAttributes[$DumpSave, HoldRest];

	(* save Matchete session *)
	With[{symb= mySymbols}, (* needds a dummy function to remove Hold *)
		tmp= $DumpSave[file,symb]/.Hold[arg_]:>arg;
	];
	Echo[tmp];
	tmp/.$DumpSave->DumpSave;

	Print["Saved the current Matchete session to the file: ", Style[file,"Code"],". \[Rule] Symbols included: ",HoldForm[symbols],"."
	];
	Print["It can be loaded again at a later point using: ",Style["Get["<>file<>"];","Code"],". This overrides the definitions of any Matchete session that might be active at this later point."];
]*)


(* ::Section:: *)
(*Deprecated symbols*)


(* ::Subsubsection::Closed:: *)
(*Functions/symbols that used to be used *)


DefineGroupRepresentation[___]:= Block[{},
	Message[DefineGroupRepresentation::usage];
	Abort[];
]


DefineGroup[___]:= Block[{},
	Message[DefineGroup::usage];
	Abort[];
]
