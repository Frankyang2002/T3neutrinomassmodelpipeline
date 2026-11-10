(* ::Package:: *)

(* ::Title:: *)
(*Initialization [Matchete`]*)


(* Generate warning for Mathematica versions earlier than 11.0.0 *)
If[$VersionNumber < 11.0, 
CellPrint[{
TextCell["Matchete` was developed for Mathematica 11.0.0 and later. 
Your current Mathematica version [" <> ToString@$Version <> "] might not be compatible. 
In case you are experiencing problems with Matchete`, please update your Mathematica version.",
"Text",Background->LightRed]
}]
]


(*Fix notation input for NCM*)
(*Notation`AutoLoadNotationPalette= False;
Needs["Notation`"]*)
(*InfixNotation[ParsedBoxWrapper["**"], Matchete`NCM];*)
$PreRead= ReplaceAll[#, "**"-> "\[CenterDot]"]&;


(*Matchete Logo is now dynamic and changes if you use ThermalFieldTheory*)
Dynamic[Global`$logoMatchete];
Global`$logoMatchete= Import[FileNameJoin[{DirectoryName[$InputFileName],"LogoMatchete.eps"}]];


(* Load the package *)
Module[{packageDirectory},
	(* Avoid double loading the package *)
	If[MemberQ[$Packages,"Matchete`"],
		Print[Style["The package Matchete` is already loaded. Please restart the kernel to reload it.",RGBColor[.8,.4706,0.2573]]]
	];
	
	(* Set directory of Matchete package *)
	packageDirectory= DirectoryName[$InputFileName, 2]<>"/Package";
	
	(*Matchete version*)
	$MatcheteVersion= StringTrim[First@ StringSplit[
		Import@ FileNameJoin[{DirectoryName@ packageDirectory, "version"}], "\n"] ];
	
	(*The whole package contents are loaded here*)
	Check[
		(*Ensure Bar updates as first thing*)
		Get[FileNameJoin[{packageDirectory, "NCM.m"}]];
		,
		Print[Style["Loading failed!",RGBColor[.6,.0706,0.1373]]];
		(*Abort[]*) (* This call of Abort[] causes issues when using Matchete without the Wolfram GUI. The warning above should be enough. *)
	];

	(*Print the logo*)
	CellPrint@ ExpressionCell[
		Row[{(*Dynamic Logo*)
			DynamicModule[{logo= Dynamic@Style[Global`$logoMatchete, Magnification-> 0.5]}, 
                Dynamic[logo]
                ,
                Initialization :> (If[!ValueQ[logo],   logo = Style[Global`$logoMatchete, Magnification -> 0.5] ])
                ,
                UnsavedVariables :> {},
                SaveDefinitions -> True
            ]
			,
			Row[{TextCell["v"<>$MatcheteVersion, Bold, 
				FontColor-> RGBColor[0, 0.3168,0.5928], 
				FontSize-> 15,
				FontFamily-> "Source Code Pro"]}, BaselinePosition->Scaled@.5]
		}, Spacer[7]]
	, CellMargins->{{70,5},{5,5}} ];

	Print[
		"by Javier Fuentes-Mart\[IAcute]n, Matthias K\[ODoubleDot]nig, Julie Pag\[EGrave]s, Anders Eller Thomsen, and Felix Wilsch \n",
		"Contributors:   Lukas Born, Luis Hourtz, Tamilarasan Ketheeswaran, Javier L\[OAcute]pez Miras, Adri\[AAcute]n Moreno-S\[AAcute]nchez, and Philip Weber \n",
		"Reference:      ", Hyperlink["arXiv:2212.04510", "https://arxiv.org/abs/2212.04510"],"\n",
		"Documentation:  ", Hyperlink["open", "paclet:Matchete/guide/Matchete"],"\n",
		"Website:        ", Hyperlink["https://matchete.gitlab.io","https://matchete.gitlab.io"],"\n",
		"Model database: ", Hyperlink["https://gitlab.com/matchete/model-database","https://gitlab.com/matchete/model-database"],"\n",
		"Repository:     ", Hyperlink["https://gitlab.com/matchete/matchete","https://gitlab.com/matchete/matchete"]
	];
];


Matchete`Matchete`PackagePrivate`AddAutoCompletion["LoadModel"][
	FileBaseName/@ FileNameTake/@ FileNames["*.m", FileNameJoin[{$MatchetePath, "Models"}], 2]];


(* Protect the symbols used in the package *)
(*
SetAttributes[
  Evaluate @ Flatten[Names /@ {"Matchete`*"}],
  {Protected, ReadProtected}
]
*)
