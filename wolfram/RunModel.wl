(* RunModel_TwoScalar.wl
   CLI: output-dir EFT-order loop-order T3-class alpha
   Example: wolframscript -file RunModel.wl output 5 1 B -1
*)
ClearAll["Global`*"];
scriptDirectory=DirectoryName@ExpandFileName[$InputFileName];
args=Rest[$ScriptCommandLine];
outputDirectory=If[Length[args]>=1,ExpandFileName@args[[1]],FileNameJoin@{scriptDirectory,"output","T3"}];
eftOrder=If[Length[args]>=2,ToExpression@args[[2]],5];
loopOrder=If[Length[args]>=3,ToExpression@args[[3]],1];
modelClass=If[Length[args]>=4,ToUpperCase@args[[4]],"B"];
ParseSignedCLI[s_String] := Which[
  StringMatchQ[s, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[s, 1]],
  StringMatchQ[s, "p" ~~ DigitCharacter ..],  ToExpression[StringDrop[s, 1]],
  True, ToExpression[s]
];
alpha=If[Length[args]>=5,ParseSignedCLI[args[[5]]],-1];
If[!DirectoryQ[outputDirectory],CreateDirectory[outputDirectory,CreateIntermediateDirectories->True]];
Print["Loading Matchete..."];
matcheteLoadResult = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoadResult],
  Print["ERROR: Matchete failed to load."];
  Exit[6]
];
Print["Matchete loaded successfully."];
Get[FileNameJoin@{scriptDirectory,"PhysicsLaTeX.wl"}];
Get[FileNameJoin@{scriptDirectory,"T3ModelCatalog.wl"}];
Get[FileNameJoin@{scriptDirectory,"LagrangianBuilder.wl"}];
Get[FileNameJoin@{scriptDirectory,"RunMatching.wl"}];
model=T3ModelFromClass[modelClass,alpha];
If[model===$Failed,Print["ERROR: unknown T3 class."];Exit[2]];
Print["Model: ",model["Class"],", alpha = ",alpha];
Print["Scalar1 (d,Y) = ", {model["Scalar1","SU2"],model["Scalar1","Y"]}];
Print["Scalar2 (d,Y) = ", {model["Scalar2","SU2"],model["Scalar2","Y"]}];
Print["Fermion (d,Y) = ", {model["Fermion","SU2"],model["Fermion","Y"]}];
result=CheckAbort[BuildT3Lagrangian[model],$Aborted];
If[!AssociationQ[result]||result===$Aborted||result===$Failed,Print["ERROR: build failed."];Exit[10]];
If[result["Status"] =!= "Success",Print["ERROR: full UV validation failed."];Exit[11]];
Print["Accepted interactions: ",result["AllowedInteractions"]];
Print["T3 ingredients present: ",result["T3IngredientsPresent"]];
matching=CheckAbort[RunT3Matching[result["LUV"],eftOrder,loopOrder],$Aborted];
If[!AssociationQ[matching]||Lookup[matching,"Status",""]=!="Success",Print["ERROR: matching failed."];Exit[12]];
matched=Lookup[matching,"MatchedEFT",Lookup[matching,"LoopEFT",0]];
uv=ExpressionToLaTeX[result["LUV"]]; bsm=ExpressionToLaTeX[result["LBSM"]]; eft=ExpressionToLaTeX[matched];

(* Detect Weinberg robustly, then separately isolate its raw GammaCC contributions. *)
weinbergData=ExtractWeinbergCoefficient[matched];
weinberg=TrueQ[Lookup[weinbergData,"Present",False]];
weinbergSector=Lookup[weinbergData,"Sector",0];
holomorphicWeinbergSector=Lookup[weinbergData,"HolomorphicSector",0];
weinbergCoefficient=Lookup[weinbergData,"Coefficient",Missing["PendingCanonicalisation"]];
weinbergTerms=Lookup[weinbergData,"Terms",{}];
holomorphicWeinbergTerms=Lookup[weinbergData,"HolomorphicTerms",{}];
conjugateWeinbergTerms=Lookup[weinbergData,"ConjugateTerms",{}];
weinbergSectorTeX=ExpressionToLaTeX[weinbergSector];
holomorphicWeinbergSectorTeX=ExpressionToLaTeX[holomorphicWeinbergSector];
weinbergCoefficientTeX=If[MissingQ[weinbergCoefficient],
  <|"Success"->False,"LaTeX"->""|>,
  ExpressionToLaTeX[weinbergCoefficient]
];
summary=<|
 "ModelClass"->model["Class"],"Alpha"->alpha,
 "Scalar1SU2"->model["Scalar1","SU2"],"Scalar1Hypercharge"->ToString@InputForm@model["Scalar1","Y"],
 "Scalar2SU2"->model["Scalar2","SU2"],"Scalar2Hypercharge"->ToString@InputForm@model["Scalar2","Y"],
 "FermionSU2"->model["Fermion","SU2"],"FermionHypercharge"->ToString@InputForm@model["Fermion","Y"],
 "BuildStatus"->result["Status"],"MatchingStatus"->matching["Status"],
 "AcceptedInteractions"->result["AllowedInteractions"],"RejectedInteractions"->result["RejectedInteractions"],
 "T3IngredientsPresent"->result["T3IngredientsPresent"],"WeinbergOperatorPresent"->TrueQ[weinberg],
 "WeinbergExtractionStatus"->Lookup[weinbergData,"Status","Unknown"],
 "WeinbergTermCount"->Lookup[weinbergData,"TermCount",0],
 "WeinbergHolomorphicTermCount"->Lookup[weinbergData,"HolomorphicTermCount",0],
 "WeinbergConjugateTermCount"->Lookup[weinbergData,"ConjugateTermCount",0],
 "WeinbergRawFile"->If[Length[weinbergTerms]>0,"c5_raw.txt",""],
 "WeinbergRawLaTeXFile"->If[Length[weinbergTerms]>0,"c5_raw.tex",""],
 "WeinbergCoefficientFile"->If[!MissingQ[weinbergCoefficient],"c5_coefficient.txt",""],
 "WeinbergCoefficientLaTeXFile"->If[TrueQ[Lookup[weinbergCoefficientTeX,"Success",False]],"c5_coefficient.tex",""],
 "WeinbergCoefficientInputForm"->If[MissingQ[weinbergCoefficient],
    "Pending exact Matchete operator canonicalisation",
    ToString[weinbergCoefficient,InputForm]],
 "WeinbergCoefficientLaTeX"->Lookup[weinbergCoefficientTeX,"LaTeX",""],
 "WeinbergSectorConversionSuccess"->TrueQ[Lookup[weinbergSectorTeX,"Success",False]],
 "WeinbergCoefficientConversionSuccess"->TrueQ[Lookup[weinbergCoefficientTeX,"Success",False]],
 "UVConversionSuccess"->TrueQ[uv["Success"]],"BSMUVConversionSuccess"->TrueQ[bsm["Success"]],"EFTConversionSuccess"->TrueQ[eft["Success"]],
 "UVLagrangianLaTeX"->uv["LaTeX"],"BSMUVLagrangianLaTeX"->bsm["LaTeX"],"EFTLagrangianLaTeX"->eft["LaTeX"]
|>;
(* Detailed Weinberg diagnostics belong in files, not normal terminal output. *)
If[Length[weinbergTerms]>0,
  Export[
    FileNameJoin@{outputDirectory,"c5_raw.txt"},
    StringRiffle[
      MapIndexed[
        "--- Weinberg term "<>ToString[First[#2]]<>" ---\n"<>ToString[#1,InputForm]&,
        weinbergTerms
      ],
      "\n\n"
    ],
    "Text"
  ];
  Export[
    FileNameJoin@{outputDirectory,"c5_raw.tex"},
    Lookup[weinbergSectorTeX,"LaTeX",""],
    "Text"
  ];
];


If[!MissingQ[weinbergCoefficient],
  Export[
    FileNameJoin@{outputDirectory,"c5_coefficient.txt"},
    ToString[weinbergCoefficient,InputForm],
    "Text"
  ];
  If[TrueQ[Lookup[weinbergCoefficientTeX,"Success",False]],
    Export[
      FileNameJoin@{outputDirectory,"c5_coefficient.tex"},
      Lookup[weinbergCoefficientTeX,"LaTeX",""],
      "Text"
    ];
  ];
];

Export[FileNameJoin@{outputDirectory,"comparison_summary.json"},summary,"RawJSON"];
Print["Weinberg operator present: ",summary["WeinbergOperatorPresent"]];
If[weinberg,
  If[summary["WeinbergTermCount"] > 0,
    Print["Weinberg contributions isolated: ",summary["WeinbergTermCount"]],
    Print["Weinberg present; contribution isolation pending"]
  ];
  If[Lookup[weinbergData,"Status",""] === "Success",
    Print["C5 extraction: Success (holomorphic terms: ",
      Lookup[weinbergData,"HolomorphicTermCount",0],
      "; HC terms: ",Lookup[weinbergData,"ConjugateTermCount",0],")"],
    Print["C5 canonicalisation: pending"]
  ];
];
Exit[0];
