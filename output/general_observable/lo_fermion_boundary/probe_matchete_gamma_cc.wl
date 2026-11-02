(* Read-only v2 diagnostic: print the installed Matchete GammaCC conventions and
   nearby source definitions. Redirect STDOUT into this worker's output folder.
   Does not modify Matchete, the T3 repository or any released boundary. *)

Print["T3_MATCHeTE_GAMMACC_DIAGNOSTIC_V2"];
Print["WOLFRAM_VERSION=", $Version];
Print["WOLFRAM_USER_APPLICATIONS=",
  FileNameJoin[{$UserBaseDirectory, "Applications"}]];
Print["WOLFRAM_SYSTEM_APPLICATIONS=",
  FileNameJoin[{$BaseDirectory, "Applications"}]];
Print["MATCHeTE_FIND_BEFORE=", ToString[InputForm[
  Quiet[Check[FindFile["Matchete`"], $Failed]]]]];

(* Mirror the frontend-aware initialization in the project's successful
   Lagrangian/RunGeneralT3FermionFirstEFT1Tree.wl.
   Deliberately do not suppress loading messages: they diagnose missing
   dependencies and FrontEnd errors. A message alone does not imply failure. *)
load = CheckAbort[
  UsingFrontEnd[Needs["Matchete`"]; True],
  $Aborted
];
file = Quiet[Check[FindFile["Matchete`"], $Failed]];
gammaNames = Names["*`GammaCC"];
packageNames = Names["Matchete`*"];
Print["FRONTEND_LOAD_RESULT=", ToString[InputForm[load]]];
Print["GAMMACC_SYMBOLS_AFTER_LOAD=", ToString[InputForm[gammaNames]]];
Print["MATCHeTE_PUBLIC_SYMBOL_COUNT=", Length[packageNames]];
Print["MATCHeTE_FILE=", ToString[InputForm[file]]];
If[load === $Aborted || gammaNames === {} || Length[packageNames] === 0,
  Print["ERROR=Matchete package initialization did not expose its GammaCC symbol"];
  Print["STATUS=PackageUnavailable"];
  Print["NEXT=Check preceding Needs/UsingFrontEnd messages and the installed Matchete application path"];
  Exit[2]
];
Print["MATCHeTE_VERSION=", ToString[InputForm[Quiet[Check[Matchete`$MatcheteVersion, Missing["Unavailable"]]]]]];

names = DeleteDuplicates[Join[
  Names["*`GammaCC"], Names["*`CC"], Names["*`CConj"], Names["*`Proj"]
]];
Print["MATCHING_SYMBOL_NAMES=", ToString[InputForm[names]]];
Do[
  With[{s = Quiet[Check[ToExpression[name], $Failed]]},
    Print["SYMBOL=", name];
    Print["ATTRIBUTES=", ToString[InputForm[Quiet[Check[Attributes[s], $Failed]]]]];
    Print["OWNVALUES=", ToString[InputForm[Quiet[Check[OwnValues[s], $Failed]]]]];
    Print["DOWNVALUES=", ToString[InputForm[Quiet[Check[DownValues[s], $Failed]]]]];
    Print["UPVALUES=", ToString[InputForm[Quiet[Check[UpValues[s], $Failed]]]]];
  ],
  {name, names}
];

If[StringQ[file] && FileExistsQ[file],
  root = DirectoryName[file];
  If[FileNameTake[root] === "Kernel", root = DirectoryName[root]];
  paths = Quiet[Check[FileNames[{"*.m", "*.wl"}, root, Infinity], {}]];
  Print["PACKAGE_SOURCE_ROOT=", root];
  Print["PACKAGE_SOURCE_FILES_SCANNED=", Length[paths]];
  total = 0;
  Do[
    content = Quiet[Check[Import[path, "Text"], ""]];
    If[StringQ[content] && content =!= "",
      occurrences = StringPosition[content, "GammaCC"];
      Do[
        If[total >= 70, Break[]];
        Print["SOURCE_FILE=", path];
        Print["GAMMACC_SOURCE_EXCERPT=", StringTake[content,
          {Max[1, pos[[1]]-220], Min[StringLength[content], pos[[2]]+320]}]];
        total++,
        {pos, Take[occurrences, UpTo[10]]}
      ];
    ];
    If[total >= 70, Break[]],
    {path, paths}
  ];
  Print["EXCERPTS_PRINTED=", total];
];
Print["STATUS=EvidenceOnly_NoPhysicalSignInferred"];
