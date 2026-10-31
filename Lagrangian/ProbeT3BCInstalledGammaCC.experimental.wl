(* Read-only Matchete GammaCC source locator, v10.1.
   Unlike v10, this does NOT require Needs["Matchete`"] to succeed:
   the source is found using FindFile / local PacletFind / $Path first.
   An optional front-end-aware Needs is used ONLY if source discovery fails.
   No EFT model is loaded, no Match is run, and no operator factor is inferred. *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[!MemberQ[{1, 2}, Length[args]],
  Print["Usage: wolframscript -file Lagrangian/ProbeT3BCInstalledGammaCC.experimental.wl OUTPUT.json [MATCHeTE_SOURCE_ROOT]"];
  Exit[2]
];
outPath = ExpandFileName[First[args]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath], CreateIntermediateDirectories -> True]
];
fmt = "t3-bc-installed-gammacc-source-probe-v1";
writeReport[status_, payload_, code_] := Module[{report, exported},
  report = Join[<|
    "format" -> fmt, "status" -> status,
    "Mathematica_version" -> ToString[$Version],
    "source_only" -> True, "matched_any_EFT" -> False,
    "absolute_GammaCC_to_Warsaw_factor_certified" -> False,
    "physical_C5" -> Null
  |>, payload];
  exported = Quiet @ Check[Export[outPath, report, "RawJSON"], $Failed];
  If[exported === $Failed,
    Print["T3_GAMMACC: JSON export failed: ", outPath]; Exit[8]
  ];
  Print["T3_GAMMACC: ", status, "; output=", outPath];
  Exit[code]
];

(* Limit automatic discovery to Matchete-specific directories.  In
   particular, do not recursively scan the entire Applications directory. *)
packageRootQ[path_] := StringQ[path] && DirectoryQ[path] &&
  StringContainsQ[ToLowerCase[FileNameTake[path]], "matchete"];
parentCandidates[file_] := Module[{dirs},
  If[!StringQ[file] || !FileExistsQ[file], Return[{}]];
  dirs = NestList[DirectoryName, DirectoryName[ExpandFileName[file]], 4];
  Select[dirs, packageRootQ]
];
pathCandidates[] := Module[{ps = Select[$Path, StringQ]},
  Select[Flatten[{
    ps,
    FileNameJoin[{#, "Matchete"}] & /@ ps
  }], packageRootQ]
];
pacletCandidates[] := Module[{found, dirs},
  (* PacletFind is a local metadata lookup, not a network fetch.  It is
     optional because older kernels need not expose the same paclet API. *)
  found = Quiet @ Check[PacletFind["Matchete"], {}];
  If[!ListQ[found], Return[{}]];
  dirs = Quiet @ Check[(#["Location"] & /@ found), {}];
  Select[dirs, packageRootQ]
];
sourceFiles[root_] := DeleteDuplicates @ Join[
  FileNames["*.wl", root, Infinity],
  FileNames["*.m", root, Infinity]
];
validRoots[candidates_] := Select[
  DeleteDuplicates[Select[candidates, packageRootQ]],
  Function[dir, Length[sourceFiles[dir]] > 0]
];

explicitRoot = If[Length[args] == 2, ExpandFileName[args[[2]]], None];
If[explicitRoot =!= None && !DirectoryQ[explicitRoot],
  writeReport["ExplicitSourceRootDoesNotExist", <|
    "specified_source_root" -> explicitRoot|>, 3]
];
Print["T3_GAMMACC_STAGE: locate-installed-source-without-loading-package"];
findFile = Quiet @ Check[FindFile["Matchete`"], $Failed];
roots = If[explicitRoot =!= None,
  {explicitRoot},
  Join[parentCandidates[findFile], pacletCandidates[], pathCandidates[]]
];
If[explicitRoot === None, roots = validRoots[roots]];
loadAttempted = False;
loadSucceeded = False;

(* The project's existing successful Matchete initialization wraps Needs
   in UsingFrontEnd.  Only invoke it if passive discovery found nothing. *)
If[roots === {} && explicitRoot === None,
  Print["T3_GAMMACC_STAGE: source-not-visible; frontend-aware-discovery-fallback"];
  loadAttempted = True;
  load = CheckAbort[
    Quiet @ Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed],
    $Aborted
  ];
  loadSucceeded = (load === True);
  findFile = Quiet @ Check[FindFile["Matchete`"], $Failed];
  pathVar = Quiet @ Check[Matchete`$MatchetePath, $Failed];
  roots = validRoots[Join[
    parentCandidates[findFile],
    If[StringQ[pathVar], {pathVar}, {}],
    pacletCandidates[], pathCandidates[]
  ]];
];
If[roots === {}, writeReport["MatcheteSourceRootNotFound", <|
  "FindFile_Matchete" -> ToString[findFile, InputForm],
  "frontend_loader_attempted" -> loadAttempted,
  "frontend_loader_succeeded" -> loadSucceeded,
  "hint" -> "Pass the exact installed Matchete package source folder as the optional second argument; it must contain .wl or .m source files."
|>, 4]];

(* Use the most complete Matchete-specific source root; explicit roots
   may have arbitrary names.  Always reject pathological source trees. *)
rootData = Table[<|"root" -> r, "files" -> sourceFiles[r]|>, {r, roots}];
rootData = Select[rootData, Length[# ["files"]] > 0 &];
If[rootData === {}, writeReport["MatcheteWLSourceNotFound", <|
  "candidate_roots" -> roots|>, 5]];
rootData = SortBy[rootData, -Length[# ["files"]] &];
root = First[rootData]["root"];
allFiles = First[rootData]["files"];
If[Length[allFiles] > 1500,
  writeReport["UnexpectedlyLargeSourceTree", <|
    "source_root" -> root, "files_found" -> Length[allFiles]|>, 6]
];
files = Select[allFiles, FileByteCount[#] <= 2000000 &];
Print["T3_GAMMACC_STAGE: scanning ", Length[files],
  " files under ", root];

needles = {"GammaCC", "CConj", "Transp", "DiracProduct", "Proj"};
extractHits[content_, needle_] := Module[{pos},
  pos = Take[StringPosition[content, needle], UpTo[20]];
  Map[Function[place,
    With[{start = Max[1, place[[1]] - 130],
      stop = Min[StringLength[content], place[[2]] + 160]},
      <|"start_offset" -> place[[1]],
        "excerpt" -> StringTake[content, {start, stop}]|>
    ]], pos]
];
(* Import every source file at most once.  Keep GammaCC hits first so a
   generous count of other symbols cannot truncate the relevant evidence. *)
allFound = Reap[
  Do[
    content = Quiet @ Check[Import[file, "Text"], $Failed];
    If[StringQ[content],
      Do[
        hits = extractHits[content, needle];
        If[hits =!= {}, Sow[<|
          "file" -> StringReplace[file, root -> "."],
          "symbol" -> needle,
          "source_snippets" -> hits
        |>]],
        {needle, needles}
      ]
    ],
    {file, files}
  ]
][[2]];
allFound = If[allFound === {}, {}, First[allFound]];
gamma = Select[allFound, # ["symbol"] === "GammaCC" &];
nonGamma = Select[allFound, # ["symbol"] =!= "GammaCC" &];
found = Join[Take[gamma, UpTo[100]],
  Take[nonGamma, UpTo[Max[0, 150 - Min[100, Length[gamma]]]]]];
Print["T3_GAMMACC: GammaCC source entries=", Length[gamma]];
If[gamma === {}, writeReport["GammaCCSourceNotLocated", <|
  "source_root" -> root,
  "files_scanned" -> Length[files],
  "frontend_loader_attempted" -> loadAttempted,
  "frontend_loader_succeeded" -> loadSucceeded,
  "source_occurrences" -> found
|>, 7]];
writeReport["InstalledGammaCCSourceLocated", <|
  "source_root" -> root,
  "files_scanned" -> Length[files],
  "frontend_loader_attempted" -> loadAttempted,
  "frontend_loader_succeeded" -> loadSucceeded,
  "GammaCC_source_files" -> DeleteDuplicates[Lookup[gamma, "file"]],
  "source_occurrences" -> found,
  "evidence_requires_human_symbol_definition_review" -> True,
  "native_NCM_GammaCC_operator_equals_Warsaw_Q_certified" -> False,
  "hbar_convention_external_reference" ->
    "https://matchete.gitlab.io/quickstart.html"
|>, 0];
