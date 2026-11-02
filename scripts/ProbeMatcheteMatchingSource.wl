(* Workstream 1: read-only diagnostic for the installed Matchete matching
   implementation. No Match[], no heavy-field integration, no EFT changes.
   This probe CANNOT certify how the already exported one-loop MatchedEFT
   handles hard/soft regions. Inspect the source evidence independently.

   wolframscript -file scripts/ProbeMatcheteMatchingSource.wl
*)
ClearAll["Global`*"];
file = ExpandFileName[$InputFileName];
root = DirectoryName[DirectoryName[file]];
out = FileNameJoin[{root, "output", "general_observable", "lo_certification_spinor"}];
If[!DirectoryQ[out], CreateDirectory[out, CreateIntermediateDirectories -> True]];
path = FileNameJoin[{out, "matchete_match_semantics_" <>
  StringReplace[CreateUUID[], "-" -> ""] <> ".json"}];
format = "t3-matchete-hard-soft-source-markers-v1";
ClearAll[finish];
finish[data_Association, exit_Integer] := Module[{rec},
  rec = Join[<|
      "format" -> format,
      "no_matching_executed" -> True,
      "previous_boundaries_unchanged" -> True,
      "hard_soft_semantics_physically_certified" -> False,
      "EFT1_soft_loop_executed" -> False,
      "native_matchete_Weyl_phase_certified" -> False,
      "physical_C5" -> Null|>, data];
  Export[path, rec, "RawJSON"];
  Print["MATCH_SOURCE_DIAGNOSTIC=", path];
  Print["STATUS=", Lookup[rec, "status", "Unknown"]];
  Exit[exit]
];
load = CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted];
If[load =!= True,
  finish[<|"status" -> "Blocked", "reason" -> "Matchete frontend load failed"|>, 2]];
matcheteFile = Quiet[Check[FindFile["Matchete`"], $Failed]];
If[!StringQ[matcheteFile] || !FileExistsQ[matcheteFile],
  finish[<|"status" -> "Blocked", "reason" -> "Matchete package file not located"|>, 2]];
packageRoot = DirectoryName[matcheteFile];
If[FileNameTake[packageRoot] === "Kernel", packageRoot = DirectoryName[packageRoot]];
(* Target the actual installed package directory, not any project outputs. *)
files = Sort[DeleteDuplicates[FileNames[{"*.m", "*.wl"}, packageRoot, Infinity]]];
truncated = Length[files] > 400;
files = Take[files, UpTo[400]];
searchTerms = {"hard region", "hardregion", "soft region", "softregion",
  "subtract", "matching", "looporder", "counterterm", "region", "eom"};
excerpts = {};
scanned = 0;
Do[
  If[FileByteCount[f] <= 2000000,
    content = Quiet[Check[Import[f, "Text"], $Failed]];
    If[StringQ[content],
      scanned++;
      lower = ToLowerCase[content];
      Do[
        positions = Take[StringPosition[lower, term], UpTo[2]];
        Do[
          If[Length[excerpts] < 100,
            AppendTo[excerpts, <|
              "file_relative_to_installed_package" -> FileNameDrop[f,
                FileNameDepth[packageRoot]],
              "search_term" -> term,
              "source_excerpt" -> StringTake[content, {
                Max[1, p[[1]] - 180], Min[StringLength[content], p[[2]] + 280]}]
              |>]],
          {p, positions}],
        {term, searchTerms}]
    ]
  ], {f, files}];
matchSymbols = Select[Names["*`Match"], StringQ];
(* Symbol definitions are diagnostics only: no Match[...] call is made. *)
definitions = Map[Function[name, Module[{sym, values, txt},
  sym = Quiet[Check[ToExpression[name], $Failed]];
  values = If[sym === $Failed, {}, Quiet[Check[DownValues[sym], {}]]];
  txt = ToString[InputForm[values]];
  <|"name" -> name, "downvalues_count" -> Length[values],
    "downvalues_sha256" -> Hash[txt, "SHA256", "HexString"],
    "downvalues_excerpt" -> StringTake[txt, UpTo[2500]]|>
]], matchSymbols];
finish[<|
  "status" -> "NativeMatcheteMatchingSourceMarkersCollected_NotSemanticsCertified",
  "installed_matchete_main_file" -> matcheteFile,
  "package_root" -> packageRoot,
  "source_file_count_scanned" -> scanned,
  "source_file_list_truncated" -> truncated,
  "matching_symbol_downvalues" -> definitions,
  "source_excerpts" -> excerpts,
  "interpretation_warning" -> "Text matches, even if they mention hard or soft, are NOT proof of Match hard/soft subtraction. Manually inspect actual implementation and EFT one-loop counterterms."
|>, 0];
