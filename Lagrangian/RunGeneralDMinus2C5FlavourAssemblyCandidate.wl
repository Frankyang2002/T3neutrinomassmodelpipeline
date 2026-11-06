(*
  General D(-2): read-only, explicitly PROVISIONAL flavour-mapped C5 assembly.

  Inputs: validated saved O(hbar) full-EFT1 heavy-scalar matched EFT, saved
  first-fermion EFT1 O(hbar) correction, and direct running insertion.
  Does NOT run Match, modify WXF, impose scalar mass relations, or assert
  a certified operator normalization.  Use the convention documented in
  RunMatching.wl: C5[p,q] = A5[p,q] + A5[q,p].

  Wolfram invocation:
    wolframscript -file Lagrangian/RunGeneralDMinus2C5FlavourAssemblyCandidate.wl \
      FIRST_THRESHOLD_OUTPUT_DIR DIRECT_RUNNING_INSERTION.wl
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2, Print["Usage: OUTPUT_DIR RUNNING_INSERTION.wl"]; Exit[2]];
outputDir = ExpandFileName[args[[1]]];
insertionPath = ExpandFileName[args[[2]]];
dataDir = FileNameJoin[{outputDir, "data"}];
profileDir = FileNameJoin[{dataDir, "general_d_m2_scalar_threshold", "profile_tree"}];
CreateDirectory[profileDir, CreateIntermediateDirectories -> True];
reportPath = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_result.json"}];
logPrefix = "[C5 flavour candidate] ";
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];

ClearAll[D2Text, D2Fail, D2Terms, D2HasHeavy, D2FlavourNames,
  D2LeptonFlavourNames, D2MapFlavour, D2ValidCoefficientQ, D2ExportExpression];
D2Text[x_] := ToString[x, InputForm];
D2Fail[code_Integer, reason_String, details_: <||>] := (
  Print[logPrefix, "ERROR: ", reason];
  Export[reportPath, Join[<|"Status" -> "Failed", "Reason" -> reason,
    "ExitCode" -> code, "NoNewMatchPerformed" -> True,
    "FixedOrderSMEFTCertified" -> False|>, details], "RawJSON"];
  Exit[code]
);
D2Terms[x_] := Module[{v = Expand[x]},
  If[v === 0, {}, If[Head[v] === Plus, List @@ v, {v}]]
];
D2HasHeavy[x_] := Module[{s = D2Text[x]},
  AnyTrue[{"NewScalar1", "T3S2R", "T3S2I"},
    StringContainsQ[s, "Field[" <> # <> ","] &]
];
(* All flavour labels, including those occurring inside Coupling objects.
   Every retained operator must refer to the same two internal labels. *)
D2FlavourNames[x_] := DeleteDuplicates @ Cases[x,
  Matchete`Index[i_Symbol, t_Symbol] /; SymbolName[t] === "Flavor" :>
    SymbolName[i], Infinity];
(* Only the flavour indices actually on lepton doublet fields.  The order
   is recorded separately, since it may matter to the operator convention. *)
D2LeptonFlavourNames[x_] := Flatten @ Cases[x,
  Matchete`Field[f_Symbol, t_Symbol, ind_, deriv_] /;
    SymbolName[f] === "l" && SymbolName[t] === "Fermion" :>
    Cases[ind, Matchete`Index[i_Symbol, flav_Symbol] /;
      SymbolName[flav] === "Flavor" :> SymbolName[i], Infinity],
  Infinity];
(* Map only Index[...,Flavor], not gauge or Lorentz dummy indices. *)
D2MapFlavour[x_, a_Symbol, b_Symbol] := x /. {
  Matchete`Index[i_Symbol, fl_Symbol] /;
      SymbolName[fl] === "Flavor" && SymbolName[i] === "d$$1" :>
    Matchete`Index[a, fl],
  Matchete`Index[i_Symbol, fl_Symbol] /;
      SymbolName[fl] === "Flavor" && SymbolName[i] === "d$$2" :>
    Matchete`Index[b, fl]
};
D2ValidCoefficientQ[x_] := x =!= $Failed && x =!= $Aborted &&
  !MissingQ[x] && x =!= 0 && !D2HasHeavy[x];
D2ExportExpression[name_String, expr_] := Module[{root},
  root = FileNameJoin[{profileDir, "profile_c5_flavour_candidate_" <> name}];
  Export[root <> ".wxf", expr, "WXF"];
  Export[root <> ".txt", D2Text[expr], "Text"];
  root
];

required = <|
  "PreviousStructureAudit" -> FileNameJoin[{profileDir,
    "profile_c5_flavour_structure_result.json"}],
  "ComponentAudit" -> FileNameJoin[{profileDir,
    "profile_fixed_order_components_result.json"}],
  "HardHbar" -> FileNameJoin[{profileDir,
    "profile_full_scalar_one_loop_hard_hbar_sector.wxf"}],
  "FirstLoop" -> FileNameJoin[{dataDir,
    "eft1_after_F_realbasis_loop.wxf"}],
  "RunningInsertion" -> insertionPath
|>;
missing = Select[Keys[required], !FileExistsQ[required[#]] &];
If[missing =!= {}, D2Fail[3, "Missing saved inputs",
  <|"MissingKeys" -> missing|>]];
structure = Quiet @ Check[Import[required["PreviousStructureAudit"], "RawJSON"], $Failed];
components = Quiet @ Check[Import[required["ComponentAudit"], "RawJSON"], $Failed];
If[!AssociationQ[structure] || !AssociationQ[components] ||
  Lookup[structure, "Status", ""] =!= "WeinbergStructureAudited" ||
  Lookup[components, "Status", ""] =!= "ComponentsAudited" ||
  Lookup[structure, "AllCandidateTermsHaveTwoLeptonsTwoHiggs", False] =!= True ||
  Lookup[structure, "AllCandidateTermsHaveUnambiguousChirality", False] =!= True ||
  Lookup[structure, "HardC5IsLinearInHbar", False] =!= True ||
  Lookup[structure, "InheritedC5IsLinearInHbar", False] =!= True,
  D2Fail[4, "Prerequisite coefficient/operator audits not passed"]
];

Print[logPrefix, "Loading Matchete/extractor and saved WXF; NO Match..."];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], $Aborted] =!= True,
  D2Fail[5, "Matchete did not load"]
];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
{hard, firstLoop} = Quiet @ Check[Import[#, "WXF"], $Failed] & /@
  {required["HardHbar"], required["FirstLoop"]};
If[MemberQ[{hard, firstLoop}, $Failed], D2Fail[6, "Saved WXF import failed"]];

(* Because L0 has no heavy-field linear source, the O(hbar) inherited action
   contributes at S=0.  As in the preceding successful audits, never run the
   O(hbar) inherited correction through another one-loop Match. *)
inherited = Expand[Total[Select[D2Terms[firstLoop], !D2HasHeavy[#] &]]];
inheritedHbar = Expand[hbar Coefficient[inherited, hbar, 1]];
If[!TrueQ[Expand[inherited - inheritedHbar] === 0],
  D2Fail[7, "Inherited pure-SM sector not strictly linear in hbar"]
];

hardExt = CheckAbort[ExtractWeinbergCoefficient[hard], $Aborted];
inhExt = CheckAbort[ExtractWeinbergCoefficient[inheritedHbar], $Aborted];
If[!AssociationQ[hardExt] || !AssociationQ[inhExt] ||
    Lookup[hardExt, "Status", ""] =!= "Success" ||
    Lookup[inhExt, "Status", ""] =!= "Success" ||
    Lookup[hardExt, "TermCount", 0] =!= 68 ||
    Lookup[inhExt, "TermCount", 0] =!= 12 ||
    Lookup[hardExt, "HolomorphicTermCount", 0] =!= 34 ||
    Lookup[inhExt, "HolomorphicTermCount", 0] =!= 6,
  D2Fail[8, "Weinberg extractor count/status changed"]
];
hardHol = Lookup[hardExt, "HolomorphicTerms", {}];
inhHol = Lookup[inhExt, "HolomorphicTerms", {}];
allHol = Join[hardHol, inhHol];

(* Strict per-operator check: don't assume a global index inventory implies
   the same flavour pair on every term.  Save diagnostics on failure. *)
leptonPairs = D2LeptonFlavourNames /@ allHol;
badPairs = Flatten @ Position[leptonPairs,
  pair_ /; Sort[pair] =!= {"d$$1", "d$$2"}, {1}];
If[badPairs =!= {},
  D2Fail[9, "Some holomorphic LLHH terms lack two distinct external flavours",
    <|"BadTermPositions" -> badPairs,
      "BadTermLeptonLabels" -> (leptonPairs[[#]] & /@ Take[badPairs, UpTo[12]]),
      "BadTermPreviews" -> (StringTake[D2Text[allHol[[#]]], UpTo[800]] & /@
        Take[badPairs, UpTo[12]])|>]
];

hardOrdered = Lookup[hardExt, "Coefficient", $Failed];
inhOrdered = Lookup[inhExt, "Coefficient", $Failed];
If[!D2ValidCoefficientQ[hardOrdered] || !D2ValidCoefficientQ[inhOrdered] ||
  !TrueQ[Expand[hardOrdered - hbar Coefficient[hardOrdered, hbar, 1]] === 0] ||
  !TrueQ[Expand[inhOrdered - hbar Coefficient[inhOrdered, hbar, 1]] === 0],
  D2Fail[10, "Extracted ordered coefficients are invalid or not O(hbar)"]
];
(* Work with loop-factor-free A5; all final candidates carry one hbar. *)
aHard = Coefficient[Expand[hardOrdered], hbar, 1];
aInh = Coefficient[Expand[inhOrdered], hbar, 1];
allowed = Sort[{"d$$1", "d$$2"}];
labelsHard = Sort[D2FlavourNames[aHard]];
labelsInh = Sort[D2FlavourNames[aInh]];
If[labelsHard =!= allowed || labelsInh =!= allowed,
  D2Fail[11, "Ordered coefficient has unexpected flavour-index labels",
    <|"HardFlavourLabels" -> labelsHard,
      "InheritedFlavourLabels" -> labelsInh|>]
];

pName = T3D2ExternalFlavorP;
qName = T3D2ExternalFlavorQ;
hardPQ = D2MapFlavour[aHard, pName, qName];
hardQP = D2MapFlavour[aHard, qName, pName];
inhPQ = D2MapFlavour[aInh, pName, qName];
inhQP = D2MapFlavour[aInh, qName, pName];
If[AnyTrue[{hardPQ, hardQP, inhPQ, inhQP},
  MemberQ[D2FlavourNames[#], "d$$1" | "d$$2"] &],
  D2Fail[12, "Internal flavour indices survived replacement"]
];
(* From RunMatching.wl WeinbergOperatorConvention[]: C5=A5(p,q)+A5(q,p).
   Keep the convention explicit; do not apply it twice to the direct RGE. *)
hardSymNoHbar = Expand[hardPQ + hardQP];
inhSymNoHbar = Expand[inhPQ + inhQP];

Print[logPrefix, "Loading independent full-flavour running insertion..."];
Get[insertionPath];
If[DownValues[EFT1D2WilsonRunningWeinbergCoefficient] === {} ||
   !ValueQ[EFT1D2WilsonRunningHeavyInsertion] ||
   EFT1D2WilsonRunningHeavyInsertion =!= 0,
  D2Fail[13, "Direct running contract changed"]
];
p = Matchete`Index[pName, Flavor];
q = Matchete`Index[qName, Flavor];
runPQ = EFT1D2WilsonRunningWeinbergCoefficient[p, q];
runQP = EFT1D2WilsonRunningWeinbergCoefficient[q, p];
If[FreeQ[runPQ, Log] || !FreeQ[runPQ, hbar],
  D2Fail[14, "Direct running must contain Log and no explicit hbar"]
];
(* Simple swap check, not a physics proof if additional tensor conventions
   are built into the running formula.  Never silently symmetrize it here. *)
runSymmetric = TrueQ[Expand[runPQ - runQP] === 0];

candidateHard = hbar hardSymNoHbar;
candidateInherited = hbar inhSymNoHbar;
candidateRunning = hbar runPQ;
candidateTotal = Expand[candidateHard + candidateInherited + candidateRunning];
If[!TrueQ[Expand[candidateTotal - hbar Coefficient[candidateTotal, hbar, 1]] === 0] ||
  D2HasHeavy[candidateTotal],
  D2Fail[15, "Candidate sum is not one-loop-only or contains heavy fields"]
];

roots = <|
  "HardSymmetric" -> D2ExportExpression["hard_symmetric", candidateHard],
  "InheritedSymmetric" -> D2ExportExpression["inherited_symmetric", candidateInherited],
  "RunningAsSupplied" -> D2ExportExpression["running_as_supplied", candidateRunning],
  "TotalProvisional" -> D2ExportExpression["total_provisional", candidateTotal]
|>;

report = <|
  "Status" -> "C5FlavourCandidateExported",
  "NoNewMatchPerformed" -> True,
  "HardHolomorphicTermCount" -> Length[hardHol],
  "InheritedHolomorphicTermCount" -> Length[inhHol],
  "EachHolomorphicTermHasTwoDifferentLeptonFlavourLabels" -> True,
  "OriginalFlavourLabels" -> allowed,
  "MappedExternalFlavourLabels" -> {"T3D2ExternalFlavorP", "T3D2ExternalFlavorQ"},
  "ProjectConvention" -> WeinbergOperatorConvention[],
  "HardSymmetricByConstruction" -> True,
  "InheritedSymmetricByConstruction" -> True,
  "RunningLeftUnsymmetrized" -> True,
  "RunningSymmetricUnderExchangeByExpand" -> runSymmetric,
  "HardCoefficientHasNativeHbar" -> True,
  "InheritedCoefficientHasNativeHbar" -> True,
  "RunningInsertedWithHbar" -> True,
  "TotalCandidateLinearInHbar" -> True,
  "SourceFiles" -> required,
  "OutputRoots" -> roots,
  "FixedOrderC5CandidateAssembled" -> True,
  "OperatorNormalizationIndependentlyCertified" -> False,
  "InheritedSMCanonicalizationCertified" -> False,
  "PhysicalScalarMassRelationsImposed" -> False,
  "FixedOrderSMEFTCertified" -> False
|>;
Export[reportPath, report, "RawJSON"];
Print[logPrefix, report["Status"], "; running symmetric: ", runSymmetric,
  "; candidate total: ", roots["TotalProvisional"], ".wxf"];
Exit[0];
