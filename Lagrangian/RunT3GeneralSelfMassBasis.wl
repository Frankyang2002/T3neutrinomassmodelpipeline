(*
  One-scale, full-general Matchete C5 for the six T3 self-bilinear models.

  Usage: wolframscript -file Lagrangian/RunT3GeneralSelfMassBasis.wl OUTPUT_DIR CLASS ALPHA
  Supported: A,D,E at alpha=-2 or 0; B(-1)/C(-1) are NOT handled here.

  The single real-neutral complex BSM scalar S = (R + I Phi)/Sqrt[2]
  is diagonalised BEFORE calling Matchete. For the source Lagrangian sign
  +mu S^dag S^dag + h.c. with real positive mu, R has mass squared
  mS^2-2 mu and Phi has mS^2+2 mu. The independent mass parameters
  T3{1,2}MinusMass / T3{1,2}PlusMass are mapped to the EFT1 running
  mS^2 and mu in GeneralT3MultiModelNeutrinos.mass_params.

  All renormalisable general Yukawas and quartics are rotated; only the
  original quadratic insertion is absorbed into the split free masses.
  A general complex mu phase is NOT handled by this real-phase pilot.
  No claim of complete sequential Wilson matching or basis certification.
*)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
projectRoot = DirectoryName[scriptDir];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3, Print["Usage: OUTPUT_DIR A|D|E m2|p0"]; Exit[2]];
out = ExpandFileName[args[[1]]];
family = ToUpperCase[args[[2]]];
alpha = Switch[args[[3]], "m2" | "-2", -2, "p0" | "0" | "+0", 0, _, 999];
If[!MemberQ[{"A", "D", "E"}, family] || !MemberQ[{-2, 0}, alpha],
  Print["ERROR: not one of the six self-bilinear T3 models"]; Exit[2]
];
If[Environment["T3_EXPERIMENTAL_GENERAL_FERMIONS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUARTICS"] =!= "1" ||
   Environment["T3_EXPERIMENTAL_GENERAL_QUADRATICS"] =!= "1",
  Print["ERROR: all three full-general flags required"]; Exit[2]
];
If[!DirectoryQ[out], CreateDirectory[out, CreateIntermediateDirectories -> True]];
printStatus[code_, msg_] := (Print["ERROR: ", msg]; Exit[code]);
Print["Loading Matchete for general self-bilinear mass-basis C5..."];
If[!TrueQ[CheckAbort[UsingFrontEnd[Needs["Matchete`"]; True], False]],
   printStatus[3, "Matchete not loaded"]];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];

model = T3ModelFromClass[family, alpha];
If[model === $Failed, printStatus[4, "model catalogue rejected point"]];
If[SetT3HeavyFields[{"F", "S1", "S2"}] === $Failed,
   printStatus[4, "heavy fields could not be assigned"]];
build = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[build] || Lookup[build,"Status",""] =!= "Success",
   printStatus[5, "full-general original-basis Lagrangian failed validation"]];
If[Lookup[build,"AdditionalQuadraticsIncluded",0] =!= 1 ||
   Lookup[build,"AdditionalYukawasIncluded",0] =!= 1 ||
   Lookup[build,"AdditionalQuarticsIncluded",-1] =!=
     Lookup[build,"AdditionalQuarticsExpected",-2],
   printStatus[5, "general interaction inventory incomplete"]];

(* Select the quadratic operator by its registered coupling name, NEVER
   by the names of its scalar fields (which would also remove quartics). *)
expanded = Expand[build["LInt"]];
allTerms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
hasExtraMassQ[term_] := !FreeQ[term,
  HoldPattern[Coupling[s_Symbol, ___]] /;
    StringStartsQ[SymbolName[s], "lambdaT3ExtraM"]];
extraMassTerms = Select[allTerms, hasExtraMassQ];
If[Length[extraMassTerms] < 2,
  printStatus[6, "could not isolate both holomorphic mass insertion and h.c."]];
interactions = Expand[Total[Select[allTerms, !hasExtraMassQ[#] &]]];
If[hasExtraMassQ[interactions],
   printStatus[6, "self-bilinear interaction still present after mass-basis extraction"]];

species = If[alpha === -2, 2, 1];
original = If[species === 2, NewScalar2, NewScalar1];
plus = If[species === 2, T3Scalar2Plus, T3Scalar1Plus];
minus = If[species === 2, T3Scalar2Minus, T3Scalar1Minus];
plusMass = If[species === 2, T32PlusMass, T31PlusMass];
minusMass = If[species === 2, T32MinusMass, T31MinusMass];
repDim = model[If[species === 2, "Scalar2", "Scalar1"], "SU2"];
If[!MemberQ[{1, 3}, repDim] ||
  model[If[species === 2, "Scalar2", "Scalar1"], "Y"] =!= 0,
  printStatus[6, "self-bilinear scalar is not a neutral real SU(2) irrep"]];
indexType = BSMIndexType[repDim];
indices = If[indexType === None, {}, {indexType}];

(* The original field must not remain in Matchete's registry: it has a
   non-diagonal holomorphic mass. Replace by TWO distinct real heavy fields.
   The mass symbols encode the nonzero bilinear; they are not zeroed. *)
If[species === 2, RemoveField[NewScalar2], RemoveField[NewScalar1]];
With[{field = plus, mass = plusMass, ix = indices},
  DefineField[field, Scalar, SelfConjugate -> True,
              Indices -> ix, Mass -> {Heavy, mass}]
];
With[{field = minus, mass = minusMass, ix = indices},
  DefineField[field, Scalar, SelfConjugate -> True,
              Indices -> ix, Mass -> {Heavy, mass}]
];

(* For real positive mu, the + mass belongs to Im[S] and the - mass to
   Re[S]. Preserve the original SU(2) indices and ALL flavour couplings.

   IMPORTANT: Matchete has ALREADY evaluated NewScalar2[i] to
     Field[NewScalar2, Scalar, {Index[i, T3BSMd3]}, {}]
   by the time build["LInt"] is returned. A rule for NewScalar2[i]
   alone cannot rotate this expression. The Bar[Field[...]] rules must
   precede the Field[...] rules, otherwise a conjugated field might be
   substituted with the wrong sign for the imaginary component.

   Keep the bare-field rules as a fallback for nonexpanded source forms.
   No self-bilinear coupling is set to zero here: its effect is encoded in
   T31/32PlusMass and T31/32MinusMass defined above. *)
rotationRules = {
    (* Matchete-evaluated conjugate fields: observed in live A(-2) log. *)
    HoldPattern[Bar[Field[NewScalar1, Scalar, ix_, opts_]]] /;
       species === 1 :>
       (Bar[Field[T3Scalar1Minus, Scalar, ix, opts]] -
        I Bar[Field[T3Scalar1Plus, Scalar, ix, opts]])/Sqrt[2],
    HoldPattern[Bar[Field[NewScalar2, Scalar, ix_, opts_]]] /;
       species === 2 :>
       (Bar[Field[T3Scalar2Minus, Scalar, ix, opts]] -
        I Bar[Field[T3Scalar2Plus, Scalar, ix, opts]])/Sqrt[2],
    (* Matchete-evaluated unbarred fields. *)
    HoldPattern[Field[NewScalar1, Scalar, ix_, opts_]] /;
       species === 1 :>
       (Field[T3Scalar1Minus, Scalar, ix, opts] +
        I Field[T3Scalar1Plus, Scalar, ix, opts])/Sqrt[2],
    HoldPattern[Field[NewScalar2, Scalar, ix_, opts_]] /;
       species === 2 :>
       (Field[T3Scalar2Minus, Scalar, ix, opts] +
        I Field[T3Scalar2Plus, Scalar, ix, opts])/Sqrt[2],
    (* Unevaluated field forms, retained for any source expressions. *)
    HoldPattern[Bar[NewScalar1[ss___]]] /; species === 1 :>
       (Bar[T3Scalar1Minus[ss]] - I Bar[T3Scalar1Plus[ss]])/Sqrt[2],
    HoldPattern[Bar[NewScalar2[ss___]]] /; species === 2 :>
       (Bar[T3Scalar2Minus[ss]] - I Bar[T3Scalar2Plus[ss]])/Sqrt[2],
    HoldPattern[Bar[NewScalar1][ss___]] /; species === 1 :>
       (Bar[T3Scalar1Minus][ss] - I Bar[T3Scalar1Plus][ss])/Sqrt[2],
    HoldPattern[Bar[NewScalar2][ss___]] /; species === 2 :>
       (Bar[T3Scalar2Minus][ss] - I Bar[T3Scalar2Plus][ss])/Sqrt[2],
    HoldPattern[NewScalar1[ss___]] /; species === 1 :>
       (T3Scalar1Minus[ss] + I T3Scalar1Plus[ss])/Sqrt[2],
    HoldPattern[NewScalar2[ss___]] /; species === 2 :>
       (T3Scalar2Minus[ss] + I T3Scalar2Plus[ss])/Sqrt[2]
};
(* Clear a diagnostic from a prior attempt; a successful rotation must not
   leave a stale failure report in the same per-model output folder. *)
rotationDebugPath = FileNameJoin[{out, "self_mass_basis_rotation_debug.json"}];
If[FileExistsQ[rotationDebugPath], DeleteFile[rotationDebugPath]];
rotatedRaw = Expand[interactions /. rotationRules];
rotated = RelabelIndices[rotatedRaw];

(* Keep the fail-closed mass-basis check, but make it diagnostically useful.
   A mass term or original field in the interaction cannot be discarded. *)
listTerms[expr_] := Module[{e = Expand[expr]},
  If[Head[e] === Plus, List @@ e, {e}]
];
shortString[x_] := Module[{str = ToString[x, InputForm]},
  StringTake[str, Min[1400, StringLength[str]]]
];
(* Check for the exact Matchete Field[old, Scalar, ...] wrapper as well as
   old bare symbols: the live failure was 23 such evaluated field terms. *)
wrappedBefore = Cases[interactions,
  HoldPattern[Field[s_Symbol, Scalar, _, _]] /; s === original,
  Infinity];
wrappedAfter = Cases[rotatedRaw,
  HoldPattern[Field[s_Symbol, Scalar, _, _]] /; s === original,
  Infinity];
Print["ROTATION FIELD WRAPPERS: before=", Length[wrappedBefore],
      "; after=", Length[wrappedAfter]];
If[wrappedBefore === {},
   printStatus[6, "no evaluated original Field[...] objects found to rotate"]];
If[FreeQ[rotatedRaw, plus] || FreeQ[rotatedRaw, minus],
   printStatus[6, "the rotation did not retain BOTH scalar mass eigenfields"]];
oldRaw = Select[listTerms[rotatedRaw], !FreeQ[#, original] &];
oldFinal = Select[listTerms[rotated], !FreeQ[#, original] &];
massRaw = Select[listTerms[rotatedRaw], hasExtraMassQ];
massFinal = Select[listTerms[rotated], hasExtraMassQ];
Print["ROTATION CHECK: old field before RelabelIndices=", Length[oldRaw],
  "; after=", Length[oldFinal],
  "; extra mass before=", Length[massRaw],
  "; after=", Length[massFinal]];
If[oldFinal =!= {} || massFinal =!= {},
  rotationDiagnostic = <|
    "model" -> family <> "_" <> If[alpha === -2, "-2", "+0"],
    "status" -> "UnrotatedFieldOrMassDetected",
    "original_field" -> SymbolName[original],
    "old_field_terms_before_relabel_count" -> Length[oldRaw],
    "old_field_terms_after_relabel_count" -> Length[oldFinal],
    "mass_terms_before_relabel_count" -> Length[massRaw],
    "mass_terms_after_relabel_count" -> Length[massFinal],
    "old_field_terms_before_relabel" -> shortString /@ Take[oldRaw, UpTo[6]],
    "old_field_terms_after_relabel" -> shortString /@ Take[oldFinal, UpTo[6]],
    "mass_terms_before_relabel" -> shortString /@ Take[massRaw, UpTo[6]],
    "mass_terms_after_relabel" -> shortString /@ Take[massFinal, UpTo[6]]
  |>;
  Export[FileNameJoin[{out, "self_mass_basis_rotation_debug.json"}],
    rotationDiagnostic, "RawJSON"];
  Print["ROTATION DEBUG saved to ",
    FileNameJoin[{out, "self_mass_basis_rotation_debug.json"}]];
  If[oldFinal =!= {},
    Print["First unrotated old-field term: ", shortString[First[oldFinal]]]];
  If[massFinal =!= {},
    Print["First remaining holomorphic mass term: ",
      shortString[First[massFinal]]]];
  printStatus[6, "unrotated field or holomorphic mass survives; see diagnostic JSON"]
];
remaining = If[species === 2, NewScalar1, NewScalar2];
smLag = Expand[build["LUV"] - build["LBSM"]];
freeEigen = With[{rest = remaining, a = minus, b = plus},
  FreeLag[NewFermion, rest, a, b]
];
newUV = Expand[smLag + freeEigen + rotated] // RelabelIndices;
Print["Self-bilinear: species=S",species,
      ", rep dimension=",repDim,
      "; bilinear absorbed in masses ",minusMass," / ",plusMass];
Print["Retained ordinary + additional Yukawas and ALL original general quartics."];

(* Fail closed if field redefinitions have not produced a valid canonical
   theory. This is a required live check, not a claim it has already passed. *)
valid = CheckAbort[UsingFrontEnd[CheckLagrangian[newUV]], $Aborted];
Print["Rotated CheckLagrangian result: ",InputForm[valid]];
If[valid =!= True,
  Export[FileNameJoin[{out,"self_mass_basis_debug.txt"}],
         ToString[Short[newUV,25],InputForm],"Text"];
  printStatus[7, "rotated Lagrangian not validated by Matchete"]
];

(* One-scale direct matching in the mass eigenbasis; do not silently invoke
   the original non-diagonal RunModel.wl path. *)
matching = CheckAbort[UsingFrontEnd[RunT3Matching[newUV,5,1]], $Aborted];
If[!AssociationQ[matching] || Lookup[matching,"Status",""] =!= "Success",
  Print["Matching status: ",InputForm[If[AssociationQ[matching],Lookup[matching,"Status"],matching]]];
  printStatus[8, "mass-eigenfield general one-loop matching failed"]
];
matched = matching["MatchedEFT"];
weinberg = ExtractWeinbergCoefficient[matched];
coefficient = Lookup[weinberg,"Coefficient",Missing["NoWeinbergCoefficient"]];
If[MissingQ[coefficient] || coefficient === 0 || coefficient === $Failed,
  printStatus[9, "no Weinberg coefficient from the general mass-basis match"]
];
sourceText = ToString[coefficient,InputForm];
If[!StringContainsQ[sourceText,"T31PlusMass" | "T32PlusMass"] ||
   !StringContainsQ[sourceText,"T31MinusMass" | "T32MinusMass"],
   printStatus[9, "C5 lost at least one self-bilinear mass eigenstate"]];
If[!StringContainsQ[sourceText,"y1"] || !StringContainsQ[sourceText,"y2"] ||
   !StringContainsQ[sourceText,If[alpha === -2,"y2Conj","y1Conj"]],
   printStatus[9, "C5 lacks a required full-general Yukawa channel"]];
Export[FileNameJoin[{out,"c5_coefficient.txt"}],sourceText,"Text"];
Export[FileNameJoin[{out,"self_mass_basis_metadata.json"}],<|
    "status"->"CandidateRealSelfBilinearGeneralOneScaleMatched",
    "model"->family<>"_"<>If[alpha===-2,"-2","+0"],
    "split_species"->species,
    "mass_symbol_plus"->ToString[plusMass,InputForm],
    "mass_symbol_minus"->ToString[minusMass,InputForm],
    "original_general_quartics"->build["AdditionalQuarticsIncluded"],
    "additional_yukawas"->build["AdditionalYukawasIncluded"],
    "self_bilinear_nonzero_by_mass_splitting"->True,
    "complex_bilinear_phase_supported"->False,
    "sequential_eft_wilson_matching_performed"->False,
    "physical_coupling_map_certified"->False
    |>,"RawJSON"];
Print["GENERAL SELF-BILINEAR MASS-BASIS C5: SUCCESS"];
Exit[0];
