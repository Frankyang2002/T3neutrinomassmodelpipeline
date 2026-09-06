(* Low-dimensional compatibility conventions for the generalized T3 basis.

   The generalized Lagrangian uses invariant tensors returned by GroupMagic /
   InvariantTensors. Those tensors are valid SU(2) invariants, but their
   normalization and phase need not equal Matchete's historical low-dimensional
   eps/gen/fStruct convention.

   For the five A-E regression benchmarks we measured

       R = C5_general / C5_legacy

   when the same numerical symbols y1, y2 and lambdaT3 are used.

   We choose the compatibility convention that leaves y1 and y2 unchanged and
   absorbs the complete basis conversion into lambdaT3:

       lambdaT3_general = lambdaT3_legacy / R.

   This map is ONLY a bridge to the historical/literature A-E convention.
   Arbitrary representations remain in the canonical generalized
   InvariantTensors convention and require no legacy conversion.
*)

ClearAll[
  T3GeneratedOverLegacyC5Ratio,
  T3LambdaT3GeneralFromLegacyFactor,
  T3LambdaT3LegacyFromGeneralFactor,
  T3LegacyConventionKnownQ
];

T3GeneratedOverLegacyC5Ratio[{1, 3, 2}] := -4/Sqrt[3];
T3GeneratedOverLegacyC5Ratio[{2, 2, 1}] := -2/Sqrt[3];
T3GeneratedOverLegacyC5Ratio[{2, 2, 3}] := 8/3;
T3GeneratedOverLegacyC5Ratio[{3, 1, 2}] := -4/Sqrt[3];
T3GeneratedOverLegacyC5Ratio[{3, 3, 2}] := -4 I Sqrt[2/3];

T3GeneratedOverLegacyC5Ratio[_List] := Missing["NoLegacyConvention"];

T3LegacyConventionKnownQ[dims_List] :=
  !MissingQ[T3GeneratedOverLegacyC5Ratio[dims]];

T3LambdaT3GeneralFromLegacyFactor[dims_List] := Module[{ratio},
  ratio = T3GeneratedOverLegacyC5Ratio[dims];
  If[MissingQ[ratio], ratio, FullSimplify[1/ratio]]
];

T3LambdaT3LegacyFromGeneralFactor[dims_List] :=
  T3GeneratedOverLegacyC5Ratio[dims];
