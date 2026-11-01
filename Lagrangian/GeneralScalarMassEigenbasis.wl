(* Physical mass-basis utilities for a neutral complex scalar in a real SU(2) representation.
   A neutral scalar with holomorphic bilinear kappa S.S + h.c. splits into real
   modes. This file prepares exact symbolic rotations and quadratic masses.

   Does not modify Matchete field registry or claim matching success.
   Matchete must receive newly defined mass-eigenstate fields with rotated
   kinetic, gauge, Yukawa, quartic, and quadratic terms before Match[].
*)
ClearAll[
  T3NeutralQuadraticRealMatrix,
  T3NeutralQuadraticSpectrum,
  T3NeutralQuadraticRotation,
  T3NeutralQuadraticFieldRules,
  T3NeutralQuadraticExport
];

(* Potential convention:
   V = mSq Conjugate[S] S + kappa S^2 + Conjugate[kappa] Conjugate[S]^2.
   Here the holomorphic term has NO 1/2; compare carefully to PlusHc
   normalization in the model builder. This convention is different from
   V= ... + (B S^2 +h.c.)/2, where B=2 kappa. *)
T3NeutralQuadraticRealMatrix[mSq_, kappa_] := {
  {mSq + 2 Re[kappa], -2 Im[kappa]},
  {-2 Im[kappa], mSq - 2 Re[kappa]}
};

T3NeutralQuadraticSpectrum[mSq_, kappa_] := {
  mSq - 2 Sqrt[Re[kappa]^2 + Im[kappa]^2],
  mSq + 2 Sqrt[Re[kappa]^2 + Im[kappa]^2]
};

(* Symmetric real rotation U with columns corresponding to ascending
   eigenvalues. Suitable for exact symbolic complex kappa. For kappa==0,
   use the identity to avoid undefined Arg and preserve the original basis. *)
T3NeutralQuadraticRotation[kappa_] := Module[{phase, theta},
  If[TrueQ[PossibleZeroQ[kappa]], Return[IdentityMatrix[2]]];
  phase = Arg[kappa];
  theta = -phase/2;
  {
    {-Sin[theta], Cos[theta]},
    {Cos[theta], Sin[theta]}
  }
];

(* Return symbolic substitutions that express the original scalar and its
   complex conjugate in terms of two real mass eigenfields rMinus, rPlus.
   The gauge transformation and field definitions must be handled separately. *)
T3NeutralQuadraticFieldRules[scalar_, scalarConjugate_, rMinus_, rPlus_,
    kappa_] := Module[{rot, fields, z},
  rot = T3NeutralQuadraticRotation[kappa];
  fields = rot . {rMinus, rPlus};
  z = (fields[[1]] + I fields[[2]])/Sqrt[2];
  {scalar -> z, scalarConjugate -> ComplexExpand[Conjugate[z]]}
];

T3NeutralQuadraticExport[filename_String, mSq_, kappa_] := Module[
  {matrix, masses, rotation, payload},
  matrix = T3NeutralQuadraticRealMatrix[mSq, kappa];
  masses = T3NeutralQuadraticSpectrum[mSq, kappa];
  rotation = T3NeutralQuadraticRotation[kappa];
  payload = <|
    "status" -> "SymbolicMassBasisOnly_NoMatcheteMatching",
    "convention" -> "V=mSq Sdag S+kappa S^2+conjugate(kappa) Sdag^2",
    "quadratic_matrix_inputform" -> ToString[matrix, InputForm],
    "mass_squared_inputform" -> (ToString[#, InputForm] & /@ masses),
    "rotation_inputform" -> ToString[rotation, InputForm],
    "matchete_field_redefinition_performed" -> False
  |>;
  Export[filename, payload, "RawJSON"];
  payload
];
