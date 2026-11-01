(* Analyze a Matchete InputForm C5 file without starting Matchete.
   Usage: wolframscript -file Lagrangian/AnalyzeGeneralC5.wl INPUT OUTPUT_DIR
   Output kernels retain hbar, phases, and physical mass symbols.
   No physical normalization or numerical matching is assumed. *)
ClearAll["Global`*"];
If[Length[Rest[$ScriptCommandLine]] != 2,
  Print["Usage: INPUT_C5_TXT OUTPUT_DIR"]; Exit[2]];
T3SourcePath = ExpandFileName[Rest[$ScriptCommandLine][[1]]];
T3OutputDir = ExpandFileName[Rest[$ScriptCommandLine][[2]]];
If[!FileExistsQ[T3SourcePath], Print["Input not found: ", T3SourcePath]; Exit[3]];
CreateDirectory[T3OutputDir, CreateIntermediateDirectories -> True];
T3Source = StringTrim[Import[T3SourcePath, "Text"]];
(* The export is a single InputForm expression, not Mathematica source code from users. *)
T3Expr = Quiet[Check[ToExpression[T3Source, InputForm], $Failed]];
If[T3Expr === $Failed,
  Print["Could not parse exported InputForm coefficient"]; Exit[4]];
T3Names = Sort@DeleteDuplicates[ToString[#, InputForm] & /@
   Cases[T3Expr, HoldPattern[Coupling[n_Symbol, ___]] :> n, Infinity]];
T3Known = {"y1", "y2", "y2Conj", "lambdaT3", "lambdaT3ExtraQ020101I1",
  "MF", "MS1", "T32PlusMass", "T32MinusMass", "T3NeutralPhase"};
Print["Coupling names: ", T3Names];
(* Preserve the original contraction structure in the separate raw export.
   The replacements below only factor the explicit Yukawa monomials. *)
T3FlavourExpression = T3Expr /. {
  HoldPattern[Bar[Coupling[y1, ___]]] :> T3Y1Bar,
  HoldPattern[Bar[Coupling[y2, ___]]] :> T3Y2Bar,
  HoldPattern[Bar[Coupling[y2Conj, ___]]] :> T3Y2ConjBar,
  HoldPattern[Coupling[y1, ___]] :> T3Y1,
  HoldPattern[Coupling[y2, ___]] :> T3Y2,
  HoldPattern[Coupling[y2Conj, ___]] :> T3Y2Conj
};
T3Products = {T3Y1Bar*T3Y2Bar, T3Y1Bar*T3Y2ConjBar,
              T3Y1*T3Y2, T3Y1*T3Y2Conj};
T3Labels = {"y1bar_y2bar", "y1bar_y2conjbar", "y1_y2", "y1_y2conj"};
T3Kernels = Table[
  Coefficient[Expand[T3FlavourExpression], T3Products[[j]]],
  {j, Length[T3Products]}];
T3Reconstruction = Total[MapThread[Times, {T3Products, T3Kernels}]];
T3Remainder = Expand[T3FlavourExpression - T3Reconstruction];
If[T3Remainder =!= 0,
  Export[FileNameJoin[{T3OutputDir, "unresolved_remainder.txt"}],
    ToString[T3Remainder, InputForm], "Text"];
  Print["ERROR: unexpected flavour structure; unresolved remainder exported"];
  Exit[5]];
T3KernelStatus = Table[
  Module[{T3Raw = T3Kernels[[j]], T3Simplified},
    T3Simplified = TimeConstrained[
      Factor[Cancel[Together[T3Raw]]], 60, $Aborted];
    If[T3Simplified === $Aborted, T3Simplified = T3Raw];
    Export[FileNameJoin[{T3OutputDir, T3Labels[[j]] <> "_kernel.txt"}],
      ToString[T3Simplified, InputForm], "Text"];
    <|"channel" -> T3Labels[[j]], "nonzero" -> (T3Raw =!= 0),
      "input_characters" -> StringLength[ToString[T3Raw, InputForm]],
      "kernel_characters" -> StringLength[ToString[T3Simplified, InputForm]],
      "simplification_completed" -> (T3Simplified =!= T3Raw)|>
  ], {j, Length[T3Products]}];
Export[FileNameJoin[{T3OutputDir, "channel_summary.json"}],
  <|"source" -> T3SourcePath,
    "coupling_symbols" -> T3Names,
    "channels" -> T3KernelStatus,
    "reconstruction_exact" -> True,
    "mass_basis_relation" -> "Mplus^2 = MS2^2 + 2 Abs[kappa]; Mminus^2 = MS2^2 - 2 Abs[kappa] (convention to verify)",
    "phase_relation" -> "phase^2 kappa = Abs[kappa] (convention to verify)",
    "normalization_certified" -> False,
    "numeric_C5_available" -> False|>, "RawJSON"];
Print["C5 decomposed; outputs in: ", T3OutputDir];
Exit[0];
