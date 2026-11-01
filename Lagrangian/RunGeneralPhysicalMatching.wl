(* Full-general T3 matching attempt in a diagonal neutral-scalar basis.
   Usage: wolframscript -file Lagrangian/RunGeneralPhysicalMatching.wl
             OUTPUT CLASS alpha
   Uses the project's actual general Lagrangian and existing Matchete
   matching/Weinberg extraction. Does not claim success unless C5 exists.
   A neutral self-bilinear is represented by two real multiplets.
   Mixed S1/S2 bilinears are explicitly unsupported in this implementation.
*)
ClearAll["Global`*"];
scriptDir = DirectoryName @ ExpandFileName[$InputFileName];
args = Rest[$ScriptCommandLine];
If[Length[args] != 3,
  Print["Usage: OUTPUT CLASS alpha"]; Exit[2]];
output = ExpandFileName[args[[1]]];
class = ToUpperCase[args[[2]]];
parseAlpha[x_String] := Which[
  StringMatchQ[x, "m" ~~ DigitCharacter ..], -ToExpression[StringDrop[x,1]],
  StringMatchQ[x, "p" ~~ DigitCharacter ..], ToExpression[StringDrop[x,1]],
  StringMatchQ[x, "-" ~~ DigitCharacter .. | DigitCharacter ..], ToExpression[x],
  True, $Failed
];
alpha = parseAlpha[args[[3]]];
If[alpha === $Failed, Print["Invalid alpha"]; Exit[2]];
CreateDirectory[output, CreateIntermediateDirectories -> True];
Needs["Matchete`"];
Get[FileNameJoin[{scriptDir, "PhysicsNotation.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "T3ModelCatalog.wl"}]];
Get[FileNameJoin[{scriptDir, "model", "LagrangianBuilder.wl"}]];
Get[FileNameJoin[{scriptDir, "RunMatching.wl"}]];
model = T3ModelFromClass[class, alpha];
If[model === $Failed, Print["Unsupported model"]; Exit[2]];
If[SetT3HeavyFields[{"F","S1","S2"}] === $Failed, Exit[2]];
Print["General physical-mass matching: ", class, "_", alpha];
build = UsingFrontEnd[BuildT3Lagrangian[model]];
If[!AssociationQ[build], Print["Lagrangian build failed"]; Exit[10]];
If[Lookup[build,"Status",""] === "RequiresMixedScalarMassDiagonalization",
  Print["Mixed-multiplet scalar diagonalization not implemented here"];
  Exit[18]
];
If[Lookup[build,"Status",""] =!= "Success",
  Print["Lagrangian build failed: ", build["Status"]]; Exit[10]];

(* Identify a holomorphic self-bilinear using the actual accepted
   interaction inventory. The coupling symbol is taken from its source
   structure and is not guessed from a generic mass matrix. *)
quadratics = Select[
  T3ExtendedExtraQuadraticStructures[model],
  Total[#] == 2 &
];
selfBilinears = Select[quadratics,
  (#[[3]]+#[[4]]==2 && #[[5]]+#[[6]]==0) ||
  (#[[5]]+#[[6]]==2 && #[[3]]+#[[4]]==0) &
];
If[Length[selfBilinears] > 1,
  Print["Multiple self bilinears require enlarged mass basis"]; Exit[18]];
If[Length[selfBilinears] == 0,
  Print["No self-bilinear: invoking existing matching without field change"];
  matching = RunSequentialT3Matching[
    build["LUV"], model, {{"F","S1","S2"}}, 5, 1, output, class],
  (* We choose field number from the two identical fields. *)
  counts = First[selfBilinears];
  which = If[counts[[3]]+counts[[4]] == 2, 1, 2];
  name = If[which == 1, NewScalar1, NewScalar2];
  mass = If[which == 1, MS1, MS2];
  scalarDimension = model[If[which==1,"Scalar1","Scalar2"],"SU2"];
  scalarHypercharge = model[If[which==1,"Scalar1","Scalar2"],"Y"];
  If[scalarHypercharge =!= 0 || EvenQ[scalarDimension],
    Print["Only neutral real SU(2) representations can split into real multiplets"];
    Exit[18]];
  couplingName = Symbol["lambdaT3ExtraM"<>StringJoin[ToString /@ counts]<>"I1"];
  kappa = Coupling[couplingName, {}, 0];
  plusMass = Symbol["T3"<>ToString[which]<>"PlusMass"];
  minusMass = Symbol["T3"<>ToString[which]<>"MinusMass"];
  (* DefineCoupling expects a bare symbol, not an applied expression. *)
  If[Check[DefineCoupling[T3NeutralPhase, SelfConjugate -> False], $Failed] === $Failed,
    Print["Unable to register neutral-scalar phase coupling"]; Exit[19]];
  phase = T3NeutralPhase[];
  plusField = If[which==1,T3S1Plus,T3S2Plus];
  minusField = If[which==1,T3S1Minus,T3S2Minus];
  (* The phase convention is:
     original S = phase (Splus + I Sminus)/Sqrt[2],
     phase^2 kappa = |kappa|, with real masses squared
     Mplus^2 = MS2sq + 2|kappa|, Mminus^2 = MS2sq - 2|kappa|.
     Matching treats phase, Mplus, Mminus as independent couplings;
     the output JSON records their physical relation. *)
  (* Cache the original scalar field encoding before resetting its registry. *)
  originalScalarFields = DeleteDuplicates @ Cases[build["LUV"], _Field, Infinity];
  Print["Original field registry encodings: ", InputForm[originalScalarFields]];
  RemoveField[name];
  indexType = BSMIndexType[scalarDimension];
  indices = If[indexType === None, {}, {indexType}];
  DefineField[plusField, Scalar, SelfConjugate -> True,
    Indices -> indices, Mass -> {Heavy,plusMass}];
  DefineField[minusField, Scalar, SelfConjugate -> True,
    Indices -> indices, Mass -> {Heavy,minusMass}];
  (* Matchete's evaluated Lagrangian contains Field[NewScalar2,Scalar,...].
     Match directly on its structural representation. Do not use Cases on
     Bar[_Field] alone: Bar is not necessarily a literal head after evaluation.
     Apply explicit field rules to all scalar factors, including derivatives;
     the new real multiplets were registered above. *)
  lagInt = build["LInt"];
  oldFieldsInUV = DeleteDuplicates @ Cases[
    build["LUV"],
    HoldPattern[Field[NewScalar2, Scalar, __]] :> "S2", Infinity];
  If[which == 1,
    oldFieldsInUV = DeleteDuplicates @ Cases[
      build["LUV"],
      HoldPattern[Field[NewScalar1, Scalar, __]] :> "S1", Infinity]
  ];
  Print["Old scalar Field objects in UV: ", Length[oldFieldsInUV]];
  If[Length[oldFieldsInUV] == 0,
    Print["ERROR: no original Field objects in UV expression."]; Exit[19]
  ];
  (* Use a rule defined with a literal symbol in each branch. The pattern
     does not depend on Mathematica's evaluation of a symbolic variable. *)
  replacements = If[which == 1,
    {
      HoldPattern[Bar[Field[NewScalar1, Scalar, ix_, dx_]]] :>
        Bar[phase] (Field[T3S1Plus, Scalar, ix, dx] -
          I Field[T3S1Minus, Scalar, ix, dx])/Sqrt[2],
      HoldPattern[Field[NewScalar1, Scalar, ix_, dx_]] :>
        phase (Field[T3S1Plus, Scalar, ix, dx] +
          I Field[T3S1Minus, Scalar, ix, dx])/Sqrt[2]
    },
    {
      HoldPattern[Bar[Field[NewScalar2, Scalar, ix_, dx_]]] :>
        Bar[phase] (Field[T3S2Plus, Scalar, ix, dx] -
          I Field[T3S2Minus, Scalar, ix, dx])/Sqrt[2],
      HoldPattern[Field[NewScalar2, Scalar, ix_, dx_]] :>
        phase (Field[T3S2Plus, Scalar, ix, dx] +
          I Field[T3S2Minus, Scalar, ix, dx])/Sqrt[2]
    }
  ];
  (* Do not confuse explicit Mass -> {Heavy,...} labels with interactions. *)
  (* The old holomorphic mass insertion is now represented exactly
     by the split masses in FreeLag; do not count it twice. *)
  terms = If[Head[Expand[lagInt]] === Plus,
    List@@Expand[lagInt], {Expand[lagInt]}];
  kept = Select[terms, FreeQ[HoldComplete[#], couplingName]&];
  Print["Neutral self-bilinear terms removed from interaction basis: ",
        Length[terms]-Length[kept]];
  transformed = Expand[Total[kept] /. replacements];
  remaining = Cases[transformed,
    HoldPattern[Field[NewScalar2, Scalar, __]] :> "S2", Infinity];
  If[which == 1, remaining = Cases[transformed,
    HoldPattern[Field[NewScalar1, Scalar, __]] :> "S1", Infinity]];
  Print["Remaining original scalar factors: ", Length[remaining]];
  Print["Transformed interaction terms: ",
    Length[If[Head[transformed] === Plus, List @@ transformed, {transformed}]]];
  If[Length[remaining] != 0,
    Print["Unconverted original scalar field remains"]; Exit[19]];
  fields = If[which==1,
    {NewFermion,plusField,minusField,NewScalar2},
    {NewFermion,NewScalar1,plusField,minusField}];
  (* Build the new free theory using the registered mass eigenfields. *)
  physicalFree = Apply[FreeLag, fields] // RelabelIndices;
  Print["Checking rebuilt free plus SM theory..."];
  baseCheck = CheckAbort[CheckLagrangian[(build["LSM"] + physicalFree) // RelabelIndices], $Aborted];
  Print["  rebuilt free + SM check: ", InputForm[baseCheck]];
  If[baseCheck =!= True,
    Export[FileNameJoin[{output,"failed_physical_free.txt"}],
      ToString[physicalFree, InputForm], "Text"];
    Print["ERROR: mass-eigenfield registration/free kinetic sector invalid."];
    Exit[21]
  ];
  (* Explicit post-transformation index canonicalisation. *)
  transformed = RelabelIndices[Expand[transformed]];
  physicalLag = RelabelIndices[Expand[build["LSM"] + physicalFree + transformed]];
  (* Do not feed an invalid theory to Match: it produces a misleading Times
     return followed by GreensSimplify index errors. *)
  Print["Checking physical Lagrangian after field replacement..."];
  (* Keep Matchete messages visible. Quiet was suppressing the specific
     non-canonical / gauge-charge failure underlying a bare False. *)
  physicalCheck = CheckAbort[CheckLagrangian[physicalLag], $Aborted];
  If[physicalCheck =!= True,
    Export[FileNameJoin[{output,"failed_physical_lagrangian.txt"}],
      ToString[physicalLag,InputForm],"Text"];
    Print["ERROR: transformed Lagrangian is not valid: ", InputForm[physicalCheck]];
    Print["Saved transformed expression for debugging."];
    Exit[20]
  ];
  Export[FileNameJoin[{output,"scalar_mass_basis.json"}],
    <|"model"->class<>"_"<>ToString[alpha],
      "split_field"->ToString[name,InputForm],
      "holomorphic_coupling"->ToString[couplingName,InputForm],
      "positive_mass_squared"->ToString[HoldForm[mass^2+2 Abs[kappa]],InputForm],
      "negative_mass_squared"->ToString[HoldForm[mass^2-2 Abs[kappa]],InputForm],
      "phase_constraint"->"T3NeutralPhase^2 * kappa = Abs[kappa]",
      "status"->"SymbolicMassBasisPrepared_MatchingNotYetCertified"|>,
      "RawJSON"];
  Print["NOTE: triplet real-basis CG and phase normalization require validation."];
  Print["Running Matchete with split real scalar multiplets"];
  matching = RunT3Matching[physicalLag,5,1]
];
If[!AssociationQ[matching] || Lookup[matching,"Status",""] =!= "Success",
  Print["Matchete matching failed: ",Lookup[matching,"Status","InvalidResult"]];
  Exit[12]];
eft = Lookup[matching,"MatchedEFT",Lookup[matching,"LoopEFT",0]];
weinberg = ExtractWeinbergCoefficient[eft];
coefficient = Lookup[weinberg,"Coefficient",Missing["Unavailable"]];
If[MissingQ[coefficient],
  Print["Matching succeeded but physical C5 extraction is pending"];
  Exit[14]];
Export[FileNameJoin[{output,"c5_coefficient.txt"}],
  ToString[coefficient,InputForm],"Text"];
Print["Full-general symbolic C5 coefficient exported"];
Exit[0];
