(* Focused probe of RGBeta's internal group/tensor reduction machinery.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaReducerDefinitions.wl
*)

ClearAll["Global`*"];

Needs["RGBeta`"];

modelFile = FileNameJoin[{
    Directory[],
    "RGE", "running", "wolfram", "T3RGBetaModel.wl"
}];
If[FileExistsQ[modelFile], Get[modelFile]];

Print["=== Focused RGBeta reducer definitions ==="];

names = {
    "RGBeta`PackageScope`RefineGroupStructures",
    "RGBeta`BetaSimplify",
    "RGBeta`PackageScope`BetaSimplify",
    "RGBeta`GroupsAndIndices`PackagePrivate`twoIndexRepDelta",
    "RGBeta`GroupsAndIndices`PackagePrivate`PerformGeneratorTraces",
    "RGBeta`GroupsAndIndices`PackagePrivate`PerformAdjAlg",
    "RGBeta`GroupsAndIndices`PackagePrivate`CanonizeAdjTrace",
    "RGBeta`GroupsAndIndices`PackagePrivate`AdjContraction",
    "RGBeta`TensorCalculations`PackagePrivate`CompleteReplace",
    "RGBeta`PackageScope`$scalarContraction",
    "RGBeta`GroupsAndIndices`PackagePrivate`replace"
};

Do[
    Print["\n>>> ", name];
    If[
        NameQ[name],
        sym = Symbol[name];
        Print["Context: ", Context[sym]];
        Print["Attributes: ", InputForm[Attributes[sym]]];
        Print["OwnValues: ", InputForm[OwnValues[sym]]];
        Print["DownValues: ", InputForm[DownValues[sym]]];
        Print["UpValues: ", InputForm[UpValues[sym]]];
        Print["SubValues: ", InputForm[SubValues[sym]]],
        Print["NOT FOUND"]
    ],
    {name, names}
];

Print["\n=== Minimal symbolic reduction experiments ==="];

repRules = If[
    NameQ["RGBeta`GroupsAndIndices`PackagePrivate`replace"],
    Symbol["RGBeta`GroupsAndIndices`PackagePrivate`replace"],
    {}
];

expr1 = delS2[SU2L, a, i, j];
expr2 = tGen[SU2L @ S2, A, a, b];
expr3 = del[SU2L @ S2, a, b];

Print["delS2 /. replace = ", InputForm[expr1 /. repRules]];
Print["tGen[S2] /. replace = ", InputForm[expr2 /. repRules]];
Print["del[S2] /. replace = ", InputForm[expr3 /. repRules]];

If[NameQ["RGBeta`PackageScope`RefineGroupStructures"],
    refine = Symbol["RGBeta`PackageScope`RefineGroupStructures"];
    Print["RefineGroupStructures[delS2] = ",
        InputForm @ Quiet @ Check[refine[expr1], $Failed]
    ];
    Print["RefineGroupStructures[tGenS2] = ",
        InputForm @ Quiet @ Check[refine[expr2], $Failed]
    ];
];

If[NameQ["RGBeta`BetaSimplify"],
    bs = Symbol["RGBeta`BetaSimplify"];
    Print["BetaSimplify[delS2] = ",
        InputForm @ Quiet @ Check[bs[expr1], $Failed]
    ];
];

Print["\nDONE"];
Exit[0];
