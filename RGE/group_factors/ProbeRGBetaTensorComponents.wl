(* Probe exact RGBeta SU(2) tensor component conventions.

Run from the repository root with:
  wolframscript -file RGE/group_factors/ProbeRGBetaTensorComponents.wl

The script intentionally uses the RGBeta tensor heads delS2 and tGen directly,
after loading the project's T3RGBetaModel.wl.  It prints explicit component
tables if those heads evaluate numerically for integer indices.  These tables
are the input needed for a component-level Ward-identity test that uses the
same conventions as RGBeta itself.
*)

ClearAll["Global`*"];

modelFile = FileNameJoin[{
    Directory[],
    "RGE", "running", "wolfram", "T3RGBetaModel.wl"
}];

If[!FileExistsQ[modelFile],
    Print["ERROR: cannot find ", modelFile];
    Exit[1];
];

Get[modelFile];

Print["=== RGBeta tensor component probe ==="];
Print["Model file: ", modelFile];

Print["\n--- delS2[SU2L,a,i,j], a=1..3, i,j=1..2 ---"];
delTable = Table[
    Quiet @ Check[
        FullSimplify[delS2[SU2L, a, i, j]],
        $Failed
    ],
    {a, 1, 3}, {i, 1, 2}, {j, 1, 2}
];
Print[InputForm[delTable]];

Print["\n--- tGen[SU2L@fund,A,i,j] ---"];
fundGen = Table[
    Quiet @ Check[
        FullSimplify[tGen[SU2L @ fund, A, i, j]],
        $Failed
    ],
    {A, 1, 3}, {i, 1, 2}, {j, 1, 2}
];
Print[InputForm[fundGen]];

Print["\n--- tGen[SU2L@S2,A,i,j] ---"];
tripGen = Table[
    Quiet @ Check[
        FullSimplify[tGen[SU2L @ S2, A, i, j]],
        $Failed
    ],
    {A, 1, 3}, {i, 1, 3}, {j, 1, 3}
];
Print[InputForm[tripGen]];

Print["\n--- T3RGBetaS2TripleInvariant[a,b,c] ---"];
triple = Table[
    Quiet @ Check[
        FullSimplify[T3RGBetaS2TripleInvariant[a, b, c]],
        $Failed
    ],
    {a, 1, 3}, {b, 1, 3}, {c, 1, 3}
];
Print[InputForm[triple]];

Print["\n--- Basic evaluation diagnostics ---"];
Print[
    "delS2 residual symbolic heads = ",
    Length @ Cases[delTable, _delS2 | _eps | _tGen, Infinity]
];
Print[
    "fund tGen residual symbolic heads = ",
    Length @ Cases[fundGen, _tGen, Infinity]
];
Print[
    "triplet tGen residual symbolic heads = ",
    Length @ Cases[tripGen, _tGen, Infinity]
];
Print[
    "triple residual symbolic heads = ",
    Length @ Cases[triple, _delS2 | _eps | _tGen, Infinity]
];

If[
    Or[
        Length @ Cases[delTable, _delS2 | _eps | _tGen, Infinity] > 0,
        Length @ Cases[fundGen, _tGen, Infinity] > 0,
        Length @ Cases[tripGen, _tGen, Infinity] > 0,
        Length @ Cases[triple, _delS2 | _eps | _tGen, Infinity] > 0
    ],
    Print[
        "\nNOTE: one or more RGBeta tensor heads remained symbolic. ",
        "If so, the next Ward checker must use RGBeta's internal tensor ",
        "reduction functions rather than direct integer component evaluation."
    ],
    Print[
        "\nSUCCESS: all requested tensor heads evaluated componentwise. ",
        "These exact tables can be used directly in the Ward-identity checker."
    ]
];

Exit[0];
