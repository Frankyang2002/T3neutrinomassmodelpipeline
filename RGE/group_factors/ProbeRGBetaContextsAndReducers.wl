(* Probe RGBeta contexts and internal tensor/group reducers.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaContextsAndReducers.wl

Unlike the previous probe, this explicitly loads RGBeta itself before loading
the T3 model file, then searches every RGBeta-related context rather than
assuming RGBeta`PackageScope`.
*)

ClearAll["Global`*"];

Print["=== RGBeta context/reducer probe ==="];

load = Quiet @ Check[Needs["RGBeta`"]; True, False];
Print["Needs[\"RGBeta`\"] success = ", load];

If[!TrueQ[load],
    Print["ERROR: RGBeta package could not be loaded with Needs[\"RGBeta`\"]."];
    Print["$Path = ", InputForm[$Path]];
    Exit[1];
];

modelFile = FileNameJoin[{
    Directory[],
    "RGE", "running", "wolfram", "T3RGBetaModel.wl"
}];

If[FileExistsQ[modelFile],
    Get[modelFile],
    Print["WARNING: T3 model file not found: ", modelFile]
];

Print["$Context = ", $Context];
Print["$ContextPath = ", InputForm[$ContextPath]];

contexts = Sort @ Select[Contexts[], StringStartsQ[#, "RGBeta`"] &];

Print["\n--- RGBeta-related contexts ---"];
Scan[Print, contexts];

patterns = {
    "Tensor", "Group", "Invariant", "Simpl", "Contract", "Reduce",
    "Index", "Quartic", "Upsilon", "Trace", "Rep", "Expand", "Evaluate"
};

allNames = DeleteDuplicates @ Flatten[
    Names[# <> "*"] & /@ contexts
];

selected = Select[
    allNames,
    Function[name,
        AnyTrue[patterns, StringContainsQ[name, #, IgnoreCase -> True] &]
    ]
];

Print["\n--- Candidate symbols across all RGBeta contexts ---"];
Scan[Print, Sort[selected]];

Print["\n--- Candidates whose definitions mention delS2/tGen/eps/del ---"];
Do[
    sym = Symbol[name];
    defs = Join[OwnValues[sym], DownValues[sym], UpValues[sym], SubValues[sym]];
    If[
        Length @ Cases[defs, _delS2 | _tGen | _eps | _del, Infinity] > 0,
        Print["\n>>> ", name];
        Print[InputForm[defs]];
    ],
    {name, selected}
];

Print["\n--- Short definitions for likely reducer names ---"];
likely = Select[
    selected,
    Function[name,
        AnyTrue[
            {"Simpl", "Contract", "Reduce", "Tensor", "Group"},
            StringContainsQ[name, #, IgnoreCase -> True] &
        ]
    ]
];

Do[
    sym = Symbol[name];
    dvs = DownValues[sym];
    If[Length[dvs] > 0,
        Print["\n>>> ", name];
        Print["Attributes: ", InputForm[Attributes[sym]]];
        Print["DownValues: ", InputForm[Take[dvs, UpTo[8]]]];
    ],
    {name, Sort[likely]}
];

Print["\nDONE"];
Exit[0];
