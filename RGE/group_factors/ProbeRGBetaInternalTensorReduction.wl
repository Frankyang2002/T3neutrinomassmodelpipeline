(* Discover RGBeta internal tensor-reduction symbols relevant to delS2/tGen.

Run from project root:
  wolframscript -file RGE/group_factors/ProbeRGBetaInternalTensorReduction.wl

This script does not modify the model.  It loads the same T3RGBetaModel.wl
used by the project, then inspects RGBeta`PackageScope` symbols whose names
suggest tensor/group/invariant reduction.  It also prints DownValues for the
most relevant candidates when present.

The goal is to identify the internal simplifier used by RGBeta itself so that
the next Ward-identity test can be performed in exactly RGBeta's convention,
rather than by imposing external component tables.
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

Print["=== RGBeta internal tensor-reduction probe ==="];
Print["Model file: ", modelFile];

all = Names["RGBeta`PackageScope`*"];

patterns = {
    "Tensor", "Group", "Invariant", "Simpl", "Contract",
    "Reduce", "Index", "Quartic", "Upsilon", "Trace", "Rep"
};

selected = Select[
    all,
    Function[name,
        AnyTrue[patterns, StringContainsQ[name, #, IgnoreCase -> True] &]
    ]
];

Print["\n--- Candidate RGBeta`PackageScope` symbols ---"];
Scan[Print, Sort[selected]];

preferredNames = {
    "RGBeta`PackageScope`QuarticTensors",
    "RGBeta`PackageScope`UpsilonQuarticTensors",
    "RGBeta`PackageScope`SimplifyGroup",
    "RGBeta`PackageScope`GroupSimplify",
    "RGBeta`PackageScope`SimplifyGroupStructure",
    "RGBeta`PackageScope`GroupStructureSimplify",
    "RGBeta`PackageScope`TensorSimplify",
    "RGBeta`PackageScope`ContractGroupIndices"
};

Print["\n--- Definitions for preferred candidates ---"];
Do[
    If[NameQ[name],
        sym = Symbol[name];
        Print["\n>>> ", name];
        Print["OwnValues: ", InputForm[OwnValues[sym]]];
        Print["DownValues: ", InputForm[DownValues[sym]]];
        Print["UpValues: ", InputForm[UpValues[sym]]];
        Print["Attributes: ", InputForm[Attributes[sym]]];
    ],
    {name, preferredNames}
];

Print["\n--- Search candidate definitions containing delS2/tGen ---"];
Do[
    sym = Symbol[name];
    defs = Join[OwnValues[sym], DownValues[sym], UpValues[sym]];
    If[
        Length @ Cases[defs, _delS2 | _tGen | _eps | _del, Infinity] > 0,
        Print["\n>>> ", name];
        Print[InputForm[defs]];
    ],
    {name, selected}
];

Print["\nDONE"];
Exit[0];
