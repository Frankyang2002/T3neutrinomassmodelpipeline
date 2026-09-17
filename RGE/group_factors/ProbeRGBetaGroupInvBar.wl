(* Probe RGBeta's conjugation map used to build quartic projectors.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaGroupInvBar.wl

The quartic UpdateProjectors branch builds the dual projector from
GroupInvBar[Invariant[...]] evaluated on Bar /@ Fields, then normalizes it by
its overlap with Lam.  This probe inspects exactly that conjugation machinery.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

Print["=== RGBeta GroupInvBar / conjugation probe ==="];

symbols = {
    "RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar",
    "RGBeta`PackageScope`GroupInvBar",
    "RGBeta`FieldsAndCouplings`PackagePrivate`Bar",
    "RGBeta`PackageScope`Bar",
    "RGBeta`FieldsAndCouplings`PackagePrivate`TStructure",
    "RGBeta`PackageScope`TStructure",
    "RGBeta`PackageScope`TsSym4"
};

Do[
    Print["\n>>> ", name];
    If[
        NameQ[name],
        s = ToExpression[name, InputForm, HoldComplete];
        sym = First[s];
        Print["OwnValues: ", InputForm[OwnValues[sym]]];
        Print["DownValues: ", InputForm[DownValues[sym]]];
        Print["UpValues: ", InputForm[UpValues[sym]]];
        Print["SubValues: ", InputForm[SubValues[sym]]],
        Print["NOT FOUND"]
    ],
    {name, symbols}
];

Print["\n=== Names containing GroupInvBar / InvBar / Conjug ==="];
ctx = Select[Contexts[], StringStartsQ[#, "RGBeta`"] &];
all = DeleteDuplicates @ Flatten[Names[# <> "*"] & /@ ctx];
sel = Select[
  all,
  Function[name,
    AnyTrue[
      {"GroupInvBar", "InvBar", "Conjug", "Bar"},
      StringContainsQ[name, #, IgnoreCase -> True] &
    ]
  ]
];
Scan[Print, Sort[sel]];

Print["\nDONE"];
Exit[0];
