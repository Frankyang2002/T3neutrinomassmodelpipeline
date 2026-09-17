(* Probe RGBeta's tensor-contraction and quartic-projector machinery.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaTtimesProjector.wl

This targets the machinery actually used by QuarticTensors:
  Ttimes[$quartics[coupling, Projector][...], BetaTensor[...]]
rather than the unused/empty $scalarContraction symbol.
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

Print["=== RGBeta Ttimes / Projector probe ==="];

Print["\n=== Ttimes ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`Ttimes]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`Ttimes]]];
Print[InputForm[SubValues[RGBeta`PackageScope`Ttimes]]];

Print["\n=== Projector ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`Projector]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`Projector]]];
Print[InputForm[SubValues[RGBeta`PackageScope`Projector]]];

Print["\n=== $quartics ==="];
Print[InputForm[OwnValues[RGBeta`PackageScope`$quartics]]];
Print[InputForm[DownValues[RGBeta`PackageScope`$quartics]]];

Print["\n=== Lam ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`Lam]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`Lam]]];

Print["\n=== TsSym4 ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`TsSym4]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`TsSym4]]];

Print["\n=== BetaTensor ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`TensorCalculations`PackagePrivate`BetaTensor
    ]
  ]
];

Print["\n=== ScalarFieldProjector ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`TensorCalculations`PackagePrivate`ScalarFieldProjector
    ]
  ]
];

Print["\n=== Names containing Ttimes/Projector/Contraction ==="];
ctx = Select[Contexts[], StringStartsQ[#, "RGBeta`"] &];
all = DeleteDuplicates @ Flatten[Names[# <> "*"] & /@ ctx];
sel = Select[
  all,
  Function[name,
    AnyTrue[
      {"Ttimes", "Projector", "Contraction", "Contract", "Scalar"},
      StringContainsQ[name, #, IgnoreCase -> True] &
    ]
  ]
];
Scan[Print, Sort[sel]];

Print["\nDONE"];
Exit[0];
