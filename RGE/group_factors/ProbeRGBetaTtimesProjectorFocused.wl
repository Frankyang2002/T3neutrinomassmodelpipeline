(* Focused RGBeta contraction/projector probe.
   Avoids printing the enormous BetaTensor definition table.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaTtimesProjectorFocused.wl
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

Print["=== Focused RGBeta Ttimes / Projector probe ==="];

Print["\n=== Ttimes ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`Ttimes]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`Ttimes]]];
Print[InputForm[SubValues[RGBeta`PackageScope`Ttimes]]];

Print["\n=== Projector ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`Projector]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`Projector]]];
Print[InputForm[SubValues[RGBeta`PackageScope`Projector]]];

Print["\n=== UpdateProjectors ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`UpdateProjectors]]];

Print["\n=== $quartics own/down values ==="];
Print[InputForm[OwnValues[RGBeta`PackageScope`$quartics]]];
Print[InputForm[DownValues[RGBeta`PackageScope`$quartics]]];

Print["\n=== ScalarFieldProjector ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`TensorCalculations`PackagePrivate`ScalarFieldProjector
    ]
  ]
];

Print["\n=== FieldsAndCouplings private projector symbols ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`FieldsAndCouplings`PackagePrivate`projector
    ]
  ]
];
Print[
  InputForm[
    OwnValues[
      RGBeta`FieldsAndCouplings`PackagePrivate`projector
    ]
  ]
];

Print["\nDONE"];
Exit[0];
