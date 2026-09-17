(* Direct probe of RGBeta internal definitions.

Run from the repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaDefinitionsDirect.wl
*)

ClearAll["Global`*"];

Needs["RGBeta`"];

Print["=== RGBeta direct internal-definition probe ==="];

Print["\n=== RefineGroupStructures ==="];
Print[
  InputForm[
    DownValues[RGBeta`PackageScope`RefineGroupStructures]
  ]
];

Print["\n=== $scalarContraction ==="];
Print[
  InputForm[
    OwnValues[RGBeta`PackageScope`$scalarContraction]
  ]
];

Print["\n=== QuarticTensors ==="];
Print[
  InputForm[
    DownValues[RGBeta`PackageScope`QuarticTensors]
  ]
];

Print["\n=== UpsilonQuarticTensors ==="];
Print[
  InputForm[
    DownValues[RGBeta`PackageScope`UpsilonQuarticTensors]
  ]
];

Print["\n=== twoIndexRepDelta ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`GroupsAndIndices`PackagePrivate`twoIndexRepDelta
    ]
  ]
];

Print["\n=== CompleteReplace ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`TensorCalculations`PackagePrivate`CompleteReplace
    ]
  ]
];

Print["\n=== PerformGeneratorTraces ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`GroupsAndIndices`PackagePrivate`PerformGeneratorTraces
    ]
  ]
];

Print["\n=== PerformAdjAlg ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`GroupsAndIndices`PackagePrivate`PerformAdjAlg
    ]
  ]
];

Print["\n=== AdjContraction ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`GroupsAndIndices`PackagePrivate`AdjContraction
    ]
  ]
];

Print["\n=== replace rule table ==="];
Print[
  InputForm[
    OwnValues[
      RGBeta`GroupsAndIndices`PackagePrivate`replace
    ]
  ]
];

Print["\nDONE"];
Exit[0];
