(* Direct probe of RGBeta conjugation machinery.
   This version avoids helper-symbol evaluation bugs by referencing each
   symbol explicitly.

Run from repository root:
  wolframscript -file RGE/group_factors/ProbeRGBetaGroupInvBarDirect.wl
*)

ClearAll["Global`*"];
Needs["RGBeta`"];

Print["=== RGBeta direct GroupInvBar / Bar probe ==="];

Print["\n=== Contexts ==="];
Print["Context[Bar] = ", Context[Bar]];
Print["Context[BarToConjugate] = ", Context[BarToConjugate]];
Print["Context[RGBeta`PackageScope`CouplingBar] = ",
      Context[RGBeta`PackageScope`CouplingBar]];
Print["Context[RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar] = ",
      Context[RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar]];

Print["\n=== GroupInvBar ==="];
Print[
  InputForm[
    DownValues[
      RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar
    ]
  ]
];
Print[
  InputForm[
    OwnValues[
      RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar
    ]
  ]
];

Print["\n=== Bar ==="];
Print[InputForm[DownValues[Bar]]];
Print[InputForm[OwnValues[Bar]]];
Print[InputForm[UpValues[Bar]]];

Print["\n=== BarToConjugate ==="];
Print[InputForm[DownValues[BarToConjugate]]];
Print[InputForm[OwnValues[BarToConjugate]]];
Print[InputForm[UpValues[BarToConjugate]]];

Print["\n=== CouplingBar ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`CouplingBar]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`CouplingBar]]];
Print[InputForm[UpValues[RGBeta`PackageScope`CouplingBar]]];

Print["\n=== TStructure ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`TStructure]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`TStructure]]];

Print["\n=== TsSym4 ==="];
Print[InputForm[DownValues[RGBeta`PackageScope`TsSym4]]];
Print[InputForm[OwnValues[RGBeta`PackageScope`TsSym4]]];

Print["\n=== Attributes ==="];
Print["Attributes[Bar] = ", InputForm[Attributes[Bar]]];
Print["Attributes[BarToConjugate] = ", InputForm[Attributes[BarToConjugate]]];
Print["Attributes[GroupInvBar] = ",
      InputForm[
        Attributes[
          RGBeta`FieldsAndCouplings`PackagePrivate`GroupInvBar
        ]
      ]];

Print["\nDONE"];
Exit[0];
