matcheteLoadResult = UsingFrontEnd[Needs["Matchete`"]; True];
If[matcheteLoadResult =!= True,
  Print["Matchete load failed."];
  Exit[1];
];

ResetAll[];
LSM = LoadModel["SM"];

Print["GetGaugeGroups[] = ", InputForm@GetGaugeGroups[]];
Print["GetGroups[] = ", InputForm@GetGroups[]];
Print["GetRepresentations[] BEFORE = ", InputForm@GetRepresentations[]];
Print["Options[DefineRepresentation] = ", InputForm@Options[DefineRepresentation]];

Print["Trying DefineRepresentation[T3Probe4, SU2L, {3}] ..."];
status = Check[
  DefineRepresentation[
    T3Probe4,
    SU2L,
    {3},
    IndexAlphabet -> {"u","v","w","x"}
  ];
  "CALL_RETURNED",
  $Failed
];
Print["status = ", InputForm@status];
Print["GetRepresentations[] AFTER = ", InputForm@GetRepresentations[]];

Print["Trying a field with SU2L[T3Probe4] ..."];
fieldStatus = Check[
  DefineField[
    T3ProbeF,
    Fermion,
    Indices -> {SU2L[T3Probe4]},
    Charges -> {U1Y -> 0},
    Mass -> MProbe
  ];
  "FIELD_CALL_RETURNED",
  $Failed
];
Print["fieldStatus = ", InputForm@fieldStatus];
Print["GetFields[] = ", InputForm@GetFields[]];

