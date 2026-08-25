
load = UsingFrontEnd[Needs["Matchete`"]; True];
If[load =!= True, Print["Matchete load failed"]; Exit[1]];

ResetAll[];
LSM = LoadModel["SM"];

DefineRepresentation[T3Probe3, SU2L, {2}, IndexAlphabet -> {"a","b","c","d"}];
DefineRepresentation[T3Probe4, SU2L, {3}, IndexAlphabet -> {"u","v","w","x"}];

Print["Registered custom reps:"];
Print[InputForm@GetRepresentations[]];

algVariants = {
  "plain" -> {{1}, {3}, {2}},
  "C2"    -> {CRep[{1}], {3}, {2}},
  "C4"    -> {{1}, CRep[{3}], {2}},
  "C2C4"  -> {CRep[{1}], CRep[{3}], {2}}
};

cgVariants = {
  "plain" -> {SU2L[fund], T3Probe4, T3Probe3},
  "bar2"  -> {Bar[SU2L[fund]], T3Probe4, T3Probe3},
  "bar4"  -> {SU2L[fund], Bar[T3Probe4], T3Probe3},
  "bar2bar4" -> {Bar[SU2L[fund]], Bar[T3Probe4], T3Probe3}
};

counter = 0;
passes = {};

Do[
  algName = First[algRule];
  algReps = Last[algRule];

  Print[""];
  Print["ALG ", algName, " = ", InputForm[algReps]];

  tensors = CheckAbort[
    Check[InvariantTensors[SU[2], algReps], $Failed],
    $Failed
  ];
  If[tensors === $Failed,
    Print["  InvariantTensors FAILED"];
    Continue[]
  ];

  Print["  invariant count = ", Length[tensors]];

  Do[
    cgName = First[cgRule];
    cgReps = Last[cgRule];
    counter++;
    sym = Symbol["ProbeCG" <> ToString[counter]];

    Print["  TRY ", algName, " / ", cgName,
      "  reps=", InputForm[cgReps]];

    status = CheckAbort[
      Check[
        UsingFrontEnd[
          DefineCG[sym, cgReps, First[tensors]];
          True
        ],
        False
      ],
      False
    ];

    Print["    DefineCG -> ", status];

    If[TrueQ[status],
      AppendTo[passes, {algName, cgName, sym}];
      Print["    GetCGTensor = ", InputForm@GetCGTensor[sym]];
    ];
  , {cgRule, cgVariants}]
, {algRule, algVariants}];

Print[""];
Print["============================================================"];
Print["PASSING CONVENTIONS"];
Print["============================================================"];
If[passes === {},
  Print["NONE"],
  Scan[Print[InputForm[#]] &, passes]
];
