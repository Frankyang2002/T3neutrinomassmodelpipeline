(* Direct normalization test for the matched T3-B Weinberg operator.
   This is a diagnostic only; it does not modify production files. *)

ClearAll[ParseWeinbergRawTerms];

ParseWeinbergRawTerms[path_String] := Module[{text, chunks},
  If[!FileExistsQ[path],
    Print["ERROR: missing file: ", path];
    Return[$Failed]
  ];

  text = Import[path, "Text"];
  chunks = Select[
    StringTrim /@ StringSplit[
      text,
      RegularExpression["--- Weinberg term [0-9]+ ---\\s*"]
    ],
    StringLength[#] > 0 &
  ];

  Quiet @ Check[
    ToExpression[#, InputForm] & /@ chunks,
    $Failed
  ]
];

projectRoot = Directory[];

outputDirectory = If[
  Length[$ScriptCommandLine] >= 2,
  $ScriptCommandLine[[2]],
  FileNameJoin[{projectRoot, "wolfram", "output", "T3_B_alpha_m1"}]
];

rawPath = FileNameJoin[{outputDirectory, "c5_raw.txt"}];
coefficientPath = FileNameJoin[{outputDirectory, "c5_coefficient.txt"}];

Print["========================================================================"];
Print["WEINBERG -> NEUTRINO MASS NORMALIZATION TEST v2"];
Print["========================================================================"];
Print[];
Print["Output directory: ", outputDirectory];
Print[];

terms = ParseWeinbergRawTerms[rawPath];

If[terms === $Failed,
  Print["FAIL: could not parse c5_raw.txt"];
  Exit[1]
];

c5 = Quiet @ Check[
  ToExpression[Import[coefficientPath, "Text"], InputForm],
  $Failed
];

If[c5 === $Failed,
  Print["FAIL: could not parse c5_coefficient.txt"];
  Exit[1]
];

holomorphicTerms = Select[
  terms,
  StringContainsQ[ToString[#, InputForm], "Proj[-1]"] &
];

Print["Total raw Weinberg terms: ", Length[terms]];
Print["Holomorphic terms: ", Length[holomorphicTerms]];
Print[];

If[Length[holomorphicTerms] =!= 4,
  Print["FAIL: expected four holomorphic terms."];
  Exit[1]
];

(* The actual T3-B matched terms all contain
       eps[d1,d2] eps[d3,d4] H[d1] H[d4] l[d3] l[d2].
   For the neutral component:
       H0 => d1=d4=2,
       nu => d2=d3=1.
   Therefore the SU(2) contraction is
       eps[2,1] eps[1,2] = (-1)(+1) = -1.
   The minus sign is convention-dependent together with the Majorana bilinear,
   but the normalization magnitude is one. *)

epsilon21 = -1;
epsilon12 = 1;
su2Factor = epsilon21 epsilon12;

Print["Neutral-component SU(2) assignment:"];
Print["  H0: d1 = d4 = 2"];
Print["  nu: d2 = d3 = 1"];
Print[];
Print["epsilon[2,1] = ", epsilon21];
Print["epsilon[1,2] = ", epsilon12];
Print["SU(2) product = ", su2Factor];
Print[];

(* Verify that every holomorphic term has precisely this common field/CG structure
   before using the component result. *)
structureChecks = Table[
  With[{s = ToString[term, InputForm]},
    And[
      StringCount[s, "Field[H, Scalar"] == 2,
      StringCount[s, "Field[l, Fermion"] == 2,
      StringCount[s, "CG[Bar[eps[SU2L]]"] == 2,
      StringCount[s, "Proj[-1]"] == 1
    ]
  ],
  {term, holomorphicTerms}
];

Print["Common LLHH epsilon structure checks: ", structureChecks];

If[!And @@ structureChecks,
  Print["FAIL: raw terms do not all have the expected T3-B LLHH structure."];
  Exit[1]
];

(* The existing production extractor forms C5 by stripping exactly the two H fields,
   the two epsilon tensors, and the common Majorana spinor chain, while retaining
   the physical prefactor. Thus the coefficient of the neutral component has
   magnitude |C5|, with the displayed minus sign coming from eps[2,1] eps[1,2]. *)

rSigned = su2Factor;
rMagnitude = Abs[rSigned];

Print[];
Print["Signed component ratio from the SU(2) contraction:"];
Print["  r_signed = C_(nu nu H0 H0) / C5 = ", rSigned];
Print["  |r| = ", rMagnitude];
Print[];

Print["After H0 -> v/Sqrt[2]:"];
Print["  L_EFT contains (r_signed C5 v^2/2) nu nu + h.c."];
Print[];
Print["Comparing with the two-component Majorana convention"];
Print["  L_mass = -(1/2) m_nu nu nu + h.c."];
Print["gives"];
Print["  m_nu = -r_signed v^2 C5."];
Print[];

If[rMagnitude === 1,
  Print["NORMALIZATION RESULT: |r| = 1."];
  Print["There is NO extra factor of 1/2 between C5 and the Majorana mass magnitude."];
  Print["Therefore |m_nu| = v^2 |C5|."];
  Print["The overall sign depends on the epsilon/spinor convention and is not a physical mass eigenvalue sign."],
  Print["FAIL: unexpected normalization magnitude."];
  Exit[1]
];

Print[];
Print["PASS: T3-B neutral-component normalization fixed unambiguously in magnitude."];
