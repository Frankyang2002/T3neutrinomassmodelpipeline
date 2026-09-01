(* ProbeT3MultiInvariantPortal.wl

   Proof-of-concept that Matchete can simultaneously keep multiple independent
   SU(2) invariant tensors for one scalar-potential field structure.

   Intended location:
     wolfram/tests/ProbeT3MultiInvariantPortal.wl

   Run from the project root with:
     wolframscript -file wolfram/tests/ProbeT3MultiInvariantPortal.wl
*)

ClearAll["Global`*"];

Print["Loading Matchete..."];
matcheteLoaded = UsingFrontEnd[Needs["Matchete`"]; True];
If[!TrueQ[matcheteLoaded],
  Print["ERROR: Matchete failed to load."];
  Exit[1]
];
Print["Matchete loaded successfully."];

LSM = LoadModel["SM"];
If[LSM === $Failed,
  Print["ERROR: SM model failed to load."];
  Exit[2]
];

(* A separate complex scalar doublet is enough to test the generic
   H^dagger H S^dagger S invariant space. *)
DefineField[
  PortalProbeScalar,
  Scalar,
  SelfConjugate -> False,
  Indices -> {SU2L[fund]},
  Charges -> {U1Y[0]},
  Mass -> {Heavy, MPortalProbe}
];

DefineCoupling[lambdaPortal1, SelfConjugate -> True];
DefineCoupling[lambdaPortal2, SelfConjugate -> True];

(* Field order:
     H^dagger, H, S^dagger, S

   For pseudoreal SU(2) doublets, use the same orientation convention as the
   generalised T3 builder.
*)
algebraReps = {
  CRep[{1}],
  {1},
  CRep[{1}],
  {1}
};

cgReps = {
  SU2L[fund],
  Bar[SU2L[fund]],
  SU2L[fund],
  Bar[SU2L[fund]]
};

Print[""];
Print["Requesting invariant tensors for H^dagger H S^dagger S..."];
tensors = Check[
  InvariantTensors[SU[2], algebraReps],
  $Failed
];

If[tensors === $Failed || !ListQ[tensors],
  Print["ERROR: InvariantTensors failed."];
  Exit[3]
];

Print["Invariant tensor count: ", Length[tensors]];

If[Length[tensors] < 2,
  Print["ERROR: Expected at least two invariant tensors for the doublet portal."];
  Exit[4]
];

Print["Defining PortalProbeCG1 and PortalProbeCG2..."];

define1 = Check[
  DefineCG[PortalProbeCG1, cgReps, tensors[[1]]];
  True,
  $Failed
];

define2 = Check[
  DefineCG[PortalProbeCG2, cgReps, tensors[[2]]];
  True,
  $Failed
];

Print["  CG1 definition: ", InputForm[define1]];
Print["  CG2 definition: ", InputForm[define2]];

If[MemberQ[{define1, define2}, $Failed],
  Print["ERROR: At least one CG definition failed."];
  Exit[5]
];

portal1 = lambdaPortal1[] *
  Bar[H[i]] H[j] Bar[PortalProbeScalar[a]] PortalProbeScalar[b] *
  CG[PortalProbeCG1, {i, j, a, b}];

portal2 = lambdaPortal2[] *
  Bar[H[i]] H[j] Bar[PortalProbeScalar[a]] PortalProbeScalar[b] *
  CG[PortalProbeCG2, {i, j, a, b}];

LFree = FreeLag[PortalProbeScalar] // RelabelIndices;

Print[""];
Print["Checking first invariant alone..."];
check1 = CheckAbort[
  Check[
    CheckLagrangian[(LSM + LFree + portal1) // RelabelIndices],
    $Failed
  ],
  $Aborted
];
Print["  result: ", InputForm[check1]];

Print["Checking second invariant alone..."];
check2 = CheckAbort[
  Check[
    CheckLagrangian[(LSM + LFree + portal2) // RelabelIndices],
    $Failed
  ],
  $Aborted
];
Print["  result: ", InputForm[check2]];

Print["Checking both independent invariants simultaneously..."];
checkBoth = CheckAbort[
  Check[
    CheckLagrangian[(LSM + LFree + portal1 + portal2) // RelabelIndices],
    $Failed
  ],
  $Aborted
];
Print["  result: ", InputForm[checkBoth]];

Print[""];
If[TrueQ[check1] && TrueQ[check2] && TrueQ[checkBoth],
  Print["PASS: Matchete accepts two independent portal invariant tensors with separate couplings."];
  Print["This mechanism can be promoted into the general scalar-potential builder."],
  Print["FAIL: The multi-invariant portal proof did not validate completely."];
  Print["Do not modify the production builder until this orientation/CG issue is resolved."];
  Exit[6]
];
