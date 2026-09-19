(* TLDR: 
- Take in SU2 and U1 representations and change into Matchete fields. 
- Create all gauge invariant interactions
- Return UV Lagrangian

These are all helper functions for the run model

Note that we use invariant tensor function from group magic in matchete. 
We can explicitly tell it the symmetric fields so that we get the tensor that we can SU2 contracts fields correctly
It returns explicit CG tensors in a list so like HH with 2x2=3+1 would be only the CG for the symmetric part {H1,H2,H3} 
*)
ClearAll[
  ValidateT3Candidate,
  SelectFirstValidByGroup,
  T3IngredientsPresentQ,
  BuildT3Lagrangian
];

(* SU(2) representation and topology-CG helpers. *)
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName], "SU2Invariants.wl"}]];

(* We define couplings for matchete *)


(* We convert into matchete fields 
1. We define our 2 scalar and 1 fermion fields to be used as NewFermion and NewScalar when called*)


(* 1.We fill our scalar field value for our 2 scalars *)


(* Scalar-potential invariant and quartic-candidate helpers. *)
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName], "ScalarPotential.wl"}]];

(* T3 BSM fields and couplings. *)
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName], "T3Fields.wl"}]];


(* Now H H S1 S2^dagger + h.c. vertex  The invariant tensor is
   generated for the exact requested representations and is symmetric in the
   two Higgs indices. 
   1. Create our new fields and conjugate the correct scalar
   2. Label indices if needed
   3. Put into CG, we name the coefficient and give it indices*)


(* T3 topology interaction candidates. *)
Get[FileNameJoin[{DirectoryName @ ExpandFileName[$InputFileName], "T3Topology.wl"}]];

(* We check with Matchete if the lagrangian with the interaction is valid 
1. use CheckLagrangian to see if Lagrangian is valid with the interaction and indices*)
(* ------------------------------ *)
(* Candidate validation           *)
(* ------------------------------ *)

ValidateT3Candidate[c_, LSM_, LFree_] := Module[{res},
  Print["Checking candidate: ", c["Name"]];
  (* Relabel indices is to prevent dummy index clashes, 
  like for A_i B_i + C_i D_i  are independent despite both using i 
  We prevent reusing index labels and throughout the file*)
  
  res = CheckAbort[
    Check[CheckLagrangian[(LSM + LFree + c["Expression"]) // RelabelIndices], $Failed],
    $Aborted
  ];
  Print["  result: ", InputForm[res]];
  Association[c, "Valid" -> TrueQ[res], "ValidationResult" -> res]
];

(* Keep the first valid contraction in each alternative group.  Candidates
   without a group are independent interactions and are all retained. 
   1. Select first valid contraction of many valid contractions to represent these fields*)
SelectFirstValidByGroup[list_List] := Module[{seen = <||>},
  Select[
    list,
    Function[candidate,
      Module[{group = Lookup[candidate, "AlternativeGroup", None]},
        If[group === None,
          True,
          If[KeyExistsQ[seen, group], False, seen[group] = True; True]
        ]
      ]
    ]
  ]
];

(* Assemble all interaction candidates without validating them.
1. Register our 2 fields as they both are valid for yukawa
2. Find all interactions vertices that come from this model, matching our field with all other fields
3. The build functions try to contract the fields correctly and check if they are valid in the Lagrangian
4. We add these interactions into a list as candidates for all candidates*)


(* A genuine T3 contribution requires both lepton-fermion-scalar Yukawas and
   the quartic that connects S1 and S2 to the two Higgs legs. 
   No weinberg without these 3 elements 
   1. Check if all 3 vertices exist*)
T3IngredientsPresentQ[valid_List] := And @@ (
  Function[class,
    AnyTrue[valid, Lookup[#, "Class", ""] === class &]
  ] /@ {"Yukawa1", "Yukawa2", "T3ScalarMix"}
);

(* We do everything here  
1. We load SM fields
2. We Define our BSM fields
3. We Define couplings
4. We define our Clebsch Gordon coefficients
5. We get our free Lagrangian
6. We get all interactions and validate them, if they are valid, we choose the first contraction that works for that combination of fields 
7. *)
(* ------------------------------ *)
(* UV model assembly              *)
(* ------------------------------ *)

BuildT3Lagrangian[model_Association] := Module[
  {LSM, LFree, cgStatus, candidates, validated, valid, rejected,
   LInt, LBSM, LUV, fullValidation, ingredientsPresent},

  ResetAll[];
Print["Build stage 1/4: load SM"];
  LSM = LoadModel["SM"];
  Print["  LoadModel result: ", If[LSM === $Failed, "$Failed", "OK"]];
  If[LSM === $Failed, Return[$Failed]];

  Print["Build stage 2/4: define BSM fields"];
  If[DefineT3Fields[model] === $Failed,
    Print["  field definition result: $Failed"];
    Return[$Failed],
    Print["  field definition result: True"]
  ];

  Print["Build stage 3/4: define couplings"];
  If[Check[DefineT3Couplings[], $Failed] === $Failed,
    Print["  coupling definition result: $Failed"];
    Return[$Failed],
    Print["  coupling definition result: True"]
  ];
  Print["Build stage 4/4: define SU(2) CG tensors"];
  cgStatus = DefineT3InvariantCGs[model];
  Print["  CG definition result: ", InputForm[cgStatus]];
  If[cgStatus === $Failed, Return[$Failed]];


  (* Get free Lagrangian *)
  LFree = FreeLag[NewFermion, NewScalar1, NewScalar2] // RelabelIndices;

  candidates = BuildT3InteractionCandidates[model];
  candidates = Association[#, "Expression" -> RelabelIndices[#["Expression"]]] & /@ candidates;

  (* We validate each interaction individually, seeing if they are compatible *)
  validated = ValidateT3Candidate[#, LSM, LFree] & /@ candidates;
  valid = SelectFirstValidByGroup @ Select[validated, TrueQ[#["Valid"]] &];
  rejected = Select[validated, !TrueQ[#["Valid"]] &];

  (* Build interaction lagrangian and then the UV *)
  LInt = Total[Lookup[valid, "Expression", {}]] // Expand // RelabelIndices;
  LBSM = (LFree + LInt) // Expand // RelabelIndices;
  LUV = (LSM + LBSM) // Expand // RelabelIndices;
  fullValidation = CheckAbort[Check[CheckLagrangian[LUV], $Failed], $Aborted];

  (* Check if it has all the vertices for T3 diagram *)
  ingredientsPresent = T3IngredientsPresentQ[valid];

  <|
    "Status" -> If[TrueQ[fullValidation], "Success", "FullValidationFailed"],
    "Model" -> model,
    "LSM" -> LSM,
    "LFree" -> LFree,
    "LInt" -> LInt,
    "LBSM" -> LBSM,
    "LUV" -> LUV,
    "AllowedInteractions" -> Lookup[valid, "Name", {}],
    "RejectedInteractions" -> Lookup[rejected, "Name", {}],
    "T3IngredientsPresent" -> ingredientsPresent,
    "WeinbergIngredientsPresent" -> ingredientsPresent,
    "FullValidation" -> fullValidation
  |>
];
