(* Generic exact/component export helpers for T3 RGE tensor exchange. *)

T3RGEExactString[value_] := ToString[
  InputForm[FullSimplify[value]],
  CharacterEncoding -> "ASCII"
];

T3RGEFactor[name_String, component_Integer, conjugated_] := <|
  "scalar_name" -> name,
  "component" -> component,
  "conjugated" -> TrueQ[conjugated]
|>;


(* Expected physical SU(2) invariant multiplicities for the T3 scalar
   potential.  SU(2) tensor products are multiplicity-free, while identical
   bosons restrict the self-quartic to the symmetric pair channels. *)

T3RGEIndexedCouplingNames[base_String, count_Integer?NonNegative] :=
  Table[base <> "Inv" <> ToString[index], {index, count}];

T3RGEInvariantFamilyMetadata[model_Association] := Module[
  {d1, d2, scalarSelf1, scalarSelf2, portal1, portal2, cross},

  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  scalarSelf1 = Ceiling[d1/2];
  scalarSelf2 = Ceiling[d2/2];
  portal1 = Min[2, d1];
  portal2 = Min[2, d2];
  cross = Min[d1, d2];

  {
    <|"family" -> "HiggsSelf", "kind" -> "quartic",
      "expected_physical_count" -> 1, "couplings" -> {"lambdaH"}|>,
    <|"family" -> "Scalar1Self", "kind" -> "quartic",
      "expected_physical_count" -> scalarSelf1,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaS1", scalarSelf1]|>,
    <|"family" -> "Scalar2Self", "kind" -> "quartic",
      "expected_physical_count" -> scalarSelf2,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaS2", scalarSelf2]|>,
    <|"family" -> "HiggsPortal1", "kind" -> "quartic",
      "expected_physical_count" -> portal1,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaH1", portal1]|>,
    <|"family" -> "HiggsPortal2", "kind" -> "quartic",
      "expected_physical_count" -> portal2,
      "couplings" -> T3RGEIndexedCouplingNames["lambdaH2", portal2]|>,
    <|"family" -> "ScalarCross", "kind" -> "quartic",
      "expected_physical_count" -> cross,
      "couplings" -> T3RGEIndexedCouplingNames["lambda12", cross]|>,
    <|"family" -> "T3Mixing", "kind" -> "quartic",
      "expected_physical_count" -> 1, "couplings" -> {"lambdaT3"}|>,
    <|"family" -> "Yukawa1", "kind" -> "yukawa",
      "expected_physical_count" -> 1, "couplings" -> {"y1"}|>,
    <|"family" -> "Yukawa2", "kind" -> "yukawa",
      "expected_physical_count" -> 1, "couplings" -> {"y2"}|>
  }
];


(* ------------------------------------------------------------------------- *)
(* Tensor component enumeration                                              *)
(* ------------------------------------------------------------------------- *)

(* ArrayRules appends a default rule such as {_,_,_}->0.  Keep only rules
   whose left-hand side is an explicit integer tuple. *)

T3RGEComponentEntries[tensor_] := Module[{rules},
  If[!ArrayQ[tensor],
    Return[{{} -> tensor}]
  ];

  rules = ArrayRules[Normal[tensor]];

  Cases[
    rules,
    (indices_List -> value_) /;
      VectorQ[indices, IntegerQ] && value =!= 0 :>
        (indices -> value)
  ]
];


(* Reinsert the component value 1 for representation positions that were
   dropped because the field is an SU(2) singlet. *)

T3RGEFullComponentTuple[
  keptComponents_List,
  keep_List,
  numberOfFields_Integer
] := Module[{full, position},
  full = ConstantArray[1, numberOfFields];

  Do[
    position = keep[[slot]];
    full[[position]] = keptComponents[[slot]],
    {slot, Length[keep]}
  ];

  full
];


(* ------------------------------------------------------------------------- *)
(* Generic quartic component export                                          *)
(* ------------------------------------------------------------------------- *)

T3RGEQuarticEntriesFromBasis[
  basisTensor_,
  keep_List,
  factorSpecifications_List,
  coefficient_
] := Module[
  {rules, numberOfFields, fullComponents},

  numberOfFields = Length[factorSpecifications];
  rules = T3RGEComponentEntries[basisTensor];

  Map[
    Function[rule,
      fullComponents = T3RGEFullComponentTuple[
        First[rule],
        keep,
        numberOfFields
      ];

      <|
        "coefficient" -> T3RGEExactString[coefficient Last[rule]],
        "factors" -> MapThread[
          T3RGEFactor[
            #1[[1]],
            #2,
            #1[[2]]
          ] &,
          {factorSpecifications, fullComponents}
        ]
      |>
    ],
    rules
  ]
];

T3RGEHermitianQuarticEntriesFromBasis[
  basisTensor_,
  keep_List,
  factorSpecifications_List,
  coefficient_
] := Module[{forward, backwardSpecifications, backward},
  forward = T3RGEQuarticEntriesFromBasis[
    basisTensor,
    keep,
    factorSpecifications,
    coefficient
  ];

  backwardSpecifications = {
    #[[1]],
    !TrueQ[#[[2]]]
  } & /@ factorSpecifications;

  backward = T3RGEQuarticEntriesFromBasis[
    Conjugate[basisTensor],
    keep,
    backwardSpecifications,
    Conjugate[coefficient]
  ];

  Join[forward, backward]
];


(* ------------------------------------------------------------------------- *)
(* Scalar-potential sectors                                                  *)
(* ------------------------------------------------------------------------- *)
