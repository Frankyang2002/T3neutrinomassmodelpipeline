(* T3ModelCatalog.wl
   Known T3 representation classes in the convention Q = T3 + Y.
   Restrepo et al. use a hypercharge label that is twice this Y, so their
   alpha is converted as alpha/2 below.
*)
ClearAll[T3ModelFromClass];

(* Matchete flavor indices require dimension > 1.  We use three heavy-fermion
   generations by default, which is also the standard phenomenological choice
   for obtaining a full neutrino-mass matrix. *)

T3ModelFromClass[class_String, alpha_Integer, multiplicity_: 3] := Module[
  {dims, y1, y2, yf},
  dims = Switch[
    ToUpperCase[class],
    "A", {1, 3, 2},
    "B", {2, 2, 1},
    "C", {2, 2, 3},
    "D", {3, 1, 2},
    "E", {3, 3, 2},
    _, Return[$Failed]
  ];
  y1 = alpha/2;
  y2 = (alpha + 2)/2;
  yf = (alpha + 1)/2;
  <|
    "Class" -> "T3-" <> ToUpperCase[class],
    "Alpha" -> alpha,
    "Lepton" -> <|"Type" -> "Fermion", "SU2" -> 2, "Y" -> -1/2|>,
    "Higgs" -> <|"Type" -> "ComplexScalar", "SU2" -> 2, "Y" -> 1/2|>,
    "Scalar1" -> <|"Type" -> "ComplexScalar", "SU2" -> dims[[1]], "Y" -> y1, "MassSymbol" -> MS1|>,
    "Scalar2" -> <|"Type" -> "ComplexScalar", "SU2" -> dims[[2]], "Y" -> y2, "MassSymbol" -> MS2|>,
    "Fermion" -> <|"Type" -> "Fermion", "SU2" -> dims[[3]], "Y" -> yf, "Multiplicity" -> multiplicity, "MassSymbol" -> MF|>
  |>
];
