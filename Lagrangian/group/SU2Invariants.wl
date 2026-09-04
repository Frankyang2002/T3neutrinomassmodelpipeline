(* SU(2) representation registration and topology invariant CG machinery.
   Extracted from LagrangianBuilder.wl; physics behaviour is unchanged. *)

(* ------------------------------------------------- *)
(* Getting our SU2 and BSM representations           *)
(* ------------------------------------------------- *)

(* For SU(2), the irrep of dimension d has highest-weight Dynkin label {d-1}*)
SU2DynkinLabel[d_Integer?Positive] := {d - 1};

(* For matchete our SM fields remain normal but BSM fields obtain its own representation to make CG work, 
we have 2 representations, one for SM doublet basis and another for BSM representation basis
This is also true for d=2 and d=3*)


(* Get BSMrep name *)
BSMRepresentationName[1] := None;
BSMRepresentationName[d_Integer?Positive] /; d >= 2 :=
  Symbol["T3BSMd" <> ToString[d]];

(* We log our representation with its respective Dynkin label *)
EnsureBSMRepresentation[1] := True;
EnsureBSMRepresentation[d_Integer?Positive] /; d >= 2 := Module[
  {rep, before, after, defineStatus},

  rep = BSMRepresentationName[d];

  before = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[before] && KeyExistsQ[before, rep], Return[True]];

  defineStatus = Check[
    DefineRepresentation[
      rep,
      SU2L,
      SU2DynkinLabel[d],
      IndexAlphabet -> {"u", "v", "w", "x", "y", "z"}
    ];
    True,
    $Failed
  ];

  If[defineStatus === $Failed,
    Print["  FAILED to define SU(2) representation d=", d,
      " (Dynkin ", SU2DynkinLabel[d], ")"];
    Return[$Failed]
  ];

  (* Never claim success unless Matchete actually registered the irrep. *)
  after = Quiet@Check[GetRepresentations[], <||>];
  If[AssociationQ[after] && KeyExistsQ[after, rep],
    Print["  Defined BSM SU(2) representation ", rep,
      " with d=", d, ", Dynkin=", SU2DynkinLabel[d]];
    True,
    Print["  ERROR: DefineRepresentation returned without registering ", rep];
    $Failed
  ]
];

(* Dimension is now an index that relates our dimension to our logged dynkin labels *)
BSMIndexType[1] := None;
BSMIndexType[d_Integer?Positive] := Module[{ok = EnsureBSMRepresentation[d]},
  If[ok === $Failed, Return[Missing["UnsupportedSU2", d]]];
  BSMRepresentationName[d]
];

(* Same thing, but for generic SU2, like for d=2 we have fund, for d=3 we have adj *)




(* Given field and SU(2) rep, get invariant Clebsch Gordon tensor to contract indices into SU(2) singlet from Matchete 
Our input is the Clebsch Gordon name, like T3Y1CG, whether fields in terms such as LFS is conjugated, so True,True,False would mean L is conj, F is conj, and S isnt
symmetricPositions tells us if we have symmetry in term, like for HHSS, the 2 S terms are symmetric, thus removing antisymmetric spaces
smPositions tells us what position a field is a standard model, where we can use in built Matchete representation for them instead of BSM
*)

(*Essentially: Dimensions -> Dynkin Reps -> Invariant Tensors -> CG that can be used in Matchete
1. We *)
DefineSU2InvariantCG[
  cgName_Symbol,
  dims_List,
  objectConjugated_List,
  symmetricPositions_: {},
  smPositions_: {}
] := Module[
  {keep, keptDims, keptConj, algebraReps, cgReps, tensors, symmetryAfterDrop, positionMap},


  If[Length[dims] =!= Length[objectConjugated], Return[$Failed]]; (* For each field we need to say if it is conjugated or not *)
  If[!AllTrue[dims, IntegerQ[#] && Positive[#] &], Return[$Failed]]; (* Only integer positive dimensions *)

  (* For each field we have, we need to make sure we have a representation for its dimension, also deleting repeated *)
  Scan[EnsureBSMRepresentation, DeleteDuplicates[dims]]; 

  (* Remove singlets as they dont have indices, create new arrays with what is left over *)
  keep = Flatten@Position[dims, _?(# > 1 &)]; 
  keptDims = dims[[keep]];
  keptConj = objectConjugated[[keep]];

  (* No CG if we have a singlet *)
  If[keptDims === {}, Return[None]];

  (* Matchete distinguishes the orientation of pseudoreal SU(2) indices.
    Note pseudoreal means that our representation is equivalent to its conjugate,
    but its indices still need to be distringuished as the conjugate has a different CG
    We need to use CRep to represent the conjugation of a representation
       - a conjugated pseudoreal field is represented by CRep[label] in
         InvariantTensors and by an UNBARRED representation in DefineCG;
       - an unconjugated pseudoreal field uses the plain Dynkin label in
         InvariantTensors and the BARRED representation in DefineCG.
     For odd-dimensional full integer isospin SU(2) irreps the representation is real, so no CRep/Bar distinction is required.
    Objectconjugated/keptConj tells us if the field object is gauge-conjugated *)

  (* The following uses this, where unconjugated becomes {3} and conjugated will be CRep[{3}], this is what we need to do for Machete *)
  algebraReps = MapThread[
    Function[{d, conjugated},
      Module[{label = SU2DynkinLabel[d]},
        If[EvenQ[d] && TrueQ[conjugated],
          CRep[label],
          label
        ]
      ]
    ],
    {keptDims, keptConj}
  ];

  (* Matchete's Invariant Tensors returns tensors with indices that transform in conjugate representation 
    So tensor index transforms opposite to the algebra field
  *)
  cgReps = MapThread[
    Function[{originalPosition, d, conjugated},
      Module[{baseRep},
        baseRep = If[MemberQ[smPositions, originalPosition],
          (* SM doublets use Matchete's built-in fundamental representation. *)
          SU2L[fund],
          BSMRepresentationName[d]
        ];

        If[EvenQ[d] && !TrueQ[conjugated],
          Bar[baseRep],
          baseRep
        ]
      ]
    ],
    {keep, keptDims, keptConj}
  ];

  (* Association thread creates a dictionary of our keep with its the integer position of everything in keep
  so if keep = {a,b,c,d}, range,length gives {1,2,3,4}. This gives us a position map, as we map a->1, b->2 etc
  with association thread*)
  positionMap = AssociationThread[keep -> Range[Length[keep]]];


  (* We want to prevent antisymmetric tensor products when we have symmetric fields *)
  symmetryAfterDrop = Select[
    Lookup[positionMap, #, Missing["Dropped"]] & /@ symmetricPositions,
    IntegerQ
  ];

  Print["  CG ", SymbolName[cgName], ": dims=", dims,
    ", oriented algebra reps=", InputForm[algebraReps],
    ", CG reps=", InputForm[cgReps],
    ", symmetric positions=", symmetryAfterDrop];

  (* We obtain Invariant tensors, noting that symmetry is tracked *)
  tensors = Check[
    If[symmetryAfterDrop === {},
      InvariantTensors[SU[2], algebraReps],
      InvariantTensors[SU[2], algebraReps, SymmetricIndices -> symmetryAfterDrop]
    ],
    $Failed
  ];

  Print["    InvariantTensors result head: ", Head[tensors],
    If[ListQ[tensors], ", count=" <> ToString[Length[tensors]], ""]];

  (* If no tensors we failed *)
  If[tensors === $Failed || !ListQ[tensors] || Length[tensors] === 0,
    Print["    FAILED while generating invariant tensor for ", SymbolName[cgName]];
    Return[$Failed]
  ];

  (* A topology vertex is defined by one named coupling and therefore must
     have exactly one physical invariant.  Never discard or merge additional
     invariant tensors silently.  Scalar-potential sectors with genuine
     multiplicity use DefineSU2InvariantFamily below instead. *)
  If[Length[tensors] > 1,
    Print["ERROR: ", SymbolName[cgName], " has ", Length[tensors],
      " invariant tensors but this topology vertex expects exactly one."];
    Return[$Failed]
  ];

  Module[{defined},
    defined = Check[
      DefineCG[cgName, cgReps, First[tensors]]; (* < We get the CG coefficients from invariant tensors *)
      cgName,
      $Failed
    ];
    Print["    DefineCG result: ", InputForm[defined]];
    defined
  ]
];

(* Register the three topology-defining invariant tensors for our 3 fields in the model.  Their field order is the same order used in the interaction. 
 We get invariant CG from each our interactions vertices, and define it for the model *)
DefineT3InvariantCGs[model_Association] := Module[
  {d1, d2, dF, y1cg, y2cg, mixcg},
  
  d1 = model["Scalar1", "SU2"];
  d2 = model["Scalar2", "SU2"];
  dF = model["Fermion", "SU2"];

  (* Bar[L] F^c S1 *)
  y1cg = DefineSU2InvariantCG[
    T3Y1CG, {2, dF, d1}, {True, True, False}, {}, {1}
  ];

  (* Bar[L] F S2^dagger *)
  y2cg = DefineSU2InvariantCG[
    T3Y2CG, {2, dF, d2}, {True, False, True}, {}, {1}
  ];

  (* H H S1 S2^dagger.  The two H fields are identical bosons, so force the
     symmetric H-H channel. *)
  mixcg = DefineSU2InvariantCG[
    T3MixCG,
    {2, 2, d1, d2},
    {False, False, False, True},
    {1, 2},
    {1, 2}
  ];

  If[MemberQ[{y1cg, y2cg, mixcg}, $Failed], Return[$Failed]];
  <|"Yukawa1" -> y1cg, "Yukawa2" -> y2cg, "Mixing" -> mixcg|>
];
