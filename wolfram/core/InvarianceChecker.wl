(*
  InvarianceChecker.wl

  Higher-level interaction invariance checks.  Electroweak U(1)_Y and SU(2)
  calculations are provided by GaugeInvariance.wl; this file adds the Z2
  condition used by the radiative-neutrino-mass models.

  Loading this file preserves the previous public behaviour:
      TemplateInvariantQ = U(1)_Y AND SU(2) AND Z2.
*)

(* Load the canonical gauge-invariance implementation from the same folder. *)
If[
  Length[DownValues[GaugeInvarianceReport]] == 0,
  Get[FileNameJoin[{DirectoryName[$InputFileName], "GaugeInvariance.wl"}]]
];

ClearAll[
  Z2InvariantQ,
  TemplateInvarianceReport,
  TemplateInvariantQ
];

(* ============================================================ *)
(* Z2                                                           *)
(* ============================================================ *)

Z2InvariantQ[fields_List] := Module[
  {charges = Lookup[fields, "Z2", Missing["Z2Charge"]]},
  If[AnyTrue[charges, MissingQ], Return[False]];
  TrueQ[Times @@ charges == 1]
];

(* ============================================================ *)
(* Full interaction report                                      *)
(* ============================================================ *)

(*
  Extend the gauge report rather than recomputing hypercharges or SU(2)
  products.  An interaction of Z2 eigenfields is allowed when the product of
  their Z2 parities is +1.
*)
TemplateInvarianceReport[
  template_Association,
  fieldData_Association
] := Module[
  {
    gaugeReport,
    resolvedFields,
    z2Charges,
    z2Product,
    z2Invariant,
    fullInvariant
  },

  gaugeReport = GaugeInvarianceReport[template, fieldData];
  resolvedFields = Lookup[gaugeReport, "ResolvedFields", {}];

  z2Charges = Lookup[resolvedFields, "Z2", Missing["Z2Charge"]];
  z2Product = If[
    resolvedFields === {} || AnyTrue[z2Charges, MissingQ],
    Missing["NotAvailable"],
    Times @@ z2Charges
  ];
  z2Invariant = !MissingQ[z2Product] && TrueQ[z2Product == 1];
  fullInvariant = TrueQ[Lookup[gaugeReport, "Invariant", False]] && z2Invariant;

  (* Keep the historical report schema/order for downstream consumers. *)
  <|
    "Name" -> Lookup[gaugeReport, "Name", "<unnamed>"],
    "ResolvedFields" -> resolvedFields,
    "HyperchargeSum" -> Lookup[gaugeReport, "HyperchargeSum", Missing["NotAvailable"]],
    "Z2Product" -> z2Product,
    "SU2ProductRepresentations" -> Lookup[gaugeReport, "SU2ProductRepresentations", {}],
    "U1Invariant" -> TrueQ[Lookup[gaugeReport, "U1Invariant", False]],
    "Z2Invariant" -> z2Invariant,
    "SU2Invariant" -> TrueQ[Lookup[gaugeReport, "SU2Invariant", False]],
    "Invariant" -> fullInvariant
  |>
];

TemplateInvariantQ[
  template_Association,
  fieldData_Association
] := TrueQ[
  TemplateInvarianceReport[template, fieldData]["Invariant"]
];
