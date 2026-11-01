(* Evaluate Matchete's exported InputForm C5 without loading Matchete.
   In this numerical context Coupling[...] and Bar[...] are inert syntax;
   loading Matchete can introduce premature evaluation and FE messages.
   CLI: coefficient.txt parameter_points.json result.json *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3, Print["Usage: coefficient.txt points.json result.json"]; Exit[2]];
raw = Import[args[[1]], "Text"];
(* Do NOT Needs["Matchete`"] here: the coefficient is already matched.
   Explicitly inert symbolic heads are enough for the numerical replacement. *)
expr = Quiet @ Check[ToExpression[raw, InputForm], $Failed];
If[expr === $Failed, Print["ERROR: cannot parse input C5"]; Exit[3]];
payload = Import[args[[2]], "RawJSON"];
If[!AssociationQ[payload] || !ListQ[Lookup[payload,"points",None]],
 Print["ERROR: malformed parameter points"]; Exit[3]];
complexValue[x_] := If[AssociationQ[x], N[x["real"] + I x["imag"]], N[x]];
(* Never use N[...] with Check[...] as a success gate: some harmless
   numerical precision messages make Check return $Failed. *)
results = Table[
  Module[{params, channel, resolve, evaluated, numeric, valid, badCouplings,
          badSymbols, otherHeads, index},
    params = AssociationMap[complexValue, point["couplings"]];
    channel = point["channel"]; index = point["index"];
    resolve[s_Symbol] := Module[{key = SymbolName[s]},
      Which[
        key === "y1", 1,
        key === "y2", If[channel === "y2", 1, 0],
        key === "y2Conj", If[channel === "y2Conj", 1, 0],
        KeyExistsQ[params, key], params[key],
        True, Missing["UnmappedCoupling", key]
      ]
    ];
    evaluated = expr /. {
      HoldPattern[Bar[Coupling[s_Symbol, ___]]] :> Conjugate[resolve[s]],
      HoldPattern[Coupling[s_Symbol, ___]] :> resolve[s]
    };
    evaluated = evaluated /. (s_Symbol /; SymbolName[s] === "hbar") :>
      payload["hbar_value"];
    badCouplings = DeleteDuplicates @ Cases[evaluated, _Coupling, Infinity];
    badSymbols = DeleteDuplicates @ Cases[evaluated,
      s_Symbol /; Context[s] =!= "System`" &&
        !MemberQ[{True,False}, s] :> ToString[s, InputForm], Infinity];
    If[badCouplings =!= {} || Cases[evaluated, _Missing, Infinity] =!= {},
      Print["ERROR: unresolved parameter in ", channel, " / ", index,
        "; couplings=", InputForm[badCouplings],
        "; missing=", InputForm[Cases[evaluated,_Missing,Infinity]]];
      <|"ok"->False,"channel"->channel,"index"->index|>,
      numeric = Quiet[N[evaluated, MachinePrecision]];
      (* A pre-encoded complex number may arrive as a RawJSON association.
         Decode it BEFORE NumberQ, rather than rejecting valid numerical
         values simply because their head is Association. *)
      If[AssociationQ[numeric] &&
          KeyExistsQ[numeric, "real"] && KeyExistsQ[numeric, "imag"],
        numeric = complexValue[numeric];
      ];
      valid = NumberQ[numeric] && FreeQ[numeric, Indeterminate|ComplexInfinity|DirectedInfinity[___]] &&
        TrueQ[Abs[Re[numeric]] < Infinity] && TrueQ[Abs[Im[numeric]] < Infinity];
      If[!TrueQ[valid],
        Print["ERROR: nonnumeric kernel in ",channel," / ",index,
          "; head=",Head[numeric],
          "; expression preview=",ToString[Short[numeric, 3],InputForm],
          "; symbols=",InputForm[badSymbols]];
        <|"ok"->False,"channel"->channel,"index"->index|>,
        <|"ok"->True,"channel"->channel,"index"->index,
          "real"->N[Re[numeric]],"imag"->N[Im[numeric]]|>
      ]
    ]
  ], {point, payload["points"]}
];
Export[args[[3]], <|"results"->results|>, "RawJSON"];
If[!AllTrue[results, TrueQ[Lookup[#,"ok",False]]&], Exit[4]];
Print["Evaluated ",Length[results]," Matchete C5 kernels."];
Exit[0];
