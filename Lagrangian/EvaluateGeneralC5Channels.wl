(* Evaluate all quadratic Yukawa monomials in an exported Matchete C5.
   Does not load Matchete. Usage: coefficient.txt parameter_points.json output.json *)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 3, Print["Expected coefficient, points JSON, output JSON"]; Exit[2]];
expr = Quiet[ToExpression[Import[args[[1]], "Text"], InputForm]];
payload = Import[args[[2]], "RawJSON"];
If[!AssociationQ[payload] || !ListQ[Lookup[payload, "points", None]], Exit[3]];
complexValue[x_] := If[AssociationQ[x], N[x["real"] + I x["imag"]], N[x]];
yukawaNames = {"y1", "y2", "y1Conj", "y2Conj"};
(* Use unique dummy symbols for bare and conjugated Yukawa factors. *)
marker[n_String, bar_] := Symbol["Global`T3YUK" <> n <> If[bar, "Bar", "Plain"]];
withY = expr /. {
  HoldPattern[Bar[Coupling[s_Symbol, ___]]] /; MemberQ[yukawaNames, SymbolName[s]] :>
    marker[SymbolName[s], True],
  HoldPattern[Coupling[s_Symbol, ___]] /; MemberQ[yukawaNames, SymbolName[s]] :>
    marker[SymbolName[s], False]
};
vars = DeleteDuplicates @ Cases[withY, x_Symbol /; StringStartsQ[SymbolName[x], "T3YUK"], Infinity];
vars = SortBy[vars, SymbolName];
Print["Yukawa markers: ", InputForm[vars]];
If[Length[vars] < 2, Print["ERROR: fewer than two Yukawa markers"]; Exit[4]];
(* CoefficientRules ensures every Yukawa monomial is kept separately;
   rejecting other degrees prevents silently discarding contributions. *)
rules = CoefficientRules[withY, vars];
If[AnyTrue[rules, Total[First[#]] != 2 &],
 Print["ERROR: C5 contains Yukawa-degree other than two: ",
   InputForm[DeleteDuplicates[Total[First[#]] & /@ rules]]]; Exit[5]
];
channels = Select[rules, Total[First[#]] == 2 &];
If[channels === {}, Print["ERROR: no quadratic Yukawa channels"]; Exit[5]];
channelName[exps_List] := StringRiffle[Flatten[Table[
  ConstantArray[SymbolName[vars[[i]]], exps[[i]]], {i, Length[vars]}]], "__"];
results = Flatten[Table[
  Module[{pt = point, index, p, val, remaining, invalid, exp, ch},
    index = pt["index"];
    p = Association @ KeyValueMap[(#1 -> complexValue[#2]) &, pt["couplings"]];
    Table[
      exp = First[rule]; ch = channelName[exp];
      val = Last[rule] /. {
        HoldPattern[Bar[Coupling[s_Symbol, ___]]] :>
          If[KeyExistsQ[p, SymbolName[s]], Conjugate[p[SymbolName[s]]],
             Missing["Unmapped", SymbolName[s]]],
        HoldPattern[Coupling[s_Symbol, ___]] :>
          If[KeyExistsQ[p, SymbolName[s]], p[SymbolName[s]],
             Missing["Unmapped", SymbolName[s]]]
      };
      (* Materialize the JSON scalar before constructing the replacement rule.
         The old RuleDelayed retained payload[hbar_value] as unevaluated code. *)
      With[{loopFactorValue = N[Lookup[payload, "hbar_value", 1.0]]},
        val = val /. (s_Symbol /; SymbolName[s] === "hbar") :> loopFactorValue
      ];
      (* Coefficient data and coupling values must be scalar, not JSON associations. *)
      If[AssociationQ[val] && KeyExistsQ[val, "real"] && KeyExistsQ[val, "imag"],
        val = complexValue[val]
      ];
      remaining = DeleteDuplicates @ Cases[val, _Missing | _Coupling, Infinity];
      invalid = remaining =!= {};
      If[!invalid, val = Quiet[N[val, MachinePrecision]];
        invalid = !NumberQ[val] || !FreeQ[val, Indeterminate | ComplexInfinity | DirectedInfinity[___]]];
      If[invalid,
        Print["ERROR index=", index, " channel=", ch, " remaining=", InputForm[remaining],
           " preview=", InputForm[Short[val, 2]]];
        <|"ok"->False,"index"->index,"channel"->ch|>,
        <|"ok"->True,"index"->index,"channel"->ch,
          "factors"->Flatten[Table[ConstantArray[SymbolName[vars[[i]]], exp[[i]]],
                         {i,Length[vars]}]],
          "real"->N[Re[val]],"imag"->N[Im[val]]|>
      ], {rule, channels}]
  ], {point, payload["points"]}], 1];
Export[args[[3]], <|"results"->results, "channels"->(channelName[First[#]]& /@ channels)|>, "RawJSON"];
If[!AllTrue[results, TrueQ[Lookup[#, "ok", False]] &], Exit[6]];
Print["Evaluated ", Length[results], " numerical C5 kernels."];
