(* Shared, fail-closed access to RGBeta's $quartics native internal records.
   RGBeta's Association keys are Symbol objects in a package context,
   not necessarily Global`Fields / Global`Invariant / etc. Never use
   info[Fields] from a separate Wolfram context.

   This file is loaded by the two Stage22 native evidence exporters.
*)
ClearAll[
  T3Stage22NativeError,
  T3Stage22NativeKeyName,
  T3Stage22NativeGet,
  T3Stage22NativeInfo,
  T3Stage22NativeString,
  T3Stage22NativeRecord
];

T3Stage22NativeError[msg_String] := (
  Print["ERROR Stage22 RGBeta native access: ", msg];
  Throw[$Failed, "T3Stage22Native"]
);

(* For example, RGBeta`Private`Fields -> "Fields". *)
T3Stage22NativeKeyName[key_] := Module[{s},
  If[StringQ[key], Return[key]];
  s = ToString[key, InputForm];
  Last[StringSplit[s, "`"]]
];

T3Stage22NativeGet[info_Association, name_String] := Module[
  {matches, value, available},
  available = T3Stage22NativeKeyName /@ Keys[info];
  matches = Select[Keys[info], T3Stage22NativeKeyName[#] === name &];
  If[Length[matches] =!= 1,
    T3Stage22NativeError[
      "Expected exactly one native key '" <> name <>
      "'; found " <> ToString[Length[matches]] <>
      "; native keys: " <> StringRiffle[available, ", "]
    ]
  ];
  value = Lookup[info, First[matches], Missing["Stage22MissingNativeKey", name]];
  If[MissingQ[value],
    T3Stage22NativeError["Native value is Missing for key " <> name]
  ];
  value
];

T3Stage22NativeInfo[symbol_Symbol] := Module[{info},
  If[!KeyExistsQ[$quartics, symbol],
    T3Stage22NativeError["Unregistered RGBeta quartic " <> SymbolName[symbol]]
  ];
  info = $quartics[symbol];
  If[!AssociationQ[info],
    T3Stage22NativeError["Quartic registration is not an Association for " <> SymbolName[symbol]]
  ];
  Do[T3Stage22NativeGet[info, name],
    {name, {"Fields", "Invariant", "Coupling", "UniqueArrangements", "SelfConjugate"}}
  ];
  info
];

T3Stage22NativeString[value_] := Module[{s},
  If[MissingQ[value] || value === $Failed || value === $Aborted,
    T3Stage22NativeError["Attempted to serialize absent native value"]
  ];
  s = ToString[value, InputForm];
  If[StringContainsQ[s, "Missing["] || StringContainsQ[s, "$Failed"] ||
     StringContainsQ[s, "$Aborted"],
    T3Stage22NativeError["Serialized native value contains an unresolved failure: " <> s]
  ];
  s
];

T3Stage22NativeRecord[symbol_Symbol] := Module[{info, f, inv, coup, perms, sc},
  info = T3Stage22NativeInfo[symbol];
  f = T3Stage22NativeGet[info, "Fields"];
  inv = T3Stage22NativeGet[info, "Invariant"];
  coup = T3Stage22NativeGet[info, "Coupling"];
  perms = T3Stage22NativeGet[info, "UniqueArrangements"];
  sc = T3Stage22NativeGet[info, "SelfConjugate"];
  If[sc =!= False,
    T3Stage22NativeError["Unexpected self-conjugate extra quartic: " <> SymbolName[symbol]]
  ];
  If[!ListQ[f] || Length[f] =!= 4,
    T3Stage22NativeError["Expected four native fields for " <> SymbolName[symbol]]
  ];
  <|
    "fields_inputform" -> T3Stage22NativeString[f],
    "invariant_inputform" -> T3Stage22NativeString[inv],
    "coupling_factor_inputform" -> T3Stage22NativeString[coup],
    "unique_arrangements_inputform" -> T3Stage22NativeString[perms],
    "self_conjugate" -> False,
    "_native_invariant" -> inv
  |>
];
