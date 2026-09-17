path = "output/full/dimensions/T3_A_alpha_p0/debug/fresh_kernel_stage_2/input_cg_registry.wxf";

If[!FileExistsQ[path],
  Print["FILE NOT FOUND: ", AbsoluteFileName[path]];
  Exit[1];
];

Print["Reading: ", AbsoluteFileName[path]];

x = Quiet @ Check[
  Import[path],
  $Failed
];

If[x === $Failed,
  Print["Plain Import failed; trying BinaryDeserialize[Import[..., \"ByteArray\"]]."];
  x = Quiet @ Check[
    BinaryDeserialize[Import[path, "ByteArray"]],
    $Failed
  ];
];

If[x === $Failed,
  Print["FAILED: could not deserialize registry."];
  Exit[2];
];

Print["Head: ", Head[x]];

names = DeleteDuplicates @ Cases[
  x,
  s_String /; StringContainsQ[s, "T3" ~~ __ ~~ "CG"] :> s,
  Infinity
];

Print["Detected CG names: ", InputForm[names]];
Print["Contains T3Y1CG: ", MemberQ[names, "T3Y1CG"]];
Print["Contains T3Y2CG: ", MemberQ[names, "T3Y2CG"]];
Print["Contains T3MixCG: ", MemberQ[names, "T3MixCG"]];

Print["\nFull InputForm:"];
Print[InputForm[x]];
