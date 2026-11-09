(* A v6: live native RGBeta SU(2) fundamental Fierz / S2 structural witness.
   This is a symbolic identity audit. tGen[SU2L[S2]] remains ABSTRACT.
   wolframscript -file scripts/RunT3ParallelAExportNativeFierz.wl OUTPUT.json
   No production code, beta functions or previous outputs are changed.
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: RunT3ParallelAExportNativeFierz.wl OUTPUT.json"]; Exit[64]
];
out = ExpandFileName[$ScriptCommandLine[[2]]];
If[FileExistsQ[out], Print["Refusing to overwrite ", out]; Exit[6]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
sources = {"RGE/running/rgbeta/T3RGBetaModel.wl",
  "scripts/RunT3ParallelAExportNativeFierz.wl"};
If[!AllTrue[sources, FileExistsQ[FileNameJoin[{root,#}]] &],
  Print["A v6 missing RGBeta/model source"]; Exit[2]
];
loaded = CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted];
If[loaded =!= True, Print["A v6 RGBeta unavailable"]; Exit[3]];
Get[FileNameJoin[{root, "RGE", "running", "rgbeta", "T3RGBetaModel.wl"}]];
T3A6Fail[message_String] := (
  Print["A v6 native Fierz BLOCKED: ", message];
  Throw[$Failed, "T3A6Fierz"]
);
T3A6Probe[] := Module[{build, fundamental, triplet, a, b, c, d, ga,
    lhs, rhs, bar, tripletGen, structure, version},
  build = CheckAbort[T3RGBetaBuildEFT1[2,2,1,-1], $Aborted];
  If[!AssociationQ[build],
    T3A6Fail["Native EFT1 base-model registration failed"]
  ];
  fundamental = SU2L[fund];
  triplet = SU2L[S2];
  If[!TrueQ[Dim[fundamental] === 2] ||
     !TrueQ[TraceNormalization[fundamental] === 1/2] ||
     !TrueQ[Casimir2[fundamental] === 3/4],
    T3A6Fail["Native SU2 fundamental representation normalization changed"]
  ];
  If[!TrueQ[Dim[triplet] === 3] ||
     !TrueQ[TraceNormalization[triplet] === 2] ||
     !TrueQ[Casimir2[triplet] === 2],
    T3A6Fail["Native SU2 S2 normalization changed"]
  ];
  (* Keep every free index DISTINCT. In RGBeta del[R,i,i] denotes an
     Einstein TRACE and becomes Dim[R], not a numerical component! *)
  {a,b,c,d,ga} = {t6a,t6b,t6c,t6d,t6A};
  lhs = tGen[fundamental,ga,a,b] tGen[fundamental,ga,c,d];
  rhs = (del[fundamental,a,d] del[fundamental,c,b] -
    del[fundamental,a,b] del[fundamental,c,d]/2)/2;
  If[!TrueQ[FullSimplify[lhs-rhs] === 0],
    Print["Native Fierz lhs: ", ToString[lhs,InputForm]];
    Print["Expected rhs: ", ToString[rhs,InputForm]];
    T3A6Fail["SU2 fundamental Fierz identity differs"]
  ];
  bar = tGen[Bar[triplet],ga,a,b];
  If[!TrueQ[bar === -tGen[triplet,ga,b,a]],
    T3A6Fail["S2 conjugate-generator minus-transpose rule changed"]
  ];
  tripletGen = tGen[triplet,ga,a,b];
  If[Head[tripletGen] =!= tGen,
    T3A6Fail["Unexpected numeric S2 generator, needs independent audit"]
  ];
  structure = delS2[SU2L,t6r,t6i,t6j] delS2[SU2L,t6s,t6i,t6j];
  If[!TrueQ[structure === del[triplet,t6r,t6s]],
    T3A6Fail["S2 symmetric-doublet isometry changed"]
  ];
  <|
    "format" -> "t3-parallel-a-native-rgbeta-fierz-witness-v6",
    "status" -> "NativeSU2FundamentalFierzExact_S2GeneratorsStillAbstract",
    "native_fundamental_Fierz_exact" -> True,
    "native_symmetric_S2_rep_dimension_three" -> True,
    "native_S2_trace_normalization_two" -> True,
    "native_triplet_generators_symbolic" -> True,
    "native_barred_generator_minus_transpose" -> True,
    "native_delS2_isometry" -> True,
    "native_Fierz_definition" ->
      "tGen[fund,A,i,j] tGen[fund,A,k,l] = (delta[i,l]delta[k,j]-delta[i,j]delta[k,l]/2)/2",
    "native_Fierz_evaluated_inputform" -> ToString[lhs,InputForm],
    "native_S2_generator_inputform" -> ToString[tripletGen,InputForm],
    "native_S2_component_intertwiner_certified" -> False,
    "physical_portal_coupling_dictionary_certified" -> False
  |>
];
report = CheckAbort[Check[Catch[T3A6Probe[], "T3A6Fierz"], $Failed], $Aborted];
If[!AssociationQ[report],
  Print["A v6 native Fierz check FAILED; no JSON written"]; Exit[5]
];
sourceHash = Association@Table[
  path -> IntegerString[FileHash[FileNameJoin[{root,path}], "SHA256"],16,64],
  {path,sources}
];
report = Join[report, <|"source_sha256" -> sourceHash|>];
If[!DirectoryQ[DirectoryName[out]],
  CreateDirectory[DirectoryName[out],CreateIntermediateDirectories -> True]
];
write = Quiet @ Check[Export[out,report,"RawJSON"],$Failed];
If[write === $Failed || !FileExistsQ[out],
  Print["A v6 Fierz report write FAILED"]; Exit[7]
];
Print["A v6 native RGBeta SU2 fundamental Fierz: EXACT; S2 remains ABSTRACT -> ",out];
