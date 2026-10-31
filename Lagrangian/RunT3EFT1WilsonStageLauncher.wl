(* Experimental local-only child-kernel launcher. Copied from
   RunMatching.wl / RunFreshThresholdStage, changing ONLY the fresh stage
   script name. Must be loaded AFTER the authoritative RunMatching.wl.
   Does not mutate the production RunFreshThresholdStage symbol. *)
ClearAll[RunT3EFT1WilsonFreshThresholdStage];
RunT3EFT1WilsonFreshThresholdStage[

  currentTree_,

  currentLoop_,

  transitionTree_,

  transitionFull_,

  model_Association,

  group_List,

  active_List,

  level_Integer,

  eftOrder_Integer,

  loopOrder_Integer,

  outputDirectory_String

] := Module[

  {

    stageScript, transitionDir, treePath, loopPath,

    transitionTreePath, transitionFullPath, cgPath,

    resultPath, logPath, command, process, result, dS1, dS2, dF, alpha,

    cgNames, cgRegistry, normalizeCGName

  },



  stageScript = FileNameJoin @ {

    DirectoryName @ ExpandFileName[$InputFileName],

    "RunThresholdStageEFT1Wilson.wl"

  };



  transitionDir = FileNameJoin @ {

    outputDirectory,

    "debug",

    "fresh_kernel_stage_" <> ToString[level]

  };

  Quiet @ CreateDirectory[

    transitionDir,

    CreateIntermediateDirectories -> True

  ];



  treePath = FileNameJoin[{transitionDir, "input_tree.wxf"}];

  loopPath = FileNameJoin[{transitionDir, "input_loop.wxf"}];

  transitionTreePath =

    FileNameJoin[{transitionDir, "input_transition_tree.wxf"}];

  transitionFullPath =

    FileNameJoin[{transitionDir, "input_transition_full.wxf"}];

  cgPath = FileNameJoin[{transitionDir, "input_cg_registry.wxf"}];

  resultPath = FileNameJoin[{transitionDir, "result.wxf"}];

  logPath = FileNameJoin[{transitionDir, "stdout.log"}];



  (* The explicit tree/loop split is retained for reporting and diagnostics.

     The separate transition pair is what the next Match actually consumes:

       transitionTree = canonical EOM EFT at O(hbar^0)

       transitionFull = canonical EOM EFT through O(hbar)

     Crucially, effective couplings remain unexpanded in these transition

     expressions. *)

  Export[treePath, currentTree, "WXF"];

  Export[loopPath, currentLoop, "WXF"];

  Export[transitionTreePath, transitionTree, "WXF"];

  Export[transitionFullPath, transitionFull, "WXF"];



  Print[

    "  Exporting canonical transition EFTs without expanding effective ",

    "couplings."

  ];



  normalizeCGName[name_] := Replace[name, Bar[x_] :> x];



  cgNames = DeleteDuplicates @ Cases[

    HoldComplete[

      currentTree,

      currentLoop,

      transitionTree,

      transitionFull

    ],

    CG[name_, ___] :> normalizeCGName[name],

    Infinity

  ];



  cgRegistry = DeleteCases[

    Map[

      Function[name,

        Module[{reps, tensor},

          reps = Quiet @ Check[

            Matchete`PackageScope`$CGproperties[name, Indices],

            $Failed

          ];

          tensor = Quiet @ Check[

            Matchete`CGManipulations`PackagePrivate`$CGtensors[name],

            $Failed

          ];

          If[

            reps === $Failed || tensor === $Failed ||

              MissingQ[reps] || MissingQ[tensor],

            Nothing,

            <|"Name" -> name, "Reps" -> reps, "Tensor" -> tensor|>

          ]

        ]

      ],

      cgNames

    ],

    Nothing

  ];



  Export[cgPath, cgRegistry, "WXF"];



  Print[

    "  Exported fresh-kernel CG registry: ",

    Length[cgRegistry],

    " definition(s)."

  ];



  dS1 = model["Scalar1", "SU2"];

  dS2 = model["Scalar2", "SU2"];

  dF = model["Fermion", "SU2"];

  alpha = model["Alpha"];



  Print[

    "\n  Launching fresh Wolfram kernel for threshold stage ",

    level,

    " (", StringRiffle[group, ", "], ")..."

  ];



  command = {

    "wolframscript",

    "-file",

    stageScript,

    resultPath,

    ToString[eftOrder],

    ToString[loopOrder],

    ToString[dS1],

    ToString[dS2],

    ToString[dF],

    If[alpha < 0, "m" <> ToString[Abs[alpha]], "p" <> ToString[alpha]],

    StringRiffle[active, ","],

    StringRiffle[group, ","],

    treePath,

    loopPath,

    transitionTreePath,

    transitionFullPath,

    cgPath

  };



  Print[

    "  Fresh-kernel CLI payload count (after -file script): ",

    Length[command] - 3

  ];



  process = RunProcess[command];



  Export[

    logPath,

    Lookup[process, "StandardOutput", ""] <>

      Lookup[process, "StandardError", ""],

    "Text"

  ];



  PrintFreshKernelSummary[

    Lookup[process, "StandardOutput", ""],

    Lookup[process, "StandardError", ""]

  ];



  If[

    Lookup[process, "ExitCode", 1] =!= 0 ||

      !FileExistsQ[resultPath],

    Print[

      "ERROR: fresh-kernel threshold stage ", level,

      " failed. See ", logPath

    ];

    Return[$Failed]

  ];



  result = Quiet @ Check[Import[resultPath, "WXF"], $Failed];



  If[

    !AssociationQ[result] ||

      Lookup[result, "Status", ""] =!= "Success",

    Print[

      "ERROR: fresh-kernel threshold stage ", level,

      " returned an invalid result. See ", logPath

    ];

    Return[$Failed]

  ];



  result

];



