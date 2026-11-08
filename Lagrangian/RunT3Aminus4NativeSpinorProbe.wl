(* A(-4) live Matchete action and ABSTRACT charge-conjugation probe.
   Uses existing full-general A(-4) UV builder with F Heavy and S1,S2 Light.
   No RGE or scalar matching. It DOES NOT set GammaCC = any 4x4 matrix.

   From the repository root, set the four experimental T3 environment flags
   to 1 (see PowerShell commands in delivery), then:
     wolframscript -file Lagrangian/RunT3Aminus4NativeSpinorProbe.wl OUTPUT.json

   A successful status means the UV action and native abstract identities
   were evaluated in a live Matchete kernel. It does NOT certify the phase,
   the 4x4 charge-conjugation matrix, or the master Weyl basis.
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 1, Print["Usage: OUTPUT.json"]; Exit[2]];
out = ExpandFileName[First[args]];
sourceDir = DirectoryName @ ExpandFileName[$InputFileName];
flags = {"T3_EXPERIMENTAL_GENERAL_FERMIONS",
  "T3_EXPERIMENTAL_GENERAL_QUARTICS",
  "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
  "T3_EXPERIMENTAL_FLAVOR_MASSES"};
If[!And @@ (Environment[#] === "1" & /@ flags),
  Print["All four general F-threshold flags must be 1."]; Exit[2]];
ClearAll[asString, writeResult, fail];
asString[x_] := ToString[InputForm[x]];
writeResult[obj_Association] := Module[{parent = DirectoryName[out]},
  If[!DirectoryQ[parent], CreateDirectory[parent, CreateIntermediateDirectories -> True]];
  If[Quiet@Check[Export[out, obj, "RawJSON"], $Failed] === $Failed,
    Print["Cannot export ", out]; Exit[7]];
];
fail[reason_] := (writeResult[<|
  "format" -> "t3-a-minus4-native-spinor-action-probe-v1",
  "model" -> "A_-4", "status" -> "BLOCKED",
  "reason" -> reason, "physical_certification" -> False|>];
  Print["Native A(-4) probe blocked: ", reason]; Exit[4]);

load = CheckAbort[Check[UsingFrontEnd[Needs["Matchete`"]; True], $Failed], $Aborted];
If[load =!= True, fail["Could not load native Matchete in a fresh kernel"]];
Get[FileNameJoin[{sourceDir,"model","T3ModelCatalog.wl"}]];
Get[FileNameJoin[{sourceDir,"model","LagrangianBuilder.wl"}]];
model = T3ModelFromClass["A", -4, 3];
If[!AssociationQ[model] || model["Fermion", "SU2"] =!= 2 ||
   model["Scalar1", "SU2"] =!= 1 || model["Scalar2", "SU2"] =!= 3,
  fail["Wrong A(-4) SU(2) representations"]];
If[SetT3HeavyFields[{"F"}] === $Failed, fail["F-only heavy configuration failed"]];
built = CheckAbort[UsingFrontEnd[BuildT3Lagrangian[model]], $Aborted];
If[!AssociationQ[built] || Lookup[built,"Status",""] =!= "Success" ||
   !TrueQ[Lookup[built,"FullValidation",False]] ||
   !TrueQ[Lookup[built,"FermionYukawaInventoryAssembled",False]] ||
   !TrueQ[Lookup[built,"AllScalarExtensionsRequested",False]] ||
   Lookup[built,"AdditionalQuadraticsExpected",-1] =!= Lookup[built,"AdditionalQuadraticsIncluded",-2] ||
   Lookup[built,"AdditionalQuarticsExpected",-1] =!= Lookup[built,"AdditionalQuarticsIncluded",-2],
   fail["Native A(-4) full-general UV action was not validated"]];

(* Source-defined topology, preserving the exact CConj and field orientation.
   These are native Lagrangian terms, not prose flags from a saved report. *)
allTerms = Module[{x = Expand[built["LInt"]]}, If[Head[x] === Plus, List @@ x, {x}]];
selectTerms[name_] := Select[allTerms,
   StringContainsQ[asString[#], "Coupling[" <> name <> ","] &];
y1Terms = selectTerms["y1"];
y2Terms = selectTerms["y2"];
heavyFree = asString[FreeLag[NewFermion]];
checksAction = <|
  "native_validated_UV_action" -> TrueQ[Lookup[built,"FullValidation",False]],
  "only_F_heavy" -> ($T3HeavyFields === {"F"}),
  "three_heavy_flavours" -> (model["Fermion","Multiplicity"] === 3),
  "both_native_topology_yukawas_present" ->
    (Length[y1Terms] > 0 && Length[y2Terms] > 0),
  "y1_native_action_contains_charge_conjugation" ->
    AnyTrue[y1Terms, StringContainsQ[asString[#],"GammaCC"] &],
  "heavy_free_action_has_MF" -> StringContainsQ[heavyFree,"MF"],
  "heavy_free_action_has_NewFermion" -> StringContainsQ[heavyFree,"NewFermion"]
|>;

(* Native abstract identities only. These identities cannot fix C vs -C. *)
lf = l[i,p];
hf = NewFermion[j,r];
spin = NCM[Bar[lf], CConj[hf]];
explicit = NCM[Bar[lf], DiracProduct[GammaCC], Transp[Bar[hf]]];
ccSquared = DiracProduct[GammaCC, GammaCC];
ccTranspose = Transp[GammaCC];
chiralL = CConj[NCM[PL,hf]];
chiralR = NCM[PR,CConj[hf]];
checksSpinor = <|
  "native_CC_square_minus_one" -> TrueQ[FullSimplify[ccSquared == -1]],
  "native_CC_transpose_antisymmetric" -> TrueQ[FullSimplify[ccTranspose == -GammaCC]],
  "native_charge_conjugated_bilinear" -> TrueQ[FullSimplify[spin == explicit]],
  "native_chiral_projector_flip" -> TrueQ[FullSimplify[chiralL == chiralR]]
|>;
passed = AllTrue[Values[checksAction],TrueQ] &&
  AllTrue[Values[checksSpinor],TrueQ];
record = <|
  "format" -> "t3-a-minus4-native-spinor-action-probe-v1",
  "model" -> "A_-4",
  "status" -> If[passed,"NativeAminus4UVActionAndAbstractCCVerified_PhaseUnmeasured","BLOCKED"],
  "checks_action" -> checksAction,
  "checks_abstract_spinor" -> checksSpinor,
  "native_heavy_fields" -> $T3HeavyFields,
  "native_light_fields" -> {"S1","S2"},
  "native_yukawa_1_terms_inputform" -> (asString /@ y1Terms),
  "native_yukawa_2_terms_inputform" -> (asString /@ y2Terms),
  "native_heavy_free_inputform" -> heavyFree,
  "native_spinor_bilinear_inputform" -> asString[spin],
  "native_spinor_bilinear_explicit_inputform" -> asString[explicit],
  "native_cc_squared_inputform" -> asString[ccSquared],
  "native_cc_transpose_inputform" -> asString[ccTranspose],
  "native_concrete_4x4_cc_matrix" -> Null,
  "native_4x4_CC_to_Weyl_epsilon_ratio" -> Null,
  "native_connected_UV_amplitude_calculated" -> False,
  "native_source_matched_contact_amplitude_calculated" -> False,
  "full_F_matching_performed_by_this_probe" -> False,
  "physical_certification" -> False,
  "note" -> "Matchete GammaCC is an abstract Dirac-algebra symbol in this probe. C and -C obey the same identities. A signed external-amplitude comparison and master basis mapping remain required."
|>;
writeResult[record];
Print["Native A(-4) spinor-action status: ",record["status"]];
If[!passed, Exit[4]];
Exit[0];
