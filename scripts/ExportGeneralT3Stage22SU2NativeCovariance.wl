(* Live Stage 22 SU(2) gauge-covariance audit for the B_-1 / C_-1
   110101 quartic pair. Read the *registered* RGBeta GroupInvariant functions
   from $quartics, not copied candidate formulae. Native-runtime diagnostic ONLY.

   wolframscript -file scripts/ExportGeneralT3Stage22SU2NativeCovariance.wl \
     output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3Stage22SU2NativeCovariance.wl OUTPUT.json"];
  Exit[64]
];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[FileExistsQ[outPath], DeleteFile[outPath]];
load = CheckAbort[Check[Needs["RGBeta`"]; True, $Failed], $Aborted];
If[load =!= True, Print["RGBeta could not be loaded"]; Exit[3]];
Get[FileNameJoin[{root, "RGE", "running", "rgbeta", "T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root, "scripts", "Stage22RGBetaNativeAccess.wl"}]];
T3S22GaugeFail[msg_] := (Print[msg]; Throw[$Failed, "T3S22Gauge"]);

(* Convert del and eps to inert component heads BEFORE inserting integer
   indices! RGBeta has del[rep,a_,a_] = Dim[rep], a TRACE convention.
   Feeding a numeric value into the original del first would turn each
   diagonal matrix entry into Dim[fund]=2 rather than KroneckerDelta=1.
   The protected heads implement component evaluation independently. *)
T3S22GaugeDelta[i_Integer, j_Integer] := KroneckerDelta[i,j];
T3S22GaugeEpsilon[i_Integer, j_Integer] := Signature[{i,j}];
T3S22GaugeComponent[expr_, inds_List] := Module[{v},
  v = expr /. {
    HoldPattern[del[SU2L[fund], i_, j_]] :> T3S22GaugeDelta[i,j],
    HoldPattern[eps[SU2L[fund], i_, j_]] :> T3S22GaugeEpsilon[i,j]
  };
  v = v /. Thread[{t3i1, t3i2, t3i3, t3i4} -> inds];
  v = FullSimplify[v];
  If[!FreeQ[v, del|eps|delS2|_Missing|Indeterminate|$Failed] ||
     !NumericQ[v] || !TrueQ[FullSimplify[Im[v] == 0]],
    T3S22GaugeFail["110101 expression has an unresolved/nonreal component: "<>
      ToString[v,InputForm]]
  ];
  v
];
T3S22GaugeTensor[expr_] := Table[
  T3S22GaugeComponent[expr,{i1,i2,i3,i4}],
  {i1,1,2},{i2,1,2},{i3,1,2},{i4,1,2}
];

(* Infinitesimal H -> U H; every Bar[field] -> Conjugate[U] Bar[field].
   Terms are contracted T_i G_{i j}, so the Lie derivative is transposed
   on each tensor index as appropriate. Three Pauli generators span su(2). *)
T3S22GaugeResiduals[ten_, n_Integer] := Module[
  {g=PauliMatrix[n]/2, gs, idx},
  gs = {g,-Conjugate[g],-Conjugate[g],-Conjugate[g]};
  Flatten[Table[
    idx={i1,i2,i3,i4};
    FullSimplify[Sum[
      Sum[gs[[axis,b,idx[[axis]]]] * Extract[ten,ReplacePart[idx,axis->b]],
        {b,1,2}], {axis,1,4}]],
    {i1,1,2},{i2,1,2},{i3,1,2},{i4,1,2}]]
];
T3S22GaugeAudit[expr_] := Module[{ten, res, nz},
  ten=T3S22GaugeTensor[expr];
  res=Table[T3S22GaugeResiduals[ten,n],{n,1,3}];
  nz=Count[#,Except[0]]& /@ res;
  <|"su2_covariant" -> (nz === {0,0,0}),
    "nonzero_residual_components_by_generator" -> nz,
    "component_tensor_flat_inputform" -> (ToString[#,InputForm]& /@ Flatten[ten])|>
];

T3S22GaugeRun[] := Module[{models=<||>,dims,build,ext,names,records,
                             exprs,audits,label,i},
  Do[
    label=spec[[1]];dims=spec[[2]];
    Print["Checking native RGBeta SU2 covariance: ",label];
    build=CheckAbort[Check[Apply[T3RGBetaBuild,dims],$Failed],$Aborted];
    If[!AssociationQ[build],T3S22GaugeFail["RGBeta model build failed: "<>label]];
    T3GeneralQuarticAddedSymbols={};
    ext=FileNameJoin[{root,"RGE","general","generated",label<>"_extra_quartics.wl"}];
    If[CheckAbort[Check[Get[ext];True,$Failed],$Aborted]=!=True,
      T3S22GaugeFail["Generated RGBeta extension could not be loaded: "<>label]];
    names={"lambdaT3GeneralQ110101I1","lambdaHHdagS1barS2barCross"};
    records=T3Stage22NativeRecord /@ (Symbol /@ names);
    If[!And@@((StringReplace[Lookup[#,"fields_inputform",""],WhitespaceCharacter..->""]===
        "{H,Bar[H],Bar[S1],Bar[S2]}")& /@ records),
      T3S22GaugeFail["Unexpected native field order: "<>label]];
    exprs=Table[records[[i]]["_native_invariant"][t3i1,t3i2,t3i3,t3i4],
      {i,Length[records]}];
    audits=T3S22GaugeAudit /@ exprs;
    AssociateTo[models,label-><|
      "native_registered_symbols"->names,
      "native_expression_inputform"->(ToString[#,InputForm]& /@ exprs),
      "native_covariance_audits"->audits
    |>],
    {spec,{{"B_-1",{2,2,1,-1}},{"C_-1",{2,2,3,-1}}}}
  ];
  <|
    "format"->"stage22-live-rgbeta-su2-covariance-v1",
    "status"->"NativeWolframSourceCovarianceDiagnosticNotPhysicalCouplingCertificate",
    "field_convention"->"H:U, three Bar fields:Conjugate[U]",
    "models"->models,
    "potential_prefactors_verified"->False,
    "physical_coupling_mapping_certified"->False,
    "production_modified"->False
  |>
];
result=CheckAbort[Check[Catch[T3S22GaugeRun[],"T3S22Gauge"],$Failed],$Aborted];
If[!AssociationQ[result] || Length[Keys[Lookup[result,"models",<||>]]] =!= 2,
  Print["Stage22 native SU2 covariance audit FAILED"];Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[outPath,result,"RawJSON"],$Failed];
If[write===$Failed || !FileExistsQ[outPath],
  Print["Stage22 native SU2 covariance JSON write FAILED"];Exit[6]];
Print["Native SU2 covariance evidence: ",outPath];
Print["This is NOT a physical Matchete-RGBeta coupling certificate"];
