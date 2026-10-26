(* Stage 25: Live RGBeta quartic ASSEMBLY check, not a physical sign certificate.
   Inspect the runtime's CouplingTensorMap, scalar field index map, and
   complex-scalar contraction. Use 16-model Stage22 registry.

   wolframscript -file scripts/ExportGeneralT3Stage25NativeQuarticAssembly.wl \
     output/general_t3_stage22/stage25_live_native_quartic_assembly.json
*)
ClearAll["Global`*"];
If[Length[$ScriptCommandLine] =!= 2,
  Print["Usage: ExportGeneralT3Stage25NativeQuarticAssembly.wl OUTPUT.json"]; Exit[64]];
outPath = ExpandFileName[$ScriptCommandLine[[2]]];
root = DirectoryName @ DirectoryName @ ExpandFileName[$InputFileName];
If[FileExistsQ[outPath],DeleteFile[outPath]];
registryPath = FileNameJoin[{root,"RGE","general","generated","general_rgbeta_quartic_registry.json"}];
registry = Quiet@Check[Import[registryPath,"RawJSON"],$Failed];
If[!AssociationQ[registry] || Length[Keys[Lookup[registry,"models",<||>]]] =!= 16,
  Print["Stage25: Incomplete model registry"]; Exit[2]];
load = CheckAbort[Check[Needs["RGBeta`"]; True,$Failed],$Aborted];
If[load =!= True,Print["Stage25: RGBeta unavailable"];Exit[3]];
Get[FileNameJoin[{root,"RGE","running","rgbeta","T3RGBetaModel.wl"}]];
Get[FileNameJoin[{root,"scripts","Stage22RGBetaNativeAccess.wl"}]];
T3S25Fail[msg_String] := (Print["Stage25 ERROR: ",msg];Throw[$Failed,"T3Stage25"]);
(* RGBeta internal identifiers do not necessarily live directly in RGBeta`.
   In particular, CouplingTensorMap has no PackageScope declaration in
   FieldsAndCouplings.m, so it may reside in RGBeta`Private`.  Names["RGBeta`*X"]
   is NOT a reliable search for the symbol X in an arbitrary subcontext.
   Search fully qualified symbol names, and require exactly one matching
   RGBeta-owned symbol.  Do not create symbols speculatively. *)
T3S25Resolve[suffix_String] := Module[{patterns, candidates},
  patterns = {
    "RGBeta`" <> suffix,
    "RGBeta`Private`" <> suffix,
    "RGBeta`*`" <> suffix,
    "*`" <> suffix
  };
  candidates = DeleteDuplicates @ Flatten[Names /@ patterns];
  candidates = Select[candidates,
    StringStartsQ[#, "RGBeta`"] && Last[StringSplit[#, "`"]] === suffix &
  ];
  If[Length[candidates] =!= 1,
    T3S25Fail["Expected one RGBeta runtime symbol named " <> suffix <>
      "; found " <> ToString[candidates, InputForm]]
  ];
  Symbol[First[candidates]]
];
(* These are genuine RGBeta implementation functions/data, not the Stage23
   Python reconstruction of the 1/N permutation weights.  Resolve them
   inside the catch-wrapped model run to avoid uncaught Throw on failure. *)
(* Stage25 compares the output of CouplingTensorMap to an independent
   reconstruction from the same registered invariant.  A nontrivial E-series
   invariant calls T3RGBetaS2TripleInvariant, which has Module[{i,j,k,l,m,n},...].
   Every independent evaluation gives new names i$NNN, j$NNN, etc.

   The live evaluated E_+0 RHS fails both previous occurrence-count tests.
   The prior Matchete/RGBeta evidence records the uncontracted invariant,
   which does satisfy term-wise pairing. Without the live failed RHS, the
   intervening simplification cannot be established. In any event, post-
   evaluation index counts are not a valid precondition for testing whether
   the native permutation assembly equals its independent reconstruction.

   For THIS ASSEMBLY-ONLY comparison, alpha-rename only the six Module-local
   names verified from T3RGBetaS2TripleInvariant in T3RGBetaModel.wl. Require
   a one-to-one base-name mapping within each expression; never rename an
   arbitrary generated symbol or any free t25i index. After renaming, compare
   the WHOLE symbolic RHS (including all numeric factors), not just indices.
   This does not certify SU(2) invariance, free-index structure, or the physical
   action. Those are distinct fail-closed checks elsewhere in the pipeline.
*)
T3S25ModuleIndexBases={"i","j","k","l","m","n"};
T3S25GeneratedLocalQ[x_Symbol] := StringMatchQ[
  SymbolName[x], RegularExpression["[A-Za-z][A-Za-z0-9]*\\$[0-9]+"]
];
T3S25AlphaNormalize[expr_] := Module[
  {generated,bases,replacements},
  generated=DeleteDuplicates@Cases[expr,
    x_Symbol /; T3S25GeneratedLocalQ[x],{0,Infinity}];
  If[generated==={},Return[expr]];
  bases=StringReplace[SymbolName[#],RegularExpression["\\$[0-9]+$"]->""]& /@ generated;
  If[!AllTrue[bases,MemberQ[T3S25ModuleIndexBases,#]&],
    T3S25Fail["Unrecognised generated local symbol in native comparison: "<>
      ToString[Thread[generated->bases],InputForm]]
  ];
  If[!DuplicateFreeQ[bases],
    T3S25Fail["More than one generated symbol shares a Module index base: "<>
      ToString[Thread[generated->bases],InputForm]]
  ];
  replacements=Thread[generated->(T3S25BoundDummy /@ bases)];
  expr /. replacements
];
T3S25SameRule[actual_,expected_] := Module[
  {lhsA,lhsB,rhsA,rhsB,canonA,canonB},
  If[Head[actual]=!=Rule || Head[expected]=!=Rule,Return[False]];
  lhsA=actual[[1]];lhsB=expected[[1]];
  If[!TrueQ[lhsA===lhsB],Return[False]];
  rhsA=actual[[2]];rhsB=expected[[2]];
  If[TrueQ[rhsA===rhsB],Return[True]];
  canonA=T3S25AlphaNormalize[rhsA];
  canonB=T3S25AlphaNormalize[rhsB];
  TrueQ[canonA===canonB] || TrueQ[Simplify[canonA-canonB]===0]
];
(* Self-test covers the exact source-defined E-series triple invariant,
   plus a generic two-branch Module expression.  Numerical coefficients,
   tensor-field labels and external indices remain part of the comparison.
   Unexpected generated symbols remain fail-closed inside AlphaNormalize. *)
T3S25AlphaSelfTest[] := Module[{u,v,e1,e2},
  u=Module[{i,j},
    T3S25TestF[i,j] T3S25TestG[i,j]+T3S25TestF[i,j] T3S25TestH[i,j]];
  v=Module[{i,j},
    T3S25TestF[i,j] T3S25TestG[i,j]+T3S25TestF[i,j] T3S25TestH[i,j]];
  e1=delS2[SU2L,a,t25i1,t25i2]*
    T3RGBetaS2TripleInvariant[a,t25i3,t25i4];
  e2=delS2[SU2L,a,t25i1,t25i2]*
    T3RGBetaS2TripleInvariant[a,t25i3,t25i4];
  If[TrueQ[u===v] ||
     !TrueQ[T3S25SameRule[{1,2,3,4}->u,{1,2,3,4}->v]] ||
     !TrueQ[T3S25SameRule[{1,2,3,4}->e1,{1,2,3,4}->e2]] ||
     TrueQ[T3S25SameRule[{1,2,3,4}->e1,{1,2,3,4}->(2e2)]] ||
     TrueQ[T3S25SameRule[{1,2,3,4}->e1,{1,2,3,4}->
       (e2 /. t25i4->t25i5)]] ||
     TrueQ[T3S25SameRule[{1,2,3,4}->e1,{1,2,3,5}->e2]],
    T3S25Fail["E-series Module-index equivalence self-test failed"]
  ];
  True
];
T3S25Inspect[sym_Symbol] := Module[
  {info,fields,inv,coup,perms,sc,idx,raw,expected,ns,
   pos,subs,coeff,complexMetricOK=True,base,phi,barPhi,indexMap,scalarContraction,
   comparisons,badIndex},
  indexMap=T3S25Resolve["$fieldIndexMap"];
  scalarContraction=T3S25Resolve["$scalarContraction"];
  info=T3Stage22NativeInfo[sym];
  fields=T3Stage22NativeGet[info,"Fields"];
  inv=T3Stage22NativeGet[info,"Invariant"];
  coup=T3Stage22NativeGet[info,"Coupling"];
  perms=T3Stage22NativeGet[info,"UniqueArrangements"];
  sc=T3Stage22NativeGet[info,"SelfConjugate"];
  If[sc =!= False || Length[fields] =!= 4 || !ListQ[perms] ||
     Length[perms]<1 || Length[perms]>24,
    T3S25Fail["Invalid native quartic field/permutation metadata for "<>SymbolName[sym]]];
  If[!AllTrue[perms,Sort[#]==={1,2,3,4}&] || Length[DeleteDuplicates[perms]]=!=Length[perms],
    T3S25Fail["Invalid native permutation list for "<>SymbolName[sym]]];
  idx={t25i1,t25i2,t25i3,t25i4};
  pos=indexMap["All"] /@ fields;
  If[!AllTrue[pos,IntegerQ] || Length[pos]=!=4,
    T3S25Fail["Field index map incomplete for "<>SymbolName[sym]]];
  ns=Length[perms];
  base=(inv@@idx)*(coup@@(T3Stage22NativeGet[info,"Indices"]@@idx));
  raw=With[{f=nativeMap},f[info,fields,idx,idx,perms]];
  If[!ListQ[raw] || Length[raw]=!=ns,
    T3S25Fail["Runtime CouplingTensorMap count differs for "<>SymbolName[sym]]];
  expected=Table[
    subs=Thread[idx->idx[[p]]];
    (pos[[Ordering[p]]] -> (base/ns /.subs)),
    {p,perms}];
  comparisons=MapThread[T3S25SameRule,{raw,expected}];
  If[!And@@comparisons,
    badIndex=First@FirstPosition[comparisons,False];
    Print["Stage25 first mismatching permutation: ",perms[[badIndex]]];
    Print["Stage25 native rule: ",ToString[raw[[badIndex]],InputForm]];
    Print["Stage25 reconstructed rule: ",ToString[expected[[badIndex]],InputForm]];
    T3S25Fail["Native quartic assembly differs from independently reconstructed tensor rule for "<>SymbolName[sym]]
  ];
  (* Index-map entries represent phi and Bar[phi], not Cartesian real/imag
     fields. Verify this convention directly from the runtime metric. *)
  Do[
    phi=If[Head[f]===Bar,f[[1]],f];barPhi=Bar[phi];
    If[!KeyExistsQ[indexMap["Scalars"],phi] || !KeyExistsQ[indexMap["Scalars"],barPhi],
      T3S25Fail["Missing complex conjugate pair in native scalar index map"]];
    If[!TrueQ[scalarContraction[[indexMap["Scalars",phi],indexMap["Scalars",barPhi]]]==2],
      complexMetricOK=False],
    {f,fields}];
  If[!complexMetricOK,T3S25Fail["Native complex-scalar metric is not factor 2"]];
  (* RGBeta's GroupsAndIndices.m defines Matrix[m_][] = m.  The stored
     Coupling pure function is 6*Matrix[sym][##1]&, but with zero flavour
     indices its evaluated VALUE is 6*sym, not 6*Matrix[sym][].
     Compare exact Wolfram expressions rather than a pre-evaluation string. *)
  coeff=coup@@(T3Stage22NativeGet[info,"Indices"]@@idx);
  If[!TrueQ[coeff===6 sym],
    T3S25Fail["Unexpected evaluated native coupling factor in "<>
      SymbolName[sym]<>": "<>ToString[coeff,InputForm]]];
  <|"symbol"->SymbolName[sym],"field_species_positions"->pos,
    "native_permutations"->ns,"native_assembly_rules_checked"->Length[raw],
    "complex_conjugate_metric_two_verified"->complexMetricOK,
    "native_coupling_value_inputform"->ToString[coeff,InputForm],
    "assembly_equals_independent_permutation_reconstruction"->True|>
];
(* Read the native declaration inventory for exact B/C reused symbols and
   all other 50 slots; its model and symbol coverage is checked fail-closed. *)
declarationsPath=FileNameJoin[{root,"output","general_t3_stage22","rgbeta_native_declarations.json"}];
declarations=Quiet@Check[Import[declarationsPath,"RawJSON"],$Failed];
If[!AssociationQ[declarations] || Length[Keys[Lookup[declarations,"models",<||>]]]=!=16,
  Print["Stage25: Need refreshed Stage22 RGBeta declarations"];Exit[4]];
T3S25RunFull[] := Module[{rows=<||>,model,meta,build,ext,names,modelRows,records,total=0,k,record,name},
  T3S25AlphaSelfTest[];
  nativeMap = T3S25Resolve["CouplingTensorMap"];
  Print["Stage25: resolved runtime quartic assembly function: ",
    ToString[nativeMap, InputForm]];
  Do[
    Print["Stage25 native quartic assembly: ",model];
    meta=Join[registry["models",model,"dimensions"],{registry["models",model,"alpha"]}];
    build=CheckAbort[Check[Apply[T3RGBetaBuild,meta],$Failed],$Aborted];
    If[!AssociationQ[build],T3S25Fail["Failed to build "<>model]];
    T3GeneralQuarticAddedSymbols={};
    ext=FileNameJoin[{root,"RGE","general","generated",model<>"_extra_quartics.wl"}];
    If[CheckAbort[Check[Get[ext];True,$Failed],$Aborted]=!=True,
      T3S25Fail["Quartic extension failed for "<>model]];
    names=Lookup[Lookup[declarations["models",model],"declarations",{}],"symbol",{}];
    If[!ListQ[names] || Length[names]=!=registry["models",model,"total_additional_quartics"] ||
       Length[DeleteDuplicates[names]]=!=Length[names],
      T3S25Fail["Invalid slot inventory for "<>model]];
    records=Lookup[declarations["models",model],"declarations",{}];
    Do[
      name=names[[k]];record=records[[k]];
      With[{current=T3Stage22NativeRecord[Symbol[name]]},
        If[!And@@(StringReplace[Lookup[record,#],WhitespaceCharacter..->""]===
                     StringReplace[Lookup[current,#],WhitespaceCharacter..->""]& /@
                     {"fields_inputform","invariant_inputform","coupling_factor_inputform","unique_arrangements_inputform"}),
          T3S25Fail["Stale Stage22 declaration: "<>model<>"/"<>name]]],
      {k,Length[names]}];
    modelRows=T3S25Inspect[Symbol[#]]& /@ names;
    AssociateTo[rows,model->modelRows];total+=Length[modelRows],
    {model,Sort[Keys[declarations["models"]]]}];
  If[total=!=50,T3S25Fail["Expected 50 native quartics, got "<>ToString[total]]];
  <|"format"->"stage25-live-rgbeta-quartic-assembly-v1",
    "status"->"NativeComplexQuarticAssemblyCheckedNotPhysicalPotentialSignCertificate",
    "source"->"RGBeta runtime CouplingTensorMap and scalar index/contraction metric",
    "models"->rows,"slots_checked"->total,"models_checked"->16,
    "all_native_permutation_reconstructions_exact"->True,
    "all_complex_conjugate_metric_two_verified"->True,
    "native_complex_species_to_cartesian_real_intertwiner_verified"->False,
    "matchete_relative_potential_sign_verified"->False,
    "physical_coupling_mapping_certified"->False,"production_modified"->False|>
];
result=CheckAbort[Check[Catch[T3S25RunFull[],"T3Stage25"],$Failed],$Aborted];
If[!AssociationQ[result] || Lookup[result,"slots_checked",0]=!=50,
  Print["Stage25 native quartic assembly FAILED"];Exit[5]];
If[!DirectoryQ[DirectoryName[outPath]],
  CreateDirectory[DirectoryName[outPath],CreateIntermediateDirectories->True]];
write=Quiet@Check[Export[outPath,result,"RawJSON"],$Failed];
If[write===$Failed || !FileExistsQ[outPath],
  Print["Stage25 native JSON export FAILED"];Exit[6]];
Print["Stage25 native complex-field assembly evidence: ",outPath];
Print["Physical sign and Cartesian real-field basis: UNCERTIFIED"];
