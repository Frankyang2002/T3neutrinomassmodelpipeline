(* Worker C, D(-2): RESTRICTED saved one-loop C5 diagnostic.
   READ ONLY inputs. No new Match[], no production integration.
   Terms with independent y2Conj are set to zero at the saved Wilson
   coefficient level. This is NOT a replacement for matching the
   restricted UV Lagrangian and NOT an independent scalar correction.

   Usage (repo root):
   wolframscript -file Lagrangian/RunWorkerCRestrictedDMinus2Finite.wl \
     output/general_sequential_pilot/D_m2_stage1 \
     output/audits/worker_c_restricted_D_m2_UNIQUE
*)
ClearAll["Global`*"];
args = Rest[$ScriptCommandLine];
If[Length[args] =!= 2, Print["Usage: FIRST_THRESHOLD_ROOT FRESH_OUTPUT_DIR"]; Exit[2]];
root = ExpandFileName[args[[1]]];
out = ExpandFileName[args[[2]]];
profile = FileNameJoin[{root,"data","general_d_m2_scalar_threshold","profile_tree"}];
If[DirectoryQ[out] || StringStartsQ[out, root <> $PathnameSeparator],
  Print["Output exists or lies inside preserved first-threshold root"]; Exit[3]];
files = <|
 "hard" -> FileNameJoin[{profile,"profile_c5_flavour_candidate_hard_symmetric.wxf"}],
 "inherited" -> FileNameJoin[{profile,"profile_c5_flavour_candidate_inherited_symmetric.wxf"}],
 "running" -> FileNameJoin[{profile,"profile_c5_flavour_candidate_running_as_supplied.wxf"}],
 "total" -> FileNameJoin[{profile,"profile_c5_flavour_candidate_total_provisional.wxf"}],
 "flavour" -> FileNameJoin[{profile,"profile_c5_flavour_candidate_result.json"}],
 "pole" -> FileNameJoin[{profile,"profile_c5_pole_cancellation_result.json"}]
|>;
If[!AllTrue[Values[files], FileExistsQ],
 Print["Missing saved source files: ", Select[Values[files], !FileExistsQ[#]&]]; Exit[4]];
flavour = Quiet@Check[Import[files["flavour"],"RawJSON"],$Failed];
poleAudit = Quiet@Check[Import[files["pole"],"RawJSON"],$Failed];
If[!AssociationQ[flavour] || !AssociationQ[poleAudit] ||
 Lookup[flavour,"Status",""] =!= "C5FlavourCandidateExported" ||
 !TrueQ[Lookup[poleAudit,"ExactPoleCancellationVerifiedForNF1",False]],
 Print["Saved flavour/pole audit prerequisites not satisfied"]; Exit[5]];
If[CheckAbort[UsingFrontEnd[Needs["Matchete`"];True],$Aborted] =!= True,
 Print["Matchete unavailable"]; Exit[6]];
ClearAll[wctext,wcname,wcconj,wcrestrict,wcNF1,wcPoleSplit,wcreduce,wcstatus,wcexport];
wctext[x_] := ToString[x,InputForm];
wcname[s_Symbol] := SymbolName[Unevaluated[s]];
wcname[x_] := wctext[x];
wcconj[x_] := !FreeQ[x, HoldPattern[Matchete`Coupling[n_Symbol,___]] /; wcname[n] === "y2Conj"];
wcrestrict[x_] := x /. (HoldPattern[Matchete`Coupling[n_Symbol,idx_List,ord_]] /; wcname[n] === "y2Conj") :> 0;
wcNF1[x_] := Module[{y},
 y=x /. (s_Symbol /; wcname[s] === "NF") :> 1;
 y=FixedPoint[Function[z,z /. HoldPattern[Sum[body_,{r_Symbol,lo_Integer,hi_Integer}]] :>
     If[hi<lo,0,Total[Table[body /. r->k,{k,lo,hi}]]]],y,12];
 y /. (HoldPattern[Matchete`Index[i_,rep_]] /; wcname[rep] === "NFlavor") :> Matchete`Index[1,rep]
];
(* Only supports a Laurent expression affine in the simple pole 1/epsilon.
   Never silently set epsilon to a numerical value. *)
wcPoleSplit[x_] := Module[{tag,eps,other,power,res,fin},
 eps = DeleteDuplicates@Cases[x,s_Symbol /; MemberQ[{"ϵ","ε","epsilon","Epsilon"},wcname[s]] :> s,Infinity];
 If[Length[eps] > 1,Return[$Failed]];
 If[eps === {},Return[<|"residue"->0,"finite"->x,"poleCount"->0|>]];
 power=Cases[x,Power[s_Symbol,k_Integer] /; MemberQ[eps,s] && k =!= -1 :> k,Infinity];
 If[power =!= {},Return[$Failed]];
 tag=Unique["wcPole"];
 other=x /. (Power[s_Symbol,-1] /; MemberQ[eps,s]) :> tag;
 If[!FreeQ[other,Alternatives@@eps],Return[$Failed]];
 other=Expand[other];
 If[!TrueQ[Coefficient[other,tag,2]===0],Return[$Failed]];
 res=Coefficient[other,tag,1]; fin=other /. tag->0;
 <|"residue"->res,"finite"->fin,"poleCount"->Length[Cases[x,Power[s_Symbol,-1] /; MemberQ[eps,s],Infinity]]|>
];
wcreduce[x_] := Module[{v=Expand[wcNF1[x]]}, If[v===0,0,Quiet@Cancel[Together[v]]]];
inputs=Map[Quiet@Check[Import[#,"WXF"],$Failed]&,KeyTake[files,{"hard","inherited","running","total"}]];
If[MemberQ[Values[inputs],$Failed|$Aborted],Print["WXF import failed"];Exit[7]];
If[!wcconj[inputs["total"]] || !wcconj[inputs["hard"]],
 Print["Saved full source lacks expected y2Conj; mismatched profile"]; Exit[8]];
restricted=Map[wcrestrict,inputs];
If[AnyTrue[Values[restricted],wcconj],Print["y2Conj remained after restriction"];Exit[9]];
(* Restrict all four components independently. Do NOT add running twice.
   The saved total is already hard + inherited + running. *)
parts=Map[wcPoleSplit,restricted];
If[AnyTrue[Values[parts],#===$Failed&],Print["Unsupported regulator structure"];Exit[10]];
residue=TimeConstrained[CheckAbort[wcreduce[parts["total"]["residue"]],$Aborted],120,$Aborted];
If[residue===$Aborted,Print["Restricted pole reduction timed out"];Exit[11]];
If[!TrueQ[residue===0],
 CreateDirectory[out,CreateIntermediateDirectories->True];
 Export[FileNameJoin[{out,"result.json"}], <|
  "status"->"RestrictedNF1PoleNotCancelled", "model"->"D_-2", "NF"->1,
  "restricted_total_pole_zero_exact_NF1"->False,
  "pole_residue_preview"->StringTake[wctext[residue],UpTo[1200]],
  "independent_nonbubble_remainder_certified"->False,
  "safe_to_add_to_worker_c_joint"->False|>,"RawJSON"];
 Print["Restricted assembled candidate has nonzero pole. Not exporting finite C5."];Exit[12]];
(* Finite coefficient is a candidate in native Matchete hbar convention.
   Never subtract the numerical scalar bubble here; their normalization
   and overlap remain unverified. *)
CreateDirectory[out,CreateIntermediateDirectories->True];
finite=parts["total"]["finite"];
Export[FileNameJoin[{out,"restricted_total_finite_candidate.wxf"}],finite,"WXF"];
Export[FileNameJoin[{out,"restricted_total_finite_candidate.txt"}],wctext[finite],"Text"];
Scan[Function[k,
  Export[FileNameJoin[{out,"restricted_"<>k<>"_finite_candidate.wxf"}],parts[k]["finite"],"WXF"]],
  {"hard","inherited","running"}];
sourceSha = Map[FileHash[#,"SHA256","HexString"]&,files];
report=<|
 "status"->"RestrictedNF1PoleFreeCandidateExported",
 "model"->"D_-2", "NF"->1,
 "no_new_matching"->True,
 "restricted_saved_matching_not_recomputed"->True,
 "independent_y2Conj_set_zero_in_saved_outputs"->True,
 "restricted_total_pole_zero_exact_NF1"->True,
 "source_sha256"->sourceSha,
 "hard_pole_occurrences"->parts["hard"]["poleCount"],
 "total_pole_occurrences"->parts["total"]["poleCount"],
 "full_restricted_theory_matched"->False,
 "finite_renormalization_certified"->False,
 "weinberg_operator_normalization_certified"->False,
 "physical_mass_relations_imposed"->False,
 "independent_nonbubble_remainder_certified"->False,
 "safe_to_add_to_worker_c_joint"->False,
 "note"->"Finite part of existing assembled symbolic candidate, NOT a new nonbubble threshold correction."
|>;
Export[FileNameJoin[{out,"result.json"}],report,"RawJSON"];
Print["Worker C: pole-free restricted NF=1 candidate exported to ",out];
Print["Independent nonbubble remainder safe to add: False"];
Exit[0];
