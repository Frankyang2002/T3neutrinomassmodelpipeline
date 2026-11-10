(* ::Package:: *)

Package["Matchete`"] 


(* ::Title:: *)
(*GroupMagic`*)


(* ::Subtitle:: *)
(*Package for Lie algebra computations *)


(* ::Text:: *)
(*Important references include:*)
(*	R. Cahn, "Semi-simple Lie algebras  and their representations" *)
(*	R. Fonseca, "Renormalization in supersymmetric models," [1310.1296]*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section::Closed:: *)
(*Scoping & usage definitions*)


PackageScope["OptionsCheck"]


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Global definitions*)


(* ::Subsection::Closed:: *)
(*Option checker*)


(* ::Text:: *)
(*General function for checking Options of functions*)


(*SetAttributes[OptionsCheck, HoldFirst];
OptionsCheck @ func_[___, opts : OptionsPattern[]] := And[
	And@@ (Message[General::invalidopt, #1, func] &)@@@ FilterRules[List@ opts, Except@ Options@ func],
	And@@ (OptionTest[func, #1][#2] || OptionMessage[#1, func, #2] &)@@@ FilterRules[List@ opts, Options@ func]
];*)


(*General::invalidopt = "Invalid option `1` given for function `2`.";
General::invalidarg = "Option `1` for function `2` received invalid value `3`.";
General::optexpectsval = "Option `1` for function `2` received invalid value `3`. A `4` is expected.";
OptionMessage[opt_, func_, val_] := Message[General::invalidarg, opt, func, val];*)


(* ::Subsubsection:: *)
(*Messages*)


(*Messages for specific options/functions*)
OptionMessage[CGsToReplace, func_, val_] := Message[General::optexpectsval, CGsToReplace, func, val, "list of specific CGs, All, or None"];
OptionMessage[Levels, func_, val_] := Message[General::optexpectsval, Levels, func, val, "Boolean"];
OptionMessage[Normalization, func_, val_] := Message[General::optexpectsval, Normalization, func, val, "positive, real number or Default"];


(* ::Subsubsection:: *)
(*Tests*)


(*Tests for specific options. Form expected is OptionTest[function, options] *)
OptionTest[_, AntisymmetricIndices]= MatchQ@ _List;
OptionTest[_, Levels]= BooleanQ;
OptionTest[_, Normalization]= Positive@ # || MatchQ[#, Default]&;
OptionTest[_, SymmetricIndices]= MatchQ@ _List;
OptionTest[_, CGsToReplace]= MatchQ@ Alternatives[All, None, _List];
