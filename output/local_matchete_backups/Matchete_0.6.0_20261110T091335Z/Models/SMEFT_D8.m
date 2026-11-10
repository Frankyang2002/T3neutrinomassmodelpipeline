(* ::Package:: *)

(* ::Title:: *)
(*SMEFT D=8 model file*)


(* ::Subtitle:: *)
(*Modified version of the D=8 Basis of Murphy http://arxiv.org/abs/2005.00059v6 *)


ParentModel["SMEFT_Warsaw"]


(* ::Subsubsection:: *)
(*Couplings*)


(* Dimension 7 *)
If[TrueQ@ BLViolation,
	(* 1 *)
	DefineCoupling[cl2WH2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)\)]\)"];
	DefineCoupling[cl2BH2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)\)]\)"];
	(* 2 *)
	DefineCoupling[cl2H4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\)]\)"];
	(* 3 *)
	DefineCoupling[cl3eH,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(3\)]\)eH\)]\)"];
	DefineCoupling[cleudH,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(leudH\)]\)"];
	DefineCoupling[cl2qdH$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdH\), \((1)\)]\)"];
	DefineCoupling[cl2qdH$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdH\), \((2)\)]\)"];
	DefineCoupling[cl2quH,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quH\)]\)"];
	(* 3 B-violating *)
	DefineCoupling[clud2H,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lud\), \(2\)]\)H\)]\)"];
	DefineCoupling[clq2dH,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)dH\)]\)"];
	DefineCoupling[cld3H,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(ld\), \(3\)]\)H\)]\)"];
	DefineCoupling[ceqd2H,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eqd\), \(2\)]\)H\)]\)"];
	(* 4 *)
	DefineCoupling[cleH3D,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)D\)]\)"];
	(* 5 *)
	DefineCoupling[cl2udD,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)udD\)]\)"];
	(* 5 B-violating*)
	DefineCoupling[clqd2D,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqd\), \(2\)]\)D\)]\)"];
	DefineCoupling[ced3D,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(ed\), \(3\)]\)D\)]\)"];
	(* 6 *)
	DefineCoupling[cl2H2D2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
	DefineCoupling[cl2H2D2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
];


(* Dimension 8 *)


(* 1 *)
DefineCoupling[cG4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cG4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cG4$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((3)\)]\)"];
DefineCoupling[cG4$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((4)\)]\)"];
DefineCoupling[cG4$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((5)\)]\)"];
DefineCoupling[cG4$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((6)\)]\)"];
DefineCoupling[cG4$7,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((7)\)]\)"];
DefineCoupling[cG4$8,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((8)\)]\)"];
DefineCoupling[cG4$9,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(4\)]\)\), \((9)\)]\)"];
DefineCoupling[cW4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cW4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cW4$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((3)\)]\)"];
DefineCoupling[cW4$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((4)\)]\)"];
DefineCoupling[cW4$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((5)\)]\)"];
DefineCoupling[cW4$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(4\)]\)\), \((6)\)]\)"];
DefineCoupling[cB4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cB4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cB4$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(4\)]\)\), \((3)\)]\)"];
DefineCoupling[cG3B$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)B\), \((1)\)]\)"];
DefineCoupling[cG3B$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)B\), \((2)\)]\)"];
DefineCoupling[cG3B$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)B\), \((3)\)]\)"];
DefineCoupling[cG3B$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)B\), \((4)\)]\)"];
DefineCoupling[cG2W2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cG2W2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cG2W2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cG2W2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cG2W2$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cG2W2$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cG2W2$7,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)\), \((7)\)]\)"];
DefineCoupling[cG2B2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cG2B2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cG2B2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cG2B2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cG2B2$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cG2B2$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cG2B2$7,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((7)\)]\)"];
DefineCoupling[cW2B2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cW2B2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cW2B2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cW2B2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cW2B2$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cW2B2$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cW2B2$7,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)\), \((7)\)]\)"];


(* 2 *)
DefineCoupling[cH8,SelfConjugate->True,NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(8\)]\)\)]\)"];


(* 3 *)
DefineCoupling[cH6$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(6\)]\)\), \((1)\)]\)"];
DefineCoupling[cH6$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(6\)]\)\), \((2)\)]\)"];


(* 4 *)
DefineCoupling[cH4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cH4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cH4$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((3)\)]\)"];


(* 5 *)
DefineCoupling[cG3H2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cG3H2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cW3H2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cW3H2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cW2BH2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cW2BH2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)\), \((2)\)]\)"];


(* 6 *)
DefineCoupling[cG2H4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cG2H4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cW2H4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cW2H4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cW2H4$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((3)\)]\)"];
DefineCoupling[cW2H4$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((4)\)]\)"];
DefineCoupling[cWBH4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cWBH4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(4\)]\)\), \((2)\)]\)"];
DefineCoupling[cB2H4$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((1)\)]\)"];
DefineCoupling[cB2H4$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)\), \((2)\)]\)"];


(* 7 *)
DefineCoupling[cG2H2D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cG2H2D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cG2H2D2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(G\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cW2H2D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cW2H2D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cW2H2D2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cW2H2D2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cW2H2D2$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cW2H2D2$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(W\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cWBH2D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cWBH2D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cWBH2D2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cWBH2D2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cWBH2D2$5,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cWBH2D2$6,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WBH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cB2H2D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cB2H2D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cB2H2D2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(B\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];


(* 8 *)
DefineCoupling[cWH4D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cWH4D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cWH4D2$3,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cWH4D2$4,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(WH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cBH4D2$1,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(BH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cBH4D2$2,SelfConjugate->True,NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(BH\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];


(* 9 *)
DefineCoupling[cleG2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leG\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cleG2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leG\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cleW2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leW\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cleW2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leW\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cleW2H$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leW\), \(2\)]\)H\), \((3)\)]\)"];
DefineCoupling[cquG2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quG\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cquG2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quG\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cquG2H$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quG\), \(2\)]\)H\), \((3)\)]\)"];
DefineCoupling[cquG2H$4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quG\), \(2\)]\)H\), \((4)\)]\)"];
DefineCoupling[cquG2H$5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quG\), \(2\)]\)H\), \((5)\)]\)"];
DefineCoupling[cquGWH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGWH\), \((1)\)]\)"];
DefineCoupling[cquGWH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGWH\), \((2)\)]\)"];
DefineCoupling[cquGWH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGWH\), \((3)\)]\)"];
DefineCoupling[cquGBH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGBH\), \((1)\)]\)"];
DefineCoupling[cquGBH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGBH\), \((2)\)]\)"];
DefineCoupling[cquGBH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quGBH\), \((3)\)]\)"];
DefineCoupling[cquW2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quW\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cquW2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quW\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cquW2H$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quW\), \(2\)]\)H\), \((3)\)]\)"];
DefineCoupling[cquWBH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quWBH\), \((1)\)]\)"];
DefineCoupling[cquWBH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quWBH\), \((2)\)]\)"];
DefineCoupling[cquWBH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quWBH\), \((3)\)]\)"];
DefineCoupling[cquB2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quB\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cquB2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quB\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cleWBH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leWBH\), \((1)\)]\)"];
DefineCoupling[cleWBH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leWBH\), \((2)\)]\)"];
DefineCoupling[cleWBH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leWBH\), \((3)\)]\)"];
DefineCoupling[cleB2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leB\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cleB2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leB\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cqdG2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdG\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cqdG2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdG\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cqdG2H$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdG\), \(2\)]\)H\), \((3)\)]\)"];
DefineCoupling[cqdG2H$4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdG\), \(2\)]\)H\), \((4)\)]\)"];
DefineCoupling[cqdG2H$5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdG\), \(2\)]\)H\), \((5)\)]\)"];
DefineCoupling[cqdGWH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGWH\), \((1)\)]\)"];
DefineCoupling[cqdGWH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGWH\), \((2)\)]\)"];
DefineCoupling[cqdGWH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGWH\), \((3)\)]\)"];
DefineCoupling[cqdGBH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGBH\), \((1)\)]\)"];
DefineCoupling[cqdGBH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGBH\), \((2)\)]\)"];
DefineCoupling[cqdGBH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdGBH\), \((3)\)]\)"];
DefineCoupling[cqdW2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdW\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cqdW2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdW\), \(2\)]\)H\), \((2)\)]\)"];
DefineCoupling[cqdW2H$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdW\), \(2\)]\)H\), \((3)\)]\)"];
DefineCoupling[cqdWBH$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdWBH\), \((1)\)]\)"];
DefineCoupling[cqdWBH$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdWBH\), \((2)\)]\)"];
DefineCoupling[cqdWBH$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qdWBH\), \((3)\)]\)"];
DefineCoupling[cqdB2H$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdB\), \(2\)]\)H\), \((1)\)]\)"];
DefineCoupling[cqdB2H$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdB\), \(2\)]\)H\), \((2)\)]\)"];



(* 10 *)
DefineCoupling[cleWH3$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leWH\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cleWH3$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leWH\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cleBH3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leBH\), \(3\)]\)\)]\)"];
DefineCoupling[cquGH3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quGH\), \(3\)]\)\)]\)"];
DefineCoupling[cquWH3$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quWH\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cquWH3$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quWH\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cquBH3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quBH\), \(3\)]\)\)]\)"];
DefineCoupling[cqdGH3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdGH\), \(3\)]\)\)]\)"];
DefineCoupling[cqdWH3$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdWH\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cqdWH3$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdWH\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cqdBH3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdBH\), \(3\)]\)\)]\)"];



(* 11 *)
DefineCoupling[cl2H2D3$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2H2D3$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2H2D3$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((3)\)]\)"];
DefineCoupling[cl2H2D3$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((4)\)]\)"];
DefineCoupling[ce2H2D3$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[ce2H2D3$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2H2D3$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2H2D3$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2H2D3$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2H2D3$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((4)\)]\)"];
DefineCoupling[cu2H2D3$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cu2H2D3$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cd2H2D3$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((1)\)]\)"];
DefineCoupling[cd2H2D3$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\), \((2)\)]\)"];
DefineCoupling[cudH2D3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(3\)]\)\)]\)"];


(* 12 *)
DefineCoupling[cleH5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(5\)]\)\)]\)"];
DefineCoupling[cquH5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(5\)]\)\)]\)"];
DefineCoupling[cqdH5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(5\)]\)\)]\)"];


(* 13 *)
DefineCoupling[cl2H4D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((1)\)]\)"];
DefineCoupling[cl2H4D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((2)\)]\)"];
DefineCoupling[cl2H4D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((3)\)]\)"];
DefineCoupling[cl2H4D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((4)\)]\)"];
DefineCoupling[ce2H4D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\)]\)"];
DefineCoupling[cq2H4D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2H4D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2H4D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2H4D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\), \((4)\)]\)"];
DefineCoupling[cu2H4D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\)]\)"];
DefineCoupling[cd2H4D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(4\)]\)D\)]\)"];
DefineCoupling[cudH4D,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udH\), \(4\)]\)D\)]\)"];


(* 14 *)
DefineCoupling[cq2G2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2G2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2G2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2G2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cq2G2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cq2W2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2W2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2W2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2W2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cq2B2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)D\)]\)"];
DefineCoupling[cu2G2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cu2G2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cu2G2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cu2G2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cu2G2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cu2W2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\)]\)"];
DefineCoupling[cu2B2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)D\)]\)"];
DefineCoupling[cd2G2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cd2G2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cd2G2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cd2G2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cd2G2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cd2W2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\)]\)"];
DefineCoupling[cd2B2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)D\)]\)"];
DefineCoupling[cq2GWD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GWD\), \((1)\)]\)"];
DefineCoupling[cq2GWD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GWD\), \((2)\)]\)"];
DefineCoupling[cq2GWD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GWD\), \((3)\)]\)"];
DefineCoupling[cq2GWD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GWD\), \((4)\)]\)"];
DefineCoupling[cq2GBD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GBD\), \((1)\)]\)"];
DefineCoupling[cq2GBD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GBD\), \((2)\)]\)"];
DefineCoupling[cq2GBD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GBD\), \((3)\)]\)"];
DefineCoupling[cq2GBD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)GBD\), \((4)\)]\)"];
DefineCoupling[cq2WBD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)WBD\), \((1)\)]\)"];
DefineCoupling[cq2WBD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)WBD\), \((2)\)]\)"];
DefineCoupling[cq2WBD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)WBD\), \((3)\)]\)"];
DefineCoupling[cq2WBD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)WBD\), \((4)\)]\)"];
DefineCoupling[cu2GBD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)GBD\), \((1)\)]\)"];
DefineCoupling[cu2GBD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)GBD\), \((2)\)]\)"];
DefineCoupling[cu2GBD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)GBD\), \((3)\)]\)"];
DefineCoupling[cu2GBD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)GBD\), \((4)\)]\)"];
DefineCoupling[cd2GBD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)GBD\), \((1)\)]\)"];
DefineCoupling[cd2GBD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)GBD\), \((2)\)]\)"];
DefineCoupling[cd2GBD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)GBD\), \((3)\)]\)"];
DefineCoupling[cd2GBD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)GBD\), \((4)\)]\)"];
DefineCoupling[cl2G2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\)]\)"];
DefineCoupling[cl2W2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cl2W2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cl2W2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cl2W2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cl2B2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)D\)]\)"];
DefineCoupling[ce2G2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(G\), \(2\)]\)D\)]\)"];
DefineCoupling[ce2W2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(W\), \(2\)]\)D\)]\)"];
DefineCoupling[ce2B2D,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(B\), \(2\)]\)D\)]\)"];
DefineCoupling[cl2WBD$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)WBD\), \((1)\)]\)"];
DefineCoupling[cl2WBD$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)WBD\), \((2)\)]\)"];
DefineCoupling[cl2WBD$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)WBD\), \((3)\)]\)"];
DefineCoupling[cl2WBD$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)WBD\), \((4)\)]\)"];


(* 15 RR *)
DefineCoupling[ce2WH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[ce2WH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[ce2WH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[ce2WH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[ce2BH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[ce2BH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[ce2BH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[ce2BH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cu2GH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cu2GH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cu2GH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cu2GH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cu2WH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cu2WH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cu2WH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cu2WH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cu2BH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cu2BH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cu2BH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cu2BH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cd2GH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cd2GH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cd2GH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cd2GH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cd2WH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cd2WH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cd2WH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cd2WH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cd2BH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cd2BH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cd2BH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cd2BH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cudGH2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udGH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cudGH2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udGH\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cudWH2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udWH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cudWH2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udWH\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cudBH2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udBH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cudBH2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(udBH\), \(2\)]\)\), \((2)\)]\)"];


(* 15 LL *)
DefineCoupling[cl2WH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cl2WH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cl2WH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cl2WH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cl2WH2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cl2WH2D$6,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((6)\)]\)"];
DefineCoupling[cl2WH2D$7,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((7)\)]\)"];
DefineCoupling[cl2WH2D$8,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((8)\)]\)"];
DefineCoupling[cl2WH2D$9,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((9)\)]\)"];
DefineCoupling[cl2WH2D$10,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((10)\)]\)"];
DefineCoupling[cl2WH2D$11,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((11)\)]\)"];
DefineCoupling[cl2WH2D$12,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((12)\)]\)"];
DefineCoupling[cl2BH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cl2BH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cl2BH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cl2BH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cl2BH2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cl2BH2D$6,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((6)\)]\)"];
DefineCoupling[cl2BH2D$7,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((7)\)]\)"];
DefineCoupling[cl2BH2D$8,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((8)\)]\)"];
DefineCoupling[cq2GH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2GH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2GH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2GH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cq2GH2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cq2GH2D$6,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((6)\)]\)"];
DefineCoupling[cq2GH2D$7,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((7)\)]\)"];
DefineCoupling[cq2GH2D$8,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(GH\), \(2\)]\)D\), \((8)\)]\)"];
DefineCoupling[cq2WH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2WH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2WH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2WH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cq2WH2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cq2WH2D$6,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((6)\)]\)"];
DefineCoupling[cq2WH2D$7,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((7)\)]\)"];
DefineCoupling[cq2WH2D$8,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((8)\)]\)"];
DefineCoupling[cq2WH2D$9,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((9)\)]\)"];
DefineCoupling[cq2WH2D$10,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((10)\)]\)"];
DefineCoupling[cq2WH2D$11,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((11)\)]\)"];
DefineCoupling[cq2WH2D$12,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(WH\), \(2\)]\)D\), \((12)\)]\)"];
DefineCoupling[cq2BH2D$1,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((1)\)]\)"];
DefineCoupling[cq2BH2D$2,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((2)\)]\)"];
DefineCoupling[cq2BH2D$3,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((3)\)]\)"];
DefineCoupling[cq2BH2D$4,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((4)\)]\)"];
DefineCoupling[cq2BH2D$5,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((5)\)]\)"];
DefineCoupling[cq2BH2D$6,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((6)\)]\)"];
DefineCoupling[cq2BH2D$7,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((7)\)]\)"];
DefineCoupling[cq2BH2D$8,Indices->{Flavor,Flavor},SelfConjugate->{2,1},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(BH\), \(2\)]\)D\), \((8)\)]\)"];


(* 16 *)
DefineCoupling[cleWHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leWHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cleWHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leWHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cleWHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leWHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cleBHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leBHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cleBHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leBHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cleBHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leBHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cquGHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quGHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cquGHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quGHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cquGHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quGHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cquWHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quWHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cquWHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quWHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cquWHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quWHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cquBHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quBHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cquBHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quBHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cquBHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quBHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cqdGHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdGHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cqdGHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdGHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cqdGHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdGHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cqdWHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdWHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cqdWHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdWHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cqdWHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdWHD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cqdBHD2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdBHD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cqdBHD2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdBHD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cqdBHD2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdBHD\), \(2\)]\)\), \((3)\)]\)"];


(* 17 *)
DefineCoupling[cleH3D2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cleH3D2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cleH3D2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cleH3D2$4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cleH3D2$5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cleH3D2$6,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cquH3D2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cquH3D2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cquH3D2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cquH3D2$4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cquH3D2$5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cquH3D2$6,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(quH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cqdH3D2$1,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cqdH3D2$2,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cqdH3D2$3,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cqdH3D2$4,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cqdH3D2$5,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cqdH3D2$6,Indices->{Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qdH\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((6)\)]\)"];


(* 18 LRRLH2 *)
DefineCoupling[cleqdH2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cleqdH2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdH\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2udH2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\)]\)"];
DefineCoupling[clequH2$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequH\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cq2udH2$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cq2udH2$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((6)\)]\)"];


(* 18 B-violating *)
If[TrueQ@ BLViolation,
	DefineCoupling[clqudH2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqudH\), \(2\)]\)\), \((1)\)]\)"];
	DefineCoupling[clqudH2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqudH\), \(2\)]\)\), \((2)\)]\)"];
	DefineCoupling[ceq2uH2,Indices->{Flavor,Flavor,Flavor,Flavor},Symmetries->{SymmetricPermutation[2,1,3,4]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)\!\(\*SuperscriptBox[\(uH\), \(2\)]\)\)]\)"];
	DefineCoupling[clq3H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"(* more symmetries *)];
	DefineCoupling[clq3H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"(* more symmetries *)];
	DefineCoupling[clq3H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"(* more symmetries *)];
	DefineCoupling[ceu2dH2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)\!\(\*SuperscriptBox[\(dH\), \(2\)]\)\)]\)"];
	DefineCoupling[clqu2H2,Indices->{Flavor,Flavor,Flavor,Flavor},Symmetries->{SymmetricPermutation[1,2,4,3]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqu\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
	DefineCoupling[clqd2H2,Indices->{Flavor,Flavor,Flavor,Flavor},Symmetries->{SymmetricPermutation[1,2,4,3]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqd\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
	DefineCoupling[ceq2dH2,Indices->{Flavor,Flavor,Flavor,Flavor},Symmetries->{SymmetricPermutation[1,2,4,3]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)\!\(\*SuperscriptBox[\(dH\), \(2\)]\)\)]\)"];
];


(* 18 LLLLH2 *)
DefineCoupling[cl4H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl4H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq4H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq4H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq4H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"];
(* The label cq4H2$4 is not used as it is redundant. See Eq. (3.20) of https://arxiv.org/pdf/2005.00059 for details *)
DefineCoupling[cq4H2$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cl2q2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2q2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2q2H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cl2q2H2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cl2q2H2$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((5)\)]\)"];


(* 18 RRRRH2 *)
DefineCoupling[ce4H2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2],SymmetricPermutation[1,4,3,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
DefineCoupling[cu4H2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
DefineCoupling[cd4H2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
DefineCoupling[ce2u2H2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
DefineCoupling[ce2d2H2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\)]\)"];
DefineCoupling[cu2d2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cu2d2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];


(* 18 LLRRH2 *)
DefineCoupling[cl2e2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2e2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2u2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2u2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2d2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2d2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2e2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2e2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2u2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2u2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2u2H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2u2H2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cq2d2H2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2d2H2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2d2H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2d2H2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((4)\)]\)"];


(* 18 LRLRH2 *)
DefineCoupling[cq2udH2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2udH2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2udH2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2udH2$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udH\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[clequH2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequH\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[clequH2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequH\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[clequH2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequH\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[clequH2$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequH\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cl2e2H2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cleqdH2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdH\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cleqdH2$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdH\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cq2u2H2$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cq2u2H2$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((6)\)]\)"];
DefineCoupling[cq2d2H2$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((5)\)]\)"];
DefineCoupling[cq2d2H2$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(H\), \(2\)]\)\), \((6)\)]\)"];


(* 19 LLLLX *)
DefineCoupling[cl4W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)W\), \((1)\)]\)"];
DefineCoupling[cl4W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)W\), \((2)\)]\)"];
DefineCoupling[cq4G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)G\), \((1)\)]\)"];
DefineCoupling[cq4G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)G\), \((2)\)]\)"];
DefineCoupling[cq4G$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)G\), \((3)\)]\)"];
DefineCoupling[cq4G$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)G\), \((4)\)]\)"];
DefineCoupling[cq4W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)W\), \((1)\)]\)"];
DefineCoupling[cq4W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)W\), \((2)\)]\)"];
DefineCoupling[cq4W$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)W\), \((3)\)]\)"];
DefineCoupling[cq4W$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)W\), \((4)\)]\)"];
DefineCoupling[cl2q2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cl2q2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cl2q2G$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)G\), \((3)\)]\)"];
DefineCoupling[cl2q2G$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)G\), \((4)\)]\)"];
DefineCoupling[cl2q2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cl2q2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cl2q2W$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((3)\)]\)"];
DefineCoupling[cl2q2W$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((4)\)]\)"];
DefineCoupling[cl2q2W$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((5)\)]\)"];
DefineCoupling[cl2q2W$6,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)W\), \((6)\)]\)"];
DefineCoupling[cl2q2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cl2q2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cl2q2B$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)B\), \((3)\)]\)"];
DefineCoupling[cl2q2B$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)B\), \((4)\)]\)"];
DefineCoupling[cl4B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)B\), \((1)\)]\)"];
DefineCoupling[cl4B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq4B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)B\), \((1)\)]\)"];
DefineCoupling[cq4B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq4B$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)B\), \((3)\)]\)"];
DefineCoupling[cq4B$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)B\), \((4)\)]\)"];


(* 19 RRRRX *)
DefineCoupling[cu4G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)G\), \((1)\)]\)"];
DefineCoupling[cu4G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)G\), \((2)\)]\)"];
DefineCoupling[cd4G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)G\), \((1)\)]\)"];
DefineCoupling[cd4G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)G\), \((2)\)]\)"];
DefineCoupling[ce2u2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[ce2u2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[ce2u2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[ce2u2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[ce2d2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[ce2d2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[ce2d2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[ce2d2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cu2d2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cu2d2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cu2d2G$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((3)\)]\)"];
DefineCoupling[cu2d2G$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((4)\)]\)"];
DefineCoupling[cu2d2G$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((5)\)]\)"];
DefineCoupling[cu2d2G$6,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((6)\)]\)"];
DefineCoupling[cu2d2G$7,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((7)\)]\)"];
DefineCoupling[cu2d2G$8,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((8)\)]\)"];
DefineCoupling[cu2d2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cu2d2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cu2d2B$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((3)\)]\)"];
DefineCoupling[cu2d2B$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((4)\)]\)"];
DefineCoupling[ce4B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(4\)]\)B\), \((1)\)]\)"(*more symmetries ?*)];
DefineCoupling[ce4B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(4\)]\)B\), \((2)\)]\)"(*more symmetries ?*)];
DefineCoupling[cu4B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)B\), \((1)\)]\)"];
DefineCoupling[cu4B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)B\), \((2)\)]\)"];
DefineCoupling[cd4B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)B\), \((1)\)]\)"];
DefineCoupling[cd4B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)B\), \((2)\)]\)"];


(* 19 LLRRX *)
DefineCoupling[cl2e2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cl2e2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cl2e2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cl2e2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cl2u2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cl2u2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cl2u2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cl2u2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cl2u2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cl2u2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cl2d2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cl2d2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cl2d2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cl2d2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cl2d2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cl2d2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq2e2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cq2e2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cq2e2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cq2e2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cq2e2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cq2e2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq2u2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cq2u2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cq2u2G$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((3)\)]\)"];
DefineCoupling[cq2u2G$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((4)\)]\)"];
DefineCoupling[cq2u2G$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((5)\)]\)"];
DefineCoupling[cq2u2G$6,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((6)\)]\)"];
DefineCoupling[cq2u2G$7,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((7)\)]\)"];
DefineCoupling[cq2u2G$8,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)G\), \((8)\)]\)"];
DefineCoupling[cq2u2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cq2u2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cq2u2W$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((3)\)]\)"];
DefineCoupling[cq2u2W$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)W\), \((4)\)]\)"];
DefineCoupling[cq2u2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cq2u2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq2u2B$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((3)\)]\)"];
DefineCoupling[cq2u2B$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)B\), \((4)\)]\)"];
DefineCoupling[cq2d2G$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((1)\)]\)"];
DefineCoupling[cq2d2G$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((2)\)]\)"];
DefineCoupling[cq2d2G$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((3)\)]\)"];
DefineCoupling[cq2d2G$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((4)\)]\)"];
DefineCoupling[cq2d2G$5,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((5)\)]\)"];
DefineCoupling[cq2d2G$6,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((6)\)]\)"];
DefineCoupling[cq2d2G$7,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((7)\)]\)"];
DefineCoupling[cq2d2G$8,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)G\), \((8)\)]\)"];
DefineCoupling[cq2d2W$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((1)\)]\)"];
DefineCoupling[cq2d2W$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((2)\)]\)"];
DefineCoupling[cq2d2W$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((3)\)]\)"];
DefineCoupling[cq2d2W$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)W\), \((4)\)]\)"];
DefineCoupling[cq2d2B$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((1)\)]\)"];
DefineCoupling[cq2d2B$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((2)\)]\)"];
DefineCoupling[cq2d2B$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((3)\)]\)"];
DefineCoupling[cq2d2B$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)B\), \((4)\)]\)"];


(* 19 LRRLX *)
DefineCoupling[cledqG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqG\), \((1)\)]\)"];
DefineCoupling[cledqG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqG\), \((2)\)]\)"];
DefineCoupling[cledqW$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqW\), \((1)\)]\)"];
DefineCoupling[cledqW$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqW\), \((2)\)]\)"];
DefineCoupling[cledqB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqB\), \((1)\)]\)"];
DefineCoupling[cledqB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ledqB\), \((2)\)]\)"];


(* 19 LRLRX *)
DefineCoupling[cq2udG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((1)\)]\)"];
DefineCoupling[cq2udG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((2)\)]\)"];
DefineCoupling[cq2udG$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((3)\)]\)"];
DefineCoupling[cq2udG$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((4)\)]\)"];
DefineCoupling[cq2udG$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((5)\)]\)"];
DefineCoupling[cq2udG$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udG\), \((6)\)]\)"];
DefineCoupling[cq2udW$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udW\), \((1)\)]\)"];
DefineCoupling[cq2udW$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udW\), \((2)\)]\)"];
DefineCoupling[cq2udW$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udW\), \((3)\)]\)"];
DefineCoupling[cq2udB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udB\), \((1)\)]\)"];
DefineCoupling[cq2udB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udB\), \((2)\)]\)"];
DefineCoupling[cq2udB$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)udB\), \((3)\)]\)"];
DefineCoupling[clequG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequG\), \((1)\)]\)"];
DefineCoupling[clequG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequG\), \((2)\)]\)"];
DefineCoupling[clequG$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequG\), \((3)\)]\)"];
DefineCoupling[clequW$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequW\), \((1)\)]\)"];
DefineCoupling[clequW$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequW\), \((2)\)]\)"];
DefineCoupling[clequW$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequW\), \((3)\)]\)"];
DefineCoupling[clequB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequB\), \((1)\)]\)"];
DefineCoupling[clequB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequB\), \((2)\)]\)"];
DefineCoupling[clequB$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequB\), \((3)\)]\)"];


(* 19 B-violating *)
If[TrueQ@ BLViolation,
	DefineCoupling[clqudG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudG\), \((1)\)]\)"];
	DefineCoupling[clqudG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudG\), \((2)\)]\)"];
	DefineCoupling[clqudG$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudG\), \((3)\)]\)"];
	DefineCoupling[clqudG$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudG\), \((4)\)]\)"];
	DefineCoupling[clqudW$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudW\), \((1)\)]\)"];
	DefineCoupling[clqudW$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudW\), \((2)\)]\)"];
	DefineCoupling[clqudB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudB\), \((1)\)]\)"];
	DefineCoupling[clqudB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lqudB\), \((2)\)]\)"];
	DefineCoupling[ceq2uG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uG\), \((1)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceq2uG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uG\), \((2)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceq2uW$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uW\), \((1)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceq2uW$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uW\), \((2)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceq2uB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uB\), \((1)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceq2uB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)uB\), \((2)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[clq3G$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)G\), \((1)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3G$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)G\), \((2)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3G$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)G\), \((3)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3G$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)G\), \((4)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3W$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)W\), \((1)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3W$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)W\), \((2)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3W$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)W\), \((3)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3B$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)B\), \((1)\)]\)"(*more symmetries*)];
	DefineCoupling[clq3B$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)B\), \((2)\)]\)"(*more symmetries*)];
	DefineCoupling[ceu2dG$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)dG\), \((1)\)]\)"];
	DefineCoupling[ceu2dG$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)dG\), \((2)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceu2dG$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)dG\), \((3)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
	DefineCoupling[ceu2dB$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)dB\), \((1)\)]\)"];
	DefineCoupling[ceu2dB$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)dB\), \((2)\)]\)"(*,Symmetries->{SymmetricPermutation[2,1,3,4]}*)];
];


(* 20 *)
DefineCoupling[cl3eHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(3\)]\)eHD\), \((1)\)]\)"];
DefineCoupling[cl3eHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(3\)]\)eHD\), \((2)\)]\)"];
DefineCoupling[cl3eHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(3\)]\)eHD\), \((3)\)]\)"];
DefineCoupling[cle3HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(le\), \(3\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cle3HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(le\), \(3\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cleq2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cleq2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cleq2HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((3)\)]\)"];
DefineCoupling[cleq2HD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((4)\)]\)"];
DefineCoupling[cleq2HD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((5)\)]\)"];
DefineCoupling[cleq2HD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leq\), \(2\)]\)HD\), \((6)\)]\)"];
DefineCoupling[cleu2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leu\), \(2\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cleu2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leu\), \(2\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cleu2HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leu\), \(2\)]\)HD\), \((3)\)]\)"];
DefineCoupling[cled2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(led\), \(2\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cled2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(led\), \(2\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cled2HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(led\), \(2\)]\)HD\), \((3)\)]\)"];
DefineCoupling[cleudHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leudHD\), \((1)\)]\)"];
DefineCoupling[cleudHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leudHD\), \((2)\)]\)"];
DefineCoupling[cleudHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(leudHD\), \((3)\)]\)"];
DefineCoupling[cl2quHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((1)\)]\)"];
DefineCoupling[cl2quHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((2)\)]\)"];
DefineCoupling[cl2quHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((3)\)]\)"];
DefineCoupling[cl2quHD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((4)\)]\)"];
DefineCoupling[cl2quHD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((5)\)]\)"];
DefineCoupling[cl2quHD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)quHD\), \((6)\)]\)"];
DefineCoupling[ce2quHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)quHD\), \((1)\)]\)"];
DefineCoupling[ce2quHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)quHD\), \((2)\)]\)"];
DefineCoupling[ce2quHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)quHD\), \((3)\)]\)"];
DefineCoupling[cq3uHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((1)\)]\)"];
DefineCoupling[cq3uHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((2)\)]\)"];
DefineCoupling[cq3uHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((3)\)]\)"];
DefineCoupling[cq3uHD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((4)\)]\)"];
DefineCoupling[cq3uHD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((5)\)]\)"];
DefineCoupling[cq3uHD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)uHD\), \((6)\)]\)"];
DefineCoupling[cqu3HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(3\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cqu3HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(3\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cqu3HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(3\)]\)HD\), \((3)\)]\)"];
DefineCoupling[cqud2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cqud2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cqud2HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((3)\)]\)"];
DefineCoupling[cqud2HD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((4)\)]\)"];
DefineCoupling[cqud2HD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((5)\)]\)"];
DefineCoupling[cqud2HD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qud\), \(2\)]\)HD\), \((6)\)]\)"];
DefineCoupling[cl2qdHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((1)\)]\)"];
DefineCoupling[cl2qdHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((2)\)]\)"];
DefineCoupling[cl2qdHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((3)\)]\)"];
DefineCoupling[cl2qdHD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((4)\)]\)"];
DefineCoupling[cl2qdHD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((5)\)]\)"];
DefineCoupling[cl2qdHD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)qdHD\), \((6)\)]\)"];
DefineCoupling[ce2qdHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)qdHD\), \((1)\)]\)"];
DefineCoupling[ce2qdHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)qdHD\), \((2)\)]\)"];
DefineCoupling[ce2qdHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)qdHD\), \((3)\)]\)"];
DefineCoupling[cq3dHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((1)\)]\)"];
DefineCoupling[cq3dHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((2)\)]\)"];
DefineCoupling[cq3dHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((3)\)]\)"];
DefineCoupling[cq3dHD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((4)\)]\)"];
DefineCoupling[cq3dHD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((5)\)]\)"];
DefineCoupling[cq3dHD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(3\)]\)dHD\), \((6)\)]\)"];
DefineCoupling[cqu2dHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((1)\)]\)"];
DefineCoupling[cqu2dHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((2)\)]\)"];
DefineCoupling[cqu2dHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((3)\)]\)"];
DefineCoupling[cqu2dHD$4,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((4)\)]\)"];
DefineCoupling[cqu2dHD$5,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((5)\)]\)"];
DefineCoupling[cqu2dHD$6,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qu\), \(2\)]\)dHD\), \((6)\)]\)"];
DefineCoupling[cqd3HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qd\), \(3\)]\)HD\), \((1)\)]\)"];
DefineCoupling[cqd3HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qd\), \(3\)]\)HD\), \((2)\)]\)"];
DefineCoupling[cqd3HD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(qd\), \(3\)]\)HD\), \((3)\)]\)"];


(* 20 B-violating *)
If[TrueQ@ BLViolation,
	DefineCoupling[clu2dHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lu\), \(2\)]\)dHD\), \((1)\)]\)"];
	DefineCoupling[clu2dHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lu\), \(2\)]\)dHD\), \((2)\)]\)"];
	DefineCoupling[clud2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lud\), \(2\)]\)HD\), \((1)\)]\)"];
	DefineCoupling[clud2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lud\), \(2\)]\)HD\), \((2)\)]\)"];
	DefineCoupling[clq2uHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)uHD\), \((1)\)]\)"];
	DefineCoupling[clq2uHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)uHD\), \((2)\)]\)"];
	DefineCoupling[clq2uHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)uHD\), \((3)\)]\)"];
	DefineCoupling[clq2dHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)dHD\), \((1)\)]\)"];
	DefineCoupling[clq2dHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)dHD\), \((2)\)]\)"];
	DefineCoupling[clq2dHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(2\)]\)dHD\), \((3)\)]\)"];
	DefineCoupling[ceq3HD,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(3\)]\)HD\)]\)"];
	DefineCoupling[cequ2HD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(equ\), \(2\)]\)HD\), \((1)\)]\)"];
	DefineCoupling[cequ2HD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(equ\), \(2\)]\)HD\), \((2)\)]\)"];
	DefineCoupling[cequdHD$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(equdHD\), \((1)\)]\)"];
	DefineCoupling[cequdHD$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(equdHD\), \((2)\)]\)"];
	DefineCoupling[cequdHD$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(equdHD\), \((3)\)]\)"];
];


(* 21 LLLLD2 *)
DefineCoupling[cl4D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl4D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq4D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq4D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq4D2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq4D2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cl2q2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2q2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2q2D2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cl2q2D2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];


(* 21 RRRRD2*)
DefineCoupling[ce4D2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\)]\)"];
DefineCoupling[cu4D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cu4D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cd4D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cd4D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(d\), \(4\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[ce2u2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[ce2u2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[ce2d2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[ce2d2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cu2d2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cu2d2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cu2d2D2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cu2d2D2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];


(* 21 LLRRD2 *)
DefineCoupling[cl2e2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2e2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2u2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2u2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cl2d2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cl2d2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(l\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2e2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2e2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(e\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2u2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2u2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2u2D2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2u2D2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(u\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];
DefineCoupling[cq2d2D2$1,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2d2D2$2,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2d2D2$3,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[cq2d2D2$4,Indices->{Flavor,Flavor,Flavor,Flavor},SelfConjugate->{2,1,4,3},Symmetries->{SymmetricPermutation[3,4,1,2]},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(d\), \(2\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\), \((4)\)]\)"];


(* 21 LRLRD2 *)
DefineCoupling[cq2udD2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cq2udD2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[cq2udD2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(q\), \(2\)]\)\!\(\*SuperscriptBox[\(udD\), \(2\)]\)\), \((3)\)]\)"];
DefineCoupling[clequD2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[clequD2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequD\), \(2\)]\)\), \((2)\)]\)"];
DefineCoupling[clequD2$3,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lequD\), \(2\)]\)\), \((3)\)]\)"];


(* 21 LRRLD2 *)
DefineCoupling[cleqdD2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdD\), \(2\)]\)\), \((1)\)]\)"];
DefineCoupling[cleqdD2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(leqdD\), \(2\)]\)\), \((2)\)]\)"];


(* 21 B-violating *)
If[TrueQ@ BLViolation,
	DefineCoupling[clqudD2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqudD\), \(2\)]\)\), \((1)\)]\)"];
	DefineCoupling[clqudD2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lqudD\), \(2\)]\)\), \((2)\)]\)"];
	DefineCoupling[ceq2uD2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eq\), \(2\)]\)\!\(\*SuperscriptBox[\(uD\), \(2\)]\)\)]\)"];
	DefineCoupling[clq3D2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(lq\), \(3\)]\)\!\(\*SuperscriptBox[\(D\), \(2\)]\)\)]\)"];
	DefineCoupling[ceu2dD2$1,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)\!\(\*SuperscriptBox[\(dD\), \(2\)]\)\), \((1)\)]\)"];
	DefineCoupling[ceu2dD2$2,Indices->{Flavor,Flavor,Flavor,Flavor},NiceForm->"\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(\!\(\*SuperscriptBox[\(eu\), \(2\)]\)\!\(\*SuperscriptBox[\(dD\), \(2\)]\)\), \((2)\)]\)"];
];


(* ::Subsubsection:: *)
(*Lagrangian*)


Module[{p,r,s,t,i,j,k,m,n,o,J,K,L,\[Alpha],\[Beta],\[Delta],A,C,D,\[Nu],\[Rho],\[Theta],\[Eta],\[Kappa],
		\[Tau]SU2L,fSU2L,\[Epsilon]SU2L,TSU3c,fSU3c,\[Epsilon]SU3c,dSU3c,HermitianCD,
		L7$1,L7$2,L7$3,L7$4,L7$5,L7$6,
		L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21,
		\[Mu]1,\[Mu]2,\[Mu]3,\[Mu]4,\[Mu]5,\[Mu]6,
		i1,i2,i3,i4,i5,i6,i7,i8,I1,I2,I3,
		a1,a2,a3,a4,a5,a6,a7,a8,A1,A2,A3,A4,A5,
		dFS,dH,
		symCD,
		Dim7Lag,
		BosonicLag, TwoFermionLag, FourFermionLag
		},

	(* generators and other CG objects *)
	\[Tau]SU2L[Jadj_,ifund_,jfund_]:= 2 CG[gen[SU2L[fund]],{Jadj,ifund,jfund}];
	fSU2L[Iadj_,Jadj_,Kadj_]:=CG[fStruct[SU2L],{Iadj,Jadj,Kadj}];
	\[Epsilon]SU2L[ifund_,jfund_]:= CG[eps[SU2L],{ifund,jfund}];
	TSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_]:=CG[gen[SU3c[fund]],{Aadj,\[Alpha]fund,\[Beta]fund}];
	fSU3c[Aadj_,Badj_,Cadj_]:=CG[fStruct[SU3c],{Aadj,Badj,Cadj}];
	\[Epsilon]SU3c[\[Alpha]fund_,\[Beta]fund_,\[Gamma]fund_]:= CG[eps[SU3c],{\[Alpha]fund,\[Beta]fund,\[Gamma]fund}];
	dSU3c[Aadj_,Badj_,Cadj_]:= CG[dSym[SU3c],{Aadj,Badj,Cadj}];
	(* Hermitian derivative: A factor of \[ImaginaryI] has been added w.r.t. Murphy's paper to make the derivative hermitian *)
	HermitianCD[ind_,term1_,term2_]:= I term1**CD[ind,term2] - I CD[ind,term1]**term2;
	dFS[l_,\[Nu]1_,\[Nu]2_,X_]:= RelabelIndices[1/2 LCTensor[\[Nu]1,\[Nu]2,\[Eta],\[Kappa]]FS[l,\[Eta],\[Kappa],X],Unique->True];
	dFS[l_,\[Nu]1_,\[Nu]2_]:= RelabelIndices[1/2 LCTensor[\[Nu]1,\[Nu]2,\[Eta],\[Kappa]]FS[l,\[Eta],\[Kappa]],Unique->True];
	dH[i_]:= RelabelIndices[\[Epsilon]SU2L[i,j]Bar[H[j]],Unique->True];
	symCD[inds_List,arg_]:=(1/Length[inds])Plus@@(CD[#,arg]&/@ Permutations[inds]);
	(*asymCD[inds_List,arg_]:=Plus@@((Signature[#]*Signature[inds]*CD[#,arg])&/@ Permutations[inds]);*)

	(* Lagrangian terms *)

	(* Dimension 7 *)

	L7$1 = If[TrueQ@ BLViolation,
		PlusHc[
			cl2WH2[p,r] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] \[Tau]SU2L[I1,i4,i5] Transp[l[i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i3,r] H[i2] H[i5] FS[W,\[Mu]1,\[Mu]2,I1] +
			cl2BH2[p,r] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Transp[l[i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i3,r] H[i2] H[i4] FS[B,\[Mu]1,\[Mu]2]
		]
	,
		0
	];
	
	L7$2 = If[TrueQ@ BLViolation,
		PlusHc[
			cl2H4[p,r] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Transp[l[i1,p]]**CC**l[i3,r] H[i2] H[i4] Bar[H[i5]] H[i5]
		]
	,
		0
	];
	
	L7$3 = If[TrueQ@ BLViolation,
		PlusHc[
			cl3eH[p,r,s,t] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Bar[e[p]]**l[i1,r] Transp[l[i2,s]]**CC**l[i3,t] H[i4] +
			cleudH[p,r,s,t] Bar@\[Epsilon]SU2L[i1,i2] Bar[d[a1,p]]**l[i1,r] Transp[u[a1,s]]**CC**e[t] H[i2] +
			cl2qdH$1[p,r,s,t] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Bar[d[a1,p]]**l[i1,r] Transp[q[a1,i2,s]]**CC**l[i3,t] H[i4] +
			cl2qdH$2[p,r,s,t] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Bar[d[a1,p]]**l[i1,r] Transp[q[a1,i3,s]]**CC**l[i2,t] H[i4] +
			cl2quH[p,r,s,t] Bar@\[Epsilon]SU2L[i1,i2] Bar[q[a1,i3,p]]**u[a1,r] Transp[l[i3,s]]**CC**l[i1,t] H[i2] +
	
			clud2H[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[l[i1,p]]**d[a1,r] Transp[u[a2,s]]**CC**d[a3,t] dH[i1] +
			clq2dH[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[l[i1,p]]**d[a1,r] Transp[q[a2,i1,s]]**CC**q[a3,i2,t] Bar[H[i2]] +
			cld3H[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[l[i1,p]]**d[a1,r] Transp[d[a2,s]]**CC**d[a3,t] H[i1] +
			ceqd2H[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[e[p]]**q[a1,i1,r] Transp[d[a2,s]]**CC**d[a3,t] Bar[H[i1]]
		]
	,
		0
	];

	L7$4 = If[TrueQ@ BLViolation,
		PlusHc[
			cleH3D[p,r] I Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Transp[l[i1,p]]**CC**\[Gamma][\[Mu]1]**e[r] H[i2] H[i3] CD[\[Mu]1,H[i4]]
		]
	,
		0
	];
	
	L7$5 = If[TrueQ@ BLViolation,
		PlusHc[
			cl2udD[p,r,s,t] I Bar@\[Epsilon]SU2L[i1,i2] Bar[d[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r] Transp[l[i1,s]]**CC**CD[\[Mu]1,l[i2,t]] +
	
			clqd2D[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[l[i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r] Transp[d[a2,s]]**CC**CD[\[Mu]1,d[a3,t]] +
			ced3D[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[e[p]]**\[Gamma][\[Mu]1]**d[a1,r] Transp[d[a2,s]]**CC**CD[\[Mu]1,d[a3,t]]
		]
	,
		0
	];
	
	L7$6 = If[TrueQ@ BLViolation,
		PlusHc[
			cl2H2D2$1[p,r] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Transp[l[i1,p]]**CC**CD[\[Mu]1,l[i2,r]] H[i3] CD[\[Mu]1,H[i4]] +
			cl2H2D2$2[p,r] Bar@\[Epsilon]SU2L[i1,i3] Bar@\[Epsilon]SU2L[i2,i4] Transp[l[i1,p]]**CC**CD[\[Mu]1,l[i2,r]] H[i3] CD[\[Mu]1,H[i4]]
		]
	,
		0
	];
	

	(* Dimension 8 *)

	(* Bosonic sector *)

	L1 =
		cG4$1[] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]FS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$2[] (FS[G,\[Mu]1,\[Mu]2,A1]dFS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]dFS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$3[] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A1]FS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$4[] (FS[G,\[Mu]1,\[Mu]2,A1]dFS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$5[] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]dFS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$6[] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A2]) +
		cG4$7[] dSU3c[A1,A2,A5] dSU3c[A3,A4,A5] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A3]FS[G,\[Mu]3,\[Mu]4,A4]) +
		cG4$8[] dSU3c[A1,A2,A5] dSU3c[A3,A4,A5] (FS[G,\[Mu]1,\[Mu]2,A1]dFS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A3]dFS[G,\[Mu]3,\[Mu]4,A4]) +
		cG4$9[] dSU3c[A1,A2,A5] dSU3c[A3,A4,A5] (FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A2]) (FS[G,\[Mu]3,\[Mu]4,A3]dFS[G,\[Mu]3,\[Mu]4,A4]) +

		cW4$1[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[W,\[Mu]3,\[Mu]4,I2]FS[W,\[Mu]3,\[Mu]4,I2]) +
		cW4$2[] (FS[W,\[Mu]1,\[Mu]2,I1]dFS[W,\[Mu]1,\[Mu]2,I1]) (FS[W,\[Mu]3,\[Mu]4,I2]dFS[W,\[Mu]3,\[Mu]4,I2]) +
		cW4$3[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I2]) (FS[W,\[Mu]3,\[Mu]4,I1]FS[W,\[Mu]3,\[Mu]4,I2]) +
		cW4$4[] (FS[W,\[Mu]1,\[Mu]2,I1]dFS[W,\[Mu]1,\[Mu]2,I2]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[W,\[Mu]3,\[Mu]4,I2]) +
		cW4$5[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[W,\[Mu]3,\[Mu]4,I2]dFS[W,\[Mu]3,\[Mu]4,I2]) +
		cW4$6[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I2]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[W,\[Mu]3,\[Mu]4,I2]) +

		cB4$1[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[B,\[Mu]3,\[Mu]4]FS[B,\[Mu]3,\[Mu]4]) +
		cB4$2[] (FS[B,\[Mu]1,\[Mu]2]dFS[B,\[Mu]1,\[Mu]2]) (FS[B,\[Mu]3,\[Mu]4]dFS[B,\[Mu]3,\[Mu]4]) +
		cB4$3[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[B,\[Mu]3,\[Mu]4]dFS[B,\[Mu]3,\[Mu]4]) +

		cG3B$1[] dSU3c[A1,A2,A3] (FS[B,\[Mu]1,\[Mu]2]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]FS[G,\[Mu]3,\[Mu]4,A3]) +
		cG3B$2[] dSU3c[A1,A2,A3] (FS[B,\[Mu]1,\[Mu]2]dFS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]dFS[G,\[Mu]3,\[Mu]4,A3]) +
		cG3B$3[] dSU3c[A1,A2,A3] (FS[B,\[Mu]1,\[Mu]2]dFS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]FS[G,\[Mu]3,\[Mu]4,A3]) +
		cG3B$4[] dSU3c[A1,A2,A3] (FS[B,\[Mu]1,\[Mu]2]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[G,\[Mu]3,\[Mu]4,A2]dFS[G,\[Mu]3,\[Mu]4,A3]) +

		cG2W2$1[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[G,\[Mu]3,\[Mu]4,A1]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$2[] (FS[W,\[Mu]1,\[Mu]2,I1]dFS[W,\[Mu]1,\[Mu]2,I1]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$3[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[W,\[Mu]3,\[Mu]4,I1]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$4[] (FS[W,\[Mu]1,\[Mu]2,I1]dFS[G,\[Mu]1,\[Mu]2,A1]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$5[] (FS[W,\[Mu]1,\[Mu]2,I1]dFS[W,\[Mu]1,\[Mu]2,I1]) (FS[G,\[Mu]3,\[Mu]4,A1]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$6[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2W2$7[] (FS[W,\[Mu]1,\[Mu]2,I1]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[G,\[Mu]3,\[Mu]4,A1]) +

		cG2B2$1[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[G,\[Mu]3,\[Mu]4,A1]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$2[] (FS[B,\[Mu]1,\[Mu]2]dFS[B,\[Mu]1,\[Mu]2]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$3[] (FS[B,\[Mu]1,\[Mu]2]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[B,\[Mu]3,\[Mu]4]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$4[] (FS[B,\[Mu]1,\[Mu]2]dFS[G,\[Mu]1,\[Mu]2,A1]) (FS[B,\[Mu]3,\[Mu]4]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$5[] (FS[B,\[Mu]1,\[Mu]2]dFS[B,\[Mu]1,\[Mu]2]) (FS[G,\[Mu]3,\[Mu]4,A1]FS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$6[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[G,\[Mu]3,\[Mu]4,A1]dFS[G,\[Mu]3,\[Mu]4,A1]) +
		cG2B2$7[] (FS[B,\[Mu]1,\[Mu]2]FS[G,\[Mu]1,\[Mu]2,A1]) (FS[B,\[Mu]3,\[Mu]4]dFS[G,\[Mu]3,\[Mu]4,A1]) +

		cW2B2$1[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[W,\[Mu]3,\[Mu]4,I1]FS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$2[] (FS[B,\[Mu]1,\[Mu]2]dFS[B,\[Mu]1,\[Mu]2]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$3[] (FS[B,\[Mu]1,\[Mu]2]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[B,\[Mu]3,\[Mu]4]FS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$4[] (FS[B,\[Mu]1,\[Mu]2]dFS[W,\[Mu]1,\[Mu]2,I1]) (FS[B,\[Mu]3,\[Mu]4]dFS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$5[] (FS[B,\[Mu]1,\[Mu]2]dFS[B,\[Mu]1,\[Mu]2]) (FS[W,\[Mu]3,\[Mu]4,I1]FS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$6[] (FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]) (FS[W,\[Mu]3,\[Mu]4,I1]dFS[W,\[Mu]3,\[Mu]4,I1]) +
		cW2B2$7[] (FS[B,\[Mu]1,\[Mu]2]FS[W,\[Mu]1,\[Mu]2,I1]) (FS[B,\[Mu]3,\[Mu]4]dFS[W,\[Mu]3,\[Mu]4,I1])
		;

	L2 =
		cH8[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) (Bar[H[i3]] H[i3]) (Bar[H[i4]] H[i4])
		;

	L3 =
		cH6$1[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) (CD[\[Mu]1,Bar[H[i3]]] CD[\[Mu]1,H[i3]]) +
		cH6$2[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) (CD[\[Mu]1,Bar[H[i4]]] \[Tau]SU2L[I1,i4,i5] CD[\[Mu]1,H[i5]])
		;

	L4 =
		cH4$1[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) (CD[\[Mu]2,Bar[H[i2]]] CD[\[Mu]1,H[i2]]) +
		cH4$2[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) (CD[\[Mu]1,Bar[H[i2]]] CD[\[Mu]2,H[i2]]) +
		cH4$3[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) (CD[\[Mu]2,Bar[H[i2]]] CD[\[Mu]2,H[i2]])
		;

	L5 =
		cG3H2$1[] fSU3c[A1,A2,A3] (Bar[H[i1]] H[i1]) FS[G,\[Mu]1,\[Mu]2,A1] FS[G,\[Mu]2,\[Mu]3,A2] FS[G,\[Mu]3,\[Mu]1,A3] +
		cG3H2$2[] fSU3c[A1,A2,A3] (Bar[H[i1]] H[i1]) FS[G,\[Mu]1,\[Mu]2,A1] FS[G,\[Mu]2,\[Mu]3,A2] dFS[G,\[Mu]3,\[Mu]1,A3] +

		cW3H2$1[] fSU2L[I1,I2,I3] (Bar[H[i1]] H[i1]) FS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]2,\[Mu]3,I2] FS[W,\[Mu]3,\[Mu]1,I3] +
		cW3H2$2[] fSU2L[I1,I2,I3] (Bar[H[i1]] H[i1]) FS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]2,\[Mu]3,I2] dFS[W,\[Mu]3,\[Mu]1,I3] +

		cW2BH2$1[] fSU2L[I1,I2,I3] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) FS[B,\[Mu]1,\[Mu]2] FS[W,\[Mu]2,\[Mu]3,I2] FS[W,\[Mu]3,\[Mu]1,I3] +
		cW2BH2$2[] fSU2L[I1,I2,I3] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) (dFS[B,\[Mu]1,\[Mu]2] FS[W,\[Mu]2,\[Mu]3,I2] FS[W,\[Mu]1,\[Mu]3,I3] + FS[B,\[Mu]1,\[Mu]2] FS[W,\[Mu]2,\[Mu]3,I2] dFS[W,\[Mu]1,\[Mu]3,I3])
		;

	L6 =
		cG2H4$1[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) FS[G,\[Mu]1,\[Mu]2,A1] FS[G,\[Mu]1,\[Mu]2,A1] +
		cG2H4$2[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) dFS[G,\[Mu]1,\[Mu]2,A1] FS[G,\[Mu]1,\[Mu]2,A1] +

		cW2H4$1[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) FS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]1,\[Mu]2,I1] +
		cW2H4$2[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) dFS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]1,\[Mu]2,I1] +
		cW2H4$3[] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) (Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]) FS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]1,\[Mu]2,I2] +
		cW2H4$4[] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) (Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]) dFS[W,\[Mu]1,\[Mu]2,I1] FS[W,\[Mu]1,\[Mu]2,I2] +

		cWBH4$1[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] FS[B,\[Mu]1,\[Mu]2] +
		cWBH4$2[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) dFS[W,\[Mu]1,\[Mu]2,I1] FS[B,\[Mu]1,\[Mu]2] +

		cB2H4$1[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) FS[B,\[Mu]1,\[Mu]2] FS[B,\[Mu]1,\[Mu]2] +
		cB2H4$2[] (Bar[H[i1]] H[i1]) (Bar[H[i2]] H[i2]) dFS[B,\[Mu]1,\[Mu]2] FS[B,\[Mu]1,\[Mu]2]
		;

	L7 =
		cG2H2D2$1[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) FS[G,\[Mu]1,\[Mu]3,A1] FS[G,\[Mu]2,\[Mu]3,A1] +
		cG2H2D2$2[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[G,\[Mu]2,\[Mu]3,A1] FS[G,\[Mu]2,\[Mu]3,A1] +
		cG2H2D2$3[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[G,\[Mu]2,\[Mu]3,A1] dFS[G,\[Mu]2,\[Mu]3,A1] +

		cW2H2D2$1[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) FS[W,\[Mu]1,\[Mu]3,I1] FS[W,\[Mu]2,\[Mu]3,I1] +
		cW2H2D2$2[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[W,\[Mu]2,\[Mu]3,I1] FS[W,\[Mu]2,\[Mu]3,I1] +
		cW2H2D2$3[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[W,\[Mu]2,\[Mu]3,I1] dFS[W,\[Mu]2,\[Mu]3,I1] +
		cW2H2D2$4[] I fSU2L[I1,I2,I3] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) FS[W,\[Mu]1,\[Mu]3,I2] FS[W,\[Mu]2,\[Mu]3,I3] +
		cW2H2D2$5[] fSU2L[I1,I2,I3] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[W,\[Mu]1,\[Mu]3,I2] dFS[W,\[Mu]2,\[Mu]3,I3] - dFS[W,\[Mu]1,\[Mu]3,I2] FS[W,\[Mu]2,\[Mu]3,I3]) +
		cW2H2D2$6[] I fSU2L[I1,I2,I3] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[W,\[Mu]1,\[Mu]3,I2] dFS[W,\[Mu]2,\[Mu]3,I3] + dFS[W,\[Mu]1,\[Mu]3,I2] FS[W,\[Mu]2,\[Mu]3,I3]) +

		cWBH2D2$1[] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,H[i2]]) FS[B,\[Mu]2,\[Mu]3] FS[W,\[Mu]2,\[Mu]3,I1] +
		cWBH2D2$2[] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,H[i2]]) FS[B,\[Mu]2,\[Mu]3] dFS[W,\[Mu]2,\[Mu]3,I1] +
		cWBH2D2$3[] I (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[B,\[Mu]1,\[Mu]3] FS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3] FS[W,\[Mu]1,\[Mu]3,I1]) +
		cWBH2D2$4[] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[B,\[Mu]1,\[Mu]3] FS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3] FS[W,\[Mu]1,\[Mu]3,I1]) +
		cWBH2D2$5[] I (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[B,\[Mu]1,\[Mu]3] dFS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3] dFS[W,\[Mu]1,\[Mu]3,I1]) +
		cWBH2D2$6[] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (FS[B,\[Mu]1,\[Mu]3] dFS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3] dFS[W,\[Mu]1,\[Mu]3,I1]) +

		cB2H2D2$1[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) FS[B,\[Mu]1,\[Mu]3] FS[B,\[Mu]2,\[Mu]3] +
		cB2H2D2$2[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[B,\[Mu]2,\[Mu]3] FS[B,\[Mu]2,\[Mu]3] +
		cB2H2D2$3[] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) FS[B,\[Mu]2,\[Mu]3] dFS[B,\[Mu]2,\[Mu]3]
		;

	(* Factors of \[ImaginaryI] have been added below (w.r.t. Murphy's paper) to make the Wilson coefficient real *)
	L8 =
		cWH4D2$1[] I (Bar[H[i1]] H[i1]) (CD[\[Mu]1,Bar[H[i2]]] \[Tau]SU2L[I1,i2,i3] CD[\[Mu]2,H[i3]]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cWH4D2$2[] I (Bar[H[i1]] H[i1]) (CD[\[Mu]1,Bar[H[i2]]] \[Tau]SU2L[I1,i2,i3] CD[\[Mu]2,H[i3]]) dFS[W,\[Mu]1,\[Mu]2,I1] +
		cWH4D2$3[] I fSU2L[I1,I2,I3] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) (CD[\[Mu]1,Bar[H[i3]]] \[Tau]SU2L[I2,i3,i4] CD[\[Mu]2,H[i4]]) FS[W,\[Mu]1,\[Mu]2,I3] +
		cWH4D2$4[] I fSU2L[I1,I2,I3] (Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]) (CD[\[Mu]1,Bar[H[i3]]] \[Tau]SU2L[I2,i3,i4] CD[\[Mu]2,H[i4]]) dFS[W,\[Mu]1,\[Mu]2,I3] +

		cBH4D2$1[] I (Bar[H[i1]] H[i1]) (CD[\[Mu]1,Bar[H[i2]]] CD[\[Mu]2,H[i2]]) FS[B,\[Mu]1,\[Mu]2] +
		cBH4D2$2[] I (Bar[H[i1]] H[i1]) (CD[\[Mu]1,Bar[H[i2]]] CD[\[Mu]2,H[i2]]) dFS[B,\[Mu]1,\[Mu]2]
		;

	(* 2-fermion sector *)

	L9 = PlusHc[
		cleG2H$1[p,r] (Bar[l[i1,p]]**e[r]) H[i1] FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +
		cleG2H$2[p,r] (Bar[l[i1,p]]**e[r]) H[i1] dFS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +

		cleW2H$1[p,r] (Bar[l[i1,p]]**e[r]) H[i1] FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cleW2H$2[p,r] (Bar[l[i1,p]]**e[r]) H[i1] dFS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cleW2H$3[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3] +

		cquG2H$1[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +
		cquG2H$2[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] dFS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +
		cquG2H$3[p,r] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] FS[G,\[Mu]1,\[Mu]2,A2]FS[G,\[Mu]1,\[Mu]2,A3] +
		cquG2H$4[p,r] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] dFS[G,\[Mu]1,\[Mu]2,A2]FS[G,\[Mu]1,\[Mu]2,A3] +
		cquG2H$5[p,r] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +

		cquGWH$1[p,r] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] dH[i2] FS[G,\[Mu]1,\[Mu]2,A1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cquGWH$2[p,r] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] dH[i2] dFS[G,\[Mu]1,\[Mu]2,A1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cquGWH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] dH[i2] FS[G,\[Mu]1,\[Mu]3,A1]FS[W,\[Mu]2,\[Mu]3,I1] +

		cquGBH$1[p,r] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] FS[G,\[Mu]1,\[Mu]2,A1]FS[B,\[Mu]1,\[Mu]2] +
		cquGBH$2[p,r] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] dFS[G,\[Mu]1,\[Mu]2,A1]FS[B,\[Mu]1,\[Mu]2] +
		cquGBH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] FS[G,\[Mu]1,\[Mu]3,A1]FS[B,\[Mu]2,\[Mu]3] +

		cquW2H$1[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cquW2H$2[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] dFS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cquW2H$3[p,r] fSU2L[I1,I2,I3](Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] dH[i2] FS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3] +

		cquWBH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] dH[i2] FS[W,\[Mu]1,\[Mu]3,I1]FS[B,\[Mu]2,\[Mu]3] +
		cquWBH$1[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] dH[i2] FS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +
		cquWBH$2[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] dH[i2] dFS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +

		cquB2H$1[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2] +
		cquB2H$2[p,r] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] dFS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2] +

		cleWBH$1[p,r] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +
		cleWBH$2[p,r] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] H[i2] dFS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +
		cleWBH$3[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]3,I1]FS[B,\[Mu]2,\[Mu]3] +

		cleB2H$1[p,r] (Bar[l[i1,p]]**e[r]) H[i1] FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2] +
		cleB2H$2[p,r] (Bar[l[i1,p]]**e[r]) H[i1] dFS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2] +

		cqdG2H$1[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] FS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +
		cqdG2H$2[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] dFS[G,\[Mu]1,\[Mu]2,A1]FS[G,\[Mu]1,\[Mu]2,A1] +
		cqdG2H$3[p,r] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] FS[G,\[Mu]1,\[Mu]2,A2]FS[G,\[Mu]1,\[Mu]2,A3] +
		cqdG2H$4[p,r] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] dFS[G,\[Mu]1,\[Mu]2,A2]FS[G,\[Mu]1,\[Mu]2,A3] +
		cqdG2H$5[p,r] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +

		cqdGWH$1[p,r] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] H[i2] FS[G,\[Mu]1,\[Mu]2,A1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cqdGWH$2[p,r] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] H[i2] dFS[G,\[Mu]1,\[Mu]2,A1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cqdGWH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] H[i2] FS[G,\[Mu]1,\[Mu]3,A1]FS[W,\[Mu]2,\[Mu]3,I1] +

		cqdGBH$1[p,r] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] FS[G,\[Mu]1,\[Mu]2,A1]FS[B,\[Mu]1,\[Mu]2] +
		cqdGBH$2[p,r] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] dFS[G,\[Mu]1,\[Mu]2,A1]FS[B,\[Mu]1,\[Mu]2] +
		cqdGBH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] FS[G,\[Mu]1,\[Mu]3,A1]FS[B,\[Mu]2,\[Mu]3] +

		cqdW2H$1[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] FS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cqdW2H$2[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] dFS[W,\[Mu]1,\[Mu]2,I1]FS[W,\[Mu]1,\[Mu]2,I1] +
		cqdW2H$3[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3] +

		cqdWBH$1[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +
		cqdWBH$2[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] H[i2] dFS[W,\[Mu]1,\[Mu]2,I1]FS[B,\[Mu]1,\[Mu]2] +
		cqdWBH$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] H[i2] FS[W,\[Mu]1,\[Mu]3,I1]FS[B,\[Mu]2,\[Mu]3] +

		cqdB2H$1[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] FS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2] +
		cqdB2H$2[p,r] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] dFS[B,\[Mu]1,\[Mu]2]FS[B,\[Mu]1,\[Mu]2]
	];

	L10 = PlusHc[
		cleWH3$1[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] H[i2] (Bar[H[i3]]H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cleWH3$2[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) H[i1] (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +

		cleBH3[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) H[i1] (Bar[H[i2]]H[i2]) FS[B,\[Mu]1,\[Mu]2] +

		cquGH3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] dH[i1] (Bar[H[i2]]H[i2]) FS[G,\[Mu]1,\[Mu]2,A1] +

		cquWH3$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] dH[i2] (Bar[H[i3]]H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cquWH3$2[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) dH[i1] (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +

		cquBH3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) dH[i1] (Bar[H[i2]]H[i2]) FS[B,\[Mu]1,\[Mu]2] +

		cqdGH3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] (Bar[H[i2]]H[i2]) FS[G,\[Mu]1,\[Mu]2,A1] +

		cqdWH3$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] H[i2] (Bar[H[i3]]H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cqdWH3$2[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) H[i1] (Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]) FS[W,\[Mu]1,\[Mu]2,I1] +

		cqdBH3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) H[i1] (Bar[H[i2]]H[i2]) FS[B,\[Mu]1,\[Mu]2]
	];

	L11 = PlusHc[ (* Contrary to Murphys's paper, these terms are non-Hermitian and the Hermitian conjugate must be added *)
		cl2H2D3$1[p,r] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,l[i1,r]]) (symCD[{\[Mu]1,\[Mu]2},Bar[H[i2]]] H[i2]) +
		cl2H2D3$2[p,r] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,l[i1,r]]) (Bar[H[i2]] symCD[{\[Mu]1,\[Mu]2},H[i2]]) +
		cl2H2D3$3[p,r] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,l[i2,r]]) \[Tau]SU2L[I1,i1,i2] (symCD[{\[Mu]1,\[Mu]2},Bar[H[i3]]] \[Tau]SU2L[I1,i3,i4] H[i4]) +
		cl2H2D3$4[p,r] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,l[i2,r]]) \[Tau]SU2L[I1,i1,i2] (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] symCD[{\[Mu]1,\[Mu]2},H[i4]]) +

		ce2H2D3$1[p,r]  (Bar[e[p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,e[r]]) (symCD[{\[Mu]1,\[Mu]2},Bar[H[i1]]]H[i1]) +
		ce2H2D3$2[p,r] I (Bar[e[p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,e[r]]) (Bar[H[i1]] symCD[{\[Mu]1,\[Mu]2},H[i1]]) +

		cq2H2D3$1[p,r] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,q[a1,i1,r]]) (symCD[{\[Mu]1,\[Mu]2},Bar[H[i2]]]H[i2]) +
		cq2H2D3$2[p,r] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,q[a1,i1,r]]) (Bar[H[i2]] symCD[{\[Mu]1,\[Mu]2},H[i2]]) +
		cq2H2D3$3[p,r] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,q[a1,i2,r]]) \[Tau]SU2L[I1,i1,i2] (symCD[{\[Mu]1,\[Mu]2},Bar[H[i3]]] \[Tau]SU2L[I1,i3,i4] H[i4]) +
		cq2H2D3$4[p,r] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,q[a1,i2,r]]) \[Tau]SU2L[I1,i1,i2] (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] symCD[{\[Mu]1,\[Mu]2},H[i4]]) +

		cu2H2D3$1[p,r] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,u[a1,r]]) (symCD[{\[Mu]1,\[Mu]2},Bar[H[i1]]]H[i1]) +
		cu2H2D3$2[p,r] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,u[a1,r]]) (Bar[H[i1]] symCD[{\[Mu]1,\[Mu]2},H[i1]]) +

		cd2H2D3$1[p,r] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,d[a1,r]]) (symCD[{\[Mu]1,\[Mu]2},Bar[H[i1]]]H[i1]) +
		cd2H2D3$2[p,r] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,d[a1,r]]) (Bar[H[i1]] symCD[{\[Mu]1,\[Mu]2},H[i1]]) +

		cudH2D3[p,r] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**CD[\[Mu]2,d[a1,r]]) (Bar[dH[i1]] symCD[{\[Mu]1,\[Mu]2},H[i1]])
	];

	L12 = PlusHc[
		cleH5[p,r] (Bar[H[i1]]H[i1]) (Bar[H[i2]]H[i2]) (Bar[l[i3,p]]**e[r]) H[i3] +
		cquH5[p,r] (Bar[H[i1]]H[i1]) (Bar[H[i2]]H[i2]) (Bar[q[a1,i3,p]]**u[a1,r]) dH[i3] +
		cqdH5[p,r] (Bar[H[i1]]H[i1]) (Bar[H[i2]]H[i2]) (Bar[q[a1,i3,p]]**d[a1,r]) H[i3]
	];

	L13 =
		cl2H4D$1[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]],H[i2]] (Bar[H[i3]]H[i3]) +
		cl2H4D$2[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (HermitianCD[\[Mu]1, Bar[H[i3]],H[i4]] \[Tau]SU2L[I1,i3,i4] (Bar[H[i5]]H[i5])+HermitianCD[\[Mu]1, Bar[H[i3]],H[i3]] \[Tau]SU2L[I1,i4,i5] (Bar[H[i4]]H[i5])) +
		cl2H4D$3[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]],H[i4]] \[Tau]SU2L[I2,i3,i4] (Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6]) +
		cl2H4D$4[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]) CD[\[Mu]1,(Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6])] +

		ce2H4D[p,r] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) HermitianCD[\[Mu]1, Bar[H[i1]],H[i1]] (Bar[H[i2]]H[i2]) +

		cq2H4D$1[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]],H[i2]] (Bar[H[i3]]H[i3]) +
		cq2H4D$2[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (HermitianCD[\[Mu]1, Bar[H[i3]],H[i4]] \[Tau]SU2L[I1,i3,i4] (Bar[H[i5]]H[i5])+HermitianCD[\[Mu]1, Bar[H[i3]],H[i3]] \[Tau]SU2L[I1,i4,i5] (Bar[H[i4]]H[i5])) +
		cq2H4D$3[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]],H[i4]] \[Tau]SU2L[I2,i3,i4] (Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6]) +
		cq2H4D$4[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]) CD[\[Mu]1,(Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6])] +

		cu2H4D[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]],H[i1]] (Bar[H[i2]]H[i2]) +

		cd2H4D[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]],H[i1]] (Bar[H[i2]]H[i2]) +

		PlusHc[cudH4D[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) HermitianCD[\[Mu]1, Bar[dH[i1]],H[i1]] (Bar[H[i2]]H[i2])]
	;

	L14 =
		cq2G2D$1[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i1,r]] FS[G,\[Mu]1,\[Mu]3,A1]FS[G,\[Mu]2,\[Mu]3,A1] +
		cq2G2D$2[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +
		cq2G2D$3[p,r] dSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +

		cq2W2D$1[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i1,r]] FS[W,\[Mu]1,\[Mu]3,I1]FS[W,\[Mu]2,\[Mu]3,I1] +
		cq2W2D$2[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3] +

		cq2B2D[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i1,r]] FS[B,\[Mu]1,\[Mu]3]FS[B,\[Mu]2,\[Mu]3] +

		cu2G2D$1[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a1,r]] FS[G,\[Mu]1,\[Mu]3,A1]FS[G,\[Mu]2,\[Mu]3,A1] +
		cu2G2D$2[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +
		cu2G2D$3[p,r] dSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +

		cu2W2D[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a1,r]] FS[W,\[Mu]1,\[Mu]3,I1]FS[W,\[Mu]2,\[Mu]3,I1] +

		cu2B2D[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a1,r]] FS[B,\[Mu]1,\[Mu]3]FS[B,\[Mu]2,\[Mu]3] +

		cd2G2D$1[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a1,r]] FS[G,\[Mu]1,\[Mu]3,A1]FS[G,\[Mu]2,\[Mu]3,A1] +
		cd2G2D$2[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +
		cd2G2D$3[p,r] dSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3] +

		cd2W2D[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a1,r]] FS[W,\[Mu]1,\[Mu]3,I1]FS[W,\[Mu]2,\[Mu]3,I1] +

		cd2B2D[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a1,r]] FS[B,\[Mu]1,\[Mu]3]FS[B,\[Mu]2,\[Mu]3] +

		cq2G2D$4[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] - dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +
		cq2G2D$5[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] + dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +

		cq2GWD$1[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i2,r]] TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] (FS[G,\[Mu]1,\[Mu]3,A1]FS[W,\[Mu]2,\[Mu]3,I1] - FS[G,\[Mu]2,\[Mu]3,A1]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2GWD$2[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i2,r]] TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] (FS[G,\[Mu]1,\[Mu]3,A1]FS[W,\[Mu]2,\[Mu]3,I1] + FS[G,\[Mu]2,\[Mu]3,A1]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2GWD$3[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i2,r]] TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] (FS[G,\[Mu]1,\[Mu]3,A1]dFS[W,\[Mu]2,\[Mu]3,I1] - FS[G,\[Mu]2,\[Mu]3,A1]dFS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2GWD$4[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i2,r]] TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] (FS[G,\[Mu]1,\[Mu]3,A1]dFS[W,\[Mu]2,\[Mu]3,I1] + FS[G,\[Mu]2,\[Mu]3,A1]dFS[W,\[Mu]1,\[Mu]3,I1]) +

		cq2GBD$1[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cq2GBD$2[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cq2GBD$3[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +
		cq2GBD$4[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a2,i1,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +

		cq2W2D$3[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[W,\[Mu]1,\[Mu]3,I2]dFS[W,\[Mu]2,\[Mu]3,I3] - dFS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3]) +
		cq2W2D$4[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[W,\[Mu]1,\[Mu]3,I2]dFS[W,\[Mu]2,\[Mu]3,I3] + dFS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3]) +

		cq2WBD$1[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]FS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2WBD$2[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]FS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2WBD$3[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]dFS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3]dFS[W,\[Mu]1,\[Mu]3,I1]) +
		cq2WBD$4[p,r] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]],\[Gamma][\[Mu]1]**q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]dFS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3]dFS[W,\[Mu]1,\[Mu]3,I1]) +

		cu2G2D$4[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] - dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +
		cu2G2D$5[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] + dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +

		cu2GBD$1[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cu2GBD$2[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cu2GBD$3[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +
		cu2GBD$4[p,r] HermitianCD[\[Mu]2,Bar[u[a1,p]],\[Gamma][\[Mu]1]**u[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +

		cd2G2D$4[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] - dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +
		cd2G2D$5[p,r] fSU3c[A1,A2,A3] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[G,\[Mu]1,\[Mu]3,A2]dFS[G,\[Mu]2,\[Mu]3,A3] + dFS[G,\[Mu]1,\[Mu]3,A2]FS[G,\[Mu]2,\[Mu]3,A3]) +

		cd2GBD$1[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cd2GBD$2[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]FS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]FS[G,\[Mu]1,\[Mu]3,A1]) +
		cd2GBD$3[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] - FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +
		cd2GBD$4[p,r] HermitianCD[\[Mu]2,Bar[d[a1,p]],\[Gamma][\[Mu]1]**d[a2,r]] TSU3c[A1,a1,a2] (FS[B,\[Mu]1,\[Mu]3]dFS[G,\[Mu]2,\[Mu]3,A1] + FS[B,\[Mu]2,\[Mu]3]dFS[G,\[Mu]1,\[Mu]3,A1]) +

		cl2G2D[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i1,r]] FS[G,\[Mu]1,\[Mu]3,A1]FS[G,\[Mu]2,\[Mu]3,A1] +

		cl2W2D$1[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i1,r]] FS[W,\[Mu]1,\[Mu]3,I1]FS[W,\[Mu]2,\[Mu]3,I1] +
		cl2W2D$2[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3] +

		cl2B2D[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i1,r]] FS[B,\[Mu]1,\[Mu]3]FS[B,\[Mu]2,\[Mu]3] +

		ce2G2D[p,r] HermitianCD[\[Mu]2,Bar[e[p]],\[Gamma][\[Mu]1]**e[r]] FS[G,\[Mu]1,\[Mu]3,A1]FS[G,\[Mu]2,\[Mu]3,A1] +

		ce2W2D[p,r] HermitianCD[\[Mu]2,Bar[e[p]],\[Gamma][\[Mu]1]**e[r]] FS[W,\[Mu]1,\[Mu]3,I1]FS[W,\[Mu]2,\[Mu]3,I1] +

		ce2B2D[p,r] HermitianCD[\[Mu]2,Bar[e[p]],\[Gamma][\[Mu]1]**e[r]] FS[B,\[Mu]1,\[Mu]3]FS[B,\[Mu]2,\[Mu]3] +

		cl2W2D$3[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[W,\[Mu]1,\[Mu]3,I2]dFS[W,\[Mu]2,\[Mu]3,I3] - dFS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3]) +
		cl2W2D$4[p,r] fSU2L[I1,I2,I3] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[W,\[Mu]1,\[Mu]3,I2]dFS[W,\[Mu]2,\[Mu]3,I3] + dFS[W,\[Mu]1,\[Mu]3,I2]FS[W,\[Mu]2,\[Mu]3,I3]) +

		cl2WBD$1[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]FS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cl2WBD$2[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]FS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3]FS[W,\[Mu]1,\[Mu]3,I1]) +
		cl2WBD$3[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]dFS[W,\[Mu]2,\[Mu]3,I1] - FS[B,\[Mu]2,\[Mu]3]dFS[W,\[Mu]1,\[Mu]3,I1]) +
		cl2WBD$4[p,r] HermitianCD[\[Mu]2,Bar[l[i1,p]],\[Gamma][\[Mu]1]**l[i2,r]] \[Tau]SU2L[I1,i1,i2] (FS[B,\[Mu]1,\[Mu]3]dFS[W,\[Mu]2,\[Mu]3,I1] + FS[B,\[Mu]2,\[Mu]3]dFS[W,\[Mu]1,\[Mu]3,I1])
	;

	(* The operators containing a hermitian derivative have an extra factor of \[ImaginaryI] w.r.t Murphy's basis to make the Wilson coefficients real *)
	L15 =
		ce2WH2D$1[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] FS[W,\[Mu]1,\[Mu]2,I1] +
		ce2WH2D$2[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		ce2WH2D$3[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]2,I1] +
		ce2WH2D$4[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] dFS[W,\[Mu]1,\[Mu]2,I1] +

		ce2BH2D$1[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		ce2BH2D$2[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] dFS[B,\[Mu]1,\[Mu]2] +
		ce2BH2D$3[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		ce2BH2D$4[p,r] (Bar[e[p]]**\[Gamma][\[Mu]2]**e[r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] dFS[B,\[Mu]1,\[Mu]2] +

		cu2GH2D$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i1]] H[i1]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cu2GH2D$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i1]] H[i1]] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cu2GH2D$3[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cu2GH2D$4[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cu2WH2D$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cu2WH2D$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cu2WH2D$3[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]2,I1] +
		cu2WH2D$4[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] dFS[W,\[Mu]1,\[Mu]2,I1] +

		cu2BH2D$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		cu2BH2D$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] dFS[B,\[Mu]1,\[Mu]2] +
		cu2BH2D$3[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		cu2BH2D$4[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**u[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] dFS[B,\[Mu]1,\[Mu]2] +

		cd2GH2D$1[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i1]] H[i1]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cd2GH2D$2[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i1]] H[i1]] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cd2GH2D$3[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cd2GH2D$4[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cd2WH2D$1[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cd2WH2D$2[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) CD[\[Mu]1,Bar[H[i1]] \[Tau]SU2L[I1,i1,i2] H[i2]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cd2WH2D$3[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]2,I1] +
		cd2WH2D$4[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] dFS[W,\[Mu]1,\[Mu]2,I1] +

		cd2BH2D$1[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		cd2BH2D$2[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) CD[\[Mu]1,Bar[H[i1]] H[i1]] dFS[B,\[Mu]1,\[Mu]2] +
		cd2BH2D$3[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		cd2BH2D$4[p,r] (Bar[d[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[H[i1]], H[i1]] dFS[B,\[Mu]1,\[Mu]2] +

		PlusHc[
		cudGH2$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[dH[i1]], H[i1]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cudGH2$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[dH[i1]], H[i1]] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cudWH2$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[dH[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]2,I1] +
		cudWH2$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[dH[i1]], H[i2]] \[Tau]SU2L[I1,i1,i2] dFS[W,\[Mu]1,\[Mu]2,I1] +

		cudBH2$1[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[dH[i1]], H[i1]] FS[B,\[Mu]1,\[Mu]2] +
		cudBH2$2[p,r] (Bar[u[a1,p]]**\[Gamma][\[Mu]2]**d[a1,r]) HermitianCD[\[Mu]1, Bar[dH[i1]], H[i1]] dFS[B,\[Mu]1,\[Mu]2]
		] +

		cl2WH2D$1[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) CD[\[Mu]1,Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$2[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) CD[\[Mu]1,Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$3[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i3]] \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$4[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i3]] \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$5[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$6[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$7[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$8[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2WH2D$9[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]] FS[W,\[Mu]1,\[Mu]2,I3] +
		cl2WH2D$10[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]] dFS[W,\[Mu]1,\[Mu]2,I3] +
		cl2WH2D$11[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I2,i3,i4] FS[W,\[Mu]1,\[Mu]2,I3] +
		cl2WH2D$12[p,r] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I2,i3,i4] dFS[W,\[Mu]1,\[Mu]2,I3] +

		cl2BH2D$1[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] FS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$2[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] dFS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$3[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] FS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$4[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] dFS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$5[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) CD[\[Mu]1,Bar[H[i2]] H[i2]] FS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$6[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) CD[\[Mu]1,Bar[H[i2]] H[i2]] dFS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$7[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] FS[B,\[Mu]1,\[Mu]2] +
		cl2BH2D$8[p,r] (Bar[l[i1,p]]**\[Gamma][\[Mu]2]**l[i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] dFS[B,\[Mu]1,\[Mu]2] +

		cq2GH2D$1[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$2[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$3[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$4[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i2,r]) TSU3c[A1,a1,a2] \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$5[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i1,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i2]] H[i2]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$6[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i1,r]) TSU3c[A1,a1,a2] CD[\[Mu]1,Bar[H[i2]] H[i2]] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$7[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i1,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2GH2D$8[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a2,i1,r]) TSU3c[A1,a1,a2] HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cq2WH2D$1[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) CD[\[Mu]1,Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$2[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) CD[\[Mu]1,Bar[H[i2]] \[Tau]SU2L[I1,i2,i3] H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$3[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i3]] \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$4[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i3]] \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$5[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$6[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$7[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i3]] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$8[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i3]] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2WH2D$9[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]] FS[W,\[Mu]1,\[Mu]2,I3] +
		cq2WH2D$10[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I2,i3,i4] H[i4]] dFS[W,\[Mu]1,\[Mu]2,I3] +
		cq2WH2D$11[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I2,i3,i4] FS[W,\[Mu]1,\[Mu]2,I3] +
		cq2WH2D$12[p,r] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I2,i3,i4] dFS[W,\[Mu]1,\[Mu]2,I3] +

		cq2BH2D$1[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] FS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$2[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]] dFS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$3[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] FS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$4[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] HermitianCD[\[Mu]1, Bar[H[i3]], H[i4]] \[Tau]SU2L[I1,i3,i4] dFS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$5[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) CD[\[Mu]1,Bar[H[i2]] H[i2]] FS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$6[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) CD[\[Mu]1,Bar[H[i2]] H[i2]] dFS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$7[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] FS[B,\[Mu]1,\[Mu]2] +
		cq2BH2D$8[p,r] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]2]**q[a1,i1,r]) HermitianCD[\[Mu]1, Bar[H[i2]], H[i2]] dFS[B,\[Mu]1,\[Mu]2]
	;

	L16 = PlusHc[
		cleWHD2$1[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,e[r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]] FS[W,\[Mu]3,\[Mu]1,I1] +
		cleWHD2$2[p,r] (Bar[l[i1,p]]**CD[\[Mu]3,e[r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]] dFS[W,\[Mu]3,\[Mu]2,I1] +
		cleWHD2$3[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]3,H[i2]] CD[\[Mu]3,FS[W,\[Mu]1,\[Mu]2,I1]] +

		cleBHD2$1[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,e[r]]) CD[\[Mu]2,H[i1]] FS[B,\[Mu]3,\[Mu]1] +
		cleBHD2$2[p,r] (Bar[l[i1,p]]**CD[\[Mu]3,e[r]]) CD[\[Mu]2,H[i1]] dFS[B,\[Mu]3,\[Mu]2] +
		cleBHD2$3[p,r] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) CD[\[Mu]3,H[i1]] CD[\[Mu]3,FS[B,\[Mu]1,\[Mu]2]] +

		cquGHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,u[a2,r]]) TSU3c[A1,a1,a2] CD[\[Mu]2,dH[i1]] FS[G,\[Mu]3,\[Mu]1,A1] +
		cquGHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,u[a2,r]]) TSU3c[A1,a1,a2] CD[\[Mu]2,dH[i1]] dFS[G,\[Mu]3,\[Mu]2,A1] +
		cquGHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]3,dH[i1]] CD[\[Mu]3,FS[G,\[Mu]1,\[Mu]2,A1]] +

		cquWHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,u[a1,r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,dH[i2]] FS[W,\[Mu]3,\[Mu]1,I1] +
		cquWHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,u[a1,r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,dH[i2]] dFS[W,\[Mu]3,\[Mu]2,I1] +
		cquWHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]3,dH[i2]] CD[\[Mu]3,FS[W,\[Mu]1,\[Mu]2,I1]] +

		cquBHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,u[a1,r]]) CD[\[Mu]2,dH[i1]] FS[B,\[Mu]3,\[Mu]1] +
		cquBHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,u[a1,r]]) CD[\[Mu]2,dH[i1]] dFS[B,\[Mu]3,\[Mu]2] +
		cquBHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) CD[\[Mu]3,dH[i1]] CD[\[Mu]3,FS[B,\[Mu]1,\[Mu]2]] +

		cqdGHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,d[a2,r]]) TSU3c[A1,a1,a2] CD[\[Mu]2,H[i1]] FS[G,\[Mu]3,\[Mu]1,A1] +
		cqdGHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,d[a2,r]]) TSU3c[A1,a1,a2] CD[\[Mu]2,H[i1]] dFS[G,\[Mu]3,\[Mu]2,A1] +
		cqdGHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,r]) TSU3c[A1,a1,a2] CD[\[Mu]3,H[i1]] CD[\[Mu]3,FS[G,\[Mu]1,\[Mu]2,A1]] +

		cqdWHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,d[a1,r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]] FS[W,\[Mu]3,\[Mu]1,I1] +
		cqdWHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,d[a1,r]]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]] dFS[W,\[Mu]3,\[Mu]2,I1] +
		cqdWHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) \[Tau]SU2L[I1,i1,i2] CD[\[Mu]3,H[i2]] CD[\[Mu]3,FS[W,\[Mu]1,\[Mu]2,I1]] +

		cqdBHD2$1[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**CD[\[Mu]3,d[a1,r]]) CD[\[Mu]2,H[i1]] FS[B,\[Mu]3,\[Mu]1] +
		cqdBHD2$2[p,r] (Bar[q[a1,i1,p]]**CD[\[Mu]3,d[a1,r]]) CD[\[Mu]2,H[i1]] dFS[B,\[Mu]3,\[Mu]2] +
		cqdBHD2$3[p,r] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) CD[\[Mu]3,H[i1]] CD[\[Mu]3,FS[B,\[Mu]1,\[Mu]2]]
	];

	L17 = PlusHc[
		cleH3D2$1[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) (Bar[l[i2,p]]**e[r]) H[i2] +
		cleH3D2$2[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,H[i2]]) (Bar[l[i3,p]]**e[r]) \[Tau]SU2L[I1,i3,i4] H[i4] +
		cleH3D2$3[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) (Bar[l[i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) H[i2] +
		cleH3D2$4[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (Bar[l[i3,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i3,i4] H[i4] +
		cleH3D2$5[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[l[i2,p]]**e[r]) CD[\[Mu]1,H[i2]] +
		cleH3D2$6[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[l[i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) CD[\[Mu]2,H[i2]] +

		cquH3D2$1[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**u[a1,r]) dH[i2] +
		cquH3D2$2[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,H[i2]]) (Bar[q[a1,i3,p]]**u[a1,r]) \[Tau]SU2L[I1,i3,i4] dH[i4] +
		cquH3D2$3[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) (Bar[q[a1,i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) dH[i2] +
		cquH3D2$4[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (Bar[q[a1,i3,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i3,i4] dH[i4] +
		cquH3D2$5[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**u[a1,r]) CD[\[Mu]1,dH[i2]] +
		cquH3D2$6[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) CD[\[Mu]2,dH[i2]] +

		cqdH3D2$1[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**d[a1,r]) H[i2] +
		cqdH3D2$2[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]1,H[i2]]) (Bar[q[a1,i3,p]]**d[a1,r]) \[Tau]SU2L[I1,i3,i4] H[i4] +
		cqdH3D2$3[p,r] (CD[\[Mu]1,Bar[H[i1]]] CD[\[Mu]2,H[i1]]) (Bar[q[a1,i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) H[i2] +
		cqdH3D2$4[p,r] (CD[\[Mu]1,Bar[H[i1]]] \[Tau]SU2L[I1,i1,i2] CD[\[Mu]2,H[i2]]) (Bar[q[a1,i3,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) \[Tau]SU2L[I1,i3,i4] H[i4] +
		cqdH3D2$5[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**d[a1,r]) CD[\[Mu]1,H[i2]] +
		cqdH3D2$6[p,r] (Bar[H[i1]] CD[\[Mu]1,H[i1]]) (Bar[q[a1,i2,p]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,r]) CD[\[Mu]2,H[i2]]
	];

	L18 = PlusHc[
		(* (LR)(RL)HH *)
		cleqdH2$1[p,r,s,t] (Bar[l[i1,p]]**e[r]) (Bar[d[a1,s]]**q[a1,i1,t]) (Bar[H[i2]] H[i2]) +
		cleqdH2$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**q[a1,i2,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cl2udH2[p,r,s,t] (Bar[l[i1,p]]**d[a1,r]) H[i1] (Bar[dH[i2]]Bar[u[a1,s]]**l[i2,t]) +

		clequH2$5[p,r,s,t] (Bar[l[i1,p]]**e[r]) H[i1] (Bar[dH[i2]]Bar[u[a1,s]]**q[a1,i2,t]) +

		cq2udH2$5[p,r,s,t] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] (Bar[dH[i2]]Bar[u[a2,s]]**q[a2,i2,t]) +
		cq2udH2$6[p,r,s,t] (Bar[q[a1,i1,p]]**d[a2,r]) TSU3c[A1,a1,a2] H[i1] (Bar[dH[i2]]Bar[u[a3,s]]**q[a4,i2,t]) TSU3c[A1,a3,a4]
	] +

	If[TrueQ@ BLViolation,
		PlusHc[(* B-violating *) (* The last 2 operators, Q_lqd2H2 & Q_eq2dH2, are missing a charge conjugation matrix, CC, in Murphy's paper *)
			clqudH2$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i2,t]) (Bar[H[i3]] H[i3]) +
			clqudH2$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
	
			ceq2uH2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[u[a3,s]]**CC**e[t]) (Bar[H[i3]] H[i2]) +
	
			clq3H2$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**l[i2,t]) (Bar[H[i5]] H[i5]) +
			clq3H2$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] Bar@\[Epsilon]SU2L[i4,i5] (Transp[q[a1,i1,p]]**CC**q[a2,i4,r]) (Transp[q[a3,i5,s]]**CC**l[i3,t]) (Bar[H[i6]] \[Tau]SU2L[I1,i6,i7] H[i7]) +
			clq3H2$3[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] \[Tau]SU2L[I1,i4,i5] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i5,s]]**CC**l[i2,t]) (Bar[H[i6]] \[Tau]SU2L[I1,i6,i7] H[i7]) +
	
			ceu2dH2[p,r,s,t]  Bar@\[Epsilon]SU3c[a1,a2,a3] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[u[a3,s]]**CC**e[t]) (Bar[H[i1]] H[i1]) +
	
			clqu2H2[p,r,s,t]  Bar@\[Epsilon]SU3c[a1,a2,a3]  Bar@\[Epsilon]SU2L[i1,i2]  Bar@\[Epsilon]SU2L[i3,i4] (Transp[l[i1,p]]**CC**q[a1,i3,r]) (Transp[u[a2,s]]**CC**u[a3,t]) (dH[i2] dH[i4]) +
	
			clqd2H2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[l[i1,p]]**CC**q[a1,i3,r]) (Transp[d[a2,s]]**CC**d[a3,t]) (H[i2] H[i4]) +
	
			ceq2dH2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[e[p]]**CC**d[a1,r]) (Transp[q[a2,i1,s]]**CC**q[a3,i3,t]) (H[i2] H[i4])
		]
	,
		0
	] +

		(* (LL)(LL)HH *)
		cl4H2$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]1]**l[i2,t]) (Bar[H[i3]] H[i3]) +
		cl4H2$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]1]**l[i3,t]) \[Tau]SU2L[I1,i2,i3] (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +

		cq4H2$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1]**q[a2,i2,t]) (Bar[H[i3]] H[i3]) +
		cq4H2$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1]**q[a2,i3,t]) \[Tau]SU2L[I1,i2,i3] (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
		cq4H2$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i4,t]) (Bar[H[i5]] H[i5]) +

		cl2q2H2$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1]**q[a2,i2,t]) (Bar[H[i3]] H[i3]) +
		cl2q2H2$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i3,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
		cl2q2H2$3[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i4,t]) (Bar[H[i5]] H[i5]) +
		cl2q2H2$4[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) \[Tau]SU2L[I1,i2,i3] (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1]**q[a2,i3,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
		cl2q2H2$5[p,r,s,t] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I2,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i4,t]) (Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6]) +

		cq4H2$5[p,r,s,t] fSU2L[I1,I2,I3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I2,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i4,t]) (Bar[H[i5]] \[Tau]SU2L[I3,i5,i6] H[i6]) +

		(* (RR)(RR)HH *)
		ce4H2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[e[s]]**\[Gamma][\[Mu]1]**e[t]) (Bar[H[i1]] H[i1]) +

		cu4H2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]1]**u[a2,t]) (Bar[H[i1]] H[i1]) +

		cd4H2[p,r,s,t] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t]) (Bar[H[i1]] H[i1]) +

		ce2u2H2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]1]**u[a1,t]) (Bar[H[i1]] H[i1]) +

		ce2d2H2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]1]**d[a1,t]) (Bar[H[i1]] H[i1]) +

		cu2d2H2$1[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t]) (Bar[H[i1]] H[i1]) +
		cu2d2H2$2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]1]**d[a4,t]) (Bar[H[i1]] H[i1]) +

		(* (LL)(RR)HH *)
		cl2e2H2$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]1]**e[t]) (Bar[H[i2]] H[i2]) +
		cl2e2H2$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]1]**e[t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cl2u2H2$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]1]**u[a1,t]) (Bar[H[i2]] H[i2]) +
		cl2u2H2$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a1,s]]**\[Gamma][\[Mu]1]**u[a1,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cl2d2H2$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]1]**d[a1,t]) (Bar[H[i2]] H[i2]) +
		cl2d2H2$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**\[Gamma][\[Mu]1]**d[a1,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cq2e2H2$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]1]**e[t]) (Bar[H[i2]] H[i2]) +
		cq2e2H2$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]1]**e[t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cq2u2H2$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]1]**u[a2,t]) (Bar[H[i2]] H[i2]) +
		cq2u2H2$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a2,s]]**\[Gamma][\[Mu]1]**u[a2,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +
		cq2u2H2$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]1]**u[a4,t]) (Bar[H[i2]] H[i2]) +
		cq2u2H2$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]1]**u[a4,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

		cq2d2H2$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t]) (Bar[H[i2]] H[i2]) +
		cq2d2H2$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +
		cq2d2H2$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]1]**d[a4,t]) (Bar[H[i2]] H[i2]) +
		cq2d2H2$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]1]**d[a4,t]) (Bar[H[i3]] \[Tau]SU2L[I1,i3,i4] H[i4]) +

	PlusHc[(* (LR)(LR)HH *)
		cq2udH2$1[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**d[a2,t]) (Bar[H[i3]] H[i3]) +
		cq2udH2$2[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**d[a2,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
		cq2udH2$3[p,r,s,t] (Bar[q[a1,i1,p]]**u[a2,r]) \[Epsilon]SU2L[i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**d[a4,t]) (Bar[H[i3]] H[i3]) +
		cq2udH2$4[p,r,s,t] (Bar[q[a1,i1,p]]**u[a2,r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i3,s]]**d[a4,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +

		clequH2$1[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**u[a2,t]) (Bar[H[i3]] H[i3]) +
		clequH2$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**u[a2,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +
		clequH2$3[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,t]) (Bar[H[i3]] H[i3]) +
		clequH2$4[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,t]) (Bar[H[i4]] \[Tau]SU2L[I1,i4,i5] H[i5]) +

		cl2e2H2$3[p,r,s,t] (Bar[l[i1,p]]**e[r]) H[i1] (Bar[l[i2,s]]**e[t]) H[i2] +

		cleqdH2$3[p,r,s,t] (Bar[l[i1,p]]**e[r]) H[i1] (Bar[q[a1,i2,s]]**d[a1,t]) H[i2] +
		cleqdH2$4[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) H[i1] (Bar[q[a1,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a1,t]) H[i2] +

		cq2u2H2$5[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) dH[i1] (Bar[q[a2,i2,s]]**u[a2,t]) dH[i2] +
		cq2u2H2$6[p,r,s,t] (Bar[q[a1,i1,p]]**u[a2,r]) dH[i1] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**u[a4,t]) dH[i2] +

		cq2d2H2$5[p,r,s,t] (Bar[q[a1,i1,p]]**d[a1,r]) H[i1] (Bar[q[a2,i2,s]]**d[a2,t]) H[i2] +
		cq2d2H2$6[p,r,s,t] (Bar[q[a1,i1,p]]**d[a2,r]) H[i1] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**d[a4,t]) H[i2]
	];

	L19 =
		(* (LL)(LL)X *)
		cl4W$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]2]**l[i3,t]) \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl4W$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]2]**l[i3,t]) \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +

		cq4G$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a3,i2,t]) TSU3c[A1,a2,a3] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq4G$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a3,i2,t]) TSU3c[A1,a2,a3] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq4G$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2]\[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]2]**q[a3,i4,t]) TSU3c[A1,a2,a3] FS[G,\[Mu]1,\[Mu]2,A1] +
		cq4G$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2]\[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]2]**q[a3,i4,t]) TSU3c[A1,a2,a3] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cq4W$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i3,t]) \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq4W$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i3,t]) \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq4W$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**\[Gamma][\[Mu]2]**q[a4,i3,t]) \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cq4W$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**\[Gamma][\[Mu]2]**q[a4,i3,t]) \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +

		cl2q2G$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i2,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		cl2q2G$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i2,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cl2q2G$3[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2]\[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a2,i4,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		cl2q2G$4[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2]\[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a2,i4,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cl2q2W$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a1,i3,t]) \[Tau]SU2L[I1,i2,i3] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2q2W$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a1,i3,t]) \[Tau]SU2L[I1,i2,i3] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2q2W$3[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i3,t]) \[Tau]SU2L[I1,i1,i2] FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2q2W$4[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i3,t]) \[Tau]SU2L[I1,i1,i2] dFS[W,\[Mu]1,\[Mu]2,I1] +
		cl2q2W$5[p,r,s,t] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I2,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i4,t]) FS[W,\[Mu]1,\[Mu]2,I3] +
		cl2q2W$6[p,r,s,t] fSU2L[I1,I2,I3] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I2,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i4,t]) dFS[W,\[Mu]1,\[Mu]2,I3] +

		cl2q2B$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a1,i2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cl2q2B$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**\[Gamma][\[Mu]2]**q[a1,i2,t]) dFS[B,\[Mu]1,\[Mu]2] +
		cl2q2B$3[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i4,t]) FS[B,\[Mu]1,\[Mu]2] +
		cl2q2B$4[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**\[Gamma][\[Mu]2]**q[a1,i4,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cl4B$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]2]**l[i2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cl4B$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**\[Gamma][\[Mu]2]**l[i2,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cq4B$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq4B$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**\[Gamma][\[Mu]2]**q[a2,i2,t]) dFS[B,\[Mu]1,\[Mu]2] +
		cq4B$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]2]**q[a2,i4,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq4B$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**\[Gamma][\[Mu]2]**q[a2,i4,t]) dFS[B,\[Mu]1,\[Mu]2] +

		(* (RR)(RR)X *)
		cu4G$1[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a3,t]) TSU3c[A1,a2,a3] FS[G,\[Mu]1,\[Mu]2,A1] +
		cu4G$2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a3,t]) TSU3c[A1,a2,a3] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cd4G$1[p,r,s,t] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a2,a3] FS[G,\[Mu]1,\[Mu]2,A1] +
		cd4G$2[p,r,s,t] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a2,a3] dFS[G,\[Mu]1,\[Mu]2,A1] +

		ce2u2G$1[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a2,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		ce2u2G$2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a2,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +

		ce2u2B$1[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		ce2u2B$2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) dFS[B,\[Mu]1,\[Mu]2] +

		ce2d2G$1[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a2,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		ce2d2G$2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a2,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +

		ce2d2B$1[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		ce2d2B$2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cu2d2G$1[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a2,a3] FS[G,\[Mu]1,\[Mu]2,A1] +
		cu2d2G$2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a2,a3] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cu2d2G$3[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		cu2d2G$4[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a3,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +
		cu2d2G$5[p,r,s,t] fSU3c[A1,A2,A3] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cu2d2G$6[p,r,s,t] fSU3c[A1,A2,A3] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +
		cu2d2G$7[p,r,s,t] dSU3c[A1,A2,A3] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cu2d2G$8[p,r,s,t] dSU3c[A1,A2,A3] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +

		cu2d2B$1[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cu2d2B$2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) dFS[B,\[Mu]1,\[Mu]2] +
		cu2d2B$3[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[B,\[Mu]1,\[Mu]2] +
		cu2d2B$4[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[B,\[Mu]1,\[Mu]2] +

		ce4B$1[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[B,\[Mu]1,\[Mu]2] +
		ce4B$2[p,r,s,t] (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[B,\[Mu]1,\[Mu]2] +

		cu4B$1[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cu4B$2[p,r,s,t] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cd4B$1[p,r,s,t] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cd4B$2[p,r,s,t] (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) dFS[B,\[Mu]1,\[Mu]2] +

		(* (LL)(RR)X *)
		cl2e2W$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2e2W$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cl2e2B$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[B,\[Mu]1,\[Mu]2] +
		cl2e2B$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[B,\[Mu]1,\[Mu]2] +

		cl2u2G$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a2,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		cl2u2G$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a2,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cl2u2W$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2u2W$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cl2u2B$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		cl2u2B$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[u[a1,s]]**\[Gamma][\[Mu]2]**u[a1,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cl2d2G$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a2,t]) TSU3c[A1,a1,a2] FS[G,\[Mu]1,\[Mu]2,A1] +
		cl2d2G$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a2,t]) TSU3c[A1,a1,a2] dFS[G,\[Mu]1,\[Mu]2,A1] +

		cl2d2W$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cl2d2W$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cl2d2B$1[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		cl2d2B$2[p,r,s,t] (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[d[a1,s]]**\[Gamma][\[Mu]2]**d[a1,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cq2e2G$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2e2G$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[G,\[Mu]1,\[Mu]2,A1] +

		cq2e2W$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2e2W$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cq2e2B$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2e2B$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[e[s]]**\[Gamma][\[Mu]2]**e[t]) dFS[B,\[Mu]1,\[Mu]2] +

		cq2u2G$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) TSU3c[A1,a2,a3] (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2u2G$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) TSU3c[A1,a2,a3] (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a3,t]) dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq2u2G$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2u2G$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2u2G$5[p,r,s,t] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cq2u2G$6[p,r,s,t] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +
		cq2u2G$7[p,r,s,t] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cq2u2G$8[p,r,s,t] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +

		cq2u2W$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2u2W$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2u2W$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2u2W$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cq2u2B$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2u2B$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[u[a2,s]]**\[Gamma][\[Mu]2]**u[a2,t]) dFS[B,\[Mu]1,\[Mu]2] +
		cq2u2B$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2u2B$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[u[a3,s]]**\[Gamma][\[Mu]2]**u[a4,t]) dFS[B,\[Mu]1,\[Mu]2] +

		cq2d2G$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) TSU3c[A1,a2,a3] (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2d2G$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) TSU3c[A1,a2,a3] (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a3,t]) dFS[G,\[Mu]1,\[Mu]2,A1] +
		cq2d2G$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2d2G$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2d2G$5[p,r,s,t] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cq2d2G$6[p,r,s,t] fSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +
		cq2d2G$7[p,r,s,t] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[G,\[Mu]1,\[Mu]2,A3] +
		cq2d2G$8[p,r,s,t] dSU3c[A1,A2,A3] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A2,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[G,\[Mu]1,\[Mu]2,A3] +

		cq2d2W$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2d2W$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +
		cq2d2W$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2d2W$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) \[Tau]SU2L[I1,i1,i2] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[W,\[Mu]1,\[Mu]2,I1] +

		cq2d2B$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2d2B$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[d[a2,s]]**\[Gamma][\[Mu]2]**d[a2,t]) dFS[B,\[Mu]1,\[Mu]2] +
		cq2d2B$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2d2B$4[p,r,s,t] (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[d[a3,s]]**\[Gamma][\[Mu]2]**d[a4,t]) dFS[B,\[Mu]1,\[Mu]2] +

	PlusHc[(* (LR)(RL)X *)
		cledqG$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) TSU3c[A1,a1,a2] (Bar[d[a1,s]]**q[a2,i1,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cledqG$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) TSU3c[A1,a1,a2] (Bar[d[a1,s]]**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i1,t]) FS[G,\[Mu]1,\[Mu]2,A1] +

		cledqW$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**q[a1,i2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cledqW$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] (Bar[d[a1,s]]**\[Sigma][\[Mu]1,\[Mu]2]**q[a1,i2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +

		cledqB$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) (Bar[d[a1,s]]**q[a1,i1,t]) FS[B,\[Mu]1,\[Mu]2] +
		cledqB$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) (Bar[d[a1,s]]**\[Sigma][\[Mu]1,\[Mu]2]**q[a1,i1,t]) FS[B,\[Mu]1,\[Mu]2]

	] +

	PlusHc[ (* (LR)(LR)X *)
		cq2udG$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a3,i2,s]]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2udG$2[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) TSU3c[A1,a2,a3] \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2udG$3[p,r,s,t] (Bar[q[a1,i1,p]]**u[a2,r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a3,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2udG$4[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) TSU3c[A1,a2,a3] \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2udG$5[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**u[a2,r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a3,i2,s]]**\[Sigma][\[Mu]3,\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		cq2udG$6[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**u[a1,r]) TSU3c[A1,a2,a3] \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**\[Sigma][\[Mu]3,\[Mu]2]**d[a3,t]) FS[G,\[Mu]1,\[Mu]2,A1] +

		cq2udW$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**d[a2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2udW$2[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		cq2udW$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**u[a1,r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a2,i3,s]]**\[Sigma][\[Mu]3,\[Mu]2]**d[a2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +

		cq2udB$1[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2udB$2[p,r,s,t] (Bar[q[a1,i1,p]]**u[a1,r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +
		cq2udB$3[p,r,s,t] (Bar[q[a1,i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**u[a1,r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a2,i2,s]]**\[Sigma][\[Mu]3,\[Mu]2]**d[a2,t]) FS[B,\[Mu]1,\[Mu]2] +

		clequG$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**u[a2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		clequG$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
		clequG$3[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**e[r]) TSU3c[A1,a1,a2] \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**\[Sigma][\[Mu]3,\[Mu]2]**u[a2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +

		clequW$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a1,i3,s]]**u[a1,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		clequW$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a1,i3,s]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
		clequW$3[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**e[r]) \[Tau]SU2L[I1,i1,i2] \[Epsilon]SU2L[i2,i3] (Bar[q[a1,i3,s]]**\[Sigma][\[Mu]3,\[Mu]2]**u[a1,t]) FS[W,\[Mu]1,\[Mu]2,I1] +

		clequB$1[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]2]**e[r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**u[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		clequB$2[p,r,s,t] (Bar[l[i1,p]]**e[r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**\[Sigma][\[Mu]1,\[Mu]2]**u[a1,t]) FS[B,\[Mu]1,\[Mu]2] +
		clequB$3[p,r,s,t] (Bar[l[i1,p]]**\[Sigma][\[Mu]1,\[Mu]3]**e[r]) \[Epsilon]SU2L[i1,i2] (Bar[q[a1,i2,s]]**\[Sigma][\[Mu]3,\[Mu]2]**u[a1,t]) FS[B,\[Mu]1,\[Mu]2]
	] +

	If[TrueQ@ BLViolation,
		PlusHc[ (* B-violating *)
			clqudG$1[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			clqudG$2[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			clqudG$3[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4]+TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			clqudG$4[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4]+TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
	
			clqudW$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
			clqudW$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i3,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
	
			clqudB$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i2,t]) FS[B,\[Mu]1,\[Mu]2] +
			clqudB$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**u[a2,r]) (Transp[q[a3,i1,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[B,\[Mu]1,\[Mu]2] +
	
			ceq2uG$1[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i2,r]) (Transp[u[a3,s]]**CC**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			ceq2uG$2[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4] + TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**q[a2,i2,r]) (Transp[u[a3,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
	
			ceq2uW$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i3,r]) (Transp[u[a3,s]]**CC**e[t]) FS[W,\[Mu]1,\[Mu]2,I1] +
	
			ceq2uB$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**q[a2,i2,r]) (Transp[u[a3,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**e[t]) FS[B,\[Mu]1,\[Mu]2] +
	
			clq3G$1[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			clq3G$2[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4] + TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
	
			clq3W$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] Bar@\[Epsilon]SU2L[i4,i5] (Transp[q[a1,i1,p]]**CC**q[a2,i4,r]) (Transp[q[a3,i5,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i3,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
			clq3W$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] Bar@\[Epsilon]SU2L[i4,i5] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**l[i5,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
	
			clq3B$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[B,\[Mu]1,\[Mu]2] +
	
			ceu2dG$1[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[u[a3,s]]**CC**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			ceu2dG$2[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] (Transp[u[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[d[a3,s]]**CC**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			ceu2dG$3[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4] + TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) (Transp[u[a1,p]]**CC**u[a2,r]) (Transp[d[a3,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**e[t]) FS[G,\[Mu]1,\[Mu]2,A1] +
	
			ceu2dB$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] (Transp[d[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[u[a3,s]]**CC**e[t]) FS[B,\[Mu]1,\[Mu]2] +
			ceu2dB$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] (Transp[u[a1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**u[a2,r]) (Transp[d[a3,s]]**CC**e[t]) FS[B,\[Mu]1,\[Mu]2] +
	
			ceq2uW$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] \[Tau]SU2L[I1,i2,i3] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[u[a3,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**e[t]) FS[W,\[Mu]1,\[Mu]2,I1] +
	
			ceq2uB$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i2,r]) (Transp[u[a3,s]]**CC**e[t]) FS[B,\[Mu]1,\[Mu]2] +
	
			clq3G$3[p,r,s,t] (TSU3c[A1,a4,a1] Bar@\[Epsilon]SU3c[a2,a3,a4] + TSU3c[A1,a4,a2] Bar@\[Epsilon]SU3c[a1,a3,a4]) Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
			clq3G$4[p,r,s,t] TSU3c[A1,a4,a3] Bar@\[Epsilon]SU3c[a4,a1,a2] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[G,\[Mu]1,\[Mu]2,A1] +
	
			clq3W$3[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] \[Tau]SU2L[I1,i4,i5] (Transp[q[a1,i1,p]]**CC**q[a2,i3,r]) (Transp[q[a3,i5,s]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**l[i2,t]) FS[W,\[Mu]1,\[Mu]2,I1] +
	
			clq3B$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a2,i3,r]) (Transp[q[a3,i4,s]]**CC**l[i2,t]) FS[B,\[Mu]1,\[Mu]2]
		]
		,
		0
	];

	L20 = PlusHc[
		cl3eHD$1[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**e[t]) CD[\[Mu]1,H[i2]] +
		cl3eHD$2[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[l[i3,s]]**e[t]) CD[\[Mu]1,H[i4]] +
		cl3eHD$3[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[CD[\[Mu]1,l[i2,s]]]**e[t]) H[i2] +

		cle3HD$1[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[l[i1,s]]**CD[\[Mu]1,e[t]]) H[i1] +

		cleq2HD$1[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[l[i2,s]]**e[t]) CD[\[Mu]1,H[i2]] +
		cleq2HD$2[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a1,i2,s]]**e[t]) CD[\[Mu]1,H[i2]] +
		cleq2HD$3[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[l[i3,s]]**e[t]) CD[\[Mu]1,H[i4]] +
		cleq2HD$4[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**e[t]) CD[\[Mu]1,H[i4]] +
		cleq2HD$5[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[l[i2,s]]**CD[\[Mu]1,e[t]]) H[i2] +
		cleq2HD$6[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4]  (Bar[l[i3,s]]**CD[\[Mu]1,e[t]]) H[i4] +

		cleu2HD$1[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[l[i1,s]]**e[t]) CD[\[Mu]1,H[i1]] +
		cleu2HD$2[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[l[i1,s]]**u[a1,t]) CD[\[Mu]1,H[i1]] +
		cleu2HD$3[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[CD[\[Mu]1,l[i1,s]]]**e[t]) H[i1] +

		cled2HD$1[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[l[i1,s]]**e[t]) CD[\[Mu]1,H[i1]] +
		cled2HD$2[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[l[i1,s]]**d[a1,t]) CD[\[Mu]1,H[i1]] +
		cled2HD$3[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[CD[\[Mu]1,l[i1,s]]]**e[t]) H[i1] +

		cleudHD$1[p,r,s,t] I Bar@\[Epsilon]SU2L[i1,i2] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[e[s]]**l[i1,t]) CD[\[Mu]1,H[i2]] +
		cleudHD$2[p,r,s,t] I Bar@\[Epsilon]SU2L[i1,i2] (Bar[e[p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[u[a1,s]]**l[i1,t]) CD[\[Mu]1,H[i2]] +
		cleudHD$3[p,r,s,t] I Bar@\[Epsilon]SU2L[i1,i2] (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[e[s]]**CD[\[Mu]1,l[i1,t]]) H[i2] +

		cle3HD$2[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[l[i1,s]]**e[t]) CD[\[Mu]1,H[i1]] +

		cl2quHD$1[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**u[a1,t]) CD[\[Mu]1,dH[i2]] +
		cl2quHD$2[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**u[a1,t]) CD[\[Mu]1,dH[i2]] +
		cl2quHD$3[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**u[a1,t]) CD[\[Mu]1,dH[i4]] +
		cl2quHD$4[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[l[i3,s]]**u[a1,t]) CD[\[Mu]1,dH[i4]] +
		cl2quHD$5[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**CD[\[Mu]1,u[a1,t]]) dH[i2] +
		cl2quHD$6[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**CD[\[Mu]1,u[a1,t]]) dH[i4] +

		ce2quHD$1[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[q[a1,i1,s]]**u[a1,t]) CD[\[Mu]1,dH[i1]] +
		ce2quHD$2[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[q[a1,i1,s]]**e[t]) CD[\[Mu]1,dH[i1]] +
		ce2quHD$3[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[CD[\[Mu]1,q[a1,i1,s]]]**u[a1,t]) dH[i1] +

		cq3uHD$1[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**u[a2,t]) CD[\[Mu]1,dH[i2]] +
		cq3uHD$2[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**u[a2,t]) CD[\[Mu]1,dH[i4]] +
		cq3uHD$3[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**u[a4,t]) CD[\[Mu]1,dH[i2]] +
		cq3uHD$4[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a3,i3,s]]**u[a4,t]) CD[\[Mu]1,dH[i4]] +
		cq3uHD$5[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[CD[\[Mu]1,q[a2,i2,s]]]**u[a2,t]) dH[i2] +
		cq3uHD$6[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[CD[\[Mu]1,q[a2,i3,s]]]**u[a2,t]) dH[i4] +

		cqu3HD$1[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[q[a2,i1,s]]**u[a2,t]) CD[\[Mu]1,dH[i1]] +
		cqu3HD$2[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**u[a4,t]) CD[\[Mu]1,dH[i1]] +
		cqu3HD$3[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[q[a2,i1,s]]**CD[\[Mu]1,u[a2,t]]) dH[i1] +

		cqud2HD$1[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[q[a2,i1,s]]**u[a2,t]) CD[\[Mu]1,dH[i1]] +
		cqud2HD$2[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[q[a2,i1,s]]**d[a2,t]) CD[\[Mu]1,dH[i1]] +
		cqud2HD$3[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**u[a4,t]) CD[\[Mu]1,dH[i1]] +
		cqud2HD$4[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**d[a4,t]) CD[\[Mu]1,dH[i1]] +
		cqud2HD$5[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[CD[\[Mu]1,q[a2,i1,s]]]**u[a2,t]) dH[i1] +
		cqud2HD$6[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[CD[\[Mu]1,q[a3,i1,s]]]**u[a4,t]) dH[i1] +

		cl2qdHD$1[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**d[a1,t]) CD[\[Mu]1,H[i2]] +
		cl2qdHD$2[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[l[i2,s]]**d[a1,t]) CD[\[Mu]1,H[i2]] +
		cl2qdHD$3[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**d[a1,t]) CD[\[Mu]1,H[i4]] +
		cl2qdHD$4[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[l[i3,s]]**d[a1,t]) CD[\[Mu]1,H[i4]] +
		cl2qdHD$5[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r]) (Bar[q[a1,i2,s]]**CD[\[Mu]1,d[a1,t]]) H[i2] +
		cl2qdHD$6[p,r,s,t] I (Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a1,i3,s]]**CD[\[Mu]1,d[a1,t]]) H[i4] +

		ce2qdHD$1[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[q[a1,i1,s]]**d[a1,t]) CD[\[Mu]1,H[i1]] +
		ce2qdHD$2[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[q[a1,i1,s]]**e[t]) CD[\[Mu]1,H[i1]] +
		ce2qdHD$3[p,r,s,t] I (Bar[e[p]]**\[Gamma][\[Mu]1]**e[r]) (Bar[CD[\[Mu]1,q[a1,i1,s]]]**d[a1,t]) H[i1] +

		cq3dHD$1[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[q[a2,i2,s]]**d[a2,t]) CD[\[Mu]1,H[i2]] +
		cq3dHD$2[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a2,i3,s]]**d[a2,t]) CD[\[Mu]1,H[i4]] +
		cq3dHD$3[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i2,s]]**d[a4,t]) CD[\[Mu]1,H[i2]] +
		cq3dHD$4[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[q[a3,i3,s]]**d[a4,t]) CD[\[Mu]1,H[i4]] +
		cq3dHD$5[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r]) (Bar[CD[\[Mu]1,q[a2,i2,s]]]**d[a2,t]) H[i2] +
		cq3dHD$6[p,r,s,t] I (Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r]) \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] (Bar[CD[\[Mu]1,q[a2,i3,s]]]**d[a2,t]) H[i4] +

		cqu2dHD$1[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[q[a2,i1,s]]**d[a2,t]) CD[\[Mu]1,H[i1]] +
		cqu2dHD$2[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[q[a2,i1,s]]**u[a2,t]) CD[\[Mu]1,H[i1]] +
		cqu2dHD$3[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**d[a4,t]) CD[\[Mu]1,H[i1]] +
		cqu2dHD$4[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**d[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**u[a4,t]) CD[\[Mu]1,H[i1]] +
		cqu2dHD$5[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r]) (Bar[CD[\[Mu]1,q[a2,i1,s]]]**d[a2,t]) H[i1] +
		cqu2dHD$6[p,r,s,t] I (Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[CD[\[Mu]1,q[a3,i1,s]]]**d[a4,t]) H[i1] +

		cqd3HD$1[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[q[a2,i1,s]]**d[a2,t]) CD[\[Mu]1,H[i1]] +
		cqd3HD$2[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a2,r]) TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] (Bar[q[a3,i1,s]]**d[a4,t]) CD[\[Mu]1,H[i1]] +
		cqd3HD$3[p,r,s,t] I (Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r]) (Bar[q[a2,i1,s]]**CD[\[Mu]1,d[a2,t]]) H[i1] +

		(* B-violating *)
		If[TrueQ@ BLViolation,
			clu2dHD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[CD[\[Mu]1,H[i1]]] (Transp[u[a1,p]]**CC**\[Gamma][\[Mu]1]**l[i1,r]) (Transp[u[a2,s]]**CC**d[a3,t]) +
			clu2dHD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[H[i1]] (Transp[u[a1,p]]**CC**\[Gamma][\[Mu]1]**l[i1,r]) (Transp[CD[\[Mu]1,u[a2,s]]]**CC**d[a3,t]) +
	
			clud2HD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**\[Gamma][\[Mu]1]**l[i1,r]) (Transp[d[a2,s]]**CC**u[a3,t]) CD[\[Mu]1,H[i2]] +
			clud2HD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[d[a1,p]]**CC**\[Gamma][\[Mu]1]**l[i1,r]) (Transp[CD[\[Mu]1,d[a2,s]]]**CC**u[a3,t]) H[i2] +
	
			(*
			clq2uHD$1[p,r,s,t] \[ImaginaryI] Bar@\[Epsilon]SU3c[a1,a2,a3] \[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Bar[CD[\[Mu]1,H[i2]]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) +
			clq2uHD$2[p,r,s,t] \[ImaginaryI] Bar@\[Epsilon]SU3c[a1,a2,a3] \[Epsilon]SU2L[i3,i2] Bar@\[Epsilon]SU2L[i1,i4] Bar[CD[\[Mu]1,H[i2]]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) +
			clq2uHD$3[p,r,s,t] \[ImaginaryI] Bar@\[Epsilon]SU3c[a1,a2,a3] \[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] Bar[H[i2]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[CD[\[Mu]1,q[a3,i1,s]]]**CC**l[i3,t]) +
			*) (* The terms above from the Murphy paper are not gauge invariant due to the incorrectly placed epsilon tensor *)
			clq2uHD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i3,i4] Bar[CD[\[Mu]1,H[i1]]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) +
			clq2uHD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i4] Bar[CD[\[Mu]1,H[i3]]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) +
			clq2uHD$3[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i3,i4] Bar[H[i1]] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[CD[\[Mu]1,q[a3,i1,s]]]**CC**l[i3,t]) +
	
			clq2dHD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**d[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) CD[\[Mu]1,H[i2]] +
			clq2dHD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i3,i2] Bar@\[Epsilon]SU2L[i1,i4] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**d[a2,r]) (Transp[q[a3,i1,s]]**CC**l[i3,t]) CD[\[Mu]1,H[i2]] +
			clq2dHD$3[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i4,p]]**CC**\[Gamma][\[Mu]1]**d[a2,r]) (Transp[CD[\[Mu]1,q[a3,i1,s]]]**CC**l[i3,t]) H[i2] +
	
			ceq3HD[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i4,i2] Bar@\[Epsilon]SU2L[i1,i3] (Transp[q[a1,i1,p]]**CC**\[Gamma][\[Mu]1]**e[r]) (Transp[CD[\[Mu]1,q[a2,i3,s]]]**CC**q[a3,i4,t]) H[i2] +
	
			cequ2HD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[CD[\[Mu]1,H[i1]]] (Transp[u[a1,p]]**CC**\[Gamma][\[Mu]1]**q[a2,i1,r]) (Transp[u[a3,s]]**CC**e[t]) +
			cequ2HD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar[H[i1]] (Transp[u[a1,p]]**CC**\[Gamma][\[Mu]1]**q[a2,i1,r]) (Transp[CD[\[Mu]1,u[a3,s]]]**CC**e[t]) +
	
			cequdHD$1[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[d[a3,s]]**CC**e[t]) CD[\[Mu]1,H[i2]] +
			cequdHD$2[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**\[Gamma][\[Mu]1]**d[a2,r]) (Transp[u[a3,s]]**CC**e[t]) CD[\[Mu]1,H[i2]] +
			cequdHD$3[p,r,s,t] I Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**\[Gamma][\[Mu]1]**u[a2,r]) (Transp[d[a3,s]]**CC**CD[\[Mu]1,e[t]]) H[i2]
		,
			0
		]
	];

	L21 = PlusHc[
		(* (LR)(RL)DD *)
		cleqdD2$1[p,r,s,t] CD[\[Mu]1,(Bar[l[i1,p]]**e[r])] CD[\[Mu]1,(Bar[d[a1,s]]**q[a1,i1,t])] +
		cleqdD2$2[p,r,s,t] HermitianCD[\[Mu]1,Bar[l[i1,p]],e[r]] HermitianCD[\[Mu]1,Bar[d[a1,s]],q[a1,i1,t]]
	] +

	If[TrueQ@ BLViolation,
		PlusHc[
			(* B-violating *)
			clqudD2$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Transp[d[a1,p]]**CC**u[a2,r])] CD[\[Mu]1,(Transp[q[a3,i1,s]]**CC**l[i2,t])] +
			(*clqudD2$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Transp[d[a1,p]]**CC**q[a2,i1,r])] CD[\[Mu]1,(Transp[u[a3,s]]**CC**l[i2,t])] +*)
			(* The operator QlqudD2$2 from the operator tables in Murphy's paper vanishes because spinors with different chiralities are directly contracted.
				This operator is replaced with the operator below which is taken from eq. 4.77 in https://arxiv.org/pdf/2005.00008 *)
			clqudD2$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Transp[l[i1,p]]**CC**\[Sigma][\[Mu]1,\[Mu]2]**q[a1,i2,r] Transp[CD[\[Mu]1,d[a2,s]]]**CC**CD[\[Mu]2,u[a3,t]] +
	
			ceq2uD2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] (Transp[q[a1,i1,p]]**CC**CD[\[Mu]1,q[a2,i2,r]]) CD[\[Mu]1,(Transp[u[a3,s]]**CC**e[t])] +
	
			clq3D2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] Bar@\[Epsilon]SU2L[i1,i2] Bar@\[Epsilon]SU2L[i3,i4] (Transp[q[a1,i3,p]]**CC**CD[\[Mu]1,q[a2,i1,r]]) CD[\[Mu]1,(Transp[q[a3,i2,s]]**CC**l[i4,t])] +
	
			ceu2dD2$1[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] (Transp[u[a1,p]]**CC**CD[\[Mu]1,u[a2,r]]) CD[\[Mu]1,(Transp[d[a3,s]]**CC**e[t])] +
			ceu2dD2$2[p,r,s,t] Bar@\[Epsilon]SU3c[a1,a2,a3] (Transp[u[a1,p]]**CC**u[a2,r]) (Transp[CD[\[Mu]1,d[a3,s]]]**CC**CD[\[Mu]1,e[t]])
		]
	,
		0
	] +

	(* (LL)(LL)DD *)
	cl4D2$1[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r])] CD[\[Mu]2,(Bar[l[i2,s]]**\[Gamma][\[Mu]1]**l[i2,t])] +
	cl4D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i1,r]] HermitianCD[\[Mu]2,Bar[l[i2,s]]**\[Gamma][\[Mu]1],l[i2,t]] +

	cq4D2$1[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r])] CD[\[Mu]2,(Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1]**q[a2,i2,t])] +
	cq4D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a1,i1,r]] HermitianCD[\[Mu]2,Bar[q[a2,i2,s]]**\[Gamma][\[Mu]1],q[a2,i2,t]] +
	cq4D2$3[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i2,r])] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] CD[\[Mu]2,(Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1]**q[a2,i4,t])] +
	cq4D2$4[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a1,i2,r]] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] HermitianCD[\[Mu]2,Bar[q[a2,i3,s]]**\[Gamma][\[Mu]1],q[a2,i4,t]] +

	cl2q2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r])] CD[\[Mu]2,(Bar[q[a1,i2,s]]**\[Gamma][\[Mu]1]**q[a1,i2,t])] +
	cl2q2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i1,r]] HermitianCD[\[Mu]2,Bar[q[a1,i2,s]]**\[Gamma][\[Mu]1],q[a1,i2,t]] +
	cl2q2D2$3[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i2,r])] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] CD[\[Mu]2,(Bar[q[a1,i3,s]]**\[Gamma][\[Mu]1]**q[a1,i4,t])] +
	cl2q2D2$4[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i2,r]] \[Tau]SU2L[I1,i1,i2] \[Tau]SU2L[I1,i3,i4] HermitianCD[\[Mu]2,Bar[q[a1,i3,s]]**\[Gamma][\[Mu]1],q[a1,i4,t]] +

	(* (RR)(RR)DD *)
	ce4D2[p,r,s,t] CD[\[Mu]2,(Bar[e[p]]**\[Gamma][\[Mu]1]**e[r])] CD[\[Mu]2,(Bar[e[s]]**\[Gamma][\[Mu]1]**e[t])] +

	cu4D2$1[p,r,s,t] CD[\[Mu]2,(Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r])] CD[\[Mu]2,(Bar[u[a2,s]]**\[Gamma][\[Mu]1]**u[a2,t])] +
	cu4D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[u[a1,p]]**\[Gamma][\[Mu]1],u[a1,r]] HermitianCD[\[Mu]2,Bar[u[a2,s]]**\[Gamma][\[Mu]1],u[a2,t]] +

	cd4D2$1[p,r,s,t] CD[\[Mu]2,(Bar[d[a1,p]]**\[Gamma][\[Mu]1]**d[a1,r])] CD[\[Mu]2,(Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t])] +
	cd4D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[d[a1,p]]**\[Gamma][\[Mu]1],d[a1,r]] HermitianCD[\[Mu]2,Bar[d[a2,s]]**\[Gamma][\[Mu]1],d[a2,t]] +

	ce2u2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[e[p]]**\[Gamma][\[Mu]1]**e[r])] CD[\[Mu]2,(Bar[u[a1,s]]**\[Gamma][\[Mu]1]**u[a1,t])] +
	ce2u2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[e[p]]**\[Gamma][\[Mu]1],e[r]] HermitianCD[\[Mu]2,Bar[u[a1,s]]**\[Gamma][\[Mu]1],u[a1,t]] +

	ce2d2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[e[p]]**\[Gamma][\[Mu]1]**e[r])] CD[\[Mu]2,(Bar[d[a1,s]]**\[Gamma][\[Mu]1]**d[a1,t])] +
	ce2d2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[e[p]]**\[Gamma][\[Mu]1],e[r]] HermitianCD[\[Mu]2,Bar[d[a1,s]]**\[Gamma][\[Mu]1],d[a1,t]] +

	cu2d2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a1,r])] CD[\[Mu]2,(Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t])] +
	cu2d2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[u[a1,p]]**\[Gamma][\[Mu]1],u[a1,r]] HermitianCD[\[Mu]2,Bar[d[a2,s]]**\[Gamma][\[Mu]1],d[a2,t]] +
	cu2d2D2$3[p,r,s,t] CD[\[Mu]2,(Bar[u[a1,p]]**\[Gamma][\[Mu]1]**u[a2,r])] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] CD[\[Mu]2,(Bar[d[a3,s]]**\[Gamma][\[Mu]1]**d[a4,t])] +
	cu2d2D2$4[p,r,s,t] HermitianCD[\[Mu]2,Bar[u[a1,p]]**\[Gamma][\[Mu]1],u[a2,r]] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] HermitianCD[\[Mu]2,Bar[d[a3,s]]**\[Gamma][\[Mu]1],d[a4,t]] +

	(* (LL)(RR)DD *)
	cl2e2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r])] CD[\[Mu]2,(Bar[e[s]]**\[Gamma][\[Mu]1]**e[t])] +
	cl2e2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i1,r]] HermitianCD[\[Mu]2,Bar[e[s]]**\[Gamma][\[Mu]1],e[t]] +

	cl2u2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r])] CD[\[Mu]2,(Bar[u[a1,s]]**\[Gamma][\[Mu]1]**u[a1,t])] +
	cl2u2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i1,r]] HermitianCD[\[Mu]2,Bar[u[a1,s]]**\[Gamma][\[Mu]1],u[a1,t]] +

	cl2d2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[l[i1,p]]**\[Gamma][\[Mu]1]**l[i1,r])] CD[\[Mu]2,(Bar[d[a1,s]]**\[Gamma][\[Mu]1]**d[a1,t])] +
	cl2d2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[l[i1,p]]**\[Gamma][\[Mu]1],l[i1,r]] HermitianCD[\[Mu]2,Bar[d[a1,s]]**\[Gamma][\[Mu]1],d[a1,t]] +

	cq2e2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r])] CD[\[Mu]2,(Bar[e[s]]**\[Gamma][\[Mu]1]**e[t])] +
	cq2e2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a1,i1,r]] HermitianCD[\[Mu]2,Bar[e[s]]**\[Gamma][\[Mu]1],e[t]] +

	cq2u2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r])] CD[\[Mu]2,(Bar[u[a2,s]]**\[Gamma][\[Mu]1]**u[a2,t])] +
	cq2u2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a1,i1,r]] HermitianCD[\[Mu]2,Bar[u[a2,s]]**\[Gamma][\[Mu]1],u[a2,t]] +
	cq2u2D2$3[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r])] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] CD[\[Mu]2,(Bar[u[a3,s]]**\[Gamma][\[Mu]1]**u[a4,t])] +
	cq2u2D2$4[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a2,i1,r]] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] HermitianCD[\[Mu]2,Bar[u[a3,s]]**\[Gamma][\[Mu]1],u[a4,t]] +

	cq2d2D2$1[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a1,i1,r])] CD[\[Mu]2,(Bar[d[a2,s]]**\[Gamma][\[Mu]1]**d[a2,t])] +
	cq2d2D2$2[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a1,i1,r]] HermitianCD[\[Mu]2,Bar[d[a2,s]]**\[Gamma][\[Mu]1],d[a2,t]] +
	cq2d2D2$3[p,r,s,t] CD[\[Mu]2,(Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1]**q[a2,i1,r])] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] CD[\[Mu]2,(Bar[d[a3,s]]**\[Gamma][\[Mu]1]**d[a4,t])] +
	cq2d2D2$4[p,r,s,t] HermitianCD[\[Mu]2,Bar[q[a1,i1,p]]**\[Gamma][\[Mu]1],q[a2,i1,r]] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] HermitianCD[\[Mu]2,Bar[d[a3,s]]**\[Gamma][\[Mu]1],d[a4,t]] +

	PlusHc[
		(* (LR)(LR)DD *)
		cq2udD2$1[p,r,s,t] CD[\[Mu]1,(Bar[q[a1,i1,p]]**u[a1,r])] \[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Bar[q[a2,i2,s]]**d[a2,t])] +
		cq2udD2$2[p,r,s,t] CD[\[Mu]1,(Bar[q[a1,i1,p]]**u[a2,r])] TSU3c[A1,a1,a2] TSU3c[A1,a3,a4] \[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Bar[q[a3,i2,s]]**d[a4,t])] +
		cq2udD2$3[p,r,s,t] HermitianCD[\[Mu]1,Bar[q[a1,i1,p]],u[a1,r]] \[Epsilon]SU2L[i1,i2] HermitianCD[\[Mu]1,Bar[q[a2,i2,s]],d[a2,t]] +

		clequD2$1[p,r,s,t] CD[\[Mu]1,(Bar[l[i1,p]]**e[r])] \[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Bar[q[a1,i2,s]]**u[a1,t])] +
		clequD2$2[p,r,s,t] CD[\[Mu]1,(Bar[l[i1,p]]**u[a1,r])] \[Epsilon]SU2L[i1,i2] CD[\[Mu]1,(Bar[q[a1,i2,s]]**e[t])] +
		clequD2$3[p,r,s,t] HermitianCD[\[Mu]1,Bar[l[i1,p]],e[r]] \[Epsilon]SU2L[i1,i2] HermitianCD[\[Mu]1,Bar[q[a1,i2,s]],u[a1,t]]
	];

	Dim7Lag = L7$1 + L7$2 + L7$3 + L7$4 + L7$5 + L7$6;
	BosonicLag     = L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8;
	TwoFermionLag  = L9 + L10 + L11 + L12 + L13 + L14 + L15 + L16 + L17;
	FourFermionLag = L18 + L19 + L20 + L21;

	RelabelIndices[
		(Dim7Lag + BosonicLag + TwoFermionLag + FourFermionLag)/.{
			FSG:FieldStrength[G,___]:>FSG/-gs[],
			FSW:FieldStrength[W,___]:>FSW/-gL[],
			FSB:FieldStrength[B,___]:>FSB/-gY[]
		} (* map from Matchete gauge coupling conventions to conventions of [arXiv:2005.00059] *)
	]
]
