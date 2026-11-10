(* ::Package:: *)

(* ::Title:: *)
(*EW sector of the SMEFT model file*)


(* ::Text:: *)
(*Based  on  the  results  of [2503.20016] .*)


(* ::Text:: *)
(*Gauge  groups*)


DefineGaugeGroup[SU2L, SU@2, gL, W, 
	FundAlphabet -> CharacterRange["\[Alpha]", "\[Lambda]"], 
	AdjAlphabet -> CharacterRange["I", "N"], 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(L\)]\)", Default}];
DefineGaugeGroup[U1Y, U1, gY, B, 
	NiceForm->{"\!\(\*SubscriptBox[\(g\), \(Y\)]\)", Default}];


(* ::Text:: *)
(*Flavor  index*)


ParameterDefault[Nf-> 3]
DefineFlavorIndex[Flavor, Nf,
	IndexAlphabet-> {"p", "r", "s", "t", "u", "v"}];


(* ::Text:: *)
(*Matter  fields*)


DefineField[q, Fermion, Indices-> {SU2L[fund], Flavor}, 
	Charges->{U1Y[1/6]}, Chiral-> LeftHanded, Mass-> 0];
DefineField[u,Fermion,Indices-> {Flavor},
	Charges->{U1Y[2/3]}, Chiral-> RightHanded, Mass-> 0];
DefineField[d,Fermion,Indices-> {Flavor},
	Charges->{U1Y[-1/3]}, Chiral-> RightHanded, Mass-> 0];
DefineField[l,Fermion,Indices-> {SU2L[fund], Flavor},
	Charges->{U1Y[-1/2]}, Chiral-> LeftHanded, Mass-> 0];
DefineField[e,Fermion,Indices-> {Flavor},
	Charges->{U1Y[-1]}, Chiral-> RightHanded, Mass-> 0];
DefineField[\[Phi], Scalar, Indices-> {SU2L[fund]}, 
	Charges-> {U1Y[1/2]}, Mass-> 0];


(* ::Text:: *)
(*Couplings*)


DefineCoupling[Yu, Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubscriptBox[\(Y\), \(u\)]\)"];
DefineCoupling[Yd, Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubscriptBox[\(Y\), \(d\)]\)"];
DefineCoupling[Ye, Indices-> {Flavor, Flavor},
	NiceForm-> "\!\(\*SubscriptBox[\(Y\), \(e\)]\)"];
(* Mass added as an interaction because it is tachyonic *)			
DefineCoupling[\[Mu]2, SelfConjugate-> True, EFTOrder-> 2,
	NiceForm-> "\[Mu]2"];
DefineCoupling[\[Lambda], SelfConjugate->True, EFTOrder-> 0];


(* ::Text:: *)
(*Effective Couplings*)


DefineCoupling[c\[Phi], SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]\)]\)"];
DefineCoupling[c\[Phi]Box, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi] \[Square]\)]\)"];
DefineCoupling[c\[Phi]D, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]D\)]\)"];


(*DefineCoupling[c\[Phi]B, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]B\)]\)"];
DefineCoupling[c\[Phi]Btilde, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi] Overscript[B, ~]\)]\)"];
DefineCoupling[c\[Phi]W, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]W\)]\)"];
DefineCoupling[c\[Phi]Wtilde, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi] Overscript[W, ~]\)]\)"];
DefineCoupling[c\[Phi]WB, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]WB\)]\)"];
DefineCoupling[c\[Phi]WtildeB, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi] Overscript[W, ~] B\)]\)"];*)


(*DefineCoupling[c3W, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(3W\)]\)"];
DefineCoupling[c3Wtilde, SelfConjugate->True, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(3 Overscript[W, ~]\)]\)"];*)


DefineCoupling[c\[Phi]q1, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(c\), \(\[Phi]q\), \((1)\)]\)"];
DefineCoupling[c\[Phi]q3, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(c\), \(\[Phi]q\), \((3)\)]\)"];
DefineCoupling[c\[Phi]l1, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(c\), \(\[Phi]l\), \((1)\)]\)"];
DefineCoupling[c\[Phi]l3, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(c\), \(\[Phi]l\), \((3)\)]\)"];


DefineCoupling[c\[Phi]u, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]u\)]\)"];
DefineCoupling[c\[Phi]d, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]d\)]\)"];
DefineCoupling[c\[Phi]e, SelfConjugate->{2,1}, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]e\)]\)"];


DefineCoupling[c\[Phi]ud, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(\[Phi]ud\)]\)"];
DefineCoupling[ce\[Phi], Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(e\[Phi]\)]\)"];
DefineCoupling[cu\[Phi], Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(u\[Phi]\)]\)"];
DefineCoupling[cd\[Phi], Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(c\), \(d\[Phi]\)]\)"];


(* ::Text:: *)
(*Lagrangian*)


Module[{\[Alpha],\[Beta],\[Delta],\[Lambda]1,\[Mu],\[Nu],p,r,\[CapitalIota],
		\[Tau]SU2L,HermitianCD,HermitianCD\[Tau],
		YukawaL,HiggsPotential,
		Lagd6\[Phi],Lagd6\[Phi]B,Lagd6\[Phi]W,Lagd6\[Phi]WB,Lagd6W,
		Lagd6\[Phi]q,Lagd6\[Phi]l,Lagd6\[Phi]u,Lagd6\[Phi]d,Lagd6\[Phi]e,
		Lagd6\[Phi]ud,Lagd6e\[Phi],Lagd6d\[Phi],Lagd6u\[Phi]},

	\[Tau]SU2L[Jadj_,ifund_,jfund_]:= 2 gen[SU2L[fund]][Jadj,ifund,jfund];
	HermitianCD[ind_,term1_,term2_]:= I term1 CD[ind,term2] - I CD[ind,term1] term2;
	
	YukawaL =  Yu[p,r] Bar@q[\[Alpha],p]**u[r] CG[eps[SU2L],{\[Alpha],\[Beta]}]Bar@\[Phi][\[Beta]] + Yd[p,r] Bar@q[\[Alpha],p]**d[r]\[Phi][\[Alpha]] + Ye[p,r] Bar@l[\[Alpha],p]**e[r]\[Phi][\[Alpha]];
	HiggsPotential = -\[Mu]2[] Bar@\[Phi][\[Alpha]]\[Phi][\[Alpha]] + \[Lambda][] Bar@\[Phi][\[Alpha]]\[Phi][\[Alpha]]Bar@\[Phi][\[Beta]]\[Phi][\[Beta]];
	
	Lagd6\[Phi]= c\[Phi][](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(Bar@\[Phi][\[Beta]]*\[Phi][\[Beta]])(Bar@\[Phi][\[Delta]]*\[Phi][\[Delta]])+c\[Phi]Box[](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])*CD[{\[Mu],\[Mu]},(Bar@\[Phi][\[Beta]]*\[Phi][\[Beta]])]+c\[Phi]D[]*(Bar@\[Phi][\[Alpha]]*CD[{\[Mu]},\[Phi][\[Alpha]]])*Bar[(Bar@\[Phi][\[Beta]]*CD[{\[Mu]},\[Phi][\[Beta]]])];
	
	(*Lagd6\[Phi]B= c\[Phi]B[] FS[B,\[Mu],\[Nu]]FS[B,\[Mu],\[Nu]](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(*+c\[Phi]Btilde[] FS[B,\[Mu],\[Nu]]FS[B,\[Mu],\[Nu]](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])*);
	Lagd6\[Phi]W= c\[Phi]W[] FS[W,\[Mu],\[Nu],\[CapitalIota]]FS[W,\[Mu],\[Nu],\[CapitalIota]](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(*+c\[Phi]Wtilde[] FS[W,\[Mu],\[Nu],\[CapitalIota]]FS[W,\[Mu],\[Nu],\[CapitalIota]](Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])*);
	Lagd6\[Phi]WB= c\[Phi]WB[] FS[W,\[Mu],\[Nu],\[CapitalIota]]FS[B,\[Mu],\[Nu]](Bar@\[Phi][\[Alpha]]*\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]]*\[Phi][\[Beta]])(*+c\[Phi]WtildeB[] FS[W,\[Mu],\[Nu],\[CapitalIota]]FS[B,\[Mu],\[Nu]](Bar@\[Phi][\[Alpha]]*\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]]*\[Phi][\[Beta]])*);
	
	Lagd6W=0;*)
	
	Lagd6\[Phi]q= c\[Phi]q1[p,r]*HermitianCD[\[Mu],Bar@\[Phi][\[Alpha]],\[Phi][\[Alpha]]]*(Bar@q[\[Beta],p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]q[\[Beta],r])+c\[Phi]q3[p,r](I Bar@\[Phi][\[Alpha]] CD[\[Mu],\[Phi][\[Beta]]]\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]] - I CD[\[Mu],Bar@\[Phi][\[Alpha]]]\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]] \[Phi][\[Beta]])(Bar@q[\[Delta],p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]q[\[Lambda]1,r])\[Tau]SU2L[\[CapitalIota],\[Delta],\[Lambda]1];
	Lagd6\[Phi]l= c\[Phi]l1[p,r]*HermitianCD[\[Mu],Bar@\[Phi][\[Alpha]],\[Phi][\[Alpha]]]*(Bar@l[\[Beta],p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]l[\[Beta],r])+c\[Phi]l3[p,r](I Bar@\[Phi][\[Alpha]] CD[\[Mu],\[Phi][\[Beta]]]\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]] - I CD[\[Mu],Bar@\[Phi][\[Alpha]]]\[Tau]SU2L[\[CapitalIota],\[Alpha],\[Beta]] \[Phi][\[Beta]])(Bar@l[\[Delta],p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]l[\[Lambda]1,r])\[Tau]SU2L[\[CapitalIota],\[Delta],\[Lambda]1];
	
	Lagd6\[Phi]u= c\[Phi]u[p,r]*HermitianCD[\[Mu],Bar@\[Phi][\[Alpha]],\[Phi][\[Alpha]]]*(Bar@u[p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]u[r]);
	Lagd6\[Phi]d= c\[Phi]d[p,r]*HermitianCD[\[Mu],Bar@\[Phi][\[Alpha]],\[Phi][\[Alpha]]]*(Bar@d[p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]d[r]);
	Lagd6\[Phi]e= c\[Phi]e[p,r]*HermitianCD[\[Mu],Bar@\[Phi][\[Alpha]],\[Phi][\[Alpha]]]*(Bar@e[p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]e[r]);
	
	Lagd6\[Phi]ud= c\[Phi]ud[p,r]*Bar@(CG[eps[SU2L],{\[Alpha],\[Beta]}]Bar@\[Phi][\[Beta]])*I*CD[\[Mu],\[Phi][\[Alpha]]](Bar@u[p]\[CenterDot]\[Gamma][\[Mu]]\[CenterDot]d[r]);
	Lagd6e\[Phi]= ce\[Phi][p,r]*(Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(Bar@l[\[Beta],p]\[CenterDot]e[r])*\[Phi][\[Beta]];
	Lagd6d\[Phi]= cd\[Phi][p,r]*(Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(Bar@q[\[Beta],p]\[CenterDot]d[r])*\[Phi][\[Beta]];
	Lagd6u\[Phi]= cu\[Phi][p,r]*(Bar@\[Phi][\[Alpha]]*\[Phi][\[Alpha]])(Bar@q[\[Beta],p]\[CenterDot]u[r])*CG[eps[SU2L],{\[Beta],\[Delta]}]Bar@\[Phi][\[Delta]];
	
	
	FreeLag[q, u, d, l, e, \[Phi], W, B] - PlusHc[YukawaL] - HiggsPotential + Lagd6\[Phi] (*+ Lagd6\[Phi]B + Lagd6\[Phi]W + Lagd6\[Phi]WB + Lagd6W*) + Lagd6\[Phi]q + Lagd6\[Phi]l + Lagd6\[Phi]u + Lagd6\[Phi]d + Lagd6\[Phi]e + PlusHc[Lagd6\[Phi]ud + Lagd6e\[Phi] +Lagd6d\[Phi] +Lagd6u\[Phi]]//RelabelIndices
]
