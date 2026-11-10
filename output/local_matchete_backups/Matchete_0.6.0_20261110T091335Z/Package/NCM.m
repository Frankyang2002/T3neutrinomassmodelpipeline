(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`NCM`*)


(* ::Subtitle:: *)
(*Contains the definitions of the NCM*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["Bar"]
PackageExport["NCM"]
PackageExport["GrassmannProduct"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["NonBarableQ"]
PackageScope["SetBarable"]
PackageScope["UnsetBarable"]


PackageScope["CommutativeQ"]
PackageScope["Commutative"]


PackageScope["FermionTrace"]
PackageScope["CanonizeSpinorLines"]


PackageScope["OpenSpinChainQ"]
PackageScope["VectorLikeSpinChainQ"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


NCM::usage = "NCM[a, b, ...] (alternatively a\[CenterDot] b\[CenterDot] ...) denotes the Dirac product of spinors, Dirac matrices, etc.";


GrassmannProduct::usage = "GrassmannProduct[a, b, ...] denotes the ordered product of Grassmann-valued fields (ghosts and fermions). It preserves ordering for correct Grassmann signs in functional derivatives.";


Bar::usage   = "Bar[field] returns the bar of a fermion field or the conjugate otherwise. Bar can also be applied to conjugate indices and charges.";


(* ::Subsubsection::Closed:: *)
(*Internal*)


NonBarableQ::usage  = "Specifies whether a term is trivial under bar or not.";
SetBarable::usage   = "Function for setting complex objects.";
UnsetBarable::usage = "Function for unsetting complex objects.";


CommutativeQ::usage      = "CommutativeQ[obj] returns True if obj is a comutative object and false otherwise.";
Commutative::usage       = "Commutative[obj] is an auxiliary head that can be used to treat the object obj temporarily as commutative.";


FermionTrace::usage        = "FermionTrace[expr] closes fermion lines in expr and performs Dirac gamma traces.";
CanonizeSpinorLines::usage = "CanonizeSpinorLines[expr] separates out nested spinor lines, breaking up the NCM accordingly.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Bar*)


(* ::Text:: *)
(*Bar denotes complex conjugation except for on fermion lines, where they are the Bar of the Dirac algebra*)


(* ::Subsubsection::Closed:: *)
(*Set rules for what is trivial under a Bar *)


(* ::Text:: *)
(*Default is no action*)


NonBarableQ[f_?NonBarableQ[x___]]:= AllTrue[{x}, NonBarableQ](*And@@ NonBarableQ/@ {x}*);
NonBarableQ[f_@ x___]:= False;

NonBarableQ[_]:= True (* this leads to wrong results in combinatino with Solve -> MapEffectiveCoupling now replaces all variables by a NonBarable[dummy]=False dummy *)
(* bug fix for solve with potential performance problems *)
(*
NonBarableQ[x_Symbol]:= If[StringMatchQ[Context[x], {"Solve`"~~___, "System`ReduceDump`"~~___, "Reduce`"~~___}], False, True]
NonBarableQ[x:Except[_Symbol]]:= True
*)

NonBarableQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold]= False;


(* ::Text:: *)
(*Objects with non-trivial behavior under Bar *)


NonBarableQ@ Alternatives[Complex, CG, Index]= False;
NonBarableQ@ _SparseArray= False;


NonBarableQ@ Alternatives[Coupling, DiracProduct, Field, FieldStrength,
	NCM, Transp, EvaOp]= False;


(* ::Subsubsection::Closed:: *)
(*General properties of Bar*)


Bar@ Bar@ x_:= x;


Bar@ x_?NonBarableQ:= x;


Bar[expr:Alternatives[_Plus, _Times, _Power, _Log, _List, _Operator]]:= Bar/@ expr;


Bar@ x_SparseArray:= Conjugate@ x;


Bar@ x_Complex:= Conjugate@ x;


(* ::Subsubsection::Closed:: *)
(*Action on Matchete symbols*)


Bar[expr_NCM]:= Bar/@ Reverse@ expr;
Bar@ x_DiracProduct:= Reverse@ x/. {Gamma5-> -Gamma5, Proj@ s_:> Proj[-s], g_GammaM:> Reverse@ g};


(*Go through list, as auto-ordering of GrassmanProduct can trigger immediately otherwise*)
Bar[expr_GrassmannProduct]:= GrassmannProduct@@ Bar/@ Reverse[List@@ expr];


Bar@ f:Field[label_, Scalar|Vector[_]|Ghost|AntiGhost, __]/; $FieldAssociation[label, SelfConjugate]:=f;
Bar@ f:Field[label_, Fermion, __]/; $FieldAssociation[label, SelfConjugate]:= 
	Transp[f]\[CenterDot] DiracProduct[GammaCC]; (* Bar for Majorana Fermions *)
Bar@ FieldStrength[label_, lind_, ginds_, cdinds_]/; $FieldAssociation[label, SelfConjugate]:= 
	FieldStrength[label, lind, Bar/@ ginds, cdinds];
Bar@ c:Coupling[label_, __]/; $CouplingAssociation[label, SelfConjugate]:= c;
Bar@ Coupling[label_, ind_List,x__]/; (Head[$CouplingAssociation[label, SelfConjugate]]===List):= Coupling[label, Permute[ind, $CouplingAssociation[label, SelfConjugate]],x];


(*Special case of charge conjugation before Bar*)
Transp/: Bar@ Transp@ f_:= Transp@ Bar@ f;


LF/: Bar@ lf_LF:= lf;
LCTensor/: Bar@ lc_LCTensor:= lc;


Bar@ Index[x_, Lorentz]:= Index[x, Lorentz];


(* ::Text:: *)
(*Properties dedicated to Finite Temperature*)


Bar[Index[\[Mu]_,Spatial]]:=Index[\[Mu],Spatial];


(* ::Subsubsection::Closed:: *)
(*For setting objects complex*)


SetBarable[x_, y__]  := (SetBarable@ x; SetBarable@ y;);
SetBarable[x_] := (NonBarableQ@ x ^= False;);


(* ::Text:: *)
(*And unseting them for removal*)


UnsetBarable[x_Symbol] := (x/: NonBarableQ@ x =.);
UnsetBarable[x:_Symbol@_] := (Evaluate@ Head@ x/: NonBarableQ@ x =.);


(* ::Section:: *)
(*NonCommutativeMultiply functionality*)


(* ::Subsubsection::Closed:: *)
(*Commutative check*)


(* ::Text:: *)
(*Hard-code commutative properties of objects to prevent iterative calls to all subparts: greatly improves performance *)


CommutativeQ@ _DiracProduct= False;
CommutativeQ@ Field[_, Fermion, __]= False;


CommutativeQ@ Field[_, Scalar| _Vector| Ghost| AntiGhost, __]= True;
CommutativeQ@ _FieldStrength= True;
CommutativeQ@ _Coupling= True;
CommutativeQ@ _CG= True;
CommutativeQ@ _LF= True;
CommutativeQ@ _MomExt= True;


(* ::Text:: *)
(*Default assumption is that everything is commutative*)


(* Memoization speeds up DeriveSubstitutions by a factor of ~2 or so *)
(* Memoization seems to have very poor performance in EOM and GreensSimplify creating 100,000s of Downvalues *)
(*CommutativeQ[x_NCM]:= (CommutativeQ[x]= ClosedSpinChainQ[x]);
CommutativeQ[f_?CommutativeQ[x___]]:= (CommutativeQ[f[x]]= And@@ CommutativeQ/@ {x});*)

CommutativeQ@ x_NCM:= ClosedSpinChainQ@ x;
CommutativeQ@ _GrassmannProduct= False;
CommutativeQ[f_?CommutativeQ[x___]] := AllTrue[{x}, CommutativeQ](*And@@ CommutativeQ/@ {x}*);
CommutativeQ[f_[x___]] := False;
CommutativeQ[_] := True;


(* ::Text:: *)
(*Define the non-commutative objects*)


CommutativeQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold] = False;


(* ::Text:: *)
(*Commutative is a head that can be used to temporarily treat an object as commutative*)


CommutativeQ@ Commutative@ _ ^= True;


(* ::Subsection:: *)
(*Properties of NCM*)


(* ::Subsubsection::Closed:: *)
(*Highjack CenterDot for NCM shorthand*)


CenterDot= NCM;
(*Protect@ CenterDot;*) (* This can cause the loading to crash! *)


(* ::Text:: *)
(*Set definitions for the non-commutative product*)


(* ::Subsubsection::Closed:: *)
(*General properties*)


(* ::Text:: *)
(*NCM contractions*)


NCM[a___, NCM[b__], c___]:= NCM[a, b, c] /;FreeQ[List@b,Fermion]
NCM[] = 1;
NCM@NCM@x___:=NCM@x;


(* Extracting commuting objects *)
NCM[a___, b_?CommutativeQ, c___]/; FreeQ[b,_Ext,All]:= b * NCM[a, c];
NCM[a___, b_?CommutativeQ * x_, c___]/; FreeQ[b,Field[_Ext,___],All]:= b * NCM[a, x, c];
(* Distributivity *)
NCM[a___, b_Plus, c___]:= NCM[a, #, c] & /@ b 


(* ::Subsubsection::Closed:: *)
(*Contraction of spin chains*)


NCM[a___, x:NCM@ b___, c___]:=
	If[ClosedSpinChainQ[x] && EvenQ@Count[x, _Ext, All],
		x NCM[a, c],
		NCM[a, b, c]
	];


(* ::Subsubsection::Closed:: *)
(*Canonical ordering of operators with transpose*)


(* ::Text:: *)
(*N.b. fermion patterns defined locally to ensure evaluation at read time of function definition (fails if defined in another file with lower read priority). *)


FermionL= Bar@ Field[_, Fermion, __] | Transp@ Field[_, Fermion, __] | EoM@ Bar@ Field[_, Fermion, __] | EoM@ Transp@ Field[_, Fermion, __];


FermionR= Field[_, Fermion, __] | Transp@ Bar@ Field[_, Fermion, __] | EoM@ Field[_, Fermion, __] | EoM@ Transp@ Bar@ Field[_, Fermion, __];


(* ::Text:: *)
(*Canonical ordering of fermion lines (excluding external legs -> this would break the UFO export)*)


NCM[f1:FermionL, g___DiracProduct, f2:FermionR]/; (FreeQ[{f1,f2}, _Ext, All] && SwitchFermionsQ[f1, f2]):= - NCM[Transp@ f2, Transp@ g, Transp@ f1];


(* ::Text:: *)
(*Determines if the spinor line should be transposed depending on the spinors.*)
(*N.b. Block over Module; it gets called frequently.*)


SwitchFermionsQ[f1_, f2_]:= Block[{l1, d1, l2, d2, r1, r2,
		transp1= !FreeQ[f1, Transp], transp2= !FreeQ[f2, Transp]},
	(*Always transpose if both spinors are transposed; never if none are*)
	If[Xnor[transp1, transp2],
		If[transp1 && transp2, 
			Return[True];
		,
			Return[False];
		];
	];
	
	(*Labels, derivatives, and selfconjugate*)
	{l1, d1}= FirstCase[f1, Field[l_, __, d_]:> {l, Length@ d}, {}, All];
	{l2, d2}= FirstCase[f2, Field[l_, __, d_]:> {l, Length@ d}, {}, All];
	{r1, r2}= {$FieldAssociation[l1, SelfConjugate], $FieldAssociation[l2, SelfConjugate]};
	
	Which[
	 !r1 && r2, (*Complex fermion w/o transp if possible*)
		transp1
	,!r2 && r1, 
		transp2
	,Head@ f2 === EoM && Head@ f1 =!= EoM, (*EoM on the right*)
		False
	,Head@ f1 === EoM && Head@ f2 =!= EoM,
		True
	,d2 > d1, (*More derivatives on the right*)
		False
	,d1 > d2,
		True
	,True, (*Defaults to canonical label ordering*)
		!OrderedQ@ {l1, l2}
	]
];


(* ::Subsubsection::Closed:: *)
(*Non-commutative symbols used in Matchete should be declared here*)


(* ::Text:: *)
(*From DiracAlgebra.m*)


CommutativeQ@ Alternatives[GammaM, Gamma5, GammaCC, Proj]= False;


(* ::Text:: *)
(*From FunctionalTools.m*)


(* SetNonCommutative[TransposeThisSpinChain] *)
(* For some reason the lie above does not work, but the one below does... *)
CommutativeQ[TransposeThisSpinChain]=False


(* ::Subsection:: *)
(*Properties of GrassmannProduct*)


(* ::Subsubsection::Closed:: *)
(*GrassmannCommuteQ \[LongDash] commutative predicate for use inside GrassmannProduct*)


(* ::Text:: *)
(*Mirrors CommutativeQ and FuncCommuteQ. Returns False for all Grassmann-valued fields (Fermion, Ghost, AntiGhost), open NCM spin chains, and pattern objects. Returns True for closed spin chains, couplings, scalars, etc.*)


GrassmannCommuteQ@ Field[Except[_Ext], Fermion|Ghost|AntiGhost, __]= False;
GrassmannCommuteQ@ Field[_Ext, __]= True;
GrassmannCommuteQ@ x_NCM:= ClosedSpinChainQ[x] && EvenQ@Count[x, _Ext, All];
GrassmannCommuteQ@ Alternatives[Pattern, Blank, BlankSequence, BlankNullSequence, Except, Hold]= False;
GrassmannCommuteQ@ Commutative@ _ ^= True;
GrassmannCommuteQ[f_?GrassmannCommuteQ[x___]]:= And@@ GrassmannCommuteQ/@ {x};
GrassmannCommuteQ[f_[x___]]:= False;
GrassmannCommuteQ[_]:= True;


(* ::Subsubsection::Closed:: *)
(*General properties*)


(* Flatness *)
GrassmannProduct[a___, GrassmannProduct[b__], c___]:= GrassmannProduct[a, b, c]
GrassmannProduct[]:= 1;

(* Extract GrassmannCommuteQ-commutative objects (scalars, closed spin chains, ...) *)
GrassmannProduct[a___, b_?GrassmannCommuteQ, c___]:= b * GrassmannProduct[a, c];
GrassmannProduct[a___, b_?GrassmannCommuteQ * x_, c___]:= b * GrassmannProduct[a, x, c];

(* Distributivity *)
GrassmannProduct[a___, b_Plus, c___]:= GrassmannProduct[a, #, c] &/@ b

(* Single commuting argument collapses; lone Grassmann field stays wrapped *)
GrassmannProduct[x_?GrassmannCommuteQ]:= x


(* ::Subsubsection::Closed:: *)
(*Canonical ordering: AntiGhost left, Ghost right*)


(* When a Ghost appears to the left of an AntiGhost, swap them with a minus sign. 
	The pattern is recursive but terminates: after the AntiGhost always moves left. *)
GrassmannProduct[a___, gh: Field[_, Ghost, ___]|Bar@ Field[_, Ghost, ___], 
	agh: Field[_, AntiGhost, ___]|Bar@ Field[_, AntiGhost, ___], d___]:=
	-GrassmannProduct[a, agh, gh, d];


(* ::Section:: *)
(*Active manipulation*)


(* ::Subsection:: *)
(*Organize fermion lines*)


(* ::Subsubsection::Closed:: *)
(*Fermion trace*)


(* ::Text:: *)
(*Function for closing fermion lines in traces as appear in some supertraces. *)


FermionTrace::cantclose = "The fermion line cannot be closed";
FermionTrace::multopens = "Multiple open fermion lines encountered in one term.";


OpenSpinChainQ= Not@* ClosedSpinChainQ;


FermionTrace@ expr_:= Module[{out= LagrangianExpand@ expr},	
	(*Trace terms in the sums separately*)
	If[Head@ out === Plus, Return[FermionTrace/@ out];]; 
	
	Switch[Count[{out}, _NCM? OpenSpinChainQ, All]
		,0 ,
			4 out
		,1 ,
			out/. x_NCM? OpenSpinChainQ:> TraceSpin@ x// ContractMetricSingleTerm
		,_ ,
			Message[FermionTrace::multopens];
			Abort[];
	]
];


TraceSpin@ spinChain_NCM:= Module[{spinCount= CumulativeSpin@ spinChain},
	If[Last@ spinCount =!= 0,
		Message[FermionTrace::cantclose];
		Abort[];	
	];
	(*Check if the fermion line can be closed through the trace *)
	If[MemberQ[spinCount, -1],
		(*Sign is from Grassmanian variables*)
		Return[-RotateLeft[spinChain, FirstPosition[spinCount, -1]]];
	];
	(*Otherwise do a Dirac trace*)
	DiracTrace@ spinChain
];


(* ::Subsubsection::Closed:: *)
(*Spinor line canonization*)


(* ::Text:: *)
(*Function for separating out nested spinor lines, breaking up the NCM accordingly*)


CanonizeSpinorLines::lopens = "The Spinor line has >1 L-open spinors.";
CanonizeSpinorLines::ropens = "The Spinor line has >1 R-open spinors.";


CanonizeSpinorLines@ fline_NCM:= Module[
		{cumSpinCount, max, min, l, r, expr},
	cumSpinCount= CumulativeSpin@ fline;
	min= Min@ cumSpinCount;
	If[min < -1,
		Message[CanonizeSpinorLines::lopens];
		Abort[];
	];
	If[Last@ cumSpinCount> 1,
		Message[CanonizeSpinorLines::ropens];
		Abort[];
	];
	
	If[min ===-1,
		cumSpinCount[[;; First@ FirstPosition[cumSpinCount, -1] -1]]--;
	]; 
	max= Max@ cumSpinCount;
		
	(*If a closed spinor line appear it is separated out*)
	{l, r}= First@ SequencePosition[cumSpinCount, {max..}, Overlaps-> False]+ {0, 1};
	If[r <= Length@ fline,
		Return[fline[[l;;r]] CanonizeSpinorLines[fline[[;; l-1]]\[CenterDot] fline[[r+ 1;;]] ] ];
	];

	fline
];


CanonizeSpinorLines@ expr_:= expr//. x_NCM:> CanonizeSpinorLines@ x;


(* ::Subsubsection::Closed:: *)
(*Vector-like spin chain*)


(* ::Text:: *)
(*Tests whether a spin chain is vector-like, i.e., contains two vector-like fermions. (If even one fermion is chiral a projector should always be present).*)


VectorLikeSpinChainQ@ NCM[
	Bar@ Field[l1_, Fermion, __], ___, Field[l2_, Fermion, __]]:= VectorLikesQ[l1, l2];
VectorLikeSpinChainQ@ NCM[
	Transp@ Field[l1_, Fermion, __], ___, Field[l2_, Fermion, __]]:= VectorLikesQ[l1, l2];
VectorLikeSpinChainQ@ NCM[
	Bar@ Field[l1_, Fermion, __], ___, Transp@ Bar@ Field[l2_, Fermion, __]]:= VectorLikesQ[l1, l2];
VectorLikeSpinChainQ@ NCM[
	Transp@ Field[l1_, Fermion, __], ___, Transp@ Bar@ Field[l2_, Fermion, __]]:= VectorLikesQ[l1, l2];
VectorLikeSpinChainQ@ ncm_:= If[FreeQ[ncm, EoM], 
		False
	, 
		VectorLikeSpinChainQ[ncm/.EoM-> Identity]
	];


VectorLikesQ[l1_, l2_]:= $FieldAssociation[l1, Chiral] === False && $FieldAssociation[l2, Chiral] === False;
