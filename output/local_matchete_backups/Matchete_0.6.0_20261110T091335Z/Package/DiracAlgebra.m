(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`DiracAlgebra`*)


(* ::Subtitle:: *)
(*Implementation of the Dirac algebra used in Matchete.*)


(* ::Text:: *)
(*\[Gamma] matrices are all decomposed into the standard form (see e.g. https://arxiv.org/abs/physics/0703214, eq.  4.15) *)
(*DiracProduct[...,Gamma5,...] is interpreted as ...\[Gamma]_5... *)
(*DiracProduct[...,GammaM[Index[\[Alpha],Lorentz]],...] is interpreted as ...\[Gamma]_\[Alpha]...*)
(*GammaM[\[Mu]1, ..., \[Mu]n] denotes \[CapitalGamma]_n^{\[Mu]1, ..., \[Mu]n} = 1/n! (\[Gamma]_\[Mu]1 \[CenterDot] ... \[CenterDot] \[Gamma]_\[Mu]n \[PlusMinus] ... \[PlusMinus] \[Gamma]_\[Mu]n \[CenterDot] ... \[CenterDot] \[Gamma]_\[Mu]1 ), the normalized anti-symmetric sum of n gamma matrices.  *)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping:*)


(* ::Subsubsection::Closed:: *)
(*Input variables*)


PackageExport["\[Gamma]"]
PackageExport["\[Sigma]"]


PackageExport["PL"]
PackageExport["PR"]


PackageExport["CC"]


PackageExport["Transp"]
PackageExport["CConj"]


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["DiracProduct"]


PackageExport["GammaM"]
PackageExport["Gamma5"]
PackageExport["GammaCC"]
PackageExport["Proj"]


PackageExport["DiracTrace"]


PackageExport["GammaReduction"]
PackageExport["Fierz"]
PackageExport["RefineDiracProducts"]
PackageExport["ExpandGammas"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["FermionQ"]
PackageScope["ClosedSpinChainQ"]
PackageScope["LOpenSpinChainQ"]
PackageScope["ROpenSpinChainQ"]
PackageScope["MajoranaQ"]
PackageScope["CumulativeSpin"]


PackageScope["ProjExpand"]


PackageScope["CollectGammaMatrices"]


PackageScope["ASymGammaExpand"]


PackageScope["\[Gamma]CommuteQ"]


PackageScope["FierzScore"]
PackageScope["FbasisQ"]


PackageScope["LC2Gamma5"]
PackageScope["Gamma52LC"]
PackageScope["LC2Delta"]
PackageScope["SchoutenID"]


PackageScope["GammaReductionInternal"]


(* ::Section:: *)
(*Usage definitions*)


(* ::Subsubsection::Closed:: *)
(*Input variables*)


\[Gamma]::usage = "\[Gamma][\[Mu]] is the shortcut for a Dirac matrix with Lorentz index \[Mu].
\[Gamma][5] is the shortcut for the Dirac matrix \!\(\*FormBox[SubscriptBox[\(\[Gamma]\), \(5\)],TraditionalForm]\).
\[Gamma][\!\(\*SubscriptBox[\(\[Mu]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Mu]\), \(n\)]\)] is the shortcut for the anti-symmetrized product of n Dirac matrices \!\(\*SubscriptBox[\(\[Gamma]\), \([\*SubscriptBox[\(\[Mu]\), \(1,  ... , \)] \*SubscriptBox[\(\[Mu]\), \(n\)]]\)]\)";
\[Sigma]::usage = "\[Sigma][\[Mu],\[Nu]] is the shortcut for the hermitian anti-symmetric Dirac matrix combination \[Sigma][\[Mu],\[Nu]] = \!\(\*FractionBox[\(\[ImaginaryI]\), \(2\)]\) [\[Gamma][\[Mu]],\[Gamma][\[Nu]]].";


PL::usage = "PL is the shortcut for the chiral left projector \!\(\*FormBox[SubscriptBox[\(P\), \(L\)],TraditionalForm]\) in a Dirac product.";
PR::usage = "PR is the shortcut for the chiral right projector \!\(\*FormBox[SubscriptBox[\(P\), \(R\)],TraditionalForm]\) in a Dirac product.";


CC::usage = "CC is the shortcut for the charge conjugation matrix in a Dirac product."


Transp::usage = "Transp[arg] is used to denote the transposed of the object arg in the Dirac algebra.";
CConj::usage  = "CConj[field] is the charge conjugation operator and can only be used to act on fermion fields, for example CConj[\[Psi]]= \!\(\*SuperscriptBox[\(\[Psi]\), \(c\)]\)= \!\(\*TemplateBox[{\"C\"},\n\"RowDefault\"]\) \!\(\*SuperscriptBox[OverscriptBox[\(\[Psi]\), \(_\)], \(T\)]\).";


(* ::Subsubsection::Closed:: *)
(*Exported*)


DiracProduct::usage = "DiracProduct[seq] represents the (non-commutative) product of Dirac matrices, charge conjugation matrices and chiral projectors in the sequence seq.";


GammaM::usage  = "GammaM[\[Mu]] is a Dirac matrix with Lorentz index \[Mu] as written inside DiracProduct.
GammaM[\!\(\*FormBox[SubscriptBox[\(\[Mu]\), \(1\)],TraditionalForm]\), ...,\!\(\*FormBox[SubscriptBox[\(\[Mu]\), \(n\)],TraditionalForm]\)] is interpreted as  \!\(\*FormBox[SubscriptBox[SuperscriptBox[\(\[CapitalGamma]\), \(n\)], \(\(\*SubscriptBox[\(\[Mu]\), \(1\)] ... \) \*SubscriptBox[\(\[Mu]\), \(n\)]\)],TraditionalForm]\)=  \!\(\*FormBox[\(\(\*SubscriptBox[\(\[Gamma]\), \([\(\)\*SubscriptBox[\(\[Mu]\), \(1\)]\)] ... \)\\\ \*SubscriptBox[\(\[Gamma]\), \(\*SubscriptBox[\(\[Mu]\), \(n\)]\(]\)\)]\),TraditionalForm]\)  (full antisymmetrization with normalization factor 1/n!).
The associated shortcut is \[Gamma][\[Mu]] or \[Gamma][\!\(\*SubscriptBox[\(\[Mu]\), \(1\)]\),...,\!\(\*SubscriptBox[\(\[Mu]\), \(n\)]\)].";(*GammaM[LoopMom] represents p slash.*)
Gamma5::usage  = "Gamma5 is the Dirac matrix \!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\) as written inside DiracProduct.
In the NDR scheme, anti-commutativity with other gamma matrices is preserved but cyclicity of the Dirac trace is abandonned, adding a consistency requirement in the reading point choice.
The associated shortcut is \[Gamma][5].";
GammaCC::usage = "GammaCC is the charge conjugation matrix as written inside DiracProduct.
The associated shortcut is CC."
Proj::usage    = "Proj[+/-1] is the chiral right/left projectors \!\(\*SubscriptBox[\(P\), \(L/R\)]\) = \!\(\*FractionBox[StyleBox[RowBox[{\"1\", \" \", \"+\", RowBox[{\"/\", RowBox[{\"-\", \" \", SubscriptBox[\"\[Gamma]\", \"5\"]}]}]}], \"Code\"], \(2\)]\) as written inside DiracProduct.
The associated shortcut is PR/PL.";


DiracTrace::usage = "DiracTrace[expr] evaluates the Dirac trace of expr.";


RefineDiracProducts::usage = "RefineDiracProducts[expr] simplifies all Dirac product and matches it into an expression to the basis of completely anti-symmetrized products.";
ExpandGammas::usage = "ExpandGammas[expr] Expands all anti-symmetric gamma matrices in an expression into individual singe-index gamma matrices.";


GammaReduction::usage            = "GammaReduction[SP1 * SP2] returns the simplification of SP1 and SP2 current after a projection of the Dirac structures to the 4D basis.";
Fierz::usage                     = "Fierz[SP1 * SP2] returns the fierzing of the product of SP1 and SP2 by exchanging the second spinor in each closed spin chain."


(* ::Subsubsection::Closed:: *)
(*Internal*)


FermionQ::usage         = "FermionQ[expr] returns True if a expr is a fermion or False otherwise.";
ClosedSpinChainQ::usage = "ClosedSpinChainQ[expr] returns True if a expr is a closed spin chain or False otherwise.";
LOpenSpinChainQ::usage  = "LOpenSpinChainQ[expr] returns True if a expr is a spin chain that is open from the left or False otherwise.";
ROpenSpinChainQ::usage  = "ROpenSpinChainQ[expr] returns True if a expr is a spin chain that is open from the right or False otherwise.";
MajoranaQ::usage        = "FermionQ[expr] returns True if a expr is a Majorana fermion or False otherwise.";
CumulativeSpin::usage   = "Returns the cumulative count of fermion spins along a NCM product. Counts +1 for R-open fermions and -1 for L-open.";


CollectGammaMatrices::usage = "CollectGammaMatrices[expr] collects all \[Gamma]-matrices in a DiracProduct and and moves them to the left inside of a NCM.";


ProjExpand::usage      = "Expand the left/right handed projector into (1-/+\!\(\*SubscriptBox[\(\[Gamma]\), \(5\)]\))/2."


ASymGammaExpand::usage = "ASymGammaExpand[DiracProduct[...]] performs the anti-symmetric expansion of the \[Gamma]-matrices inside the DiracProduct";


\[Gamma]CommuteQ::usage = "\[Gamma]CommuteQ[obj] checks if obj commutes with gamma matrices.";


LC2Gamma5::usage= "LC2Gamma5[expr] combines \[Epsilon]^\[Mu]\[Nu]\[Rho]\[Sigma] with gamma matrices in expr to obtain Subscript[\[Gamma], 5], using the relation Subscript[\[Gamma], 5] = -I/4! \[Epsilon]^\[Mu]\[Nu]\[Rho]\[Sigma] \[Gamma][\[Mu]]\[CenterDot] \[Gamma][\[Nu]]\[CenterDot] \[Gamma][\[Rho]]\[CenterDot] \[Gamma][\[Sigma]].";
Gamma52LC::usage= "Gamma52LC[expr] decomposes Subscript[\[Gamma], 5] in expr as Subscript[\[Gamma], 5] = -I/4! \[Epsilon]^\[Mu]\[Nu]\[Rho]\[Sigma] \[Gamma][\[Mu]]\[CenterDot] \[Gamma][\[Nu]]\[CenterDot] \[Gamma][\[Rho]]\[CenterDot] \[Gamma][\[Sigma]].";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*SNDR scheme*)


(* ::Subsubsection::Closed:: *)
(*Convenient notation*)


PL= DiracProduct@ Proj@ -1;
PR= DiracProduct@ Proj@ 1;


\[Gamma]@ Sequence[\[Mu]_Symbol,\[Nu]___]:= DiracProduct@ GammaM@ Sequence[Index[\[Mu], Lorentz], Sequence@@(Index[#,Lorentz]&/@List[\[Nu]])];
\[Gamma]@ Sequence[\[Mu]: Index[_Symbol, Lorentz],\[Nu]___]:= DiracProduct@ GammaM@ Sequence[\[Mu],\[Nu]];
\[Gamma]@ 5:= DiracProduct@ Gamma5;
\[Sigma]@ Sequence[\[Mu]_Symbol,\[Nu]_Symbol]:= I DiracProduct@ GammaM@ Sequence[Index[\[Mu], Lorentz],Index[\[Nu], Lorentz]];
\[Sigma]@ Sequence[\[Mu]: Index[_Symbol, Lorentz],\[Nu]: Index[_Symbol, Lorentz]]:= I DiracProduct@ GammaM@ Sequence[\[Mu],\[Nu]];


CC = DiracProduct@ GammaCC;


(* ::Subsection:: *)
(*Properties of DiracProduct*)


(* ::Subsubsection::Closed:: *)
(*General properties*)


DiracProduct/: NCM[a___, x_DiracProduct, y_DiracProduct, b___]:=
	NCM[a, Join[x, y], b];


DiracProduct/: MomNCM[a___, x_DiracProduct, y_DiracProduct, b___]:=
	MomNCM[a, Join[x, y], b];


DiracProduct[a___, GammaM[] | Transp@ GammaM[], b___] := DiracProduct[a, b];
DiracProduct[a___, c_?NumberQ d_, b___]:= c DiracProduct[a, d, b];
DiracProduct[___,0,___]= 0;
DiracProduct[]= 1;


DiracProduct[A___,1,B___]:=DiracProduct[A,B];


(* ::Subsubsection::Closed:: *)
(*Properties of antisymmetrized \[Gamma] matrices*)


(* ::Text:: *)
(*Anti-symmetric in the indices.*)


GammaM@ OrderlessPatternSequence[\[Mu]_, \[Mu]_, ___]:= 0;


GammaM[a_, b__]/; !OrderedQ@ {a, b}:= Signature@{a, b} GammaM @@ Sort@{a, b}


(* ::Subsubsection::Closed:: *)
(*SNDR specific properties*)


(* ::Text:: *)
(*Subscript[\[Gamma], 5] : contract pairs (or with projectors) and move singles to the right.*)


DiracProduct[a___, Gamma5 | Transp@Gamma5 , b___/; FreeQ[List@ b, Gamma5| Transp@Gamma5| GammaCC| Proj| Transp@Proj], Gamma5 | Transp@Gamma5 , c___]:=
	(-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List] * DiracProduct[a, b, c];
DiracProduct[a___, Gamma5 | Transp@Gamma5, b___/; FreeQ[List@ b, Gamma5 | Transp@Gamma5| GammaCC| Proj| Transp@Proj], p: Proj@s_ | Transp@Proj@s_ , c___]:=
	s (-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List] * DiracProduct[a, b, p, c];
DiracProduct[a___, p: Proj@s_ | Transp@Proj@s_, b___/; FreeQ[List@ b, Gamma5 | Transp@Gamma5| GammaCC| Proj| Transp@Proj], Gamma5 | Transp@Gamma5, c___]:=
	s (-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List] * DiracProduct[a, p, b, c];
DiracProduct[a___, g5: Gamma5 | Transp@Gamma5, b__/; FreeQ[List@ b, Gamma5 | Transp@Gamma5| GammaCC| Proj| Transp@Proj]]:=
	(-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List] * DiracProduct[a, b, g5];


(* ::Text:: *)
(*Projectors: contract pairs and move singles to the right.*)


DiracProduct[a___, Proj@ r_ , b___/; FreeQ[List@ b, Gamma5| Transp@ Gamma5| GammaCC| Proj],  Proj@s_, c___]:=
	If[r s === (-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List], DiracProduct[a, b, Proj@ s, c], 0];
DiracProduct[a___, Proj@ s_, b__/; FreeQ[List@ b, Gamma5| Transp@ Gamma5| GammaCC| Proj]]:=
	DiracProduct[a, b, Proj[s (-1)^Length@ Flatten[List@ b/. Transp-> Identity/. GammaM-> List] ] ];


(* ::Text:: *)
(*Subscript[\[Gamma], \[Mu]] : Contract adjacent gamma matrices with common indices.*)


DiracProduct[a___, GammaM@ \[Mu]_Index , GammaM@\[Mu]_Index, b___]:= \[ScriptD] DiracProduct[a, b];


DiracProduct[a___, Transp@GammaM@ \[Mu]_Index, Transp@GammaM@ \[Mu]_Index, b___]:= \[ScriptD] DiracProduct[a, b];


(* ::Text:: *)
(*p slash :  combine two into a propagator (using use {\[Gamma]_\[Mu], \[Gamma]_\[Nu]}= 2 g_{\[Mu]\[Nu]} if needed).*)


DiracProduct[a___, GammaM@ LoopMom, GammaM@ LoopMom, b___]:= Power[Prop[0], -1] DiracProduct[a, b];
DiracProduct[a___, gL: GammaM@ LoopMom, gn: GammaM@ \[Nu]_Index, b___, GammaM@ LoopMom, c___]:=
	2 LoopMom@ \[Nu] * DiracProduct[a, b, gL, c] - DiracProduct[a, gn, gL, b, gL, c];


DiracProduct[a___, Transp@ GammaM@ LoopMom, Transp@ GammaM@ LoopMom, b___]:= Power[Prop[0], -1] DiracProduct[a, b];
DiracProduct[a___, gL: Transp@ GammaM@ LoopMom, gn: Transp@ GammaM@ \[Nu]_Index, b___, Transp@ GammaM@ LoopMom, c___]:=
	2 LoopMom@ \[Nu] * DiracProduct[a, b, gL, c] - DiracProduct[a, gn, gL, b, gL, c];


(* ::Subsubsection::Closed:: *)
(*4D function to remove Levi-Civita tensor*)


(* ::Text:: *)
(*Contract Levi-Civita tensor into the Dirac algebra to obtain a \[Gamma]5:*)
(*	\[CurlyEpsilon]^{\[Mu]\[Nu]\[Rho]\[Sigma]} \[Gamma]_\[Sigma] = -i \[CapitalGamma]^{\[Mu]\[Nu]\[Rho]} \[Gamma]5*)
(*This agrees with the conventions in Peskin (3.68). *)
(*The identity generalizes to anti-symmetrized gammas by (note anti symmetrization of the \[Alpha]i indices enclosed by [...])*)
(*	\[CurlyEpsilon]_{\[Mu]\[Nu]\[Rho]\[Sigma]} \[CapitalGamma]^{\[Sigma] \[Alpha]1... \[Alpha]n} = -i/(n+1) \[Sum]_{k=0}^{n} (-1)^k \[CapitalGamma]^{ [\[Alpha]1... \[Alpha]k } \[CapitalGamma]_{\[Mu]\[Nu]\[Rho]} \[Gamma]5 \[CapitalGamma]^{ \[Alpha](k+1)...\[Alpha]n] }*)


LC2Gamma5[expr_Plus]:=LC2Gamma5/@expr;


LC2Gamma5@ expr_:= Module[
		{gammaInds, gammaPos, gammaPositions, LCinds, LCtens, remainingInds, repInd, setk, setnk, sign, out, temp, n, k, TranspQ},
	(*Find the indices of the LC tensor*)
	LCtens= Cases[expr, _LCTensor, All];
	If[Length@ LCtens =!= 1, Return@ expr; ];
	LCinds= List@@ First@ LCtens;
	(*Find the positions and indices of candidate gammas*)
	gammaPositions= Position[expr, g_GammaM/; IntersectingQ[List@@ g, LCinds], All];
	If[gammaPositions === {}, Return@ expr; ];

	(*Only the gammas with the highest number of shared indices*)
	gammaPositions= MaximalBy[gammaPositions, Length@ Intersection[LCinds, List@@ Extract[expr, #]] &];
	out= 1/Length[gammaPositions]* Sum[
			gammaInds= List@@ Extract[expr, gammaPos];
			(*Find (a) repeated index and determine the sign
			(note indices of LCTensor + GammaM are automatically canonized)*)
			repInd= First@ Intersection[LCinds, gammaInds];
			sign= - Power[-1, Total@ Position[{LCinds, gammaInds}, repInd][[;;, 2]]];
			gammaInds= DeleteCases[gammaInds, repInd];
			remainingInds= DeleteCases[LCinds, repInd];
			n= Length@ gammaInds;

			(*Performance will likely be improved by explicit handling of situation with multiple contracted indices*)
			(*Perform reduction*)
			-I* sign/ (n+ 1)* Sum[
					setnk= Complement[gammaInds, setk];
					temp= expr;
					TranspQ= (Head@temp[[Sequence@@ Drop[gammaPos,-1]]]===Transp);
					If[!TranspQ,
						temp[[Sequence@@ gammaPos]]= Sequence[GammaM@@ setk, GammaM@@ remainingInds, Gamma5, GammaM@@ setnk],

						temp[[Sequence@@ Drop[gammaPos,-1]]]=  Sequence[Transp@* GammaM@@ setnk, Transp@ Gamma5, Transp@* GammaM@@ remainingInds, Transp@* GammaM@@ setk]];
					Power[-1, k]/ Binomial[n, k]* Signature@ setk* Signature@ setnk* Signature@ Join[setk, setnk]* temp
				, {k, 0, n}, {setk, Subsets[gammaInds, {k}]}]/. _LCTensor-> 1// RefineDiracProducts// ContractMetric
		, {gammaPos, gammaPositions}];

	(*The identities are four-dimensional so we should set \[ScriptD]-> 4*)
	out/. \[ScriptD]-> 4// RelabelIndices
]


(* ::Text:: *)
(*Write \[Gamma]5 with a Levi-Cevita tensor*)


Gamma52LC[exp_]:=exp/.DiracProduct[dp___,Gamma5]:> -I/4! LCTensor[\[Alpha],\[Beta],\[Delta],\[Eta]]DiracProduct[dp]\[CenterDot] \[Gamma][\[Alpha]]\[CenterDot] \[Gamma][\[Beta]]\[CenterDot] \[Gamma][\[Delta]]\[CenterDot] \[Gamma][\[Eta]] //RelabelIndices ;


(* ::Text:: *)
(*Expand product of two Levi-Civita tensors into antisymmetric product of metric*)


LC2Delta[exp_]:= exp/.{
		LCTensor[a:Sequence[_Index..]] * LCTensor[b:Sequence[_Index..]] :> If[FreeQ[exp, Lorentz], 1 , -1] Det[Outer[Metric,List@a,List@b]],
		_LCTensor^2-> If[FreeQ[exp, Lorentz], ($Dimensions-1)! , - $Dimensions!]
	};


(* ::Text:: *)
(*Perform cyclic permutations of indices on the product of a Levi-Civita tensor with a metric*)


SchoutenID[exp_]:=Module[{tempind, Tind, LCt, LCind, LMind, newLind, out},
	If[Length[LCt=Cases[exp,_LCTensor ,Infinity]]=!=1,
		exp,
		LCt=First@LCt;
		LCind=List@@LCt;
		LMind=Cases[exp,Index[_,Lorentz],Infinity];
		LMind=Part[Select[Tally@LMind,Part[#,2]===2&],All,1];
		LMind=Complement[LMind, LCind];
		If[LMind==={},
			exp,
			(*If several indices, another than First might be prefered*)
			LMind=First@LMind;
			out= ReleaseHold[exp/.LMind->Index[HoldForm@Unique[],tempind]];
			Tind=Cases[out, Index[_,tempind],Infinity];
			newLind=Index[Unique[d], Lorentz];
			out= out/.Tind[[1]]->LMind/.Tind[[2]]->newLind;
			out= - Metric[LCind[[1]], LMind] (out /.{LCt->LCTensor[LCind[[2]],LCind[[3]],LCind[[4]],newLind]})
				 - Metric[LCind[[2]], LMind] (out /.{LCt->LCTensor[LCind[[3]],LCind[[4]],newLind,LCind[[1]]]})
				 - Metric[LCind[[3]], LMind] (out /.{LCt->LCTensor[LCind[[4]],newLind,LCind[[1]],LCind[[2]]]})
				 - Metric[LCind[[4]], LMind] (out /.{LCt->LCTensor[newLind, LCind[[1]],LCind[[2]],LCind[[3]]]});
			out //RelabelIndices//Contract
			]
		]
	]


(* ::Subsection:: *)
(*Dirac trace*)


(* ::Subsubsection::Closed:: *)
(*Trace of an expression*)


(* ::Text:: *)
(*Move the DiracTrace from acting on any expression to act on the DiracProduct inside the expression*)


Options@ DiracTrace= {Dimensions->\[ScriptD]};


DiracTrace[expr_, opt:OptionsPattern[]]:= Module[
		{allterms= TermsToList@ expr, term},
	Sum[
		Switch[Count[term, DiracProduct, Infinity, Heads-> True],
			0, 4 term,
			(*Prevent the result of the Diractrace from being inserted inside NCM,
				which is an expensive operation.*)
			1, term/.{NCM@ x_DiracProduct:> DiracTraceProduct@ x,
				x_DiracProduct:> DiracTraceProduct@ x},
			_, Message[DiracTrace::canteval]; term
		],
		{term, allterms}
	]/.\[ScriptD]-> OptionValue@ Dimensions
];


(* ::Text:: *)
(*Error Messages*)


DiracTrace::canteval= "Dirac trace encountered term with multiple Dirac products. No trace is applied.";


(* ::Subsubsection::Closed:: *)
(*Split of anti-symmetrized \[Gamma] matrices *)


(* ::Text:: *)
(*Anti-symmetric expansion of the \[Gamma]-matrices inside the DiracProduct*)


ASymGammaExpand@ DiracProduct[a___, GammaM[\[Mu]___]/; Length@{\[Mu]}>1, opt___]:=Module[{permutations= Permutations@{\[Mu]}},
	Signature@ {\[Mu]}/Length @permutations *
		Plus@@((Signature@ # * ASymGammaExpand@ DiracProduct[a, Sequence@@ (GammaM/@ #), opt]&)/@ permutations)//
	LagrangianExpand];
ASymGammaExpand@ x_= x;


(* ::Text:: *)
(*Exported expansion function for all Gamma matrices in an expression*)


ExpandGammas@ expr_:= expr/. d_DiracProduct:> ASymGammaExpand@ d;


(* ::Subsubsection::Closed:: *)
(*Dirac trace*)


(* ::Text:: *)
(*Expand projections inside DiracTraces*)


DiracTraceProduct@ 1= 4;


DiracTraceProduct[DiracProduct[a___, Proj@ s_]]:=
	DiracTraceProduct@ DiracProduct@ a/2+ s DiracTraceProduct@ DiracProduct[a, Gamma5]/2 ;


(* ::Text:: *)
(*Caching of traces*)


DiracTraceProduct@ prod:DiracProduct[___, _GammaM| Gamma5]:= Module[
		{len=Total[Length/@ List@@ prod]},
	If[OddQ@ len, Return@ 0];

	(*Cache expressions with less than 7 open indices*)
	If[2* Length@ DeleteDuplicates@ Cases[prod, _Index, {2}] - len < 7,
		Return@ With[{xPats= prod/. Index[a_, Lorentz]:> Index[Pattern[#, Blank[]]&@ a, Lorentz]},
			DiracTraceProduct@ xPats= ContractMetric@ DiracTraceEvaluate@ prod
		];
	];

	DiracTraceEvaluate@ prod
];


(* ::Text:: *)
(*Traces w/o gamma5*)


DiracTraceEvaluate@ prod:DiracProduct[_GammaM...]:= Module[
		{indPos, pairCombs},
	indPos= Flatten[Table[{m, n}, {m, Length@ prod}, {n, Length@ prod[[m]]}], 1];
	(*If[OddQ@ Length@ indPos, Return@ 0;];*)

	(*No metric with indices from the same GammaM: these indices are anti-symmetric*)
	pairCombs= NonOverlappingPairs@ indPos;
	If[!FreeQ[prod, GammaM[_, __]],
		pairCombs= DeleteCases[pairCombs,
			List@ OrderlessPatternSequence[{{a_, _}, {a_, _}}, ___]];
	];

	Sum[
		4* Signature[Join@@ pairs]* Product[Metric@@ Extract[p]@ prod, {p, pairs}]
	, {pairs, pairCombs}]
]


(* ::Text:: *)
(*Traces w/ gamma5*)
(*Inspired by https://inspirehep.net/literature/318268, we substitute the Gamma5 for Levi-Civita contracted with \[CapitalGamma]_4^{mu1,... mu4}, then refine the product of all remaining gammas. The Clifford algebra basis of \[CapitalGamma]_n is orthogonal in the trace. *)


DiracTraceEvaluate@ prod:DiracProduct[_GammaM..., Gamma5]:= Module[{aSym, indPos, len, pairCombs},
	indPos= Flatten[Table[{m, n}, {m, Length@ prod- 1}, {n, Length@ prod[[m]]}], 1];
	len= Length@ indPos;
	Which[(*OddQ@ len ||*) len < 4,
		Return@ 0;
	,len === 4,
		Return[-4I * LCTensor@@ Extract[indPos]@ prod];
	];

	Sum[
		aSym= Complement[indPos, set];
		(*No metric with indices from the same GammaM: these indices are anti-symmetric*)
		pairCombs= DeleteCases[NonOverlappingPairs@ set,
			List@ OrderlessPatternSequence[{{a_, _}, {a_, _}}, ___]];
		Sum[
			-4I* Signature@ Join[Sequence@@ pairs, aSym]* Product[Metric@@ Extract[p]@ prod, {p, pairs}]*
			LCTensor@@ Extract[aSym]@ prod
		, {pairs, pairCombs}]
	, {set, Subsets[indPos, {len - 4}]}]
]


(* ::Subsection:: *)
(*Put Dirac Products to basis*)


(* ::Subsubsection::Closed:: *)
(*RefineDiracProducts*)


(* ::Text:: *)
(*Anti-symmetrize all the DiracProducts inside DiracTraces*)


RefineDiracProducts@ expr_:= expr/. x_DiracProduct:> RefineDiracProductsInternal@ x// ContractMetric;


(* ::Text:: *)
(*Cannot simplify Dirac products with a single factor*)


RefineDiracProductsInternal@ d:DiracProduct[_]:= d;


(* ::Text:: *)
(*Charge conjugation matrices, chiral projectors, and Subscript[\[Gamma], 5] are left  alone  *)


RefineDiracProductsInternal@ DiracProduct[GammaCC, x__]:=
	DiracProduct@ GammaCC\[CenterDot] RefineDiracProductsInternal@ DiracProduct@x;


RefineDiracProductsInternal@ DiracProduct[x__, chiral: (Gamma5| _Proj)]:=
	RefineDiracProductsInternal@ DiracProduct@ x \[CenterDot] DiracProduct@ chiral;


(* ::Text:: *)
(*Expedite evaluation when multiple Lorentz indices are contracted (performance might be improved by including the more advanced methods of [1905.00429]. *)


RefineDiracProductsInternal[DiracProduct[a___, gm: GammaM@ \[Mu]_Index, gn: GammaM@\[Nu]_Index, b___, GammaM@ \[Mu]_Index, c___]/; FreeQ[List[b, c], \[Nu]]]:=
	RefineDiracProducts[2 Metric[\[Mu], \[Nu]] DiracProduct[a, b, gm, c]- DiracProduct[a, gn, gm, b, gm, c]];


(* ::Text:: *)
(*Expansion of Dirac product unto basis of anti-symmetric gamma matrices with the techniques of [1905.00429, Sec 2.9.2].*)


RefineDiracProductsInternal[product:DiracProduct[_GammaM, _GammaM..]]:= Module[
		{aSym, indPos, pairCombs, prod= product},
	(* indPos is ordered by construction: its signature is +1 *)
	indPos= Flatten[Table[{m, n}, {m, Length@ prod}, {n, Length@ prod[[m]]}], 1];
	Sum[
		aSym= Complement[indPos, set];
		(*No metric with indices from the same GammaM: these indices are anti-symmetric*)
		pairCombs= DeleteCases[NonOverlappingPairs@ set,
			List@ OrderlessPatternSequence[{{a_, _}, {a_, _}}, ___]];
		Sum[
			Signature@ Join[Sequence@@ pairs, aSym]* Product[Metric@@ Extract[prod, p], {p, pairs}]*
			DiracProduct[GammaM@@ Extract[prod, aSym]]
		, {pairs, pairCombs}]
	, {k, 0, Length@ indPos, 2}, {set, Subsets[indPos, {k}]}]
];


(* ::Text:: *)
(*Dealing with transposed matrices*)


RefineDiracProducts::transp= "Some, but not all, matrices in `1` was transp. Could not match to basis";
RefineDiracProductsInternal@ x:DiracProduct[_Transp, ___]:= Module[{},
	If[!MatchQ[x, DiracProduct[_Transp..]],
		Message[RefineDiracProducts::transp, x];
		Return@ x;
	];

	(*Use RefineDircaProducts in case of a factor of -1 (e.g. from anti-commuting \[Gamma]5)*)
	RefineDiracProducts[Transp[x]]/. d_DiracProduct:> Transp@ d
];


(* ::Subsection:: *)
(*Collecting and manipulating the gamma matrices*)


(* ::Subsubsection::Closed:: *)
(*\[Gamma]CommuteQ*)


(* ::Text:: *)
(*Check if something commutes with gamma matrices: *)


\[Gamma]CommuteQ[Xop| DiracProduct| MomNCM]= False;
\[Gamma]CommuteQ@ Field[_, Fermion, ___]= False;
\[Gamma]CommuteQ@ TransposeThisSpinChain= False;
\[Gamma]CommuteQ@x_NCM:=Nor[LOpenSpinChainQ@x, ROpenSpinChainQ@x];
\[Gamma]CommuteQ[f_?\[Gamma]CommuteQ[x___]]:= AllTrue[{x}, \[Gamma]CommuteQ] (*And@@ \[Gamma]CommuteQ/@ {x}*);
\[Gamma]CommuteQ[f:Except[NCM][x___]]:=  False;
\[Gamma]CommuteQ@ Alternatives[Pattern, _Blank, _BlankSequence, _BlankNullSequence, _Except] = False;
\[Gamma]CommuteQ[_]:= True;


(* ::Subsubsection::Closed:: *)
(*Collect gamma matrices*)


(* ::Text:: *)
(*Collect all gamma matrices in Dirac products on the left of the non-commutative multiply:*)


CollectGammaMatrices@expr_:=expr/.{
		x_DiracProduct :> Hold@ x
	}//.{
		MomNCM[a___, x_?\[Gamma]CommuteQ, Hold@ g_DiracProduct, b___]:>
		MomNCM[a, Hold@g, x, b]
	}//ReleaseHold;


(* ::Section:: *)
(*Operations*)


(* ::Subsection:: *)
(*Fermion test functions*)


(* ::Subsubsection::Closed:: *)
(*Type of fermions*)


FermionL= Bar@ Field[_, Fermion, __] | Transp@ Field[_, Fermion, __] | EoM@ Bar@ Field[_, Fermion, __] | EoM@ Transp@ Field[_, Fermion, __];


FermionR= Field[_, Fermion, __] | Transp@ Bar@ Field[_, Fermion, __] | EoM@ Field[_, Fermion, __] | EoM@ Transp@ Bar@ Field[_, Fermion, __];


(* ::Subsubsection::Closed:: *)
(*FermionQ*)


FermionQ[FermionL | FermionR]:=True;
FermionQ[expr_]:=False;


(* ::Subsubsection::Closed:: *)
(*L/R OpenSpinChainQ*)


LOpenSpinChainQ@x_NCM:=LOpenSpinChainQ@ First@x
LOpenSpinChainQ[FermionL]:=False
LOpenSpinChainQ@x_:=!\[Gamma]CommuteQ@x


ROpenSpinChainQ@x_NCM:=ROpenSpinChainQ@ Last@x
ROpenSpinChainQ[FermionR]:=False
ROpenSpinChainQ@x_:=!\[Gamma]CommuteQ@x


(* ::Subsubsection::Closed:: *)
(*ClosedSpinChainQ*)


ClosedSpinChainQ@ x_NCM:= !FreeQ[x, Fermion
	] && Nor[
		LOpenSpinChainQ@ x,
		ROpenSpinChainQ@ x
	] && Last@ CumulativeSpin@ x === 0;
ClosedSpinChainQ@ _:= False;


(* ::Subsubsection::Closed:: *)
(*MajoranaQ *)


MajoranaQ@x_:=!FreeQ[x, Field[label_,Fermion,___]/;$FieldAssociation[label][SelfConjugate]]


(* ::Subsubsection::Closed:: *)
(*CumulativeSpin count*)


(* ::Text:: *)
(*Returns the cumulative count of fermion spins along an NCM product. Counts +1 for R-open fermions and -1 for L-open.*)


CumulativeSpin@ op:(_NCM|_MomNCM):= Module[{temp},
	temp= FermionCount/@ (List@@ op);
	FoldList[Plus, temp]
];


FermionCount@ FermionL= 1;
FermionCount@ FermionR= -1;
FermionCount@ x:(_NCM| _MomNCM| _Times):=
	Plus@@ FermionCount/@ (List@@ x);
FermionCount@ x:_Plus:= FermionCount@ First@ x;
FermionCount@_ = 0;


(* ::Subsection:: *)
(*Dirac structures*)


(* ::Subsubsection::Closed:: *)
(*GammaCC*)


DiracProduct[a___, GammaCC, GammaCC, b___]:= - DiracProduct[a, b];
DiracProduct[a___, s: _Proj | Gamma5, GammaCC, b___]:= DiracProduct[a, GammaCC, Transp@ s, b];
DiracProduct[a___, Transp@s: _Proj | Gamma5, GammaCC, b___]:= DiracProduct[a, GammaCC, s, b];
DiracProduct[a___, g_GammaM, GammaCC, b___]:= Power[-1, Ceiling[Length@ g/2]] DiracProduct[a, GammaCC, Transp@ g, b];
DiracProduct[a___, Transp@g_GammaM, GammaCC, b___]:= Power[-1, Ceiling[Length@ g/2]] DiracProduct[a, GammaCC, g, b];

(* For UFO purposes *)


(* ::Subsubsection::Closed:: *)
(*Transp*)


(* ::Text:: *)
(*Transp should work on closed NCM lines for convenience, although they are commutative*)


(* external fermions do not anti-commute *)
Transp@ expr_NCM:= (-1)^Floor[Count[expr, Field[Except[_Ext],Fermion,__], All]/2] Transp/@ Reverse@ expr;


(*Transp@ expr_NCM:= (-1)^Floor[(Count[expr, _?FermionQ]-Count[expr, _Ext])/2] Transp/@ Reverse@ expr;*)


Transp@ expr_?CommutativeQ:= expr/. ncm_NCM:> Transp@ ncm;


Transp[expr:(_Plus|_Times|_Operator)]:= Transp/@ expr;


Transp@ Transp@ x_:= x;
Transp@ EoM@ x_:= EoM@ Transp@ x;
Transp@ GammaCC= -GammaCC;
Transp@ Gamma5= Gamma5;
Transp@ p_Proj:= p;
Transp[]:= 1;
Transp@ x_DiracProduct:= Transp/@ x/. y_DiracProduct:> Reverse@ y; (*Transp and Reverse can cause multiplication with -1; hence, replacement rather than function call.*)


(* ::Text:: *)
(*Transp with multiple arguments (can occur in NCM rules for automatically transposing spinor lines).*)


Transp[x_, y__]:= Reverse[NCM[x, y]/. d_DiracProduct:> Transp@ d];


(* ::Subsubsection::Closed:: *)
(*CConj*)


(* ::Text:: *)
(*Properties of charge conjugation working on fermion fields *)


CConj::notfermion= "CConj is only supposed to act on a fermion field."


CConj@f:Field[label_,Fermion,_,_]/;$FieldAssociation[label][SelfConjugate]:= f;
CConj@f:Bar[Field[label_,Fermion,_,_]]/;$FieldAssociation[label][SelfConjugate]:= f;
CConj@f:Bar[Field[_,Fermion,_,_]]:= Transp@Bar@f \[CenterDot] DiracProduct@GammaCC;
CConj@f:Field[_,Fermion,_,_]:= DiracProduct@ GammaCC \[CenterDot] Transp@Bar@f;
CConj[DiracProduct[Proj[s_]]\[CenterDot] (f:Field[_,Fermion,_,_])]:= DiracProduct[Proj[-s]]\[CenterDot] CConj[f];
CConj[(f:Bar[Field[_,Fermion,_,_]])\[CenterDot] DiracProduct[Proj[s_]]]:= CConj[f]\[CenterDot] DiracProduct[Proj[-s]];
CConj@_:=(Message[CConj::notfermion]; Abort[]);


(* ::Subsubsection::Closed:: *)
(*Helper functions*)


ProjExpand[exp_]:=exp/. DiracProduct[b___,Proj[s_]]:> (DiracProduct[b] + s  DiracProduct[b,Gamma5])/2 //Expand;


(* ::Subsection:: *)
(*Fierzing*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


Fierz::notFierzable= "The Fierz function can only be used on the product of two closed spin chains with head NCM containing Dirac structures with head DiracProduct (or nothing for scalar currents). ";
Fierz::Order = "Fierz function only support Fierz Order {1,3,4,2} or {1,4,3,2}."


(* ::Subsubsection::Closed:: *)
(*Property*)


Options[Fierz]={Order->{1,4,3,2}, Evanescent-> True};
Options[Fierz4D]={Order->{1,4,3,2}};


(* ::Subsubsection::Closed:: *)
(*Helper functions*)


FieldLabel[f:_Field|Bar@_Field|Transp@_Field|Transp@Bar@_Field]:=First@Cases[{f}, Field[l_,_, _,_]:>l,Infinity];
FieldIndices[f:_Field|Bar@_Field|Transp@_Field|Transp@Bar@_Field]:=First@Cases[{f}, Field[_,_,ind_,_]:>ind,Infinity];


(* ::Subsubsection::Closed:: *)
(*4D Fierz function*)


Fierz4D[(SP1:NCM[field1_,G12___,field2_]),(SP2:NCM[field3_,G34___,field4_]), opt:OptionsPattern[]]:=
Fierz4D[SP1,SP2, opt]=
Module[
{Gbasisup, Gbasisdown, ind1, ind2, Gcoeff, i,j, C1=1, C2=1,
G1=G12*1, G2=G34*1, P1=1, P2=1, P3=1, P4=1, f1=field1, f2=field2, f3=field3, f4=field4, sgnO=1},

	ind1={Unique[],Unique[],Unique[],Unique[]};
	ind2={Unique[],Unique[],Unique[],Unique[]};

	Gbasisup[a_,b_,c_,d_]:={PL,PR,PR\[CenterDot] \[Gamma][a],PL\[CenterDot] \[Gamma][b],\[Sigma][c,d]/2}/.NCM[x_]:>x;
	Gbasisdown[a_,b_,c_,d_]:={PL,PR,PL\[CenterDot] \[Gamma][a],PR\[CenterDot] \[Gamma][b],\[Sigma][c,d]/2}/.NCM[x_]:>x;

	If[OptionValue@Order === {1,2,3,4}, Return[SP1*SP2]];
	If[!MemberQ[{{1,3,4,2},{1,4,3,2},{1,2,3,4}},OptionValue@Order], Message[Fierz::Order];Abort[]];
	If[OptionValue@Order==={1,3,4,2},
			f3=Transp@field4;
			f4=Transp@field3;
			sgnO=-1;
			If[G2=!=1,
				G2=Transp@G2;
				If[First@G2===-1, G2=-G2; sgnO=+1;]];
	];

	If[!ClosedSpinChainQ@SP1 || !ClosedSpinChainQ@SP2, Message[Fierz::notFierzable]; Abort[]];
	If[(Head[G1]=!= DiracProduct && G1=!=1) || (Head[G2]=!= DiracProduct  &&  G2=!=1), Message[Fierz::notFierzable]; Abort[]];

	If[Head[f1]===Transp, G1=(-CC)\[CenterDot] G1; f1=f1\[CenterDot] CC];
	If[Head[f2]===Transp, G1=G1\[CenterDot] (-CC); f2=CC\[CenterDot] f2];

	If[Head[f4]===Transp, G2=G2\[CenterDot] (-CC); f4=CC\[CenterDot] f4];
	If[Head[f3]===Transp, G2=(-CC)\[CenterDot] G2; f3=f3\[CenterDot] CC];

	If[!FreeQ[G1, _Proj],
		P1=DiracProduct@First@Cases[Transp@G1,_Proj,Infinity];
		P2=DiracProduct@First@Cases[G1,_Proj,Infinity];
	];

	If[!FreeQ[G2, _Proj],
		P3=DiracProduct@First@Cases[Transp@G2,_Proj,Infinity];
		P4=DiracProduct@First@Cases[G2,_Proj,Infinity];
	];

	G1 = G1 /._Proj->1;
	G2 = G2 /._Proj->1;

	If[!FreeQ[G1, x_GammaM/;(Length[x]===2)],
		sgnO*=-I;
		G1*=I;
	];
	If[!FreeQ[G2, x_GammaM/;(Length[x]===2)],
		sgnO*=-I;
		G2*=I;
	];

	Gcoeff=IdentityMatrix[5];
	For[j=1,j<6,++j,
		For[i=1,i<6,++i,
			Gcoeff[[i,j]]= -(1/4) DiracTrace[G1\[CenterDot] (Gbasisdown@@ind1)[[i]]\[CenterDot] G2\[CenterDot] (Gbasisdown@@ind2)[[j]],Dimensions->4] //Contract;
		]
	];

	sgnO*(Sum[
			Gcoeff[[i,j]]
				* (f1\[CenterDot] P1\[CenterDot] (Gbasisup@@ind2)[[j]]\[CenterDot] P4\[CenterDot] f4) * (f3\[CenterDot] P3\[CenterDot] (Gbasisup@@ind1)[[i]]\[CenterDot] P2\[CenterDot] f2)
		,{i,1,5},{j,1,5}] //Expand//Contract //RelabelIndices //LC2Gamma5) /.GammaM[a_,b_]/;(!OrderedQ[{a,b}]):> -GammaM[b,a]
]


(* ::Subsubsection::Closed:: *)
(*Fierz score*)


FbasisQ[SP1_,SP2_]:=FreeQ[OperatorToNormalForm@Operator[SP1,SP2]//.Join[{_Field->1, GammaCC->1,Proj[1]->1,Proj[-1]->1},Thread[NCM/@{DiracProduct@GammaM[x_],DiracProduct@Transp@GammaM[x_],DiracProduct@GammaM[x_, y_],
DiracProduct@Transp@GammaM[x_, y_]}->1]],NCM];


FierzScore[op_Operator]:=FierzScore[OperatorToNormalForm@op];


FierzScore[(SP1 : NCM[field1_, G12___, field2_] ),(SP2 : NCM[field3_, G34___, field4_]),rest___] :=
 Module[
  {result = 0, fields = {field1, field2, field3, field4}, order = OptionValue@Order,
   GaugeGroups, gind, pos, posmax, epsrep, G1=G12*1,G2=G34*1},
      (*Replace SU(2) epsilon by delta for counting*)
      epsrep=Cases[List@rest, CG[eps[_] | Bar[eps[_]] ,ind_/;(Length[ind]===2)]:>ind, Infinity]/.Bar[x_]:>x;
      fields=fields/.((#[[1]]->#[[2]])& /@epsrep);
  	(*Default rules: same fields, higher groups together*)
  	result = Plus[0,
    	    (*same labels favored*)
    		If[FieldLabel[fields[[1]]] === FieldLabel[fields[[2]]], +1 , -1],
    		If[FieldLabel[fields[[3]]] === FieldLabel[fields[[4]]], +1 , -1],
    		(*many transposed fields penalized*)
    		If[Xor[Head[fields[[1]]] === Transp, Head[fields[[2]]] === Transp] && Xor[Head[fields[[3]]] === Transp, Head[fields[[4]]] === Transp] , -6 , +.5],
    		(*tensor currents  penalized*)
    		If[! FreeQ[G1, GammaM[a___] /; (Length[List@a] > 1)] && ! FreeQ[G2, GammaM[b___] /; (Length[List@b] > 1)], -7, 0],
    		(*scalar currents  preferred over vector currents*)
    		If[FreeQ[G1, _GammaM] && FreeQ[G2, _GammaM], 0, -0.06],
    		(*Dirac structures not part of the basis penalized*)
    		If[FbasisQ[SP1, SP2], 0, -10],
    		(*more projectors penalized*)
    		Count[SP1*SP2,_Proj,All]*(-0.04),
    		(*contracted indices in bilinears encouraged*)
    		Total[If[MemberQ[FieldIndices@ fields[[2]], #], +0.01* GroupDimension@ $Groups@ GroupFromInd@ #, 0]&/@ FieldIndices@ fields[[1]]],
    		Total[If[MemberQ[FieldIndices@ fields[[4]], #], +0.01* GroupDimension@ $Groups@ GroupFromInd@ #, 0]&/@ FieldIndices@ fields[[3]]]
    		];

  	(*same group in only two fields favored, starting from biggest group*)
  	GaugeGroups = Reverse@SortBy[Keys@$GaugeGroups, GroupDimension@* $GaugeGroups[#, Group] &];
  	gind = Reverse[SortBy[Select[FieldIndices@#, (MemberQ[GaugeGroups, GroupFromRep@ #[[2]]] &)], GroupDimension@* $GaugeGroups[#[[2]], Group] & ]] & /@ fields;
      pos = Position[gind /. Index[_, g_] :> GroupFromRep@g, #] & /@ GaugeGroups /. {} -> Nothing;
    	If[pos =!= {} && pos =!= Nothing,
    		posmax = First@Transpose@First@pos;
    		If[Length[posmax]=!=2 && Length@pos>1,posmax=First@Transpose[pos[[2]]]];
     		result += If[MemberQ[{{1, 2}, {3, 4}}, posmax], +1 , 0];
     	];
    	result
    ]


FierzScore[expr_]:=(*FierzScore[expr]=*)(* cached identites can worsen performance significantly if expr is very complicated *)
Module[{NCMs,SP1,SP2, rest},
		NCMs=Cases[RemovePower@ expr, _NCM,All];
		If[Length@NCMs === 2,
			{SP1,SP2}=NCMs;
			rest=expr/._NCM->1;
			FierzScore[SP1,SP2,rest]
			,
			0
		]
	]


(* ::Subsubsection::Closed:: *)
(*Fierz function*)


Fierz[arg_, opt:OptionsPattern[]]:=
Module[
{rest, NCMs, SP1, SP2, result, groupStruct, fieldsinop, coeff, evaOperator, iniOperator, expr=RelabelIndices[arg,Unique->True]},
	NCMs=Cases[RemovePower@ expr, _NCM,All];
	If[Length@NCMs =!= 2,  Message[Fierz::notFierzable]; Abort[]];

	{SP1,SP2}=NCMs;

	rest=expr/._NCM->1;

	iniOperator=rest *SP1*SP2//RefineDiracProducts//RelabelIndices;

	If[Head@iniOperator===Plus, Fierz[#,opt]&/@iniOperator,
	result=Fierz4D[SP1,SP2,Order->OptionValue@Order];

	rest=Flatten@{rest/.Times->List};
	groupStruct=Times@@Cases[rest,_Delta|_CG];
	fieldsinop=Times@@Cases[rest,_Field|_FieldStrength|Bar[_Field]];
	rest=Times@@rest;

	evaOperator=If[OptionValue@Evanescent,
			DefineEvanescentOperator[groupStruct*fieldsinop*(SP1*SP2),groupStruct*fieldsinop*result,(*Exponent[rest*1,hbar]+1,*)Fierz]*rest/groupStruct/fieldsinop,
			0];
	evaOperator+result*rest//ReleaseHold//RelabelIndices//ContractCGs//Contract
	]
]


(* ::Subsection:: *)
(*Gamma Products Reduction*)


(* ::Subsubsection::Closed:: *)
(*4D GammaReduction function*)


(* ::Text:: *)
(*Projects two contracted spinor chains to the physical basis with gamma reduction. *)


Basis4D[]:=Module[{\[Mu],\[Nu]},
	{{PL,PL},{PR,PR},{PL,PR},{PR,PL},
	{\[Gamma]@\[Mu]\[CenterDot] PL,\[Gamma]@\[Mu]\[CenterDot] PL},{\[Gamma]@\[Mu]\[CenterDot] PR,\[Gamma]@\[Mu]\[CenterDot] PR},{\[Gamma]@\[Mu]\[CenterDot] PL,\[Gamma]@\[Mu]\[CenterDot] PR},{\[Gamma]@\[Mu]\[CenterDot] PR,\[Gamma]@\[Mu]\[CenterDot] PL},
	{\[Sigma][\[Mu],\[Nu]]\[CenterDot] PL,\[Sigma][\[Mu],\[Nu]]\[CenterDot] PL},{\[Sigma][\[Mu],\[Nu]]\[CenterDot] PR,\[Sigma][\[Mu],\[Nu]]\[CenterDot] PR}}
	];


InverseBasisTrace[]:= InverseBasisTrace[]= Simplify@ Inverse@ Table[
		bi[[1]]\[CenterDot] bj[[1]]\[CenterDot] bi[[2]]\[CenterDot] bj[[2]]// DiracTrace// ContractMetric
	, {bi, Basis4D[]}, {bj, Basis4D[]}];


GammaReduction4D[(SP1:TempNCM[field1_,G12_,field2_]),
	(SP2:TempNCM[field3_,G34_,field4_])]:= GammaReduction4D[SP1, SP2]=
Module[{Bcoeff,basis,f1=field1,f2=field2,f3=field3,f4=field4, G1= G12, G2= G34},

	If[Head[f1]===Transp, G1=(-CC)\[CenterDot] G1; f1=f1\[CenterDot] CC];
	If[Head[f2]===Transp, G1=G1\[CenterDot] (-CC); f2=CC\[CenterDot] f2];

	If[Head[f4]===Transp, G2=G2\[CenterDot] (-CC); f4=CC\[CenterDot] f4];
	If[Head[f3]===Transp, G2=(-CC)\[CenterDot] G2; f3=f3\[CenterDot] CC];

	(*Check if in the basis*)
	If[MatchQ[{G1, G2}/. Times[_Integer, ncm_NCM]-> ncm/.
		NCM-> Identity, $basis4DPattern],
		Return[SP1* SP2/. TempNCM-> NCM];
	];

	basis=Basis4D[];

	Bcoeff= InverseBasisTrace[] . Table[
		DiracTrace[bj[[1]]\[CenterDot] G1\[CenterDot] bj[[2]]\[CenterDot] G2]//ContractMetric
	,{bj,basis}]//Simplify;

	Sum[Bcoeff[[i]]*(f1\[CenterDot] basis[[i,1]]\[CenterDot] f2) * (f3\[CenterDot] basis[[i,2]]\[CenterDot] f4) , {i,1,Length@Bcoeff}]
]


(* ::Text:: *)
(*Patterns used to discern basis elements of the spinor chains*)


$basis4DPattern= Module[{\[Mu],\[Nu]},
	Alternatives[
		{PL, PL}, {PR, PR}, {PL, PR}, {PR, PL},
		{\[Gamma]@ \[Mu]\[CenterDot] PL, \[Gamma]@ \[Mu]\[CenterDot] PL}, {\[Gamma]@ \[Mu]\[CenterDot] PR, \[Gamma]@ \[Mu]\[CenterDot] PR},
		{\[Gamma]@ \[Mu]\[CenterDot] PL, \[Gamma]@ \[Mu]\[CenterDot] PR}, {\[Gamma]@ \[Mu]\[CenterDot] PR, \[Gamma]@ \[Mu]\[CenterDot] PL},
		{\[Gamma][\[Mu], \[Nu]]\[CenterDot] PL, \[Gamma][\[Mu], \[Nu]]\[CenterDot] PL}, {\[Gamma][\[Mu], \[Nu]]\[CenterDot] PR, \[Gamma][\[Mu], \[Nu]]\[CenterDot] PR}
	]/. {\[Mu]-> \[Mu]_, \[Nu]-> \[Nu]_}/. NCM-> Identity
];


(* ::Subsubsection::Closed:: *)
(*Internal GammaReduction function*)


(* ::Text:: *)
(*Applies gamma reduction identities to any operator with exactly two spinor-chains*)
(*Internal version with OperatorForm input -> Operator Form output*)


GammaReductionInternal[expr_, opt__]:= expr/. op_Operator:> GammaReductionInternal[op, opt];


GammaReductionInternal[op_Operator, produceEva_? BooleanQ, fourDim_? BooleanQ]:=
Module[{rest, NCMs, SP1, SP2, result, groupStruct, evaOperator, iniOperator= op},
	If[!Reducible4dOpQ@ iniOperator, Return@ op];
	
	NCMs= Cases[iniOperator, _NCM, All];
	If[Length@ NCMs =!= 2, Return@ op; ];
	{SP1, SP2}= NCMs;
	(*TempNCM is a dummy head. It prevents the automatic ordering of transposed fermion fields.*)
	NCMs= TempNCM@@@ NCMs;
	rest= iniOperator/._NCM-> 1// OperatorToNormalForm;
	
	Which[
	fourDim, (*purely four dimensional*)
		result= GammaReduction4D@@ NCMs;
		result= EpsExpand@ result;
	,FreeQ[NCMs, GammaCC], (*No charge conjugation, but d-dimensional*)
		result= GammaReduction4D@@ NCMs;
	,True, (*Average over transposed versions of spin lines with charge conjugation*)
		NCMs= NCMs/. ncm_TempNCM?(!FreeQ[#, GammaCC]&):> AverageSpinLine@ ncm;
		(*Expanding out sums and prefactors*)
		result= Threading@@ NCMs/. Threading-> GammaReduction4D
	];

	evaOperator= If[produceEva,
			groupStruct=Times@@Cases[{rest},_Delta|_CG,Infinity];
			DefineEvanescentOperator[groupStruct* SP1* SP2, groupStruct* result, GammaReduction]*
				rest/ groupStruct
		,
			0
		];
	evaOperator+ result* rest// ContractCGs// ContractDelta //Contract// NormalToOperatorForm
];


AverageSpinLine@ ncm_TempNCM:=Module[{out},
	out= ncm -Transp/@ Reverse@ ncm; (*sign from Grassmanian exchange*)
	out/ 2/. TempNCM[f1_, Times[d_DiracProduct, rest___], f2_]:> Times@ rest *TempNCM[f1, d, f2]//Expand
]


Threading[a___, p_Plus, b___]:= Threading[a, #, b]&/@ p;
Threading[a___, Times[x_, ncm_TempNCM], b___]:= x Threading[a, ncm, b];


(* ::Subsubsection::Closed:: *)
(*Public GammaReduction function*)


(* ::Text:: *)
(*Public version with *)
(*NormalForm input -> NormalForm output*)


Options[GammaReduction]= {FourDimensional-> False, Evanescent-> True};


GammaReduction[expr_, OptionsPattern[]](*? OptionsCheck*):= Block[{out},
	out= GammaReductionInternal[Operator@ ContractMetric@ RefineDiracProducts@ expr, OptionValue@ Evanescent, False];
	out= out// OperatorToNormalForm// ContractDelta// RelabelIndices;
	If[OptionValue@ FourDimensional,
		EpsExpand@ out
	,
		out
	]
];
