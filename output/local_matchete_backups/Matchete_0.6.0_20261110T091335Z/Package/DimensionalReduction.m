(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`DimensionalReduction`*)


(* ::Subtitle:: *)
(*Tools for computations of finite-temperature dimensional reduction*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["DRMatch"]
PackageExport["RemoveConstantTerms"]


PackageExport["DefineZeroModes"]


PackageExport["\[CapitalTau]"]
PackageExport["\[CapitalLambda]\[CapitalTau]"]
PackageExport["Thermalhbar"]


PackageExport["ZeroMode"]
PackageExport["ZeroModeNames"]


PackageExport["Spatial"]


PackageExport["SeriesThermalEFT"]
PackageExport["ThermalPowerCounting"]
PackageExport["SuperSoftMatch"]


PackageExport["\[CapitalLambda]bar2"]


PackageExport["ZeroModeMixing"]


PackageExport["EuclideanConventions"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["CollectMomentaSpatial"]
PackageScope["ExtractMomenta"]


PackageScope["WilsonToStatic"]


PackageScope["Temporal"]
PackageScope["TemporalProjector"]
PackageScope["TemporalExpand"]
PackageScope["TemporalWilson"]
PackageScope["TemporalField"]


PackageScope["IntroduceDebyeMasses"]
PackageScope["IntroduceEffectiveGaugeCouplings"]
PackageScope["RestoreDebyeMasses"]
PackageScope["RestoreEffectiveGaugeCouplings"]


PackageScope["LagrangianToStatic"]


PackageScope["CheckZeroModeMixing"]
PackageScope["DefinedZMQ"]


PackageScope["SeparateOutZeroModes"]
PackageScope["FreeOfZeroModes"]


PackageScope["LoopSumIntegrate"]
PackageScope["$FSumInt"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


DRMatch::usage = "DRMatch[Lag] returns the static EFT Lagrangian resulting from integrating out the massive Matsubara modes of Lag. Its options are EFTOrder."
RemoveConstantTerms::usage = "An option for DRMatch. If set to True, the final result removes the renormalization of the vacuum."
SuperSoftMatch::usage = "SuperSoftMatch[Lag] integrates out the temporal components of the Gauge bosons, assuming that the Debye masses are heavy in the SuperSoft limit."


\[CapitalLambda]\[CapitalTau]::usage = "Represents a factor 4\[Pi]\[CapitalTau], with T>0 being the temperature, appearing in finite-temperature computations."
\[CapitalTau]::usage = "Represents the temperature, with \[CapitalTau]>0, appearing in finite-temperature computations."
ZeroMode::usage = "Fields with this properties being True are the Matsubara zero mode of a field."


Spatial::usage = "Index type for the spacial components of a Lorentz index."


SeriesThermalEFT::usage = "SeriesThermalEFT[expr,EFTOrder->6] truncates expr in the order given by EFTOrder but takin into account the thermal power counting associated to the different couplings.  "
ThermalPowerCounting::usage = "ThermalPowerCounting is an option for various routine specifying the targeted order in terms of thermal EFT order. It also defines how does it scales with temperature the couplings."
Thermalhbar::usage = "\!\(\*SubscriptBox[\(\[HBar]\), \(T\)]\) is the power-counting parameter for the loop-order of dimensional reduction in the supersoft-scale."


\[CapitalLambda]bar2::usage = "Variable denoting the square of the matching scale in MSbar for dimensional reduction."


ZeroModeMixing::usage = "Option for the CheckLagrangian function. It takes the values True/False to determine whether the check is performed."
ZeroModeNames::usage = "Option for giving the names (and printing) of the zero modes of the different fields. It is given as an association <|\[Phi]->{\[Phi]0, \!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\)}|> or <|\[Phi]->\[Phi]0|>."


(* ::Subsubsection::Closed:: *)
(*Internal*)


DefineZeroModes::usage = "Set the Matsubara zero modes of the already defined fields."


Temporal::usage = "Index type for the temporal component of a Lorentz index."
TemporalProjector::usage = "TemporalProjector[\[Mu]] represents the four-dimensional vector \!\(\*SubscriptBox[\(\[DoubleStruckCapitalP]\), \(\[Mu]\)]\)=(1,0,0,0). In case of various entries, it represents the product of various vectors."


LoopSumIntegrate::usage = "Performs the scalar loop sum-integral over the loop momentum propagators."
$FSumInt::usage = "Represents the fermionic factor of difference between a bosonic and a fermionic sum-integral."


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Setting the Lagrangian to its static limit*)


(* ::Subsection:: *)
(*Zero Modes*)


(* ::Subsubsection::Closed:: *)
(*Useful definitions for spacial gauge group*)


NonSpGaugeGroups[]:= First@GaugeGroupByProperty[Field->#]&/@Select[FieldByProperty[ZeroMode->False],GaugeFieldQ]
SpGaugeGroups[]:= First@GaugeGroupByProperty[Field->#]&/@Select[FieldByProperty[ZeroMode->True],GaugeFieldQ]
ToSpGaugeGroup[]:= #->Symbol[ToString@#<>"Sp"]&/@ NonSpGaugeGroups[]
ToGaugeGroup[]:= Symbol[ToString@#<>"Sp"]->#&/@ NonSpGaugeGroups[]


NonSpGaugeCouplings[]:= $GaugeGroups[#,Coupling]&/@ NonSpGaugeGroups[]
ToSpGaugeCouplings[]:=#->Symbol[ToString@#<>"Sp"]&/@NonSpGaugeCouplings[]


(* ::Subsubsection::Closed:: *)
(*Define zero modes*)


Options@DefineZeroModes = {ZeroModeNames-> Automatic};


DefineZeroModes::definedlabel= "The label `1`  is already defined in other context for the zero mode of `2`. Please select another label with the option ZeroModeNames. ";


(* ::Text:: *)
(*Function that defines the corresponding zero modes of any already defined field and the scalar temporal component of every gauge group.*)


DefineZeroModes[OptionsPattern[]]:=Module[{gaugefields, scalars, vectors, properties, zeromodenames= OptionValue@ZeroModeNames,label, printing},

	(*Define the temporal and spatial components of the gauge bosons*)
	gaugefields= {#,$GaugeGroups[#][Field]}&/@NonSpGaugeGroups[];
	MapThread[
		Block[{},
		 If[!TrueQ[DefinedZMQ[#1]],
			DefineGaugeGroup[Symbol[ToString[#1]<>"Sp"],$GaugeGroups[#1,Group],Symbol[ToString[$GaugeGroups[#1,Coupling]]<>"Sp"],Symbol[ToString[$GaugeGroups[#1,Field]]<>"Sp"],
				ZeroMode->True,
				NiceForm->{ToString[OverHat[ToString@$GaugeGroups[#1,Coupling]],StandardForm],ToString[OverVector[ToString@$GaugeGroups[#1,Field]],StandardForm]},
				FundAlphabet-> If[MemberQ[Keys@$IndexAlphabets,#1[fund]], ToString/@DeleteCases[Values@$IndexAlphabets[#1[fund]],_SubscriptBox], None],
				AdjAlphabet-> If[MemberQ[Keys@$IndexAlphabets,#1[adj]], ToString/@DeleteCases[Values@$IndexAlphabets[#1[adj]],_SubscriptBox], None]
			];
			
			If[zeromodenames===Automatic || FreeQ[Keys@zeromodenames, #2[[2]]], 
				(*If no given name to the zeromode of B, associate it B0*)
				label= Symbol[ToString[#2[[2]]]<>"0"];
				printing= ToString[Subscript[ToString[#2[[2]]],"0"],StandardForm];
			,
				If[Head[zeromodenames[#2[[2]]]]===List, 
					(*Name given and also a printing*)
					label= zeromodenames[#2[[2]]][[1]]; printing= ToString[zeromodenames[#2[[2]]][[2]],StandardForm];
					,
					(*Name given but no printing*)
					label= zeromodenames[#2[[2]]]; printing=  ToString[label,StandardForm];
				];
			];
			
			(*If the label is already defined, ask the user to give another*)
			If[Defined[label], Message[DefineZeroModes::definedlabel, label, #1]; Abort[];];
			
			(*Define the temporal gauge boson*)
			DefineField[label,Scalar,SelfConjugate->True,Indices->If[$GaugeGroups[#1,Group]===U1,{},Symbol[ToString@#1<>"Sp"][adj]],
				Mass->0,
				ZeroMode->True,
				NiceForm->printing
			];
			DefinedZMQ[#1]= True;
			ZM[#1]= label;
		 ];
		]&

	,
		{NonSpGaugeGroups[],gaugefields}
	];
	
	(*Define zero modes for non-gauge vectors*)
	vectors=KeyDrop[$FieldAssociation,Select[Keys@$FieldAssociation,GaugeFieldQ]];
	vectors=Select[vectors,(#[Type]===Vector&&#[ZeroMode]===False)&]//Normal;
	properties=DeleteCases[Normal@Values@vectors, Rule[Type,_]|Rule[Heavy,_]|Rule[Mass,_]|Rule[ZeroMode,_],All];
	MapThread[
		Block[{},
		 If[!TrueQ[DefinedZMQ[#1]],
		 
			If[zeromodenames===Automatic || FreeQ[Keys@zeromodenames, #2[[2]]], 
				(*If no given name to the zeromode of B, associate it B0*)
				label= Symbol[ToString[#2[[2]]]<>"0"];
				printing= ToString[Subscript[ToString[#2[[2]]],"0"],StandardForm];
			,
				If[Head[zeromodenames[#2[[2]]]]===List, 
					(*Name given and also a printing*)
					label= zeromodenames[#2[[2]]][[1]]; printing= ToString[zeromodenames[#2[[2]]][[2]],StandardForm];
					,
					(*Name given but no printing*)
					label= zeromodenames[#2[[2]]]; printing=  ToString[label,StandardForm];
				];
			];
			
			(*If the label is already defined, ask the user to give another*)
			If[Defined[label], Message[DefineZeroModes::definedlabel, label, #1]; Abort[];];
			
			(DefineField[label,Vector,##2/.ToSpGaugeGroup[],Mass->Light,ZeroMode->True,NiceForm->printing]);
			DefinedZMQ[#1]=True;
			ZM[#1]= label;
		 ];	
		]&
	,
		{Keys@vectors,properties}
	];
	
	(*Define zero modes for scalars*)
	scalars=Select[$FieldAssociation,(#[Type]===Scalar&&#[ZeroMode]===False)&]//Normal;
	properties=DeleteCases[Normal@Values@scalars, Rule[Type,_]|Rule[Heavy,_]|Rule[Mass,_]|Rule[ZeroMode,_],All];
	MapThread[
		Block[{},
		  If[!TrueQ[DefinedZMQ[#1]],
			If[zeromodenames===Automatic || FreeQ[Keys@zeromodenames, #1], 
				(*If no given name to the zeromode of \[Phi], associate it \[Phi]0*)
				label= Symbol[ToString[#1]<>"0"];
				printing= ToString[Subscript[ToString[#1],"0"],StandardForm];
			,
				If[Head[zeromodenames[#1]]===List, 
					(*Name given and also a printing*)
					label= zeromodenames[#1][[1]]; printing= ToString[zeromodenames[#1][[2]],StandardForm];
					,
					(*Name given but no printing*)
					label= zeromodenames[#1]; printing=  ToString[label,StandardForm];
				];
			];
			
			(*If the label is already defined, ask the user to give another*)
			If[Defined[label], Message[DefineZeroModes::definedlabel, label, #1]; Abort[];];
			
			(*Define the zero mode field*)
			(DefineField[label,Scalar,##2/.ToSpGaugeGroup[],Mass->0,ZeroMode->True,NiceForm-> printing]);
			
			DefinedZMQ[#1]=True;
			ZM[#1]= label;
		 ];	
		]&
	,
		{Keys@scalars,properties}
	];
]


(* ::Subsubsection::Closed:: *)
(*Zero mode mixing check*)


(* ::Text:: *)
(*A term in the Lagrangian cannot have a mixing of zero modes and four dimensional fields*)


CheckLagrangian::zeromodemixing= "The following term(s) in the expression is a mixing of static and dynamic fields: `1`";


ZeroModeMixingQ@ term_:= Module[{zeromodesQ, nonzeromodesQ},
	zeromodesQ= AnyTrue[OccurringFields[term],($FieldAssociation[#][ZeroMode])&];
	nonzeromodesQ= AnyTrue[OccurringFields[term],(Not[$FieldAssociation[#][ZeroMode]])&];
	And[zeromodesQ,nonzeromodesQ]
];


CheckZeroModeMixing[Lag_, detoutput_?BooleanQ]:=Module[{ZMmixingterms=0, mZeroModeMixing},
	Do[
		If[ZeroModeMixingQ@term,
		ZMmixingterms+=term;
		];
	,
	{term, TermsToList@Lag}];
	
	If[(mZeroModeMixing= (ZMmixingterms=!=0)), 
		Message[CheckLagrangian::zeromodemixing,Format[ZMmixingterms,NiceForm]];
	];
	{mZeroModeMixing, ZMmixingterms}
]


(* ::Subsubsection::Closed:: *)
(*FreeLag for zero modes*)


FreeLag[ZeroMode]:=Sum[FreeLag@ f, {f, FieldByProperty[BackgroundField-> False,ZeroMode->True]}];


(* ::Subsubsection::Closed:: *)
(*Derivatives on zero modes*)


Field[\[Phi]_?($FieldAssociation[#,ZeroMode]===True&),type_,ind_,{x___,Index[\[Mu]_,Lorentz],y___}]:=Field[\[Phi],type,ind,{x,Index[\[Mu],Spatial],y}]
FieldStrength[A_?($FieldAssociation[#,ZeroMode]===True&),lor_,ind_,{x___,Index[\[Mu]_,Lorentz],y___}]:=FieldStrength[A,lor,ind,{x,Index[\[Mu],Spatial],y}]


(* ::Subsubsection::Closed:: *)
(*FieldStrenght and vector field of zero modes*)


FieldStrength[A_?($FieldAssociation[#,ZeroMode]===True&),{Index[mu_,Lorentz],Index[nu_,Lorentz]},ind_,der_]:= FieldStrength[A,{Index[mu,Spatial],Index[nu,Spatial]},ind,der]
Field[A_?($FieldAssociation[#,ZeroMode]===True&),Vector@Index[mu_,Lorentz],ind_,der_]:= Field[A,Vector@Index[mu,Spatial],ind,der]


(* ::Subsubsection::Closed:: *)
(*Free of Zero Modes*)


FreeOfZeroModes[expr_]:= And@@(!$FieldAssociation[#][ZeroMode]&/@OccurringFields[expr] && FreeQ[expr,Spatial])


(* ::Subsection:: *)
(*Static limit*)


(* ::Subsubsection::Closed:: *)
(*Lagrangian to static*)


LagrangianToStatic[expr_]:=Module[{ZeroModePart, out=expr},
	out=RemovePower@ out;
	
	(*Lorentz derivatives*)
	out=out//.Field[field_,inds__,{x___,Index[\[Mu]_,Lorentz],y___}]:>(Field[field,inds,{x,Index[\[Mu],Spatial],y}] + Field[field,inds,{x,Index[\[Mu],Temporal],y}]);
	out=out//.FieldStrength[field_,inds__,{x___,Index[\[Mu]_,Lorentz],y___}]:>(FieldStrength[field,inds,{x,Index[\[Mu],Spatial],y}] + FieldStrength[field,inds,{x,Index[\[Mu],Temporal],y}]);
	
	(*Gauge indices in CGs*)
	out=out/. cg_CG :> (cg/.ToSpGaugeGroup[]);
	
	(* Gauge coupling *)
	out=out/. ToSpGaugeCouplings[];
	
	(*Decomposition of field-strenght tensors*)
	out=out/.FieldStrength[field_,{Index[\[Mu]_,\[Mu]Type_],Index[\[Nu]_,\[Nu]Type_]},inds_,devs_]:> Module[{Spinds=inds/.ToSpGaugeGroup[], Spfield=Symbol[ToString[field]<>"Sp"]},
		Switch[{\[Mu]Type,\[Nu]Type},
			{Lorentz,Lorentz},
				FieldStrength[Spfield,{Index[\[Mu],Spatial],Index[\[Nu],Spatial]},Spinds,devs] + FieldStrength[Spfield,{Index[\[Mu],Temporal],Index[\[Nu],Spatial]},Spinds,devs] + FieldStrength[Spfield,{Index[\[Mu],Spatial],Index[\[Nu],Temporal]},Spinds,devs]
			,
			{Lorentz,Spatial|Temporal},
				FieldStrength[Spfield,{Index[\[Mu],Spatial],Index[\[Nu],\[Nu]Type]},Spinds,devs] + FieldStrength[Spfield,{Index[\[Mu],Temporal],Index[\[Nu],\[Nu]Type]},Spinds,devs]
			,
			{Spatial|Temporal,Lorentz},
				FieldStrength[Spfield,{Index[\[Mu],\[Mu]Type],Index[\[Nu],Spatial]},Spinds,devs] + FieldStrength[Spfield,{Index[\[Mu],\[Mu]Type],Index[\[Nu],Temporal]},Spinds,devs]
			,
			{Spatial|Temporal,Spatial|Temporal},
				FieldStrength[Spfield,{Index[\[Mu],\[Mu]Type],Index[\[Nu],\[Nu]Type]},Spinds,devs]
		]
	];
		
	out=out//ReleaseHold//CanonizeSpatial;
	
	(*Non Abelian FS*)
	out=out/.{fs:FieldStrength[field_,{\[Mu]:Index[_,\[Mu]Type_],\[Nu]:Index[_,\[Nu]Type_]},inds_/;(Length@inds>0),devs_]:>Switch[{\[Mu]Type,\[Nu]Type},
		{Spatial,Temporal},
			TemporalProjector[\[Nu]]CD[Append[devs,\[Mu]],TemporalComponent[GroupFromInd[First@inds],First@inds]]
		,
		{Temporal,Spatial},
			-TemporalProjector[\[Mu]]CD[Append[devs,\[Nu]],TemporalComponent[GroupFromInd[First@inds],First@inds]]
		,
		{Spatial,Spatial},
			fs
		,
		{Temporal,Temporal},
			0
		]};
		
	(*Abelian FS*)
	out=out/.Table[With[{tcomponent=TemporalComponent[group,{}]},
		fs:FieldStrength[$GaugeGroups[group,Field],{\[Mu]:Index[_,\[Mu]Type_],\[Nu]:Index[_,\[Nu]Type_]},{},devs_]:>Switch[{\[Mu]Type,\[Nu]Type},
			{Spatial,Temporal},
				TemporalProjector[\[Nu]]CD[Append[devs,\[Mu]],tcomponent]
			,
			{Temporal,Spatial},
				-TemporalProjector[\[Mu]]CD[Append[devs,\[Nu]],tcomponent]
			,
			{Spatial,Spatial},
				fs
			,
			{Temporal,Temporal},
				0	
		]]
	,{group,GaugeGroupByProperty[Abelian->True]}];
		
	out= out//CanonizeSpatial//TemporalExpand;
	out= out/.(Field[#,Scalar,inds__]:>Field[ZM[#],Scalar,inds]&/@FieldByProperty[ZeroMode->False])/.ToSpGaugeGroup[];	
	out//BetterContractCGs//RelabelIndices
]


(* ::Subsubsection::Closed:: *)
(*Wilson to static*)


WilsonToStatic[expr_Plus]:= WilsonToStatic/@expr;


WilsonToStatic[expr_]:= Block[{out},
	out= expr/.WilsonTerm-> WilsonSplit;
	out= out //.{WilsonSplit[field_,inds_,{x___,Index[\[Mu]_,Lorentz],y___}]:>
		WilsonSplit[field,inds,{x,Index[\[Mu],Spatial],y}] + WilsonSplit[field,inds,{x,Index[\[Mu],Temporal],y}]};
	out= CanonizeSpatial[out]/.WilsonSplit-> WilsonTerm
]


(* ::Subsubsection::Closed:: *)
(*Temporal Component*)


(* ::Text:: *)
(*Function to obtain the scalar temporal component of a Gauge group.*)


TemporalComponent[gaugegroup_]:= ZM[gaugegroup/.ToGaugeGroup[]]
TemporalComponent[gaugegroup_,index_]:= ZM[gaugegroup/.ToGaugeGroup[]][index]*$GaugeGroups[gaugegroup,Coupling][] /. ToSpGaugeCouplings[]
TemporalComponent[gaugegroup_,{}]:= ZM[gaugegroup/.ToGaugeGroup[]][]*$GaugeGroups[gaugegroup,Coupling][]/. ToSpGaugeCouplings[]


(* ::Subsubsection::Closed:: *)
(*Temporal Projector*)


TemporalProjector[x_Symbol]:=TemporalProjector[Index[x,Temporal]];
TemporalProjector[Index[x_,Lorentz]]:=TemporalProjector[Index[x,Temporal]];
TemporalProjector[Index[x_,Spatial]]:=0;


Bar@x_TemporalProjector:=x


TemporalProjector/:TemporalProjector[x_]^2=1


(* ::Subsection:: *)
(*Temporal Expand*)


(* ::Subsubsection::Closed:: *)
(*Temporal Wilson*)


(* ::Text:: *)
(*Temporal derivatives acting on a Wilson Term*)


WilsonTerm[field_,{i_,j_},{Index[\[Mu]_,Temporal],rest___}]:=Module[{x},TemporalProjector[\[Mu]] TemporalWilson[field,{i,x},{}]*WilsonTerm[field,{x,j},{rest}]];


WilsonTerm[field_,{i_,j_},{inds__/;FreeQ[{inds},Temporal|Lorentz],Index[\[Mu]_,Temporal],rest___}]:=Module[
	{x,subsets,indsUnique}
,
	indsUnique=Map[Index[Unique[],#[[2]]]&,{inds}];
	subsets=Subsets@indsUnique;
	Sum[TemporalProjector[\[Mu]]TemporalWilson[field,{i,x},Complement[indsUnique,s]]*WilsonTerm[field,{x,j},Join[s,{rest}]],{s,subsets}]/.Thread[indsUnique->{inds}]
]


(* ::Text:: *)
(*Giving gauge indices to TemporalWilson*)


TemporalWilson[field_,{i_Symbol,j_Symbol},devs_]:=Module[
		{conj,fieldLabel,indices,lorentz, 
		gaugeIndices, flavorIndices, flavorDeltas},

	fieldLabel= If[(conj= MatchQ[field, _Conj]), First@ field, field];
	lorentz= If[$FieldAssociation[fieldLabel, Type] === Vector, Metric[i, j], 1];
	indices= If[conj, Bar, Identity]@ Map[{Index[i, #1], Bar@ Index[j, #1]} &,
		$FieldAssociation[fieldLabel, Indices]];
		
	gaugeIndices= Cases[indices, _? (MemberQ[Keys@ $GaugeGroups, GroupFromInd@ First@ #] &)];
	flavorIndices= Complement[indices, gaugeIndices];
	flavorDeltas= Times@@ Delta@@@ flavorIndices;
		
	lorentz* flavorDeltas*TemporalWilson[field,gaugeIndices,devs]
]


(* ::Text:: *)
(*For singlets TemporalWilson is vanishing*)


TemporalWilson[\[Phi]_/;HoldPattern@GaugeSingletQ[\[Phi]/.Conj->Identity],rest__]:=0;


(* ::Text:: *)
(*TemporalWilson objects accept derivatives on them*)


Matchete`CovariantDerivative`PackagePrivate`DZeroQ[TemporalWilson]=False;
TakeDev[newDevs_, HoldPattern@TemporalWilson[field_,{i_Symbol,j_Symbol},devs_]]:= TemporalWilson[field,{i,j},Join[newDevs, devs]];


(* ::Subsubsection::Closed:: *)
(*Temporal Field*)


(* ::Text:: *)
(*Properties of TemporalField*)


NonBarableQ@TemporalField=False;
Bar@TemporalField[field_,rest__]:=TemporalField[Conj@field,rest];


(* ::Text:: *)
(*Temporal derivatives acting on a Field (it already has gauge indices)*)


Field[field_,type_,inds_,{Index[\[Mu]_,Temporal],rest___}]:=Module[{x},TemporalProjector[\[Mu]] TemporalField[field,{inds,inds/.Index[_,rep_]->Index[x,rep]},{}]Field[field,type,inds/.Index[_,rep_]->Index[x,rep],{rest}]]


Field[field_,type_,inds_,{devs__/;FreeQ[{devs},Temporal|Lorentz],Index[\[Mu]_,Temporal],rest___}]:=Module[
	{x,subsets,devsUnique}
,
	devsUnique=Map[Index[Unique[],#[[2]]]&,{devs}];
	subsets=Subsets@devsUnique;
	Sum[TemporalProjector[\[Mu]] TemporalField[field,{inds,inds/.Index[_,rep_]->Index[x,rep]},Complement[devsUnique,s]]Field[field,type,inds/.Index[_,rep_]->Index[x,rep],Join[s,{rest}]],
	{s,subsets}]/.Thread[devsUnique->{devs}]
]


(* ::Text:: *)
(*For singlets TemporalField is vanishing*)


TemporalField[\[Phi]_/;HoldPattern@GaugeSingletQ[\[Phi]/.Conj->Identity],rest__]:=0;


(* ::Text:: *)
(*TemporalField objects accept derivatives on them*)


Matchete`CovariantDerivative`PackagePrivate`DZeroQ[TemporalField]=False;
TakeDev[newDevs_, HoldPattern@TemporalField[field_,{inds1_,inds2_},devs_]]:= TemporalField[field,{inds1,inds2},Join[newDevs, devs]];


(* ::Subsubsection::Closed:: *)
(*Temporal Field Strength*)


(* ::Text:: *)
(*Temporal  derivatives  acting  on  a  FieldStrength  (it  already  has  gauge  indices)*)


FieldStrength[field_,lorinds_,inds_,{Index[\[Mu]_,Temporal],rest___}]:= Module[{x},TemporalProjector[\[Mu]]TemporalField[field,{inds,inds/.Index[_,rep_]->Index[x,rep]},{}]FieldStrength[field,lorinds,inds/.Index[_,rep_]->Index[x,rep],{rest}]]


FieldStrength[field_,lorinds_,inds_,{devs__/;FreeQ[{devs},Temporal|Lorentz],Index[\[Mu]_,Temporal],rest___}]:=Module[
	{x,subsets,devsUnique}
,
	devsUnique=Map[Index[Unique[],#[[2]]]&,{devs}];
	subsets=Subsets@devsUnique;
	Sum[TemporalProjector[\[Mu]] TemporalField[field,{inds,inds/.Index[_,rep_]->Index[x,rep]},Complement[devsUnique,s]]FieldStrength[field,lorinds,inds/.Index[_,rep_]->Index[x,rep],Join[s,{rest}]],
	{s,subsets}]/.Thread[devsUnique->{devs}]
]


(* ::Subsubsection::Closed:: *)
(*Temporal Expand*)


(* ::Text:: *)
(*Gives the Temporal Wilson in terms of temporal gauge bosons*)


TemporalWilsonExpand[field_,gaugeIndices_,devInds_]:= Module[{conj, fieldCharges, fieldLabel, temporalfields},
	fieldLabel= If[(conj= MatchQ[field, _Conj]), First@ field, field];

	fieldCharges= If[conj, MapAt[Minus, #, {All, 1}]&, Identity]@ $FieldAssociation[fieldLabel, Charges];

	temporalfields=If[conj,I,-I]*Sum[Module[{b,A},
			CD[devInds,TemporalComponent[GroupFromInd[First@gind],A]]*If[conj, Bar, Identity]@FieldGenerators[field/.Conj->Identity,GroupFromInd@First@gind,{A,##}&@@First/@(gind/.Bar->Identity)]Times@@(Delta@@@Complement[gaugeIndices,{gind}])],{gind,gaugeIndices}]
			-I*Times@@(Delta@@@gaugeIndices)Sum[First@chr*CD[devInds,TemporalComponent[Head@chr,{}]],{chr,fieldCharges}];

	(* The gauge indices in the CGs should be in the Spatial gauge *)
	temporalfields/. cg_CG :> (cg/. ToSpGaugeGroup[])
]


(* ::Text:: *)
(*Gives the Temporal Wilson in terms of temporal gauge bosons*)


TemporalFieldExpand[field_,{inds1_,inds2_},devInds_]:= Module[{conj, fieldCharges, indices, fieldLabel, flavorIndices, gaugeIndices, flavorDeltas, temporalfields},
	fieldLabel= If[(conj= MatchQ[field, _Conj]), First@ field, field];
	indices= If[conj, Bar, Identity]/@ Transpose[{inds1,inds2}];	

	(*Determines gauge and flavor indices*)
	gaugeIndices= Cases[indices, _? (MemberQ[Keys@ $GaugeGroups, GroupFromInd@ First@ #] &)];
	flavorIndices= Complement[indices, gaugeIndices];
	fieldCharges= If[conj, MapAt[Minus, #, {All, 1}]&, Identity]@ $FieldAssociation[fieldLabel, Charges];
	flavorDeltas= Times@@ Delta@@@ flavorIndices;

	temporalfields=If[conj,I,-I]*
		Sum[Module[{b,A},CD[devInds,TemporalComponent[GroupFromInd[First@gind],A]]* If[conj, Bar, Identity]@FieldGenerators[field/.Conj->Identity,GroupFromInd@First@gind,{A,First@(gind[[1]]/.Bar->Identity),First@(gind[[2]]/.Bar->Identity)}]Times@@(Delta@@@Complement[gaugeIndices,{gind}])],{gind,gaugeIndices}]
		-I*Times@@(Delta@@@gaugeIndices)* Sum[First@chr*CD[devInds,TemporalComponent[Head@chr,{}]],{chr,fieldCharges}];

	(* The gauge indices in the CGs should be in the Spatial gauge *)
	flavorDeltas* temporalfields/. cg_CG :> (cg/. ToSpGaugeGroup[])
]


(* ::Text:: *)
(*Main Function*)


TemporalExpand[expr_]:=expr//. TemporalWilson->TemporalWilsonExpand //. TemporalField->TemporalFieldExpand


(* ::Section:: *)
(*Dimensional Reduction*)


(* ::Subsection:: *)
(*Dimensional Reduction*)


(* ::Text:: *)
(*Main function that performs the dimensional reduction of a Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Options*)


Options@ DRMatch = { LoopOrder->1, EFTOrder -> 6 , Verbose -> True, RemoveConstantTerms -> True, ZeroModeNames -> Automatic };


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DRMatch::orderArgument= "The option EFTOrder-> `1`  is neither an integer nor Automatic and cannot be used as an order specification.";


DRMatch::heavyvec= "Lagrangian contains massive vector(s) `1`. Massive vectors are not yet supported for loop-level matching. Stay tuned for future updates\[Ellipsis]";


(* ::Subsubsection::Closed:: *)
(*Convert Heavy Masses to Light*)


(* ::Text:: *)
(*Function to treat temporarly heavy masses as light ones*)


ConvertHeavyMassestoLight[lag_]:= Module[{heavyfields,masses,changeofmasses,  out},
	heavyfields=Select[OccurringFields@lag,$FieldAssociation[#,Heavy]&];
	If[heavyfields==={},out=lag,
	Print["WARNING: Matchete does not currently support heavy fields in dimensional reduction. The mass of the fields "<> ToString[heavyfields]<>" has been changed to light."];
	masses= {$FieldAssociation[#,Mass],$FieldAssociation[#,Type]}&/@heavyfields;
	changeofmasses= (Coupling[#[[1]],{},0]->Coupling[#[[1]],{},If[#[[2]]==Fermion,1,2]])&/@masses;
	out= lag/.changeofmasses;
	];	
	out
]


(* ::Subsubsection::Closed:: *)
(*Main function*)


DRMatch[expr_,OptionsPattern[]]:=Module[{loopOrder=OptionValue@LoopOrder,order,
													lag=expr, ZeroModePart, TreeLevelLag, OneLoopLag,
													powerTraces,logTraces,traceResults, ptr,
													heavyVectors, LDR,field,i, time},
	(* Suggest the relevant literature and change logo to thermal *)
	AddToBibliography["Matchotter","Using dimensional reduction"];
	ChangeMatcheteLogo["Matchotter"];
	
	time=First@Timing[
		(* Determine order argument and cut the Lagrangian to the given order *)
		{order, lag}= TruncateLagrangian[OptionValue@ EFTOrder, lag];
		
		(*Heavy masses are treated as light*)
		lag= ConvertHeavyMassestoLight[lag];
		
		(* Substract off terms with zero modes *)
		{ZeroModePart,lag}= SeparateOutZeroModes@lag;
		ZeroModePart=Which[
						loopOrder===0, ZeroModePart/. hbar->0
					,
						loopOrder==={1}, hbar Coefficient[ZeroModePart, hbar]
					,
						loopOrder===1, ZeroModePart
					];		
	
		(* Define the zero modes of the present fields *)
		DefineZeroModes[ZeroModeNames -> OptionValue@ZeroModeNames];
	
		(* At tree level, we just need to remove fermions and take the static limit (see later) *)
		TreeLevelLag = If[loopOrder==={1}, 0, lag/.Field[_,Fermion,__]->0];
	
		(* One-loop piece *)
		If[MatchQ[loopOrder, 1|{1}],
			(*Check that no heavy vectors appear in loop-level matching*)
			If[MatchQ[loopOrder, 1|{1}],
				heavyVectors= MassiveVectorsInLag@ lag;
				If[Length@ heavyVectors =!= 0,
					Message[DRMatch::heavyvec, heavyVectors];
					Abort[];
				];
			];
		
			(* set the expanded Lagrangian and determine X terms *)
			SetCurrentLagrangian[lag, 1, order, Mode -> DimReduction];
	
			(* determine the Log traces *)
			$MonitorString="";
			logTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
						Sum[
							If[Or[
								Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field, Charges-> {__}] > 0,
								Length@ FieldByProperty[Sequence@@ Normal@ $FieldTypes@ field,
								Indices-> inds_/; IntersectingQ[GroupFromRep/@ inds, Keys@ $GaugeGroups]] > 0
							],
								$MonitorString = StringForm["Evaluating log-type supertrace: `1`", field/. fieldFormat];
									LogTypeSTr[field, order,  Mode -> DimReduction]
							,
								0
							]
						,{field,{hScalar,hFermion,hVector,hGhost,hAntiGhost,lScalar,lFermion,lVector,lGhost,lAntiGhost}}
						],
						$MonitorString					
					];
		
			(* list the possible power traces, including light fields *)
			ptr = ListPowerTypeTraces[ order, True ];
		
			powerTraces = OptionalMonitor[TrueQ@ OptionValue@ Verbose,
								Table[
									$MonitorString = StringForm["Evaluating power-type supertrace: `1` \t (`2` / `3`)", ptr[[i]] /. fieldFormat, i, Length@ ptr];
									PowerTypeSTr[ptr[[i]],order,Mode->DimReduction],			
								{i, Length@ptr}]
							,
								$MonitorString];
			traceResults= Riffle[ptr,powerTraces];
			powerTraces=Total@powerTraces;
		
			(* This is the ine-Loop dimensionally reduced Lagrangian *)
			OneLoopLag =  I (logTraces + powerTraces);
		,
			OneLoopLag=0	
		];
	];
	
	(* Collect tree and loop pieces, take the static limit and apply Euclidean transformation *)
	OptionalMonitor[TrueQ@ OptionValue@ Verbose,
		LDR = TreeLevelLag + OneLoopLag// LagrangianToStatic//EuclideanConventions;
	,
		StringForm["Taking the static limit..."]
	];
	
	(*Collect all the pieces and simplify*)
	LDR = LDR + ZeroModePart;
			
	(* Substract off constants unless specified otherwise *)
	LDR = If[OptionValue@RemoveConstantTerms===True, LDR-(LDR/.{_Field->0, _FieldStrength->0}), LDR];

		
	(* This is only performed if the result of the matching should be saved for validation *)
	If[$ValidationRun,
		QuietEcho@ SaveValidationResults[traceResults, LDR, time, lag, "Directory"->"DRResults"];
	];
	
	LDR
]


(* ::Subsubsection::Closed:: *)
(*Separate out zero modes *)


(* ::Text:: *)
(*Splits zero-mode terms off from a Lagrangian. *)
(*	Return: {zero-mode terms, remaining expression}*)


SeparateOutZeroModes@ expr_:= Module[{const, terms, ZeroModeTerms, rest},
	(*substract off constants*)
	{const, terms}= SeparateOutConstants@ expr;
	
	(*zero mode terms do not mix with non-zero modes*)
	ZeroModeTerms= terms/. (#->0&/@FieldByProperty[ZeroMode->False]);
	rest=terms/. (#->0&/@FieldByProperty[ZeroMode->True]);
	
	(*add back the constants to the non-zero-mode terms*)
	{ZeroModeTerms, rest+ const}
]


(* ::Subsection:: *)
(*Change of metric signature *)


(* ::Subsubsection::Closed:: *)
(*Euclidean conventions*)


EuclideanConventions[expr_]:=Module[{out,GaugeGroups, GaugeCouplingReplacement, $C2Euc},
	(* The scalar from the 0th component of a gauge field picks a factor of I *)
	out = expr/.(x:Field[#,a__]:>I*Field[#,a]&/@(TemporalComponent/@NonSpGaugeGroups[]));

	(* Every derivative adds a factor of I (FieldStrength counts as two derivatives) *)
	out= out /.{f:Field[__,devs_]:> Power[$C2Euc,Length@devs] f, fs:FieldStrength[__,devs_]:> Power[$C2Euc,Length@devs+2] fs} /.$C2Euc->I;
	
	(* Every spatial LCTensor adds a factor of -I *)
	out= out /.lc_LCTensor :> -$C2Euc lc /.$C2Euc->I;

	(* Recover the factor of \[CapitalTau] *)
	out= 1/\[CapitalTau] out /.f_Field :> Sqrt[\[CapitalTau]] f;
	
	(* Reescale of the gauge coupling and flip the sign of g^2 to have regular gauge kinetic terms *)
	out= out  /.($GaugeGroups[#,Coupling][]-> I Sqrt[1/\[CapitalTau]] $GaugeGroups[#,Coupling][]&/@SpGaugeGroups[]);
	
	(* Flip overall sign to match conventions  *)
	-out 
]


(* ::Section:: *)
(*SeriesThermalEFT*)


(* ::Subsubsection::Closed:: *)
(*DetermineThermalPowerCounting*)


(* ::Text:: *)
(*Function to determine the thermal power counting of an effective coupling if it is not determined already*)


DetermineThermalPowerCounting[label_]:= Module[{effexpr, thermalcount= $CouplingAssociation[label,ThermalPowerCounting],\[Lambda]}, 
	
	(*Determine if it does not have an associated value*)
	If[!NumberQ@thermalcount,
	effexpr= Cases[$TempCouplingRules, Rule[Coupling[label,__], term_]->term, All];
	thermalcount= effexpr/. {g: Coupling[x_,__]:> Power[\[Lambda],DetermineThermalPowerCounting[x]]*g};
	
	
	If[FreeQ[thermalcount, \[Lambda]], 
		thermalcount= 0,
		thermalcount= Assuming[\[Lambda]>0, thermalcount/. \[Lambda]-> (\[Lambda]+O[\[Lambda]]^2)/.sd_SeriesData:> sd[[4]]/sd[[6]]][[1]];
	];
	(*Give the value*)
	ThermalPowerCounting[label]:= thermalcount;
	];
	
	thermalcount	
]


(* ::Subsubsection::Closed:: *)
(*Error messages*)


SeriesThermalEFT::notarationale = "The order for the EFT series expansion specified `1` is neither an integer nor a rationale.";
SeriesThermalEFT::fatalpwrcounting = "The ThermalPowerCounting expecified `1` is not a list of the form {g->1,\[Lambda]->2,...}";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[SeriesThermalEFT]={ThermalPowerCounting -> 6};


(* ::Subsubsection::Closed:: *)
(*Main function*)


SeriesThermalEFT[arg_,OptionsPattern[]]:=Module[
	{
		\[Lambda],
		expr  = LagrangianExpand[arg],
		order= OptionValue@ThermalPowerCounting,
		result
	},
		
	expr= expr//RemovePower;
	
	(*Look if it is necessary to determine the thermal power counting*)
	expr= expr/.g: Coupling[label_,__]:> Power[\[Lambda],DetermineThermalPowerCounting[label]]*g;
	expr=ReleaseHold[expr];	
	
	Switch[order,
		(* select terms up to the given order *)
		_Integer | _Rational,
			If[Head[expr]===Plus,
				result = BetterSeries[#,{\[Lambda],0,order}]&/@expr
				,
				result = BetterSeries[expr,{\[Lambda],0,order}]
			],
		(* select terms of the exact order *)
		_,
			Message[SeriesThermalEFT::notarationale, order];
			Abort[]
	];
	result/.\[Lambda]->1
]


(* ::Section:: *)
(*SuperSoft Matching*)


(* ::Text:: *)
(*Thermalhbar is introduced to retain formal two-loop effects in SuperSoft matching (one-loop DR +  one loop SuperSoft).  Only one-loop DR terms are retained, i.e., Thermalhbar^2 ->0.*)


Thermalhbar^x_/;x>=2 ^:= 0


(* ::Subsection:: *)
(*Introduce and Restore Effective Gauge Couplings*)


$4DGaugeCouplings={};


(* ::Subsubsection::Closed:: *)
(*Options*)


Options @ IntroduceEffectiveGaugeCouplings = {EffectiveCouplingSymbol -> "\[ScriptG]", OverrideDuplicateCouplingCheck -> True};


(* ::Subsubsection::Closed:: *)
(*Main function*)


IntroduceEffectiveGaugeCouplings[L_ , OptionsPattern[]] := Module[{Lkin, Lrest,effcoups},
	(* take only gauge-kinetic terms *)
	Lkin = SelectOperatorClass[L,{},4];
	Lrest = L-Lkin//GreensSimplify;
	
	Lkin = HcSimplifyInternal[ Lkin, InternalOpRepresentation-> True];
	Lkin = If[Head@Lkin===Plus, List@@Lkin,{Lkin}];
	
	(*Introduce effective couplings*)
	Lkin  = Lkin  /. (x:(_AtomicOp|_CompOp)) * c_ :> x ToEffectiveCoupling[c/.hbar->Thermalhbar, x, Superleading -> False, Internal -> False, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol], OverrideDuplicateCouplingCheck -> OptionValue[OverrideDuplicateCouplingCheck] ];
	
	(*Return to normal form*)
	Lkin=Total@RelabelIndices@AtomicToNormalForm @ Lkin;
	
	(*List of effective gauge couplings with its associated gauge group*)
	effcoups=Cases[Lkin,-(FieldStrength[B_,__]^2/(4 Coupling[g_,{},0]^2))->{g,B},All];
	
	(*Keep track of the original 4D gauge couplings for later*)
	$4DGaugeCouplings={#,$GaugeGroups[#,Coupling]}&/@Keys@$GaugeGroups;
	
	(*Make the new coupling a dynamical gauge coupling for the soft match*)
	Do[
		With[{g=coups[[1]],B=coups[[2]]},
			group=First@*Keys@Select[$GaugeGroups,#[Field]===B &];
			$GaugeGroups[group,Coupling]=g;
		]
	,{coups,effcoups}];
	
	 Lkin+Lrest 
]


(* ::Subsubsection::Closed:: *)
(*RestoreEffectiveGaugeCouplings*)


RestoreEffectiveGaugeCouplings[]:= Do[
	$GaugeGroups[gcouplings[[1]],Coupling]=gcouplings[[2]]; 
	,
	{gcouplings,$4DGaugeCouplings}
]


(* ::Subsection:: *)
(*Introduce Debye Masses*)


(* ::Subsubsection::Closed:: *)
(*Options*)


Options @ IntroduceDebyeMasses = {EffectiveCouplingSymbol -> "M", OverrideDuplicateCouplingCheck -> True};


(* ::Subsubsection::Closed:: *)
(*Main function*)


IntroduceDebyeMasses[expr_ , OptionsPattern[]] := Module[{L=expr//ReplaceEffectiveCouplings,gaugescalars,Lkin,LB0kin,MB0,Lrest},
	gaugescalars=TemporalComponent/@SpGaugeGroups[];
	Lrest=L;
	Lkin=0;
	
	Do[
		$FieldAssociation[B0,Heavy]=True;
		LB0kin= SelectOperatorClass[L,{B0,B0},0];
		Lrest= Lrest-LB0kin;
		LB0kin= InternalSimplify[LB0kin,InternalOpRepresentation->True]/. (x:(_AtomicOp|_CompOp))*c_ :>  -(1/2)*x Sow@( ToEffectiveCoupling[Sqrt[-2*c/.$TempCouplingRules/.hbar->Thermalhbar], x, Superleading -> True, Internal -> False, EffectiveCouplingSymbol -> OptionValue[EffectiveCouplingSymbol], OverrideDuplicateCouplingCheck -> False, FactorOutHbar -> True ]^2);
		MB0= First@Cases[LB0kin,Coupling[M_,__]^2->M];
		DetermineThermalPowerCounting[MB0];
	
		$FieldAssociation[B0,Mass]=MB0;
	
		Echo["An effective Debye mass: "<>ToString[MB0[]//NiceForm]<>" has been introduced for the temporal gauge field "<>ToString[B0]<>"\n" ];
	
		LB0kin= RelabelIndices@AtomicToNormalForm @LB0kin;
	
		Lkin= Lkin+LB0kin;	
	,{B0,gaugescalars}];
	
	Lkin+Lrest
]


(* ::Subsubsection::Closed:: *)
(*Restore Debye Masses*)


RestoreDebyeMasses[]:=Module[{gaugescalars=TemporalComponent/@SpGaugeGroups[]},
	Do[
		$FieldAssociation[B0,Heavy]=False;
		$FieldAssociation[B0,Mass]=0;
		$currentFieldAssociation[B0,Heavy]=False;
		$currentFieldAssociation[B0,Mass]=0;
	,{B0,gaugescalars}];
]


(* ::Subsection:: *)
(*SuperSoft Match*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


SuperSoftMatch::NotZeroModes= "The Lagrangian contains Lorentz indices and/or fields that are not zero modes. Please run DRMatch first."


(* ::Subsubsection::Closed:: *)
(*Options*)


Options@ SuperSoftMatch = { EFTOrder -> 6 , Verbose -> True, ThermalPowerCounting-> 6};


(* ::Subsubsection::Closed:: *)
(*Main function*)


SuperSoftMatch[expr_, OptionsPattern[]]:=Module[{out=LagrangianExpand@expr,constants},
	(*Check that the input consists only of zero modes.*)
	If[!FreeQ[expr, Lorentz] || Complement[OccurringFields[out], FieldByProperty[ZeroMode->True]]=!={},
		Message[SuperSoftMatch::NotZeroModes];
		Abort[]
	];
	
	(*Suggest the relevant literature*)
	AddToBibliography["Matchotter","Using supersoftMatch"];
	
	(*Canonically normalize the kinetic terms*)
	{constants, out}= SeparateOutConstants@out;
	out=out/.\[Epsilon]^-1->0//EOMSimplify;
	
	(*Introduce Debye masses*)
	out= out//IntroduceDebyeMasses;
	
	(*Introduce effective gauge couplings*)
	out= out//IntroduceEffectiveGaugeCouplings;	
	
	out=out/.hbar-> Thermalhbar;
	
	(*Perform the matching with the usual code (see hack below)*)
	out= Match[out, EFTOrder-> OptionValue@ EFTOrder]//ReplaceEffectiveCouplings;
	RestoreEffectiveGaugeCouplings[];
	RestoreDebyeMasses[];
	
	SeriesThermalEFT[out + constants, ThermalPowerCounting-> OptionValue@ThermalPowerCounting]//CollectOperators
]


(* ::Section:: *)
(*Loop (sum)integrals*)


(* ::Subsection:: *)
(*SumIntegrals*)


(* ::Subsubsection::Closed:: *)
(*General loop integration function*)


(* ::Text:: *)
(*Performs  the  1 - loop  sumintegral  over  an  expression . Assumes  that  all  loop - momenta  occurs  in ' Prop' ,' PropTp' and 'PropSp'  objects . *)


LoopSumIntegrate::onescale="There is more than one scale in the LogTerm integration."


Options@ LoopSumIntegrate= {LogTerm-> False};


LoopSumIntegrate[expr_, opt:OptionsPattern[]]:=Module[{out},
	out= expr/.{Power[_Prop|_PropTp|_PropSp,n_/;IntegerQ[n]===False]-> 0};
	out= LagrangianExpand[out* sumintegralType[0, 0, 0]];
	out= out/. {$FSumInt x_sumintegralType:> SolveSumInt[Sequence@@x, $Dimensions, Type->Fermion]}/. {x_sumintegralType:> SolveSumInt[Sequence@@ x, $Dimensions]};
	out= out
];


(* ::Subsubsection::Closed:: *)
(*Identifying propagators in the expressions*)


(* ::Text:: *)
(*mprops  is  an  association  of  mass  propagators  and  their  corresponding  power, for  example  a  denominator  1/((p^2 - m1^2)^n1 (p^2 - m1^2)^n2)  reads  <|m1 -> n1, m2 -> n2|>  in  mprops .*)
(*  \[Alpha]  is  the  power  of  (1/p^2)^\[Alpha] .*)


sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] Prop@ 0:= sumintegralType[ \[Alpha]+ 1, \[Beta], \[Gamma]];
sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] Power[Prop@ 0, n_]:= sumintegralType[\[Alpha]+ n, \[Beta], \[Gamma]];
sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] PropTp@ 0:= sumintegralType[\[Alpha], \[Beta]+1, \[Gamma]];
sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] Power[PropTp@ 0, n_]:= sumintegralType[\[Alpha], \[Beta]+n, \[Gamma]];
sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] PropSp@ 0:= sumintegralType[\[Alpha], \[Beta], \[Gamma]+1];
sumintegralType/: sumintegralType[\[Alpha]_, \[Beta]_, \[Gamma]_] Power[PropSp@ 0, n_]:= sumintegralType[\[Alpha], \[Beta], \[Gamma]+ n];


(* ::Subsubsection::Closed:: *)
(*Solve SumIntegrals*)


Options@ SolveSumInt= {Type->Scalar};


SolveSumInt[\[Alpha]_,\[Beta]_,\[Gamma]_,dim_,opts:OptionsPattern[]]:= SolveSumInt[\[Alpha],\[Beta],\[Gamma],dim,opts]=
Module[{PreFact, FerPrefact, FerTerm, GammaPreFact, LogCoeff, d=dim, mpFact, type=OptionValue@Type},
	PreFact=(-1)^(\[Alpha]+\[Beta]+\[Gamma]) (4\[Pi])^(2-(d-1)/2) 2\[CapitalTau] (\[CapitalLambda]\[CapitalTau]/2)^((d-1)-2(\[Alpha]+\[Beta]+\[Gamma]))/(Gamma[\[Alpha]] Gamma[(d-1)/2]);
	GammaPreFact=Gamma[(d-1)/2-\[Gamma]]Gamma[\[Alpha]+\[Gamma]-(d-1)/2];
	
	FerPrefact=If[type===Fermion, 2^(2(\[Alpha]+\[Beta]+\[Gamma])-(d-1))-1, 1];
	FerTerm=If[type===Fermion, 2^(2(\[Alpha]+\[Beta]+\[Gamma])-(d-2))Log[2], 0];
	LogCoeff=If[type===Fermion, \[CapitalLambda]\[CapitalTau]^2/16, \[CapitalLambda]\[CapitalTau]^2];
	
	If[EvenQ[d],
		If[2(\[Alpha]+\[Beta]+\[Gamma])-(d-1)===1, (* Expand around \[Zeta] pole*)
			Return[ PreFact GammaPreFact 1/2(1/\[Epsilon]+Log[\[CapitalLambda]bar2/LogCoeff]+FunctionExpand[Log[4]+HarmonicNumber[(d-3)/2 ]-HarmonicNumber[(d-3)/2- \[Gamma]]+HarmonicNumber[\[Alpha]+\[Gamma]-(d+1)/2]]) ]
		]
	,
		If[\[Alpha]+\[Gamma]<= (d-1)/2, (* Expand around Gamma pole *) 
			If[Zeta[2(\[Alpha]+\[Beta]+\[Gamma])-(d-1)]===0,		
				Return[ PreFact FerPrefact Gamma[(d-1)/2- \[Gamma]](-1)^(\[Alpha]+\[Gamma]-(d-1)/2)/((d-1)/2-\[Alpha]-\[Gamma])! 2 FunctionExpand@Derivative[1][Zeta][2 (\[Alpha]+\[Beta]+\[Gamma])-(d-1)]]
			,						
				Return[ PreFact Gamma[(d-1)/2- \[Gamma]](-1)^(\[Alpha]+\[Gamma]-(d-1)/2)/((d-1)/2-\[Alpha]-\[Gamma])! *
					Zeta[2(\[Alpha]+\[Beta]+\[Gamma])-(d-1)] (FerPrefact (1/\[Epsilon]+ Log[\[CapitalLambda]bar2/\[CapitalLambda]\[CapitalTau]^2]+Simplify[-Log[16\[Pi]^2]+HarmonicNumber[(d-3)/2]-HarmonicNumber[(d-3)/2-\[Gamma]]+Log[4]+2Derivative[1][Zeta][2 (\[Alpha]+\[Beta]+\[Gamma])-(d-1)]/Zeta[2(\[Alpha]+\[Beta]+\[Gamma])-(d-1)]]) + FerTerm) ]	
			]
		]
	];

	Return[ PreFact FerPrefact GammaPreFact Zeta[2(\[Alpha]+\[Beta]+\[Gamma])-(d-1)] ]
	
(* Remove factors of \[Pi] *)
]/. \[Pi]-> \[CapitalLambda]\[CapitalTau]/(4\[CapitalTau]);
