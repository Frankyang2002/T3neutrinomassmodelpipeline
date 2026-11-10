(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`ExportUFO`*)


(* ::Subtitle:: *)
(*Create UFO model for Matchete.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["ExportUFO"]


PackageExport["DefaultParamCard"]


PackageExport["InputFile"]
PackageExport["OutputDirectory"]


PackageExport["DefineCouplingOrder"]
PackageExport["GetCouplingOrders"]
PackageExport["RemoveCouplingOrder"]
PackageExport["ResetCouplingOrders"]


PackageExport["UFO"]


(* ::Subsubsection::Closed:: *)
(*Scoped*)


PackageScope["ExplicitEinsteinSums"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


ExportUFO::usage= "ExportUFO[lag] computes the Feynman rules for the given Lagrangian lag and generates the corresponding UFO files.";


DefaultParamCard::usage= "DefaultParamCard[lag] creates a template for the parameter card for the given Lagrangian lag, which is saved as a JSON file in the current working directory or, if specified, in the directory dir given as option OutputDirectory \[Rule] dir.";


InputFile::usage= "InputFile is an Option of ExportUFO that allows to specify the path to a parameter input file using the JSON format.";
OutputDirectory::usage= "OutputDirectory is an Option of ExportUFO that allows to specify the output directory where the UFO files should be saved.
OutputDirectory is an Option of DefaultParamCard that allows to specify the output directory where the parameter card should be saved.";


DefineCouplingOrder::usage= "DefineCouplingOrder[name, couplings, expansionOrder, hierarchy, opts] defines a coupling order based on the given coupling(s). The expansion order and hierarchy should be stated for numerical reasons.";
GetCouplingOrders::usage= "GetCouplingOrders[] returns the association of all couplings orders that are defined.
GetCouplingOrders[ord] returns the information on the couplings order ord.";
RemoveCouplingOrder::usage= "RemoveCouplingOrder[ord] removes the definitions of the coupling order ord.";
ResetCouplingOrders::usage= "ResetCouplingOrders[] removes all previously defined coupling orders.";


UFO::usage= "UFO is an option of DefineField, DefineCoupling, and DefineGaugeGroup that can be used to provide auxiliary information that is used by DefaultParamCard and ExportUFO but not needed for the core Matchete functionality. If the UFO option is not provided, default values will be automatically determined based on the filed and coupling labels used in Matchete. The only exception are PDG codes for particles, which are mandatory for generating UFO files.";


(* ::Chapter:: *)
(*Private:*)


(* ::Section:: *)
(*Main UFO functions*)


(* ::Subsection:: *)
(*ExportUFO*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


(* move to Matchete.m ? *)
(*OptionMessage[Mode, ExportUFO, val_]:= Message[General::optexpectsval, Mode, ExportUFO, val, "\"Standard\" or \"MadGraph\""];*)


ExportUFO::nofile      = "Parameter file `1` not found.";
ExportUFO::missingdata = "Definitions for the following parameters are missing: `1`. MATCHETE only checks for parameters of the lagrangian (couplings, etc.)";
ExportUFO::missingpdg  = "pdg code for particle `1` is missing. Please define fields with DefineField[..., UFO -> <|\"pdg\" -> ...|>] and gauge fields with DefineGaugeGroup[..., UFO -> <|\"pdg\" -> ...|>].";


ExportUFO::dir = "Failed to create the UFO directory at `1`."


ExportUFO::multGroups = "The following non-Abelian groups are present in the given Lagrangian: `1`. Only a single SU(3) gauge group is allowed, representing QCD.";
ExportUFO::noQCD      = "The non-Abelian group `1` present in the given Lagrangian has the algebra `2`, which does not correspond to a SU(3) group (i.e. QCD) for which we expect `3`.";
ExportUFO::gluonPDG   = "The gauge boson `3` of the SU(3) group `1` present in the given Lagrangian has the PDG code `2`, which does not correspond to the gluon PDG code of 21.";
ExportUFO::noGroupID  = "Could not itentify the `1` gauge group!"


ExportUFO::missingPDG = "The PDG codes for the following fields have not been defined: `1`.";


ExportUFO::nongf= "It looks like the input Lagrangian has not been gauge-fixed despite using Gauge-> Manual option. This is likely an error."


ExportUFO::defaultParams = "Warning: Using the default parameter card which initialized all Lagrangian parameters to zero and does not provide a specific parameterization. Use DefaultParamCard[lagrangian, OutputDirectory\[Rule]\"path\"] to obtain a template of this parameter card which you can modify as desired and then used by calling ExportUFO[lagrangian, InputFile\[Rule]\"path/model_parameters.json\"].";


ExportUFO::numext = "External parameters must be numbers. The parameter `1` was given the value `2`."


ExportUFO::attr = "The parameter `1` does not have all required attributes. The required attributes are: \"name\", \"nature\", \"type\", \"value\", \"texname\", \"lhablock\", and \"lhacode\". The last two are only required for external parameters."
ExportUFO::extAttr = "The external parameter `1` does not have all required attributes. The required attributes are: \"name\", \"nature\", \"type\", \"value\", \"texname\", \"lhablock\", and \"lhacode\"."
ExportUFO::intAttr = "The internal parameter `1` does not have all required attributes. The required attributes are: \"name\", \"nature\", \"type\", \"value\", and \"texname\"."


ExportUFO::paramSort = "Sorting of the parameter card failed. Most likely the param card is inconsistent.";


ExportUFO::paramRemove = "Removed parameters `1` that were defined in the input file but are not used in the model.";


ExportUFO::chargeconj = "The following 4-fermion vertices with charge conjugation matrices were found: `1`. Such vertices might not be fully supported by MadGraph.";


ExportUFO::massBasis = "Field `1` was not defined in the mass basis. Nevertheless, ExportUFO assumes it to be in the mass basis."


ExportUFO::defaultValues = "Ignoring the default SM input values. DefaultValues\[Rule]True is only allowed when InputFile\[Rule]None is given."


ExportUFO::ghosts = "Ghost fields detected in the input. These might have been introduced by a call of GaugeFixLagrangian. Currently only tree-level UFO files can be generated with Matchete. Dropping all ghost terms\[Ellipsis]"


ExportUFO::frgauge = "Matchete is not applying any gauge fixing when the Feynman rules are provided manually as the gauge fixing has to be performed before deriving the Feynman rules. Consider using the option Gauge\[Rule]Manual for ExportUFO to suppress this message."


ExportUFO::ui = "The UI for modifying the parameter card was closed. UFO files are not generated.";


ExportUFO::vanishingD4 = "The following terms vanish in D=4 spacetime dimensions and are removed from the UFO files: `1`";


ExportUFO::unphysgamma = "Fermionic currents with more than two Gamma matrices detected. Consider reducing them to a physical four-dimensional basis before exporting to the UFO format.";


ExportUFO::rxiwidth= "The current implementation of widths is not fully consistent in \!\(\*SubscriptBox[\(R\), \(\[Xi]\)]\) gauges and can break gauge invariance (particularly for large \[Xi] values)."


ExportUFO::flavinvcard= "The option FlavorInvariants is ignored whenever a parameter card is provided through the option InputFile.";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[ExportUFO] = {
	DefaultValues   -> False,
	FeynmanRules    -> None, 
	FlavorInvariants-> <||>,
	Gauge           -> Automatic, (* Manual | Unitary | R\[Xi] *)
	Information     -> <|
		"author"        -> "Matchete",
		"model_version" -> "0.0",
		"arxiv"         -> "0000.00000"
	|>,
	InputFile       -> None,
	Legs            -> All,
	(*Mode            -> "MadGraph", (* "MadGraph" | "Standard" *)*)
	OutputDirectory :> Automatic,
	R\[Xi]              -> 1, (* General *)
	Verbose         -> True
	,
	"ExternalMomentum"   -> "outgoing",
	"FermionFlow"        -> "DEHK",
	"IncludeD4Vanishing" -> False,
	"\[Sigma]Normalization"     -> I/4
};


OptionTest[ExportUFO, "ExternalMomentum"]   = MatchQ["incoming"|"outgoing"];
OptionTest[ExportUFO, "FermionFlow"]        = MatchQ["Standard"|"DEHK"];
OptionTest[ExportUFO, "IncludeD4Vanishing"] = BooleanQ;
OptionTest[ExportUFO, "\[Sigma]Normalization"]     = NumberQ;


OptionMessage["ExternalMomentum", ExportUFO, val_]   := Message[General::optexpectsval, "ExternalMomentum", ExportUFO, val, "\"incoming\" or \"outgoing\""];
OptionMessage["FermionFlow", ExportUFO, val_]        := Message[General::optexpectsval, "FermionFlow", ExportUFO, val, "\"Standard\" or \"DEHK\""];
OptionMessage["IncludeD4Vanishing", ExportUFO, val_] := Message[General::optexpectsval, "IncludeD4Vanishing", ExportUFO, val, "boolean"];
OptionMessage["\[Sigma]Normalization", ExportUFO, val_]     := Message[General::optexpectsval, "\[Sigma]Normalization", ExportUFO, val, "number"];


(* ::Subsubsection::Closed:: *)
(*ExportUFO*)


(* ::Text:: *)
(*Storage of the current QCD and QED group names*)


$QCD$UFO= Nothing
$QED$UFO= Nothing


ExportUFO[lagrangian_, opt:OptionsPattern[]] ? OptionsCheck := Module[
	{
		vertices,
		allFields,
		vertexData = {},
		dir,
		fields,
		fieldTypes,
		expandedVertex,
		indexTypes,
		lag= lagrangian,
		ghostLag,
		lorentzStructure,
		groupFlavorStructure,
		nonGFed,
		structure,
		couplings,
		lorentzData,
		couplingData,
		defaultCardParameters,
		parametersNeeded,
		parametersExternal,
		parametersFile = OptionValue[InputFile],
		fieldTypeTranslation = <|"V" -> 3, "S" -> 1, "F" -> 2, "U" -> -1|>,
		lorentzRules,
		brokenSymmetries,
		couplingRules,
		compFields,
		rulesCompFields,
		paramsZERO, zeroCouplingsRule,
		defaultCard,
		\[Sigma]Rescaling= (I/2)/OptionValue["\[Sigma]Normalization"], (* rescale \[Sigma][\[Mu],\[Nu]] by a factor of 2 for MadGraph conventions, i.e. the default for OptionValue["\[Sigma]Normalization"] is \[ImaginaryI]/4 *)
		(*\[Sigma]Rescaling= If[OptionValue[Mode]==="MadGraph", 2, 1],*) (* rescale \[Sigma][\[Mu],\[Nu]] by a factor of 2 for MadGraph conventions *)
		metadata,
		info= OptionValue[Information],
		initFile,
		countCC, 
		fieldLabels, 
		ghostMessage= False,
		vanishingDiracStructures,
		propInfo,
		gaugeInfo, gaugedata,
		massiveVectors, droppedGBs
	}
	,
	(* using FlavorInvariants is only possible without a param card *)
	If[(OptionValue[FlavorInvariants]=!=<||>) && (OptionValue[InputFile]=!=None),
		Message[ExportUFO::flavinvcard]
	];
	
	(*Check the Lagrangian*)
	If[OptionValue@ FeynmanRules === None && !CheckLagrangian[lag, 
				CanonicallyNormalized -> False, 
				GaugeAnomalies        -> False,
				GaugeInvariant        -> False,
				HeavyMassBasis        -> False,
				Hermiticity           -> False,
				UndefinedObject       -> False
			], 
		Abort[]
	];
	
	(* check and (if required) include gauge-fixing terms *)
	{lag, ghostLag, gaugeInfo}= GaugeFixingUFO[lag, OptionValue[Gauge], OptionValue[R\[Xi]]];
	
	(* print warning for R\[Xi] gauges *)
	If[OptionValue[R\[Xi]]===General,
		Echo[ExportUFO::rxiwidth]
	];
	
	(*
	(* throw warning if the Lagrangian is not gauge fixed yet *)
	gaugeInfo= NonGaugeFixedLagQ[lag, DetailedOutput-> True];
	If[First@ gaugeInfo,
		Message[DefaultParamCard::gaugefix];
		Abort[]
	];
	*)
	
	(*
	(* check for ghost fields and remove them *)
	If[!FreeQ[lag, Ghost|AntiGhost],
		(*Message[ExportUFO::ghosts];*)
		Echo["Dropping ghost terms for tree-level UFO..."];
		ghostMessage= True; (* just a flag whether this message has been printed already *)
		(* drop the ghosts since we work at tree level for now *)
		lag= lag/. Field[_,Ghost|AntiGhost, ___]-> 0;
	];
	*)
	
	(* check that pdg codes of all particles have been defined *)
	CheckPDGcodes[lag];
	
	(* determine QCD and QED group *)
	IdentifyQCDandQED[lag];
	
	(* Throw a warning if a field was not defined in its mass basis *)
	fieldLabels= Intersection[
		DeleteDuplicates@ Cases[lag, Field[l_,___]->l, All],
		GetFieldsByProperty[Mass->Except[0]]
	];
	Do[
		(* check that the mass is either flavorless or fully diagonal *)
		If[!MatchQ[GetCouplings[GetFields[f,Mass], DiagonalCoupling], {}|{True..}],
			Message[ExportUFO::massBasis, f]
		]
		,
		{f, fieldLabels}
	];
	
	(*
	(* Gauge must be fixed before now *)
	(*Check if Lagrangian is or should be gauge-fixed*)	
	If[OptionValue@ FeynmanRules === None,
		nonGFed= NonGaugeFixedLagQ@ lag;
		Switch[{OptionValue@ Gauge, nonGFed}
		,{Manual, True},
			Null;
		,{Manual, False},
			Message[ExportUFO::nongf];
		,{_, True},
			lag= GaugeFixLagrangian[lag, Gauge-> OptionValue@ Gauge];
		,{_, False},
			Abort[];
		]
		,
		If[OptionValue@ Gauge =!= Manual,
			Message[ExportUFO::frgauge]
		]
	];
	*)
	
	(*
	(* drop the ghosts since we work at tree level for now *)
	lag= lag/. Field[_,Ghost|AntiGhost, ___]-> 0;
	*)
	
	(* check for charge conjugation matrices in (n\[GreaterEqual]4)-fermion operators *)
	countCC= DeleteDuplicates@ Cases[Operator@ lagrangian, _Operator, All];
	countCC= Cases[countCC, x_/; (Count[x, _NCM, All]>=2 && !FreeQ[x,GammaCC,All])];
	If[Length[countCC]>0,
		Message[ExportUFO::chargeconj, Format[countCC, NiceForm]]
	];
	
	(* check for Dirac structures vanishing in D=4 dimensions *)
	If[!TrueQ[OptionValue["IncludeD4Vanishing"]] && OptionValue[FeynmanRules]===None,
		vanishingDiracStructures= Cases[
			lag, 
			x_/; !FreeQ[x, NCM[Bar[_Field],DiracProduct[GammaM[Index[\[Mu]_,Lorentz],Index[\[Nu]_,Lorentz]],Proj[-1]],_Field]NCM[Bar[_Field],DiracProduct[GammaM[Index[\[Mu]_,Lorentz],Index[\[Nu]_,Lorentz]],Proj[+1]],_Field]],
			1
		];
		If[Length[vanishingDiracStructures]>0,
			Message[ExportUFO::vanishingD4, Format[Short[Total@ vanishingDiracStructures, 30], NiceForm]];
			(* 
			it might be better to use:
				lag= CollectOperators@ GammaReduction[lag, FourDimensional-> True, Evanescent-> False]
			though it should only be applied to the problematic terms
			*)
			lag= lag/. {
				NCM[Bar[_Field],DiracProduct[GammaM[Index[\[Mu]_,Lorentz],Index[\[Nu]_,Lorentz]],Proj[-1]],_Field]NCM[Bar[_Field],DiracProduct[GammaM[Index[\[Mu]_,Lorentz],Index[\[Nu]_,Lorentz]],Proj[+1]],_Field]-> 0
			}
		];
	];
	
	(* create the UFO output directory, copy default files there, and set as working directory *)
	If[OptionValue[OutputDirectory]===Automatic,
		dir= FileNameJoin@ {Directory[], "UFO_output"},
		dir= OptionValue[OutputDirectory]
	];
	dir= CreateNewDir[dir];
	If[dir===$Failed,
		Message[ExportUFO::dir, dir];
		Abort[]
	];
	CopyDefaultUFOFiles[dir];
	SetDirectory[dir];
	
	(* this is the metadata string for __init__.py *)
	metadata= StringJoin[
	  "__author__         = \"", "author"/.info, "\"\n",
	  "__date__           = \"", DateString[{"Day",".","Month",".","Year"}], "\"\n",
	  "__model_version__  = \"", "model_version"/.info, "\"\n",
	  "__arxiv__          = \"", "arxiv"/.info, "\"\n",
	  "__UFO_version__    = \"2.0\"\n",
	  "__python_version__ = 3\n",
	  "\n"
	];
	
	(* determine the gauges with which this UFO can be used and add the to __init__.py *)
	massiveVectors= MassiveVectorsInLag[lag, Light-> True];
	gaugedata= If[ContainsOnly[Last[gaugeInfo]/@ massiveVectors, {Unitary}],
		"# 0: unitary only
gauge = [0]"
	,
		If[OptionValue[R\[Xi]]===1,
			"# 0: unitary (MadGraph drops goldstones and uses the default unitary propagators internally);
# 1: feynman
gauge = [0, 1]"
		,
			"# 1: feynman only -> avoids deletion of goldstones while propagators have to be overwritten so feynman gauge is effectively not used; general Rxi gauge cannot be reduced to unitary gauge (for xi=!=1)
gauge = [1]"
		]
	];
	gaugedata= "\n\n"<>gaugedata<>"\n";

	
	initFile= Import["__init__.py", "Text"];
	initFile= metadata <> initFile <> gaugedata;
	Export["__init__.py", initFile, "Text"];
	
	(* import parameter file if specified, otherwise use default param card*)
	If[parametersFile === None,
		parametersExternal= DefaultParamCard[ghostLag (*lag is not gauge invariant*), "Internal"->True, DefaultValues-> OptionValue[DefaultValues], Gauge-> Manual, R\[Xi]-> OptionValue[R\[Xi]], FlavorInvariants-> OptionValue[FlavorInvariants]]/. Association-> List;
		parametersExternal= parametersExternal/. ("lhacode" -> {x_}) :> ("lhacode" -> x);
		defaultCard= True; (* flag for whether DefaultParamCard is used *)
		(*Message[ExportUFO::defaultParams];*)
		,
		If[OptionValue[DefaultValues], Message[ExportUFO::defaultValues]];
		If[FileExistsQ[parametersFile],
			parametersExternal= Import[parametersFile]/. ("lhacode" -> {x_}) :> ("lhacode" -> x);
			defaultCard= False; (* flag for whether DefaultParamCard is used *)
			,
			Message[ExportUFO::nofile, parametersFile]; 
			Abort[]
		]
	];
	
	(* determine minimal set of required input parameters *)
	{defaultCardParameters, droppedGBs}= Reap[DefaultParamCard[ghostLag (*lag is not gauge invariant*), "Internal"->True, Gauge-> Manual, R\[Xi]-> OptionValue[R\[Xi]],FlavorInvariants-> OptionValue[FlavorInvariants]]["parameters"]];
	(* dropped GBs is the list of GBs and their gauge parameters that require a custom propagtor in R\[Xi] *)
	
	(*Write parameters.py; returnes the list of parameters that should be set to zero entirely*)
	OptionalMonitor[OptionValue@Verbose,
		paramsZERO= WritePyParameters[defaultCardParameters, parametersExternal, defaultCard];
	,"Writing the parameters.py \[Ellipsis]"];
	
	(*Replace leftover CGs*)
	lag = lag// ReplaceCGs;

	(* calculate Feynman rules *)
	Module[{unphysGamma},
	unphysGamma= Reap@ If[OptionValue[FeynmanRules] === None,
		OptionalMonitor[OptionValue@Verbose,
			(*  
				- flip sign of all Levi-Civitas to map from Matchete to UFO convention 
				- expand anti-sym. GammaM with more than 2 indices since these cannot be handled by the UFO format
				- in "MadGraph" mode rescale GammaM[_,_] (i.e. \[Sigma][\[Mu],\[Nu]]) by a factor of 2 to match MadGraph conventions 
			*)
			vertices = FeynmanRules[
				lag/. (lct_LCTensor:> -lct)/. (dp_DiracProduct/; !FreeQ[dp, HoldPattern@GammaM[_,_,__], All]:> (Sow[dp]; ASymGammaExpand[dp]))/. (sigma:(GammaM[\[Mu]_,\[Nu]_]):> \[Sigma]Rescaling*sigma),
				Legs               -> OptionValue@ Legs,
				"ExternalMomentum" -> OptionValue@ "ExternalMomentum",
				"FermionFlow"      -> OptionValue@ "FermionFlow"
			];
		,"Deriving Feynman rules \[Ellipsis]"];
	,
		Echo["Ignoring options \"ExternalMomentum\", \"FermionFlow\", and \"IncludeD4Vanishing\" when using the 'FeynmanRules' option."];
		vertices= OptionValue[FeynmanRules]/. lct_LCTensor:> -lct/. (dp_DiracProduct/; !FreeQ[dp, HoldPattern@GammaM[_,_,__], All]:> (Sow[dp]; ASymGammaExpand[dp]))/. (sigma:(GammaM[\[Mu]_,\[Nu]_]):> \[Sigma]Rescaling*sigma);
		(* check for ghost fields and remove them *)
		If[!FreeQ[vertices, Ghost|AntiGhost],
			If[!ghostMessage, 
				Echo["Dropping ghost terms for tree-level UFO..."];
				(*Message[ExportUFO::ghosts]; *)
				ghostMessage= True;
			];
			(* drop the ghosts since we work at tree level for now *)
			vertices= DeleteCases[vertices, x_/; !FreeQ[x, Ghost|AntiGhost, All]];
		]
	];
	If[Length[Last@ unphysGamma]> 0, Echo@ StringForm[ExportUFO::unphysgamma, Format[Short[Last@ unphysGamma, 30], NiceForm]]];
	];
	
	(* keep only FR w/o keys *)
	If[AssociationQ[vertices],
		vertices= Values[vertices]
	];

	(*Contract Flavor deltas with external legs for better performance*)
	vertices = ContractFlavorDelta/@ vertices;

	(* introduce explicit flavor indices *)
	OptionalMonitor[OptionValue@Verbose,
		vertices= DeleteCases[
			Flatten[ExplicitEinsteinSums[LagrangianExpand[#], List->True, Internal->False]&/@ vertices]
			,0,1
		];
		(* perform Einstein summation over flavor inidices *)
		(*vertices= DeleteCases[ExplicitEinsteinSums[LagrangianExpand[#], List->True], 0, 1]&/@ vertices;*)
		(* define components fields for flavor multiplets and replace the latter by the former *)
		{compFields,rulesCompFields}= IntroduceFlavorComponentFields[vertices];
		vertices= vertices/. rulesCompFields;
	,"Introducing explicit flavor indices in Feynman rules \[Ellipsis]"];
	
	(* drop ZERO couplings *)
	zeroCouplingsRule= Map[
		(#->(TranslateToUFO[#]/.((#->0)&/@paramsZERO)))&,
		(*DeleteDuplicates@Cases[vertices,c_Coupling,All]*)
		DeleteDuplicates@Cases[
			ExplicitEinsteinSums[#, List->False, Internal->True]& /@ Map[
				If[Head[#] === Plus, (Times@@Cases[#, c_Coupling, All]&) /@ #,
				Times@@Cases[#, c_Coupling, All]]&,
				vertices
				]
			, c_Coupling, All]
	];
	zeroCouplingsRule= Dispatch@Select[zeroCouplingsRule, #[[2]]===0&];
	vertices= DeleteCases[vertices/. zeroCouplingsRule, 0]; 
	
	
	(* get labels of all fields *)
	allFields = DeleteDuplicates@Cases[vertices, Field[Ext[x_,_], ___] :> (x/.(*OverBar|*)Conj->Identity), All];
	
	(* list of local & global symmetry groups as well as flavor index types defined *)
	indexTypes = DeleteDuplicates@ Replace[Flatten[$FieldAssociation[#][Indices]&/@ Keys[$FieldAssociation]], s_Symbol[___] :> s, {1}];

	(* loop over all verices *)
	Do[
		(* sorted list of all fields of given vertex *)
		fields         = SortBy[DeleteDuplicates@GetFieldsInTerm[vertex], Cases[#, Ext[_, int_Integer]:> int, \[Infinity]] &];
		(* list of corresponding field types *)
		fieldTypes     = StringTake[SymbolName[If[AtomQ[#], #, Head[#]]], 1]&/@ (fields /. Bar[x_] -> x)[[All, 2]] /. {"G"->"U", "A"->"U"};
		(* expanded form of the vertex *)
		expandedVertex = TermsToList[vertex, Expand->True];
		(* extract lorentz/spinor structures *)
		lorentzStructure = (Times@@ Join[GetDiracStructureSingleTerm[ReleasePseudoTimes@ #], Cases[#, _Metric|_ExtMom|_LCTensor, All]]&)/@ PseudoTimes/@ expandedVertex;
		expandedVertex   = expandedVertex/. _Field|_Metric|_ExtMom|_LCTensor|_DiracProduct -> 1;
		(* extract all remaining structures (gauge groups, flavor) *)
		(* loop over groups and associate factors related to them *)
		groupFlavorStructure = Association@ Table[
			ToString[indexType] -> Map[
				Function[
					term, Module[
						{factors},
						(* break the product into its factors if it is a Times product *)
						factors = If[Head[term] === Times, List@@term, {term}];
						(* recombine only those factors that contain the current index type *)
						Times@@ Select[factors, Not[FreeQ[#, indexType]] &]
					]
				],
				expandedVertex
			],
			{indexType, indexTypes}
		];
		(* replace dummy indices by negative integers *)
		groupFlavorStructure= Map[RenameDummyIndicesToUFO, groupFlavorStructure, {2}];
		
		(* combine structures *)
		structure = Append[groupFlavorStructure, "Lorentz"-> lorentzStructure];
		
		(* the remaining part of the Feynman rule is the coupling *)
		couplings = ExplicitEinsteinSums[#, List->False, Internal->True]& /@ (expandedVertex /. (_CG|_Delta)-> 1);
		couplings = couplings /. zeroCouplingsRule;
		AppendTo[structure, "Couplings"-> couplings];
		(* The line below only deletes couplings but not the corresponding entries in the other associations *)
		(*couplings = DeleteCases[couplings, 0];*)
		(* skip this structure if all couplings vanish [empty lists cannot be handled by the code] *)
		If[MatchQ[couplings, {0..}], Continue[]];
		(* remove all entries in structure where the corresponding couplings vanish *)
		structure= ReplacePart[#, Position[couplings, 0, 1]-> Nothing]&/@ structure; (* this has to be done before the additional info is added below, also ReplaceAt does not work *)
		
		(* add fields, ... to structure *)
		AppendTo[structure, "FieldsContent" -> StringJoin[fieldTypes]];
		AppendTo[structure, "Fields"        -> fields//. Field[Ext[x_,_], ___]:> x];
		AppendTo[structure, "FieldsSorted"  -> Replace[fieldTypes, fieldTypeTranslation, {1}]];
		AppendTo[structure, "Feynmanrule"   -> vertex];
		
		(* add final result to vertexData *)
		AppendTo[vertexData, structure]
		,
		{vertex, vertices}
	];
	
	(* drop field labels from momenta *)
	vertexData = vertexData /. ExtMom[\[Mu]_, Ext[_, n_]] :> ExtMom[\[Mu], n]; 

	vertexData = AddRedundantLorentzStructures /@ vertexData;

	(* create a list of unique couplings *)
	couplingData = DeleteDuplicates@ Flatten[
		Table[
			\[ScriptL],
			{ver,vertexData},
			{\[ScriptL], ver["Couplings"]}
		]
	];
	
	(* create a list of unique Lorentz structures *)
	lorentzData = DeleteDuplicates@ Flatten[
		Table[
			<|"Lorentz" -> \[ScriptL], "FieldsContent" -> ver["FieldsContent"], "FieldsSorted" -> ver["FieldsSorted"]|>
			,
			{ver, vertexData},(* loop over each association *)
			{\[ScriptL], ver["Lorentz"]}(* loop over each element of its Lorentz list *)
		]
	];
	
	(* add names to Lorentz structures *)
	lorentzData = Flatten@ KeyValueMap[
		Function[
			{fieldCounts, group},
			MapIndexed[
				Append[#, "Name" -> fieldCounts <> ToString[First[#2]]] &,
				group
			]
		],
		GroupBy[lorentzData, #["FieldsContent"]&]
	];
	
	(* replace each Lorentz structure in vertexData with its corresponding name *)
	lorentzRules = Association[
		Table[
			{l["Lorentz"], l["FieldsContent"]} -> l["Name"]
			,
			{l, lorentzData}
		]
	];
	vertexData = Map[
		Function[v,
			Module[
				{fieldsContent = v["FieldsContent"], oldL = v["Lorentz"], newL},
				(* produce a new list of names, one per old lorentz structure *)
				newL = Lookup[
					lorentzRules,
					Thread[{oldL, ConstantArray[fieldsContent, Length@oldL]}],
					oldL  (* fallback: keep original if no rule found *)
				];
				(* return a new association which is like v but with "Lorentz"->newL *)
				ReplacePart[v, "Lorentz" -> newL]
			]
		],
		vertexData
	];
	
	(* drop any broken symmetry, i.e., where no symmetry structures are present *)
	brokenSymmetries = Select[
		Keys[First[vertexData]],
		Function[symmetry, AllTrue[vertexData, MatchQ[# [symmetry], {1 ..}] &]]
	];
	vertexData = KeyDrop[#, brokenSymmetries] & /@ vertexData;
	
	(* edit coupling data *)
	(* The resulting list is of the form {{couplings},{names},{UFO names},{coupling orders}} *)
	couplingData = List[
		couplingData,
		"GC_" <> ToString[#]& /@ Range[Length@couplingData],
		TranslateToUFO /@ couplingData,
		StringJoin @ {
			"{", 
			StringRiffle[
				"'" <>
				#1 <> 
				"':" <> 
				ToString[#2, InputForm] & @@@ Select[Normal[ClassifyCouplingOrder[#]], #[[2]] =!= 0 &],
				", "
			],
			"}"
		}& /@ couplingData
	];
	
	(* replace each Coupling in the vertex data with its corresponding name *)
	couplingRules = AssociationThread[couplingData[[1]] -> couplingData[[2]]];
	vertexData    = Map[MapAt[Lookup[couplingRules, #, #] &, "Couplings"] @ # &, vertexData];

	(* translate from MATCHETE convention to UFO format *)

	lorentzData = MapAt[RenameDummyIndicesToUFO, #, Key["Lorentz"]]& /@ lorentzData;
	lorentzData = Map[MapAt[TranslateToUFO       , #, Key["Lorentz"]] &, lorentzData /. (index_Symbol /; StringMatchQ[SymbolName[index], "d$$"~~__] :> ToExpression @ StringDrop[SymbolName[index], 3])];
	vertexData  = Map[MapAt[TranslateToUFO/@ # & , #, Key[ToString@$QCD$UFO]] &   , vertexData /. (index_Symbol /; StringMatchQ[SymbolName[index], "d$$"~~__] :> ToExpression @ StringDrop[SymbolName[index], 3])];
	
	(*Remove redundancies in color structure*)
	vertexData = Map[RemoveColorRedundancies, vertexData];
	
	
	(*** WRITING REMAINING OUTPUT FILES ***)
	(*Write propagators.py*)
	propInfo= WritePYPropagators[lag, gaugeInfo, droppedGBs, Switch[OptionValue[R\[Xi]], 1, True, General, False, _, Abort[]]];
	(*Write particles.py*)
	WritePYParticles[allFields, lag, rulesCompFields, paramsZERO, propInfo];
	(*Write vertices.py*)
	WritePyVertices[vertexData];
	(*Write lorentz.py*)
	WritePyLorentz[lorentzData];
	(*Write couplings.py*)
	WritePyCouplings[couplingData];
	(*Write coupling_orders.py*)
	WritePyCouplingOrders[GetCouplingOrders[]];
	(*Write*)
	ResetDirectory[];
	
	(* remove all internally defined flavor component fields *)
	RemoveField/@ compFields;
	
	dir
]


(* ::Subsubsection::Closed:: *)
(*Auxiliary files*)


ExportUFO::noQCDord= "No \"QCD\" coupling order was defined. This can cause event generators such as MadGraph to crash. Please define the coupling order using, for example: 'DefineCouplingOrder[\"QCD\", gs, 1]' where 'gs' is the QCD gauge coupling of the broken phase, controlling gluon emissions.";


(* ::Text:: *)
(*Function that determines the labels of the QCD and QED gauge groups in the given Lagrangian and sets some global flags with their names:*)


(* ::Text:: *)
(*This function also checks whether a QCD coupling order was defined.*)


IdentifyQCDandQED[lagrangian_]:=Module[{gaugeGroups},
	(* determine QCD group *)
	$QCD$UFO= Nothing; (* reset *)
	(* list of all non-Abelian groups *)
	gaugeGroups= DeleteDuplicates@Cases[lagrangian,Index[_,rep:Except[Lorentz]/;!MemberQ[Keys@$FlavorIndices,rep]]:>Head@rep,All];
	(* check non-Abelian groups *)
	Switch[Length[gaugeGroups]
	,1, (* identify QCD group *)
		If[!MatchQ[$Groups[First@gaugeGroups],SU[3]],
			Message[ExportUFO::noQCD, First[gaugeGroups], $Groups@First[gaugeGroups], SU[3]]; 
			Abort[];
		,
			(* check that the gauge-boson pdg code matches the gluon code *)
			If[MatchQ[$FieldAssociation[$GaugeGroups[First[gaugeGroups]][Field]][UFO,"pdg"], 21],
				$QCD$UFO= First[gaugeGroups]
			,
				Message[ExportUFO::gluonPDG, First[gaugeGroups], 
					$FieldAssociation[$GaugeGroups[First[gaugeGroups]][Field]][UFO,"pdg"], 
					$FieldAssociation[$GaugeGroups[First[gaugeGroups]][Field]] ]; 
				Abort[];
			];
		];
		(* check for QCD coupling order *)
		If[!MemberQ[Keys@$CouplingOrderAssociation, "QCD"],
			Message[ExportUFO::noQCDord]
		]
	,0, (* no QCD group found *)
		Message[ExportUFO::noGroupID, "QCD"];
	,_, (* at most one non-Abelian group can be present *)
		Message[ExportUFO::multGroups, gaugeGroups]; 
		Abort[]
	];
	
	(* determine QED group *)
	$QED$UFO= Nothing; (* reset *)
	Do[
		If[KeyExistsQ[$FieldAssociation[$GaugeGroups[u1Group,Field]],UFO],
			If[KeyExistsQ[$FieldAssociation[$GaugeGroups[u1Group,Field]][UFO],"pdg"],
				If[MatchQ[$FieldAssociation[$GaugeGroups[u1Group,Field]][UFO,"pdg"], 22],
					$QED$UFO= u1Group; Break[]
				]
			]
		]
		,
		{u1Group, Keys@Select[$GaugeGroups,!FreeQ[#,U1]&]}
	];
	If[$QED$UFO === Nothing, Message[ExportUFO::noGroupID, "QED"]];
];


(* ::Text:: *)
(*Check that all PDG codes have been defined:*)


CheckPDGcodes[lagrangian_]:=If[
	!AllTrue[
		DeleteDuplicates@Cases[lagrangian,field:(_Field|_FieldStrength):>First[field],All],
		KeyExistsQ[$FieldAssociation[#,UFO],"pdg"]&
	]
	,
	Message[
		ExportUFO::missingPDG, 
		Cases[
			DeleteDuplicates@Cases[lagrangian,field:(_Field|_FieldStrength):>First[field],All],
			x_/;!KeyExistsQ[$FieldAssociation[x,UFO],"pdg"]
		]
	]; 
	Abort[]
];


(* ::Subsection:: *)
(*Default UFO files*)


(* ::Subsubsection::Closed:: *)
(*File names and directory*)


(* ::Text:: *)
(*List of generic files present in all UFO models:*)


$GenericUFOFiles = {
	"object_library.py",
	"__init__.py",
	"function_library.py",
	"write_param_card.py"
	(*,"propagators.py"*) (* auto generated now *)
};


(* ::Text:: *)
(*Directory of these files*)


$GenericUFOFilesDirectory= FileNameJoin[{$MatchetePath, "Package", "Export", "UFO", "StandardFilesUFO"}];


(* ::Subsubsection::Closed:: *)
(*Copying default files*)


(* ::Text:: *)
(*Copy all default files to given directory*)


CopyDefaultUFOFiles[dir_]:= CopyFile[
	FileNameJoin[{$GenericUFOFilesDirectory, #}],
	FileNameJoin[{dir,#}]
]&/@ $GenericUFOFiles; 


(* ::Subsection:: *)
(*Writing *)


(* ::Subsubsection::Closed:: *)
(*Parameters*)


WritePyParameters[defaultCardParameters_, defaultData_List, UI_] := Module[
	{
		outfile,
		defaultParams,
		params,
		outputString,
		headers,
		ZEROparams
	}
	,
	outfile = "parameters.py";

	defaultParams = DeleteCases[
		Keys[defaultCardParameters],
		x_/;(StringMatchQ[x, "Re_"~~___] || StringMatchQ[x, "Im_"~~___])
	];
	
	(*Get all relevant parameters via UI (need to introduce standard parameters for SM and SMEFT)*)
	headers = {"name", "nature", "type", "value", "texname", "lhablock", "lhacode"};
	params = If[UI,
		(* prepare parameter table w/o UI if concrete ParamCard was given *)
		ParamsTableUI[headers, defaultData, defaultParams],
		(* prepare parameter table w/ UI for modifications if no concrete ParamCard was given and the default is used *)
		ParamsTable[headers, defaultData, defaultParams]
	];
	If[params===$Canceled, Message[ExportUFO::ui]; Abort[]];
	params = Select[params, !AllTrue[#, # == "" &] &]; (*Remove empty sublists*)
	params = DeleteCases[#, "-" | ""] & /@ params; (*Delete all "-" to ensure correct formatting of internal and external degrees of freedom*)
	
	(* determine all parameters that should be entirely removed from the UFO *)
	ZEROparams= Cases[params, {x_,_,_,"ZERO",___}->x, 1];
	params= DeleteCases[params, x_/;MemberQ[ZEROparams, First[x]], 1];
	(* set the dependent parameters to 0 *)
	params= MapAt[
		(*StringReplace[#, (left___~~(Alternatives@@ZEROparams)~~right___) :> (left<>"0"<>right)]&,*)
		StringReplace[#, (WordBoundary~~(Alternatives@@ZEROparams)~~WordBoundary) :> "0"]&,
		params,
		{;;,4}
	];
	
	(* check that all params are properly defined with all required attributes *)
	Do[
		If[Length@param<2, Message[ExportUFO::attr, param[[1]]]; Abort[]];
		Switch[param[[2]],
			"external",
				If[Length[param]<7, Message[ExportUFO::extAttr, param[[1]]]; Abort[]],
			"internal",
				If[Length[param]<5, Message[ExportUFO::intAttr, param[[1]]]; Abort[]]
		]
		,
		{param,params}
	];
	
	(* remove parameters that are contained in the input file but unused *)
	params= RemoveUnusedParams[defaultParams, params, Values@Values@defaultCardParameters];
	(* sort parameters: external first, then internal in an ordering such that they only depend on previously defined params *)
	params= SortParams[params, Values@Values@defaultCardParameters];
	
	outputString = StringJoin@ Table[
		StringJoin[
			paramValues[[1]],
			" = Parameter( \n", 
			StringJoin[
				Table[
					StringJoin[
						"   ",
						headers[[i]],
						Switch[headers[[i]],
							"value",
								If[
									paramValues[[2]] === "external",
									(* check for non-numeric values of external parameters *)
									If[!NumberQ@ ToExpression[StringReplace[paramValues[[i]], "e"->"*^"]], Message[ExportUFO::numext,paramValues[[1]],paramValues[[i]]]; Abort[]];
									" = " <> paramValues[[i]] <> ",\n",
									" = '" <> paramValues[[i]] <> "',\n"
								],
							"lhacode",
								" = [ " <> paramValues[[i]] <> " ],\n",
							_,
								" = '" <> paramValues[[i]] <> "',\n"
						]
					],
					{i, Length[paramValues]}
				],
				")\n"
			]
		],
		{paramValues, params}
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from object_library import all_parameters, Parameter\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n\n"];
	WriteString[outfile, StringJoin[
		"ZERO = Parameter(\n",
		"   name = 'ZERO',\n",
		"   nature = 'internal',\n",
		"   type = 'real',\n",
		"   value = '0.0',\n",
		"   texname = '0'\n",
		")\n"
		]
	];
	WriteString[outfile, outputString];
	Close[outfile];
	
	(* return labels of parameters with value ZERO *)
	ZEROparams
]


(* ::Text:: *)
(*Tries to sort the parameters: external first, then internal parameters in such a way that the value of each parameters only depends on previous params.*)


(*Optimized for large number of default parameters*)
SortParams[params_List, defaultCardParameters_List]:= Module[
	{
		paramsExternal,
		paramsInternal,
		paramsDefaultInternal,
		paramsRest,
		res = {},
		pattern,
		checkLen,
		paramsMissed
	},
	(* external before internal *)
	paramsExternal = Select[params, Part[#,2]==="external"&];
	paramsInternal = Select[params, Part[#,2]==="internal"&];
	(* Special treatment for default complex parameters as their dependence on external parameters is known
		and they potentially account for a large number of the parameters *)
	paramsDefaultInternal = Intersection[paramsInternal, defaultCardParameters];
	
	(* These have to be sorted *)
	paramsRest= Complement[paramsInternal, paramsDefaultInternal];
	(* add params that were potentially overlooked *)
	paramsMissed= Cases[Complement[paramsInternal, paramsRest], x_/; !StringFreeQ[x[[4]], paramsRest[[;;,1]]], 1];
	paramsRest= Join[paramsRest, paramsMissed];
	paramsDefaultInternal= Complement[paramsDefaultInternal, paramsMissed];
	
	(* as long as there are dependent params *)
	While[Length[paramsRest]>0,
		checkLen= Length[paramsRest];
		pattern= paramsRest[[;;,1]];
		
		(* loop over params *)
		Do[
			(* if param is indpenedent of other params add it to result and delete it from params *)
			If[StringFreeQ[paramsRest[[i,4]], Drop[pattern,{i}]], 
				AppendTo[res, paramsRest[[i]]];
				(* mark parameters which should be deleted afterwards *)
				paramsRest[[i]]= nothing
			]
			,
			(* drop the found params only now to have correct part specifications above *)
			{i,Length[paramsRest]}
		];
		
		paramsRest= paramsRest/. nothing->Nothing;
		If[Length[paramsRest]==checkLen, 
			Message[ExportUFO::paramSort];
			res= res~Join~paramsRest;
			Break[]
		]
	];
	
	Join[paramsExternal, paramsDefaultInternal, res]
]



(* ::Text:: *)
(*Removes parameters that are defined in the input file but unused (optimized for large numbers of parameters).*)


RemoveUnusedParams[defaultParams_List, params_List, defaultCardParameters_List] := Module[
	{
		paramsDefault,
		paramsDefaultInternal,
		paramsRestInternal,
		paramsRest,
		paramsReduced = params,
		length
	},

	(* These are needed *)
	paramsDefault = Intersection[paramsReduced, defaultCardParameters];
	paramsDefaultInternal = Select[paramsDefault, Part[#,2]==="internal"&];

	(* These have to be checked *)
	paramsRest = Complement[paramsReduced, paramsDefault];
	paramsRestInternal = Select[paramsRest, Part[#,2]==="internal"&];

	(* Reduce the number of parameters iteratively as the dependencies can be nested *)
	length = Length@paramsRest+1;
	While[Length@paramsRest < length,
		length = Length@paramsRest;
		(* Delete parameters that are not contained in the Lagrangian and not required by other parameters *)
		paramsRest = Select[
			paramsRest
			,
			Function[{param},
				Or[MemberQ[defaultParams, First[param]], Not@AllTrue[Part[paramsRestInternal,;;,4], StringFreeQ[First@param]], Not@AllTrue[Part[paramsDefaultInternal,;;,4], StringFreeQ[First@param]]]
			]
		];
		paramsRestInternal = Select[paramsRest, Part[#,2]==="internal"&];
	];

	paramsReduced = Join[paramsDefault, paramsRest];

	If[Not[Length@paramsReduced === Length@params],
		Echo@ StringForm[ExportUFO::paramRemove, Complement[params, paramsReduced][[;;,1]]];
	];

	paramsReduced
];


(* ::Subsubsection::Closed:: *)
(*Particles*)


WritePYParticles[fields_, lagrangian_, rulesCompFields_, paramsZERO_, propInfo_] := Module[
	{
		outfile         = "particles.py",
		massAssociation = GetFieldsAndMasses[lagrangian, rulesCompFields],
		masses          = FindMasses[lagrangian],
		goldstones      = GetFieldsByProperty[GoldstoneBoson-> True],
		spinAssociation = <|Scalar->"1", Fermion->"2", Vector->"3", Ghost->"-1", AntiGhost->"-1"|>,
		colorAssociation= <|fund->"3", adj->"8"|>,
		listKeys        = {"pdg_code", "name", "antiname", "spin", "color", "mass", "width", "propagator", "texname", "antitexname", "charge", "goldstone"},
		fieldAttributesList,
		outputString,
		lab, labProp,
		nonCanonicalMasses, newParam, newPCEntry, newPCEntries= {},
		\[Lambda]pc,
		stream
	}
	,	
	(* check masses *)
	(* Taylor expand masses *)
	masses= {Normal[Series[First[#]/.c_Coupling:> \[Lambda]pc^Last[c]* c, {\[Lambda]pc, 0, GetMaxEFTOrder[lagrangian, False]-3}]]/. \[Lambda]pc-> 1, Last[#]}&/@ masses;
	(* get the non-canonical mass terms *)
	nonCanonicalMasses= Cases[Values[masses]/. _\[Xi]Gauge-> 1, {Except[_Coupling], l_}:> l, 1];
	
	(* get all relevant data for all fields *)          
	fieldAttributesList = KeyValueMap[
		{
			(* check that PDG codes have been defined *)
			If[KeyExistsQ[#2[UFO], "pdg"], 
				ToString[#2[UFO]["pdg"]], 
				Message[ExportUFO::missingpdg, #1]; Abort[]
			],
			(* check if UFO names have been defined *)
			If[KeyExistsQ[#2[UFO], "name"], 
				ToString @ #2[UFO]["name"], 
				TranslateToUFO[#1]
			],
			(* add conjugate fields *)
			If[#2[SelfConjugate], 
				If[KeyExistsQ[#2[UFO], "name"], 
					ToString @ #2[UFO]["name"], 
					TranslateToUFO[#1]
				], 
				If[KeyExistsQ[#2[UFO], "name"], 
					ToString @ #2[UFO]["name"], 
					TranslateToUFO[#1]
				] <> "~"
			], 
			(* field type *)
			spinAssociation[#2[Type]], 
			(* QCD indices *)
			If[FreeQ[#2[Indices], $QCD$UFO], 
				"1", 
				(* conjugate representations get a minus *)
				ToString[First[SignedQCDRepresentationDimension/@ Select[#2[Indices], !FreeQ[#,$QCD$UFO]&]]]
			], 
			(* mass *)
			If[massAssociation[#1], 
				(* massive particles *)
				If[MemberQ[goldstones, #1] && MemberQ[Values[KeySelect[masses, !MemberQ[goldstones,#]&]][[;;,1]], (First[masses[#1]]/. _\[Xi]Gauge-> 1)],
					(* goldstone with mass of gauge boson *)
					lab= First@ Keys@ Select[
						KeySelect[masses, !MemberQ[goldstones,#]&],
						MatchQ[{(First[masses[#1]]/. _\[Xi]Gauge-> 1), _}]
					];
					labProp= GetFields[lab];
					"Param."<>If[KeyExistsQ[labProp[UFO], "mass"], 
						(* use mass label if given *)
						ToString @ labProp[UFO]["mass"], 
						If[KeyExistsQ[labProp[UFO], "name"],
							(* otherwise use name if given *)
							"M"<>ToString @ labProp[UFO]["name"],
							(* use label if nothing is given *)
							TranslateToUFO[labProp[Mass]]
						]
						(* add '_mod' for non-canonical mass term *)
					] <> If[MemberQ[nonCanonicalMasses,$FieldAssociation[lab,Mass]], "_mod", ""]
				,
					(* any other massive particle *)
					If[!MemberQ[nonCanonicalMasses,$FieldAssociation[#1,Mass]],
						"Param." <> If[KeyExistsQ[masses,#1],
							TranslateToUFO[First@ masses[#1]],
							If[KeyExistsQ[#2[UFO], "mass"], 
								(* use mass label if given *)
								ToString @ #2[UFO]["mass"], 
								If[KeyExistsQ[#2[UFO], "name"],
									(* otherwise use name if given *)
									"M"<>ToString @ #2[UFO]["name"],
									(* use label if nothing is given *)
									TranslateToUFO[#2[Mass]]
								]
							]
						]
					,
						(* get new parameter name *)
						newParam=  If[KeyExistsQ[#2[UFO], "mass"], 
							(* use mass label if given *)
							ToString @ #2[UFO]["mass"], 
							If[KeyExistsQ[#2[UFO], "name"],
								(* otherwise use name if given *)
								"M"<>ToString @ #2[UFO]["name"],
								(* use label if nothing is given *)
								TranslateToUFO[#2[Mass]]
							]
						] <> "_mod";
						(* write parameter card entry for that new param *)
						newPCEntry= newParam <> " = Parameter(
   name = '" <> newParam <> "',
   nature = 'internal',
   type = 'real',
   value = '(" <> TranslateToUFO@First[masses[#1]] <> ").real',
   texname = '" <> newParam <> "'
)";
						AppendTo[newPCEntries, newPCEntry];
						"Param." <> newParam
					]
				]
			, 
				(* massless particles *)
				"Param.ZERO"
			], 
			(* width *)
			If[massAssociation[#1],
				(* massive particles *)
				If[MemberQ[goldstones, #1] && MemberQ[Values[KeySelect[masses, !MemberQ[goldstones,#]&]][[;;,1]], (First[masses[#1]]/. _\[Xi]Gauge-> 1)],
					(* goldstone with mass/width of gauge boson *)
					lab= First@ Keys@ Select[
						KeySelect[masses, !MemberQ[goldstones,#]&],
						MatchQ[{(First[masses[#1]]/. _\[Xi]Gauge-> 1), _}]
					];
					labProp= GetFields[lab];
					"Param."<>If[KeyExistsQ[labProp[UFO], "width"], 
						(* use width label if given *)
						ToString @ labProp[UFO]["width"], 
						If[KeyExistsQ[labProp[UFO], "name"],
							(* otherwise use name if given *)
							"W"<>ToString @ labProp[UFO]["name"],
							(* use label if nothing is given *)
							"W"<>TranslateToUFO[lab]
						]
					]
				,
					(* any other massive particle *)
					"Param."<>If[KeyExistsQ[#2[UFO], "width"], 
						(* use width label if given *)
						ToString @ #2[UFO]["width"], 
						If[KeyExistsQ[#2[UFO], "name"],
							(* otherwise use name if given *)
							"W"<>ToString @ #2[UFO]["name"],
							(* use label if nothing is given *)
							"W"<>TranslateToUFO[#1]
						]
					]
				]
			, 
				(* massless particles *)
				"Param.ZERO"
			], 
			(* propagator *)
			If[KeyExistsQ[propInfo, #1],
				"Prop."<>propInfo[#1],
				"Default"
			],
			(* texname *)
			If[KeyExistsQ[#2[UFO], "name"], 
				ToString @ #2[UFO]["name"], 
				TranslateToUFO[#1]
			],
			(* conjugate texname *)
			If[#2[SelfConjugate] == True, 
				If[KeyExistsQ[#2[UFO], "name"], 
					ToString @ #2[UFO]["name"], 
					TranslateToUFO[#1]
				], 
				If[KeyExistsQ[#2[UFO], "name"], 
					ToString @ #2[UFO]["name"], 
					TranslateToUFO[#1]
				] <> "~"],
			(* electric charges *)
			If[FreeQ[#2[Charges], $QED$UFO], 
				"0", 
				ToString[First[Cases[#2[Charges], $QED$UFO[x_] :> x]], InputForm]
			]
			,
			(* Goldstone *)
			(* Syntax different from the description in (arXiv:2304.09883) *)
			If[#2[GoldstoneBoson],
				"True",
				"False"
			]
			,
			(* selfconjugate *)
			#2[SelfConjugate]
		}&
		, 
		KeyTake[GetFields(*Updated*)[], fields]
	];
	
	(* set specified masses and widths to Param.ZERO *)
	fieldAttributesList= fieldAttributesList/. Map[("Param."<>#->"Param.ZERO")&, paramsZERO];
	
	(* create the output string *)
	outputString = StringJoin@ Table[
		StringJoin[
			fieldAttributes[[2]],
			" = Particle( \n", 
			StringJoin[
				Table[
					If[(listKeys[[i]]==="goldstone" && MatchQ[fieldAttributes[[i]], "'false'"|"False"|"'False'"|False]) || (listKeys[[i]]==="propagator" && MatchQ[fieldAttributes[[i]], "Default"]),
						(* drop goldstone property for non-Goldstones *)
						Nothing
					,
						StringJoin[
							"   ",
							listKeys[[i]],
							If[
								StringContainsQ[listKeys[[i]], "name"],
								" = '" <> fieldAttributes[[i]] <> "',\n",
								" = "  <> fieldAttributes[[i]] <> ",\n"
							]
						]
					],
					{i, Length[fieldAttributes] - 1}
				],
				")\n\n",
				If[
					!(Last @ fieldAttributes),
					StringJoin[
						fieldAttributes[[2]],
						"__tilde__ = ",
						fieldAttributes[[2]],
						".anti()\n\n"
					],
					""
				]
			]
		],
		{fieldAttributes, fieldAttributesList}
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from __future__ import division\n"];
	WriteString[outfile, "from object_library import all_particles, Particle\n"];
	WriteString[outfile, "import parameters as Param\n\n"];
	WriteString[outfile, "import propagators as Prop\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
	
	(* append new parameters to parameters.py *)
	stream = OpenAppend["parameters.py", CharacterEncoding -> "UTF-8"];
	WriteString[stream, # <> "\n"]&/@ newPCEntries;  (* leading \n guards against a missing trailing newline *)
	Close[stream];
]


SignedQCDRepresentationDimension[expr_]:= Module[
	{
		rep,
		isBared,
		representations = GetRepresentations[]
	},
	isBared = Head[expr]===Bar;
	If[isBared, rep = First@expr, rep = expr];

	(*Conjugate representation if second Dynkin coefficient greater than the first one*)
	(-1)^(Boole@isBared + Boole[(Last[#] > First[#])&@ representations[rep][DynkinCoefficients]]) representations[rep][Dimension]
];


(* ::Subsubsection::Closed:: *)
(*Couplings*)


WritePyCouplings[data_] := Module[
	{
	outfile = "couplings.py",
	outputString
	},

	outputString = StringJoin @ Table[
		StringJoin[
			coupling[[2]],
			" = Coupling(\n   name = '",
			coupling[[2]],
			"',\n   value = '",
			coupling[[3]],
			"',\n   order = ",
			coupling[[4]],
			"\n)\n\n"
		],
		{coupling, Transpose[data]}
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from object_library import all_couplings, Coupling\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Vertices*)


WritePyVertices[data_] := Module[
	{
	outfile = "vertices.py",
	headers,
	outputString
	},
	
	outputString = StringJoin @ MapIndexed[
		Module[
			{vertex = #1, i = First@#2},
			
			StringJoin[
				"V_",
				ToString[i],
				" = Vertex(\n   name = 'V_",
				ToString[i],
				"',\n   particles = ",
				"[ ",
				Riffle["P." <> If[KeyExistsQ[$FieldAssociation[# /. (*OverBar|*)Conj[x_] :> x][UFO], "name"], ToString @ $FieldAssociation[# /. (*OverBar|*)Conj[x_] :> x][UFO]["name"], TranslateToUFO[# /. (*OverBar|*)Conj[x_] :> x]] <> If[!MatchQ[#, (*OverBar|*)Conj[___]] && $FieldAssociation[#][SelfConjugate]===False, "__tilde__", ""] & /@ vertex["Fields"], ", "],
				" ],\n   color = [ ",
				Riffle["'" <> ToString@# <> "'" & /@ vertex[ToString@$QCD$UFO], ", "],
				" ],\n   lorentz = [ ",
				Riffle["L." <> StringReplace[SpokenString[#], " " -> ""] & /@ vertex["Lorentz"], ", "],
				" ],\n   couplings = {",
				StringRiffle[
					MapIndexed[
						Function[
							{structure, idx}, 
							"(" <> ToString[structure[[2]]] <> "," <> ToString[idx[[1]] - 1] <> "):C." <> ToString[structure[[1]]]
						],
						Transpose[{vertex["Couplings"], vertex["Colormap"]}]
					],
					", "
				],
				
				"}\n)\n"
			]
		] &,
		data
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from object_library import all_vertices, Vertex\n"];
	WriteString[outfile, "import particles as P\n"];
	WriteString[outfile, "import couplings as C\n"];
	WriteString[outfile, "import lorentz as L\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Lorentz*)


WritePyLorentz[data_] := Module[
	{
	outfile = "lorentz.py",
	outputString
	},
	 
	outputString = StringJoin @ Table[
		StringJoin[
			structure["Name"],
			" = Lorentz(\n   name = '",
			structure["Name"],
			"',\n   spins = [ ",
			StringJoin[Riffle[ToString /@ structure["FieldsSorted"], ", "]],
			" ],\n   structure = '",
			structure["Lorentz"],
			"'\n)\n\n"
		],
		{structure, data}
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from object_library import all_lorentz, Lorentz\n"];
	WriteString[outfile, "from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot\n"];
	WriteString[outfile, "try:\n   import form_factors as ForFac\nexcept ImportError:\n   pass\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
];


(* ::Subsubsection::Closed:: *)
(*Coupling Orders*)


WritePyCouplingOrders[data_] := Module[
	{
	outfile = "coupling_orders.py",
	outputString
	},
	
	outputString = StringJoin @ Table[
		StringJoin[
			ToString@couplingOrder["Name"],
			" = CouplingOrder(\n   name = '",
			ToString@couplingOrder["Name"],
			"',\n   expansion_order = ",
			ToString@couplingOrder["ExpansionOrder"],
			",\n   hierarchy = ",
			ToString@couplingOrder["Hierarchy"],
			(* only for NLO *)
			(*If[
				couplingOrder["PerturbativeExps"] =!= None,
				",\n   perturbative_expansion = " <> ToString[couplingOrder["PerturbativeExps"]],
				""
			],*)
			"\n)\n\n"
		],
		{couplingOrder, data}
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "from object_library import all_orders, CouplingOrder \n\n\n"];
	WriteString[outfile, outputString];
	Close[outfile];
]


(* ::Subsubsection::Closed:: *)
(*Header*)


WriteHeader[outfile_]:=Module[{},
	WriteString[outfile, "# This file was created:\n"];
	WriteString[outfile, "#  - on "<>ToString[DateString[]]<>"\n"];
	WriteString[outfile, "#  - using Matchete v"<>ToString[Global`$MatcheteVersion]<>"\n"];
	WriteString[outfile, "#  - running on Mathematica "<>ToString[$Version]<>".\n\n\n"];
]


(* ::Subsubsection::Closed:: *)
(*Propagators*)


ExportUFO::unknowngauge= "The gauge of the field `1` could not be determined.";
ExportUFO::notgaugefixed= "The given Lagrangian appears to be not gauge fixed.";


WritePYPropagators[lagrangian_, gauges_, droppedGBs_, feynmanGauge_: True] := Module[
	{
		outfile= "propagators.py",
		prop, denominator,
		outputString,
		(* determine gauge fixing of every field *)
		(*gaugeFixing= NonGaugeFixedLagQ[lagrangian, DetailedOutput-> True],*)
		gaugeFixing= gauges,
		(* find all massive vector fields *)
		massiveVectros= MassiveVectorsInLag[lagrangian, Light-> True],
		RxiCounter= 3,
		vecPropAssoc,
		gb
	}
	,
	(* check gauge fixing *)
	If[First@ gaugeFixing,
		Message[ExportUFO::notgaugefixed]; Abort[]
	,
		(* convert to list *)
		gaugeFixing= Normal@ Last[gaugeFixing]
	];
	
	(* predefined denominator *)
	(*
	denominator= "# define only once the denominator since this is always the same
denominator = \"P('mu', id)**2 - Mass(id)**2 + complex(0,1) * Mass(id) * Width(id)\"";
*)
	denominator= "# define only once the denominator since this is always the same
denominator = \"P(-1, id)**2 - Mass(id)**2 + complex(0,1) * Mass(id) * Width(id)\"";
	
	(* the various types  of propagators *)
	prop["Scalar"]= "# propagator for the scalar [currently not in use]
S = Propagator(name = \"S\",
               numerator = \"1\",
               denominator = denominator
               )";
               
	prop["Fermion"]= "# propagator for the incoming fermion [currently not in use]
# the one for the outcomming is computed on the flight 
F = Propagator(name = \"F\",
                numerator = \"(Gamma(-1, 1, 2) * P(-1, id) + Mass(id) * Identity(1, 2))\",
                denominator = denominator
               )";
               
	prop["Unitary"]= "# massive vector in the unitary gauge, cannot be use for massless particles [currently only used when Unitary and Rxi gauges get mixed]
V1 = Propagator(name = \"V1\",
                numerator = \"(-1 * Metric(1, 2) + P(1, id) * P(2, id) / Mass(id)**2)\",
                denominator = denominator
               )";
               
	prop["Feynman"]= "# vector in Feynman gauge [currently not in use]
V2 = Propagator(name = \"V2\",
                numerator = \"-1 * Metric(1, 2)\",
                denominator =  denominator
               )";
               
	prop["Rxi"][xi_, c_]:= "# massless vector in R_xi gauge
V"<>ToString[c]<>" = Propagator(name = \"V"<>ToString[c]<>"\",
                numerator = \"(-Metric(1, 2) + (1 - "<>ToString[xi]<>") * P(1, id) * P(2, id) / (P(-1, id)**2))\",
                denominator = \"P(-1, id)**2\"
               )";
               
	prop["Rxi-massive"][xi_, c_]:= "# massive vector in R_xi gauge
V"<>ToString[c]<>" = Propagator(name = \"V"<>ToString[c]<>"\",
                numerator = \"(-Metric(1, 2) + (1 - "<>ToString[xi]<>") * P(1, id) * P(2, id) / (P(-1, id)**2 - "<>ToString[xi]<>" * (Mass(id)**2 - complex(0,1) * Mass(id) * Width(id))))\",
                denominator = denominator
               )";
	prop["Rxi-GB"][xi_, c_]:= "# Goldstone in R_xi gauge
# WARNING: this is gauge invariant only in the limit of vanishing widths
S"<>ToString[c]<>" = Propagator(name = \"S"<>ToString[c]<>"\",
                numerator = \"1\",
                denominator = \"P(-1, id)**2 - "<>ToString[xi]<>" * Mass(id)**2 + complex(0,1)* "<>ToString[xi]<>" * Mass(id) * Width(id)\"
               )";
	
	(* create the output string containing model independent propagators [sclar, fermion, unitary vector, Feynman vector] *)
	outputString = StringJoin@ {"\n", denominator, "\n\n", prop["Scalar"], "\n\n", prop["Fermion"], "\n\n", prop["Unitary"], "\n\n", prop["Feynman"], "\n\n"};
	
	(* if everything is unitary or everything is Feynman use MadGraph default propagators -> propagators.py file only needed so that import commands can be kept in other files *)
	If[feynmanGauge && Or[
		(ContainsOnly[Last[gauges]/@ massiveVectros, {Unitary}] && ContainsOnly[Values@ KeyDrop[Last@ gauges, massiveVectros], {R\[Xi]}]),
		ContainsOnly[Values@ Last@ gauges, {R\[Xi]}]
		]
	,
		(* only default propagators *)
		vecPropAssoc= <||>
	,
		(* add custom propagators *)
		(* associate propagator types to vector fields and add model-dependent R\[Xi] propagators if required *)
		vecPropAssoc= Reap@ Do[
			Switch[Last@ vec,
				(* R\[Xi] gauges *)
				R\[Xi],
					If[feynmanGauge,
						(* Feynman gauge *)
						(*Sow[First@vec-> "V2"]*)
						Nothing (* use default propagator in Feynman gauge *)
					,
						(* general \[Xi] *)
						Sow[First@vec-> "V"<>ToString[RxiCounter]];
						(* add new propagator for every field *)
						If[MemberQ[massiveVectros, First@ vec], (*distinguish massless and massive vectors *)
							(* massive *)
							outputString= StringJoin@ {outputString, prop["Rxi-massive"][TranslateToUFO[\[Xi]Gauge[First@vec]], RxiCounter], "\n\n"};
							(* if associated GB is in GBs requiring custom propagators *)
							gb= FirstCase[droppedGBs[[1,1]], x_/; !FreeQ[x, \[Xi]Gauge[First@vec], All], None, 1];
							If[gb=!=None,
								Sow[First@gb-> "S"<>ToString[RxiCounter]];
								outputString= StringJoin@ {outputString, prop["Rxi-GB"][TranslateToUFO[\[Xi]Gauge[First@vec]], RxiCounter], "\n\n"};
							]
							
						,
							(* massless *)
							outputString= StringJoin@ {outputString, prop["Rxi"][TranslateToUFO[\[Xi]Gauge[First@vec]], RxiCounter], "\n\n"}
						];
						RxiCounter++
					],
				(* unitary gauge *)
				Unitary,
					Sow[First@vec-> "V1"] (* this is only reached if mixed gauges are present and in this case it is cruial to restrict to unitary propagators *),
				(* no other gauges currently supported *)
				_,
					Message[ExportUFO::unknowngauge, First@vec];
					Abort[]
			];
			
			,
			{vec, gaugeFixing}
		];
		vecPropAssoc= vecPropAssoc[[2,1]];
	];
	
	WriteHeader[outfile];
	WriteString[outfile, "# NOTE: the propagators below are ignored when using regular Unitary or Feynman gauge; Only the Rxi propgators (if present) are used \n\n"];
	WriteString[outfile, "from object_library import all_propagators, Propagator\n\n"];
	(*WriteString[outfile, "import parameters as Param\n\n"];*)
	WriteString[outfile, outputString];
	Close[outfile];
	
	(* return association of vectors and their propagators *)
	Association[vecPropAssoc]
]


(* ::Section:: *)
(*Parameter card*)


(* ::Subsection:: *)
(*Create default parameter card*)


(* ::Subsubsection::Closed:: *)
(*DefaultParamCard*)


(* ::Text:: *)
(*Function providing a default parameter card for use in UFO *)


DefaultParamCard::massbasis= "The mass `1` carries more than one flavor index. Only Lagrangians in the mass eigenbasis are supported, where at most one flavor index should be carried by masses.";
DefaultParamCard::pdg= "The fields `1` have no associated pdg value. This needs to be defined with the fields.";
DefaultParamCard::unknownreim= "Could not determine whether the flavor structure `1` for the coupling `2` is a real or complex. The given flavor structure might not match the coupling symmetries."
DefaultParamCard::gaugefix= "It looks like the provided Lagrangian is not yet gauge fixed. Gauge fixing terms can be included using, for example, the GaugeFixLagrangian routine."
DefaultParamCard::noGPinGBmass= "No gauge fixing parameter could be found in the Goldstonbe boson mass terms."


Options@ DefaultParamCard= {
	DefaultValues   -> False,
	FlavorInvariants-> <||>,
	Gauge           -> Automatic, (* Manual | Unitary | R\[Xi] *)
	R\[Xi]              -> 1, (* General *)
	OutputDirectory -> Automatic
	, 
	"Internal"      -> False
}


(* Do not cash this function! It will crash if flavor symmetries change between different runs *)
DefaultParamCard[L_, opt:OptionsPattern[]] ? OptionsCheck:= Module[
	{
		invalidFields,
		fieldUFO,
		p,
		params,
		paramAssoc,
		fullJSON,
		lag= L,
		ghostLag,
		massiveFields, massTerms,
		(* To use GetFieldsUpdated one would have to run SetCurrentLagrangian first, which leads to some other issues here right now... *)
		fieldAssoc= $FieldAssociation(*GetFieldsUpdated[]*),
		sortingIndex= 0,
		masses,
		counterLHA= 1,
		dir,
		paramNonRed,
		default= OptionValue[DefaultValues],
		flavorCounter,
		currentPDG,
		inputRelationsSM,
		gaugeInfo,
		goldstones= GetFieldsByProperty[GoldstoneBoson-> True],
		droppedGBs,
		requiredNonZeroCouplings, nonZeroLabels
	}
	,
	(* check and (if required) include gauge-fixing terms *)
	{lag, ghostLag, gaugeInfo}= GaugeFixingUFO[lag, OptionValue[Gauge], OptionValue[R\[Xi]], Verbose-> !OptionValue["Internal"]];
	
	(*
	(* throw warning if the Lagrangian is not gauge fixed yet *)
	If[NonGaugeFixedLagQ[lag] && !OptionValue["Internal"],
		Message[DefaultParamCard::gaugefix];
		Abort[]
	];
	*)
	
	(*
	(* check for ghost fields and remove them *)
	If[!FreeQ[lag, Ghost|AntiGhost]  && !OptionValue["Internal"],
		Echo["Dropping ghost terms for tree-level UFO..."];
		(* drop the ghosts since we work at tree level for now *)
		lag= lag/. Field[_,Ghost|AntiGhost, ___]-> 0;
	];
	*)
	
	(* determine all parameters of the Lagrangian *)
	(* this misses all widths and (possibly) some of the masses *)
	params= Union[
		(* get all symbols *)
		Flatten@Cases[
			DeleteCases[Hold@lag, _Coupling|_Field|_FieldStrength|_LCTensor|_CG|_Index, Infinity],
			s_Symbol /; Context[s] === "Global`",
			Infinity
		],
		(* get all couplings *)
		DeleteDuplicates@ Cases[
			(* this line drops mass terms which are treated separately *)
			DeleteCases[lag, x_/; (Count[RemovePower@x, Field[_,_,_,{}], All]==2 && Count[RemovePower@x, _Field, All]==2 && FreeQ[x, _FieldStrength, All]), 1], 
			Coupling[c_, __]:> c, 
			Infinity
		]
	];
	(* determine which params could lead to a divergence if assigned default value 0.0 -> this assumes rational coefficients (no sums in denominators) *)
	requiredNonZeroCouplings= DeleteDuplicates@ Cases[
		Denominator/@ DeleteCases[TermsToList[lag], x_/; Count[RemovePower[x], _Field|_FieldStrength, All]<=2, 1],
		Alternatives@@ DeleteDuplicates[Join[params, DeleteDuplicates@ Cases[lag, Coupling[l_,___]->l, All]]],
		All
	];

	(* determine all massive fields from Lagrangian *)
	massiveFields= Keys@ Select[GetFieldsAndMasses[lag,{}], MatchQ[True]];
	
	(*Check that all fields have a valid pdg-code*)
	invalidFields= Query[Key/@ massiveFields, Key@ UFO, KeyExistsQ["pdg"]]@ fieldAssoc;
	(* Note: ghosts generally require PDG codes, however we ignore them right now at tree level *)
	(*invalidFields= KeySelect[invalidFields, !MatchQ[fieldAssoc[#, Type], Ghost|AntiGhost]&];*)
	invalidFields= Keys@ Select[invalidFields, Not];
	If[Length@ invalidFields > 0,
		Message[DefaultParamCard::pdg, invalidFields];
		Abort[];
	];

(*	
	(* get the mass labels *)
	masses= If[fieldAssoc[#,Mass]===0,
		(* get new unique mass label *)
		Module[{newLabel= ToString["M"]<>If[KeyExistsQ[GetFields[#, UFO], "name"], GetFields[#, UFO]["name"], ToString[#]]},
			newLabel= If[Defined[ToExpression@ newLabel],
				Unique[newLabel]
				,
				ToExpression@ newLabel
			];
			Echo["The field "<>ToString[#]<>" was defined as massless, but a mass term was found in the Lagrangian. Created new label \""<>ToString[newLabel]<>"\" for that mass."];
			newLabel
		]
		,
		fieldAssoc[#,Mass]
	]&/@ massiveFields;
*)

	(* find massive fields *)
	massTerms= FindMasses[lag];
	(* delete goldstone masses if expressed through gauge masses *)
	droppedGBs= Keys[massTerms];
	massTerms= KeySelect[
		massTerms, 
		Or[
			!MemberQ[goldstones, #],
			!MemberQ[
				Values[KeySelect[massTerms, !MemberQ[goldstones,#]&]][[;;,1]], 
				First[massTerms[#]]/. _\[Xi]Gauge-> 1
			]
		]&
	];
	(* GB labels dropped from param card *)
	droppedGBs= Complement[droppedGBs, Keys[massTerms]];
	(* Sow the dropped GBs with their gauge parameters in Internal mode -> needed for R\[Xi] custom propagators *)
	If[OptionValue["Internal"] && OptionValue[R\[Xi]]===General, 
		Sow[{#, FirstCase[
			First[FindMasses[lag][#]], 
			_\[Xi]Gauge, 
			Message[DefaultParamCard::noGPinGBmass]; Abort[], 
			All
		]}&/@ droppedGBs];
	];
	(* get only labels *)
	masses= Normal[massTerms][[;;,2,2]];
	
	(* prepend gauge parameters and masses to remaining parameters  *)
	params= DeleteDuplicates@ Join[
		DeleteDuplicates@ Cases[lag, _\[Xi]Gauge, All],
		masses,
		params		
	];
	
	(* determine SM input relations *)
	If[default,
		inputRelationsSM= DetermineSMInputReplacements[lag]
	];
	
	(*Create template association covering all the parameters*)
	counterLHA= 1000;
	{paramAssoc, nonZeroLabels}= Reap@ Table[
			paramNonRed= NonRedundantUFOParams[param, OptionValue@ FlavorInvariants];
			flavorCounter= 0;
			{Table[
				flavorCounter++; (* increase counter for flavor multiplets to find the right PDG code *)
				If[MemberQ[masses, param], 
					(* if this is a mass, get the PDG code of the particle *)
					currentPDG= ($FieldAssociation[First@ Keys@ Select[massTerms, !FreeQ[#[[-1]], param]&,1]][UFO,"pdg"])/. l_List:> l[[flavorCounter]],
					currentPDG= None
				];
				First@ p -> <|
						"name"    -> p[[1]],
						"nature"  -> If[MemberQ[masses, param] && !MatchQ[First[Values@ Select[massTerms, !FreeQ[#[[-1]], param]&, 1], False], {Coupling[l_,___],l_}],
							(* if masses have not been replaced -> include their definitions *)
							p[[3, 2]]= "internal"
							,
							If[default && MemberQ[masses, param] && 24===currentPDG,
								p[[3, 2]]= "internal" (* make the W-mass internal when using the defaults *)
								,
								If[default && MemberQ[inputRelationsSM[[;;,1]], param],
									p[[3, 2]]= "internal" (* make the parameters assigned to SM inputs internal when using the defaults *)
									,
									p[[3, 2]]
								]
							]
						],
						"type"    -> p[[3, 1]],
						"value"   -> If[!default && MemberQ[masses, param] && !MatchQ[First[Values@ Select[massTerms, !FreeQ[#[[-1]], param]&, 1], False], {Coupling[l_,___],l_}],
							(* if masses have not been replaced -> include their definitions *)
							(* need to check flavor indices *)
							If[Length[paramNonRed] >1,
								TranslateToUFO[(First@ First@ Values@ Select[massTerms, !FreeQ[#[[-1]], param]&, 1])/. Index[Global`d$$1,_]-> (First@ First@ Position[paramNonRed[[;;,1]], p[[1]]])]
							,
								(* flavorless field *)
								TranslateToUFO@ First@ First@ Values@ Select[massTerms, !FreeQ[#[[-1]], param]&, 1]
							]
							,
							If[default && MemberQ[masses, param],
								(* for default SM masses *)
								If[MemberQ[($FieldAssociation[#,UFO,"pdg"]&/@ GetFieldsByProperty[GoldstoneBoson->True]), currentPDG],
									(* goldstones *)
									TranslateToUFO@ First@ First@ Values@ Select[massTerms, !FreeQ[#[[-1]], param]&, 1]
									,
									(* rest *)
									DefaultMassesSM[currentPDG]
								]
								,
								If[default && MemberQ[inputRelationsSM[[;;,1]], param],
									param/. inputRelationsSM
								,
									(* change default value away from zero if it leads to a divergence *)
									If[MemberQ[requiredNonZeroCouplings, param] && p[[3,2]]==="external",
										p[[2]]/. ("0"| "0.0"| 0| 0.):> (If[!TrueQ[OptionValue["Internal"]], Sow[p[[1]]]]; "1")
									,
										p[[2]]	
									]
								]
							]
						],
						"texname" -> p[[1]],
						If[MatchQ[p[[3, 2]],"external"],
							(*"lhablock" -> If[MemberQ[masses, param], "MASS", "COUPLINGS"]*)
							"lhablock" -> Switch[param,
								Alternatives@@ masses, "MASS",
								_\[Xi]Gauge, "GAUGE",
								_, "COUPLINGS"]
							,
							Nothing
						],
						If[MatchQ[p[[3, 2]],"external"],
							"lhacode" -> If[MemberQ[masses, param], 
								(* for masses pick the PDG code *)
								{currentPDG}, 
								counterLHA++
							]
							,
							Nothing
						]
					|>
			, {p, paramNonRed}], sortingIndex++}
		, {param, params}];
	
	If[Length[First[nonZeroLabels, {}]]>0,
		Echo["The parameters '" <> ToString[First@ nonZeroLabels] <> "' were assigned default value '1' since setting them to '0' potentially leads to divergent coefficients."];
	];
	
	(*... and widths*)
	paramAssoc= Join[paramAssoc, Table[
			{WidthUFOParams[p, DefaultValues-> default], sortingIndex++}
		, {p, Complement[massiveFields, droppedGBs]}] ];
	
	(*Put the sorting index by each parameter*)
	paramAssoc= Flatten[Thread/@ paramAssoc, 1];

	(*Put all internal parameters to the end, and subsequently sort by sorting index*)
	paramAssoc= SortBy[paramAssoc, {If[#[[1, 2]]["nature"] === "internal", 1, 0]&, Last}]; 
	paramAssoc= Association@@ paramAssoc[[;;, 1]];
	
	(* include default SM inputs if required *)
	If[default,
		paramAssoc= PrependTo[paramAssoc, DefaultSMinputs]
	];

	(* full template for export *)
	fullJSON = <|
		"description" -> "Minimal JSON with all parameters", 
		"modelname"   -> "Your model name", 
		"info"        -> <|
							"creation date"       -> ToString@DateString[],
							"Matchete version"    -> ToString@Global`$MatcheteVersion,
							"Mathematica version" -> ToString@$Version
						|>,
		"parameters"  -> paramAssoc
	|>;
	
	(* for internal use return only the association *)
	If[OptionValue["Internal"], Return[fullJSON]];
	
	(* determine file path and potentially check if it should be overwritten *)
	dir= If[OptionValue[OutputDirectory]===Automatic, 
		Directory[], 
		OptionValue[OutputDirectory]
	];
	dir= FileNameJoin@{dir, "model_parameters.json"};
	(* open interactive dialog if frontend is available *)
	If[And[FileExistsQ[dir], $FrontEnd =!= Null],
		If[
			!DialogInput@ Column@ {
				Style["The file \"" <> ToString[dir] <> "\" already exists.", Bold, 14],
				Spacer[10],
				Style["Do you want to overwrite it?", Bold, 14],
				Spacer[10],
				Row@ {
					DefaultButton["Yes", DialogReturn[True]],
					Spacer[20],
					CancelButton["No", DialogReturn[False]]
				}
			}
			,
			Abort[]
		]
	];
	 
	(* export as JSON file *)
	Export[
		dir,
		fullJSON,
		"JSON",
		"Indentation" -> 2
	]
]


(* ::Text:: *)
(*Default parameters template for the decay widths associated with a field*)


Options@ WidthUFOParams= {DefaultValues-> False}


WidthUFOParams[field_Symbol, OptionsPattern[]]:= Module[{fieldUFO, len, n, pdgs, ufoNames, default= OptionValue[DefaultValues]},
	fieldUFO= GetFields[field, UFO]; (* don't use GetFieldsUpdated here! *)
	pdgs= If[Head@ # === List, #, {#}]&@ fieldUFO["pdg"];
	
	(*Checks which if any information can be used to provide default UFO names*)
	ufoNames= Which[
		KeyExistsQ[fieldUFO, "width"],
			EnforceList@ fieldUFO["width"],
		KeyExistsQ[fieldUFO, "name"],
			"W"<> #&/@ EnforceList@ fieldUFO["name"],
		(len= Length@ pdgs) > 1, 
			"W"<> TranslateToUFO@ field<> "_"<> ToString@ #&/@ Range@ len,
		True,
			{"W"<> TranslateToUFO@ field}
		];
	
	(*Loops over each flavor*)
	Table[
		ufoNames[[n]] -> <|
			"name"     -> ufoNames[[n]],
			"nature"   -> "external", (* widths are not internal *)
			"type"     -> "real",
			"value"    -> If[default, DefaultWidthsSM[pdgs[[n]]], 0.0],
			"texname"  -> ufoNames[[n]],
			"lhablock" -> "DECAY",
			"lhacode"  -> {pdgs[[n]]}
		|>
	, {n, Length@ pdgs}]
]
EnforceList@ l_List:= l;
EnforceList@ x_:= {x};


(* ::Subsubsection::Closed:: *)
(*Default SM inputs*)


(* ::Text:: *)
(*This is a set of default SM input parameters (Gf, aEM, aS) for the 5-flavor scheme and the masses (mZ, mh, mt, mW) as well as all widths are below*)


DefaultSMinputs= <|
	"aEM" -> <|
		"name"        -> "aEM",
		"nature"      -> "external",
		"type"        -> "real",
		"value"       -> 1./127.930, (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.15] *)
		"texname"     -> "\\alpha_{\\text{EM}}",
		"lhablock"    -> "SMINPUTS",
		"lhacode"     -> {1},
		"description" -> "Electromagnetic coupling (\[Alpha]EM(MZ)=e^2/4\[Pi])"
	|>,
	"Gf" -> <|
		"name"        -> "Gf",
		"nature"      -> "external",
		"type"        -> "real",
		"value"       -> 1.1663785*10^-5, (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.17] *)
		"texname"     -> "G_f",
		"lhablock"    -> "SMINPUTS",
		"lhacode"     -> {2},
		"description" -> "Fermi's constant"
	|>,
	"aS" -> <|
		"name"        -> "aS",
		"nature"      -> "external",
		"type"        -> "real",
		"value"       -> 0.1180, (* https://pdg.lbl.gov/2025/reviews/contents_sports.html *)
		"texname"     -> "\\alpha_s",
		"lhablock"    -> "SMINPUTS",
		"lhacode"     -> {3},
		"description" -> "QCD coupling (aS = gs^2/4\[Pi])"
	|>
|>;


(* ::Text:: *)
(*This is a list with default values for the masses of SM particles. [5-flavor scheme]*)


UFONames[$FieldAssociation[FieldFromPDG[23],Mass]]


DefaultMassesSM[pdg_]:= Switch[pdg,
	1|2|3|4|5|11|13|15 , "ZERO",  (* d,u,s,c,b,e,\[Mu],\[Tau] *)
	6                  , 172.56,  (* t *) (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.21] *)
	23                 , 91.1880, (* Z *) (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.19] *)
	24                 , With[{mZ= UFONames[$FieldAssociation[FieldFromPDG[23],Mass]]}, (* W *)
							"cmath.sqrt("<>mZ<>"**2/2. * (1 + cmath.sqrt(1 - (4.*cmath.pi*aEM) * (1./(cmath.sqrt(2.)*Gf)) / "<>mZ<>"**2)))"
						],    
	25                 , 125.08,  (* Higgs *) (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.20] *)
	_                  , 0.0      (* unknown PDG *)
] 


(* ::Text:: *)
(*This is a list with default values for the widths of SM particles. [5-flavor scheme]*)


DefaultWidthsSM[pdg_]:= Switch[pdg,
	1|2|3|4|5|11|13|15 , "ZERO", (* d,u,s,c,b,e,\[Mu],\[Tau] *)
	6                  , 1.42,   (* t *)     (* PDG average *)
	23                 , 2.4955, (* Z *)     (* https://pdg.lbl.gov/2025/reviews/contents_sports.html [p.19] *)
	24                 , 2.14,   (* W *)     (* PDG average *)
	25                 , 0.0037, (* Higgs *) (* PDG average *)
	_                  , 0.0     (* non-SM particle *)
] 


(* ::Text:: *)
(*This function determines how to assign the SM input values to the parameters of the model*)


DefaultParamCard::SMinput= "Matchete encountered issues associating the default SM input values to this model. Please make sure that the resulting parameter card is actually consistent.";


DetermineSMInputReplacements[lag_]:= Module[
	{
		(* get field labels from their PDG codes *)
		higgs, gluon, photon, WBosons,
		QCD, QED, higgsTrilinear, higgsQuartic, WWhh,
		L= lag/. hbar-> 0,
		sol, res,
		higgsMassLabel, WMassLabel
	}
	,
Check[	
	(* get field labels *)
	higgs   = FieldFromPDG[25];
	gluon   = FieldFromPDG[21];
	photon  = FieldFromPDG[22];
	WBosons = FieldFromPDG[24];
	(* get W and Higgs UFO mass labels *)
	higgsMassLabel = UFONames[$FieldAssociation[higgs,Mass]];
	WMassLabel     = UFONames[$FieldAssociation[WBosons,Mass]];
	
	(* get Higgs trilinear and quartic interction terms *)
	higgsTrilinear = SelectOperatorClass[L, {higgs,higgs,higgs}, 0]/. _Field-> 1;
	higgsQuartic   = SelectOperatorClass[L, {higgs,higgs,higgs,higgs}, 0]/. _Field-> 1;
	(* get gluon and photon kinetic terms *)
	QCD = PowerExpand[(-4*SelectOperatorClass[L, {gluon,gluon}, 0]/. _FieldStrength-> 1)^(-1/2)];
	QED = PowerExpand[(-4*SelectOperatorClass[L, {photon,photon}, 0]/. _FieldStrength-> 1)^(-1/2)];
	(*  get WWhh interaction *)
	WWhh = SelectOperatorClass[L, {WBosons,Bar@WBosons,higgs,higgs}, 0]/. _Field-> 1;
	
	(* find couplings/variables to solve for *)
	sol = DeleteDuplicates@ Cases[{QCD,QED,higgsTrilinear,higgsQuartic,WWhh}, _Coupling, All];
	
	(* solve teh system of equations *)
	sol = Last@ Solve[{
		-(1/2)"vev"*"lam"==higgsTrilinear,
		-(1/8)"lam"==higgsQuartic,
		"gQED"==QED,
		"gQCD"==QCD,
		1/4 ("gQED")^2/("sw")^2==WWhh
	}, sol];
	
	(* determine replacement rule for internal variables with input params *)
	res= sol//. {
		"gQED" -> Sqrt["4*\[Pi]""aEM"],
		"gQCD" -> Sqrt["4*\[Pi]""aS"],
		"vev"  -> 1/Sqrt[Sqrt["2"]"Gf"],
		"lam"  -> higgsMassLabel^2/("vev")^2,
		"sw"   -> ("gQED"*"vev")/(2*WMassLabel)
	};
	
	(* add cw - this assumes cw was defined through DefineCoSinePair *)
	AppendTo[res, 
		First@Flatten[
			Cases[$trigPairs, x_/; !FreeQ[x, #]:> (FirstCase[x, y_/;(y=!=#), 1]->Sqrt[1-UFONames[First@#]^2]),1]&/@ res[[;;,1]]
		]
	];
	
	(* print replacements *)
	Echo["The following SM input relations have been determined:"]
	Print[Format[res, NiceForm]];
	
	(* some final string manipulations *)
	res = res/. HoldPattern@Rule[a_,b_]:> Rule[a, ToString[b, FortranForm]];
	res = res/. str_String:> StringReplace[str, {"\""->"", "\[Pi]"->"cmath.pi", "Sqrt"->"cmath.sqrt"}];
	res = res/. Coupling[l_,__]->l
,
	Message[DefaultParamCard::SMinput];
	res
]]


(* ::Subsection:: *)
(*Reading input parameter card*)


(* ::Subsubsection::Closed:: *)
(*Reading input parameter card w/o UI for custom param card*)


(* ::Text:: *)
(*Re-formats the data from the input ParamCard without UI (used if a custom ParamCard is used)*)


ParamsTable[headers_List, defaultData_List : {"parameters" -> {}}, defaultParamsSymbols_List : {}] := Module[
    {
	    defaultParamsData,
	    table,
	    paramsWOData,
	    optionsNature,
	    optionsType,
	    addRow
    }
    ,
    (* Extract data out of defaultData*)
    defaultParamsData = Lookup[defaultData, "parameters"];
    (*
    table= Table[
		ToString@ Lookup[Lookup[defaultParamsData, paramName], header, ""],
		{paramName, Keys[defaultParamsData]},
		{header, headers}
    ];
    *)
    (* performance improved version of code above *)
    defaultParamsData= defaultParamsData/. List->Association;
	table= Values@ Map[
		ToString/@ Lookup[#, headers, ""]&, 
		defaultParamsData
	];
	
    (* Check if all MATCHETE parameters get a value from the default data if applicable *)
    paramsWOData = Complement[defaultParamsSymbols, Keys @ defaultParamsData];
    
    If[
		defaultParamsData =!= {} && paramsWOData =!= {},
		Message[ExportUFO::missingdata, paramsWOData];
		Abort[]
    ];
	
	(* Initialize table from parameters out of the lagrangian *)
    table = Join[
		table,
		PadRight[List@ #, Length[headers], ""]& /@ paramsWOData
	];
	
	table
]


(* ::Subsubsection::Closed:: *)
(*Reading input parameter card w/ UI for default param card*)


(* ::Text:: *)
(*Re-formats the data from the input ParamCard with UI (called if no InputFile is given and the DefaultParamCard is used instead)*)


ParamsTableUI[headers_List, defaultData_List : {"parameters" -> {}}, defaultParamsSymbols_List : {}] := DialogInput[DynamicModule[
    {
    defaultParamsData,
    table,
    paramsWOData,
    optionsNature,
    optionsType,
    addRow
    },
    
    (* Extract data out of defaultData*)
    defaultParamsData = Lookup[defaultData, "parameters"]/. ("lhacode" -> {x_}) :> ("lhacode" -> x);
    (*
    table= Table[
		ToString@ Lookup[Lookup[defaultParamsData, paramName], header, ""],
		{paramName, Keys[defaultParamsData]},
		{header, headers}
    ];
    *)
    (* performance improved version of code above *)
    defaultParamsData= defaultParamsData/. List->Association;
	table= Values@ Map[
		ToString/@ Lookup[#, headers, ""]&, 
		defaultParamsData
	];
	
    (* Check if all MATCHETE parameters get a value from the default data if applicable *)
    paramsWOData = Complement[defaultParamsSymbols, Keys @ defaultParamsData];
    
    If[
		defaultParamsData =!= {} && paramsWOData =!= {},
		Message[ExportUFO::missingdata, paramsWOData];
		Abort[]
    ];
	
	(* Initialize table from parameters out of the lagrangian *)
    table = Join[
		table,
		PadRight[List @ #, Length[headers], ""] & /@ paramsWOData
	];
	
    (* Dropdown options *)
    optionsNature = {"external", "internal"};
    optionsType = {"real", "complex"};

    (* Add Row *)
    addRow[] := (table = Append[table, ConstantArray["", Length[headers]]]);

    (* UI *)
    Column[
        {Pane[Dynamic[Grid[
            Prepend[
                Table[
                    With[
                        {ii = i, jj = j},
                        DynamicModule[{},
                            Switch[
                                jj,
                                2, (* Nature Column *)
                                PopupMenu[
                                    Dynamic[
                                        table[[ii, jj]],
                                        (table[[ii, jj]] = #;
                                         Which[
											# === "external", 
											table[[ii, 3]] = "real",
											# === "internal",
											table[[ii, 6]] = "-";
											table[[ii, 7]] = "-"]) & (* Fix type to real *)
                                    ],
                                    optionsNature
                                ],
                                3, (* Type Column *)
                                PopupMenu[
                                    Dynamic[
                                        table[[ii, jj]],
                                        (If[table[[ii, 2]] === "external", table[[ii, jj]] = "real", table[[ii, jj]] = #]) &
                                    ],
                                    optionsType,
                                    Enabled -> (table[[ii, 2]] =!= "external") (* Disable type selection if external *)
                                ],
                                6 | 7, (* lhablock and lhacode Columns *)
                                InputField[
                                    Dynamic[table[[ii, jj]]],
                                    String,
                                    Enabled -> (table[[ii, 2]] =!= "internal") (* Disable when internal *)
                                ],
                                _, (* Default case *)
                                InputField[Dynamic[table[[ii, jj]]], String]
                            ]
                        ]
                    ],
                    {i, Length[table]},
                    {j, Length[headers]}
                ],
                headers
            ],
            Frame -> All
        ]],
		{Automatic, 300},
		Scrollbars -> True
		],
        Row[
            {
            Button["Add Row", addRow[]],
            Spacer[10],
            DefaultButton["Done", DialogReturn[table]]
            }
        ]}
    ]
]];


(* ::Section:: *)
(*Coupling Orders*)


(* ::Subsection:: *)
(*Helper functions*)


(* ::Subsubsection::Closed:: *)
(*Extracting powers of couplings*)


ExtractCouplingPower[expr_, couplingName_] := Module[
	{
	dummy,
	expanded
	},
	
	dummy = Unique["dummy"];
	expanded = BetterExpand[expr /. Coupling[couplingName, ___] -> dummy];
	Max[Exponent[#, dummy] & @ expanded]
];


(* ::Subsubsection::Closed:: *)
(*Extracting labels of couplings*)


(* ::Text:: *)
(*Function that extracts the labels of a coupling from some expression*)


ExtractCouplingLabel[expr_] := Module[
	{head},
	head = Head[expr];
	Which[
		head === Symbol,
		expr,
		head === Power || head === Times,
		ExtractCouplingLabel[First[expr]],
		head === Coupling,
		First[expr],
		True,
		Missing["UnknownForm",expr]
	]
]


(* ::Subsubsection::Closed:: *)
(*Powers*)


(* ::Text:: *)
(*Counts the powers of symbol sym in the given expression expr.*)


PowerOf[expr_,sym_] := Module[
	{
		factors,
		exp
	}
	,
	(*1) Turn a product into a list of factors*)
	factors = If[Head[expr] === Times, List @@ expr, {expr}];
	(*2) Sum up all explicit Powers of sym*)
	exp = Total @ Cases[factors, Power[sym, n_] :> n];
	(*3) Count any bare sym\[CloseCurlyQuote]s (each is+1)*)
	exp += Count[factors, sym];
	exp
]


(* ::Subsubsection::Closed:: *)
(*ClassifyCouplingOrder*)


ClassifyCouplingOrder[term_]:= Module[
	{
		couplingOrders= $CouplingOrderAssociation,
		(*dummy, couplingOrderAssoc,*)
		relevantCouplingOrders
	}
	,
	relevantCouplingOrders= Select[Keys@ couplingOrders, (!FreeQ[term, Alternatives@@ Keys[couplingOrders[#]["Couplings"]], All])&];
	
	AssociationMap[
		Function[
			order,
			(*
			couplingOrderAssoc = Power[dummy, #]& /@ KeyMap[Coupling[#, __]&, couplingOrders[order]["Couplings"]];
			Exponent[term/. Normal[couplingOrderAssoc], dummy, Min]
			*)
			(* not having to re-build couplingOrderAssoc for every term saves significant time *)
			Exponent[term/. $couplingOrderRules[order], "dummy_CO_"<>order, Min]
		]
		,
		relevantCouplingOrders
	]
]


(* ::Text:: *)
(*This is global association of rules for counting the coupling order of terms. Not having to reconstruct this rule for every term improves performance significantly.*)


$couplingOrderRules= <||>;


SetCouplingOrderRules[]:= $couplingOrderRules= Association@@ Table[
	coupOrd-> Normal[Power["dummy_CO_"<>coupOrd, #]&/@ KeyMap[Coupling[#, __]&, $CouplingOrderAssociation[coupOrd]["Couplings"]]]
,
	{coupOrd, Keys@ $CouplingOrderAssociation}
]


(* ::Subsection:: *)
(*Main functions*)


(* ::Subsubsection::Closed:: *)
(*Association*)


(* ::Text:: *)
(*Association where all information on coupling orders are saved.*)


$CouplingOrderAssociation = <||>;


(* ::Subsubsection::Closed:: *)
(*GetCouplingOrders*)


GetCouplingOrders[names___]:= Return @ $CouplingOrderAssociation[names];


(* ::Subsubsection::Closed:: *)
(*DefineCouplingOrder*)


DefineCouplingOrder::name      = "The label `1` specified for the coupling order must be a string that is not yet used for any other coupling order.";
DefineCouplingOrder::couplings = "The couplings `1` have not been defined before.";
DefineCouplingOrder::format    = "Inconsistent format provided for the coupling `1`. Please provide a symbol, a power of a symbol, or a list of a symbol paired with an integer. Only integer orders are allowed, so when using powers only inverse integer exponents are allowed.";


(* only required for NLO *)
(*Options[DefineCouplingOrder] = {"PerturbativeExps" -> None};*)


DefineCouplingOrder[name_, couplings:Except[_List|_Symbol|Power[_Symbol, _?IntegerQ]], hierarchy_, expansionOrder_:99]:= (Message[DefineCouplingOrder::format,couplings]; Abort[])


DefineCouplingOrder[name_, couplings:{_,_?IntegerQ}, hierarchy_, expansionOrder_:99]:= DefineCouplingOrder[name, {couplings}, hierarchy, expansionOrder]


DefineCouplingOrder[name_, Power[coupling_, n_?IntegerQ], hierarchy_, expansionOrder_:99]:= DefineCouplingOrder[name, {{coupling, 1/n}}, hierarchy, expansionOrder]


DefineCouplingOrder[name_, coupling_Symbol, hierarchy_, expansionOrder_:99]:= DefineCouplingOrder[name, {{coupling, 1}}, hierarchy, expansionOrder]


DefineCouplingOrder[name_, couplings:Except[{{_,_?NumberQ}..}, _List], hierarchy_, expansionOrder_:99]:= DefineCouplingOrder[
	name, 
	Table[
		Switch[coup,
			_Symbol,
				{coup, 1},
			Power[_Symbol, n_/; IntegerQ[1/n]],
				{First@coup, 1/Last@coup},
			{_Symbol, _?IntegerQ},
				coup,
			_,
				Message[DefineCouplingOrder::format, coup]; Abort[]
		]
	,
		{coup, couplings}
	],
	hierarchy,
	expansionOrder
]


DefineCouplingOrder[name_, couplings:{{_,_?NumberQ}..}, hierarchy_, expansionOrder_:99]:= 
DefineCouplingOrder[name, couplings, hierarchy, expansionOrder]= Module[
	{
		(* get labels of all couplings *)
		(*
		couplingNames= ExtractCouplingLabel[#]&/@ If[!ListQ[couplings], List@ couplings, couplings],
		*)
		couplingNames= couplings[[;;,1]],
		couplingCount
	}
	,
	(* check that coupling order symbol is not yet defiend *)
	If[(!StringQ[name]) || MemberQ[Keys@$CouplingOrderAssociation, name],
		Message[DefineCouplingOrder::name, name];
		Abort[]
	];
	
	(* make couplingNames to list *)
	If[!ListQ[couplingNames],
		couplingNames= List@ couplingNames
	];
	
	(* check that all couplings have been defined *)
	If[
		!AllTrue[couplingNames, KeyExistsQ[GetCouplings[], #] & ],
		Message[DefineCouplingOrder::couplings, StringRiffle[ToString/@ Cases[couplingNames, x_/; !KeyExistsQ[GetCouplings[],x], 1], ", "]];
		Abort[]
	];
	
	(*
	couplingCount= AssociationThread[couplingNames, MapThread[PowerOf, {If[!ListQ[couplings], List@ couplings, couplings]//. Coupling[label_, ___]:> label, couplingNames}]];
	*)
	couplingCount= Association[Rule@@@ couplings];
	
	(* add to global coupling orders *)
	AppendTo[
		$CouplingOrderAssociation,
		name-> <|
			"Name"             -> name,
			"Couplings"        -> couplingCount,
			"ExpansionOrder"   -> expansionOrder,
			"Hierarchy"        -> hierarchy
			(* only for NLO *)
			(*, PerturbativeExps" -> OptionValue@ "PerturbativeExps"*)
		|>
	];
	
	(* rebuild some association *)
	SetCouplingOrderRules[];
	
	(*name::usage= "Coupling order " <> ToString@ name <> " for the couplings " <> StringRiffle[ToString/@ couplingNames, ", "];*)
];


(* ::Subsubsection::Closed:: *)
(*RemoveCouplingOrder*)


RemoveCouplingOrder::noco= "The coupling order \"`1`\" does not match any of the defined coupling orders. The currently defined coupling orders are `2`."


RemoveCouplingOrder[name_]:= If[KeyExistsQ[$CouplingOrderAssociation, name],
	DownValues[DefineCouplingOrder]= DeleteCases[DownValues[DefineCouplingOrder], x_/;!FreeQ[x, name, All]];
	KeyDropFrom[$CouplingOrderAssociation, name];
	
	ClearAll @ name;
	(* rebuild rules for counting coupling orders *)
	SetCouplingOrderRules[];
,
	Message[RemoveCouplingOrder::noco, name, Short[Keys@ $CouplingOrderAssociation, 10]]
];


(* ::Subsubsection::Closed:: *)
(*ResetCouplingOrders*)


ResetCouplingOrders[] := (RemoveCouplingOrder/@ Keys[$CouplingOrderAssociation];);


(* ::Section:: *)
(*Help Functions*)


(* ::Subsection:: *)
(*Miscellaneous*)


(* ::Subsubsection::Closed:: *)
(*GetFieldFromPDG*)


(* ::Text:: *)
(*This function gives the field associated to a PDG code:*)


FieldFromPDG[n_Integer]:= First[
	Keys@ KeySelect[$FieldAssociation, ($FieldAssociation[#,UFO,"pdg"]===n)&], 
	None
]


(* ::Subsubsection::Closed:: *)
(*GetFieldsAndMasses*)


(* ::Text:: *)
(*Check which fields of the Lagrangian have mass terms*)


GetFieldsAndMasses[lag_, rulesCompFields_] := Module[
	{
		terms,
		(* rules to split flavor mutiplets into components *)
		rules= DeleteDuplicates[rulesCompFields/. (*OverBar|*)Conj-> Identity/. Ext[lab_,_]->lab]
	}
	,
	(* extract mass terms and split them by flavor indices *)
	terms= FlavorExpandedLagrangian@ IsolateMassTerms[lag];
	terms= Collect[RelabelIndices@ terms, _NCM| Power[_Field, _]| _Field];
	
	(* use component fields *)
	terms= terms/. rules;
	
	(* determine massive fields *)
	terms= DeleteDuplicates@ Cases[terms,Field[l_,___]:>l,All];
	
	(* get rule for labels only *)
	rules= (#[[1,1]]->#[[2,1]])&/@ rules;
	
	(* association of massive and massles fields *)
	Association@ Table[
		f-> MemberQ[terms, f]
		,
		{
			f, 
			Complement[
				(* get list of all component labels *)
				DeleteDuplicates@ Flatten@ Map[
					(#/.(List/@rules))&,
					(* get all field labels on flavor-multiplet level *)
					DeleteDuplicates@Cases[lag, x:(_Field|_FieldStrength):>First[x], All]
				]
				,
				rules[[;;,1]] (* this drops the multiplet labels *)
			]
		}
	]
]


(* ::Subsubsection::Closed:: *)
(*FindMasses*)


(* ::Text:: *)
(*This function takes a Lagrangian as input and returns an association where the keys are the labels of the massive fields and the values are lists where the first entry is the expression for their mass couplings in the Lagrangian and the second entry is the canonical mass label assigned.*)


DefaultParamCard::massmix= "Mass mixing detected in the term `1`. Please provide the Lagrangia in the mass basis."
DefaultParamCard::unknownmass= "The following massterm cannot be handled: `1`."
DefaultParamCard::massused= "The newly defined mass parameter \"`1`\"  might be already used elsewhere. Please check for naming overlap."


ExportUFO::missflavrot= "The masses `1` are not in the mass basis and require a falvor rotation to reach it. Off-diagonal entries are simply ignored for the masses in the UFO, which can lead to inconsitencies if the off-diagonal entries also appear in other couplings where they are not set to zero. Please use SetSSBReplacements to properly rotate to the mass basis.";


FindMasses[lagrangian_]:= FindMasses[lagrangian]= Module[
	{
		expr= LagrangianExpand@ lagrangian,
		LHcTerms, LHTerms,
		field, prop,
		fieldAssoc= $FieldAssociation,
		getMassLabel,
		res, res2, resGB,
		flavInd= Alternatives@@ Keys@ $FlavorIndices,
		nonMassBasisTerms,
		goldstones= GetFieldsByProperty[GoldstoneBoson-> True]
	}
	,
	(* extrac mass terms *)
	expr= RelabelIndices@ IsolateMassTerms[expr];
	
	(* find terms not in the mass basis due to missing flavor rotations *)
	nonMassBasisTerms= Cases[expr, x_/; !FreeQ[x, Except[Index[Global`d$$1, flavInd], Index[_, flavInd]]], 1];
	If[Length[nonMassBasisTerms]>0,
		Message[ExportUFO::missflavrot, Format[Total@ nonMassBasisTerms, NiceForm]];
	];

	(* make lists with atomics *)
	expr= TermsToList[HcSimplifyInternal[expr, InternalOpRepresentation-> True], Expand-> False];

	(* separate Hermitian and non-Hermitian parts *)
	LHcTerms= Total@ Cases[expr, _HcTerms]/. HcTerms-> Identity;
	LHTerms = Total[expr/. _HcTerms-> 0];

	Switch[Head[LHTerms],
		Plus,  LHTerms= List@@ LHTerms,
		Times, LHTerms= {LHTerms}
	];
	Switch[Head[LHcTerms],
		Plus,  LHcTerms= List@@ LHcTerms,
		Times, LHcTerms= {LHcTerms}
	];

	(* non-Hermitian masses *)
	LHcTerms= Table[
		(*find fields*)
		field= DeleteDuplicates[Cases[term, AtomicOp[{{lab_,lab_},_},___]:>lab, All]];
		If[Length@field != 1,
			Message[DefaultParamCard::massmix, Format[AtomicToNormalForm@term,NiceForm]]; 
			Abort[] (*mass mixing*)
			,
			field= First@ field;
			prop= $FieldAssociation[field];
			(* only chiral fermions have non-hermitian mass terms *)
			If[prop[Type]===Fermion && (!prop[SelfConjugate]) && MatchQ[prop[Chiral], LeftHanded|RightHanded],
				field-> -2*term/. _AtomicOp-> 1,
				Message[DefaultParamCard::unknownmass, Format[AtomicToNormalForm@term,NiceForm]];
				Abort[]
			]
		]
		,
		{term, LHcTerms}
	];

	(* Hermitian masses *)
	LHTerms= Table[
		(*find fields*)
		field = DeleteDuplicates[Cases[term, AtomicOp[{{lab_|Conj[lab_],lab_|Conj[lab_]},_},___]:>lab, All]];
		If[Length@field != 1,
			Message[DefaultParamCard::massmix, Format[AtomicToNormalForm@term,NiceForm]]; 
			Abort[] (*mass mixing*)
			,
			field= First@ field;
			prop= $FieldAssociation[field];
			Switch[prop[Type],
				Scalar,
					If[prop[SelfConjugate],
						field-> Sqrt[-2*term/. _AtomicOp-> 1], (*real*)
						field-> Sqrt[-term/. _AtomicOp-> 1]    (*complex*)
					],
				Vector,
					If[prop[SelfConjugate],
						field-> Sqrt[2*term/. _AtomicOp-> 1], (*real*)
						field-> Sqrt[term/. _AtomicOp-> 1]    (*complex*)
					],
				Fermion,
					If[prop[SelfConjugate],
						field-> -2*term/. _AtomicOp-> 1, (*real*)
						field-> -term/. _AtomicOp-> 1    (*complex*)
					],
				_ , 
					Message[DefaultParamCard::unknownmass, Format[AtomicToNormalForm@term,NiceForm]];
					Abort[]
			]
		]
		,
		{term, LHTerms}
	];
	
	(* turn into association and power expand *)
	res= PowerExpand/@ (Association@@ Join[LHTerms, LHcTerms]);
	
	(* go to mas basis *)
	res2= RelabelIndices/@ (res/. Index[Except[Global`d$$1], type: flavInd]:> Index[Global`d$$1, type]);
	
	(* print the affected couplings *)
	nonMassBasisTerms= Cases[Complement[res, res2], Coupling[_,{_,__},_], All];
	If[Length[nonMassBasisTerms] > 0,
		Echo[StringForm["Neglecting off-diagonal contributions from the couplings `1` to masses.", Format[nonMassBasisTerms, NiceForm]]]
	];
	res= res2;
	
	getMassLabel[l_]:= If[fieldAssoc[l, Mass]===0,
		(* get new unique mass label *)
		Module[{newLabel, flavCount= 1},
			newLabel= If[KeyExistsQ[fieldAssoc[l,UFO],"name"], 
				(* if UFO name for field has been specified *)
				If[Head[fieldAssoc[l,UFO]["name"]]===List,
					(* flavor multiplet *)
					(ToString["M"]<>#)&/@ fieldAssoc[l,UFO]["name"]
				,
					(* flavor singlet *)
					ToString["M"]<>fieldAssoc[l,UFO]["name"]
				]
			,
				(* if no UFO name has been specified *)
				If[Length[Cases[$FieldAssociation[l,Indices], flavInd]]>0,
					(* flavor multiplet *)
					ToString["M"]<>ToString[l]<>ToString[flavCount++]
				,
					(* flavor singlet *)
					ToString["M"]<>ToString[l]
				]
			];
			newLabel= ToExpression@ newLabel;
			
			(* Sow for printing of new masses *)
			Sow[{l,newLabel}];
			(* check if label has been defined before *)
			If[Defined[ToExpression@ #],
				Message[DefaultParamCard::massused, ToString[newLabel]]
			]&/@ Flatten[{newLabel}];
			
			newLabel
		]
	,
		fieldAssoc[l,Mass]
	];
	
	(* separate out goldstones *)
	resGB= KeySelect[res, MemberQ[goldstones,#]&];
	res=   KeySelect[res, !MemberQ[goldstones,#]&];
	
	(* get new mass labels *)
	res= Reap[
		res=   Association@@ KeyValueMap[(#1-> {#2, getMassLabel[#1]})&, res];
		(* for GBs check if mass matches a vector mass *)
		resGB= Association@@ KeyValueMap[(#1-> {#2, FirstCase[Values[res], {#2/. _\[Xi]Gauge-> 1, x_}-> x, getMassLabel[#1], 1]})&, resGB];
		Join[res, resGB]
	];
	
	(* print warning if new masses have been defined *)
	If[Length[res[[-1]]]>0,
		Echo["The fields "<>ToString[Transpose[res[[-1,1]]][[1]]]<>" have been defined as massless, but a mass term was found in the Lagrangian. Created new mass labels "<>ToString[Transpose[res[[-1,1]]][[2]]]<>"."];
	];
	
	res[[1]]
]


(* ::Subsubsection::Closed:: *)
(*Creating output directory*)


CreateNewDir[dir_] := If[DirectoryQ[dir],
	(* EXISTING DIRECTROY *)
	If[Or[$FrontEnd === Null, (* check if frontend is available *)
		(* open interactive dialog for confirming overwriting *)
		DialogInput@ Column@ {
			Style["The directory \"" <> ToString[dir] <> "\" already exists.", Bold, 14],
			Spacer[10],
			Style["Do you want to overwrite it?", Bold, 14],
			Spacer[10],
			Row@ {
				DefaultButton["Yes", DialogReturn[True]],
				Spacer[20],
				CancelButton["No", DialogReturn[False]]
			}
		}]
	,
		(* overwriting *)
		(*Echo["Overwriting the directory: "<> dir];*)
		DeleteDirectory[dir, DeleteContents->True];
		CreateDirectory[dir]
	,
		(* NOT overwriting *)
		$Failed
	]
	,
	(* NEW DIRECTROY *)
	(*Echo["Creating the directory: "<> dir];*)
	CreateDirectory[dir]
]


(*
CreateNewDir[] := Module[
	{
		parent, 
		newDirName, 
		newDirPath
	},
	
	(* Let the user choose the parent directory *)
	parent = SystemDialogInput["Directory"];
	If[parent === $Canceled, Abort[]];
	
	(* Prompt the user for the new directory name *)
	newDirName = InputString["Enter the name for the new directory:"];
	If[!StringQ[newDirName] || newDirName === "", 
    Return[Message[InputString::nval, newDirName]]
	];
	
	(* Create the new directory if it doesn't already exist *)
	newDirPath = FileNameJoin[{parent, newDirName}];
	If[DirectoryQ[newDirPath],
		Print["Directory already exists at: ", newDirPath];
		Abort[],
		CreateDirectory[newDirPath];
		Print["Directory created at: ", newDirPath]
	];
	newDirPath
];
*)


(* ::Subsubsection::Closed:: *)
(*Determining Dirac structures*)


(* Contraction fermion*fermion *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] == {} && !FreeQ[term, Fermion, {0, Infinity}] := Module[
	{
		spinorStructure,
		spinorFields
	}
	,
	spinorStructure = Cases[term, _NCM, {0, Infinity}];
	(*
	spinorFields    = Cases[spinorStructure, Field[___], {0, Infinity}] /. Bar[x_] :> x;
	
	{Delta[Index[spinorFields[[1, 1, 1, 2]], Spinor], Index[spinorFields[[2, 1, 1, 2]], Spinor]]}
	*)
	
	Table[
		spinorFields= Cases[spinChain, Field[___], {0, Infinity}] /. Bar[x_] :> x;
		Delta[Index[spinorFields[[1, 1, 2]], Spinor], Index[spinorFields[[2, 1, 2]], Spinor]]
		,
		{spinChain, spinorStructure}
	]
]


(* Contraction fermion*spinorstructure*fermion *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] =!= {} := Module[
	{
		diracStructures,
		positionsDirac,
		diracStructureFields,
		listInd = {},
		dummyIndexNumber = 1,
		result,
		dummyrules,
		spinorStructure, spinorFields
	},
	(* find Dirac structures and their positions *)
	diracStructures = Cases[term, DiracProduct[___], {0, Infinity}];
	positionsDirac  = Position[term, _?(MemberQ[diracStructures, #] &), Infinity];
	
	(* loop over all spin chains in the given term *)
	Do[
		(* get all fermion fields *)
		diracStructureFields = Cases[Extract[term, Most@dirac], Field[___], {0, Infinity}] /. Bar[x_] :> x;
		(* Add UFO compatible spinor-index pairs depending on how many Dirac structures are in one spin chain *)
		Switch[Length@ Extract[term, dirac],
			1,
				AppendTo[listInd, {Index[diracStructureFields[[1, 1, 2]], Spinor], Index[diracStructureFields[[2, 1, 2]], Spinor]}],
			2,
				listInd = Join[
					listInd,
					{{Index[diracStructureFields[[1, 1, 2]], Spinor], Index[-dummyIndexNumber, Spinor]}, {Index[-dummyIndexNumber, Spinor], Index[diracStructureFields[[2, 1, 2]], Spinor]}}
				];
				dummyIndexNumber++,
			_,
				AppendTo[listInd, {Index[diracStructureFields[[1, 1, 2]], Spinor], Index[-dummyIndexNumber, Spinor]}];
				Do[
					AppendTo[listInd, {Index[-dummyIndexNumber, Spinor],Index[-dummyIndexNumber-1, Spinor]}];
					dummyIndexNumber++,
					{i, Length[Extract[term, dirac]] - 2}
				];
				AppendTo[listInd, {Index[-dummyIndexNumber, Spinor], Index[diracStructureFields[[2, 1, 2]], Spinor]}];
				dummyIndexNumber++
		],
		{dirac, positionsDirac}
	];
	
	(* Add spinor-index pair just determined to the dirac structures *)
	result = MapThread[
		Function[
			{expr,idxList},
			If[Head[expr] === Symbol,
				expr[idxList],
				(*Transposed objects get reversed spinor indices*)
				If[Head[expr] === Transp,
					If[Head[First@expr] === Symbol,
						First[expr][Reverse[idxList]],
						Head[First[expr]][Sequence @@ (List @@ First[expr]), Reverse[idxList]]
					],
					Head[expr][Sequence @@ (List @@ expr), idxList]
				]
			]
		],
		{Flatten[List @@@ diracStructures], listInd}
	];
	
	(* create replacement rule potential Lorentz dummies following UFO convention *)
	dummyrules = MapIndexed[
		#1 -> Index[-(dummyIndexNumber + First[#2]) - 1, Lorentz]&, 
		Cases[FindDummyIndices[result], Index[_, Lorentz]]
	];
	
	result= result /. dummyrules;
	
	(* include Delta for ferm*ferm currents *)
	spinorStructure = Cases[term, ncm_NCM/;FreeQ[ncm,_DiracProduct,All], {0, Infinity}];
	result= result~Join~Table[
		spinorFields= Cases[spinChain, Field[___], {0, Infinity}] /. Bar[x_] :> x;
		Delta[Index[spinorFields[[1, 1, 2]], Spinor], Index[spinorFields[[2, 1, 2]], Spinor]]
		,
		{spinChain, spinorStructure}
	];
	
	result
]


(* Exception *)
GetDiracStructureSingleTerm[term_] /; Cases[term, DiracProduct[___], {0, Infinity}] == {} := {};


(* ::Subsubsection::Closed:: *)
(*Get all fields in one term*)


(* ::Text:: *)
(*Obtains a list of all fields and bared fields in the given expression (neglecting derivatives)*)


GetFieldsInTerm[expr_Plus] := Flatten[GetFieldsInTerm /@ List@@expr];

GetFieldsInTerm[expr_]:= Join[
	Cases[
		DeleteCases[
			Hold@ expr, 
			Bar[___], 
			{0,Infinity}
		], 
		Field[___], 
		{0,Infinity}
	], 
	Cases[
		expr, 
		Bar[Field[___]], 
		{0,Infinity}
	]
]/. Field[label_, type_, indices_, ___] :> Field[label, type, indices, {}]		


(* ::Subsubsection::Closed:: *)
(*Get new and unused PDG code*)


(* ::Text:: *)
(*This function returns a list of all PDG codes currently in use.*)


GetAllDefinedPDGCodes[]:= Sort@ Flatten[#["pdg"]&/@ ((#[UFO]&/@ Values[$FieldAssociation])/. <||>-> Nothing)]


(* ::Text:: *)
(*This function returns a new PDG code not yet used. By default it will take the next largest integer equal or greater than 50k that is not in use. Optionally, one can also provide an integer as argument and the next greater integer from this number will be used.*)


$autoGeneratedPDGCodes= {}; (* this ensures that GetNewPDGCode never returns the same number twice *)


GetNewPDGCode[n_Integer:50000]:= Module[{currentPDG= $autoGeneratedPDGCodes~Join~GetAllDefinedPDGCodes[], newPDG= n},
	(* increase newPDG until an unused integer is found *)
	While[MemberQ[currentPDG, newPDG],
		newPDG++
	];
	AppendTo[$autoGeneratedPDGCodes, newPDG]; (* this remebers not to reuse this number *)
	newPDG
]


(* ::Subsubsection::Closed:: *)
(*GaugeFixingUFO*)


(* ::Text:: *)
(*This function check the gauge fixing terms (and includes them if not present yet) for the UFO functions (ExportUFO and DefaultParamCard). *)
(*INPUT: a Lagrangian, as well as the option values received for Gauge and R\[Xi].*)
(*OUTPUT: {Lagrangian, gauge fixing information}*)


GaugeFixingUFO::partialGF = "The given Lagrangian appears to be already partially gauge fixed using non-unitary gauges: `1`. When providing partially gauge-fixed Lagrangians only unitary gauge is allowed. For R\[Xi] gauge only fully fixed R\[Xi]-gauge Lagrangian are allowed."
GaugeFixingUFO::noGF      = "The given Lagrangian is not fully gauge fixed yet. Use the function 'GaugeFixLagrangian' first or set the option 'Gauge' to 'Automatic', 'Unitary', or 'R\[Xi]'."
GaugeFixingUFO::wrongGF   = "The given Lagrangian cannot be fixed to `1` gauge. The determined gauges are: `2`."
GaugeFixingUFO::unknownGF = "The gauge fixing option `1` is unknown."


Options@ GaugeFixingUFO= {Verbose-> True}


GaugeFixingUFO[lagrangian_, gauge_, r\[Xi]_, OptionsPattern[]]:= Module[
	{
		lag= lagrangian,
		initialGF,
		finalGF,
		ghostLag,
		massiveVectors,
		echo= If[OptionValue[Verbose], Echo, Identity](* suppress Echo if needed *)
	}
	,
	(* get all massive vectors *)
	massiveVectors= MassiveVectorsInLag[lag, Light-> True];
	
	(* determine GF status of given Lagrangian *)
	initialGF= Quiet[NonGaugeFixedLagQ[lag, DetailedOutput-> True], {NonGaugeFixedLagQ::nongaugekinfix, NonGaugeFixedLagQ::noghostkin}];
	
	(* check for and include GF depending on the OptionValue of Gauge *)
	Switch[gauge,
		Automatic,
			If[initialGF[[1]],
				(* GF in R\[Xi] gauge if no GF present at all *)
				lag= GaugeFixLagrangian[lag, Gauge-> R\[Xi], SubstituteMasses-> True];
				Switch[r\[Xi],
					1,       echo["Gauge fixing using Feynman gauge. Note: this UFO will still be usable in unitary gauge. For MadGraph use 'set gauge unitary' to run in unitary gauge (default) and 'set gauge Feynman' to actually use Feynman gauge."],
					General, echo["Gauge fixing using R\[Xi] gauge."]
				]
			,
				If[ContainsNone[initialGF[[2]], {Indeterminate}],
					(* fully gauge fixed already *)
					echo@ StringForm["The input Lagrangian is already fully gauge fixed using the gauges: `1`.", initialGF[[2]]]
				,
					If[ContainsOnly[initialGF[[2]]/@ massiveVectors, {Unitary, Indeterminate}] && ContainsOnly[KeyDrop[initialGF[[2]], massiveVectors], {Indeterminate}],
						lag= GaugeFixLagrangian[lag, Full-> False, Gauge-> Unitary, SubstituteMasses-> True];
						Switch[r\[Xi],
							1,       echo["Gauge fixing using unitary gauge (massless vectors in Feynman gauge)."],
							General, echo["Gauge fixing using unitary gauge (massless vectors in R\[Xi] gauge)."]
						]
					,
						Message[GaugeFixingUFO::partialGF, initialGF[[2]]];
						Abort[]
					]
				]
			],
		Manual,
			(* if lag is not yet fully GF abort *)
			If[initialGF[[1]] || ContainsAny[initialGF[[2]], {Indeterminate}],
				Message[GaugeFixingUFO::noGF];
				Abort[]
			],
		R\[Xi],
			If[initialGF[[1]],
				(* no previous gauge fixing at all *)
				lag= GaugeFixLagrangian[lag, Gauge-> R\[Xi], SubstituteMasses-> True];
				Switch[r\[Xi],
					1,       echo["Gauge fixing using Feynman gauge. Note: this UFO will still be usable in unitary gauge. For MadGraph use 'set gauge unitary' to run in unitary gauge (default) and 'set gauge Feynman' to actually use Feynman gauge."],
					General, echo["Gauge fixing using R\[Xi] gauge."]
				]
			,
				(* if previous GF present it must be all R\[Xi] *)
				If[!ContainsOnly[initialGF[[2]],{R\[Xi]}],
					(* abort if it is not fully R\[Xi] *)
					Message[GaugeFixingUFO::wrongGF, "R\[Xi]", initialGF[[2]]];
					Abort[]
				]
			],
		Unitary,
			If[initialGF[[1]],
				(* no previous gauge fixing at all *)
				lag= GaugeFixLagrangian[lag, Gauge-> Unitary, SubstituteMasses-> True];
				Switch[r\[Xi],
					1,       echo["Gauge fixing using unitary gauge (massless vectors in Feynman gauge)."],
					General, echo["Gauge fixing using unitary gauge (massless vectors in R\[Xi] gauge)."]
				]
			,
				(* when previous GF present ensure massive vectors are Unitary or Indeterminate *)
				If[ContainsOnly[initialGF[[2]]/@ massiveVectors, {Unitary, Indeterminate}],
					If[ContainsOnly[KeyDrop[initialGF[[2]], massiveVectors], {Indeterminate}],
						(* if massless vectors are Indeterminate gauge fix again *)
						lag= GaugeFixLagrangian[lag, Full-> False, Gauge-> Unitary, SubstituteMasses-> True];
						Switch[r\[Xi],
							1,       echo["Gauge fixing using unitary gauge (massless vectors in Feynman gauge)."],
							General, echo["Gauge fixing using unitary gauge (massless vectors in R\[Xi] gauge)."]
						]
					,
						(* If not full GF already throw an abort *)
						If[!(ContainsOnly[initialGF[[2]]/@ massiveVectors, {Unitary}] && ContainsOnly[KeyDrop[initialGF[[2]], massiveVectors], {R\[Xi]}]),
							Message[GaugeFixingUFO::wrongGF, "unitary", initialGF[[2]]];
							Abort[]
						]
					];
					
				,
					(* abort if it is not Unitary or Indeterminate *)
					Message[GaugeFixingUFO::wrongGF, "unitary", initialGF[[2]]];
					Abort[]
				]
			],
		_,
			Message[GaugeFixingUFO::unknownGF, gauge];
			Abort[]
	];
	
	(* determine GF of the output Lagrangian *)
	
	finalGF= Quiet[NonGaugeFixedLagQ[lag, DetailedOutput-> True], {NonGaugeFixedLagQ::nongaugekinfix, NonGaugeFixedLagQ::noghostkin}];
	
	(* save the ghost Lagrangian since the lag obtained below by setting ghosts to zero is not gauge invariant *)
	ghostLag= lag;
	
	(* drop the ghosts since we work at tree level for now *)
	lag= lag/. Field[_, Ghost|AntiGhost, ___]-> 0;
	(* print a warning only if ghosts were present in the input already *)
	If[!FreeQ[lagrangian, Ghost|AntiGhost, All],
		echo["Dropping ghost terms for tree-level UFO..."]
	];
	
	(* for Feynman gauge replace all \[Xi]Gauge[_] by unity *)
	If[r\[Xi]===1,
		lag= lag/. \[Xi]Gauge[_]-> 1;
		If[!initialGF[[1]],
			echo["Employing Feynman gauge and setting \[Xi]=1 for all gauge parameters. Use the option 'R\[Xi]\[Rule]General' to keep generic gauge parameters (\[Xi])."];
		];
	];
	
	(* return lag and information about its GF terms *)
	{lag, ghostLag, finalGF}
]


(* ::Subsection:: *)
(*UFO translation functions*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


TranslateToUFO::unidentified = "Structure cannot be identified: `1`";
TranslateToUFO::noSextet     = "The generator `1` cannot be exported to UFO. It does not correspond to the fundamental or adjoint representation, neither could a sextet be identified.";
TranslateToUFO::genUnknown   = "The generator `1` cannot be exported to UFO. It does not correspond to the fundamental, adjoint, or sextet representation.";
TranslateToUFO::unknownCG    = "The CG coefficient `1` cannot be exported to UFO.";


(* ::Subsubsection::Closed:: *)
(*TranslateToUFO*)


(* Metrics *)
TranslateToUFO[expr_Metric] := Module[
	{inds},
	inds = expr /. Metric[Index[d_, Lorentz], Index[e_, Lorentz]] :> {d, e};
	"Metric(" <> StringJoin[Riffle[ToString /@ inds, ","]] <> ")"
];

(* Gamma matrices *)
TranslateToUFO[expr:HoldPattern@GammaM[_Index,_List]] := Module[
	{inds}, 
	inds = expr /. GammaM[Index[d_, Lorentz], {Index[e_, Spinor], Index[f_, Spinor]}] :> {d, e, f};
	"Gamma(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* \[Sigma] matrices using \[Sigma][a,b] = \[ImaginaryI] GammaM[a,b] *)
TranslateToUFO[expr:HoldPattern@GammaM[_Index, _Index, _List]] := Module[
	{inds}, 
	inds = expr /. GammaM[Index[d1_, Lorentz], Index[d2_, Lorentz], {Index[e_, Spinor], Index[f_, Spinor]}] :> {d1, d2, e, f};
	"complex(0,-1)*Sigma(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
	(* rescaling by factor of 2 for "MadGraph" mode is done in ExportUFO on Lagrangian level *)
];

(* Gamma5 matrix *)
TranslateToUFO[expr_Gamma5] := Module[
	{inds}, 
	inds = expr /. Gamma5[{Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"Gamma5(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Charge conjugation operator *)
TranslateToUFO[expr_GammaCC] := Module[
	{inds}, 
	inds = expr /. GammaCC[{Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"C(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Projectors *)
TranslateToUFO[expr_Proj] /; MatchQ[expr, Proj[-1, _]] := Module[
	{inds}, 
	inds = expr /. Proj[_, {Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"ProjM(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];
TranslateToUFO[expr_Proj] /; MatchQ[expr, Proj[1, _]] := Module[
	{inds}, 
	inds = expr /. Proj[_, {Index[e_, Spinor], Index[f_, Spinor]}] :> {e, f};
	"ProjP(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Momenta *)
TranslateToUFO[expr_ExtMom] := Module[
	{inds}, 
	inds = expr /. ExtMom[Index[e_, Lorentz], f_] :> {e, f};
	"P(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* Deltas *)
TranslateToUFO[Delta[Index[a_, _], Index[b_, _]]] := "Identity(" <> ToString@a <> "," <> ToString@b <> ")";

(* Levi Civita Tensors *)
TranslateToUFO[LCTensor[Index[a_, _], Index[b_, _], Index[c_, _], Index[d_, _]]] := "Epsilon(" <> StringJoin[Riffle[ToString /@ {a, b, c, d}, ","]] <> ")";

(* anti-symmetric structure functions *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[fStruct[_], ___]] := Module[
	{inds},
	inds = expr /. CG[fStruct[_], {Index[d_, _], Index[e_, _], Index[f_, _]}] :> {d, e, f};
	"f(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* fully symmetric CG functions *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[dSym[_], ___]] := Module[
	{inds},
	inds = expr /. CG[dSym[_], {Index[d_, _], Index[e_, _], Index[f_, _]}] :> {d, e, f};
	"d(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* fundamental generators of SU(n) *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[gen[_[fund]], ___]] := Module[
	{inds},
	inds = expr /. CG[gen[_[fund]], {Index[d_, _], Index[e_, _], Bar[Index[f_, _]]}] :> {d, e, f};
	"T(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* anti-symmetric color tensor [fundamental and anti-fundamental] *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[eps[_], ___]] := Module[
	{inds},
	inds = expr /. CG[eps[_], {Index[d_, _], Index[e_, _], Index[f_, _]}] :> {d, e, f};
	"Epsilon(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[Bar@eps[_], ___]] := Module[
	{inds},
	inds = expr /. CG[Bar@eps[_], {Bar@Index[d_, _], Bar@Index[e_, _], Bar@Index[f_, _]}] :> {d, e, f};
	"EpsilonBar(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* del CG *)
(* should never be used since this is capured by the Deltas rule above *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[del[_], ___]] := Module[
	{inds},
	inds = expr /. CG[del[_], {Index[d_, _], Bar@Index[e_, _]}] :> {e, d};
	"Identity(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
];

(* sextet generator *)
TranslateToUFO[expr_CG] /; MatchQ[expr, CG[gen[Except[_@fund|_@adj]], ___]] := Module[
	{inds, reps, sextet}
	,
	(* list of QCD reps *)
	reps=   GetGaugeGroups[$QCD$UFO, Representations];
	reps=   KeySelect[GetRepresentations[], MatchQ[Alternatives@@reps]];
	(* sextet rep *)
	sextet= Select[reps, MatchQ[KeyValuePattern[DynkinCoefficients->{2,0}]]];
	sextet= First[Keys@sextet, Message[TranslateToUFO::noSextet, expr];Abort[]];
	
	(* if given CG is the sextet generator *)
	If[MatchQ[expr, CG[gen[sextet],___]],
		inds = expr /. CG[gen[sextet], {Index[d_, _], Index[e_, _], Bar@Index[f_, _]}] :> {d, f, e};
		"T6(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"
		,
		Message[TranslateToUFO::genUnknown, expr]; Abort[]
	]
];

(* CG contracting sextet and 2 funds *)
TranslateToUFO[expr_CG]:= Module[
	{inds, reps, sextet, cgs, k6, k6anti}
	,
	(* list of QCD reps *)
	reps=   GetGaugeGroups[$QCD$UFO, Representations];
	reps=   KeySelect[GetRepresentations[], MatchQ[Alternatives@@reps]];
	(* sextet rep *)
	sextet= Select[reps, MatchQ[KeyValuePattern[DynkinCoefficients->{2,0}]]];
	sextet= First[Keys@sextet, Message[TranslateToUFO::noSextet, expr];Abort[]];
	(* get all QCD CGs *)
	cgs=    GetGaugeGroups[$QCD$UFO, ClebschGordanCoefficients];
	(* get the CG properties *)
	cgs=    KeySelect[$CGproperties, MatchQ[Alternatives@@cgs]];
	(* select CG contracting the sextet with two anti-fundamentals *)
	k6=     Select[cgs, MatchQ[KeyValuePattern[Indices->{OrderlessPatternSequence[sextet, Bar@$QCD$UFO[fund], Bar@$QCD$UFO[fund]]}]]];
	(* if there are multiple select the one in the basis *)
	If[Length[k6]>1,
		k6= Select[k6, MatchQ[KeyValuePattern[InBasis->True]]];
	];
	k6= Keys@k6;
	(* K6 generator *)
	If[1==Length@k6,
		If[MatchQ[expr, CG[k6,___]],
			inds = expr /. CG[k6, {OrderlessPatternSequence[Index[d_, sextet], Bar@Index[e_, $QCD$UFO[fund]], Bar@Index[f_, $QCD$UFO[fund]]]}] :> {d, e, f};
			Return["K6(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"]
			,
			If[MatchQ[expr, CG[Bar@k6,___]],
				inds = expr /. CG[Bar@k6, {OrderlessPatternSequence[Bar@Index[d_, sextet], Index[e_, $QCD$UFO[fund]], Index[f_, $QCD$UFO[fund]]]}] :> {d, e, f};
				Return["K6Bar(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"]
				,
				Message[TranslateToUFO::unknownCG, expr]; Abort[]
			]
		]
		,
		(* in case K6 was defined through its bar *)
		k6anti= Select[cgs, MatchQ[KeyValuePattern[Indices->{OrderlessPatternSequence[Bar@sextet, $QCD$UFO[fund], $QCD$UFO[fund]]}]]];
		(* if there are multiple select the one in the basis *)
		If[Length[k6anti]>1,
			k6anti= Select[k6anti, MatchQ[KeyValuePattern[InBasis->True]]];
		];
		k6anti= Keys@k6anti;
		If[1==Length@k6anti,
			If[MatchQ[expr, CG[k6anti,___]],
				inds = expr /. CG[k6anti, {OrderlessPatternSequence[Bar@Index[d_, sextet], Index[e_, $QCD$UFO[fund]], Index[f_, $QCD$UFO[fund]]]}] :> {d, e, f};
				Return["K6Bar(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"]
				,
				If[MatchQ[expr, CG[Bar@k6anti,___]],
					inds = expr /. CG[Bar@k6anti, {OrderlessPatternSequence[Index[d_, sextet], Index[e_, Bar@$QCD$UFO[fund]], Bar@Index[f_, $QCD$UFO[fund]]]}] :> {d, e, f};
					Return["K6(" <> StringJoin[Riffle[ToString/@inds,","]] <> ")"]
					,
					Message[TranslateToUFO::unknownCG, expr]; Abort[]
				]
			]
			,
			Message[TranslateToUFO::unknownCG, expr]; Abort[]
		]
	]
]

(* Couplings (placeholder)*)
(*TranslateToUFO[c_Coupling] := StringJoin[
	TranslateToUFO @ c[[1]],
	Sequence @@ (ToString /@ c[[2]])
];*)
TranslateToUFO[c_Coupling]:= Block[{name},
	If[KeyExistsQ[$CouplingAssociation[First@c, UFO],"name"],
		$CouplingAssociation[First@c, UFO, "name"][[Sequence@@(c[[2]])]]
	,
		name= TranslateToUFO @ c[[1]];
		If[Length@ c[[2]] > 0, 
			StringJoin[name, "_",
				Sequence @@ (ToString /@ c[[2]])
			]
		,
			name
		]
	]
];

(* Operations *)
TranslateToUFO[Plus[a_, b__]] := StringJoin @ Riffle[TranslateToUFO /@ {a, b}, "+"];
TranslateToUFO[Times[a_, b__]] := StringJoin @ Riffle[TranslateToUFO /@ {a, b}, "*"];
TranslateToUFO[DiracProduct[a___]] := StringJoin @ Riffle[TranslateToUFO /@ {a}, "*"];
TranslateToUFO[Rational[a_, b_]] := StringJoin @ Riffle[TranslateToUFO /@ {a, b}, "/"];
TranslateToUFO[Power[a_, 1/2]] := "cmath.sqrt(" <> TranslateToUFO[a] <> ")";
TranslateToUFO[Power[a_, -1/2]] := "1/cmath.sqrt(" <> TranslateToUFO[a] <> ")";
TranslateToUFO[Power[a_, 2]] := If[MatchQ[Head[a], Coupling|ExtMom],
	TranslateToUFO[a] <> "**2"
,
	"(" <> TranslateToUFO[a] <> ")**2"
];
TranslateToUFO[Power[a_, exponent_]] := "(" <> TranslateToUFO[a] <> ")**(" <> TranslateToUFO[exponent] <> ")";
TranslateToUFO[Bar[a_]] := "complexconjugate(" <> TranslateToUFO@a <> ")";
TranslateToUFO[Conjugate[a_]] := "complexconjugate(" <> TranslateToUFO@a <> ")";
TranslateToUFO[TempSymb[a_String]] := a;

(* Standard expressions *)
TranslateToUFO[1] := "1";
TranslateToUFO[0] := "0";
TranslateToUFO[n_Integer] := ToString @ n <> ".";
TranslateToUFO[n_Symbol] := StringReplace[SpokenString[n], {" " -> "", "$" -> ""}];
TranslateToUFO[z_Complex] := "complex("<>TranslateToUFO[Re @ z]<>","<>TranslateToUFO[Im @ z]<>")";
TranslateToUFO[a_String] := a;

(* gauge parameters *)
TranslateToUFO[\[Xi]Gauge[f_]] := "xi"<>If[KeyExistsQ[$FieldAssociation[f, UFO], "name"], $FieldAssociation[f, UFO]["name"], ToString[f]];

(* Error *)
TranslateToUFO[expr_] := Message[TranslateToUFO::unidentified, expr];


(*NonRedundantUFOParams[\[Xi]Gauge[f_],___]:= {{
	"xi"<>If[KeyExistsQ[GetFields[f, UFO], "name"], GetFields[f, UFO]["name"], ToString[f]], 
	1, (* default to Feynman gauge *)
	{"real", "external"}
}}*)


(* ::Subsubsection::Closed:: *)
(*Renaming of dummy indices to UFO conventions*)


(* ::Text:: *)
(*Replace dummy indices by negative integers following UFO convention*)


RenameDummyIndicesToUFO[sum_Plus] := RenameDummyIndicesToUFO /@ sum;

RenameDummyIndicesToUFO[term_] := Module[
	{
		flavors = GetFlavorIndices[],
		dummies,
		replacements,
		counter = -1
	}
	,
	dummies = FindDummyIndices[term /. Index[_, Alternatives @@ Keys @ flavors] -> 1];
	replacements = Association@ Table[
		sym /. Index[sym_, grp_] :> Index[sym, grp] -> Index[ToString[counter--], grp],
		{sym, dummies}
	];
	
	term /. Index[sym_, grp_] :> Lookup[replacements, Index[sym, grp], Index[sym, grp]]
];


(* ::Subsection:: *)
(*ExplicitFlavorStructure*)


(* ::Subsubsection::Closed:: *)
(*Caching of flavor expansions*)


FlavorExpandedLagrangian[lag_]:= FlavorExpandedLagrangian[lag]= ExplicitEinsteinSums@LagrangianExpand[lag]


(* ::Subsubsection::Closed:: *)
(*Explicitly carry out Einstein summation over flavor indices*)


(* ::Text:: *)
(*This function performs the explicit summation over flavor indices.*)


(*Options@ExplicitEinsteinSums= {List->False}


ExplicitEinsteinSums[sum_, OptionsPattern[]]:= If[OptionValue[List],
	Flatten[ExplicitEinsteinSumTerm/@ TermsToList[sum, Expand->True]],
	Plus@@ Flatten[ExplicitEinsteinSumTerm/@ TermsToList[sum, Expand->True]]
]


ExplicitEinsteinSumTerm[term_, OptionsPattern[]]:= Module[
	{
		repeatedInds,
		flavors = GetFlavorIndices[],
		res
	},
	
	(* determine repeated indices *)
	repeatedInds = FindDummyIndices[term];
	(* filter for flavor indices *)
	repeatedInds = DeleteCases[repeatedInds, Index[_, Except[Alternatives @@ Keys @ flavors]]];
	(* ensure that repeated indices are ordered in the same way in different terms *)
	repeatedInds = Sort@repeatedInds;
	(* determine dimensions of these indices *)
	repeatedInds = Table[
		{ind, flavors[Last[ind]][IndexDimension]}
		,
		{ind, repeatedInds}
	];
	(* perform explicit sums over flavor indices *)
	res= If[Length[repeatedInds]>0, 
		Flatten@ Table[term, Evaluate[Sequence@@ repeatedInds]], 
		term
	]/. FlavorSum[_?NumericQ]->1;
	
	(* replace flavor deltas *)
	res/. Delta[n_Integer,m_Integer] :> If[n==m,1,0]
]*)


(* ::Text:: *)
(*Old code version below:*)
(*This old verion assumes its argument is already fully expanded. It will give wrong results (double counting of some terms) if that is not the case!*)


(*ExplicitEinsteinSums[sum_Plus, OptionsPattern[]]:= If[OptionValue[List],
	Flatten[ExplicitEinsteinSums[#, List->OptionValue[List]]&/@ (List@@sum)],
	Plus@@ Flatten[ExplicitEinsteinSums[#, List->OptionValue[List]]&/@ (List@@sum)]
]*)


(*ExplicitEinsteinSums[term:Except[_Plus], OptionsPattern[]]:= Module[
	{
		repeatedInds,
		flavors = GetFlavorIndices[],
		res
	},
	
	(* determine repeated indices *)
	repeatedInds = FindDummyIndices[term];
	(* filter for flavor indices *)
	repeatedInds = DeleteCases[repeatedInds, Index[_, Except[Alternatives @@ Keys @ flavors]]];
	(* ensure that repeated indices are ordered in the same way in different terms *)
	repeatedInds = Sort@repeatedInds;
	(* determine dimensions of these indices *)
	repeatedInds = Table[
		{ind, flavors[Last[ind]][IndexDimension]}
		,
		{ind, repeatedInds}
	];
	(* perform explicit sums over flavor indices *)
	res= If[OptionValue[List],
		If[Length[repeatedInds] > 0, Flatten@Table[term, Evaluate[Sequence @@ repeatedInds]], term]/.FlavorSum[_?NumericQ]->1
		,
		If[Length[repeatedInds] > 0, Sum[term, Evaluate[Sequence @@ repeatedInds]], term]/.FlavorSum[_?NumericQ]->1
	];
	
	(* replace flavor deltas *)
	res/. Delta[n_Integer,m_Integer] :> If[n==m,1,0]
]*)


Options@ExplicitEinsteinSums= {List->False, Internal->True}


ExplicitEinsteinSums[sum_Plus, OptionsPattern[]]:= ExplicitEinsteinSums[#, List->OptionValue[List], Internal->OptionValue@Internal]&/@ sum


ExplicitEinsteinSums[term:Except[_Plus], OptionsPattern[]]:= Module[
	{
		repeatedInds,
		flavors = GetFlavorIndices[],
		res,
		fields, barFields, 
		fieldInds, repeatedIntInds, repeatedExtInds
	},
	
	(* determine repeated indices *)
	repeatedInds = FindDummyIndices[term];
	(* filter for flavor indices *)
	repeatedInds = DeleteCases[repeatedInds, Index[_, Except[Alternatives@@ Keys@ flavors]]];
	(* find out the dummy indices carried by external fields *)
	fields = Cases[
		DeleteCases[
			Hold@ term, 
			Bar[___], 
			{0,Infinity}
		],
		Field[Ext[label_, _], rest__] -> Field[label, rest],
		{0,Infinity}
	];
	barFields = Cases[
		term, 
		Bar[Field[Ext[(*OverBar|*)Conj[label_], _], rest__]] -> Bar[Field[label, rest]], 
		{0,Infinity}
	];
	If[Length[Join[fields,barFields]]>0,
		(*Get tuples of index and corresponding field*)
		fieldInds = Tuples/@ If[!FreeQ[Join[fields, barFields], Index[_, Alternatives@@ Keys@ flavors]],
			Join[
				({List@ReplacePart[#, 3->{}], Cases[#[[3]], Index[_, Alternatives@@ Keys@ flavors]]}&)/@ fields,
				({List@ReplacePart[#, {1,3}->{}], Cases[First[#][[3]], Index[_, Alternatives@@ Keys@ flavors]]}&)/@ barFields
			],
			{}
		];
		fieldInds = Flatten[fieldInds, 1];
		(*Sort so that that indices are summed in the same order in all terms of a sum*)
		fieldInds = SortBy[fieldInds, Part[#, 2]&];

		(* divide into external and internal dummy indices *)
		repeatedIntInds = Complement[repeatedInds, Last/@ fieldInds];
		(* determine dimensions of these indices *)
		repeatedIntInds = Table[
			{ind, flavors[Last[ind]][IndexDimension]}
			,
			{ind, repeatedIntInds}
		];
		(*repeatedExtInds = Table[
			{ind, flavors[Last[ind]][IndexDimension]}
			,
			{ind, Last /@ fieldInds}
		];*)
		(*Make sure that vertices with the same external fields are not counted multiple times*)
		fieldInds = Map[Last, GatherBy[fieldInds, {First@#, Last@Last@#}&], {2}];

		repeatedExtInds = Sequence@@@ Map[
			Function[inds,
				Transpose@List[inds, Join[{1}, Most@inds], ConstantArray[flavors[Last[First@inds]][IndexDimension], Length@inds]]
			],
			fieldInds
		];
		,
		repeatedIntInds = Table[
			{ind, flavors[Last[ind]][IndexDimension]}
			,
			{ind, repeatedInds}
		];
		repeatedExtInds = {};
	];

	(* perform explicit sums over flavor indices *)
	res= If[OptionValue[List],
		If[Length[repeatedExtInds]>0 && Length[repeatedIntInds]>0, 
			Table[
				If[OptionValue@Internal, 
					Sum[term, Evaluate[Sequence@@ repeatedIntInds]], 
					term
				],
				Evaluate[Sequence@@ repeatedExtInds]
			],
			If[Length[repeatedExtInds]>0,
				Table[
					term,
					Evaluate[Sequence@@ repeatedExtInds]
				],
				If[Length[repeatedIntInds]>0,
					If[OptionValue@Internal, 
						Sum[term, Evaluate[Sequence@@ repeatedIntInds]], 
						term
					],
					term
				]
			]
		]/. FlavorSum[_?NumericQ]-> 1
		,
		If[Length[repeatedExtInds]>0 && Length[repeatedIntInds]>0, 
			Sum[
				If[OptionValue@Internal, 
					Sum[term, Evaluate[Sequence@@ repeatedIntInds]], 
					term
				],
				Evaluate[Sequence@@ repeatedExtInds]
			],
			If[Length[repeatedExtInds]>0,
				Sum[
					term,
					Evaluate[Sequence@@ repeatedExtInds]
				],
				If[Length[repeatedIntInds]>0,
					If[OptionValue@Internal, 
						Sum[term, Evaluate[Sequence@@ repeatedIntInds]], 
						term
					],
					term
				]
			]
		]/. FlavorSum[_?NumericQ]-> 1
	];
	
	(* replace flavor deltas *)
	res/. Delta[n_Integer, m_Integer]:> If[n==m, 1, 0]
]


(* ::Subsubsection::Closed:: *)
(*Introduce component fields for flavor multiplets*)


(* ::Text:: *)
(*This function introduces component fields for all flavor multiplets appearing in the given Lagrangian.*)
(*It returns a list of the labels of the defined fields and the replacement rules to obtain them from the flavor multiplets.*)


IntroduceFlavorComponentFields[lag_]:=Module[
	{
		replRules={},
		flavorComponentFields={},
		(* get all fields of the Lagrangian *)
		presentFields=KeySelect[$FieldAssociation,!FreeQ[lag,#,All]&]
	}
	,
	(* loop over all fields *)
	Do[
		(* only consider flavored fields *)
		If[!FreeQ[presentFields[field,Indices],Alternatives@@Keys@GetFlavorIndices[]],
			Module[
				{
					fieldLabel,
					(* flavor index types of the fiels *)
					flavs=KeySelect[GetFlavorIndices[],MemberQ[presentFields[field,Indices],#]&],
					(* get properties of flavor multiplet parent field *)
					fieldProp=GetFields(*Updated*)[field],
					counters,counterLabels,
					pdgNumber,pdgName,pdgMass,pdgWidth,
					(* group indices *)
					groupInd
				}
				,
				(* loop over all possible flavor index combinations *)
				Do[
					(* determine courent flavor index combination *)
					counters=ToExpression["i"<>ToString[#]&/@Keys[flavs]];
					(* make to $-separated string *)
					counterLabels=StringJoin/@ToString/@Riffle[counters,"$"];
					(* create label for component field of the flavor multiplet field *)
					fieldLabel=ToExpression[ToString[field]<>counterLabels];
					AppendTo[flavorComponentFields,fieldLabel];
					(* determine pdg number and name if specified *)
					pdgNumber=If[KeyExistsQ[fieldProp[UFO],"pdg"],
						 "pdg"-> fieldProp[UFO,"pdg"][[Sequence@@counters]],
						 Nothing
					];
					pdgName=If[KeyExistsQ[fieldProp[UFO],"name"],
						"name"-> fieldProp[UFO,"name"][[Sequence@@counters]],
						Nothing
					];
					pdgMass=If[KeyExistsQ[fieldProp[UFO],"mass"],
						"mass"-> fieldProp[UFO,"mass"][[Sequence@@counters]],
						Nothing
					];
					pdgWidth=If[KeyExistsQ[fieldProp[UFO],"width"],
						"width"-> fieldProp[UFO,"width"][[Sequence@@counters]],
						Nothing
					];
					
					(* define the flavor component field *)
					DefineField[fieldLabel,fieldProp[Type],
						(* drop flavor idices *)
						Indices->DeleteCases[fieldProp[Indices],x_/;MemberQ[Keys@GetFlavorIndices[],x]],
						Charges->fieldProp[Charges],
						SelfConjugate->fieldProp[SelfConjugate],
						Chiral->fieldProp[Chiral],
						Mass->If[fieldProp[Mass]===0,
							0
							,
							{If[fieldProp[Heavy],Heavy,Light],ToExpression[ToString[fieldProp[Mass]]<>counterLabels]}
						],
						GoldstoneBoson->fieldProp[GoldstoneBoson],
						BackgroundField->fieldProp[BackgroundField],
						UFO-> Association@{pdgNumber,pdgName,pdgMass,pdgWidth}
						(*UFO->Switch[{pdgNumber,pdgName},
							{Nothing,Nothing},<||>,
							{_,Nothing},<|"pdg"->pdgNumber|>,
							{Nothing,_},<|"name"->pdgName|>,
							{_,_},<|"pdg"->pdgNumber,"name"->pdgName|>
						]*)
					];
					
					(* get group indices *)
					groupInd=Length@DeleteCases[fieldProp[Indices],x_/;MemberQ[Keys@GetFlavorIndices[],x]];
					groupInd=Table[ToExpression["i"<>ToString[j]],{j,groupInd}];
					
					(* replacement rules for the field and the bared field *)
					AppendTo[replRules,
						Rule[
							Field[Ext[field,sub_],     fieldProp[Type], Join[Map[Pattern[#,_]&,groupInd],counters], devs_],
							Field[Ext[fieldLabel,sub], fieldProp[Type], groupInd,                         devs]
						]
					];
					AppendTo[replRules,
						Rule[
							Field[Ext[(*OverBar*)Conj@field,sub_],     fieldProp[Type], Join[Map[Pattern[#,_]&,groupInd],counters], devs_],
							Field[Ext[(*OverBar*)Conj@fieldLabel,sub], fieldProp[Type], groupInd,                         devs]
						]
					];
					(* safe field names for later removal and create replacement rule *)
					,
					Evaluate[Sequence@@Table[{ToExpression["i"<>ToString[fl]],1,flavs[fl,IndexDimension]},{fl,Keys[flavs]}]]
				]
			]
		]
		,
		{field,Keys@presentFields}
	];
	
	{flavorComponentFields,replRules}
]


(* ::Subsubsection::Closed:: *)
(*Expand flavored couplings*)


(* ::Text:: *)
(*Unflavored couplings*)


ExpandCouplingFlavor[coup:Coupling[_,{},_]]:= coup


(* ::Text:: *)
(*Flavored particles*)


ExpandCouplingFlavor[coup:Coupling[l_,{__},n_]]:= Module[
	{
		(* get all flavor indices *)
		inds= $CouplingAssociation[l, Indices],
		dims,
		myList
	}
	,
	(* get dimension of all indices *)
	dims= $FlavorIndices[#, IndexDimension]&/@ inds;
	(* values all flavor indices can assume *)
	dims= Range/@ dims;
	(* all possible combinations of indices *)
	inds= Flatten[Outer[myList, Sequence@@ dims]]/. myList->List;
	(* get all coupling combinations *)
	Coupling[l,#,n]&/@ inds
]


(* ::Subsubsection::Closed:: *)
(*Non-redundant coupling parameterization*)


(* ::Text:: *)
(*Takes a gauge fixing parameter and returns a list of parameter information that can be used in parameter cards.*)


NonRedundantUFOParams[\[Xi]Gauge[f_],___]:= {{
	"xi"<>If[KeyExistsQ[GetFields[f, UFO], "name"], GetFields[f, UFO]["name"], ToString[f]], 
	1, (* default to Feynman gauge *)
	{"real", "external"}
}}


(* ::Text:: *)
(*Takes a coupling and provides a non-redundant parameterization compatible with symmetries and reality assumptions.*)
(*Returns a list of parameter information that can be used in parameter cards. Decomposes complex coefficients into real and imaginary parts.*)


(* this is needed, e.g., if neutrino masses are generated by the Weinberg operator after SSB but the fields have been defined massless *)
NonRedundantUFOParams[list_List, flavSym_:<||>]:= Flatten[NonRedundantUFOParams[#, flavSym]&/@ list, 1]


NonRedundantUFOParams[name_Symbol, flavSym_:<||>]:= Module[
		{additionalSymbols, baseName, conj, im, inds, newSymb, permutations, re, 
			set, symbols, syms, types, values, flavorStructures= flavSym, structureCounter= 1, struct, conjProps},
	(*Not a coupling*)
	If[!KeyExistsQ[$CouplingAssociation, name], 
		Return@ {{TranslateToUFO@ name, "0", {"real", "external"}}};
	];
	
	(*Properties of the coupling*)
	{conj, inds, syms}= Lookup[$CouplingAssociation@ name, {SelfConjugate, Indices, Symmetries}];
	baseName= UFONames@ name;
	(* special case: field masses are assumed diagonal even if defined non-diagonal *)
	If[Length@(*GetFieldsUpdatedByProperty*)GetFieldsByProperty[Mass-> name]>0,
		{conj, inds, syms}= {True, DeleteDuplicates@inds, <|Range[Length@DeleteDuplicates@inds]->1|>}
	];
	
	(*Return if scalar coupling*)
	If[Length@ inds === 0,
		If[$CouplingAssociation[name, SelfConjugate],
			Return@ {{baseName, "0", {"real", "external"}}};	
		,
			{re, im}= {"Re_"<> baseName, "Im_"<> baseName};
			Return@ {
				{baseName, TranslateToUFO[TempSymb@ re + I TempSymb@ im], {"complex", "internal"}},
				{re, "0", {"real", "external"}},
				{im, "0", {"real", "external"}}
			};
		];
	];
	
	(*If real, the conjugate is the same index permutation*)
	If[conj, conj= Range@ Length@ inds];
	inds= $FlavorIndices[#, IndexDimension]&/@ inds;
	values= ConstantArray[None, inds];
	types= ConstantArray[{}, inds];
	inds= Tuples[Range/@ inds];
	
	(* START WIP *)
	(* if a flavor structure/symmetry is defined *)
	additionalSymbols= Reap[If[KeyExistsQ[flavorStructures, name],
		(* get flavor structures for this coupling *)
		flavorStructures= flavorStructures[name];
		(* handle couplings forbidden by symmetry *)
		If[flavorStructures==={},
			(* no allowed flavor structures *)
			Do[
				Part[values, Sequence@@ set]= "ZERO"; (* remove this coupling *)
				Part[types,  Sequence@@ set]= {"real", "internal"}; (* the couplings might not be real, but since they are removed, it does not matter *)
				,
				{set, inds}
			]
			,
			(* remove all entries that are forbidden *)
			Do[
				Part[values, Sequence@@ set]= "ZERO"; (* remove this coupling *)
				Part[types,  Sequence@@ set]= {"real", "internal"}; (* the couplings might not be real, but since they are removed, it does not matter *)
				,
				{set, Complement[inds, Sequence@@ (Keys/@ (Association@ DeleteCases[ArrayRules@ #, HoldPattern[_-> 0]]&/@ Flatten[flavorStructures,1]))]}
			];
			(* loop over allowed flavor structures and populate non-empty entries *)
			structureCounter= 1;
			Do[
				(* convert SparseArray to association w/o 0 entries *)
				struct= Association@ DeleteCases[ArrayRules@ Total@ structSA, HoldPattern[_-> 0]];
				
				(*Wrap the string in a dymmy head to allow UFO translation of formula*)
				
				(* determine label for coupling: use only ONE of the following options *)
				(* OPTION 1: find index combination unique to this structure *)
				(*
				newSymb= Reap@ Extract[
					baseName,
					First[
						Complement[Keys@ First[structSA], Sequence@@ (Keys[First[#]]&/@ flavorStructures)],
						Sow[1]; First@ Keys@ First[structSA]
					]
				];
				If[Last[newSymb]==={},
					newSymb= First[newSymb],
					newSymb= First[newSymb]<>"_"<>ToString[structureCounter++] (* if there is no unique index combination add a counter for unique name *)
				];
				*)
				(* OPTION 2: simply take the lowest index combination and always add a counter *)
				newSymb= Extract[baseName, First@ Keys@ struct]<>"_"<>ToString[structureCounter++];
				(*newSymb= If[Length[flavorStructures]>1,
					Extract[baseName, First@ Keys@ struct]<>"_"<>ToString[structureCounter++],
					Extract[baseName, First@ Keys@ struct] (* this leads to self-referenced parameters *)
				];*)
				
				(* OPTION 3: [BROKEN] ignore indices and simply take the bare coupling label and append counter [this ignores UFO\[Rule]\[LeftAssociation]"name"\[Rule] \[Ellipsis]\[RightAssociation]]] *)
				(*
				newSymb= TranslateToUFO[name]<>"_"<>ToString[structureCounter++];
				*)
				
				(* distinguish intrinsically real and complex couplings *)
				If[conj === False,
					(* complex couplings *)
					{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
					Do[
						Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/.None-> 0) + struct@ p* (TempSymb@ re + I TempSymb@ im);
						Part[types,  Sequence@@ p]= {"complex", "internal"};
						,
						{p, Keys@ struct}
					];
					(*Ensure that unitary matrices default to identity matrix*)
					If[$CouplingAssociation[name, Unitary] && SameQ@@ First[Keys@ struct], 
						Sow@ {re, "1", {"real", "external"}};
						,
						Sow@ {re, "0", {"real", "external"}};
					];
					Sow@ {im, "0", {"real", "external"}};
					,
					(* not intrinsically complex couplings can have real/imaginary/complex constituents *)
					(* determine whether the coupling associated to this flavor structure is Re, Im, or ReIm *)
					conjProps= DetermineConjugationPropertiesOfFlavorStructure[Total@ structSA, conj, syms];
					Switch[First[conjProps],
						(* real *)
						Re,
							Do[
								Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/. None-> 0) + struct@ p* (TempSymb@ newSymb);
								Part[types,  Sequence@@ p]= {"real", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {newSymb, "0", {"real", "external"}},
						(* imaginary *)
						Im, 
							im= "Im_"<> newSymb;
							Do[
								Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/. None-> 0) + struct@ p* (TempSymb@ im); (* no factor of \[ImaginaryI] needed since struct[p] contains it *)
								Part[types,  Sequence@@ p]= {"complex", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {im, "0", {"real", "external"}},
						(* complex *)
						ReIm, 
							If[Length[structSA]=!=2, Message[DetermineConjugationPropertiesOfFlavorStructure::missconj, structSA, name]; Abort[]];
							{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
							Do[
								If[MemberQ[Last[structSA]["NonzeroPositions"], p] (*OLD->*)(*MemberQ[conjProps[[2,2]],p]*), 
									Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/. None-> 0) + struct@ p* (TempSymb@ re -  I TempSymb@ im),
									Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/. None-> 0) + struct@ p* (TempSymb@ re +  I TempSymb@ im)
								];
								
								Part[types,  Sequence@@ p]= {"complex", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {re, "0", {"real", "external"}};
							Sow@ {im, "0", {"real", "external"}}
					]
					(*
					Switch[structSA,
						(* real *)
						Transpose[structSA, conj],
							Do[
								Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/.None-> 0) + struct@ p* (TempSymb@ newSymb);
								Part[types, Sequence@@ p]= {"real", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {newSymb, "0", {"real", "external"}};
							,
						(* complex selfconjugate *)
						ConjugateTranspose[structSA, conj],
							{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
							Do[
								Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/.None-> 0) + struct@ p* (TempSymb@ re + (*this sign should flip if struct[p] contains an \[ImaginaryI]*) I TempSymb@ im);
								Part[types, Sequence@@ p]= {"complex", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {re, "0", {"real", "external"}};
							Sow@ {im, "0", {"real", "external"}};
							,
						(* should never appear *)
						_,
							Message[DefaultParamCard::unknownreim, structSA, name];
							{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
							Do[
								Part[values, Sequence@@ p]= (Part[values, Sequence@@ p]/.None-> 0) + struct@ p* (TempSymb@ re + (*this sign should flip if struct[p] contains an \[ImaginaryI]*) I TempSymb@ im);
								Part[types, Sequence@@ p]= {"complex", "internal"};
								,
								{p, Keys@ struct}
							];
							Sow@ {re, "0", {"real", "external"}};
							Sow@ {im, "0", {"real", "external"}};
					]
					*)
				]
				,
				{structSA, flavorStructures}
			]
		]
		,
		(* when no flavor structures have been defined for this coupling *)
		Do[
			(*Skip entry if already populated*)
			If[Extract[values, set] =!= None, Continue[];];
			(*Wrap the string in a dymmy head to allow UFO translation of formula*)
			newSymb= Extract[baseName, set];
			
			(*Get all signed permutations of the current indices from symmetries*)
			permutations= Merge[KeyValueMap[set[[#1]] -> #2&, syms], CombinedSym];
			Which[
				MemberQ[permutations, 0], 
					(*entries are automatically zero if a permutation has both signs: entry= - entry*)
					Do[
						Part[values, Sequence@@ p]= "ZERO";
						Part[types, Sequence@@ p]= {"real", "internal"};
					, {p, Keys@ permutations}];
				,conj === False,
					(*Coupling is inherently complex*)
					{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
					Do[
						Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re + I TempSymb@ im);
						Part[types, Sequence@@ p]= {"complex", "internal"};
					, {p, Keys@ permutations}];
					(*Ensure that unitary matrices default to identity matrix*)
					If[$CouplingAssociation[name, Unitary] && SameQ@@ set, 
						Sow@ {re, "1", {"real", "external"}};
					,
						Sow@ {re, "0", {"real", "external"}};
					];
					Sow@ {im, "0", {"real", "external"}};
				,True,
					Switch[permutations@ set[[conj]]
					, _Missing, 
						(*The conjugated entry is different from the symmetry-related entries*)
						{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
						Do[
							Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re + I TempSymb@ im);
							Part[types, Sequence@@ p]= {"complex", "internal"};
						, {p, Keys@ permutations}];
						permutations= Merge[KeyValueMap[set[[conj]][[#1]] -> #2&, syms], CombinedSym];
						Do[
							Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re - I TempSymb@ im);
							Part[types, Sequence@@ p]= {"complex", "internal"};
						, {p, Keys@ permutations}];
						Sow@ {re, "0", {"real", "external"}};
						Sow@ {im, "0", {"real", "external"}};
					, +1,
						(*The conjugated entry equals a symmetric permutation*)
						Do[
							Part[values, Sequence@@ p]= permutations@ p* TempSymb@ newSymb;
							Part[types, Sequence@@ p]= {"real", "internal"};
						, {p, Keys@ permutations}];
						Part[values, Sequence@@ set]= 0;
						Part[types, Sequence@@ set]= {"real", "external"};
					, -1, 
						(*The conjugated entry equals an antisymmetric permutation*)
						im= "Im_"<> newSymb;
						Do[
							Part[values, Sequence@@ p]= permutations@ p* I * TempSymb@ im;
							Part[types, Sequence@@ p]= {"complex", "internal"};
						, {p, Keys@ permutations}];
						Sow@ {im, "0", {"real", "external"}};
					, _,
						(*Should never reach*)
						Abort[]
					];
			];
			,
			{set, inds}
		]
	]][[2]];
	(* END WIP *)
	
	(*
	(*Loop over indices in the tensor array*)
	additionalSymbols= Reap[Do[
			(*Skip entry if already populated*)
			If[Extract[values, set] =!= None, Continue[];];
			(*Wrap the string in a dymmy head to allow UFO translation of formula*)
			newSymb= Extract[baseName, set];
			
			(*Get all signed permutations of the current indices from symmetries*)
			permutations= Merge[KeyValueMap[set[[#1]] -> #2&, syms], CombinedSym];
			Which[
			MemberQ[permutations, 0], 
				(*entries are automatically zero if a permutation has both signs: entry= - entry*)
				Do[
					Part[values, Sequence@@ p]= 0; (* replace by "ZERO" ? *)
					Part[types, Sequence@@ p]= {"real", "internal"};
				, {p, Keys@ permutations}];
			,conj === False,
				(*Coupling is inherently complex*)
				{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
				Do[
					Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re + I TempSymb@ im);
					Part[types, Sequence@@ p]= {"complex", "internal"};
				, {p, Keys@ permutations}];
				(*Ensure that unitary matrices default to identity matrix*)
				If[$CouplingAssociation[name, Unitary] && SameQ@@ set, 
					Sow@ {re, "1", {"real", "external"}};
				,
					Sow@ {re, "0", {"real", "external"}};
				];
				Sow@ {im, "0", {"real", "external"}};
			,True,
				Switch[permutations@ set[[conj]]
				, _Missing, 
					(*The conjugated entry is different from the symmetry-related entries*)
					{re, im}= {"Re_"<> newSymb, "Im_"<> newSymb};
					Do[
						Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re + I TempSymb@ im);
						Part[types, Sequence@@ p]= {"complex", "internal"};
					, {p, Keys@ permutations}];
					permutations= Merge[KeyValueMap[set[[conj]][[#1]] -> #2&, syms], CombinedSym];
					Do[
						Part[values, Sequence@@ p]= permutations@ p* (TempSymb@ re - I TempSymb@ im);
						Part[types, Sequence@@ p]= {"complex", "internal"};
					, {p, Keys@ permutations}];
					Sow@ {re, "0", {"real", "external"}};
					Sow@ {im, "0", {"real", "external"}};
				, +1,
					(*The conjugated entry equals a symmetric permutation*)
					Do[
						Part[values, Sequence@@ p]= permutations@ p* TempSymb@ newSymb;
						Part[types, Sequence@@ p]= {"real", "internal"};
					, {p, Keys@ permutations}];
					Part[values, Sequence@@ set]= 0;
					Part[types, Sequence@@ set]= {"real", "external"};
				, -1, 
					(*The conjugated entry equals an antisymmetric permutation*)
					im= "Im_"<> newSymb;
					Do[
						Part[values, Sequence@@ p]= permutations@ p* I * TempSymb@ im;
						Part[types, Sequence@@ p]= {"complex", "internal"};
					, {p, Keys@ permutations}];
					Sow@ {im, "0", {"real", "external"}};
				, _,
					(*Should never reach*)
					Abort[]
				];
			];
		, {set, inds}]][[2]];
	*)
	
	Join[
		Table[
			{Extract[baseName, set], TranslateToUFO@ Extract[values, set], Extract[types, set]}
			,
			{set, inds}
		], 
		Flatten[additionalSymbols, 1]
	]
]


(* ::Text:: *)
(*Converts mixed lists to 0, and pure +/-1 lists to single number*)


CombinedSym[signs_List]:= Switch[signs
	,{1..}, 1
	,{-1..}, -1
	,_, 0
];


(* ::Text:: *)
(*Determine the UFO name of a coupling or produces an array of UFO names for every single index in a coupling.*)


UFONames[name_Symbol]:= Module[{array, field, inds},
	(*Check if a field is associated with a mass *)
	field= (*GetFieldsUpdatedByProperty*)GetFieldsByProperty[Mass-> name];
	Which[
	(*A UFO name was provided*)
	KeyExistsQ[$CouplingAssociation[name, UFO], "name"],
		$CouplingAssociation[name, UFO, "name"]
	(*Masses can take name from the associated field*)
	,Length@ field> 0 && KeyExistsQ[GetFields(*Updated*)[First@ field, UFO], "name"],
		"M"<> #&/@ GetFields(*Updated*)[First@ field, UFO]["name"]
	,True,
		inds= $CouplingAssociation[name, Indices];
		If[Length@ inds === 0,
			Return@ TranslateToUFO[name];
		];
		
		inds= $FlavorIndices[#, IndexDimension]&/@ inds;
		array= ConstantArray[0, inds];
		inds= Tuples[Range/@ inds];
		Do[
			Part[array, Sequence@@ ind]= TranslateToUFO[name] <> "_" <> StringJoin@@ ToString/@ ind; 
		, {ind, inds}];
		array		
	]
];



(* ::Text:: *)
(*Export a python file that creates a restriction file for MadGraph*)


WriteRestrictionFile[parametersFile_String, restrictedParams_List, outputDirectory_String, flavSym_:<||>] := Module[{card, params, value, parametersExternal, parametersNeeded, dir},
	card= Import[parametersFile]/. {("lhacode" -> {x_}) :> ("lhacode" -> x), List -> Association};
	params= First/@ Join@@(NonRedundantUFOParams[#, flavSym]&/@ restrictedParams);
	Do[
		If[KeyExistsQ[card["parameters"], param],
			value= card["parameters", param, "value"];
			If[value == "0", card["parameters", param, "value"] = "1e-99"];
		];
		,
		{param, params}
	];

	parametersExternal= card/. Association->List;
	parametersNeeded= DeleteCases[
		Keys[card["parameters"]],
		x_/;(StringMatchQ[x, "Re_"~~___] || StringMatchQ[x, "Im_"~~___])
	];

	dir= CreateNewDir[outputDirectory];
	SetDirectory[dir]; (* ResetDirectory missing, maybe better specify the full path *)
	WritePyParameters[parametersNeeded, parametersExternal, False];
];


RemoveColorRedundancies[vertexAssoc_] := Module[
	{
		lorentz = vertexAssoc["Lorentz"],
		color = vertexAssoc[ToString@$QCD$UFO],
		lorentzUnique,
		colorUnique,
		mappingLorentz,
		mappingColor,
		result
	},

	colorUnique = DeleteDuplicates@color;

	mappingColor = Lookup[AssociationThread[colorUnique -> Range[Length[colorUnique]]], color] - 1;

	result = Join[vertexAssoc, <|"Colormap" -> mappingColor|>];
	result[ToString@$QCD$UFO] = colorUnique;

	result
]


(* ::Text:: *)
(*Add Lorentz structures with the same coupling and color structure, excluding multifermion rules*)


AddRedundantLorentzStructures[vertexAssoc_] := Module[
	{
		result = vertexAssoc,
		structures,
		keys
	},
	(*Apply only if vertex contains less than 4 fermions*)
	If[Count[Characters@vertexAssoc["FieldsContent"], "F"] < 4,
		(*Get relevant keys from the vertex association*)
		keys = Complement[Keys@vertexAssoc, {"FieldsContent", "Fields", "FieldsSorted", "Feynmanrule", "Lorentz"}];
		keys = Join[keys, {"Lorentz"}];
		(*Get corresponding vertex structures*)
		structures = Transpose[Lookup[vertexAssoc, keys]];
		(*Gather by same color structures + couplings and sum*)
		structures = Map[Join[Most@First[#], {Total[#[[All, -1]]]}]&, GatherBy[structures, Most[#]&]];
		(*Update the association wih the summed structures*)
		result = ReplacePart[
			result,
			AssociationThread[Key /@ keys -> Transpose@structures]
		];
	];

	result
]


(* ::Text:: *)
(**)


(*Contracts Flavor deltas with external legs in Feynman rules.
Keeps always the first of the indices in a Flavor delta to have the correct index names for ExplicitEinsteinSums*)
ContractFlavorDelta[expr_]:= Module[{result=LagrangianExpand@expr, flavors=GetFlavorIndices[], indicesDelta, term},
	If[Head@result === Plus,
		Map[
			Function[
				{arg},
				term = arg;
				indicesDelta = Cases[term, Delta[Index[a_,rep: Alternatives@@ Keys@ flavors],Index[b_,rep: Alternatives@@ Keys@ flavors]]:>{Index[a,rep],Index[b,rep]},All];
				indicesEXT = Cases[result, f_Field :> f[[3]], All];
				Do[
					(*First index is the field index*)
					If[FreeQ[indicesEXT, First@ind] || FreeQ[indicesEXT, Last@ind],
						term = (term /. Delta[First@ind,Last@ind]->1)/.Last@ind->First@ind
					];
					,
					{ind,indicesDelta}
				];
				term
			]
			,
			result
		]
		,
		indicesDelta = Cases[result, Delta[Index[a_,rep: Alternatives@@ Keys@ flavors],Index[b_,rep: Alternatives@@ Keys@ flavors]]:>{Index[a,rep],Index[b,rep]},All];
		indicesEXT = Cases[result, f_Field :> f[[3]], All];
		Do[
			If[FreeQ[indicesEXT, First@ind] || FreeQ[indicesEXT, Last@ind],
				result = (result /. Delta[First@ind,Last@ind]->1)/.Last@ind->First@ind
			];
			,
			{ind,indicesDelta}
		];
		result
	]
];
