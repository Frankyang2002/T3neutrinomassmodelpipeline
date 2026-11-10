(* ::Package:: *)

Package["Matchete`"]


(* ::Title:: *)
(*Matchete`Definitions`*)


(* ::Subtitle:: *)
(*Field and coupling definitions for Matchete paclet.*)


(* ::Chapter:: *)
(*Public:*)


(* ::Section:: *)
(*Scoping*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


PackageExport["Field"]
PackageExport["DefineField"]
PackageExport["RemoveField"]
PackageExport["ResetFields"]
PackageExport["GetFields"]
PackageExport["GetFieldsByProperty"]


PackageExport["FieldGenerators"]


PackageExport["Scalar"]
PackageExport["Fermion"]
PackageExport["Vector"]
PackageExport["Ghost"]
PackageExport["AntiGhost"]


PackageExport["Indices"]
PackageExport["Charges"]
PackageExport["SelfConjugate"]
PackageExport["Mass"]
PackageExport["Heavy"]
PackageExport["Light"]
PackageExport["Type"]
PackageExport["Chiral"]
PackageExport["LeftHanded"]
PackageExport["RightHanded"]
PackageExport["GoldstoneBoson"]
PackageExport["BackgroundField"]
PackageExport["GhostFields"]
PackageExport["Unitary"]


PackageExport["Symmetries"]
PackageExport["SymmetricPermutation"]
PackageExport["AntisymmetricPermutation"]


PackageExport["Coupling"]
PackageExport["DefineCoupling"]
PackageExport["RemoveCoupling"]
PackageExport["GetCouplings"]
PackageExport["ResetCouplings"]
PackageExport["DiagonalCoupling"]


PackageExport["GetFlavorIndices"]
PackageExport["IndexDimension"]
PackageExport["DefineFlavorIndex"]
PackageExport["RemoveFlavorIndex"]
PackageExport["ResetFlavorIndices"]


PackageExport["IndexAlphabet"]


PackageExport["Abelian"]
PackageExport["Representations"]
PackageExport["ClebschGordanCoefficients"]
PackageExport["ClebschGordanIndices"]


PackageExport["U1"]
PackageExport["fund"]
PackageExport["tFundf"]
PackageExport["adj"]
PackageExport["gen"]
PackageExport["fStruct"]
PackageExport["eps"]
PackageExport["dSym"]


PackageExport["FieldStrength"]
PackageExport["FS"]


PackageExport["GetGaugeGroups"]
PackageExport["FundAlphabet"]
PackageExport["AdjAlphabet"]
PackageExport["DefineGaugeGroup"]
PackageExport["RemoveGaugeGroup"]
PackageExport["ResetGaugeGroups"]


PackageExport["GetGlobalGroups"]
PackageExport["DefineGlobalGroup"]
PackageExport["RemoveGlobalGroup"]
PackageExport["ResetGlobalGroups"]


PackageExport["ResetAll"]


PackageExport["PlusHc"]


PackageExport["FreeLag"]
PackageExport["LoadModel"]
PackageExport["ParentModel"]
PackageExport["ParameterDefault"]
PackageExport["SetModelOptions"]
PackageExport["GetModels"]
PackageExport["ModelParameters"]


(* ::Subsubsection::Closed:: *)
(*Internal*)


PackageScope["SymmetryOverride"]
PackageScope["CouplingAssociationEntry"]


PackageScope["GroupFromRep"]
PackageScope["GroupFromInd"]
PackageScope["RepFromInd"]
PackageScope["FieldsFromGroup"]
PackageScope["CouplingsFromGroup"]


PackageScope["$FieldAssociation"]
PackageScope["$GhostFields"]
PackageScope["FieldByProperty"]
PackageScope["$CouplingAssociation"]
PackageScope["$FlavorIndices"]
PackageScope["$GaugeGroups"]
PackageScope["GaugeGroupByProperty"]
PackageScope["$GlobalGroups"]
PackageScope["GetCouplingByProperty"]


PackageScope["KinOpLagrangian"]


PackageScope["$DropDiagonalCouplings"]


PackageScope["ReadFileExpression"]


PackageScope["SortOptions"]


(* ::Section:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Exported*)


Field::usage               = "Field[label,type,{indices},{CDerivs}] denotes a field with given label, type and list of {indices}. The argument {CDerivs} is the list of Lorentz indices of the covariant derivatives acting on the field.";
DefineField::usage         = "DefineField[label, type] defines a new field of given type and label. The type can be Scalar,Fermion or Vector. The default options for this function are: Indices->{}, Charges->{}, SelfConjugate->False, Chiral -> False, Mass->Heavy, NiceForm-> Default, GhostFields-> Automatic, UFO -> <||>.";
RemoveField::usage         = "RemoveField[alias] undefines the field with label alias together with its corresponding mass coupling.";
ResetFields::usage         = "Undefines all fields (excluding gauge fields) and field masses."
GetFields::usage           = "GetFields[] returns an association of all fields that have been defined by DefineField and their properties. GetFields[field] returns an association with all the properties of an already defined field.";
GetFieldsByProperty::usage = "GetFieldsByProperty[properties] returns all fields matching the property, list of properties or sequence of properties determined by the argument.";


FieldGenerators::usage      = "FieldGenerators[Field,GaugeGroup,Indices] returns the generator of the given GaugeGroup in the same representation as the given Field, with given Indices."


Scalar::usage  =    "Field[\[Phi],Scalar,{indices},{CDerivs}] denotes a scalar field \[Phi].";
Vector::usage  =    "Field[A,Vector[Index[\[Mu],Lorentz]],{indices},{CDerivs}] denotes a vector field \!\(\*SubscriptBox[\(A\), \(\[Mu]\)]\).";
Fermion::usage =    "Field[\[Psi],Fermion,{indices},{CDerivs}] denotes a fermion field \[Psi].";


Indices::usage       = "Indices is an option for the routine DefineField that specifies a list of representations under which the field transforms. By default, Indices is an empty list {}.";
Charges::usage       = "Charges is an option for the routine DefineField that specifies a list of U(1) groups and corresponding charges of the field. By default, Charges is an empty list {}.";
SelfConjugate::usage = "SelfConjugate is an option for the routine DefineField that specifies whether the field is self-conjugated. By default, SelfConjugate is False.";
Mass::usage          = "Mass is an option for the routine DefineField that specifies the mass and EFT scaling of the field. The allowed values are either Heavy|Light|0, {Light,0}, {Heavy|Light,mlabel}, or {Heavy|Light,mlabel,{Indices}} where mlabel is a custom label for the mass and Indices must be a subset of the field indices. The default value is Mass->Heavy. If no mlabel is given, the default value Mflabel, with flabel being the label of the corresponding field, is used.";
Heavy::usage         = "Heavy specifies that the mass of a field is at the UV scale.";
Light::usage         = "Light specifies that the mass of a field is at the IR scale.";
Type::usage          = "Type specifies the type of a field. It can be Scalar, Vector or Fermion.";
Chiral::usage        = "Chiral is an option for the routine DefineField that specifies if a fermion field is LeftHanded or RightHanded.";
LeftHanded::usage    = "LeftHanded specifies that a field is of left-handed chirality."
RightHanded::usage   = "RightHanded specifies that a field is of right-handed chirality."
GoldstoneBoson::usage= "GoldstoneBoson indicates if the field is a GB associated with a heavy vector field."
BackgroundField::usage= "Is an option used in Matchete."
Unitary::usage       = "Is an option used in Matchete."
GhostFields::usage = "Is an option used in Matchete."


Symmetries::usage               = "Symmetries is an option that takes a list of symmetries. The possible option values are: SymmetricIndices[n1, n2, ...], AntisymmetricIndices[n1, n2, ...], SymmetricPermutation[n1, n2, ...] or AntisymmetricPermutation[n1, n2, ...]." ;
SymmetricPermutation::usage     = "SymmetricPermutation is a possible value for the option Symmetries. SymmetricPermutation[n1, n2, ...] denotes symmetry under a particular permutation of indices." ;
AntisymmetricPermutation::usage = "AntisymmetricPermutation is a possible value for the option Symmetries. AntisymmetricPermutation[n1, n2, ...] denotes antisymmetry under a particular permutation of indices." ;


Coupling::usage         = "Coupling[label,{indices},EFTOrder] denotes the coupling label with a set of indices and an EFTOrder."
DefineCoupling::usage   = "DefineCoupling[label] defines a new coupling with a given label. The default options for this function are: EFTOrder -> 0, Indices->{}, SelfConjugate->False, Symmetries-> {}, DiagonalCoupling -> False, , NiceForm-> Default.";
RemoveCoupling::usage   = "RemoveCoupling[alias] undefines the coupling with label alias.";
GetCouplings::usage     = "GetCouplings[] returns an association of all couplings that have been defined by DefineCoupling and their properties.GetCouplings[coupling] returns an association with all the properties of an already defined coupling.";
ResetCouplings::usage   = "Undefines all couplings (excluding the field masses and the gauge couplings)."
DiagonalCoupling::usage = "DiagonalCoupling is an option for the routine DefineCoupling. Allowed values are a list of booleans or an empty list. The n-th boolean specifies whether the n-th flavor index is a diagonal in flavor space. By default all flavor indices are assumed to be non-diagonal.
DiagonalCoupling is also the corresponding key in the association returned by GetCouplings[]."


GetFlavorIndices::usage   = "GetFlavorIndices[] returns an association of all flavor indices that have been defined by DefineFlavorIndex and their dimension. GetFlavorIndices[flavor] returns an association with the dimension of an already defined flavor index.";
IndexDimension::usage     = "IndexDimension is a key for the association returned by GetFlavorIndices.";
DefineFlavorIndex::usage  = "DefineFlavorIndex[IndexName, Dimension] defines a flavor index with given IndexName and Dimension. This routine has one optional argument, IndexAlphabet, to define the printing alphabet of the representation.";
RemoveFlavorIndex::usage  = "RemoveFlavorIndex[alias] undefines the flavor index with given alias.";
ResetFlavorIndices::usage = "Undefines all flavor indices.";


IndexAlphabet::usage             = "IndexAlphabet is an option for the routine DefineRepresentation and DefineFlavorIndex that specifies a printing alphabet for the representation or flavor indices.";


Abelian::usage                   = "Abelian is a key for the association returned by GetGaugeGroups.";
Representations::usage           = "Representations is a key for the association returned by GetGaugeGroups.";
ClebschGordanCoefficients::usage = "ClebschGordanCoefficients is a key for the association returned by GetGaugeGroups.";
ClebschGordanIndices::usage      = "ClebschGordanIndices[CGname] return the representation of the indices (ordered) for the Clebsch-Gordan coefficient CGname."


U1::usage      = "Standard name for Abelian U(1) group.";
fund::usage    = "Fundamental representations of groups are referred to by group@ fund."
tFundf::usage  = "Standard head for the Clebsch-Gordan coefficient, CG[tFundf@ group, {i, j, A, B}] associated with the commonly occuring combination \!\(\*SubscriptBox[SuperscriptBox[\(T\), \(Ci\)], \(j\)]\)\!\(\*SuperscriptBox[\(f\), \(CAB\)]\), where T is the generator of the fundamental representation of the group.";
adj::usage     = "Adjoint representations of groups are refered to by group@ adj."
gen::usage     = "Standard head for generator Clebsch-Gordan coefficients. E.g., CG[gen@ rep, {A, i, j}] for \!\(\*SubscriptBox[SuperscriptBox[\(T\), \(Ai\)], \(j\)]\), where A is an adjoint index and i an index of the representation rep.";
fStruct::usage = "Standard head for structure constant Clebsch-Gordan coefficients. E.g., CG[fStruct@ group, {A, B, C}] for \!\(\*SuperscriptBox[\(f\), \(ABC\)]\)";
eps::usage     = "Standard head for Levi-Civita type Clebsch-Gordan coefficients for SU(N) groups. E.g., CG[eps@ group, {a, b, ...}] for \!\(\*SuperscriptBox[\(\[Epsilon]\), \(ab ... \)]\), with indices in the fundamental representation.";
dSym::usage    = "Standard head for the totally symmetric d Clebsch-Gordan coefficients for SU(N\[GreaterEqual]3) groups. E.g., CG[dSym@ group, {A, B, C}] for \!\(\*SuperscriptBox[\(d\), \(ABC\)]\).";


FieldStrength::usage = "FieldStrength[V, {Index[\[Mu],Lorentz], Index[\[Nu],Lorentz]}, {indices}, {CDerivs}] denotes the field-strength tensor \!\(\*SubscriptBox[\(V\), \(\[Mu]\[Nu]\)]\) of the vector field V."
FS::usage            = "FS[V,\[Mu],\[Nu],indices] denotes the field-strength tensor of an already defined vector field V with Lorentz indices \[Mu] and \[Nu], and gauge indices as the last argument."


GetGaugeGroups::usage   = "GetGaugeGroups[] returns an association of all gauge groups that have been defined by DefineGaugeGroup and their properties. GetGaugeGroups[GaugeGroup] returns an association with all the properties of an already defined gauge group.";
FundAlphabet::usage     = "FundAlphabet is an option for the routine DefineGaugeGroup that specifies a printing alphabet for the indices of the fundamental representation.";
AdjAlphabet::usage      = "AdjAlphabet is an option for the routine DefineGaugeGroup that specifies a printing alphabet for the indices of the adjoint representation.";
DefineGaugeGroup::usage = "DefineGaugeGroup[groupName, LieGroup, coupling, fieldName] defines a gauge group labelled groupName from a given simple LieGroup, together with its corresponding coupling and gauge bosons labelled fieldName. DefineGaugeGroup also defines the fundamental and adjoint representations of the corresponding group, with some basic Clebsch-Gordan coefficients. This routine has two optional arguments, FundAlphabet and AdjAlphabet, to define respectively the printing alphabet of the fundamental and adjoint representations.";
RemoveGaugeGroup::usage = "RemoveGaugeGroup[alias] undefines the gauge group with given alias, together with its associated representations, Clebsch-Gordan coefficients, gauge coupling and bosons.";
ResetGaugeGroups::usage = "Undefines all gauge groups and associated bosons, couplings, representations, and Clebsch-Gordan coefficients."


GetGlobalGroups::usage   = "GetGlobalGroups[] returns an association of all global groups that have been defined by DefineGlobalGroup and their properties. GetGlobalGroups[GlobalGroup] returns an association with all the properties of an already defined global group.";
DefineGlobalGroup::usage = "DefineGlobalGroup[groupName, LieGroup] defines a global group labelled groupName from a given simple LieGroup. DefineGlobalGroup also defines the fundamental and adjoint representations of the corresponding group, together with some basic Clebsch-Gordan coefficients. This routine has two optional arguments, FundAlphabet and AdjAlphabet, to define respectively the printing alphabet of the fundamental and adjoint representations.";
RemoveGlobalGroup::usage = "RemoveGlobalGroup[alias] undefines the global group with given alias together with its associated Clebsch-Gordan coefficients.";
ResetGlobalGroups::usage = "Undefines all global groups and associated representations, and Clebsch-Gordan coefficients."


ResetAll::usage = "Undefines all fields, couplings, flavor indices, global and gauge groups."


PlusHc::usage = "PlusHc[expr] returns expr + Bar[expr] and can be used for convenience when writing the Lagrangian."


FreeLag::usage         = "FreeLag[f] returns the free Lagrangian for the fields specified by the argument f, which can be either a single label or a sequence of labels. If the function is called without arguments, it returns the free Lagrangian of all fields which have been defined.";
LoadModel::usage       = "LoadModel[ModelName] loads the model 'ModelName' previously defined in the Models folder in Matchete. The optional arguments lets the user change the name/value of any of the model parameters or add/modify options of fields, couplings, or groups.";
GetModels::usage       = "Returns an association with all available model files that can be used with the LoadModel[...] command.";
ModelParameters::usage = "ModelParameters is an option of the LoadModel routine."
SetModelOptions::usage = "SetModelOptions is an option of the LoadModel routine."
ParentModel::usage     = "ParentModel[<model file>] is used to indicate that a model file builds directly on another."
ParameterDefault::usage= "ParameterDefault[<parameter> -> <default value>] is used to indicate a default value for a parameter in a model file if no other value is given during loading."


(* ::Subsubsection::Closed:: *)
(*Internal*)


$FieldAssociation::usage    = "$FieldAssociation is an association containing the complete information on all defined fields."
$CouplingAssociation::usage = "$CouplingAssociation is an association containing the complete information on all defined couplings."
$FlavorIndices::usage       = "$FlavorIndices is an association containing the complete information on all defined flavor indices."
$GaugeGroups::usage         = "$GaugeGroups is an association containing information on all defined gauge groups. It does not contain neither representations nor Clebsch-Gordan coefficients associated to the group."
$GlobalGroups::usage        = "$GlobalGroups is an association containing information on all defined global groups. It does not contain neither representations nor Clebsch-Gordan coefficients associated to the group."


KinOpLagrangian::usage      = "KinOpLagrangian[field] returns the  the part of the free Lagrangian that reproduces the kinetic piece of the fluctuation operator (see [2012.08506] eq. (2.3))."


GroupFromRep::usage         = "Returns the name of the group associtated to a given representation."
FieldsFromGroup::usage      = "Returns a list of fields that are charged under a given group."
CouplingsFromGroup::usage   = "Returns a list of couplings that are charged under a given group."


Ghost::usage   =    "Field[c,Ghost,{indices},{CDerivs}] denotes a ghost field c.";
AntiGhost::usage =    "Field[c,AntiGhost,{indices},{CDerivs}] denotes an antighost field c.";


(* ::Chapter:: *)
(*Private:*)


(* ::Subsubsection::Closed:: *)
(*Sorting of Options before caching *)


SortOptions::usage= "SortOptions[func, opts] sorts the options (opts) for a given function (func) and includes the missing options with default values."


SortOptions[func_, opts___]:= Module[
	{op={opts}}
	,
	Sequence@@Sort@Join[
		FilterRules[Options@func,Except[op[[;;,1]]]],
		op
	]
]


(* ::Section:: *)
(*Field definitions*)


(* ::Subsubsection::Closed:: *)
(*Association of all fields*)


$FieldAssociation=<||>;
$GhostFields=<||>;


FieldByProperty[]:= $FieldAssociation
FieldByProperty[propsSeq__]:= FieldByProperty[List[propsSeq]]
FieldByProperty[props_Association]:= FieldByProperty@ Normal@ props;
FieldByProperty[props:_List|_Rule]:=
	Keys@ Select[$FieldAssociation, MatchQ[#, KeyValuePattern[props]]&]


GetFields[FieldName___]:= 
KeyDrop[$FieldAssociation,FieldByProperty[Type->Ghost|AntiGhost]][FieldName]


GetFieldsByProperty[]:=GetFields[]

GetFieldsByProperty[propsSeq__]:=GetFieldsByProperty[List[propsSeq]]

GetFieldsByProperty[props_Association]:= GetFieldsByProperty@ Normal@ props

GetFieldsByProperty[props:_List|_Rule]:=
	Keys@ Select[GetFields[], MatchQ[#, KeyValuePattern[props]]&]


(* ::Subsection:: *)
(*Defining new fields*)


(* ::Subsubsection::Closed:: *)
(*Field properties*)


Field[0, __]= 0;


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineField::bkgMass                = "Background fields do not have a mass. They are Formally taken to have mass 0. Mass\[Rule] 0 is the only valid option for background fields.";
DefineField::FieldLabel             = "The label '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another label.";
DefineField::FieldType              = "The field type '`1`' is invalid. Please use either Scalar, Fermion or Vector.";
DefineField::Heavy0                 = "The mass of the '`1`' cannot be simultaneously `Heavy` and `0`.";
DefineField::LorentzIndex           = "Lorentz cannot be part of the Indices list.";
DefineField::masslessVector         = "A massless vector field can only be defined as the gauge boson of a gauge group via DefineGaugeGroup. Define the standalone vector '`1`' with a non-zero mass instead.";
DefineField::MassIndices            = "The indices in the Mass option must be a subset of the field flavor Indices";
DefineField::Chiral                 = "The option Chiral can only be different from False for fields of the Fermion type with SelfConjugate->False";
DefineField::ChiralMass             = "A chiral fermion with U(1) charges or complex indices can only be created with a vanishing mass, i.e. by adding the option Mass->0."
DefineField::Complex                = "A real field cannot transform as a complex representation or be charged under an U(1) group."
DefineField::LightIndices           = "DefineField does not support light scalar/vector masses with flavor indices. Light mass should be set to 0 in DefineField and the squared mass should be defined separately with DefineCoupling with EFTOrder->2. ";
DefineField::DoubleGaugeCharges     = "The field contains multiple charges or indices assigned to the same group, only a single charge and a single index is allowed. For non-Abelian groups, define a new representation to have a field transforming under representations other than adjoint and fundamental.";
DefineField::GBscalar               = "Only a scalar field can be a Goldstone boson."
DefineField::GhostFieldsVector      = "The GhostFields option is valid only for Vector fields.";
DefineField::GhostLabel             = "The ghost label '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another label.";
DefineField::nondiagheavymass       = "A heavy field cannot have a non-diagonal mass.";
DefineField::pseudorealmassterm     = "The chiral fermion is in a pseudoreal representation: a mass term is possible only if it has flavor indices. If no mass term is desired, please provide the option \"Mass -> 0\".";
DefineField::pseudorealdiagonalmass = "The chiral fermion is in a pseudoreal representation and a flavor-diagonal mass term is impossible. Anti-symmetric heavy masses are not presently supported in Matchete.";
DefineField::PDGcodeInteger         = "Received invalid PDG code `1` for field `2`. Only integers or Lists/tensors of integers are allowed.";
DefineField::PDGcodeFlavor          = "Received invalid PDG code `1` for flavored field `2`. The Dimension of the Value of \"pdg\" must match the dimensions of the flavor indices of the field, i.e., for every flavor a seperate PDG code has to be given. The PDG code Dimension is `4` whereas the flavor index dimensions are `3`.";
DefineField::PDGnameString          = "Received invalid PDG name `1` for field `2`. Only strings or Lists/tensors of strings are allowed.";
DefineField::PDGnameFlavor          = "Received invalid PDG name `1` for flavored field `2`. The Dimension of the Value of \"name\" must match the dimensions of the flavor indices of the field, i.e., for every flavor a seperate PDG name has to be given. The PDG name Dimension is `4` whereas the flavor index dimensions are `3`.";
DefineField::noPDGcode              = "No PDG code was given to UFO for the field `1`. If the OptionValue of UFO is a non-empty Association it must contain the Key \"pdg\" to uniquely idetify the field for the ExportUFO function.";


DefineField::UFO$OptionType    = "Received invalid UFO$Option \"`3`\"\[Rule]`1` for field `2`. Only `4` or lists/tensors of `4` are allowed.";
DefineField::UFO$OptionFlavor  = "Received invalid UFO$Option \"`5`\"\[Rule]`1` for flavored field `2`. The dimension of the value of \"`5`\" must match the dimensions of the flavor indices of the field, i.e., for every flavor a seperate values must be given. The dimensions of \"`5`\" are `4` whereas the flavor index dimensions are `3`.";


(* ::Text:: *)
(*Thermal Field Theory*)


DefineField::ZeroModeFermion = "Fermion fields do not have massless Matsubara modes so you cannot define a Fermion field with ZeroMode->True.";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[DefineField]={
	BackgroundField-> False,
	Charges-> {},
	Chiral-> False,
	GhostFields-> Automatic,
	GoldstoneBoson-> False,
	Indices->{},
	Mass-> Heavy,
	NiceForm-> Default,
	SelfConjugate-> False,
	UFO -> <||>,
	ZeroMode-> False
};
Options[DefineFieldInternal]= Options[DefineField]


OptionTest[DefineField, BackgroundField]= BooleanQ;
OptionTest[DefineField, Charges]= If[ListQ@ #, And@@ ChargeQ/@ #, ChargeQ@ #]&;
OptionTest[DefineField, Indices]:= If[ListQ@ #,
		And@@ fieldIndQ/@# && DuplicateFreeQ[GroupFromRep/@ Cases[#, Not@* MemberQ[Keys@ $FlavorIndices]]]
	,
		fieldIndQ@ #
	]&;
OptionTest[DefineField, GhostFields]= MatchQ[
	Automatic |
	{Automatic | {_Symbol, _Symbol},
	 Automatic | {_String, _String}}
];
OptionTest[DefineField, NiceForm]= MatchQ[_String | Default | {_String | Default , _String | Default}];
OptionTest[DefineField, UFO] := AssociationQ;


OptionTest[DefineField, Mass]= MatchQ[Heavy| Light| 0| {Light, 0}| {Heavy| Light, _}|
	{Heavy|Light, _, _? (SubsetQ[Keys@ $FlavorIndices[], #]&)}];


OptionTest[DefineField, ZeroMode]= MatchQ[True| False];


fieldIndQ= MemberQ[Keys@ $FlavorIndices, #] || MemberQ[
	Join[Keys@ $GlobalGroups, Keys@ $GaugeGroups], GroupFromRep@ #]&


ChargeQ@ symb_Symbol[_Integer | _Rational | _Symbol]:= MemberQ[GaugeGroupByProperty[Group-> U1], symb];
ChargeQ@ symb_Symbol[a_Times | a_Plus]:= AllTrue[List@@a, ChargeQ@symb[#]&](*And@@(ChargeQ@symb[#]&/@(List@@a))*);
ChargeQ@ _= False;


OptionMessage[BackgroundField, DefineField, val_]:= Message[General::optexpectsval, BackgroundField, DefineField, val, "Boolean"];
OptionMessage[GhostFields, DefineField, val_]:= Message[General::optexpectsval, GhostFields, DefineField, val, "Automatic, or a list {labelsSpec, niceFormSpec} where each spec is either Automatic or a pair {ghostValue, antighostValue}"];
OptionMessage[Charges, DefineField, val_]:= Message[General::optexpectsval, Charges, DefineField, val, "(list of) valid group charge(s)"];
OptionMessage[Indices, DefineField, val_]:= Message[General::optexpectsval, Indices, DefineField, val, "(list of) already defined group representation(s) (only one per group) and/or flavor indice(s)"];
OptionMessage[Mass, DefineField, val_]:= Message[General::optexpectsval, Mass, DefineField, val, "value Heavy, Light, 0, {Light,0}, {Heavy,MassLabel}, {Light,MassLabel}, {Heavy,MassLabel,{FlavorIndex, \[Ellipsis]}} or {Light,MassLabel,{FlavorIndex, \[Ellipsis]}}, with FlavorIndex among the flavor indices of the field, "];
OptionMessage[UFO, DefineField, val_]:=Message[General::optexpectsval, UFO, DefineField, val, "association of necessary UFO parameters, e.g. <|\"pdgcode\" -> ..., ...|>"];


(* ::Subsubsection::Closed:: *)
(*DefineField*)


DefineField[fieldLabel_, type_, opts:OptionsPattern[]]? OptionsCheck:=
DefineFieldInternal[fieldLabel, type, SortOptions[DefineField, opts]]

DefineFieldInternal[fieldLabel_, type_, opts:OptionsPattern[]]:=
DefineFieldInternal[fieldLabel, type, opts]= Module[
	{
		massInfo= OptionValue@ Mass,
		scale,
		flavorInds,
		gh, ghbar,
		ghNF, ghbarNF,
		groupInds,
		massIndDiags,
		massIndExchange,
		isPseudoreal,
		massInds,
		massIndCount,
		massLabel,
		fieldInds,
		fieldCharges,
		tmpAssociation,
		ufoMassOption,
		zeromode= OptionValue@ ZeroMode
	},

	fieldInds= If[Head@ # === List, #, {#}]&@ OptionValue@ Indices;
	fieldCharges= If[Head@ # === List, #, {#}]&@ OptionValue@ Charges;
	fieldCharges= fieldCharges/. ( _[0|0.] -> Nothing ); (* remove charges explicitly specified to vanish *)

	(*Check that the symbol 'label' is not yet in use.*)
	If[Defined[fieldLabel],
			Message[DefineField::FieldLabel,fieldLabel];
			Abort[]
	];

	(*Check that 'type' is valid.*)
	If[!MatchQ[type,Scalar|Vector|Fermion|Ghost|AntiGhost],
		Message[DefineField::FieldType,type];
		Abort[]
	];

	(*Check Lorentz indices are not part on the internal indices *)
	If[MemberQ[fieldInds,Lorentz],
		Message[DefineField::LorentzIndex];
		Abort[]
	];

	(*A massless vector is only consistent as a gauge boson*)
	If[type === Vector && MatchQ[massInfo, 0| {Light, 0}] && !GaugeFieldQ@ fieldLabel && !OptionValue@ BackgroundField,
		Message[DefineField::masslessVector, fieldLabel];
		Abort[]
	];

	(*Check that only complex fermion fields are chiral *)
	If[OptionValue@Chiral=!=False&&(type=!=Fermion||OptionValue@SelfConjugate===True),
		Message[DefineField::Chiral];
		Abort[]
	];

	(* Check that charges appear only once per group *)
	If[!DuplicateFreeQ[Head/@(fieldCharges)]||
			!DuplicateFreeQ[DeleteCases[GroupFromRep/@ fieldInds, None]],
		Message[DefineField::DoubleGaugeCharges];
		Abort[];
	];

	(*Ensure that background fields have zero mass*)
	If[OptionValue@ BackgroundField,
		If[FilterRules[{opts}, Mass] =!= {} && !MatchQ[massInfo, 0| {Light, 0}],
			Message[DefineField::bkgMass];
			Abort[];
		];
		massInfo= 0;
	];

	(* Check that chiral fermion mass is zero if field is charged under any U(1) or has complex indices *)
	If[(type===Fermion && OptionValue@Chiral=!=False
			&& (Length[fieldCharges]=!=0 || Or@@(MemberQ[Select[Keys@$Representations,$Representations[#][Reality]===0&],#]&/@fieldInds)) )
			&&  !MatchQ[massInfo,0|{Light,0}],
		Message[DefineField::ChiralMass];
		Abort[];
	];

	(* Check that indices are not complex when SelfConjugate -> True*)
	If[(OptionValue@SelfConjugate===True
			&& (Length[fieldCharges]=!=0 || Or@@(MemberQ[Select[Keys@$Representations,$Representations[#][Reality]===0&],#]&/@fieldInds))),
		Message[DefineField::Complex];
		Abort[];
	];

	(* Check only scalar fields are given Goldstone boson option*)
	If[OptionValue@ GoldstoneBoson && type =!= Scalar,
		Message[DefineField::GBscalar];
		Abort[];
	];

	(* GhostFields only valid for Vector *)
	If[OptionValue@ GhostFields =!= Automatic && type =!= Vector,
		Message[DefineField::GhostFieldsVector];
		Abort[]
	];
	
	(*For Vetor fields ensure that ghost and anti-ghost labels are valid *)
	If[type === Vector,
		(*Define associated (anti-)ghost fields*)
		{gh, ghbar}= Replace[First[OptionValue@ GhostFields, Automatic],
			Automatic :> {Symbol["gh" <> ToString@ fieldLabel], Symbol["antigh" <> ToString@ fieldLabel]}
		];
		If[Defined@ gh,
			Message[DefineField::GhostLabel, gh];
			Abort[]
		];
		If[Defined@ ghbar,
			Message[DefineField::GhostLabel, ghbar];
			Abort[]
		];
	];

	(* Check UFO *)
	If[Length[OptionValue@ UFO]>0,
		(* check that PDG code is always give if any UFO is given *)
		If[!KeyExistsQ[OptionValue@UFO,"pdg"] && !MatchQ[Keys@OptionValue@UFO,{"mass"}],
			Message[DefineField::noPDGcode,fieldLabel];
			Abort[]
		];
		(* check correct dimensions and formats *)
		If[KeyExistsQ[OptionValue@UFO,#],
			If[FreeQ[OptionValue@Indices, Alternatives@@Keys[GetFlavorIndices[]], All],
				(* non-flavored particles *)
				If[!MatchQ[OptionValue[UFO][#], If[#==="pdg",_Integer,_String]],
					Message[DefineField::UFO$OptionType, OptionValue[UFO][#], fieldLabel,#,If[#==="pdg","integers","strings"]];
					Abort[]
				]
				,
				(* flavored particles *)
				If[GetFlavorIndices[#,IndexDimension]&/@Cases[Flatten@{OptionValue@Indices}, Alternatives@@Keys[GetFlavorIndices[]]]=!=Dimensions[OptionValue[UFO][#]],
					Message[DefineField::UFO$OptionFlavor, OptionValue[UFO][#], fieldLabel, GetFlavorIndices[#,IndexDimension]&/@Cases[OptionValue@Indices, Alternatives@@Keys[GetFlavorIndices[]]],Dimensions[OptionValue[UFO][#]],#];
					Abort[];
				];
				If[!MatchQ[Flatten[OptionValue[UFO][#]], {If[#==="pdg",_Integer,_String]..}],
					Message[DefineField::UFO$OptionType, OptionValue[UFO][#], fieldLabel,#,If[#==="pdg","integers","strings"]];
					Abort[]
				]
			]
		]&/@ {"pdg", "name", "mass", "width"}
	];
	
	(*No fermion zero modes*)
	If[type === Fermion && zeromode,
		Message[DefineField::ZeroModeFermion];
		Abort[]
	];

	(*Check if mass term should have be symmetric or anti-symmetric in flavor indices depending on
		whether the chiral fermions are in a pseudoreal representation or not*)
	isPseudoreal= If[type === Fermion && OptionValue@ Chiral =!= False,
			flavorInds= Intersection[Keys@ $FlavorIndices, fieldInds];
			groupInds= Complement[fieldInds, flavorInds];
			OddQ[Plus@@ (If[$Representations[#, Reality] === -1, 1, 0]&)/@ groupInds]
		,
			False
		];

	(*Interpret mass option*)
	{scale, massLabel, massInds}= Switch[massInfo
		,Heavy|Light,
			If[isPseudoreal,
				Message[DefineField::pseudorealmassterm];
				Abort[];
			];
			{massInfo, Symbol[If[massInfo === Heavy, "M", "m"] <> ToString@ fieldLabel], {}}
		,0| {Light, 0},
			{Light, 0, {}}
		,{Heavy, 0},
			Message[DefineField::Heavy0, fieldLabel];
			Abort[];
		,{_, Except@ 0},
			Append[massInfo, {}]
		,{_, _, _},
			(* Check that mass indices are consistent with field indices *)
			If[!SubsetQ[fieldInds, Last@ massInfo] || (*only field indices*) 
					!FreeQ[Tally[Last@ massInfo][[;;, 2]]-2, _Integer?Positive], (*no index repeated more than twice*)
				Message[DefineField::MassIndices];
				Abort[];
			];
						
			(*Do not allow heavy scalar and vector masses with flavor indices*)
			If[First@ massInfo === Heavy && !DuplicateFreeQ@ Last@ massInfo,
				Message[DefineField::nondiagheavymass];
				Abort[];
			];
			
			(* Purely diagonal pseudoreal fields is not supported *)
			If[isPseudoreal && DuplicateFreeQ@ Last@ massInfo,
				Message[DefineField::pseudorealdiagonalmass];
				Abort[];
			];
			
			massInfo
		];

	(*Define mass coupling*)
	(*Catch/Throw determines the control flow*)
	Catch[
		If[massLabel === 0, Throw[Null];];
		(*Booleans for which indices are diagonal*)
		massIndDiags= Count[massInds, #] === 1 &/@ massInds;
		
		If[scale === Heavy,
			DefineCoupling[massLabel, 
				EFTOrder-> 0,
				SelfConjugate-> True,
				Indices-> massInds,
				DiagonalCoupling-> massIndDiags,
				NiceForm-> Last[OptionValue[NiceForm], Default],
				UFO-> If[KeyExistsQ[OptionValue@UFO, "mass"], <|"name"-> OptionValue[UFO]["mass"]|>, <||>]
			];
			Throw[Null];
		];
		
		ufoMassOption= If[KeyExistsQ[OptionValue@ UFO, "mass"] && 
			$FlavorIndices[#, IndexDimension]&/@ massInds === Dimensions@ OptionValue[UFO]["mass"],
				<|"name"->OptionValue[UFO]["mass"]|>
			,
				<||>
			];
		(*Determines the permutation needed to exchange all non-diagonal indices with their parneters*)
		massIndExchange= Module[{temp}, Table[
				temp= Flatten@ Position[massInds, massInds[[i]], {1}];
				First@ If[Length@ temp > 1, Complement[temp, {i}], temp]
			, {i, Length@ massInds}] ];
		
		(*Chiral fermion*)
		If[type === Fermion && OptionValue@ Chiral =!= False,
			DefineCoupling[massLabel, 
				EFTOrder-> 1,
				SelfConjugate-> False,
				Indices-> massInds,
				DiagonalCoupling-> massIndDiags,
				Symmetries-> {If[isPseudoreal, AntisymmetricPermutation, SymmetricPermutation]@@
					massIndExchange},
				NiceForm-> Last[OptionValue[NiceForm], Default],
				UFO-> ufoMassOption
			];
			Throw[Null];
		];
		
		(*Purely diagonal coupling*)
		If[And@@ massIndDiags,
			DefineCoupling[massLabel, 
				EFTOrder-> 1,
				SelfConjugate-> True,
				Indices-> massInds,
				DiagonalCoupling-> massIndDiags,
				NiceForm-> Last[OptionValue[NiceForm], Default],
				UFO-> ufoMassOption
			];
			Throw[Null];
		];
		
		(*Default*)
		DefineCoupling[massLabel, 
			EFTOrder-> If[type === Fermion, 1, 2],
			SelfConjugate-> If[MatchQ[type, Scalar|Vector] && OptionValue@ SelfConjugate, 
				True, massIndExchange],
			Indices-> massInds,
			DiagonalCoupling-> massIndDiags,
			NiceForm-> Last[OptionValue[NiceForm], Default],
			UFO-> ufoMassOption
		];
	];

	(*Add field to the list of fields*)
	AppendTo[$FieldAssociation,fieldLabel->
		<|
		Type-> type,
		Indices-> fieldInds,
		Charges-> fieldCharges,
		SelfConjugate-> OptionValue@SelfConjugate,
		Chiral-> OptionValue@Chiral,
		Mass-> massLabel,
		Heavy-> scale === Heavy,
		ZeroMode-> zeromode,
		GoldstoneBoson-> OptionValue@ GoldstoneBoson,
		BackgroundField-> OptionValue@ BackgroundField,
		UFO-> OptionValue@ UFO

		|>
	];
	
	(* safe NiceForm printing info *)
	If[!MatchQ[OptionValue[NiceForm], Default|{Default,_}],
		AppendTo[$LabelsNiceForm[Field], fieldLabel-> First[OptionValue[NiceForm], OptionValue[NiceForm]]];
	];

	(*In case of vectors, add Lorentz to indices list *)
	If[type===Vector, PrependTo[fieldInds,Lorentz]];

	(* Create the usage message for the new field *)
	If[Length[fieldInds]=== 0,
		fieldLabel::usage= ToString[fieldLabel]<>"[]: Gives a " <> ToString[type] <> " field with label "<>ToString[fieldLabel]<>".";
	,
		fieldLabel::usage= ToString[fieldLabel]<>"[indices]: Gives a " <> ToString[type] <> " field with label "<>ToString[fieldLabel]<>". The argument must be a sequence of indices following the order: "<> StringJoin[ToString/@Riffle[fieldInds,","]] <>".";
	];
	If[Length[fieldCharges]=!= 0,
		fieldLabel::usage= fieldLabel::usage<>" The charges of the field are "<> ToString[fieldCharges,InputForm] <> "."
	];
	fieldLabel::IndexNumber= "Incorrect number of indices specified. Expected `2` indices, but given `1`.";

	(* Define the new field *)
	With[{fi= fieldInds},
		fieldLabel[indexlist___]:= Module[{i= fi, l=fieldLabel, t=type, inputInd= List[indexlist]},
	
			(*check the length of indices*)
			If[Length[i] =!= Length[inputInd],
				Message[fieldLabel::IndexNumber,Length[inputInd],Length[i]];
				Abort[]
			];
			
			(*Extract Lorentz index for vectors*)
			If[t===Vector,
				t = Vector@Index[Quiet[First@inputInd,First::normal],Lorentz];
				inputInd = Drop[inputInd,1];
				i = Drop[i,1];
			];
			
			If[OptionValue@Chiral===False,
				Field[l, t, Thread@Index[inputInd,i], {}]
			,
				If[OptionValue@Chiral===LeftHanded,PL,PR]\[CenterDot] Field[l, t, Thread@Index[inputInd,i], {}]
			]
		];
	];
	
	(*In case of vectors*)
	If[type === Vector,
		(*Define field-strength tensor*)
		With[{i= Drop[fieldInds, 1]},
			FS[fieldLabel, \[Mu]_, \[Nu]_, indexlist___]:= Module[{inputInd= List@ indexlist},
	
				(*check the length of indices*)
				If[Length[i] =!= Length[inputInd],
					Message[fieldLabel::IndexNumber, Length@ inputInd, Length@ i];
					Abort[];
				];
	
				FieldStrength[fieldLabel, {Index[\[Mu], Lorentz], Index[\[Nu], Lorentz]}, Thread@ Index[inputInd, i], {}]
			];
		];

		{ghNF, ghbarNF}= Replace[Last[OptionValue@ GhostFields, Automatic],
			Automatic :> {"\!\(\*SubscriptBox[\(\[Omega]\), \(" <> ToString@ fieldLabel <> "\)]\)",
			              "\!\(\*SubscriptBox[OverscriptBox[\(\[Omega]\), \(^\)], \(" <> ToString@ fieldLabel <> "\)]\)"}
		];

		fieldInds= fieldInds[[2;;]];
		DefineField[gh, Ghost,
			Indices-> fieldInds,
			Charges-> fieldCharges,
			Mass-> 0,
			SelfConjugate-> OptionValue@ SelfConjugate,
			NiceForm-> ghNF(*,
			UFO-> <|"pdg"-> $pdgCounterGhosts++|>*)
		];
		DefineField[ghbar, AntiGhost,
			Indices-> Bar/@ fieldInds,
			Charges-> MapAt[Minus, fieldCharges, {All, 1}],
			Mass-> 0,
			SelfConjugate-> OptionValue@ SelfConjugate,
			NiceForm-> ghbarNF(*,
			UFO-> <|"pdg"-> $pdgCounterGhosts++|>*)
		];
		$GhostFields[fieldLabel]= {gh, ghbar};
	];

	(*Setup of Conj (used to clasify operator types)*)
	If[OptionValue@ SelfConjugate,
		Conj@ fieldLabel= fieldLabel;
	];
]


(* ::Text:: *)
(*The pdg codes for ghosts need to be handled more carefully once this extended to one loop.*)


(* ::Text:: *)
(*This is a counter for the PDG codes of ghosts:*)


(*$pdgCounterGhosts= 9000001;*)


(* ::Subsection:: *)
(*Remove fields*)


(* ::Subsubsection::Closed:: *)
(*Function*)


RemoveField::rmgaugefield="A gauge fields must be removed by removing the associated gauge group with RemoveGaugeGroup[group] or by resetting all gauge groups with ResetGaugeGroups[]. Continuing without removing the gauge field.";
RemoveField::unkwn= "The label `1` is not a recognized field label.";
RemoveField::gh= "Cannot remove the (anti-)ghost field `1`. Remove the associated vector instead.";


(* ::Text:: *)
(*Public function for removing field definitions*)


RemoveField[label_]:= Module[{l},
	If[!KeyExistsQ[$FieldAssociation, label],
		Message[RemoveField::unkwn, label];
		Abort[];
	];
	If[GaugeFieldQ[label],
		Message[RemoveField::rmgaugefield];
		Abort[];
	];
	
	Switch[$FieldAssociation[label, Type]
	,Ghost|AntiGhost,
		Message[RemoveField::gh, label];
		Abort[];
	,Vector, 
		RemoveFieldInternal@ label;
		Do[
			RemoveFieldInternal@ l;
		, {l, $GhostFields@ label}];
		KeyDropFrom[$GhostFields, label];
	,_,
		RemoveFieldInternal@ label;
	];
];


(* ::Text:: *)
(*Only valid fields are passed to RemoveField Internal*)


RemoveFieldInternal[alias_]:= Module[{mass= $FieldAssociation[alias, Mass]},
	If[KeyExistsQ[$LabelsNiceForm[Field], alias], KeyDropFrom[$LabelsNiceForm[Field], alias]];
	KeyDropFrom[$FieldAssociation, alias];
	(* delete some auxiliary information *)
	If[Head@$currentFieldAssociation === Association && KeyExistsQ[$currentFieldAssociation, alias],
		KeyDropFrom[$currentFieldAssociation, alias];
	];
	If[KeyExistsQ[$currentHeavyDims, alias],
		KeyDropFrom[$currentHeavyDims, alias];
	];
	$currentEOMs= DeleteCases[$currentEOMs, KeyValuePattern[Field[alias|{alias,_},___]:>_], All];
	Quiet[Conj@ alias=.];
	If[mass=!=0 && !MemberQ[Lookup[Mass]@ Values@ $FieldAssociation, mass],RemoveCoupling[mass]];
	(* ClearAllValues[alias] is rather slow here, instead clear specific places *)
	ClearAllValues[alias, {
		(*"DefineField",*)
		"Matchete`Definitions`PackagePrivate`DefineFieldInternal",
		"Matchete`PackageScope`ExpandVectorFluctuations",
		"Matchete`FunctionalTools`PackagePrivate`VarDraw"
	}];
	(*Zero mode is not defined anymore*)
	DefinedZMQ[alias]=False;

	(* unset $currentLagrangian if necessary *)
	If[!FreeQ[$currentLagrangian, alias, All],
		Quiet[$currentLagrangian=.]
	];
	ClearAll@ alias;
];


(* ::Text:: *)
(*Remove all non-gauge fields *)


ResetFields[]:= Module[{fields},
	fields= FieldByProperty[Type-> Scalar|Fermion|Vector];
	fields= Select[fields, Not@* GaugeFieldQ];
	RemoveField/@ fields;
];


(* ::Section:: *)
(*Coupling definitions*)


(* ::Subsection:: *)
(*Coupling setup*)


(* ::Subsubsection::Closed:: *)
(*Properties*)


Coupling[0, __]= 0;


(* ::Subsubsection::Closed:: *)
(*Association of all couplings*)


$CouplingAssociation=<||>;
GetCouplings[CouplingName___]:=$CouplingAssociation[CouplingName]//Map[ReplaceAll[<|x___,Symmetries-> val_,y___|>:><|x,Symmetries->If[Length[List@@val]===1,{}, val],y|>]]


GetCouplingByProperty[propsSeq__]:= GetCouplingByProperty[List[propsSeq]]

GetCouplingByProperty[props_Association]:= GetCouplingByProperty@ Normal@ props;

GetCouplingByProperty[props:_List|_Rule]:=
	Keys@ Select[$CouplingAssociation, MatchQ[#, KeyValuePattern[props]]&];


(* ::Subsubsection::Closed:: *)
(*Remove an already defined coupling*)


RemoveCoupling::rmgaugecoupling="A gauge coupling must be removed by removing the associated gauge group with RemoveGaugeGroup[group] or by resetting all gauge groups with ResetGaugeGroups[]. Continuing without removing the gauge coupling.";


RemoveCoupling[alias_,OptionsPattern[{"gauge"->False}]]:= (
	If[!OptionValue["gauge"] && GaugeCouplingLabelQ[alias],
		Message[RemoveCoupling::rmgaugecoupling];
		Return[]
	];
	If[KeyExistsQ[$LabelsNiceForm[Coupling],alias],KeyDropFrom[$LabelsNiceForm[Coupling], alias]];
	If[KeyExistsQ[$CouplingAssociation,alias],
		KeyDropFrom[$CouplingAssociation, alias];
		(* ClearAllValues[alias] is rather slow here, hence specify where to clean up *)
		ClearAllValues[alias, {
			(*"DefineCoupling",*)
			"Matchete`Definitions`PackagePrivate`DefineCouplingInternal",
			(*"DefineGaugeGroup",*)
			"Matchete`Definitions`PackagePrivate`DefineGaugeGroupInternal",
			"EvaluateLoopFunctions",
			"ThermalPowerCounting"
		}];
		(* unset $currentLagrangian if necessary *)
		If[!FreeQ[$currentLagrangian,alias,All],
			Quiet[$currentLagrangian=.]
		];
		ClearAll@alias;
		(*regenerates the replacement rule to drop diagonal flavor indices on couplings*)
		UpdateDropDiagonalCouplings[];
	];
)


ResetCouplings[]:=Module[{mList1,mList2,cList},
					mList1=$FieldAssociation[#][Mass]&/@(Keys@$FieldAssociation);
					mList2=If[$GaugeGroups===<||>,{},Values@Transpose[$GaugeGroups][Coupling]];
					cList=Keys@$CouplingAssociation;
					RemoveCoupling/@Complement[cList,Join[mList1,mList2]];
				  ];


GaugeCouplingLabelQ[alias_]:=If[MemberQ[$GaugeGroups, KeyValuePattern[Coupling->alias]],
	True,
	False
]


(* ::Subsection:: *)
(*Defining new couplings*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineCoupling::CouplingLabel=    "Cannot use the label '`1`' to define a new coupling. The label '`1`' is either not a Symbol, or it is already used in some loaded context, or it already has some definitions that does not match the given definition exactly. Please use another label or check the properties of the coupling.";
DefineCoupling::LorentzIndex=     "Lorentz cannot be part of the Indices list.";
DefineCoupling::DiagonalCoupling= "The Option DiagonalCoupling must be (list of) boolean(s) of the same length as the number of indices carried by the coupling.";
DefineCoupling::nonmatunitary=    "Only a matrix coupling can be unitary.";
DefineCoupling::unitaryord=       "Unitary couplings must have EFTOrder 0.";

DefineCoupling::PDGnameString  = "Received invalid UFO name `1` for coupling `2`. Only strings or Lists/tensors of strings are allowed.";
DefineCoupling::PDGnameFlavor  = "Received invalid UFO name `1` for flavored coupling `2`. The Dimension of the Value of \"name\" must match the dimensions of the flavor indices of the coupling, i.e., for every flavor index combination a seperate UFO name has to be given. The UFO name Dimension is `4` whereas the flavor index dimensions are `3`.";


(* ::Subsubsection::Closed:: *)
(*Options*)


Options[DefineCoupling]= {
		DiagonalCoupling-> Default,
		EFTOrder-> 0,
		Indices-> {},
		NiceForm-> Default,
		SelfConjugate-> False,
		Symmetries-> {},
		ThermalPowerCounting-> 1,
		UFO-> <||>,
		Unitary-> False
	};
Options[DefineCouplingInternal]= Options[DefineCoupling]


(* ::Text:: *)
(*Option checks and messages  *)


OptionTest[DefineCoupling, DiagonalCoupling]= MatchQ[{_?BooleanQ...}| _?BooleanQ | Default];
OptionTest[DefineCoupling, EFTOrder]= MatchQ[_Integer? NonNegative];
OptionTest[DefineCoupling, Indices]= flavorIndQ@ # || (ListQ[#]&& And@@ flavorIndQ/@#)&;
OptionTest[DefineCoupling, NiceForm]= MatchQ[_String | Default];
OptionTest[DefineCoupling, Symmetries]= (ListQ[#] && And@@((MatchQ[Head[#1],SymmetricIndices|AntisymmetricIndices|SymmetricPermutation|AntisymmetricPermutation] && VectorQ[List@@#1,Positive])&)/@#) || (MatchQ[#,_SymmetryOverride]) &;
OptionTest[DefineCoupling, Unitary]= BooleanQ;
OptionTest[DefineCoupling, ThermalPowerCounting]= MatchQ[_Integer|{_Integer}|_ThermalPowerCounting];


flavorIndQ= MemberQ[Keys@ $FlavorIndices, #]|| MemberQ[Keys@$GlobalGroups, GroupFromRep@ #]&


OptionMessage[DiagonalCoupling, func:DefineCoupling, val_]:= Message[General::optexpectsval, DiagonalCoupling, func, val, "(list of) boolean(s) or Default"];
OptionMessage[EFTOrder, DefineCoupling, val_]:= Message[General::optexpectsval, EFTOrder, DefineCoupling, val, "non-negative integer"];
OptionMessage[Indices, func:DefineCoupling, val_]:= Message[General::optexpectsval, Indices, func, val, "(list of) already defined flavor or global indice(s)"];
OptionMessage[Symmetries, func:DefineCoupling, val_]:= Message[General::optexpectsval, Symmetries, func, val, "a list of terms of the form SymmetricIndices[n1,n2,..], AntisymmetricIndices[n1,n2,..], SymmetricPermutation[n1,n2,..], or AntisymmetricPermutation[n1,n2,..], with n1,n2,... being positive integers indicating the index positions, or an empty list (in case of no symmetries)"];
OptionMessage[Unitary, func:DefineCoupling, val_]:= Message[General::optexpectsval, Unitary, func, val, "Boolean"];


OptionMessage[ThermalPowerCounting, DefineCoupling, val_]:= Message[General::optexpectsval, ThermalPowerCounting, DefineCoupling, val, "integer"];


(* ::Subsubsection::Closed:: *)
(*DefineCoupling*)


(* this creates the coupling association entry both used when defining a new coupling, but also for the duplicate check in IntroduceEffectiveCouplings *)
CouplingAssociationEntry[order_, indInternal_, sc_, symmetries_, dc_, thpwrcounting_, uni_:False, ufo_:<||>]:=
	<|
		EFTOrder-> order,
		Indices-> indInternal,
		SelfConjugate-> sc,
		Symmetries-> symmetries,
		DiagonalCoupling-> dc,
		ThermalPowerCounting-> thpwrcounting,
		Unitary-> uni,
		UFO-> ufo
	|>


(* ::Text:: *)
(*Thread over labels list*)


DefineCoupling[labelList_List, opts:OptionsPattern[]]? OptionsCheck:= 
	(DefineCoupling[#, opts]&/@ labelList;);


(* ::Text:: *)
(*Main function *)


DefineCoupling[label_, opts:OptionsPattern[]]? OptionsCheck:=
DefineCouplingInternal[label, SortOptions[DefineCoupling, opts]]

DefineCouplingInternal[label_, opts:OptionsPattern[]]:=
DefineCouplingInternal[label, opts]= Module[
	{
		order = OptionValue@ EFTOrder,
		indInternal,
		symmetries,
		tmpAssociation,
		diagCoup= OptionValue@ DiagonalCoupling,
		thpwrcounting= OptionValue@ ThermalPowerCounting
	},

	indInternal= If[Head@ # === List, #, {#}]&@ OptionValue@ Indices;

	(*Check that the symbol 'label' is not yet in use.*)
	If[Defined[label],
		Message[DefineCoupling::CouplingLabel,label];
		Abort[]
	];

	(*Check that DiagonalCoupling and Indices have the same length*)
	diagCoup= Switch[diagCoup
		,Default, ConstantArray[False, Length@ indInternal]
		,_?BooleanQ, {diagCoup}
		,_, diagCoup
		];
	If[Length[diagCoup] =!= Length[indInternal],
		Message[DefineCoupling::DiagonalCoupling]; Abort[]
	];

	(*Check Lorentz indices are not part on the internal indices *)
	If[MemberQ[indInternal,Lorentz],
		Message[DefineCoupling::LorentzIndex];
		Abort[]
	];

	(*Perform check of coupling input*)
	symmetries= CouplingSymmetries[Length@ indInternal, OptionValue@ Symmetries];
	
	(*Check that only matrices can be unitary*)
	If[OptionValue@ Unitary && !MatchQ[indInternal, {t_, t_}],
		Message[DefineCoupling::nonmatunitary];
		Abort[];
	];
	(*Check that unitary couplings don't have any power counting*)
	If[OptionValue@ Unitary && order > 0,
		Message[DefineCoupling::unitaryord];
		Abort[];
	];
	
	(* Check UFO *)
	If[Length[OptionValue@ UFO]>0,
		If[KeyExistsQ[OptionValue@UFO,"name"],
			If[FreeQ[OptionValue@Indices, Alternatives@@Keys[GetFlavorIndices[]], All],
				(* non-flavored coupling *)
				If[!MatchQ[OptionValue[UFO]["name"],_String],
					Message[DefineCoupling::PDGnameString, OptionValue[UFO]["name"], label];
					Abort[]
				]
				,
				(* flavored coupling *)
				If[GetFlavorIndices[#,IndexDimension]&/@Cases[OptionValue@Indices, Alternatives@@Keys[GetFlavorIndices[]]]=!=Dimensions[OptionValue[UFO]["name"]],
					Message[DefineCoupling::PDGnameFlavor, OptionValue[UFO]["name"], label, GetFlavorIndices[#,IndexDimension]&/@Cases[OptionValue@Indices, Alternatives@@Keys[GetFlavorIndices[]]],Dimensions[OptionValue[UFO]["name"]]];
					Abort[];
				];
				If[!MatchQ[Flatten[OptionValue[UFO]["name"]],{_String..}],
					Message[DefineCoupling::PDGnameString, OptionValue[UFO]["name"], label];
					Abort[]
				]
			]
		];
	];
	
	(*Add coupling to the list of couplings*)
	$CouplingAssociation@ label= CouplingAssociationEntry[order,
		indInternal, OptionValue@ SelfConjugate, symmetries, diagCoup, thpwrcounting, OptionValue@ Unitary, OptionValue@ UFO];

	(* safe NiceForm printing info *)
	If[OptionValue[NiceForm]=!=Default,
		AppendTo[$LabelsNiceForm[Coupling],label->OptionValue[NiceForm]];
	];

	(* Create the usage message for the new coupling *)
	If[Length[indInternal] === 0,
	label::usage=ToString[label]<>"[]: Gives a coupling with label "<>ToString[label]<>" and EFTOrder "<>ToString[order]<>".";
	,
	label::usage=ToString[label]<>"[indices]: Gives a coupling with label "<>ToString[label]<>" and EFTOrder "<>ToString[order]<>". The argument must be a sequence of indices following the order: "<> StringJoin[ToString/@Riffle[indInternal,","]] <>".";
	];
	label::IndexNumber="Incorrect number of indices specified. Expected `2` indices, but given `1`.";

	(* Define the new coupling *)
	label[indexlist___]:=Module[
		{
			l=label,
			ind=List[indexlist],
			i=indInternal
		},
		(*check the length of indices*)
		If[Length[i] =!= Length[ind],
			Message[label::IndexNumber,Length[ind],Length[i]];
			Abort[]
		];

		Coupling[l, Thread@ Index[ind,i],order]
	];

	(*regenerates the replacement rule to drop diagonal flavor indices on couplings*)
	UpdateDropDiagonalCouplings[];
]


(* ::Subsubsection::Closed:: *)
(*Determine all symmetric permutation of coupling indices *)


(* ::Text:: *)
(*Returns an association with viable coupling permutations with their sign*)


CouplingSymmetries[noInds_, symmetries_]:=Module[{couplingPermutations, evenGens, perm, oddGens, newPerm},
	couplingPermutations=<|Range@ noInds-> 1|>;

	evenGens={}; oddGens={};
	Do[
		Switch[Head@ perm
		,SymmetricIndices,
			evenGens= evenGens~ Join~ {Cycles@ List[List@@ perm], Cycles@ List[List@@ perm[[{1, 2}]]]};
		,AntisymmetricIndices,
			AppendTo[oddGens, Cycles@ List[List@@ perm[[{1, 2}]]]];
			(*oddGens= oddGens~ Join~ {Cycles@ List[List@@ perm], Cycles@ List[List@@ perm[[{1, 2}]]]};*)
			If[OddQ@ Length@ perm,
				AppendTo[evenGens, Cycles@ List[List@@ perm] ];,
				AppendTo[oddGens, Cycles@ List[List@@ perm] ];
			]
		,SymmetricPermutation,
			AppendTo[evenGens, FindPermutation@ perm];
		,AntisymmetricPermutation,
			AppendTo[oddGens, FindPermutation@ perm];
		]
	,{perm, symmetries}];

	(*Compose the symmetries *)
	newPerm= couplingPermutations;
	While[Length@ newPerm > 0,
		newPerm= Association@@ Join[
				KeyValueMap[Function[{key, val}, Permute[key, #]-> val], newPerm]&/@ evenGens,
				KeyValueMap[Function[{key, val}, Permute[key, #]-> -val], newPerm]&/@ oddGens
			];
		KeyDropFrom[newPerm, Keys@ couplingPermutations];
		couplingPermutations= Join[couplingPermutations, newPerm];
	];
	couplingPermutations
]


(* ::Text:: *)
(*Override to directly give the the symmetries in the internal form*)


CouplingSymmetries[ _ , SymmetryOverride[x_]]:= Module[{},
	Switch[Head@x,
		List, Association @@ x,
		Association, x,
		Rule, Association @ x,
		_, Association[ {1,2,3,4}->1 ]
	]
]


(* ::Subsubsection::Closed:: *)
(*Rule for dropping diagonal indices on Couplings*)


(* ::Text:: *)
(*Rule to drop diagonal flavor indices from couplings when relabeling indices. Defined here once globally for better performance.*)


UpdateDropDiagonalCouplings[]:= Block[{},
	$DropDiagonalCouplings= Join[{_FlavorSum->1}, Table[
			With[{tmp= Position[$CouplingAssociation[x][DiagonalCoupling],True], lab= x},
				If[tmp==={},
					Nothing,
					Coupling[x,ind_List,ord_]:> Coupling[lab, ReplacePart[ind,tmp:>Unique[]],ord] (* Unique ensures that 1/(m[p]-m[r]) yields 1/(m[]-m[])=ComplexInfinity *)
				]
			]
			,
			{x,Keys@$CouplingAssociation}
		]
	];
];


UpdateDropDiagonalCouplings[]


(* ::Section:: *)
(*Flavor indices*)


(* ::Subsection:: *)
(*Flavor indices setup*)


(* ::Subsubsection::Closed:: *)
(*Association of all flavor indices*)


$FlavorIndices= <||>;
GetFlavorIndices[FlavorIndex___]:=Return@$FlavorIndices[FlavorIndex];


(* ::Subsubsection::Closed:: *)
(*Removing  flavor indices*)


RemoveFlavorIndex[alias_]:= Module[{},
					       If[KeyExistsQ[$FlavorIndices,alias],
					          RemoveCoupling/@CouplingsFromGroup[alias];
					          RemoveField/@FieldsFromGroup[alias];
					          KeyDropFrom[$IndexAlphabets, alias];
							  KeyDropFrom[$FlavorIndices, alias];
							  GenerateDeltaContractionPatterns[]; (* this sets up the new Delta contraction patterns *)
							  Bar[Index[ind_,alias]]=.;
							  ClearAllValues[alias];
							  ClearAll@alias;
					       ]
					      ]


ResetFlavorIndices[]:=(RemoveFlavorIndex/@Keys[$FlavorIndices[]];);


(* ::Subsection:: *)
(*Defining flavor indices*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineFlavorIndex::FlavorName     = "The flavor index name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another flavor index name.";
DefineFlavorIndex::IndexDimension = "The flavor index dimension '`1`' is not a positive integer larger than 1 or an undefined Symbol."


(* ::Subsubsection::Closed:: *)
(*Options*)


(* Define the option patterns for DefineFlavorIndices *)
Options[DefineFlavorIndex]={IndexAlphabet -> None};
Options[DefineFlavorIndexInternal]= Options[DefineFlavorIndex]


(* ::Subsubsection::Closed:: *)
(*DefineFlavorIndices*)


(* ::Text:: *)
(*Function to define a flavor indices*)


DefineFlavorIndex[flavorName_, indexDim_, opts:OptionsPattern[]]? OptionsCheck:=
DefineFlavorIndexInternal[flavorName, indexDim, SortOptions[DefineFlavorIndex, opts]]

DefineFlavorIndexInternal[flavorName_, indexDim_, opts:OptionsPattern[]]:= 
DefineFlavorIndexInternal[flavorName, indexDim, opts]= Module[{},
	(*Check if 'Name' is already defined as a flavor index.*)
	If[Defined[flavorName]||KeyExistsQ[$FlavorIndices,flavorName]||KeyExistsQ[$GaugeGroups,flavorName]||KeyExistsQ[$GlobalGroups,flavorName],
		Message[DefineFlavorIndex::FlavorName,flavorName];
		Abort[]
	];

	(*Check if 'IndexDimension' is an undefined symbol or a positive integer larger than 1.*)
	If[Defined[indexDim]&&!(IntegerQ[indexDim] && indexDim > 1),
		Message[DefineFlavorIndex::IndexDimension,indexDim];
		Abort[]
	];

	(*Update list of flavor indices*)
	AppendTo[$FlavorIndices, flavorName-> <|
			IndexDimension-> indexDim
	|>];
	(* this sets up the new delta contraction patterns *)
	GenerateDeltaContractionPatterns[];

	(* If provided, add index alphabet for the representation *)
	If[OptionValue@IndexAlphabet=!=None,AppendTo[$IndexAlphabets,flavorName -> BuildIndexAssoc[OptionValue@IndexAlphabet]]];

	(* Define a usage message for this object *)
	flavorName::usage=ToString[flavorName]<>" is a label for a flavor index.";

	(* Define flavor indices to be neutral under Bar[]*)
	Bar@Index[ind_,flavorName]:=Index[ind,flavorName];
]


(* ::Section:: *)
(*Gauge and global groups*)


(* ::Subsubsection::Closed:: *)
(*Representation and  Clebsch Gordan info*)


(* ::Text:: *)
(*Shows the representation indices for a given Clebsch Gordan*)


ClebschGordanIndices::ClebschGordanName   = "The Clebsch-Gordan coefficient has not been defined. Please use another name or use DefineCG to define it."


ClebschGordanIndices[CGname_]:= Module[{},
								(*Check that the Clebsch-Gordan is defined.*)
								If[!KeyExistsQ[$CGproperties,CGname],
									Message[ClebschGordanIndices::ClebschGordanName,CGname];
									Abort[]
								];
							    $CGproperties[CGname][Indices]
							  ]


(* ::Subsubsection::Closed:: *)
(*GroupFromRep*)


(* ::Text:: *)
(*For finding the group associated to a representation*)


GroupFromRep@ Bar@ rep_:= GroupFromRep@ rep;
GroupFromRep@ rep_:= Module[{gr},
	gr= Lookup[$Representations, rep, None]@ Group;
	If[KeyExistsQ[$GaugeGroups, gr]||KeyExistsQ[$GlobalGroups, gr], gr, None]
];


(* ::Text:: *)
(*Returns the group associated to an index*)


GroupFromInd@ Bar@ ind_:= GroupFromInd@ ind;
GroupFromInd@ Index[_, rep_]:= GroupFromRep@ rep;


(* ::Text:: *)
(*Returns the representation associated to an index*)


RepFromInd@ Bar@ ind_:= RepFromInd@ ind;
RepFromInd@ Index[_, rep_]:= rep;


(* ::Subsubsection::Closed:: *)
(*FieldsFromGroup*)


(* ::Text:: *)
(*Gives all fields associated with the group: fields charged under the group or the associated vector and ghost fields*)


FieldsFromGroup[group_]:= Module[{fields, vec},
	If[MemberQ[Keys@ $FlavorIndices, group],
		Return@ Select[Keys@ $FieldAssociation, MemberQ[Transpose[$FieldAssociation][Indices][#], group]&];
	];
	fields= If[MatchQ[Join[$GaugeGroups,$GlobalGroups][group][Group], U1],
			Select[Keys@$FieldAssociation,MemberQ[Head/@Transpose[$FieldAssociation][Charges][#],group]&]
		,
			Select[Keys@$FieldAssociation,MemberQ[GroupFromRep/@Transpose[$FieldAssociation][Indices][#],group]&]
		];
	If[MemberQ[Keys@ $GaugeGroups, group],
		vec= $GaugeGroups[group, Field];
		fields= DeleteDuplicates@ Join[fields, Append[$GhostFields@ vec, vec]];
	];
	fields
];


(* ::Subsubsection::Closed:: *)
(*CouplingsFromGroup*)


CouplingsFromGroup[group_]:= Module[{},
	If[MemberQ[Keys@$FlavorIndices,group],
		Return[Select[Keys@$CouplingAssociation,MemberQ[Transpose[$CouplingAssociation][Indices][#],group]&]];
		,
		If[MatchQ[Join[$GaugeGroups,$GlobalGroups][group][Group], U1],
			Return[Select[Keys@$CouplingAssociation,MemberQ[Head/@Transpose[$CouplingAssociation][Charges][#],group]&]];
			,
			Return[Select[Keys@$CouplingAssociation,MemberQ[GroupFromRep/@Transpose[$CouplingAssociation][Indices][#],group]&]];
		]
	]
]


(* ::Subsection:: *)
(*Get field generators*)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


FieldGenerators::Field       = "The field '`1`' has not been defined. Please define the field using DefineField or use another field name.";
FieldGenerators::GroupName   = "The group '`1`' has not been defined. Please define the group using DefineGaugeGroup or DefineGlobalGroup, or use another group name.";
FieldGenerators::Abelian     = "The group '`1`' is Abelian so no indices should be provided. Please use FieldGenerators[`2`, `1`] instead or use a different group name."
FieldGenerators::NotAbelian  = "The group '`1`' is not Abelian so generator indices should be provided. Please use FieldGenerators[`2`, `1`,Indices] instead or use a different group name."
FieldGenerators::Charged     = "The field '`1`' is not charged under the group '`2`'. Please provide another group name."
FieldGenerators::Indices     = "The indices argument '`1`' is incorrect. Please provide a list of three Symbols";


(* ::Subsubsection::Closed:: *)
(*FieldGenerators*)


(* ::Text:: *)
(*For finding the generators associated to a given field and gauge group*)
(*With Hermitian generators, the conjugate representation has generators -T**)


FieldGenerators[Bar@field_,GroupName_,indices_]:= -Bar@ FieldGenerators[field,GroupName,indices]
FieldGenerators[field_,GroupName_,indices_]:= Module[{rep,generators},

	(*Check that the field is already defined*)
	If[!KeyExistsQ[$FieldAssociation,field],
		Message[FieldGenerators::Field,field];
		Abort[]
	];

	(*Check that the group is already defined*)
	If[!KeyExistsQ[$GaugeGroups,GroupName] && !KeyExistsQ[$GlobalGroups,GroupName],
		Message[FieldGenerators::GroupName,GroupName];
		Abort[]
	];

	(* Check that the group is not Abelian *)
	If[$GaugeGroups[GroupName][Abelian] || $GlobalGroups[GroupName][Abelian],
		Message[FieldGenerators::Abelian,GroupName,field];
		Abort[]
	];

	(*Check that the field is charged under the given group*)
	If[!MemberQ[DeleteCases[GroupFromRep/@$FieldAssociation[field][Indices],None],GroupName],
		Message[FieldGenerators::Charged,field,GroupName];
		Abort[]
	];

	(*Check that indices is a list of three indices*)
	If[!MatchQ[indices,{_Symbol,_Symbol,_Symbol}],
		Message[FieldGenerators::Indices,indices];
		Abort[]
	];

	rep=First@Select[$FieldAssociation[field][Indices],GroupFromRep[#]===GroupName &];
	If[Head@rep === Bar,
		-Bar@ CG[gen[rep[[1]]], indices]
	,
		CG[gen[rep],indices]
	]
];


(* ::Text:: *)
(*In the case of Abelian fields it returns the charge *)


FieldGenerators[Bar@field_,GroupName_]:=Bar@FieldGenerators[field,GroupName]
FieldGenerators[field_,GroupName_]:= Module[{generators},

	(*Check that the field is already defined*)
	If[!KeyExistsQ[$FieldAssociation,field],
		Message[FieldGenerators::Field,field];
		Abort[]
	];

	(*Check that the group is already defined*)
	If[!KeyExistsQ[$GaugeGroups,GroupName] && !KeyExistsQ[$GlobalGroups,GroupName],
		Message[FieldGenerators::GroupName,GroupName];
		Abort[]
	];

	(* Check that the group is Abelian *)
	If[!$GaugeGroups[GroupName, Abelian] || !$GlobalGroups[GroupName, Abelian],
		Message[FieldGenerators::NotAbelian,GroupName,field];
		Abort[]
	];

	(*Check that the field is charged under the given group*)
	If[!MemberQ[Head/@$FieldAssociation[field][Charges],GroupName],
		Message[FieldGenerators::Charged,field,GroupName];
		Abort[]
	];

	Select[$FieldAssociation[field][Charges],Head[#]==GroupName&][[1]][[1]]
];


(* ::Section:: *)
(*Gauge groups*)


(* ::Subsection:: *)
(*Gauge setup*)


(* ::Subsubsection::Closed:: *)
(*Field-strength and structure constant properties*)


(* ::Text:: *)
(*Field strength*)


FieldStrength[0, __]:= 0;
FieldStrength[_, {\[Mu]_, \[Mu]_}, __]:= 0;
FieldStrength[_, _, {a_, a_}, __]:= 0;
FieldStrength[_, _, {a_, Bar@ a_}, __]:= 0;
FieldStrength[_, _, {Bar@ a_, a_}, __]:= 0;


(* Order the Lorentz indices canonically *)
FieldStrength[label_, {\[Mu]_,\[Nu]_}, rest___] := -FieldStrength[label, {\[Nu],\[Mu]}, rest] /; !OrderedQ[{\[Mu],\[Nu]}]


(* ::Text:: *)
(*Order indices of structure constants*)


CG[fStruct[group_], indices_List]/; !OrderedQ[indices]:=
	Signature[indices]CG[fStruct[group], Sort@ indices]


(* ::Subsubsection::Closed:: *)
(*Association of all gauge groups*)


(* ::Text:: *)
(*Basic gauge group information*)


$GaugeGroups= <||>;
GetGaugeGroups[GaugeGroupName___]:=Module[{$tmpGauge},
	$tmpGauge=$GaugeGroups;
	Table[
		AppendTo[$tmpGauge[i], ClebschGordanCoefficients -> Keys@Select[$CGproperties,#[Group]===i&]];
		AppendTo[$tmpGauge[i], Representations -> Keys@Select[$Representations,#[Group]===i&]];
	,{i,Keys@$GaugeGroups}];
	$tmpGauge[GaugeGroupName]
];


GaugeGroupByProperty[]:= Keys@ $GaugeGroups;
GaugeGroupByProperty[propsSeq__]:= GaugeGroupByProperty@ List@ propsSeq;
GaugeGroupByProperty[props_Association]:= GaugeGroupByProperty@ Normal@ props;
GaugeGroupByProperty[props: _List| _Rule]:=
	Keys@ Select[$GaugeGroups, MatchQ[#, KeyValuePattern[props]]&];


(* ::Subsubsection::Closed:: *)
(*Removing gauge group*)


RemoveGaugeGroup::notgg= "`1` is not a previously defined gauge group.";


(* ::Text:: *)
(*Public function for removing gauge groups*)


RemoveGaugeGroup[groupName_]:= Block[{},
	If[!KeyExistsQ[$GaugeGroups, groupName],
		Message[RemoveGaugeGroup::notgg, groupName];
		Abort[];
	];
	
	RemoveCoupling[$GaugeGroups[groupName, Coupling], "gauge"-> True];
	RemoveFieldInternal/@ FieldsFromGroup[groupName];
	KeyDropFrom[$GhostFields, $GaugeGroups[groupName, Field]];
	
	RemoveCG/@ GetGaugeGroups[groupName][ClebschGordanCoefficients];
	RemoveRepresentation/@ GetGaugeGroups[groupName][Representations];
	KeyDropFrom[$IndexAlphabets, Keys@ Select[$Representations, #[Group] === groupName &]];
	If[!MatchQ[$GaugeGroups[groupName, Group], U1], KeyDropFrom[$Groups, groupName];];
	KeyDropFrom[$GaugeGroups, groupName];
	ClearAllValues@ groupName;
	ClearAll@ groupName;
	DefinedZMQ[groupName]=False;
];


ResetGaugeGroups[]:=(RemoveGaugeGroup/@ Keys[$GaugeGroups[]];);


(* ::Subsection:: *)
(*Defining gauge group *)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineGaugeGroup::CouplingLabel   = "The coupling '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another symbol for the gauge coupling.";
DefineGaugeGroup::GroupName       = "The group name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another group name.";
DefineGaugeGroup::GaugeFieldLabel = "The gauge field name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another symbol for the gauge field.";


(* ::Subsubsection::Closed:: *)
(*DefineGaugeGroup*)


(* ::Text:: *)
(*Function to define a gauge group together with some basic representations and Clebsch-Gordan coefficients*)


Options[DefineGaugeGroup]= {AdjAlphabet-> None, FundAlphabet-> None, GhostFields-> Automatic, NiceForm-> Default, UFO -> <||>, ZeroMode-> False}
Options[DefineGaugeGroupInternal]= Options[DefineGaugeGroup]


OptionTest[DefineGaugeGroup, GhostFields]= MatchQ[
	Automatic |
	{Automatic | {_Symbol, _Symbol},
	 Automatic | {_String, _String}}
];
OptionTest[DefineGaugeGroup, NiceForm]= MatchQ[Default | {_String | Default, _String | Default}];
OptionTest[DefineGaugeGroup, UFO] := AssociationQ;
OptionTest[DefineGaugeGroup, ZeroMode]= MatchQ[True | False];


OptionMessage[GhostFields, DefineGaugeGroup, val_]:= Message[General::optexpectsval, GhostFields, DefineGaugeGroup, val, "Automatic, or a list {labelsSpec, niceFormSpec} where each spec is either Automatic or a pair {ghostValue, antighostValue}"];
OptionMessage[NiceForm, DefineGaugeGroup, val_]:= Message[General::optexpectsval, NiceForm, DefineGaugeGroup, val, "value Default or a list with two elements each of which is either Default or a string"];
OptionMessage[UFO, DefineGaugeGroup, val_]:=Message[General::optexpectsval, UFO, DefineGaugeGroup, val, "association of necessary UFO parameters, e.g. <|\"pdgcode\" -> ..., ...|>"];


DefineGaugeGroup[grName_,lieGroup_, coupling_, gaugeField_, opts:OptionsPattern[]]? OptionsCheck:=
DefineGaugeGroupInternal[grName,lieGroup, coupling, gaugeField, SortOptions[DefineGaugeGroup, opts]]

(* This splitting is necessary so that the cached results do not depend on the ordering of the options *)

DefineGaugeGroupInternal[grName_,lieGroup_, coupling_, gaugeField_, opts:OptionsPattern[]]:=
DefineGaugeGroupInternal[grName,lieGroup, coupling, gaugeField, opts]= Module[{fundRep, adjRep, rank, nfCoupling, nfField},
	(*Check that the symbol 'coupling' is not yet in use.*)
	If[Defined[coupling] && !KeyExistsQ[$GaugeGroups,coupling],
		Message[DefineGaugeGroup::CouplingLabel,coupling];
		Abort[]
	];

	(*Check that the symbol 'grName' is not yet in use.*)
	If[Defined[grName]||KeyExistsQ[$FlavorIndices,grName]||KeyExistsQ[$GaugeGroups,grName]||KeyExistsQ[$GlobalGroups,grName],
		Message[DefineGaugeGroup::GroupName,grName];
		Abort[]
	];

	(*Check that the symbol 'gaugeField' is not yet in use.*)
	If[Defined[gaugeField],
		Message[DefineGaugeGroup::GaugeFieldLabel, gaugeField];
		Abort[]
	];

	(*Check that the ghost labels are not yet in use.*)
	With[{ghLabels= Replace[First[OptionValue@ GhostFields, Automatic],
		Automatic :> {Symbol["gh" <> ToString@ gaugeField], Symbol["antigh" <> ToString@ gaugeField]}
	]},
		If[Defined@ First@ ghLabels, Message[DefineField::GhostLabel, First@ ghLabels]; Abort[]];
		If[Defined@ Last@ ghLabels, Message[DefineField::GhostLabel, Last@ ghLabels]; Abort[]]
	];

	(*Check that the symbol 'lieGroup' is a valid simple Lie group.*)
	If[!MatchQ[lieGroup, U1], TestAlg@ lieGroup];
	
	(*Update list of gauge groups*)
	AppendTo[$GaugeGroups, grName-> <|
			Group-> lieGroup,
			Coupling-> coupling,
			Field-> gaugeField,
			Abelian-> lieGroup === U1
	|>];

	(*Define usage message*)
	grName::usage=ToString[grName]<>" is the label of a gauge group with gauge field "<>ToString[gaugeField]<> " and gauge coupling "<>ToString[coupling]<>".";

	(*Set up Lie group and CG definitions for non-Ablelian gauge groups*)
	If[lieGroup=!= U1,
		AddNonAbelianGroupDefinitions[grName, lieGroup, FilterOptions[opts, AddNonAbelianGroupDefinitions] ];
	];
	
	nfCoupling= First[OptionValue[NiceForm], Default];
	nfField= Last[OptionValue[NiceForm], Default];

	(*Add gauge coupling*)
	DefineCoupling[coupling, SelfConjugate -> True, NiceForm->nfCoupling, UFO-> If[KeyExistsQ[OptionValue@UFO, "coupling"],<|"name"->OptionValue[UFO]["coupling"]|>,<||>]];

	(*Add gauge field*)
	DefineField[gaugeField, Vector,
		Indices-> If[lieGroup =!= U1, {grName@ adj}, {}],
		Mass-> 0, (*This is allowed only because the field has been already been added to $GaugeGroups*)
		SelfConjugate-> True,
		NiceForm-> nfField,
		GhostFields-> OptionValue@ GhostFields,
		UFO-> KeyDrop[OptionValue@ UFO, "coupling"],
		ZeroMode-> OptionValue@ ZeroMode
	];
];


(* ::Text:: *)
(*Makes all the standard CG definitions for non-Abelian groups*)


Options@ AddNonAbelianGroupDefinitions= {AdjAlphabet-> None, FundAlphabet -> None}; 


AddNonAbelianGroupDefinitions[grName_, lieGroup_, OptionsPattern[]]:= Module[{adjRep, fundRep, rank, $a, $i, $j, $b, $c},
	(*Add group*)
	DefineCGGroup[grName, lieGroup];
	adjRep= AdjointRepresentation@ lieGroup;
	fundRep= FundamentalRepresentation@ lieGroup;
	rank= AlgebraRank@ lieGroup;
	
	(*Add adjoint representation*)
	DefineRepresentation[grName@ adj, grName, adjRep, IndexAlphabet-> OptionValue@ AdjAlphabet];
	(*Structure constants and identification with adj generators*)
	DefineCG[fStruct@ grName, {grName@ adj, grName@ adj, grName@ adj}, StructureConstants@ lieGroup];
	CG[gen[grName[adj]], indices_]:= -I CG[fStruct@ grName, indices];
	
	(*Add fundamental representaion*)
	If[!MatchQ[lieGroup, Alg["E", 8]|Alg["B", 1]], (*For B1 and E8, the fundamental conincide with the Adjoint*)
		DefineRepresentation[grName@ fund, grName, fundRep,
			IndexAlphabet -> OptionValue@ FundAlphabet];
		(*Define the commonly occuring t^ai_j f^{abc}*)
		DefineCompositeCG[tFundf@ grName, {gen@ grName@ fund, fStruct@ grName}, {{$a, $i, $j}, {$a, $b, $c}}];
	];

	(*Add LC CGs for SU(N) algebras (limited to N\[LessEqual]5 for performance reasons)*)
	If[MatchQ[lieGroup, Alg["A", _?(# > 0 && # < 5 &)]],
		DefineCG[eps@grName, ConstantArray[grName@ fund, rank+ 1],
			First@InvariantTensors[lieGroup, Table[fundRep, rank+ 1],
				AntisymmetricIndices-> Range[rank+ 1],
				Normalization-> (rank+ 1)!]];
	];
	(*Add symmetric 3-index adjoint tensor for SU(N>2)*)
	If[MatchQ[lieGroup, Alg["A", _?(# > 1 &)]],
		DefineCG[dSym@ grName, {grName@ adj, grName@ adj, grName@ adj},
			First@ InvariantTensors[lieGroup, {adjRep, adjRep, adjRep},
				SymmetricIndices-> {1, 2, 3},
				Normalization-> ((rank+1)^2-1)((rank+1)^2-4)/(rank+1)]];
	];

	(*Add two-index invariant CGs for Sp(N) algebras*)
	If[MatchQ[lieGroup, Alg["C", _]],
		DefineCG[eps@ grName, Table[grName@ fund, 2],
			First@ InvariantTensors[lieGroup, Table[fundRep, 2],
				AntisymmetricIndices-> {1, 2}, Normalization-> 2* lieGroup[[2]]]];
	];

]


(* ::Section:: *)
(*Global groups*)


(* ::Subsection:: *)
(*Global groups setup*)


(* ::Subsubsection::Closed:: *)
(*Association of all global groups*)


(* ::Text:: *)
(*Basic global group information*)


$GlobalGroups= <||>;
GetGlobalGroups[GlobalGroupName___]:=Module[{$tmpGlobal},
	$tmpGlobal=$GlobalGroups;
	Table[
		AppendTo[$tmpGlobal[i], ClebschGordanCoefficients -> Keys@Select[$CGproperties,#[Group]===i&]];
		AppendTo[$tmpGlobal[i], Representations -> Keys@Select[$Representations,#[Group]===i&]];
	,{i,Keys@$GlobalGroups}];
	$tmpGlobal[GlobalGroupName]
];


(* ::Subsubsection::Closed:: *)
(*Removing global group*)


RemoveGlobalGroup[groupName_]:=
	If[KeyExistsQ[$GlobalGroups, groupName],
		RemoveCG/@ GetGlobalGroups[groupName][ClebschGordanCoefficients];
		RemoveRepresentation/@ GetGlobalGroups[groupName][Representations];
		KeyDropFrom[$IndexAlphabets, Keys@ Select[$Representations, #[Group] === groupName &]];
		If[!MatchQ[$GlobalGroups[groupName, Group], U1], KeyDropFrom[$Groups, groupName];];
		KeyDropFrom[$GlobalGroups, groupName];
		ClearAllValues@ groupName;
		ClearAll@ groupName;
	];


ResetGlobalGroups[]:=(RemoveGlobalGroup/@Keys[$GlobalGroups[]];);


(* ::Subsection:: *)
(*Defining global group *)


(* ::Subsubsection::Closed:: *)
(*Error messages*)


DefineGlobalGroup::GroupName = "The group name '`1`' is not a Symbol, it is already used in some loaded context or already has some definitions. Please use another group name.";


(* ::Subsubsection::Closed:: *)
(*Options*)


(* Define the option patterns for DefineGlobalGroup *)
Options[DefineGlobalGroup]= {FundAlphabet-> None, AdjAlphabet-> None}
Options[DefineGlobalGroupInternal]= Options[DefineGlobalGroup]


(* ::Subsubsection::Closed:: *)
(*DefineGlobalGroup*)


(* ::Text:: *)
(*Function to define a global group together with some basic representations and Clebsch-Gordan coefficients*)


DefineGlobalGroup[grName_, lieGroup_, opts:OptionsPattern[]]? OptionsCheck:=
DefineGlobalGroupInternal[grName, lieGroup, SortOptions[DefineGlobalGroup, opts]]


DefineGlobalGroupInternal[grName_, lieGroup_, opts:OptionsPattern[]]:=
DefineGlobalGroupInternal[grName, lieGroup, opts]= Module[{},
	(*Check that the symbol 'grName' is not yet in use.*)
	If[Defined[grName]||KeyExistsQ[$FlavorIndices,grName]||KeyExistsQ[$GaugeGroups,grName]||KeyExistsQ[$GlobalGroups,grName],
		Message[DefineGlobalGroup::GroupName,grName];
		Abort[]
	];

	(*Check that the symbol 'lieGroup' is a valid simple Lie group.*)
	If[!MatchQ[lieGroup, U1], TestAlg@ lieGroup];

	(*Update list of global groups*)
	AppendTo[$GlobalGroups, grName-> <|
			Group-> lieGroup,
			Abelian-> lieGroup === U1
	|>];

	(*Define a usage message*)
	grName::usage=ToString[grName]<> " is the label of a global symmetry.";

	(*Set up Lie group and CG definitions for non-Ablelian gauge groups*)
	If[lieGroup=!= U1,
		AddNonAbelianGroupDefinitions[grName, lieGroup, FilterOptions[opts, AddNonAbelianGroupDefinitions] ];
	];
];


(* ::Section:: *)
(*Reset all*)


(* ::Subsubsection::Closed:: *)
(*Remove all global data*)


ResetAll[]:= Block[{},
		Quiet[$currentLagrangian=.];
		RemoveAssociatedDownValues@ LoadModel[_String, ___]; (* Clear LoadModel cache*)
		ResetAuxiliaryMatchingInformation[];
		ResetGaugeGroups[];
		ResetCouplings[];
		ResetFields[];
		ResetFlavorIndices[];
		ResetGlobalGroups[];
		ResetIndexAlphabets[];
		ResetInternalCouplings[];
		ResetOperatorAssociations[];
		ResetSymmetryBreaking[];
		ResetSubstitutionInformation[];
		ResetTempCouplings[];
		ResetTrigonometrics[];
		ResetEvanescentOperators[];
		
		RemoveFlavorStructures[]; (* removes all defiend flavor invariants/symmetries *)
		ResetCouplingOrders[];    (* resets couplings orders for UFO export *)

		ClearGroups[];(*needed to reset the CG numbering*)
		
		SetDimensions[4];
	];


(* ::Section:: *)
(*Removing values associated with all Matchete` symbols*)


(* ::Subsubsection::Closed:: *)
(*Function for removing all Up-, Down-, and SubValues off all functions that contain a given symbol*)


PackageScope["ClearAllValues"]


ClearAllValues[var_]:=Do[
	If[NameQ[function],
		With[{func=function},
			If[!FreeQ[DownValues[func],var,All],
				DownValues[func]=DeleteCases[DownValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			If[!FreeQ[UpValues[func],var,All],
				UpValues[func]=DeleteCases[UpValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			If[!FreeQ[SubValues[func],var,All],
				SubValues[func]=DeleteCases[SubValues[func],x_/;!FreeQ[First@x,var,All]];
			];
			(* OwnValues cannot contain var on the LHS *)
			(*OwnValues[func]=DeleteCases[OwnValues[func],x_/;!FreeQ[First@x,var,All]];*)
		]
	];
	,
	{function,Names["Matchete`"~~___]}
]


(* ::Subsubsection::Closed:: *)
(*Function clearing all Up-, Down-, and SubValues off a specific functions that contain a given symbol*)


(* ::Text:: *)
(*Note: func must be given as a String including the full Context*)


ClearAllValues[var_, func_]:=If[NameQ[func],
	DownValues[func] = DeleteCases[DownValues[func], x_/;!FreeQ[First@x,var,All]];
	UpValues[func]   = DeleteCases[UpValues[func],   x_/;!FreeQ[First@x,var,All]];
	SubValues[func]  = DeleteCases[SubValues[func],  x_/;!FreeQ[First@x,var,All]];
	(* OwnValues cannot contain var on the LHS *)
	(*OwnValues[func]=DeleteCases[OwnValues[func],x_/;!FreeQ[First@x,var,All]];*)
]


(* ::Text:: *)
(*Note: functions must be given as a List of Strings including the full Context*)


ClearAllValues[var_, functions_List]:=Do[
	If[NameQ[func],
		With[{f=func},
			DownValues[f] = DeleteCases[DownValues[f], x_/;!FreeQ[First@x,var,All]];
			UpValues[f]   = DeleteCases[UpValues[f],   x_/;!FreeQ[First@x,var,All]];
			SubValues[f]  = DeleteCases[SubValues[f],  x_/;!FreeQ[First@x,var,All]];
			(* OwnValues cannot contain var on the LHS *)
			(*OwnValues[f]=DeleteCases[OwnValues[f],x_/;!FreeQ[First@x,var,All]];*)
		]
	]
	,
	{func,functions}
]


(* ::Subsubsection::Closed:: *)
(*Function for checking in which definitions a symbol appears*)


PackageScope["GetAllValues"]

GetAllValues::usage="GetAllValues[var] lists all DownValues, UpValues, SubValues, and OwnValues of any symbol inside the Matchete` context or its subcontexts that contains the symbol var.";

GetAllValues[var_]:=Do[
	If[NameQ[function],
		With[{func=function},
			Cases[DownValues[func], x_/;!FreeQ[x,var,All] :> Echo[x, "DownValue"] ];
			Cases[UpValues[func],   x_/;!FreeQ[x,var,All] :> Echo[x, "UpValue"]   ];
			Cases[SubValues[func],  x_/;!FreeQ[x,var,All] :> Echo[x, "SubValue"]  ];
			Cases[OwnValues[func],  x_/;!FreeQ[x,var,All] :> Echo[x, "OwnValue"]  ];
		]
	]
	,
	{function,Names["Matchete`"~~___]}
]


(* ::Section:: *)
(*Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Hermitian conjugate*)


(* ::Text:: *)
(*Ease of life function for adding '+H.c.' terms to the Lagrangian*)


PlusHc[x_]:=x+Bar@x


(* ::Subsubsection::Closed:: *)
(*Free Lagrangian*)


FreeLag::bkgfield = "The field '`1`' is a background field and has no associated kinetic term/free Lagrangian.";
FreeLag::UndefinedField = "The field '`1`' has not been defined.";


FreeLag[field_Symbol]:= Module[
	{fieldType, isChiralFermion, massDiags, massInd, normalization, fieldInd, indK, ind1, ind2, indM, m, i, j, a, \[Mu], \[Nu]},

	(* Check if the field has been defined *)
	If[!KeyExistsQ[$FieldAssociation, field], Message[FreeLag::UndefinedField, field]; Abort[]];
	If[$FieldAssociation[field, BackgroundField], Message[FreeLag::bkgfield, field]; Abort[]];

	isChiralFermion= $FieldAssociation[field, Chiral] =!= False;

	massInd= GetCouplings[$FieldAssociation[field, Mass]][Indices];
	massDiags= GetCouplings[$FieldAssociation[field, Mass]][DiagonalCoupling];
	fieldInd= $FieldAssociation[field, Indices];
	fieldType= $FieldAssociation[field, Type];

	indK = ConstantArray[a, Length@ fieldInd]/. List-> Sequence;

	If[$FieldAssociation[field, Mass] =!= 0,
		ind1= Sequence@@ ConstantArray[i, Length@ fieldInd];
		ind2= Sequence@@ ConstantArray[j, Length@ fieldInd];
		m= $FieldAssociation[field, Mass]@@ Module[{encountered= {}},
				If[MemberQ[encountered,#], j, AppendTo[encountered, #]; i]&/@ massInd
			];
		(*Light bosons masses with purely diagonal indices are squared*)
		If[$FieldAssociation[field, Type] =!= Fermion && And@@ massDiags, m= m^2];
		(*Insert the group structure to make the mass term invariant. Also includes delts for the diagonal mass indices*)
		m*= Times@@ Table[
				If[!isChiralFermion || MemberQ[Keys@ $FlavorIndices, indType] ||
					$Representations[indType/. Bar-> Identity, Reality] =!= -1
				,
					Delta[Index[i, indType], Index[j, indType]]
				,
					If[Head@ indType === Bar, Identity, Bar]@ CG[eps@ GroupFromRep@ indType, {i, j}]
				]
			, {indType, Complement[fieldInd, Pick[massInd, massDiags, False]]}];
	,
		ind1= indK;
		ind2= indK;
		m= 0;
	];

	If[$FieldAssociation[field, SelfConjugate] === True, 1/2, 1]* Switch[fieldType,
		Scalar,
			(Bar[CD[\[Mu],field[indK]]] CD[\[Mu],field[indK]] - m Bar[field[ind1]] field[ind2]),
		Fermion,
			If[$FieldAssociation[field][Chiral] === False,
				(I Bar[field[indK]]\[CenterDot] \[Gamma][\[Mu]]\[CenterDot] CD[\[Mu], field[indK]] - m Bar[field[ind1]]\[CenterDot] field[ind2])
			,
				(I Bar[field[indK]]\[CenterDot] \[Gamma][\[Mu]]\[CenterDot] CD[\[Mu], field[indK]] - (m/2 Bar[CConj[field[ind1]]]\[CenterDot] field[ind2]+ Bar@ m/2 Bar[field[ind1]]\[CenterDot] CConj[field[ind2]]))
			],
		Vector,
			(*Gauge fields are normalized with their couplings*)
			normalization= First[Query[Select[#@ Field === field&], Key@ Coupling]@ $GaugeGroups, 1&][]^-2;
			(-1/2 normalization* Bar[FS[field,\[Mu],\[Nu],indK]]FS[field,\[Mu],\[Nu],indK] + m Bar[field[\[Mu],ind1]] field[\[Mu],ind2])
		,Ghost|AntiGhost,
			0
	]//Contract//RelabelIndices
];


FreeLag[field_, fields__]:= Plus@@ FreeLag/@ List[field, fields];

(*Defaults to all defined (non-bkg) matter and vector fields*)
FreeLag[]:= Sum[FreeLag@ f, {f, FieldByProperty[BackgroundField-> False, Type-> (Scalar|Fermion|Vector), ZeroMode-> False]}];


(* ::Subsubsection::Closed:: *)
(*Kinetic Lagrangian*)


(* ::Text:: *)
(*This is the part of the free Lagrangian that reproduces the Kinetic piece, \[CapitalDelta], of the fluctuation operator (see [2012.08506] eq. (2.3)).*)


KinOpLagrangian[field_, fields__]:= Plus@@ KinOpLagrangian/@ {field, fields};


KinOpLagrangian[]= 0;
KinOpLagrangian[field_Symbol]:=Module[
		{FieldType, MassInd, FieldInd, indK, ind1, ind2, indM, m, i, j, a, \[Mu], \[Nu],fluctuation},
	FieldInd= $FieldAssociation[field, Indices];
	FieldType= $FieldAssociation[field, Type];

	indK= Sequence@@ Table[a, {n, Length@ FieldInd}];

	If[GetFieldsUpdated[field, Heavy],
		MassInd= DeleteDuplicates@ GetCouplings[GetFieldsUpdated[field, Mass], Indices];
		ind1 = Sequence@@ Table[If[MemberQ[MassInd, FieldInd[[n]]], i, a], {n, Length@ FieldInd}];
		indM = Sequence@@ Table[i, Length@ MassInd];
		ind2 = ind1;
		m= GetFieldsUpdated[field, Mass]@ indM;
	,
		ind1= ind2= indK;
		m= 0;
	];

	If[$FieldAssociation[field, SelfConjugate], 1/2, 1]Switch[FieldType
		,Scalar,
			Bar[BackgroundCD[\[Mu],field[indK]]] BackgroundCD[\[Mu],field[indK]] - m^2 Bar[field[ind1]] field[ind2]
		,Fermion,
			If[$FieldAssociation[field][Chiral]===False,
				I Bar[field[indK]]\[CenterDot] \[Gamma][\[Mu]]\[CenterDot] BackgroundCD[\[Mu],field[indK]] - m Bar[field[ind1]]\[CenterDot] field[ind2]
			,
				I Bar[field[indK]]\[CenterDot] \[Gamma][\[Mu]]\[CenterDot] BackgroundCD[\[Mu],field[indK]] - m/2 (Bar[CConj[field[ind1]]]\[CenterDot] field[ind2]+ Bar[field[ind1]]\[CenterDot] CConj[field[ind2]])
			]
		,Vector,
			If[GaugeFieldQ@ field,
				-GetFieldsUpdated[field, Coupling]  BackgroundCD[\[Mu], field[\[Nu], indK]]^2
			,
				-Bar@ BackgroundCD[\[Mu], field[\[Nu], indK]]  BackgroundCD[\[Mu], field[\[Nu], indK]] + m^2 Bar[field[\[Mu],ind1]]  field[\[Mu],ind2]
			]
		,Ghost|AntiGhost,
			0
	]//Expand//RelabelIndices
];


(* ::Subsection:: *)
(*Routine to load Lagrangians*)


(* ::Subsubsection::Closed:: *)
(*Options*)


OptionTest[LoadModel, SetModelOptions]=
	MatchQ[{{_Symbol, _Symbol, _Symbol-> _}...}| {_Symbol, _Symbol, _Symbol-> _}];


OptionMessage[SetModelOptions, LoadModel, val_]:=
	Message[General::optexpectsval, SetModelOptions, LoadModel, val, "list of sets of the kind {<DefinitionFunction>, <Model symbol>, <Option> -> <New value>}"];


(* ::Subsubsection::Closed:: *)
(*Load model*)


LoadModel::str= "The model name must be a string."


Options[LoadModel]={ModelParameters -> {}, SetModelOptions-> {}, Verbose-> False}


LoadModel[fileName_,opts:OptionsPattern[]]? OptionsCheck:= LoadModel[fileName, opts]= Module[{modelPath},

	(* Check that the file name is a string *)
	If[Head@ fileName =!= String, Message[LoadModel::str]; Abort[]; ];
	
	(* determine path to the model file *)
	modelPath= DetermineModelPath[fileName, OptionValue@ Verbose];
	
	(* loading the model file *)
	ExecuteModelDefinition[modelPath, FilterOptions[opts, ExecuteModelDefinition] ]
]


(* ::Subsubsection::Closed:: *)
(*Determine model path*)


(* ::Text:: *)
(*Determine which model file is indicated by a string and return full path *)


LoadModel::fileext= "The extension \"`1`\" of model name \"`2`\" is not valid. Use .m or .wl (or leave balnk).";
LoadModel::nofiles= "No files matches the model name \"`1`\".";
(*LoadModel::multfiles= "Multiple files `1` matches the model name \"`2`\". Please specify an extension and/or a valid model path.";*)
LoadModel::multfiles= "Multiple files `1` matche the model name \"`2`\". Loading the file `3`. Use the full path if you want to load another file.";


DetermineModelPath[fileName_String, verbose_]:= Block[
		{dir, extension, matches, models, modName, modPossibilities},
	{dir, modName}= Through@ {DirectoryName, FileNameTake}@ fileName;

	modName= StringReplace[modName, RegularExpression["\\.$"]-> ""]; (*Remove any trailing '.'*)
	modPossibilities= Switch[extension= FileExtension@ modName
		,"",
			{modName <> ".m", modName <> ".wl"}
		,"m"|"wl",
			{modName}
		,_,
			Message[LoadModel::fileext, extension, fileName];
			Abort[];
	];

	(*If directory is specified*)
	If[dir =!= "" && !StringStartsQ[dir, "https"],
		matches= FileNames[modPossibilities, dir];
		Switch[Length@ matches
		,0,
			Message[LoadModel::nofiles, fileName];
			Abort[];
		,1,
			If[verbose, Print["Loading model ", First@ matches]];
			Return@ First@ matches;
		,_,
			Message[LoadModel::multfiles, matches, fileName];
			Abort[];
		];
	];

	(*Full file match*)
	models= Normal@ GetModels[Path-> True];
	matches= Position[models, Alternatives@@ modPossibilities];
	If[Length@ matches === 0,
		(*Partial match*)
		matches= Position[models, mod_String/; StringContainsQ[mod, fileName]];
	];

	matches=  FileNameJoin/@ Transpose@
		{Extract[models, {First@ #, 1}&/@ matches], Extract[models, matches]};
	Switch[Length@ matches
	,0,
		If[CheckGitLabAccess[],
			(* if model database is available *)
			DetermineModelPathOnline[fileName, modPossibilities, verbose]
			,
			(* if model database is not available *)
			Message[GetModels::GitLab]; 
			Message[LoadModel::nofiles, fileName];
			Abort[];
		]
		
	,1,
		If[verbose, Print["Loading model ", First@ matches]];
		First@ matches
	,_,
		Echo@ StringForm[LoadModel::multfiles, matches, fileName, Last@ matches];
		If[verbose, Print["Loading model ", Last@ matches]];
		Last@ matches
		(*Abort[];*)
	]

]


DetermineModelPathOnline[fileName_,modPossibilities_,verbose_]:= Module[
	{models,matches}
	,
	models= Normal@ GetModels[Path-> True, "Database"->True];
	matches= Position[models, Alternatives@@ modPossibilities];

	If[Length@ matches === 0,
		(*Partial match*)
		matches= Position[models, mod_String/; StringContainsQ[mod, fileName]];
	];

	matches=  FileNameJoin/@ Transpose@
		{Extract[models, {First@ #, 1}&/@ matches], Extract[models, matches]};
		
	matches= StringReplace[matches,"https:/gitlab"->"https://gitlab"]; (* due to black magic aforwards slash is disappears in the code block above *)
		
	Switch[Length@ matches
	,0,
		Message[LoadModel::nofiles, fileName];
		Abort[];
	,1,
		Print["Loading model from database:\n \[Rule] ",First@ matches];
		First@ matches
	,_,
		Message[LoadModel::multfiles, matches, fileName];
		Abort[];
	]
]


(* ::Text:: *)
(*Get models*)


GetDirectories[]:= Block[{dir, modelDir, notebookDir, subDirs},
	(*Returns $Failed if the notebook is unsaved*)
	notebookDir= Quiet[NotebookDirectory[], NotebookDirectory::nosv];
	modelDir= FileNameJoin@{$MatchetePath, "Models"};
	subDirs= Select[FileNames[RegularExpression@"\\w*", modelDir], DirectoryQ]; (*Subdirectories from Models*)
	
	DeleteDuplicatesBy[Association@@ Flatten@ {
		(* the ordering determines which models are red when multiple matches exist *)
		FileNameJoin@ {"Matchete", "Models"}-> modelDir,
		Table[
			FileNameJoin@ {"Matchete", FileNameTake[dir, -2]}-> dir
		, {dir, subDirs}],
		"Directory"-> Directory[],
		If[notebookDir =!= $Failed, "NotebookDirectory"-> NotebookDirectory[], Nothing]
	}, FileNameDrop[#, 0]&]
];


Options@ GetModels= {Path-> False, "Database"-> False, Verbose->True};


GetModels[OptionsPattern[]]? OptionsCheck:= Block[{models, dirs},
	If[OptionValue["Database"],
		(* online *)
		OptionalMonitor[OptionValue@Verbose,
			models= GetModelDatabase[];
			dirs= (#->("https://gitlab.com/matchete/model-database/-/raw/master/"<>#))&/@Keys[models],
		"Searching the model database at https://gitlab.com/matchete/model-database."]
	,
		(* offline *)
		dirs= GetDirectories[];
		models= DeleteCases[(FileNameTake/@ FileNames[{"*.m", "*.wl"}, #]&)/@ dirs, {}];
	];
	If[OptionValue@ Path,
		models= KeyMap[#/. dirs&, models];
	];
	models
]


Quiet[NotebookDirectory[],{NotebookDirectory::nosv}]


(* ::Text:: *)
(*Get models from database*)


(* cach results since api access is rather slow *)
GetModelDatabase[]:= GetModelDatabase[]= Module[
	{
		repoURL= "https://gitlab.com/matchete/model-database",
		encodedPath, projectInfo, projectID, ListFiles
	}
	,
	(* check internet and GitLab access  *)
	If[!CheckGitLabAccess[],
		Message[GetModels::GitLab]; Abort[]
	];
	
	(* Convert repository path to the form expected by the GitLab API *)
	encodedPath= URLEncode@StringDrop[repoURL, StringLength["https://gitlab.com/"]];
	
	(* --- GET PROJECT ID --- *)
	projectInfo= Import["https://gitlab.com/api/v4/projects/" <> encodedPath, "RawJSON"];
	projectID= projectInfo["id"];
	
	(* --- FUNCTION TO RECURSE THROUGH REPO TREE --- *)
	ListFiles[path_]:= Module[
		{tree,subdirs,files},
		tree=    Import["https://gitlab.com/api/v4/projects/"<>ToString[projectID]<>"/repository/tree?path="<>URLEncode[path]<>"&per_page=100","RawJSON"];
		subdirs= Select[tree,#["type"]==="tree"&];
		files=   Select[tree,#["type"]==="blob"&];
		Flatten@Join[
			{
				(*recurse into subdirectories*)
				ListFiles[FileNameJoin[{path,#["name"]}/.""->Nothing]]&/@subdirs,
				(*return matching .m or.wl files*)
				path-> (FileNameJoin[{#["name"]}]&/@Select[files,StringMatchQ[#["name"],___~~(".m"|".wl")]&])
			}/.{}->Nothing/.(_->Nothing)->Nothing
		]
	];
	
	Association@@ListFiles[""]
]


(* Returns True if GitLab is accessible and false otherwise *)
CheckGitLabAccess[]:= Quiet@Check[URLRead["https://gitlab.com/matchete/model-database"],$Failed]=!=$Failed


GetModels::GitLab= "Could not connect to https://gitlab.com/matchete/model-database.";


(* ::Subsubsection::Closed:: *)
(*ReadFileExpression*)


ReadFileExpression@ file_:= Module[{out},
	out= ReadList[file, Hold[Expression]];
		(*Catches Wolfram NonCommutativeMultiply for backwards compatibility*)
	out= out/. HoldPattern@ NonCommutativeMultiply-> NCM;
	DeleteCases[out, Hold@ Null]
]


(* ::Subsubsection::Closed:: *)
(*Read and execute model file*)


(* ::Text:: *)
(*Read model file, change parameters and options accordingly, and execute definition.*)


LoadModel::defsymb      = "The symbol `1` is already in use and cannot be used as a new model parameter.";
LoadModel::defaultpar   = "The LHS of the default value `1` from the (parent) model file must be a Symbol.";
LoadModel::modelpar     = "The parameter `1` does not appear in the model file and cannot be renamed.";
LoadModel::localvar     = "The redefined symbol(s) \"`1`\" appear as local variables in Module or Block in the (parent) model. This redefinition is likely to cause an error and is not allowed.";
LoadModel::nulllag      = "The Lagrangian returned by the model file `1` is Null (Parent file level `2`). Ensure that the Lagrangian expression is the last statement in the model file and that it DOES NOT terminate with a semicolon ' ; '.";
LoadModel::optscall     = "The error happens with the Option call SpecifyOptions\[Rule] {\[Ellipsis], `1` ,\[Ellipsis]}.";
LoadModel::parcall      = "The error happens with the Option call ModelParameters\[Rule] {\[Ellipsis], `1` ,\[Ellipsis]}.";
LoadModel::unkwnfunc    = "The model file does not define `1`[`2`, \[Ellipsis]], which you are attempting to change options for. Perhaps you changed the model parameter `2`.";
LoadModel::unkwnopt     = "Attempt to modify the option `1` for a definition with `2`. This option does not exist. ";
LoadModel::unkwnoptfunc = "Cannot change options associated with function `1`. SetModelOptions should only be used to change options associated with functions `2`.";


Options@ ExecuteModelDefinition= {
	ModelParameters-> {},
	SetModelOptions-> {}
};


ExecuteModelDefinition[file_String, OptionsPattern[]]:= Module[
		{change, changedPars, childLag, defaultParams, lag, modelDef, modelDir, modelOpts, modelLocalVars, n, overlap,
		parameterSubs, symb},
	modelDir= DirectoryName@ file;

	(*Load held Model definition*)
	modelDef= ReadFileExpression@ file;

	modelDef= ReadParentModels[modelDef, modelDir];

	(*Determine which are the local variables in the model file(s)*)
	modelLocalVars= Join@@ Cases[modelDef, HoldPattern@ (Module|Block|With)[pars_List, _]:> pars, All]// DeleteDuplicates;

	(*Check parameter changes*)
	parameterSubs= OptionValue[ModelParameters];
	Do[
		If[FreeQ[modelDef, symb[[1]]],
			Message[LoadModel::modelpar, symb[[1]]];
			Message[LoadModel::parcall, symb];
			Abort[];
		];
		(* the checks below are problematic in light of parameter name overloading done consistently between UV and EFT *)
		(*
		If[MatchQ[symb[[2]], _Symbol] && Defined@ symb[[2]],
			Message[LoadModel::defsymb, symb[[2]]];
			Message[LoadModel::parcall, symb];
			Abort[];
		];
		*)

		overlap= Intersection[modelLocalVars, List@@ symb];
		If[Length@ overlap =!= 0,
			Message[LoadModel::localvar, StringJoin@@ Riffle[ToString/@ overlap, ", "] ];
			Message[LoadModel::parcall, symb];
		];
	, {symb, parameterSubs}];

	(*Remove Niceforms where parameter changes are performed*)
	changedPars= parameterSubs[[;;, 1]];
	modelDef= modelDef/. {
			HoldPattern@ DefineCoupling[p_, l___, Rule[NiceForm, _], r___]/; MemberQ[changedPars, p]:>
				DefineCoupling[p, l, r],
			HoldPattern@ DefineField[p_, l___, Rule[NiceForm, _], r___]/; MemberQ[changedPars, p]:>
				DefineField[p, l, r],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, coup_, field_], l___, Rule[NiceForm, _], r___]/; SubsetQ[changedPars, {coup, field}]:>
				DefineGaugeGroup[s, l, r],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, _, field_], l___, Rule[NiceForm, {nfCoup_, _}], r___]/; MemberQ[changedPars, field]:>
				DefineGaugeGroup[s, l, r, NiceForm-> {nfCoup, Default}],
			HoldPattern@ DefineGaugeGroup[s:PatternSequence[_, _, coup_, _], l___, Rule[NiceForm, {_, nfField_}], r___]/; MemberQ[changedPars, coup]:>
				DefineGaugeGroup[s, l, r, NiceForm-> {Default, nfField}]
		};

	(*Incorporate defualt parameters*)
	defaultParams= Cases[modelDef, ParameterDefault[par_-> val_]:>
		If[MatchQ[par, _Symbol],
			par-> val
		,
			Message[LoadModel::defaultpar, FullForm@ par-> val];
			Abort[];
		]
		, All];
	modelDef= DeleteCases[modelDef, _ParameterDefault, All];

	(*Use value from ModelParameters or the value from the most junior model that provides a value*)
	parameterSubs= DeleteDuplicatesBy[Join[parameterSubs, Reverse@ defaultParams], First];

	(*Change parameter names*)
	modelDef= modelDef/. parameterSubs;

	(*Set/change specific options*)
	modelOpts= OptionValue@ SetModelOptions;
	If[MatchQ[modelOpts, {_Symbol, __}], modelOpts= {modelOpts}];

	Do[
		With[{func= change[[1]], var= change[[2]], opt= change[[3, 1]], newVal= change[[3, 2]]},
			If[!MemberQ[$setModelOptionsFuncs, func],
				Message[LoadModel::unkwnoptfunc, func, $setModelOptionsFuncs];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			If[FreeQ[modelDef, HoldPattern@ func[var, ___]],
				Message[LoadModel::unkwnfunc, func, var];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			If[!MemberQ[Options[func][[;;, 1]], opt],
				Message[LoadModel::unkwnopt, opt, func];
				Message[LoadModel::optscall, change];
				Abort[];
			];

			modelDef= If[FreeQ[modelDef, HoldPattern@ func[var, ___, Rule[opt, _], ___]],
					modelDef/. HoldPattern@ func[var, rest___]:> func[var, rest, opt-> newVal]
				,
					modelDef/. HoldPattern@ func[var, rest1___, Rule[opt, _], rest2___]:>
						func[var, rest1, opt-> newVal, rest2]
				];
		]
	, {change, modelOpts}];

	(*Run model definition over consecutive child models*)
	lag= 0;
	
	Do[
		childLag= Last[ ReleaseHold@ modelDef[[n]], 0];
		If[childLag === Null,
		Message[LoadModel::nulllag, FileBaseName@ file, Length@ modelDef- n];
			Abort[];
		];
		lag+= childLag;
	, {n, Length@ modelDef}];
	
	lag
]


(* ::Text:: *)
(*Valid functions to change with SetModelOptions*)


$setModelOptionsFuncs= {DefineCoupling, DefineField, DefineFlavorIndex, DefineGaugeGroup, DefineGlobalGroup};


(* ::Subsubsection::Closed:: *)
(*Determine parent model*)


(* ::Text:: *)
(*Recursively inserts the model definition from ParentModel files into the Model definition *)


LoadModel::findparent= "The parentmodel '`1`' could not be found.";
LoadModel::parentarg= "The argument/option `1` in `2` is invalid.";
LoadModel::parentopt= "The option \"UseLagrangian\" in `1` takes a boolean value.";
LoadModel::parfile= "The model in ParentModel[`1`] must be a String.";
LoadModel::parmoduse= "ParentModel[<model>] must appear as the first expression in the model file.";


(*This is a faux option given that ParentModel never evaluates. It nevertheless is interpreted as an option in LoadModel*)
Options@ ParentModel= {
	"UseLagrangian"-> True
};


ReadParentModels[modelDef_List, modelDir_String, includeParentLag_:True]:= Module[
		{childDef, parentDef, parentFile, parentOps, parentPath, useParentLag= includeParentLag},
	If[FreeQ[modelDef, _ParentModel],
		Return@ {modelDef};
	];

	parentFile= First@ modelDef;
	childDef= Rest@ modelDef;

	{parentFile, parentOps}= Switch[parentFile
	, Hold[ParentModel[_, ___] ],
		{parentFile[[1, 1]], List@@ parentFile[[1, 2;;]]}
	, Hold[CompoundExpression[ParentModel[_, ___], Null]], (* If file reads ParentModel[...]; *)
		{parentFile[[1, 1, 1]], List@@ parentFile[[1, 1, 2;;]]}
	, _,
		Message[LoadModel::parmoduse];
		Abort[];
	];

	If[!MatchQ[parentFile, _String],
		Message[LoadModel::parfile, parentFile];
		Abort[];
	];

	(*Check parent options*)
	Do[
		Switch[opt
		, "UseLagrangian"->  _?BooleanQ,
			useParentLag= useParentLag && opt[[2]];
		, "UseLagrangian"->  _,
			Message[LoadModel::parentopt, ParentModel[parentFile, Sequence@@ parentOps]];
			Abort[];
		, _,
			Message[LoadModel::parentarg, opt, ParentModel[parentFile, Sequence@@ parentOps]];
			Abort[];
		];
	, {opt, parentOps}];

	(*Read parent model definitions*)
	(*Also check for model path with the directory of the child model*)
	parentPath= CheckAbort[Quiet@DetermineModelPath[FileNameJoin@ {modelDir, parentFile}, False],
			CheckAbort[Quiet@DetermineModelPath[parentFile, False],
				Message[LoadModel::findparent, parentFile];
				Abort[];
			]
		];

	parentDef= ReadFileExpression@ parentPath;

	(*Return 0 if Lagrangian should not be used. This way maintains all definitions that may
		be performed in the last file expression, e.g., in a module.*)
	If[!useParentLag,
		AppendTo[parentDef, 0];
	];

	(*Recursively check for parents*)
	parentDef= ReadParentModels[parentDef, parentPath, useParentLag (*don't include lags if disabled further down the tree*)]; 
	childDef= ReadParentModels[childDef, modelDir]; (*Catches additional parrent files from the child*)

	Join[parentDef, childDef]
]
