(* ::Package:: *)

InstallMatchete::version="Warning: `1` has only been tested on Mathematica versions \[GreaterEqual] `2` and you are using Mathematica `3`";


(* in-house solution for copying directories since Mathematica does not allow overwriting *)
BetterCreateDirectory[dir_]:=If[FileExistsQ[dir], dir, CreateDirectory[dir]]

BetterCopyDirectory[from_,to_]:=Module[{source,relativeDirs,out, files, targetfiles},

	source= FileNames[All,from,Infinity];
	relativeDirs = FileNameDrop[#,FileNameDepth[from]]&/@Select[source,DirectoryQ];
	relativeDirs = FileNameJoin[{to,#}]&/@relativeDirs;

	(* create the directory structure *)
	out = BetterCreateDirectory/@PrependTo[relativeDirs,to];
	If[MemberQ[out,$Failed],Return[$Failed]];

	(* get file names and create target *)
	files = Select[source,!DirectoryQ[#]&];
	targetfiles=FileNameJoin[{to,#}]&/@(FileNameDrop[#,FileNameDepth[from]]&/@ files);

	out=CopyFile[#[[1]],#[[2]],OverwriteTarget->True]&/@( {files,targetfiles}\[Transpose]);
	If[MemberQ[out,$Failed],Return[$Failed],Return[out]]
]


InstallMatchete[]:=Module[{yn,packageName,packageDir,MinVersion,MatcheteLink,QuestionOverwrite,tmpFile,unzipDir,zipDir},

	(* Definitions *)
	packageName="Matchete";
	
	If[Head[Matchete`$MatchetePath]===String,
		(* if this is run from Matchete's update function, retrieve the path so we install in the right location *)
		packageDir = Matchete`$MatchetePath
	,
		packageDir = FileNameJoin[{$UserBaseDirectory,"Applications","Matchete"}]
	];
	
	MinVersion=11.0;
	MatcheteLink="https://gitlab.com/matchete/matchete/-/archive/master/matchete-master.zip";

	(* Messages *)
	QuestionOverwrite="Matchete is already installed. Do you want to replace the content of "<>packageDir<>" with a newly downloaded version?";
	
	(* Check Mathematica version *)
	If[$VersionNumber<MinVersion,
		Message[InstallMatchete::version,packageName,ToString[MinVersion],$VersionNumber];
	];
	
	(* Check if Matchete has already been installed *)
	If[
		DirectoryQ[packageDir],
		If[
			ChoiceDialog[QuestionOverwrite,{"Yes"->True,"No"->False},WindowFloating->True,WindowTitle->"Matchete installation detected"],
			Null,
			Abort[]
		];
	];

	(* Download Matchete *)
	Print["Downloading Matchete from ",MatcheteLink];

	tmpFile=Quiet@URLSave[MatcheteLink];

	If[tmpFile===$Failed,
		Print["Failed to download Matchete.\nInstallation aborted!"];
		Abort[]
	];

	(* Unzip Matchete file *)
	Print["Extracting Matchete zip file"];

	unzipDir=tmpFile<>".dir/";
	ExtractArchive[tmpFile,unzipDir];

	(* Move files to the Mathematica packages folder *)
	Print["Copying Matchete to "<>packageDir];

	zipDir=FileNames["README.md",unzipDir,Infinity];
	Echo@zipDir;
	If[
		BetterCopyDirectory[DirectoryName[zipDir[[1]]],packageDir]===$Failed, 
		"The installation failed. Make sure you can access the target directory."
	];

	(* Delete the extracted archive *)
	Quiet@DeleteDirectory[unzipDir,DeleteContents->True];
	(* link documentation *)
	LinkDocu[packageDir];
	
	(*Remove development folders*)
	Quiet@ DeleteDirectory[FileNameJoin[{packageDir, "DocumentationSource"}],   DeleteContents-> True];
	Quiet@ DeleteDirectory[FileNameJoin[{packageDir, "DocumentationPreBuilt"}], DeleteContents-> True];
	Quiet@ DeleteDirectory[FileNameJoin[{packageDir, "Validation"}],            DeleteContents-> True];
	Quiet@ DeleteDirectory[FileNameJoin[{packageDir, "Package", "DevTools"}],   DeleteContents-> True];
	
	Print["Installation complete!"];
];


InstallMatchete::noDocu= "No pre-built documentation was found.";


LinkDocu[packageDir_]:=Module[{},
	(* check that pre-built documentation is present *)
	If[!DirectoryQ[FileNameJoin@{packageDir,"DocumentationPreBuilt"}],
		Message[InstallMatchete::noDocu];
		Return[]
	];
	
	(* if already present delete current documentation directory *)
	Quiet[
		DeleteDirectory[FileNameJoin@{packageDir,"Documentation"}, DeleteContents-> True]
		,{DeleteDirectory::nodir}
	];
	
	(* copy pre-built documentation *)
	CopyDirectory[
		FileNameJoin@{packageDir,"DocumentationPreBuilt"},
		FileNameJoin@{packageDir,"Documentation"}
	];
	
	(* load matchete *)
	Get[FileNameJoin[{packageDir,"Kernel","init.m"}]];
	
	(* link the new documentation *)
	Matchete`PackageScope`LinkDocumentation[];
]


InstallMatchete[];
