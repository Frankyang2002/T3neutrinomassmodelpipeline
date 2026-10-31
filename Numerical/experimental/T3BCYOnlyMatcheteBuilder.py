"""Build a complete experimental Matchete Y-only matching program from the installed native probe.

Does not modify or evaluate the authoritative RunGeneralT3EFT1ScalarThresholdProbe.wl.
The generated standalone program resides alongside RunMatching.wl. Every source
edit is anchored in the inspected native version, with fail-closed drift guards.

IMPORTANT: A source-level adaptation, not a verified live Matchete calculation.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path

BASE_NAME = 'RunGeneralT3EFT1ScalarThresholdProbe.wl'
GENERATED_NAME = 'RunT3BCYOnlyMatching.experimental.wl'
FORMAT = 't3-bc-y-only-matchete-source-adaptation-v1'


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError('Y-only Matchete source: ' + message)


def replace_once(source: str, old: str, new: str, label: str) -> str:
    n = source.count(old)
    require(n == 1, f'{label}: expected one anchor, found {n}')
    return source.replace(old, new, 1)


def adapt(source: str) -> str:
    require(len(source) > 80000, 'native program unexpectedly short')
    require('T3MELeadingC5OneLoopProjection' in source and
            'RunT3Matching[scalarEFT,5,1]' in source and
            'T3MEVerifiedRenormalizableReconstruction' in source,
            'missing native matching and verification machinery')
    s = source
    s = replace_once(s,
        'checkHeavyX=Quiet@Check[\n   DefineField[NewScalarX,Scalar,SelfConjugate->False,\n     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Heavy,MSX}];True,$Failed];',
        'checkHeavyX=Quiet@Check[\n   DefineField[NewScalarX,Scalar,SelfConjugate->False,\n     Indices->{T3BSMd2},Charges->{U1Y[-1/2]},Mass->{Light,MSX}];True,$Failed];',
        'retain-X-as-light')
    # The source already removes and redefines both fields, and forms their
    # free Lagrangian under the new field registry. Only Y remains heavy.
    s = replace_once(s,
        'StringContainsQ[completeText,"Field[NewScalarX"] ||\n    StringContainsQ[completeText,"Field[NewScalarY"] ||',
        '!StringContainsQ[completeText,"Field[NewScalarX"] ||\n    StringContainsQ[completeText,"Field[NewScalarY"] ||',
        'matched-EFT-has-X-no-Y')
    s = replace_once(s,
        '{"Field[NewScalarX","Field[NewScalarY","Field[NewFermion"}',
        '{"Field[NewScalarY","Field[NewFermion"}',
        'diagnostic-heavy-field-list')
    s = replace_once(s,
        '"heavy_scalar_fields"->{"NewScalarX","NewScalarY"}',
        '"heavy_scalar_fields"->{"NewScalarY"},\n    "light_scalar_fields"->{"H","NewScalarX"}',
        'report-field-inventory')
    s = replace_once(s,
        '"no_scalar_fields_in_matched_eft"->True,',
        '"heavy_Y_absent_and_light_X_retained_in_matched_eft"->True,',
        'report-output-field-inventory')
    s = replace_once(s,
        '"format"->"t3-mass-eigen-eft1-scalar-threshold-probe-v1",',
        '"format"->"t3-bc-y-only-matchete-eigenfield-probe-v1",',
        'report-format')
    # Two occurrences: returned report and final status check. Exact count is important.
    old_status='ScalarThresholdOneLoopRawSectorsExtractedNoEFT1Running'
    require(s.count(old_status)==2,
            f'unexpected native stage-status count {s.count(old_status)}')
    s=s.replace(old_status,'YOnlyOneLoopRawEFT2SectorsExtractedNotPhysicalMatching')
    # In the inherited probe only the opt-in leading-C5 subtheory explicitly
    # re-based the light scalar free terms before matching. For a Y-only FULL
    # EFT run we must also re-base the free Lagrangian (Y changed classification),
    # then require CheckLagrangian to pass before invoking Matchete.
    s = replace_once(s,
        ' T3MECheckpoint["scalar-threshold-EFT1-lagrangian-check-started"];',
        '  If[!TrueQ[$T3MELeadingC5Only],\n'
        '    T3MECheckpoint["Y-only-full-EFT-free-rebase-started"];\n'
        '    scalarEFT=Expand[scalarEFT-lightEigenFree+heavyEigenFree] // RelabelIndices;\n'
        '    If[scalarEFT===$Failed || scalarEFT===$Aborted ||\n'
        '       FreeQ[scalarEFT,NewScalarX] || FreeQ[scalarEFT,NewScalarY] ||\n'
        '       FreeQ[scalarEFT,GammaCC],\n'
        '      T3MEFail["YOnlyFullEFTFreeLagRebaseInvalid"]];\n'
        '    T3MECheckpoint["Y-only-full-EFT-free-rebase-returned"]];\n'
        ' T3MECheckpoint["scalar-threshold-EFT1-lagrangian-check-started"];',
        'full-EFT-heavy-Y-free-lagrangian-rebase')
    # Program's existing LLHH filtering is valid as an isolated operator selection.
    # Additionally extract LLXX and LLHX operators that remain in EFT2.
    anchor=' <|\n   "format"->"t3-bc-y-only-matchete-eigenfield-probe-v1",'
    addition=r'''  (* Additional EFT2 operator inventory after native Y-only matching.
     Only terms containing two light leptons and exactly two H/X scalars
     are exported; the complete matched EFT2 is separately preserved above.
     These are raw one-loop operators and are not renormalized Wilson seeds. *)
  $T3BCYOnlySectorCounts=<||>;
  Module[{baseTerms, labels, eligible, part, nH, nX, filePath, label},
    baseTerms=Select[deltaTerms,Function[term,
      With[{s=ToString[InputForm[term]]},
        StringContainsQ[s,"GammaCC"] &&
        StringCount[s,"Field[l, Fermion"]===2 &&
        !StringContainsQ[s,"Field[NewScalarY"] &&
        !StringContainsQ[s,"Field[NewFermion"]]]];
    labels={{"LLHH",2,0},{"LLHX",1,1},{"LLXX",0,2}};
    Do[
      label=entry[[1]];nH=entry[[2]];nX=entry[[3]];
      eligible=Select[baseTerms,Function[term,
        With[{s=ToString[InputForm[term]]},
          StringCount[s,"Field[H, Scalar"]===nH &&
          StringCount[s,"Field[NewScalarX, Scalar"]===nX]]];
      part=If[eligible==={},0,Total[eligible]];
      filePath=FileNameJoin[{DirectoryName[outPath],
        FileBaseName[outPath]<>"_raw_"<>ToLowerCase[label]<>".inputform.txt"}];
      If[Quiet@Check[Export[filePath,ToString[InputForm[part]],"Text"],$Failed]===$Failed,
        T3MEFail["YOnlyEFT2OperatorSectorExportFailed",label]];
      $T3BCYOnlySectorCounts[label]=Length[eligible],
      {entry,labels}];
    Print["Y-only EFT2 raw Wilson sector counts: ",$T3BCYOnlySectorCounts]
  ];
'''
    s = replace_once(s, anchor, addition+anchor, 'additional-EFT2-operator-extraction')
    s = replace_once(s,
        '"raw_scalar_llhh_term_count"->Length[scalarCandidates],',
        '"raw_scalar_llhh_term_count"->Length[scalarCandidates],\n    "raw_EFT2_operator_sector_counts"->$T3BCYOnlySectorCounts,\n    "Y_only_eigenfield_matching_experimental"->True,\n    "renormalized_sequential_matching_certified"->False,',
        'report-operator-counts')
    # Avoid falsely describing generated EFT2 files as SMEFT. Keeping the
    # filename and content-compatible for downstream debugging is intentional.
    require('Mass->{Heavy,MSY}' in s and
            'Mass->{Heavy,MSX}];True,$Failed]' not in s and
            'Mass->{Light,MSX}];True,$Failed]' in s,
            'incorrect heavy/light field reclassification')
    require('!StringContainsQ[completeText,"Field[NewScalarX"]' in s and
            '"heavy_scalar_fields"->{"NewScalarY"}' in s,
            'missing EFT2 survivor checks')
    return ('(* GENERATED EXPERIMENTAL Y-ONLY MATCHING FROM INSTALLED NATIVE PROBE.\n'
            '   This source is NOT a certified finite threshold coefficient.\n'
            '   Source fidelity and exact edits are checked by the Python builder. *)\n') + s


def build(source_path: Path, target_path: Path, *, overwrite: bool=False) -> dict:
    source_path=Path(source_path).resolve()
    target_path=Path(target_path).resolve()
    require(source_path.is_file(), 'missing native source '+str(source_path))
    require(source_path.name==BASE_NAME, 'unexpected source basename')
    require(target_path.name==GENERATED_NAME, 'unexpected generated basename')
    require(source_path.parent==target_path.parent,
            'generated program must be next to original and RunMatching.wl')
    require(source_path!=target_path, 'never replace authoritative Matchete source')
    code=adapt(source_path.read_text(encoding='utf-8-sig'))
    if target_path.exists():
        require(target_path.read_text(encoding='utf-8')==code,
                'existing generated file differs; resolve before regeneration')
    else:
        target_path.write_text(code,encoding='utf-8')
    return {
        'format':FORMAT,
        'native_source_sha256':sha256(source_path.read_bytes()).hexdigest(),
        'generated_source_sha256':sha256(target_path.read_bytes()).hexdigest(),
        'source':str(source_path),'generated':str(target_path),
        'source_edits_fail_closed':True,'physical_matching_certified':False,
        'full_EFT1_to_EFT2_one_loop_numeric_matching_completed':False,
    }
