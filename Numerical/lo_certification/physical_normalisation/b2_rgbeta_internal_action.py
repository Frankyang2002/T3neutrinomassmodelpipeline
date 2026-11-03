"""B2: validate a fresh RGBeta *internal tensor* prefactor witness.

This is not a native component-axis witness, action density, Feynman rule,
Matchete-to-RGBeta coupling map, or a physical-source certificate.
All source files and native JSONs are read-only. A conditional native probe
is required; no analytical reference JSON can impersonate live evidence.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any
import json
import re

from b2_gauge_coupling import AuditBlocked, CLASSES, SOURCE_BLOBS, git_result, models, require
from b2_rgbeta_ward_witness import source_provenance

FORMAT = 'b2-rgbeta-internal-normalisation-native-v1'
OPERATOR_ORDER = ('y1', 'y2', 'MF', 'lambdaT3', 'lambdaH1', 'lambdaH2')
EXPECTED_RATIOS = {'y1': 'Sqrt[2]', 'y2': 'Sqrt[2]', 'MF': '2',
                   'lambdaT3': '6', 'lambdaH1': '6', 'lambdaH2': '6'}
INTERNAL_SEMANTICS = {
    'y1': 'AddYukawa complex scalar real-component convention',
    'y2': 'AddYukawa complex scalar real-component convention',
    'MF': 'AddFermionMass symmetrised fermion tensor',
    'lambdaT3': 'AddQuartic four complex scalar real-component convention',
    'lambdaH1': 'AddQuartic four complex scalar real-component convention',
    'lambdaH2': 'AddQuartic four complex scalar real-component convention',
}
ALLOWED_B2_PATHS = (
    'RGE/running/rgbeta/T3RGBetaModel.wl',
    'Lagrangian/interactions/T3Topology.wl',
    'Lagrangian/model/SU2Invariants.wl',
)


def read_report(path: Path) -> dict:
    require(path.is_file(), f'missing native output {path}')
    require(path.stat().st_size > 0,
            f'empty native output {path}; the Wolfram export did not complete')
    try:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise AuditBlocked(f'malformed/incomplete native JSON {path}: {exc}') from exc
    require(isinstance(data, dict), 'native output must be JSON association')
    return data


def native_snapshot_guard(data: dict, root: Path, provenance: dict) -> None:
    """Bind a report to current Git HEAD and canonical tracked source bytes.

    This is consistency of native capture with current checkout, NOT a
    cryptographic attestation of code execution or installed RGBeta version.
    """
    require(data.get('git_head_at_capture') == provenance['repository_head_measured'],
            'native report was captured at a different Git HEAD')
    blobs=data.get('git_source_blobs_at_capture')
    require(isinstance(blobs,dict) and set(blobs)==set(ALLOWED_B2_PATHS),
            'native JSON missing or unexpected source blob inventory')
    for rel in ALLOWED_B2_PATHS:
        require(blobs[rel]==SOURCE_BLOBS[rel],
                'stale/mismatched native source blob: '+rel)
        rc,current=git_result(root,'hash-object',f'--path={rel}','--',rel)
        require(rc==0 and current==blobs[rel],
                'native source snapshot does not match local checkout: '+rel)


def expected_fields(family: str, coupling: str, majorana: bool) -> str:
    if coupling=='y1':return '{Bar[S1], l, F}'
    if coupling=='y2':return '{S2, l, '+('F' if majorana else 'Fc')+'}'
    if coupling=='MF':return '{F, '+('F' if majorana else 'Fc')+'}'
    if coupling=='lambdaT3':return '{H, H, S1, Bar[S2]}'
    if coupling=='lambdaH1':return '{Bar[H], H, Bar[S1], S1}'
    if coupling=='lambdaH2':return '{Bar[H], H, Bar[S2], S2}'
    raise AuditBlocked(f'unknown coupling {coupling}')


def normalized_expr(value: str) -> str:
    require(isinstance(value,str) and bool(value.strip()), 'missing native expression')
    # Do not evaluate untrusted Wolfram expressions. Only remove whitespace.
    return re.sub(r'\s+','',value)


def check_row(row: dict, family: str, coupling: str, maj: bool) -> dict:
    require(isinstance(row,dict) and row.get('coupling')==coupling,
            f'{family}: missing/reordered native {coupling}')
    kind=('yukawa' if coupling in ('y1','y2') else
          'mass' if coupling=='MF' else 'quartic')
    require(row.get('kind')==kind and row.get('status')=='NativeRegistryRatioMeasured',
            f'{family} {coupling}: native registry ratio not measured')
    require(normalized_expr(row.get('registered_fields',''))==
            normalized_expr(expected_fields(family,coupling,maj)) and
            row.get('actual_expected_fields_match') is True,
            f'{family} {coupling}: F/Fc, barred field or order changed')
    require(row.get('chirality')==('not_applicable' if kind=='quartic' else 'Right'),
            f'{family} {coupling}: unexpected chirality')
    require(row.get('self_conjugate')==
            ((coupling!='lambdaT3') if kind=='quartic' else None),
            f'{family} {coupling}: unexpected self-conjugacy or h.c. convention')
    require(row.get('native_exact_ratio_and_metadata_verified') is True and
            row.get('physical_Lagrangian_coefficient_verified') is False,
            f'{family} {coupling}: no exact native ratio / false physical claim')
    exp=EXPECTED_RATIOS[coupling]
    require(normalized_expr(row.get('expected_internal_ratio_inputform','')) == exp,
            f'{family} {coupling}: reference internal factor changed')
    # Native evaluated ratio strings are diagnostic. Admit exactly the simple
    # canonical InputForm constants; do NOT sympify arbitrary Wolfram input.
    for field in ('sample_Coupling_over_Matrix_inputform',
                  'sample_CouplingBar_over_Matrix_inputform'):
        require(normalized_expr(row.get(field,''))==exp,
                f'{family} {coupling}: incorrect measured {field}')
    for field in ('sample_Coupling_inputform','sample_CouplingBar_inputform',
                  'registered_invariant_inputform'):
        require(isinstance(row.get(field),str) and row[field].strip() and
                '$Failed' not in row[field] and 'Missing[' not in row[field],
                f'{family} {coupling}: malformed native sample {field}')
    # RGBeta stores flavour-indexed Yukawas and masses as Matrix[...] but
    # zero-flavour-index quartics are scalar symbols.  Check the exact native
    # InputForm seen for each operator, rather than requiring Matrix[ for all.
    # This checks both coupling orientation and the explicit internal factor;
    # the ratio checks above remain independent guards.
    if kind == 'quartic':
        expected_coupling = f'{exp}*{coupling}'
        expected_bar = (f'{exp}*Bar[{coupling}]' if coupling == 'lambdaT3'
                        else expected_coupling)
    else:
        expected_coupling = f'{exp}*Matrix[Bar[{coupling}]][b4p,b4r]'
        expected_bar = f'{exp}*Matrix[{coupling}][b4p,b4r]'
    require(normalized_expr(row['sample_Coupling_inputform']) == expected_coupling and
            normalized_expr(row['sample_CouplingBar_inputform']) == expected_bar,
            f'{family} {coupling}: unexpected native coupling sample or conjugation')
    return {'coupling':coupling,'kind':kind,'internal_tensor_prefactor':exp,
            'sample_Coupling_inputform':row['sample_Coupling_inputform'],
            'sample_CouplingBar_inputform':row['sample_CouplingBar_inputform'],
            'registered_fields':row['registered_fields'],
            'source_of_ratio':'native-evaluated RGBeta coupling registry functions',
            'physical_action_factor_from_this_ratio':None,
            'physical_action_verified':False}


def validate_native_data(data: dict) -> dict:
    """Pure JSON validator: can be tested with adversarial fixtures."""
    require(data.get('format')==FORMAT and
            data.get('status')=='NativeInternalTensorPrefactorsCapturedPhysicalActionPending',
            'wrong native internal-witness format/status')
    require(data.get('internal_prefactor_is_physical_action_factor') is False and
            data.get('numeric_component_axes_verified') is False and
            data.get('physical_source_certified') is False,
            'false physical-source, component, or action claim')
    classes=data.get('classes')
    require(isinstance(classes,list) and len(classes)==5,'not all five native classes')
    output={}
    for row in classes:
        require(isinstance(row,dict),'nonrecord native family')
        family=row.get('family')
        require(family in CLASSES and family not in output,
                'duplicate/unknown native family')
        d1,d2,df,alphas=CLASSES[family]
        dim=row.get('dimensions_alpha')
        require(isinstance(dim,list) and dim[:3]==[d1,d2,df] and
                dim[3] in alphas, f'{family}: incorrect representative')
        maj=bool(dim[3]==-1 and df%2==1)
        require(row.get('majorana') is maj and
                row.get('status')=='AllSixRegistryRatiosMeasured' and
                row.get('physical_action_certified') is False,
                f'{family}: incomplete/overclaimed native family')
        operators=row.get('operators')
        require(isinstance(operators,list) and len(operators)==6,
                f'{family}: expected six internal operator rows')
        extracted=[check_row(op,family,name,maj) for op,name in
                   zip(operators,OPERATOR_ORDER)]
        output[family]={'representative':dim,'majorana_representative':maj,
                        'operator_internal_conventions':extracted,
                        'physical_action_certified':False}
    require(set(output)==set(CLASSES),'native family coverage incomplete')
    return {'families':output, 'families_checked':5,
            'native_registry_ratios_measured':30,
            'physical_coupling_ratios_deducible_from_registry_only':False,
            'physical_source_certified':False}


def expand_all16(validated: dict) -> list[dict]:
    """Transport *representation-level internal factors* only.

    Nothing here transports a physical CG or loop coefficient in an alpha
    dependent mass/rephasing convention.
    """
    rows=[]
    for m in models():
        cls=m['family']
        ops=validated['families'][cls]['operator_internal_conventions']
        rows.append({'model':m['label'], 'family':cls, 'alpha':m['alpha'],
                     'internal_tensor_factor_by_operator':{
                         r['coupling']:r['internal_tensor_prefactor'] for r in ops},
                     'field_ordering_has_alpha_dependent_Fc':m['alpha']!=-1,
                     'SU2_ward_conditional':True,
                     'physical_action_verified':False,
                     'physical_source_certified':False})
    require(len(rows)==16 and len({r['model'] for r in rows})==16,
            'wrong all-16 model enumeration')
    return rows


def analyse(root: Path, native: Path) -> dict:
    provenance=source_provenance(root)
    data=read_report(native)
    native_snapshot_guard(data,root,provenance)
    validated=validate_native_data(data)
    sixteen=expand_all16(validated)
    return {'format':'b2-internal-action-crosspackage-boundary-v1',
            'status':'NATIVE_RGBETA_INTERNAL_RATIOS_30_OF_30_ACTION_CONVERSION_PENDING',
            'source_provenance':provenance,
            'native_file':native.relative_to(root).as_posix() if native.is_relative_to(root) else str(native),
            'native_package_path_from_report':data.get('rgbeta_package_file'),
            'GroupInvBar_symbols_from_report':data.get('GroupInvBar_symbols'),
            'GroupInvBar_downvalues_from_report':data.get('GroupInvBar_downvalues_inputform'),
            'native':validated,'all16_internal_factor_transport':sixteen,
            'verified_native_internal_registration':True,
            'physical_source_certified':False,
            'physical_crosspackage_complete':False,
            'action_prefactor_equated_to_registry_factor':False,
            'native_numeric_RGBeta_component_axes_verified':False,
            'signed_y1_y2_MF_lambdaT3_map':None,
            'minimal_missing_native_evidence':[
                'Explicit RGBeta barred SU2 index map and delS2/eps component realization in Matchete weight axes',
                'Independent physical Lagrangian coefficient/Feynman-rule witness for y1, y2, MF, lambdaT3',
                'B1 four-spinor to left-Weyl relative conjugation phase eta',
            ]}


def write_outputs(output_dir:Path, report:dict) -> None:
    require(not output_dir.exists(),f'refusing existing output directory {output_dir}')
    output_dir.mkdir(parents=True)
    (output_dir/'internal_action_boundary.json').write_text(
        json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    (output_dir/'HANDOFF.md').write_text(
        '# B2 native RGBeta internal action-boundary handoff\n\n'
        f"Status: `{report['status']}`.\n\n"
        f"Measured local Git revision: `{report['source_provenance']['repository_head_measured']}`.\n\n"
        'Source guard: all six upstream-inspected tracked blobs verified via Git clean filters.\n\n'
        'Native checks: five families, six registry couplings each = 30 exact native '
        'internal tensor prefactors. Yukawas `Sqrt[2]`, fermion mass `2`, quartics `6`. '
        'These numbers are **not physical Lagrangian or Feynman-rule factors**.\n\n'
        'Tested orders: `y1:{Bar[S1],l,F}`, `y2:{S2,l,F/Fc}`, '
        '`MF:{F,F/Fc}`, `lambdaT3:{H,H,S1,Bar[S2]}`. '
        'Mixing coupling is nonselfconjugate; portals are real.\n\n'
        'All-16 internal representation-level factors are expanded in the JSON. '
        'Physical source certificates: **0/16**. Not safe to integrate as a physical coupling map.\n\n'
        'Still required: actual RGBeta SU2 index realization, action-level factors, '
        'and B1 spinor-phase interface. The B/C literal Ward mismatch remains conditional.\n',
        encoding='utf-8')
