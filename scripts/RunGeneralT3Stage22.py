"""Stage 22: fail-closed quartic tensor certification across 16 T3 models.

Static checks alone NEVER certify Matchete/RGBeta conventions. Exact physical
maps require two independent native tensor exports and explicit potential
prefactors. No numeric fitting, tolerance thresholds, or identity fallback.

Run: python -m scripts.RunGeneralT3Stage22 --static-only
     python -m scripts.RunGeneralT3Stage22
"""
from __future__ import annotations

import argparse
import ast
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from scripts.Stage22TensorAlgebra import (TensorCertificateError, exact, potential_map,
                                          same_species_symmetry, solve_basis, tensor_rank, normalized_record)

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'
MANIFEST = ROOT/'docs/general_scalar_trial_operator_manifest.json'
GENERATED = ROOT/'RGE/general/generated'
RGBETA_WL = ROOT/'RGE/running/rgbeta/T3RGBetaModel.wl'
MAT_WL = ROOT/'Lagrangian/interactions/GeneralScalarOperators.wl'
VERIFIED_WL = ROOT/'Lagrangian/interactions/GeneralScalarVerifiedQuartics.wl'
REPORT = ROOT/'output/general_t3_stage22/stage22_summary.json'
M_TENSORS = ROOT/'output/general_t3_stage22/matchete_tensors.json'
R_TENSORS = ROOT/'output/general_t3_stage22/rgbeta_native_tensors.json'
R_DECLARATIONS = ROOT/'output/general_t3_stage22/rgbeta_native_declarations.json'
R_INVARIANTS = ROOT/'output/general_t3_stage22/rgbeta_native_invariants.json'
EXPECTED = {'matchete_extra_quartics':50,'rgbeta_preexisting_structural':12,
            'rgbeta_new_vertices':38,'quadratic_couplings_blocked':8}
SPECIES = ('H','Hdag','S1','S1dag','S2','S2dag')
REP = {'A':(1,3,2),'B':(2,2,1),'C':(2,2,3),'D':(3,1,2),'E':(3,3,2)}


def ensure(condition: bool, message: str) -> None:
    if not condition: raise TensorCertificateError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    ensure(isinstance(obj, dict),f'{path}: expected JSON object')
    return obj


def expected_symbols(entry: dict) -> list[tuple[tuple[int,...],str]]:
    out = [(tuple(x['counts']),x['name']) for x in entry['new_rgbeta_symbols']]
    out += [(tuple(x['counts']),x['existing_name']) for x in entry['preexisting_rgbeta_symbols']]
    return out


def static_audit() -> tuple[dict, dict]:
    for f in (REGISTRY,MANIFEST,RGBETA_WL,MAT_WL,VERIFIED_WL):
        ensure(f.is_file(),f'Missing source: {f}')
    reg, man = load_json(REGISTRY), load_json(MANIFEST)
    ensure(reg.get('status')=='GeneratedCandidateBasisNotWolframVerified',
           'Registry changed: inspect before using Stage 22')
    ensure(reg.get('totals')==EXPECTED, 'Stage 20 slot totals changed')
    models = reg.get('models')
    ensure(isinstance(models,dict) and len(models)==16, 'Registry does not contain exactly 16 models')
    ensure(man.get('model_count')==16 and set(man.get('models',{}))==set(models),
           'Matchete structural manifest does not match RGBeta registry')
    native_src=RGBETA_WL.read_text(encoding='utf-8-sig')
    mat_src=MAT_WL.read_text(encoding='utf-8-sig')
    ensure('T3ExtendedBosonicTensorBasis' in mat_src and
           'PhysicalSU2InvariantBasis' in mat_src and 'T3ExtendedSymmetrizeGroup' in mat_src,
           'Matchete source construction unexpectedly changed')
    verified_table={}
    for k, payload in re.findall(r'"(\d+_\d+_-?\d+)"\s*->\s*(\{[^\n]*\})',
                                 VERIFIED_WL.read_text(encoding='utf-8-sig')):
        verified_table[k] = {tuple(row[0]):int(row[1])
                             for row in ast.literal_eval(payload.replace('{','[').replace('}',']'))}
    ensure(len(verified_table)==13, 'Verified Matchete multiplicity table unexpectedly changed')
    stats={'models':16,'new_vertices':0,'reused_vertices':0,'total_slots':0,
           'quadratics_still_blocked':0,'blocks':0}
    blocks = {}
    for model,entry in sorted(models.items()):
        cls, a = model.split('_')
        ensure(cls in REP and int(a)==entry['alpha'] and tuple(entry['dimensions'])==REP[cls],
               f'{model}: inconsistent model assignment')
        manifest = man['models'][model]
        expected_counts = {tuple(x['counts']):int(x['invariants'])
                           for x in manifest['extra_quartic_structures']}
        key = f'{entry["dimensions"][0]}_{entry["dimensions"][1]}_{entry["alpha"]}'
        ensure(key in verified_table and verified_table[key] == expected_counts,
               f'{model}: independent verified Wolfram multiplicity table differs')
        grouped=defaultdict(list)
        for counts,symbol in expected_symbols(entry): grouped[counts].append(symbol)
        ensure({c:len(v) for c,v in grouped.items()}==expected_counts,
               f'{model}: Matchete multiplicities disagree with symbol registry')
        ext=GENERATED/f'{model}_extra_quartics.wl'
        ensure(ext.is_file(),f'Missing generated RGBeta extension: {ext}')
        source=ext.read_text(encoding='utf-8-sig')
        ensure(re.search(r'T3GeneralQuarticExpectedCount\s*=\s*'+str(entry['total_additional_quartics'])+r'\s*;',source),
               f'{model}: incorrect registered total')
        for v in entry['new_rgbeta_symbols']:
            s=v['name']
            ensure(len(v['counts'])==6 and sum(v['counts'])==4,
                   f'{model}/{s}: bad quartic counts')
            ensure(len(v['fields'])==4 and isinstance(v['invariant'],str),
                   f'{model}/{s}: incomplete invariant')
            ensure(re.search(r'AddQuartic\[\s*'+re.escape(s)+r'\s*,',source) is not None,
                   f'{model}/{s}: absent AddQuartic source')
            ensure(re.search(r'AppendTo\[\s*T3GeneralQuarticAddedSymbols\s*,\s*'+re.escape(s)+r'\s*\]',source) is not None,
                   f'{model}/{s}: absent added-symbol registration')
            ensure(''.join(v['invariant'].split()) in ''.join(source.split()),
                   f'{model}/{s}: source invariant changed relative to registry')
        for v in entry['preexisting_rgbeta_symbols']:
            name=v['existing_name']
            ensure(re.search(r'AddQuartic\[\s*'+re.escape(name)+r'\s*,',native_src) is not None,
                   f'{model}/{name}: reused vertex absent from native RGBeta source')
        stats['new_vertices']+=len(entry['new_rgbeta_symbols'])
        stats['reused_vertices']+=len(entry['preexisting_rgbeta_symbols'])
        stats['total_slots']+=sum(map(len,grouped.values()))
        stats['quadratics_still_blocked']+=int(entry['quadratics_not_implemented'])
        stats['blocks']+=len(grouped)
        blocks[model] = {''.join(map(str,c)):names for c,names in sorted(grouped.items())}
    ensure((stats['new_vertices'],stats['reused_vertices'],stats['total_slots'],
            stats['quadratics_still_blocked'])==(38,12,50,8), 'Broken 16-model totals')
    files=[REGISTRY,MANIFEST,RGBETA_WL,MAT_WL,VERIFIED_WL]
    files += [GENERATED/f'{k}_extra_quartics.wl' for k in models]
    return reg,{'status':'SourceStructureOnly','totals':stats,
                'hashes':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files},
                'expected_rgbeta_symbols_by_block':blocks}


def tensor_index(payload: dict, model: str, expected_basis: str,
                 explicit_conventions: bool) -> dict[tuple[int,...],list[dict]]:
    ensure(payload.get('format')=='stage22-exact-tensors-v1',
           f'{expected_basis}: wrong tensor-export format')
    ensure(payload.get('basis')==expected_basis,
           f'Expected native {expected_basis} tensors')
    ensure(payload.get('source')==('MatcheteInvariantTensors' if expected_basis=='Matchete'
                                  else 'RGBetaNativeQuartic'),
           f'{expected_basis}: cannot use candidate/reconstructed tensor as native proof')
    if explicit_conventions:
        ensure(payload.get('native_invariant_definitions_verified') is True and
               payload.get('native_potential_prefactors_verified') is True and
               payload.get('field_component_alignment_verified') is True,
               'RGBeta native tensor/potential/convention evidence is missing')
    models=payload.get('models',{})
    ensure(model in models, f'{expected_basis}/{model}: native export missing')
    tensors=models[model].get('tensors')
    ensure(isinstance(tensors,list),f'{expected_basis}/{model}: missing tensors')
    out=defaultdict(list)
    for record in tensors:
        same_species_symmetry(record)
        counts=tuple(record['counts'])
        expected_species=[x for i,x in enumerate(SPECIES) for _ in range(record['counts'][i])]
        ensure(record.get('species')==expected_species,
               f'{expected_basis}/{model}/{record.get("symbol")}: component order unverified')
        ensure(record.get('potential_prefactor') is not None,
               f'{expected_basis}/{model}/{record.get("symbol")}: missing Lagrangian prefactor')
        ensure(exact(record['potential_prefactor'])!=0,
               f'{expected_basis}/{model}/{record.get("symbol")}: zero Lagrangian prefactor')
        out[counts].append(record)
    return out


def certify(reg: dict, mat: dict, rgb: dict) -> tuple[dict, bool]:
    result={}; failures=0; ncert=0
    for model,entry in sorted(reg['models'].items()):
        expected=defaultdict(list)
        for c,name in expected_symbols(entry): expected[c].append(name)
        model_rows=[]
        try:
            m=tensor_index(mat,model,'Matchete',False)
            r=tensor_index(rgb,model,'RGBeta',True)
        except (TensorCertificateError,KeyError,TypeError) as e:
            m=r={}
            global_error=str(e)
        else: global_error=None
        if set(m)!=set(expected) or set(r)!=set(expected):
            if global_error is None:
                global_error='Native exports have missing/extra quartic field-count blocks'
        for counts,names in sorted(expected.items()):
            key=''.join(map(str,counts))
            try:
                if global_error: raise TensorCertificateError(global_error)
                mm=m[counts]; rr=r[counts]
                expected_m=['lambdaT3ExtraQ'+key+'I'+str(i+1) for i in range(len(names))]
                ensure([x['symbol'] for x in mm]==expected_m,
                       f'{model}/{key}: Matchete I# basis order mismatch')
                ensure([x['symbol'] for x in rr]==names,
                       f'{model}/{key}: native RGBeta symbols/order mismatch')
                # The exporter must align conjugation orientation explicitly.
                orient=rgb['models'][model].get('holomorphic_orientation')
                ensure(orient in ('direct','conjugate'),
                       f'{model}/{key}: unknown complex conjugation orientation')
                transformation=solve_basis(mm,rr)
                physics=potential_map(transformation,
                                      [x['potential_prefactor'] for x in mm],
                                      [x['potential_prefactor'] for x in rr],orient)
                model_rows.append({'counts':list(counts),'status':'ExactAlgebraPassed',
                                   'tensor_map':transformation,'potential_map':physics,
                                   'certification_scope':'ConditionalOnIndependentNativeExportProvenance'})
                ncert+=len(names)
            except (TensorCertificateError,ValueError,KeyError,TypeError) as e:
                model_rows.append({'counts':list(counts),'status':'Uncertified',
                                   'symbols':names,'reason':str(e)})
                failures+=len(names)
        result[model]={'slots':len(expected_symbols(entry)),
                       'blocks':model_rows,'status':('ExactAlgebraPassed'
                       if all(x['status']=='ExactAlgebraPassed' for x in model_rows)
                       and not global_error else 'Uncertified')}
    return {'models':result,'exact_algebra_slots':ncert,'uncertified_slots':failures}, failures==0



def uncertified_inventory(reg:dict,reason:str) -> dict:
    models={}
    for model,entry in sorted(reg['models'].items()):
        groups=defaultdict(list)
        for counts,symbol in expected_symbols(entry):groups[counts].append(symbol)
        models[model]={'status':'Uncertified', 'slots':len(expected_symbols(entry)),
                       'blocks':[{'counts':list(c),'symbols':v, 'status':'Uncertified',
                                  'reason':reason} for c,v in sorted(groups.items())]}
    return {'models':models,'exact_algebra_slots':0,
            'uncertified_slots':sum(x['slots'] for x in models.values())}


def _native_inputform(value: Any, context: str, field: str) -> str:
    """Reject Wolfram's serialized Missing[] (a nonempty string!)."""
    ensure(isinstance(value, str) and value.strip(),
           f'{context}: absent native {field}')
    bad = ('Missing[', '$Failed', '$Aborted', 'Indeterminate',
           'TemplateBox[', 'Unevaluated[')
    ensure(not any(marker in value for marker in bad),
           f'{context}: unresolved native {field}: {value[:130]}')
    return value


def audit_native_declarations(reg:dict,path:Path) -> dict:
    payload=load_json(path)
    ensure(payload.get('format')=='stage22-native-rgbeta-declarations-v1'
           and payload.get('source')=='RGBetaRegisteredQuarticAssociation'
           and payload.get('status')=='NativeDeclarationMetadataOnly',
           'Unexpected RGBeta native-declaration format')
    ensure(payload.get('matchete_rgbeta_normalisation_verified') is False,
           'Declarations incorrectly claim a tensor certificate')
    actual=payload.get('models',{})
    ensure(set(actual)==set(reg['models']),
           'Native RGBeta declarations missing one or more models')
    seen=0
    for model,entry in reg['models'].items():
        expect=[name for _,name in expected_symbols(entry)]
        found=actual[model].get('declarations')
        ensure(isinstance(found,list) and [x.get('symbol') for x in found]==expect,
               f'{model}: native RGBeta declarations disagree with candidate registry')
        for x in found:
            context=f'{model}/{x.get("symbol")}'
            ensure(x.get('self_conjugate') is False,
                   f'{context}: unexpected self-conjugate native RGBeta quartic')
            for field in ('fields_inputform','invariant_inputform',
                          'coupling_factor_inputform','unique_arrangements_inputform'):
                _native_inputform(x.get(field),context,field)
            field_list=x['fields_inputform'].strip()
            ensure(field_list.startswith('{') and field_list.endswith('}'),
                   f'{context}: native Fields must be an explicit four-field list')
        seen+=len(found)
    ensure(seen==50,'RGBeta native declaration audit incomplete')
    return {'status':'NativeRGBetaDeclarationsVerifiedNoCommonBasisComponents',
            'coupling_slots_checked':seen,'sha256':sha(path),
            'component_normalisation_certified':False}


def audit_native_matchete(reg:dict,path:Path) -> dict:
    """Independently validate all 50 Matchete components without RGBeta inputs.

    In particular, accept SU(2) singlets only as singleton component axes;
    test Bose permutations and full rank of every individual field block.
    """
    payload=load_json(path)
    ensure(payload.get('format')=='stage22-exact-tensors-v1'
           and payload.get('source')=='MatcheteInvariantTensors'
           and payload.get('basis')=='Matchete',
           'Unrecognised Matchete tensor provenance')
    ensure(set(payload.get('models',{})) == set(reg['models']),
           'Matchete export must contain exactly the 16 registry models')
    slots=0
    blocks=0
    for model, entry in sorted(reg['models'].items()):
        native=tensor_index(payload, model,'Matchete',False)
        expected=defaultdict(list)
        for counts,name in expected_symbols(entry):
            expected[counts].append(name)
        ensure(set(native)==set(expected),
               f'Matchete/{model}: missing or extra field-content block')
        for counts, rgbeta_names in sorted(expected.items()):
            mats=native[counts]
            key=''.join(map(str,counts))
            wanted=['lambdaT3ExtraQ'+key+'I'+str(i+1)
                    for i in range(len(rgbeta_names))]
            ensure([r['symbol'] for r in mats]==wanted,
                   f'Matchete/{model}/{key}: unexpected independent basis order')
            species=mats[0]['species']
            dims=[2 if name.startswith('H') else
                  entry['dimensions'][0] if name.startswith('S1') else
                  entry['dimensions'][1] for name in species]
            ensure(all(row['dimensions']==dims for row in mats),
                   f'Matchete/{model}/{key}: incorrect representation dimensions')
            cols=[normalized_record(row, row['symbol'])['components'] for row in mats]
            ensure(tensor_rank(cols)==len(mats),
                   f'Matchete/{model}/{key}: rank deficient physical tensors')
            slots+=len(mats)
            blocks+=1
    ensure(slots==50 and blocks==46,
           f'Matchete native audit expects 50 slots/46 blocks, found {slots}/{blocks}')
    return {'status':'MatcheteNativeComponentsAndRanksVerified',
            'models_checked':len(reg['models']),
            'tensor_slots_checked':slots,'field_content_blocks_checked':blocks,
            'rgbeta_tensor_alignment_verified':False,
            'matchete_rgbeta_normalisation_certified':False,
            'sha256':sha(path)}


def audit_native_rgbeta_invariants(reg:dict,path:Path) -> dict:
    """Audit native stored RGBeta invariant expressions, NOT components."""
    payload=load_json(path)
    ensure(payload.get('format')=='stage22-native-rgbeta-invariants-v1'
           and payload.get('source')=='RGBetaRegisteredQuarticInvariantFunctions'
           and payload.get('status')=='NativeInvariantExpressionsNotCommonBasisTensors',
           'Unrecognised RGBeta native invariant-expression provenance')
    ensure(payload.get('native_invariant_component_alignment_verified') is False
           and payload.get('native_potential_prefactors_verified') is False,
           'Diagnostic RGBeta expressions may not assert physical certification')
    ensure(set(payload.get('models',{}))==set(reg['models']),
           'RGBeta native invariant expressions missing some models')
    count=0
    for model,entry in reg['models'].items():
        expected=[name for _,name in expected_symbols(entry)]
        found=payload['models'][model].get('invariants')
        ensure(isinstance(found,list)
               and [x.get('symbol') for x in found]==expected,
               f'RGBeta/{model}: native expression symbol order mismatch')
        for row in found:
            ensure(row.get('self_conjugate') is False
                   and row.get('component_tensor_exported') is False
                   and row.get('verified_potential_prefactor') is None,
                   f'RGBeta/{model}/{row.get("symbol")}: diagnostic metadata not fail-closed')
            for field in ('native_fields_inputform','native_invariant_function_inputform',
                          'native_invariant_expression_inputform',
                          'native_coupling_function_inputform',
                          'native_unique_arrangements_inputform'):
                _native_inputform(row.get(field),
                                  f'RGBeta/{model}/{row.get("symbol")}',field)
            field_list=row['native_fields_inputform'].strip()
            ensure(field_list.startswith('{') and field_list.endswith('}'),
                   f'RGBeta/{model}/{row.get("symbol")}: native Fields must be a list')
        count+=len(found)
    ensure(count==50, f'RGBeta native expression audit expects 50 slots, got {count}')
    return {'status':'RGBetaNativeInvariantExpressionsVerifiedNotComponents',
            'models_checked':16,'quartic_slots_checked':count,
            'component_tensors_exported':False,
            'matchete_rgbeta_normalisation_certified':False,
            'sha256':sha(path)}


def audit_native_rgbeta_coherence(declarations_path: Path, invariants_path: Path) -> dict:
    """Require the two independent diagnostic exports to agree in every slot.

    This checks native register/field/invariant-function/coupling metadata only.
    It does not equate either to Matchete tensors or establish RGE conventions.
    """
    declarations=load_json(declarations_path).get('models',{})
    invariants=load_json(invariants_path).get('models',{})
    ensure(set(declarations)==set(invariants),
           'RGBeta native declaration and invariant exports cover different models')
    pairs=(('fields_inputform','native_fields_inputform'),
           ('invariant_inputform','native_invariant_function_inputform'),
           ('coupling_factor_inputform','native_coupling_function_inputform'),
           ('unique_arrangements_inputform','native_unique_arrangements_inputform'))
    n=0
    for model in sorted(declarations):
        left=declarations[model].get('declarations',[])
        right=invariants[model].get('invariants',[])
        ensure([x.get('symbol') for x in left] == [x.get('symbol') for x in right],
               f'{model}: native RGBeta evidence order/symbols disagree')
        for a,b in zip(left,right):
            symbol=a['symbol']
            for f1,f2 in pairs:
                _native_inputform(a.get(f1),f'{model}/{symbol}',f1)
                _native_inputform(b.get(f2),f'{model}/{symbol}',f2)
                ensure(a[f1]==b[f2],
                       f'{model}/{symbol}: independent native exports disagree on {f1}')
            ensure(a.get('self_conjugate') is False and b.get('self_conjugate') is False,
                   f'{model}/{symbol}: unexpected native quartic conjugation')
            n+=1
    ensure(n==50 and len(declarations)==16,
           f'Native RGBeta evidence coherence incomplete: {len(declarations)} models, {n} slots')
    return {'status':'TwoNativeRGBetaDiagnosticChannelsAgree',
            'coupling_slots_checked':n,'native_tensor_components_certified':False,
            'matchete_rgbeta_coupling_mapping_certified':False}


def execute(out:Path=REPORT, mat_path:Path=M_TENSORS, rgb_path:Path=R_TENSORS,
            static_only:bool=False, declarations_path:Path=R_DECLARATIONS,
            invariants_path:Path=R_INVARIANTS, inspect_native:bool=False) -> dict:
    output={'stage':22,'certification_status':'Uncertified','production_modified':False,
            'normalisation_certified':False,'errors':[]}
    try:
        reg, audit=static_audit()
        output['static_audit']=audit
        # A static audit MUST NOT read stale or partly written native files.
        # This is essential when the Matchete exporter is being repaired.
        if not static_only or inspect_native:
            if declarations_path.is_file():
                output['native_rgbeta_declarations']=audit_native_declarations(reg,declarations_path)
            if mat_path.is_file():
                output['native_matchete']=audit_native_matchete(reg,mat_path)
            if invariants_path.is_file():
                output['native_rgbeta_invariant_expressions']=audit_native_rgbeta_invariants(reg,invariants_path)
            if declarations_path.is_file() and invariants_path.is_file():
                output['native_rgbeta_coherence']=audit_native_rgbeta_coherence(
                    declarations_path,invariants_path)
        if static_only:
            output['certification']=uncertified_inventory(reg,'Static data only; exact native tensors not supplied')
            output['errors'].append('Static audit is NOT a coupling normalisation proof')
        else:
            missing=[str(p) for p in (mat_path,rgb_path) if not p.is_file()]
            if missing:
                reason='Missing independent native tensor export(s): '+', '.join(missing)
                output['certification']=uncertified_inventory(reg,reason)
                output['errors'].append(reason)
            else:
                m=load_json(mat_path);r=load_json(rgb_path)
                output['native_export_hashes']={str(mat_path):sha(mat_path),str(rgb_path):sha(rgb_path)}
                report,passed=certify(reg,m,r)
                output['certification']=report
                if passed:
                    output['certification_status']='ExactTensorAndExplicitPrefactorAlgebraPassed'
                    output['normalisation_certified']=True
                else:
                    output['errors'].append(f"{report['uncertified_slots']} of 50 coupling slots remain uncertified")
    except (TensorCertificateError,ValueError,KeyError,FileNotFoundError) as e:
        output['errors'].append(str(e))
        if 'static_audit' in output and 'certification' not in output:
            output['certification']=uncertified_inventory(
                reg, 'Invalid or incomplete native evidence: '+str(e))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return output


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--static-only',action='store_true')
    parser.add_argument('--matchete-tensors',type=Path,default=M_TENSORS)
    parser.add_argument('--rgbeta-tensors',type=Path,default=R_TENSORS)
    parser.add_argument('--rgbeta-declarations',type=Path,default=R_DECLARATIONS)
    parser.add_argument('--rgbeta-invariants',type=Path,default=R_INVARIANTS)
    parser.add_argument('--inspect-native',action='store_true',
                        help='Check both diagnostic native exports, without claiming a physical map')
    parser.add_argument('--output',type=Path,default=REPORT)
    args=parser.parse_args()
    output=execute(args.output,args.matchete_tensors,args.rgbeta_tensors,
                   args.static_only or args.inspect_native, args.rgbeta_declarations,
                   args.rgbeta_invariants, args.inspect_native)
    print('Stage 22:',output['certification_status'])
    if 'static_audit' in output:
        print('Static slots:',output['static_audit']['totals'])
    for e in output['errors']: print('UNCERTIFIED:',e)
    print('Report:',args.output)
    if args.inspect_native:
        inspected=('native_matchete' in output and
                   'native_rgbeta_invariant_expressions' in output and
                   'native_rgbeta_coherence' in output and
                   'static_audit' in output and
                   output['errors']==['Static audit is NOT a coupling normalisation proof'])
        return 0 if inspected else 1
    if args.static_only:
        return 0 if ('static_audit' in output and
                     output['errors']==['Static audit is NOT a coupling normalisation proof']) else 1
    return 0 if output['normalisation_certified'] else 1


if __name__=='__main__':
    raise SystemExit(main())
