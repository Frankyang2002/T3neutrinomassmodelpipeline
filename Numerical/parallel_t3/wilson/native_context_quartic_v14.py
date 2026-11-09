"""Read-only source-guarded v14 RGBeta PackageScope Lam / quartic projector audit.

Unlike the v10/v13 Global` calls, the v14 exporter invokes the installed
RGBeta`PackageScope`Lam and RGBeta`PackageScope`RefineGroupStructures.
No native component frame or physical Wilson normalization is certified.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any
import numpy as np

from .native_quartic_projectors_v10 import (
    DIMENSIONS, EXPECTED_MODELS, expected_symbols, compare_with_v9,
    V10Blocked,
)
from .portal_real_action_v9 import audit_from_paths as recheck_v9

FORMAT = 't3-parallel-a-native-rgbeta-quartic-projectors-v14'
SMOKE = 'NativeRGBetaQuarticProjectorSmoke3ModelsPassed'
FULL = 'NativeRGBetaRealQuarticLamProjectorResponsesObserved_Conditional'
BLOCKED = 'Blocked_NativeRGBetaQuarticProjectorResponses'
STATUS = 'NativePackageScopeQuartic150Diagonal110PortalResponsesVerified_Complex28ExploratoryConditional'
SOURCE_FILES = (
    'RGE/running/rgbeta/T3RGBetaModel.wl',
    'scripts/RGBetaNativeAccess.wl',
    'scripts/RunT3ParallelANativeQuarticProjectorsV14.wl',
)
API = {
    'Lam': ['RGBeta`PackageScope`Lam'],
    'RefineGroupStructures': ['RGBeta`PackageScope`RefineGroupStructures'],
    'TStructure': ['RGBeta`PackageScope`TStructure'],
}
SMOKE_MODELS = {'A_-4', 'B_-1', 'E_-2'}


class V14Blocked(ValueError):
    pass


def require(cond: bool, why: str) -> None:
    if not cond:
        raise V14Blocked('v14 native RGBeta: '+why)


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8-sig'), parse_constant=lambda v: (_ for _ in ()).throw(ValueError('Nonfinite JSON: '+v)))
    require(isinstance(obj, dict), 'not a JSON object: '+str(path))
    return obj


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def response(row: dict, what: str) -> float:
    require(isinstance(row, dict) and row.get('ok') is True,
            what+': native contraction is blocked: '+str(row)[:180])
    val = row.get('value')
    require(type(val) in (int, float) and math.isfinite(val), what+': nonfinite/non-numeric response')
    return float(val)


def check_provenance(raw: dict, root: Path) -> None:
    require(raw.get('format') == FORMAT, 'unexpected native format')
    require(raw.get('native_package_scoped_api_verified') is True,
            'PackageScope APIs were not verified by installed kernel')
    require(raw.get('resolved_native_package_symbols') == API,
            'native RGBeta PackageScope symbol names changed')
    hashes = raw.get('source_sha256')
    require(isinstance(hashes, dict) and set(hashes) == set(SOURCE_FILES),
            'missing or unexpected v14 source hashes')
    for p in SOURCE_FILES:
        path = root / p
        require(path.is_file() and hashes[p] == file_hash(path),
                'v14 native source changed: '+p)
    installed = raw.get('installed_RGBeta_main_file')
    require(isinstance(installed, str) and Path(installed).is_file() and
            file_hash(Path(installed)) == raw.get('installed_RGBeta_main_file_sha256'),
            'installed RGBeta entrypoint has changed')
    for key in ('physical_Wilson_evolution_authorized', 'full_scalar_potential_dictionary_certified',
                'native_numeric_S2_frame_measured', 'native_complex_projector_28_channels_certified',
                'native_Matchete_scalar_CG_evaluated'):
        require(raw.get(key) is False, 'unsupported physics certificate: '+key)
    require(raw.get('EFT_order') ==
            'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
            'EFT threshold order changed')


def audit_witness(raw: dict, root: Path, *, allow_blocked: bool = False) -> dict:
    """Independently validate full/smoke counts, native responses and complex observations.

    Blocked evidence is inspectable only when allow_blocked=True; it is NEVER
    returned as certified, nor accepted as a successful full audit.
    """
    root = root.resolve()
    check_provenance(raw, root)
    smoke = raw.get('smoke_only')
    require(type(smoke) is bool, 'missing smoke-only flag')
    models = raw.get('models')
    expected_models = SMOKE_MODELS if smoke else EXPECTED_MODELS
    require(isinstance(models, dict) and set(models) == expected_models and
            raw.get('registered_model_count') == len(expected_models),
            'native model census incorrect')
    status = raw.get('status')
    accepted_status = SMOKE if smoke else FULL
    require(status in ({accepted_status, BLOCKED} if allow_blocked else {accepted_status}),
            'native status is blocked or unexpected: '+str(status))
    if status == accepted_status:
        require(raw.get('native_Lam_evaluated') is True and
                raw.get('native_quartic_projectors_evaluated') is True,
                'claimed pass without complete Lam/projector evaluation')
    else:
        require(raw.get('native_Lam_evaluated') is False and
                raw.get('native_quartic_projectors_evaluated') is False,
                'blocked native witness incorrectly claims all projectors evaluated')
    counts = Counter()
    rows = {}
    failures = []
    complex_counts = Counter()
    for tag in sorted(expected_models):
        model = models[tag]
        require(isinstance(model, dict), tag+': missing model')
        sym = expected_symbols(tag)
        d1, d2, df = DIMENSIONS[tag[0]]
        require(model.get('dimensions') == [d1, d2, df] and
                model.get('alpha') == int(tag.split('_')[1]) and
                model.get('registered_quartics_count') == len(sym) and
                model.get('registry_exact_source_inventory') is True,
                tag+': source inventory mismatch')
        inventory = model.get('inventory')
        require(isinstance(inventory, list) and len(inventory) == len(sym) and
                {x.get('symbol') for x in inventory} == set(sym),
                tag+': malformed quartic inventory')
        for record in inventory:
            name = record['symbol']
            require(record.get('native_self_conjugate') is sym[name] and
                    record.get('has_installed_RGBeta_projector') is True,
                    tag+'/'+name+': projector registration/self conjugacy wrong')
            perms = record.get('native_unique_arrangement_count')
            require(type(perms) is int and 1 <= perms <= 24,
                    tag+'/'+name+': invalid permutation count')
            raw_coup = ''.join(record.get('native_coupling_inputform', '').split())
            require(raw_coup == f'6*Matrix[{name}][##1]&',
                    tag+'/'+name+': RGBeta complex-field factor != 6*Matrix')
        diagonal = model.get('native_real_self_projectors')
        require(isinstance(diagonal, list) and len(diagonal) == sum(sym.values()) and
                {x.get('symbol') for x in diagonal} == {x for x, real in sym.items() if real},
                tag+': missing real diagonal contractions')
        # Diagonal and portal checks may be blocked: record each native failure exactly.
        for row in diagonal:
            counts['real_diagonal_attempted'] += 1
            r = row.get('native_response')
            if (isinstance(r, dict) and r.get('ok') is True and
                    type(r.get('value')) in (float, int) and math.isfinite(r['value']) and
                    abs(r['value']-1) < 2.e-8):
                counts['real_diagonal_unit'] += 1
            else:
                failures.append({'model': tag, 'type': 'real_diagonal', 'symbol': row['symbol'],
                                 'native_response': r})
        blocks = model.get('native_portal_response_blocks')
        require(isinstance(blocks, list) and len(blocks) == 2,
                tag+': missing portal response blocks')
        sectors = []
        for which, dim in ((1, d1), (2, d2)):
            block = blocks[which-1]
            names = [f'lambdaH{which}'] + ([f'lambdaH{which}Adj'] if dim > 1 else [])
            vals = block.get('native_response_matrix')
            require(block.get('which') == which and block.get('dimension') == dim and
                    block.get('symbols') == names and
                    isinstance(vals, list) and len(vals) == len(names) and
                    all(isinstance(a, list) and len(a) == len(names) for a in vals),
                    tag+f'/S{which}: bad native portal response block')
            numeric = True
            R = []
            for i, line in enumerate(vals):
                obs = []
                for j, x in enumerate(line):
                    counts['portal_response_attempted'] += 1
                    if (isinstance(x, dict) and x.get('ok') is True and
                            type(x.get('value')) in (int, float) and math.isfinite(x['value'])):
                        obs.append(float(x['value']))
                        counts['portal_response_numeric'] += 1
                    else:
                        obs.append(0.)
                        numeric = False
                        failures.append({'model':tag,'type':'portal_response',
                                         'target':names[i],'source':names[j],
                                         'native_response':x})
                R.append(obs)
            if numeric:
                mat = np.array(R, dtype=float)
                eig = np.linalg.svd(mat, compute_uv=False)
                if np.max(abs(np.diag(mat)-1)) >= 2e-8 or eig[-1] < 1.e-8 or eig[0]/eig[-1] >= 1.e8:
                    failures.append({'model':tag, 'type':'portal_response_invalid',
                                     'which':which, 'matrix':R})
            sectors.append({'which': which,'dimension':dim,'symbols':names,
                            'native_projector_response_matrix':R if numeric else None})
            counts['portal_blocks_attempted'] += 1
        complex_rows = model.get('native_complex_conjugate_projector_probes')
        require(isinstance(complex_rows, list) and len(complex_rows) == sum(not v for v in sym.values()) and
                {x.get('symbol') for x in complex_rows} == {x for x, real in sym.items() if not real},
                tag+': missing native complex source/conjugate probe inventory')
        for cr in complex_rows:
            x = cr.get('native_complex_conjugate_response')
            complex_counts['complex_attempted'] += 1
            good = (isinstance(x, dict) and x.get('ok') is True and
                    x.get('unit_direct_zero_conjugate') is True and
                    type(x.get('direct_channel')) in (int,float) and
                    type(x.get('conjugate_channel')) in (int,float) and
                    all(math.isfinite(v) for v in (x['direct_channel'],x['conjugate_channel'])) and
                    abs(x['direct_channel']-1)<2e-8 and abs(x['conjugate_channel'])<2e-8)
            if good:
                complex_counts['complex_direct_unit_conj_zero'] += 1
            else:
                failures.append({'model':tag,'type':'complex_exploratory',
                                 'symbol':cr['symbol'],'native_response':x})
        counts['quartic_registered'] += len(sym)
        counts['real_registered'] += sum(sym.values())
        counts['complex_registered'] += sum(not v for v in sym.values())
        rows[tag] = {'portal_sectors':sectors,
                     'real_diagonals_attempted':len(diagonal),
                     'complex_channels_attempted':len(complex_rows)}
    require((counts['quartic_registered'], counts['real_registered'],
             counts['complex_registered'], counts['real_diagonal_attempted'],
             counts['portal_response_attempted'], complex_counts['complex_attempted']) ==
            ((38,29,9,29,21,9) if smoke else (178,150,28,150,110,28)),
            'unexpected native result census: '+str(dict(counts)))
    real_success = (counts['real_diagonal_unit'] == counts['real_diagonal_attempted'] and
                    counts['portal_response_numeric'] == counts['portal_response_attempted'] and
                    not any(f['type']=='portal_response_invalid' for f in failures))
    require(real_success == (status == accepted_status),
            'claimed native success differs from individual contraction evidence')
    require(all(model.get('ok') is real_success for model in models.values()) if real_success
            else any(model.get('ok') is False for model in models.values()),
            'model-level native status inconsistent')
    report = {'format':'t3-parallel-a-rgbeta-native-context-quartic-audit-v14',
              'status':('NativeRGBetaRealProjectorsSmokeVerified_ComplexExploratory'
                        if smoke and real_success else
                        STATUS if real_success else 'Blocked_NativeRGBetaRealProjectors_PackageScopeProbe'),
              'smoke_only':smoke,'source_native_sha256':raw['source_sha256'],
              'source_census':dict(counts),
              'complex_probe_census':dict(complex_counts),
              'complex_28_certified':False,
              'real_projector_conditionally_verified':real_success,
              'native_complex_exploratory_complete':complex_counts['complex_direct_unit_conj_zero']==complex_counts['complex_attempted'],
              'native_numeric_S2_frame_verified':False,
              'full_scalar_potential_dictionary_certified':False,
              'physical_Wilson_evolution_authorized':False,
              'models':rows,'failures':failures}
    if not real_success and not allow_blocked:
        raise V14Blocked('real native response blocked: '+str(failures[:3])[:500])
    return report


def audit_from_paths(root: Path, native_file: Path, v9_file: Path) -> dict:
    """Fresh full audit, including source-grounded v9 portal tensor revalidation."""
    root = root.resolve()
    native = read(native_file)
    report = audit_witness(native,root)
    require(not report['smoke_only'],'full native all-16 witness required')
    base=root/'output/parallel_t3/wilson'
    v9inputs=[base/x for x in (
        'matchete_portals_native_v5.json','rgbeta_portals_native_v5.json',
        'portal_dictionary_conditional_v5.json','rgbeta_fundamental_fierz_native_v6.json',
        'portal_orientation_conditional_v6.json','s2_intertwiner_conditional_v7.json',
        'rgbeta_s2_reduction_native_v8.json','s2_spectrum_conditional_v8.json')]
    v9_fresh=recheck_v9(*v9inputs,root,samples=3)
    saved=read(v9_file)
    require(saved.get('format') == v9_fresh['format'] and
            saved.get('status') == v9_fresh['status'] and
            saved.get('input_sha256') == v9_fresh['input_sha256'],
            'saved v9 portal-action source provenance stale')
    for tag in EXPECTED_MODELS:
        old = saved['models'][tag]
        new = v9_fresh['models'][tag]
        require(len(old)==len(new)==2,tag+': v9 portal sectors changed')
        for a,b in zip(old,new):
            require(a['which']==b['which'] and
                    np.allclose(a['conditional_rgbeta_from_matchete_real_matrix'],
                                b['conditional_rgbeta_from_matchete_real_matrix'],atol=2e-10),
                    tag+': v9 Matchete coupling matrix changed')
    # Transfer observation basis only. Does NOT identify physical CG components.
    bridge=compare_with_v9({'models':report['models']},v9_fresh)
    report['v9_native_portal_bridge']=bridge
    report['input_sha256']={'native_v14':file_hash(native_file),
                            'v9':file_hash(v9_file),
                            **{p.name:file_hash(p) for p in v9inputs}}
    return report
