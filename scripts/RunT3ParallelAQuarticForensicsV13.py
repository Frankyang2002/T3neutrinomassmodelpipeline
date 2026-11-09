"""Read-only forensic analysis of the *existing* failed native v10 smoke JSON.

No numerical RGBeta certification, no mutations to the saved witness.
Writes a separate source-hash-guarded report.  Exit 0 means the failure was
classified, NOT that native projectors passed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import (
    FORMAT, SOURCE_FILES, expected_symbols, file_hash, DIMENSIONS,
)

EXPECTED = ('A_-4', 'B_-1', 'E_-2')
BLOCKED = 'Blocked_NativeRGBetaQuarticProjectorResponses'
STATUS = 'NativeQuarticSmokeFailureClassified_NotCertified'


class ForensicBlocked(ValueError):
    pass


def require(cond: bool, why: str) -> None:
    if not cond:
        raise ForensicBlocked(why)


def analyze(smoke: dict[str, Any], root: Path, witness_sha: str) -> dict:
    """Strict provenance and row-oriented failure classification."""
    require(smoke.get('format') == FORMAT, 'not the original native v10 format')
    require(smoke.get('status') == BLOCKED, 'expected a blocked native smoke witness')
    require(smoke.get('smoke_only') is True and
            smoke.get('registered_model_count') == 3, 'not a three-model smoke witness')
    require(smoke.get('physical_Wilson_evolution_authorized') is False,
            'native witness falsely authorizes physics')
    require(smoke.get('native_numeric_S2_frame_measured') is False and
            smoke.get('full_scalar_potential_dictionary_certified') is False,
            'native witness falsely certifies component/physical normalization')
    hashes = smoke.get('source_sha256')
    require(isinstance(hashes, dict) and set(hashes) == set(SOURCE_FILES),
            'unexpected native source inventory')
    for name in SOURCE_FILES:
        path = root / name
        require(path.is_file() and hashes[name] == file_hash(path),
                'saved smoke source hash differs from local file: '+name)
    ep = smoke.get('installed_RGBeta_main_file')
    require(isinstance(ep, str) and Path(ep).is_file() and
            smoke.get('installed_RGBeta_main_file_sha256') == file_hash(Path(ep)),
            'installed RGBeta entrypoint differs from saved smoke')
    models = smoke.get('models')
    require(isinstance(models, dict) and set(models) == set(EXPECTED),
            'smoke model list changed')
    failures: list[dict] = []
    counts: Counter = Counter()
    totals = Counter()
    selection: list[dict] = []
    model_rows: dict[str, dict] = {}
    for tag in EXPECTED:
        row = models[tag]
        require(isinstance(row, dict), f'{tag}: model witness missing')
        expected = expected_symbols(tag)
        require(row.get('registered_quartics_count') == len(expected),
                f'{tag}: source quartic inventory count differs')
        inventory = row.get('inventory')
        require(isinstance(inventory, list) and
                {x.get('symbol') for x in inventory} == set(expected) and
                len(inventory) == len(expected), f'{tag}: native inventory differs')
        require(row.get('registry_exact_source_inventory') is True,
                f'{tag}: native inventory did not pass')
        real = {k for k,v in expected.items() if v}
        diag = row.get('native_real_self_projectors')
        require(isinstance(diag, list) and
                {x.get('symbol') for x in diag} == real and
                len(diag) == len(real), f'{tag}: native real diagonal coverage changed')
        model_counts = Counter()
        def classify(target: str, source: str, raw: Any, area: str) -> None:
            totals[area] += 1
            if not isinstance(raw, dict):
                category='malformed_response'; reason='not_object'; val=None
            elif raw.get('ok') is not True:
                category='native_evaluation_failed'
                reason=str(raw.get('reason','missing_failure_reason'))
                val=None
            else:
                val=raw.get('value')
                if not isinstance(val,(int,float)) or isinstance(val,bool):
                    category='malformed_response';reason='nonnumeric_reported_value'
                elif not (-1e100<float(val)<1e100):
                    category='malformed_response';reason='nonfinite_reported_value'
                elif target==source and abs(float(val)-1)>2e-8:
                    category='non_unit_self_response';reason='diagonal_differs_from_one'
                else:
                    category='observed';reason='numeric_native_observation'
            counts[category]+=1; model_counts[category]+=1
            if category!='observed':
                failures.append({'model':tag,'target':target,'source':source,
                    'area':area,'class':category,'reason':reason,
                    'value':val,'native_excerpt':{
                        k:str(raw[k])[:500] for k in ('expanded_inputform','coefficient_inputform')
                        if isinstance(raw,dict) and k in raw},})
                if not any(x['model']==tag for x in selection):
                    selection.append({'model':tag,'target':target,'source':source,
                                      'reason':reason,'area':area})
        for dr in diag:
            s=dr.get('symbol')
            classify(s,s,dr.get('native_response'),'diagonal')
        portal = row.get('native_portal_response_blocks')
        require(isinstance(portal,list) and len(portal)==2, f'{tag}: portal block count')
        for which, sector in enumerate(portal,1):
            d=DIMENSIONS[tag[0]][which-1]
            ss=[f'lambdaH{which}']+([f'lambdaH{which}Adj'] if d>1 else [])
            require(sector.get('which')==which and sector.get('symbols')==ss,
                    f'{tag}/S{which}: native portal symbols changed')
            matrix=sector.get('native_response_matrix')
            require(isinstance(matrix,list) and len(matrix)==len(ss) and
                    all(isinstance(z,list) and len(z)==len(ss) for z in matrix),
                    f'{tag}/S{which}: malformed native portal matrix')
            for i,target in enumerate(ss):
                for j,source in enumerate(ss):
                    classify(target,source,matrix[i][j],'portal')
            if all(isinstance(matrix[i][j],dict) and
                   matrix[i][j].get('ok') is True and
                   isinstance(matrix[i][j].get('value'),(int,float)) and
                   np.isfinite(matrix[i][j]['value'])
                   for i in range(len(ss)) for j in range(len(ss))):
                mat=np.array([[matrix[i][j]['value'] for j in range(len(ss))]
                             for i in range(len(ss))],dtype=float)
                sv=np.linalg.svd(mat,compute_uv=False)
                if sv[-1] <= 1e-8 or sv[0]/sv[-1] >= 1e8:
                    failures.append({'model':tag,'area':'portal','which':which,
                        'class':'singular_response_block',
                        'reason':'native_projector_response_matrix_is_rank_deficient',
                        'singular_values':sv.tolist()})
                    counts['singular_response_block']+=1
                    if not any(x['model']==tag for x in selection):
                        selection.append({'model':tag,'target':ss[0],
                          'source':ss[-1],
                          'reason':'singular_portal_response_matrix',
                          'area':'portal'})
        if (model_counts.get('native_evaluation_failed',0)==0 and
                model_counts.get('non_unit_self_response',0)==0 and
                row.get('ok') is not True):
            failures.append({'model':tag,'area':'model','class':'acceptance_flag_blocked',
                             'reason':'model rejected but individual rows appear numeric'})
        model_rows[tag]={
            'registry_exact_source_inventory':True,
            'native_projector_unit_flag':row.get('native_real_self_projectors_all_unit'),
            'native_portal_block_flag':row.get('native_portal_response_blocks_evaluated'),
            'response_counts':dict(sorted(model_counts.items())),
        }
    require(bool(failures), 'blocked witness contains no identifiable failure')
    return {
        'format':'t3-parallel-a-native-rgbeta-quartic-forensics-v13',
        'status':STATUS,
        'input_smoke_sha256':witness_sha,
        'input_smoke_source_sha256':hashes,
        'installed_RGBeta_main_file_sha256':smoke['installed_RGBeta_main_file_sha256'],
        'attempted_response_count':int(sum(totals.values())),
        'attempted_by_type':dict(sorted(totals.items())),
        'response_class_counts':dict(sorted(counts.items())),
        'models':model_rows,
        'failures':failures,
        'minimal_reprobe_targets':selection,
        'native_Lam_projector_certified':False,
        'numerical_S2_frame_certified':False,
        'physical_Wilson_evolution_authorized':False,
        'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT',
    }


def main(argv: list[str]|None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('.'))
    p.add_argument('--smoke',type=Path,default=Path('output/parallel_t3/wilson/rgbeta_quartic_projectors_native_v10_smoke.json'))
    p.add_argument('--output',type=Path,default=None)
    a=p.parse_args(argv)
    try:
        root=a.project_root.resolve()
        smoke_path=a.smoke if a.smoke.is_absolute() else root/a.smoke
        raw=smoke_path.read_bytes(); digest=sha256(raw).hexdigest()
        record=json.loads(raw.decode('utf-8-sig'))
        report=analyze(record,root,digest)
        output=(a.output if a.output is not None else
                root/'output/parallel_t3/wilson'/f'rgbeta_quartic_forensics_v13_{digest[:12]}.json')
        if not output.is_absolute():output=root/output
        if output.exists():
            old=json.loads(output.read_text(encoding='utf-8-sig'))
            require(old==report,'existing forensic output differs; refusing overwrite')
            print('v13: identical forensic output already exists; verified, not overwritten')
        else:
            output.parent.mkdir(parents=True,exist_ok=True)
            tmp=output.with_suffix(output.suffix+'.partial')
            require(not tmp.exists(),'refusing to overwrite partial diagnostic')
            tmp.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
            tmp.rename(output)
        print(f'v13: classified {len(report["failures"])} failure records across '
              f'{report["attempted_response_count"]} response evaluations')
        print('Breakdown: '+json.dumps(report['response_class_counts'],sort_keys=True))
        for x in report['failures'][:18]:
            print(f'  {x["model"]} {x.get("area")} '
                  f'{x.get("target","?")} <- {x.get("source","?")}: '
                  f'{x["class"]} / {x["reason"]}')
            if x.get('native_excerpt'):
                print('    '+json.dumps(x['native_excerpt'],sort_keys=True)[:350])
        print('Forensic output: '+str(output))
        print('Native v10 still BLOCKED. Classification success is not projector certification.')
        return 0
    except (ValueError,KeyError,OSError,TypeError,IndexError,ForensicBlocked,
            UnicodeError,ArithmeticError) as exc:
        print('v13 forensic BLOCKED: '+str(exc),file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
