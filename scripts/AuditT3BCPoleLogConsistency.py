"""Audit formal Matchete scalar-threshold pole/log identities and finite quartics.

Consumes completed ExtractT3BCFiniteScalarTerms v1 outputs, not live Matchete.
Does not certify MS-bar counterterms, full EFT matching, or neutrino masses.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from Numerical.experimental.T3BCPoleLogConsistency import (
    ConsistencyError, compare_residue_and_log, finite_groups, LAMBDA, terms,
)

ROOT=Path(__file__).resolve().parents[1]
INPUT=Path('output/audits/t3_bc_scalar_finite_term_extraction_v1')
OUTPUT=Path('output/audits/t3_bc_pole_log_consistency_v2')
MODELS=('B_-1','C_-1')
SECTORS=(('Y','LLHH'),('Y','LLHX'),('Y','LLXX'),('X','LLHH'))


def load(path:Path)->str:
    if not path.is_file():
        raise ConsistencyError('Missing audit evidence: '+str(path))
    return path.read_text(encoding='utf8').strip()


def _sha(data:str)->str:
    return sha256(data.encode('utf8')).hexdigest()


def run(input_dir:Path,output_dir:Path, models:tuple[str,...],root:Path|None=None)->dict:
    source=json.loads(load(input_dir/'summary.json'))
    if (source.get('format')!='t3-bc-native-local-matchete-pole-log-split-v1' or
        source.get('status')!='SymbolicLocalAmplitudePartitionedNoRenormalizedPhysicalMatching'):
        raise ConsistencyError('Unexpected finite-term audit format/status')
    rows=source['sector_results']
    by_index={(r['model'],r['threshold'],r['sector']):r for r in rows}
    if len(by_index)!=len(rows):
        raise ConsistencyError('Duplicate sector in source summary')
    planned=[]
    for model in models:
        if model not in MODELS:
            raise ConsistencyError('Unrecognized T3 benchmark model')
        for stage,sector in SECTORS:
            key=(model,stage,sector)
            if key not in by_index:
                raise ConsistencyError('Missing source summary sector: '+str(key))
            ref=by_index[key]
            for tag in ('pole_residue','pole_free_unrenormalized',
                        'at_mu_equal_heavy_mass_unrenormalized'):
                rel=Path(ref[tag+'_file'])
                if rel.is_absolute() or '..' in rel.parts:
                    raise ConsistencyError('Unsafe source relative path')
            residue=load(input_dir/ref['pole_residue_file'])
            pole_free=load(input_dir/ref['pole_free_unrenormalized_file'])
            constant=load(input_dir/ref['at_mu_equal_heavy_mass_unrenormalized_file'])
            for data,k in ((residue,'pole_residue_sha256'),
                           (pole_free,'pole_free_sha256'),
                           (constant,'scale_equal_mass_sha256')):
                if _sha(data)!=ref[k]:
                    raise ConsistencyError('InputForm SHA256 mismatch: '+str(key)+' '+k)
            if root is not None:
                source_path=root/ref['source_file']
                if _sha(load(source_path))!=ref['source_sha256']:
                    raise ConsistencyError('Original raw Matchete source hash changed: '+str(key))
            audit=compare_residue_and_log(residue,pole_free,stage)
            if (audit['pole_terms']!=ref['pole_terms'] or
                audit['log_terms']!=ref['log_terms'] or
                audit['constant_terms']!=ref['constant_terms']):
                raise ConsistencyError('Terms differ from prior counts')
            if len(terms(constant))!=ref['constant_terms']:
                raise ConsistencyError('Finite constant inventory changed')
            groups=finite_groups(constant,sector) if stage=='Y' and sector!='LLHX' else {}
            if sum(g['raw_additive_terms'] for g in groups.values()) != ref['constant_terms']:
                if ref['constant_terms']!=0:
                    raise ConsistencyError('Finite quartic partition incomplete')
            fragments={}
            for label in groups:
                selected=[t for t in terms(constant) if label in LAMBDA.findall(t)]
                if len(selected)!=groups[label]['raw_additive_terms']:
                    raise ConsistencyError('Finite quartic group cardinality mismatch')
                fragments[label]=' '.join(selected)
            planned.append((model,stage,sector,{**audit,'quartic_groups':groups,
                'model':model,'threshold':stage,'sector':sector,
                'source_sha256':ref['source_sha256'],
                'physical_matching_certified':False},fragments))
    # All inputs and models checked before any output is written.
    output_dir.mkdir(parents=True,exist_ok=True)
    detail=[]
    for model,stage,sector,audit,fragments in planned:
        sub=output_dir/model.replace('_-','_m')
        sub.mkdir(exist_ok=True)
        for label,fragment in fragments.items():
            (sub/f'{stage}_{sector}_{label}_finite_at_mass.inputform.txt').write_text(
                fragment+'\n',encoding='utf8')
        (sub/f'{stage}_{sector}_pole_log_audit.json').write_text(
            json.dumps(audit,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
        detail.append(audit)
        print(f'{model} {stage}/{sector}: pole=log coefficient PASS; '
              f'finite={audit["constant_terms"]}; quartic groups={len(fragments)}; '
              'physical=NOT CERTIFIED',flush=True)
    with (output_dir/'summary.csv').open('w',encoding='utf8',newline='') as f:
        fields=['model','threshold','sector','pole_terms','log_terms','constant_terms',
                'nonzero_unique_opaque_pole_monomials',
                'pole_equals_log_coefficient_in_opaque_factor_algebra',
                'physical_matching_certified']
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(detail)
    result={'format':'t3-bc-opaque-factor-pole-log-consistency-v2',
       'status':'FormalPoleLogIdentityCheckedNoPhysicalCountertermOrRGE',
       'models':list(models),'sector_results':detail,
       'pole_log_identity_is_exact_for_opaque_matchete_factors':True,
       'renormalized_matching_coefficient_obtained':False,
       'rge_threshold_decoupling_certified':False,
       'original_fermion_first_eomeft_equivalence_proven':False,
       'scalar_bubble_full_mass_resummation_certified':False,
       'numerical_neutrino_masses_changed':False,
       'physical_sequential_matching_certified':False}
    (output_dir/'summary.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    return result


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--all',action='store_true')
    g.add_argument('--model',choices=MODELS)
    ap.add_argument('--input-dir',type=Path,default=INPUT)
    ap.add_argument('--output-dir',type=Path,default=OUTPUT)
    ap.add_argument('--skip-original-source-validation',action='store_true',
                    help='Use only validated v1 split hashes, for archival replay when raw sources are unavailable')
    opt=ap.parse_args()
    inp=(ROOT/opt.input_dir).resolve()
    out=(ROOT/opt.output_dir).resolve()
    if not out.is_relative_to(ROOT):
        raise ConsistencyError('Output must remain under the project root')
    run(inp,out,MODELS if opt.all else (opt.model,),
        root=None if opt.skip_original_source_validation else ROOT)
    print('Pole-log consistency summary:',out/'summary.csv')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
