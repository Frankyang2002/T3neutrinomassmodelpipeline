"""Read-only source-inventory and physical-matching gate audit.

Intended project path: scripts/AuditT3MatchingEvidenceArchive.py
Usage: python -X utf8 scripts/AuditT3MatchingEvidenceArchive.py --archive PATH --out PATH
This never modifies numerical output or asserts physical matching completeness.
"""
from __future__ import annotations
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import zipfile

MODELS = ('A_-4', 'D_-2')
SUFFIXES = {'A_-4': 'A_m4', 'D_-2': 'D_m2'}
ROOT = 'output/parallel_t3/'


def read_json(z, path):
    return json.loads(z.read(path).decode('utf-8-sig'))


def audit_archive(archive):
    result = {'format':'t3-batch2-source-evidence-gate-audit-v1',
              'archive_sha256':hashlib.sha256(Path(archive).read_bytes()).hexdigest(),
              'models': {}, 'physical_matching_certified':False}
    with zipfile.ZipFile(archive) as z:
        for model in MODELS:
            tag=SUFFIXES[model]
            prefix=ROOT+f'full_general_batch2/thesis_full_general_batch2_raw_v1/{model}/'
            ledger=read_json(z,prefix+'one_loop_operator_ledger.json')
            laurent=read_json(z,ROOT+f'full_general_batch2/thesis_full_general_batch2_laurent_v2/{model}/laurent_audit.json')
            ct=read_json(z,ROOT+f'eft1_msbar_counterterms_all16_batch6_v1/{tag}/counterterm_and_scale_ledger.json')
            pole=read_json(z,ROOT+f'native_pole_closure_all16_batch7_v1/{tag}/source_pole_and_matching_ledger.json')
            sectors=ledger['sector_counts']
            assert sum(sectors.values()) == ledger['term_count'], f'{model}: sector total mismatch'
            assert len(ledger['terms']) == ledger['term_count'], f'{model}: ledger length mismatch'
            assert sum(bool(t['contains_hbar']) for t in ledger['terms']) == laurent['raw_hbar1_term_count'], f'{model}: hbar count mismatch'
            assert ledger['term_count'] == laurent['original_source_term_count'], f'{model}: source count mismatch'
            assert laurent['uv_hard_minus_eft1_soft_subtraction_performed'] is False
            assert laurent['local_MSbar_counterterms_added'] is False
            assert ct['complete_one_loop_UV_EFT1_matching'] is False
            assert pole['complete_one_loop_UV_to_EFT1_matching'] is False
            gates=pole['physical_matching_certification']['gates']
            result['models'][model]={
                'operator_total':ledger['term_count'],
                'hbar_one_count':laurent['raw_hbar1_term_count'],
                'source_sectors':sectors,
                'pole_residue_term_count':laurent['raw_pole_residue_term_count'],
                'laurent_status':laurent['status'],
                'eft1_scalar_bubble_ct_cancelled':ct['EFT1_scalar_bubble_pole_cancelled_by_its_own_CT'],
                'eft1_scalar_pole_projection_residual':ct['EFT1_scalar_pole_projection_residual'],
                'native_raw_pole_vs_explicit_eft1_ct_relative':ct['native_raw_pole_vs_explicit_EFT1_CT_relative'],
                'physical_matching_gates':gates,
                'unresolved_gates':[k for k,v in gates.items() if not v],
                'status':'SOURCE_AUDITED_MATCHING_UNCERTIFIED',
                'matchete_ledger_sha256':hashlib.sha256(z.read(prefix+'one_loop_operator_ledger.json')).hexdigest(),
            }
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise SystemExit(f'Refusing to overwrite: {args.out}')
    result=audit_archive(args.archive)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('PASS: both source ledgers reconcile with Laurent counts')
    for model,r in result['models'].items():
        print(model, 'terms=',r['operator_total'],'one-loop=',r['hbar_one_count'],
              'poles=',r['pole_residue_term_count'],'unresolved=',len(r['unresolved_gates']))
    print('Matching certification: BLOCKED (no physical Wilson boundary injected)')
    print('Output:',args.out)

if __name__=='__main__':
    main()
