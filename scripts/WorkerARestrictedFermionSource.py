"""Source-exact, all-operator restricted-Yukawa projection at F threshold.

Matchete raw UV->EFT1 source; NOT UV-minus-EFT1, a renormalised Wilson
boundary, or a claim that counterterms/operator normalisation are known.
Only imposes independent y1Conj=y2Conj=0 on complete raw expression.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import zipfile

SCHEMA = 'worker-a-restricted-all-operator-first-threshold-v1'
MODELS = ('A_-4', 'D_-2')
FORBIDDEN = frozenset(('y1Conj', 'y2Conj'))
TARGET = ('CandidateLeptonNumberViolatingLLScalarScalar', 'CandidateLeptonNumberViolatingLLHH')
COUPLING_TOKEN = re.compile(r'Coupling\[([A-Za-z_][A-Za-z0-9_]*)\s*,')


def sha(data: bytes) -> str:
    return sha256(data).hexdigest()


def require(ok: bool, message: str):
    if not ok:
        raise ValueError(message)


def find_unique(members: list[str], ending: str) -> str:
    found = [p for p in members if p.endswith(ending)]
    require(len(found) == 1, f'Expected one {ending}; found {len(found)}')
    return found[0]


def audit_model(z: zipfile.ZipFile, model: str) -> tuple[dict, dict[str, bytes]]:
    names = z.namelist()
    rawname = find_unique(names, f'/thesis_full_general_batch2_raw_v1/{model}/native_matching/data/eft1_after_F_one_loop_inputform.txt')
    ledgername = find_unique(names, f'/thesis_full_general_batch2_raw_v1/{model}/one_loop_operator_ledger.json')
    statusname = find_unique(names, f'/thesis_full_general_batch2_raw_v1/{model}/native_matching/native_status.json')
    rawbytes = z.read(rawname)
    raw = rawbytes.decode('utf-8-sig').strip()
    ledger = json.loads(z.read(ledgername))
    native = json.loads(z.read(statusname))
    require(ledger.get('source_sha256') == sha(rawbytes), 'Raw expression hash differs from ledger')
    require(ledger.get('status') == 'UnrenormalisedOneLoopTermsEnumerated', 'Unexpected ledger status')
    require(ledger.get('exact_reconstruction_by_expansion') is True, 'Ledger lacks exact reconstruction audit')
    require(ledger.get('only_fermions_integrated_out') is True, 'Non-fermionic first threshold')
    require(native.get('integrated_fields') == ['F'] and native.get('light_bsm_fields') == ['S1','S2'], 'Wrong heavy/light inventory')
    terms = ledger['terms']
    require(len(terms) == ledger['term_count'], 'Term count mismatch')
    require([t['index'] for t in terms] == list(range(1,len(terms)+1)), 'Missing source indices')
    kept, excluded = [], []
    reasons = Counter()
    for term in terms:
        expression = term['raw_inputform']
        lexical = set(COUPLING_TOKEN.findall(expression))
        recorded = set(term['coupling_symbols'])
        require(lexical == recorded, f'{model} term {term["index"]}: symbol inventory mismatch: {lexical ^ recorded}')
        absent = lexical & FORBIDDEN
        item = {'source_index': term['index'], 'sector':term['sector'],
                'source_term_sha256':sha(expression.encode('utf-8')),
                'has_hbar':term['contains_hbar'], 'has_explicit_pole':term['contains_dim_reg_pole'],
                'symbols':sorted(lexical)}
        if absent:
            item['set_to_zero'] = sorted(absent)
            excluded.append(item)
            for symbol in absent: reasons[symbol] += 1
        else:
            kept.append((item, expression))
    require(len(kept)+len(excluded) == len(terms), 'Projection partition incomplete')
    require(not any(set(item['symbols']) & FORBIDDEN for item,_ in kept), 'Forbidden coupling survived')
    # Each input term is an independently additive Matchete ledger entry.
    # Parenthesise each term so the sum is unambiguous irrespective of leading sign.
    expr = ' +\n'.join('('+s+')' for _,s in kept) + '\n'
    selected = [item for item,_ in kept if item['sector'] in TARGET]
    sector_counts = Counter(item['sector'] for item,_ in kept)
    excluded_counts = Counter(item['sector'] for item in excluded)
    metadata = {
        'schema':SCHEMA, 'model':model, 'scope':'all-operator unrenormalised F-only Matchete expression',
        'field_content':{'integrated':['F'], 'retained':['S1','S2']},
        'substitution':{'independent_y1Conj':0, 'independent_y2Conj':0},
        'restriction_method':'remove additive ledger terms containing exact Coupling[y1Conj,...] or Coupling[y2Conj,...] tokens',
        'original':{'ledger_source_path':rawname,'source_sha256':sha(rawbytes), 'ledger_sha256':sha(z.read(ledgername)), 'term_count':len(terms)},
        'restricted':{'expression_sha256':sha(expr.encode()), 'total_terms':len(kept), 'sector_counts':dict(sorted(sector_counts.items())),
                      'candidate_LNV_terms':len(selected),'candidate_LNV_by_sector':{k:sector_counts[k] for k in TARGET},
                      'hbar0_terms':sum(not item['has_hbar'] for item,_ in kept),
                      'explicit_pole_terms':sum(item['has_explicit_pole'] for item,_ in kept)},
        'excluded':{'total_terms':len(excluded),'sector_counts':dict(sorted(excluded_counts.items())), 'by_forbidden_symbol':dict(sorted(reasons.items()))},
        'excluded_terms':excluded, 'retained_term_metadata':[item for item,_ in kept],
        'native_status':native.get('status'),
        'uv_minus_eft1_soft_subtraction_performed':False,
        'counterterms_added':False,
        'loop_function_evaluation_performed':False,
        'operator_normalization_certified':False,
        'finite_msbar_wilson_coefficients_computed':False,
    }
    out = {'restricted_raw_inputform.txt':expr.encode(), 'restricted_manifest.json':(json.dumps(metadata,indent=2,allow_nan=False)+'\n').encode()}
    return {'model':model,'original':len(terms),'retained':len(kept),'excluded':len(excluded), 'retained_sectors':dict(sector_counts),'excluded_sectors':dict(excluded_counts)},out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources-zip',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True,help='Output ZIP (different from source ZIP)')
    args=p.parse_args()
    require(args.sources_zip.is_file(),'Missing source ZIP')
    require(args.sources_zip.resolve()!=args.output.resolve(),'Output may not replace source')
    results={}
    with zipfile.ZipFile(args.sources_zip) as src:
        require(src.testzip() is None,'Corrupt source archive')
        artifacts={}
        for model in MODELS:
            report, payload=audit_model(src,model)
            results[model]=report
            artifacts.update({model+'/'+name:data for name,data in payload.items()})
    artifacts['summary.json']=(json.dumps({'schema':SCHEMA, 'models':results,
        'finite_matching_ready':False},indent=2)+'\n').encode()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(args.output,'w',compression=zipfile.ZIP_DEFLATED) as dest:
        for name,data in artifacts.items():dest.writestr(name,data)
    for model,report in results.items():
        print(f'{model}: all sectors {report["original"]} => {report["retained"]} retained, {report["excluded"]} set to zero; retained LLSS={report["retained_sectors"].get(TARGET[0],0)}, LLHH={report["retained_sectors"].get(TARGET[1],0)}')
    print('Source-restricted raw Matchete EFT1, NOT renormalised one-loop Wilson matching')
    print('Wrote',args.output)

if __name__=='__main__':main()
