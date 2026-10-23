"""Audit symbolic Wilson-beta couplings against numerical intermediate-EFT names.

Requires the actual exported Matchete quartic seed and RGBeta payload.
Does not change the physics pipeline. ``λ -> lambdaH`` is allowed only after
checking the exported scalar-quartic term and the numerical state convention.
"""
from __future__ import annotations
import argparse
import json
import re
from pathlib import Path
import sympy as sp

NUMERIC_SYMBOLS = {
    'gY','g2','g3','yu','yd','ye','lambdaH',
    'lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12','lambdaT3',
    'lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj',
    'lambda12Adj','lambda12Cross',
    'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
    'lambdaS1bar2S2bar2','lambdaS1barS2S2bar2',
    'lambdaS1S1bar2S2bar','lambdaHHdagS1barS2barCross',
}

def symbols_in_matrix(validation):
    names=set()
    for expr in validation['anomalous_dimension_matrix'].values():
        names.update(s.name for s in sp.sympify(expr,locals={'conjugate':sp.conjugate}).free_symbols)
    return sorted(names)

def collect_quartic_coupling_terms(seed):
    result={}
    for term in seed.get('ScalarQuarticTerms', []):
        raw=term.get('TermInputForm','')
        for name in re.findall(r'Coupling\[([^,\]]+),',raw):
            result.setdefault(name,[]).append({'name':term.get('Name'), 'term':raw})
    return result

def validate(validation,quartic_seed,rgbeta):
    symbolic=symbols_in_matrix(validation)
    couplings=collect_quartic_coupling_terms(quartic_seed)
    available=set(rgbeta.get('report_betas',{}))
    if not available:
        raise ValueError('RGBeta payload has no report_betas')
    # This alias is intentionally checked for evidence rather than asserted.
    lambda_terms=couplings.get('λ',[])
    lambda_h_report='lambdaH' in available
    verified_higgs_alias=bool(lambda_terms and lambda_h_report and all(
        ('Field[H,' in t['term'] or 'Field[h,' in t['term'] or
         'Field[Phi,' in t['term']) for t in lambda_terms))
    # The raw seed may use other Higgs field names: fail closed rather than guess.
    mapping={'λ':'lambdaH'} if verified_higgs_alias else {}
    unknown=[x for x in symbolic if x not in NUMERIC_SYMBOLS and x not in mapping]
    possible_missing=[x for x in symbolic if x in NUMERIC_SYMBOLS and x not in available and x not in {'gY','g2','g3','yu','yd','ye'}]
    return {
        'status': ('CandidateMappingNeedsNormalizationCheck' if mapping and not unknown and not possible_missing else ('ReadyForNumericalEvaluator' if not unknown and not possible_missing else 'NeedsMappingReview')),
        'symbolic_symbols':symbolic,
        'quartic_seed_names':sorted(couplings),
        'rgbeta_report_beta_names':sorted(available),
        'candidate_mapping':mapping,
        'lambda_source_term_count':len(lambda_terms),
        'lambda_source_terms':lambda_terms[:4],
        'unknown_symbol_names':unknown,
        'missing_scalar_running_keys':possible_missing,
        'notes':[
            'The lambda alias requires confirmation of identical potential normalization, not merely field labels.',
            'Other quartic invariant names require an explicit mapping to RGBeta state variables.',
            'Three-generation Yukawa traces must be used when integrating full-flavor Wilson tensors.',
        ]
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--quartic-seed',type=Path,required=True)
    p.add_argument('--rgbeta',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    j=validate(*[json.loads(path.read_text(encoding='utf8')) for path in (a.validation,a.quartic_seed,a.rgbeta)])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(j,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps({k:v for k,v in j.items() if k not in ('lambda_source_terms','rgbeta_report_beta_names','quartic_seed_names')},indent=2,ensure_ascii=False))
    return 0 if j['status']=='ReadyForNumericalEvaluator' else 2
if __name__=='__main__':
    raise SystemExit(main())
