"""Symbolic leading-log C12 -> C5 -> C12 -> C5 feedback diagnostic.

Uses the *same* Matchete matched-C12 basis for both off-diagonal
anomalous dimensions. Inputs are previously exported symbolic EFT1 RGE
and the ProjectC5ToC12.py reverse-projection result.

This is NOT a full coupled EFT RGE solution. It isolates the off-diagonal
mixing feedback for constant anomalous dimensions over a log interval.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'output/full/comparison/smallY_smallL'
REVERSE = ROOT / 'Reports/output/analytical/GroupFactors/C12_reverse/c5_to_c12.json'
OUT = ROOT / 'Reports/output/analytical/GroupFactors/C12_feedback'


def expr(x):
    return sp.sympify(str(x), locals={'conjugate': sp.conjugate, 'sqrt': sp.sqrt, 'I': sp.I})


def clean(x):
    return sp.factor(sp.simplify(sp.expand(x)))


def as_text(x):
    return sp.sstr(clean(x))


def unit_matched_seed(x):
    """The reverse projection defines C12 basis by y1=y2=MF=1.

    Require only the matched seed symbols, not arbitrary replacements.
    """
    value = expr(x)
    replacements = {s: sp.S.One for s in value.free_symbols if str(s) in {'y1', 'y2', 'MF'}}
    value = clean(value.subs(replacements))
    if any(str(s) in {'y1', 'y2', 'MF'} for s in value.free_symbols):
        raise ValueError('Failed to remove the matched-seed prefactor')
    return value


def locate(model):
    folder = DATA / model
    hits = list(folder.rglob('eft1_wilson_rge.json'))
    if len(hits) != 1:
        raise FileNotFoundError(f'{model}: expected exactly one eft1_wilson_rge.json; found {len(hits)}')
    return hits[0]


def run(reverse_file, output, log_mu_low_over_high):
    payload = json.loads(reverse_file.read_text(encoding='utf-8-sig'))
    if not isinstance(payload, list):
        raise ValueError('Reverse projection JSON must be a list of model results')
    h = 1 / (16 * sp.pi**2)
    logarithm = sp.log(sp.sympify(log_mu_low_over_high)) if '/' in log_mu_low_over_high else sp.sympify(log_mu_low_over_high)
    results = []
    for record in payload:
        model = record['model']
        if int(record['residual_component_count']) != 0:
            status = 'Not closed: reverse projection has residuals'
        else:
            status = 'Single-C12 projection closed (reverse sector only)'
        rge_file = locate(model)
        rge = json.loads(rge_file.read_text(encoding='utf-8-sig'))
        check = rge.get('weinberg_subspace_validation', {})
        if check.get('matches_weinberg_subspace') is not True:
            raise ValueError(f'{model}: forward Weinberg projection not closed')
        gf = unit_matched_seed(check['beta_kappa_16pi2'])
        gr = clean(expr(record['beta_C12_from_C5']))
        # With initial C5=0 and constant off-diagonal gamma, the first
        # feedback term in C5 is h^3 L^3 gf^2 gr C12 / 3!.
        # Its ratio to the leading h L gf C12 is h^2 L^2 gf gr / 6.
        ratio = clean(h**2 * logarithm**2 * gf * gr / 6)
        results.append({
            'model': model,
            'gamma_5_from_12': as_text(gf),
            'gamma_12_from_5': as_text(gr),
            'feedback_to_leading_C5_ratio': as_text(ratio),
            'reverse_residual_count': int(record['residual_component_count']),
            'basis_status': status,
            'forward_source': str(rge_file.relative_to(ROOT)),
        })
        print(model, '| forward:', as_text(gf), '| reverse:', as_text(gr),
              '| residuals:', record['residual_component_count'], flush=True)
    output.mkdir(parents=True, exist_ok=True)
    (output/'leading_log_feedback.json').write_text(json.dumps({
        'method': 'constant off-diagonal anomalous dimensions, C5(mu_F)=0, leading feedback at O(hbar^3)',
        'log_mu_low_over_high': as_text(logarithm),
        'hbar': '1/(16*pi**2)',
        'warning': 'A nonzero reverse residual means a full multioperator basis is required; projected estimate is not a complete physical result.',
        'results': results,
    }, indent=2)+'\n', encoding='utf-8')
    with (output/'leading_log_feedback.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]) if results else ['model'])
        writer.writeheader()
        writer.writerows(results)
    print(f'Saved {len(results)} rows to {output}', flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reverse-json', type=Path, default=REVERSE)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--mu-high', type=sp.sympify, default=sp.Integer(100000))
    parser.add_argument('--mu-low', type=sp.sympify, default=sp.Integer(1000))
    args=parser.parse_args()
    if args.mu_high <= 0 or args.mu_low <= 0 or args.mu_low >= args.mu_high:
        parser.error('Require 0 < mu-low < mu-high')
    logval = str(sp.log(args.mu_low/args.mu_high))
    run(args.reverse_json,args.output,logval)

if __name__=='__main__':
    main()
