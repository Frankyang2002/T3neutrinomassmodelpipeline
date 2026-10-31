#!/usr/bin/env python3
"""Verify native extended C5 expressions; optionally evaluate explicit Matchete-basis inputs."""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.T3ExtendedC5ExpressionReduction import MODELS, verify_model, self_test, evaluate_conditional_c5


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--zip',type=Path,help='Input uploaded expression archive; otherwise local native export directory')
    parser.add_argument('--input',type=Path,help='Optional explicitly declared Matchete numerical input JSON')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'Reports/output/workstream_a/extended_c5_expression_reduction')
    args=parser.parse_args()
    if args.self_test:
        self_test();return 0
    if not args.all and not args.input:
        parser.error('Choose --all, --input, or --self-test')
    source=args.zip if args.zip else ROOT/'Reports/output/workstream_a/extended_c5_native/native'
    args.output_dir.mkdir(parents=True,exist_ok=True)
    models=list(MODELS) if args.all else []
    obj=None
    if args.input:
        obj=json.loads(args.input.read_text(encoding='utf-8-sig'))
        model=obj.get('model')
        if model not in models:models.append(model)
    rows=[]
    for model in models:
        audit=verify_model(source,model)
        rows.append(audit)
        print(f"PASS {model}: group={audit['common_group_prefactor_exact']} | extra={audit['extra_yukawa']} | quartic={audit['extra_quartic_matchete']}",flush=True)
    summary={
        'status':'ExactOneStepMatcheteOrderedC5Reduction6of6' if len(rows)==6 else 'ExactOneStepMatcheteOrderedC5Reduction',
        'analysed_models':len(rows),'results':rows,
        'unindexed_MF_equal_indexed_MF_per_generation_assumed':True,
        'cross_package_physical_dictionary_certified':False,
        'hierarchical_matching_certified':False,
        'physical_neutrino_mass_prediction_certified':False,
    }
    (args.output_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if obj is not None:
        computed=evaluate_conditional_c5(obj['model'],obj)
        dest=args.output_dir/f"{obj['model']}_conditional_c5.json"
        dest.write_text(json.dumps(computed,indent=2)+'\n',encoding='utf-8')
        print(f'Conditional one-step Matchete-basis numerical result: {dest}',flush=True)
    print(f'Native one-step expression reduction: {len(rows)}/{len(models)} verified, report={args.output_dir / "summary.json"}',flush=True)
    return 0

if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as exc:
        print('BLOCKED:',type(exc).__name__+':',exc,file=sys.stderr)
        raise SystemExit(4)
