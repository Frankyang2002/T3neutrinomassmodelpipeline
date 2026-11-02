"""Final Workstream3 conditional HHSS unit-sensitivity report (read-only inputs)."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from Numerical.lo_certification.completeness.response_closeout import (
    DEFAULT_CAPS, CloseoutError, inspect_run)

OWNED=ROOT/'output/general_observable/lo_certification_completeness'
SOURCE=OWNED/'portal_completion_02'
OUTPUT=OWNED/'final_sensitivity_closeout_01'


def inside_owned(path:Path)->Path:
    p=path.resolve()
    base=OWNED.resolve()
    if p==base or base not in p.parents:
        raise CloseoutError(f'Refusing path outside Workstream3: {p}')
    return p


def write_json_new(path:Path, obj:dict)->None:
    with path.open('x',encoding='utf-8') as f:
        json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False)
        f.write('\n')


def write_csv_new(path:Path,fields:tuple[str,...], rows:list[dict])->None:
    with path.open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        for row in rows:writer.writerow(row)


def render_report(report:dict)->str:
    models=report['models']
    caps=report['illustrative_unit_coefficient_caps']
    top=report['reference_sensitivity_ranking'][:8]
    lines=[
        '# Workstream 3 — final conditional numerical closeout',
        '',
        f"**Status:** `{report['status']}`",
        '',
        'This report describes **exploratory source-backed numerical responses**, not the',
        'complete physical one-loop SMEFT Weinberg coefficient. No physical',
        '`GammaCC` factor or extra HHSS named coefficient is assigned.',
        '',
        '## Source-backed completion',
        '',
        f"- Baseline completed: {report['baseline_completed_models']}/16 models.",
        f"- Missing baseline portal directions incorporated: {report['baseline_added_portal_directions']}/14.",
        f"- Extra named source symbols: {report['extra_named_symbols_in_inventory']}/24.",
        f"- Extra unit-reference directions: {report['extra_reference_directions']} total; "
        f"{report['valid_unit_reference_directions']} valid conditional projected responses, "
        f"{report['ward_rejected_reference_directions']} Ward-rejected, "
        f"{report['other_blocked_reference_directions']} other blocked.",
        '',
        '## Definition of the sensitivity numbers',
        '',
        r'For the saved baseline matrix $C_{5,0}$ and a unit-reference invariant $r$,',
        r'the code records $A_r=\partial C_5/\partial\operatorname{Re}c_r$',
        r'and $B_r=\partial C_5/\partial\operatorname{Im}c_r$ in the declared V3 convention.',
        r'For hypothetical $|c_r|\leq\kappa$,',
        '',
        r'$$\max_{|c_r|\leq\kappa}\frac{\|\operatorname{Re}c_r A_r+\operatorname{Im}c_r B_r\|_F}{\|C_{5,0}\|_F}',
        r'=\kappa\frac{\sqrt{\lambda_{\max}(G_r)}}{\|C_{5,0}\|_F},$$',
        '',
        r'where $(G_r)_{ij}=\operatorname{Re}\operatorname{Tr}(M_i^\dagger M_j)$ for $M=(A_r,B_r)$.',
        'The sum over reference directions is a triangle-inequality upper envelope for',
        'independent hypothetical unit-reference coefficients of that size. It is',
        '**not a physical upper limit** and need not use the correct Matchete/RGBeta',
        'phase or coupling basis. A missing/rejected direction makes the full-sector',
        'envelope unavailable rather than equal to zero.',
        '',
        '## Model comparison',
        '',
        '| Model | Valid/total extra refs | Largest per-unit response / baseline |',
    ]
    for cap in caps:
        lines[-1]+=f' Envelope $|c|\\leq {cap:g}$ |'
    lines.append('|---|---:|---:|'+'---:|'*len(caps))
    for m in models:
        valid=m['valid_reference_directions'];count=m['reference_directions']
        pieces=[f"| {m['model']} | {valid}/{count} | "
                f"{m['max_single_direction_relative_response_per_unit']:.4g} |"]
        env=m['illustrative_max_relative_frobenius_shift_at_caps']
        for cap in caps:
            value=env[f'{cap:g}']
            pieces.append(('not available' if value is None else f'{value:.4g}')+' |')
        lines.append(' '.join(pieces))
    lines.extend(['', 'These ratios are **dimensionless reference-basis sensitivities**;',
        'a value of 1 does not imply a 100% physical correction.', '',
        '## Largest valid single-reference sensitivities', '',
        '| Model | Reference | Unassigned source names | Per-unit relative norm |'
        ' Reference $|c|$ for 10% baseline response |',
        '|---|---|---|---:|---:|'])
    for row in top:
        ref=row['reference_id'];names=', '.join(row['unassigned_source_symbols'])
        ratio=row['max_relative_to_baseline_per_unit_coefficient']
        c=row['single_reference_coefficient_for_10pct_baseline']
        lines.append(f"| {row['model']} | {ref} {row['sector_counts']} | {names} | "
                     f"{ratio:.4g} | {f'{c:.4g}' if c is not None else 'n/a'} |")
    if not top:
        lines.append('| — | — | No valid nonzero unit directions | — | — |')
    lines.extend([
        '', '## Conclusion and limits', '',
        'The all-16 completed **baseline** portal calculation is an adequate',
        'exploratory numerical deliverable in its stated scope. Sensitivities above',
        'describe how strongly the saved unit-reference directions could contribute',
        '**if independently specified coefficients were supplied**. They cannot',
        'determine the actual relative physical impact without those coefficients.',
        '', '**Unresolved physical work (outside this numerical closeout):**', '',
    ])
    lines.extend('- '+x for x in report['missing_physics_not_resolved'])
    lines.extend(['', '**EFT ordering:** '+report['eft_sequence'], '',
        '**Final classification:** numerical Workstream 3 reported with explicit',
        'limitations; neither physical $C_5$ nor complete one-loop matching certified.',
        ''])
    return '\n'.join(lines)


def main(argv:list[str]|None=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=SOURCE,
                        help='Existing portal_completion_02 output (never modified)')
    parser.add_argument('--output-dir',type=Path,default=OUTPUT,
                        help='NEW Workstream3 child directory')
    parser.add_argument('--caps',default=','.join(f'{x:g}' for x in DEFAULT_CAPS),
                        help='Illustrative reference-coefficient magnitude caps, comma separated')
    parser.add_argument('--check-only',action='store_true',
                        help='Audit and print results without writing files')
    args=parser.parse_args(argv)
    try:
        src=inside_owned(args.input_dir)
        out=inside_owned(args.output_dir)
        caps=tuple(float(x) for x in args.caps.split(','))
        if not args.check_only:
            if out.exists():raise CloseoutError(f'Output already exists (no overwrite): {out}')
            if src==out or src in out.parents or out in src.parents:
                raise CloseoutError('Input and output cannot contain one another')
        report=inspect_run(src,caps)
        print(json.dumps({k:report[k] for k in (
            'status','model_count','baseline_completed_models','extra_named_symbols_in_inventory',
            'extra_reference_directions','valid_unit_reference_directions',
            'ward_rejected_reference_directions','other_blocked_reference_directions',
            'physical_C5_certified')},indent=2))
        if args.check_only:
            return 0
        out.mkdir(parents=True,exist_ok=False)
        write_json_new(out/'summary.json',report)
        fields=('model','status','baseline_dim5_frobenius_per_unknown_N_gamma_gev_inverse',
                'reference_directions','named_source_symbols','valid_reference_directions',
                'ward_rejected_reference_directions','other_blocked_reference_directions',
                'max_single_direction_relative_response_per_unit',
                'sum_single_direction_relative_response_per_unit')
        envelope_fields=tuple('illustrative_upper_relative_for_abs_c_le_'+f'{x:g}'.replace('.','p')
                              for x in caps)
        rows=[]
        for m in report['models']:
            row={k:m.get(k) for k in fields}
            for cap,key in zip(caps,envelope_fields):
                row[key]=m['illustrative_max_relative_frobenius_shift_at_caps'][f'{cap:g}']
            rows.append(row)
        write_csv_new(out/'model_sensitivity.csv',fields+envelope_fields,rows)
        ref_fields=('model','sector_counts','reference_id','unassigned_source_symbols',
                    'status','re_status','im_status','re_dim5_frobenius_per_unit',
                    'im_dim5_frobenius_per_unit','max_dim5_frobenius_per_unit_complex_coefficient',
                    'max_relative_to_baseline_per_unit_coefficient',
                    'single_reference_coefficient_for_1pct_baseline',
                    'single_reference_coefficient_for_10pct_baseline',
                    'single_reference_coefficient_for_100pct_baseline',
                    'max_saved_projection_residual')
        refs=[]
        for model in report['models']:
            for r in model['references']:
                row={k:r.get(k) for k in ref_fields}
                row['sector_counts']=':'.join(str(x) for x in r['sector_counts'])
                row['unassigned_source_symbols']=';'.join(r['unassigned_source_symbols'])
                refs.append(row)
        write_csv_new(out/'reference_sensitivity.csv',ref_fields,refs)
        (out/'WORKSTREAM3_CLOSEOUT.md').write_text(render_report(report),encoding='utf-8')
        print('SAVED: '+str(out))
        return 0
    except (ValueError,OSError,KeyError,TypeError,OverflowError,ArithmeticError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
