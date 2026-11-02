"""All 16 source-backed T3 lambdaT3-only conditional points, including 6 self-bilinears.

Physics is explicitly NOT certified: UV->all F(muF)->EFT1(S1,S2)->both S(muS)->SMEFT.
New output only in output/general_observable/lo_eft1_mixing/, never modifies v1.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))

from RGE.running.intermediate.All16PreliminaryNumerics import (
    SPLIT_MODELS, evaluate_model, summary_row, tag,
)
from RGE.running.intermediate.SixSelfBilinearPreliminary import evaluate_split_model, FORMAT
from Numerical.joint_scalar_matching.adapter import ALL_MODELS,normalize_model

SAFE=ROOT/'output/general_observable/lo_eft1_mixing'
SUMMARY_FORMAT='t3-all16-lambdaT3-conditional-self-bilinear-v2'


def build_summary(rows):
    passing=sum(r.get('status')=='ConditionalWeinbergLambdaT3OneChannel' for r in rows)
    nonward=sum(r.get('status')=='NonWeinbergTensorBlockedMassInterpretation' for r in rows)
    failed=sum(r.get('status')=='InputError' for r in rows)
    return {'format':SUMMARY_FORMAT,
        'status':'All16ConditionalOneChannel_NotPhysical' if passing==len(rows) else
                 'ConditionalPartialBatch_NotPhysical',
        'models_selected':len(rows),
        'models_with_Ward_valid_one_channel_results':passing,
        'models_non_Weinberg':nonward,'models_input_errors':failed,
        'models_blocked_special_scalar_basis':0,
        'six_self_bilinear_cartesian_hypotheses':sum(
            r.get('model') in SPLIT_MODELS and r.get('status')=='ConditionalWeinbergLambdaT3OneChannel'
            for r in rows),
        'source_to_Weyl_factor':None,'physical_neutrino_masses_computed':False,
        'full_seven_term_EFT1_running_included':False,
        'other_quartic_channels_included':False,
        'input_class_warning':('Eight ordinary models use illustrative non-synchronised scalar RG trajectories; '
             'B_-1 and C_-1 use recorded mixed-doublet bases; six self-bilinear models use a '
             'DECLARED native-Cartesian minus/plus map with assumed RGBeta/Matchete '
             'bilinear and lambdaT3 coupling normalization. All numerical masses are per '
             '|N_GammaCC| and are NOT certified physical neutrino masses.'),
        'EFT_sequence':'UV -> ALL F at muF -> EFT1(S1,S2) -> ALL scalars at ONE muS -> SMEFT',
        'models':rows}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',help='One model (default all 16)')
    p.add_argument('--output-root',type=Path,default=SAFE/'all16_lambdaT3_self_v2')
    p.add_argument('--check-only',action='store_true',help='Run input/numerical checks with no writes')
    p.add_argument('--require-all',action='store_true',help='Exit code 2 unless every model is Ward-valid')
    args=p.parse_args(argv)
    dest=args.output_root.resolve();base=SAFE.resolve()
    if dest==base or base not in dest.parents:
        p.error(f'Output must be a NEW directory under {base}')
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        p.error('Nonempty output directory: do not overwrite previous numerical results')
    models=[normalize_model(args.model)] if args.model else list(ALL_MODELS)
    docs=[];rows=[]
    for model in models:
        try:
            doc=evaluate_split_model(ROOT,model) if model in SPLIT_MODELS else evaluate_model(ROOT,model)
            row=summary_row(doc)
            if model in SPLIT_MODELS:
                row['self_bilinear_source_map_certified']=False
                row['self_bilinear_mass_bilinear_normalization_certified']=False
            docs.append(doc)
        except (OSError,ValueError,KeyError,IndexError,TypeError,RuntimeError) as exc:
            row={'model':model,'status':'InputError',
                 'reason':f'{type(exc).__name__}: {exc}',
                 'physical_mass_certified':False,'masses_ev_per_abs_NGamma':None}
        rows.append(row)
        print(f'{model}: {row["status"]}'+(f' — {row["reason"]}' if 'reason' in row else ''),flush=True)
    summary=build_summary(rows)
    if args.check_only:
        print(json.dumps(summary,indent=2,allow_nan=False))
    else:
        dest.mkdir(parents=True,exist_ok=True)
        for doc in docs:
            path=dest/tag(doc['model'])/'lambdaT3_one_channel.json'
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x',encoding='utf-8') as f:
                json.dump(doc,f,indent=2,allow_nan=False); f.write('\n')
        with (dest/'summary.json').open('x',encoding='utf-8') as f:
            json.dump(summary,f,indent=2,allow_nan=False);f.write('\n')
        with (dest/'summary.csv').open('x',encoding='utf-8',newline='') as f:
            w=csv.writer(f);w.writerow(('model','status','input_class','m1_ev_per_abs_NGamma',
                       'm2_ev_per_abs_NGamma','m3_ev_per_abs_NGamma','SU2_beta_residual',
                       'matching_scale_cancellation','physical_mass_certified'))
            for r in rows:
                m=r.get('masses_ev_per_abs_NGamma') or (None,None,None)
                w.writerow((r['model'],r['status'],r.get('input_class'),*m,
                    r.get('SU2_beta_residual'),r.get('matching_scale_cancellation'),False))
        print('Summary:',dest/'summary.json',flush=True)
    return 2 if summary['models_input_errors'] or (args.require_all and
                summary['models_with_Ward_valid_one_channel_results']!=len(models)) else 0


if __name__=='__main__':
    raise SystemExit(main())
