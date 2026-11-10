#!/usr/bin/env python3
"""Generate updated thesis Appendix B.1--B.4 from 32 fresh T3 EFT runs.

Use two verified result.json Yukawa benchmarks, transfer each unchanged across
all 16 models, and evaluate actual UV -> all F -> EFT1(S1,S2) -> all S -> SMEFT.

This uses the installed *exploratory* v3 projection and B/C portal convention.
It makes no claims of certified gauge-invariant one-loop predictions.
No physics result is invented or copied from thesis PDF.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys
import traceback
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.all16_rough_yukawa_fit import (
    BENCHMARK, YUKAWA_SCALE, FreshAll16Evaluator, all_models,
    physical_yukawas, preflight, read_json, save_json,
)
from Numerical.parallel_t3.background.contracts import tag
from Numerical.parallel_t3.rough_inputs import builder as builder_module
from Numerical.parallel_t3.threshold import preflight as preflight_module
from Numerical.parallel_t3.threshold import rough_numerics as numerics_module
from Numerical.parallel_t3.wilson import workflow, rough_stream2
from scripts.RunT3All16ExploratoryFit import (
    corrected_quartic_unit_tensors, consistent_project_c5,
    exploratory_preflight, exploratory_rough_case,
)

OBS_NAMES = (
    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)
FIELDS = ('reference','model','status','chi2', 'm1_eV','m2_eV','m3_eV',
          'dm21_1e-5_eV2','dm3l_1e-3_eV2', 'sin2_theta12','sin2_theta13',
          'sin2_theta23','theta12_deg','theta13_deg','theta23_deg',
          'epsilon_C5','k_re','k_im','r_weinberg_source',
          'r_weinberg_bubble','error')


def load_reference(model: str, directories: list[Path]) -> dict:
    """Choose lowest verified best chi2 from the explicitly supplied outputs."""
    candidates = []
    for directory in directories:
        f = directory / model / 'result.json'
        if not f.is_file():
            continue
        doc = read_json(f)
        best = doc.get('best')
        if not isinstance(best, dict):
            continue
        check = doc.get('verification', {})
        if check.get('passed') is not True:
            continue
        try:
            y1 = np.asarray(best['UV_y1'], float)
            y2 = np.asarray(best['UV_y2'], float)
            chi2 = float(best['chi2'])
        except (KeyError, TypeError, ValueError):
            continue
        if y1.shape != (3,3) or y2.shape != (3,3) or not (
            np.isfinite(y1).all() and np.isfinite(y2).all() and
            math.isfinite(chi2) and chi2 >= 0):
            continue
        candidates.append((chi2, f, y1, y2))
    if not candidates:
        raise ValueError(f'No verified {model} best point found in supplied directories')
    chi2, f, y1, y2 = min(candidates, key=lambda v: v[0])
    x = np.r_[y1.ravel(),y2.ravel()] / YUKAWA_SCALE
    if x.shape != (18,) or not np.isfinite(x).all():
        raise ValueError(f'Invalid 18-real benchmark for {model}')
    return {'x': x, 'y1': y1, 'y2': y2, 'source': str(f),
            'original_chi2': chi2}


def _cleanup_scratch(root: Path, token: str):
    if not (token.startswith('fit16_') and len(token) == 16 and
            token[6:].isalnum()):
        raise RuntimeError('Refusing to clean scratch with malformed token')
    base = root / 'output' / 'parallel_t3'
    for part in ('background','rough_inputs','rough_wilson','rough_results',
                 'fit16_adapted'):
        folder = base / part / token
        if folder.is_dir() and folder.name == token:
            shutil.rmtree(folder)


def _extract_numerics(root: Path, model: str, token: str) -> dict:
    folder = root / 'output/parallel_t3/rough_results' / token / tag(model)
    low_path = folder / 'C5_smeft_low.npz'
    at_ms_path = folder / 'C5_smeft_at_MS.npz'
    if not low_path.is_file() or not at_ms_path.is_file():
        raise FileNotFoundError(f'Missing actual SMEFT Wilson output: {low_path}')
    with np.load(low_path, allow_pickle=False) as z:
        c5 = np.asarray(z['C5_SMEFT_low'], np.complex128).copy()
    if c5.shape != (3,3) or not np.isfinite(c5).all():
        raise ValueError('Low-scale C5 must be finite complex 3x3')
    with np.load(at_ms_path, allow_pickle=False) as z:
        c5_eft1 = np.asarray(z['C5_EFT1_MS']).copy()
        c5_bubble = np.asarray(z['C5_scalar_bubble_MS']).copy()
        c5_final = np.asarray(z['C5_SMEFT_MS']).copy()
    obs = read_json(folder / 'observables.json')
    return {'c5': c5, 'c5_eft1_norm':float(np.linalg.norm(c5_eft1)),
            'c5_scalar_bubble_norm':float(np.linalg.norm(c5_bubble)),
            'c5_smeft_ms_norm':float(np.linalg.norm(c5_final)),
            'weinberg_source_residual':obs.get('source_Weinberg_relative_projection_residual'),
            'weinberg_bubble_residual':obs.get('scalar_bubble_Weinberg_relative_projection_residual'),
            'observables':obs}


def _complex_matrix(z: np.ndarray) -> dict:
    return {'real':z.real.tolist(), 'imag':z.imag.tolist()}


def _k_epsilon(current: np.ndarray, reference: np.ndarray) -> tuple[complex,float]:
    denom = np.vdot(reference,reference)
    if abs(denom) == 0:
        raise ValueError('Reference C5 is zero; proportionality undefined')
    k = np.vdot(reference,current) / denom
    norm = np.linalg.norm(current)
    if norm == 0:
        return complex(k), float('nan')
    return complex(k), float(np.linalg.norm(current-k*reference)/norm)


def _s(v, ndigits=9):
    if v is None:
        return '--'
    try:
        if not math.isfinite(float(v)):
            return '--'
        return f'{float(v):.{ndigits}g}'
    except (TypeError,ValueError):
        return '--'


def _tex_table(path: Path, caption: str, label: str, headers: list[str],
               values: list[list[str]]):
    n = len(headers)
    lines = [r'\begin{table}[htbp]',r'\centering',r'\small',
             r'\setlength{\tabcolsep}{3.4pt}',
             r'\begin{tabular}{' + 'l' + 'r'*(n-1) + '}',r'\hline',
             ' & '.join(headers)+r' \\',r'\hline']
    lines += [' & '.join(row)+r' \\' for row in values]
    lines += [r'\hline',r'\end{tabular}',
              r'\caption{' + caption + '}',r'\label{' + label + '}',
              r'\end{table}','']
    path.write_text('\n'.join(lines), encoding='utf-8')


def _export_tables(out: Path, results: dict, *, mu_low: float):
    rows=[]
    for ref in ('D_-2','A_-4'):
        for model in all_models():
            r=results.get(ref,{}).get(model,{})
            obs=r.get('observables',{})
            masses=obs.get('masses_eV') or [None,None,None]
            if len(masses)!=3:
                masses=[None,None,None]
            rows.append({'reference':ref, 'model':model,
               'status':r.get('status','NOT_RUN'), 'chi2':r.get('chi2'),
               'm1_eV':masses[0], 'm2_eV':masses[1], 'm3_eV':masses[2],
               'dm21_1e-5_eV2':_div(obs.get('dm21_sq_eV2'),1e-5),
               'dm3l_1e-3_eV2':_div(obs.get('dm31_sq_eV2'),1e-3),
               'sin2_theta12':obs.get('sin2_theta12'),
               'sin2_theta13':obs.get('sin2_theta13'),
               'sin2_theta23':obs.get('sin2_theta23'),
               'theta12_deg':obs.get('theta12_deg'),
               'theta13_deg':obs.get('theta13_deg'),
               'theta23_deg':obs.get('theta23_deg'),
               'epsilon_C5':r.get('epsilon_C5'),
               'k_re':r.get('k_re'),'k_im':r.get('k_im'),
               'r_weinberg_source':r.get('weinberg_source_residual'),
               'r_weinberg_bubble':r.get('weinberg_bubble_residual'),
               'error':r.get('error')})
    with (out/'all32_observables.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS)
        writer.writeheader();writer.writerows(rows)

    for ref,tab in (('D_-2','B1'),('A_-4','B2')):
        selection=[r for r in rows if r['reference']==ref]
        cols=['Model',r'$m_1$ (eV)',r'$m_2$ (eV)',r'$m_3$ (eV)',
              r'$\Delta m^2_{21}$ ($10^{-5}$ eV$^2$)',
              r'$\Delta m^2_{31}$ ($10^{-3}$ eV$^2$)',r'$\chi^2$',r'$\epsilon_{C_5}$']
        body=[[r['model'].replace('_','(')+')', *[_s(r[k]) for k in
               ('m1_eV','m2_eV','m3_eV','dm21_1e-5_eV2','dm3l_1e-3_eV2',
                'chi2','epsilon_C5')]] for r in selection]
        _tex_table(out/f'Table_{tab}.tex',
            f'New exploratory full-chain fixed-Yukawa transfer with the {ref.replace("_", "(") + ")"} UV benchmark. '
            f'All three heavy fermions are integrated out before simultaneous scalar matching; '
            f'SMEFT coefficients are evolved to {mu_low:g} GeV. '
            r'$\epsilon_{C_5}$ is measured against the reference model. '
            'Non-Weinberg projection and other conditional approximations apply.',
            f'tab:generalised_{tab.lower()}',cols,body)
    cols=['Reference','Model',r'$\sin^2\theta_{12}$',
          r'$\sin^2\theta_{13}$',r'$\sin^2\theta_{23}$',r'$\chi^2$']
    body=[[ref.replace('_','(')+')',r['model'].replace('_','(')+')',
           *[_s(r[k]) for k in ('sin2_theta12','sin2_theta13','sin2_theta23','chi2')]]
          for r in rows for ref in [r['reference']]]
    _tex_table(out/'Table_B3.tex',
       'PMNS mixing angles in the 32 full-chain fixed-Yukawa transfers; Yukawas '
       'are not reoptimised for the target representation. Exploratory matching.',
       'tab:generalised_b3',cols,body)
    cols=['Reference','Model',r'$\mathrm{Re}\,k$',r'$\mathrm{Im}\,k$',
          r'$\epsilon_{C_5}$']
    body=[[r['reference'].replace('_','(')+')',r['model'].replace('_','(')+')',
           _s(r['k_re']),_s(r['k_im']),_s(r['epsilon_C5'])] for r in rows]
    _tex_table(out/'Table_B4.tex',
       f'Complex Frobenius proportionality coefficients and residuals for the '
       f'low-scale $C_5$ matrices at {mu_low:g} GeV. Exploratory matching.',
       'tab:generalised_b4',cols,body)
    return rows


def _div(x,y):
    return None if x is None else float(x)/y


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--reference-results',type=Path,nargs='+',required=True,
       help='One or more EXISTING result-root directories, each with A_-4/result.json and/or D_-2/result.json')
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--batch1',default='output/batch1_general_contract/batch1_native_v1')
    ap.add_argument('--mu-low',type=float,default=100.)
    ap.add_argument('--wilson-solver',choices=['RK2','DOP853'],default='RK2')
    ap.add_argument('--wilson-steps',type=int,default=2)
    ap.add_argument('--mass-bilinears',choices=['zero','nonzero-exploratory'],default='zero')
    ap.add_argument('--gamma-factor',type=float,default=1.)
    ap.add_argument('--keep-trial-files',action='store_true')
    ap.add_argument('--debug',action='store_true')
    args=ap.parse_args(argv)
    root=ROOT.resolve()
    output=(root/args.output).resolve()
    allowed=(root/'output/parallel_t3').resolve()
    if allowed not in output.parents or output.exists():
        ap.error('--output must be a NEW subdirectory of output/parallel_t3')
    refdirs=[(root/p).resolve() for p in args.reference_results]
    for d in refdirs:
        if not d.is_dir():
            ap.error('Reference results directory not found: '+str(d))
    if args.wilson_steps<2 or args.mu_low<=0:
        ap.error('Invalid numerical controls')
    references={ref:load_reference(ref,refdirs) for ref in ('D_-2','A_-4')}
    cfg=read_json(root/BENCHMARK)
    missing={m:v['missing'] for m in all_models()
             for v in [preflight(root,m,batch1=args.batch1)] if not v['input_ready']}
    if missing:
        ap.error('Missing prerequisite source inputs: '+str(missing))
    output.mkdir(parents=True)
    results:dict[str,dict]={ref:{} for ref in references}
    provenance={'format':'t3-exploratory-generalised-thesis-tables-v1',
       'reference_benchmarks':{k:{'result_json':v['source'],'original_chi2':v['original_chi2'],
                              'UV_y1':v['y1'].tolist(),'UV_y2':v['y2'].tolist()}
                               for k,v in references.items()},
       'pipeline':'fresh UV RGBeta -> all heavy F at MF -> EFT1 two S -> one common MS -> SMEFT',
       'conditional_not_physical':True,
       'B_C_lambdaT3_portal_exploratory_correction':True,
       'L2_nonWeinberg_projection_exploratory':True,
       'finite_F_hard_C5_zero_benchmark':True,
       'reduced_scalar_potential':True,
       'mass_bilinears':args.mass_bilinears,
       'gamma_factor_assumed':args.gamma_factor,
       'wilson_solver':args.wilson_solver,'wilson_steps':args.wilson_steps,
       'mu_low_gev':args.mu_low,'runs_requested':32,
       'low_scale_C5_basis':'same Weyl/real-scalar basis as the repository',
       'ordering':'ascending singular masses; dm31_sq used for NO target'}
    save_json(output/'provenance.json',provenance)
    with patch.object(builder_module,'quartic_unit_tensors',corrected_quartic_unit_tensors),\
         patch.object(workflow,'project_c5',consistent_project_c5),\
         patch.object(rough_stream2,'project_c5',consistent_project_c5),\
         patch.object(preflight_module,'preflight_pair',exploratory_preflight),\
         patch.object(numerics_module,'run_rough_case',exploratory_rough_case):
        for ref in ('D_-2','A_-4'):
            for index,model in enumerate(all_models(),1):
                print(f'RUN {ref} -> {model} ({index}/16)',flush=True)
                evaluator=FreshAll16Evaluator(root,model,cfg,batch1=args.batch1,
                   gamma_factor=args.gamma_factor,mu_low=args.mu_low,
                   vev=float(cfg['vev_GeV']),reduced_potential=True,
                   mass_bilinears=args.mass_bilinears,
                   wilson_solver=args.wilson_solver,wilson_steps=args.wilson_steps,
                   keep_trial_files=True)
                token=None
                try:
                    trial=evaluator(references[ref]['x'].copy())
                    token=trial.meta.get('scratch_token')
                    if not token:
                        raise RuntimeError('Pipeline did not return scratch token for C5 extraction')
                    num=_extract_numerics(root,model,token)
                    obs=trial.observables
                    r={'status':'VALID_EXPLORATORY','model':model,
                       'reference':ref,'chi2':trial.chi2,
                       'residuals':trial.residuals,'observables':obs,
                       'c5_low':_complex_matrix(num['c5']),
                       'c5_eft1_ms_norm':num['c5_eft1_norm'],
                       'c5_scalar_bubble_ms_norm':num['c5_scalar_bubble_norm'],
                       'c5_smeft_ms_norm':num['c5_smeft_ms_norm'],
                       'weinberg_source_residual':num['weinberg_source_residual'],
                       'weinberg_bubble_residual':num['weinberg_bubble_residual'],
                       'meta':trial.meta}
                    if not args.keep_trial_files:
                        _cleanup_scratch(root,token)
                except Exception as exc:
                    r={'status':'BLOCKED','model':model,'reference':ref,
                       'error':f'{type(exc).__name__}: {exc}'}
                    if args.debug:
                        r['traceback']=traceback.format_exc(limit=14)
                    print('BLOCKED',ref,model,r['error'],flush=True)
                    # On pipeline failure token may be unknown; intentionally retain
                    # its scratch files rather than deleting unverified directories.
                results[ref][model]=r
                save_json(output/ref/f'{model}.json',r)
                save_json(output/'results.json',results)
                print('DONE',ref,model,r['status'],'chi2',r.get('chi2'),flush=True)
    for ref in ('D_-2','A_-4'):
        reference=results[ref].get(ref)
        if not reference or reference['status']!='VALID_EXPLORATORY':
            print(f'WARNING: no valid {ref} reference C5; proportionality unavailable',flush=True)
            continue
        a=reference['c5_low']; cref=np.asarray(a['real'])+1j*np.asarray(a['imag'])
        for model,r in results[ref].items():
            if r['status']!='VALID_EXPLORATORY':
                continue
            c=r['c5_low']; current=np.asarray(c['real'])+1j*np.asarray(c['imag'])
            try:
                k,eps=_k_epsilon(current,cref)
                r.update(k_re=float(k.real),k_im=float(k.imag),epsilon_C5=eps)
            except ValueError as exc:
                r['proportionality_error']=str(exc)
            save_json(output/ref/f'{model}.json',r)
    save_json(output/'results.json',results)
    table_rows=_export_tables(output,results,mu_low=args.mu_low)
    count=sum(r['status']=='VALID_EXPLORATORY' for ref in results.values() for r in ref.values())
    save_json(output/'summary.json',{'format':provenance['format'],'valid_runs':count,
       'requested_runs':32,'conditional_not_physical':True,
       'rows':[{'reference':r['reference'],'model':r['model'],'status':r['status'],
                'chi2':r['chi2'],'error':r['error']} for r in table_rows]})
    print(f'Finished {count}/32. Tables: {output}',flush=True)
    return 0 if count==32 else 1

if __name__=='__main__':
    raise SystemExit(main())
