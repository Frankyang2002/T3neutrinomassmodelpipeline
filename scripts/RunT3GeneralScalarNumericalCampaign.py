"""Non-destructive all-study general-scalar *readiness and B/C refit* campaign.

Scope is intentionally NOT an all-16 general-potential C5 fitter. The current
mass-mixed Matchete evaluator supports only B_-1 and C_-1. All other model-case
cells are explicitly BLOCKED. Never substitute the historical restricted C5.

  python -m scripts.RunT3GeneralScalarNumericalCampaign --plan
  python -m scripts.RunT3GeneralScalarNumericalCampaign --run --campaign general_20261028
  python -m scripts.RunT3GeneralScalarNumericalCampaign --run --campaign general_fresh_1 --fresh

The runner does not touch historical output directories or pre-existing B/C
pilot, fit or continuation directories. Every new fit is a fresh run ID; reuse
of previously successful *exploratory* fits is explicit in the campaign report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('output/studies/general_scalar_extended')
CAMPAIGN_FORMAT = 't3-general-scalar-numerical-campaign-v1'
CUTOFF = 5.86
STUDIES = {
    'small_uv_running_5e-9': ('output/studies/lambda_5e-9_refit', True, 5e-9),
    'small_uv_frozen_5e-9': ('output/studies/fermion_start_lambda_comparison/5e-9', False, 5e-9),
    'large_uv_frozen_0p5': ('output/studies/fermion_start_lambda_comparison/0p5', False, 0.5),
}
MODEL_REP = {
    'A': (1,3,2), 'B': (2,2,1), 'C': (2,2,3),
    'D': (3,1,2), 'E': (3,3,2),
}
MODEL_NAMES = (
    'A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
    'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
    'E_+0','E_+2','E_-2','E_-4',
)
SUPPORTED = frozenset(('B_-1','C_-1'))
BLOCKER = ('No model-specific general-scalar C5 matching/fitting adapter. '
           'The historical restricted C5 must not be substituted and labelled general.')


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('General-scalar campaign: '+message)


def run_id(value: str) -> str:
    require(isinstance(value, str) and bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,34}', value)),
            'campaign ID must contain 1..35 safe characters')
    return value


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relpath(path: Path, root: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def historical_stem(model: str) -> str:
    require(model in MODEL_NAMES, 'unknown model')
    cls, alp = model.split('_')
    d1,d2,df = MODEL_REP[cls]
    a = int(alp)
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{"p" if a>=0 else "m"}{abs(a)}'


def validate_registry(root: Path) -> dict:
    path = root/'RGE/general/generated/general_rgbeta_quartic_registry.json'
    raw = json.loads(path.read_text(encoding='utf-8'))
    require(set(raw.get('models',{}))==set(MODEL_NAMES),'model inventory must be exactly 16')
    totals=raw.get('totals',{})
    require(totals.get('matchete_extra_quartics')==50 and
            totals.get('rgbeta_preexisting_structural')==12 and
            totals.get('rgbeta_new_vertices')==38 and
            totals.get('quadratic_couplings_blocked')==8,
            'general-scalar inventory unexpectedly changed')
    return {'path':relpath(path,root),'sha256':digest(path),'totals':totals}


def historical_row(root: Path, case: str, model: str) -> dict:
    folder, uv_on, lam = STUDIES[case]
    stem=historical_stem(model)
    cfg=root/folder/'fitted'/f'{stem}.json'
    fit=root/folder/'fits'/f'{stem}.json'
    row={'case':case,'model':model,'historical_uv_running':uv_on,
         'historical_lambdaT3':lam/2 if model.startswith('E_') else lam,
         'historical_fit_path':str(fit.relative_to(root).as_posix()),
         'historical_config_path':str(cfg.relative_to(root).as_posix()),
         'supported_general_c5':model in SUPPORTED,
         'status':'BlockedGeneralMatching' if model not in SUPPORTED else 'ReadyForExploratoryFit',
         'reason':BLOCKER if model not in SUPPORTED else
             'Conditional B/C candidate RGBeta + mixed one-scale Matchete C5 only; not physical.'}
    if not cfg.is_file() or not fit.is_file():
        row.update(status='MissingHistoricalInput',reason='Missing historical fit or fitted config')
        return row
    fitted=json.loads(fit.read_text(encoding='utf-8'))
    config=json.loads(cfg.read_text(encoding='utf-8'))
    require(fitted.get('model')==stem,'historical fit model mismatch: '+stem)
    rep=config.get('representation',{})
    expected={'d_s1':MODEL_REP[model[0]][0],
              'd_s2':MODEL_REP[model[0]][1],
              'd_f':MODEL_REP[model[0]][2],
              'alpha':int(model.split('_')[1]),'shared_scalar':False}
    require(rep==expected,'historical representation mismatch: '+stem)
    lt=config.get('base_state',{}).get('ordinary',{}).get('lambdaT3')
    require(isinstance(lt,dict) and
            math.isclose(float(lt.get('real',float('nan'))),row['historical_lambdaT3'],rel_tol=1e-12) and
            float(lt.get('imag',float('nan')))==0,
            'historical lambdaT3 mismatch: '+stem)
    chi=float(fitted['chi2_fit'])
    require(math.isfinite(chi) and chi>=0,'invalid historical chi2: '+stem)
    row.update(historical_chi2=chi,
               historical_fit_screen_pass=(fitted.get('status')=='Success' and chi<=CUTOFF),
               historical_fit_sha256=digest(fit),historical_config_sha256=digest(cfg))
    return row


def available_fits(root: Path, case: str, model: str, *,
                   quartic_fraction: float=0.1, mu2_fraction: float=0.0) -> list[dict]:
    """Only saved accepted candidate fits with the same scalar assumptions.

    This does not call the numerical evaluator or claim re-verification of chi2.
    """
    base=root/OUT/case
    tag=model.replace('-','m')
    records=[]
    if not base.is_dir():
        return records
    for path in sorted(base.glob(f'*/{tag}_exploratory_fit.json')):
        try:
            item=json.loads(path.read_text(encoding='utf-8'))
            if (item.get('model')!=model or item.get('case')!=case or
                item.get('status')!='ExploratoryAcceptedFit' or
                item.get('passes_approx_chi2_5p86') is not True or
                item.get('physical_matchete_rgbeta_map_certified') is not False or
                item.get('full_sequential_general_scalar_matching_performed') is not False):
                continue
            q=item.get('input_extra_quartic_fraction_of_abs_lambdaT3')
            u=item.get('input_mu2_fraction_of_min_scalar_m2')
            if q is None or u is None:
                source=item.get('source_manifest',{})
                q=source.get('quartic_fraction')
                u=source.get('mu2_fraction')
            if q is None or u is None or not math.isclose(float(q),quartic_fraction,abs_tol=1e-15) or not math.isclose(float(u),mu2_fraction,abs_tol=1e-15):
                continue
            provenance=item.get('source_manifest',{})
            hashes=provenance.get('files_sha256') or provenance.get('prior_input_file_sha256')
            if not isinstance(hashes,dict) or len(hashes)<1:
                continue
            valid_hashes=True
            for relative, expected_hash in hashes.items():
                if not isinstance(relative,str) or not isinstance(expected_hash,str):
                    valid_hashes=False; break
                source_path=(root/relative).resolve()
                if not source_path.is_relative_to(root.resolve()) or not source_path.is_file() or digest(source_path)!=expected_hash:
                    valid_hashes=False; break
            if not valid_hashes:
                continue
            chi=float(item['chi2_final'])
            if not math.isfinite(chi) or chi<0 or chi>CUTOFF:
                continue
            config=path.with_name(f'{tag}_exploratory_fitted_config.json')
            if not config.is_file():
                continue
            spec=json.loads(config.read_text(encoding='utf-8'))
            meta=spec.get('exploratory_general_scalar_fit',{})
            if meta.get('model')!=model or meta.get('case')!=case or meta.get('run_id')!=path.parent.name:
                continue
            records.append({'chi2':chi,'run_id':path.parent.name,
                            'fit_path':relpath(path,root), 'fit_sha256':digest(path),
                            'config_path':relpath(config,root),'config_sha256':digest(config)})
        except (OSError,ValueError,TypeError,KeyError,json.JSONDecodeError):
            # Never claim malformed/partial historical output is accepted.
            continue
    return sorted(records,key=lambda x:(x['chi2'],x['run_id']))


def inventory(root: Path, *, quartic_fraction: float=0.1,
              mu2_fraction: float=0.0) -> dict:
    root=Path(root).resolve()
    registry=validate_registry(root)
    rows=[historical_row(root,case,model) for case in STUDIES for model in MODEL_NAMES]
    for row in rows:
        if row['status']=='ReadyForExploratoryFit':
            recs=available_fits(root,row['case'],row['model'],
                                  quartic_fraction=quartic_fraction,mu2_fraction=mu2_fraction)
            if recs:
                row.update(status='ExistingExploratoryAccepted',
                           existing_fit=recs[0],
                           reason='Saved accepted conditional candidate result; not recomputed in inventory.')
    require(len(rows)==48 and sum(r['supported_general_c5'] for r in rows)==6,
            'unexpected 3-study × 16-model inventory')
    return {'format':CAMPAIGN_FORMAT,'status':'InventoryOnly',
            'registry':registry,'cases':list(STUDIES),
            'model_count':16,'study_count':3,'cell_count':48,
            'supported_cell_count':6,'general_c5_not_implemented_cell_count':42,
            'scalar_boundary':{'quartic_fraction_of_abs_lambdaT3':quartic_fraction,
                               'mu2_fraction_of_min_m2':mu2_fraction},
            'physical_coupling_map_certified':False,
            'sequential_general_scalar_matching_performed':False,
            'rows':rows}


def new_ids(campaign: str, model: str) -> tuple[str,str]:
    short=model.replace('-','m').replace('+','p')
    return f'{campaign}_{short}_fit',f'{campaign}_{short}_cont'


def validate_outputs_free(root: Path, campaign: str, inv: dict, *, fresh: bool) -> None:
    report=root/OUT/'campaigns'/campaign
    require(not report.exists(),'campaign report path already exists: '+str(report))
    for row in inv['rows']:
        if row['status']=='MissingHistoricalInput':
            require(False,'missing historical input: '+row['case']+' '+row['model'])
        if row['model'] not in SUPPORTED or (not fresh and row['status']=='ExistingExploratoryAccepted'):
            continue
        for name in new_ids(campaign,row['model']):
            dest=root/OUT/row['case']/name
            require(not dest.exists(),'new fit output path already exists: '+str(dest))


def write_report(root: Path, campaign: str, data: dict) -> None:
    out=root/OUT/'campaigns'/campaign
    require(out.is_dir(),'campaign directory not created')
    json_path=out/'campaign.json'; md_path=out/'summary.md'
    stats={}
    for row in data['rows']:
        stats[row['status']]=stats.get(row['status'],0)+1
    data['status_counts']=stats
    incomplete=any(r['status']=='ReadyForExploratoryFit' for r in data['rows'])
    bad=any(r['status'] in ('FitFailed','MissingHistoricalInput','NewExploratoryPoorFit') for r in data['rows'])
    data['status']='Running' if incomplete else ('CompletedWithFailuresOrPoorFits' if bad else 'CompletedWithMatchingBlocks')
    temp=json_path.with_name('campaign.json.tmp')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    temp.replace(json_path)
    lines=['# General-scalar numerical campaign', '',
           f'Campaign: `{campaign}`',
           'Three historical studies × 16 models = 48 model-study cells.',
           '**Only B(-1)/C(-1) (6 cells) have the current conditional mixed-scalar C5 fitter.**',
           'Other 42 cells remain explicitly blocked; this report is not 48 refits.',
           'Existing or newly accepted fits use candidate one-scale matching, not certified physical/sequential C5.', '',
           '| Study | Model | Status | Best exploratory chi2 | Source |',
           '|---|---|---|---:|---|']
    for row in data['rows']:
        accepted=row.get('new_fit') or row.get('existing_fit')
        chi=f'{accepted["chi2"]:.8g}' if accepted else '—'
        path=accepted.get('fit_path','—') if accepted else '—'
        lines.append(f'| {row["case"]} | {row["model"]} | {row["status"]} | {chi} | `{path}` |')
    lines+=['','## Status counts','']
    lines.extend(f'- {k}: {v}' for k,v in sorted(stats.items()))
    lines+=['','## Physics restrictions','',
            '- Input extra quartics are illustrative; the signed Matchete/RGBeta map is conditional.',
            '- New fits are at 1 TeV using the saved one-scale C5 coefficient, then Weinberg running to 100 GeV.',
            '- No newly certified 16-model general-potential matching, quadratic-operator basis or vacuum-stability result.',
            '- Neither historical studies nor earlier exploratory run directories are modified.','']
    md_path.write_text('\n'.join(lines),encoding='utf-8')


def run_campaign(root: Path,campaign: str,*,fresh: bool=False,
                 quartic_fraction:float=0.1,mu2_fraction:float=0.0,
                 sobol_points:int=8,max_nfev:int=60,fail_fast:bool=False) -> dict:
    """Execute only the six supported B/C candidate refits, with 42 explicit blocks."""
    root=Path(root).resolve()
    campaign=run_id(campaign)
    require(0 < quartic_fraction<=1 and 0 <= mu2_fraction<=0.2,
            'invalid scalar fractions')
    require(sobol_points>=2 and sobol_points&(sobol_points-1)==0,
            'Sobol sample must be a power of two >=2')
    require(10<=max_nfev<=500,'invalid max-nfev')
    inv=inventory(root,quartic_fraction=quartic_fraction,mu2_fraction=mu2_fraction)
    validate_outputs_free(root,campaign,inv,fresh=fresh)
    from scripts import RunT3GeneralScalarBCExploratoryFit as runner
    from scripts import ContinueT3GeneralScalarBCExploratoryFit as continue_fit
    require(root==runner.ROOT.resolve(),'campaign must run from installed project checkout')
    output=root/OUT/'campaigns'/campaign
    output.mkdir(parents=True,exist_ok=False)
    inv.update(status='Running',campaign=campaign,
               run_mode='fresh' if fresh else 'reuse_existing_accepted',
               physics_fitter='Batch11 + Batch12 (conditional B/C only)',
               historical_outputs_overwritten=False)
    write_report(root,campaign,inv)
    for row in inv['rows']:
        if row['model'] not in SUPPORTED:
            continue
        if row['status']=='ExistingExploratoryAccepted' and not fresh:
            print(f"{row['case']} {row['model']}: reuse accepted chi2={row['existing_fit']['chi2']:.8g}",flush=True)
            continue
        case,model=row['case'],row['model']
        fit_id,cont_id=new_ids(campaign,model)
        try:
            print(f'{case} {model}: starting new 18-real-Yukawa candidate fit',flush=True)
            first=runner.new_run(case,model,fit_id,
                quartic_fraction=quartic_fraction,mu2_fraction=mu2_fraction,
                sobol_points=sobol_points,seed=20261028,
                bound=3.0,sobol_fraction=0.5,sobol_floor=1e-6,
                seed_count=3,fit=True,root=root)
            result=first
            if not first['passes_approx_chi2_5p86']:
                print(f'{case} {model}: stage 1 chi2={first["chi2_final"]:.8g}; continue all 18',flush=True)
                result=continue_fit.continue_run(case,model,fit_id,cont_id,
                          max_nfev=max_nfev,root=root)
            tag=model.replace('-','m')
            output_id=result['run_id']
            path=root/OUT/case/output_id/f'{tag}_exploratory_fit.json'
            conf=path.with_name(f'{tag}_exploratory_fitted_config.json')
            require(path.is_file() and conf.is_file(),'fitter returned no saved result')
            score=float(result['chi2_final']); accepted=bool(result['passes_approx_chi2_5p86'])
            row.update(status='NewExploratoryAccepted' if accepted else 'NewExploratoryPoorFit',
                       reason='Candidate fit output, physical matching unverified',
                       new_fit={'chi2':score,'run_id':output_id,
                                'fit_path':relpath(path,root),'fit_sha256':digest(path),
                                'config_path':relpath(conf,root),'config_sha256':digest(conf)})
            print(f'{case} {model}: {row["status"]} chi2={score:.8g}',flush=True)
        except Exception as exc:
            row.update(status='FitFailed',reason=f'{type(exc).__name__}: {exc}')
            print(f'{case} {model}: FAILED: {exc}',flush=True)
            write_report(root,campaign,inv)
            if fail_fast:
                raise
        write_report(root,campaign,inv)
    return inv


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    action=p.add_mutually_exclusive_group(required=True)
    action.add_argument('--plan',action='store_true',help='Read-only inventory; no numerical solver')
    action.add_argument('--run',action='store_true',help='Reuse accepted B/C or refit missing; all other models blocked')
    p.add_argument('--campaign',help='New report/run identifier for --run')
    p.add_argument('--fresh',action='store_true',help='Force new fits of all six B/C, even if accepted fit exists')
    p.add_argument('--quartic-fraction',type=float,default=0.1)
    p.add_argument('--mu2-fraction',type=float,default=0.0)
    p.add_argument('--sobol-points',type=int,default=8)
    p.add_argument('--max-nfev',type=int,default=60)
    p.add_argument('--fail-fast',action='store_true')
    a=p.parse_args()
    if a.plan:
        require(not a.campaign and not a.fresh,'--plan cannot create or force a campaign')
        inv=inventory(ROOT,quartic_fraction=a.quartic_fraction,mu2_fraction=a.mu2_fraction)
        from collections import Counter
        stats=Counter(r['status'] for r in inv['rows'])
        for case in STUDIES:
            subset=[r for r in inv['rows'] if r['case']==case]
            print(f"{case}: " + ', '.join(f'{k}={v}' for k,v in sorted(Counter(r['status'] for r in subset).items())))
        print('All model-study cells: 48; supported conditional B/C C5: 6; general-C5 unsupported: 42')
        print('No files written; physical C5 certification: false')
        return
    require(bool(a.campaign),'--run requires --campaign')
    inv=run_campaign(ROOT,a.campaign,fresh=a.fresh,
            quartic_fraction=a.quartic_fraction,mu2_fraction=a.mu2_fraction,
            sobol_points=a.sobol_points,max_nfev=a.max_nfev,fail_fast=a.fail_fast)
    stats={}
    for row in inv['rows']:
        stats[row['status']]=stats.get(row['status'],0)+1
    print('Final campaign:',', '.join(f'{k}={v}' for k,v in sorted(stats.items())))
    print('Saved:',ROOT/OUT/'campaigns'/a.campaign/'summary.md')


if __name__=='__main__':
    main()
