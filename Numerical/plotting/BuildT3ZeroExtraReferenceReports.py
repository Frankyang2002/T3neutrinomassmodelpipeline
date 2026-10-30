"""Exact existing optimised_rge.html canvas + original oval renderer, zero extras.

No scalar-RGE plotting, no visual redesign, no old fitted/64-case data substitution.
The existing frozen-64 reference-template adapter owns the HTML/CSS/JS.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from scripts import T3Frozen64Sources as policy
from scripts import T3ZeroExtraReportSources as zero
from Numerical.plotting import BuildT3FullScalarFrozen64HTML as old
from scripts import GenerateT3ZeroExtraWeinbergRunning as weinberg
from scripts import GenerateT3ZeroExtraGamma as gamma_runner

OUTPUT = zero.ROOT / 'Reports/output/display_zero_extra_all16'


def validated_weinberg(run_id: str, mode: str, src: dict, output_id: str,
                       summary: dict):
    path=zero.run_root(run_id)/output_id/mode/(policy.tag(src['cell']['model'])+'.json')
    if not path.is_file():
        return None
    d=policy.read_json(path)
    proof=d.get('source',{})
    policy.require(d.get('status')=='VerifiedSMWeinbergEndPointParity' and
                   d.get('mode')==mode and d.get('model')==src['cell']['model'] and
                   d.get('benchmark')==summary['benchmark'] and
                   proof.get('trajectory_sha256')==src['cell']['trajectory_sha256'] and
                   proof.get('neutrino_sha256')==src['cell']['neutrino_sha256'] and
                   proof.get('control_summary_sha256')==summary['summary_sha256'] and
                   proof.get('control_settings_sha256')==summary['settings_sha256'] and
                   proof.get('adapter_generator_sha256')==policy.digest(Path(weinberg.__file__)) and
                   proof.get('weinberg_solver_adapter_sha256')==policy.digest(Path(weinberg.old.__file__)),
                   f'{mode}/{src["cell"]["model"]}: stale/foreign Weinberg running')
    points=d.get('points',[])
    policy.require(len(points)>=3, 'too few numerical C5 running points')
    last=points[-1]
    policy.require(math.isclose(float(points[0]['mu_gev']),1e3) and
                   math.isclose(float(last['mu_gev']),1e2) and
                   all(float(a['mu_gev'])>float(b['mu_gev']) for a,b in zip(points,points[1:])),
                   'incorrect Weinberg scale interval')
    for key in policy.FIELDS:
        policy.require(math.isclose(float(last['observables'][key]),
                         float(src['neutrino']['oscillation_observables_normal_ordering'][key]),
                         rel_tol=3e-4,abs_tol=1e-9), 'final neutrino observable parity mismatch')
    return d


def validated_scan(run_id: str, mode: str, src: dict, scan_id: str, summary: dict):
    path=zero.run_root(run_id)/scan_id/mode/(policy.tag(src['cell']['model'])+'.csv')
    meta=path.with_suffix('.json')
    if not path.exists() and not meta.exists():
        return None
    policy.require(path.is_file() and meta.is_file(), 'partial gamma data: '+str(path))
    doc=policy.read_json(meta)
    cell=src['cell']
    expect={'model':cell['model'],'mode':mode,'benchmark':summary['benchmark'],
            'gamma':list(zero.NARROW_GAMMA),
            'trajectory_sha256':cell['trajectory_sha256'],
            'neutrino_sha256':cell['neutrino_sha256'],
            'source_summary_sha256':summary['summary_sha256'],
            'source_settings_sha256':summary['settings_sha256'],
            'generator_sha256':policy.digest(Path(gamma_runner.__file__)),
            'frozen_gamma_calculator_sha256':policy.digest(
                zero.ROOT/'scripts/GenerateT3FrozenFullScalarGamma.py')}
    policy.require(doc.get('source')==expect and
                   doc.get('csv_sha256')==policy.digest(path) and
                   doc.get('rows')==33 and doc.get('failures')==0,
                   'gamma source hash/grid mismatch; do not reuse another scan')
    with path.open(newline='',encoding='utf-8-sig') as stream:
        rows=list(csv.DictReader(stream))
    policy.require(len(rows)==33,'gamma CSV must have 33 physical evaluations')
    directions={d:[] for d in policy.DIRECTIONS}
    seen=set()
    for row in rows:
        d=row['direction']; g=float(row['gamma'])
        policy.require(d in directions and g in zero.NARROW_GAMMA and
                       (d,g) not in seen and row['model']==cell['model'] and
                       row['benchmark']==summary['benchmark'] and
                       row['status']=='NeutrinoEvaluated',
                       'duplicate/foreign/failed gamma point')
        seen.add((d,g))
        fields={k:float(row[k]) for k in policy.FIELDS}
        policy.require(all(math.isfinite(v) for v in fields.values()),
                       'nonfinite observable in gamma CSV')
        directions[d].append({'gamma':g,**fields})
    centre=[]
    for d in policy.DIRECTIONS:
        directions[d].sort(key=lambda p:p['gamma'])
        policy.require([p['gamma'] for p in directions[d]]==list(zero.NARROW_GAMMA),
                       'incorrect gamma grid (must be 0.95 through 1.05)')
        centre.append(next(x for x in directions[d] if x['gamma']==1.))
    saved=src['neutrino']['oscillation_observables_normal_ordering']
    for k in policy.FIELDS:
        for c in centre:
            policy.require(math.isclose(c[k],float(saved[k]),rel_tol=1e-6,abs_tol=1e-10),
                           'gamma=1 does not reproduce saved zero-extra central neutrino')
    return {'model':cell['model'], 'label':old.reference_model_name(cell['model']),
            'status':src['neutrino']['status'],
            'fit_stage':'historical_fitted_seed_zero_extra_'+mode,
            'centre':centre[0], 'directions':directions}


def sensitivity_html(mode: str, entries: list[dict], scan_id: str, summary: dict,
                     run_id: str):
    """Reuse original BuildFittedModelEnvelopes.HTML; no alternative canvas."""
    from Numerical.plotting import BuildFittedModelEnvelopes as original
    ref=policy.read_json(zero.ROOT/'docs/FROZEN_LAMBDA_TARGET_REFERENCE.json')
    target=ref['target']
    models=[];skipped=[]
    for src in entries:
        item=validated_scan(run_id,mode,src,scan_id,summary)
        if item is None:
            skipped.append({'model':src['cell']['model'],'status':'MissingNarrowGammaScan'})
        else:
            c=item['centre']
            item['chi2']=sum((c[k]-target['central_values'][i])**2/target['covariance'][i][i]
                             for i,k in enumerate(target['observable_names']))
            models.append(item)
    data={'models':models,'skipped':skipped,'target':target,
          'labels':{k:original.LABELS[k] for k in policy.FIELDS},
          'observables':list(policy.FIELDS),
          'source':f'zero-extra all16 {run_id}/{mode}/{scan_id}: new 0.95–1.05 T1/T2/T3 scans',
          'future_sigma':ref['future_sigma'],'future_source':ref['future_source'],
          'envelope_definition':'2x empirical sample covariance about saved gamma=1; descriptive, not confidence regions'}
    policy.require(original.HTML.count('__DATA__')==1, 'original sensitivity template changed')
    html=original.HTML.replace('__DATA__',json.dumps(data,allow_nan=False).replace('<','\\u003c'))
    html=html.replace('Individually fitted T3 model sensitivity',
                      'T3 zero-additional-scalar sensitivity — '+mode)
    html=html.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                      'Saved optimised fitted Yukawas; 33 new gamma evaluations per model in '
                      '0.95–1.05. '+mode+' control. No full sequential Wilson matching.')
    html=html.replace('Fitted centres', 'Candidate centres')
    return html,len(models)


def build(run_id='zero_extra_all16_low_v1', output_id='weinberg_running_v1',
          scan_id='gamma_095_105_v1', dest: Path|None=None,
          mode='boundary_zero'):
    for name in (output_id,scan_id):
        policy.require(bool(policy.RUN_ID.fullmatch(name)), 'invalid output/scan ID')
    policy.require(mode in (*zero.MODES, 'all'), 'invalid --mode')
    selected_modes = zero.MODES if mode == 'all' else (mode,)
    summary=zero.load_summary(run_id)
    destination=Path(dest) if dest else OUTPUT/run_id
    destination.mkdir(parents=True,exist_ok=True)
    counts={};coverage=[]
    nav=[]
    for mode in selected_modes:
        entries=[zero.sources(run_id,summary,mode,m) for m in policy.MODELS]
        records=[]
        for src in entries:
            w=validated_weinberg(run_id,mode,src,output_id,summary)
            records.append((src['cell'],{'neutrino':src['neutrino']},w))
            coverage.append({'model':src['cell']['model'],'mode':mode,
                'central_neutrino':True,'weinberg_running_verified':w is not None,
                'c12_intermediate_running':False,
                'sensitivity_verified':validated_scan(run_id,mode,src,scan_id,summary) is not None,
                'mass_scale_check':old.endpoint_screen(src['neutrino']),
                'matching_status':src['neutrino']['status']})
        folder=destination/mode
        folder.mkdir(parents=True,exist_ok=True)
        # EXACT user reference HTML/CSS/canvas from the existing adapter.
        rge=old.operator_html(mode+' '+summary['benchmark'],records)
        (folder/'interactive_rge.html').write_text(rge,encoding='utf-8')
        sensitivity,count=sensitivity_html(mode,entries,scan_id,summary,run_id)
        (folder/'sensitivity_ovals.html').write_text(sensitivity,encoding='utf-8')
        th,trs,columns,rows=old.table_rows(records)
        (folder/'all_results.html').write_text(old.make_page(
            mode+' — zero-extra central neutrino results',
            '<nav><a href="../index.html">Both modes</a>'
            '<a href="interactive_rge.html">Interactive RGE</a>'
            '<a href="sensitivity_ovals.html">Sensitivity</a></nav>'
            f'<div class="card scroller"><table><thead>{th}</thead><tbody>{trs}</tbody></table></div>'),encoding='utf-8')
        with (folder/'all_results.csv').open('w',encoding='utf-8',newline='') as f:
            wr=csv.writer(f);wr.writerow(columns);wr.writerows(rows)
        verified=sum(w is not None for _,_,w in records)
        counts[mode]={'weinberg':verified,'gamma':count,'models':16}
        nav.append(f'<div class="card"><h2>{old.esc(mode)}</h2><nav>'
                   f'<a href="{mode}/interactive_rge.html">Interactive RGE</a>'
                   f'<a href="{mode}/sensitivity_ovals.html">Sensitivity ovals</a>'
                   f'<a href="{mode}/all_results.html">Neutrino results</a>'
                   f'</nav><p>{verified}/16 C₅ running; {count}/16 narrow gamma scans</p></div>')
    with (destination/'coverage.csv').open('w',encoding='utf-8',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(coverage[0]));wr.writeheader();wr.writerows(coverage)
    proof={'source_run':run_id,'benchmark':summary['benchmark'],
           'source_summary_sha256':summary['summary_sha256'],
           'source_settings_sha256':summary['settings_sha256'],
           'mode_selection':list(selected_modes),
           'numerical_weinberg_output_id':output_id,'gamma_scan_id':scan_id,
           'reference_rge_template':'Numerical/plotting/templates/T3OptimisedReferenceRGE.html',
           'counts':counts,'intermediate_c12_running':False,
           'full_sequential_matching_certified':False}
    (destination/'provenance.json').write_text(json.dumps(proof,indent=2)+'\n',encoding='utf-8')
    intro=('<p>16 fitted T3 models in '+old.esc(', '.join(selected_modes))+
           '. The original topology coupling λT3 is nonzero. '
           'Uses the user-provided interactive RGE canvas and original sensitivity ovals.</p>'
           '<div class="card warning">Additional couplings are zero at 100 TeV. '
           '<b>boundary_zero</b> allows radiative generation; '
           '<b>held_zero</b> artificially removes their beta sources. '
           'Only the selected mode is displayed. Intermediate C₁₂(μ) not solved. '
           'Neutrino matching is candidate/conditional; no new Yukawa refit. '
           'Sensitivity γ = 0.95–1.05, Δγ = 0.01.</div>')
    (destination/'index.html').write_text(old.make_page(
        'T3 zero additional scalar couplings — interactive neutrino reports',
        intro+''.join(nav)+'<p><a href="coverage.csv">Coverage</a> · '
        '<a href="provenance.json">Provenance</a></p>'),encoding='utf-8')
    return proof


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='zero_extra_all16_low_v1')
    p.add_argument('--weinberg-output-id',default='weinberg_running_v1')
    p.add_argument('--scan-id',default='gamma_095_105_v1')
    p.add_argument('--output',type=Path)
    p.add_argument('--mode',choices=(*zero.MODES,'all'),default='boundary_zero')
    a=p.parse_args(argv)
    print(json.dumps(build(a.run_id,a.weinberg_output_id,a.scan_id,a.output,a.mode),indent=2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
