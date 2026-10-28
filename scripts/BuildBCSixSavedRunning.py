"""Presentation ONLY for the SIX saved, replay-validated B(-1)/C(-1) cases.

Reads the already existing Batch18 trajectory_samples.csv and validation report.
Produces interactive UV/EFT1 scalar-running HTML, vector PDF figure, and
source-verified Yukawa matrices for the six imported B/C fits. Also collects
B(+1)'s two accepted saved configurations when their chi2 matches the frozen CSV.
No RGE integration, matching, fitting, or historical output changes.

PowerShell, repository root:
    python -m scripts.BuildBCSixSavedRunning
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import html
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = Path('output/studies/general_scalar_extended/bulk_48/bulk48_v1/campaign_frozen.csv')
VALIDATION = Path('output/studies/general_scalar_extended/validation/bc_six_validation_v1')
EXPECTED = 'B3B9FE0EF2F3730AE49504B51CDD7AED4D4232E6812CBBFE0D16D78CD23ED3AF'
SCENARIOS = ['small_uv_frozen_5e-9','large_uv_frozen_0p5','small_uv_running_5e-9']


def read(root: Path):
    csvpath = root / FROZEN
    if hashlib.sha256(csvpath.read_bytes()).hexdigest().upper() != EXPECTED:
        raise ValueError('Frozen 48-study checksum mismatch')
    with csvpath.open(newline='',encoding='utf-8-sig') as f:
        frozen = {(r['model'],r['study']):r for r in csv.DictReader(f)}
    report = json.loads((root/VALIDATION/'report.json').read_text(encoding='utf-8'))
    if report.get('status')!='ReplayValidated' or report.get('candidate_count')!=6 or not report.get('perturbative_coordinate_screens_all_pass'):
        raise ValueError('Six-fit validation registry is not in verified state')
    with (root/VALIDATION/'trajectory_samples.csv').open(newline='',encoding='utf-8-sig') as f:
        samples = list(csv.DictReader(f))
    expected = {(m,s) for m in ('B_-1','C_-1') for s in SCENARIOS}
    fits = {(r['model'],r['case']):r for r in report['results']}
    if set(fits)!=expected or set((r['model'],r['case']) for r in samples)!=expected:
        raise ValueError('Wrong fit identities in six-case validated samples')
    for row in samples:
        if row['phase'] not in ('UV','EFT1'):
            raise ValueError('Unexpected trajectory phase')
        for key in ['mu_gev','scalar_light_tev','scalar_heavy_tev','rhoF','max_abs_quartic','max_yukawa_singular','quadratic_mu2_abs_gev2']:
            row[key]=float(row[key])
    for key, fit in fits.items():
        orig = frozen[key]
        if not math.isclose(float(fit['chi2']),float(orig['chi2']),rel_tol=1e-10,abs_tol=1e-14):
            raise ValueError('Validation and frozen fit chi2 differ: '+str(key))
        ends=[r for r in samples if (r['model'],r['case'],r['phase'])==(key[0],key[1],'EFT1')]
        if not ends:raise ValueError('EFT1 samples missing for '+str(key))
        last=min(ends,key=lambda z:z['mu_gev'])
        import ast
        masses=ast.literal_eval(orig['scalar_masses_muS_tev'])
        if any(not math.isclose(x,y,rel_tol=1e-9) for x,y in zip(
            (last['scalar_light_tev'],last['scalar_heavy_tev']),masses)):
            raise ValueError('EFT1 last sample fails frozen endpoint comparison: '+str(key))
    return frozen,fits,samples


def write_csv(out: Path, records: list[dict]):
    fields=['model','study','tensor','i','j','real','imag','source','provenance']
    with out.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\r\n')
        w.writeheader();w.writerows(records)


def yukawa_matrix(root: Path, row: dict, source: Path, provenance: str) -> list[dict]:
    data=json.loads(source.read_text(encoding='utf-8'))
    ordinary=data['base_state']['ordinary']
    out=[]
    for tensor in ('y1','y2'):
        mat=ordinary[tensor]
        for i in range(3):
            for j in range(3):
                out.append({'model':row['model'],'study':row['study'],'tensor':tensor,'i':i,'j':j,
                            'real':mat['real'][i][j],'imag':mat['imag'][i][j],
                            'source':source.relative_to(root).as_posix(),'provenance':provenance})
    return out


def collect_yukawas(root: Path, frozen: dict, six: dict) -> tuple[list[dict],list[dict]]:
    collected=[];meta=[]
    for key,r in sorted(six.items()):
        cfg=root/Path(r['cfg_path'].replace('\\','/'))
        fit=root/Path(r['fit_path'].replace('\\','/'))
        for path,sha in [(cfg,r['config_sha256']),(fit,r['fit_sha256'])]:
            if not path.is_file():raise FileNotFoundError('Validated historical source missing: '+str(path))
            if hashlib.sha256(path.read_bytes()).hexdigest().lower()!=sha.lower():
                raise ValueError('Validated historical source hash mismatch: '+str(path))
        collected.extend(yukawa_matrix(root,frozen[key],cfg,'Batch18 replay validated; source SHA-256 verified'))
        meta.append({'model':key[0],'study':key[1],'source':cfg.relative_to(root).as_posix()})
    for scenario in SCENARIOS[:2]:
        key=('B_+1',scenario)
        folder=root/'output/studies/general_scalar_extended'/scenario/'b_p1_three_fit_v1'
        cfg,fit=folder/'fitted_config.json',folder/'fit.json'
        if not cfg.is_file() or not fit.is_file():continue
        j=json.loads(fit.read_text(encoding='utf-8'))
        chi=(j.get('best_mass_safe') or {}).get('chi2')
        if j.get('status')!='CandidateConstrainedAccepted' or (j.get('model'),j.get('case'))!=key or chi is None or not math.isclose(float(chi),float(frozen[key]['chi2']),rel_tol=1e-10,abs_tol=1e-14):
            raise ValueError('B(+1) saved fit does not match frozen endpoint')
        collected.extend(yukawa_matrix(root,frozen[key],cfg,'Batch21 saved fit; chi2 verified against frozen'))
        meta.append({'model':key[0],'study':key[1],'source':cfg.relative_to(root).as_posix()})
    return collected,meta


def make_figure(samples:list[dict], out:Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes=plt.subplots(2,3,figsize=(12,6.8))
    for j,model in enumerate(('B_-1','C_-1')):
        for i,scenario in enumerate(SCENARIOS):
            ax=axes[j][i]
            for phase,color in [('UV','#ae764e'),('EFT1','#245a86')]:
                part=sorted([r for r in samples if (r['model'],r['case'],r['phase'])==(model,scenario,phase)],key=lambda r:r['mu_gev'])
                if not part:continue
                ax.plot([r['mu_gev'] for r in part],[r['scalar_light_tev'] for r in part],'-',lw=1.6,c=color,label=phase+' light')
                ax.plot([r['mu_gev'] for r in part],[r['scalar_heavy_tev'] for r in part],'--',lw=1.6,c=color,label=phase+' heavy')
            ax.set_xscale('log')
            ax.grid(alpha=.25)
            ax.set_xlabel('mu [GeV]',fontsize=9)
            if i==0:ax.set_ylabel('scalar quadratic mass [TeV]',fontsize=9)
            ax.set_title(model+' · '+['small frozen','large frozen','small UV running'][i],fontsize=10)
    handles,labels=axes[-1][-1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.025),ncol=4,fontsize=9)
    fig.suptitle('Saved candidate general-quartic scalar running — six B/C fits',y=1.085,fontsize=14)
    fig.text(.09,.015,'Source: Batch18 independently replay-validated trajectory_samples.csv. Only 6/48 cases; not pole masses or complete EFT.',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.955))
    fig.savefig(out.with_suffix('.pdf'),bbox_inches='tight')
    fig.savefig(out.with_suffix('.png'),dpi=170,bbox_inches='tight')
    plt.close(fig)


def make_html(samples:list[dict],out:Path):
    payload=json.dumps(samples,ensure_ascii=True,separators=(',',':')).replace('</','<\\/')
    prefix='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>T3 saved B/C running</title><style>
body{font:15px/1.55 system-ui,Arial;color:#16243a;background:#f5f7fa;margin:0}main{max-width:1130px;margin:auto;padding:24px}h1{font-size:26px}a{color:#1d558b}
.panel{padding:22px;border:1px solid #d7e0e8;border-radius:12px;background:#fff}label{display:inline-block;margin:6px 18px 12px 0}select{font:inherit;padding:7px;border:1px solid #a8b9ca;border-radius:7px}svg{width:100%;max-height:500px}p{max-width:100ch}small{color:#627181}
</style></head><body><main><p><a href="index.html">48-case overview</a></p><h1>Candidate general-quartic RGE replay: B(−1) and C(−1)</h1>
<p>Genuine saved Batch18 UV/EFT1 scalar running from six independently validated fits. This is <strong>6 of 48 cases</strong>; the frozen campaign does not contain matching scale-by-scale data for the other 42. Masses are running scalar quadratic eigenmasses, not pole masses. The general scalar coupling map and the one-scale C5 are conditional.</p>
<div class="panel"><label>Model <select id="model"><option>B_-1</option><option>C_-1</option></select></label>
<label>Scenario <select id="scenario"><option value="small_uv_frozen_5e-9">Small / UV frozen</option><option value="large_uv_frozen_0p5">Large / UV frozen</option><option value="small_uv_running_5e-9">Small / UV running</option></select></label>
<label>Quantity <select id="qty"><option value="scalar_light_tev">Light scalar mass [TeV]</option><option value="scalar_heavy_tev">Heavy scalar mass [TeV]</option><option value="max_abs_quartic">Max absolute quartic</option><option value="rhoF">rhoF</option></select></label><div id="chart"></div>
<small>Brown = UV, blue = EFT1. Frozen UV cases intentionally have no UV sample. Only stored points are plotted, joined by straight lines in log-scale.</small></div>
<p><a href="figures/bc6_saved_scalar_running.pdf">Vector-PDF six-case figure</a> · <a href="tables/bc6_saved_running.csv">Exact CSV samples</a></p>
'''
    js=r'''<script>const D=DATAJSON;const $=id=>document.getElementById(id);
function draw(){let m=$('model').value,s=$('scenario').value,k=$('qty').value;
 let phases=['UV','EFT1'].map(p=>D.filter(r=>r.model===m&&r.case===s&&r.phase===p).sort((a,b)=>a.mu_gev-b.mu_gev));
 let d=phases.flat(),svg=document.createElementNS('http://www.w3.org/2000/svg','svg');svg.setAttribute('viewBox','0 0 950 480');svg.setAttribute('role','img');svg.setAttribute('aria-label','Existing six-fit scalar running');
 function elem(tag,attrs={},t){let e=document.createElementNS('http://www.w3.org/2000/svg',tag);for(let x in attrs)e.setAttribute(x,attrs[x]);if(t!==undefined)e.textContent=t;return e}
 let ys=d.map(r=>r[k]),low=Math.min(...ys),high=Math.max(...ys),span=Math.max(high-low,Math.max(Math.abs(low),1)*.04);low-=span*.2;high+=span*.2;
 const x=v=>70+820*(Math.log10(v)-3)/4,y=v=>400-330*(v-low)/(high-low);
 svg.append(elem('line',{x1:70,x2:920,y1:400,y2:400,stroke:'#6b8198'}));svg.append(elem('line',{x1:70,x2:70,y1:25,y2:400,stroke:'#6b8198'}));
 for(let i=0;i<5;i++){let val=low+(high-low)*i/4,yy=y(val);svg.append(elem('text',{x:65,y:yy+4,'text-anchor':'end','font-size':11},val.toPrecision(4)));svg.append(elem('line',{x1:70,y1:yy,x2:920,y2:yy,stroke:'#e6edf5'}))}
 for(let i=3;i<=7;i++)svg.append(elem('text',{x:x(10**i),y:425,'text-anchor':'middle','font-size':12},'10^'+i));
 phases.forEach((p,i)=>{if(!p.length)return;let points=p.map(r=>x(r.mu_gev)+','+y(r[k])).join(' ');svg.append(elem('polyline',{points,fill:'none',stroke:i?'#245a86':'#ae764e','stroke-width':3}));
 p.forEach(r=>{let e=elem('circle',{cx:x(r.mu_gev),cy:y(r[k]),r:4,fill:i?'#245a86':'#ae764e'});e.append(elem('title',{},r.phase+' | '+r.mu_gev+' GeV | '+k+'='+r[k]));svg.append(e)})});
 svg.append(elem('text',{x:475,y:465,'text-anchor':'middle','font-size':13},'mu [GeV], log scale'));$('chart').replaceChildren(svg)}
 for(let id of ['model','scenario','qty'])$(id).addEventListener('change',draw);draw();</script></main></body></html>'''
    out.write_text(prefix+js.replace('DATAJSON',payload),encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();root=args.root.resolve()
    out=(args.output or root/'Reports/output/general_scalar_48_final').resolve()
    frozen,fits,samples=read(root)
    imported,meta=collect_yukawas(root,frozen,fits)
    (out/'figures').mkdir(parents=True,exist_ok=True)
    (out/'tables').mkdir(parents=True,exist_ok=True)
    make_figure(samples,out/'figures/bc6_saved_scalar_running')
    make_html(samples,out/'bc6_saved_running.html')
    with (out/'tables/bc6_saved_running.csv').open('w',encoding='utf-8',newline='') as f:
        keys=['case','model','run_id','phase','mu_gev','scalar_light_tev','scalar_heavy_tev','min_fermion_tev','rhoF','max_abs_quartic','max_yukawa_singular','quadratic_mu2_abs_gev2']
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\r\n');w.writeheader();w.writerows(samples)
    write_csv(out/'tables/imported_saved_yukawas.csv',imported)
    index = out/'index.html'
    if index.is_file():
        source_html=index.read_text(encoding='utf-8')
        anchor='<h2>New genuine RGE / sensitivity data not in the frozen CSV</h2>'
        addition=('<h2>Available new scalar-sector running subset (saved B/C validation)</h2>'
                  '<div class="panel"><p><a href="bc6_saved_running.html">Interactive B/C six-case UV/EFT1 running</a> · '
                  '<a href="figures/bc6_saved_scalar_running.pdf">Vector PDF</a> · '
                  '<a href="tables/imported_saved_yukawas.csv">Verified imported Yukawa matrices</a></p>'
                  '<p>Six independently replay-validated B(-1)/C(-1) benchmark samples only. Not 48/48; '
                  'no new sensitivity ovals or final Weinberg curves.</p></div>')
        if anchor in source_html and 'Available new scalar-sector running subset' not in source_html:
            index.write_text(source_html.replace(anchor,addition+anchor),encoding='utf-8')
    (out/'bc6_reporting_provenance.json').write_text(json.dumps({'input_frozen_sha256':EXPECTED,'validated_case_count':6,
        'saved_points':len(samples),'imported_yukawa_cases':meta,'physics_rerun':False,
        'scope':'6 validated B/C cases, UV/EFT1 saved scalar trajectories; B(+1) accepted Yukawa sources where available'},indent=2)+'\n',encoding='utf-8')
    print(f'Saved validated B/C RGE figure and interactive view: {out} | points: {len(samples)} | imported Yukawa case records: {len(meta)}')

if __name__=='__main__':main()
