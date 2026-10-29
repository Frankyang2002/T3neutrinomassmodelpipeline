"""Overnight conditional 16 x 3 x 2 general-scalar Sobol / shared-Yukawa campaign.

Runs the repository's real 48-case Sobol + least-squares campaign, then the
existing shared A(-2) evaluator and a read-only HTML/CSV result report.
No claims of certified full general-scalar matching or gamma sensitivity.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BULK = ROOT / 'output/studies/general_scalar_extended/bulk_48'
SIX = ROOT / 'output/studies/general_scalar_extended/six_benchmarks'
MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3','C_+1','C_-1','C_-3',
          'D_+0','D_+2','D_-2','E_+0','E_+2','E_-2','E_-4')
CASES = ('small_uv_frozen_5e-9','large_uv_frozen_0p5','small_uv_running_5e-9')
FIELDS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')


def check_id(name: str) -> str:
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', name):
        raise ValueError('Invalid run identifier')
    return name


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copied_betas(source: Path, dest: Path) -> dict:
    """Re-use actual compiled native EFT beta sources, not a fresh Wolfram export."""
    inventory = read(ROOT / 'RGE/general/generated/general_rgbeta_quartic_registry.json')
    if len(inventory.get('models', {})) != 16:
        raise ValueError('Expected exactly 16 models in registry')
    required = [m for m in inventory['models'] if m not in ('B_-1','C_-1')]
    hashes = {}
    for model in required:
        stem = model.replace('+','p').replace('-','m')
        path = source / 'betas' / f'{stem}.json'
        if not path.is_file():
            raise FileNotFoundError('Missing saved native EFT beta: ' + str(path))
        data = read(path)
        if data.get('model') != model or data.get('status') != 'Success':
            raise ValueError('Incompatible saved native beta for ' + model)
        hashes[model] = sha(path)
    dest.mkdir(parents=True, exist_ok=True)
    for model in required:
        stem = model.replace('+','p').replace('-','m')
        original = source / 'betas' / f'{stem}.json'
        target = dest / f'{stem}.json'
        if target.exists():
            if sha(target) != hashes[model]:
                raise ValueError('Existing target beta differs: ' + str(target))
        else:
            shutil.copy2(original, target)
    return hashes


def run_logged(name: str, command: list[str], log: Path) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    print('STAGE', name, ' '.join(command), flush=True)
    with log.open('a', encoding='utf-8') as handle:
        handle.write('\n=== ' + datetime.now(timezone.utc).isoformat() + ' ' + name + ' ===\n')
        handle.flush()
        process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True,
                                   encoding='utf-8', errors='replace', bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end='', flush=True)
            handle.write(line)
            handle.flush()
        code = process.wait()
        handle.write('=== EXIT ' + str(code) + ' ===\n')
        handle.flush()
    return code


def render_report(csv_path: Path, html_path: Path) -> dict:
    with csv_path.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 96 or {(r['case'], r['model'], r['yukawa_mode']) for r in rows} != {
            (c,m,v) for c in CASES for m in MODELS for v in ('optimised','shared_a_m2')}:
        raise ValueError('Not a full 96-row six-benchmark summary')
    counts = {}
    for r in rows:
        key = r['yukawa_mode'] + ': ' + r['status']
        counts[key] = counts.get(key, 0) + 1
    headings = ('case','model','yukawa_mode','status','chi2','donor_status',*FIELDS)
    encoded = json.dumps([{k:r.get(k,'') for k in headings} for r in rows]).replace('<','\\u003c')
    summary = html.escape(json.dumps(counts, indent=2))
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>General-scalar six benchmark Sobol comparison</title>
<style>body{font:15px system-ui;background:#101725;color:#eff4fc;margin:2rem}a{color:#8bd9ff}input,select{padding:.45rem;background:#202e45;color:white;border:1px solid #56677f;border-radius:5px}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}td,th{border-bottom:1px solid #34425a;padding:.45rem;text-align:left;white-space:nowrap}th{cursor:pointer;position:sticky;top:0;background:#253550}tr:hover{background:#21334b}.scroll{overflow:auto;max-height:73vh}p{max-width:920px;line-height:1.55}pre{font-size:12px}</style>
<h1>General-scalar six benchmarks — Sobol / Yukawa comparison</h1>
<p>16 models × 3 scalar scenarios × 2 Yukawa prescriptions (96 records). Optimised results originate in the native 48-case Sobol + local fitting campaign; shared results insert same-scenario A(-2) Yukawa matrices into other models. An unavailable donor or failed evaluation is NOT a fit.</p>
<p><strong>Scientific limitation:</strong> these are conditional candidate general-quartic calculations, not certified full general-scalar matching. The 48-case source omits additional holomorphic quadratic operators; full sequential Wilson evolution and matched T1/T2/T3 sensitivity scans are not provided here. A poor-fit A(-2) donor, if explicitly allowed, is labelled accordingly.</p>
<p><a href="comparison.csv">CSV table</a> | <a href="comparison.json">JSON results</a> | <a href="a_m2_donors.json">Donor provenance</a></p>
<label>Filter <input id="filter" placeholder="model / scenario / status"></label>
<label>Mode <select id="mode"><option value="">Both</option><option>optimised</option><option>shared_a_m2</option></select></label>
<p id="count"></p><div class="scroll"><table><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div>
<h2>Counts</h2><pre>''' + summary + '''</pre>
<script>const H=''' + json.dumps(headings) + ''', D=''' + encoded + ''';let sort='case',reverse=false;
const head=document.getElementById('head'),body=document.getElementById('body');
for(const k of H){let th=document.createElement('th');th.textContent=k;th.onclick=()=>{reverse=(sort===k&&!reverse);sort=k;draw()};head.appendChild(th)}
function draw(){const q=document.getElementById('filter').value.toLowerCase(),m=document.getElementById('mode').value;let x=D.filter(r=>(!m||r.yukawa_mode===m)&&Object.values(r).join(' ').toLowerCase().includes(q));x.sort((a,b)=>{let n=Number(a[sort]),v=Number(b[sort]);let cmp=a[sort]!==''&&b[sort]!==''&&Number.isFinite(n)&&Number.isFinite(v)?n-v:String(a[sort]).localeCompare(String(b[sort]));return reverse?-cmp:cmp});body.replaceChildren();for(const r of x){let tr=document.createElement('tr');for(const k of H){let td=document.createElement('td');td.textContent=r[k]||'';tr.appendChild(td)}body.appendChild(tr)}document.getElementById('count').textContent=x.length+' records'}
document.getElementById('filter').oninput=draw;document.getElementById('mode').onchange=draw;draw();</script></html>'''
    html_path.write_text(page, encoding='utf-8')
    return counts


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=True, type=check_id)
    p.add_argument('--source-campaign', default='bulk48_v1', type=check_id)
    p.add_argument('--sobol-points', type=int, default=128)
    p.add_argument('--max-nfev', type=int, default=100)
    p.add_argument('--seed', type=int, default=20261030)
    p.add_argument('--extra-real', type=float, default=-0.01)
    p.add_argument('--allow-poor-a2', action='store_true', help='Evaluate poor-fit A(-2) donor as an explicitly labelled diagnostic')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--skip-fit', action='store_true', help='Run shared stage only from an existing completed bulk campaign')
    p.add_argument('--prepare-only', action='store_true', help='Check inputs and prepare saved native EFT betas without launching fits')
    args = p.parse_args(argv)
    if args.sobol_points < 2 or args.sobol_points & (args.sobol_points-1):
        p.error('--sobol-points must be a power of two >= 2')
    if not 0 <= args.max_nfev <= 1000:
        p.error('--max-nfev must be 0..1000')
    bulk_id = args.run_id + '_bulk48'
    six_id = args.run_id + '_six'
    base = BULK / bulk_id
    result = SIX / six_id
    log = SIX / (args.run_id + '_overnight.log')
    source = BULK / args.source_campaign
    if not args.skip_fit:
        if base.exists() and not args.resume:
            p.error('Destination already exists: use --resume or a new --run-id: ' + str(base))
        if not (source/'campaign.json').is_file():
            p.error('Source campaign unavailable: ' + str(source))
        hashes = copied_betas(source, base/'betas')
        signature={'source_campaign':args.source_campaign,'source_campaign_sha256':sha(source/'campaign.json'),
                   'beta_sha256':hashes,'sobol_points':args.sobol_points,
                   'max_nfev':args.max_nfev,'seed':args.seed,'extra_real':args.extra_real}
        sigfile=base/'overnight_settings.json'
        if sigfile.exists() and read(sigfile)!=signature:
            p.error('Saved overnight settings differ: use a NEW --run-id')
        if not sigfile.exists(): atomic(sigfile,signature)
        print('Prepared saved native EFT beta inputs for',len(hashes),'models',flush=True)
        if args.prepare_only: return 0
        code=run_logged('SOBOL_48', [sys.executable,'-u','-X','utf8','-m','scripts.RunT3GeneralScalar48Campaign',
            '--run-id',bulk_id,'--phase','fit','--resume',
            '--sobol-points',str(args.sobol_points),'--max-nfev',str(args.max_nfev),
            '--seed',str(args.seed),'--extra-real',str(args.extra_real)],log)
        if code: print('Bulk fit returned',code,'; inspect log and campaign summary.',flush=True)
        if not (base/'campaign.json').is_file(): return 2
    else:
        if args.prepare_only: p.error('--prepare-only and --skip-fit are incompatible')
        if not (base/'campaign.json').is_file(): p.error('Missing generated bulk campaign: '+str(base))
    cmd=[sys.executable,'-u','-X','utf8','-m','scripts.CompareT3GeneralScalarSharedA2',
         '--campaign',bulk_id,'--run-id',six_id,'--extra-real',str(args.extra_real)]
    if args.allow_poor_a2:cmd.append('--allow-poor-donor')
    if result.exists():
        if not args.resume: p.error('Shared result directory exists; use --resume')
        cmd.append('--resume')
    code=run_logged('SHARED_A2_48',cmd,log)
    if code:
        print('Shared stage failed; run with --resume to continue.',flush=True)
        return code
    csvfile=result/'comparison.csv'
    counts=render_report(csvfile,result/'index.html')
    atomic(result/'overnight_summary.json',{'run_id':args.run_id,'optimised_campaign':str(base.relative_to(ROOT)),
        'shared_report':str(result.relative_to(ROOT)),'status_counts':counts,
        'provenance':{'campaign_sha256':sha(base/'campaign.json')},
        'full_general_scalar_matching_certified':False,
        't1_t2_t3_sensitivity_calculated':False})
    print('DONE:',result/'index.html',flush=True)
    print('STATUS:',counts,flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
