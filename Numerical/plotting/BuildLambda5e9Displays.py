"""Build linked/copyable HTML displays from saved 5e-9 comparison outputs.

Presentation-only: no fit, matching, ODE or sensitivity evaluations.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
STUDY = PROJECT_ROOT / 'output' / 'studies' / 'lambda_5e-9_refit' / 'comparison'
DEFAULT_DEST = PROJECT_ROOT / 'Reports' / 'output' / 'display_lambda_5e9'
MODES = {'optimised': 'Individually optimised Yukawas', 'shared_a_m2': 'Shared A(-2) Yukawas'}
IMAGE_EXTENSIONS = {'.png', '.jpg', '.jpeg', '.svg', '.webp'}


def esc(x):
    return html.escape(str(x), quote=True)


def page(title, contents):
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title><style>
:root{{font-family:system-ui,Arial,sans-serif;color-scheme:light dark}}
body{{margin:0 auto;max-width:1380px;padding:20px 26px;line-height:1.5}}
a{{color:#2987c6}} header{{border-bottom:1px solid #8885;margin-bottom:20px}}
nav a{{margin-right:20px}} .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:18px}}
figure{{margin:0;border:1px solid #8885;border-radius:9px;overflow:hidden;padding:9px}}
figure img{{width:100%;height:auto;display:block}} figcaption{{font-size:.9em;padding:6px}}
section{{margin:24px 0}} table{{border-collapse:collapse;font-size:.9em;width:100%}}
th,td{{border-bottom:1px solid #8885;padding:7px;text-align:left}}
small,.muted{{opacity:.75}} .alert{{padding:12px;border-left:3px solid #c99;background:#8881}}
</style></head><body><header><h1>{esc(title)}</h1><nav><a href="index.html">Study overview</a><a href="optimised.html">Optimised</a><a href="shared_a_m2.html">Shared A(-2)</a></nav></header>{contents}</body></html>'''


def read_json(path):
    with path.open(encoding='utf-8') as f:
        return json.load(f)


def copy_asset(src, output, relative):
    dst = output / relative
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return relative.as_posix()


def images_as_gallery(paths, source, output, mode):
    figures = []
    for src in paths:
        if src.suffix.lower() not in IMAGE_EXTENSIONS or not src.is_file():
            continue
        relative = Path('assets') / mode / src.relative_to(source)
        url = copy_asset(src, output, relative)
        figures.append(f'<figure><a href="{esc(url)}"><img loading="lazy" src="{esc(url)}" alt="{esc(src.stem)}"></a><figcaption>{esc(src.stem.replace("_", " "))}</figcaption></figure>')
    return '<div class="grid">' + ''.join(figures) + '</div>' if figures else '<p class="muted">No saved figures found.</p>'


def mode_page(mode, output):
    source = STUDY / mode
    if not source.exists():
        return '<div class="alert">Missing input directory: '+esc(source)+'</div>', []
    summary_path = source / 'summary.json'
    summary = read_json(summary_path) if summary_path.is_file() else {}
    rows = summary.get('models', [])
    records = []
    for r in rows:
        if isinstance(r, dict):
            records.append((r.get('model', ''), r.get('status', ''), r.get('chi2')))
    text = f'<p>Source: <code>{esc(source)}</code>. This display uses saved study data only; no new physics was evaluated.</p>'
    text += '<h2>All-model figures</h2>'
    top_images = sorted(p for p in source.glob('*.png'))
    text += images_as_gallery(top_images, source, output, mode)
    envelope = source / 'envelopes.html'
    if envelope.is_file():
        url = copy_asset(envelope, output, Path('assets')/mode/'envelopes.html')
        text += f'<p><strong><a href="{esc(url)}">Open interactive T1 / T2 / T3 envelope dashboard</a></strong></p>'
    else:
        text += '<p class="alert">Envelope HTML missing. Run BuildLambda5e9Envelopes display for this mode after its scans are complete.</p>'
    text += '<h2>Per-model trajectories and saved figures</h2>'
    model_paths = sorted(p for p in (source/'models').iterdir() if p.is_dir()) if (source/'models').is_dir() else []
    record_by_name = {r['model']:r for r in rows if isinstance(r,dict) and 'model' in r}
    for model_dir in model_paths:
        r = record_by_name.get(model_dir.name, {})
        chi2 = r.get('chi2')
        chi_text = f', χ²={chi2:.6g}' if isinstance(chi2, (int,float)) else ''
        imgs = sorted(p for p in (model_dir/'figures').rglob('*') if p.suffix.lower() in IMAGE_EXTENSIONS) if (model_dir/'figures').exists() else []
        text += f'<section id="{esc(model_dir.name)}"><h3>{esc(model_dir.name)} <small>({esc(r.get("status", "unknown"))}{esc(chi_text)})</small></h3>'
        text += images_as_gallery(imgs, source, output, mode)
        text += '</section>'
    if not records:
        text += '<p class="alert">No summary models found: check whether comparisons completed.</p>'
    return text, records


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=DEFAULT_DEST)
    args=ap.parse_args()
    output=args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    mode_rows={}
    for mode,label in MODES.items():
        body,rows=mode_page(mode,output)
        mode_rows[mode]=rows
        (output/(mode+'.html')).write_text(page(label,body),encoding='utf-8')
        print(f'{mode}: {len(rows)} saved model records -> {output/(mode+".html")}')
    names=sorted({r[0] for rows in mode_rows.values() for r in rows})
    by_mode={mode:{name:(status,score) for name,status,score in rows} for mode,rows in mode_rows.items()}
    content='''<p>Two UV-input comparisons at λ<sub>T3</sub>=5×10<sup>−9</sup> for classes A–D (class E: 2.5×10<sup>−9</sup>). This is a saved-data display, not a rerun.</p>'''
    content+='<div class="grid"><figure><h2><a href="optimised.html">Individually optimised</a></h2><p>Each target model uses its own fitted Yukawa matrices.</p></figure><figure><h2><a href="shared_a_m2.html">Shared A(-2)</a></h2><p>All 16 models use the fitted A(-2) Yukawa matrices.</p></figure></div>'
    content+='<h2>Cross-study summary</h2><table><thead><tr><th>Model</th><th>Optimised χ²</th><th>Shared A(-2) χ²</th></tr></thead><tbody>'
    for name in names:
        cells=[]
        for mode in MODES:
            status,value=by_mode[mode].get(name,('Missing',None))
            cells.append(f'{value:.6g}' if status=='Success' and isinstance(value,(int,float)) else status)
        content+='<tr><td>'+esc(name)+'</td>'+''.join('<td>'+esc(x)+'</td>' for x in cells)+'</tr>'
    content+='</tbody></table><p class="muted">χ² values refer to each target model’s five-observable NuFIT evaluation; these are not combined goodness-of-fit measures.</p>'
    (output/'index.html').write_text(page('T3 λ=5×10⁻⁹: all-16 displays',content),encoding='utf-8')
    print('Display index:',output/'index.html')

if __name__=='__main__':
    main()
