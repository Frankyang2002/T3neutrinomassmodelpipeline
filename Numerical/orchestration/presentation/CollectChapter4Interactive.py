"""Build a read-only Chapter 4 catalogue of EXISTING local T3 interactive/PNG results.

No physics is recalculated, no original files are moved or altered.  HTML links
point to the existing repository files; keep the catalogue in the repository.

PowerShell, from repository root:
    python -m Numerical.orchestration.presentation.CollectChapter4Interactive
    python -m Numerical.orchestration.presentation.CollectChapter4Interactive --strict
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import html
import json
import os
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[3]
SCENARIO = 'smallY_smallL'
BASE = Path('Reports/output/full')
DISPLAY = Path('Reports/output/display')
OUTPUT = Path('Reports/thesis_chapter4/interactive')
MODEL_DIMS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3), 'D': (3, 1, 2), 'E': (3, 3, 2)}
SELECTED = ('A-m2', 'A-p0', 'B-m3', 'B-p1', 'D-m2', 'D-p0')


def model_dir(short: str) -> str:
    cls, signed = short.split('-')
    d1, d2, df = MODEL_DIMS[cls]
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{signed}'


@dataclass(frozen=True)
class Item:
    section: str
    title: str
    path: Path
    kind: str
    purpose: str
    scope: str
    recommendation: str


def candidates() -> list[Item]:
    items = [
        Item('4.1–4.2', 'Common-input RGE dashboard (16 models)', DISPLAY / 'common_input_rge_smallY_smallL.html',
             'interactive', 'Model selectors and saved scale-dependent observables; use as a companion to Figures 4.1–4.2.', 'common-input', 'primary-interactive'),
        Item('4.1–4.2', 'Common-input model sensitivity envelopes', DISPLAY / 'common_input_model_envelopes_all.html',
             'interactive', 'Explore common-input prediction-space trajectories; distinguish sensitivity scans from baseline comparison.', 'common-input', 'supplementary'),
        Item('4.3', 'Saved heavy-fermion scan: solar splitting', BASE / 'common_input_fermion_scale' / SCENARIO / 'plots/delta_m21_sq_ev2_running_vs_frozen.png',
             'image', 'Full sixteen-model running and frozen mass-scale scan; reference for selecting six curves.', 'common-input', 'reference-only'),
        Item('4.3', 'Saved heavy-fermion scan: atmospheric splitting', BASE / 'common_input_fermion_scale' / SCENARIO / 'plots/delta_m3l_sq_ev2_running_vs_frozen.png',
             'image', 'Existing all-model atmospheric comparison.', 'common-input', 'reference-only'),
        Item('4.4', 'Saved Yukawa scan: solar splitting', BASE / 'common_input_coupling_scan' / SCENARIO / 'plots/yukawa_delta_m21_sq_ev2_running_vs_frozen.png',
             'image', 'Existing five-factor running/frozen Yukawa response (all models).', 'common-input', 'reference-only'),
        Item('4.4', 'Saved Yukawa scan: atmospheric splitting', BASE / 'common_input_coupling_scan' / SCENARIO / 'plots/yukawa_delta_m3l_sq_ev2_running_vs_frozen.png',
             'image', 'Existing five-factor atmospheric response (all models).', 'common-input', 'reference-only'),
        Item('4.4', 'Saved scalar-coupling scan: solar splitting', BASE / 'common_input_coupling_scan' / SCENARIO / 'plots/lambda_delta_m21_sq_ev2_running_vs_frozen.png',
             'image', 'Existing five-factor lambdaT3 response (all models).', 'common-input', 'reference-only'),
        Item('4.4', 'Saved scalar-coupling scan: atmospheric splitting', BASE / 'common_input_coupling_scan' / SCENARIO / 'plots/lambda_delta_m3l_sq_ev2_running_vs_frozen.png',
             'image', 'Existing five-factor atmospheric response (all models).', 'common-input', 'reference-only'),
        Item('Supplement (not 4.1–4.5)', 'Independently fitted RGE dashboard', DISPLAY / 'fitted_rge_smallY_smallL.html',
             'interactive', 'Fitted results have model-dependent inputs; never combine their numerical values with common-input tables.', 'fitted', 'supplementary-only'),
        Item('Supplement (not 4.1–4.5)', 'Fitted running-stage comparison', DISPLAY / 'study_comparisons/smallY_smallL/StudyA/interactive_comparison.html',
             'interactive', 'Exploratory fitted comparison, not a common-input representation scan.', 'fitted', 'supplementary-only'),
        Item('Supplement (not 4.1–4.5)', 'Fitted comparison dashboard', DISPLAY / 'comparison/smallY_smallL/interactive_comparison.html',
             'interactive', 'Different study and input prescription: use for context only.', 'fitted', 'supplementary-only'),
        Item('Supplement (not 4.1–4.5)', 'Fitted parameter correlation dashboard', DISPLAY / 'interactive_correlations.html',
             'interactive', 'Sensitivity/correlation visualisation; outside the Chapter 4 common-input figure series.', 'fitted-or-other', 'supplementary-only'),
    ]
    base = BASE / 'common_input_rge' / SCENARIO
    for short in SELECTED:
        d = base / model_dir(short) / 'figures'
        for section, name, purpose in (
            ('4.5', 'c5_threshold_contributions.png', 'Existing magnitudes |C5_ij| for hard, direct and combined; not signed interference.'),
            ('4.5', 'intermediate_direct_weinberg_running.png', 'Direct intermediate-EFT Weinberg generation, saved trace.'),
            ('4.1–4.2', 'neutrino_mass_splitting_running.png', 'Low-energy running of both mass-squared splittings.'),
            ('4.1–4.2', 'c5_running.png', 'Final-EFT |C5_ij| running, useful explanatory supplementary plot.'),
        ):
            items.append(Item(section, f'{short}: {name.removesuffix(".png").replace("_", " ")}', d / name,
                              'image', purpose, 'common-input', 'selected-six-supplement'))
    return items


def relative_url(target: Path, output: Path) -> str:
    rel = os.path.relpath(target, start=output.parent)
    return '/'.join(quote(part) for part in Path(rel).parts)


def render(items: list[dict], out: Path) -> str:
    cards = []
    for i, x in enumerate(items):
        url = relative_url(x['absolute'], out)
        if x['exists']:
            preview = (f'<img loading="lazy" src="{html.escape(url)}" alt="{html.escape(x["title"])}">'
                       if x['kind'] == 'image' else '<div class="html-icon">Interactive HTML<br><small>Open or preview below</small></div>')
            action = (f'<a href="{html.escape(url)}" target="_blank" rel="noopener">Open original ↗</a>'
                      f'<button type="button" onclick="preview({i})">Preview here</button>')
        else:
            preview = '<div class="missing">Not found in local repository</div>'
            action = '<span class="missing">Missing; no substitute generated</span>'
        cards.append(f'''<article class="card" data-scope="{html.escape(x['scope'])}" data-section="{html.escape(x['section'])}">
        {preview}<div class="content"><div class="meta">{html.escape(x['section'])} · {html.escape(x['scope'])} · {html.escape(x['recommendation'])}</div>
        <h3>{html.escape(x['title'])}</h3><p>{html.escape(x['purpose'])}</p>
        <div class="actions">{action}</div><small class="path">{html.escape(x['path'])}</small></div></article>''')
    jsdata = json.dumps([{'url': relative_url(x['absolute'], out), 'kind': x['kind'], 'title': x['title'], 'exists': x['exists']} for x in items]).replace('</', '<\\/')
    return '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Chapter 4 | existing interactive figure catalogue</title>
    <style>body{background:#f7f8fa;color:#182433;font:15px/1.5 system-ui,Arial;margin:0}header{background:#152a45;color:white;padding:28px max(22px,5vw)}main{max-width:1320px;margin:auto;padding:20px}.warn{background:#fff1d6;padding:12px 16px;border-left:4px solid #a96a05} .controls{display:flex;gap:12px;flex-wrap:wrap;margin:20px 0}.controls label{display:flex;flex-direction:column;font-weight:600}select{padding:9px;min-width:190px}#cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(305px,1fr));gap:16px}.card{border:1px solid #d6dde7;background:white;border-radius:10px;overflow:hidden;display:flex;flex-direction:column}.card img{width:100%;height:175px;object-fit:contain;background:white}.content{padding:15px}.meta{font-size:12px;color:#4d6075}.content h3{margin:7px 0}.content p{margin:4px 0 12px}.actions{display:flex;gap:12px;align-items:center;flex-wrap:wrap}.actions a,.actions button{font:inherit;color:#084e9a;cursor:pointer;background:transparent;border:0;padding:0;text-decoration:underline}.path{display:block;overflow-wrap:anywhere;color:#65768a;margin-top:12px}.html-icon,.missing{padding:55px 10px;text-align:center;background:#edf2f7;color:#345}.missing{background:#fff3ee;color:#8b4425}#viewer{display:none;margin:20px 0;background:white;padding:16px;border:1px solid #cad2db;border-radius:10px}#viewer iframe{width:100%;height:690px;border:1px solid #ddd}#viewer img{max-width:100%}button.close{float:right;cursor:pointer}footer{color:#596b7d;margin:28px 0}</style></head><body>
    <header><h1>Chapter 4 · existing figure catalogue</h1><p>Read-only collection of saved T3 dashboards and diagnostic figures. Interactive HTML opens from its original location.</p></header><main>
    <p class="warn"><strong>Provenance:</strong> Only common-input items are candidates for Chapter 4 comparisons. Fitted/other dashboards are segregated for supplementary browsing. Existing 16-model scans are not substitutes for the required six-model thesis figures. Absolute C5 entries are not signed interference terms.</p>
    <div class="controls"><label>Input class<select id="scope"><option value="common-input">Common input (Chapter 4)</option><option value="all">All including supplement</option><option value="fitted">Fitted only</option></select></label><label>Section<select id="section"><option value="all">All sections</option><option value="4.1–4.2">4.1–4.2 Baseline / RGE</option><option value="4.3">4.3 Fermion scale</option><option value="4.4">4.4 Couplings</option><option value="4.5">4.5 Thresholds</option><option value="Supplement (not 4.1–4.5)">Supplement</option></select></label></div>
    <section id="viewer"><button class="close" onclick="closePreview()">Close preview</button><h2 id="viewertitle"></h2><div id="viewerbody"></div></section><section id="cards">''' + ''.join(cards) + '''</section><footer>This catalogue does not regenerate calculations or edit original outputs. Open from its project location to keep relative links valid.</footer></main>
    <script>const assets=''' + jsdata + ''';function filter(){const s=document.getElementById('scope').value,t=document.getElementById('section').value;document.querySelectorAll('.card').forEach(e=>{const sc=e.dataset.scope; e.style.display=( (s==='all'||(s==='fitted'?sc!=='common-input':sc==='common-input'))&&(t==='all'||e.dataset.section===t))?'':'none'})}document.querySelectorAll('select').forEach(e=>e.addEventListener('change',filter));function preview(i){const a=assets[i];if(!a.exists)return;const el=document.getElementById('viewer'),body=document.getElementById('viewerbody');body.replaceChildren();document.getElementById('viewertitle').textContent=a.title;const node=document.createElement(a.kind==='interactive'?'iframe':'img');node.src=a.url;if(a.kind==='interactive')node.title=a.title;body.appendChild(node);el.style.display='block';el.scrollIntoView({behavior:'smooth'})}function closePreview(){document.getElementById('viewer').style.display='none';document.getElementById('viewerbody').replaceChildren()}filter();</script></body></html>'''


def build(root: Path, strict: bool = False) -> int:
    output_dir = root / OUTPUT
    output_dir.mkdir(parents=True, exist_ok=True)
    items = []
    for candidate in candidates():
        absolute = root / candidate.path
        exists = absolute.is_file() and absolute.stat().st_size > 0
        items.append(dict(section=candidate.section, title=candidate.title, path=candidate.path.as_posix(),
                          absolute=absolute, kind=candidate.kind, purpose=candidate.purpose,
                          scope=candidate.scope, recommendation=candidate.recommendation, exists=exists,
                          bytes=absolute.stat().st_size if exists else 0))
    page = output_dir / 'index.html'
    page.write_text(render(items, page), encoding='utf-8')
    with (output_dir / 'inventory.csv').open('w', encoding='utf-8', newline='') as stream:
        fields = ('section', 'title', 'path', 'kind', 'purpose', 'scope', 'recommendation', 'exists', 'bytes')
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for item in items:
            writer.writerow({k: item[k] for k in fields})
    lines = ['# Chapter 4 interactive and existing-figure index', '',
             'Generated from files already saved in this local checkout; no physics rerun or original output change.',
             '', 'Open [interactive/index.html](index.html) for the clickable gallery.', '',
             '| Section | Saved item | Provenance | Status | Proposed use |',
             '|---|---|---|---|---|']
    for item in items:
        status = 'verified local file' if item['exists'] else '**missing locally**'
        lines.append(f"| {item['section']} | `{item['path']}` | {item['scope']} | {status} | {item['recommendation']} |")
    lines.extend(['', '## Scientific restrictions', '',
                  '- Common-input studies and independently fitted benchmark studies must not be combined numerically.',
                  '- The overview scans include sixteen models; Figures 4.3–4.6 still require their six-model thesis versions.',
                  '- `c5_threshold_contributions.png` compares absolute matrix entries; no sign or interference can be inferred from bars.',
                  '- The common-input RGE dashboard can display a threshold-anchored visual continuation; it is not necessarily the true intermediate-EFT C5 coefficient.',
                  '- This gallery links to original files, and will not work as a standalone relocated HTML without the repository.', ''])
    (output_dir / 'INTERACTIVE_INDEX.md').write_text('\n'.join(lines), encoding='utf-8')
    found = sum(x['exists'] for x in items)
    print(f'Indexed {found}/{len(items)} existing files; missing files are explicitly labelled.')
    for item in items:
        if not item['exists']:
            print('MISSING:', item['path'])
    print('Gallery:', page)
    print('Index:', output_dir / 'INTERACTIVE_INDEX.md')
    return 1 if strict and found != len(items) else 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT, help='Local repository root')
    p.add_argument('--strict', action='store_true', help='Fail on missing listed assets')
    args = p.parse_args()
    return build(args.root.resolve(), args.strict)


if __name__ == '__main__':
    raise SystemExit(main())
