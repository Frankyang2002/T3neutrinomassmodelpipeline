"""Read-only presentation builder for the frozen T3 16 x 3 general-quartic campaign.

No RGE integration, fitting, matching, or experimental calculations are called.
The exact frozen CSV SHA-256 is mandatory. All physics source paths are read-only.
Figures depict saved endpoints and screening decisions, NOT new RGE trajectories
or UV-sensitivity ellipses. Historical analytical and dashboard links are labelled.

PowerShell from repository root:
    python -m scripts.BuildGeneralScalar48FinalReports --verify-only
    python -m scripts.BuildGeneralScalar48FinalReports
Dependencies for PDFs: matplotlib, reportlab. HTML/CSVs require only stdlib.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
from typing import Any

REPO = Path(__file__).resolve().parents[1]
FROZEN_REL = Path('output/studies/general_scalar_extended/bulk_48/bulk48_v1')
EXPECTED_SHA = 'B3B9FE0EF2F3730AE49504B51CDD7AED4D4232E6812CBBFE0D16D78CD23ED3AF'
GITHUB = 'https://github.com/Frankyang2002/T3neutrinomassmodelpipeline/blob/main/'
SCENARIOS = [
    'small_uv_frozen_5e-9',
    'large_uv_frozen_0p5',
    'small_uv_running_5e-9',
]
TITLES = {
    SCENARIOS[0]: 'Small topology quartic / frozen UV',
    SCENARIOS[1]: 'Large topology quartic / frozen UV',
    SCENARIOS[2]: 'Small topology quartic / running UV',
}
DIMS = {'A': (1, 3, 2), 'B': (2, 2, 1), 'C': (2, 2, 3),
        'D': (3, 1, 2), 'E': (3, 3, 2)}
# Verified against RGE/general/generated/general_rgbeta_quartic_registry.json
# on main, tree bfcf03408cd868a5a8f6cc4b2a2f5c67d36cd5fc.
EXTRA_QUARTICS = {
    'A_+0': 5, 'A_-2': 5, 'A_-4': 0,
    'B_+1': 1, 'B_-1': 7, 'B_-3': 1,
    'C_+1': 1, 'C_-1': 7, 'C_-3': 1,
    'D_+0': 5, 'D_+2': 0, 'D_-2': 5,
    'E_+0': 6, 'E_+2': 0, 'E_-2': 6, 'E_-4': 0,
}
ACCEPTED = {'CandidateConstrainedAccepted', 'ExistingConditionalAccepted'}
POOR = 'CandidateConstrainedPoorFit'
SCREEN = 'NoMassSafeNeutrinoFit'
METRICS = ['delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
           'sin2_theta12', 'sin2_theta13', 'sin2_theta23']


def model_sort(name: str) -> tuple[str, int]:
    a, b = name.split('_')
    return (a, int(b))


def model_slug(name: str) -> str:
    a, b = name.split('_')
    return f'{a}_{"p" if b.startswith("+") else "m"}{abs(int(b))}'


def model_report_folder(name: str) -> str:
    a, b = name.split('_')
    return f'T3-{a}-{ "p" if b.startswith("+") else "m"}{abs(int(b))}'


def relative_link(root: Path, report: Path, relative_file: str) -> str:
    file = root / relative_file
    if file.is_file():
        return Path(os.path.relpath(file, report.parent)).as_posix()
    return GITHUB + relative_file.replace(' ', '%20')


def fnum(x: Any, precision: int = 5) -> str:
    if x is None or x == '':
        return 'unavailable'
    try:
        v = float(x)
        return f'{v:.{precision}g}' if math.isfinite(v) else 'unavailable'
    except (TypeError, ValueError):
        return 'unavailable'


def safe_float(x: Any) -> float | None:
    if x is None or x == '':
        return None
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def read_rows(root: Path) -> tuple[list[dict], str]:
    source = root / FROZEN_REL / 'campaign_frozen.csv'
    if not source.is_file():
        raise FileNotFoundError(f'Frozen input missing: {source}')
    digest = hashlib.sha256(source.read_bytes()).hexdigest().upper()
    if digest != EXPECTED_SHA:
        raise ValueError(f'Frozen CSV SHA-256 mismatch: {digest}; expected {EXPECTED_SHA}. '
                         'Do not regenerate or silently re-encode the frozen source.')
    with source.open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    models = sorted({r['model'] for r in rows}, key=model_sort)
    if len(rows) != 48 or len(models) != 16 or set(models) != set(EXTRA_QUARTICS):
        raise ValueError(f'Expected 48 rows / 16 known models, found {len(rows)} / {len(models)}')
    pairs = [(r['model'], r['study']) for r in rows]
    if len(set(pairs)) != 48 or {r['study'] for r in rows} != set(SCENARIOS):
        raise ValueError('Frozen dataset must have one row per model/scenario')
    counts = Counter(r['status'] for r in rows)
    if (sum(counts[s] for s in ACCEPTED), counts[POOR], counts[SCREEN]) != (20, 6, 22):
        raise ValueError(f'Unexpected frozen status distribution: {counts}')
    if sum(r['boundary_relaxed'].lower() == 'true' for r in rows) != 2:
        raise ValueError('Expected exactly two modified-boundary cases')
    for r in rows:
        r['chi2'] = safe_float(r['chi2'])
        for key in METRICS + ['rhoF', 'lambdaT3_input']:
            r[key] = safe_float(r[key])
        if r['scalar_masses_muS_tev']:
            masses = ast.literal_eval(r['scalar_masses_muS_tev'])
            if not isinstance(masses, list) or len(masses) != 2:
                raise ValueError('Expected two scalar mass entries')
            r['scalar_masses_muS_tev'] = [float(m) for m in masses]
        else:
            r['scalar_masses_muS_tev'] = None
        r['boundary_relaxed'] = r['boundary_relaxed'].lower() == 'true'
        r['available_observables'] = all(r[x] is not None for x in METRICS)
        r['class'] = r['model'].split('_')[0]
        r['dimensions'] = DIMS[r['class']]
        r['alpha'] = int(r['model'].split('_')[1])
        r['extra_quartic_count_registry'] = EXTRA_QUARTICS[r['model']]
    return sorted(rows, key=lambda x: (model_sort(x['model']), SCENARIOS.index(x['study']))), digest


def read_registry(root: Path) -> tuple[str, dict[str, Any] | None]:
    path = root / 'RGE/general/generated/general_rgbeta_quartic_registry.json'
    if not path.is_file():
        return 'not available in current local source', None
    payload = json.loads(path.read_text(encoding='utf-8'))
    models = payload['models']
    for name, count in EXTRA_QUARTICS.items():
        if int(models[name]['total_additional_quartics']) != count or tuple(models[name]['dimensions']) != DIMS[name[0]]:
            raise ValueError(f'Registry changed unexpectedly for {name}: inspect before reporting')
    return str(payload.get('status', 'unknown')), payload


def extract_saved_yukawas(root: Path, row: dict) -> tuple[list[dict], str]:
    """Extract explicit saved y1/y2 entries only. Never use guessed/historical seeds."""
    folder = root / FROZEN_REL / row['study'] / model_slug(row['model'])
    fit = folder / 'fit.json'
    cfg = folder / 'fitted_config.json'
    values: list[dict] = []
    source = ''
    if cfg.is_file():
        d = json.loads(cfg.read_text(encoding='utf-8'))
        ordinary = d.get('base_state', {}).get('ordinary', {})
        try:
            for tensor in ('y1', 'y2'):
                matrix = ordinary[tensor]
                real = matrix['real']
                imag = matrix.get('imag')
                if len(real) != 3 or any(len(x) != 3 for x in real):
                    raise ValueError('not 3x3')
                for i in range(3):
                    for j in range(3):
                        values.append({'tensor': tensor, 'i': i, 'j': j,
                                       'real': float(real[i][j]),
                                       'imag': float(imag[i][j]) if imag is not None else 0.0})
            source = cfg.relative_to(root).as_posix()
        except (KeyError, IndexError, TypeError, ValueError):
            values = []
    if not values and fit.is_file():
        d = json.loads(fit.read_text(encoding='utf-8'))
        params = (d.get('best_mass_safe') or {}).get('parameters', {})
        try:
            for tensor in ('y1', 'y2'):
                for i in range(3):
                    for j in range(3):
                        values.append({'tensor': tensor, 'i': i, 'j': j,
                                       'real': float(params[f'{tensor}_{i}{j}']), 'imag': 0.0})
            source = fit.relative_to(root).as_posix() + ' (real parameters)'
        except (KeyError, TypeError, ValueError):
            values = []
    return values, source


def make_matrix_figure(rows: list[dict], dest: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    models = sorted({r['model'] for r in rows}, key=model_sort)
    by = {(r['model'], r['study']): r for r in rows}
    values = []
    labels = []
    for name in models:
        vrow = []
        lrow = []
        for scenario in SCENARIOS:
            r = by[name, scenario]
            if r['boundary_relaxed']:
                vrow.append(3); lrow.append('A*')
            elif r['status'] in ACCEPTED:
                vrow.append(2); lrow.append('A')
            elif r['status'] == POOR:
                vrow.append(1); lrow.append('P')
            else:
                vrow.append(0); lrow.append('X')
        values.append(vrow); labels.append(lrow)
    fig, ax = plt.subplots(figsize=(9, 8))
    ax.imshow(values, vmin=-0.5, vmax=3.5, aspect='auto', cmap=ListedColormap(['#e9e9e9', '#eed9b2', '#a6d5c4', '#82a9df']))
    ax.set_xticks(range(3), ['Small / UV frozen', 'Large / UV frozen', 'Small / UV running'], fontsize=10)
    ax.set_yticks(range(16), [n.replace('_', '(') + ')' for n in models], fontsize=10)
    for i in range(16):
        for j in range(3):
            ax.text(j, i, labels[i][j], ha='center', va='center', weight='bold', fontsize=10)
    ax.set_xticks([-.5, .5, 1.5, 2.5], minor=True)
    ax.set_yticks([i-.5 for i in range(17)], minor=True)
    ax.grid(which='minor', color='white', linewidth=1.8)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.set_title('Frozen T3 general-quartic benchmark outcomes', fontsize=14, pad=17)
    fig.text(.13, .018, 'A = conditional accepted; A* = modified boundary; P = poor fit; X = no mass-safe fit.', fontsize=9)
    fig.text(.13, .001, 'Saved endpoint classification only. No new RGE trajectories or sensitivity scans.', fontsize=8)
    fig.tight_layout(rect=(0, .043, 1, 1))
    fig.savefig(dest.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(dest.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_chi2_figure(rows: list[dict], dest: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected = [r for r in rows if r['status'] in ACCEPTED and r['chi2'] is not None and r['chi2'] > 0]
    selected.sort(key=lambda r: (SCENARIOS.index(r['study']), model_sort(r['model'])))
    fig, ax = plt.subplots(figsize=(12, 6))
    shades = {SCENARIOS[0]: '#397d8e', SCENARIOS[1]: '#6c6b93', SCENARIOS[2]: '#ba7836'}
    for i, r in enumerate(selected):
        ax.scatter(i, math.log10(r['chi2']), c=shades[r['study']], s=65,
                   marker='D' if r['boundary_relaxed'] else 'o', zorder=3)
    ax.axhline(math.log10(5.86), linestyle='--', color='#777777', linewidth=1.3, label='chi2 = 5.86 screen')
    ax.set_xticks(range(len(selected)), [r['model'] for r in selected], rotation=70, fontsize=8)
    ax.set_ylabel('log10(chi2)')
    ax.set_title('20 accepted conditional fits — saved 48-study endpoints')
    ax.grid(axis='y', alpha=.2)
    fig.text(.12, .018, 'Symbols with diamonds have modified topology-coupling UV boundaries; other points retain original boundaries.', fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, 1))
    fig.savefig(dest.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(dest.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_observable_fig(rows: list[dict], dest: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    available = [r for r in rows if r['status'] in ACCEPTED and r['available_observables']]
    # Five independently scaled axes: no covariance/sensitivity ellipse is implied.
    fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True)
    keys = METRICS
    names = ['Delta m21^2 [eV^2]', 'Delta m3l^2 [eV^2]', 'sin^2 theta12', 'sin^2 theta13', 'sin^2 theta23']
    models = sorted({r['model'] for r in rows}, key=model_sort)
    dx = {s: (i-1)*.23 for i, s in enumerate(SCENARIOS)}
    colors = {SCENARIOS[0]: '#397d8e', SCENARIOS[1]: '#6c6b93', SCENARIOS[2]: '#ba7836'}
    for ax, key, lab in zip(axes, keys, names):
        for s in SCENARIOS:
            for r in available:
                if r['study'] != s: continue
                x = models.index(r['model']) + dx[s]
                ax.scatter(x, r[key], s=28, color=colors[s], marker='D' if r['boundary_relaxed'] else 'o')
        ax.set_ylabel(lab, fontsize=9)
        ax.grid(axis='y', alpha=.2)
    axes[0].set_title('Saved oscillation observables for accepted conditional fits')
    axes[-1].set_xticks(range(16), models, rotation=45, ha='right')
    fig.text(.1, .016, 'Points are 100 GeV endpoint predictions (not RGE trajectories). No covariance or model sensitivity contours shown.', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, 1))
    fig.savefig(dest.with_suffix('.pdf'), bbox_inches='tight')
    fig.savefig(dest.with_suffix('.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def csv_write(path: Path, header: list[str], records: list[dict]) -> None:
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=header, extrasaction='ignore', lineterminator='\r\n')
        w.writeheader()
        for rec in records:
            w.writerow({k: '' if rec.get(k) is None else rec.get(k) for k in header})


def make_pdf(rows: list[dict], root: Path, out: Path, registry_status: str, yukawa_sources: dict[str, str]) -> None:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle, Image, KeepTogether
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='AtlasTitle', parent=styles['Title'], fontSize=19, leading=22, spaceAfter=12))
    styles.add(ParagraphStyle(name='AtlasSub', parent=styles['Heading2'], spaceBefore=12, spaceAfter=5))
    styles.add(ParagraphStyle(name='AtlasSmall', parent=styles['Normal'], fontSize=8.2, leading=11))
    styles.add(ParagraphStyle(name='AtlasTiny', parent=styles['Normal'], fontSize=7, leading=9))
    document = SimpleDocTemplate(str(out), pagesize=A4, rightMargin=19*mm, leftMargin=19*mm,
                                 topMargin=17*mm, bottomMargin=17*mm, title='Frozen T3 general-quartic benchmarks')
    story: list[Any] = []
    add = lambda text: story.append(Paragraph(html.escape(text), styles['AtlasSmall']))
    story.append(Paragraph('T3 general-quartic benchmark atlas', styles['AtlasTitle']))
    add('Frozen 16 models x 3 scenarios | read-only presentation of campaign_frozen.csv')
    story.append(Spacer(1, 7*mm))
    add('SOURCE INTEGRITY: SHA-256 ' + EXPECTED_SHA)
    add('Accepted conditional: 20 (18 original-boundary + 2 modified-boundary); poor fits: 6; screened out: 22.')
    add('All five oscillation observables and scalar masses below come directly from the frozen CSV. Never interpret a screened result as an evaluated oscillation failure.')
    story.append(Paragraph('Scientific status', styles['AtlasSub']))
    for item in [
        'Candidate/conditional general-quartic comparison, not a full general-potential or fully sequential matched EFT.',
        'Additional renormalisable UV/EFT1 scalar quartics are used where exported. Holomorphic quadratic extensions, full general-scalar Wilson matching, and physical Matchete/RGBeta normalisation certification are absent.',
        'Four extreme hypercharge cases use legacy common-scale C5 matching sources. B(-1)/C(-1) are imported, independently validated historical conditional fits; their running-UV cases used modified topology-coupling boundaries.',
        'Earlier Reports/output analytical equations, RGE dashboards and sensitivity ovals are retained as historical references, not recalculated general-quartic figures.',
        'No complete frozen 48-study RGE trajectories or local sensitivity covariance/scan outputs are in the campaign table. They are not invented in this atlas.',
    ]:
        story.append(Paragraph('&#8226; ' + html.escape(item), styles['AtlasSmall']))
    story.append(Paragraph('Source and scope', styles['AtlasSub']))
    add('Dataset: ' + (FROZEN_REL / 'campaign_frozen.csv').as_posix())
    add('Additional quartic registry status: ' + registry_status)
    add('The exact three scenario names are reproduced in the following pages. E-class lambdaT3 inputs are half-normalised (2.5e-9 / 0.25) as frozen.')
    story.append(PageBreak())
    story.append(Paragraph('Model representations and analytical references', styles['Heading1']))
    add('Representation tuple is (SU(2) dimension of S1, S2, F); alpha is the model hypercharge parameter. These are not complete field quantum-number tables.')
    models = sorted({r['model'] for r in rows}, key=model_sort)
    tbl = [['Model', '(dS1,dS2,dF)', 'alpha', 'Registered extra quartics', 'Analytical source folder']]
    for name in models:
        tbl.append([name, str(DIMS[name[0]]), str(int(name.split('_')[1])), str(EXTRA_QUARTICS[name]),
                    'Reports/output/models/' + model_report_folder(name)])
    make_table(tbl, story, widths=[37, 94, 43, 90, 208], fs=7.0)
    story.append(Paragraph('Analytical inventory', styles['AtlasSub']))
    add('Each model-specific source folder contains UV and EFT1 Lagrangian, RGE and group-factor LaTeX, plus EFT2/Weinberg-stage files. Existing master PDFs/TeX are in Reports/output/analytical/. Generated additional-quartic symbolic sources are in RGE/general/generated/. These analytical files have not been relabelled as 48-study equations.')
    story.append(PageBreak())
    story.append(Paragraph('16 x 3 screening matrix', styles['Heading1']))
    add('A = accepted original boundary; A* = accepted modified boundary; P = poor fit; X = no mass-safe fit.')
    by = {(r['model'], r['study']): r for r in rows}
    tbl = [['Model', 'Small / frozen', 'Large / frozen', 'Small / running']]
    for name in models:
        statuses = []
        for scenario in SCENARIOS:
            r = by[name, scenario]
            statuses.append('A*' if r['boundary_relaxed'] else 'A' if r['status'] in ACCEPTED else 'P' if r['status'] == POOR else 'X')
        tbl.append([name] + statuses)
    make_table(tbl, story, widths=[85, 130, 130, 130], fs=8)
    story.append(Paragraph('Fit quality by scenario', styles['AtlasSub']))
    for s in SCENARIOS:
        rs = [r for r in rows if r['study'] == s]
        add(TITLES[s] + ': ' + str(sum(r['status'] in ACCEPTED for r in rs)) + ' accepted; ' +
            str(sum(r['status'] == POOR for r in rs)) + ' poor; ' +
            str(sum(r['status'] == SCREEN for r in rs)) + ' screened out.')
    for s in SCENARIOS:
        story.append(PageBreak())
        story.append(Paragraph(TITLES[s], styles['Heading1']))
        add('Saved low-energy endpoints; masses at muS (TeV). Missing entries are explicitly unavailable. B/C existing imported fits are not new bulk optimisations.')
        tbl = [['Model', 'status', 'chi2', 'mS1', 'mS2', 'lambdaT3']]
        for r in rows:
            if r['study'] == s:
                masses = r['scalar_masses_muS_tev']
                tbl.append([r['model'], 'Accepted*' if r['boundary_relaxed'] else 'Accepted' if r['status'] in ACCEPTED else 'Poor' if r['status'] == POOR else 'Screened',
                            fnum(r['chi2'], 6), fnum(masses[0], 5) if masses else '—',
                            fnum(masses[1], 5) if masses else '—', fnum(r['lambdaT3_input'], 4)])
        make_table(tbl, story, widths=[59, 85, 102, 80, 80, 80], fs=7)
        story.append(Paragraph('Matching-source accounting', styles['AtlasSub']))
        src = Counter(r['c5_source_kind'] for r in rows if r['study'] == s)
        add('; '.join(f'{k}: {v}' for k, v in sorted(src.items())))
        story.append(PageBreak())
        story.append(Paragraph('Neutrino observables — '+TITLES[s], styles['Heading1']))
        add('Numerical predictions are recorded only where the frozen campaign reached the neutrino evaluator. Blank results were screened out, not assigned artificial residuals.')
        tbl = [['Model', 'Delta m21^2', 'Delta m3l^2', 'sin2 th12', 'sin2 th13', 'sin2 th23']]
        for r in rows:
            if r['study'] == s:
                tbl.append([r['model']] + [fnum(r[k], 5) if r[k] is not None else '—' for k in METRICS])
        make_table(tbl, story, widths=[59, 85, 90, 86, 86, 80], fs=7)
    story.append(PageBreak())
    story.append(Paragraph('Saved Yukawa information and missing outputs', styles['Heading1']))
    add('Matrix entries are exported separately in tables/saved_yukawas.csv when an explicit local fitted_config.json or fit.json contains all 18 entries. The code does not infer Yukawas from low-energy observables, or substitute generic defaults for missing fit records.')
    add('Yukawa sources located for this reporting environment: ' + str(len(yukawa_sources)) + '/48 cases. Unavailable cases remain blank in tables/yukawa_availability.csv.')
    story.append(Paragraph('RGE and oval figures', styles['AtlasSub']))
    add('The saved CSV contains only endpoint diagnostics. It is insufficient to construct genuine UV/EFT1/final-C5 running functions of scale or the UV-parameter response required for sensitivity ovals. The existing historical HTML visualisations remain accessible via the interactive index with a historical label. No numerical integration was performed during this reporting build.')
    story.append(Paragraph('Reproducibility', styles['AtlasSub']))
    add('The source and missing-data manifest is in provenance.json. This PDF and its plots are generated by scripts/BuildGeneralScalar48FinalReports.py using only existing frozen data. No fitting, matching or running is executed.')
    def footer(c, d):
        c.saveState(); c.setFont('Helvetica', 8)
        c.setFillColor(colors.HexColor('#666666'))
        c.drawString(19*mm, 11*mm, 'T3 / frozen general-quartic / conditional, not fully matched')
        c.drawRightString(A4[0]-19*mm, 11*mm, f'{d.page}')
        c.restoreState()
    document.build(story, onFirstPage=footer, onLaterPages=footer)


def make_table(table: list[list], story: list, widths: list[int], fs: float) -> None:
    from reportlab.lib import colors
    from reportlab.platypus import Table, TableStyle
    obj = Table(table, colWidths=widths, repeatRows=1, hAlign='LEFT')
    obj.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#e6edf5')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor('#182638')),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f5f7fa')]),
        ('LINEBELOW', (0, 0), (-1, 0), .6, colors.HexColor('#8da0b8')),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('LEFTPADDING', (0, 0), (-1, -1), 5), ('RIGHTPADDING', (0, 0), (-1, -1), 4),
        ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ('FONTSIZE', (0, 0), (-1, -1), fs),
    ]))
    story.append(obj)


def fmt_html(x: Any) -> str:
    if isinstance(x, list): return ', '.join(fnum(t, 6) for t in x)
    if isinstance(x, (int, float)): return fnum(x, 7)
    return str(x) if x else 'unavailable'


def html_head(title: str) -> str:
    return '''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>''' + html.escape(title) + '''</title><style>
:root{--ink:#16243a;--muted:#576779;--line:#d8e0e9;--surface:#fff;--bg:#f4f7fb;--blue:#254f86}
*{box-sizing:border-box}body{margin:0;background:var(--bg);font:15px/1.55 system-ui,Segoe UI,sans-serif;color:var(--ink)}
main{max-width:1260px;margin:auto;padding:32px 22px 70px}h1{font-size:2.1rem;line-height:1.16;margin:0 0 10px}h2{margin:32px 0 12px;font-size:1.3rem}p{max-width:100ch}a{color:#1c5592}small,.subtle{color:var(--muted)}
.panel{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px;margin:18px 0;overflow-x:auto}
.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px}.metric{background:#fff;border:1px solid var(--line);border-radius:10px;padding:14px}.metric b{display:block;font-size:27px}
select,button{font:inherit;padding:7px 10px;border:1px solid #9eb0c4;border-radius:8px;background:#fff;color:var(--ink)}button{cursor:pointer}
.chip{padding:2px 7px;border-radius:5px;background:#eaf0f7;font-size:12px;font-weight:600}table{width:100%;border-collapse:collapse;font-size:13px}th,td{text-align:left;padding:7px;border-bottom:1px solid #e4e9f0}th{background:#f2f5f9}td.num{font-variant-numeric:tabular-nums}tr:hover td{background:#f3f8fd}
.legend span{margin-right:16px;white-space:nowrap}svg{max-width:100%;height:auto}details{margin-top:14px}summary{cursor:pointer;font-weight:600}
footer{color:#68798c;font-size:12px;margin-top:30px}a:focus,button:focus,select:focus{outline:2px solid #4d80ad}
</style></head><body><main>'''


def make_model_pages(rows: list[dict], root: Path, out: Path, availability: dict[str, str]) -> None:
    folder = out / 'models'
    folder.mkdir(parents=True, exist_ok=True)
    for model in sorted({r['model'] for r in rows}, key=model_sort):
        rp = folder / (model_slug(model) + '.html')
        these = [r for r in rows if r['model'] == model]
        part = [html_head('T3 '+model+' — saved benchmark report'), '<a href="../index.html">Back to 48-study dashboard</a>',
                '<h1>T3 ' + html.escape(model) + '</h1>',
                f'<p>SU(2) dimensions (S1,S2,F): {DIMS[model[0]]}; alpha = {int(model.split("_")[1]):+d}. '
                f'Registered additional quartics: {EXTRA_QUARTICS[model]} (candidate-basis count, not certification).</p>',
                '<h2>Frozen endpoint predictions</h2><div class="panel"><table><thead><tr><th>Scenario</th><th>Status</th><th>chi2</th><th>mS1, mS2 [TeV]</th><th>Matching source</th></tr></thead><tbody>']
        for r in these:
            part.append('<tr>' + ''.join('<td>'+html.escape(str(v))+'</td>' for v in [TITLES[r['study']], r['status']+(' (modified boundary)' if r['boundary_relaxed'] else ''),
                         fnum(r['chi2'], 7), fmt_html(r['scalar_masses_muS_tev']),r['c5_source_kind']])+'</tr>')
        part.append('</tbody></table></div><h2>Five neutrino observables</h2><div class="panel"><table><thead><tr><th>Scenario</th>'+
                    ''.join(f'<th>{html.escape(k)}</th>' for k in METRICS)+'</tr></thead><tbody>')
        for r in these:
            part.append('<tr><td>'+html.escape(TITLES[r['study']])+'</td>'+ ''.join('<td>'+fnum(r[k], 7)+'</td>' for k in METRICS)+'</tr>')
        part.append('</tbody></table></div><h2>Saved Yukawa availability</h2><div class="panel">')
        for r in these:
            key = r['model']+'|'+r['study']
            part.append('<p>'+html.escape(TITLES[r['study']])+': '+html.escape(availability.get(key) or 'unavailable in 48-study local fit/config files')+'</p>')
        part.append('</div><h2>Analytical material (historical source)</h2><div class="panel">')
        base = 'Reports/output/models/'+model_report_folder(model)
        for sector, names in [('Lagrangian', ['UV.tex','EFT_1_after_F.tex','EFT_2_after_S1_S2.tex','C5.tex']),
                              ('RGE',['UV.tex','EFT_1_after_F.tex','EFT_2_after_S1_S2.tex']),
                              ('GroupFactors',['GF_UV.tex','GF_EFT_1_after_F.tex','GF_EFT_2_after_S1_S2.tex'])]:
            part.append(f'<h3>{sector}</h3><p>')
            for fname in names:
                rel=f'{base}/{sector}/{fname}'
                url=relative_link(root,rp,rel)
                part.append(f'<a href="{html.escape(url,quote=True)}">{html.escape(fname)}</a> &nbsp; ')
            part.append('</p>')
        alpha=int(model.split('_')[1]); expr=f'RGE/general/generated/{model[0]}_{alpha:+d}_extra_quartics.wl'
        part.append('<p>New extra-quartic symbolic source: <a href="'+html.escape(relative_link(root,rp,expr),quote=True)+'">'
                    +html.escape(expr)+'</a>. Candidate basis; not a certified general-scalar EFT RGE report.</p>')
        part.append('</div><p><strong>Not available:</strong> 48-study UV/EFT1/Weinberg trajectory arrays and newly evaluated oval sensitivity responses. Earlier dashboards are historical only.</p>')
        part.append('<footer>Frozen campaign SHA-256 '+EXPECTED_SHA+'<br>Presentation only; no physics recalculated.</footer></main></body></html>')
        rp.write_text(''.join(part),encoding='utf-8')


def make_index(rows: list[dict], root: Path, out: Path) -> None:
    target=out/'index.html'
    by = sorted(rows,key=lambda r:(model_sort(r['model']),SCENARIOS.index(r['study'])))
    # Minimal rows with exact stored floats; no simulated points.
    payload=[{k:v for k,v in r.items() if k in ('model','study','status','chi2','rhoF',
             'scalar_masses_muS_tev','boundary_relaxed','lambdaT3_input','c5_source_kind','reason',
             'delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')}
             for r in by]
    jsdata=json.dumps(payload,ensure_ascii=True,separators=(',',':')).replace('</','<\\/')
    h=[html_head('T3 frozen general-quartic 48-study dashboard'),
       '<p class="subtle">FROZEN GENERAL-QUARTIC / CONDITIONAL REPORT / 16 MODELS × 3 BENCHMARKS</p>',
       '<h1>Benchmark outcomes and neutrino endpoints</h1>',
       '<p>This is a presentation of saved endpoint data, not a new RGE integration or a complete general-scalar EFT. '
       'Older analytical equations and dashboards are linked separately and are <strong>not</strong> recalculated with the expanded quartic terms.</p>',
       '<div class="summary"><div class="metric"><small>Accepted conditional</small><b>20 / 48</b><small>18 original; 2 modified boundary</small></div>'
       '<div class="metric"><small>Poor fits</small><b>6</b><small>Oscillation residual evaluated</small></div>'
       '<div class="metric"><small>Screened out</small><b>22</b><small>No mass-safe candidate</small></div></div>',
       '<div class="panel"><strong>Physics scope:</strong> registered general-quartic UV/EFT1 running where implemented; '
       'no holomorphic quadratic completion, no complete sequential general-potential Wilson matching, '
       'and no certified Matchete/RGBeta coupling map. Four extreme-hypercharge cases use legacy C5 sources. '
       'E-class frozen topology-coupling inputs are half-normalised.</div>',
       '<h2>Interactive 16 × 3 status matrix</h2><div class="panel">',
       '<p class="legend"><span>A: accepted</span><span>A*: accepted, modified boundary</span><span>P: poor fit</span><span>X: screened out</span></p>',
       '<div id="heatmap"></div><p class="subtle">Select a cell to view its exact frozen observables, masses and matching source.</p></div>',
       '<h2>Selected result</h2><div class="panel"><label>Model <select id="model"></select></label> &nbsp; '
       '<label>Scenario <select id="scenario"></select></label><div id="detail"></div></div>',
       '<h2>Endpoint observable comparison</h2><div class="panel"><label>Observable <select id="observable">'
       '<option value="sin2_theta12">sin² θ12</option><option value="sin2_theta13">sin² θ13</option><option value="sin2_theta23">sin² θ23</option>'
       '<option value="delta_m21_sq_ev2">Δm²21</option><option value="delta_m3l_sq_ev2">Δm²3l</option></select></label>'
       '<div id="scatter"></div><p class="subtle">Saved accepted endpoint values only; not RG trajectories or sensitivity contours.</p></div>',
       '<h2>Model-specific reports and historical analytical equations</h2><div class="panel"><div id="modellinks"></div></div>',
       '<h2>Historical RGE dashboards and sensitivity ovals</h2><div class="panel"><p>These visualisations come from earlier benchmarks and have not been updated to the frozen 48-study general-quartic equations.</p>']
    old=[('Historical fermion-start / 5e-9 / optimised RGE', 'Reports/output/display_fermion_start/5e-9/reused/optimised/interactive_rge.html'),
         ('Historical fermion-start / 5e-9 / optimised ovals', 'Reports/output/display_fermion_start/5e-9/reused/optimised/sensitivity_ovals.html'),
         ('Historical lambda 5e-9 / optimised RGE', 'Reports/output/display_lambda_5e9/optimised_rge.html'),
         ('Historical thesis-ready figures','Reports/output/thesis_assets/figures/fig08_t3b_weinberg_mass_splitting_running.pdf')]
    for label, path in old:
        h.append('<p><a href="'+html.escape(relative_link(root,target,path),quote=True)+'">'+html.escape(label)+'</a></p>')
    h.append('</div><h2>New genuine RGE / sensitivity data not in the frozen CSV</h2>'
             '<div class="panel"><p>Numerical scale-by-scale UV, EFT1, and final Weinberg trajectories and parameter-variation/covariance samples '
             'are unavailable in this frozen 48-study CSV and its bulk fit records. No replacement curves have been synthesised. '
             'To update these graphs, a separate <em>read-only trajectory export from saved fitted inputs</em> would be required; '
             'that step is deliberately not performed by this presentation builder.</p></div>')
    h.append('<h2>Downloads</h2><div class="panel"><p>'
             '<a href="benchmark_atlas.pdf">Benchmark atlas (PDF)</a> · '
             '<a href="figures/status_matrix.pdf">Status matrix (PDF)</a> · '
             '<a href="figures/chi2_accepted.pdf">Fit quality (PDF)</a> · '
             '<a href="figures/accepted_observables.pdf">Oscillation observables (PDF)</a> · '
             '<a href="tables/benchmark_comparison.csv">All 48 results (CSV)</a> · '
             '<a href="tables/model_inventory.csv">Representations (CSV)</a> · '
             '<a href="tables/yukawa_availability.csv">Yukawa provenance (CSV)</a> · '
             '<a href="provenance.json">Reproducibility (JSON)</a></p></div>')
    h.append('<footer>Frozen input SHA-256 '+EXPECTED_SHA+'. Presentation only. Conditional, incomplete general scalar sector.</footer>')
    h.append('<script>const D='+jsdata+';const S='+json.dumps(SCENARIOS)+';const labels='+json.dumps(TITLES)+';'
             'const M='+json.dumps(sorted({r['model'] for r in rows},key=model_sort))+';\n')
    h.append(r'''const $=id=>document.getElementById(id); const mk=(tag,attrs={},text='')=>{let e=document.createElementNS('http://www.w3.org/2000/svg',tag); for(const k in attrs)e.setAttribute(k,attrs[k]);e.textContent=text;return e};
const lookup=(m,s)=>D.find(r=>r.model===m&&r.study===s);
const status=r=>r.boundary_relaxed?'A*':r.status.includes('Accepted')?'A':r.status.includes('PoorFit')?'P':'X';
const fill=r=>status(r)==='A*'?'#82a9df':status(r)==='A'?'#a6d5c4':status(r)==='P'?'#eed9b2':'#e9e9e9';
const num=(v,p=6)=>v===null||v===undefined?'unavailable':Number(v).toPrecision(p);
for (let m of M){let o=document.createElement('option');o.value=m;o.textContent=m;$('model').appendChild(o)}
for(let s of S){let o=document.createElement('option');o.value=s;o.textContent=labels[s];$('scenario').appendChild(o)}
$('scenario').value=S[0];$('model').value=M[0];
const H=770,W=950,step=38;
let svg=mk('svg',{viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':'16 by 3 frozen benchmark status matrix'});
svg.append(mk('title',{},'Select a model and scenario cell'));
S.forEach((s,j)=>svg.append(mk('text',{x:180+j*240,y:24,'font-size':14,'font-weight':650},labels[s].replace(' topology quartic',''))));
M.forEach((m,i)=>{
 const y=44+i*43;svg.append(mk('text',{x:70,y:y+27,'font-size':15,'text-anchor':'end'},m));
 S.forEach((s,j)=>{const r=lookup(m,s), x=160+j*240;
   let g=mk('g',{'tabindex':'0',role:'button','aria-label':m+' '+labels[s]+' '+r.status,style:'cursor:pointer'});
   g.append(mk('rect',{x,y,width:200,height:38,rx:5,fill:fill(r),stroke:'#c0cbd8'}));
   g.append(mk('text',{x:x+100,y:y+26,'text-anchor':'middle','font-size':15,'font-weight':650},status(r)));
   g.append(mk('title',{},m+' — '+labels[s]+' — '+r.status+'; chi2='+num(r.chi2)));
   const select=()=>{$('model').value=m;$('scenario').value=s;render()};g.addEventListener('click',select);
   g.addEventListener('keydown',ev=>{if(ev.key==='Enter'||ev.key===' '){ev.preventDefault();select()}});svg.append(g)
 })
});$('heatmap').appendChild(svg);
const esc=t=>String(t).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function render(){let m=$('model').value,s=$('scenario').value,r=lookup(m,s);
 let fields=[['Status',r.status],['chi²',num(r.chi2,9)],['Scalar masses at muS (TeV)',r.scalar_masses_muS_tev?r.scalar_masses_muS_tev.map(v=>num(v,8)).join(', '):'unavailable'],
 ['Δm²21 (eV²)',num(r.delta_m21_sq_ev2,9)],['Δm²3l (eV²)',num(r.delta_m3l_sq_ev2,9)],
 ['sin² θ12',num(r.sin2_theta12,9)],['sin² θ13',num(r.sin2_theta13,9)],['sin² θ23',num(r.sin2_theta23,9)],
 ['rhoF',num(r.rhoF)],['topology coupling input',num(r.lambdaT3_input)],['C5 source',r.c5_source_kind],['modified boundary',r.boundary_relaxed?'YES':'no'],['note',r.reason||'—']];
 $('detail').innerHTML='<h3>'+esc(m)+' / '+esc(labels[s])+'</h3><table><tbody>'+fields.map(([k,v])=>'<tr><th>'+esc(k)+'</th><td class="num">'+esc(v)+'</td></tr>').join('')+'</tbody></table>'+
 '<p><a href="models/'+m.replace('_+','_p').replace('_-','_m').replace('p0','p0')+'.html">Model analytical-source and fit page</a></p>';
}
$('model').addEventListener('change',render);$('scenario').addEventListener('change',render);
const plot=()=>{let s=$('scenario').value,key=$('observable').value,rs=D.filter(r=>r.study===s&&r.status.includes('Accepted')&&r[key]!==null);
  const plotSVG=mk('svg',{viewBox:'0 0 960 310',role:'img','aria-label':'Saved endpoint observable values'}), vals=rs.map(r=>r[key]);
  let lo=Math.min(...vals),hi=Math.max(...vals);if(!Number.isFinite(lo)){ $('scatter').textContent='No saved accepted endpoints';return}
  const span=Math.max(hi-lo,Math.max(Math.abs(lo),1e-12)*0.002,1e-12);lo-=span*.2;hi+=span*.2;
  plotSVG.append(mk('line',{x1:75,y1:12,x2:75,y2:265,stroke:'#8194a9'}));plotSVG.append(mk('line',{x1:75,y1:265,x2:940,y2:265,stroke:'#8194a9'}));
  for(let z=0;z<=4;z++){let y=265-z*60;plotSVG.append(mk('text',{x:66,y:y+4,'text-anchor':'end','font-size':11},num(lo+(hi-lo)*z/4,4)));plotSVG.append(mk('line',{x1:75,y1:y,x2:940,y2:y,stroke:'#ecf0f4'}))}
  M.forEach((m,i)=>{let x=92+i*53;plotSVG.append(mk('text',{x,y:289,'text-anchor':'middle','font-size':11},m))});
  rs.forEach(r=>{let x=92+M.indexOf(r.model)*53, y=265-(r[key]-lo)/(hi-lo)*240;const c=mk('circle',{cx:x,cy:y,r:5.5,fill:r.boundary_relaxed?'#ba7836':'#397d8e'});
     c.append(mk('title',{},r.model+' : '+num(r[key],9)+'; chi2='+num(r.chi2,5)));plotSVG.append(c)});
  $('scatter').replaceChildren(plotSVG);
};$('observable').addEventListener('change',plot);$('scenario').addEventListener('change',plot);
$('modellinks').innerHTML=M.map(m=>'<a style="display:inline-block;min-width:90px;padding:6px" href="models/'+m.replace('_+','_p').replace('_-','_m')+'.html">'+esc(m)+'</a>').join('');
render();plot();
</script></main></body></html>''')
    target.write_text(''.join(h),encoding='utf-8')


def build(root: Path, out: Path, verify_only: bool) -> None:
    rows, digest = read_rows(root)
    registry_status, registry = read_registry(root)
    counts = Counter(r['status'] for r in rows)
    print('Frozen source SHA-256:', digest)
    print('48 rows / 16 models / 3 scenarios:', dict(counts))
    print('Quartic registry status:', registry_status)
    if verify_only:
        print('VERIFIED. No files written.')
        return
    # Output must not overlap the frozen input or source report directories.
    if out.resolve() == (root / FROZEN_REL).resolve() or (root / FROZEN_REL).resolve() in out.resolve().parents:
        raise ValueError('Do not write reporting outputs into the frozen source directory')
    out.mkdir(parents=True,exist_ok=True)
    (out/'figures').mkdir(exist_ok=True)
    (out/'tables').mkdir(exist_ok=True)
    (out/'models').mkdir(exist_ok=True)
    models = sorted({r['model'] for r in rows},key=model_sort)
    avail: dict[str,str] = {}
    saved: list[dict] = []
    avail_rows: list[dict] = []
    for row in rows:
        y, src = extract_saved_yukawas(root, row)
        key=row['model']+'|'+row['study']
        avail[key]=src
        avail_rows.append({'model':row['model'],'study':row['study'], 'available': bool(y), 'source':src or 'UNAVAILABLE'})
        for cell in y:
            saved.append({'model':row['model'],'study':row['study'], 'status':row['status'], **cell,'source':src})
    comparison=[]
    for r in rows:
        m=r['scalar_masses_muS_tev'] or [None,None]
        comparison.append({**r,'scalar_mass_1_muS_tev':m[0],'scalar_mass_2_muS_tev':m[1],
                           'representation':str(r['dimensions']),
                           'screen_decision':'accepted_modified_boundary' if r['boundary_relaxed'] else 'accepted_original_boundary' if r['status'] in ACCEPTED else 'poor_fit' if r['status']==POOR else 'screened_out'})
    head = ['model','representation','alpha','study','status','screen_decision','chi2',
            *METRICS,'scalar_mass_1_muS_tev','scalar_mass_2_muS_tev','rhoF',
            'lambdaT3_input','c5_source_kind','extra_quartic_count_registry','reason']
    csv_write(out/'tables/benchmark_comparison.csv',head,comparison)
    csv_write(out/'tables/model_inventory.csv', ['model','dS1','dS2','dF','alpha','extra_quartic_count_registry','analytical_reference'],
              [dict(model=m,dS1=DIMS[m[0]][0],dS2=DIMS[m[0]][1],dF=DIMS[m[0]][2],alpha=int(m.split('_')[1]),
                    extra_quartic_count_registry=EXTRA_QUARTICS[m],
                    analytical_reference='Reports/output/models/'+model_report_folder(m)) for m in models])
    csv_write(out/'tables/yukawa_availability.csv', ['model','study','available','source'], avail_rows)
    csv_write(out/'tables/saved_yukawas.csv', ['model','study','status','tensor','i','j','real','imag','source'], saved)
    shutil.copyfile(root/FROZEN_REL/'campaign_frozen.csv',out/'tables/frozen_source_sha_verified.csv')
    make_matrix_figure(rows,out/'figures/status_matrix')
    make_chi2_figure(rows,out/'figures/chi2_accepted')
    make_observable_fig(rows,out/'figures/accepted_observables')
    make_pdf(rows,root,out/'benchmark_atlas.pdf',registry_status,{x:y for x,y in avail.items() if y})
    make_model_pages(rows,root,out,avail)
    make_index(rows,root,out)
    new_paths = sorted(p.relative_to(out).as_posix() for p in out.rglob('*') if p.is_file() and p.name!='provenance.json')
    provenance={
        'format':'t3-frozen-general-quartic-48-presentation-v1',
        'input':'output/studies/general_scalar_extended/bulk_48/bulk48_v1/campaign_frozen.csv',
        'input_sha256':digest,
        'source_tree_sha_verified':'bfcf03408cd868a5a8f6cc4b2a2f5c67d36cd5fc',
        'read_only_physics':True,'fits_performed':0,'ode_integrations_performed':0,
        'matching_performed':False,'full_general_scalar_eft':False,
        'quartic_registry_status':registry_status,
        'model_count':16,'scenario_count':3,'row_count':48,'accepted':20,
        'accepted_original_boundary':18,'accepted_modified_boundary':2,
        'poor_fits':6,'screened_out':22,
        'saved_yukawa_cases_available_in_local_bulk_folder':sum(bool(x) for x in avail.values()),
        'missing_48study_full_scale_trajectories':True,
        'missing_48study_sensitivity_oval_scans':True,
        'legacy_figure_link_policy':'Historical links only; never labelled general-quartic recalculations',
        'limitations':[
            'Full general scalar holomorphic quadratic sector not implemented in campaign',
            'Complete model-dependent general-quartic sequential Wilson matching not implemented',
            'Matchete/RGBeta map not physically certified',
            'Four extreme-hypercharge model sources use legacy C5',
            'B(-1)/C(-1) original historical fits reused, not reoptimised in bulk',
            'No saved 48-study RGE scale arrays or oval response scans',
        ],
        'generated_files':new_paths,
    }
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    (out/'README.md').write_text(
        '# Frozen general-quartic 48-study reports\n\n'
        'Open `index.html` for the offline interactive dashboard.\n\n'
        '`benchmark_atlas.pdf` is the PDF report; `figures/` contains thesis-ready vector PDF figures, '
        '`tables/` contains machine-readable 48-case data and available saved Yukawa matrices, '
        'and `models/` links all 16 analytical references.\n\n'
        '**Scientific scope:** conditional general-quartic benchmark only; not full scalar potential or '
        'fully matched EFT. Previously generated RGE dashboards and oval plots are *historical*, '
        'not new general-quartic outputs. The frozen data do not contain the trajectories/sensitivity '
        'scan points required to produce their genuine 48-study replacements. No physics was rerun.\n\n'
        'Frozen CSV SHA-256: `' + EXPECTED_SHA + '`.\n',encoding='utf-8')
    print('WROTE',out,'with',len(new_paths)+2,'files; Yukawa cases available:',sum(bool(x) for x in avail.values()))


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=REPO,help='Local repository root')
    p.add_argument('--output',type=Path,help='Dedicated reporting destination; default Reports/output/general_scalar_48_final')
    p.add_argument('--verify-only',action='store_true',help='Validate frozen data and local registry; do not write')
    a=p.parse_args()
    root=a.root.resolve()
    out=(a.output or root/'Reports/output/general_scalar_48_final').resolve()
    build(root,out,a.verify_only)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
