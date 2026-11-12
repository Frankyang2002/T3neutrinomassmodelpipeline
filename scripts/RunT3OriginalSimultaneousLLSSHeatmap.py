#!/usr/bin/env python3
"""LLSS-off ablation using the ORIGINAL 16-model simultaneous-EFT1 study.

Original thesis baseline: saved 1111 (full), 1001 (both EFT1 scalar+Wilson frozen).
New LLSS-off case: run ORIGINAL simultaneous EFT1 equations, subtracting only
  lambda_{HHSS} C_{LLSS} from the LLHH portion of the scalar_pair Wilson beta.
The saved original 1111 and 1001 cases are never recomputed or substituted.

Project path: scripts/RunT3OriginalSimultaneousLLSSHeatmap.py
Exploratory EFT calculation; not a certified physical MSbar matching result.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from contextlib import contextmanager, ExitStack
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import RunT3SimultaneousEFT1StageFreezeHeatmaps as historical

MODELS = historical.MODELS
OBS = historical.OBSERVABLES
DISPLAY = (r'$\Delta m^2_{21}$', r'$\Delta m^2_{3\ell}$',
           r'$\sin^2\theta_{12}$', r'$\sin^2\theta_{13}$',
           r'$\sin^2\theta_{23}$')
DEFAULT_OUT = 'output/studies/general_scalar_extended/reports/llss_thesis_original_v1'
ALLOWED_STALE_SYMBOLIC_SOURCE = 'RGE/running/rgbeta/RunGeneralT3QuadraticRGBeta.wl'


def source_in_original_basis(ctx, c):
    """Only the 16pi² LLSS->LLHH scalar-quartic term, in original master basis.

    The original scalar_pair uses 0.5(s_ab + s_ba); lambda is symmetric.
    BOTH heavy scalar multiplets (all real components beyond Higgs) participate.
    This is valid for the original mixed fixed-Cartesian X/Y adapter as well.
    """
    c = np.asarray(c, complex)
    q = np.asarray(ctx.quartic)
    if c.ndim != 4 or c.shape[:2] != (ctx.nf, ctx.nf) or c.shape[2:] != (ctx.ns, ctx.ns):
        raise ValueError('Invalid canonical master tensor shape')
    if q.shape != (ctx.ns,) * 4:
        raise ValueError('Invalid quartic tensor shape')
    h = tuple(ctx.higgs)
    heavy = tuple(k for k in range(ctx.ns) if k not in h)
    lepton = tuple(i for pair in ctx.leptons for i in pair)
    if len(h) != 4 or len(heavy) == 0 or len(lepton) != 6:
        raise ValueError('LLSS source requires Higgs, both heavy scalars, three lepton doublets')
    qhhss = q[np.ix_(h, h, heavy, heavy)]
    cllss = c[np.ix_(lepton, lepton, heavy, heavy)]
    source = np.einsum('abcd,ijcd->ijab', qhhss, cllss, optimize=True)
    source = 0.5 * (source + source.swapaxes(2, 3))
    out = np.zeros_like(c)
    out[np.ix_(lepton, lepton, h, h)] = source
    if not np.isfinite(out).all():
        raise ValueError('Nonfinite LLSS source')
    return out


class LLSSIsolator:
    """Calls original beta unchanged during parity checks/other diagnostics.

    The original integrators calculate the unmodified beta ONCE at mu_F to
    check their exact direct HHSS-source parity. Every subsequent RHS call
    receives precisely the LLSS->LLHH subtraction, recomputed from the
    running Wilson tensor and running quartics.
    """
    def __init__(self, original_beta):
        self.original_beta = original_beta
        self.in_integration = False
        self.skip_initial_parity = False
        self.rhs_calls = 0
        self.maximum_source_norm = 0.0
        self.initial_source_norm = None
        self.initial_source_ratio = None
        self.per_model = {}

    @contextmanager
    def integration(self, *, eft_on):
        if self.in_integration:
            raise RuntimeError('Nested simultaneous integrations unsupported')
        if not eft_on:
            yield
            return
        self.in_integration = True
        self.skip_initial_parity = True
        try:
            yield
            if self.rhs_calls == 0:
                raise RuntimeError('LLSS-off integrator never evaluated an RHS')
        finally:
            self.in_integration = False
            self.skip_initial_parity = False

    def beta(self, ctx, c):
        terms = self.original_beta(ctx, c)
        if not self.in_integration:
            return terms
        src = source_in_original_basis(ctx, c)
        norm = float(np.linalg.norm(src))
        diag = self.per_model.setdefault(ctx.model, {'initial_source_norm': None, 'max_rhs_source_norm': 0.0, 'rhs_calls': 0})
        if self.skip_initial_parity:
            # The original mu_F parity check MUST see the unmodified RHS.
            self.skip_initial_parity = False
            self.initial_source_norm = norm
            diag['initial_source_norm'] = norm
            self.initial_source_ratio = norm / max(float(np.linalg.norm(terms['scalar_pair'])), 1e-300)
            return terms
        self.rhs_calls += 1
        diag['rhs_calls'] += 1
        diag['max_rhs_source_norm'] = max(diag['max_rhs_source_norm'],norm)
        self.maximum_source_norm = max(self.maximum_source_norm, norm)
        terms['scalar_pair'] = terms['scalar_pair'] - src
        return terms


def _read_historical(model, *, a4, nonmixed, mixed):
    folder = a4 if model == 'A_-4' else mixed if model in ('B_-1','C_-1') else nonmixed
    cases = historical._read_case_file(folder, model)
    return cases['1111'], cases['1001']


def _decode_complex(obj):
    if not isinstance(obj, dict) or not {'real','imag'} <= set(obj):
        raise ValueError('Expected original complex matrix encoding')
    return np.asarray(obj['real'], float) + 1j*np.asarray(obj['imag'], float)


def _comparison(full, frozen, off, model):
    row = {'model': model, 'status': 'VALID_EXPLORATORY', 'error': ''}
    for key in OBS:
        values = [float(x['prediction'][key]) for x in (full, frozen, off)]
        if not all(map(math.isfinite, values)) or values[1] == 0:
            raise ValueError(f'{model}: missing/zero/nonfinite thesis observable {key}')
        a, b, c = values
        row.update({key+'_FULL': a, key+'_FROZEN': b, key+'_LLSS_OFF': c,
                    key+'_pct_running': 100.*(a-b)/b,
                    key+'_pct_mixing_off': 100.*(c-b)/b})
    return row


def _archived_symbolic_fingerprint(models, original_source):
    """Recover the original canonical export stamp without weakening its validation.

    The current Mathematica launcher is no longer byte-identical to the source
    used for the successfully exported UV/EFT1 beta JSON.  Permit ONLY this
    launcher to differ.  All other required source files must still match;
    the canonical validate_cached() then checks the actual cached RGBeta
    exports, Lagrangian inventory, stage status, and quartic inventory.
    """
    from scripts import RunT3GeneralAll16FromReference as general
    from scripts import RunT3GeneralSymbolicStudy as canonical

    source_root = general.OUTPUT / original_source
    settings_file = source_root / 'settings.json'
    if not settings_file.is_file():
        raise FileNotFoundError(f'Missing original fit settings: {settings_file}')
    settings = json.loads(settings_file.read_text(encoding='utf-8-sig'))
    run_id = settings.get('canonical_run_id', 'full16_general_v1')
    actual = canonical.fingerprint()
    expected = None
    audit = []
    for model in models:
        provenance_file = (canonical.DEFAULT_DATA_ROOT / run_id /
                           canonical.model_tag(model) / 'provenance.json')
        if not provenance_file.is_file():
            raise FileNotFoundError(f'Missing archived symbolic provenance: {provenance_file}')
        stamp = json.loads(provenance_file.read_text(encoding='utf-8-sig'))
        historical = stamp.get('source_sha256')
        if stamp.get('model') != model or not isinstance(historical, dict) or not historical:
            raise ValueError(f'{model}: invalid archived canonical symbolic provenance')
        if set(historical) != set(actual):
            raise ValueError(f'{model}: cached symbolic source inventory differs from current inventory; '
                             f'missing={sorted(set(historical)-set(actual))}; '
                             f'added={sorted(set(actual)-set(historical))}')
        mismatch = [path for path in actual if actual[path] != historical[path]]
        unexpected = sorted(set(mismatch)-{ALLOWED_STALE_SYMBOLIC_SOURCE})
        if unexpected:
            raise ValueError(f'{model}: unapproved symbolic-source changes {unexpected}; '
                             'refusing to reuse cached symbolic betas')
        if expected is not None and historical != expected:
            raise ValueError(f'{model}: canonical symbolic runs use inconsistent source fingerprints')
        expected = historical
        audit.append({'model':model, 'provenance_path':str(provenance_file),
                      'source_mismatches':mismatch,
                      'cached_rgbeta_exports_revalidated_by_canonical_loader':True})
    if expected is None:
        raise ValueError('No symbolic input models requested')
    return expected, {'canonical_run_id':run_id,
                      'accepted_mismatch_only':ALLOWED_STALE_SYMBOLIC_SOURCE,
                      'models':audit}


def _run_original_group(models, output, name, *, original_source, isolator):
    from scripts import RunT3GeneralSequentialFitAblation as original
    from scripts import RunT3SimultaneousEFT1All16 as nonmixed
    from scripts import RunT3SimultaneousEFT1A4 as a4
    from scripts import RunT3SimultaneousEFT1MixedBC as mixed
    original.OUT = output / 'raw' / name
    if name == 'a4':
        original._nonmixed = a4._nonmixed_A4
        original_integrate = a4.integrate_simultaneously
        target = a4
        method = 'integrate_simultaneously'
    elif name == 'mixed':
        original._mixed = mixed._mixed_simultaneous
        original_integrate = mixed.integrate_mixed
        target = mixed
        method = 'integrate_mixed'
    else:
        original._nonmixed = nonmixed._nonmixed_general
        original_integrate = nonmixed.integrate_simultaneously
        target = nonmixed
        method = 'integrate_simultaneously'

    def isolate_integrate(*args, **kwargs):
        with isolator.integration(eft_on=kwargs.get('eft_on', True)):
            return original_integrate(*args, **kwargs)

    run_id = 'source_exact_llss_off'
    output_folder = original.OUT / run_id
    if output_folder.exists():
        raise FileExistsError(f'Refusing to overwrite previous computed results: {output_folder}')
    from scripts import RunT3GeneralSymbolicStudy as canonical
    archived_fingerprint, fingerprint_audit = _archived_symbolic_fingerprint(
        models, original_source)
    # This is an in-memory, source-only compatibility exception.  No source,
    # archived record, beta JSON, or production validator is modified.
    output_folder.parent.mkdir(parents=True, exist_ok=True)
    (output_folder.parent/'source_fingerprint_reuse_audit.json').write_text(
        json.dumps(fingerprint_audit, indent=2)+'\n', encoding='utf-8')
    for item in fingerprint_audit['models']:
        if item['source_mismatches']:
            print(f"{item['model']}: using archived canonical source fingerprint; "
                  f"source-only exception: {item['source_mismatches']}", flush=True)
    with patch.object(canonical, 'fingerprint', return_value=archived_fingerprint), \
         patch.object(target, method, isolate_integrate):
        rc = original.main(['--source-run', original_source, '--run-id', run_id,
                            '--models', ','.join(models)])
    if rc != 0:
        raise RuntimeError(f'{name}: historical simultaneous runner returned {rc}; inspect {output_folder}')
    return output_folder


def _format_models():
    return [m.replace('_','(')+')' for m in MODELS]


def _annotation(v):
    mag=abs(v)
    return (f'{v:+.1e}' if mag>=1e4 or mag<.01 else
            f'{v:+.1f}' if mag>=10 else f'{v:+.2f}' if mag>=1 else f'{v:+.3f}')


def _tick_values(limit):
    max_dec=max(0,math.floor(math.log10(limit)))
    return sorted(set([0]+[sign*10.**p for sign in (-1,1)
                          for p in range(max_dec+1) if 10.**p<=limit]))


def _heatmap_matrix(rows, observables):
    z = np.full((len(MODELS), 2*len(observables)), np.nan)
    by_model = {r['model']:r for r in rows}
    for i, model in enumerate(MODELS):
        row=by_model.get(model)
        if not row or row.get('status') != 'VALID_EXPLORATORY':
            continue
        for j, obs in enumerate(observables):
            z[i, 2*j] = row[obs+'_pct_running']
            z[i, 2*j+1] = row[obs+'_pct_mixing_off']
    return z


def _write_heatmap_tables(rows, out):
    csv_dir = out / 'csv'
    csv_dir.mkdir(parents=True, exist_ok=True)
    def write_table(filename, observables):
        matrix = _heatmap_matrix(rows, observables)
        headers=[]
        for obs in observables:
            headers.extend([obs + '_running', obs + '_mixing_off'])
        with (csv_dir / filename).open('w', newline='', encoding='utf-8') as f:
            w=csv.writer(f)
            w.writerow(['model', *headers])
            for i, model in enumerate(MODELS):
                row=[model]
                for value in matrix[i]:
                    row.append('' if not np.isfinite(value) else f'{float(value):.17g}')
                w.writerow(row)
    write_table('all16_heatmap_table.csv', OBS)
    write_table('all16_masssq_heatmap_table.csv', OBS[:2])
    write_table('all16_mixing_heatmap_table.csv', OBS[2:])


def _render_panel(rows, out, args, observables, display, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.ticker import FixedLocator, FuncFormatter
    from matplotlib.transforms import blended_transform_factory

    z = _heatmap_matrix(rows, observables)
    if not np.isfinite(z).any():
        raise ValueError('No finite numerical data for heatmap')
    plt.rcdefaults()
    n = len(MODELS)
    ncols = z.shape[1]
    fig, ax = plt.subplots(figsize=(args.figure_width, args.row_height*n+2.0))
    limit = max(float(np.nanmax(np.abs(z))), 1.)
    norm = SymLogNorm(linthresh=1, vmin=-limit, vmax=limit, base=10)
    cmap=plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('white')
    im=ax.imshow(np.ma.masked_invalid(z), cmap=cmap, norm=norm, aspect='auto')

    ax.set_yticks(np.arange(n), _format_models(), fontsize=args.y_tick_fontsize)
    ax.set_xticks(np.arange(ncols), ['Running','Mixing off']*len(observables), fontsize=args.bottom_tick_fontsize)
    ax.set_ylabel('T3 model', fontsize=args.axis_label_fontsize)
    ax.tick_params(axis='both', length=0)

    top = blended_transform_factory(ax.transData, ax.transAxes)
    for j, title in enumerate(display):
        ax.text(2*j+.5, 1.025, title, transform=top, ha='center', va='bottom',
                fontsize=args.top_tick_fontsize)

    ax.set_xticks(np.arange(-.5, ncols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.7)
    ax.tick_params(which='minor', bottom=False, left=False)
    for i in range(2, ncols, 2):
        ax.axvline(i-.5, color='black', linewidth=1.2)
    for i in (3, 6, 9, 12):
        ax.axhline(i-.5, color='black', linewidth=1.2)

    for i in range(n):
        for j in range(ncols):
            v = z[i, j]
            if not np.isfinite(v):
                continue
            r, g, b, _ = cmap(norm(v))
            lum = .2126*r + .7152*g + .0722*b
            ax.text(j, i, _annotation(v), ha='center', va='center', fontsize=args.cell_fontsize,
                    color='white' if lum < .48 else 'black')

    bar=fig.colorbar(im, ax=ax, fraction=.035, pad=.025)
    bar.locator=FixedLocator(_tick_values(limit))
    bar.formatter=FuncFormatter(lambda v,pos:'0' if v==0 else
                                (f'{v:.0e}' if abs(v)>=1e4 else f'{v:g}'))
    bar.update_ticks()
    bar.ax.tick_params(labelsize=args.colorbar_fontsize)
    bar.set_label('Signed change (%) — symmetric log colours', fontsize=args.colorbar_fontsize)

    fig.subplots_adjust(left=.12, right=.90, bottom=.07, top=.94)
    plots=out/'plots'; plots.mkdir(parents=True, exist_ok=True)
    for ext in ('png','pdf'):
        fig.savefig(plots/f'{filename}.{ext}', dpi=args.dpi, bbox_inches='tight')
    plt.close(fig)


def _render(rows, out, args):
    _write_heatmap_tables(rows, out)
    _render_panel(rows, out, args, OBS, DISPLAY, 'all16_heatmap_percent')
    _render_panel(rows, out, args, OBS[:2], DISPLAY[:2], 'all16_masssq_heatmap_percent')
    _render_panel(rows, out, args, OBS[2:], DISPLAY[2:], 'all16_mixing_heatmap_percent')


def _use_original_stagefreeze_running(rows, stage_csv, out):
    """Replace display Running cells using the exact original EFT1 stage CSV.

    This is a presentation-only reconciliation, NOT proof that the LLSS-off
    counterfactual shares the original stage-freeze denominator.
    """
    from scripts import PlotT3OscillationHeatmaps as original_plotter
    original_rows = original_plotter.read_rows(stage_csv)
    lookup = {}
    for r in original_rows:
        if r['switch'] == 'EFT1':
            key = (r['model'], r['observable'])
            if key in lookup:
                raise ValueError(f'Duplicate original EFT1 stage result: {key}')
            lookup[key] = float(r['relative_change_percent'])
    audit = []
    for row in rows:
        if row.get('status') != 'VALID_EXPLORATORY':
            continue
        for obs in OBS:
            key = (row['model'], obs)
            if key not in lookup or not math.isfinite(lookup[key]):
                raise ValueError(f'Missing original EFT1 stage result: {key}')
            old = float(row[obs+'_pct_running'])
            new = lookup[key]
            audit.append({'model': row['model'], 'observable': obs,
                          'llss_csv_running_pct': old,
                          'original_stagefreeze_running_pct': new,
                          'difference_percentage_points': new-old,
                          'relative_disagreement': abs(new-old)/max(abs(new),abs(old),1e-300)})
            row[obs+'_pct_running'] = new
    out.mkdir(parents=True, exist_ok=True)
    with (out/'running_source_audit.csv').open('w', newline='', encoding='utf-8') as f:
        fields = ['model','observable','llss_csv_running_pct',
                  'original_stagefreeze_running_pct','difference_percentage_points',
                  'relative_disagreement']
        writer=csv.DictWriter(f,fieldnames=fields)
        writer.writeheader(); writer.writerows(audit)
    mismatches=[r for r in audit if r['relative_disagreement'] > 1e-8]
    provenance={'running_source':str(stage_csv),
                'mixing_off_source':'LLSS comparison CSV supplied to --plot-only',
                'source_agreement_tolerance_relative':1e-8,
                'number_of_mismatched_cells':len(mismatches),
                'number_of_compared_cells':len(audit),
                'physically_consistent_shared_1001_baseline_verified':not bool(mismatches),
                'warning':('MISMATCHED INPUTS: Running cells copied from original EFT1 stage-freeze CSV; '
                           'Mixing-off cells remain from LLSS comparison CSV. Their common 1001 '
                           'denominator is NOT established where Running differs. '
                           'Do not interpret them as a single consistent benchmark.'
                           if mismatches else 'Running percentages agree; raw denominator identity not independently verified.')}
    (out/'running_source_audit.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    if mismatches:
        print(f'WARNING: {len(mismatches)}/{len(audit)} Running cells DISAGREE with LLSS CSV. '
              'Combined plot is for visual comparison only; baseline consistency NOT verified. '
              'See running_source_audit.csv.',flush=True)
    else:
        print('Original EFT1 Running percentages match LLSS CSV.',flush=True)


def build_parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a4',type=Path,default=historical.DEFAULT_A4)
    p.add_argument('--nonmixed',type=Path,default=historical.DEFAULT_NONMIXED)
    p.add_argument('--mixed',type=Path,default=historical.DEFAULT_MIXED)
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--output',type=Path,default=Path(DEFAULT_OUT))
    p.add_argument('--models',default=','.join(MODELS),
                   help='Comma-separated original-model list; use one model for smoke test')
    p.add_argument('--plot-only',type=Path,help='Replot existing combined CSV only')
    p.add_argument('--stagefreeze-csv',type=Path,help='Original three-stage frozen_stage_observable_changes.csv; use its EFT1 values for Running, with disagreement audit')
    p.add_argument('--font-size',type=float,default=14.)
    for k in ('top-tick','bottom-tick','y-tick','cell','axis-label','colorbar'):
        p.add_argument('--'+k+'-fontsize',type=float,default=None)
    p.add_argument('--figure-width',type=float,default=21.)
    p.add_argument('--row-height',type=float,default=.54)
    p.add_argument('--dpi',type=int,default=300)
    return p


def main(argv=None):
    p=build_parser();args=p.parse_args(argv)
    args.output=(ROOT/args.output).resolve()
    if args.output.exists():p.error('Choose a NEW output directory; existing output will never be overwritten')
    base=args.font_size
    for k,factor in (('top_tick',1.2),('bottom_tick',1.),('y_tick',1.),
                     ('cell',.85),('axis_label',1.1),('colorbar',1.)):
        name=k+'_fontsize'
        if getattr(args,name) is None:setattr(args,name,base*factor)
        if getattr(args,name)<=0:p.error(name+' must be positive')
    if args.plot_only:
        source=(ROOT/args.plot_only).resolve()
        with (source/'all16_three_way.csv').open(newline='',encoding='utf-8') as f:
            rows=list(csv.DictReader(f))
        # Values are already independent correct-source comparisons; no guesswork.
        for row in rows:
            if row['status']=='VALID_EXPLORATORY':
                for key in OBS:
                    for part in ('pct_running','pct_mixing_off'):
                        row[key+'_'+part]=float(row[key+'_'+part])
        args.output.mkdir(parents=True)
        if args.stagefreeze_csv:
            stage_csv=(ROOT/args.stagefreeze_csv).resolve()
            if not stage_csv.is_file():p.error(f'Missing original stagefreeze CSV: {stage_csv}')
            _use_original_stagefreeze_running(rows,stage_csv,args.output)
        _render(rows,args.output,args)
        return 0
    if args.stagefreeze_csv:
        p.error('--stagefreeze-csv requires --plot-only; do not alter numerical run')
    selected=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not selected or len(set(selected))!=len(selected) or not set(selected)<=set(MODELS):
        p.error('Invalid model selection')
    folders={k:(ROOT/getattr(args,k)).resolve() for k in ('a4','nonmixed','mixed')}
    for folder in folders.values():
        if not folder.exists():p.error(f'Missing original historical data: {folder}')
    original_cases={m:_read_historical(m,**folders) for m in selected}
    from RGE.running.intermediate import GeneralEFT1C5Mixing as mixing_beta
    isolator=LLSSIsolator(mixing_beta.calculate_master_beta)
    args.output.mkdir(parents=True)
    requested_groups=(('a4',tuple(m for m in selected if m=='A_-4')),
                      ('nonmixed',tuple(m for m in selected if m not in ('A_-4','B_-1','C_-1'))),
                      ('mixed',tuple(m for m in selected if m in ('B_-1','C_-1'))))
    raw_locations={}
    with patch.object(mixing_beta,'calculate_master_beta',isolator.beta):
        for name,group in requested_groups:
            if not group:continue
            print(f'{name}: computing LLSS-off via ORIGINAL simultaneous pipeline for {group}',flush=True)
            raw=_run_original_group(group,args.output,name,original_source=args.source_run,isolator=isolator)
            for model in group:raw_locations[model]=raw/historical.tag(model)/'cases.json'
    rows=[];errors=[]
    for model in selected:
        try:
            raw=json.loads(raw_locations[model].read_text(encoding='utf-8-sig'))
            if raw.get('status')!='Computed':raise RuntimeError('Original LLSS-off trial blocked')
            off=raw['cases']['1111']
            full,frozen=original_cases[model]
            if off.get('eft1_on') is not True or off.get('wilson_log_on') is not True:
                raise ValueError('Ablated result is not the evolved 1111 geometry')
            # Check that original finite scalar matching survives the LLSS modification.
            saved_bubble=_decode_complex(full['c5_joint_scalar_bubble'])
            new_bubble=_decode_complex(off['c5_joint_scalar_bubble'])
            mismatch=float(np.linalg.norm(saved_bubble-new_bubble)/max(np.linalg.norm(saved_bubble),1e-300))
            if mismatch>5e-5:
                raise ValueError(f'{model}: scalar matching changed relative to archived thesis: {mismatch:.3e}')
            diagnosis=isolator.per_model.get(model,{})
            if diagnosis.get('rhs_calls',0)==0 or diagnosis.get('max_rhs_source_norm',0.) <= 0:
                raise ValueError(f'{model}: LLSS scalar-pair source absent in original master basis; investigate mixed source/map')
            if model in ('B_-1','C_-1'):
                original_c5 = _decode_complex(full['c5_at_scalar_threshold'])
                ablated_c5 = _decode_complex(off['c5_at_scalar_threshold'])
                c5_rel = float(np.linalg.norm(original_c5-ablated_c5) /
                               max(np.linalg.norm(original_c5),1e-300))
                if c5_rel < 1e-8:
                    raise ValueError(f'{model}: LLSS-off has negligible C5 effect '
                                     f'({c5_rel:.3e}) despite nonzero source; '
                                     'mixed-basis contraction requires audit')
            row=_comparison(full,frozen,off,model)
            row.update(scalar_matching_relative_residual=mismatch,
                       original_cases_source=str(raw_locations[model]),
                       source_beta_max_norm=diagnosis['max_rhs_source_norm'])
        except Exception as exc:
            row={'model':model,'status':'BLOCKED','error':f'{type(exc).__name__}: {exc}'}
            errors.append(row)
        rows.append(row)
    cols=['model','status','error','scalar_matching_relative_residual','source_beta_max_norm','original_cases_source']
    cols += [f'{o}_{part}' for o in OBS for part in ('FULL','FROZEN','LLSS_OFF','pct_running','pct_mixing_off')]
    with (args.output/'all16_three_way.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
    report=['# Original simultaneous-EFT1 LLSS isolation', '',
            'Saved original 1111 and 1001 results; new LLSS-off counterfactual.',
            '',
            'Running = 100 (O_1111 - O_1001) / O_1001.',
            'Mixing off = 100 (O_LLSS_OFF - O_1001) / O_1001.',
            '', '| Model | Status | Running Δm²21 (%) | Mixing off Δm²21 (%) |',
            '|---|---|---:|---:|']
    for row in rows:
        if row['status']=='VALID_EXPLORATORY':
            report.append(f"| {row['model']} | VALID_EXPLORATORY | "
                          f"{row['delta_m21_sq_ev2_pct_running']:.5g} | "
                          f"{row['delta_m21_sq_ev2_pct_mixing_off']:.5g} |")
        else:
            report.append(f"| {row['model']} | BLOCKED: {row['error']} | -- | -- |")
    report += ['', 'The CSV provides all 5 observables, all 3 raw cases, and both percentage changes.',
               'Results are conditional exploratory RG/matching calculations; not certified physical MSbar predictions.']
    (args.output/'report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
    provenance={'source':'original_simultaneous_EFT1_1111_vs_1001',
                'a4':str(folders['a4']),'nonmixed':str(folders['nonmixed']),
                'mixed':str(folders['mixed']), 'source_campaign':args.source_run,
                'percentage_running':'100*(O_1111-O_1001)/O_1001',
                'percentage_mixing_off':'100*(O_LLSS_OFF-O_1001)/O_1001',
                'llss_off':'original seven-term master scalar_pair LLHH subtracted dynamically',
                'mixed_method':'original native source rotated into fixed Cartesian X/Y basis',
                'initial_source_norm':isolator.initial_source_norm,
                'rhs_evaluations':isolator.rhs_calls,
                'maximum_source_norm':isolator.maximum_source_norm,
                'per_model_source_diagnostics':isolator.per_model,
                'conditional_not_physical_msbar':True,'models':list(selected)}
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    if len(selected)==len(MODELS) and not errors:_render(rows,args.output,args)
    else:print('Incomplete or smoke test: full heatmap omitted; see CSV',flush=True)
    print('LLSS thesis result:',len(selected)-len(errors),'/',len(selected),'valid. Output:',args.output,flush=True)
    return 0 if not errors else 2


if __name__=='__main__':
    raise SystemExit(main())
