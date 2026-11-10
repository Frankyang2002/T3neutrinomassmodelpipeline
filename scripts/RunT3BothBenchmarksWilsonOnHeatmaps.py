
#!/usr/bin/env python3
"""Two fixed-Yukawa donors x 16 T3 models: original-axis four heatmaps.

Stage ablations:
  UV    freezes UV renormalisable RG only; EFT1 Wilson and SMEFT run.
  EFT1  freezes EFT1 RENORMALISABLE scalar RG only; FULL Wilson beta RUNS.
  SMEFT freezes final SM+Weinberg running; EFT1 Wilson still runs.
Fermion matching and joint scalar matching are retained in all scenarios.

This reuses the repository's exploratory 32-case full-running baseline. The
results are not the separately corrected D(-2) Ward-verified model and are
not a certified general one-loop matching calculation. Missing results remain
NaN/blank. No Yukawa refit and no production-code modifications.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
import subprocess
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from scripts import RunT3GeneralisedRunningComparisons as general
from Numerical.parallel_t3.wilson import workflow
from scripts.BuildT3GeneralisedThesisTables import load_reference

DEFAULT_TABLES=Path('output/parallel_t3/thesis_generalised_tables_v1')
DEFAULT_REFS=Path('output/parallel_t3/A4_Dm2_targeted_fit01')
DEFAULT_RUNNING=Path('output/parallel_t3/all16_both_donors_wilson_on_v1')
DEFAULT_FIGURES=Path('output/parallel_t3/all16_both_donors_wilson_on_four_heatmaps_v1')
REFERENCE=('D_-2','A_-4')
STAGES=('UV','EFT1','SMEFT')
MODELS=16

def ensure_baseline(tables:Path, reference_results:Path):
    if not (tables/'results.json').is_file() or not (tables/'provenance.json').is_file():
        raise FileNotFoundError('Existing 32-case baseline results/provenance not found: '+str(tables))
    provenance=json.loads((tables/'provenance.json').read_text(encoding='utf-8-sig'))
    if not isinstance(provenance.get('reference_benchmarks'),dict):
        raise ValueError('Baseline provenance missing fixed-Yukawa donor references')
    # The original driver itself enforces EXACT source paths. Here we check
    # donor numerical arrays, avoiding Windows path-case mistakes.
    import numpy as np
    for label in REFERENCE:
        loaded=load_reference(label,[reference_results])
        saved=provenance['reference_benchmarks'].get(label,{})
        for name,saved_name in [('y1','UV_y1'),('y2','UV_y2')]:
            a=np.asarray(loaded[name],dtype=float)
            b=np.asarray(saved.get(saved_name),dtype=float)
            if a.shape!=(3,3) or b.shape!=(3,3) or not np.array_equal(a,b):
                raise ValueError(f'{label} {name}: current donor differs from baseline')
        print(f'Validated donor {label}: {loaded["source"]}',flush=True)
    return provenance

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables',type=Path,default=DEFAULT_TABLES)
    parser.add_argument('--reference-results',type=Path,default=DEFAULT_REFS)
    parser.add_argument('--running-output',type=Path,default=DEFAULT_RUNNING)
    parser.add_argument('--figures-output',type=Path,default=DEFAULT_FIGURES)
    parser.add_argument('--plot-only',action='store_true')
    parser.add_argument('--keep-trial-files',action='store_true')
    parser.add_argument('--debug',action='store_true')
    args=parser.parse_args()
    tables=(ROOT/args.tables).resolve()
    refs=(ROOT/args.reference_results).resolve()
    running=(ROOT/args.running_output).resolve()
    figures=(ROOT/args.figures_output).resolve()
    baseline=ensure_baseline(tables,refs)
    solver=baseline['wilson_solver']
    steps=int(baseline['wilson_steps'])
    if not args.plot_only:
        if running.exists():
            parser.error('Output exists; choose a NEW --running-output: '+str(running))
        # Crucial behavior: the original runner temporarily replaces
        # workflow.beta with general._frozen_wilson_beta when EFT1 is off.
        # Replace that function variable with the ORIGINAL fully coupled
        # Wilson beta, so the patch is a no-op for Wilson evolution.
        original_beta=workflow.beta
        if original_beta is general._frozen_wilson_beta:
            raise RuntimeError('Unexpected identity: production Wilson beta already frozen')
        general._frozen_wilson_beta=original_beta
        print('EFT1-off: freeze IntermediateScalarRunner only; full Wilson beta remains ON',flush=True)
        print(f'Reusing original FULL 32-case baselines; evaluating up to '
              f'{len(REFERENCE)*MODELS*len(STAGES)} counterfactual cases.',flush=True)
        arguments=[
            '--tables',str(tables),
            '--reference-results',str(refs),
            '--reference','both',
            '--output',str(running),
            '--wilson-solver',solver,
            '--wilson-steps',str(steps)]
        if args.keep_trial_files:arguments.append('--keep-trial-files')
        if args.debug:arguments.append('--debug')
        rc=general.main(arguments)
        if rc != 0:
            print(f'Numerical driver returned {rc}; blocked cases remain blank. '
                  'Proceeding with available cases.',flush=True)
        if not (running/'running_comparisons.csv').is_file():
            raise FileNotFoundError('Driver did not write its comparison CSV')
        # The reused original script's provenance says EFT1 Wilson frozen;
        # immediately correct this metadata in this *new* output directory.
        pp=running/'provenance.json'
        doc=json.loads(pp.read_text(encoding='utf-8'))
        doc.update({
            'format':'t3-generalised-stage-freeze-Wilson-ON-v1',
            'eft1_freeze_includes_full_Wilson_running':False,
            'eft1_renormalisable_background_frozen':True,
            'eft1_wilson_evolution_always_active':True,
            'uv_freeze_only_uv_renormalisable_rg':True,
            'smeft_stage_off_disables_final_Weinberg_running':True,
            'baseline_32case_reused':str(tables),
            'generalised_exploratory_matching_not_certified':True,
            'not_the_new_Ward_corrected_Dm2_single_benchmark':True,
        })
        pp.write_text(json.dumps(doc,indent=2)+'\n',encoding='utf-8')
        # Reused analysis file descriptions would otherwise misstate the
        # ablation as turning off Wilson. Correct those descriptions too.
        for ref in REFERENCE:
            md=running/f'Section_5_2_results_{ref}.md'
            if md.is_file():
                s=md.read_text(encoding='utf-8')
                s=s.replace('EFT1-off means both scalar RGBeta evolution AND full dim-5 Wilson evolution are frozen.',
                            'EFT1-off means scalar RGBeta evolution is frozen; full dimension-five Wilson evolution REMAINS ON.')
                md.write_text(s,encoding='utf-8')
            tex=running/f'Section_5_2_results_{ref}.tex'
            if tex.is_file():
                s=tex.read_text(encoding='utf-8')
                s=s.replace('The EFT1-off calculation freezes both the renormalisable scalar background and',
                            'The EFT1-off calculation freezes only the renormalisable scalar background;')
                s=s.replace('the coupled dimension-five Wilson coefficients.',
                            'the coupled dimension-five Wilson coefficients continue evolving.')
                tex.write_text(s,encoding='utf-8')
    else:
        if not (running/'running_comparisons.csv').is_file():
            parser.error('Plot-only requires an existing completed running comparison CSV')

    if figures.exists() and any(figures.iterdir()):
        parser.error('Figures output exists and is nonempty; choose a NEW --figures-output')
    with (running/'running_comparisons.csv').open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    count=sum(r.get('status')=='VALID_EXPLORATORY' for r in rows)
    print(f'Valid counterfactuals: {count}/{len(rows)}; blocked cells will be blank',flush=True)
    rc=subprocess.run(
        [sys.executable,str(ROOT/'scripts/PlotT3GeneralisedOriginalHeatmaps.py'),
         '--running',str(running),'--reference','both','--output',str(figures)],
        cwd=ROOT,check=False).returncode
    if rc:raise RuntimeError('Original heatmap renderer returned '+str(rc))
    outputs=[]
    for label in REFERENCE:
        for kind in ('masssq','mixing'):
            image=figures/label/'plots'/f'oscillation_{kind}_heatmap_percent.png'
            if not image.is_file():raise FileNotFoundError('Missing figure: '+str(image))
            outputs.append(str(image))
    status={
        'schema':'t3-both-donors-wilson-on-heatmaps-v1',
        'reference_benchmarks':list(REFERENCE),
        'requested_counterfactuals':2*16*3,
        'computed_counterfactuals':count,
        'wilson_running_enabled_during_eft1_background_freeze':True,
        'smeft_running_frozen_only_for_smeft_ablation':True,
        'full_running_baselines_reused':str(tables),
        'physics_status':'exploratory; canonical D(-2) Ward fixes not propagated to 16 models',
        'figures':outputs}
    figures.mkdir(parents=True,exist_ok=True)
    (figures/'run_summary.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print('\nFour original-axis figures:',flush=True)
    for f in outputs:print('  '+f,flush=True)
    print('No Yukawa refitting, no changes to production files.',flush=True)

if __name__=='__main__':
    main()
