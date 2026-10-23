"""Rebuild existing fermion-start sensitivity-oval HTML with stable display framing.

Does not modify sensitivities, covariance ovals, NuFIT/JUNO inputs, or the
original BuildFittedModelEnvelopes.py. Presentation-only wrapper.

Run:
 python -m Numerical.plotting.BuildFermionStartOvalsReadable --case all --mode all
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Numerical.plotting import BuildFittedModelEnvelopes as original

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'output/studies/fermion_start_lambda_comparison'
REPORTS = ROOT / 'Reports/output/display_fermion_start'


def adjusted_template(template: str, *, case: str, mode: str) -> str:
    """Retain the original renderer but prevent numerical-noise-level framing.

    The minimum axis span is a small fraction of the marginal NuFIT sigma.
    Extremely small trajectories then appear nearly stationary, rather than
    being magnified to machine precision. Explicit zoom controls still work.
    """
    needle = 'let px=Number($(\'xPad\').value),py=Number($(\'yPad\').value);'
    assert template.count(needle) == 1, 'Upstream template changed: framing insertion point not found'
    # Bound the *initial automatic* frame but never alter the saved data or
    # empirical oval covariance. nuFIT covariance is in the current units.
    stable = r'''// Display-only stability floor: do not autozoom to roundoff noise.
const targetNames=(D.target&&D.target.observable_names)||[];
const targetCov=(D.target&&D.target.covariance)||[];
function marginalSigma(key){const k=targetNames.indexOf(key);return k>=0&&targetCov[k]&&Number.isFinite(targetCov[k][k])&&targetCov[k][k]>0?Math.sqrt(targetCov[k][k]):null}
const sx=marginalSigma(x),sy=marginalSigma(y);
// Four percent of a 1-sigma width: retains meaningful sensitivity but stops
// microscopic variations from filling the entire viewport.
const floorX=sx?0.04*sx:1e-12*Math.max(1,Math.abs((xmin+xmax)/2));
const floorY=sy?0.04*sy:1e-12*Math.max(1,Math.abs((ymin+ymax)/2));
if(xmax-xmin<floorX){const mid=(xmin+xmax)/2;xmin=mid-floorX/2;xmax=mid+floorX/2;dx=floorX}
if(ymax-ymin<floorY){const mid=(ymin+ymax)/2;ymin=mid-floorY/2;ymax=mid+floorY/2;dy=floorY}
'''
    template = template.replace(needle, stable + '\n' + needle)
    # The original title says UV sensitivity even for this UV-frozen study.
    template = template.replace('Individually fitted T3 model sensitivity',
                                f'Fermion Start — UV Off: {case} / {mode} sensitivity')
    template = template.replace('saved UV sensitivity points',
                                'saved fermion-start sensitivity points')
    template = template.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters.',
                                'Fermion boundary at μF; UV running disabled. Optimised or shared parameters as labelled.')
    # Graph still plots every actual observation. The note avoids implying
    # precision beyond the floating point resolution of saved predictions.
    template = template.replace('descriptive envelopes, NOT model confidence regions.',
                                'descriptive envelopes, NOT model confidence regions. Automatic axes use a minimum span of 4% of marginal NuFIT sigma to avoid magnifying roundoff noise.')
    return template


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', choices=('5e-9', '5e-10', 'all'), default='all')
    p.add_argument('--mode', choices=('optimised', 'shared_a_m2', 'all'), default='all')
    args = p.parse_args()
    cases = ('5e-9', '5e-10') if args.case == 'all' else (args.case,)
    modes = ('optimised', 'shared_a_m2') if args.mode == 'all' else (args.mode,)
    old_html = original.HTML
    try:
        for case in cases:
            for mode in modes:
                source = STUDY / case / 'reports_existing' / mode / 'sensitivity'
                destination = REPORTS / case / 'reused' / mode / 'sensitivity_ovals.html'
                original.HTML = adjusted_template(old_html, case=case, mode=mode)
                info = original.build(source, destination)
                if len(info['models']) != 16:
                    raise RuntimeError(f'{case} {mode}: expected 16 models, found {len(info["models"])}')
                print(f'Validated {case} {mode}: 16/16 models, {len(info["skipped"])} excluded')
    finally:
        original.HTML = old_html


if __name__ == '__main__':
    main()
