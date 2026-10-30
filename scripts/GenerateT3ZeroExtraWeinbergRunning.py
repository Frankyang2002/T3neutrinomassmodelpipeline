"""Generate verified 1 TeV -> 100 GeV C5 and neutrino-observable curves.

Uses the *existing* frozen-64 Weinberg ODE calculator on the separate saved
zero-extra control. No intermediate C12 trajectory is invented.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import T3Frozen64Sources as policy
from scripts import T3ZeroExtraReportSources as zero
from scripts import GenerateT3Frozen64WeinbergRunning as old


def generate(run_id='zero_extra_all16_low_v1', output_id='weinberg_running_v1',
             points=41, resume=False):
    policy.require(bool(policy.RUN_ID.fullmatch(output_id)) and
                   output_id not in zero.MODES, 'unsafe output ID')
    meta = zero.load_summary(run_id)
    base = zero.run_root(run_id)
    generator_hash = policy.digest(Path(__file__))
    old_generator_hash = policy.digest(Path(old.__file__))
    n = 0
    for mode in zero.MODES:
        for model in policy.MODELS:
            src = zero.sources(run_id, meta, mode, model)
            out = base / output_id / mode / (policy.tag(model) + '.json')
            identity = {'control_summary_sha256': meta['summary_sha256'],
                        'control_settings_sha256': meta['settings_sha256'],
                        'trajectory_sha256': src['cell']['trajectory_sha256'],
                        'neutrino_sha256': src['cell']['neutrino_sha256'],
                        'adapter_generator_sha256': generator_hash,
                        'weinberg_solver_adapter_sha256': old_generator_hash,
                        'points_requested': points}
            if out.exists():
                policy.require(resume, f'previous output exists; specify --resume: {out}')
                saved = policy.read_json(out)
                policy.require(saved.get('source') == identity and
                               saved.get('status') == 'VerifiedSMWeinbergEndPointParity' and
                               saved.get('model') == model and saved.get('mode') == mode,
                               'cached Weinberg data have changed; use another output ID')
                n += 1
                continue
            value = old.calculate(src['trajectory'], src['neutrino'], points=points)
            value.update(mode=mode, benchmark=meta['benchmark'], source=identity)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            n += 1
            print(f'{mode}/{model}: {points}-point C5 and observables validated', flush=True)
    summary = {'run_id': run_id, 'output_id': output_id,
               'benchmark': meta['benchmark'], 'verified': n, 'expected': 32,
               'c12_intermediate_numerical_trajectory_available': False,
               'physical_sequential_matching_certified': False}
    (base / output_id / 'summary.json').write_text(
        json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', default='zero_extra_all16_low_v1')
    p.add_argument('--output-id', default='weinberg_running_v1')
    p.add_argument('--points', type=int, default=41)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args(argv)
    print(json.dumps(generate(a.run_id, a.output_id, a.points, a.resume), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
