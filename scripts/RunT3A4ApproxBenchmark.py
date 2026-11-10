#!/usr/bin/env python3
"""Executable A(-4) rough benchmark: F-hard + S-soft -> SM EFT -> PMNS.

Reads existing F and S numerical endpoints, avoids re-running Matchete/RGBeta.
Results are provisional because full one-loop matching is not yet implemented.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.a4_rough_benchmark import (
    RoughBenchmarkError, estimate_boundary, run_smeft_observables,
)
from Numerical.parallel_t3.background.general_batch1_api import read_json, encode_native, decode_native

B1 = ROOT/'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/A_-4'
DIRECT = ROOT/'output/parallel_t3/a4_direct_hard_soft_bridge_v1/A_-4/A4_direct_hard_soft_formal_audit.json'
RG = ROOT/'output/parallel_t3/a4_native_wilson_rg_v1/A_-4/A4_formal_tree_master_RG.json'
CONFIG = ROOT/'configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json'
OUT = ROOT/'output/parallel_t3/a4_rough_neutrino_benchmark_v1/A_-4'


def _source(paths: dict, name: str, path: Path) -> dict:
    if not path.is_file():
        raise RoughBenchmarkError('Missing input '+str(path))
    paths[name] = {'file': str(path), 'sha256': sha256(path.read_bytes()).hexdigest()}
    return read_json(path)


def compute(batch1_dir: Path, config_path: Path, *, operator_factor: float = 1.,
            mu_low: float = 100., direct_report: Path | None = DIRECT,
            rg_report: Path | None = RG) -> tuple[dict, dict[str,np.ndarray]]:
    sources: dict[str,dict] = {}
    f = _source(sources, 'F_numerical_endpoint', batch1_dir/'fermion_matching_input.json')
    s = _source(sources, 'S_numerical_endpoint', batch1_dir/'scalar_matching_input.json')
    cfg = _source(sources, 'pilot_configuration', config_path)
    if cfg.get('schema') != 't3-thesis-benchmark-two-model-pilot-v1':
        raise RoughBenchmarkError('Wrong thesis benchmark configuration')
    if cfg.get('scales_GeV') != {'UV':1e7,'F':1e5,'S':1e3,'low':100}:
        raise RoughBenchmarkError('Unexpected EFT decoupling scale hierarchy')
    target_file = ROOT/cfg['fitting']['target']
    target = _source(sources, 'NuFIT_target', target_file)
    if target.get('ordering') != 'NO':
        raise RoughBenchmarkError('This rough benchmark uses NO targets only')

    boundary = estimate_boundary(f, s, operator_factor=operator_factor)
    direct = None
    if direct_report is not None and Path(direct_report).is_file():
        direct = _source(sources, 'prior_direct_source_audit', Path(direct_report))
        if direct.get('model') != 'A_-4' or direct.get('source_audit',{}).get('formal_ordered_pole_cancelled') is not True:
            raise RoughBenchmarkError('Prior direct-source audit not in expected state')
        saved = np.asarray(decode_native(direct['numerical_formal_channel']['ordered_formal_finite_remainder_GeVinv']), complex)
        if not np.allclose(boundary.ordered,saved,rtol=1e-11,atol=1e-26):
            raise RoughBenchmarkError('Fresh direct finite remainder differs from prior source audit')

    rged = None
    if rg_report is not None and Path(rg_report).is_file():
        rged = _source(sources,'prior_formal_Wilson_RG_diagnostic',Path(rg_report))
        if rged.get('model') != 'A_-4' or rged.get('status') != 'FormalTreeMasterRGComputed_NotPhysicalC5':
            raise RoughBenchmarkError('Prior Wilson RG audit not in expected state')

    result, arrays = run_smeft_observables(
        boundary, s, target, vev_gev=float(cfg['vev_GeV']), mu_low_gev=mu_low)
    result['format'] = 't3-a4-rough-neutrino-benchmark-v1'
    result['source_files'] = sources
    result['prior_direct_channel_audit_reused'] = direct is not None
    result['prior_Wilson_RG_diagnostic_reused_as_reference_only'] = rged is not None
    result['formal_Wilson_RG_LLHH_norm_GeVinv'] = (
        rged.get('C5_LLHH_muS_formal_frobenius_GeVinv') if rged else None)
    result['ordered_direct_formal_norm_GeVinv'] = float(np.linalg.norm(boundary.ordered))
    result['no_EFT1_leading_log_double_count'] = True
    result['source_terms_included'] = [
        'three-flavour F direct LLHH one-loop hard channel',
        'joint S1/S2 one-loop scalar-soft partner',
        'existing SM + Weinberg one-loop running from 1000 GeV to low scale',
        'existing Takagi/PMNS and NuFIT Gaussian five-observable comparison',
    ]
    result['source_terms_omitted'] = [
        'other F-hard dimension-five matching contributions and local MSbar counterterms',
        'complete EFT1 Wilson anomalous-dimension evolution consistently partitioned with thresholds',
        'remaining finite one-loop two-scalar matching beyond the paired direct channel',
        'independent GammaCC/Weyl normalisation and Matchete/RGBeta quartic dictionary',
    ]
    return result, arrays


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch1-dir', type=Path, default=B1)
    p.add_argument('--config', type=Path, default=CONFIG)
    p.add_argument('--out', type=Path, default=OUT)
    p.add_argument('--operator-factor', type=float, default=1.,
                   help='Assumed signed GammaCC/Weyl factor; not source-certified')
    p.add_argument('--mu-low', type=float, default=100.)
    args=p.parse_args(argv)
    if args.out.exists():
        p.error('Output directory already exists; use --out for a new run')
    try:
        result, matrices = compute(args.batch1_dir,args.config,
                                   operator_factor=args.operator_factor,
                                   mu_low=args.mu_low)
        args.out.mkdir(parents=True, exist_ok=False)
        (args.out/'benchmark.json').write_text(
            json.dumps(result, indent=2, allow_nan=False)+'\n',encoding='utf-8')
        np.savez_compressed(args.out/'numerical_matrices.npz', **matrices)
        obs = result['observables']
        print('A(-4) approximate two-threshold numerical benchmark: COMPUTED')
        print('Assumed GammaCC/Weyl operator factor:',args.operator_factor)
        print('C5 estimate at 1000 GeV norm [GeV^-1]:',
              f"{result['C5_MS_approx_frobenius_GeVinv']:.9e}")
        print('Neutrino masses at 100 GeV [eV]:',obs['masses_eV'])
        print('Delta m21^2, m31^2 [eV^2]:',obs['dm21_sq_eV2'],obs['dm31_sq_eV2'])
        if obs['mixing_defined']:
            print('sin^2(theta12,theta13,theta23):',obs['sin2_theta12'],
                  obs['sin2_theta13'],obs['sin2_theta23'])
            print('NuFIT Gaussian chi2:',
                  result['nufit_gaussian_approximation']['chi2_gaussian_uncorrelated'])
        else:
            print('PMNS mixing unavailable:',obs.get('mixing_blocker'))
        print('Complete renormalised physical matching: NOT COMPUTED')
        print('Output:',args.out/'benchmark.json')
        return 0
    except (OSError,ValueError,RuntimeError,KeyError,TypeError) as err:
        print('BLOCKED:',err,file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
