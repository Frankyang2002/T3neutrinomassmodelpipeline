"""Rough T3 fermion-first EFT1 -> joint-scalar SMEFT -> low-energy neutrinos.

Project path: output/general_observable/lo_fermion_boundary/run_rough_smeft_neutrinos.py

Consumes the EXISTING all-16 exploratory v6 joint-scalar output (source-resolved
LLSS at muF, fixed-order EFT1 C5 leading log, simultaneous scalar bubble at
muS). Assumes a user-chosen signed overall GammaCC-to-Weyl factor, default +1.
Runs the repository's SM+Weinberg RGE from muS to mu_low and computes the
Majorana neutrino masses in the charged-lepton mass basis. No certification.
Nothing in the existing release, scalar matching or running code is modified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_SOURCE = ROOT / 'output/general_observable/lo_joint_scalar_matching/exploratory_all16_v6_native_fix'
DEFAULT_RUNNING = ROOT / 'output/general_observable/all16_v1'
DEFAULT_OUTPUT = HERE / 'rough_smeft_pilot_v1'
BENCHMARKS = ('A_-4', 'B_-1', 'C_-1')
FORMAT = 't3-rough-fermion-first-smeft-neutrinos-v1'
VEV = 246.22


def slug(model: str) -> str:
    return model.replace('-', 'm').replace('+', 'p')


def read_json(path: Path) -> dict:
    with Path(path).open(encoding='utf-8-sig') as f:
        x = json.load(f)
    if not isinstance(x, dict):
        raise ValueError(f'Expected JSON object: {path}')
    return x


def hash_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decode_matrix(raw: dict, *, label: str = 'matrix') -> np.ndarray:
    if not isinstance(raw, dict) or not {'real', 'imag'} <= raw.keys():
        raise ValueError(f'{label}: expected real/imag matrices')
    try:
        re = np.asarray(raw['real'], dtype=float)
        im = np.asarray(raw['imag'], dtype=float)
        z = re + 1j*im
    except (ValueError, TypeError) as exc:
        raise ValueError(f'{label}: numeric matrix required') from exc
    if z.shape != (3, 3) or not (np.isfinite(re).all() and np.isfinite(im).all()):
        raise ValueError(f'{label}: expected finite complex 3x3')
    return z


def encode_matrix(z: np.ndarray) -> dict:
    a = np.asarray(z, complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError('Cannot encode non-finite 3x3 matrix')
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def decoded_state_matrix(raw: object, label: str) -> np.ndarray:
    """Historical all16 running state encodes each element as {real,imag}."""
    if not isinstance(raw, list) or len(raw) != 3 or any(not isinstance(r, list) or len(r) != 3 for r in raw):
        raise ValueError(f'{label}: expected 3x3 list')
    try:
        out = np.array([[complex(float(x['real']), float(x['imag']))
                         if isinstance(x, dict) else complex(x)
                         for x in row] for row in raw], dtype=complex)
    except (TypeError, ValueError, KeyError) as exc:
        raise ValueError(f'{label}: unsupported complex encoding') from exc
    if not np.isfinite(out).all():
        raise ValueError(f'{label}: non-finite entries')
    return out


def take_scalar_c5(src: dict, model: str, factor: float) -> tuple[np.ndarray, dict]:
    """Read complete fixed-order LLSS C5; never add the log or bubble twice."""
    if src.get('model') != model or src.get('format') != 't3-conditional-c5-per-unknown-gammacc-mixed-split-v6':
        raise ValueError(f'{model}: wrong source-resolved joint-scalar result')
    resp = src.get('response')
    if not isinstance(resp, dict) or resp.get('model') != model:
        raise ValueError(f'{model}: missing v6 response')
    if resp.get('status') != 'ConditionalOneLoopLLSSResponse_PerUnitUncertifiedSpinorNormalization':
        raise ValueError(f'{model}: scalar matching result is not complete')
    if not math.isfinite(factor) or factor == 0:
        raise ValueError('Choose a finite nonzero signed spinor factor')
    total = decode_matrix(resp['C5_LLSS_total_one_loop_per_N_gev_inverse'], label='total C5')
    log = decode_matrix(resp['C5_EFT1_induced_leading_log_per_N_gev_inverse'], label='EFT1 C5 log')
    bubble = decode_matrix(resp['C5_joint_scalar_finite_bubble_per_N_gev_inverse'], label='scalar bubble')
    scale = max(float(np.max(abs(total))), 1.e-30)
    if float(np.max(abs(total-log-bubble))) > 1e-10*scale + 1e-28:
        raise ValueError(f'{model}: leading log + joint bubble differs from saved total')
    if not np.allclose(total, total.T, rtol=1e-10, atol=1e-27):
        raise ValueError(f'{model}: C5 is not complex symmetric')
    if (not resp.get('mu_f_gev', 0) > resp.get('mu_s_gev', 0) > 0 or
            resp.get('renormalization_scheme') != 'MSbar'):
        raise ValueError(f'{model}: malformed EFT scales/scheme')
    return factor * total, {
        'mu_F_gev': float(resp['mu_f_gev']),
        'mu_S_gev': float(resp['mu_s_gev']),
        'C5_EFT1_leading_log_gev_inv': encode_matrix(factor*log),
        'C5_joint_scalar_finite_bubble_gev_inv': encode_matrix(factor*bubble),
        'C5_LLSS_total_gev_inv': encode_matrix(factor*total),
        'scalar_eigenmasses_gev': resp.get('scalar_eigenstate_masses_gev'),
        'selected_PL_source_count': len(resp.get('per_PL_source', [])),
        'matching_sequence': resp.get('sequence'),
    }


def compute_model(model: str, scalar_path: Path, running_path: Path,
                  *, normalization: float = 1.0, mu_low: float = 91.1876,
                  vev: float = VEV, ordering: str = 'AUTO') -> dict:
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions, evolve_sm_weinberg
    from physics.NeutrinoTrajectory import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables

    c5data = read_json(scalar_path)
    c5_s, origin = take_scalar_c5(c5data, model, normalization)
    trajectory = read_json(running_path)
    if trajectory.get('model') != model or trajectory.get('status') != 'NumericalRenormalisableRunningOnly':
        raise ValueError(f'{model}: missing numerical EFT1 running state')
    scales = trajectory['scales_gev']
    mu_s = origin['mu_S_gev']
    if not math.isclose(float(scales['S']), mu_s, rel_tol=1e-9) or not math.isclose(float(scales['F']), origin['mu_F_gev'], rel_tol=1e-9):
        raise ValueError(f'{model}: running and scalar matching scales differ')
    if not (math.isfinite(mu_low) and 0 < mu_low <= mu_s and math.isfinite(vev) and vev > 0):
        raise ValueError('Require 0 < mu_low <= mu_S and a positive finite Higgs vev')
    sm = trajectory['eft1_at_scalar_threshold']
    ini = SMWeinbergInitialConditions(
        gY=float(sm['gY']), g2=float(sm['g2']), g3=float(sm['g3']), lambdaH=float(sm['lambdaH']),
        ye=decoded_state_matrix(sm['ye'], 'ye'),
        yu=decoded_state_matrix(sm['yu'], 'yu'),
        yd=decoded_state_matrix(sm['yd'], 'yd'), K=c5_s,
    )
    evolved = evolve_sm_weinberg(ini, mu_s, mu_low) if mu_low != mu_s else ini
    c5_low = evolved.K
    ye_low = evolved.ye
    mass_flavour_gev = charged_lepton_mass_basis_matrix(
        neutrino_mass_from_c5(c5_low, vev_gev=vev), ye_low)
    obs = calculate_neutrino_observables(mass_flavour_gev, ordering=ordering)
    return {
        'format': FORMAT,
        'status': 'RoughSMEFTAndNeutrinosComputed', 'model': model,
        'assumed_signed_GammaCC_to_Weyl_factor': normalization,
        'vev_gev': vev, 'mu_low_gev': mu_low,
        'C5_at_scalar_threshold_gev_inv': encode_matrix(c5_s),
        'C5_at_low_scale_gev_inv': encode_matrix(c5_low),
        'mass_matrix_charged_lepton_basis_gev': encode_matrix(mass_flavour_gev),
        'mass_matrix_charged_lepton_basis_ev': encode_matrix(mass_flavour_gev*1e9),
        'masses_ev': obs.masses_ev.tolist(),
        'delta_m21_sq_ev2': obs.delta_m21_sq_ev2,
        'delta_m31_sq_ev2': obs.delta_m31_sq_ev2,
        'delta_m32_sq_ev2': obs.delta_m32_sq_ev2,
        'pmns_abs': obs.pmns_abs.tolist(),
        'ordering': obs.ordering,
        'takagi_residual': obs.takagi_residual,
        'scalar_matching': origin,
        'provenance': {
            'scalar_source': str(scalar_path), 'scalar_source_sha256': hash_file(scalar_path),
            'running_source': str(running_path), 'running_source_sha256': hash_file(running_path),
            'scalar_matching_worker': 'scripts/RunT3JointScalarAll16V6.py (precomputed result)',
            'low_energy_worker': 'Numerical.running.SMWeinbergEvolution.evolve_sm_weinberg',
        },
        'approximations': [
            'Illustrative original benchmark parameters; not a neutrino oscillation fit',
            'Assumed signed Matchete GammaCC-to-master-Weyl normalization (default +1)',
            'Source-to-RGBeta scalar quartic aliases as used by v6',
            'Leading-log EFT1 C5 plus one simultaneous finite scalar bubble, not full one-loop matching',
            'No additional hard heavy-fermion one-loop matching or full EFT1 Wilson evolution',
        ],
    }


def models_for(args: argparse.Namespace) -> list[str]:
    from Numerical.joint_scalar_matching.adapter import ALL_MODELS, normalize_model
    if args.all16:
        return list(ALL_MODELS)
    if args.model:
        return [normalize_model(args.model)]
    return list(BENCHMARKS)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    choice = ap.add_mutually_exclusive_group()
    choice.add_argument('--model', help='e.g. B_-1; defaults to A_-4 B_-1 C_-1')
    choice.add_argument('--all16', action='store_true')
    ap.add_argument('--source-dir', type=Path, default=DEFAULT_SOURCE,
                    help='Successful v6 exploratory scalar matching directory')
    ap.add_argument('--running-dir', type=Path, default=DEFAULT_RUNNING)
    ap.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument('--normalization', type=float, default=1.0,
                    help='ASSUMED signed GammaCC-to-Weyl factor, default +1')
    ap.add_argument('--mu-low', type=float, default=91.1876)
    ap.add_argument('--vev', type=float, default=VEV)
    ap.add_argument('--ordering', choices=['AUTO', 'NO', 'IO'], default='AUTO')
    args = ap.parse_args(argv)
    if not math.isfinite(args.normalization) or args.normalization == 0:
        ap.error('--normalization must be finite and nonzero')
    out = args.output_dir.resolve()
    # Only create new outputs under the current pilot work area, never overwrite releases.
    if not (out == HERE or HERE in out.parents):
        ap.error(f'--output-dir must be inside {HERE}')
    models = models_for(args)
    rows, results = [], []
    try:
        if (out/'summary.json').exists():
            raise FileExistsError(f'Pilot summary exists; use a NEW --output-dir: {out}')
        for model in models:
            path = args.source_dir/'C5_per_unknown_NGamma_exploratory'/f'{slug(model)}.json'
            running_path = args.running_dir/slug(model)/'running.json'
            record = compute_model(model, path, running_path,
                                   normalization=args.normalization, mu_low=args.mu_low,
                                   vev=args.vev, ordering=args.ordering)
            out.mkdir(parents=True, exist_ok=True)
            target = out/f'{slug(model)}.json'
            if target.exists():
                raise FileExistsError(f'Existing pilot model output: {target}')
            target.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            masses = record['masses_ev']
            rows.append({'model': model, 'mu_low_gev': args.mu_low,
                         'm1_ev': masses[0], 'm2_ev': masses[1], 'm3_ev': masses[2],
                         'dm21_sq_ev2': record['delta_m21_sq_ev2'],
                         'dm31_sq_ev2': record['delta_m31_sq_ev2'],
                         'ordering': record['ordering'],
                         'C5_frobenius_mu_low_gev_inv': float(np.linalg.norm(
                             decode_matrix(record['C5_at_low_scale_gev_inv']))),
                         'normalization_assumed': args.normalization})
            results.append({'model': model, 'output_json': str(target),
                            'masses_ev': masses, 'ordering': record['ordering']})
            print(f'{model}: m_nu [eV] = '+', '.join(f'{x:.6g}' for x in masses), flush=True)
        csv_path = out/'rough_neutrino_masses.csv'
        with csv_path.open('w', newline='', encoding='utf-8') as fp:
            writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        summary = {'format': FORMAT, 'status': 'RoughSMEFTPilotCompleted',
                   'models_computed': len(results), 'models': results,
                   'mu_low_gev': args.mu_low, 'assumed_spinor_factor': args.normalization,
                   'scalar_matching_recomputed': False,
                   'csv': str(csv_path),
                   'sequence': 'UV -> integrate out all fermions -> EFT1 with both scalars -> integrate out all scalars -> SMEFT -> low-energy m_nu',
                   'rough_only_not_a_physical_fit': True}
        (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
        print('CREATED:', csv_path)
        print('CREATED:', out/'summary.json')
        return 0
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
