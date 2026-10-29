"""Non-destructive 16-model validation of nonzero allowed T3 quadratic masses.

Checks the full Cartesian Hessian against independently predicted spectra,
phase invariance, the mu2->0 limit, and (for B/C(-1)) the independently
implemented complex-doublet rotation and full quadratic-tensor transport.

This does NOT certify Matchete/RGBeta couplings or physical threshold matching.

Run: python -m scripts.RunT3QuadraticMixingValidation
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    analytic_mass_squared, load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
    MixedDoubletMassBasis, quadratic_mass_tensor_transport,
)

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'docs/general_scalar_trial_operator_manifest.json'
DEFAULT_OUTPUT = ROOT / 'output/audits/t3_nonzero_quadratic_validation.json'
MIXED = {'B_-1', 'C_-1'}
EXPECTED_KINDS = {
    'A_+0': 'self_s1', 'A_-2': 'self_s2',
    'B_-1': 'mixed', 'C_-1': 'mixed',
    'D_+0': 'self_s1', 'D_-2': 'self_s2',
    'E_+0': 'self_s1', 'E_-2': 'self_s2',
}


def check(condition: bool, explanation: str) -> None:
    if not condition:
        raise ValueError('T3 quadratic validation: ' + explanation)


def close(x, y, *, rtol: float = 2e-11, atol: float = 1e-6) -> bool:
    return bool(np.allclose(x, y, rtol=rtol, atol=atol))


def validate_model(manifest: dict, model: str, m1: float, m2: float,
                   mu: complex) -> dict:
    kind = quadratic_kind(manifest, model)
    check(kind == EXPECTED_KINDS.get(model, 'none'), model + ': wrong quadratic kind')
    chosen = mu if kind != 'none' else 0j
    zero = quadratic_mass_hessian(manifest, model, m1, m2, 0j)
    result = quadratic_mass_hessian(manifest, model, m1, m2, chosen)
    spec = np.asarray(result['mass_squared_eigenvalues_gev2'], dtype=float)
    pred = analytic_mass_squared(manifest, model, m1, m2, chosen)
    base = np.asarray(zero['mass_squared_eigenvalues_gev2'], dtype=float)
    check(close(spec, pred), model + ': analytical spectrum mismatch')
    check(not bool(result['tachyonic']), model + ': tachyonic benchmark masses')
    check(bool(result['analytic_spectrum_verified']), model + ': missing Hessian check')
    check(base.shape == spec.shape, model + ': changing number of real modes')
    check(close(base, analytic_mass_squared(manifest, model, m1, m2, 0j)),
          model + ': zero coefficient spectrum mismatch')
    check(close(spec, np.asarray(quadratic_mass_hessian(
        manifest, model, m1, m2, chosen * np.exp(0.81j)
    )['mass_squared_eigenvalues_gev2'])), model + ': mass spectrum depends on mu phase')
    check(close(spec, np.asarray(quadratic_mass_hessian(
        manifest, model, m1, m2, chosen.conjugate()
    )['mass_squared_eigenvalues_gev2'])), model + ': conjugation changes spectrum')
    check(close(spec, np.asarray(quadratic_mass_hessian(
        manifest, model, m1, m2, -chosen
    )['mass_squared_eigenvalues_gev2'])), model + ': sign changes spectrum')
    check(close(base, np.asarray(quadratic_mass_hessian(
        manifest, model, m1, m2, chosen * 1e-12
    )['mass_squared_eigenvalues_gev2'])), model + ': small-mu limit incorrect')
    output = {
        'quadratic_kind': kind,
        'nonzero_mu2_tested': kind != 'none',
        'mass_squared_eigenvalues_gev2': spec.tolist(),
        'minimum_mass_squared_gev2': float(spec[0]),
        'analytic_hessian_agreement': True,
        'phase_conjugation_sign_invariance': True,
        'zero_boundary_limit_verified': True,
        'mixed_doublet_tensor_transport': 'not_applicable',
    }
    if model in MIXED:
        check(kind == 'mixed', model + ': expected mixed doublet mass')
        basis = MixedDoubletMassBasis.build(model, m1, m2, chosen)
        repeated = np.repeat(basis.eigenvalues_sq, 4)
        check(close(spec, repeated), model + ': real/complex mass basis disagreement')
        tensor = quadratic_mass_tensor_transport(basis)
        # The helper raises on invalid contraction or non-diagonal terms.
        check(isinstance(tensor, dict), model + ': no tensor transport report')
        output['mixed_doublet_tensor_transport'] = 'verified_by_existing_helper'
        output['sin_2theta'] = float(basis.mixing_sin_2theta)
        output['mass_eigenvalues_sq_gev2'] = basis.eigenvalues_sq.tolist()
    return output


def run(*, m_s1_gev: float = 1350., m_s2_gev: float = 1780.,
        mu2_real_gev2: float = 25000., mu2_imag_gev2: float = 12000.,
        manifest_path: Path = MANIFEST) -> dict:
    args = (m_s1_gev, m_s2_gev, mu2_real_gev2, mu2_imag_gev2)
    check(all(math.isfinite(x) for x in args), 'nonfinite input')
    check(m_s1_gev > 0 and m_s2_gev > 0, 'positive input masses required')
    check(abs(complex(mu2_real_gev2, mu2_imag_gev2)) > 0,
          'nonzero extra mu2 required to test mixing')
    manifest = load_manifest(manifest_path)
    m1, m2 = m_s1_gev**2, m_s2_gev**2
    rows = {}
    for model in sorted(manifest['models']):
        rows[model] = validate_model(
            manifest, model, m1, m2, complex(mu2_real_gev2, mu2_imag_gev2)
        )
    check(len(rows) == 16, 'not 16 model cases')
    check(sum(r['nonzero_mu2_tested'] for r in rows.values()) == 8,
          'expected exactly eight nonzero quadratic tests')
    check(sum(r['mixed_doublet_tensor_transport'] == 'verified_by_existing_helper'
              for r in rows.values()) == 2,
          'expected exactly two mixed doublet transport checks')
    return {
        'status': 'CandidateQuadraticMassSpectrumAndTransportValidated',
        'models': rows,
        'models_checked': 16,
        'extra_quadratic_model_slots_checked_nonzero': 8,
        'mixed_doublet_tensor_checks': 2,
        'input': {'m_s1_gev': m_s1_gev, 'm_s2_gev': m_s2_gev,
                  'mu2_real_gev2': mu2_real_gev2,
                  'mu2_imag_gev2': mu2_imag_gev2},
        'source_manifest_sha256': sha256(manifest_path.read_bytes()).hexdigest(),
        'global_matchete_rgbeta_normalisation_certified': False,
        'physical_scalar_threshold_matching_performed': False,
        'full_general_one_loop_c5_certified': False,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ms1-gev', type=float, default=1350.)
    p.add_argument('--ms2-gev', type=float, default=1780.)
    p.add_argument('--mu2-real-gev2', type=float, default=25000.)
    p.add_argument('--mu2-imag-gev2', type=float, default=12000.)
    p.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    a = p.parse_args()
    result = run(m_s1_gev=a.ms1_gev, m_s2_gev=a.ms2_gev,
                 mu2_real_gev2=a.mu2_real_gev2,
                 mu2_imag_gev2=a.mu2_imag_gev2)
    destination = a.output if a.output.is_absolute() else ROOT / a.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print('Status:', result['status'])
    print('Models:', result['models_checked'], 'nonzero quadratic slots:',
          result['extra_quadratic_model_slots_checked_nonzero'])
    print('Mixed doublet tensor checks:', result['mixed_doublet_tensor_checks'])
    print('Physical threshold matching: NOT CERTIFIED')
    print('Report:', destination)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
