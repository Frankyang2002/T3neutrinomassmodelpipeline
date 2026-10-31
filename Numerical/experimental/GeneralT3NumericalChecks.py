"""Stage 4 read-only numerical checks for conditional T3 Wilson/C5 artifacts.

The original B/C outputs are CONDITIONAL and are not physical predictions.
Numerical parity never upgrades missing Y/X finite matching to a known zero.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from Numerical.experimental.GeneralT3EFTHandoffController import (
    BlockSpec, HandoffContractError,
)
from physics.NeutrinoMass import numerical_neutrino_mass_matrix
from physics.NeutrinoObservables import takagi_factorization

EXPECTED_BC_MASSES_EV = {
    'B_-1': (0.03853362292521312, 0.01561382145082619, 0.00045848474457153946),
    'C_-1': (0.046405741900074446, 0.045923847867146636, 0.021684644219742994),
}
REPORT_FORMAT = 't3-bc-conditional-sequential-hx-to-sm-v1'
REPORT_STATUS = 'EFT2HXRunningXCrossingSMWeinbergConditionalFiniteScalarMatchingMissing'


def complex_matrix_payload(raw: Any, *, shape: tuple[int, ...]) -> np.ndarray:
    """Read the exact {real,imag} encoding in the saved B/C reports."""
    if not isinstance(raw, Mapping) or set(raw) != {'real', 'imag'}:
        raise HandoffContractError('Complex payload requires exactly real and imag arrays')
    try:
        real = np.asarray(raw['real'], dtype=float)
        imag = np.asarray(raw['imag'], dtype=float)
    except (ValueError, TypeError, OverflowError) as exc:
        raise HandoffContractError('Non-numeric complex matrix payload') from exc
    if real.shape != shape or imag.shape != shape or not (np.isfinite(real).all() and np.isfinite(imag).all()):
        raise HandoffContractError(f'Invalid complex payload, expected {shape}')
    return real + 1j*imag


def validate_c5(c5: Any) -> np.ndarray:
    """Require complex symmetric flavour C5, NOT Hermitian or diagonal C5."""
    spec = BlockSpec('C5', (3,3), 'complex', True, units='GeV^-1',
                     operator_basis='Weinberg / external-provenance-required',
                     normalization='mnu=-(v^2/2)C5 / project-convention', mass_dimension=-1)
    z = spec.validate(c5)
    scale = max(float(np.max(np.abs(z))), 1e-30)
    if float(np.max(np.abs(z-z.T))) / scale > 1e-9:
        raise HandoffContractError('C5 is not transpose-symmetric at its physical magnitude')
    return z


@dataclass(frozen=True)
class NeutrinoMatrixCheck:
    masses_ev_descending: np.ndarray
    takagi_residual: float
    svd_max_relative_error: float
    neutrino_mass_gev: np.ndarray


def check_neutrino_matrix(c5: Any, *, vev_gev: float = 246.22,
                          max_residual: float = 1e-7) -> NeutrinoMatrixCheck:
    """Call existing project mnu and Takagi routines; never derive a new C5 map."""
    matrix = validate_c5(c5)
    if not np.isfinite(vev_gev) or vev_gev <= 0:
        raise HandoffContractError('Invalid vev')
    m = numerical_neutrino_mass_matrix(matrix, vev_gev=vev_gev)
    _, takagi_masses, residual = takagi_factorization(m)
    reference_svd = np.linalg.svd(m, compute_uv=False)
    takagi_desc = takagi_masses[::-1]
    delta = float(np.max(np.abs(takagi_desc-reference_svd)) /
                  max(float(np.max(reference_svd)), 1e-30))
    if not np.isfinite(residual) or residual > max_residual or delta > max_residual:
        raise HandoffContractError('Native Takagi residual or singular masses disagree')
    return NeutrinoMatrixCheck(takagi_desc*1e9, float(residual), delta, m.copy())


def validate_bc_wilson_tensor(value: Any) -> np.ndarray:
    """Validate project B/C 6-independent-flavour-pair, real-scalar tensor.

    Array layout (3,3,N_s,N_s,6), N_s=8 (H,X) or 12 (H,S1,S2).
    The last axis follows T3BCCoupledWilsonTensor.PAIRS. This only checks
    array shape/symmetry, NOT a physical operator normalization certificate.
    """
    z = np.asarray(value, dtype=complex)
    if z.ndim != 5 or z.shape[:2] != (3,3) or z.shape[2] not in (8,12) or \
            z.shape[3] != z.shape[2] or z.shape[4] != 6 or not np.isfinite(z).all():
        raise HandoffContractError('B/C Wilson tensor expected shape (3,3,8|12,8|12,6)')
    magnitude = max(float(np.max(np.abs(z))), 1e-30)
    lf = float(np.max(np.abs(z-z.swapaxes(0,1))))/magnitude
    sf = float(np.max(np.abs(z-z.swapaxes(2,3))))/magnitude
    if max(lf,sf) > 1e-8:
        raise HandoffContractError('B/C Wilson pair or scalar exchange symmetry fails')
    return np.array(z, copy=True)


def bc_wilson_block(name: str, n_scalar: int, *, operator_basis: str,
                    normalization: str, mass_dimension: int) -> BlockSpec:
    """Declare Wilson shape and metadata; independent physical verification is external."""
    if n_scalar not in (8,12) or not isinstance(name, str) or not name.strip():
        raise HandoffContractError('Invalid B/C tensor block specification')
    if not operator_basis.strip() or not normalization.strip() or type(mass_dimension) is not int:
        raise HandoffContractError('Explicit Wilson provenance is mandatory')
    return BlockSpec(name, (3,3,n_scalar,n_scalar,6), 'complex',
                     operator_basis=operator_basis, normalization=normalization,
                     mass_dimension=mass_dimension)


def check_bc_conditional_report(report: Mapping[str, Any], *,
                                expected_model: str | None = None,
                                masses_atol_ev: float = 1e-10) -> dict[str, Any]:
    """Recalculate singular masses from SAVED C5; no RGE or matching rerun.

    Checks both report mnu and C5 -> mnu against the original saved B/C
    conditional masses. Fails on modified physical-certification flags.
    """
    if not isinstance(report, Mapping) or report.get('format') != REPORT_FORMAT or \
            report.get('status') != REPORT_STATUS:
        raise HandoffContractError('Unexpected B/C conditional report schema/status')
    model = report.get('model')
    if model not in EXPECTED_BC_MASSES_EV or (expected_model is not None and model != expected_model):
        raise HandoffContractError('Unexpected B/C model in report')
    for name in ('finite_Y_threshold_delta_C_computed','finite_X_threshold_delta_C5_computed',
                 'full_sequential_matching_performed','physical_neutrino_mass_prediction_certified'):
        if report.get(name) is not False:
            raise HandoffContractError('Missing finite matching or conditional physics status changed: '+name)
    if report.get('finite_threshold_corrections_status') != 'BOTH_MISSING_NOT_PHYSICAL_ZERO' or \
            report.get('leading_log_sequential_candidate_completed') is not True:
        raise HandoffContractError('Finite scalar matching cannot be treated as zero/certified')
    scale_y = float(report['mu_Y_gev'])
    scale_x = float(report['mu_X_corrected_gev'])
    if not (np.isfinite(scale_y) and np.isfinite(scale_x) and scale_y > scale_x > 100):
        raise HandoffContractError('B/C Y->X->100GeV threshold ordering invalid')
    c5_x = complex_matrix_payload(report['C5_at_X_tree_gev_inverse'], shape=(3,3))
    validate_c5(c5_x)
    c5 = complex_matrix_payload(report['C5_at_100gev_gev_inverse'], shape=(3,3))
    m_report_ev = complex_matrix_payload(report['mnu_at_100gev_ev'], shape=(3,3))
    if not np.allclose(m_report_ev, m_report_ev.T, atol=1e-11, rtol=0):
        raise HandoffContractError('Reported mnu not complex symmetric')
    mass_check = check_neutrino_matrix(c5)
    expected_m = numerical_neutrino_mass_matrix(c5) * 1e9
    # The report stores mnu in the charged-lepton mass basis, and hence
    # a unitary congruence rather than equal matrix elements. Compare spectra.
    report_sv = np.linalg.svd(m_report_ev, compute_uv=False)
    values = np.asarray(report['takagi_masses_ev_descending'], dtype=float)
    reference = np.asarray(EXPECTED_BC_MASSES_EV[model], dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all() or \
            not np.all(values[:-1] >= values[1:]) or np.any(values < 0):
        raise HandoffContractError('Invalid saved Takagi mass spectrum')
    for name, got in (('reported mnu singulars',report_sv),
                      ('C5 Takagi singulars',mass_check.masses_ev_descending),
                      ('saved mass spectrum',values)):
        if not np.allclose(got, reference, rtol=0, atol=masses_atol_ev):
            raise HandoffContractError(f'{model}: {name} disagrees with saved reference')
    if not np.allclose(np.linalg.svd(expected_m, compute_uv=False), values,
                       rtol=0, atol=masses_atol_ev):
        raise HandoffContractError('C5 -> mnu conversion differs from saved convention')
    return {
        'model': model, 'status': 'ConditionalSpectrumRecomputedFromSavedC5',
        'mu_Y_gev': scale_y, 'mu_X_gev': scale_x,
        'masses_ev_descending': values.tolist(),
        'recomputed_c5_takagi_masses_ev_descending': mass_check.masses_ev_descending.tolist(),
        'c5_takagi_residual': mass_check.takagi_residual,
        'physical_neutrino_mass_prediction_certified': False,
        'finite_Y_X_matching_missing': True,
    }
