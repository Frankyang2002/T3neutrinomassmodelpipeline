"""Batch 6: genuine EFT1 scalar-bubble MSbar pole counterterm and RG closure.

At the single F matching threshold BOTH heavy scalar multiplets are still
propagating. The dimension-five tree LLSS operator generates an LLHH 1PI
scalar bubble with the project's exact scalar/Weyl action normalization:

    A_soft_bare = -[(C_LLSS * lambda_HHSS)/(2*16*pi**2)]
                      * [1/epsbar + B0_finite(mu)].
    Delta C5_CT_EFT1 = +project_master(C_LLSS * lambda_HHSS)
                      /(2*16*pi**2*epsbar).

The CT above is derived from the explicit EFT1 diagram, NOT fitted to the
native Matchete residual. Its coefficient is source-resolved by heavy fermion
and scalar quartic. The full seven-channel candidate operator MSbar pole is
also reconstructed from the existing strict ADM, with an explicitly declared
beta/(2*16*pi**2) convention, for comparison only.

No UV native counterterm, Matchete GammaCC phase, or other hard topology is
silently inferred or cancelled. Matching-scale scans keep the saved couplings
fixed: diagnostic variations of a common SUBTRACTION scale, not new thresholds.
"""
from __future__ import annotations

from math import isfinite, log, pi
from typing import Mapping, Sequence

import numpy as np

LOOP = 16.0 * pi**2
SCHEMA = 't3-eft1-explicit-scalar-MSbar-CT-native-pole-scale-closure-batch6-v1'
SCENARIOS = ('analytic_triangle', 'native_formal', 'native_as_UV_minus_soft_hypothesis')
ADM_CHANNELS = (
    'yukawa_wavefunction', 'scalar_pair', 'mixed_yukawa_gauge',
    'crossed_yukawa', 'conjugate_coefficient_yukawa',
    'scalar_anomalous_dimension', 'fermion_anomalous_dimension')
SCALE_FACTORS = (0.5, 1.0, 2.0)


def relative(a: np.ndarray, b: np.ndarray) -> float:
    x, y = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(x-y) / max(float(np.linalg.norm(y)), 1e-100))


def symmetric3(data: np.ndarray, name: str) -> np.ndarray:
    x = np.asarray(data, complex)
    if (x.shape != (3, 3) or not np.isfinite(x).all() or
            relative(x, x.T) > 1e-9):
        raise ValueError(f'{name}: expected finite complex symmetric 3x3')
    return (x + x.T) / 2


def source_scalar_counterterm(c_by_fermion: np.ndarray,
                              quartic_channels: np.ndarray,
                              names: Sequence[str]) -> dict:
    """Compute the EFT1 LLHH divergent scalar bubble and its MSbar CT.

    The pole is mass-independent. All scalar basis rotations cancel between C
    and lambda; using the saved canonical real basis is therefore equivalent
    to summing scalar mass eigenstates. No epsilon pole is formally dropped.

    Returns [quartic_channel,heavy_fermion,flavour,flavour] pole coefficients
    multiplying 1/epsbar. Every full scalar index is contracted, no S1 or S2
    projection or scalar-eigenstate decoupling occurs at MF.
    """
    from Numerical.parallel_t3.threshold.matching import master_weinberg_template

    c = np.asarray(c_by_fermion, complex)
    q = np.asarray(quartic_channels, float)
    names = tuple(names)
    if (c.ndim != 7 or c.shape[:5] != (3, 3, 3, 2, 2) or
            c.shape[-1] != c.shape[-2] or c.shape[-1] not in (8, 12)):
        raise ValueError('Need exactly three fermions and both scalar multiplets')
    n = c.shape[-1]
    if (q.ndim != 5 or q.shape != (len(names), 4, 4, n, n) or
            len(names) == 0 or len(names) != len(set(names)) or
            any(not isinstance(x, str) or not x for x in names)):
        raise ValueError('Source quartic channels must cover the same full scalar basis')
    if not all(np.isfinite(z).all() for z in (c, q)):
        raise ValueError('Nonfinite scalar counterterm input')
    if (not np.allclose(c, c.transpose(0, 2, 1, 4, 3, 5, 6), rtol=1e-9, atol=1e-23)
            or not np.allclose(c, c.swapaxes(-1, -2), rtol=1e-9, atol=1e-23)):
        raise ValueError('Tree LLSS lepton/scalar pair exchange fails')
    if (not np.allclose(q, q.swapaxes(1, 2), rtol=1e-10, atol=1e-12)
            or not np.allclose(q, q.swapaxes(-1, -2), rtol=1e-10, atol=1e-12)):
        raise ValueError('Quartic HHSS pair exchange fails')
    t = master_weinberg_template()
    # [channel, F, flavour, flavour, Lweak, Lweak, Hreal, Hreal]
    amps = np.einsum('kfgijab,xABab->xkfgijAB', c, q, optimize=True)
    amp_scale = 1. / (2. * LOOP)
    c5_ct = amp_scale * np.einsum('xkfgijAB,ijAB->xkfg',
              amps, t.conj(), optimize=True) / np.vdot(t, t)
    if not np.isfinite(c5_ct).all():
        raise ValueError('Counterterm contraction nonfinite')
    c5 = c5_ct.sum(axis=(0, 1))
    if relative(c5, c5.T) > 1e-8:
        raise ValueError('Total C5 counterterm is not flavour symmetric')
    c5 = symmetric3(c5, 'total EFT1 scalar C5 CT')
    # Per-source projections may be non-singlets. Their sum, not a forced
    # projection of the full Wilson, is the supplied physical tensor test.
    residue = amps.sum(axis=(0, 1)) * amp_scale
    fitted = np.einsum('fg,ijAB->fgijAB', c5, t, optimize=True)
    pure = relative(residue, fitted)
    return {
        'CT_per_channel_per_heavy': c5_ct,
        'CT_per_channel': c5_ct.sum(axis=1),
        'CT_per_heavy': c5_ct.sum(axis=0),
        'CT_scalar_LLHH': c5,
        'EFT1_soft_bare_pole': -c5,
        'EFT1_MSbar_counterterm_pole': c5,
        'renormalized_soft_pole': np.zeros((3, 3), complex),
        'projected_tensor_purity_residual': pure,
        'channel_names': names,
        'scalar_components_active_at_MF': n,
        'fermions_removed_together_at_MF': 3,
        'counterterm_origin': 'actual_EFT1_LLSS_HHSS_scalar_bubble_B0_UV_pole',
        'sign_convention': 'EFT1 bubble pole -C*lambda/(2*16pi2); MSbar CT +same',
    }


def candidate_adm_counterterms(strict: Mapping[str, np.ndarray],
                               log_ms_over_mf: float) -> dict:
    """Compare EFT1 source scalar CT to one-loop *candidate* ADM residues.

    For the declared d=4-2epsilon convention, the residue of a one-loop
    operator MSbar counterterm from a beta channel is beta_C/(2*16pi**2).
    The saved strict evolution is Delta C = beta_C * ln(MS/MF)/(16pi**2),
    so the residue is Delta C/[2 ln(MS/MF)]. NOT a UV native CT.
    """
    L = float(log_ms_over_mf)
    if not isfinite(L) or L >= -1e-5:
        raise ValueError('Require actual descending EFT1 scale interval')
    pieces = {}
    for name in ADM_CHANNELS:
        key = 'C5_EFT1_strict_running_' + name
        if key not in strict:
            raise ValueError('Missing saved full-tensor ADM source channel ' + key)
        pieces[name] = symmetric3(np.asarray(strict[key], complex) / (2*L), key)
    combined = sum(pieces.values(), np.zeros((3, 3), complex))
    actual = symmetric3(strict['C5_EFT1_strict_running_delta'], 'running delta')
    if relative(combined, actual / (2*L)) > 1e-9:
        raise ValueError('Seven-channel candidate ADM residue reassembly failed')
    return {'CT_per_ADM_channel': pieces,
            'CT_candidate_total_7channel': symmetric3(combined, 'candidate ADM total CT'),
            'CT_candidate_scalar_pair': pieces['scalar_pair'],
            'counterterm_origin': 'candidate_ADM_beta_over_2epsilonbar_from_saved_Batch2',
            'native_UV_counterterm_calculated': False}


def native_pole_by_quartic(kernel_constant: np.ndarray,
                           quartic_labels: Sequence[str]) -> dict[str, np.ndarray]:
    """Do not infer a quartic map: retain original Matchete source labels."""
    arr = np.asarray(kernel_constant, complex)
    labels = tuple(quartic_labels)
    if (arr.ndim != 4 or arr.shape[1:] != (3, 3, 3) or
            len(labels) != len(arr) or not np.isfinite(arr).all()):
        raise ValueError('Native Laurent kernel/source metadata shape mismatch')
    out = {}
    for label, coeff in zip(labels, arr):
        if not isinstance(label, str) or not label.startswith('lambda'):
            raise ValueError('Native source quartic name absent')
        out.setdefault(label, np.zeros((3, 3), complex))
        out[label] += coeff.sum(axis=0)
    for k in out:
        out[k] = symmetric3(out[k], k)
    if relative(sum(out.values(), np.zeros((3, 3), complex)), arr.sum(axis=(0, 1))) > 1e-10:
        raise ValueError('Native Laurent grouping by original quartic omitted a source')
    return out


def matching_scale_scenarios(*,
    native_constant_by_kernel_fermion: np.ndarray,
    masses_gev: np.ndarray, c5_analytic_hard_mf: np.ndarray,
    c5_native_formal_mf: np.ndarray, c5_running_delta_mf_to_ms: np.ndarray,
    c5_scalar_bubble_ms: np.ndarray,
    soft_finite: Mapping[float, np.ndarray],
    mu_f: float = 1e5, mu_s: float = 1e3,
    factors: Sequence[float] = SCALE_FACTORS) -> dict:
    """Compose three transparent strict one-loop subtractions, not new thresholds.

    Match-scale variation holds the numerical UV couplings fixed. At mu=mu_F
    analytic_triangle reproduces Batch 2, native_formal reproduces Batch 4,
    and native_as_UV_minus_soft_hypothesis reproduces Batch 5. Only the
    source-triangle hard+soft scale relation is a derived loop identity.
    """
    k = np.asarray(native_constant_by_kernel_fermion, complex)
    mf, ms = float(mu_f), float(mu_s)
    m = np.asarray(masses_gev, float)
    factors = tuple(float(f) for f in factors)
    if (k.ndim != 4 or k.shape[1:] != (3, 3, 3) or m.shape != (3,) or
            np.any(m <= 0) or not np.isfinite(k).all() or not 0 < ms < mf or
            not factors or len(set(factors)) != len(factors) or 1. not in factors or
            any(not isfinite(f) or mf*f <= ms for f in factors)):
        raise ValueError('Matching-scale scan requires one MF, one MS, all three F')
    ah = symmetric3(c5_analytic_hard_mf, 'analytic MF hard')
    nh = symmetric3(c5_native_formal_mf, 'native formal MF hard')
    running = symmetric3(c5_running_delta_mf_to_ms, 'strict EFT1 running')
    bubble = symmetric3(c5_scalar_bubble_ms, 'MS scalar bubble')
    logarithms = np.log(mf / m)
    n0 = np.einsum('kfij,f->ij', k, (1. + 2*logarithms), optimize=True)
    if relative(n0, nh) > 1e-10:
        raise ValueError('Native MF finite source differs from Laurent coefficients')
    root_log = log(ms/mf)
    rows = {}
    for f in factors:
        sf = symmetric3(soft_finite[f], 'soft finite matching scale')
        dlog = log(f)
        n = symmetric3(nh + 2*dlog*k.sum(axis=(0,1)), 'native matching-scale log')
        # The analytic slope is fixed by the actual scalar-bubble MSbar CT
        # (passed through its identity with the analytic matching amplitude).
        a = symmetric3(ah + 2*dlog*analytic_ct_from_soft(soft_finite, factors),
                       'analytic matching-scale log')
        r = symmetric3(running * (log(ms/(mf*f))/root_log), 'EFT1 running at scan scale')
        choices = {
            'analytic_triangle': a,
            'native_formal': n,
            'native_as_UV_minus_soft_hypothesis': n - sf,
        }
        rows[f] = {name: {
            'hard_boundary_MF_renormalization_scale': symmetric3(z, name),
            'EFT1_running_strict_at_MS': r,
            'scalar_bubble_at_MS': bubble,
            'C5_SMEFT_MS': symmetric3(z+r+bubble, name + ' final MS'),
        } for name,z in choices.items()}
    return {'scales': rows, 'native_one_loop_slope': symmetric3(2*k.sum(axis=(0,1)), 'native slope'),
            'strict_EFT1_running_slope': symmetric3(running/root_log, 'EFT1 beta')}


def analytic_ct_from_soft(soft: Mapping[float, np.ndarray],
                          factors: Sequence[float]) -> np.ndarray:
    """Recover analytic hard log slope solely from the *evaluated* B0 finite.

    As B0bar(mu) changes by +2ln f, C5_soft(mu) changes by -2ln f*CT.
    Thus CT = -[soft(f)-soft(1)]/[2ln(f)], entirely independent of native.
    """
    ss = symmetric3(soft[1.], 'soft at MF')
    estimates = []
    for f in factors:
        if f == 1.: continue
        x = symmetric3(soft[f], 'soft matching-scale sample')
        estimates.append(-(x-ss)/(2*log(f)))
    if not estimates:
        raise ValueError('Need matching-scale samples away from MF')
    base = estimates[0]
    if any(relative(v, base)>1e-8 for v in estimates[1:]):
        raise ValueError('EFT1 finite scalar B0 scale slope inconsistent')
    return symmetric3(base, 'analytical B0-derived EFT1 CT')
