"""Provisional A(-4) coupled master psi^2 phi^2 transport using *running* EFT1 tensors.

The ordinary A(-4) physical field content and 45-Weyl/12-real-scalar
renormalizable tensor snapshots are taken from Aminus4EFT1Couplings.  The
fermion-tree boundary remains FORMAL until its native Matchete GammaCC/Weyl
normalization and scalar-CG-to-RGBeta basis change are independently fixed.

No scalar-threshold matching, SMEFT running or neutrino mass is supplied.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping
import json
import math

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.running.T3EFT1MasterSourceBridge import (
    SourceBridgeError, embed_formal_source, load_formal_source,
)

LOOP = 1.0 / (16.0 * math.pi**2)


@dataclass
class Aminus4Trajectory:
    """Existing sampled *running* renormalizable EFT1 tensors."""
    t: np.ndarray
    mu: np.ndarray
    y: np.ndarray
    quartic: np.ndarray
    fermion_labels: tuple[str, ...]
    scalar_labels: tuple[str, ...]
    couplings: tuple[dict, ...]
    gauge: dict[str, tuple[tuple[np.ndarray, np.ndarray], ...]]

    @classmethod
    def load(cls, npz_path: str | Path, json_path: str | Path | None = None):
        p = Path(npz_path)
        doc = json.loads(Path(json_path or p.with_name('A_m4_eft1_tensor_basis.json'))
                         .read_text(encoding='utf-8-sig'))
        if doc.get('model') != 'A_-4' or doc.get('both_heavy_scalars_active_until_mu_s') is not True:
            raise SourceBridgeError('Expected A(-4) EFT1 with both heavy scalars present')
        if doc.get('matchete_scalar_and_Wilson_basis_map_certified') is not False:
            raise SourceBridgeError('Unexpected change to native basis certification: re-audit')
        from Numerical.running.Aminus4EFT1Couplings import F_LABELS, S_LABELS
        flabels = tuple(doc.get('fermion_basis', []))
        slabels = tuple(doc.get('real_scalar_basis', []))
        if flabels != F_LABELS or slabels != S_LABELS:
            raise SourceBridgeError('Running tensor point has incompatible label order')
        with np.load(p, allow_pickle=False) as data:
            mu = np.asarray(data['mu_gev'], dtype=float)
            Y = np.asarray(data['Y'], dtype=np.complex128)
            L = np.asarray(data['lambda4'], dtype=np.complex128)
            gauge = {}
            for kind, n in (('SU2', 3), ('U1', 1), ('SU3', 8)):
                f = np.asarray(data[f'generator_{kind}_fermion'], dtype=np.complex128)
                s = np.asarray(data[f'generator_{kind}_scalar'], dtype=np.complex128)
                if f.shape != (n, 45, 45) or s.shape != (n, 12, 12):
                    raise SourceBridgeError(f'Invalid {kind} generator dimensions')
                gauge[kind] = tuple((f[i], s[i]) for i in range(n))
        if mu.ndim != 1 or len(mu) < 3 or np.any(mu <= 0) or not np.isfinite(mu).all():
            raise SourceBridgeError('Invalid EFT1 sample scales')
        if Y.shape != (len(mu), 45, 45, 12) or L.shape != (len(mu), 12, 12, 12, 12):
            raise SourceBridgeError('Running Yukawa/quartic tensor shapes mismatch')
        if not np.isfinite(Y).all() or not np.isfinite(L).all():
            raise SourceBridgeError('Nonfinite EFT1 tensor')
        rows = doc.get('running_couplings')
        if not isinstance(rows, list) or len(rows) != len(mu):
            raise SourceBridgeError('Missing running EFT1 gauge coupling rows')
        order = np.argsort(np.log(mu))
        if np.any(np.diff(mu[order]) <= 0):
            raise SourceBridgeError('Duplicate EFT1 sample scales')
        for row in rows:
            for key in ('g2', 'gY', 'g3'):
                if key not in row or not np.isfinite(row[key]) or row[key] <= 0:
                    raise SourceBridgeError('Invalid running gauge coupling ' + key)
        return cls(np.log(mu[order]), mu[order], Y[order], L[order],
                   flabels, slabels, tuple(rows[int(i)] for i in order), gauge)

    def context(self, mu: float, mu_f: float, mu_s: float):
        """Interpolate full renormalizable tensors in t=ln(mu); no extrapolation."""
        from RGE.running.intermediate.GeneralEFT1C5Mixing import Context
        t = math.log(float(mu))
        eps = 1e-12
        if t < self.t[0] - eps or t > self.t[-1] + eps:
            raise SourceBridgeError('Wilson query lies outside EFT1 running interval')
        t = float(np.clip(t, self.t[0], self.t[-1]))
        hi = int(np.searchsorted(self.t, t, side='right'))
        if hi == 0:
            lo = hi = 0
        elif hi == len(self.t):
            lo = hi = len(self.t) - 1
        else:
            lo = hi - 1
        f = 0.0 if lo == hi else (t - self.t[lo]) / (self.t[hi] - self.t[lo])
        Y = (1-f)*self.y[lo] + f*self.y[hi]
        L = (1-f)*self.quartic[lo] + f*self.quartic[hi]
        coupling = {key: (1-f)*self.couplings[lo][key] + f*self.couplings[hi][key]
                    for key in ('gY','g2','g3')}
        sectors = tuple((coupling[name], self.gauge[kind]) for name, kind in
                        (('g2','SU2'),('gY','U1'),('g3','SU3')))
        fm = {x:i for i,x in enumerate(self.fermion_labels)}
        leptons = tuple((fm[f'L[{p},1]'], fm[f'L[{p},2]']) for p in (1,2,3))
        return Context('A_-4', 45, 12, leptons, (0,1,2,3),
                       mu_f, mu_s, Y, L, sectors, (), 'synthetic', {})



def integrate_master_tensor(initial: np.ndarray, context_at_mu, mu_f: float,
                            mu_s: float, *, rtol: float = 1e-6,
                            atol: float = 1e-18, max_step_log: float = 0.75
                            ) -> tuple[np.ndarray, dict]:
    """Integrate the FULL master tensor with running context, no fixed-order log.

    Generic low-level solver: inputs are not declared physical by this function.
    The tensor is transported all at once, so C_LLSS self-running and generated
    C_LLHH are included within the chosen psi^2 phi^2 operator class.
    """
    initial = np.asarray(initial, dtype=np.complex128)
    if initial.ndim != 4 or initial.shape[0] != initial.shape[1] or initial.shape[2] != initial.shape[3]:
        raise SourceBridgeError('Full master tensor must be (nf,nf,ns,ns)')
    if not np.isfinite(initial).all() or not 0 < mu_s < mu_f:
        raise SourceBridgeError('Nonfinite initial tensor or invalid threshold ordering')
    if not (0 < rtol < 1 and 0 < atol < 1 and max_step_log > 0):
        raise SourceBridgeError('Invalid integration tolerances/step')
    from RGE.running.intermediate.GeneralEFT1C5Mixing import calculate_master_beta
    shape = initial.shape
    calls = 0

    def rhs(t, packed):
        nonlocal calls
        calls += 1
        ctx = context_at_mu(math.exp(t))
        terms = calculate_master_beta(ctx, packed.reshape(shape))
        if len(terms) != 7:
            raise SourceBridgeError('Expected exactly seven master beta terms')
        val = LOOP * sum(terms.values())
        if val.shape != shape or not np.isfinite(val).all():
            raise SourceBridgeError('Malformed/nonfinite master beta')
        return val.reshape(-1)

    sol = solve_ivp(rhs, (math.log(mu_f), math.log(mu_s)), initial.ravel(),
                    method='DOP853', rtol=rtol, atol=atol,
                    max_step=max_step_log)
    if not sol.success or not np.isfinite(sol.y[:, -1]).all():
        raise SourceBridgeError(f'Wilson ODE failed: {sol.message}')
    return sol.y[:, -1].reshape(shape), {
        'solver':'DOP853_complex', 'rhs_calls':calls, 'nfev':int(sol.nfev),
        'relative_tolerance':rtol,'absolute_tolerance':atol,
        'one_loop_beta_terms':7, 'renormalizable_couplings_running':True,
    }


def transport(
    source_npz: str | Path,
    couplings_npz: str | Path,
    *, couplings_json: str | Path | None = None,
    heavy_source_from_rgbeta: np.ndarray | None = None,
    spinor_signed_factor: float | None = None,
    allow_formal_identity: bool = False,
    rtol: float = 1e-6, atol: float = 1e-18,
    max_step_log: float = 0.75,
) -> tuple[dict, np.ndarray]:
    """Integrate all master psi^2 phi^2 components to common scalar threshold.

    No implied physical normalization: all results are PROVISIONAL even when
    a map/factor is supplied, pending independent evidence and operator closure.
    """
    if heavy_source_from_rgbeta is None and not allow_formal_identity:
        raise SourceBridgeError('Explicit native-to-RGBeta scalar map required; '
                                'use allow_formal_identity ONLY for a diagnostic')
    if spinor_signed_factor is None and not allow_formal_identity:
        raise SourceBridgeError('Explicit GammaCC/Weyl signed factor required')
    if not (0 < rtol < 1 and 0 < atol < 1 and max_step_log > 0):
        raise SourceBridgeError('Invalid integration tolerances/step')
    src = load_formal_source(source_npz, model='A_-4')
    tr = Aminus4Trajectory.load(couplings_npz, couplings_json)
    mu_f = src.mu_f_gev
    mu_s = float(tr.mu[0])
    if abs(math.log(mu_f) - tr.t[-1]) > 1e-7 or not 0 < mu_s < mu_f:
        raise SourceBridgeError('F boundary and renormalizable EFT1 trajectory do not meet')
    initial, emb = embed_formal_source(
        src.raw, model='A_-4', fermion_labels=tr.fermion_labels,
        scalar_labels=tr.scalar_labels,
        heavy_source_from_target=heavy_source_from_rgbeta,
        spinor_signed_factor=spinor_signed_factor,
    )
    from RGE.running.intermediate.GeneralEFT1C5Mixing import project_c5
    shape = initial.shape
    final, solver = integrate_master_tensor(
        initial, lambda mu: tr.context(mu, mu_f, mu_s), mu_f, mu_s,
        rtol=rtol, atol=atol, max_step_log=max_step_log)
    fermion_sym = float(np.max(abs(final - final.swapaxes(0,1))))
    scalar_sym = float(np.max(abs(final - final.swapaxes(2,3))))
    ctx = tr.context(mu_s, mu_f, mu_s)
    c5, ward, flavour_sym = project_c5(final, ctx)
    output = {
        'status': 'ProvisionalFullMasterTensorODE_NotPhysicalPrediction',
        'model': 'A_-4', 'mu_F_gev': mu_f, 'mu_S_gev': mu_s,
        'hierarchy': 'UV -> all fermions out -> EFT1(S1,S2) -> [all scalars later] -> SMEFT',
        'success': True, **solver,
        'master_tensor_shape': list(shape),
        'formal_source_bridge': emb,
        'C5_LLHH_formal_gev_inverse': {'real': c5.real.tolist(), 'imag': c5.imag.tolist()},
        'weinberg_LLHH_absolute_ward_residual': float(ward),
        'flavour_symmetry_absolute_residual': float(flavour_sym),
        'master_fermion_pair_absolute_residual': fermion_sym,
        'master_scalar_pair_absolute_residual': scalar_sym,
        'gamma_cc_to_weyl_normalization_certified': False,
        'matchete_to_rgbeta_scalar_map_certified': False,
        'full_dimension_five_operator_closure_verified': False,
        'finite_joint_scalar_matching_included': False,
        'SMEFT_running_included': False,
        'physical_neutrino_mass_claimed': False,
        'leading_log_must_not_be_added_again': True,
        'warning': 'Formal transport only. Missing independent GammaCC/Weyl factor and native-to-RGBeta SU(2) CG/basis verification.',
    }
    return output, final
