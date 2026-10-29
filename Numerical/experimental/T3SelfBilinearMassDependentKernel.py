"""Mass-sensitive *one-scale* T3 A/D/E self-bilinear C5 spectral continuation.

Use only a mass-dependent source kernel. A hierarchical Wilson coefficient with
no MS1/MS2 dependence CANNOT be split by numerical mass substitution.

This is a conditional one-scale calculation, not sequential EFT matching.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import math
import numpy as np

from Numerical.experimental.T3SelfBilinearAnomalousKernel import (
    MODELS, anomalous_from_unit_kernels, model_alpha,
)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 mass-dependent self-bilinear C5: ' + message)


def source_path(root: Path, model: str) -> Path:
    require(model in MODELS, 'unsupported model')
    group, alpha = model.split('_')
    tag = ('m' if alpha.startswith('-') else 'p') + alpha[1:]
    return Path(root) / f'output/hypercharge/T3_{group}_alpha_{tag}/data/c5_coefficient.txt'


def verify_legacy_source(path: Path) -> dict:
    """Check exact structural mass dependence before parsing with project code."""
    path = Path(path)
    require(path.is_file(), 'missing mass-dependent Matchete source: ' + str(path))
    raw = path.read_text(encoding='utf-8-sig')
    require(bool(raw.strip()) and len(raw) > 120, 'empty/unexpected one-scale C5 source')
    required = ('Coupling[MF,', 'Coupling[MS1,', 'Coupling[MS2,',
                'Coupling[lambdaT3,', 'Coupling[y1,', 'Coupling[y2,')
    require(all(s in raw for s in required),
            'source does not contain both scalar masses and the T3 topology')
    require(not any(s in raw for s in ('$Failed', 'Missing[', 'BetaTerm[')),
            'source contains unevaluated/failure expression')
    return {'status':'ExplicitMassDependentLegacyOneScaleSource',
            'scalar_mass_symbols':['MS1','MS2'], 'source_length':len(raw)}


def checked_c5(value: np.ndarray, label: str) -> np.ndarray:
    arr = np.asarray(value, dtype=complex)
    require(arr.shape == (3, 3) and np.isfinite(arr.real).all()
            and np.isfinite(arr.imag).all(), label + ': nonfinite 3x3 C5')
    require(np.allclose(arr, arr.T, rtol=1e-10, atol=1e-25),
            label + ': nonsymmetric C5')
    return 0.5 * (arr + arr.T)


def split_unit_kernels(kernel, uv, eft, basis) -> tuple[np.ndarray, np.ndarray, float]:
    """Vary the actual self-scalar mass in an independently mass-sensitive C5.

    λT3 is set to 1 for unit kernels, to avoid dividing by small fitted λT3.
    The result must change with m_±² if the split is nonzero.
    """
    require(basis.species in ('S1', 'S2'), 'unsupported self-scalar species')
    sq = tuple(float(x) for x in basis.eigen_masses_squared)
    require(len(sq) == 2 and all(math.isfinite(x) and x > 0 for x in sq),
            'invalid split masses')
    key = 'mS1Sq' if basis.species == 'S1' else 'mS2Sq'
    unit = []
    for mass_sq in sq:
        candidate = replace(eft, **{key: mass_sq, 'lambdaT3': 1+0j}).validated()
        unit.append(checked_c5(kernel.evaluate(uv, candidate), 'unit one-scale C5'))
    plus, minus = unit
    scale = max(float(np.linalg.norm(plus)), float(np.linalg.norm(minus)), 1e-300)
    sensitivity = float(np.linalg.norm(plus - minus) / scale)
    require(math.isfinite(sensitivity), 'nonfinite scalar-mass response')
    if not math.isclose(sq[0], sq[1], rel_tol=0, abs_tol=0):
        require(sensitivity > 1e-12,
                'one-scale kernel is insensitive to the split scalar mass; '
                'refuse spurious zero anomalous matching')
    # Verify λT3-linear topology against direct nonunit numerical evaluation.
    lam = complex(eft.lambdaT3)
    check = checked_c5(kernel.evaluate(uv, eft), 'unsplit one-scale C5')
    center = replace(eft, lambdaT3=1+0j).validated()
    base = checked_c5(kernel.evaluate(uv, center), 'unsplit unit C5')
    parity = float(np.linalg.norm(check-lam*base) /
                   max(float(np.linalg.norm(check)), float(np.linalg.norm(lam*base)), 1e-300))
    require(parity < 3e-8, 'legacy source is not linear in lambdaT3')
    return plus, minus, sensitivity


def spectral_c5(plus: np.ndarray, minus: np.ndarray, *,
                lam_t3: complex, mu_sq: complex, alpha: int,
                extra_matchete: complex, tensor_factor: complex) -> dict:
    plus = checked_c5(plus, 'plus kernel')
    minus = checked_c5(minus, 'minus kernel')
    lam = complex(lam_t3)
    require(np.isfinite([lam.real, lam.imag]).all(), 'invalid topology coupling')
    normal = checked_c5(0.5 * lam * (plus+minus), 'normal channel')
    anomalous = checked_c5(anomalous_from_unit_kernels(
        plus, minus, mu_sq=mu_sq, alpha=alpha,
        matchete_extra=extra_matchete, tensor_factor=tensor_factor), 'anomalous channel')
    combined = checked_c5(normal+anomalous, 'combined channel')
    denominator = float(np.linalg.norm(normal))
    return {'normal': normal, 'anomalous': anomalous, 'combined': combined,
            'anomalous_over_normal': (float(np.linalg.norm(anomalous)/denominator)
                                      if denominator > 0 else None)}
