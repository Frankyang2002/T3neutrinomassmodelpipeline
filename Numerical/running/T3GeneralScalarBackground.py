"""Worker B: strict, source-injected EFT1 real-scalar background.

This module DOES NOT manufacture unverified SU(2) invariant tensors. For each
model the caller must supply the complete direct quartic basis (or its checked
affine columns) and all quadratic directions from the local source inventory.
It preserves the repository's 4 Higgs + 8 BSM real-coordinate convention and
keeps both S1/S2 until the joint scalar threshold.

V2 = 1/2 phi_a M2_ab phi_b; V4 = 1/4! Q_abcd phi_a phi_b phi_c phi_d.
Physical eigenfields eta satisfy phi = U eta, U.T M2 U = diag(eigenvalues).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Callable, Mapping
import numpy as np

MODEL_DIMS = {f'{f}_{a:+d}': dims for f, alphas, dims in (
    ('A', (-4, -2, 0), (1, 3)),
    ('B', (-3, -1, 1), (2, 2)),
    ('C', (-3, -1, 1), (2, 2)),
    ('D', (-2, 0, 2), (3, 1)),
    ('E', (-4, -2, 0, 2), (3, 3)),
) for a in alphas}


def _complex(x: Any) -> complex:
    if isinstance(x, Mapping):
        if set(x) != {'real', 'imag'}:
            raise ValueError('complex mapping needs exactly real,imag')
        x = complex(float(x['real']), float(x['imag']))
    z = complex(x)
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError('nonfinite coupling')
    return z


def _real_array(x: Any, shape: tuple[int, ...], label: str) -> np.ndarray:
    z = np.asarray(x, dtype=complex)
    if z.shape != shape or not np.isfinite(z).all():
        raise ValueError(f'{label}: expected finite shape {shape}')
    if np.max(np.abs(z.imag)) > 1e-11 * max(1.0, float(np.max(np.abs(z.real)))):
        raise ValueError(f'{label}: expected real tensor')
    return np.array(z.real, dtype=float, copy=True)


def _check_symmetric(q: np.ndarray, rank: int, tol: float = 1e-9) -> None:
    if rank == 2:
        checks = [q.T]
    elif rank == 4:
        checks = [q.transpose(1, 0, 2, 3), q.transpose(0, 2, 1, 3),
                  q.transpose(0, 1, 3, 2)]
    else:
        raise ValueError('unsupported rank')
    scale = max(float(np.linalg.norm(q)), 1e-100)
    if max(float(np.linalg.norm(q - x)) / scale for x in checks) > tol:
        raise ValueError('bosonic tensor symmetry violated')


def replace_source_hhss(q: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Replace the six H,H,S,S index placements; never accumulate or average."""
    q = _real_array(q, (12,)*4, 'Q')
    h = _real_array(h, (4, 4, 8, 8), 'HHSS')
    if not np.allclose(h, h.swapaxes(0, 1), rtol=1e-9, atol=1e-12) or not np.allclose(
            h, h.swapaxes(2, 3), rtol=1e-9, atol=1e-12):
        raise ValueError('HHSS exchange symmetry violated')
    q[:4,:4,4:,4:] = h
    q[:4,4:,:4,4:] = h.transpose(0,2,1,3)
    q[:4,4:,4:,:4] = h.transpose(0,2,3,1)
    q[4:,:4,:4,4:] = h.transpose(2,0,1,3)
    q[4:,:4,4:,:4] = h.transpose(2,0,3,1)
    q[4:,4:,:4,:4] = h.transpose(2,3,0,1)
    _check_symmetric(q, 4)
    return q


def _deltas(reference: Mapping[str, Any], state: Mapping[str, Any],
            directions: tuple[tuple[str, str], ...]) -> np.ndarray:
    out = []
    for name, axis in directions:
        if name not in reference or name not in state:
            raise ValueError(f'missing running coupling {name}')
        if axis not in ('real', 'imag'):
            raise ValueError(f'invalid coordinate {axis}')
        d = _complex(state[name]) - _complex(reference[name])
        out.append(d.real if axis == 'real' else d.imag)
    return np.asarray(out, float)


@dataclass(frozen=True)
class AffineTensor:
    """Validated cached basis; coefficients are differences from reference state."""
    reference: Mapping[str, Any]
    directions: tuple[tuple[str, str], ...]
    value_at_reference: np.ndarray
    columns: np.ndarray
    rank: int

    def __post_init__(self):
        shape = (12,) * self.rank
        if self.rank not in (2,4):
            raise ValueError('rank must be 2 or 4')
        base = _real_array(self.value_at_reference, shape, 'base')
        cols = _real_array(self.columns, (len(self.directions),) + shape, 'columns')
        _check_symmetric(base, self.rank)
        for col in cols: _check_symmetric(col, self.rank)
        for name, axis in self.directions:
            if name not in self.reference or axis not in ('real', 'imag'):
                raise ValueError('invalid response direction')
        if len(set(self.directions)) != len(self.directions):
            raise ValueError('duplicate response direction')
        object.__setattr__(self, 'value_at_reference', base)
        object.__setattr__(self, 'columns', cols)

    def evaluate(self, state: Mapping[str, Any]) -> np.ndarray:
        delta = _deltas(self.reference, state, self.directions)
        result = self.value_at_reference + np.tensordot(delta, self.columns, axes=(0,0))
        _check_symmetric(result, self.rank)
        return result


def construct_affine(reference: Mapping[str, Any], names: tuple[str, ...],
                     complex_names: frozenset[str], direct: Callable[[Mapping[str,Any]], Any],
                     rank: int, *, tolerance: float = 2e-9) -> AffineTensor:
    """Probe an explicitly complete direct evaluator; verify a multi-coupling point.

    A passing affine test establishes interpolation of supplied source only,
    NOT whether that source contains every invariant or correct normalization.
    """
    if len(names) != len(set(names)) or not set(complex_names) <= set(names):
        raise ValueError('invalid inventory')
    shape = (12,)*rank
    r = dict(reference)
    base = _real_array(direct(r), shape, 'direct reference')
    dirs, cols = [], []
    for name in names:
        if name not in r: raise ValueError(f'missing reference {name}')
        z = _complex(r[name])
        if name not in complex_names and abs(z.imag)>1e-12:
            raise ValueError(f'{name}: real coupling has imaginary part')
        for axis in (('real','imag') if name in complex_names else ('real',)):
            trial = dict(r)
            zz = z + (1j if axis == 'imag' else 1.)
            trial[name] = {'real':zz.real,'imag':zz.imag} if name in complex_names else zz.real
            cols.append(_real_array(direct(trial), shape, name) - base)
            dirs.append((name,axis))
    obj = AffineTensor(r, tuple(dirs), base, np.stack(cols) if cols else
                       np.zeros((0,)+shape), rank)
    trial = dict(r)
    for j,name in enumerate(names):
        z = _complex(r[name]) + 0.137*np.sin(j+1) + (0.083j*np.cos(j+1) if name in complex_names else 0)
        trial[name] = {'real':z.real,'imag':z.imag} if name in complex_names else z.real
    exact = _real_array(direct(trial), shape, 'multicoupling test')
    err = np.linalg.norm(exact - obj.evaluate(trial)) / max(np.linalg.norm(exact), 1e-100)
    if err > tolerance: raise ValueError(f'nonlinear/incomplete affine reconstruction {err:g}')
    return obj


@dataclass(frozen=True)
class ScalarBackground:
    model: str
    scale_gev: float
    q_fixed: np.ndarray
    m2_fixed: np.ndarray
    rotation_original_from_mass: np.ndarray
    eigenvalues_gev2: np.ndarray
    q_mass: np.ndarray
    assumptions: tuple[str, ...]

    @property
    def scalar_mass_matrix(self) -> np.ndarray:
        return self.m2_fixed[4:,4:]

    @property
    def scalar_rotation(self) -> np.ndarray:
        return self.rotation_original_from_mass[4:,4:]


def build_scalar_background(model: str, scale_gev: float, state: Mapping[str, Any],
                            quartic: AffineTensor, quadratic: AffineTensor, *,
                            hhss_source: np.ndarray | None = None,
                            expected_quartic_names: set[str] | None = None,
                            expected_quadratic_names: set[str] | None = None,
                            basis_rotation: np.ndarray | None = None,
                            certificate: str = 'candidate') -> ScalarBackground:
    """Build both scalar multiplets' joint 8x8 mass basis at the specified EFT1 scale.

    For B_-1/C_-1 require a supplied dedicated physical 8x8 basis transformation;
    its provenance must come from the LOCAL mixed-scalar implementation.
    `quadratic` must include all model-specific bilinears (including D_-2).
    """
    if model not in MODEL_DIMS: raise ValueError('unknown T3 model')
    if not np.isfinite(scale_gev) or scale_gev <= 0: raise ValueError('invalid scale')
    if quartic.rank != 4 or quadratic.rank != 2: raise ValueError('tensor ranks')
    qnames = {name for name,_ in quartic.directions}
    mnames = {name for name,_ in quadratic.directions}
    if expected_quartic_names is None or expected_quadratic_names is None:
        raise ValueError('explicit complete inventories are required for both tensors')
    if qnames != set(expected_quartic_names) or mnames != set(expected_quadratic_names):
        raise ValueError(f'inventory mismatch: quartic missing {set(expected_quartic_names)-qnames}, '
                         f'quadratic missing {set(expected_quadratic_names)-mnames}')
    q = quartic.evaluate(state)
    if hhss_source is not None: q = replace_source_hhss(q, hhss_source)
    m2 = quadratic.evaluate(state)
    if not np.allclose(m2[:4,4:],0,atol=1e-10,rtol=1e-10):
        raise ValueError('Higgs/BSM mixing incompatible with requested EFT1 basis')
    heavy = m2[4:,4:]
    # Orthogonal input transformation from source component basis to dedicated
    # mixed-model coordinates, followed by an exact diagonalization.
    if model in ('B_-1','C_-1') and basis_rotation is None:
        raise ValueError('B/C mixed scalar convention requires dedicated basis_rotation')
    v = np.eye(8) if basis_rotation is None else _real_array(basis_rotation,(8,8),'basis')
    if not np.allclose(v.T@v,np.eye(8),atol=1e-9): raise ValueError('nonorthogonal basis')
    vals, u = np.linalg.eigh(v.T@heavy@v)
    R = np.eye(12)
    R[4:,4:] = v@u
    if not np.allclose(R.T@m2@R, np.diag(np.diag(R.T@m2@R)), atol=1e-7,
                       rtol=1e-10):
        # Higgs block need not be diagonal in the unbroken-field convention.
        if not np.allclose(R[4:,4:].T@heavy@R[4:,4:],np.diag(vals),atol=1e-7,rtol=1e-10):
            raise ValueError('heavy scalar diagonalization failed')
    q_mass = np.einsum('abcd,aA,bB,cC,dD->ABCD',q,R,R,R,R,optimize=True)
    assumptions = (
        'All 8 BSM real scalar components retained and diagonalized jointly.',
        'Input EFT1 running state and complete invariant inventories are caller-supplied.',
        'Affine tests verify interpolation, not independent physical normalization.',
        'HHSS source replacement' if hhss_source is not None else 'HHSS from candidate input; not source-certified',
        'Mixed B/C rotation is caller-supplied' if model in ('B_-1','C_-1') else 'No dedicated B/C rotation required',
        'Negative eigenvalues are returned, not silently repaired.',
        f'Provenance status: {certificate}',
    )
    return ScalarBackground(model, float(scale_gev), q, m2, R, vals, q_mass, assumptions)
