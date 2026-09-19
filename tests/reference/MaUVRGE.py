from __future__ import annotations

"""One-loop ultraviolet RGEs for the minimal Ma/scotogenic model.

Conventions follow Merle and Platscher, arXiv:1507.06314:

    D = 16 pi^2 d/d ln(mu)
    V = mH2 H^dag H + mEta2 eta^dag eta
        + lambda1/2 (H^dag H)^2
        + lambda2/2 (eta^dag eta)^2
        + lambda3 (H^dag H)(eta^dag eta)
        + lambda4 (H^dag eta)(eta^dag H)
        + lambda5/2 [(H^dag eta)^2 + h.c.]

The matrix convention is the one used in the paper: the heavy-neutrino
Yukawa matrix h has shape (n_heavy, 3), and the SM Yukawa matrices multiply
their left-handed flavour space from the right.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


LOOP = 16.0 * np.pi**2
MATRIX_SIZE = 3
SCALAR_COUNT = 10
MATRIX_NAMES = ("Ye", "Yu", "Yd", "h", "M")


def _matrix(value, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=complex)
    if result.shape != (MATRIX_SIZE, MATRIX_SIZE):
        raise ValueError(f"{name} must be a 3x3 matrix.")
    return result


def _real(value: complex, name: str, tolerance: float = 1e-9) -> float:
    value = complex(value)
    if abs(value.imag) > tolerance * max(1.0, abs(value.real)):
        raise ValueError(f"{name} unexpectedly acquired an imaginary part: {value}")
    return float(value.real)


def _dagger(matrix: np.ndarray) -> np.ndarray:
    return matrix.conj().T


@dataclass(frozen=True)
class MaUVInitialConditions:
    gY: float
    g2: float
    g3: float
    lambda1: float
    lambda2: float
    lambda3: float
    lambda4: float
    lambda5: float
    mH2: float
    mEta2: float
    Ye: np.ndarray
    Yu: np.ndarray
    Yd: np.ndarray
    h: np.ndarray
    M: np.ndarray

    def validated(self) -> "MaUVInitialConditions":
        M = _matrix(self.M, "M")
        if not np.allclose(M, M.T, rtol=1e-9, atol=1e-12):
            raise ValueError("The Majorana mass matrix M must be symmetric.")
        if float(self.mEta2) <= 0.0:
            raise ValueError("mEta2 must be positive in the inert phase.")

        return MaUVInitialConditions(
            gY=float(self.gY),
            g2=float(self.g2),
            g3=float(self.g3),
            lambda1=float(self.lambda1),
            lambda2=float(self.lambda2),
            lambda3=float(self.lambda3),
            lambda4=float(self.lambda4),
            lambda5=float(self.lambda5),
            mH2=float(self.mH2),
            mEta2=float(self.mEta2),
            Ye=_matrix(self.Ye, "Ye"),
            Yu=_matrix(self.Yu, "Yu"),
            Yd=_matrix(self.Yd, "Yd"),
            h=_matrix(self.h, "h"),
            M=M,
        )


@dataclass(frozen=True)
class MaUVRGEResult:
    mu_initial: float
    mu_final: float
    state: MaUVInitialConditions
    solver_success: bool
    solver_message: str
    nfev: int


@dataclass(frozen=True)
class MaUVTrajectoryResult:
    scales_gev: np.ndarray
    states: tuple[MaUVInitialConditions, ...]
    solver_success: bool
    solver_message: str
    nfev: int


def _pack(state: MaUVInitialConditions) -> np.ndarray:
    scalars = np.array(
        [
            state.gY,
            state.g2,
            state.g3,
            state.lambda1,
            state.lambda2,
            state.lambda3,
            state.lambda4,
            state.lambda5,
            state.mH2,
            state.mEta2,
        ],
        dtype=float,
    )
    matrices = [state.Ye, state.Yu, state.Yd, state.h, state.M]
    real_parts = np.concatenate([matrix.real.reshape(-1) for matrix in matrices])
    imag_parts = np.concatenate([matrix.imag.reshape(-1) for matrix in matrices])
    return np.concatenate([scalars, real_parts, imag_parts])


def _unpack(vector: np.ndarray) -> MaUVInitialConditions:
    scalars = vector[:SCALAR_COUNT]
    number_of_matrix_entries = len(MATRIX_NAMES) * MATRIX_SIZE * MATRIX_SIZE
    real_start = SCALAR_COUNT
    imag_start = real_start + number_of_matrix_entries
    real_data = vector[real_start:imag_start]
    imag_data = vector[imag_start:imag_start + number_of_matrix_entries]

    matrices: list[np.ndarray] = []
    block_size = MATRIX_SIZE * MATRIX_SIZE
    for index in range(len(MATRIX_NAMES)):
        start = index * block_size
        stop = start + block_size
        matrices.append(
            (real_data[start:stop] + 1j * imag_data[start:stop]).reshape(
                MATRIX_SIZE,
                MATRIX_SIZE,
            )
        )

    return MaUVInitialConditions(
        gY=scalars[0],
        g2=scalars[1],
        g3=scalars[2],
        lambda1=scalars[3],
        lambda2=scalars[4],
        lambda3=scalars[5],
        lambda4=scalars[6],
        lambda5=scalars[7],
        mH2=scalars[8],
        mEta2=scalars[9],
        Ye=matrices[0],
        Yu=matrices[1],
        Yd=matrices[2],
        h=matrices[3],
        M=matrices[4],
    )


def ma_uv_beta(state: MaUVInitialConditions) -> MaUVInitialConditions:
    """Return the one-loop beta coefficients before division by 16 pi^2."""

    s = state
    I3 = np.eye(MATRIX_SIZE, dtype=complex)

    Ye2 = _dagger(s.Ye) @ s.Ye
    Yu2 = _dagger(s.Yu) @ s.Yu
    Yd2 = _dagger(s.Yd) @ s.Yd
    h2 = _dagger(s.h) @ s.h
    hh_dagger = s.h @ _dagger(s.h)

    T = _real(np.trace(Ye2 + 3.0 * Yu2 + 3.0 * Yd2), "T")
    Tnu = _real(np.trace(h2), "Tnu")
    T4 = _real(
        np.trace(Ye2 @ Ye2 + 3.0 * Yu2 @ Yu2 + 3.0 * Yd2 @ Yd2),
        "T4",
    )
    T4nu = _real(np.trace(h2 @ h2), "T4nu")
    Tnue = _real(np.trace(h2 @ Ye2), "Tnue")

    beta_gY = 7.0 * s.gY**3
    beta_g2 = -3.0 * s.g2**3
    beta_g3 = -7.0 * s.g3**3

    beta_Yu = s.Yu @ (
        1.5 * (Yu2 - Yd2)
        + (T - 17.0 / 12.0 * s.gY**2 - 9.0 / 4.0 * s.g2**2 - 8.0 * s.g3**2) * I3
    )
    beta_Yd = s.Yd @ (
        1.5 * (Yd2 - Yu2)
        + (T - 5.0 / 12.0 * s.gY**2 - 9.0 / 4.0 * s.g2**2 - 8.0 * s.g3**2) * I3
    )
    beta_Ye = s.Ye @ (
        1.5 * Ye2
        + 0.5 * h2
        + (T - 15.0 / 4.0 * s.gY**2 - 9.0 / 4.0 * s.g2**2) * I3
    )
    beta_h = s.h @ (
        1.5 * h2
        + 0.5 * Ye2
        + (Tnu - 3.0 / 4.0 * s.gY**2 - 9.0 / 4.0 * s.g2**2) * I3
    )
    beta_M = hh_dagger @ s.M + s.M @ hh_dagger.conj()

    gauge_quartic_same = 0.75 * (
        s.gY**4 + 2.0 * s.gY**2 * s.g2**2 + 3.0 * s.g2**4
    )
    gauge_quartic_mixed = 0.75 * (
        s.gY**4 - 2.0 * s.gY**2 * s.g2**2 + 3.0 * s.g2**4
    )
    gauge_linear = s.gY**2 + 3.0 * s.g2**2

    beta_lambda1 = (
        12.0 * s.lambda1**2
        + 4.0 * s.lambda3**2
        + 4.0 * s.lambda3 * s.lambda4
        + 2.0 * s.lambda4**2
        + 2.0 * s.lambda5**2
        + gauge_quartic_same
        - 3.0 * s.lambda1 * gauge_linear
        + 4.0 * s.lambda1 * T
        - 4.0 * T4
    )
    beta_lambda2 = (
        12.0 * s.lambda2**2
        + 4.0 * s.lambda3**2
        + 4.0 * s.lambda3 * s.lambda4
        + 2.0 * s.lambda4**2
        + 2.0 * s.lambda5**2
        + gauge_quartic_same
        - 3.0 * s.lambda2 * gauge_linear
        + 4.0 * s.lambda2 * Tnu
        - 4.0 * T4nu
    )
    beta_lambda3 = (
        2.0 * (s.lambda1 + s.lambda2) * (3.0 * s.lambda3 + s.lambda4)
        + 4.0 * s.lambda3**2
        + 2.0 * s.lambda4**2
        + 2.0 * s.lambda5**2
        + gauge_quartic_mixed
        - 3.0 * s.lambda3 * gauge_linear
        + 2.0 * s.lambda3 * (T + Tnu)
        - 4.0 * Tnue
    )
    beta_lambda4 = (
        2.0 * (s.lambda1 + s.lambda2) * s.lambda4
        + 8.0 * s.lambda3 * s.lambda4
        + 4.0 * s.lambda4**2
        + 8.0 * s.lambda5**2
        + 3.0 * s.gY**2 * s.g2**2
        - 3.0 * s.lambda4 * gauge_linear
        + 2.0 * s.lambda4 * (T + Tnu)
        + 4.0 * Tnue
    )
    beta_lambda5 = s.lambda5 * (
        2.0 * (s.lambda1 + s.lambda2)
        + 8.0 * s.lambda3
        + 12.0 * s.lambda4
        - 3.0 * gauge_linear
        + 2.0 * (T + Tnu)
    )

    beta_mH2 = (
        6.0 * s.lambda1 * s.mH2
        + 2.0 * (2.0 * s.lambda3 + s.lambda4) * s.mEta2
        + s.mH2 * (2.0 * T - 1.5 * gauge_linear)
    )
    heavy_mass_trace = _real(
        # With U.T M U = D and h_mass = U.T h, the invariant form of
        # sum_i M_i^2 (h_mass h_mass^dagger)_ii is
        # Tr[M M^dagger h h^dagger].
        np.trace(s.M @ _dagger(s.M) @ hh_dagger),
        "heavy mass trace",
    )
    beta_mEta2 = (
        6.0 * s.lambda2 * s.mEta2
        + 2.0 * (2.0 * s.lambda3 + s.lambda4) * s.mH2
        + s.mEta2 * (2.0 * Tnu - 1.5 * gauge_linear)
        - 4.0 * heavy_mass_trace
    )

    return MaUVInitialConditions(
        gY=beta_gY,
        g2=beta_g2,
        g3=beta_g3,
        lambda1=beta_lambda1,
        lambda2=beta_lambda2,
        lambda3=beta_lambda3,
        lambda4=beta_lambda4,
        lambda5=beta_lambda5,
        mH2=beta_mH2,
        mEta2=beta_mEta2,
        Ye=beta_Ye,
        Yu=beta_Yu,
        Yd=beta_Yd,
        h=beta_h,
        M=beta_M,
    )


def _ode(_t: float, vector: np.ndarray) -> np.ndarray:
    return _pack(ma_uv_beta(_unpack(vector))) / LOOP


def evolve_ma_uv(
    initial: MaUVInitialConditions,
    mu_initial: float,
    mu_final: float,
    *,
    rtol: float = 1e-8,
    atol: float = 1e-11,
) -> MaUVRGEResult:
    """Run the full Ma theory between two scales where all fields are active."""

    if mu_initial <= 0.0 or mu_final <= 0.0:
        raise ValueError("RGE scales must be positive.")

    initial = initial.validated()
    solution = solve_ivp(
        _ode,
        (np.log(mu_initial), np.log(mu_final)),
        _pack(initial),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Ma UV RGE integration failed: {solution.message}")

    final = _unpack(solution.y[:, -1])
    final = MaUVInitialConditions(
        **{
            **final.__dict__,
            "M": 0.5 * (final.M + final.M.T),
        }
    ).validated()

    return MaUVRGEResult(
        mu_initial=float(mu_initial),
        mu_final=float(mu_final),
        state=final,
        solver_success=True,
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
    )


def sample_ma_uv_trajectory(
    initial: MaUVInitialConditions,
    mu_initial: float,
    mu_final: float,
    *,
    samples: int = 64,
    rtol: float = 1e-8,
    atol: float = 1e-11,
) -> MaUVTrajectoryResult:
    """Return logarithmically spaced states along one UV RGE integration."""

    if mu_initial <= 0.0 or mu_final <= 0.0:
        raise ValueError("RGE scales must be positive.")
    if samples < 2:
        raise ValueError("At least two trajectory samples are required.")

    initial = initial.validated()
    scales = np.geomspace(mu_initial, mu_final, num=samples)
    times = np.log(scales)
    solution = solve_ivp(
        _ode,
        (np.log(mu_initial), np.log(mu_final)),
        _pack(initial),
        t_eval=times,
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Ma UV trajectory integration failed: {solution.message}")

    states: list[MaUVInitialConditions] = []
    for vector in solution.y.T:
        state = _unpack(vector)
        states.append(
            MaUVInitialConditions(
                **{
                    **state.__dict__,
                    "M": 0.5 * (state.M + state.M.T),
                }
            )
        )
    return MaUVTrajectoryResult(
        scales_gev=scales,
        states=tuple(states),
        solver_success=True,
        solver_message=str(solution.message),
        nfev=int(solution.nfev),
    )


def takagi_majorana(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return U,masses satisfying U.T @ M @ U = diag(masses)."""

    M = _matrix(M, "M")
    if not np.allclose(M, M.T, rtol=1e-8, atol=1e-11):
        raise ValueError("Takagi input must be complex symmetric.")

    eigenvalues, U = np.linalg.eigh(_dagger(M) @ M)
    order = np.argsort(np.maximum(eigenvalues, 0.0))
    U = U[:, order]
    diagonal = np.diag(U.T @ M @ U)
    phases = np.exp(-0.5j * np.angle(diagonal))
    U = U @ np.diag(phases)
    transformed = U.T @ M @ U
    masses = np.real(np.diag(transformed))

    scale = max(1.0, float(np.linalg.norm(M)))
    if np.linalg.norm(transformed - np.diag(masses)) > 1e-7 * scale:
        raise RuntimeError("Takagi factorisation residual is too large.")
    if np.any(masses <= 0.0):
        raise ValueError("All heavy Majorana masses must be positive.")

    return U, masses


def _loop_f_over_mass(majorana_mass: float, inert_mass_sq: float) -> float:
    """Return f(M,m_eta)/M from Eq. (6) of arXiv:1507.06314."""

    M = float(majorana_mass)
    r = float(inert_mass_sq) / M**2
    epsilon = r - 1.0
    if abs(epsilon) < 1e-5:
        f = -0.5 + epsilon / 3.0 - epsilon**2 / 4.0
    else:
        f = np.log(r) / (1.0 - r) ** 2 + 1.0 / (1.0 - r)
    return float(f / M)


def match_scotogenic_c5(state: MaUVInitialConditions) -> dict:
    """Match the small-lambda5 Ma model onto the SM Weinberg coefficient.

    The returned matrix obeys the project convention m_nu = -v^2 C5.
    ``lambda5`` and ``h`` are the paper's scalar-potential and Yukawa
    couplings.  The bridge to the existing Matchete T3-B convention is
    ``lambdaT3 = -lambda5`` and ``y1 = y2 = conjugate(h)``.  The inert mass
    argument is the unbroken-phase parameter ``mEta2``; electroweak-suppressed
    corrections to the physical neutral-scalar average are not included.
    """

    state = state.validated()
    U, masses = takagi_majorana(state.M)
    h_mass = U.T @ state.h
    weights = np.array(
        [_loop_f_over_mass(mass, state.mEta2) for mass in masses],
        dtype=float,
    )
    C5 = state.lambda5 / LOOP * h_mass.T @ np.diag(weights) @ h_mass
    C5 = 0.5 * (C5 + C5.T)
    return {
        "C5": C5,
        "heavy_masses": masses,
        "h_mass_basis": h_mass,
        "loop_weights": weights,
        "takagi_U": U,
    }
