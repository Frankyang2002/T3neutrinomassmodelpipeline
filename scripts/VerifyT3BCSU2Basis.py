"""Workstream B: independent SU(2) covariance/generator audit for B(-1)/C(-1).

Consumes the real Matchete generator arrays, v1 native Matchete CG,
RGBeta v5 native mass records, and the v6 conditional bilinear report.
A relative U(1) phase of the two scalar fields CANNOT be inferred from SU(2)
generator data. Hence no physical cross-tool phase certificate is produced.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts import VerifyT3BCMixedBilinearMap as previous

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("B_-1", "C_-1")
FORMAT = "t3-bc-live-matchete-su2-generators-v1"
EPS = np.array([[0, 1], [-1, 0]], dtype=complex)
PAULI = np.array([
    [[0, 1], [1, 0]],
    [[0, -1j], [1j, 0]],
    [[1, 0], [0, -1]],
], dtype=complex) / 2.0
IDENTITY = np.eye(2, dtype=complex)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError("B/C SU(2) basis audit: " + message)


def load(path: Path) -> dict:
    require(path.is_file(), "Missing input: " + str(path))
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(data, dict), "Expected JSON object: " + str(path))
    return data


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(gens_real: object, gens_imag: object) -> np.ndarray:
    try:
        real = np.asarray(gens_real, dtype=float)
        imag = np.asarray(gens_imag, dtype=float)
    except (ValueError, TypeError) as exc:
        raise ValueError("B/C SU(2) basis audit: generator components not numeric") from exc
    require(real.shape == (3, 2, 2) and imag.shape == real.shape,
            "Unexpected generator-array dimensions")
    require(np.all(np.isfinite(real)) and np.all(np.isfinite(imag)),
            "Nonfinite native generator component")
    return real + 1j * imag


def generators_check(gens: np.ndarray, *, tolerance: float = 1e-10) -> dict:
    """Compare actual doublet matrices with a canonical Pauli realisation.

    O changes the ADJOINT generator labels; it is not an inferred scalar
    doublet field rephasing or a unique cross-tool fundamental intertwiner.
    """
    gens = np.asarray(gens, dtype=complex)
    require(gens.shape == (3, 2, 2), "Wrong generator dimensions")
    require(np.all(np.isfinite(gens)), "Nonfinite generator")
    h = max(np.max(abs(g - g.conj().T)) for g in gens)
    traces = max(abs(np.trace(g)) for g in gens)
    trnorm = max(abs(np.trace(gens[a] @ gens[b]) - (0.5 if a == b else 0.0))
                 for a in range(3) for b in range(3))
    require(h < tolerance and traces < tolerance and trnorm < tolerance,
            "Non-Hermitian, non-traceless or incorrectly normalized generators")
    o_complex = np.array([[2 * np.trace(gens[a] @ PAULI[b])
                           for b in range(3)] for a in range(3)])
    o = o_complex.real
    o_imag = float(np.max(abs(o_complex.imag)))
    reconstruction = float(np.max(abs(gens - np.einsum("ab,bij->aij", o, PAULI))))
    orthogonality = float(np.max(abs(o @ o.T - np.eye(3))))
    det = float(np.linalg.det(o))
    require(o_imag < tolerance and reconstruction < tolerance and
            orthogonality < tolerance and abs(abs(det) - 1) < tolerance,
            "Generator basis does not admit an orthogonal adjoint map")
    # The handedness of the ordered native generator basis is measured, not
    # guessed. An orientation-reversed basis has det(O)=-1 and reversed f_abc.
    comm_residual = 0.0
    for a in range(3):
        for b in range(3):
            comm = gens[a] @ gens[b] - gens[b] @ gens[a]
            ref = 1j * np.sign(det) * sum(
                levi(a, b, c) * gens[c] for c in range(3))
            comm_residual = max(comm_residual, float(np.max(abs(comm - ref))))
    require(comm_residual < tolerance,
            "Native SU(2) commutators incompatible with adjoint rotation")
    return {
        "max_hermiticity_residual": float(h),
        "max_traceless_residual": float(traces),
        "max_trace_normalization_residual": float(trnorm),
        "max_adjoint_rotation_imaginary": o_imag,
        "max_generator_reconstruction_residual": reconstruction,
        "max_adjoint_orthogonality_residual": orthogonality,
        "max_su2_commutator_residual": comm_residual,
        "adjoint_rotation_matchete_to_canonical": o.tolist(),
        "adjoint_rotation_determinant": det,
    }


def levi(a: int, b: int, c: int) -> int:
    if len({a, b, c}) != 3:
        return 0
    triples = {(0, 1, 2), (1, 2, 0), (2, 0, 1)}
    return 1 if (a, b, c) in triples else -1


def invariant_check(gens: np.ndarray, cg: np.ndarray, *, trials: int = 32) -> dict:
    """Live native generator/CG compatibility and finite gauge transformations."""
    gens = np.asarray(gens, dtype=complex)
    cg = np.asarray(cg, dtype=complex)
    require(cg.shape == (2, 2), "Native mixed CG is not 2x2")
    k = cg[0, 1]
    require(abs(k) > 1e-12 and np.max(abs(cg - k * EPS)) < 1e-12,
            "Native bilinear CG is not a nonzero epsilon multiple")
    infinitesimal = float(max(np.max(abs(g.T @ cg + cg @ g)) for g in gens))
    require(infinitesimal < 1e-10, "Native CG is not SU(2) invariant")
    rng = np.random.default_rng(20261031)
    finite = 0.0
    action = 0.0
    phase_action = 0.0
    phase_shift_detectable = 0.0
    mu_m = 0.67 - 0.29j
    alpha1, alpha2 = 0.31, -0.47
    for _ in range(trials):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = float(rng.uniform(-3, 3))
        # exp(2 i angle n_a G_a), since (2n.G)^2=1 for fundamental SU2.
        u = np.cos(angle)*IDENTITY + 2j*np.sin(angle)*np.einsum("a,aij->ij", axis, gens)
        finite = max(finite, float(np.max(abs(u.T @ cg @ u - cg))),
                     float(np.max(abs(u.conj().T @ u - IDENTITY))),
                     float(abs(np.linalg.det(u) - 1)))
        s1 = rng.normal(size=2) + 1j*rng.normal(size=2)
        s2 = rng.normal(size=2) + 1j*rng.normal(size=2)
        m_poly = np.conjugate(s1) @ cg @ np.conjugate(s2)
        m_rot = np.conjugate(u @ s1) @ cg @ np.conjugate(u @ s2)
        action = max(action, abs(m_poly-m_rot))
        # Rephase both Matchete fields relative to fixed RGBeta component labels.
        m_s1 = np.exp(1j*alpha1)*s1
        m_s2 = np.exp(1j*alpha2)*s2
        lm = 2*np.real(mu_m * (np.conjugate(m_s1) @ cg @ np.conjugate(m_s2)))
        mu_r = - k * np.exp(-1j*(alpha1+alpha2)) * mu_m
        vr = 2*np.real(mu_r * (np.conjugate(s1) @ EPS @ np.conjugate(s2)))
        phase_action=max(phase_action,abs(lm+vr))
        wrong_mu_r = -k*mu_m
        wrong_vr = 2*np.real(wrong_mu_r * (np.conjugate(s1) @ EPS @ np.conjugate(s2)))
        phase_shift_detectable=max(phase_shift_detectable,abs(lm+wrong_vr))
    require(finite < 1e-10 and action < 1e-10 and phase_action < 1e-10,
            "Gauge covariance or phase-compensated action equality failed")
    require(phase_shift_detectable > 0.01,
            "Independent field rephasing negative control was insensitive")
    return {
        "matchete_over_canonical_epsilon": [float(k.real),float(k.imag)],
        "infinitesimal_epsilon_covariance_residual": infinitesimal,
        "max_finite_su2_covariance_residual": float(finite),
        "max_gauge_rotated_bilinear_residual": float(action),
        "max_rephased_cross_tool_action_residual": float(phase_action),
        "unrephased_coupling_negative_control": float(phase_shift_detectable),
        "field_phase_example_radians": [alpha1,alpha2],
        "canonical_ratio_if_same_component_phase": [-float(k.real),-float(k.imag)],
        "independent_scalar_u1_phase_fixed_by_generators": False,
    }


def check_generator_document(data: dict) -> dict:
    require(data.get("format") == FORMAT and
            data.get("physical_mapping_certified") is False and
            data.get("relative_field_phase_certified") is False and
            data.get("rgbeta_numeric_generators_extracted") is False,
            "Native generator report format or flags invalid")
    rows = data.get("models")
    require(isinstance(rows,dict) and set(rows) == set(MODELS),
            "Missing native B/C generator results")
    return rows


def audit(generator_report: dict, mat: dict, rgb: dict,
          bilinear_report: dict, *, root: Path = ROOT,
          live_paths: tuple[Path,Path] | None = None) -> dict:
    rows = check_generator_document(generator_report)
    verified = previous.audit(mat,rgb,root)
    require(bilinear_report.get("format") == verified["format"] and
            bilinear_report.get("status") == verified["status"] and
            bilinear_report.get("physical_mapping_certified") is False and
            bilinear_report.get("live_matchete_cg_tensor_checked") is True and
            bilinear_report.get("live_rgbeta_massive_lam_symbol_detected") is True,
            "Prior conditional map report is not valid")
    old = bilinear_report.get("models")
    require(isinstance(old,list) and len(old) == 2 and
            [x.get("model") for x in old] == list(MODELS),
            "Wrong prior model ordering")
    for a,b in zip(old,verified["models"]):
        require(a.get("conditional_mu_rgbeta_over_mu_matchete_exact") ==
                b["conditional_mu_rgbeta_over_mu_matchete_exact"] and
                a.get("native_matchete_registered_cg") == b["native_matchete_registered_cg"],
                "Prior conditional model result differs from native source")
    if live_paths is not None:
        saved = bilinear_report.get("live_evidence_sha256")
        require(isinstance(saved,dict) and
                sorted(saved.values()) == sorted(digest(p) for p in live_paths),
                "Native evidence changed after the prior bilinear report")
    result=[]
    for m,old_row in zip(MODELS,old):
        row=rows[m]
        require(isinstance(row,dict) and row.get("model") == m and
                row.get("status") == "LiveMatcheteSU2GeneratorsCaptured" and
                row.get("accessor") == "Matchete`Generators[SU[2], {1}]" and
                row.get("generator_dimensions") == [3,2,2] and
                row.get("relative_field_phase_certified") is False and
                row.get("physical_mapping_certified") is False,
                m + ": native generator probe not successful")
        gens=matrix_payload(row.get("generator_real"),row.get("generator_imag"))
        gencheck=generators_check(gens)
        cg=np.array([[complex(previous.parse_component(v)) for v in line]
                     for line in mat["models"][m]["cg_components_inputform"]],dtype=complex)
        invcheck=invariant_check(gens,cg)
        expected_ratio=-cg[0,1]
        recorded=complex(previous.parse_component(
            old_row["conditional_mu_rgbeta_over_mu_matchete_exact"]))
        require(abs(expected_ratio-recorded) < 1e-12,
                m + ": prior conditional ratio inconsistent with native generators/CG")
        result.append({
            "model":m,
            "native_matchete_su2_generators_checked":True,
            "matchete_generator_and_mixed_cg_covariance_checked":True,
            "rgbeta_canonical_pauli_realization_used":True,
            "rgbeta_native_epsilon_symbol_registered":True,
            "rgbeta_numeric_generator_basis_independently_extracted":False,
            "relative_scalar_field_phase_fixed":False,
            "conditional_mu_rgbeta_over_mu_matchete":old_row[
                "conditional_mu_rgbeta_over_mu_matchete_exact"],
            "generators":gencheck,
            "invariant":invcheck,
            "physical_mapping_certified":False,
        })
    return {
        "format":"t3-bc-independent-su2-native-generator-audit-v1",
        "status":"MatcheteGeneratorsAndCanonicalSU2CovarianceCheckedRelativeFieldPhaseUndetermined",
        "models":result,
        "why_not_physical_certificate":(
            "The RGBeta eps[SU2L@fund,...] is abstract and is realized here "
            "with canonical epsilon/Pauli components. The SU(2) algebra "
            "and invariant tensors do not determine independent U(1) "
            "rephasings of S1 and S2 between the two packages. Under "
            "S1_M=e^(i alpha1) S1_R, S2_M=e^(i alpha2) S2_R, "
            "mu_R=-k exp[-i(alpha1+alpha2)] mu_M."
        ),
        "physical_mapping_certified":False,
        "production_modified":False,
    }


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--generators",type=Path,required=True)
    p.add_argument("--matchete",type=Path,required=True)
    p.add_argument("--rgbeta",type=Path,required=True)
    p.add_argument("--bilinear",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--root",type=Path,default=ROOT)
    args=p.parse_args(argv)
    require(not args.output.exists(),"Refusing to overwrite previous SU2 output")
    require(args.output not in (args.generators,args.matchete,args.rgbeta,args.bilinear),
            "Output path conflicts with an input")
    paths=(args.generators,args.matchete,args.rgbeta,args.bilinear)
    payload=[load(x) for x in paths]
    result=audit(*payload,root=args.root,live_paths=(args.matchete,args.rgbeta))
    result["input_sha256"]={str(path):digest(path) for path in paths}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for x in result["models"]:
        print(x["model"],"native SU2 covariance = PASS; adjoint determinant =",
              round(x["generators"]["adjoint_rotation_determinant"],6))
    print(result["status"])
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
