"""T3 A_-4 Batch 1: source-pinned restricted benchmark and fail-closed spinor audit.

This is NOT a new matching solver. It reuses the original frozen benchmark,
provisional source-resolved Matchete F-tree boundary, and existing spinor
validator; it cannot certify Matchete GammaCC from abstract identities.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

MODEL = "A_-4"
REPO_COMMIT = "e581a434c4114ea1b100c50378c085529fb85be7"
INPUT_REL = "Numerical/lo_certification/benchmark/A_m4_input.json"
OWNED_OUTPUT = "output/general_observable/lo_certification_benchmark"

# Git blob hashes from a read-only inspection of the main branch on 2026-11-03.
# A mismatch is diagnosed, not silently accepted as current source truth.
GIT_BLOBS = {
    "Lagrangian/model/T3ModelCatalog.wl": "9a6c47e10fc83d9480f08353682151180c9fbb4a",
    "Lagrangian/interactions/T3Topology.wl": "a45a099a8b766cbad0114a0559cb0bcb80895649",
    "Numerical/lo_certification/spinor/spinors.py": "52575a0cd1b90a2f1f3db01587c65ac730b5b2b7",
    "Numerical/lo_certification/spinor/boundary.py": "9ebc6eca3b0713341d31405d4c18374fa602ec33",
    "Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_m4.json": "431f38262d06042877d009932271dcfde074c73f",
    "output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f/A_m4/boundary.json": "bf5c08d2a2d7563167c58605cca3df4ffc1686e3",
    "output/general_observable/lo_certification_couplings/ws2_extra_source_v1/A_m4/V6_source_evidence.json": "c447e652c7f6f56ff81e45f96a96e5bde9871680",
}


class AuditError(RuntimeError):
    """Missing, stale, mismatched or physically inconsistent benchmark source."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def read_json(path: Path) -> dict:
    require(path.is_file(), f"Missing required local file: {path}")
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj, dict), f"Expected JSON object: {path}")
    return obj


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_blob_sha(path: Path) -> str | None:
    if not path.is_file():
        return None
    raw = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def drift_audit(root: Path) -> dict[str, dict[str, Any]]:
    out = {}
    for rel, expected in GIT_BLOBS.items():
        observed = git_blob_sha(root / rel)
        out[rel] = {
            "github_blob_sha": expected, "local_blob_sha": observed,
            "match": observed == expected,
            "status": "MATCH" if observed == expected else
                      ("MISSING" if observed is None else "SOURCE_DRIFT"),
        }
    return out


def validate_matrix(data: Any, name: str) -> None:
    require(isinstance(data, list) and len(data) == 3, f"{name}: expected 3x3")
    for row in data:
        require(isinstance(row, list) and len(row) == 3, f"{name}: expected 3x3")
        require(all(isinstance(x, (int, float)) and math.isfinite(x) for x in row),
                f"{name}: non-finite/non-numeric")


def compare_frozen(contract: dict, frozen: dict) -> dict:
    """Require exact numeric identity with the prior frozen source, not a new fit."""
    require(contract.get("format") == "t3-restricted-a-m4-batch1-input-v1", "Wrong contract format")
    require(contract.get("model") == MODEL, "Model mismatch")
    r = contract["representation"]
    require((r["scalar1"]["su2_dimension"], r["scalar2"]["su2_dimension"],
             r["heavy_fermion"]["su2_dimension"]) == (1, 3, 2), "Incorrect A representation")
    require((r["scalar1"]["hypercharge"], r["scalar2"]["hypercharge"],
             r["heavy_fermion"]["hypercharge"]) == (-2.0, -1.0, -1.5),
            "Incorrect alpha=-4 hypercharge assignment")
    require(r["heavy_fermion"]["mass_type"] == "Dirac", "A_-4 is a complex F doublet")
    fr = frozen["representation"]
    require([fr[k] for k in ("d_s1", "d_s2", "d_f", "alpha")] == [1, 3, 2, -4],
            "Frozen source representation differs")
    scales = contract["matching"]
    sc = frozen["scales"]
    require(sc["mu_fermion_threshold_gev"] == scales["mu_F_gev"] and
            sc["mu_scalar_threshold_gev"] == scales["mu_S_gev"] and
            sc["mu_low_gev"] == scales["mu_low_gev"] and
            frozen["base_state"]["mu_gev"] == scales["mu_UV_gev"],
            "Frozen scale mismatch")
    require(scales["mu_UV_gev"] > scales["mu_F_gev"] > scales["mu_S_gev"] > scales["mu_low_gev"] > 0,
            "Invalid two-threshold ordering")
    saved = frozen["base_state"]["ordinary"]
    v = contract["frozen_values_at_mu_UV"]
    for name in ("y1", "y2"):
        real = v[name + "_real"]
        validate_matrix(real, name)
        require(real == saved[name]["real"], f"{name} real matrix differs from frozen source")
        require(v[name + "_imag"] == "all exactly zero in source" and
                all(x == 0 for row in saved[name]["imag"] for x in row),
                f"{name} imaginary components differ from source")
    mf = saved["MF"]
    validate_matrix(mf["real"], "MF")
    masses = v["M_F_diagonal_gev"]
    require(masses == [mf["real"][i][i] for i in range(3)], "Fermion masses differ")
    require(all(mf["real"][i][j] == 0 for i in range(3) for j in range(3) if i != j) and
            all(x == 0 for row in mf["imag"] for x in row), "Expected real diagonal MF")
    for name, src in (("m_S1_squared_gev2", "mS1Sq"),
                      ("m_S2_squared_gev2", "mS2Sq"),
                      ("lambda_S1", "lambdaS1"), ("lambda_S2", "lambdaS2"),
                      ("lambda_H1", "lambdaH1"), ("lambda_H2", "lambdaH2"),
                      ("lambda_12", "lambda12"),
                      ("lambda_H2Adj", "lambdaH2Adj"),
                      ("lambda_S2Adj", "lambdaS2Adj")):
        require(v[name] == saved[src], f"Frozen numerical parameter differs: {name}")
    require(v["lambda_T3"] == saved["lambdaT3"], "lambdaT3 differs")
    require(v["vev_gev"] == frozen["numerical"]["vev_gev"], "VEV differs")
    require(v["lambda_H2Adj"] == v["lambda_S2Adj"] == 0,
            "Restricted baseline zero-coupling contract changed")
    m_s = [math.sqrt(v["m_S1_squared_gev2"]), math.sqrt(v["m_S2_squared_gev2"])]
    require(all(math.isfinite(x) and x > 0 for x in m_s), "Invalid positive scalar masses")
    require(min(masses) > max(m_s), "Fermion not heavier than BOTH scalars")
    require(max(abs(x) for n in ("y1_real", "y2_real") for row in v[n] for x in row) < 1,
            "Nonperturbative Yukawa example")
    require(abs(complex(**v["lambda_T3"])) < 4 * math.pi,
            "Nonperturbative quartic example")
    return {
        "values_are_from_frozen_UV_source": True,
        "yukawa_input_scale_gev": scales["mu_UV_gev"],
        "couplings_run_to_mu_F": False,
        "fermion_mass_eigenvalues_gev": masses,
        "scalar_unbroken_mass_parameters_gev": m_s,
        "mF_min_over_mS_max": min(masses) / max(m_s),
        "mS_max_over_mF_min": max(m_s) / min(masses),
        "scalar_positive_unbroken_quadratic_parameters": True,
        "scalar_broken_phase_masses_certified": False,
        "lambda_T3_source_value": v["lambda_T3"],
        "y1_real_at_UV": v["y1_real"],
        "y2_real_at_UV": v["y2_real"],
        "omitted_couplings_zero_at": "mu_UV only; RG stability unverified",
    }


def signed_ambiguity() -> dict:
    """Show why abstract C^2 and antisymmetry cannot determine the sign."""
    import numpy as np
    z = np.zeros((2, 2), dtype=complex)
    eye = np.eye(2, dtype=complex)
    sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    gamma0 = np.block([[z, eye], [eye, z]])
    gamma2 = np.block([[z, sigma2], [-sigma2, z]])
    gamma5 = np.block([[-eye, z], [z, eye]])
    pl = (np.eye(4) - gamma5) / 2
    c = 1j * gamma2 @ gamma0
    cinv = np.linalg.inv(c)
    eps = 1j * sigma2
    a = np.r_[np.array([1 + 2j, 3 - 1j]), [0j, 0j]]
    b = np.r_[np.array([2 - 1j, -1 + 4j]), [0j, 0j]]
    w = a[:2].T @ eps @ b[:2]
    require(abs(w) > 0, "Degenerate spinor witness")
    results = []
    for label, obj, factor in (("C", c, 1), ("C_inverse", cinv, -1)):
        require(np.allclose(obj @ obj, -np.eye(4)) and
                np.allclose(obj.T, -obj), "Abstract GammaCC algebra failed")
        observed = a.T @ obj @ pl @ b
        require(np.allclose(observed, factor*w), "Signed spinor ratio mismatch")
        results.append({"candidate": label, "PL_over_Weyl_epsilon": factor,
                        "GammaCC_square_minus_I": True, "GammaCC_antisymmetric": True})
    return {
        "two_indistinguishable_abstract_choices": results,
        "abstract_GammaCC_identity_sufficient_to_fix_sign": False,
        "native_GammaCC_to_master_Weyl_factor": None,
        "required_resolution": "a native, signed, connected four-point UV/EFT1 amplitude and the master Weyl-operator prefactor with matching external-leg order",
    }


def audit(root: Path) -> dict:
    root = root.resolve()
    contract = read_json(root / INPUT_REL)
    src = contract["source"]
    frozen_path = root / src["frozen_benchmark_path"]
    frozen = read_json(frozen_path)
    require(sha256_file(frozen_path) == src["frozen_benchmark_sha256"],
            "Frozen source SHA256 differs; do not silently change the benchmark")
    numeric = compare_frozen(contract, frozen)
    v6 = read_json(root / src["v6_extra_source_path"])
    require(v6["model"] == MODEL and v6["extra_HHSS_channels"] == 0 and
            v6["RGBeta_additional_quartics"] == [],
            "A_-4 no-extra-HHSS assumption no longer supported by V6 source")
    bound_path = root / src["provisional_boundary_path"]
    bound = read_json(bound_path)
    # Reuse the existing audited component-tensor implementation; no parallel solver.
    from Numerical.lo_certification.spinor.boundary import validate_model, verify_reference_hashes
    validated = validate_model(bound, MODEL)
    verify_reference_hashes(root, MODEL, bound)
    require(bound["mu_F_gev"] == contract["matching"]["mu_F_gev"],
            "Boundary evaluated at different fermion threshold")
    require(bound["spinor_normalization"]["gamma_cc_to_master_weyl_factor"] is None,
            "Unexpected physical-phase promotion in provisional release")
    rows = bound["summed_independent_symmetric_coordinates"]
    example = next((r for r in rows if r["coordinate"] ==
                   contract["source_verification"]["sample_coordinate"]), None)
    require(example is not None, "Pinned sample LLSS coordinate disappeared")
    expected = contract["source_verification"]["sample_formal_coefficient_gev_inv"]
    require(math.isclose(example["coefficient_real"], expected["real"], rel_tol=1e-14) and
            math.isclose(example["coefficient_imag"], expected["imag"], abs_tol=1e-30),
            "Pinned provisional LLSS sample coefficient changed")
    from Numerical.lo_certification.spinor.spinors import convention_probe, exchange_probe
    old_matrix = convention_probe()
    old_exchange = exchange_probe()
    require(old_matrix["exact_matrices_and_complex_spinor_checks_passed"] and
            old_exchange["exact_grassmann_saddle_signs_verified"],
            "Existing conditional Weyl checks failed")
    ambigu = signed_ambiguity()
    drift = drift_audit(root)
    exact = all(d["match"] for d in drift.values())
    # No value of the physical phase is deduced from gamma matrix identities.
    sample = complex(example["coefficient_real"], example["coefficient_imag"])
    return {
        "format": "t3-restricted-batch1-a-m4-audit-v1",
        "model": MODEL,
        "github_read_only_commit_inspected": REPO_COMMIT,
        "local_matches_inspected_github_blobs": exact,
        "source_blobs": drift,
        "source_validation": {
            "frozen_benchmark_sha256": sha256_file(frozen_path),
            "provisional_boundary_sha256": sha256_file(bound_path),
            "V6_extra_HHSS_channels": 0,
            "A_m4_native_PL_sources": validated["verified_source_count_PL"],
            "A_m4_native_PR_conjugates_counted": 0,
            "A_m4_independent_LLSS_rows": validated["summed_coordinate_rows"],
            "A_m4_max_source_sum_relative_residual": validated["max_relative_source_sum_residual"],
            "reference_source_SHA256_checks_passed": True,
        },
        "numeric_input": numeric,
        "fermion_boundary": {
            "operator": "LLSS", "dimension_gev": -1,
            "source_numeric_status": "PROVISIONAL_FORMAL_POLYNOMIAL_ONLY",
            "sample_coordinate": example["coordinate"],
            "sample_formal_coefficient_gev_inv": expected,
            "conditional_if_GammaCC_equals_C": {"real": sample.real, "imag": sample.imag},
            "conditional_if_GammaCC_equals_C_inverse": {"real": -sample.real, "imag": -sample.imag},
            "numeric_physical_C_LLSS": None,
            "numeric_physical_C5_gev_inv": None,
            "neutrino_mass_ev": None,
        },
        "spinor_algebra": {
            "existing_spinor_checks_passed": True,
            "existing_Grassmann_checks_passed": True,
            "declared_Dirac_cross_channel_factor": old_exchange["dirac_mixed_channel_factor"],
            "native_ambiguity": ambigu,
        },
        "matching": {
            "order": "UV -> integrate out ALL heavy fermions at mu_F -> EFT1 with BOTH scalars -> integrate out BOTH scalars at mu_S -> SMEFT",
            "mu_F_gev": contract["matching"]["mu_F_gev"],
            "mu_S_gev": contract["matching"]["mu_S_gev"],
            "fermion_tree_LLSS": "provisional numeric components only",
            "hard_one_loop_fermion_matching_complete": False,
            "scalar_matching_performed": False,
        },
        "physical_certification": {
            "status": "BLOCKED_GAMMACC_AND_MATCHING_ORDER",
            "fermion_boundary_certified": False,
            "ready_for_physical_EFT1_or_scalar_matching": False,
            "required_gates": [
                "native signed GammaCC-to-master-Weyl 4-point amplitude including Dirac exchange/crossing and CG",
                "native master Wilson prefactor and external Grassmann ordering",
                "hard one-loop fermion threshold terms classified at required perturbative order",
                "UV frozen couplings consistently interpreted/evolved to mu_F",
                "scalar phase, portal and two-scale matching conventions (Batch 2)"
            ],
            "no_physical_sign_was_assigned": True,
        },
        "execution_context": "local execution on supplied project root only; not a Windows certification by author",
    }


def markdown_report(result: dict) -> str:
    s = result["source_validation"]
    n = result["numeric_input"]
    d = result["source_blobs"]
    drift_lines = [f"- `{path}`: {info['status']}" for path, info in d.items() if not info["match"]]
    drift_section = "\n".join(drift_lines) if drift_lines else "- Inspected GitHub blobs match local files."
    return f"""# T3 Batch 1 — A_-4 fermion-boundary normalization audit

## Current result

**Physical certification: BLOCKED.** The declared Weyl-basis identities establish
`C PL / epsilon = +1` and `C^-1 PL / epsilon = -1`, but the Matchete native
`GammaCC` phase and the signed master Wilson-operator prefactor are unmeasured.
Both candidate signs satisfy the same abstract conjugation identities. Neither
sign has been assigned to the physical boundary.

## Source-checked provisional data

- Native model: A_-4; SU(2) reps (S1,S2,F)=(1,3,2), Dirac F.
- V6 additional HHSS channels: {s['V6_extra_HHSS_channels']}.
- PL source terms: {s['A_m4_native_PL_sources']}; independent LLSS coordinates: {s['A_m4_independent_LLSS_rows']}.
- Formal sample C[1,1,1,1,5,7] = {result['fermion_boundary']['sample_formal_coefficient_gev_inv']['real']:.16e} GeV^-1.
- Already-applied symmetric orbit factors are NOT reapplied. PR conjugates are NOT counted.
- This is a numeric **formal** tensor; no physical overall sign/prefactor is certified.

## Frozen source benchmark (UV values, not RG-evolved to mu_F)

- mu_UV={result['numeric_input']['yukawa_input_scale_gev']:g} GeV.
- mu_F={result['matching']['mu_F_gev']:g} GeV; mu_S={result['matching']['mu_S_gev']:g} GeV.
- Dirac F masses: {n['fermion_mass_eigenvalues_gev']} GeV.
- Unbroken scalar mass parameters: {n['scalar_unbroken_mass_parameters_gev']} GeV.
- min(MF)/max(mS)={n['mF_min_over_mS_max']:.8g}; hence a perturbative fermion-heavy expansion is plausible.
- lambda_T3={n['lambda_T3_source_value']}. See JSON input for full Yukawa matrices and baseline quartics.
- Additional independent couplings fixed to zero in this **restricted benchmark at mu_UV**, not claimed RG-stable.
- Broken-phase scalar masses, physical source portal normalisation and running between mu_UV/mu_F remain uncertified.

## GitHub snapshot drift audit

Read-only main at `{REPO_COMMIT}`. Comparison to local checkout:
{drift_section}

## Next missing amplitude and prescription

Compute one signed, connected native Matchete UV four-point amplitude for
`l l S1 S2` in the Dirac-doublet A_-4 model at zero external momentum, including
both exchange orientations, native `T3MixCG`, Grassmann ordering and the operator
prefactor. Project it independently onto a **declared** two-component Weyl
`epsilon_12=+1` amplitude. Do not choose `+1` from `C=i gamma2 gamma0` alone.
Then check whether already-Matchete-matched one-loop F terms contain hard/soft
subtraction before adding any EFT1 subtraction or counterterm; avoid double counting.

No C5 or m_nu calculation is claimed by this batch.
"""


def write_immutable(path: Path, body: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == body, f"Refusing to overwrite changed file: {path}")
    else:
        with path.open("xb") as f:
            f.write(body)


def run(root: Path, run_id: str, check_only: bool = False) -> tuple[dict, Path | None]:
    require(run_id and len(run_id) <= 40 and
            all(c.isascii() and (c.isalnum() or c in "_-") for c in run_id),
            "Unsafe run ID")
    result = audit(root)
    if check_only:
        return result, None
    path = root.resolve() / OWNED_OUTPUT / f"batch1_A_m4_{run_id}"
    require(path.resolve().is_relative_to((root.resolve() / OWNED_OUTPUT).resolve()),
            "Output escaped owned workstream")
    js = (json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    md = markdown_report(result).encode("utf-8")
    # All source validation is completed before writing either report.
    write_immutable(path / "summary.json", js)
    write_immutable(path / "physical_normalization_report.md", md)
    return result, path
