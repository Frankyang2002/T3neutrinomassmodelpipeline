"""Reconstruct the ordinary-T3 charge audits from the 16 saved native CG exports.

This is a read-only consumer of the existing UV CG exports and frozen benchmarks.
The candidate representation generators are copied from the repository's unit-
test fixtures (tests/test_t3_matchete_generator_analysis.py); they are NOT
asserted to have been extracted live from Matchete. Invariance of EVERY
exported tensor under these generators must be verified before a charge-basis
coefficient can be produced.

Generated charge-resolved weights are group/CG contraction diagnostics only,
NOT complete one-loop contributions or an overall-normalisation certificate.
Neither existing Lagrangian/EFT modules nor existing numeric outputs are edited.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from Numerical.diagnostics.AnalyseMatcheteGenerators import analyse as analyse_generators
from Numerical.diagnostics.AuditT3ChargeSpectrum import audit as audit_spectrum
from Numerical.diagnostics.AuditT3ChargeChannels import audit as audit_channels
from Numerical.diagnostics.VerifyT3CGCovariance import audit as audit_covariance
from Numerical.diagnostics.ResolveT3ChargeWeights import run as resolve_weights

EXPECTED_CG_NAMES = frozenset(("T3Y1CG", "T3Y2CG", "T3MixCG"))
EXPECTED_COUNT = 16


def generator_candidate() -> dict[str, Any]:
    """Generators from checked-in test fixtures, *not* a live Matchete probe.

    Source: tests/test_t3_matchete_generator_analysis.py.
    Original ordering is G1, G2(Cartan), G3; SU(2) positive orientation uses
    G1, G3, G2 in AnalyseMatcheteGenerators and VerifyT3CGCovariance.
    """
    return {
        "status": "candidate_from_repository_test_fixtures_not_live_export",
        "candidate_provenance": "tests/test_t3_matchete_generator_analysis.py",
        "generator_basis_is_live_export": False,
        "representations": [
            {
                "dimension": 2,
                "generators_input_form": [
                    [["0", "1/2"], ["1/2", "0"]],
                    [["1/2", "0"], ["0", "-1/2"]],
                    [["0", "-I/2"], ["I/2", "0"]],
                ],
            },
            {
                "dimension": 3,
                "generators_input_form": [
                    [["0", "0", "0"], ["0", "0", "I"], ["0", "-I", "0"]],
                    [["0", "0", "-I"], ["0", "0", "0"], ["I", "0", "0"]],
                    [["0", "I", "0"], ["-I", "0", "0"], ["0", "0", "0"]],
                ],
            },
        ],
    }


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError(f"Missing required file: {path}")
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _model_key(obj: dict[str, Any], names: tuple[str, ...]) -> tuple[int, ...]:
    return tuple(int(obj[name]) for name in names)


def verify_source_inventory(
    project_root: Path, cg_root: Path, components: dict[str, Any]
) -> list[tuple[int, int, int, int]]:
    """Cross-check every export and existing component audit against frozen inputs.

    Prevents accidental acceptance of 16 duplicated/wrong-alpha outputs, even
    if the upstream audit reported 16 records and no parse failures.
    """
    benchmark_dir = project_root / "Numerical" / "benchmarks" / "fitted_reference"
    benchmarks = sorted(benchmark_dir.glob("T3_*.json"))
    if len(benchmarks) != EXPECTED_COUNT:
        raise ValueError(f"Expected 16 frozen benchmarks, found {len(benchmarks)}")
    source_keys: set[tuple[int, int, int, int]] = set()
    for path in benchmarks:
        representation = _load(path)["representation"]
        expected = _model_key(representation, ("d_s1", "d_s2", "d_f", "alpha"))
        exported = _load(cg_root / path.stem / "uv_cg_registry.json")
        if exported.get("status") != "Success":
            raise ValueError(f"Native CG export not successful: {path.stem}")
        got = _model_key(exported["metadata"], ("dS1", "dS2", "dF", "alpha"))
        if got != expected:
            raise ValueError(f"CG metadata mismatch for {path.stem}: {got} != {expected}")
        names = [str(x["Name"]).split("`")[-1] for x in exported["CGRegistry"]]
        if set(names) != EXPECTED_CG_NAMES or len(names) != 3:
            raise ValueError(f"Unexpected/duplicate native CG tensors for {path.stem}: {names}")
        if got in source_keys:
            raise ValueError(f"Duplicate native model: {got}")
        source_keys.add(got)
    if components.get("errors"):
        raise ValueError("Existing component audit contains errors")
    cases = components.get("cases")
    if not isinstance(cases, list) or len(cases) != EXPECTED_COUNT:
        raise ValueError("Existing component audit must contain exactly 16 cases")
    audit_keys = [_model_key(c["model"], ("dS1", "dS2", "dF", "alpha")) for c in cases]
    if len(set(audit_keys)) != EXPECTED_COUNT or set(audit_keys) != source_keys:
        raise ValueError("Component-audit metadata does not match the frozen 16-model input set")
    for case in cases:
        if set(case["vertices"]) != EXPECTED_CG_NAMES:
            raise ValueError("One or more native CG vertices are absent from component audit")
    return sorted(source_keys)


def _write_new_or_identical(path: Path, obj: dict[str, Any]) -> None:
    """Never silently overwrite any pre-existing source or numerical evidence."""
    contents = json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text(encoding="utf-8-sig") == contents:
            return
        raise FileExistsError(f"Refusing to overwrite a different existing audit: {path}")
    # Exclusive creation prevents accidental overwriting by concurrent runs.
    with path.open("x", encoding="utf-8") as f:
        f.write(contents)


def run(project_root: Path, cg_root: Path, audit_dir: Path) -> dict[str, Any]:
    components = _load(audit_dir / "component_support_audit.json")
    expected = verify_source_inventory(project_root, cg_root, components)
    print(f"Native CG registry metadata independently checked: {len(expected)}/16")

    spectrum = audit_spectrum(components)
    charges = audit_channels(spectrum)
    if len(spectrum["models"]) != EXPECTED_COUNT or len(charges["models"]) != EXPECTED_COUNT:
        raise ValueError("Spectrum/charge inventory incomplete")
    _write_new_or_identical(audit_dir / "charge_spectrum_audit.json", spectrum)
    _write_new_or_identical(audit_dir / "charge_channels.json", charges)
    print("Charge spectra and allowed channels: 16/16 (not physical loop multiplicities)")

    candidate = generator_candidate()
    analysis = analyse_generators(candidate)
    if any(row["status"] != "generator_algebra_verified" for row in analysis["representations"]):
        raise ValueError("Repository test-fixture generators failed SU(2) algebra checks")
    _write_new_or_identical(audit_dir / "generator_candidate_probe.json", candidate)
    _write_new_or_identical(audit_dir / "matchete_generator_analysis.json", analysis)
    covariance = audit_covariance(components, candidate)
    _write_new_or_identical(audit_dir / "cg_covariance_audit.json", covariance)
    if (covariance.get("status") != "cg_generator_covariance_verified" or
        len(covariance.get("models", [])) != EXPECTED_COUNT or
        not all(row.get("all_covariant") is True for row in covariance["models"])):
        raise ValueError(
            "Candidate generators are NOT verified against all 16 native CG tensors. "
            "See cg_covariance_audit.json. No charge weights were produced."
        )
    print("SU(2) generator-candidate covariance against source Matchete CGs: 16/16")
    weights = resolve_weights(components, analysis, covariance, charges)
    if len(weights.get("models", [])) != EXPECTED_COUNT:
        raise ValueError("Charge-resolved CG contraction inventory incomplete")
    _write_new_or_identical(audit_dir / "charge_resolved_loop_weights.json", weights)
    print("Charge-resolved CG coefficients: 16/16 (not normalized physical UV masses)")
    return {"models": len(expected), "covariant": True, "physical_normalization_certified": False}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, default=Path.cwd())
    p.add_argument("--cg-root", type=Path)
    p.add_argument("--audit-dir", type=Path)
    a = p.parse_args()
    project = a.project_root.resolve()
    cg_root = a.cg_root or project / "output/general_observable/lo_direct_one_loop/source_cg_all16"
    audit_dir = a.audit_dir or project / "output/general_observable/lo_direct_one_loop/source_audits"
    try:
        run(project, cg_root.resolve(), audit_dir.resolve())
    except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        p.exit(2, f"BLOCKED: {exc}\n")
    print("Status: source CG charge audit complete; physical UV vertex conventions pending")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
