"""Rebuild the B_-1/C_-1 F-first Matchete source in the SAME running mass basis.

Uses the most recent full-general numerical ``running.json`` at mu_F, and
reuses the repository's validated B/C scalar tensor transport, eigenbasis
Matchete validation and fermion-first *tree* matching. No historical mixed
scalar endpoint is substituted. Matchete quartic/RGBeta physical
normalisations and actual EFT1 Wilson RGEs remain uncertified/unimplemented.

PowerShell:
  python -X utf8 -m scripts.RunT3GeneralMixedFSynchronised --all --prepare-only
  python -X utf8 -m scripts.RunT3GeneralMixedFSynchronised --all

Only new output directories are written; no tracked inputs are modified.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("B_-1", "C_-1")
FORMAT = "t3-mixed-F-boundary-synchronised-to-running-v1"
OUTPUT = Path("output/general_observable/mixed_f_synchronised_v1")
RUN_ROOTS = (
    Path("output/general_observable/mixed_adapter_v3"),
    Path("output/general_observable/all16_v1"),
)
MANIFEST = Path("docs/general_scalar_trial_operator_manifest.json")
REGISTRY = Path("RGE/general/generated/general_rgbeta_quartic_registry.json")
MIXING_NAME = "muT3GeneralS1S2"


def tag(model: str) -> str:
    if model not in MODELS:
        raise ValueError(f"Only {MODELS} supported, got {model}")
    return model.replace("-", "m").replace("+", "p")


def read(path: Path) -> dict:
    obj = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object at {path}")
    return obj


def save(path: Path, obj: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def number(raw: Any, name: str) -> complex:
    if isinstance(raw, dict) and set(raw) == {"real", "imag"}:
        value = complex(float(raw["real"]), float(raw["imag"]))
    elif isinstance(raw, (int, float, complex)) and not isinstance(raw, bool):
        value = complex(raw)
    else:
        raise ValueError(f"Invalid numeric value for {name}")
    if not np.isfinite(value.real) or not np.isfinite(value.imag):
        raise ValueError(f"Nonfinite {name}")
    return value


def positive_real(raw: Any, name: str) -> float:
    z = number(raw, name)
    if abs(z.imag) > 1e-12 * max(1.0, abs(z)) or z.real <= 0:
        raise ValueError(f"{name}: positive real mass-squared required")
    return z.real


def project_path(path: Path) -> Path:
    """Use absolute repository paths for downstream Matchete certificates.

    The existing mass-eigenstate/fermion-tree modules call
    Path.relative_to(ROOT) without resolving their input paths first.
    Thus passing a relative output directory can fail *after* a valid
    Wolfram calculation, when producing its return record.
    """
    path = Path(path)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def running_file(model: str, roots: tuple[Path, ...]) -> Path:
    candidates = [ROOT / root / tag(model) / "running.json" for root in roots]
    found = [p for p in candidates if p.is_file()]
    if not found:
        raise FileNotFoundError(f"{model}: running.json absent in {candidates}")
    return found[0]


def running_state(record: Mapping, model: str) -> tuple[dict, float]:
    if record.get("model") != model or record.get("status") != "NumericalRenormalisableRunningOnly":
        raise ValueError(f"{model}: unexpected general running status/model")
    states = record.get("uv_at_fermion_threshold")
    if not isinstance(states, dict):
        raise ValueError(f"{model}: missing evolved UV state at mu_F")
    needed = ("mS1Sq", "mS2Sq", MIXING_NAME, "y1", "y2", "MF", "lambdaH")
    absent = set(needed) - set(states)
    if absent:
        raise ValueError(f"{model}: missing UV threshold parameters: {sorted(absent)}")
    scales = record.get("scales_gev", {})
    # Support the two existing general running output contracts, but never
    # assume a default matching scale when the value is absent.
    mu = scales.get("F", scales.get("mu_F", scales.get("mu_f")))
    if mu is None:
        raise ValueError(f"{model}: no explicit fermion matching scale")
    mu = float(mu)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError(f"{model}: invalid fermion matching scale")
    return states, mu


def make_mass_basis(model: str, running: Mapping, source_file: Path,
                    manifest_file: Path, registry_file: Path,
                    *, trials: int = 4) -> dict:
    """Recalculate all 17 tensor rotations and exact X/Y mass basis at mu_F.

    The physical RGBeta<->Matchete coupling normalization is NOT established
    by any of these algebraic mass/field redefinitions.
    """
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import (
        MixedDoubletMassBasis, cjson, transport_operator, verify_su2_covariance,
        verify_tensor_transport, quadratic_mass_tensor_transport,
    )
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import quadratic_mass_hessian
    from scripts import RunGeneralT3MixedScalarMassBasis as old

    if model not in MODELS or not 1 <= trials <= 100:
        raise ValueError("Invalid B/C model or tensor-transport trial count")
    state, mu_f = running_state(running, model)
    m1 = positive_real(state["mS1Sq"], "mS1Sq")
    m2 = positive_real(state["mS2Sq"], "mS2Sq")
    mu2 = number(state[MIXING_NAME], MIXING_NAME)
    basis = MixedDoubletMassBasis.build(model, m1, m2, mu2)
    if basis.tachyonic or basis.degenerate or np.any(basis.eigenvalues_sq <= 0):
        raise ValueError(f"{model}: tachyonic/degenerate F-scale mixed-scalar basis")
    # The two remaining fields must be light at the F threshold. Report
    # scalar masses above mu_F as an invalid grouped F-first hierarchy.
    masses = np.sqrt(np.asarray(basis.eigenvalues_sq, dtype=float))
    if not np.max(masses) < mu_f:
        raise ValueError(f"{model}: scalar eigenmass >= grouped F threshold")

    manifest, registry = read(manifest_file), read(registry_file)
    operators = old.quartic_operators(manifest, registry, model)
    if len(operators) != 17:
        raise ValueError(f"{model}: expected all 17 quartic tensors")
    prototype = quadratic_mass_hessian(manifest, model, m1, m2, mu2)
    expected = np.repeat(basis.eigenvalues_sq, 4)
    actual = np.asarray(prototype["mass_squared_eigenvalues_gev2"], dtype=float)
    if actual.shape != (8,) or not np.allclose(actual, expected, rtol=2e-11, atol=1e-5):
        raise ValueError(f"{model}: independent eight-real-component Hessian mismatch")

    # Reuse the exact source field-count inventory and transport functions.
    # The running JSON flattens SM couplings; old._state_value expects lambdaH
    # under the sm substate. Construct only that view, never invent couplings.
    old_state = dict(state, sm={"lambdaH": state["lambdaH"]})
    result = {
        "format": "t3-mixed-doublet-mass-basis-transport-v1",
        "status": old.STATUS,
        "model": model, "endpoint": "final", "scale_gev": mu_f,
        "basis": basis.info(),
        "quadratic_hessian_agreement": True,
        "exact_quadratic_tensor_transport": quadratic_mass_tensor_transport(basis),
        "max_su2_covariance_absolute_error": verify_su2_covariance(basis),
        "source_quartic_count": 17, "source_new_quartics": 1,
        "source_reused_quartics": 6, "source_baseline_quartics": 10,
        "quartics": {}, "yukawa_scalar_leg_transforms": {},
        "physical_mass_matching_performed": False,
        "MatcheteRGBetaPhysicalCouplingMapVerified": False,
        "MatcheteHeavyMassBasisCheckOnRedefinedLagrangian": False,
        "source_eft1_dirac_yukawas_active": False,
        "source_sha256": {relative(p): sha(p)
                          for p in (source_file, manifest_file, registry_file)},
        "running_boundary_provenance": {
            "running_file": relative(source_file),
            "running_sha256": sha(source_file),
            "boundary": "uv_at_fermion_threshold",
            "boundary_scale_gev": mu_f,
            "historical_mixed_EFT1_endpoint_reused": False,
        },
    }
    max_error = 0.0
    for name, fields in operators.items():
        z = old._state_value(old_state, name, {})
        error = float(verify_tensor_transport(fields, basis.U, trials=trials))
        max_error = max(max_error, error)
        result["quartics"][name] = {
            "source_rgbeta_fields": list(fields),
            "source_rgbeta_coupling": cjson(z),
            "independent_tensor_preserved": True,
            "source_invariant_component_tensor_supplied": False,
            "coefficient_normalization_certified": False,
            "rotation_terms": transport_operator(fields, basis.U),
            "transport_identity_max_relative_error": error,
        }
    for name, field in old.YUKAWA_SCALAR_LEGS.items():
        max_error = max(max_error, float(verify_tensor_transport((field,), basis.U, trials=trials)))
        result["yukawa_scalar_leg_transforms"][name] = {
            "rgbeta_source_scalar_field": field,
            "non_scalar_indices_and_original_cg_unchanged": True,
            "yukawa_numerical_value_not_defined_in_eft1": True,
            "rotation_terms": transport_operator((field,), basis.U),
        }
    if max_error > 3e-12:
        raise ValueError(f"{model}: numerical scalar tensor transport residual {max_error}")
    result["max_source_tensor_transport_relative_error"] = max_error
    return result


def process(model: str, roots: tuple[Path, ...], output_root: Path,
            *, prepare_only: bool = False, recover_existing: bool = False,
            wolframscript: str = "wolframscript",
            timeout_sec: int = 3600, trials: int = 4) -> dict:
    """Create new matched source, never substituting a historical U(mu_S)."""
    from scripts import RunGeneralT3MassEigenMatchete as matchete
    from scripts import RunGeneralT3FermionFirstEFT1Tree as tree
    output_root = project_path(output_root)
    work = output_root / tag(model)
    mass_dir = work / "mass_basis"
    native_dir = work / "matchete"
    tree_dir = work / "fermion_tree"
    source = running_file(model, roots)
    mass_basis = make_mass_basis(model, read(source), source,
                                 ROOT / MANIFEST, ROOT / REGISTRY, trials=trials)
    mass_path = mass_dir / f"{tag(model)}_final_mass_basis.json"
    previous_mass_basis = read(mass_path) if recover_existing and mass_path.is_file() else None
    if not recover_existing:
        save(mass_path, mass_basis)
    record = {
        "format": FORMAT,
        "model": model,
        "status": "CurrentRunningMassBasisPrepared" if prepare_only else "BlockedPendingLiveMatchete",
        "running_file": relative(source),
        "mass_basis_file": relative(mass_path),
        "mu_F_gev": mass_basis["scale_gev"],
        "scalar_eigenmasses_gev": mass_basis["basis"]["masses_gev"],
        "scalar_basis_updated_from_current_mu_F": True,
        "historical_mu_S_basis_reused": False,
        "numerical_flavor_wilson_matrix_computed": False,
        "evolved_MF_and_Yukawas_inserted_into_Matchete_tree": False,
        "full_eft1_wilson_rge_computed": False,
        "scalar_threshold_C5_computed": False,
        "physical_matchete_rgbeta_coupling_map_certified": False,
    }
    if recover_existing and not prepare_only:
        # The previous run can have completed Matchete and tree matching before
        # failing on Path.relative_to(ROOT).  Reuse it ONLY when both existing
        # native validators pass and the previous mass basis is derived from
        # this identical running.json and the identical scalar state.
        native_report = native_dir / f"{tag(model)}_final_matchete.json"
        tree_stem = tree._stem(model, "final")
        tree_report = tree_dir / f"{tree_stem}.json"
        if previous_mass_basis is None or not native_report.is_file() or not tree_report.is_file():
            raise ValueError(f"{model}: no complete previous certificates to recover; rerun without --recover-existing")
        old_source = previous_mass_basis.get("running_boundary_provenance", {})
        if (old_source.get("running_sha256") != sha(source) or
                old_source.get("boundary") != "uv_at_fermion_threshold" or
                old_source.get("boundary_scale_gev") != mass_basis["scale_gev"] or
                previous_mass_basis.get("source_sha256") != mass_basis["source_sha256"] or
                set(previous_mass_basis.get("quartics", {})) != set(mass_basis["quartics"]) or
                any(previous_mass_basis["quartics"][key].get("source_rgbeta_coupling") !=
                    mass_basis["quartics"][key]["source_rgbeta_coupling"]
                    for key in mass_basis["quartics"])):
            raise ValueError(f"{model}: previous mass-basis input differs from current running source; rerun without --recover-existing")
        old_basis = previous_mass_basis.get("basis", {})
        new_basis = mass_basis["basis"]
        for field in ("unitary_flavour_rotation_U", "mass_squared_flavour_matrix_gev2"):
            old = old_basis.get(field, {})
            new = new_basis[field]
            try:
                old_matrix = np.asarray(old["real"], dtype=float) + 1j*np.asarray(old["imag"], dtype=float)
                new_matrix = np.asarray(new["real"], dtype=float) + 1j*np.asarray(new["imag"], dtype=float)
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"{model}: missing saved {field} for recovery") from exc
            if old_matrix.shape != (2, 2) or not np.allclose(old_matrix, new_matrix,
                                                            rtol=1e-12, atol=1e-9):
                raise ValueError(f"{model}: previous {field} is inconsistent with current running")
        matchete.validate_wolfram_output(read(native_report), model,
                                         new_basis["unitary_flavour_rotation_U"])
        matched_eft, matched_llss = tree.validate_tree_report(
            read(tree_report), model, tree_dir, tree_stem)
        report_data = read(tree_report)
        record.update(status="CurrentRunningMassBasisFThresholdTreeLLSSExported",
                      native_parameters_file=relative(native_dir / f"{tag(model)}_final_parameters.json"),
                      mass_basis_matchete_certificate=relative(native_report),
                      matched_llss_file=relative(matched_llss),
                      matched_eft_file=relative(matched_eft),
                      llss_term_count=report_data["llss_term_count"],
                      source_operator_matching_at_running_mu_F=True,
                      reused_existing_validated_wolfram_outputs=True)
        save(work / "status.json", record)
        return record
    if prepare_only:
        # This also performs the strict preflight of the existing Matchete
        # source contract without requiring a Wolfram license.
        prepared = matchete.prepare_model(model, which="final", folder=mass_dir,
                                          outdir=native_dir)
        record.update(mass_basis_preflight="Passed",
                      native_parameters_file=relative(prepared["parameters"]))
    else:
        cert = matchete.run_model(model, which="final", folder=mass_dir,
                                  outdir=native_dir, wolframscript=wolframscript,
                                  timeout_sec=timeout_sec)
        # The tree routine independently checks the PRECISE fresh U matrix
        # against the prior Matchete certificate. No cached reference source.
        result = tree.run_model(model, which="final", folder=mass_dir,
                                certificate_dir=native_dir, outdir=tree_dir,
                                wolframscript=wolframscript,
                                timeout_sec=timeout_sec)
        record.update(status="CurrentRunningMassBasisFThresholdTreeLLSSExported",
                      native_parameters_file=relative(cert["parameters"]),
                      mass_basis_matchete_certificate=relative(native_dir / f"{tag(model)}_final_matchete.json"),
                      matched_llss_file=result["llss_file"],
                      matched_eft_file=result["matched_eft_file"],
                      llss_term_count=result["llss_term_count"],
                      source_operator_matching_at_running_mu_F=True)
    save(work / "status.json", record)
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--all", action="store_true")
    group.add_argument("--model", choices=MODELS)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--recover-existing", action="store_true",
                        help="Revalidate pre-existing native Matchete/tree outputs without Wolfram rerun")
    parser.add_argument("--running-roots", default=",".join(map(str, RUN_ROOTS)))
    parser.add_argument("--output-root", type=Path, default=OUTPUT)
    parser.add_argument("--wolframscript", default="wolframscript")
    parser.add_argument("--timeout-sec", type=int, default=3600)
    parser.add_argument("--trials", type=int, default=4)
    args = parser.parse_args(argv)
    if args.prepare_only and args.recover_existing:
        parser.error("--prepare-only and --recover-existing cannot be combined")
    output_root = project_path(args.output_root)
    selected = MODELS if args.all else (args.model,)
    roots = tuple(Path(x.strip()) for x in args.running_roots.split(",") if x.strip())
    results = []
    for i, model in enumerate(selected, 1):
        print(f"[{i}/{len(selected)}] {model}: F-scale mixed mass-basis refresh", flush=True)
        try:
            report = process(model, roots, output_root,
                             prepare_only=args.prepare_only, recover_existing=args.recover_existing,
                             wolframscript=args.wolframscript,
                             timeout_sec=args.timeout_sec, trials=args.trials)
            print(f"  {report['status']}; m_X,Y = {report['scalar_eigenmasses_gev']}", flush=True)
        except Exception as exc:
            report = {"model": model, "status": "Blocked", "error": f"{type(exc).__name__}: {exc}",
                      "physical_matchete_rgbeta_coupling_map_certified": False,
                      "full_eft1_wilson_rge_computed": False}
            print(f"  BLOCKED: {report['error']}", flush=True)
            save(output_root / tag(model) / "status.json", report)
        results.append(report)
    count = sum(r["status"] == "CurrentRunningMassBasisFThresholdTreeLLSSExported" for r in results)
    summary = {"format": FORMAT, "results": results, "model_count": len(results),
               "prepared_count": sum(r["status"] == "CurrentRunningMassBasisPrepared" for r in results),
               "live_tree_success_count": count,
               "full_eft1_wilson_rge_computed": False,
               "physical_matching_certified": False}
    save(output_root / "summary.json", summary)
    print(f"Synchronized B/C F-first tree: {count}/{len(results)}; EFT1 Wilson RGEs NOT performed", flush=True)
    return 0 if all(r["status"] != "Blocked" for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
