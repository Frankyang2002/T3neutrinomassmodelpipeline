"""Read-only SU(2) boundary/ADM diagnostic for retained T3 fit16 trials.

Run from repository root:
    python -X utf8 scripts/DiagnoseT3Fit16Ward.py --capture output/parallel_t3/all16_ward_capture01

No Wilson coefficients are altered; no EFT evolution or fitting is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.adm import beta, total_beta
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical
from Numerical.parallel_t3.threshold.matching import (
    extract_lepton_sectors, master_weinberg_template,
)
from Numerical.parallel_t3.threshold.rough_projection import diagnose_llhh


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def full_ward(c: np.ndarray, tf: np.ndarray, ts: np.ndarray) -> float:
    """Same index action as wilson.adm.gauge_ward_residual, Frobenius-normalised."""
    variation = (
        np.einsum("ki,kjab->ijab", tf, c, optimize=True)
        + np.einsum("kj,ikab->ijab", tf, c, optimize=True)
        + np.einsum("ca,ijcb->ijab", ts, c, optimize=True)
        + np.einsum("cb,ijac->ijab", ts, c, optimize=True)
    )
    return float(np.linalg.norm(variation) / max(np.linalg.norm(c), 1e-100))


def full_wards(c: np.ndarray, basis) -> dict:
    generators = basis.generators()
    return {
        f"{group}_{i+1}": full_ward(c, tf, ts)
        for group in ("SU2", "U1")
        for i, (tf, ts) in enumerate(generators[group])
    }


def llhh(c: np.ndarray, manifest: dict, basis) -> dict:
    hh, _ = extract_lepton_sectors(
        c,
        manifest["lepton_indices"],
        manifest["higgs_source_indices"],
        manifest["heavy_source_indices"],
        number_of_real_scalars=basis.ns,
    )
    d = diagnose_llhh(hh, master_weinberg_template(), top_components=0)
    return {
        "norm": d["source_LLHH_norm_gev_inverse"],
        "purity": d["relative_Weinberg_purity_residual"]
        if d["source_LLHH_norm_gev_inverse"] > 1e-100 else None,
        "wards": {
            k: v["relative_to_LLHH"]
            for k, v in d["ward_SU2xU1"].items()
        } if d["source_LLHH_norm_gev_inverse"] > 1e-100 else None,
    }


def diagnose(manifest_path: Path, project: Path) -> dict:
    manifest_path = manifest_path.resolve()
    token = manifest_path.parent.parent.name
    model = "A_-4"
    if not token.startswith("fit16_") or manifest_path.parent.name != "A_m4":
        raise ValueError("Expected retained fit16_*/A_m4/wilson.json")
    meta = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    if meta.get("model") != model:
        raise ValueError("Manifest model mismatch")
    if meta.get("adm_variant") != "fonseca_eq22_six_term":
        raise ValueError("Unexpected ADM variant; refusing to assume trial settings")
    npz_path = manifest_path.parent / meta["tensor_npz"]
    if meta.get("tensor_npz_sha256") and sha256(npz_path) != meta["tensor_npz_sha256"]:
        raise ValueError("Wilson NPZ checksum mismatch")
    with np.load(npz_path, allow_pickle=False) as f:
        endpoint = {key: f[key] for key in f.files}

    canonical = project / "output/parallel_t3/rough_inputs" / token / "a/A_m4"
    case = load_canonical(project, canonical, model,
                          c5_boundary="require", assume_formal_identity=True)
    basis = case.basis
    boundary = np.asarray(endpoint["C_psi2phi2_tree_mf"], complex)
    evolved = np.asarray(endpoint["C_psi2phi2_ms"], complex)
    if not np.allclose(boundary, case.initial, rtol=1e-12, atol=1e-100):
        raise ValueError("Retained canonical boundary differs from Wilson endpoint boundary")
    if not math.isclose(float(meta["mu_f_gev"]), case.mf, rel_tol=1e-10):
        raise ValueError("MF mismatch")
    y, lam, gauge = case.background.at(math.log(case.mf))
    channels, backend = beta(boundary, y, lam, gauge,
                             variant=meta["adm_variant"])
    total = total_beta(channels)
    if "beta_C_at_MF_16pi2" in endpoint:
        stored = endpoint["beta_C_at_MF_16pi2"]
        mismatch = float(np.linalg.norm(total - stored) /
                         max(np.linalg.norm(stored), 1e-100))
        if mismatch > 1e-8:
            raise ValueError(f"Recomputed MF beta disagrees with stored beta: {mismatch:g}")
    else:
        mismatch = None

    def metrics(c: np.ndarray) -> dict:
        return {
            "full_norm": float(np.linalg.norm(c)),
            "full_ward_relative": full_wards(c, basis),
            "LLHH": llhh(c, meta, basis),
        }

    return {
        "status": "READ_ONLY_DIAGNOSTIC_NOT_A_REPAIR",
        "model": model,
        "token": token,
        "wilson_manifest": str(manifest_path),
        "wilson_npz_sha256": sha256(npz_path),
        "canonical_input": str(canonical),
        "MF_GeV": case.mf,
        "MS_GeV": case.ms,
        "adm_variant": meta["adm_variant"],
        "adm_backend": backend,
        "stored_beta_relative_mismatch": mismatch,
        "boundary_MF": metrics(boundary),
        "evolved_MS": metrics(evolved),
        "beta_total_MF_16pi2": metrics(total),
        "beta_channels_MF_16pi2": {
            name: metrics(value) for name, value in channels.items()
        },
        "interpretation_note": (
            "Full-tensor Ward tests use Frobenius norms and source basis generators. "
            "A zero LLHH boundary has no meaningful relative LLHH Ward/purity. "
            "Individual beta channels need not be gauge invariant separately; "
            "their sum and the full boundary are the primary covariance tests. "
            "A nonzero beta LLHH Ward residual alone cannot distinguish a bad "
            "boundary from a bad ADM."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--capture", type=Path,
                   default=Path("output/parallel_t3/all16_ward_capture01"))
    p.add_argument("--project-root", type=Path, default=Path("."))
    args = p.parse_args()
    project = args.project_root.resolve()
    capture = (project / args.capture).resolve()
    if not capture.is_dir():
        raise FileNotFoundError(f"Capture directory missing: {capture}")
    tokens = sorted({
        path.name.split("_ward_diagnostic.json")[0]
        for path in capture.glob("fit16_*_ward_diagnostic.json")
        if path.name.endswith("_ward_diagnostic.json")
    })
    if not tokens:
        raise FileNotFoundError("No retained trial diagnostic JSONs in capture directory")
    for token in tokens:
        if not token.startswith("fit16_"):
            raise ValueError(f"Unexpected token: {token}")
        manifest = (project / "output/parallel_t3/rough_wilson" /
                    token / "A_m4/wilson.json")
        if not manifest.is_file():
            raise FileNotFoundError(f"Retained Wilson manifest missing: {manifest}")
        result = diagnose(manifest, project)
        output = capture / f"{token}_full_boundary_adm_ward.json"
        output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")
        print(f"\n{token}")
        print(f"  MF full Ward: {result['boundary_MF']['full_ward_relative']}")
        print(f"  MS full Ward: {result['evolved_MS']['full_ward_relative']}")
        print(f"  MF beta full Ward: {result['beta_total_MF_16pi2']['full_ward_relative']}")
        print(f"  MF beta LLHH: {result['beta_total_MF_16pi2']['LLHH']}")
        print(f"  MF beta channel LLHH norms:")
        for name, d in result["beta_channels_MF_16pi2"].items():
            print(f"    {name}: {d['LLHH']['norm']:.6e}")
        print(f"  Report: {output}")


if __name__ == "__main__":
    main()
