"""General T3 F-threshold Wilson *flavour-factor* handoff, all 16 models.

Uses the actual Matchete tree-EFT operator strings. For each nonmixed model,
replaces its two source Yukawa factors and one heavy mass denominator by ONE,
retaining all scalar fields, CGs, signs, phases and symmetry prefactors in the
original operator template. The corresponding ordered 3x3 flavour matrix is
computed from the evolved UV state at mu_F in the existing heavy mass basis:

    W[p,q] = sum_r Y_A[p,r]^(bar?) Y_B[q,r]^(bar?) / M_F[r].

Thus a source term is represented as operator_template * W[p,q]. This is
an *operator-level factorisation*, not a canonical real C_ijab tensor.
The tree source needs its own exact SU(2)/scalar mass-basis map before the
existing scalar-only Wilson master RGE can be applied to every model.

B_-1/C_-1 preserve their independently matched X/Y operator source but do
NOT import numerical masses/couplings from an unrelated mixed mass basis.
Neither low-energy C5 nor full sequential EFT is computed here.

From repo root:
 python -X utf8 -m scripts.RunT3GeneralWilsonBoundary --all
 python -X utf8 -m scripts.RunT3GeneralWilsonBoundary --model A_-2
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODELS = tuple(f"{f}_{n:+d}" for f, nums in {
    "A": (-4, -2, 0), "B": (-3, -1, 1), "C": (-3, -1, 1),
    "D": (-2, 0, 2), "E": (-4, -2, 0, 2),
}.items() for n in nums)
MIXED = frozenset(("B_-1", "C_-1"))
DEFAULT_F_ROOTS = (
    Path("output/general_observable/fermion_first_self_fix_v2"),
    Path("output/general_observable/fermion_first_all16_v1"),
)
DEFAULT_R_ROOTS = (
    Path("output/general_observable/mixed_adapter_v3"),
    Path("output/general_observable/vectorlike_fix_v2"),
    Path("output/general_observable/all16_v1"),
)
DEFAULT_OUT = Path("output/general_observable/general_wilson_boundary_v1")
FORMAT = "t3-source-exact-F-tree-flavour-operator-handoff-v1"
# Require explicit 3x3 flavour factors and shared internal heavy index.
_YUK = re.compile(
    r"(?P<bar>Bar\[)?Coupling\[(?P<name>y1Conj|y2Conj|y1|y2),\s*"
    r"\{Index\[(?P<external>[^\],]+),\s*Flavor\],\s*"
    r"Index\[(?P<internal>[^\],]+),\s*NFlavor\]\},\s*0\]"
    r"(?P<close>\])?"
)
_MF = re.compile(
    r"(?P<bar>Bar\[)?Coupling\[MF,\s*(?:\{\}|\{Index\[(?P<internal>[^\],]+),\s*NFlavor\]\}),\s*0\]"
    r"(?P<close>\])?"
)


def tag(model: str) -> str:
    return model.replace("+", "p").replace("-", "m")


def read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def save(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def decode(value: Any) -> Any:
    if isinstance(value, dict) and set(value) == {"real", "imag"}:
        return complex(float(value["real"]), float(value["imag"]))
    if isinstance(value, list):
        return [decode(x) for x in value]
    return value


def encode_matrix(value: np.ndarray) -> dict:
    m = np.asarray(value, dtype=complex)
    if m.shape != (3, 3) or not np.isfinite(m).all():
        raise ValueError("Expected finite 3x3 complex matrix")
    return {"real": m.real.tolist(), "imag": m.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(m))}


@dataclass(frozen=True)
class FlavorFactor:
    name: str
    conjugated: bool
    external_index: str
    heavy_index: str
    span: tuple[int, int]


def _factors(text: str) -> list[FlavorFactor]:
    out = []
    for match in _YUK.finditer(text):
        barred = match.group("bar") is not None
        if barred != (match.group("close") is not None):
            raise ValueError("Unbalanced Bar wrapper on Yukawa factor")
        out.append(FlavorFactor(match.group("name"), barred,
                                match.group("external"), match.group("internal"),
                                match.span()))
    return out


def _normalise_term(text: str) -> dict:
    """Extract only algebraically exact factorizations; fail closed otherwise."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Empty matched Wilson term")
    if "GammaCC" not in text or "Field[NewFermion" in text:
        raise ValueError("Not a first-threshold light-field Majorana Wilson term")
    factors = _factors(text)
    if len(factors) != 2:
        raise ValueError(f"Expected exactly two source Yukawa factors, found {len(factors)}")
    if len(set(f.heavy_index for f in factors)) != 1:
        raise ValueError("No common heavy-flavour Yukawa dummy")
    if len(set(f.external_index for f in factors)) != 2:
        raise ValueError("Expected two distinct external lepton flavour dummies")
    mf = list(_MF.finditer(text))
    if len(mf) != 1:
        raise ValueError(f"Expected exactly one MF symbol, got {len(mf)}")
    mass = mf[0]
    if (mass.group("bar") is not None) != (mass.group("close") is not None):
        raise ValueError("Unbalanced Bar wrapper on MF")
    if mass.group("internal") is not None and mass.group("internal") != factors[0].heavy_index:
        raise ValueError("Heavy mass and Yukawas use distinct NFlavor indices")
    # MF must occur once in the denominator. Matchete's E-family source
    # contains an exact /(Sqrt[3]*MF), whereas most families use /MF.
    # Retain the Sqrt[3] factor in the symbolic operator template.
    prefix = text[:mass.start()].rstrip()
    if prefix.endswith("/"):
        pass
    elif prefix.endswith("*"):
        scale_product = re.search(r"/\(\s*(?:Sqrt\[\d+\]|\d+(?:/\d+)?)\s*\*\s*$", prefix)
        if not scale_product or not text[mass.end():].lstrip().startswith(")"):
            raise ValueError("MF embedded in unverified mass-denominator product")
    else:
        raise ValueError("MF not in a recognizable inverse-mass tree factor")
    if len(re.findall(r"Coupling\[(?:y1Conj|y2Conj|y1|y2|MF)\b", text)) != 3:
        raise ValueError("Unrecognized additional Yukawa/heavy mass insertion")
    spans = [(f.span[0], f.span[1]) for f in factors] + [mass.span()]
    spans.sort()
    if any(a[1] > b[0] for a, b in zip(spans, spans[1:])):
        raise ValueError("Overlapping flavor-factor spans")
    template = text
    for start, stop in reversed(spans):
        template = template[:start] + "1" + template[stop:]
    if "Coupling[MF" in template or re.search(r"Coupling\[y(?:1|2)", template):
        raise ValueError("Unreplaced threshold flavor symbols")
    return {
        "source_term_inputform": text,
        "operator_template_unit_flavor_inputform": template,
        "factors": [dict(name=f.name, conjugated=f.conjugated,
                         external_index=f.external_index,
                         heavy_index=f.heavy_index) for f in factors],
        "heavy_mass_power": -1,
        "symmetric_flavour_combination_assumed": False,
    }


def _validate_seed(seed: dict, model: str) -> list[dict]:
    if seed.get("IntegratedFields") != ["F"] or seed.get("ActiveHeavyFields") != ["S1", "S2"]:
        raise ValueError(f"{model}: not an F-only Matchete seed")
    terms = seed.get("TreeWilsonTerms")
    if not isinstance(terms, list) or seed.get("TreeWilsonTermCount") != len(terms) or not terms:
        raise ValueError(f"{model}: invalid tree Wilson inventory")
    selected = [x for x in terms if x.get("Chirality") == "PL"]
    if not selected:
        raise ValueError(f"{model}: no physical PL threshold Wilson sector")
    result = []
    for term in selected:
        factorisation = _normalise_term(term["TermInputForm"])
        factorisation["source_index"] = int(term["Index"])
        factorisation["cg_names"] = term.get("CGNames", [])
        result.append(factorisation)
    return result


def flavour_matrix(row: Mapping, masses: np.ndarray,
                   yukawas: Mapping[str, np.ndarray]) -> np.ndarray:
    """Numerically lift the source's ORDERED Yukawa product; no ad hoc 1/2."""
    fs = row["factors"]
    m = np.asarray(masses, dtype=float)
    if m.shape != (3,) or not np.isfinite(m).all() or np.any(m <= 0):
        raise ValueError("Invalid positive heavy-mass eigenvalues")
    a = np.asarray(yukawas[fs[0]["name"]], dtype=complex)
    b = np.asarray(yukawas[fs[1]["name"]], dtype=complex)
    if a.shape != (3, 3) or b.shape != (3, 3) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Yukawa must be finite 3x3 complex")
    if fs[0]["conjugated"]:
        a = a.conj()
    if fs[1]["conjugated"]:
        b = b.conj()
    return (a / m[None, :]) @ b.T


def _root_paths(text: str | None, defaults: tuple[Path, ...]) -> tuple[Path, ...]:
    return tuple(Path(x.strip()) for x in text.split(",") if x.strip()) if text else defaults


def _find_model_file(model: str, roots: tuple[Path, ...], tail: str) -> Path:
    found = [ROOT / p / tag(model) / tail for p in roots]
    existing = [p for p in found if p.is_file()]
    if not existing:
        raise FileNotFoundError(f"No {tail} for {model}; searched: {found}")
    return existing[0]


def _running_basis(model: str, running_path: Path, required: set[str]):
    from Numerical.matching.GeneralT3FlavourC5 import prepare_heavy_flavour
    data = read(running_path)
    if data.get("model") != model or data.get("status") != "NumericalRenormalisableRunningOnly":
        raise ValueError(f"{model}: invalid numerical UV/EFT1 output")
    state = data["uv_at_fermion_threshold"]
    required = set(required) | {"y1", "y2"}
    if not required <= set(state):
        raise ValueError(f"{model}: missing evolved Yukawas {sorted(required - set(state))}")
    uv = SimpleNamespace(
        MF=np.asarray(decode(state["MF"]), dtype=complex),
        y1=np.asarray(decode(state["y1"]), dtype=complex),
        y2=np.asarray(decode(state["y2"]), dtype=complex),
        general_couplings={k: np.asarray(decode(state[k]), dtype=complex)
                           for k in required - {"y1", "y2"}},
        representation=SimpleNamespace(majorana_fermion=model.endswith("_-1")),
    )
    masses, yukawas, diagnostics = prepare_heavy_flavour(uv, required)
    return data, masses, yukawas, diagnostics


def _first_threshold(model: str, roots: tuple[Path, ...]) -> tuple[dict, Path, Path]:
    folder_status = _find_model_file(model, roots, "status.json")
    status = read(folder_status)
    if status.get("model") != model:
        raise ValueError(f"Wrong model in {folder_status}")
    if model in MIXED:
        if status.get("status") != "FermionThresholdMixedMassTreeLLSSExported":
            raise ValueError("Mixed model lacks Matchete X/Y tree certificate")
        raw = Path(status["llss_operator_file"])
        if not raw.is_absolute():
            raw = ROOT / raw
        if not raw.is_file():
            raise FileNotFoundError(f"Missing mixed-mass Matchete LLSS source {raw}")
        return status, raw, folder_status
    if status.get("status") != "FermionThresholdTreeWilsonExported":
        raise ValueError(f"{model}: first-threshold tree matching did not succeed")
    folder = folder_status.parent
    sources = list((folder / "data").glob("*wilson_seed.json"))
    if len(sources) != 1:
        raise ValueError(f"{model}: expected exactly one Wilson seed JSON")
    return status, sources[0], folder_status


def process(model: str, f_roots: tuple[Path, ...], r_roots: tuple[Path, ...],
            output: Path) -> dict:
    f, source, certificate = _first_threshold(model, f_roots)
    target = output / tag(model) / "wilson_boundary.json"
    source_info = dict(first_threshold_certificate=str(certificate),
                       first_threshold_sha256=digest(certificate),
                       matchete_operator_source=str(source),
                       matchete_operator_sha256=digest(source))
    if model in MIXED:
        report = {
            "format": FORMAT, "model": model,
            "status": "VerifiedMixedMassFThresholdSource_NumericalFlavourBoundaryNotConnected",
            "source": source_info,
            "operator_term_count": f["tree_wilson_term_count"],
            "mixed_scalar_mass_basis": f.get("mass_basis"),
            "reason": "X/Y source is from the separate prepared mixed-mass basis; current benchmark U(mu_S) not mapped to that Matchete source",
            "numerical_wilson_boundary_computed": False,
            "canonical_real_wilson_tensor_computed": False,
            "eft1_wilson_rge_computed": False,
            "sequential_c5_computed": False,
        }
        save(target, report)
        return {"model": model, "status": report["status"],
                "operator_count": report["operator_term_count"], "output": str(target)}
    seed = read(source)
    rows = _validate_seed(seed, model)
    running_path = _find_model_file(model, r_roots, "running.json")
    required = {f["name"] for row in rows for f in row["factors"]}
    running, masses, ys, basis_diag = _running_basis(model, running_path, required)
    for row in rows:
        row["flavour_matrix_axes_external_indices"] = [f["external_index"] for f in row["factors"]]
        row["ordered_flavour_matrix_gev_inv"] = encode_matrix(flavour_matrix(row, masses, ys))
    payload = {
        "format": FORMAT, "model": model,
        "status": "SourceFactoredOrderedFlavourMatricesAtF_NotCanonicalWilsonTensor",
        "source": {**source_info, "running_json": str(running_path),
                   "running_sha256": digest(running_path)},
        "mu_F_gev": float(running["scales_gev"]["F"]),
        "heavy_masses_at_F_gev": [float(x) for x in masses],
        "heavy_flavour_mass_basis_diagnostics": basis_diag,
        "ordered_flavour_factor_count": len(rows),
        "operator_terms": rows,
        "source_CG_registry": seed.get("CGRegistry", []),
        "operator_template_definition": "Original Matchete TermInputForm with the two Yukawa factors and inverse MF set to 1. Its product with ordered_flavour_matrix_gev_inv reconstructs only the factorised flavour dependence, retaining every remaining symbolic field/CG/normalisation factor.",
        "full_general_yukawa_flavour_channels_retained": True,
        "external_lepton_order": "Matrix axes follow the two Yukawa factors in the exact Matchete term; see flavour_matrix_axes_external_indices. No unverified p/q or Majorana symmetrisation is applied.",
        "numeric_flavour_contractions_computed": True,
        "canonical_real_wilson_tensor_computed": False,
        "eft1_wilson_rge_computed": False,
        "scalar_threshold_matching_computed": False,
        "physical_matchete_rgbeta_dictionary_certified": False,
        "sequential_c5_computed": False,
    }
    save(target, payload)
    return {"model": model, "status": payload["status"],
            "operator_count": len(rows), "output": str(target)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument("--all", action="store_true")
    sel.add_argument("--model", choices=MODELS)
    sel.add_argument("--models", help="Comma-separated models")
    ap.add_argument("--first-threshold-roots", help="Comma-separated root folders, highest priority first")
    ap.add_argument("--running-roots", help="Comma-separated root folders, highest priority first")
    ap.add_argument("--output-root", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    models = MODELS if args.all else (args.model,) if args.model else tuple(
        x.strip() for x in args.models.split(",") if x.strip())
    if not models or any(m not in MODELS for m in models) or len(set(models)) != len(models):
        ap.error("Invalid or duplicate model selection")
    f_roots = _root_paths(args.first_threshold_roots, DEFAULT_F_ROOTS)
    r_roots = _root_paths(args.running_roots, DEFAULT_R_ROOTS)
    output = ROOT / args.output_root
    results = []
    for i, model in enumerate(models, 1):
        try:
            item = process(model, f_roots, r_roots, output)
        except Exception as exc:
            item = {"model": model, "status": "Blocked", "error": f"{type(exc).__name__}: {exc}"}
            save(output / tag(model) / "error.json", item)
        results.append(item)
        print(f"[{i}/{len(models)}] {model}: {item['status']}"
              + (" — " + item["error"] if "error" in item else ""), flush=True)
        save(output / "summary.json", {
            "format": FORMAT, "results": results,
            "source_factored_numerical_count": sum(r["status"].startswith("SourceFactored") for r in results),
            "mixed_operator_source_count": sum(r["status"].startswith("VerifiedMixed") for r in results),
            "blocked_count": sum(r["status"] == "Blocked" for r in results),
            "fully_sequential_C5_certified": False,
        })
    print(f"Numerical source-factorized Wilson operators: {sum(r['status'].startswith('SourceFactored') for r in results)}/{len(models)}; mixed X/Y source-only: {sum(r['status'].startswith('VerifiedMixed') for r in results)}; full EFT1 running: not performed")
    return 2 if any(r["status"] == "Blocked" for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
