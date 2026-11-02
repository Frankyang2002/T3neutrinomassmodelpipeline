"""Read-only, source-verified bridge from the all-16 fermion-first LLSS release.

The release contains a *formal commuting-spinor polynomial*, not a physically
normalised EFT1 Wilson tensor. This module preserves its exact orbit convention
and DOES NOT supply an implicit GammaCC/Weyl factor or EFT1 RG evolution.

An independently evidenced normalization certificate and an explicit same-scale
EFT1 endpoint are required before preparing inputs for joint scalar matching.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from .adapter import ALL_MODELS, CONVENTION, InputSet, normalize_model, scalar_groups
from .core import EigenBasis, check_tensors

SOURCE_FORMAT = "t3-lo-fermion-boundary-source-resolved-provisional-v1"
SOURCE_STATUS = "NumericalLLSSBySourceVerified_SpinorNormalizationOpen"
RELEASE_STATUS = "16of16NumericalLLSSTreeVerified_PhysicalNormalizationOpen"
ENDPOINT_FORMAT = "t3-joint-scalar-eft1-endpoint-v1"
BRIDGE_FORMAT = "t3-joint-scalar-source-bridge-v1"
_ORTHO_RTOL = 2e-10


def require(ok: bool, reason: str) -> None:
    if not bool(ok):
        raise ValueError("Joint-scalar fermion bridge: " + reason)


def read_json(path: Path) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    require(isinstance(data, dict), f"JSON object required: {path}")
    return data


def digest(path: Path) -> str:
    return sha256(Path(path).read_bytes()).hexdigest()


def slug(model: str) -> str:
    return normalize_model(model).replace("-", "m").replace("+", "p")


def canonical(coordinate: tuple[int, ...]) -> tuple[int, ...]:
    i, p, j, q, a, b = coordinate
    if (i, p) > (j, q):
        i, p, j, q = j, q, i, p
    return (i, p, j, q, min(a, b), max(a, b))


def orbit(coordinate: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    i, p, j, q, a, b = coordinate
    return tuple(sorted({
        (i,p,j,q,a,b), (i,p,j,q,b,a),
        (j,q,i,p,a,b), (j,q,i,p,b,a),
    }))


def _validate_labels(doc: dict, n: int) -> tuple[str, ...]:
    offsets = doc.get("scalar_field_offsets_one_based")
    dimensions = doc.get("scalar_field_complex_dimensions")
    require(isinstance(offsets, dict) and isinstance(dimensions, dict),
            "explicit scalar basis offsets/dimensions required")
    require(offsets.keys() == dimensions.keys(), "scalar offset/dimension labels differ")
    labels: list[str | None] = [None] * n
    for name, first_raw in offsets.items():
        require(isinstance(name,str) and name, "scalar label must be nonempty")
        first = int(first_raw)
        size = int(dimensions[name])
        require(size > 0, "scalar dimension must be positive")
        # Self-bilinear split-real multiplets are already canonical real fields.
        count = size if re.fullmatch(r"T3Scalar[12](Minus|Plus)", name) else 2*size
        require(5 <= first and first+count <= 5+n, f"scalar block {name} outside heavy range")
        for global_index in range(first, first+count):
            index = global_index - 5
            require(labels[index] is None, "overlapping heavy scalar basis blocks")
            labels[index] = name
    require(all(label is not None for label in labels),
            "heavy scalar blocks do not cover every real component")
    return tuple(str(x) for x in labels)


def _row_coordinates(rows: list[dict], n: int, label: str) -> dict[tuple[int,...], complex]:
    require(isinstance(rows,list), f"{label}: coordinate list required")
    result: dict[tuple[int,...], complex] = {}
    for item in rows:
        require(isinstance(item,dict) and set(("coordinate","coefficient_real","coefficient_imag","orbit_size")) <= set(item),
                f"{label}: incomplete coordinate")
        idx = item["coordinate"]
        require(isinstance(idx,list) and len(idx)==6 and all(isinstance(k,int) and not isinstance(k,bool) for k in idx),
                f"{label}: six integer coordinates required")
        key = tuple(idx)
        i,p,j,q,a,b = key
        require(i in (1,2) and j in (1,2) and p in (1,2,3) and q in (1,2,3),
                f"{label}: lepton/generation index out of range")
        require(5<=a<5+n and 5<=b<5+n, f"{label}: scalar index outside heavy region")
        require(key == canonical(key), f"{label}: noncanonical representative")
        require(item["orbit_size"] == len(orbit(key)), f"{label}: orbit size mismatch")
        require(key not in result, f"{label}: duplicate independent coordinate")
        z = complex(float(item["coefficient_real"]),float(item["coefficient_imag"]))
        require(math.isfinite(z.real) and math.isfinite(z.imag),f"{label}: nonfinite coefficient")
        result[key] = z
    return result


def _dense(rows: dict[tuple[int,...],complex], n: int) -> np.ndarray:
    tensor = np.zeros((3,3,2,2,n,n),dtype=np.complex128)
    for key, value in rows.items():
        for i,p,j,q,a,b in orbit(key):
            tensor[p-1,q-1,i-1,j-1,a-5,b-5] = value
    # Use the combined lepton SU(2)+flavour swap, not flavour alone.
    require(np.allclose(tensor, tensor.transpose(1,0,3,2,4,5),rtol=1e-12,atol=1e-30),
            "lepton exchange symmetry violated")
    require(np.allclose(tensor,tensor.swapaxes(4,5),rtol=1e-12,atol=1e-30),
            "scalar exchange symmetry violated")
    return tensor


@dataclass(frozen=True)
class SourceBoundary:
    model: str
    path: Path
    hash: str
    mu_f_gev: float
    raw: np.ndarray
    by_source: np.ndarray
    source_indices: tuple[int,...]
    labels: tuple[str,...]
    normalization_certified: bool
    source_status: str


def load_source(path: Path, *, expected_model: str | None = None) -> SourceBoundary:
    path = Path(path)
    doc = read_json(path)
    require(doc.get("format")==SOURCE_FORMAT and doc.get("status")==SOURCE_STATUS,
            "wrong source release format/status; no silent migration")
    model = normalize_model(doc.get("model", ""))
    if expected_model is not None:
        require(model==normalize_model(expected_model), "source model mismatch")
    n = 2*{"A":4,"B":4,"C":4,"D":4,"E":6}[model[0]]
    require(doc.get("real_scalar_dimension")==4+n,
            "release dimension must include four Higgs coordinates and all heavy scalars")
    require(doc.get("index_base")==1, "expected one-based Matchete coordinates")
    require(doc.get("index_order")==["lepton_su2_i","lepton_flavour_p","lepton_su2_j",
                                       "lepton_flavour_q","real_scalar_a","real_scalar_b"],
            "release tensor axis ordering changed")
    require(doc.get("units")=="GeV^-1", "wrong Wilson mass dimension")
    require(doc.get("formal_polynomial")==
            "sum C[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b]",
            "release polynomial convention changed")
    require(doc.get("symmetric_tensor_orbit_rule")==
            "All distinct i,p<->j,q and a<->b orbit members equal representative. No extra orbit factor.",
            "formal polynomial orbit rule changed")
    fields = doc.get("field_content",{})
    require(fields.get("heavy_fermions_integrated_out")=="all at mu_F" and
            fields.get("heavy_scalar_multiplets")=="both retained in EFT1",
            "incorrect EFT1 field content")
    flags = doc.get("completeness",{})
    require(flags.get("F_tree_LLSS_included") is True and
            flags.get("distinct_PL_sources_retained") is True and
            flags.get("PR_Hermitian_conjugate_not_double_counted") is True,
            "incomplete PL/PR boundary handoff")
    require(flags.get("EFT1_running_performed") is False and
            flags.get("mu_S_scalar_matching_performed") is False,
            "source is not a muF-only boundary")
    mf = float(doc.get("mu_F_gev",float("nan")))
    require(math.isfinite(mf) and mf>0,"invalid muF")
    labels = _validate_labels(doc,n)
    summed = _row_coordinates(doc.get("summed_independent_symmetric_coordinates"),n,"summed")
    source_rows: list[np.ndarray] = []
    source_indices: list[int] = []
    sources=doc.get("sources")
    require(isinstance(sources,list) and len(sources)>0,"missing source-resolved PL operators")
    for src in sources:
        require(src.get("native_chirality")=="PL","PR or unknown operator in source list")
        number = src.get("matchete_source_index")
        require(isinstance(number,int) and number not in source_indices,"duplicate/invalid Matchete source index")
        source_indices.append(number)
        parsed = _row_coordinates(src.get("independent_symmetric_coordinates"),n,f"PL {number}")
        source_rows.append(_dense(parsed,n))
    summed_dense = _dense(summed,n)
    by_source = np.stack(source_rows,axis=0)
    total = np.sum(by_source,axis=0)
    # Relative norm against the true summed amplitudes, with a small absolute floor
    # for cancellation-protected entries. Never discard a small source term.
    scale = max(np.max(np.abs(total)),np.max(np.abs(summed_dense)),1e-30)
    discrepancy = float(np.max(abs(total-summed_dense)))
    require(discrepancy <= 1e-9*scale + 1e-28,
            f"source sum mismatches recorded total (max |difference|={discrepancy:.3e})")
    norm = doc.get("spinor_normalization",{})
    require(isinstance(norm,dict) and norm.get("certified") is False and
            norm.get("gamma_cc_to_master_weyl_factor") is None and
            norm.get("gamma_cc_to_master_weyl_sign") is None,
            "provisional release must explicitly leave GammaCC normalization unresolved")
    return SourceBoundary(model,path,digest(path),mf,summed_dense,by_source,
                          tuple(source_indices),labels,bool(norm["certified"]),doc["status"])


def load_release_manifest(release_dir: Path) -> dict[str,Path]:
    release_dir=Path(release_dir)
    report=read_json(release_dir/"summary.json")
    require(report.get("format")==SOURCE_FORMAT and report.get("status")==RELEASE_STATUS,
            "unsupported fermion-boundary release")
    require(report.get("n_models_verified")==16 and
            report.get("absolute_spinor_normalization_certified") is False,
            "source summary is not the audited provisional 16-model release")
    require(report.get("physical_eft1_workers_must_reject_this_as_certified_input") is True,
            "source summary unexpectedly changed its normalization gate")
    rows=report.get("models")
    require(isinstance(rows,list) and len(rows)==16,"expected exactly 16 model entries")
    out={}
    for item in rows:
        model=normalize_model(item["model"])
        require(model not in out,"duplicate model in source manifest")
        expected=f"{slug(model)}/boundary.json"
        require(item.get("path")==expected,"unexpected model path in source manifest")
        require(item.get("status")==SOURCE_STATUS,"source manifest model failed")
        out[model]=release_dir/expected
    require(set(out)==set(ALL_MODELS),"model inventory not all 16")
    return out


def write_provisional_source(target: Path, source: SourceBoundary) -> None:
    """Write only uncertified muF tensors; NEVER an input for physical C5."""
    target=Path(target)
    target.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(target,
         C_LLSS_formal_muF=source.raw,
         C_LLSS_by_PL_source_formal_muF=source.by_source,
         matchete_PL_source_indices=np.array(source.source_indices,dtype=np.int64))
    metadata={"format":BRIDGE_FORMAT,"model":source.model,
              "status":"FormalLLSSExported_PhysicalNormalizationUncertified",
              "source_boundary_path":str(source.path),"source_boundary_sha256":source.hash,
              "mu_f_gev":source.mu_f_gev,"n_heavy_real":len(source.labels),
              "source_scalar_real_labels":list(source.labels),
              "PL_source_indices":list(source.source_indices),
              "gamma_cc_weyl_normalization_certified":source.normalization_certified,
              "eft1_muS_running_included":False,
              "source_field_order":"1-based Higgs 1..4; heavy real 5..(4+n)",
              "source_tensor_shape":"(flavor_p,flavor_q,lepton_su2_i,lepton_su2_j,heavy_a,heavy_b)",
              "full_physical_C5_computed":False,
              "note":"No GammaCC factor or EFT1 RG evolution is applied. Only a formal Matchete-source polynomial."}
    target.with_suffix(".json").write_text(json.dumps(metadata,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def _load_certificate(meta: dict, endpoint_path: Path, source: SourceBoundary) -> float:
    cert=meta.get("normalization_certificate")
    require(isinstance(cert,dict),"independent GammaCC/Weyl normalization certificate missing")
    require(cert.get("status")=="independent_signed_amplitude_verified",
            "GammaCC/Weyl sign and factor have not been independently verified")
    require(cert.get("source_boundary_sha256")==source.hash,
            "normalization certificate is not tied to this source release")
    require(cert.get("master_operator_convention")==CONVENTION,
            "spinor factor pertains to a different Wilson normalization")
    factor = cert.get("gamma_cc_to_master_weyl_signed_factor")
    require(isinstance(factor,(int,float)) and not isinstance(factor,bool) and
            math.isfinite(float(factor)) and factor!=0,
            "explicit signed GammaCC/Weyl factor required; never assume +/-1")
    name=cert.get("evidence_file")
    require(isinstance(name,str) and Path(name).name==name and name and name not in (".",".."),
            "normalization evidence must be an adjacent file")
    p=endpoint_path.parent/name
    require(p.is_file(),"normalization evidence file absent")
    require(cert.get("evidence_sha256")==digest(p),"normalization evidence hash mismatch")
    return float(factor)


def _strict_tensors_at_muS(ep: dict, arrays: dict, source: SourceBoundary) -> InputSet:
    meta_path=ep["__path"]
    require(ep.get("format")==ENDPOINT_FORMAT,"unsupported EFT1 muS endpoint format")
    require(normalize_model(ep.get("model",""))==source.model,"endpoint model mismatch")
    require(ep.get("source_boundary_sha256")==source.hash,"endpoint boundary hash mismatch")
    require(ep.get("source_scalar_real_labels")==list(source.labels),
            "endpoint source field ordering differs; explicit map required")
    require(ep.get("renormalization_scheme")=="MSbar", "EFT1 endpoint must use MSbar")
    require(ep.get("tensor_convention")==CONVENTION,"mismatched Weinberg normalization")
    mf=float(ep.get("mu_f_gev",float("nan")))
    ms=float(ep.get("mu_s_gev",float("nan")))
    require(math.isfinite(mf) and math.isfinite(ms) and 0<ms<=mf and
            abs(mf/source.mu_f_gev-1)<1e-11,"endpoint scale mismatch")
    mass_origin=ep.get("scalar_mass_origin")
    require(mass_origin in ("unbroken_electroweak","EWSB_resummed"),
            "real endpoint must declare scalar mass/EFT truncation")
    ev=ep.get("eft1_evolution",{})
    require(isinstance(ev,dict) and ev.get("status")=="renormalized_full_EFT1_to_muS" and
            ev.get("same_muS_wilson_quartic_mass_C5") is True and
            ev.get("no_separate_scalar_threshold") is True and
            ev.get("C_LLSS_used_in_bubble_order")=="hbar0" and
            ev.get("C5_EFT1_includes_running_once") is True,
            "missing certified, same-scale EFT1 endpoint/power-counting handoff")
    factor=_load_certificate(ep,meta_path,source)
    n=len(source.labels)
    required={"source_to_matching_real_U","C_LLSS_muF_matching",
              "C_LLSS_muS_hbar0","lambda_HHSS_muS","mass_squared_gev2_muS","C5_EFT1_muS"}
    require(required.issubset(arrays),f"endpoint NPZ missing {sorted(required-set(arrays))}")
    U=np.asarray(arrays["source_to_matching_real_U"])
    require(U.shape==(n,n) and np.isrealobj(U) and np.isfinite(U).all(),
            "explicit source-to-matching real orthogonal basis map required")
    require(np.allclose(U.T@U,np.eye(n),rtol=_ORTHO_RTOL,atol=_ORTHO_RTOL),
            "source-to-matching transformation not orthogonal")
    # s_source = U @ s_matching, so C_matching = U.T @ C_source @ U.
    expected=np.einsum("fgijab,ak,bl->fgijkl",factor*source.raw,U,U,optimize=True)
    boundary=np.asarray(arrays["C_LLSS_muF_matching"],dtype=np.complex128)
    require(boundary.shape==expected.shape and np.isfinite(boundary).all(),
            "muF tensor has wrong shape/nonfinite entries")
    scale=max(float(np.max(abs(expected))),1e-30)
    error=float(np.max(abs(boundary-expected)))
    require(error<=1e-9*scale+1e-28,
            f"muF tensor differs from signed source boundary after basis rotation ({error:.3e})")
    C=np.asarray(arrays["C_LLSS_muS_hbar0"],dtype=complex)
    q=np.asarray(arrays["lambda_HHSS_muS"])
    mass2=np.asarray(arrays["mass_squared_gev2_muS"])
    c5=np.asarray(arrays["C5_EFT1_muS"],complex)
    check_tensors(C,q,n)
    EigenBasis.from_mass2(mass2)
    require(c5.shape==(3,3) and np.isfinite(c5).all() and
            np.allclose(c5,c5.T,rtol=1e-10,atol=1e-27),
            "inherited C5_EFT1(muS) must be finite, complex symmetric 3x3")
    channels=arrays.get("lambda_HHSS_channels_muS")
    names=tuple(ep.get("quartic_channel_names",[]))
    require(channels is not None and isinstance(names,tuple) and len(names)>0,
            "all independent quartic channels required for physical-input matching")
    channels=np.asarray(channels)
    require(channels.ndim==5 and channels.shape[1:]==q.shape and len(names)==len(channels)
            and len(set(names))==len(names) and all(isinstance(x,str) and x for x in names),
            "quartic channel labels/shapes mismatch")
    for ch in channels:
        check_tensors(C,ch,n)
    require(np.allclose(np.sum(channels,axis=0),q,rtol=1e-10,atol=1e-12*max(np.max(abs(q)),1e-30)),
            "independent quartic channels do not reconstruct complete portal tensor")
    ref=float(ep.get("reference_scale_gev",ms))
    require(math.isfinite(ref) and ref>0,"invalid matching reference scale")
    return InputSet(source.model,mf,ms,ref,mass2,C,q,c5,
       {"source":"existing fermion-first release + externally evidenced EFT1 muS endpoint",
        "source_boundary_file_sha256":source.hash,"signed_spinor_factor":factor,
        "external_normalization_evidence_sha256":ep["normalization_certificate"]["evidence_sha256"],
        "endpoint_npz_sha256":ep["__npz_sha256"],
        "independent_full_scalar_threshold_amplitude_certified":False},
       False,mass_origin,channels,names)


def load_eft1_endpoint(path: Path, source: SourceBoundary) -> InputSet:
    path=Path(path)
    meta=read_json(path)
    name=meta.get("tensor_npz")
    require(isinstance(name,str) and name and Path(name).name==name and name not in (".",".."),
            "endpoint tensor NPZ must be adjacent to JSON")
    npz_path=path.parent/name
    require(npz_path.is_file(),f"endpoint NPZ not found: {npz_path}")
    with np.load(npz_path,allow_pickle=False) as z:
        arr={key:np.array(z[key]) for key in z.files}
    meta["__path"]=path
    meta["__npz_sha256"]=digest(npz_path)
    return _strict_tensors_at_muS(meta,arr,source)
