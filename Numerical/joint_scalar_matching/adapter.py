"""Strict Chat-1 -> simultaneous scalar matcher NPZ+JSON adapter.

This layer deliberately does not guess a Matchete CG normalization, infer a
quartic from its coupling name, or turn a scalar mixing angle into a tensor.
Chat 1 supplies the explicit canonical real-component tensors, with all
independent SU(2) contractions already included.
"""
from __future__ import annotations
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .core import check_tensors, EigenBasis

CLASS_DIMS = {"A":(1,3,2),"B":(2,2,1),"C":(2,2,3),
              "D":(3,1,2),"E":(3,3,2)}
MODEL_ALPHAS = {"A":(-4,-2,0),"B":(-3,-1,1),
                "C":(-3,-1,1),"D":(-2,0,2),"E":(-4,-2,0,2)}
ALL_MODELS = tuple(f"{letter}_{alpha:+d}" if alpha>=0 else f"{letter}_{alpha}"
                   for letter, nums in MODEL_ALPHAS.items() for alpha in nums)
SELF_SPLIT = {"A_+0","A_-2","D_+0","D_-2","E_+0","E_-2"}
MIXED = {"B_-1","C_-1"}
SCHEMA_VERSION = 1
CONVENTION = "T3_REAL_LLSS_HHSS_WEINBERG_v1"


def normalize_model(model: str) -> str:
    """Map A_p0/A_m2, T3-A, T3_dS..., A_0 and A_-2 to A_+0/A_-2."""
    model = model.strip().replace("T3-", "")
    if model.startswith("T3_") and "_alpha_" in model:
        import re
        m = re.fullmatch(r"T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp]\d+)", model)
        if m:
            dims = tuple(int(x) for x in m.groups()[:3]);
            letter = next((key for key,val in CLASS_DIMS.items() if val==dims), None)
            if letter is None: raise ValueError("Unsupported model representation")
            model = letter+"_"+m.group(4)
    import re
    m = re.fullmatch(r"([A-E])_([mp]\d+|[+-]?\d+)", model)
    if not m: raise ValueError(f"Unknown T3 model {model!r}")
    x = m.group(2)
    a = (-int(x[1:]) if x.startswith("m") else int(x[1:]) if x.startswith("p") else int(x))
    label = f"{m.group(1)}_{a:+d}" if a >=0 else f"{m.group(1)}_{a}"
    if label not in ALL_MODELS: raise ValueError(f"Not one of the 16 models: {label}")
    return label


def scalar_groups(model: str) -> list[str]:
    model = normalize_model(model)
    d1,d2,_ = CLASS_DIMS[model[0]]
    return ["S1"]*(2*d1) + ["S2"]*(2*d2)


@dataclass(frozen=True)
class InputSet:
    model: str
    mu_f_gev: float
    mu_s_gev: float
    reference_scale_gev: float
    mass_squared: np.ndarray
    c_llss: np.ndarray
    lambda_hhss: np.ndarray
    c5_eft1: np.ndarray
    provenance: dict
    is_fixture: bool
    scalar_mass_origin: str
    quartic_channels: np.ndarray | None = None
    quartic_channel_names: tuple[str, ...] = ()


def load_input(metadata_path: Path) -> InputSet:
    """Read a metadata JSON and NPZ; no conversion from unknown tensor conventions."""
    path = Path(metadata_path)
    meta = json.loads(path.read_text(encoding="utf-8"))
    if meta.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported input schema_version")
    if meta.get("tensor_convention") != CONVENTION:
        raise ValueError("Uncertified tensor convention; explicit Chat-1 real-basis adapter required")
    if meta.get("renormalization_scheme") != "MSbar":
        raise ValueError("Only MSbar input supported, no implicit scheme conversion")
    model = normalize_model(meta["model"])
    mf,ms = float(meta["mu_f_gev"]),float(meta["mu_s_gev"])
    if not(np.isfinite([mf,ms]).all() and mf >= ms > 0):
        raise ValueError("Expect single scalar threshold 0 < mu_S <= mu_F")
    ref = float(meta.get("reference_scale_gev",ms))
    if not (np.isfinite(ref) and ref>0):raise ValueError("Invalid reference scale")
    mass_origin = meta.get("scalar_mass_origin")
    if mass_origin not in ("unbroken_electroweak", "EWSB_resummed", "synthetic_algebra_only"):
        raise ValueError("Must state scalar_mass_origin and EFT truncation implications")
    if mass_origin == "synthetic_algebra_only" and not bool(meta.get("synthetic_fixture",False)):
        raise ValueError("synthetic_algebra_only is valid only for explicitly marked fixtures")
    filename = meta.get("tensor_npz")
    if not filename or Path(filename).is_absolute() or len(Path(filename).parts)!=1:
        raise ValueError("tensor_npz must be a sibling basename")
    with np.load(path.parent/filename,allow_pickle=False) as z:
        required={"mass_squared_gev2","C_LLSS","lambda_HHSS","C5_EFT1"}
        if not required.issubset(z.files):
            raise ValueError(f"NPZ is missing {sorted(required-set(z.files))}")
        mass2 = np.array(z["mass_squared_gev2"])
        C = np.array(z["C_LLSS"])
        q = np.array(z["lambda_HHSS"])
        direct = np.array(z["C5_EFT1"])
        channels = np.array(z["lambda_HHSS_channels"]) if "lambda_HHSS_channels" in z.files else None
    groups = scalar_groups(model)
    if mass2.shape != (len(groups), len(groups)):
        raise ValueError("Scalar basis size incompatible with the model dimensions")
    EigenBasis.from_mass2(mass2)
    check_tensors(C,q,len(groups))
    if direct.shape!=(3,3) or not np.allclose(direct,direct.T,rtol=1e-10,atol=1e-28):
        raise ValueError("C5_EFT1 must be complex symmetric 3x3")
    provenance = meta.get("provenance")
    if not isinstance(provenance,dict) or "source" not in provenance:
        raise ValueError("Explicit provenance.source required")
    channel_names = tuple(meta.get("quartic_channel_names", []))
    if channels is not None:
        if channels.ndim != 5 or channels.shape[1:] != q.shape:
            raise ValueError("lambda_HHSS_channels shape should be (k,4,4,n,n)")
        if len(channel_names) != len(channels) or len(set(channel_names))!=len(channel_names) or not all(isinstance(x,str) and x.strip() for x in channel_names):
            raise ValueError("Every independent quartic channel needs a unique nonempty name")
        for tensor in channels:
            check_tensors(C,tensor,len(groups))
        if not np.allclose(np.sum(channels,axis=0),q,rtol=1e-11,
                           atol=1e-13*max(np.max(abs(q)),1e-300)):
            raise ValueError("Independent quartic channels do not sum to lambda_HHSS")
    elif channel_names:
        raise ValueError("quartic_channel_names provided without lambda_HHSS_channels")
    is_fixture = bool(meta.get("synthetic_fixture",False))
    return InputSet(model,mf,ms,ref,mass2,C,q,direct,provenance,is_fixture,mass_origin,channels,channel_names)


def save_input(destination: Path, inp: InputSet) -> Path:
    path=Path(destination)
    path.parent.mkdir(parents=True,exist_ok=True)
    npz_name=path.stem+".npz"
    tensor_payload={"mass_squared_gev2":inp.mass_squared,
                    "C_LLSS":inp.c_llss,"lambda_HHSS":inp.lambda_hhss,
                    "C5_EFT1":inp.c5_eft1}
    if inp.quartic_channels is not None:
        tensor_payload["lambda_HHSS_channels"]=inp.quartic_channels
    np.savez_compressed(path.parent/npz_name,**tensor_payload)
    meta={"schema_version":SCHEMA_VERSION,"model":inp.model,
          "tensor_convention":CONVENTION,"renormalization_scheme":"MSbar",
          "mu_f_gev":inp.mu_f_gev,"mu_s_gev":inp.mu_s_gev,
          "reference_scale_gev":inp.reference_scale_gev,
          "scalar_mass_origin":inp.scalar_mass_origin,"tensor_npz":npz_name,
          "synthetic_fixture":inp.is_fixture,"provenance":inp.provenance,
          "quartic_channel_names":list(inp.quartic_channel_names)}
    path.write_text(json.dumps(meta,indent=2)+"\n",encoding="utf-8")
    return path
