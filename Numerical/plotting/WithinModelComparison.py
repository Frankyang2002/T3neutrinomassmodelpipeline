"""Build within-model comparisons of the four fixed T3 benchmark scenarios.

Presentation-only. This module reads existing ``running_diagnostics.json``
files below ``output/full/comparison`` and never recomputes matching or RGEs.

The scenario lambda value is the common reference value. The actual initialized
mixing quartic is representation-normalized,

    lambdaT3_eff = f_R * lambdaT3_ref,

with f_R = (1, 1, sqrt(3), 1, 1/sqrt(2)) for classes A--E.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "output" / "full" / "comparison"
DEFAULT_OUTPUT = PROJECT_ROOT / "Reports" / "output" / "full" / "within_model_comparison"

SCENARIOS = {
    "smallY_smallL": (0.005, 0.1),
    "smallY_largeL": (0.005, 1.0),
    "largeY_smallL": (0.5, 0.1),
    "largeY_largeL": (0.5, 1.0),
}
SCENARIO_ORDER = tuple(SCENARIOS)

MODEL_CLASS_BY_DIMS = {
    (1, 3, 2): "A",
    (2, 2, 1): "B",
    (2, 2, 3): "C",
    (3, 1, 2): "D",
    (3, 3, 2): "E",
}
LAMBDA_T3_NORMALISATION = {
    "A": 1.0,
    "B": 1.0,
    "C": math.sqrt(3.0),
    "D": 1.0,
    "E": 1.0 / math.sqrt(2.0),
}


def _load(path: Path) -> dict[str, Any]:
    payload=json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload,dict) or payload.get("status")!="Success":
        raise ValueError("diagnostics must be a successful JSON object")
    return payload


def _model_name(path: Path,payload: dict[str,Any]) -> str:
    rep=payload.get("benchmark_search",{}).get("model")
    if isinstance(rep,dict) and rep.get("model_key"):
        return str(rep["model_key"])
    for part in reversed(path.parts):
        if part.startswith("T3_dS1_"):
            return part
    return path.parent.parent.name


def _model_class(model: str, payload: dict[str, Any] | None = None) -> str:
    if payload is not None:
        rep=payload.get("benchmark_search",{}).get("model")
        if isinstance(rep,dict):
            try:
                dims=(int(rep["d_s1"]),int(rep["d_s2"]),int(rep["d_f"]))
            except (KeyError,TypeError,ValueError):
                pass
            else:
                if dims in MODEL_CLASS_BY_DIMS:
                    return MODEL_CLASS_BY_DIMS[dims]

    match=re.search(r"T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_",model)
    if match:
        dims=tuple(int(value) for value in match.groups())
        if dims in MODEL_CLASS_BY_DIMS:
            return MODEL_CLASS_BY_DIMS[dims]
    raise ValueError(f"Cannot determine T3 model class from {model!r}.")


def _lambda_values(
    model: str,
    scenario: str,
    payload: dict[str, Any] | None = None,
) -> tuple[float,float,float,str]:
    _,reference=SCENARIOS[scenario]
    model_class=_model_class(model,payload)
    factor=LAMBDA_T3_NORMALISATION[model_class]
    return reference,reference*factor,factor,model_class


def collect(input_root: Path) -> dict[str,dict[str,dict[str,Any]]]:
    """Return model -> scenario -> successful diagnostic payload."""
    result: dict[str,dict[str,dict[str,Any]]]={}
    root=Path(input_root)
    for scenario in SCENARIO_ORDER:
        directory=root/scenario
        if not directory.is_dir():
            continue
        for path in sorted(directory.rglob("running_diagnostics.json")):
            try:
                payload=_load(path)
            except (OSError,json.JSONDecodeError,ValueError):
                continue
            model=_model_name(path,payload)
            if scenario in result.setdefault(model,{}):
                raise RuntimeError(
                    f"Duplicate successful diagnostic for {model} / {scenario}"
                )
            result[model][scenario]=payload
    return result


def _save(path: Path) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    plt.tight_layout()
    plt.savefig(path,dpi=220,bbox_inches="tight")
    plt.close()


def _label(model: str, scenario: str, payload: dict[str,Any]) -> str:
    y,_=SCENARIOS[scenario]
    reference,effective,_,_=_lambda_values(model,scenario,payload)
    return (
        rf"$Y={y:g},\ \lambda_{{T3}}^{{ref}}={reference:g},"
        rf"\ \lambda_{{T3}}^{{eff}}={effective:.4g}$"
    )


def _plot_final_scalar(
    model: str,
    scenarios: dict[str,dict[str,Any]],
    *,
    key: str,
    ylabel: str,
    title: str,
    output: Path,
    absolute: bool=False,
    logy: bool=False,
) -> None:
    plt.figure(figsize=(7.8,4.9))
    plotted=False
    all_positive=True
    for scenario in SCENARIO_ORDER:
        payload=scenarios.get(scenario)
        if payload is None:
            continue
        block=payload["final_running"]
        mu=np.asarray(block["mu_gev"],dtype=float)
        values=np.asarray(block[key],dtype=float)
        if absolute:
            values=np.abs(values)
        if values.shape!=(mu.size,):
            raise ValueError(f"final_running.{key} has unexpected shape")
        all_positive &= bool(np.all(values>0.0))
        plt.plot(mu,values,label=_label(model,scenario,payload),linewidth=1.35)
        plotted=True
    if not plotted:
        plt.close()
        return
    plt.xscale("log")
    if logy and all_positive:
        plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True,alpha=0.25)
    plt.legend(fontsize=7.5)
    _save(output)


def _plot_c5_norm(
    model: str, scenarios: dict[str,dict[str,Any]], output: Path
) -> None:
    plt.figure(figsize=(7.8,4.9))
    plotted=False
    positive=True
    for scenario in SCENARIO_ORDER:
        payload=scenarios.get(scenario)
        if payload is None:
            continue
        block=payload["final_running"]
        mu=np.asarray(block["mu_gev"],dtype=float)
        c5=np.asarray(block["c5_abs"],dtype=float)
        if c5.shape!=(mu.size,3,3):
            raise ValueError("final_running.c5_abs has unexpected shape")
        norm=np.sqrt(np.sum(c5*c5,axis=(1,2)))
        positive &= bool(np.all(norm>0.0))
        plt.plot(mu,norm,label=_label(model,scenario,payload),linewidth=1.35)
        plotted=True
    if not plotted:
        plt.close()
        return
    plt.xscale("log")
    if positive:
        plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$\|C_5\|_F$ [GeV$^{-1}$]")
    plt.title(r"Within-model Weinberg-coefficient comparison")
    plt.grid(True,alpha=0.25)
    plt.legend(fontsize=7.5)
    _save(output)


def _plot_direct_c5_norm(
    model: str, scenarios: dict[str,dict[str,Any]], output: Path
) -> None:
    plt.figure(figsize=(7.8,4.9))
    plotted=False
    for scenario in SCENARIO_ORDER:
        payload=scenarios.get(scenario)
        if payload is None:
            continue
        block=payload.get("intermediate_direct_weinberg")
        if not isinstance(block,dict):
            continue
        mu=np.asarray(block["mu_gev"],dtype=float)
        c5=np.asarray(block["delta_c5_abs"],dtype=float)
        if c5.shape!=(mu.size,3,3):
            raise ValueError(
                "intermediate_direct_weinberg.delta_c5_abs has unexpected shape"
            )
        norm=np.sqrt(np.sum(c5*c5,axis=(1,2)))
        positive=norm>0.0
        if not np.any(positive):
            continue
        plt.plot(mu[positive],norm[positive],
                 label=_label(model,scenario,payload),linewidth=1.35)
        plotted=True
    if not plotted:
        plt.close()
        return
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Renormalisation scale $\mu$ [GeV]")
    plt.ylabel(r"$\|\Delta C_5^{\rm direct}\|_F$ [GeV$^{-1}$]")
    plt.title("Within-model direct intermediate Weinberg generation")
    plt.grid(True,alpha=0.25)
    plt.legend(fontsize=7.5)
    _save(output)


def _plot_mixing(
    model: str, scenarios: dict[str,dict[str,Any]], output: Path
) -> None:
    fig,axes=plt.subplots(3,1,figsize=(7.8,9.0),sharex=True)
    keys=(
        ("sin2_theta12",r"$\sin^2\theta_{12}$"),
        ("sin2_theta13",r"$\sin^2\theta_{13}$"),
        ("sin2_theta23",r"$\sin^2\theta_{23}$"),
    )
    plotted=False
    for scenario in SCENARIO_ORDER:
        payload=scenarios.get(scenario)
        if payload is None:
            continue
        block=payload["final_running"]
        mu=np.asarray(block["mu_gev"],dtype=float)
        for ax,(key,label) in zip(axes,keys):
            values=np.asarray(block[key],dtype=float)
            ax.plot(mu,values,label=_label(model,scenario,payload),linewidth=1.2)
            ax.set_ylabel(label)
            ax.grid(True,alpha=0.25)
        plotted=True
    if not plotted:
        plt.close(fig)
        return
    axes[-1].set_xscale("log")
    axes[-1].set_xlabel(r"Renormalisation scale $\mu$ [GeV]")
    axes[0].legend(fontsize=7.5,ncol=2)
    fig.suptitle("Within-model neutrino mixing comparison")
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.tight_layout(rect=(0,0,1,0.97))
    fig.savefig(output,dpi=220,bbox_inches="tight")
    plt.close(fig)


def write_model_report(
    model: str,
    scenarios: dict[str,dict[str,Any]],
    output_dir: Path,
) -> Path:
    output=Path(output_dir)/model
    output.mkdir(parents=True,exist_ok=True)

    _plot_direct_c5_norm(model,scenarios,output/"direct_weinberg_norm_comparison.png")
    _plot_c5_norm(model,scenarios,output/"c5_norm_comparison.png")
    _plot_final_scalar(
        model,scenarios,key="delta_m21_sq_ev2",
        ylabel=r"$\Delta m_{21}^2$ [eV$^2$]",
        title="Within-model solar splitting comparison",
        output=output/"dm21_comparison.png",logy=True,
    )
    _plot_final_scalar(
        model,scenarios,key="delta_m3l_sq_ev2",
        ylabel=r"$|\Delta m_{3\ell}^2|$ [eV$^2$]",
        title="Within-model atmospheric splitting comparison",
        output=output/"dm3l_comparison.png",absolute=True,logy=True,
    )
    _plot_mixing(model,scenarios,output/"mixing_comparison.png")

    lines=[
        f"# {model} — four-benchmark comparison","",
        "Presentation-only comparison of existing diagnostics. No physics was recomputed.","",
        "The scenario value is lambda_T3^ref. The actual initialized coupling is",
        "lambda_T3^eff = f_R lambda_T3^ref, using the verified direct",
        "LLSS -> Weinberg representation-factor normalization.","",
        "## Available scenarios","",
    ]
    for scenario in SCENARIO_ORDER:
        y,_=SCENARIOS[scenario]
        payload=scenarios.get(scenario)
        status="available" if payload is not None else "missing"
        if payload is not None:
            ref,eff,factor,model_class=_lambda_values(model,scenario,payload)
            lines.append(
                f"- `{scenario}`: Y={y:g}, lambda_T3^ref={ref:g}, "
                f"lambda_T3^eff={eff:.8g}, f_R={factor:.8g}, "
                f"class={model_class} — {status}"
            )
        else:
            ref=SCENARIOS[scenario][1]
            lines.append(
                f"- `{scenario}`: Y={y:g}, lambda_T3^ref={ref:g} — {status}"
            )
    lines += [
        "",
        "The C5 norm is the Frobenius norm reconstructed from the stored absolute",
        "matrix entries. The direct intermediate norm uses the stored full-flavor",
        "Delta C5^direct trajectory. These plots do not use the symbolic LLSS",
        "self-running diagnostic and do not modify the authoritative final C5.",
    ]
    (output/"README.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return output


def build_within_model_comparisons(
    input_root: Path=DEFAULT_INPUT,
    output_root: Path=DEFAULT_OUTPUT,
) -> Path:
    data=collect(Path(input_root))
    if not data:
        raise FileNotFoundError(
            f"No successful comparison diagnostics found under {input_root}"
        )
    output=Path(output_root)
    output.mkdir(parents=True,exist_ok=True)
    for model in sorted(data):
        write_model_report(model,data[model],output)
    print(f"Within-model comparisons written to: {output}")
    print(f"Models: {len(data)}")
    return output


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=DEFAULT_INPUT)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    args=parser.parse_args()
    build_within_model_comparisons(args.input,args.output)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
