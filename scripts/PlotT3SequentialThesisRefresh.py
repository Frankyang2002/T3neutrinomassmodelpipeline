"""Regenerate thesis Figures 5.2-5.3 and Appendix B Tables B.1-B.4.

Input: completed RunT3SequentialThesisRefresh --phase fit and --phase transfer.
Creates separately labelled D(-2)-shared / A(-4)-shared heatmaps, PDFs/PNGs,
CSV and LaTeX appendix tables. Does NOT modify historical assets, Figure 5.1,
Figure 5.4 (coupling sensitivity), or the group-theory tables.

Unlike historical Section 5.2, these plots use one fixed donor Yukawa pair
per 16-model campaign, not 16 separately fitted Yukawa pairs. The figure
captions and manifest mark that distinction explicitly.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import json
import math
import re

import numpy as np

from scripts.RunT3SequentialThesisRefresh import (
    DONORS, MODEL_ORDER, FORMAT, FIT_THRESHOLD, OUT, ROOT, jread, jwrite,
    csvwrite, extract_obs,
)

MIXING = ("sin2_theta12", "sin2_theta13", "sin2_theta23")
SPLITTINGS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2")
STAGES = ("UV", "EFT1", "SMEFT")
OUTPUT = ROOT / "Reports/output/thesis_sequential_refresh"


def model_label(model: str) -> str:
    letter, charge = model.split("_")
    number = int(charge)
    return f"{letter}({number:+d})" if number > 0 else f"{letter}({number})"


def model_short(model: str) -> str:
    return model.replace("-", "m").replace("+", "p")


def stage_map() -> dict:
    return {"UV": "U0_E1_S1", "EFT1": "U1_E0_S1", "SMEFT": "U1_E1_S0"}


def matrix(value: dict) -> np.ndarray:
    r = np.asarray(value["real"], dtype=float)
    i = np.asarray(value["imag"], dtype=float)
    if r.shape != (3,3) or i.shape != (3,3):
        raise ValueError("Expected 3x3 real/imag matrix")
    m = r + 1j*i
    if not np.isfinite(m).all():
        raise ValueError("Nonfinite matrix")
    return m


def alignment(ref: np.ndarray, candidate: np.ndarray) -> tuple[complex,float]:
    """C5(candidate) = k*C5(reference) + epsilon, Frobenius metric."""
    if ref.shape != (3,3) or candidate.shape != (3,3):
        raise ValueError("Bad C5 matrix")
    nr = float(np.linalg.norm(ref))
    nt = float(np.linalg.norm(candidate))
    if nr <= 1e-35 or nt <= 1e-35:
        raise ValueError("C5 is too small for stable relative alignment")
    k = np.vdot(ref, candidate) / np.vdot(ref, ref)
    eps = float(np.linalg.norm(candidate-k*ref) / nt)
    return complex(k), eps


def gather(src: Path, *, allow_poor: bool) -> dict:
    manifest = jread(src / "manifest.json")
    if manifest.get("format") != FORMAT:
        raise ValueError("Wrong sequential thesis campaign format")
    fits, results = {}, {}
    for donor in DONORS:
        fit = jread(src / "fits" / (model_short(donor)+".json"))
        if fit.get("model") != donor or fit.get("kind") != "donor_fit":
            raise ValueError("Wrong donor-fit provenance")
        if fit["status"] != "accepted_conditional_fit" and not allow_poor:
            raise ValueError(f"{donor} donor fit is poor (chi2={fit['final_chi2']}); need --include-poor")
        fits[donor] = fit
        results[donor] = {}
        for model in MODEL_ORDER:
            p = src/"transfers"/model_short(donor)/model_short(model)/"transfer.json"
            item = jread(p)
            if (item.get("format") != FORMAT or item.get("kind") != "shared_Yukawa_transfer"
                or item.get("model") != model or item.get("reference") != donor):
                raise ValueError(f"Incorrect transfer provenance: {p}")
            if item.get("extra_yukawa_policy") != manifest.get("extra_yukawa_policy"):
                raise ValueError(f"Extra Yukawa policy mismatch: {p}")
            for case in ("U1_E1_S1", "U0_E1_S1", "U1_E0_S1", "U1_E1_S0", "U0_E0_S1"):
                if case not in item["cases"]:
                    raise ValueError(f"Missing {case}: {p}")
            results[donor][model] = item
    return {"manifest": manifest, "fits": fits, "results": results}


def heatmap_data(results: dict, donor: str, observables: tuple[str, ...]) -> tuple[np.ndarray,list[str]]:
    columns = [f"{stage}\n{obs}" for stage in STAGES for obs in observables]
    raw = np.full((len(MODEL_ORDER),len(columns)),float("nan"))
    for row, model in enumerate(MODEL_ORDER):
        item = results[donor][model]
        full = item["cases"]["U1_E1_S1"]["observables"]
        for si, stage in enumerate(STAGES):
            frozen = item["cases"][stage_map()[stage]]["observables"]
            for oi, obs in enumerate(observables):
                if frozen[obs] != 0:
                    raw[row,si*len(observables)+oi] = 100*(full[obs]-frozen[obs])/frozen[obs]
    if not np.isfinite(raw).all():
        raise ValueError("Heatmap has nonfinite percentage changes")
    return raw, columns


def draw_heatmap(out: Path, dataset: dict, donor: str,
                 observables: tuple[str,...], fig_id: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    data, cols = heatmap_data(dataset["results"], donor, observables)
    ncol = len(cols)
    width = 12.0 if ncol == 9 else 9.5
    fig, ax = plt.subplots(figsize=(width,10.0), constrained_layout=True)
    maxabs = max(0.01,float(np.max(abs(data))))
    norm = TwoSlopeNorm(vmin=-maxabs,vcenter=0,vmax=maxabs)
    im = ax.imshow(data,aspect="auto",cmap="RdBu_r", norm=norm)
    ax.set_xticks(np.arange(ncol))
    symbols = {
        "sin2_theta12":r"$\sin^2\theta_{12}$",
        "sin2_theta13":r"$\sin^2\theta_{13}$",
        "sin2_theta23":r"$\sin^2\theta_{23}$",
        "delta_m21_sq_ev2":r"$\Delta m^2_{21}$",
        "delta_m3l_sq_ev2":r"$\Delta m^2_{3\ell}$",
    }
    ticks = [stage+"\n"+symbols[obs] for stage in STAGES for obs in observables]
    ax.set_xticklabels(ticks,fontsize=8)
    ax.set_yticks(np.arange(len(MODEL_ORDER)))
    ax.set_yticklabels([model_label(m) for m in MODEL_ORDER],fontsize=9)
    ax.set_ylabel("T3 model")
    ax.set_xlabel("Frozen stage (other stages running) and observable")
    for i in range(len(MODEL_ORDER)):
        for j in range(ncol):
            z = float(data[i,j])
            number = f"{z:.3g}" if abs(z)>=0.1 else f"{z:.2g}"
            ax.text(j,i,number,ha="center",va="center",fontsize=7,
                    color="white" if abs(z)>0.66*maxabs else "black")
    ax.set_xticks(np.arange(-.5,ncol,1), minor=True)
    ax.set_yticks(np.arange(-.5,len(MODEL_ORDER),1),minor=True)
    ax.grid(which="minor",color="black",linewidth=.3)
    ax.tick_params(which="minor",bottom=False,left=False)
    for x in (len(observables)-0.5,2*len(observables)-0.5):
        ax.axvline(x,color="black",linewidth=1.5)
    for y in (3.5,6.5,8.5,11.5):
        ax.axhline(y,color="black",linewidth=1.2)
    cbar = fig.colorbar(im,ax=ax,shrink=.8)
    cbar.set_label(r"Signed percentage change $100(X_{111}-X_{\mathrm{stage\ frozen}})/X_{\mathrm{stage\ frozen}}$")
    title = ("Figure 5.2 | PMNS mixing-angle running" if fig_id == "5_2"
             else "Figure 5.3 | Neutrino mass-squared running")
    ax.set_title(f"{title}\nNew sequential EFT | fixed {model_label(donor)} Yukawa reference",fontsize=12,pad=16)
    out.mkdir(parents=True,exist_ok=True)
    basename = f"fig{fig_id}_sequential_{model_short(donor)}"
    fig.savefig(out/(basename+".pdf"),bbox_inches="tight")
    fig.savefig(out/(basename+".png"),bbox_inches="tight",dpi=220)
    plt.close(fig)


def fig_alignment(out: Path, alig: list[dict]) -> None:
    """Additional thesis-relevant comparison: proportionality residual vs donor."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.arange(len(MODEL_ORDER))
    fig, ax = plt.subplots(figsize=(12,5.3),constrained_layout=True)
    width=0.37
    for idx,donor in enumerate(DONORS):
        d = {r["model"]:r for r in alig if r["reference"]==donor}
        y=np.asarray([d[m]["epsilon_c5"] for m in MODEL_ORDER],float)
        ax.bar(x+(idx-.5)*width,y,width,label=f"Fixed {model_label(donor)} Yukawas")
    ax.set_xticks(x)
    ax.set_xticklabels([model_label(m) for m in MODEL_ORDER],rotation=55,ha="right")
    ax.set_ylabel(r"Relative flavour residual $\epsilon = \Vert C_5-kC_{5,ref}\Vert_F/\Vert C_5\Vert_F$")
    ax.set_title("Low-scale Weinberg matrix proportionality across 16 T3 representations")
    ax.set_ylim(bottom=0)
    ax.legend()
    fig.savefig(out/"fig_sequential_c5_proportionality.pdf",bbox_inches="tight")
    fig.savefig(out/"fig_sequential_c5_proportionality.png",bbox_inches="tight",dpi=220)
    plt.close(fig)


def fig_fit_quality(out: Path, fits: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(6.4,4.2),constrained_layout=True)
    names=[model_label(m) for m in DONORS]
    values=[fits[m]["final_chi2"] for m in DONORS]
    ax.bar(names,values)
    ax.axhline(FIT_THRESHOLD,color="black",linestyle="--",label=f"Project screening cut ({FIT_THRESHOLD})")
    ax.set_yscale("symlog",linthresh=1e-4)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(r"Approximate oscillation $\chi^2$")
    ax.set_title("Two sequential-EFT Yukawa reference fits")
    ax.legend()
    fig.savefig(out/"fig_sequential_two_donor_fits.pdf",bbox_inches="tight")
    fig.savefig(out/"fig_sequential_two_donor_fits.png",bbox_inches="tight",dpi=220)
    plt.close(fig)


def output_tables(output: Path, dataset: dict) -> dict:
    res=dataset["results"]
    table1=[]; table2=[]; table3=[]; table4=[]; matrix_rows=[]
    for donor in DONORS:
        ref = matrix(res[donor][donor]["cases"]["U1_E1_S1"]["C5_at_muLow"])
        for model in MODEL_ORDER:
            item=res[donor][model]
            full=item["cases"]["U1_E1_S1"]
            obs=full["observables"]
            cand=matrix(full["C5_at_muLow"])
            k, eps=alignment(ref,cand)
            r={"reference":donor,"model":model,
               "m1_ev":obs["m1_ev"],"m2_ev":obs["m2_ev"],"m3_ev":obs["m3_ev"],
               "delta_m21_sq_1e_minus5_ev2":obs["delta_m21_sq_ev2"]*1e5,
               "delta_m3l_sq_1e_minus3_ev2":obs["delta_m3l_sq_ev2"]*1e3,
               "chi2":obs["chi2"],"epsilon_c5":eps}
            (table1 if donor=="D_-2" else table2).append(r)
            table3.append({"model":model,"reference":donor,
                           **{k:obs[k] for k in MIXING}})
            table4.append({"model":model,"reference":donor,
                           "k_real":float(k.real),"k_imag":float(k.imag),
                           "epsilon_c5":eps})
            matrix_rows.append({"reference":donor,"model":model,
                "k_real":float(k.real),"k_imag":float(k.imag),"epsilon_c5":eps})
    tables=output/"tables"
    csvwrite(tables/"table_B1_Dm2_observables.csv",table1)
    csvwrite(tables/"table_B2_Am4_observables.csv",table2)
    csvwrite(tables/"table_B3_pmns_shared.csv",table3)
    csvwrite(tables/"table_B4_c5_proportionality.csv",table4)
    write_tex_tables(tables,table1,table2,table3,table4)
    return {"table1":table1,"table2":table2,"table3":table3,"table4":table4,
            "matrix_rows":matrix_rows}


def write_tex_tables(out:Path, t1:list[dict],t2:list[dict],t3:list[dict],t4:list[dict]):
    """Standalone tabular snippets; existing appendix LaTeX is NOT rewritten."""
    lb=r"\\"
    for title, rows, name in (("D(-2)",t1,"table_B1_Dm2_observables.tex"),
                              ("A(-4)",t2,"table_B2_Am4_observables.tex")):
        lines=[r"\begin{tabular}{lrrrrrrr}",r"\hline",
               r"Model & $m_1$ & $m_2$ & $m_3$ & $\Delta m^2_{21}/10^{-5}$ & $\Delta m^2_{3\ell}/10^{-3}$ & $\chi^2$ & $\epsilon_{C_5}$ "+lb,
               r"\hline"]
        for r in rows:
            lines.append(f"{model_label(r['model'])} & {r['m1_ev']:.6g} & {r['m2_ev']:.6g} & {r['m3_ev']:.6g} & {r['delta_m21_sq_1e_minus5_ev2']:.6g} & {r['delta_m3l_sq_1e_minus3_ev2']:.6g} & {r['chi2']:.5g} & {r['epsilon_c5']:.4g} "+lb)
        lines.extend([r"\hline",r"\end{tabular}",f"% Conditional sequential soft-sector fit; {title} donor; fixed Y1/Y2; no independent hard matching."])
        (out/name).write_text("\n".join(lines)+"\n",encoding="utf-8")
    lines=[r"\begin{tabular}{llrrr}",r"\hline",
           r"Model & Yukawa reference & $\sin^2\theta_{12}$ & $\sin^2\theta_{13}$ & $\sin^2\theta_{23}$ "+lb,
           r"\hline"]
    for r in t3:
        lines.append(f"{model_label(r['model'])} & {model_label(r['reference'])} & {r['sin2_theta12']:.7g} & {r['sin2_theta13']:.7g} & {r['sin2_theta23']:.7g} "+lb)
    lines.extend([r"\hline",r"\end{tabular}"])
    (out/"table_B3_pmns_shared.tex").write_text("\n".join(lines)+"\n",encoding="utf-8")
    indexed={(r["reference"],r["model"]):r for r in t4}
    lines=[r"\begin{tabular}{lrrrrrr}",r"\hline",
           r"Model & $\Re k_D$ & $\Im k_D$ & $\epsilon_D$ & $\Re k_A$ & $\Im k_A$ & $\epsilon_A$ "+lb,
           r"\hline"]
    for model in MODEL_ORDER:
        d=indexed[("D_-2",model)];a=indexed[("A_-4",model)]
        lines.append(f"{model_label(model)} & {d['k_real']:.6g} & {d['k_imag']:.3g} & {d['epsilon_c5']:.5g} & {a['k_real']:.6g} & {a['k_imag']:.3g} & {a['epsilon_c5']:.5g} "+lb)
    lines.extend([r"\hline",r"\end{tabular}"])
    (out/"table_B4_c5_proportionality.tex").write_text("\n".join(lines)+"\n",encoding="utf-8")


def captions_file(out: Path, fit_status: dict) -> None:
    lines=[
        "# Sequential thesis figure captions (separate proposed replacements)",
        "",
        "All figures are CONDITIONAL leading-log + two-scalar finite-bubble results.",
        "Complete independent hard matching, GammaCC/Weyl normalization and mixed scalar basis identification are unverified.",
        "The original Section 5.2 heatmaps used individually optimized Yukawas for each model;",
        "the new alternatives deliberately transfer ONE fixed D(-2) or A(-4) pair to 16 models.",
        "SMEFT running means the SM + Weinberg operator from 1 TeV to 100 GeV.",
        "The UV/EFT1 'off' cases freeze only renormalisable couplings, NOT the fixed-order LLSS log.",
        "",
    ]
    for donor in DONORS:
        d=model_short(donor)
        lines.extend([
            f"## Reference {model_label(donor)} ({fit_status[donor]['status']}, chi2={fit_status[donor]['final_chi2']:.5g})",
            f"**Figure 5.2 alternative (`fig5_2_sequential_{d}.pdf`):** Percentage changes in sin²(theta12), sin²(theta13), sin²(theta23) from switching UV, scalar EFT1 or SMEFT running on, compared with the corresponding one-stage-frozen result. Same fixed {model_label(donor)}-optimised Y1/Y2 in all 16 models.",
            f"**Figure 5.3 alternative (`fig5_3_sequential_{d}.pdf`):** Analogous percentage changes in solar and atmospheric mass-squared splittings, using the identical transferred Yukawa pair and three stage comparisons.",
            "",
        ])
    lines += [
        "## Unchanged content",
        "- Figure 5.1 (Feynman diagrams): analytical; unchanged.",
        "- Figure 5.4 (coupling sensitivity): intentionally NOT regenerated.",
        "- Table 5.1/5.2 (group theory): unchanged; not recomputed from neutrino fits.",
        "",
        "## Appendix B",
        "The four regenerated tabulars B.1-B.4 are in tables/*.tex and tables/*.csv.",
        "B.4 is a complex Frobenius projection of C5 at the **low scale** (100 GeV).",
        "The transfer retains each target model's original saved non-Yukawa parameters.",
        "Additional model-specific Yukawa channels follow the campaign's explicit zero/saved policy.",
    ]
    (out/"captions_and_scope.md").write_text("\n".join(lines)+"\n",encoding="utf-8")


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-id",required=True)
    p.add_argument("--input-root",type=Path,default=OUT)
    p.add_argument("--output-root",type=Path,default=OUTPUT)
    p.add_argument("--include-poor",action="store_true",help="Explicitly draw provisional plots from poor-fit donors")
    args=p.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}",args.run_id):
        raise ValueError("Unsafe run-id")
    src=Path(args.input_root)/args.run_id
    dst=Path(args.output_root)/args.run_id
    if dst.exists():
        raise FileExistsError(f"Output exists; figures are write-once: {dst}")
    data=gather(src,allow_poor=args.include_poor)
    dst.mkdir(parents=True,exist_ok=False)
    tables=output_tables(dst,data)
    for donor in DONORS:
        draw_heatmap(dst/"figures",data,donor,MIXING,"5_2")
        draw_heatmap(dst/"figures",data,donor,SPLITTINGS,"5_3")
    fig_alignment(dst/"figures",tables["matrix_rows"])
    fig_fit_quality(dst/"figures",data["fits"])
    captions_file(dst,data["fits"])
    jwrite(dst/"figure_manifest.json",{
        "format":FORMAT,"run_id":args.run_id,"source_manifest":str(src/"manifest.json"),
        "donor_fit_status":{d:{"status":x["status"],"chi2":x["final_chi2"]} for d,x in data["fits"].items()},
        "thesis_figures_redone":["Figure 5.2 (two donors)","Figure 5.3 (two donors)"],
        "extra_figures":["C5 proportionality","two donor fit quality"],
        "appendix_tables_redone":["B.1","B.2","B.3","B.4"],
        "thesis_figures_intentionally_unchanged":["Figure 5.1","Figure 5.4 (sensitivity)"],
        "hard_matching_certified":False,"physical_normalization_certified":False,
        "new_figure_methodology":"fixed donor Y1/Y2 across 16 models, differs from original individually optimised Section 5.2",
    })
    print(f"Figures: {dst/'figures'}")
    print(f"Appendix tables: {dst/'tables'}")
    print(f"Captions / scope: {dst/'captions_and_scope.md'}")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
