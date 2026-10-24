"""Read-only T3-B-m1 Yukawa-basis and C12 texture diagnostic.

Reads a single fitted configuration, interpreted as input at the *fermion*
threshold. No UV running, EFT evolution, matching, fitting or symbolic stages
are called. Its C12 is the tree-level heavy-fermion matching texture, not the
running C12(mu_S) or an independent low-energy neutrino prediction.

A full unitary SVD of Y1 is a mathematical *diagnostic* basis: it generally
rotates charged-lepton and heavy-Majorana mass matrices away from diagonal.
The separate texture-angle calculation retains the physical charged-lepton
basis and the positive diagonal heavy-mass basis of the original input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from html import escape
from pathlib import Path

import numpy as np

DEFAULT_CONFIG = Path("Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json")
DEFAULT_OUT = Path("output/studies/t3b_yukawa_basis/fitted_reference_frozen_boundary")


def _load(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return obj


def _matrix(value, key: str) -> np.ndarray:
    a = np.asarray(value["real"], float) + 1j * np.asarray(value["imag"], float)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError(f"{key}: finite complex (3,3) matrix required")
    return a


def _json_matrix(matrix: np.ndarray) -> dict:
    arr = np.asarray(matrix, complex)
    return {"real": arr.real.tolist(), "imag": arr.imag.tolist(),
            "frobenius_norm": float(np.linalg.norm(arr))}


def _shape(reference: np.ndarray, candidate: np.ndarray) -> dict:
    rr = np.vdot(reference, reference)
    cc = float(np.linalg.norm(candidate))
    if abs(rr) == 0 or cc == 0:
        raise ValueError("Cannot compare zero matrices")
    alpha = np.vdot(reference, candidate) / rr
    return {"alpha_real": float(alpha.real), "alpha_imag": float(alpha.imag),
            "residual_relative_to_candidate": float(np.linalg.norm(candidate-alpha*reference)/cc),
            "perpendicular_relative_to_reference": float(np.linalg.norm(candidate-alpha*reference)/np.linalg.norm(reference))}


def _heavy_basis(y1, y2, mf):
    """Use positive diagonal Majorana MF; refuse an unsafe arbitrary rotation.

    The original T3-B fitted configurations have exactly this form. For a
    non-diagonal MF use the repository's validated Majorana Takagi machinery
    separately and pass the already-rotated Yukawas/masses to a fitting study.
    """
    scale = max(np.linalg.norm(mf), 1e-30)
    if np.linalg.norm(mf-mf.T) / scale > 1e-10:
        raise ValueError("Majorana MF is not symmetric")
    offdiag = mf-np.diag(np.diag(mf))
    if np.linalg.norm(offdiag)/scale > 1e-12:
        raise ValueError("MF is not diagonal. This read-only diagnostic refuses a nontrivial heavy Takagi basis; use the project's validated mass-basis converter.")
    masses_complex = np.diag(mf)
    masses = np.abs(masses_complex)
    if np.any(masses <= 0) or not np.isfinite(masses).all():
        raise ValueError("Heavy masses must be finite positive")
    phase = np.exp(-0.5j*np.angle(masses_complex))
    return y1*phase[np.newaxis, :], y2*phase[np.newaxis, :], masses


def _texture(y1, y2, masses) -> tuple[np.ndarray, np.ndarray]:
    """Matched C12_r=(conj y1_r conj y2_r^T + swap)/(2 M_r)."""
    terms = []
    for r, mr in enumerate(masses):
        u = np.conj(y1[:, r]); v = np.conj(y2[:, r])
        terms.append((np.outer(u, v) + np.outer(v, u))/(2*mr))
    pieces = np.asarray(terms)
    total = pieces.sum(axis=0)
    if np.linalg.norm(total-total.T) > 1e-11 * max(np.linalg.norm(total),1e-50):
        raise RuntimeError("C12 texture lost Majorana symmetry")
    return pieces, total


def _charged_basis(matrix, ye):
    # Same left-lepton SVD basis convention as physics/NeutrinoTrajectory.py.
    left, singular, _ = np.linalg.svd(ye)
    left = left[:, np.argsort(singular)]
    return left.T @ matrix @ left


def _texture_angles(c12, ye) -> dict:
    """Takagi singular-vector proxy from H = M^dagger M.

    Absolute scale is immaterial. Assumes normal-ordering index convention,
    as in the fitted T3-B example. Reject numerical degeneracy.
    """
    m = _charged_basis(c12, ye)
    h = m.conj().T @ m
    w, u = np.linalg.eigh(h)
    scale = max(float(w[-1]), 1e-300)
    if min(w[1]-w[0],w[2]-w[1]) < 1e-10*scale:
        raise ValueError("Texture singular values nearly degenerate: angles basis-ambiguous")
    if np.any(w < -1e-10*scale):
        raise ValueError("M^dagger M negative beyond numerical rounding")
    a = np.abs(u)
    s13 = np.clip(a[0,2]**2,0,1)
    c13 = 1-s13
    if c13 < 1e-12:
        raise ValueError("theta13 near pi/2")
    s12 = np.clip(a[0,1]**2/c13,0,1)
    s23 = np.clip(a[1,2]**2/c13,0,1)
    return {"theta12_deg": float(np.degrees(np.arcsin(np.sqrt(s12)))),
            "theta13_deg": float(np.degrees(np.arcsin(np.sqrt(s13)))),
            "theta23_deg": float(np.degrees(np.arcsin(np.sqrt(s23)))),
            "pmns_abs": a.tolist(),
            "squared_singular_values_over_largest": (w/scale).tolist()}


def _figures(output, y1, y2_rot, components, total, variation):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1,2,figsize=(9.5,4),layout="constrained")
    mats = [(np.diag(np.linalg.svd(y1, compute_uv=False))/np.linalg.norm(y1), "D1 / ||Y1|| in Y1 SVD basis"),
            (np.abs(y2_rot)/np.linalg.norm(y2_rot), "|U1† Y2 V1|/||Y2||")]
    for ax,(data,label) in zip(axes,mats):
        im=ax.imshow(data,vmin=0);fig.colorbar(im,ax=ax,fraction=.047,pad=.04)
        ax.set_title(label,fontsize=10)
        ax.set_xticks(range(3));ax.set_yticks(range(3))
    fig.savefig(output/"yukawa_svd.png",dpi=170)
    plt.close(fig)

    fig, axes=plt.subplots(1,4,figsize=(13,3.6),layout="constrained")
    denom=max(float(np.linalg.norm(total)),1e-300)
    values=[*components,total]
    bound=max(float(np.max(np.abs(z.real/denom))) for z in values)
    for j,(ax,mat) in enumerate(zip(axes,values)):
        im=ax.imshow(mat.real/denom,vmin=-bound,vmax=bound,cmap="coolwarm")
        ax.set_title(f"C12 from F{j+1}" if j<3 else "Total C12",fontsize=11)
        ax.set_xticks(range(3),["e","mu","tau"])
        ax.set_yticks(range(3),["e","mu","tau"])
    fig.colorbar(im,ax=axes.tolist(),fraction=.016,pad=.017,label="Re(C12)/||C12 total||")
    fig.savefig(output/"heavy_c12_terms.png",dpi=170,bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5,4.6),layout="constrained")
    x=np.arange(3); width=0.23
    labels=["theta12_deg","theta13_deg","theta23_deg"]
    for j,(name,record) in enumerate(variation.items()):
        if name=="global_y2":continue
        ax.bar(x+(j-1)*width,[record["angle_differences_deg"][k] for k in labels],width,label=name)
    ax.axhline(0,lw=.8,color="gray")
    ax.set_xticks(x,[r"$\theta_{12}$",r"$\theta_{13}$",r"$\theta_{23}$"])
    ax.set_ylabel("Change in texture-only mixing angle (degrees)")
    ax.set_title("Y2 heavy-column perturbations, no RG evolution")
    ax.legend(title="Varied Y2 column")
    fig.savefig(output/"column_texture_sensitivity.png",dpi=170)
    plt.close(fig)


def _html_report(out: Path, r: dict):
    s=r["svd"]; h=r["heavy_contributions"]; t=r["texture_test"]
    rows="".join(f"<tr><td>F{i+1}</td><td>{q['norm_relative_total']:.5f}</td><td>{q['plus_fraction_texture_shift']:.5f}</td></tr>" for i,q in enumerate(h["per_heavy"]))
    varied="".join(f"<tr><td>{escape(k)}</td>"+"".join(f"<td>{v['angle_differences_deg'][x]:+.4f}</td>" for x in ("theta12_deg","theta13_deg","theta23_deg"))+"</tr>" for k,v in t["variations"].items())
    htm=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>T3-B Yukawa basis</title><style>body{{font:15px system-ui;max-width:1100px;margin:26px auto;padding:0 18px;line-height:1.5}}table{{border-collapse:collapse}}th,td{{text-align:left;border-bottom:1px solid #ccc;padding:8px 12px}}img{{max-width:100%}}code{{overflow-wrap:anywhere}}</style></head><body>
<h1>T3-B-m1: Yukawa basis and matched C12 texture</h1>
<p><strong>Read-only; no UV/EFT running or fitting.</strong> Input couplings are <em>interpreted at the fermion matching threshold</em>, even if their source configuration was historically a UV-start fitted reference. The texture-only angles below are not new low-energy predictions.</p>
<p>Input: <code>{escape(str(r['source']['config']))}</code>; λT3={r['source']['lambdaT3_real']:.7g}, heavy masses (GeV): {', '.join(f'{v:g}' for v in r['source']['masses_gev'])}</p>
<h2>Yukawa SVD</h2><table><tr><td>||Y2|| / ||Y1||</td><td>{s['y2_norm_over_y1_norm']:.6g}</td></tr><tr><td>||offdiag(U1†Y2V1)|| / ||Y2||</td><td><strong>{s['rotated_y2_offdiag_fraction']:.6g}</strong></td></tr><tr><td>Y2 residual to best scalar multiple of Y1</td><td>{s['y2_vs_y1_shape_residual']:.6g}</td></tr></table>
<p>SVD is a mathematical comparison, not a basis where both charged leptons and Majorana heavy masses remain diagonal.</p><img src="yukawa_svd.png" alt="Yukawa structures">
<h2>Heavy-state contributions to C12</h2><table><tr><th>State</th><th>||C12(r)|| / ||C12||</th><th>Texture change after +10% Y2 column</th></tr>{rows}</table>
<p>Sum of component norms / norm of sum: {h['component_norm_sum_over_total']:.6g}. Values above one for individual contributions signify cancellations among heavy states.</p><img src="heavy_c12_terms.png" alt="Matched C12 heavy-flavour contribution matrices">
<h2>Controlled texture perturbation</h2><p>Baseline texture angles (degrees): θ12={t['baseline']['theta12_deg']:.4f}; θ13={t['baseline']['theta13_deg']:.4f}; θ23={t['baseline']['theta23_deg']:.4f}.</p>
<p>All cases are generated by algebraically changing the indicated Yukawa entries while keeping the heavy masses and charged-lepton basis fixed. The <strong>global +10% Y2</strong> direction preserves C12's texture exactly, while <strong>column-specific changes</strong> do not.</p>
<table><tr><th>Y2 variation</th><th>Δθ12 (deg)</th><th>Δθ13 (deg)</th><th>Δθ23 (deg)</th></tr>{varied}</table>
<img src="column_texture_sensitivity.png" alt="Texture-only mixing-angle changes under heavy-column variations">
<p>These are diagnostic leading C12 texture angles, not PMNS observables after scalar matching and Weinberg running. The upstream coupled C12-C5 flow is not used here.</p><p><a href="summary.json">Full numerical summary</a></p></body></html>'''
    (out/"index.html").write_text(htm,encoding="utf-8")


def analyse(config: Path, output: Path, fraction: float):
    data=_load(config)
    if data.get("representation") != {"d_s1":2,"d_s2":2,"d_f":1,"alpha":-1,"shared_scalar":False}:
        # Permit simple excerpts without representation? No: provenance identity matters.
        raise ValueError("Expected ordinary T3-B alpha=-1 (2,2,1,-1)")
    if not np.isfinite(fraction) or fraction <= 0 or fraction >= 1:
        raise ValueError("0 < --fraction < 1 required")
    ordinary=data["base_state"]["ordinary"]
    y1=_matrix(ordinary["y1"],"y1")
    y2=_matrix(ordinary["y2"],"y2")
    mf=_matrix(ordinary["MF"],"MF")
    ye=_matrix(data["base_state"]["sm"]["ye"],"ye")
    y1,y2,masses=_heavy_basis(y1,y2,mf)
    if min(np.linalg.norm(y1),np.linalg.norm(y2))<=0:
        raise ValueError("Both Yukawas must be nonzero")
    u,s,vh=np.linalg.svd(y1)
    rotated=u.conj().T@y2@vh.conj().T
    components,total=_texture(y1,y2,masses)
    if np.linalg.norm(total)==0:
        raise ValueError("Matched total C12 vanishes: no flavour texture")
    base_angles=_texture_angles(total,ye)
    variations={}
    for key, idx in (("heavy_1",0),("heavy_2",1),("heavy_3",2),("global_y2",None)):
        if idx is None:
            yt=y2*(1+fraction)
        else:
            yt=y2.copy();yt[:,idx]*=(1+fraction)
        _,next_t=_texture(y1,yt,masses)
        angles=_texture_angles(next_t,ye)
        variations[key]={"angles_deg":{k:angles[k] for k in ("theta12_deg","theta13_deg","theta23_deg")},
                         "angle_differences_deg":{k:angles[k]-base_angles[k] for k in ("theta12_deg","theta13_deg","theta23_deg")},
                         "shape_residual_to_baseline":_shape(total,next_t)["residual_relative_to_candidate"]}
    per_heavy=[]
    for idx,term in enumerate(components):
        per_heavy.append({"state":idx+1,"mass_gev":float(masses[idx]),
                          "norm":float(np.linalg.norm(term)),
                          "norm_relative_total":float(np.linalg.norm(term)/np.linalg.norm(total)),
                          "plus_fraction_texture_shift":variations[f"heavy_{idx+1}"]["shape_residual_to_baseline"],
                          "c12_term":_json_matrix(term)})
    result={"status":"ReadOnlyT3BYukawaTextureComplete",
            "source":{"config":str(config),"sha256":hashlib.sha256(config.read_bytes()).hexdigest(),
                      "lambdaT3_real":float(ordinary["lambdaT3"]["real"]),"masses_gev":masses.tolist(),
                      "boundary_interpretation":"Input values imposed directly at mu_F; not historical UV-running endpoint"},
            "svd":{"y1_frobenius":float(np.linalg.norm(y1)),"y2_frobenius":float(np.linalg.norm(y2)),
                   "y2_norm_over_y1_norm":float(np.linalg.norm(y2)/np.linalg.norm(y1)),
                   "y1_singular_values":s.tolist(),"y2_singular_values":np.linalg.svd(y2,compute_uv=False).tolist(),
                   "rotated_y2_offdiag_fraction":float(np.linalg.norm(rotated-np.diag(np.diag(rotated)))/np.linalg.norm(rotated)),
                   "y2_vs_y1_shape_residual":_shape(y1,y2)["residual_relative_to_candidate"],
                   "rotated_y2":_json_matrix(rotated)},
            "heavy_contributions":{"c12_at_mu_f_tree":_json_matrix(total),
                                   "component_norm_sum_over_total":float(sum(np.linalg.norm(t) for t in components)/np.linalg.norm(total)),
                                   "per_heavy":per_heavy},
            "texture_test":{"fraction":fraction,
                            "baseline":base_angles,
                            "variations":variations,
                            "notice":"C12 texture-only Takagi directions in charged-lepton mass basis, not the final low-energy C5 PMNS angles"},
            "limitations":["No UV running, EFT1 Wilson RGE, scalar matching or Weinberg evolution performed",
                           "Not the exact locally saved optimised 2.29% result unless --config selects that benchmark",
                           "Y1 SVD basis is not a physical simultaneous mass basis",
                           "Heavy Majorana MF must be diagonal in input basis; positive phase rotations are allowed",
                           "The theta differences are controlled leading-texture responses, not phenomenological predictions"]}
    output.mkdir(parents=True,exist_ok=True)
    (output/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    _figures(output,y1,rotated,components,total,variations)
    _html_report(output,result)
    print(json.dumps({"output":str(output/"index.html"),
       "y2_over_y1":result["svd"]["y2_norm_over_y1_norm"],
       "y2_offdiag_fraction_in_y1_svd_basis":result["svd"]["rotated_y2_offdiag_fraction"],
       "c12_heavy_cancellation_ratio":result["heavy_contributions"]["component_norm_sum_over_total"],
       "angle_changes_deg":{k:v["angle_differences_deg"] for k,v in variations.items()}},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config",type=Path,default=DEFAULT_CONFIG)
    p.add_argument("--output",type=Path,default=DEFAULT_OUT)
    p.add_argument("--fraction",type=float,default=0.1,help="Fractional Y2 change per heavy column (default 0.1)")
    a=p.parse_args()
    analyse(a.config,a.output,a.fraction)

if __name__=="__main__":main()
