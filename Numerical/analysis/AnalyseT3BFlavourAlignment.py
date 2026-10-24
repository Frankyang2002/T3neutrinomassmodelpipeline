"""Read-only alignment diagnostics for saved T3-B-m1 C12/C5 matrices.

Requires existing exploratory coupled and saved low-energy comparison JSON.
Never runs UV/EFT RGEs, optimisation, threshold matching or symbolic code.

The mass-basis portion is restricted to real-symmetric C5 benchmarks. It checks
that their eigensystem reproduces the saved |PMNS| before interpreting it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from html import escape
from pathlib import Path

import numpy as np

DEFAULT_COUPLED = Path("output/t3b_fitted_common_scale_threshold.json")
DEFAULT_LOW = Path("output/t3b_fitted_low_energy_comparison.json")
DEFAULT_CONTROL = Path("output/t3b_fitted_fixed_order_control.json")
DEFAULT_OUT = Path("output/studies/t3b_flavour_alignment/fitted_reference_saved")


def load(path: Path) -> dict:
    with path.open("r", encoding="utf-8-sig") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"JSON object required: {path}")
    return data


def matrix(item, name: str) -> np.ndarray:
    if not isinstance(item, dict) or not {"real", "imag"} <= item.keys():
        raise ValueError(f"{name}: expected real/imag object")
    m = np.asarray(item["real"], dtype=float) + 1j * np.asarray(item["imag"], dtype=float)
    if m.shape != (3, 3) or not np.all(np.isfinite(m)):
        raise ValueError(f"{name}: finite complex 3x3 expected")
    if np.linalg.norm(m - m.T) > 1e-9 * max(np.linalg.norm(m), 1e-300):
        raise ValueError(f"{name}: Majorana matrix not symmetric")
    return m


def relative(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(b))
    if denom == 0:
        raise ValueError("Zero-norm reference matrix")
    return float(np.linalg.norm(a-b) / denom)


def projection(reference: np.ndarray, candidate: np.ndarray) -> dict:
    """Complex Frobenius projection of candidate onto reference."""
    d = np.vdot(reference, reference)
    if not np.isfinite(d) or abs(d) == 0 or np.linalg.norm(candidate) == 0:
        raise ValueError("Cannot project a zero matrix")
    alpha = np.vdot(reference, candidate) / d
    perp = candidate-alpha*reference
    return {
        "alpha_real": float(alpha.real), "alpha_imag": float(alpha.imag),
        "shape_residual_fraction_of_candidate": float(np.linalg.norm(perp)/np.linalg.norm(candidate)),
        "orthogonal_fraction_of_reference": float(np.linalg.norm(perp)/np.linalg.norm(reference)),
        "frobenius_cosine": float(abs(np.vdot(reference,candidate)) / (np.linalg.norm(reference)*np.linalg.norm(candidate))),
    }


def safe_stored_relative(stored: object, computed: float, key: str):
    if stored is None:
        return
    if not np.isclose(float(stored), computed, rtol=1e-6, atol=1e-10):
        raise ValueError(f"Saved {key}={stored} differs from reconstructed {computed}")


def real_mass_basis(fixed: np.ndarray, corrected: np.ndarray, pmns_abs) -> tuple[dict, np.ndarray]:
    """Eigen/Takagi basis for verified real-symmetric CP-conserving matrices.

    Returns signed-eigenvalue perturbative rotations; do not use this shortcut
    for complex CP-violating inputs or a nontrivial charged-lepton basis.
    """
    if max(np.linalg.norm(fixed.imag),np.linalg.norm(corrected.imag)) > 1e-10 * np.linalg.norm(fixed):
        raise ValueError("Mass-basis shortcut only supports effectively real C5; use repository Takagi/charged-lepton routine")
    eig, U = np.linalg.eigh(fixed.real)
    order = np.argsort(np.abs(eig))
    eig, U = eig[order], U[:,order]
    pmns = np.asarray(pmns_abs, float)
    if pmns.shape != (3,3):
        raise ValueError("Saved PMNS absolute matrix missing or malformed")
    pmns_error = float(np.max(np.abs(np.abs(U)-pmns)))
    if pmns_error > 1e-6:
        raise ValueError(f"Mass-basis vectors differ from saved |PMNS|: {pmns_error:.4g}. Charged lepton rotation/ordering may matter.")
    delta = (corrected-fixed).real
    rotated = U.T@delta@U
    offdiag = rotated-np.diag(np.diag(rotated))
    denom = float(np.linalg.norm(fixed))
    gaps = {}
    for i in range(3):
        for j in range(i+1,3):
            gap = eig[j]-eig[i]
            if abs(gap) <= 1e-10*max(np.linalg.norm(eig),1e-300):
                gaps[f"{i+1}{j+1}"] = {"status":"near_degenerate"}
            else:
                theta = float(rotated[i,j]/gap)
                gaps[f"{i+1}{j+1}"] = {"signed_eigenvalue_gap": float(gap),
                                       "first_order_rotation_radians": theta,
                                       "first_order_rotation_degrees": float(np.degrees(theta))}
    metrics = {
        "pmns_abs_reconstruction_max_error":pmns_error,
        "signed_fixed_eigenvalues":eig.tolist(),
        "offdiagonal_norm_relative_to_fixed":float(np.linalg.norm(offdiag)/denom),
        "offdiagonal_norm_relative_to_delta":float(np.linalg.norm(offdiag)/np.linalg.norm(delta)),
        "mass_basis_delta_over_fixed_norm":(rotated/denom).tolist(),
        "first_order_signed_eigenvector_rotations":gaps,
        "warning":"Perturbative signed-eigenvalue rotations are diagnostic only; saved exact PMNS differences take precedence.",
    }
    return metrics, rotated


def plot_heatmaps(C12, fixed, delta, rotated, fixed_low, out: Path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    scale=float(np.linalg.norm(fixed))
    alpha=np.vdot(C12,fixed)/np.vdot(C12,C12)
    matrices=[('Aligned C12 / ||C5||', (alpha*C12/scale).real),
              ('Fixed C5 / ||C5||', (fixed/scale).real),
              ('Delta C5 / ||C5||', (delta/scale).real)]
    fig,axs=plt.subplots(1,3,figsize=(10,3.6),layout='constrained')
    for ax,(heading,arr) in zip(axs,matrices):
        im=ax.imshow(arr)
        ax.set_title(heading,fontsize=10)
        ax.set_xticks(range(3),['e','mu','tau']);ax.set_yticks(range(3),['e','mu','tau'])
        fig.colorbar(im,ax=ax,fraction=0.045,pad=0.04)
    fig.suptitle('Scalar threshold: matrix shapes (same flavour basis)',fontsize=12)
    fig.savefig(out/'c12_c5_alignment.png',dpi=160)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.5,4.4),layout='constrained')
    norm_rot=rotated/float(np.linalg.norm(fixed_low))
    im=ax.imshow(norm_rot)
    fig.colorbar(im,ax=ax,label='Delta C5(low) / ||fixed C5(low)||')
    for i in range(3):
        for j in range(3):ax.text(j,i,f'{norm_rot[i,j]:+.2e}',ha='center',va='center',fontsize=9,color='black')
    ax.set_xticks(range(3),['1','2','3']);ax.set_yticks(range(3),['1','2','3'])
    ax.set_title('Low-energy correction in fixed-order mass basis')
    fig.savefig(out/'mass_basis_delta.png',dpi=170)
    plt.close(fig)


def write_report(payload: dict, out: Path):
    a=payload['ms']['delta_vs_fixed']
    low=payload['low_energy']
    mb=low['mass_basis']
    angle_rows=''.join('<tr><td>'+escape(k)+'</td><td>'+f'{v:.6g}'+'</td></tr>' for k,v in low['saved_angle_differences_deg'].items())
    html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>T3-B-m1 flavour alignment</title>
<style>body{{font:15px system-ui;margin:28px auto;max-width:1020px;padding:0 18px;line-height:1.5}}img{{max-width:100%}}table{{border-collapse:collapse}}td,th{{padding:6px 14px;border-bottom:1px solid #ccc;text-align:left}}code{{overflow-wrap:anywhere}}.subtle{{color:#555}}</style></head><body>
<h1>T3-B-m1 saved-matrix flavour alignment</h1>
<p>Read-only analysis. <strong>No UV evolution, EFT running, fitting or threshold calculation</strong> was performed by this script. Saved upstream coupled results remain exploratory, not a validated physical prediction.</p>
<p>Fitted-reference source: {escape(str(payload['source_inputs'].get('source_config','not independently verified')))}. This historical dataset was generated using its existing running prescription; this diagnostic does not convert it into a UV-frozen study.</p>
<h2>1. Scalar-threshold shape alignment</h2><table>
<tr><th>Quantity</th><th>Value</th></tr>
<tr><td>||Delta C5|| / ||C5 fixed||</td><td>{payload['ms']['delta_relative_to_fixed']:.8g}</td></tr>
<tr><td>Shape residual C12 to fixed C5</td><td>{payload['ms']['c12_vs_fixed']['shape_residual_fraction_of_candidate']:.8g}</td></tr>
<tr><td>Shape residual C12 to direct C5</td><td>{payload['ms']['c12_vs_direct']['shape_residual_fraction_of_candidate']:.8g}</td></tr>
<tr><td>Parallel correction scale Re(alpha)</td><td>{a['alpha_real']:.8g}</td></tr>
<tr><td>||perpendicular Delta C5||/||fixed C5||</td><td>{a['orthogonal_fraction_of_reference']:.8g}</td></tr></table>
<img src="c12_c5_alignment.png" alt="C12, C5 and Delta C5 shape matrices">
<h2>2. Low-energy mass eigenbasis</h2><p>PMNS eigenvector reproduction error: {mb['pmns_abs_reconstruction_max_error']:.3g}. Off-diagonal correction relative to fixed C5: {mb['offdiagonal_norm_relative_to_fixed']:.8g}.</p>
<img src="mass_basis_delta.png" alt="Correction rotated into mass eigenbasis">
<h3>Saved exact mixing-angle changes</h3><table><tr><th>Angle</th><th>Difference (degrees)</th></tr>{angle_rows}</table>
<p class="subtle">The sign and phase conventions for a Majorana Takagi decomposition must be respected. The mass-basis shortcut is permitted only when the saved real symmetric matrix reproduces |PMNS|. Full RGE/operator closure is not established by this diagnostic.</p>
<p><a href="summary.json">Machine-readable numerical results</a></p></body></html>'''
    (out/'index.html').write_text(html,encoding='utf-8')


def analyse(coupled_path: Path, low_path: Path, control_path: Path, output: Path):
    coupled, low, control=load(coupled_path),load(low_path),load(control_path)
    if coupled.get('status')!='ExploratoryCoupledMatrixComplete':
        raise ValueError('Unexpected coupled status')
    if coupled.get('threshold_scheme')!='common_scale_log_boundary':
        raise ValueError('Require common_scale_log_boundary, not legacy extrapolation')
    if control.get('status')!='FixedOrderControlComplete' or low.get('status')!='ExploratoryLowEnergyComparisonComplete':
        raise ValueError('Unexpected saved baseline/low-energy status')
    c12=matrix(coupled['c12_at_ms'],'C12(MS)')
    hard=matrix(coupled['production_hard_c5_at_ms'],'hard C5')
    direct=matrix(coupled['production_direct_c5_at_ms'],'direct C5')
    fixed=matrix(control['fixed_order_c5_at_ms'],'fixed C5')
    corrected=matrix(coupled['exploratory_total_c5_at_ms'],'coupled C5')
    if relative(hard+direct,fixed)>1e-8:
        raise ValueError('Hard+direct does not reproduce production fixed C5')
    if relative(hard+matrix(coupled['coupled_eft1_c5_at_ms'],'coupled EFT1'),corrected)>1e-8:
        raise ValueError('Hard+coupled EFT1 does not reproduce corrected C5')
    delta=corrected-fixed
    reported=relative(corrected,fixed)
    safe_stored_relative(coupled.get('relative_difference_from_fixed_order'),reported,'relative_difference_from_fixed_order')
    safe_stored_relative(low.get('relative_c5_difference_at_ms'),reported,'low MS relative C5')
    fixed_low=matrix(low['fixed_order']['c5_at_low_scale'],'fixed low C5')
    coupled_low=matrix(low['exploratory_coupled']['c5_at_low_scale'],'coupled low C5')
    low_delta=coupled_low-fixed_low
    low_relative=relative(coupled_low,fixed_low)
    safe_stored_relative(low.get('relative_c5_difference_at_low_scale'),low_relative,'low relative C5')
    mass_basis,rotated=real_mass_basis(fixed_low,coupled_low,low['fixed_order']['observables']['pmns_abs'])
    angles=low.get('differences',{}).get('mixing_angle_differences_deg',{})
    if set(angles)!={'theta12_deg','theta13_deg','theta23_deg'}:
        raise ValueError('Saved low-energy angle difference set missing')
    result={
        'status':'ReadOnlyAlignmentDiagnosticComplete',
        'source_inputs':{
            'source_config':low.get('source_inputs',{}).get('config'),
            'coupled':str(coupled_path),'low_energy':str(low_path),'control':str(control_path),
            'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (coupled_path,low_path,control_path)},
        },
        'scope':'Saved matrices only. No new dynamics or validated higher-order predictions.',
        'ms':{'mu_gev':float(control['mu_s_gev']),
              'delta_relative_to_fixed':reported,
              'c12_vs_fixed':projection(c12,fixed),
              'c12_vs_direct':projection(c12,direct),
              'c12_vs_hard':projection(c12,hard),
              'delta_vs_fixed':projection(fixed,delta)},
        'low_energy':{'mu_gev':float(low['mu_low_gev']),
                      'delta_relative_to_fixed':low_relative,
                      'delta_vs_fixed':projection(fixed_low,low_delta),
                      'mass_basis':mass_basis,
                      'saved_angle_differences_deg':angles},
        'limitations':['The historical saved reference used its original running prescription; no UV freeze is retroactively imposed',
                       'Exploratory 2-operator C12/C5 closure and higher-order matching are unverified',
                       'Real-symmetric CP-conserving eigenbasis shortcut verified against saved |PMNS|'],
    }
    output.mkdir(parents=True,exist_ok=True)
    (output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    plot_heatmaps(c12,fixed,delta,rotated,fixed_low,output)
    write_report(result,output)
    print(json.dumps({'output':str(output/'index.html'),
                      'delta_relative_to_fixed':result['ms']['delta_relative_to_fixed'],
                      'shape_residual_c12_fixed':result['ms']['c12_vs_fixed']['shape_residual_fraction_of_candidate'],
                      'perpendicular_delta_relative_fixed':result['ms']['delta_vs_fixed']['orthogonal_fraction_of_reference'],
                      'low_offdiagonal_relative_fixed':mass_basis['offdiagonal_norm_relative_to_fixed'],
                      'pmns_reconstruction_error':mass_basis['pmns_abs_reconstruction_max_error']},indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coupled',type=Path,default=DEFAULT_COUPLED)
    p.add_argument('--low-energy',type=Path,default=DEFAULT_LOW)
    p.add_argument('--control',type=Path,default=DEFAULT_CONTROL)
    p.add_argument('--output',type=Path,default=DEFAULT_OUT)
    args=p.parse_args()
    analyse(args.coupled,args.low_energy,args.control,args.output)

if __name__=='__main__':
    main()
