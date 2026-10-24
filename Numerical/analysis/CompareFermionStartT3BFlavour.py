"""Read-only T3-B-m1 flavour diagnostic for already-fitted UV-FROZEN inputs.

Reconstruct the existing 011 fixed-order pipeline from mu_F, then optionally
run the repository's opt-in EXPLORATORY C12<->C5 two-operator system. Neither
fit parameters nor production code/data are changed. No UV beta evaluation.

Run: python -m Numerical.analysis.CompareFermionStartT3BFlavour
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import traceback

import numpy as np

MODEL = 'T3_dS1_2_dS2_2_dF_1_alpha_m1'
CASES = ('5e-9', '0p5')
LAMBDAS = {'5e-9': 5e-9, '0p5': 0.5}
OUT = Path('output/studies/t3b_uv_frozen_flavour_alignment')
STUDY = Path('output/studies/fermion_start_lambda_comparison')


def read(path):
    with Path(path).open(encoding='utf-8-sig') as f:
        result = json.load(f)
    if not isinstance(result, dict):
        raise ValueError('JSON object expected: '+str(path))
    return result


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def mat(obj):
    m = np.asarray(obj['real'], float) + 1j*np.asarray(obj['imag'], float)
    if m.shape != (3, 3) or not np.isfinite(m).all():
        raise ValueError('Expected finite complex 3x3 array')
    return m


def serial(m):
    x = np.asarray(m, complex)
    return {'real': x.real.tolist(), 'imag': x.imag.tolist(),
            'frobenius_norm': float(np.linalg.norm(x))}


def proj(reference, candidate):
    """Projection onto the complex line of a nonzero Frobenius matrix."""
    nr = float(np.linalg.norm(reference)); nc = float(np.linalg.norm(candidate))
    if min(nr, nc) <= 0:
        raise ValueError('Matrix projection requires two nonzero matrices')
    a = np.vdot(reference, candidate)/np.vdot(reference, reference)
    rem = candidate-a*reference
    return {'alpha_real': float(a.real), 'alpha_imag': float(a.imag),
            'orthogonal_over_reference': float(np.linalg.norm(rem)/nr),
            'orthogonal_over_candidate': float(np.linalg.norm(rem)/nc),
            'frobenius_cosine': float(abs(np.vdot(reference, candidate))/(nr*nc))}


def angles(pmns):
    a = np.asarray(pmns, float)
    s13 = float(a[0, 2]**2); c13 = 1-s13
    if c13 <= 0:
        raise ValueError('Invalid PMNS angle convention')
    vals = [float(a[0,1]**2/c13), s13, float(a[1,2]**2/c13)]
    return dict(zip(('theta12_deg','theta13_deg','theta23_deg'),
                    [float(np.degrees(np.arcsin(np.sqrt(np.clip(v,0,1))))) for v in vals]))


def neutrino_summary(K, ye, vev):
    from physics.NeutrinoTrajectory import (neutrino_mass_from_c5,
                                            charged_lepton_mass_basis_matrix)
    from physics.NeutrinoObservables import calculate_neutrino_observables, takagi_factorization
    mass = charged_lepton_mass_basis_matrix(neutrino_mass_from_c5(K, vev_gev=vev), ye)
    obs = calculate_neutrino_observables(mass, ordering='NO')
    if obs.takagi_residual > 1e-7:
        raise RuntimeError('Takagi numerical residual too large: '+str(obs.takagi_residual))
    U, masses, residual = takagi_factorization(mass)
    if residual > 1e-7 or np.max(np.abs(abs(U)-obs.pmns_abs)) > 1e-7:
        raise RuntimeError('Takagi/PMNS conventions inconsistent')
    return ({'masses_ev': np.asarray(obs.masses_ev).tolist(),
             'theta_deg': angles(obs.pmns_abs),
             'pmns_abs': np.asarray(obs.pmns_abs).tolist(),
             'takagi_residual': float(obs.takagi_residual),
             'delta_m21_sq_ev2':float(obs.delta_m21_sq_ev2),
             'delta_m31_sq_ev2':float(obs.delta_m31_sq_ev2)}, mass, U)


def yukawa_diagnostic(y1, y2, masses, ye, fraction=0.1):
    """Yukawa SVD + C12 heavy structure; angles here are TEXTURE ONLY."""
    from Numerical.diagnostics.T3BHeavyLogAlignment import flavor_kernel_by_heavy
    from physics.NeutrinoTrajectory import charged_lepton_mass_basis_matrix
    from physics.NeutrinoObservables import calculate_neutrino_observables
    left, sv, vh = np.linalg.svd(y1)
    transformed_y2 = left.conj().T @ y2 @ vh.conj().T
    off = transformed_y2 - np.diag(np.diag(transformed_y2))
    terms = flavor_kernel_by_heavy(y1, y2, masses)
    total = terms.sum(axis=0)
    if np.linalg.norm(total) <= 0:
        raise ValueError('Vanishing C12 texture')

    def texture_angles(C):
        # Only eigenvector angles of the matched flavour texture, not physical C5.
        basis = charged_lepton_mass_basis_matrix(C, ye)
        return angles(calculate_neutrino_observables(basis, ordering='NO').pmns_abs)

    initial = texture_angles(total)
    trials = {}
    for label, idx in [('global_y2', None), ('heavy_1',0), ('heavy_2',1), ('heavy_3',2)]:
        varied = np.array(y2,copy=True)
        if idx is None:
            varied *= (1+fraction)
        else:
            varied[:,idx] *= (1+fraction)
        candidate = flavor_kernel_by_heavy(y1,varied,masses).sum(axis=0)
        trial = texture_angles(candidate)
        trials[label] = {'angle_differences_deg':{k: trial[k]-initial[k] for k in initial},
                         'shape_difference_over_candidate':proj(total,candidate)['orthogonal_over_candidate']}
    return ({'y1_frobenius':float(np.linalg.norm(y1)),
             'y2_frobenius':float(np.linalg.norm(y2)),
             'y2_over_y1':float(np.linalg.norm(y2)/np.linalg.norm(y1)),
             'y1_singular_values':sv.tolist(),
             'y2_offdiagonal_fraction_in_y1_svd_basis':float(np.linalg.norm(off)/np.linalg.norm(y2)),
             'transformed_y2':serial(transformed_y2),
             'heavy_masses_gev':np.asarray(masses).tolist(),
             'heavy_component_norm_over_total': [float(np.linalg.norm(x)/np.linalg.norm(total)) for x in terms],
             'sum_component_norms_over_total':float(sum(np.linalg.norm(x) for x in terms)/np.linalg.norm(total)),
             'texture_only_baseline_angles_deg':initial,
             'ten_percent_y2_variations':trials,
             'c12_tree_at_mu_f':serial(total)}, terms, total)


def _plot(path, C12, fixed, delta, rotated, mass_fixed_norm, y1, y2, terms):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def heatmap(name, arrays, titles, normalization=1.0, axis_labels=('e','mu','tau')):
        fig, axes = plt.subplots(1, len(arrays), figsize=(4.1*len(arrays),4),
                                 layout='constrained',squeeze=False)
        for ax,x,title in zip(axes[0], arrays, titles):
            z=np.asarray(x)/normalization
            im=ax.imshow(np.abs(z))
            ax.set_title(title, fontsize=10)
            ax.set_xticks(range(3),axis_labels)
            ax.set_yticks(range(3),axis_labels)
            fig.colorbar(im,ax=ax,fraction=.05,pad=.04)
        fig.savefig(path/name,dpi=155);plt.close(fig)
    d=float(np.linalg.norm(fixed))
    a=np.vdot(C12,fixed)/np.vdot(C12,C12)
    heatmap('c12_vs_c5.png',[a*C12,fixed,delta],['Projected C12','Fixed C5','Coupled minus fixed'],d)
    heatmap('mass_basis_delta.png',[rotated],['|U^T delta(mnu) U| / ||mnu_fixed||'],mass_fixed_norm,('1','2','3'))
    heatmap('yukawa_matrices.png',[y1,y2],['|Y1|','|Y2|'])
    heatmap('c12_heavy_terms.png',[*terms,terms.sum(axis=0)],['Heavy 1','Heavy 2','Heavy 3','Sum'],float(np.linalg.norm(terms.sum(axis=0))))


def _html_case(case, result):
    yc=result['yukawa'];a=result['alignment'];obs=result['observables'];m=result['mass_basis']
    tables=[]
    for title,metrics in [('C12 vs C5',a['c12_vs_fixed']),('Correction vs C5',a['delta_vs_fixed'])]:
        tables.append(f'<tr><td>{html.escape(title)}</td><td>{metrics["orthogonal_over_candidate"]:.6g}</td><td>{metrics["orthogonal_over_reference"]:.6g}</td></tr>')
    ats=''.join(f'<tr><td>{html.escape(k)}</td><td>{obs["fixed"]["theta_deg"][k]:.6f}</td><td>{obs["coupled"]["theta_deg"][k]:.6f}</td><td>{v:+.6f}</td></tr>' for k,v in obs['angle_changes_deg'].items())
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>T3B {html.escape(case)} UV-frozen flavour</title><style>body{{font:15px system-ui;max-width:1120px;margin:28px auto;padding:0 20px;line-height:1.5}}img{{max-width:100%}}td,th{{padding:7px 14px;border-bottom:1px solid #aaa;text-align:left}}table{{border-collapse:collapse}}code{{overflow-wrap:anywhere}}</style><body>
<h1>T3-B-m1: UV-frozen fitted {html.escape(case)} flavour diagnostics</h1><p><strong>Numerical scope:</strong> fixed-order fermion matching + true EFT1 running + hard/direct C5 matching + SM/Weinberg running. UV beta evolution is zero. An independent <strong>exploratory</strong> coupled C12-C5 LLSS/LLHH projection is compared without overwriting production C5. Full operator closure and complete higher-order matching are not established.</p>
<p>Input config: <code>{html.escape(result['input']['config'])}</code> (SHA-256 {result['input']['sha256']}).</p>
<h2>Yukawa structure at fermion matching</h2><p>||Y2||/||Y1|| = {yc['y2_over_y1']:.6g}; Y2 off-diagonal fraction in Y1 SVD basis = {yc['y2_offdiagonal_fraction_in_y1_svd_basis']:.6g}; heavy C12 cancellation ratio = {yc['sum_component_norms_over_total']:.6g}.</p>
<img src="yukawa_matrices.png" alt="Yukawa heatmaps"><img src="c12_heavy_terms.png" alt="Heavy-state C12 decomposition">
<h2>Wilson flavour alignment at scalar threshold</h2><p>||delta C5||/||fixed C5|| = {a['delta_relative_to_fixed']:.8g}. Projected complex rescaling coefficient alpha = {a['delta_vs_fixed']['alpha_real']:.8g} + i({a['delta_vs_fixed']['alpha_imag']:.8g}).</p>
<table><tr><th>Alignment</th><th>Orthogonal / candidate</th><th>Orthogonal / reference</th></tr>{''.join(tables)}</table><img src="c12_vs_c5.png" alt="Wilson alignment">
<h2>Low-energy mass-eigenbasis</h2><p>||offdiag(U^T delta M U)||/||M fixed|| = {m['offdiag_over_fixed']:.8g}; offdiagonal / ||delta M|| = {m['offdiag_over_delta']:.8g}.</p><img src="mass_basis_delta.png" alt="Mass eigenbasis correction">
<table><tr><th>Angle</th><th>Fixed</th><th>Coupled (exploratory)</th><th>Change (degrees)</th></tr>{ats}</table>
<h2>Interpretation</h2><p>Yukawa +10% column perturbations in <a href="summary.json">summary.json</a> are matched C12 texture-only; they are NOT new physical PMNS predictions. Small mass-basis offdiagonal entries are consistent with stable angles, but no new validated NLO prediction is claimed.</p></body></html>'''


def _run_one(case, config, fit_path, artifact_map, out, n_samples):
    from types import SimpleNamespace
    from Numerical.orchestration import RunFermionStartLambdaStudy as study
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.fitting.ScanCLI import build_uv_state_from_config
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
    from Numerical.running.T3Trajectory import continue_with_weinberg_running
    from Numerical.diagnostics.FinalC5ContributionDiagnostics import evaluate_final_c5_contributions
    from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs, evaluate_final_c5_on_trajectory
    from Numerical.orchestration.T3BIntegratedCoupledPipeline import _intermediate_samples
    from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import (
        calculate as coupled_calculate, fixed_order_log_from_common_boundary)
    from Numerical.diagnostics.T3BHeavyLogAlignment import matrix_from_json

    cfg=read(config); fit=read(fit_path)
    if cfg.get('representation') != dict(d_s1=2,d_s2=2,d_f=1,alpha=-1,shared_scalar=False):
        raise ValueError('Input not ordinary T3-B-m1')
    if fit.get('status')!='Success' or float(fit.get('chi2_fit',float('inf'))) >= 5.86:
        raise ValueError('Not a successfully optimised fit with chi2 < 5.86')
    marker=cfg.get('fermion_start_study',{})
    if marker.get('case')!=case or 'UV RGEs disabled' not in marker.get('boundary_definition',''):
        raise ValueError('Not a verified fermion-start frozen-UV fitted configuration')
    lam=complex(cfg['base_state']['ordinary']['lambdaT3']['real'],cfg['base_state']['ordinary']['lambdaT3'].get('imag',0))
    if abs(lam.imag)>1e-13 or not np.isclose(lam.real,LAMBDAS[case],rtol=1e-12,atol=0):
        raise ValueError(f'Wrong scalar coupling for {case}: {lam}')
    scales=cfg['scales']; mu_f=float(scales['mu_fermion_threshold_gev']); mu_s=float(scales['mu_scalar_threshold_gev']);mu_low=float(scales['mu_low_gev'])
    mu_uv=float(cfg['base_state']['mu_gev']);opt=cfg.get('numerical',{})
    rtol=float(opt.get('rtol',1e-8));atol=float(opt.get('atol',1e-11));vev=float(opt.get('vev_gev',246.22))
    if cfg.get('ordering','NO').upper()!='NO':
        raise ValueError('Current mass-basis comparison requires normal ordering')
    amap=read(artifact_map) if artifact_map else None
    artifacts=base.locate_artifacts(MODEL,amap)
    initial=build_uv_state_from_config(SimpleNamespace(raw={**cfg,'scan':cfg.get('scan',{'bindings':{}})}),{}).validated()
    uv=study.frozen_result(initial,mu_f,uv=True)
    eft_initial=build_intermediate_scalar_boundary(uv.final_state,matching_scale_gev=mu_f)
    beta=study.load_rgbeta_payload(artifacts[1])
    eft=run_intermediate_scalar_segment(eft_initial,beta,mu_s,
                                        save_scales_gev=np.geomspace(mu_f,mu_s,n_samples),
                                        rtol=rtol,atol=atol)
    traj=study.with_uv_eft_states(uv,eft,mu_uv,mu_f,mu_s)
    assert uv.nfev == 0, 'UV RGEs unexpectedly active'
    breakdown=evaluate_final_c5_contributions(artifacts[2],traj)
    fixed=breakdown.combined
    independent=evaluate_final_c5_on_trajectory(artifacts[2],traj)
    fixed_check=float(np.linalg.norm(fixed-independent)/max(np.linalg.norm(fixed),1e-300))
    if fixed_check>2e-6:
        raise ValueError(f'Fixed C5 breakdown and official evaluator disagree: {fixed_check:.4g}')
    heavy=extract_final_c5_inputs(traj)
    samples=_intermediate_samples(traj)
    predicted=fixed_order_log_from_common_boundary(samples,heavy.y1,heavy.y2,heavy.heavy_masses_gev)
    direct_check=float(np.linalg.norm(predicted-breakdown.direct_running)/max(np.linalg.norm(breakdown.direct_running),1e-300))
    if direct_check>1e-8:
        raise ValueError(f'Unfitted leading-log direct C5 normalization mismatch: {direct_check:.4g}')
    control={'status':'FixedOrderControlComplete','mu_f_gev':mu_f,'mu_s_gev':mu_s,
             'intermediate_running':{'samples':samples},'fixed_order_c5_at_ms':serial(fixed)}
    norm={'status':'UnfittedNormalizationPassed','production_direct_c5':serial(breakdown.direct_running),
          'predicted_direct_c5':serial(predicted)}
    exploratory=coupled_calculate(control,norm,heavy.y1,heavy.y2,heavy.heavy_masses_gev,
                                   diagonal=True,threshold_scheme='common_scale_log_boundary')
    candidate=matrix_from_json(exploratory['exploratory_total_c5_at_ms'])
    c12s=matrix_from_json(exploratory['c12_at_ms'])
    delta=candidate-fixed
    ysummary,terms,c12f=yukawa_diagnostic(heavy.y1,heavy.y2,heavy.heavy_masses_gev,traj.eft1_threshold_state.sm.ye)
    fixed_low=continue_with_weinberg_running(traj,fixed,mu_low,rtol=rtol,atol=atol).weinberg
    coupled_low=continue_with_weinberg_running(traj,candidate,mu_low,rtol=rtol,atol=atol).weinberg
    low_fixed,mass_fixed,U=neutrino_summary(fixed_low.K,fixed_low.ye,vev)
    low_coupled,mass_coupled,_=neutrino_summary(coupled_low.K,coupled_low.ye,vev)
    if np.linalg.norm(fixed_low.ye-coupled_low.ye)>1e-8*max(np.linalg.norm(fixed_low.ye),1e-30):
        raise ValueError('Fixed and coupled Weinberg paths have different charged-lepton backgrounds')
    dmass=mass_coupled-mass_fixed
    dmass_basis=U.T@dmass@U
    off=dmass_basis-np.diag(np.diag(dmass_basis))
    dnorm=float(np.linalg.norm(dmass))
    f_norm=float(np.linalg.norm(mass_fixed))
    result={'status':'UVFrozenT3BExploratoryFlavourComparisonComplete',
            'scope':'Actual fixed-order 011, separate opt-in exploratory two-operator C12-C5 LL comparison',
            'input':{'config':str(config),'sha256':hashlib.sha256(Path(config).read_bytes()).hexdigest(),
                     'fit':str(fit_path),'fit_chi2':float(fit['chi2_fit']),
                     'symbolic_artifacts':[str(x) for x in artifacts],
                     'case':case,'lambda_t3':lam.real,'mu_f_gev':mu_f,'mu_s_gev':mu_s,
                     'mu_low_gev':mu_low,'uv_running':False,'uv_nfev':int(uv.nfev),
                     'eft1_nfev':int(eft.nfev),'eft1_samples':int(eft.n_points)},
            'normalization_checks':{'fixed_c5_evaluator_relative':fixed_check,
                                    'direct_c5_leading_log_relative':direct_check},
            'yukawa':ysummary,
            'alignment':{'delta_relative_to_fixed':float(np.linalg.norm(delta)/np.linalg.norm(fixed)),
                         'c12_at_ms':serial(c12s),
                         'fixed_c5_at_ms':serial(fixed),
                         'coupled_exploratory_c5_at_ms':serial(candidate),
                         'delta_c5_at_ms':serial(delta),
                         'c12_vs_fixed':proj(c12s,fixed),
                         'c12_vs_direct':proj(c12s,breakdown.direct_running),
                         'delta_vs_fixed':proj(fixed,delta),
                         'c12_initial_vs_final':proj(c12f,c12s),
                         'reverse_feedback_fraction':float(exploratory['relative_reverse_feedback'])},
            'mass_basis':{'offdiag_over_fixed':float(np.linalg.norm(off)/f_norm),
                          'offdiag_over_delta':float(np.linalg.norm(off)/dnorm) if dnorm>0 else None,
                          'delta_mass_in_fixed_takagi_basis':serial(dmass_basis),
                          'charged_lepton_mass_basis':'repo physics.NeutrinoTrajectory.charged_lepton_mass_basis_matrix'},
            'observables':{'fixed':low_fixed,'coupled':low_coupled,
                           'angle_changes_deg':{k:low_coupled['theta_deg'][k]-low_fixed['theta_deg'][k] for k in low_fixed['theta_deg']}},
            'warnings':['Coupled C12-C5 LLSS/LLHH RGE is exploratory; complete operator closure, consistent NLO matching and sequential heavy thresholds are not verified.',
                        'Texture-only Yukawa perturbations are not fitted neutrino predictions.',
                        'C5 scalar matching is the existing production fixed-order result; high lambda perturbativity not established by this diagnostic.']}
    out.mkdir(parents=True,exist_ok=True)
    save(out/'summary.json',result)
    (out/'index.html').write_text(_html_case(case,result),encoding='utf-8')
    _plot(out,c12s,fixed,delta,dmass_basis,f_norm,heavy.y1,heavy.y2,terms)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',default=','.join(CASES),help='Comma-separated cases: 5e-9,0p5')
    parser.add_argument('--artifact-map',type=Path)
    parser.add_argument('--output',type=Path,default=OUT)
    parser.add_argument('--samples',type=int,default=121)
    args=parser.parse_args()
    names=args.cases.split(',')
    if not names or len(set(names))!=len(names) or not set(names)<=set(CASES):
        parser.error('Cases must be unique members of 5e-9,0p5')
    if args.samples<3:parser.error('At least 3 EFT1 background samples required')
    args.output.mkdir(parents=True,exist_ok=True)
    items=[]
    for case in names:
        cfg=STUDY/case/'fitted'/f'{MODEL}.json'
        fit=STUDY/case/'fits'/f'{MODEL}.json'
        if not cfg.is_file() or not fit.is_file():
            row={'case':case,'status':'MissingFittedInput','config':str(cfg),'fit':str(fit)}
            print(f'{case}: missing successfully fitted input; no surrogate generated')
        else:
            try:
                r=_run_one(case,cfg,fit,args.artifact_map,args.output/case,args.samples)
                a=r['alignment'];m=r['mass_basis'];o=r['observables'];
                row={'case':case,'status':'Success','lambdaT3':r['input']['lambda_t3'],
                     'chi2':r['input']['fit_chi2'],
                     'y2_over_y1':r['yukawa']['y2_over_y1'],
                     'y2_offdiag_fraction':r['yukawa']['y2_offdiagonal_fraction_in_y1_svd_basis'],
                     'c12_heavy_cancellation_ratio':r['yukawa']['sum_component_norms_over_total'],
                     'delta_c5_relative':a['delta_relative_to_fixed'],
                     'flavour_orthogonal_relative':a['delta_vs_fixed']['orthogonal_over_reference'],
                     'mass_basis_offdiag_relative':m['offdiag_over_fixed'],
                     **{'delta_'+k:o['angle_changes_deg'][k] for k in ('theta12_deg','theta13_deg','theta23_deg')},
                     'report':str(args.output/case/'index.html')}
                print(f'{case}: corrected C5={100*a["delta_relative_to_fixed"]:.6g}%, '
                      f'mass-basis offdiag={m["offdiag_over_fixed"]:.3e}')
            except Exception as exc:
                row={'case':case,'status':'Failed','error':f'{type(exc).__name__}: {exc}',
                     'traceback':traceback.format_exc()}
                print(f'{case}: {row["error"]}')
        items.append(row)
    save(args.output/'comparison.json',{'model':MODEL,'rows':items,
        'warning':'Saved fitted configs are imposed at mu_F; do not interpret the optional LL coupled result as a validated NLO physical prediction.'})
    keys=['case','status','lambdaT3','chi2','y2_over_y1','y2_offdiag_fraction',
          'c12_heavy_cancellation_ratio','delta_c5_relative','flavour_orthogonal_relative',
          'mass_basis_offdiag_relative','delta_theta12_deg','delta_theta13_deg','delta_theta23_deg','report','error']
    with (args.output/'comparison.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for row in items:w.writerow({k:row.get(k,'') for k in keys})
    rows_html=''.join('<tr>'+''.join(f'<td>{html.escape(str(x.get(k,"—")))}</td>' for k in
        ['case','status','delta_c5_relative','flavour_orthogonal_relative','mass_basis_offdiag_relative'])+'</tr>' for x in items)
    links=''.join(f'<li>{html.escape(x["case"])}: '+
          (f'<a href="{html.escape(x["case"]+"/index.html")}">Open numerical report</a>' if x['status']=='Success' else html.escape(x.get('error',x['status'])))+'</li>' for x in items)
    doc='''<!doctype html><html lang="en"><meta charset="utf-8"><title>UV-frozen T3-B-m1 flavour comparison</title><style>body{font:15px system-ui;max-width:1100px;margin:28px auto;padding:0 22px;line-height:1.5}th,td{padding:8px 15px;border-bottom:1px solid #aaa}table{border-collapse:collapse}</style><h1>UV-frozen T3-B-m1: 5e-9 versus 0.5</h1><p>Each result uses the saved Sobol-optimised input, UV running disabled, matched fixed-order C5 and the separate exploratory coupled LLSS/LLHH system. Missing fits are not replaced with seeds or guessed numerical results.</p><table><tr><th>Case</th><th>Status</th><th>Relative Delta C5</th><th>Orthogonal Delta C5 / C5</th><th>Mass-basis offdiagonal / M</th></tr>'''+rows_html+'</table><ul>'+links+'</ul><p><a href="comparison.json">JSON</a> | <a href="comparison.csv">CSV</a></p>'
    (args.output/'index.html').write_text(doc,encoding='utf-8')
    print('Comparison:',args.output/'index.html')
    failures=[x for x in items if x['status']=='Failed']
    if failures:
        raise RuntimeError(f'{len(failures)} case(s) failed. See comparison.json for tracebacks.')


if __name__=='__main__':main()
