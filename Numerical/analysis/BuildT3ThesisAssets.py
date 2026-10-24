"""Reproduce thesis-ready figures/tables for saved UV-frozen T3 benchmarks.

Reads saved repository JSON files only. Does not run RGEs or change any physics
code, fits, symbolic artifacts, or existing diagnostics.

Use --source snapshot to reproduce the bundled 2026-10-24 GitHub extraction.
Use --source local to compute from the current local saved repository data.
Use --source auto to select local only when ALL required sources exist.

The archived snapshot contains 11 exact samples of the 41-point Weinberg EFT
trajectory; local mode uses all saved points.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'Reports' / 'output' / 'thesis_assets'
SNAPSHOT = DEST / 'data' / 'saved_reference_snapshot.json'
CASE = '5e-9'
SOURCE = ROOT / 'output' / 'studies' / 'fermion_start_lambda_comparison' / CASE
REPRESENTATIVES = [
    'T3_dS1_1_dS2_3_dF_2_alpha_m2',   # A(-2)
    'T3_dS1_2_dS2_2_dF_1_alpha_m1',   # B(-1)
    'T3_dS1_2_dS2_2_dF_3_alpha_p1',   # C(+1)
    'T3_dS1_3_dS2_1_dF_2_alpha_m2',   # D(-2)
    'T3_dS1_3_dS2_3_dF_2_alpha_m2',   # E(-2)
]
T3B = REPRESENTATIVES[1]
KKEYS = ('hard', 'direct_running', 'combined')
BLUE, TEAL, ORANGE, PURPLE = '#2563A6', '#0C847E', '#C46B38', '#7553A5'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def model_label(name):
    m = re.fullmatch(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)', name)
    if not m:
        return name
    a, b, c, sign, n = m.groups()
    cls = {(1, 3, 2): 'A', (2, 2, 1): 'B', (2, 2, 3): 'C',
           (3, 1, 2): 'D', (3, 3, 2): 'E'}[(int(a), int(b), int(c))]
    alpha = ('+' if sign == 'p' and int(n) > 0 else '-' if sign == 'm' else '') + str(int(n))
    return f'{cls}({alpha})'


def full_config(m):
    return SOURCE / 'comparison' / 'optimised' / 'models' / m / 'fermion_start_config.json'


def full_diagnostic(m):
    return SOURCE / 'reports_existing' / 'optimised' / 'models' / m / 'running_diagnostics.json'


def to_complex(spec):
    if isinstance(spec, dict):
        r = np.asarray(spec['real'], dtype=float)
        im = np.asarray(spec.get('imag', np.zeros_like(r)), dtype=float)
        return r + 1j * im
    return np.asarray(spec, dtype=complex)


def real_matrix(mat, what):
    mat = np.asarray(mat, dtype=complex)
    if np.max(np.abs(mat.imag)) > 1e-10 * max(np.max(np.abs(mat.real)), 1e-30):
        raise ValueError(f'{what} is complex: this diagnostic uses a real orthogonal C5 eigenbasis; extend with Takagi first.')
    return mat.real


def calc_flavour(config, c5_matrices):
    ordinary = config['base_state']['ordinary']
    y1 = real_matrix(to_complex(ordinary['y1']), 'Y1')
    y2 = real_matrix(to_complex(ordinary['y2']), 'Y2')
    mf = real_matrix(to_complex(ordinary['MF']), 'MF')
    if not np.allclose(mf, np.diag(np.diag(mf)), rtol=0, atol=1e-8):
        raise ValueError('Fermion masses are not in a diagonal basis.')
    masses = np.diag(mf)
    if np.min(masses) <= 0:
        raise ValueError('Expected positive fermion masses.')
    W = np.diag(1. / masses)
    unsym = y1 @ W @ y2.T
    F = (unsym + unsym.T) / 2
    terms = [
        (np.outer(y1[:, r], y2[:, r]) + np.outer(y2[:, r], y1[:, r])) / (2*masses[r])
        for r in range(3)
    ]
    assert np.allclose(sum(terms), F, rtol=1e-12, atol=1e-22)
    K = real_matrix(np.asarray(c5_matrices['combined'], dtype=float), 'matched C5')
    _, U = np.linalg.eigh(K)
    eig = np.linalg.eigvalsh(K)
    U = U[:, np.argsort(np.abs(eig))]
    U = np.array(U, dtype=float)
    transform = lambda m: U.T @ m @ U
    fraction = lambda m: float(np.linalg.norm(m - np.diag(np.diag(m))) / np.linalg.norm(m))
    projected = transform(F)
    source_metrics = {}
    for key in KKEYS:
        value = real_matrix(c5_matrices[key], f'C5 {key}')
        source_metrics[key] = {
            'norm_gev_inverse': float(np.linalg.norm(value)),
            'offdiag_fraction_in_C5_basis': fraction(transform(value)),
        }
    proj_alpha = float(np.vdot(F, K).real / np.vdot(F, F).real)
    shape_residual = float(np.linalg.norm(K - proj_alpha * F) / np.linalg.norm(K))
    return {
        'C12': F, 'unsym': unsym, 'terms': terms, 'U': U, 'C12_mass_basis': projected,
        'unsym_mass_basis': transform(unsym),
        'terms_mass_basis': [transform(t) for t in terms],
        'offdiag_fraction_C12_in_C5_basis': fraction(projected),
        'offdiag_fraction_unsymmetrised_in_C5_basis': fraction(transform(unsym)),
        'offdiag_fraction_direct_in_C5_basis': source_metrics['direct_running']['offdiag_fraction_in_C5_basis'],
        'projected_shape_residual_C5_to_C12': shape_residual,
        'heavy_cancellation_ratio': float(sum(np.linalg.norm(t) for t in terms) / np.linalg.norm(F)),
        'c5_metrics': source_metrics,
    }


def load_local():
    snap = read(SNAPSHOT)
    fits = []
    for entry in snap['fits']:
        m = entry['model']
        path = SOURCE / 'comparison' / 'optimised' / 'models' / m / 'evaluation.json'
        payload = read(path)
        row = {'model': m, 'status': payload['status'], 'chi2': payload['chi2'], **payload['prediction']}
        fits.append(row)
    flavour = []
    for m in REPRESENTATIVES:
        config = read(full_config(m))
        diagnostic = read(full_diagnostic(m))
        mats = {name: np.asarray(diagnostic['c5_threshold_contributions'][name]['real']) for name in KKEYS}
        result = calc_flavour(config, mats)
        flavour.append({'model': m, **{k: result[k] for k in ('offdiag_fraction_C12_in_C5_basis',
            'offdiag_fraction_unsymmetrised_in_C5_basis', 'offdiag_fraction_direct_in_C5_basis',
            'projected_shape_residual_C5_to_C12', 'heavy_cancellation_ratio')}})
    ref = read(full_config(T3B))
    diag = read(full_diagnostic(T3B))
    mats = {k: diag['c5_threshold_contributions'][k]['real'] for k in KKEYS}
    w = diag['final_running']
    trajectories = [dict(mu_gev=float(mu), **{
        key: float(w[key][i]) for key in (
            'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
            'delta_m21_sq_ev2', 'delta_m3l_sq_ev2')})
        for i, mu in enumerate(w['mu_gev'])]
    return {'metadata': {**snap['metadata'], 'source_mode': 'local',
                         'weinberg_samples': f'all {len(trajectories)} saved points'},
            'fits': fits, 'flavour_five': flavour, 't3b_reference_config': ref,
            't3b_matched_c5': mats, 'weinberg_running_t3b': trajectories}


def plot_setup():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.titlesize': 12, 'axes.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'figure.dpi': 130, 'savefig.dpi': 200,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.grid': False,
    })


def savefig(tag, fig, description, manifest):
    folder = DEST / 'figures'
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f'{tag}.png', bbox_inches='tight', facecolor='white')
    fig.savefig(folder / f'{tag}.pdf', bbox_inches='tight', facecolor='white')
    manifest.append({'id': tag, 'png': f'figures/{tag}.png',
                     'pdf': f'figures/{tag}.pdf', 'description': description})
    plt.close(fig)


def fig_flow(manifest):
    fig, ax = plt.subplots(figsize=(10.3, 3.1))
    ax.set_xlim(0, 10); ax.set_ylim(0, 3); ax.axis('off')
    blocks = [
        (0.03, 'UV theory', '$Y_1, Y_2, M_F, \\lambda_{T3}$', BLUE),
        (2.55, 'Fermion matching', '$C_{12}(\\mu_F)$', TEAL),
        (5.07, 'Scalar EFT', '$C_{12}\\to C_5$', PURPLE),
        (7.59, 'SM + Weinberg', '$C_5(\\mu)\\to m_\\nu$', ORANGE),
    ]
    for i,(x,title,eq,color) in enumerate(blocks):
        ax.add_patch(FancyBboxPatch((x,1.05),2.28,1.2,boxstyle='round,pad=0.12,rounding_size=0.08',
                                    facecolor='white',edgecolor=color,lw=2))
        ax.text(x+1.14,1.85,title,ha='center',va='center',weight='bold',fontsize=10)
        ax.text(x+1.14,1.4,eq,ha='center',va='center',fontsize=11)
        if i<3:
            ax.annotate('',xy=(x+2.49,1.65),xytext=(x+2.32,1.65),
                        arrowprops={'arrowstyle':'->','lw':2,'color':'#475569'})
    ax.text(3.69,0.65,r'$\mu_F$: remove all $F_r$',ha='center',size=9,color='#475569')
    ax.text(6.23,0.65,r'$\mu_S$: remove scalars',ha='center',size=9,color='#475569')
    ax.text(5,2.77,'Threshold matching and running in the hierarchical T3 calculation',
            ha='center',size=13,weight='bold')
    savefig('fig01_t3_effective_theory_flow',fig,'Schematic: production EFT stages and operator mixing (not numerical data).',manifest)


def fig_scaling(manifest):
    fig, ax = plt.subplots(figsize=(7.0,4.4))
    t=np.linspace(-2,8,201)
    ax.plot(t,t,lw=2.4,color=ORANGE,label=r'$\lambda_{T3}/\lambda_{T3,0}=\gamma$')
    ax.plot(t,-t/2,lw=2.4,color=BLUE,label=r'$Y_{1,2}/Y_{1,2,0}=\gamma^{-1/2}$')
    ax.plot(t,np.zeros_like(t),lw=2.4,color=TEAL,label=r'$(\lambda_{T3}Y_1Y_2)/(\lambda_{T3}Y_1Y_2)_0=1$')
    ax.axvline(8,linestyle=':',color='#475569',alpha=.7)
    ax.set(xlabel=r'$\log_{10}\gamma$',ylabel=r'$\log_{10}$(ratio to original)',
           title='Analytic leading-order matching degeneracy',xlim=(-2,8.2))
    ax.grid(alpha=.17)
    ax.legend(loc='lower left',fontsize=8,frameon=False)
    savefig('fig02_analytic_coupling_scaling',fig,'Analytic scaling relation, not an RGE trajectory or parameter scan.',manifest)


def fig_fit(fits,manifest):
    fig,ax=plt.subplots(figsize=(11.5,4.6))
    xs=np.arange(len(fits)); vals=[max(float(x['chi2']),1e-35) for x in fits]
    colors=[ORANGE if v>=.1 else BLUE for v in vals]
    ax.bar(xs,vals,color=colors,alpha=.86,width=.78)
    ax.set_yscale('log'); ax.set_ylim(1e-28,1e3)
    ax.set_xticks(xs,[model_label(x['model']) for x in fits],rotation=55,ha='right')
    ax.set_ylabel(r'Saved fit $\chi^2$ (logarithmic)')
    ax.set_title('UV-frozen optimised neutrino fits across 16 T3 models')
    ax.grid(axis='y',alpha=.2,which='major')
    savefig('fig03_optimised_fit_quality_16_models',fig,'Saved fit objective for all 16 representations; poor fits shown without censoring.',manifest)


def fig_flavour(rows,manifest):
    fig,ax=plt.subplots(figsize=(8.4,4.7))
    x=np.arange(len(rows));w=.34
    offu=np.array([100*r['offdiag_fraction_unsymmetrised_in_C5_basis'] for r in rows]);
    offs=np.array([100*r['offdiag_fraction_C12_in_C5_basis'] for r in rows]);
    b1=ax.bar(x-w/2,offu,w,color=BLUE,label=r'Unsymmetrised $Y_1^*M_F^{-1}Y_2^\dagger$')
    b2=ax.bar(x+w/2,offs,w,color=TEAL,label=r'Symmetrised $C_{12}$')
    for bars in (b1,b2):
        for b in bars: ax.text(b.get_x()+b.get_width()/2,b.get_height()+1,f'{b.get_height():.1f}%',ha='center',fontsize=9)
    ax.set_xticks(x,[model_label(r['model']) for r in rows]);ax.set_ylim(0,102)
    ax.set_ylabel('Off-diagonal Frobenius norm fraction (%)')
    ax.set_title(r'Yukawa products in each model\'s matched $C_5$ eigenbasis')
    ax.legend(fontsize=8,frameon=False,loc='upper right')
    ax.grid(axis='y',alpha=.16)
    savefig('fig04_five_class_flavour_alignment',fig,'Five fitted class representatives; real-Yukawa basis of the corresponding matched C5.',manifest)


def fig_matrix(algebra,manifest):
    mats=[algebra['unsym_mass_basis'],algebra['C12_mass_basis']]
    mmax=max(float(np.max(np.abs(m))) for m in mats)/1e-6
    fig,axs=plt.subplots(1,2,figsize=(10.6,4.45),constrained_layout=True)
    for ax,arr,name in zip(axs,mats,[r'Unsymmetrised $A$',r'Symmetrised $C_{12}=(A+A^T)/2$']):
        arr=arr/1e-6
        im=ax.imshow(arr,cmap='RdBu_r',vmin=-mmax,vmax=mmax,interpolation='nearest')
        for i in range(3):
            for j in range(3):
                ax.text(j,i,f'{arr[i,j]:+.3f}',ha='center',va='center',fontsize=9,
                        color='white' if abs(arr[i,j])>0.6*mmax else '#15202b')
        ax.set_xticks([0,1,2],[1,2,3]);ax.set_yticks([0,1,2],[1,2,3]);
        ax.set_xlabel('Neutrino mass eigenstate index');ax.set_ylabel('Neutrino mass eigenstate index')
        ax.set_title(name)
    fig.colorbar(im,ax=axs,shrink=.78,label=r'Entries ($10^{-6}$ GeV$^{-1}$)',pad=.03)
    fig.suptitle('T3-B(-1): Yukawa-product symmetrisation in the matched C5 eigenbasis',fontsize=12)
    savefig('fig05_t3b_yukawa_product_matrices',fig,'Real fitted T3-B(-1), before/after Majorana symmetrisation, fixed matched C5 basis.',manifest)


def fig_cancellation(algebra,manifest):
    fig,ax=plt.subplots(figsize=(8.2,4.9))
    ix=[(0,1),(0,2),(1,2)]
    x=np.arange(3);w=.19
    for r,(color,name) in enumerate([(BLUE,'Heavy 1'),(TEAL,'Heavy 2'),(ORANGE,'Heavy 3')]):
        val=[algebra['terms_mass_basis'][r][i,j]/1e-6 for i,j in ix]
        ax.bar(x+(r-1.5)*w,val,w,color=color,label=name)
    sums=[algebra['C12_mass_basis'][i,j]/1e-6 for i,j in ix]
    ax.bar(x+1.5*w,sums,w,color=PURPLE,label='Sum')
    ax.axhline(0,color='#475569',lw=.9)
    ax.set_xticks(x,[r'$(12)$',r'$(13)$',r'$(23)$'])
    ax.set_ylabel(r'Off-diagonal contribution ($10^{-6}$ GeV$^{-1}$)')
    ax.set_xlabel('Flavour pair in the neutrino mass basis')
    ax.set_title('T3-B(-1): heavy-fermion terms cancel in the symmetrised flavour tensor')
    ax.legend(ncol=4,frameon=False,fontsize=9)
    ax.grid(axis='y',alpha=.2)
    savefig('fig06_t3b_heavy_state_cancellations',fig,'Three mass-weighted heavy-state terms and their sum in C5 eigenbasis.',manifest)


def fig_c5_contributions(algebra,manifest):
    data=algebra['c5_metrics']; keys=['hard','direct_running','combined'];
    vals=[data[k]['norm_gev_inverse']/1e-15 for k in keys]
    fig,ax=plt.subplots(figsize=(6.5,4.6))
    b=ax.bar(range(3),vals,color=[BLUE,TEAL,PURPLE],width=.56)
    for v in b:
        ax.text(v.get_x()+v.get_width()/2,v.get_height()+.03,f'{v.get_height():.3f}',ha='center',fontsize=10)
    ax.set_xticks(range(3),['Hard matching','Direct running','Combined'])
    ax.set_ylabel(r'Frobenius norm ($10^{-15}$ GeV$^{-1}$)')
    ax.set_title(r'T3-B(-1): contributions to matched $C_5(\mu_S)$')
    ax.grid(axis='y',alpha=.2)
    fig.text(.50,.015,'Combined norm is not the sum of norms: matrices interfere.',ha='center',fontsize=8,color='#475569')
    savefig('fig07_t3b_c5_threshold_components',fig,'Actual saved hard, direct-running and combined C5 matrices at the scalar threshold.',manifest)


def fig_weinberg(data,manifest):
    x=np.array([r['mu_gev'] for r in data],dtype=float)
    fig,ax=plt.subplots(figsize=(7.2,4.6))
    keys=['delta_m21_sq_ev2','delta_m3l_sq_ev2']
    names=[r'$\Delta m^2_{21}$',r'$\Delta m^2_{3\ell}$']
    for k,name,c in zip(keys,names,[BLUE,ORANGE]):
        y=np.array([r[k] for r in data]);ax.plot(x,100*(y/y[0]-1),'o-',color=c,lw=2,ms=3,label=name)
    ax.set_xscale('log');ax.invert_xaxis();ax.set_xlim(x.max()*1.04,x.min()*.96)
    ax.set_xlabel(r'Running scale $\mu$ (GeV)')
    ax.set_ylabel('Change from scalar-threshold value (%)')
    ax.set_title('T3-B(-1): mass-squared splittings during SM + Weinberg running')
    ax.legend(frameon=False);ax.grid(alpha=.2)
    savefig('fig08_t3b_weinberg_mass_splitting_running',fig,'Actual saved Weinberg-stage RGE solution; UV frozen, 5e-9, 1000 to 100 GeV.',manifest)

    fig,ax=plt.subplots(figsize=(7.2,4.6))
    for key,txt,color in [('sin2_theta12',r'$\theta_{12}$',BLUE),
                          ('sin2_theta13',r'$\theta_{13}$',ORANGE),
                          ('sin2_theta23',r'$\theta_{23}$',TEAL)]:
        ang=np.degrees(np.arcsin(np.sqrt([r[key] for r in data])))
        ax.plot(x,ang-ang[0],'o-',lw=2,ms=3,color=color,label=txt)
    ax.set_xscale('log');ax.invert_xaxis();ax.set_xlim(x.max()*1.04,x.min()*.96)
    ax.set_xlabel(r'Running scale $\mu$ (GeV)')
    ax.set_ylabel('Change from scalar-threshold angle (degrees)')
    ax.set_title('T3-B(-1): mixing-angle drift during SM + Weinberg running')
    ax.legend(frameon=False);ax.grid(alpha=.2)
    savefig('fig09_t3b_weinberg_angle_running',fig,'Actual saved Weinberg-stage mixing observables, not the scalar-only C5 trajectory.',manifest)


def write_csv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def sci_tex(value,digits=3):
    if value is None:return '--'
    if value == 0:return '$0$'
    val=float(value)
    e=int(np.floor(np.log10(abs(val))))
    if e>=-2 and e<=3:return f'${val:.{digits}f}$'
    coef=val/10**e
    return f'${coef:.{digits}g}\\times10^{{{e}}}$'


def latex_table(path,caption,label,columns,rows):
    # columns is a sequence of (data-key, display-label, format-function)
    header=' & '.join(col[1] for col in columns)
    body=[' & '.join(col[2](row.get(col[0])) for col in columns)+r' \\' for row in rows]
    colspec='l'+'r'*(len(columns)-1)
    text=(r'\begin{table}[htbp]'+'\n'+r'\centering'+'\n'+
          r'\caption{'+caption+'}\n'+r'\label{'+label+'}\n'+
          r'\begin{tabular}{'+colspec+'}\n'+r'\toprule'+'\n'+header+r' \\'+'\n'+
          r'\midrule'+'\n'+'\n'.join(body)+'\n'+r'\bottomrule'+'\n'+r'\end{tabular}'+'\n'+r'\end{table}'+'\n')
    path.write_text(text,encoding='utf-8')


def tables(data,algebra):
    out=DEST/'tables';out.mkdir(exist_ok=True,parents=True)
    fit=data['fits']; fl=data['flavour_five']
    fkeys=['model','chi2','delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23']
    write_csv(out/'table01_uv_frozen_fits_16_models.csv',fit,fkeys)
    latex_table(out/'table01_uv_frozen_fits_16_models.tex',
        r'UV-frozen $\lambda_{T3}=5\times10^{-9}$ optimised fit results. All 16 points report Success; small $\chi^2$ values reflect fitted targets, not model evidence.',
        'tab:t3-uv-frozen-fits',
        [('model','Model',lambda s:model_label(s)),('chi2',r'$\chi^2$',lambda v:sci_tex(v,3)),
         ('delta_m21_sq_ev2',r'$\Delta m_{21}^2$ (eV$^2$)',lambda v:sci_tex(v,4)),
         ('delta_m3l_sq_ev2',r'$\Delta m_{3\ell}^2$ (eV$^2$)',lambda v:sci_tex(v,4))],fit)
    fk=['model','offdiag_fraction_unsymmetrised_in_C5_basis','offdiag_fraction_C12_in_C5_basis',
        'projected_shape_residual_C5_to_C12','heavy_cancellation_ratio']
    write_csv(out/'table02_five_class_flavour_alignment.csv',fl,fk)
    latex_table(out/'table02_five_class_flavour_alignment.tex',
        r'Mass-basis off-diagonal norm fractions of the unsymmetrised Yukawa product and $C_{12}$ for one fitted benchmark per T3 class.',
        'tab:t3-flavour-alignment',
        [('model','Model',lambda s:model_label(s)),
         ('offdiag_fraction_unsymmetrised_in_C5_basis',r'Unsym. (\%)',lambda x:f'${100*x:.2f}$'),
         ('offdiag_fraction_C12_in_C5_basis',r'$C_{12}$ (\%)',lambda x:f'${100*x:.2f}$'),
         ('projected_shape_residual_C5_to_C12',r'Shape residual (\%)',lambda x:f'${100*x:.2f}$'),
         ('heavy_cancellation_ratio','Cancellation',lambda x:f'${x:.2f}$')],fl)
    cm=algebra['c5_metrics'];rows=[{'component':key,**cm[key]} for key in KKEYS]
    write_csv(out/'table03_t3b_c5_decomposition.csv',rows,
              ['component','norm_gev_inverse','offdiag_fraction_in_C5_basis'])
    latex_table(out/'table03_t3b_c5_decomposition.tex',
        r'Hard, direct-running, and combined matched $C_5$ for the UV-frozen T3-B(-1) benchmark.',
        'tab:t3-c5-decomposition',
        [('component','Component',lambda s:s.replace('_',' ').capitalize()),
         ('norm_gev_inverse',r'$\|C_5\|_F$ (GeV$^{-1}$)',lambda x:sci_tex(x,4)),
         ('offdiag_fraction_in_C5_basis',r'Off-diagonal (\%)',lambda x:f'${100*x:.3f}$')],rows)
    w=data['weinberg_running_t3b']; endpoints=[]
    for name in ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23'):
        a=float(w[0][name]);b=float(w[-1][name]);
        endpoints.append({'quantity':name,'muS':a,'muLow':b,'absolute_change':b-a,'fractional_change':(b-a)/a})
    write_csv(out/'table04_t3b_weinberg_running_endpoints.csv',endpoints,
              ['quantity','muS','muLow','absolute_change','fractional_change'])
    latex_table(out/'table04_t3b_weinberg_running_endpoints.tex',
        r'Saved T3-B(-1) Weinberg-stage predictions at $\mu_S=1$ TeV and $\mu_{\rm low}=100$ GeV. Angles are presented as $\sin^2\theta_{ij}$.',
        'tab:t3-weinberg-endpoints',
        [('quantity','Observable',lambda s:s.replace('_',r'\_')),
         ('muS',r'1 TeV',lambda x:sci_tex(x,6)),('muLow','100 GeV',lambda x:sci_tex(x,6)),
         ('fractional_change',r'Change (\%)',lambda x:f'${100*x:+.4f}$')],endpoints)
    cancellation=[]
    for i,j in ((0,1),(0,2),(1,2)):
        item={'index_pair':f'{i+1}{j+1}',**{f'heavy_{r+1}':float(M[i,j]) for r,M in enumerate(algebra['terms_mass_basis'])},
              'sum':float(algebra['C12_mass_basis'][i,j])}
        cancellation.append(item)
    write_csv(out/'table05_t3b_heavy_cancellation.csv',cancellation,['index_pair','heavy_1','heavy_2','heavy_3','sum'])
    latex_table(out/'table05_t3b_heavy_cancellation.tex',
        r'Individual heavy-state contributions to off-diagonal $C_{12}$ entries in the eigenbasis of the full matched $C_5$, units $10^{-6}$ GeV$^{-1}$.',
        'tab:t3-heavy-cancellation',
        [('index_pair','Entry',lambda s:s),('heavy_1','$r=1$',lambda x:f'${x*1e6:+.4f}$'),
         ('heavy_2','$r=2$',lambda x:f'${x*1e6:+.4f}$'),('heavy_3','$r=3$',lambda x:f'${x*1e6:+.4f}$'),
         ('sum','Sum',lambda x:f'${x*1e6:+.4f}$')],cancellation)
    return endpoints


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',choices=['auto','local','snapshot'],default='auto')
    args=p.parse_args(argv)
    local_ok=all([*[full_config(m).is_file() and full_diagnostic(m).is_file() for m in REPRESENTATIVES],
                  *[(SOURCE/'comparison'/'optimised'/'models'/r['model']/'evaluation.json').is_file() for r in read(SNAPSHOT)['fits']]])
    if args.source=='local' and not local_ok:
        raise SystemExit('Missing local saved JSON inputs. Run --source snapshot or restore repository saved results.')
    mode='local' if (args.source=='local' or (args.source=='auto' and local_ok)) else 'snapshot'
    data=load_local() if mode=='local' else read(SNAPSHOT)
    algebra=calc_flavour(data['t3b_reference_config'],data['t3b_matched_c5'])
    rows=data['flavour_five']
    here=next(r for r in rows if r['model']==T3B)
    if not np.isclose(algebra['offdiag_fraction_C12_in_C5_basis'],here['offdiag_fraction_C12_in_C5_basis'],rtol=1e-7,atol=1e-9):
        raise ValueError('Source inconsistency: local/snapshot fitted C12 alignment fails.')
    if not np.isclose(algebra['offdiag_fraction_unsymmetrised_in_C5_basis'],here['offdiag_fraction_unsymmetrised_in_C5_basis'],rtol=1e-7,atol=1e-9):
        raise ValueError('Source inconsistency: unsymmetrised product alignment fails.')
    if not np.allclose(sum(algebra['terms']),algebra['C12'],rtol=1e-13,atol=1e-21):
        raise AssertionError('C12 is not the sum of individual heavy states.')
    plot_setup()
    manifest=[]
    fig_flow(manifest)
    fig_scaling(manifest)
    fig_fit(data['fits'],manifest)
    fig_flavour(rows,manifest)
    fig_matrix(algebra,manifest)
    fig_cancellation(algebra,manifest)
    fig_c5_contributions(algebra,manifest)
    fig_weinberg(data['weinberg_running_t3b'],manifest)
    endpoints=tables(data,algebra)
    summary={'status':'Success','mode':mode,'input_source':data['metadata'],
             'fit_count':len(data['fits']),'representative_count':len(rows),
             't3b_offdiag_C12_in_C5_basis':algebra['offdiag_fraction_C12_in_C5_basis'],
             't3b_offdiag_unsym_in_C5_basis':algebra['offdiag_fraction_unsymmetrised_in_C5_basis'],
             't3b_heavy_cancellation_ratio':algebra['heavy_cancellation_ratio'],
             't3b_hard_and_direct_norms':algebra['c5_metrics'],
             'weinberg_sample_count':len(data['weinberg_running_t3b']),
             'weinberg_endpoints':endpoints,
             'figures':manifest,
             'warning':'These figures describe a UV-frozen fitted sample. No UV-running sensitivity or high-lambda=0.5 prediction is inferred.'}
    write_json(DEST/'data'/'summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('status','mode','fit_count','representative_count',
        't3b_offdiag_C12_in_C5_basis','t3b_heavy_cancellation_ratio','weinberg_sample_count')},indent=2))
    print(f'Generated {len(manifest)} PDF+PNG figure pairs and 5 CSV+LaTeX table pairs in {DEST}')

if __name__=='__main__':
    main()
