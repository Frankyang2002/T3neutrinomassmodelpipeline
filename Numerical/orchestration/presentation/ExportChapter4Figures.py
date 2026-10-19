"""Export Chapter 4 thesis figures/tables from EXISTING saved T3 results.

Run from the local repository root:
    python Numerical/orchestration/ExportChapter4Figures.py

Never performs RGE evolution, matching or fitting. Does not touch GitHub.
Requires: numpy, matplotlib (Python standard-library csv/json otherwise).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Patch
from matplotlib.lines import Line2D

CLASSES = {(1,3,2):'A', (2,2,1):'B', (2,2,3):'C', (3,1,2):'D', (3,3,2):'E'}
ORDER = ['A-m4','A-m2','A-p0','B-m3','B-m1','B-p1','C-m3','C-m1','C-p1',
         'D-m2','D-p0','D-p2','E-m4','E-m2','E-p0','E-p2']
SIX = ['A-m2','A-p0','B-m3','B-p1','D-m2','D-p0']
CLASS_COLORS = {'A':'#1479B0','B':'#D55E00','C':'#17A673','D':'#C67AA8','E':'#E5A000'}
SIX_COLORS = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']
OBS = [('delta_m21_sq_ev2',r'$\Delta m_{21}^{2}$',1e5,r'$10^{-5}\,\mathrm{eV}^{2}$'),
       ('delta_m3l_sq_ev2',r'$\Delta m_{31}^{2}$',1e3,r'$10^{-3}\,\mathrm{eV}^{2}$')]
PATTERN = re.compile(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)$')


def label(full):
    m = PATTERN.fullmatch(full)
    if not m:
        raise ValueError(f'Unknown model identifier: {full}')
    dims = tuple(map(int,m.group(1,2,3)))
    return f'{CLASSES[dims]}-{m.group(4)}{m.group(5)}'


def require(path):
    if not path.is_file():
        raise FileNotFoundError(f'Missing saved input: {path}')
    return path


def csv_rows(path):
    with require(path).open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))


def json_load(path):
    return json.loads(require(path).read_text(encoding='utf-8'))


def success(rows, name):
    bad = [r for r in rows if r.get('status') not in ('Success','AlreadyComplete')]
    if bad:
        raise ValueError(f'{name}: {len(bad)} non-successful rows. Refusing to silently omit them.')


def select_unique(rows, mode=None):
    ret={}
    for row in rows:
        if mode is not None and row.get('mode') != mode:
            continue
        key=(label(row['model']),float(row['factor']))
        if key in ret: raise ValueError(f'Duplicate result: {key}')
        ret[key]=row
    return ret


def require_grid(data, models, factors, what):
    for m in models:
        for f in factors:
            if (m,f) not in data:
                raise ValueError(f'Missing {what}: {m}, factor={f}')


def write_csv(path, fields, rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def save(fig, out, stem):
    fig.savefig(out / f'{stem}.png',dpi=240,bbox_inches='tight')
    fig.savefig(out / f'{stem}.pdf',bbox_inches='tight')
    plt.close(fig)


def spans(ax):
    for p in [3,6,9,12]:ax.axvline(p-0.5,lw=0.7,color='0.8',zorder=0)


def models_axis(ax):
    ax.set_xticks(range(16),ORDER,rotation=55,ha='right',fontsize=8)
    spans(ax);ax.grid(axis='y',alpha=.18)


def baseline_outputs(ferm, figs, tables):
    base={m:ferm[m,1.0] for m in ORDER}
    t=[]
    for m in ORDER:
        r=base[m]
        t.append({'model':m,'delta_m21_sq_ev2':r['running_delta_m21_sq_ev2'],
                  'delta_m31_sq_ev2_NO':r['running_delta_m3l_sq_ev2']})
    write_csv(tables/'table_4_1_baseline_16.csv',list(t[0]),t)
    fig, axs=plt.subplots(2,1,figsize=(12,8),constrained_layout=True)
    colors=[CLASS_COLORS[m[0]] for m in ORDER]
    for ax,(key,tex,scale,unit) in zip(axs,OBS):
        vals=[float(base[m]['running_'+key])*scale for m in ORDER]
        ax.bar(range(16),vals,color=colors,width=.78)
        ax.set_ylabel(f'{tex} [{unit}]');models_axis(ax)
    fig.suptitle('Common-input predictions across T3 configurations')
    save(fig,figs,'fig_4_1_representation_comparison')

    t=[]
    fig,axs=plt.subplots(2,1,figsize=(12,8),constrained_layout=True)
    for m in ORDER:
        r=base[m]
        t.append({'model':m,'rg_percent_delta_m21_sq':r['relative_percent_delta_m21_sq_ev2'],
                  'rg_percent_delta_m31_sq_NO':r['relative_percent_delta_m3l_sq_ev2']})
    write_csv(tables/'table_4_2_rg_corrections_16.csv',list(t[0]),t)
    for ax,(key,tex,_,_) in zip(axs,OBS):
        vals=[float(base[m]['relative_percent_'+key]) for m in ORDER]
        ax.bar(range(16),vals,color=colors,width=.78)
        ax.axhline(0,color='black',lw=.9)
        ax.set_ylabel(rf'$\delta_{{\mathrm{{RG}}}}({tex.strip(chr(36))})$ [\%]')
        models_axis(ax)
    fig.suptitle('Effect of parameter running relative to frozen controls')
    save(fig,figs,'fig_4_2_rg_corrections')


def scan_figures(ferm,couplings,figs):
    factors=[.5,.75,1.,1.5,2.]
    require_grid(ferm,SIX,factors,'heavy-fermion scan')
    # Figure 4.3: running observables; the fermion matching scale is also rescaled.
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.4),constrained_layout=True)
    for ax,(key,tex,scale,unit) in zip(axs,OBS):
        for m,c in zip(SIX,SIX_COLORS):
            ax.plot(factors,[float(ferm[m,f]['running_'+key])*scale for f in factors],
                    marker='o',label=m,color=c,lw=1.7,markersize=4)
        ax.axvline(1,color='0.4',lw=.9,ls='--');ax.set_xlabel(r'$r$ (scales $M_F$ and $\mu_F$)')
        ax.set_ylabel(f'{tex} [{unit}]');ax.grid(alpha=.2)
    axs[1].legend(ncol=2,fontsize=8,loc='best')
    save(fig,figs,'fig_4_3_fermion_mass_dependence')
    # Figure 4.4: signed percent differences.
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.4),constrained_layout=True)
    for ax,(key,tex,_,_) in zip(axs,OBS):
        for m,c in zip(SIX,SIX_COLORS):
            ax.plot(factors,[float(ferm[m,f]['relative_percent_'+key]) for f in factors],
                    marker='o',label=m,color=c,lw=1.7,markersize=4)
        ax.axvline(1,color='0.4',lw=.9,ls='--');ax.axhline(0,color='0.25',lw=.8)
        ax.set_xlabel(r'$r$ (scales $M_F$ and $\mu_F$)')
        ax.set_ylabel(rf'$\delta_{{\mathrm{{RG}}}}({tex.strip(chr(36))})$ [\%]');ax.grid(alpha=.2)
    axs[1].legend(ncol=2,fontsize=8,loc='best')
    save(fig,figs,'fig_4_4_fermion_rg_dependence')
    # Figure 4.5: show distinct power-law predictions separately, with 6-model points.
    fig,axs=plt.subplots(1,2,figsize=(12.5,4.6),constrained_layout=True)
    for ax,mode,power in zip(axs,('yukawa','lambda'),(4,2)):
        data=couplings[mode]
        require_grid(data,SIX,factors,mode+' scaling scan')
        for m,c in zip(SIX,SIX_COLORS):
            ref=float(data[m,1.0]['running_delta_m21_sq_ev2'])
            vals=[float(data[m,f]['running_delta_m21_sq_ev2'])/ref for f in factors]
            ax.plot(factors,vals,marker='o',markersize=3.7,color=c,lw=1.35,label=m)
        xx=np.linspace(.5,2.,150)
        ax.plot(xx,xx**power,ls='--',lw=2,color='black',label=rf'$k^{power}$')
        ax.axvline(1,color='0.5',lw=.8,ls=':')
        ax.set(xlabel=r'$k$',ylabel=r'$\Delta m_{21}^{2}(k)/\Delta m_{21}^{2}(1)$')
        ax.set_title('Simultaneous Yukawa rescaling' if mode=='yukawa' else r'$\lambda_{T3}$ rescaling')
        ax.grid(alpha=.2)
    axs[1].legend(fontsize=8,ncol=2)
    save(fig,figs,'fig_4_5_coupling_dependence')


def ellipse_patch(centre,cov,delta,ax,**style):
    values,vec=np.linalg.eigh(np.asarray(cov,dtype=float))
    if np.any(values <= 0):return
    angle=math.degrees(math.atan2(vec[1,1],vec[0,1]))
    ax.add_patch(Ellipse(centre,2*math.sqrt(values[1]*delta),2*math.sqrt(values[0]*delta),angle=angle,**style))


def correlation_figure(root,figs):
    # Use three selected common-input trajectories; they are NOT independently fitted.
    selected=['A-m2','B-m3','D-p0']
    sens=root/'Reports/output/full/common_input_sensitivity/smallY_smallL'
    index={}
    for directory in sens.iterdir():
        if directory.is_dir():
            try: index[label(directory.name)]=directory
            except ValueError: pass
    data={}
    for m in selected:
        if m not in index:raise FileNotFoundError('Missing common-input sensitivity for '+m)
        j=json_load(index[m]/'sensitivity_diagnostics.json')
        if j['status']!='Success':raise ValueError('Sensitivity not successful for '+m)
        data[m]=j
    target=data[selected[0]]['target']
    names=target['observable_names'];mu=np.array(target['central_values']);cov=np.array(target['covariance'])
    specs=[('sin2_theta12','delta_m21_sq_ev2',r'$\sin^2\theta_{12}$',r'$\Delta m_{21}^{2}$ [$10^{-5}\,\mathrm{eV}^{2}$]',1e5),
           ('sin2_theta23','delta_m3l_sq_ev2',r'$\sin^2\theta_{23}$',r'$\Delta m_{31}^{2}$ [$10^{-3}\,\mathrm{eV}^{2}$]',1e3)]
    fig,axs=plt.subplots(1,2,figsize=(13,5.3),constrained_layout=True)
    model_colors={m:SIX_COLORS[SIX.index(m)] for m in selected}
    for ax,(x,y,xlab,ylab,scale) in zip(axs,specs):
        i=names.index(x);j=names.index(y)
        centre=np.array([mu[i],mu[j]*scale]);cv=cov[np.ix_([i,j],[i,j])].copy()
        cv[1,:]*=scale;cv[:,1]*=scale
        ellipse_patch(centre,cv,11.83,ax,facecolor='#bcd5f1',edgecolor='#6a92c1',alpha=.35,lw=1.25,zorder=0)
        ellipse_patch(centre,cv,2.30,ax,facecolor='#85addb',edgecolor='#487db7',alpha=.35,lw=1.25,zorder=0)
        for m in selected:
            s=data[m];color=model_colors[m]
            for direction,ls in [('T1','-'),('T2','--'),('T3',':')]:
                pts=[p for p in s['directions'][direction]['points'] if p.get('status')=='Success' and 'prediction' in p]
                pts.sort(key=lambda p:float(p['gamma']))
                if not pts:continue
                xx=[float(p['prediction'][x]) for p in pts]
                yy=[float(p['prediction'][y])*scale for p in pts]
                ax.plot(xx,yy,ls=ls,lw=1.3,color=color,alpha=.9)
                cent=min(pts,key=lambda p:abs(float(p['gamma'])-1))
                ax.scatter([float(cent['prediction'][x])],[float(cent['prediction'][y])*scale],color=color,s=20,zorder=4)
            # Descriptive empirical covariance envelope, centred at baseline, not statistical confidence.
            all_pts=[]
            for d in ('T1','T2','T3'):
                all_pts += [p for p in s['directions'][d]['points'] if p.get('status')=='Success'
                            and 'prediction' in p and not math.isclose(float(p['gamma']),1,abs_tol=1e-9)]
            if len(all_pts)>=3:
                z=np.array([[float(p['prediction'][x]),float(p['prediction'][y])*scale] for p in all_pts])
                baseline=[]
                for d in ('T1',):
                    ps=s['directions'][d]['points'];c=min(ps,key=lambda p:abs(float(p['gamma'])-1))
                    baseline=[float(c['prediction'][x]),float(c['prediction'][y])*scale]
                co=(z-np.array(baseline)).T@(z-np.array(baseline))/(len(z)-1)
                ellipse_patch(baseline,co,4,ax,facecolor='none',edgecolor=color,linestyle='-.',lw=.9,alpha=.7)
        ax.set_xlabel(xlab);ax.set_ylabel(ylab);ax.grid(alpha=.12)
    handles=[Patch(facecolor='#85addb',alpha=.4,label='NuFIT Gaussian 1σ (approx.)'),
             Patch(facecolor='#bcd5f1',alpha=.4,label='NuFIT Gaussian 3σ (approx.)')]
    handles +=[Line2D([0],[0],color=model_colors[m],lw=2,label=m) for m in selected]
    handles +=[Line2D([0],[0],color='black',ls=ls,lw=1.5,label=direction) for direction,ls in [('T1','-'),('T2','--'),('T3',':')]]
    handles +=[Line2D([0],[0],color='black',ls='-.',lw=1,label='Descriptive model envelope (not confidence)')]
    fig.legend(handles=handles,loc='outside lower center',ncol=4,fontsize=8)
    save(fig,figs,'fig_4_6_observable_correlations')


def c5_complex(block):
    return np.asarray(block['real'],dtype=float)+1j*np.asarray(block['imag'],dtype=float)


def c5_outputs(root,figs,tables):
    directory=root/'Reports/output/full/common_input_rge/smallY_smallL'
    names={}
    for folder in directory.iterdir():
        if folder.is_dir():
            try:names[label(folder.name)]=folder
            except ValueError:pass
    rows=[]
    for model in SIX:
        j=json_load(names[model]/'running_diagnostics.json')
        a=j['c5_threshold_contributions'];h=c5_complex(a['hard']);d=c5_complex(a['direct_running']);t=c5_complex(a['combined'])
        nh=np.linalg.norm(h);nd=np.linalg.norm(d);nt=np.linalg.norm(t)
        residual=np.linalg.norm(t-(h+d))/max(nt,1e-100)
        if residual>1e-9:raise ValueError(f'Inconsistent C5 decomposition for {model}: {residual}')
        interf=float(np.real(np.vdot(h,d))/(nh*nd))
        rows.append({'model':model,'hard_norm_GeV_inv':nh,'direct_norm_GeV_inv':nd,
                     'total_norm_GeV_inv':nt,'interference_cosine':interf,'sum_residual_relative':residual})
    write_csv(tables/'table_4_3_c5_contributions_6.csv',list(rows[0]),rows)
    fig,axs=plt.subplots(1,2,figsize=(12,4.4),constrained_layout=True)
    pos=np.arange(len(SIX));w=.23
    for offs,key,color,labeltext in [(-w,'hard_norm_GeV_inv','#5279ad','Hard'),(0,'direct_norm_GeV_inv','#39967d','Direct'),(w,'total_norm_GeV_inv','#c68a37','Combined')]:
        axs[0].bar(pos+offs,[r[key]*1e15 for r in rows],width=w,color=color,label=labeltext)
    axs[0].set_ylabel(r'$\|C_5\|_F$ [$10^{-15}\,\mathrm{GeV}^{-1}$]')
    axs[0].set_xticks(pos,SIX,rotation=35,ha='right');axs[0].legend(fontsize=8);axs[0].grid(axis='y',alpha=.15)
    axs[1].bar(pos,[r['interference_cosine'] for r in rows],color=[CLASS_COLORS[m[0]] for m in SIX])
    axs[1].set_xticks(pos,SIX,rotation=35,ha='right');axs[1].set_ylim(-1.05,0.05)
    axs[1].set_ylabel(r'$\mathrm{Re\,Tr}(C_{5,h}^{\dagger}C_{5,d})/(\|C_{5,h}\|_F\|C_{5,d}\|_F)$')
    axs[1].axhline(0,color='black',lw=.7);axs[1].grid(axis='y',alpha=.15)
    save(fig,figs,'fig_4_7_c5_contributions')


def write_latex(dest):
    """Write a ready-to-include LaTeX snippet for Figure 4.2.

    Paths are relative to the thesis main .tex file; adjust the graphicspath
    or the includegraphics path according to your thesis folder layout.
    """
    target=dest/'latex'
    target.mkdir(parents=True,exist_ok=True)
    source=r"""% Equation and Figure 4.2; requires \usepackage{graphicx}.
% If the thesis has its own subsection heading, do not repeat it here.
The relative effect of renormalisation-group evolution on an observable
$O$ is defined by
\begin{equation}
\delta_{\mathrm{RG}}(O)=100\,
\frac{O_{\mathrm{running}}-O_{\mathrm{frozen}}}
     {O_{\mathrm{frozen}}},
\label{eq:rg_correction_definition}
\end{equation}
where both calculations retain radiative threshold matching.
The frozen-parameter calculation disables parameter evolution, and
is not a tree-level neutrino-mass calculation.

\begin{figure}[htbp]
    \centering
    \includegraphics[width=0.95\textwidth]{Reports/thesis_chapter4/figures/fig_4_2_rg_corrections.pdf}
    \caption{Relative corrections from renormalisation-group evolution
    to (top) $\Delta m_{21}^{2}$ and (bottom) $\Delta m_{31}^{2}$
    for sixteen common-input T3 configurations, as defined in
    Eq.~\eqref{eq:rg_correction_definition}.
    Both running and frozen controls retain radiative threshold matching.
    A--E model families use consistent colours; the zero line
    separates enhancements and suppressions. Normal ordering is assumed,
    so the stored $\Delta m_{3\ell}^{2}$ is $\Delta m_{31}^{2}$.}
    \label{fig:rg_baseline_corrections}
\end{figure}
"""
    (target/'figure_4_2_rg_corrections.tex').write_text(source,encoding='utf-8')


def write_index(dest):
    (dest/'FIGURE_INDEX.md').write_text('''# Chapter 4 thesis exports

Generated using `Numerical/orchestration/ExportChapter4Figures.py` from saved numerical data only. PNG and PDF versions are exported.

- **4.1:** `figures/fig_4_1_representation_comparison.*`, `tables/table_4_1_baseline_16.csv`. Sixteen common-input predictions, grouped A–E.
- **4.2:** `figures/fig_4_2_rg_corrections.*`, `tables/table_4_2_rg_corrections_16.csv`, `latex/figure_4_2_rg_corrections.tex`. Signed running-versus-frozen effects at factor 1; LaTeX includes the definition equation and figure caption.
- **4.3:** `figures/fig_4_3_fermion_mass_dependence.*`; six-model scan of **both** heavy mass `MF` and matching scale `muF`.
- **4.4:** `figures/fig_4_4_fermion_rg_dependence.*`; fractional effects using the same scan.
- **4.5:** `figures/fig_4_5_coupling_dependence.*`; baseline-normalised **running** delta-m21-squared for Yukawa and lambda scans against k^4 and k^2. This plot deliberately avoids redundant atmospheric panels.
- **4.6:** `figures/fig_4_6_observable_correlations.*`; A-m2, B-m3, D-p0, common-input UV sensitivity directions T1/T2/T3. **Blue regions are approximate Gaussian NuFIT ellipses, not official correlated likelihood contours. Model dash-dot ovals are empirical 2-sigma descriptive envelopes, NOT model confidence regions.** No RGE rerun.
- **4.7:** `figures/fig_4_7_c5_contributions.*`, `tables/table_4_3_c5_contributions_6.csv`. Frobenius norms and matrix-level interference. These are not additive contributions to mass-squared differences.

Labels show $\\Delta m_{31}^{2}$ for the stored normal-ordering `delta_m3l_sq_ev2`; original field names are untouched. Class E uses its distinct lambdaT3 normalisation. Figures and tables cover *common-input*, not individually refitted, points.
''',encoding='utf-8')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path.cwd(),help='Local T3 repository root (default: working directory)')
    a=p.parse_args();root=a.root.resolve()
    base=root/'Reports/output/full'
    f=csv_rows(base/'common_input_fermion_scale/smallY_smallL/scan.csv')
    c=csv_rows(base/'common_input_coupling_scan/smallY_smallL/scan.csv')
    success(f,'fermion scan');success(c,'coupling scan')
    ferm=select_unique(f)
    require_grid(ferm,ORDER,[1.],'baseline comparison')
    cp={k:select_unique(c,k) for k in ('yukawa','lambda')}
    dest=root/'Reports/thesis_chapter4';figs=dest/'figures';tables=dest/'tables'
    figs.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.titlesize':11,'figure.dpi':125,'savefig.facecolor':'white'})
    baseline_outputs(ferm,figs,tables)
    scan_figures(ferm,cp,figs)
    correlation_figure(root,figs)
    c5_outputs(root,figs,tables)
    write_index(dest)
    write_latex(dest)
    print('Created 7 figures as PNG+PDF, 3 CSV tables, a LaTeX Figure 4.2 snippet, and FIGURE_INDEX.md under:')
    print(dest)
    print('No RGE, fitting, matching or GitHub operations were performed.')

if __name__=='__main__':main()
