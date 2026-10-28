"""Render final reduced-scalar results with the original lambda-study HTML layout.

HTML: Numerical.plotting.BuildLambda5e9Displays.page / images_as_gallery
PDF:  Reports.ReportGeneration.compile_latex_document

The lambda-study HTML renderer does not itself generate PDFs; the PDF uses the
repository's existing LaTeX compilation routine. Both consume only saved data.
No fitting, running, or Wilson matching occurs here.

    python -m scripts.BuildReducedScalarLambdaStyleReports --check
    python -m scripts.BuildReducedScalarLambdaStyleReports
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import shutil

from Numerical.plotting import BuildLambda5e9Displays as lambda_html

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'output/studies/reduced_scalar_final_reports/final_reduced_v2'
DESTINATION = ROOT / 'Reports/output/reduced_scalar_lambda_style_v1'
CASES = ('small_uv_frozen_5e-9', 'large_uv_frozen_0p5')
CASE_LABELS = ('Small topology coupling', 'Large topology coupling')
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12',
       'sin2_theta13', 'sin2_theta23')
MODELS = ('A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
          'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
          'E_+0', 'E_+2', 'E_-2', 'E_-4')
CUTOFF = 5.86
FIGURES = ('frozen_chi2_comparison.png', 'running_uv_formal_mass_endpoints.png',
           'running_uv_scalar_trajectories.png')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('Reduced-scalar lambda-style reporting: ' + message)


def read_csv(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), f'missing source CSV: {path}')
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def number(value: object, label: str, *, positive: bool = False) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'Invalid {label}: {value!r}') from exc
    require(math.isfinite(x) and (not positive or x > 0), f'nonfinite or nonpositive {label}')
    return x


def flag(value: object) -> bool:
    require(str(value).lower() in ('true', 'false'), f'invalid boolean {value!r}')
    return str(value).lower() == 'true'


def filehash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(source: Path) -> dict:
    """Validate identity, coverage and basic consistency before creating reports."""
    files = {name: source / name for name in (
        'frozen_comparison.csv', 'oscillation_observable_pulls.csv',
        'running_uv_formal_endpoints.csv', 'summary.md', 'manifest.json')}
    require(all(p.is_file() for p in files.values()), 'incomplete final_reduced_v2 source')
    manifest = json.loads(files['manifest.json'].read_text(encoding='utf-8-sig'))
    require(manifest.get('format') == 't3-reduced-final-report-v1', 'unexpected manifest format')
    require(manifest.get('primary_frozen_run') == 'reduced_frozen_overnight_v2',
            'not the completed newly fitted frozen study')
    require(manifest.get('diagnostic_running_run') == 'reduced_48_v1',
            'unexpected running diagnostic study')
    frozen = read_csv(files['frozen_comparison.csv'])
    require(len(frozen) == 32, 'expected exactly 32 frozen comparison rows')
    fits = {}
    for row in frozen:
        key = (row['case'], row['model'])
        require(key[0] in CASES and key[1] in MODELS and key not in fits,
                f'unexpected/duplicate frozen case {key}')
        ind = number(row['chi2_individual'], 'individual chi2')
        ref = number(row['chi2_shared_A_m2'], 'shared chi2')
        require(ind >= 0 and ref >= 0, 'negative chi2')
        require(flag(row['individual_accepted']) == (ind <= CUTOFF),
                'individual acceptance disagrees with chi2')
        require(flag(row['shared_A_m2_accepted']) == (ref <= CUTOFF),
                'shared acceptance disagrees with chi2')
        lam = number(row['lambdaT3_input'], 'topology coupling', positive=True)
        expected = (5e-9 if key[0] == CASES[0] else .5) * (.5 if key[1].startswith('E_') else 1.)
        require(math.isclose(lam, expected, rel_tol=1e-12, abs_tol=0),
                f'topology coupling differs for {key}')
        require(bool(row['c5_source_kind'].strip()), 'missing C5 provenance')
        fits[key] = row
    require(set(fits) == {(c,m) for c in CASES for m in MODELS}, 'frozen coverage incomplete')
    accepted = (sum(number(r['chi2_individual'],'chi2') <= CUTOFF for r in frozen),
                sum(number(r['chi2_shared_A_m2'],'chi2') <= CUTOFF for r in frozen))
    require(accepted == (32,10), f'unexpected accepted totals: {accepted}')
    require(manifest.get('number_of_individual_accepted') == 32 and
            manifest.get('number_of_shared_accepted') == 10,
            'manifest count inconsistent')
    for case in CASES:
        donor = fits[(case,'A_-2')]
        require(math.isclose(number(donor['chi2_individual'],'donor'),
                             number(donor['chi2_shared_A_m2'],'shared donor'),
                             rel_tol=1e-7,abs_tol=1e-10), 'reference donor mismatch')
    running_rows = read_csv(files['running_uv_formal_endpoints.csv'])
    require(len(running_rows) == 16, 'expected exactly 16 running-UV rows')
    running = {}
    for row in running_rows:
        m = row['model']
        require(m in MODELS and m not in running, 'unexpected/duplicate UV model')
        require(not flag(row['screen_pass']), 'historical UV status inconsistent')
        for k in ('mS1_formal_muS_tev','mS2_formal_muS_tev'):
            require(number(row[k],k,positive=True)>5, 'unexpected scalar window endpoint')
        number(row['rhoF_at_muF'], 'rhoF', positive=True)
        running[m] = row
    require(set(running)==set(MODELS), 'running-UV coverage incomplete')
    require(manifest.get('number_of_running_historical_scalar_screen_pass')==0,
            'running-UV manifest status inconsistent')
    pulls_rows = read_csv(files['oscillation_observable_pulls.csv'])
    require(len(pulls_rows)==320, 'expected 320 observable records')
    pulls = {}
    for row in pulls_rows:
        k=(row['case'],row['model'],row['yukawa_choice'],row['observable'])
        require(k[:2] in fits and k[2] in ('individual','shared_A_-2') and
                k[3] in OBS and k not in pulls, f'unexpected/duplicate pull {k}')
        prediction=number(row['prediction'],'prediction')
        target=number(row['target'],'target')
        sigma=number(row['sigma'],'sigma',positive=True)
        require(math.isclose((prediction-target)/sigma,number(row['pull'],'pull'),
                             rel_tol=1e-7,abs_tol=1e-7),f'inconsistent observable pull {k}')
        expected = fits[k[:2]]['chi2_individual' if k[2]=='individual' else 'chi2_shared_A_m2']
        require(math.isclose(number(row['chi2_five_observables'],'observable chi2'),
                             number(expected,'fit chi2'),rel_tol=1e-7,abs_tol=1e-8),
                f'chi2 differs from fitted model for {k}')
        pulls[k]=row
    require(len(pulls)==320, 'missing observable records')
    assets = {name:source/name for name in FIGURES if (source/name).is_file()}
    hashes = {p.name:filehash(p) for p in files.values()}
    hashes.update({name:filehash(p) for name,p in assets.items()})
    return {'source':source,'manifest':manifest,'fits':fits,'running':running,
            'pulls':pulls,'assets':assets,'hashes':hashes}


def esc(value: object) -> str:
    return lambda_html.esc(value)


def table(headers: list[str], rows: list[list[object]]) -> str:
    return ('<table><thead><tr>'+''.join('<th>'+esc(x)+'</th>' for x in headers)+
            '</tr></thead><tbody>\n'+''.join('<tr>'+''.join('<td>'+esc(x)+'</td>' for x in row)+'</tr>\n'
                                      for row in rows)+'</tbody></table>')


def fmt(x: object) -> str:
    return f'{number(x,"value"):.6g}'


def model_file(m: str) -> str:
    return 'model_'+m.replace('+','p').replace('-','m')+'.html'


def gallery(data: dict, dest: Path, mode: str, names: tuple[str,...]) -> str:
    paths=[data['assets'][n] for n in names if n in data['assets']]
    # Call the original lambda-study gallery copier, not a custom copy/renderer.
    return lambda_html.images_as_gallery(paths,data['source'],dest,mode)


def write_page(dest: Path, filename: str, title: str, contents: str) -> None:
    # Exact page wrapper/CSS from the original lambda=5e-9 display builder.
    (dest/filename).write_text(lambda_html.page(title,contents),encoding='utf-8')


def html_reports(data: dict, dest: Path) -> None:
    f=data['fits']; r=data['running']; p=data['pulls']
    items=[]
    for m in MODELS:
        items.append([f'<a href="{model_file(m)}">{esc(m)}</a>',
            fmt(f[(CASES[0],m)]['chi2_individual']),fmt(f[(CASES[0],m)]['chi2_shared_A_m2']),
            fmt(f[(CASES[1],m)]['chi2_individual']),fmt(f[(CASES[1],m)]['chi2_shared_A_m2'])])
    body=('<p>Reduced T3 scalar sector: 16 models, two frozen-UV scenarios. '
          'Independently fitted and transferred A(-2) Yukawa matrices are compared using '
          'the five-observable Gaussianised NuFIT target.</p>'
          '<p><strong>32/32</strong> individual fits accepted; <strong>10/32</strong> shared '
          'A(-2) evaluations accepted (chi-square &le; 5.86). '
          '<strong>0/16</strong> historical running-UV seeds satisfy the original scalar window.</p>'
          '<p><a href="benchmark_report.pdf">PDF report</a> | '
          '<a href="optimised.html">Individual fits</a> | '
          '<a href="shared_a_m2.html">Shared A(-2)</a> | '
          '<a href="running_uv.html">Exploratory UV running</a></p>'
          '<h2>Frozen-UV comparisons</h2>')
    body+=table(['Model','Low: individual chi2','Low: shared chi2',
                 'High: individual chi2','High: shared chi2'],items)
    body+='<h2>Saved figures</h2>'+gallery(data,dest,'overview',('frozen_chi2_comparison.png',))
    body+=('<h2>Interpretation</h2><p>Fitted chi-square close to zero reflects flexible '
           'Yukawa optimisation (18 real parameters for five observables), not predictivity. '
           'The four extreme-hypercharge models retain different legacy C5 matching sources. '
           'Full Matchete/RGBeta normalisation and vacuum stability are not certified.</p>'
           '<p>No new running, fitting or matching was performed. '
           'Numerical sources are preserved; see provenance.json.</p>')
    write_page(dest,'index.html','T3 reduced-scalar study | lambda-style overview',body)

    for mode,label,field in [('optimised','Individually fitted Yukawas','chi2_individual'),
                             ('shared_a_m2','Transferred A(-2) Yukawas','chi2_shared_A_m2')]:
        parts=[f'<p>{esc(label)} at each frozen-UV boundary. The A(-2) donor is fitted '
               'separately in the low and high topology-coupling scenarios.</p>']
        for case,case_label in zip(CASES,CASE_LABELS):
            parts.extend([f'<section><h2>{esc(case_label)}</h2>',
                          table(['Model','chi2','Within threshold?','C5 source'],[
                              [f'<a href="{model_file(m)}">{esc(m)}</a>',
                               fmt(f[(case,m)][field]),
                               'Yes' if number(f[(case,m)][field],'chi2')<=CUTOFF else 'No',
                               f[(case,m)]['c5_source_kind']] for m in MODELS]),'</section>'])
        parts.append('<h2>Comparison plot</h2>'+gallery(data,dest,mode,('frozen_chi2_comparison.png',)))
        write_page(dest,mode+'.html',f'T3 reduced-scalar | {label}', ''.join(parts))

    uv_body=('<p>Historical high-scale scalar and Yukawa inputs at 10^7 GeV were evolved '
             'with the candidate one-loop RGEs. The nominal 1 TeV scalar-EFT endpoint '
             'masses below are formal diagnostics, not physically matched thresholds '
             'when the scalar masses are large. No neutrino-data exclusion follows.</p>')
    uv_body+=table(['Model','Formal mS1 at 1 TeV [TeV]','Formal mS2 at 1 TeV [TeV]',
                    'rhoF at 100 TeV'],[
        [m,fmt(r[m]['mS1_formal_muS_tev']),fmt(r[m]['mS2_formal_muS_tev']),
         fmt(r[m]['rhoF_at_muF'])] for m in MODELS])
    uv_body+='<h2>Saved UV figures</h2>'+gallery(data,dest,'running',
        ('running_uv_formal_mass_endpoints.png','running_uv_scalar_trajectories.png'))
    uv_body+='<p><a href="index.html">Back to overview</a></p>'
    write_page(dest,'running_uv.html','T3 reduced-scalar | exploratory UV running',uv_body)

    for m in MODELS:
        parts=[f'<p><a href="index.html">All-model comparison</a> | '
               f'<a href="running_uv.html">UV-running diagnostics</a></p>']
        for case,case_label in zip(CASES,CASE_LABELS):
            fit=f[(case,m)]
            parts.extend([f'<h2>{esc(case_label)}</h2>',
                          table(['Yukawa choice','chi2','Accepted?'],[
                              ['Independent',fmt(fit['chi2_individual']),'Yes'],
                              ['Shared A(-2)',fmt(fit['chi2_shared_A_m2']),
                               'Yes' if flag(fit['shared_A_m2_accepted']) else 'No']]),
                          '<p>C5 source: '+esc(fit['c5_source_kind'])+
                          '; topology coupling: '+esc(fmt(fit['lambdaT3_input']))+'</p>'])
            obs_rows=[]
            for obs in OBS:
                a=p[(case,m,'individual',obs)]
                b=p[(case,m,'shared_A_-2',obs)]
                obs_rows.append([obs,fmt(a['target']),fmt(a['prediction']),fmt(a['pull']),
                                 fmt(b['prediction']),fmt(b['pull'])])
            parts.append(table(['Observable','Target','Individual','Individual pull',
                                'Shared A(-2)','Shared pull'],obs_rows))
        uv=r[m]
        parts.append('<h2>Exploratory UV-running endpoint</h2>')
        parts.append(table(['mS1 [TeV]','mS2 [TeV]','rhoF'],[[
            fmt(uv['mS1_formal_muS_tev']),fmt(uv['mS2_formal_muS_tev']),
            fmt(uv['rhoF_at_muF'])]]))
        parts.append('<p>Masses are evaluated at the nominal scalar-EFT endpoint, not '
                     'rematched physical thresholds.</p>')
        write_page(dest,model_file(m),f'T3 {m} | reduced-scalar benchmarks',''.join(parts))


def tex_escape(s: object) -> str:
    chars={'\\':r'\textbackslash{}','{':r'\{','}':r'\}', '#':r'\#',
           '$':r'\$', '%':r'\%', '&':r'\&', '_':r'\_',
           '^':r'\textasciicircum{}','~':r'\textasciitilde{}'}
    return ''.join(chars.get(c,c) for c in str(s))


def latex_report(data: dict, dest: Path) -> Path:
    f=data['fits']; r=data['running']
    tex=dest/'benchmark_report.tex'
    lines=[r'\documentclass[10pt,a4paper]{article}',
           r'\usepackage[margin=19mm]{geometry}',
           r'\usepackage{booktabs,longtable,graphicx}',
           r'\usepackage[T1]{fontenc}',
           r'\usepackage{lmodern}',
           r'\usepackage[colorlinks=true,linkcolor=blue]{hyperref}',
           r'\setlength{\parskip}{5pt}',
           r'\begin{document}',
           r'\begin{center}\LARGE\bfseries T3 reduced-scalar numerical benchmarks\end{center}',
           r'\noindent\textbf{Scope:} 16 representations, two frozen-UV scenarios; '
           r'additional scalar couplings initially zero. Historical UV-running '
           r'inputs are examined as exploratory RGE diagnostics.\\',
           r'\textbf{Acceptance:} Five-observable Gaussianised NuFIT '
           r'$\chi^2\leq 5.86$. Individual fits: 32/32. '
           r'Unchanged scenario-specific A(-2) Yukawas: 10/32.',
           r'\section*{Frozen-UV Yukawa comparisons}',
           r'\small\begin{longtable}{lrrrr}',
           r'\toprule Model & Low ind. $\chi^2$ & Low A(-2) $\chi^2$ & High ind. $\chi^2$ & High A(-2) $\chi^2$ \\',
           r'\midrule\endfirsthead',
           r'\toprule Model & Low ind. $\chi^2$ & Low A(-2) $\chi^2$ & High ind. $\chi^2$ & High A(-2) $\chi^2$ \\',
           r'\midrule\endhead']
    for m in MODELS:
        a=f[(CASES[0],m)];b=f[(CASES[1],m)]
        lines.append(' & '.join([tex_escape(m),fmt(a['chi2_individual']),
             fmt(a['chi2_shared_A_m2']),fmt(b['chi2_individual']),
             fmt(b['chi2_shared_A_m2'])])+r' \\')
    lines.extend([r'\bottomrule\end{longtable}\normalsize',
         r'\noindent Class E uses $\lambda_{T3}=2.5\times10^{-9}$ (low) and '
         r'$0.25$ (high); classes A--D use $5\times10^{-9}$ and $0.5$.',
         r'\section*{Exploratory fixed-input UV-running diagnostic}',
         r'\noindent The nominal 1 TeV scalar-EFT endpoints below are '
         r'\emph{formal diagnostics}, not rematched threshold masses. '
         r'No historical seed satisfies the original 0.5--5 TeV scalar window.',
         r'\small\begin{longtable}{lrrr}',
         r'\toprule Model & Formal $m_{S_1}$ [TeV] & Formal $m_{S_2}$ [TeV] & $\rho_F$ \\',
         r'\midrule\endfirsthead',
         r'\toprule Model & Formal $m_{S_1}$ [TeV] & Formal $m_{S_2}$ [TeV] & $\rho_F$ \\',
         r'\midrule\endhead'])
    for m in MODELS:
        u=r[m]
        lines.append(' & '.join([tex_escape(m),fmt(u['mS1_formal_muS_tev']),
            fmt(u['mS2_formal_muS_tev']),fmt(u['rhoF_at_muF'])])+r' \\')
    lines.extend([r'\bottomrule\end{longtable}\normalsize',
                  r'\section*{Interpretation and limitations}',
                  r'\begin{itemize}',
                  r'\item Independent Yukawa fits are flexible: 18 real Yukawa parameters '
                  r'are fitted to five oscillation observables, so a small $\chi^2$ '
                  r'is not an independent prediction.',
                  r'\item The zero-extra-scalar-coupling boundary is a restricted '
                  r'leading-order approximation, not the complete scalar potential.',
                  r'\item Additional B/C quartics generated by one-loop RGEs are not '
                  r'inserted into the leading fixed-order $C_5$ coefficient.',
                  r'\item Four extreme-hypercharge models use nonuniform legacy $C_5$ sources. '
                  r'Full Matchete/RGBeta normalisation and general sequential Wilson '
                  r'matching have not been certified.',
                  r'\item Heavy scalar diagnostics from running UV do not justify '
                  r'using the nominal 1 TeV scalar matching as a physical threshold.',
                  r'\end{itemize}'])
    for img,caption in [('frozen_chi2_comparison.png','Frozen-UV independent versus shared fits'),
                        ('running_uv_formal_mass_endpoints.png','Formal scalar masses at nominal 1 TeV'),
                        ('running_uv_scalar_trajectories.png','UV-running scalar trajectories, 10$^7$--10$^5$ GeV')]:
        if img in data['assets']:
            lines.extend([r'\clearpage',r'\begin{figure}[p]\centering',
                          r'\includegraphics[width=\textwidth,height=0.78\textheight,keepaspectratio]{assets/'+img+'}',
                          r'\caption{'+caption+'}',r'\end{figure}',r'\clearpage'])
    lines.append(r'\end{document}')
    tex.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return tex


def protect(source: Path, destination: Path) -> None:
    s=source.resolve(); d=destination.resolve()
    require(s != d and s not in d.parents and d not in s.parents,
            'output must not overlap numerical input directory')
    require(d != ROOT.resolve() and d != (ROOT/'Reports/output').resolve(),
            'refusing broad report output path')
    require(not destination.exists(), 'output directory already exists; choose --output')


def build(data: dict, destination: Path, *, html_only: bool = False) -> Path:
    protect(data['source'],destination)
    if not html_only and not (shutil.which('pdflatex') or shutil.which('latexmk')):
        raise RuntimeError('PDF requires pdflatex/latexmk. Install a LaTeX distribution '
                           'or run with --html-only.')
    destination.mkdir(parents=True)
    # Saved figures are used in both the lambda-gallery HTML and the PDF.
    assets=destination/'assets';assets.mkdir()
    for filename,path in data['assets'].items():
        shutil.copy2(path,assets/filename)
    tex=latex_report(data,destination)
    if not html_only:
        if not (shutil.which('pdflatex') or shutil.which('latexmk')):
            raise RuntimeError('PDF requires pdflatex/latexmk. .tex was saved; '
                               'install a LaTeX distribution or use --html-only.')
        from Reports.ReportGeneration import compile_latex_document
        compile_latex_document(tex)
        require(tex.with_suffix('.pdf').is_file(),
                f'LaTeX compilation did not create PDF: {tex} (inspect .compile.log)')
    html_reports(data,destination)
    provenance={'format':'t3-reduced-lambda-style-report-v1',
                'source':str(data['source']),
                'html_renderer':'Numerical.plotting.BuildLambda5e9Displays.page/images_as_gallery',
                'pdf_compiler':'Reports.ReportGeneration.compile_latex_document',
                'pdf_status':'skipped-by-request' if html_only else 'compiled',
                'source_sha256':data['hashes'],
                'notes':['Read-only saved-data rendering. No physics rerun.',
                         'PDF generation uses the existing TeX compiler, since the original lambda HTML renderer does not produce PDFs.']}
    (destination/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    return destination


def main(argv: list[str] | None = None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--output',type=Path,default=DESTINATION)
    parser.add_argument('--check',action='store_true',help='validate saved CSVs without rendering')
    parser.add_argument('--html-only',action='store_true',help='skip PDF compilation')
    args=parser.parse_args(argv)
    data=load(args.source.resolve())
    print('Verified: 32 frozen fits; 10 shared A(-2) accepted; 16 historical UV diagnostics.')
    if args.check:
        return 0
    output=build(data,args.output.resolve(),html_only=args.html_only)
    print('HTML:',output/'index.html')
    print('PDF:',output/'benchmark_report.pdf' if not args.html_only else 'skipped')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
