"""Saved-general-output adapter for the original Reports/ report generator.

Presentation only. Does not run Matchete/RGBeta or modify beta functions.
Output: GeneralSymbolic/UV/<model>/UV.pdf, EFT_1_after_F/<model>/..., GroupFactors/...

Wolfram is used solely to typeset the saved Matchete InputForm. If its
conversion fails, this command fails rather than substituting plain-text PDFs.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / 'output/studies/general_scalar_extended/complete_scalar_symbolic'
REPORT_ROOT = ROOT / 'Reports/output'


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def read_json(path):
    require(path.is_file(), f'Missing input: {path}')
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), f'JSON object required: {path}')
    return value


def tag(model):
    return model.replace('+', 'p').replace('-', 'm')


def esc(value):
    text = str(value)
    for a, b in [('\\',r'\textbackslash{}'),('&',r'\&'),('%',r'\%'),('$',r'\$'),('#',r'\#'),('_',r'\_'),('{',r'\{'),('}',r'\}')]:
        text = text.replace(a,b)
    return text


def preamble(title, subtitle):
    from Reports.ReportGeneration import latex_document_preamble, latex_escape_text
    lines = latex_document_preamble(
        '10pt,a4paper', 'margin=20mm',
        'amsmath,amssymb,mathtools,booktabs,longtable,array,adjustbox',
        r'\usepackage[T1]{fontenc}', r'\allowdisplaybreaks',
        r'\setlength{\parindent}{0pt}',
    )
    return lines + [r'\begin{center}',
                    r'{\Large\bfseries '+latex_escape_text(title)+r'}\\[4pt]',
                    latex_escape_text(subtitle), r'\end{center}',r'\bigskip']


def compile_report(folder, stem, contents, pdflatex):
    from Reports.ReportGeneration import compile_latex_document
    folder.mkdir(parents=True, exist_ok=True)
    tex = folder / (stem + '.tex')
    tex.write_text('\n'.join(contents + [r'\end{document}', '']), encoding='utf-8')
    compile_latex_document(tex)
    pdf = tex.with_suffix('.pdf')
    require(pdf.is_file() and pdf.stat().st_size > 1000,
            f'Original report compiler did not create {pdf}; see the LaTeX log')
    return pdf


def paper_symbol(name):
    from Reports.ReportGeneration import paper_symbol_latex
    return paper_symbol_latex(str(name))


def beta_latex(name, expression):
    # RGBeta's own TeXForm export is preferred; its expressions already use
    # the same formatter as the original restricted RGE report writer.
    from Reports.RGEComparison import normalise_rgbeta_latex
    from Reports.ReportGeneration import split_latex_terms
    expr = normalise_rgbeta_latex(expression)
    terms = split_latex_terms(expr) or [expr]
    # Each independent sum term gets its own line; no single huge TeX box.
    lhs = r'16\pi^2\,\frac{d '+paper_symbol(name)+r'}{d\ln\mu}'
    block = [r'\begin{align*}',lhs+r' &='+terms[0].strip()+r' \\']
    for term in terms[1:]:
        t=term.strip()
        block.append(r'&\quad '+(t if t.startswith(('+','-')) else '+'+t)+r' \\')
    block.append(r'\end{align*}')
    return block


# Mathematica reads its own exported InputForm file and returns separate TeXForm
# expressions for the top-level Lagrangian sum terms. No algebra is performed.
WOLFRAM_CONVERTER = r'''args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
input=args[[1]]; output=args[[2]];
text=Import[input,"Text"];
expr=ToExpression[text,InputForm];
parts=If[Head[expr]===Plus,List@@expr,{expr}];
tex=Map[ToString[TeXForm[#]]&,parts];
Export[output,<|"terms"->tex|>,"RawJSON"];
'''


def matchete_terms(input_path, output_json, wolframscript):
    if output_json.exists() and output_json.stat().st_mtime >= input_path.stat().st_mtime:
        cached = read_json(output_json)
        if isinstance(cached.get('terms'),list) and cached['terms']:
            return cached['terms']
    require(shutil.which(wolframscript) is not None,
            f'Wolfram converter required for typeset BSM Lagrangian: {wolframscript}')
    script = output_json.parent/'format_matchete_input.wl'
    script.write_text(WOLFRAM_CONVERTER,encoding='utf-8')
    p=subprocess.run([wolframscript,'-file',str(script),str(input_path),str(output_json)],
                     capture_output=True,text=True,errors='replace',check=False)
    require(p.returncode==0 and output_json.is_file(),
            f'Wolfram TeXForm conversion failed for {input_path}: {p.stderr[-1200:]} {p.stdout[-1200:]}')
    terms=read_json(output_json).get('terms')
    require(isinstance(terms,list) and terms and all(isinstance(x,str) and x for x in terms),
            f'No converted LaTeX terms: {output_json}')
    return terms


# Earlier saved RGBeta files export `report_betas` as InputForm, but do not
# necessarily contain `report_beta_latex`. Obtain TeXForm from exactly those
# saved beta expressions. This is formatting, NOT a new RGBeta calculation.
RGBETA_TEX_CONVERTER = r'''args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
input=args[[1]]; output=args[[2]];
payload=Import[input,"RawJSON"];
If[!AssociationQ[payload] || !KeyExistsQ[payload,"report_betas"],Exit[3]];
raw=payload["report_betas"];
If[!AssociationQ[raw],Exit[4]];
tex=Association@KeyValueMap[Function[{name,expression},
    name->ToString[TeXForm[ToExpression[expression,InputForm]]]],raw];
Export[output,<|"report_beta_latex"->tex|>,"RawJSON"];
'''


def obtain_beta_latex(payload, source, destination, wolframscript):
    """Use existing LaTeX if complete, otherwise typeset saved InputForm."""
    betas=payload.get('report_betas')
    require(isinstance(betas,dict) and betas, f'No RGEs in {source}')
    latex=payload.get('report_beta_latex')
    if isinstance(latex,dict) and all(isinstance(latex.get(k),str) and latex[k].strip() for k in betas):
        return latex
    cache=destination/'rgbeta_beta_latex.json'
    if cache.is_file() and cache.stat().st_mtime>=source.stat().st_mtime:
        old=read_json(cache).get('report_beta_latex')
        if isinstance(old,dict) and all(isinstance(old.get(k),str) and old[k].strip() for k in betas):
            return old
    require(shutil.which(wolframscript) is not None,
            f'Cannot format saved RGBeta InputForm without {wolframscript}: {source}')
    destination.mkdir(parents=True,exist_ok=True)
    script=destination/'format_saved_rgbeta.wl'
    script.write_text(RGBETA_TEX_CONVERTER,encoding='utf-8')
    result=subprocess.run([wolframscript,'-file',str(script),str(source),str(cache)],
                          capture_output=True,text=True,errors='replace',check=False)
    require(result.returncode==0 and cache.is_file(),
            f'RGBeta TeXForm conversion failed for {source}: '
            f'{result.stderr[-1500:]} {result.stdout[-1500:]}')
    latex=read_json(cache).get('report_beta_latex')
    require(isinstance(latex,dict) and all(isinstance(latex.get(k),str) and latex[k].strip() for k in betas),
            f'Incomplete converted beta LaTeX in {cache}')
    return latex




# Presentation-only decomposition: every expanded term becomes a coefficient
# and a symbolic structure. This extracts what RGBeta already produced; it
# does not infer or verify group factors independently.
RGBETA_STRUCTURE_CONVERTER = r'''args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
payload=Import[args[[1]],"RawJSON"];
If[!AssociationQ[payload] || !KeyExistsQ[payload,"report_betas"],Exit[3]];
raw=payload["report_betas"];
If[!AssociationQ[raw],Exit[4]];
make[expression_]:=Module[{expr,terms,rows,fac,c,structure},
  expr=ToExpression[expression,InputForm];
  terms=If[Head[Expand[expr]]===Plus,List@@Expand[expr],{Expand[expr]}];
  rows=Map[Function[term,
     fac=FactorTermsList[term];
     c=First[fac];structure=Times@@Rest[fac];
     If[!NumericQ[c] || !FreeQ[c,_Symbol],c=1;structure=term];
     <|"structure"->ToString[InputForm[structure]],
       "structure_latex"->ToString[TeXForm[structure]],
       "coefficient"->ToString[InputForm[c]],
       "coefficient_latex"->ToString[TeXForm[c]]|>
  ],terms];
  rows
];
Export[args[[2]],<|"terms"->Association@KeyValueMap[
    Function[{name,expression},name->make[expression]],raw]|>,"RawJSON"];
'''


def obtain_beta_structures(payload, source, destination, wolframscript):
    betas = payload.get('report_betas')
    require(isinstance(betas, dict) and betas, f'No beta functions in {source}')
    cache = destination / 'rgbeta_structures.json'
    if cache.is_file() and cache.stat().st_mtime >= source.stat().st_mtime:
        data = read_json(cache).get('terms')
        if isinstance(data, dict) and set(data) == set(betas):
            return data
    require(shutil.which(wolframscript) is not None,
            f'Wolfram needed to format saved beta structures: {wolframscript}')
    destination.mkdir(parents=True, exist_ok=True)
    script = destination / 'format_saved_structures.wl'
    script.write_text(RGBETA_STRUCTURE_CONVERTER, encoding='utf-8')
    proc = subprocess.run([wolframscript, '-file', str(script), str(source), str(cache)],
                          capture_output=True, text=True, errors='replace')
    require(proc.returncode == 0 and cache.is_file(),
            f'Failed to format RGE structures for {source}: '
            f'{proc.stdout[-1200:]} {proc.stderr[-1200:]}')
    data = read_json(cache).get('terms')
    require(isinstance(data, dict) and set(data) == set(betas),
            f'Missing term inventory for {source}')
    for name, terms in data.items():
        require(isinstance(terms, list) and all(
            isinstance(t, dict) and all(k in t for k in
              ('structure','structure_latex','coefficient','coefficient_latex'))
            for t in terms), f'Bad term inventory: {name} in {source}')
    return data


def numeric_display(value):
    # Render rational numbers from the native coefficient string as LaTeX.
    # No physical group-factor formulas are used here.
    from fractions import Fraction
    try:
        n = Fraction(value)
        if n.denominator == 1:
            return str(n.numerator)
        return (r'\frac{' + str(n.numerator) + '}{' + str(n.denominator) + '}')
    except (ValueError, TypeError, ZeroDivisionError):
        return None


def build_factor_rows(model_stages, stage):
    """Join native RGE term structures across models without recomputing betas."""
    names = sorted({name for data in model_stages.values() for name in data[stage]})
    for name in names:
        structures = {}
        for model, data in model_stages.items():
            for t in data[stage].get(name, []):
                key = t['structure']
                item = structures.setdefault(key, {'tex': t['structure_latex'], 'values': {}})
                if model in item['values']:
                    # Equal structures in the same beta are summed as exact
                    # rational coefficients; non-rational values are displayed
                    # individually instead of inventing a numeric prefactor.
                    old = item['values'][model]
                    from fractions import Fraction
                    try:
                        item['values'][model] = str(Fraction(old) + Fraction(t['coefficient']))
                    except (ValueError, ZeroDivisionError):
                        item['values'][model] = f'({old})+({t["coefficient"]})'
                else:
                    item['values'][model] = t['coefficient']
        yield name, structures


def factor_tex(stage, model_stages, group_rows):
    """Original GF comparison tables, populated with general RGBeta terms."""
    from types import SimpleNamespace
    from fractions import Fraction
    from collections import defaultdict
    import sympy as sp
    from Reports.GroupFactorReports import _saved_rge_term_table, _run_key
    from Reports.ReportGeneration import paper_symbol_latex, paper_notation_key_lines
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import REP

    records = []
    for model in model_stages:
        cls, signed = model.split('_')
        d1, d2, df = REP[cls]
        records.append(SimpleNamespace(name=cls, alpha=int(signed), shared_scalar=False,
                      d_s1=d1, d_s2=d2, d_f=df))
    lines = preamble('T3 general group factors: '+stage,
                     'Original RGE-structure / prefactor comparison layout')
    lines += paper_notation_key_lines()
    lines += [r'\section*{Field representations}',
              r'\begin{longtable}{@{}lccc@{}}',r'\toprule',
              r'Model & $S_1:(d,Y)$ & $S_2:(d,Y)$ & $F:(d,Y)$ \\',r'\midrule']
    for model in model_stages:
        fld={g['field']:g for g in group_rows if g['model']==model}
        vals = [f"$({fld[k]['d_SU2']},{sp.latex(sp.Rational(fld[k]['Y']))})$"
                for k in ('S1','S2','F')]
        lines += [esc(model)+' & '+' & '.join(vals)+r' \\']
    lines += [r'\bottomrule',r'\end{longtable}',
              r'\section*{One-loop renormalisable RGE structure prefactors}']
    for name, structures in build_factor_rows(model_stages,stage):
        per_run = defaultdict(dict)
        for structure,item in structures.items():
            # Match original structure/prefactor table contract.
            tex = item['tex']
            for rec in records:
                model = rec.name+'_'+('+' if rec.alpha>=0 else '')+str(rec.alpha)
                value=item['values'].get(model)
                if value is None:
                    continue
                try: per_run[_run_key(rec)][tex] = sp.sympify(value)
                except (sp.SympifyError,TypeError,ValueError):
                    per_run[_run_key(rec)][tex] = sp.Symbol('unparsed')
        if per_run:
            lines += _saved_rge_term_table(
                r'\beta_{'+paper_symbol_latex(name)+r'}^{(1)}',records,per_run)
    return lines


def lagrangian_latex(terms):
    out=[r'\section*{Complete BSM UV Lagrangian}',
         r'\emph{Matchete BSM Lagrangian, excluding the pure Standard Model contribution.}',
         r'\medskip']
    for idx,term in enumerate(terms):
        # The term may be large. adjustbox prevents hard page-overflow while
        # independent terms retain real TeX maths rather than code listings.
        out += [r'\noindent\begin{adjustbox}{max width=\linewidth}',
                r'$\displaystyle '+(r'\mathcal{L}_{\rm BSM} = ' if idx==0 else ('+' if not term.lstrip().startswith('-') else ''))+term+r'$',
                r'\end{adjustbox}\par\smallskip']
    return out


def render(*,study='full16_general_symbolic',run_id='full16_general_v1',
           pdflatex='pdflatex',wolframscript='wolframscript'):
    for part in (study,run_id):
        require(bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',part)),f'Bad name: {part}')
    out=REPORT_ROOT/study/'GeneralSymbolic'
    data=DATA_ROOT/run_id
    summary=read_json(out/'summary.json')
    require(summary.get('all_requested_models_validated') is True,'Saved general study incomplete')
    with (out/'summary.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    require(rows and len(rows)==int(summary['model_count_validated']),'Incomplete model inventory')
    with (out/'group_invariants.csv').open(encoding='utf-8-sig',newline='') as f:
        groups=list(csv.DictReader(f))
    require(len(groups)==3*len(rows),'Incomplete group invariant input')
    completed=[]
    model_stages={}
    for row in rows:
        model=row['model']; mt=tag(model); saved=data/mt
        lag=read_json(saved/'lagrangian.json')
        require(lag.get('general_scalar_quadratics_requested') is True and
                lag.get('extra_yukawas_requested') is True, f'{model}: not a general Lagrangian')
        native=Path(lag['bsm_lagrangian_inputform_file'])
        require(native.is_file(),f'Missing BSM InputForm: {native}')
        model_heading=r'\section*{'+esc(model)+r'}'
        intro=[r'Additional quartics: '+esc(row['extra_quartics'])+', additional bilinears: '+esc(row['extra_bilinears'])+
               ', extra Yukawa: '+esc(row['extra_yukawa'] or 'none')+r'.\\',
               r'\emph{General matching and Wilson coefficients are not computed in this symbolic report.}']
        model_stages[model]={}
        ux=out/'UV'/mt
        ux.mkdir(parents=True,exist_ok=True)
        terms=matchete_terms(native,ux/'matchete_bsm_latex.json',wolframscript)
        require(all('\\mathcal{L}_{\\rm BSM}' not in t for t in terms),'Unexpected nested LaTeX')
        bsm=lagrangian_latex(terms)
        from Reports.ReportGeneration import paper_notation_key_lines
        bsm = paper_notation_key_lines() + bsm
        bsm += [r'\section*{Interaction inventory}',r'\begin{itemize}']
        bsm += [r'\item \texttt{'+esc(x)+'}' for x in lag['interactions']]
        bsm += [r'\end{itemize}']
        for stage,stagefolder,stem in [('UV','UV','UV'),('EFT1','EFT_1_after_F','EFT_1_after_F')]:
            payload=read_json(saved/(stage.lower()+'_rgbeta.json'))
            betas=payload.get('report_betas')
            require(isinstance(betas,dict) and betas,f'No {stage} RGEs: {model}')
            dest=out/stagefolder/mt
            latex=obtain_beta_latex(payload,saved/(stage.lower()+'_rgbeta.json'),dest,wolframscript)
            model_stages[model][stage]=obtain_beta_structures(payload,saved/(stage.lower()+'_rgbeta.json'),dest,wolframscript)
            contents=preamble(f'T3 {model}: '+('full general UV Lagrangian and RGEs' if stage=='UV' else 'general EFT1 RGEs'),
                              'Matchete/RGBeta symbolic output in mathematical notation')
            contents += [model_heading]+intro
            if stage=='UV':
                contents += bsm
            contents += [r'\clearpage',r'\section*{One-loop renormalisable beta functions}',
                         r'$16\pi^2\,dX/d\ln\mu=\beta_X$']
            for name in betas:
                contents += [r'\subsection*{$'+paper_symbol(name)+r'$}']
                contents += beta_latex(name,str(latex[name]))
            dest=out/stagefolder/mt
            pdf=compile_report(dest,stem,contents,pdflatex)
            completed.append(pdf)
            print('PDF:',pdf,flush=True)
    for stage, stagefolder, stem in [('UV','GroupFactors','GF_UV'),
                                      ('EFT1','GroupFactors','GF_EFT1')]:
        report = factor_tex(stage, model_stages, groups)
        pdf = compile_report(out / stagefolder / ('UV' if stage == 'UV' else 'EFT1'),
                             stem, report, pdflatex)
        completed.append(pdf)
        print('PDF:',pdf,flush=True)
    index=['# General symbolic LaTeX/PDF reports','',
           'UV and EFT1 PDF reports are separated by model; GF_UV and GF_EFT1 compare all models.','',
           '| Model | UV | EFT1 |','|---|---|---|']
    for r in rows:
        m=tag(r['model'])
        index += [f'| {r["model"]} | [UV](UV/{m}/UV.pdf) | [EFT1](EFT_1_after_F/{m}/EFT_1_after_F.pdf) |']
    index += ['', '[UV group-factor comparison](GroupFactors/UV/GF_UV.pdf)',
              '[EFT1 group-factor comparison](GroupFactors/EFT1/GF_EFT1.pdf)', '']
    (out/'PDF_REPORTS.md').write_text('\n'.join(index)+'\n',encoding='utf-8')
    return completed


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',default='full16_general_symbolic')
    p.add_argument('--run-id',default='full16_general_v1')
    p.add_argument('--pdflatex',default='pdflatex')
    p.add_argument('--wolframscript',default='wolframscript')
    a=p.parse_args(argv)
    try:
        render(study=a.study,run_id=a.run_id,pdflatex=a.pdflatex,wolframscript=a.wolframscript)
        return 0
    except (OSError,KeyError,ValueError,RuntimeError,json.JSONDecodeError) as e:
        print('General T3 LaTeX PDF generation failed:',e,file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
