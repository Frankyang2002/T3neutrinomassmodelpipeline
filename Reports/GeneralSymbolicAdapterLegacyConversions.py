"""Saved Wolfram InputForm -> TeXForm conversion utilities only."""
from pathlib import Path
import json, shutil, subprocess

def require(ok, message):
    if not ok:raise RuntimeError(message)

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

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




