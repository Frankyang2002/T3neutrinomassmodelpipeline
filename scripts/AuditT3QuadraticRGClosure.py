"""Batch 4: zero-boundary RG closure for the eight extra T3 quadratic terms.

Read-only analysis of *saved* UV/EFT1 RGBeta exports from the current general
T3 experiment.  Never invokes Wolfram or changes the scalar model.

It substitutes mu2=Bar[mu2]=0 exactly in the one-loop beta expression,
using a deliberately small arithmetic AST and failing on unsupported forms.
It also evaluates the resulting source at the SAVED EFT1 initial boundary and
makes a fixed-slope leading-log estimate.  This is NOT a new zero-mu2 ODE
trajectory, a physical Matchete/RGBeta mapping or a C5 error estimate.

    python -m scripts.AuditT3QuadraticRGClosure
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import operator
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'output' / 'general_t3_quadratic_eft1'
OUT = ROOT / 'output' / 'scalar_relevance'
LOOP = 16 * math.pi**2
MODELS = ('A_+0', 'A_-2', 'A_-4', 'B_+1', 'B_-1', 'B_-3',
          'C_+1', 'C_-1', 'C_-3', 'D_+0', 'D_+2', 'D_-2',
          'E_+0', 'E_+2', 'E_-2', 'E_-4')
MASS_SYMBOLS = {
    'self_s1': 'muT3GeneralS1Sq',
    'self_s2': 'muT3GeneralS2Sq',
    'mixed': 'muT3GeneralS1S2',
    'none': None,
}
EXPECTED = {'A_+0':'self_s1','A_-2':'self_s2','B_-1':'mixed',
            'C_-1':'mixed','D_+0':'self_s1','D_-2':'self_s2',
            'E_+0':'self_s1','E_-2':'self_s2'}
STATUS = 'CandidateOneLoopQuadraticZeroBoundarySourceAudit'


def die(message: str) -> None:
    raise ValueError('Batch 4 fail-closed: ' + message)


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        die(f'missing required input: {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(obj, dict):
        die(f'not a JSON object: {path}')
    return obj


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tag(model: str) -> str:
    name,alpha = model.split('_')
    return f'{name}_{"m"+alpha[1:] if alpha.startswith("-") else "p"+alpha[1:]}'


def replace_wolfram_functions(expression: str) -> str:
    """Reduce the *restricted* beta-mass expression to Python arithmetic.

    Tr[...] can appear multiplying mu2 in UV beta functions.  A distinct
    token records each trace; the mu2=0 result MUST NOT contain a trace.
    Bar[symbol] becomes BAR_symbol. No general Wolfram syntax is accepted.
    """
    expr = expression
    counter = 0
    while True:
        m = re.search(r'(?<![A-Za-z0-9_])Tr\[', expr)
        if m is None:
            break
        begin = m.start()
        k = m.end()
        depth = 1
        while k < len(expr) and depth:
            depth += (expr[k] == '[') - (expr[k] == ']')
            k += 1
        if depth:
            die('unbalanced Tr[...] in beta expression')
        expr = expr[:begin] + f'__TR_{counter}' + expr[k:]
        counter += 1
    expr = re.sub(r'\bBar\[([A-Za-z][A-Za-z0-9_]*)\]',r'BAR_\1',expr)
    if any(sym in expr for sym in ('[',']','.')):
        die('unsupported Wolfram construct in quadratic mass beta: '+expr)
    return expr.replace('^','**')


def parse_beta(expression: str) -> ast.expr:
    if not isinstance(expression,str) or not expression.strip():
        die('missing quadratic mass beta')
    try:
        tree = ast.parse(replace_wolfram_functions(expression),mode='eval').body
    except SyntaxError as exc:
        die('unparseable mass beta: '+str(exc))
    validate_ast(tree)
    return tree


def validate_ast(node: ast.AST) -> None:
    if isinstance(node, ast.Constant):
        if not isinstance(node.value,(float,int,complex)):
            die('non-numeric AST constant')
    elif isinstance(node, ast.Name):
        if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*',node.id):
            die('unexpected variable name')
    elif isinstance(node, ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):
        validate_ast(node.operand)
    elif isinstance(node, ast.BinOp) and isinstance(node.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):
        validate_ast(node.left)
        validate_ast(node.right)
    else:
        die('unsupported AST element '+type(node).__name__)


def C(v: int | float) -> ast.Constant:
    return ast.Constant(value=v)


def iszero(a: ast.expr) -> bool:
    return isinstance(a,ast.Constant) and a.value == 0


def substitute(node: ast.expr, zeros: set[str]) -> ast.expr:
    """Exact syntactic substitution + annihilation for arithmetic beta ASTs."""
    if isinstance(node,ast.Name):
        return C(0) if node.id in zeros else ast.Name(id=node.id,ctx=ast.Load())
    if isinstance(node,ast.Constant):
        return C(node.value)
    if isinstance(node,ast.UnaryOp):
        term = substitute(node.operand,zeros)
        if iszero(term):
            return C(0)
        return ast.UnaryOp(op=node.op,operand=term)
    if isinstance(node,ast.BinOp):
        lhs,rhs = substitute(node.left,zeros),substitute(node.right,zeros)
        if isinstance(node.op,ast.Mult) and (iszero(lhs) or iszero(rhs)):
            return C(0)
        if isinstance(node.op,ast.Div):
            if names(node.right) & zeros:
                die('quadratic mass coefficient appears in denominator')
            if iszero(lhs):
                return C(0)
        if isinstance(node.op,ast.Pow):
            if names(node.right) & zeros:
                die('quadratic mass coefficient appears in exponent')
            if iszero(lhs):
                if not isinstance(rhs,ast.Constant) or rhs.value <= 0:
                    die('unsupported zero-base power')
                return C(0)
        if isinstance(node.op,(ast.Add,ast.Sub)):
            if iszero(lhs) and iszero(rhs):
                return C(0)
            if iszero(rhs):
                return lhs
            if iszero(lhs):
                return ast.UnaryOp(op=ast.USub(),operand=rhs) if isinstance(node.op,ast.Sub) else rhs
        if isinstance(lhs,ast.Constant) and isinstance(rhs,ast.Constant):
            ops = {ast.Add: operator.add,ast.Sub:operator.sub,ast.Mult:operator.mul,
                   ast.Div:operator.truediv,ast.Pow:operator.pow}
            try:
                value=ops[type(node.op)](lhs.value,rhs.value)
            except (ValueError,OverflowError,ZeroDivisionError) as exc:
                die('invalid arithmetic simplification: '+str(exc))
            if not isinstance(value,(float,int,complex)) or not math.isfinite(abs(value)):
                die('non-finite arithmetic simplification')
            return C(value)
        return ast.BinOp(left=lhs,op=node.op,right=rhs)
    die('invalid node substitution')


def names(expr: ast.expr) -> set[str]:
    return {n.id for n in ast.walk(expr) if isinstance(n,ast.Name)}


def exprtext(expr: ast.expr) -> str:
    return ast.unparse(ast.fix_missing_locations(expr))


def evaluate(expr: ast.expr, values: dict[str,complex | float]) -> complex:
    syms=names(expr)
    if not syms.issubset(values):
        die('missing scalar boundary values for '+str(sorted(syms-set(values))))
    def go(n:ast.expr):
        if isinstance(n,ast.Constant): return n.value
        if isinstance(n,ast.Name): return values[n.id]
        if isinstance(n,ast.UnaryOp):
            x=go(n.operand)
            return x if isinstance(n.op,ast.UAdd) else -x
        if isinstance(n,ast.BinOp):
            lhs,rhs=go(n.left),go(n.right)
            funcs={ast.Add:operator.add,ast.Sub:operator.sub,ast.Mult:operator.mul,
                   ast.Div:operator.truediv,ast.Pow:operator.pow}
            return funcs[type(n.op)](lhs,rhs)
        die('unknown evaluation node')
    try: result=complex(go(expr))
    except (ValueError,OverflowError,ZeroDivisionError) as e:
        die('evaluation failure: '+str(e))
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        die('nonfinite numerical source')
    return result


def complex_value(raw: Any) -> complex:
    if isinstance(raw,dict) and set(raw)=={'real','imag'}:
        return complex(float(raw['real']),float(raw['imag']))
    if isinstance(raw,(float,int)):
        return complex(raw)
    die('invalid scalar coupling encoding '+repr(raw))


def source_at_zero(expression: str, symbol: str, allowed_quartics: set[str]) -> tuple[ast.expr, bool, list[str]]:
    tree=parse_beta(expression)
    zeros={symbol,'BAR_'+symbol}
    source=substitute(tree,zeros)
    remaining=names(source)
    if any(k.startswith('__TR_') for k in remaining):
        die('trace remained in zero-mu2 beta; cannot certify simple source')
    if any(k.startswith('BAR_') for k in remaining):
        die('barred symbol remained; requires extended source parser')
    unexpected=remaining-allowed_quartics-{'mS1Sq','mS2Sq'}
    if unexpected:
        die('nonquartic source structure: '+str(sorted(unexpected)))
    noextra = substitute(source,allowed_quartics)
    if not iszero(noextra):
        die('mu2 beta does not close with all additional quartics set zero')
    return source,iszero(source),sorted(remaining & allowed_quartics)


def check_live(model: str, stage: str, raw: dict[str,Any], kind: str) -> tuple[str | None, ast.expr | None, list[str]]:
    if raw.get('format')!='general-t3-quadratic-rgbeta-experimental-v1' or raw.get('status')!='SuccessCandidateMassBeta':
        die(f'{model}/{stage}: missing live candidate beta certificate')
    meta=raw.get('metadata',{})
    expected={'GeneralScalarModel':model,'GeneralScalarStage':stage,
              'QuadraticKind':kind,'QuadraticMassSymbol':MASS_SYMBOLS[kind],
              'PhysicalMatcheteRGBetaMassBasisCertified':False,
              'QuadraticPrefactorConvention':'V2 = mu2 * invariant + h.c.; no 1/2',
              'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'}
    for k,v in expected.items():
        if meta.get(k)!=v:
            die(f'{model}/{stage}: mismatched {k}: {meta.get(k)!r} != {v!r}')
    betas=raw.get('report_betas',{})
    if not isinstance(betas,dict) or not betas:
        die(f'{model}/{stage}: missing full beta inventory')
    sym=MASS_SYMBOLS[kind]
    if kind=='none':
        if any(n.startswith('muT3General') for n in betas):
            die(model+' unexpected quadratic beta')
        return None,None,[]
    beta=betas.get(sym)
    if not isinstance(beta,str):
        die(f'{model}/{stage}: missing {sym} beta')
    if 'BetaTerm[' in beta or '$Failed' in beta:
        die(f'{model}/{stage}: unevaluated beta expression')
    registry=set(meta.get('GeneralScalarNewBetaKeys',[]))
    # Existing (reused) RGBeta special quartics for B(-1)/C(-1).
    reused={'lambdaS1barS2S2bar2','lambdaS1bar2S2bar2','lambdaS1S1bar2S2bar',
            'lambdaHHdagS2S2','lambdaHHdagS1barS1bar','lambdaHHdagS1barS2barCross'} if model in ('B_-1','C_-1') else set()
    src,trivial,needed=source_at_zero(beta,sym,registry|reused)
    if trivial:
        die(f'{model}/{stage}: expected nontrivial source from extra quartics; manual check required')
    return beta,src,needed


def boundary(model: str, data: dict[str,Any], kind: str) -> tuple[dict[str,complex],float,float]:
    if data.get('model')!=model or data.get('format')!='general-t3-candidate-quadratic-eft1-v1' or data.get('quadratic_kind')!=kind:
        die(model+': inconsistent saved EFT1 running report')
    if data.get('status')!='CandidateQuadraticAndFullQuarticEFT1Integrated' or data.get('physical_matchete_rgbeta_mapping_certified') is not False or data.get('physical_neutrino_mass_prediction') is not False:
        die(model+': saved running report makes unexpected claims')
    if kind!='none' and data.get('mu2_beta_included') is not True:
        die(model+': candidate quadratic was not evolved')
    first,last=data.get('mu_initial_gev'),data.get('mu_final_gev')
    if not all(isinstance(x,(int,float)) and math.isfinite(x) and x>0 for x in (first,last)) or last>=first:
        die(model+': invalid EFT1 scale interval')
    info=data.get('initial_baseline',{})
    m1,m2=info.get('mS1Sq'),info.get('mS2Sq')
    if not all(isinstance(x,(float,int)) and math.isfinite(x) and x>0 for x in (m1,m2)):
        die(model+': invalid scalar diagonal masses')
    values={'mS1Sq':complex(m1),'mS2Sq':complex(m2)}
    for key,value in info.items():
        if key.startswith('lambda') and value is not None:
            values[key]=complex_value(value)
    for key,value in data.get('initial_new_quartics',{}).items():
        if key in values:
            die(model+': duplicate extra-quartic boundary symbol')
        values[key]=complex_value(value)
    return values,float(first),float(last)


def model_audit(model: str, directory: Path) -> tuple[dict[str,Any],dict[str,str]]:
    kind=EXPECTED.get(model,'none')
    stem=tag(model)
    files={st:directory/f'{stem}_{st.lower()}_quadratic_live_rgbeta.json' for st in ('UV','EFT1')}
    files['running']=directory/f'{stem}_eft1_quadratic_running.json'
    payloads={key:read_json(path) for key,path in files.items()}
    rows={}
    for stage in ('UV','EFT1'):
        beta,source,drivers=check_live(model,stage,payloads[stage],kind)
        rows[stage]={'beta_inputform':beta,'zero_mu2_source_expression':exprtext(source) if source is not None else None,
                     'source_quartic_symbols':drivers,'one_loop_mu2_zero_closed_with_all_extra_quartics_zero':True}
    if kind!='none':
        if rows['UV']['zero_mu2_source_expression']!=rows['EFT1']['zero_mu2_source_expression']:
            die(model+': UV/EFT1 zero-mu2 source differs; scope must be reconsidered')
    values,scale1,scale2=boundary(model,payloads['running'],kind)
    slope=None
    source_at_boundary=None
    ratio=None
    if kind!='none':
        src,_trivial,_terms=source_at_zero(rows['EFT1']['beta_inputform'], MASS_SYMBOLS[kind],
                       set(rows['EFT1']['source_quartic_symbols']))
        value=evaluate(src,values)
        shift=value*math.log(scale2/scale1)/LOOP
        baseline=min(values['mS1Sq'].real,values['mS2Sq'].real)
        source_at_boundary={'real':value.real,'imag':value.imag,'unit':'GeV^2 (16pi^2 beta coefficient)'}
        slope={'real':shift.real,'imag':shift.imag,'abs_gev2':abs(shift),
               'method':'frozen-boundary first-log slope; not an integrated trajectory'}
        ratio=abs(shift)/baseline
    result={
        'model':model,'kind':kind,'extra_quadratic_exists':kind!='none',
        'uv':rows['UV'],'eft1':rows['EFT1'],
        'mu_initial_gev':scale1,'mu_final_gev':scale2,
        'saved_nonzero_mu2_run_not_a_zero_boundary_comparison':kind!='none',
        'eft1_zero_mu2_source_at_saved_nonzero_quartic_boundary':source_at_boundary,
        'eft1_fixed_slope_estimated_generated_mu2_gev2':slope,
        'ratio_abs_generated_mu2_over_min_input_m2':ratio,
        'zero_bilinear_alone_generally_rg_closed':False if kind!='none' else None,
        'all_additional_quartics_and_bilinear_zero_one_loop_mass_source':True,
        'physical_mass_basis_and_c5_not_certified':True,
    }
    return result,{str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path):sha(path)
                   for path in files.values()}


def build(directory: Path=DATA) -> dict[str,Any]:
    results={}
    hashes={}
    for model in MODELS:
        row,h=model_audit(model,directory)
        results[model]=row
        hashes.update(h)
    active=[r for r in results.values() if r['extra_quadratic_exists']]
    if len(active)!=8 or sum(r['kind']=='mixed' for r in active)!=2:
        die('expected 8 quadratic slots (2 mixed) across 16 models')
    return {
        'status':STATUS,'basis':'candidate native RGBeta, physical coupling map uncertified',
        'model_count':16,'quadratic_slots':8,'models_with_nonzero_saved_boundary_source':
            [r['model'] for r in active if abs(complex(
                r['eft1_zero_mu2_source_at_saved_nonzero_quartic_boundary']['real'],
                r['eft1_zero_mu2_source_at_saved_nonzero_quartic_boundary']['imag'])) > 1e-8],
        'zero_quadratic_not_closed_at_generic_nonzero_additional_quartics':True,
        'source_vanishes_if_all_additional_quartics_zero_at_one_loop':True,
        'numerical_result_only_a_fixed_boundary_slope':True,
        'not_a_zero_mu2_trajectory':True,'not_a_neutrino_mass_error_estimate':True,
        'uv_and_eft1_beta_exports_checked':True,
        'source_sha256':hashes,'models':results,
    }


def write(report:dict[str,Any], directory:Path=OUT) -> None:
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'quadratic_rg_zero_boundary_audit.json').write_text(
        json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    headers=['model','kind','source_quartics','eft1_beta_zero_mu2_real_gev2',
             'eft1_beta_zero_mu2_imag_gev2','estimated_generated_mu2_abs_gev2',
             'estimated_generated_mu2_over_min_m2']
    with (directory/'quadratic_rg_zero_boundary_audit.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=headers)
        writer.writeheader()
        for model,row in report['models'].items():
            v=row['eft1_zero_mu2_source_at_saved_nonzero_quartic_boundary'] or {}
            writer.writerow({'model':model,'kind':row['kind'],
                             'source_quartics':';'.join(row['eft1']['source_quartic_symbols']),
                             'eft1_beta_zero_mu2_real_gev2':v.get('real'),
                             'eft1_beta_zero_mu2_imag_gev2':v.get('imag'),
                             'estimated_generated_mu2_abs_gev2':(row['eft1_fixed_slope_estimated_generated_mu2_gev2'] or {}).get('abs_gev2'),
                             'estimated_generated_mu2_over_min_m2':row['ratio_abs_generated_mu2_over_min_input_m2']})
    lines=['# T3 scalar relevance — Batch 4: quadratic RG source audit','',
           '- 16 model reports checked; **8** additional bilinears (6 self, 2 mixed).',
           '- Saved UV and EFT1 **candidate RGBeta beta functions** were inspected.',
           '- The eight extra quadratic beta functions have generically nonzero sources at mu2=0 when their source quartics are nonzero.',
           '- When **all additional quartics** are also set zero, each beta-mass source vanishes algebraically at this one-loop order.',
           '- Fixed-slope figures below use the saved EFT1 upper-scale boundary, with mu2 manually set zero; no independent zero-mu2 trajectory was integrated.',
           '',
           '| Model | Kind | Source quartics | Estimated |mu2|/min(m²) from one EFT1 log |',
           '|---|---|---:|---:|']
    for model,row in report['models'].items():
        v=row['ratio_abs_generated_mu2_over_min_input_m2']
        lines.append(f'| {model} | {row["kind"]} | {len(row["eft1"]["source_quartic_symbols"])} | {"—" if v is None else f"{v:.5%}"} |')
    lines += ['', '## Perturbative interpretation','',
      'For mu2(mu_F)=0 and nonzero source quartics at mu_F,',
      '', '    mu2(mu_S) ≈ beta_mu2^(1)(mu_F;mu2=0) ln(mu_S/mu_F)/(16*pi²).',
      '', 'This is an **illustrative fixed-slope leading-log term**; it neglects coupling running, feedback, and the UV segment. It is not the difference between two complete physical parameter trajectories.',
      '', 'If extra quartics are themselves set zero at the high-scale boundary and generated only at one loop, induced mu2 may first appear at higher iterated order. A full proof of RG closure requires considering the complete coupled beta system, not only beta_mu2.',
      '', '## Scope limits','',
      '- Matchete–RGBeta coupling and scalar-field normalizations are NOT certified.',
      '- The existing quadratic RG trajectories started from a nonzero mu2=0.01*min(m²), not from zero.',
      '- No physical one-loop C5 with nonzero quadratic mixing is computed here.',
      '- A small induced mass parameter is not a validated relative neutrino-mass error.',
      '- Source quartics can affect scalar RG evolution even if absent in the direct one-vertex C5 graph.',
      '']
    (directory/'quadratic_rg_zero_boundary_summary.md').write_text('\n'.join(lines),encoding='utf-8')


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,default=DATA)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    opts=parser.parse_args()
    input_dir=opts.input_dir if opts.input_dir.is_absolute() else ROOT/opts.input_dir
    output_dir=opts.output_dir if opts.output_dir.is_absolute() else ROOT/opts.output_dir
    report=build(input_dir)
    write(report,output_dir)
    print(f'Batch 4 one-loop zero-boundary audit: {report["model_count"]} models, '
          f'{report["quadratic_slots"]} quadratic slots; physical C5 NOT certified')
    print(f'Source files: {len(report["source_sha256"])} hash-verified')
    print('Output:',output_dir/'quadratic_rg_zero_boundary_summary.md')

if __name__=='__main__':
    main()
