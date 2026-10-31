"""Exact formal monomial audit of saved Matchete pole and mass-log coefficients.

The Matchete InputForm expressions are parsed as *commutative products of opaque
non-polynomial factors*. A complete NCM[...] is one indivisible factor; no
SU(2), Majorana or flavour identities are assumed. This proves an equality in
the free commutative algebra of these opaque factors, not in a physical EFT.

No eval/sympify or Mathematica execution; unsupported syntax fails closed.
"""
from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
import re

MASS_LOGS = {
    'Y': 'Log[μbar2/Coupling[MSY, {}, 0]^2]',
    'X': 'Log[μbar2/Coupling[MSX, {}, 0]^2]',
}
CALL = re.compile(r'[A-Za-z_μϵ$][A-Za-z_0-9μϵ$]*\[')
NAME = re.compile(r'[A-Za-z_μϵ$][A-Za-z_0-9μϵ$]*')
TOKENIZED = re.compile(r'[A-Za-z_0-9 +*/()^\-\s]*')
LAMBDA = re.compile(r'(?<![A-Za-z_0-9])Coupling\[(lambda[A-Za-z_0-9]+),')


class ConsistencyError(ValueError):
    pass


def terms(source: str) -> tuple[str, ...]:
    """Only top-level signs split a sum; nested Mathematica calls remain opaque."""
    source = source.strip()
    if not source:
        raise ConsistencyError('Empty InputForm')
    if source == '0':
        return ()
    pairs = {')': '(', ']': '[', '}': '{'}
    stack = []
    start = 0
    found = []
    for i, ch in enumerate(source):
        if ch in '([{':
            stack.append(ch)
        elif ch in ')]}':
            if not stack or stack.pop() != pairs[ch]:
                raise ConsistencyError('Unbalanced input')
        elif ch in '+-' and not stack and i > start:
            found.append(source[start:i].strip())
            start = i
    if stack:
        raise ConsistencyError('Unbalanced input')
    found.append(source[start:].strip())
    if any(not s or (i and s[0] not in '+-') for i,s in enumerate(found)):
        raise ConsistencyError('Unrecognized additive expression')
    return tuple(found)


def _atomize(expr: str, table: dict[str,str]) -> str:
    """Replace each entire Foo[...], including nested brackets, by one atom."""
    pieces=[]
    i=0
    while i<len(expr):
        m=CALL.match(expr, i)
        if not m:
            pieces.append(expr[i]);i+=1;continue
        start=i; j=m.end(); depth=1
        while j<len(expr) and depth:
            if expr[j]=='[':
                depth+=1
            elif expr[j]==']':
                depth-=1
            j+=1
        if depth:
            raise ConsistencyError('Unclosed function call')
        opaque=expr[start:j]
        if opaque not in table:
            table[opaque]='atom_'+str(len(table))
        pieces.append(table[opaque])
        i=j
    return ''.join(pieces)


def _merge(a: dict[str,int], b: dict[str,int], sign=1) -> dict[str,int]:
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,0)+sign*v
        if out[k]==0:
            del out[k]
    return out


def _mul(left, right):
    a,x=left; b,y=right
    return (a*b, _merge(x,y))


def _evaluate(node):
    """Exact rational Laurent monomial; reject sums inside a single term."""
    if isinstance(node,ast.Constant) and type(node.value) is int:
        return Fraction(node.value), {}
    if isinstance(node,ast.Name):
        return Fraction(1), {node.id:1}
    if isinstance(node,ast.UnaryOp) and isinstance(node.op, (ast.USub,ast.UAdd)):
        a,x=_evaluate(node.operand)
        return ((-a if isinstance(node.op,ast.USub) else a),x)
    if isinstance(node,ast.BinOp):
        if isinstance(node.op,ast.Mult):
            return _mul(_evaluate(node.left),_evaluate(node.right))
        if isinstance(node.op,ast.Div):
            a,x=_evaluate(node.left); b,y=_evaluate(node.right)
            if b==0:
                raise ConsistencyError('Division by zero')
            return a/b, _merge(x,y,-1)
        if isinstance(node.op,ast.Pow):
            a,x=_evaluate(node.left)
            if isinstance(node.right,ast.Constant) and type(node.right.value) is int:
                exponent=node.right.value
            elif (isinstance(node.right,ast.UnaryOp) and isinstance(node.right.op,ast.USub) and
                  isinstance(node.right.operand,ast.Constant) and type(node.right.operand.value) is int):
                exponent=-node.right.operand.value
            else:
                raise ConsistencyError('Non-integer exponent')
            if abs(exponent)>12 or (a==0 and exponent<0):
                raise ConsistencyError('Unsupported power')
            return a**exponent,{k:v*exponent for k,v in x.items() if v*exponent}
        # Some Matchete monomials have parenthesized integer sums.
        if isinstance(node.op,(ast.Add,ast.Sub)):
            a,x=_evaluate(node.left);b,y=_evaluate(node.right)
            if x!=y or x:
                raise ConsistencyError('Non-monomial additive factor')
            return (a+b if isinstance(node.op,ast.Add) else a-b),{}
    raise ConsistencyError('Unsupported InputForm arithmetic node: '+type(node).__name__)


def polynomial(sources: tuple[str,...], *, table: dict[str,str],
               remove_log: str|None=None) -> dict[tuple[tuple[str,int],...],Fraction]:
    """Exact signed coefficient map for a sum of opaque-factor monomials."""
    result: dict[tuple[tuple[str,int],...],Fraction]={}
    for raw in sources:
        expr=raw
        if remove_log is not None:
            if expr.count(remove_log)!=1 or expr.count('Log[')!=1:
                raise ConsistencyError('Non-unique or unexpected logarithm')
            expr=expr.replace(remove_log,'1')
        elif 'Log[' in expr:
            raise ConsistencyError('Unexpected logarithm')
        if 'ϵ' in expr:
            raise ConsistencyError('Unexpected epsilon in residue')
        expr=_atomize(expr,table)
        if not TOKENIZED.fullmatch(expr):
            raise ConsistencyError('Unsupported tokens in expression')
        # No source-execution: parse a fixed, restricted arithmetic AST.
        try:
            tree=ast.parse(expr.replace('^','**'),mode='eval')
            coeff,exponents=_evaluate(tree.body)
        except (SyntaxError,ValueError,TypeError,OverflowError) as exc:
            raise ConsistencyError('Could not parse exact scalar monomial') from exc
        key=tuple(sorted(exponents.items()))
        result[key]=result.get(key,Fraction(0))+coeff
        if result[key]==0:
            del result[key]
    return result


def compare_residue_and_log(pole_residue: str, pole_free: str, stage: str) -> dict:
    if stage not in MASS_LOGS:
        raise ConsistencyError('Unknown heavy mass')
    residue=terms(pole_residue)
    finite=terms(pole_free)
    logs=tuple(t for t in finite if 'Log[' in t)
    constants=tuple(t for t in finite if 'Log[' not in t)
    if bool(residue)!=bool(logs):
        raise ConsistencyError('Pole/log inventories differ')
    table={}
    p=polynomial(residue,table=table)
    l=polynomial(logs,table=table,remove_log=MASS_LOGS[stage])
    if p!=l:
        raise ConsistencyError('Pole and logarithmic coefficients differ')
    return {
        'pole_terms':len(residue),'log_terms':len(logs),
        'constant_terms':len(constants),
        'nonzero_unique_opaque_pole_monomials':len(p),
        'pole_equals_log_coefficient_in_opaque_factor_algebra':True,
        'formal_d_by_d_ln_mu_of_pole_free_equals_2_pole':True,
        'physical_rge_consistency_certified':False,
    }


def finite_groups(finite_mass_expression: str, sector: str) -> dict:
    grouped:dict[str,list[str]]={}
    for t in terms(finite_mass_expression):
        labels=set(LAMBDA.findall(t))
        if len(labels)!=1:
            raise ConsistencyError(f'Expected exactly one quartic in {sector}: {labels}')
        label=next(iter(labels))
        grouped.setdefault(label,[]).append(t)
        if 'Log[' in t or 'ϵ' in t or 'Coupling[MSX' in t or 'Coupling[MSY' in t:
            raise ConsistencyError('Unexpected scale, pole, or mass dependence at threshold')
    table={}
    result={}
    for label, group in sorted(grouped.items()):
        p=polynomial(tuple(group),table=table)
        result[label]={
            'raw_additive_terms':len(group),
            'independent_opaque_monomials_after_exact_combining':len(p),
            'formal_expression_nonzero_in_opaque_factor_algebra':bool(p),
            'physical_coefficient_nonzero_certified':False,
        }
    return result
