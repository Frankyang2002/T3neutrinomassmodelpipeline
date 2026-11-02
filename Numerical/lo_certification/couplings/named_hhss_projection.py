"""V4: exact audit of named RGBeta HHSS GroupInvariant expressions.

This code DOES NOT identify an RGBeta coupling with a Matchete coefficient:
RGBeta's index duality, AddQuartic internal normalization and the original
Matchete CG convention need independent validation. Uncertified projections
are diagnostic only; no physical EFT1 tensor is exported here.

Uses fixed, source-identified WL primitive vocabulary; unsupported expressions
fail closed, rather than using Mathematica evaluation or guessed CG factors.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import re
import sympy as sp

from .extra_hhss_geometry import DIMS, reference_invariants, _check_counts

FORMAT='t3-ws2-source-named-rgbeta-HHSS-projection-v4'
_TYPES={'H':0,'Bar @ H':1,'S1':2,'Bar @ S1':3,'S2':4,'Bar @ S2':5}
_DIMENSIONS={'H':2,'S1':None,'S2':None}


class SourceInvariantError(ValueError):
    pass


def _bracketed(s: str, start: int, opener: str='[', closer: str=']') -> tuple[str,int]:
    if start>=len(s) or s[start]!=opener:
        raise SourceInvariantError('Expected opening bracket')
    depth=0
    for j in range(start,len(s)):
        if s[j]==opener: depth+=1
        elif s[j]==closer:
            depth-=1
            if depth==0:return s[start+1:j],j+1
    raise SourceInvariantError('Unbalanced Wolfram source expression')


def _split_top_level(s:str)->list[str]:
    starts={'[':']','{':'}','(':')'}
    stack=[];rows=[];last=0
    for i,c in enumerate(s):
        if c in starts:stack.append(starts[c])
        elif stack and c==stack[-1]:stack.pop()
        elif c==',' and not stack:
            rows.append(s[last:i].strip());last=i+1
    if stack:raise SourceInvariantError('Unbalanced Wolfram function parameters')
    rows.append(s[last:].strip())
    return rows


def generated_declaration(source:str,name:str)->dict:
    """Extract literal AddQuartic and its GroupInvariant from local .wl.

    Compares registry invariant separately, so drift in the generated Wolfram
    source is not accepted on the strength of the JSON registry alone.
    """
    match=list(re.finditer(r'\bAddQuartic\[\s*'+re.escape(name)+r'\s*,',source))
    if len(match)!=1:raise SourceInvariantError(f'{name}: expected one AddQuartic declaration, found {len(match)}')
    op=source.index('[',match[0].start())
    body,_=_bracketed(source,op)
    parts=_split_top_level(body)
    if parts[0].strip()!=name:raise SourceInvariantError(f'{name}: declaration mismatch')
    fields=parts[1].strip()
    if not(fields.startswith('{') and fields.endswith('}')):raise SourceInvariantError('Unrecognized fields list')
    fields=_split_top_level(fields[1:-1]);fields=[re.sub(r'\s+',' ',x).strip() for x in fields]
    if any(f not in _TYPES for f in fields):raise SourceInvariantError(f'{name}: unknown field type {fields}')
    groups=[p.split('->',1)[1].strip() for p in parts[2:] if p.strip().startswith('GroupInvariant')]
    if len(groups)!=1:raise SourceInvariantError(f'{name}: expected exactly one GroupInvariant')
    expr=groups[0]
    if '&' not in expr:raise SourceInvariantError(f'{name}: GroupInvariant not a slot function')
    expr=expr.replace('&','').strip()
    expr=_strip_parens(expr)
    conj=[p for p in parts[2:] if p.strip().startswith('SelfConjugate')]
    return {'symbol':name,'fields':fields,'group_invariant_expression':expr,
            'self_conjugate_option':conj[0].split('->',1)[1].strip() if conj else 'default'}


def _strip_parens(expr:str)->str:
    s=expr.strip()
    while s.startswith('('):
        inner,end=_bracketed(s,0,'(',')')
        if end!=len(s):break
        s=inner.strip()
    return s


@lru_cache(maxsize=None)
def _literal_wl(expr:str)->str:
    """Finite whitelist translator. No evaluation of arbitrary WL/Python code."""
    expr=expr.strip()
    allowed=('delS2','del','eps','T3RGBetaS2TripleInvariant')
    rx=re.compile(r'\b('+'|'.join(allowed)+r')\s*\[')
    # Replace Wolfram contractions in non-nested arguments from right to left.
    while True:
        matches=list(rx.finditer(expr))
        if not matches:break
        m=matches[-1]
        args,end=_bracketed(expr,m.end()-1)
        parts=_split_top_level(args)
        name=m.group(1)
        if name=='del':
            if len(parts)!=3 or parts[0] not in ('SU2L @ fund','SU2L @ S2'):
                raise SourceInvariantError(f'Unsupported del structure: {args}')
            fn='D2' if parts[0]=='SU2L @ fund' else 'D3'
            out=f'{fn}({parts[1]},{parts[2]})'
        elif name=='eps':
            if len(parts)!=3 or parts[0]!='SU2L @ fund':
                raise SourceInvariantError(f'Unsupported eps structure: {args}')
            out=f'E2({parts[1]},{parts[2]})'
        elif name=='delS2':
            if len(parts)!=4 or parts[0]!='SU2L':
                raise SourceInvariantError(f'Unsupported delS2 structure: {args}')
            out=f'CG3({parts[1]},{parts[2]},{parts[3]})'
        else:
            if len(parts)!=3:raise SourceInvariantError(f'Unsupported triple invariant: {args}')
            out=f'TRI({parts[0]},{parts[1]},{parts[2]})'
        expr=expr[:m.start()]+out+expr[end:]
    expr=re.sub(r'#([1-4])\b',r'i\1',expr)
    if '[' in expr or ']' in expr or '#' in expr or '@' in expr:
        raise SourceInvariantError(f'Unsupported Wolfram source: {expr}')
    if not re.fullmatch(r'[\s\d+\-*/().,iaDGECTR234I]+',expr):
        raise SourceInvariantError(f'Unexpected characters in invariant: {expr}')
    return expr


@lru_cache(maxsize=None)
def _tri(a:int,b:int,c:int):
    """Literal T3RGBetaS2TripleInvariant definition in T3RGBetaModel.wl.

    The last eight SU(2) doublet index contractions use the six epsilons
    shown in that source; no unverified normalizing constant is inserted.
    """
    cg=_cg3;eps=_eps2
    val=sp.S.Zero
    for i,j,k,l,m,n in product(range(2),repeat=6):
        val+=cg(a,i,j)*cg(b,k,l)*cg(c,m,n)*eps(i,k)*eps(j,m)*eps(l,n)
    return sp.simplify(val)


def _eps2(i,j):return sp.Integer(1 if (i,j)==(0,1) else -1 if (i,j)==(1,0) else 0)


def _cg3(a,i,j):
    if a==0 and i==j==0:return sp.S.One
    if a==2 and i==j==1:return sp.S.One
    if a==1 and {i,j}=={0,1}:return 1/sp.sqrt(2)
    return sp.S.Zero


@lru_cache(maxsize=None)
def _evaluate_raw(expr:str, indices:tuple[int,...]):
    parsed=_literal_wl(expr)
    def evaluate(a):
        loc={'D2':lambda i,j:sp.Integer(i==j),
             'D3':lambda i,j:sp.Integer(i==j),
             'E2':_eps2,'CG3':_cg3,'TRI':_tri,
             **{f'i{k+1}':v for k,v in enumerate(indices)},'a':a}
        # SymPy's parser receives only a restricted expression and named
        # numeric index functions; no arbitrary builtins or eval permitted.
        return sp.sympify(parsed,locals=loc,evaluate=True)
    return sp.simplify(sum((evaluate(a) for a in range(3)),sp.S.Zero) if re.search(r'\ba\b',parsed) else evaluate(0))


def _symmetrize(raw:sp.Matrix,types:tuple[int,...],dims:tuple[int,...])->sp.Matrix:
    out=raw
    from .extra_hhss_geometry import _permute_tensor
    for typ in sorted(set(types)):
        positions=[i for i,t in enumerate(types) if t==typ]
        if len(positions)==2:out=(out+_permute_tensor(out,dims,*positions))/2
    return sp.simplify(out)


def _source_vector(decl:dict,model:str,types:tuple[int,...],dims:tuple[int,...])->sp.Matrix:
    fields=decl['fields']
    sdims=(2,2,*DIMS[model[0]])
    def fielddim(t:int):return (2,2,sdims[2],sdims[2],sdims[3],sdims[3])[t]
    ftypes=[_TYPES[f] for f in fields]
    if sorted(ftypes)!=sorted(types):
        raise SourceInvariantError(f'{decl["symbol"]}: field content differs from registry')
    # Ordered fields in the generated AddQuartic are in a different order
    # from the deterministic order of reference_invariants.
    source_positions={}
    for i,t in enumerate(ftypes):source_positions.setdefault(t,[]).append(i)
    original_positions={t:0 for t in ftypes}
    inverse=[]
    for t in types:
        index=original_positions[t];inverse.append(source_positions[t][index]);original_positions[t]+=1
    result=[]
    for coord in product(*(range(d) for d in dims)):
        srcidx=[0]*4
        for reference_idx,srcpos in enumerate(inverse):srcidx[srcpos]=coord[reference_idx]
        result.append(_evaluate_raw(decl['group_invariant_expression'],tuple(srcidx)))
    return _symmetrize(sp.Matrix(result),types,dims)


def _dual_matrix(d:int)->sp.Matrix:
    if d==1:return sp.eye(1)
    if d==2:return sp.Matrix([[0,1],[-1,0]])
    if d==3:return sp.Matrix([[0,0,1],[0,-1,0],[1,0,0]])
    raise SourceInvariantError('Unsupported SU2 dimension')


def _dual_variants(raw:sp.Matrix, types:tuple[int,...],dims:tuple[int,...]):
    """Search exact pseudoreality maps, but label them UNVERIFIED choices.

    A raw RGBeta tensor that is not a canonical-invariant coefficient may
    become one after a dual map. Such a map is not certified by a Ward test.
    """
    unique=[t for t in sorted(set(types)) if dims[types.index(t)]>1]
    # Shared type gets shared map, so exchange symmetry is preserved.
    for switches in product((False,True),repeat=len(unique)):
        tags=dict(zip(unique,switches))
        mats=[_dual_matrix(d) if tags.get(t,False) else sp.eye(d) for t,d in zip(types,dims)]
        tensor=sp.kronecker_product(*mats)
        yield tuple(t for t in unique if tags[t]),_symmetrize(tensor*raw,types,dims)


def projection_report(model:str,counts:tuple[int,...],decl:dict,*,max_variants:int=16)->dict:
    counts=_check_counts(counts)
    types,dims,basis=reference_invariants(*DIMS[model[0]],counts)
    if len(decl['fields'])!=4:raise SourceInvariantError('Expected four external scalar fields')
    raw=_source_vector(decl,model,types,dims)
    B=sp.Matrix.hstack(*basis)
    raw_nonzero=any(sp.simplify(x)!=0 for x in raw)
    if not raw_nonzero:
        raise SourceInvariantError('Generated GroupInvariant source is identically zero')
    variants=[]
    for tags,vector in _dual_variants(raw,types,dims):
        if len(variants)>=max_variants:raise SourceInvariantError('Too many orientation possibilities')
        try:
            solved=B.gauss_jordan_solve(vector)[0]
            if any(sp.simplify(x)!=0 for x in B*solved-vector):raise ValueError('nonzero residual')
            coefficients=[sp.sstr(sp.simplify(x)) for x in solved]
            in_span=True
        except (ValueError,sp.polys.polyerrors.PolificationFailed):
            coefficients=None;in_span=False
        variants.append({'dual_map_field_types':[
            ('H','BarH','S1','BarS1','S2','BarS2')[t] for t in tags],
            'in_exact_canonical_SU2_reference_span':in_span,
            'source_to_reference_polynomial_coefficients':coefficients})
    matches=[v for v in variants if v['in_exact_canonical_SU2_reference_span']]
    return {'rgbeta_symbol':decl['symbol'],'counts':list(counts),
            'source_fields':decl['fields'],
            'group_invariant_expression':decl['group_invariant_expression'],
            'self_conjugate_option':decl['self_conjugate_option'],
            'raw_component_tensor_nonzero':raw_nonzero,
            'reference_dimension':len(basis),
            'tested_orientation_candidates':len(variants),
            'valid_canonical_invariant_orientation_count':len(matches),
            'candidate_orientation_projection_maps':matches,
            'physical_Matchete_RGBeta_coupling_dictionary':None,
            'status':('conditional_exact_RGBeta_index_projection_orientation_ambiguous'
                      if matches else 'blocked_source_tensor_not_in_canonical_invariant_span'),
            'reason_not_certified':'SU2 Ward identity alone does not verify RGBeta/Machete CG coordinate and phase orientation, AddQuartic internal factors or coefficient RG evolution',
            'safe_physical_input':False}


def named_hhss_sources(root:Path,model:str,registry:dict)->dict:
    """Fail-closed cross-check against local Wolfram declaration and registry."""
    family=model[0]
    extra_path=root/f'RGE/general/generated/{model}_extra_quartics.wl'
    common_path=root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    sources={'new':extra_path.read_text(encoding='utf-8'),
             'preexisting':common_path.read_text(encoding='utf-8')}
    sectors=[]
    for category,entries in (('new',registry.get('new_rgbeta_symbols',[])),
                             ('preexisting',registry.get('preexisting_rgbeta_symbols',[]))):
        for item in entries:
            counts=tuple(item['counts'])
            if sum(counts[:2])!=2 or sum(counts[2:])!=2:continue
            name=item['name'] if category=='new' else item['existing_name']
            decl=generated_declaration(sources[category],name)
            if [decl['fields'].count(t) for t in ('H','Bar @ H','S1','Bar @ S1','S2','Bar @ S2')]!=list(counts):
                raise SourceInvariantError(f'{model}/{name}: generated WL fields disagree with registry counts')
            if category=='new':
                literal=re.sub(r'\s+','',_strip_parens(decl['group_invariant_expression']))
                declared=re.sub(r'\s+','',_strip_parens(item['invariant']))
                if literal!=declared:
                    raise SourceInvariantError(f'{model}/{name}: generated WL invariant differs from registry')
            sectors.append({'category':category,'counts':counts,'declaration':decl})
    return {'model':model,'source_sha256':{
              'generated_extra_wl':hashlib.sha256(sources['new'].encode()).hexdigest(),
              'core_rgbeta_wl':hashlib.sha256(sources['preexisting'].encode()).hexdigest()},
            'sectors':sectors}
