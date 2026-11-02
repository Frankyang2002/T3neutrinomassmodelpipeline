"""Strict exact Matchete SparseArray/term evidence reader (no Wolfram runtime).

Reads the actual `CGRegistry` and `ScalarQuarticTerms` in the EFT1 source
seed. Every tensor entry is an exact SymPy number, not a float or inferred
Clebsch–Gordan convention. Does not assert RGBeta's internal AddQuartic factor.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
import hashlib
import json
import re
import sympy as sp

from .named_hhss_projection import _bracketed, _split_top_level, SourceInvariantError


class MatcheteCGError(SourceInvariantError):
    pass


def _list(s: str):
    s=s.strip()
    if not (s.startswith('{') and s.endswith('}')):
        raise MatcheteCGError('Expected literal Wolfram list')
    if len(s)==2:return []
    chunks=_split_top_level(s[1:-1])
    return [_list(t) if t.strip().startswith('{') else t.strip() for t in chunks]


def _integer(s: str)->int:
    if not re.fullmatch(r'[-+]?\d+',s):
        raise MatcheteCGError('Expected integer coordinate')
    return int(s)


def _exact_number(s:str)->sp.Expr:
    """Restricted source-number grammar: rationals, I, Sqrt, powers, signs.

    Security: no arbitrary Wolfram symbols, names, functions or code allowed.
    """
    s=s.strip().replace(' ','')
    if not s:raise MatcheteCGError('Empty exact number')
    if not re.fullmatch(r'[\dIiSqrt+*/^().,\[\]-]+',s):
        raise MatcheteCGError(f'Unsupported exact numeric expression {s[:80]}')
    s=re.sub(r'\bSqrt\[', 'sqrt(',s)
    s=s.replace(']',' )').replace('^','**')
    if '[' in s or ']' in s:raise MatcheteCGError('Unsupported source-number brackets')
    try:
        num=sp.sympify(s,locals={'I':sp.I,'sqrt':sp.sqrt},evaluate=True)
    except Exception as exc:
        raise MatcheteCGError(f'Invalid exact number {s}') from exc
    if num.free_symbols or not num.is_number or num.has(sp.zoo,sp.nan,sp.oo,-sp.oo):
        raise MatcheteCGError(f'Nonfinite or symbolic coefficient {s}')
    return sp.simplify(num)


@dataclass(frozen=True)
class ExactCG:
    name: str
    reps: tuple[str,...]
    dims: tuple[int,...]
    entries: tuple[tuple[tuple[int,...],sp.Expr],...]
    source_sha256: str

    def component(self, coordinates:tuple[int,...])->sp.Expr:
        return dict(self.entries).get(coordinates,sp.S.Zero)

    def invariant_summary(self)->dict:
        return {'name':self.name,'reps':list(self.reps),
                'dims':list(self.dims),'nonzero_components':len(self.entries),
                'norm_squared_exact':sp.sstr(sp.simplify(sum((sp.conjugate(v)*v for _,v in self.entries),sp.S.Zero))),
                'contains_explicit_complex_phase':any(sp.simplify(sp.im(v))!=0 for _,v in self.entries),
                'tensor_source_sha256':self.source_sha256,
                'tensor_exact_nonzero_components':[{'coordinates_zero_based':list(k), 'value_exact':sp.sstr(v)}
                                                   for k,v in self.entries],
                'tensor_exact_component_digest_sha256':hashlib.sha256(json.dumps(
                  [(list(k),sp.sstr(v)) for k,v in self.entries],separators=(',',':')).encode()).hexdigest()}


def _rep_dimension(rep:str)->int:
    rep=rep.strip()
    if rep.startswith('Bar[') and rep.endswith(']'):
        rep=rep[4:-1].strip()
    if rep=='SU2L[fund]':return 2
    match=re.fullmatch(r'T3BSMd([123])',rep)
    if match:return int(match.group(1))
    raise MatcheteCGError(f'Unknown native CG representation {rep}')


def parse_exact_cg(item:dict)->ExactCG:
    required={'Name','RepsInputForm','TensorInputForm'}
    if not required.issubset(item):raise MatcheteCGError('Missing CG registry fields')
    name=item['Name'];source=item['TensorInputForm']
    if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z]\w*',name):
        raise MatcheteCGError('Invalid CG name')
    reps=tuple(_split_top_level(item['RepsInputForm'].strip()[1:-1]))
    if not item['RepsInputForm'].startswith('{') or not item['RepsInputForm'].endswith('}'):
        raise MatcheteCGError('Invalid CG representation list')
    if not source.startswith('SparseArray['):raise MatcheteCGError(f'{name}: only explicit SparseArray supported')
    body,end=_bracketed(source,source.index('['))
    if end!=len(source):raise MatcheteCGError(f'{name}: trailing SparseArray text')
    fields=_split_top_level(body)
    if len(fields)!=4 or fields[0]!='Automatic':
        raise MatcheteCGError(f'{name}: unsupported SparseArray structure')
    dims=tuple(_integer(x) for x in _list(fields[1]))
    if not dims or any(d<1 or d>3 for d in dims):
        raise MatcheteCGError(f'{name}: invalid SU2 dimensions')
    if _exact_number(fields[2])!=0:raise MatcheteCGError(f'{name}: nonzero SparseArray default')
    metadata=_list(fields[3])
    if len(metadata)!=3 or _integer(metadata[0])!=1 or len(metadata[1])!=2:
        raise MatcheteCGError(f'{name}: unsupported CSR representation')
    offsets=[_integer(x) for x in metadata[1][0]]
    tails=metadata[1][1]
    values=[_exact_number(x) for x in metadata[2]]
    if len(reps)!=len(dims):raise MatcheteCGError(f'{name}: rep rank mismatch')
    if tuple(_rep_dimension(r) for r in reps)!=dims:
        raise MatcheteCGError(f'{name}: CG representation dimensions differ from SparseArray')
    if len(offsets)!=dims[0]+1 or offsets[0]!=0 or offsets[-1]!=len(values):
        raise MatcheteCGError(f'{name}: invalid CSR row offsets')
    if any(a>b for a,b in zip(offsets,offsets[1:])) or len(tails)!=len(values):
        raise MatcheteCGError(f'{name}: invalid CSR column data')
    decoded={}
    for row in range(dims[0]):
        for k in range(offsets[row],offsets[row+1]):
            tail=tails[k]
            if len(tail)!=len(dims)-1:raise MatcheteCGError(f'{name}: incomplete column tuple')
            coords=(row,)+tuple(_integer(x)-1 for x in tail)
            if any(i<0 or i>=d for i,d in zip(coords,dims)):
                raise MatcheteCGError(f'{name}: out of range CG coordinate')
            if coords in decoded:raise MatcheteCGError(f'{name}: duplicate CG coordinate')
            decoded[coords]=values[k]
    if not decoded or all(v==0 for v in decoded.values()):raise MatcheteCGError(f'{name}: zero CG')
    return ExactCG(name,reps,dims,tuple(sorted((k,v) for k,v in decoded.items() if v!=0)),
                   hashlib.sha256(source.encode()).hexdigest())


def parse_all_cgs(seed:dict)->dict[str,ExactCG]:
    if not isinstance(seed.get('CGRegistry'),list):raise MatcheteCGError('No native CG registry')
    cgs={}
    for entry in seed['CGRegistry']:
        cg=parse_exact_cg(entry)
        if cg.name in cgs:raise MatcheteCGError(f'Duplicate CG name {cg.name}')
        cgs[cg.name]=cg
    return cgs


_FIELD_RX=re.compile(r'(Bar\[)?Field\[(H|NewScalar1|NewScalar2)\s*,\s*Scalar\s*,')
_COUP_RX=re.compile(r'Coupling\[([^,\]\s]+)\s*,')
_CG_RX=re.compile(r'CG\[(?:Bar\[)?(T3\w+)')
_NAMEMAP={'H':(0,1),'NewScalar1':(2,3),'NewScalar2':(4,5)}


def inspect_quartic_term(term:dict,cgs:dict[str,ExactCG])->dict:
    source=term.get('TermInputForm')
    if not isinstance(source,str) or not source:raise MatcheteCGError('Missing quartic term input')
    hits=list(_FIELD_RX.finditer(source))
    counts=[0]*6
    for hit in hits:
        t=_NAMEMAP[hit.group(2)][1 if hit.group(1) else 0]
        counts[t]+=1
    if sum(counts)!=4:
        raise MatcheteCGError(f'Term {term.get("Index")}: expected exactly four complex scalar field factors; got {counts}')
    if tuple(sum(counts[k:k+2]) for k in (0,2,4))!=tuple(term.get(x) for x in
          ('HiggsFieldCount','Scalar1FieldCount','Scalar2FieldCount')):
        raise MatcheteCGError(f'Term {term.get("Index")}: declared field counts disagree with source')
    names=tuple(dict.fromkeys(_COUP_RX.findall(source)))
    cg_names=tuple(dict.fromkeys(_CG_RX.findall(source)))
    for name in cg_names:
        if name not in cgs:raise MatcheteCGError(f'Term {term.get("Index")}: undeclared CG {name}')
    listed=set(term.get('CGNames',[]))
    if listed!=set(cg_names):
        raise MatcheteCGError(f'Term {term.get("Index")}: CGNames differs from source {cg_names} vs {listed}')
    return {'term_index':term['Index'], 'counts':counts, 'native_coupling_symbols':list(names),
            'CGNames':list(cg_names),'term_source_sha256':hashlib.sha256(source.encode()).hexdigest(),
            'has_hermitian_conjugate_coupling':bool(re.search(r'Bar\[Coupling\[',source))}


def inspect_all_terms(seed:dict,cgs:dict[str,ExactCG])->list[dict]:
    rows=seed.get('ScalarQuarticTerms')
    if not isinstance(rows,list) or len(rows)!=seed.get('ScalarQuarticTermCount'):
        raise MatcheteCGError('Incorrect source quartic term count')
    indices=[r.get('Index') for r in rows]
    if len(indices)!=len(set(indices)):
        raise MatcheteCGError('Duplicate quartic term indices')
    return [inspect_quartic_term(r,cgs) for r in rows]


def candidate_raw_tensor(cg:ExactCG,source_types:tuple[int,...],dims_by_type:tuple[int,...]):
    """Enumerate unassigned raw CG axis placements consistent with source dimensions.

    This is a PRECONDITION DIAGNOSTIC, not a physical field-to-index map.
    Rep bars/dualities are intentionally not assigned by dimension matching.
    The caller must separately check actual Matchete Field/CG index linking.
    """
    import itertools
    dims=tuple(dims_by_type[t] for t in source_types)
    nonsinglet=[i for i,d in enumerate(dims) if d>1]
    if len(nonsinglet)!=len(cg.dims):return []
    routes=[]
    for assignment in itertools.permutations(nonsinglet):
        if tuple(dims[i] for i in assignment)!=cg.dims:continue
        raw=[]
        for coords in product(*(range(d) for d in dims)):
            raw.append(cg.component(tuple(coords[i] for i in assignment)))
        routes.append({'axis_to_external_slot':list(assignment), 'raw_tensor':sp.Matrix(raw)})
    return routes
