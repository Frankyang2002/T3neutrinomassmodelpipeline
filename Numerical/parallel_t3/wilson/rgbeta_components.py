"""Conditional RGBeta SU(2) invariant component construction (not native CG certification).

RGBeta keeps group tensors symbolic; applying its pure invariant function to
*integer* slot labels is NOT a numeric component evaluation. In particular,
RGBeta's del[rep,a,a] evaluates to Dim[rep] by Einstein summation rules.

This module evaluates the *source-exported, unmodified* InputForm expressions
with four independent component slots in one declared orthonormal SU(2) frame:
  eps[fund,1,2]=+1, eps[fund,2,1]=-1
  S2[1] = e1 e1; S2[2] = (e1 e2 + e2 e1)/sqrt(2); S2[3] = e2 e2.

RGBeta's abstract triplet index does not select this basis or its phases.
The outcome is a useful conditional component reference, NOT a certified
Matchete<->RGBeta mapping, coupling normalization or RGBeta native-frame tensor.

The 50-slot declaration inventory uses del, eps, delS2 and the project-defined
T3RGBetaS2TripleInvariant. The latter is expanded using its *actual* definition
in RGE/running/rgbeta/T3RGBetaModel.wl (three delS2 and three eps).
An explicitly repeated dummy triplet label a is summed over 1..3. Parsing uses a closed AST interpreter; no eval/exec.
"""
from __future__ import annotations

import ast
from hashlib import sha256
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping


FORMAT = 't3-parallel-a-rgbeta-declared-invariants-conditional-frame-v2'
REQUIRED_MODELS = frozenset({
    'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
    'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
    'E_-4','E_-2','E_+0','E_+2',
})


class ComponentFrameError(ValueError):
    """A declaration is missing, unsafe, or inconsistent with its model."""


def delta(rep: str, i: int, j: int) -> float:
    if rep not in ('fund', 'triplet'):
        raise ComponentFrameError('Unknown delta representation: ' + str(rep))
    dim = 2 if rep == 'fund' else 3
    if not (1 <= i <= dim and 1 <= j <= dim):
        raise ComponentFrameError('Delta index out of range')
    return float(i == j)


def epsilon(rep: str, i: int, j: int) -> float:
    if rep != 'fund' or not (1 <= i <= 2 and 1 <= j <= 2):
        raise ComponentFrameError('eps requires fundamental SU(2) indices')
    return float((i, j) == (1, 2)) - float((i, j) == (2, 1))


def symmetric_cg(group: str, a: int, i: int, j: int) -> float:
    if group != 'SU2L' or not (1 <= a <= 3 and 1 <= i <= 2 and 1 <= j <= 2):
        raise ComponentFrameError('delS2 requires SU2L S2 x fund x fund')
    if a == 1:
        return float(i == j == 1)
    if a == 3:
        return float(i == j == 2)
    return float((i, j) in ((1, 2), (2, 1))) / math.sqrt(2)


@lru_cache(maxsize=None)
def triple_invariant(a: int, b: int, c: int) -> float:
    """Exact source structure of T3RGBetaS2TripleInvariant in chosen CG frame.

    Source: T3RGBetaModel.wl lines 462-470; index labels i,j,k,l,m,n
    are internal and summed, independently of the four external slots.
    """
    if not all(type(x) is int and 1<=x<=3 for x in (a,b,c)):
        raise ComponentFrameError('Triple invariant requires three triplet indices')
    result=0.0
    for i,j,k,l,m,n in product((1,2),repeat=6):
        result+=(symmetric_cg('SU2L',a,i,j)*
                 symmetric_cg('SU2L',b,k,l)*
                 symmetric_cg('SU2L',c,m,n)*
                 epsilon('fund',i,k)*epsilon('fund',j,m)*epsilon('fund',l,n))
    return result


def normalise_inputform(expr: str) -> str:
    if not isinstance(expr, str):
        raise ComponentFrameError('Invariant must be an InputForm string')
    raw = expr.strip()
    if not raw.endswith('&'):
        raise ComponentFrameError('Expected a native pure function ending in &')
    raw = raw[:-1].strip()
    # The supported real SU(2) invariant grammar was verified against the
    # source-exported Stage22 declaration list, not inferred from the symbol.
    raw = raw.replace('del[', 'delta[')
    raw = raw.replace('SU2L[fund]', 'fund')
    raw = raw.replace('SU2L[S2]', 'triplet')
    raw = re.sub(r'(?<![A-Za-z0-9_])SU2L(?![A-Za-z0-9_])', 'SU2L', raw)
    raw = re.sub(r'#([1-4])(?!\d)', r'idx\1', raw)
    raw = raw.replace('[', '(').replace(']', ')').replace('^', '**')
    return raw


def evaluate_inputform(expr: str, indices: tuple[int, int, int, int]) -> float:
    """Evaluate the *abstract invariant* in the explicit chosen component frame."""
    if len(indices) != 4 or any(type(x) is not int for x in indices):
        raise ComponentFrameError('Exactly four integer component slots required')
    source = normalise_inputform(expr)
    try:
        node = ast.parse(source, mode='eval')
    except (SyntaxError, ValueError) as exc:
        raise ComponentFrameError('Unsupported native invariant syntax') from exc
    env = {f'idx{i+1}': k for i, k in enumerate(indices)}
    env.update({'fund': 'fund', 'triplet': 'triplet', 'SU2L': 'SU2L'})
    functions = {'delta': delta, 'eps': epsilon, 'delS2': symmetric_cg,
                 'T3RGBetaS2TripleInvariant': triple_invariant}

    def run(x: ast.AST):
        if isinstance(x, ast.Constant) and type(x.value) in (int, float):
            return x.value
        if isinstance(x, ast.Name) and x.id in env:
            return env[x.id]
        if isinstance(x, ast.UnaryOp) and type(x.op) in (ast.UAdd, ast.USub):
            v = run(x.operand)
            if not isinstance(v, (int, float)):
                raise ComponentFrameError('Non-numeric unary argument')
            return v if isinstance(x.op, ast.UAdd) else -v
        if isinstance(x, ast.BinOp) and type(x.op) in (ast.Add, ast.Sub, ast.Mult, ast.Div):
            a, b = run(x.left), run(x.right)
            if type(a) not in (int, float) or type(b) not in (int, float):
                raise ComponentFrameError('Non-numeric binary operand')
            if isinstance(x.op, ast.Add): return a+b
            if isinstance(x.op, ast.Sub): return a-b
            if isinstance(x.op, ast.Mult): return a*b
            if b == 0: raise ComponentFrameError('Zero denominator')
            return a/b
        if isinstance(x, ast.Call) and isinstance(x.func, ast.Name) and x.func.id in functions:
            if x.keywords or any(isinstance(a, ast.Starred) for a in x.args):
                raise ComponentFrameError('Unexpected call arguments')
            return functions[x.func.id](*(run(a) for a in x.args))
        raise ComponentFrameError('Unsupported AST node: ' + type(x).__name__)

    try:
        # RGBeta uses Einstein contraction over a project-defined triplet
        # index 'a' in T3RGBetaMixInvariant[2,2] and [3,3].
        # A free 'a' symbol appears only as a dummy and is summed here.
        summed = any(isinstance(z, ast.Name) and z.id == 'a' for z in ast.walk(node))
        if summed:
            value=0.0
            for internal_index in (1,2,3):
                env['a']=internal_index
                value+=run(node.body)
        else:
            value=run(node.body)
    except TypeError as exc:
        raise ComponentFrameError('Wrong native invariant arity') from exc
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ComponentFrameError('Non-real or nonfinite invariant component')
    return float(value)


def species_list(native_fields: str) -> list[str]:
    if not isinstance(native_fields, str):
        raise ComponentFrameError('Missing Fields InputForm')
    m = re.fullmatch(r'\{\s*(.*?)\s*\}', native_fields.strip())
    if not m:
        raise ComponentFrameError('Fields must be a 4-slot native list')
    fields = [f.strip().replace(' ', '') for f in m.group(1).split(',')]
    if len(fields) != 4 or any(f not in (
            'H','S1','S2','Bar[H]','Bar[S1]','Bar[S2]') for f in fields):
        raise ComponentFrameError('Unsupported scalar fields')
    return fields


def model_dimensions(setting: Mapping[str, Any]) -> dict[str, int]:
    dims = setting.get('dimensions')
    if not isinstance(dims, list) or len(dims) != 3 or any(
            type(i) is not int or i not in (1,2,3) for i in dims):
        raise ComponentFrameError('Missing/invalid model dimension triple')
    return {'H': 2, 'S1': dims[0], 'S2': dims[1]}


def validate_frame() -> dict[str, float]:
    # Orthonormality of triplet tensor axes. The projector is the symmetric
    # subspace of 2 x 2 (same normalization as RGBeta's symbolic identity).
    err_orth = 0.0
    err_proj = 0.0
    for a in range(1,4):
        for b in range(1,4):
            val = sum(symmetric_cg('SU2L',a,i,j)*symmetric_cg('SU2L',b,i,j)
                      for i in (1,2) for j in (1,2))
            err_orth = max(err_orth, abs(val-float(a==b)))
    for i,j,k,l in product((1,2), repeat=4):
        lhs = sum(symmetric_cg('SU2L',a,i,j)*symmetric_cg('SU2L',a,k,l)
                  for a in (1,2,3))
        rhs = 0.5*(float(i==k and j==l)+float(i==l and j==k))
        err_proj = max(err_proj,abs(lhs-rhs))
    return {'s2_orthonormality_max_abs':err_orth, 's2_projector_max_abs':err_proj}


def reconstruct(declarations: Mapping[str, Any], registry: Mapping[str, Any],
                *, declarations_sha256: str, registry_sha256: str,
                strict_all16: bool = True) -> dict[str, Any]:
    if declarations.get('format') != 'stage22-native-rgbeta-declarations-v1' or \
       declarations.get('status') != 'NativeDeclarationMetadataOnly':
        raise ComponentFrameError('Requires actual native Stage22 declarations, NOT probe v1')
    if not isinstance(registry.get('models'), dict) or not isinstance(declarations.get('models'), dict):
        raise ComponentFrameError('Model dictionary missing')
    names = set(declarations['models'])
    if strict_all16 and (names != REQUIRED_MODELS or set(registry['models']) != REQUIRED_MODELS):
        raise ComponentFrameError('Require all 16 registered models')
    if not names.issubset(registry['models']):
        raise ComponentFrameError('Unknown model in declarations')
    all_records: dict[str, dict] = {}
    total = 0
    for name in sorted(names):
        settings = registry['models'][name]
        dimensions = model_dimensions(settings)
        block = declarations['models'][name]
        records = block.get('declarations')
        if block.get('tensors_not_exported') is not True or not isinstance(records,list):
            raise ComponentFrameError(f'{name}: expected native declaration-only block')
        expected = [x['name'] for x in settings['new_rgbeta_symbols']]
        expected += [x['existing_name'] for x in settings['preexisting_rgbeta_symbols']]
        actual = [rec.get('symbol') for rec in records]
        if len(actual) != len(set(actual)) or set(actual) != set(expected):
            raise ComponentFrameError(f'{name}: native symbol inventory mismatch')
        if settings['total_additional_quartics'] != len(records):
            raise ComponentFrameError(f'{name}: quartic slot count mismatch')
        result = []
        for rec in records:
            fields = species_list(rec['fields_inputform'])
            dims = [dimensions[f[4:-1] if f.startswith('Bar[') else f] for f in fields]
            expr = rec['invariant_inputform']
            components = [evaluate_inputform(expr,tuple(idx)) for idx in
                          product(*(range(1,d+1) for d in dims))]
            if not components or max(abs(x) for x in components) == 0:
                raise ComponentFrameError(f'{name}:{rec["symbol"]}: zero invariant')
            result.append({
                'symbol':rec['symbol'],
                'status':'Constructed_ConditionalSU2Frame_NativeNormalizationUncertified',
                'species': fields, 'dimensions':dims,
                'components':components,
                'component_convention':'orthonormal_S2_symmetric_pair_and_eps12_plus1',
                'native_invariant_inputform': expr,
                'native_coupling_function':rec['coupling_factor_inputform'],
                'native_unique_arrangements':rec['unique_arrangements_inputform'],
                'native_self_conjugate':rec['self_conjugate'],
                'prefactor_certified':False,
                'native_triplet_frame_certified':False,
                'matchete_scalar_basis_aligned':False,
                'physical_wilson_normalisation_certified':False,
            })
        total+=len(result)
        all_records[name]={'records':result,'n_records':len(result),
                           'status':'ConditionalComponents_NotPhysicalDictionary'}
    if strict_all16 and total != 50:
        raise ComponentFrameError('Expected exactly 50 additional quartic slots')
    checks=validate_frame()
    if max(checks.values())>1e-14:
        raise ComponentFrameError('Orthonormality/projector check failed')
    return {
        'format':FORMAT,'status':'ConditionalFrameConstructed_NotRGBetaNativeOrMatcheteCertified',
        'model_count':len(all_records),'quartic_slot_count':total,
        'models':all_records,'frame_checks':checks,
        'provenance':{'declarations_sha256':declarations_sha256,
                      'registry_sha256':registry_sha256},
        'interpretation':{
            'RGBeta_symbolic_invariant_provenance':True,
            'component_axes_are_source_native':False,
            'orientation_eps12_plus1_is_a_choice':True,
            'S2_orthonormal_triplet_basis_is_a_choice':True,
            'native_coupling_prefactor_applied':False,
            'potential_polynomial_normalized':False,
            'physical_matchete_to_RGBeta_dictionary_certified':False,
            'gamma_cc_to_weyl_signed_factor_certified':False,
        }}


def load_reconstruct(declarations_path: str | Path,
                     registry_path: str | Path, *, strict_all16: bool=True) -> dict[str,Any]:
    d=Path(declarations_path).read_bytes()
    r=Path(registry_path).read_bytes()
    return reconstruct(json.loads(d.decode('utf-8-sig')),
                       json.loads(r.decode('utf-8-sig')),
                       declarations_sha256=sha256(d).hexdigest(),
                       registry_sha256=sha256(r).hexdigest(),
                       strict_all16=strict_all16)
