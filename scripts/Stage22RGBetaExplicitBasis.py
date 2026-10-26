"""Diagnostic ONLY: evaluate exported RGBeta SU(2) invariants in explicit Cartesian
component labels and compare to Matchete. Never asserts physical normalisation.

This is intentionally separate from RunGeneralT3Stage22 and does not write
rgbeta_native_tensors.json, whose semantics require independently certified
component-basis identification and native potential normalisation.
"""
from __future__ import annotations
import argparse
import json
import re
import hashlib
from collections import defaultdict
from itertools import product
from pathlib import Path

from scripts.Stage22TensorAlgebra import exact, solve_basis, TensorCertificateError

SOURCE_ORDER = ('H', 'Hdag', 'S1', 'S1dag', 'S2', 'S2dag')
_FIELD_RE = re.compile(r'^(?:H|S1|S2|Bar\[(?:H|S1|S2)\])$')
_ATOM = re.compile(r'(del|eps)\[SU2L\[(fund|S2)\],\s*(t3i[1-4]),\s*(t3i[1-4])\]')
_BRIDGE = re.compile(r'delS2\[|T3RGBetaS2TripleInvariant\[')


def parse_fields(form):
    if not isinstance(form, str) or not form.startswith('{') or not form.endswith('}'):
        raise ValueError(f'Not an explicit RGBeta field list: {form!r}')
    fields = [s.strip() for s in form[1:-1].split(',')]
    if len(fields) != 4 or any(_FIELD_RE.fullmatch(s) is None for s in fields):
        raise ValueError(f'Invalid RGBeta field list: {form!r}')
    return [f'{s[4:-1]}dag' if s.startswith('Bar[') else s for s in fields]


def safe_arith(expr):
    """Exact Mathematica scalar arithmetic after delta/epsilon replacement."""
    if '[' in expr or ']' in expr or '$' in expr:
        raise ValueError('Unreduced Wolfram tensor expression')
    return exact(expr)


def eval_component(expression, index_tuple):
    """Restricted index expansion: no arbitrary Wolfram/Python evaluation."""
    if not isinstance(expression,str) or _BRIDGE.search(expression):
        raise ValueError('Requires triplet CG component basis')
    # Only 1, del, eps, Rational arithmetic in supported expressions.
    # Replace each exact named atom by Kronecker or 2D epsilon.
    def atom(match):
        typ, rep, a, b = match.groups()
        i, j = index_tuple[int(a[-1])-1], index_tuple[int(b[-1])-1]
        if typ == 'del': return '1' if i == j else '0'
        if rep != 'fund': raise ValueError('Nonfundamental epsilon')
        return '1' if (i,j)==(0,1) else '-1' if (i,j)==(1,0) else '0'
    replaced = _ATOM.sub(atom,expression)
    if re.search(r'\b(del|eps|t3i\d|SU2L)\b',replaced) or re.search(r'[A-Za-z_]',replaced):
        raise ValueError(f'Unsupported tensor primitive: {expression[:70]}')
    return safe_arith(replaced)


def reconstructed_tensor(native, mat_reference):
    field_names = parse_fields(native['native_fields_inputform'])
    if field_names != mat_reference['species']:
        raise ValueError(f'Field order mismatch: RGBeta {field_names}, Matchete {mat_reference["species"]}')
    dims = mat_reference['dimensions']
    expression = native['native_invariant_expression_inputform']
    # Explicitly not enough to identify triplet CG bases or triplet complex conjugation.
    if _BRIDGE.search(expression):
        raise ValueError('Unaligned RGBeta delS2 triplet CG convention')
    return {**{k: mat_reference[k] for k in ['counts','dimensions']},
            'symbol':native['symbol'],
            'species':field_names,
            'components':[str(eval_component(expression, i)) for i in product(*(range(d) for d in dims))]}


def verify_sources(mat_json, inv_json, declarations):
    if mat_json.get('quartic_slot_count') != 50 or mat_json.get('model_count') != 16:
        raise ValueError('Incorrect Matchete inventory')
    if inv_json.get('total_slots') != 50:
        raise ValueError('Incorrect RGBeta invariant inventory')
    if set(mat_json['models']) != set(inv_json['models']) or set(mat_json['models']) != set(declarations['models']):
        raise ValueError('Model inventories disagree')
    for model in mat_json['models']:
        mats = mat_json['models'][model]['tensors']
        invs = inv_json['models'][model]['invariants']
        decs = declarations['models'][model]['declarations']
        if len(mats) != len(invs) or len(invs) != len(decs):
            raise ValueError(f'{model}: inconsistent slot counts')
        dd = {x['symbol']:x for x in decs}
        if len(dd) != len(decs): raise ValueError(f'{model}: duplicate declarations')
        for inv in invs:
            symbol = inv['symbol']
            if symbol not in dd: raise ValueError(f'{model}: missing native declaration for {symbol}')
            d = dd[symbol]
            for a,b in [('fields_inputform','native_fields_inputform'),
                        ('invariant_inputform','native_invariant_function_inputform'),
                        ('coupling_factor_inputform','native_coupling_function_inputform'),
                        ('unique_arrangements_inputform','native_unique_arrangements_inputform')]:
                if d[a] != inv[b]:raise ValueError(f'{model}/{symbol}: incoherent {a}')
            if 'Missing[' in str(inv) or 'Missing[' in str(d):
                raise ValueError(f'{model}/{symbol}: native association key missing')
            if not re.fullmatch(r'6\*Matrix\['+re.escape(symbol)+r'\]\[##1\] &\s*',d['coupling_factor_inputform']):
                raise ValueError(f'{model}/{symbol}: unexpected RGBeta coupling normalization; inspect source')
    return True


def compare(mat_json, inv_json, declarations):
    verify_sources(mat_json, inv_json, declarations)
    output={'format':'stage22-rgbeta-explicit-geometry-proposals-v1',
            'status':'GeometryCandidatesNotCertifiedCouplingMaps',
            'physical_normalisation_verified':False,
            'native_triplet_basis_verified':False,
            'rgbeta_lagrangian_prefactor_verified':False,
            'epsilon_orientation':'epsilon[0,1]=+1 (diagnostic convention)',
            'models':{},'totals':defaultdict(int)}
    for model in sorted(mat_json['models']):
        mb=defaultdict(list)
        for t in mat_json['models'][model]['tensors']:
            mb[tuple(t['counts'])].append(t)
        rb=defaultdict(list)
        for t in inv_json['models'][model]['invariants']:
            # Reused RGBeta vertices do NOT share Matchete symbols: match by the
            # explicitly ordered native field list, not by a guessed name.
            fields = parse_fields(t['native_fields_inputform'])
            counts=tuple(fields.count(s) for s in SOURCE_ORDER)
            rb[counts].append(t)
        if set(mb)!=set(rb):raise ValueError(f'{model}: quartic block counts disagree')
        outputs=[]
        for counts,mats in sorted(mb.items()):
            invs=rb[counts]; ref=mats[0]
            code=''.join(str(x) for x in counts)
            o={'block':code, 'matchete_symbols':[x['symbol'] for x in mats],
               'rgbeta_symbols':[x['symbol'] for x in invs],
               'status':'Uncertified'}
            try:
                if any(x['species']!=ref['species'] for x in mats):raise ValueError('Inconsistent Matchete species order')
                rgb=[reconstructed_tensor(t,ref) for t in invs]
                answer=solve_basis(mats,rgb)
                o.update({'status':'ExactGeometryCandidate',
                          'basis_change_r_in_m':answer['r_in_m'],
                          'basis_change_m_in_r':answer['m_in_r'],
                          'rank':len(mats),'exact_component_residual':'0',
                          'physical_coupling_map_certified':False,
                          'warning':'Formal RGBeta contraction evaluated in chosen component labels; native component-basis and physical prefactors remain unverified.'})
            except (ValueError,TensorCertificateError) as exc:
                o['reason']=str(exc)
                # Explicitly differentiate mismatch and remaining unknown triplet basis.
                if 'Unaligned RGBeta delS2' in str(exc):o['status']='TripletBasisPending'
                else:o['status']='ExplicitGeometryMismatch'
            outputs.append(o)
            output['totals'][o['status']]+=len(mats)
            output['totals']['slots']+=len(mats)
            output['totals']['blocks']+=1
        output['models'][model]=outputs
    output['totals']=dict(output['totals'])
    if output['totals']['slots']!=50 or output['totals']['blocks']!=46:
        raise ValueError('Wrong structural totals')
    return output


def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--input-dir',type=Path,default=Path('output/general_t3_stage22'))
    pa.add_argument('--output',type=Path,default=None)
    ns=pa.parse_args()
    def load(name):return json.loads((ns.input_dir/name).read_text(encoding='utf-8-sig'))
    summary=compare(load('matchete_tensors.json'),load('rgbeta_native_invariants.json'),load('rgbeta_native_declarations.json'))
    summary['input_sha256']={name:hashlib.sha256((ns.input_dir/name).read_bytes()).hexdigest()
                             for name in ('matchete_tensors.json','rgbeta_native_invariants.json','rgbeta_native_declarations.json')}
    dest=ns.output or ns.input_dir/'stage22_explicit_geometry_candidates.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print('Stage 22 explicit SU(2) geometry:',summary['totals']); print('Physical coupling mapping remains UNCERTIFIED'); print('Report:',dest)

if __name__=='__main__':main()
