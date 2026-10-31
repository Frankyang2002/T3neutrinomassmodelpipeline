#!/usr/bin/env python3
"""Workstream C: independently verify native SU(2) arrays and close the EFT1
component-equivalence/Hermiticity handoff without asserting physical Y matching.

This diagnostic does NOT modify the Lagrangian or calculate a renormalized
C5. All evidence must be pre-existing output from the native Wolfram fork.

Usage (from repository root):
 python scripts/ResolveT3BCEFT1HermiticityAndYGate.py \
  --residual-archive T3BC_FullEFT1Residual_Outputs.zip \
  --native-dir output/audits/t3_bc_su2_component_equivalence_v1/native \
  --matching-archive xy_fresh_native_evidence_v3.zip \
  --output output/audits/t3_bc_hermiticity_y_handoff_v1/summary.json
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import zipfile

import sympy as sp

from T3BCSU2ComponentEquivalence import Evidence, Parser, audit_model

MODELS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
FORMAT = 't3-bc-eft1-component-hermiticity-y-handoff-v1'


def fail(msg: str) -> None:
    raise ValueError(msg)


def numeric_ast(node):
    """Evaluate only small arithmetic Matchete InputForm ASTs; never exec code."""
    kind = node[0]
    if kind == 'num':
        return sp.Integer(node[1])
    if kind == 'sym':
        if node[1] == 'I':
            return sp.I
        fail('unrecognized numeric symbol in tensor array: ' + str(node[1]))
    if kind == 'neg':
        return -numeric_ast(node[1])
    if kind == 'inv':
        return 1 / numeric_ast(node[1])
    if kind == 'add':
        return sp.Add(*(numeric_ast(item) for item in node[1]))
    if kind == 'mul':
        return sp.Mul(*(numeric_ast(item) for item in node[1]))
    if kind == 'pow':
        return numeric_ast(node[1]) ** numeric_ast(node[2])
    if kind == 'list':
        return [numeric_ast(item) for item in node[1]]
    if kind == 'call':
        name, args = node[1:]
        args = [numeric_ast(arg) for arg in args]
        if name == 'Sqrt' and len(args) == 1:
            return sp.sqrt(args[0])
        if name == 'Rational' and len(args) == 2:
            return args[0] / args[1]
        if name == 'Complex' and len(args) == 2:
            return args[0] + sp.I * args[1]
        fail('unrecognized numeric function in tensor array: ' + name)
    fail('unrecognized numeric AST: ' + repr(node)[:120])


def parse_native_arrays(text: str):
    arr = numeric_ast(Parser(text.strip()).parse())
    if not isinstance(arr, list) or len(arr) != 2:
        fail('native array must contain exactly two CG tensors')
    def index4(t, ijkl):
        value = t
        for k in ijkl:
            if not isinstance(value, list) or len(value) != 2:
                fail('native CG tensor must have shape 2x2x2x2')
            value = value[k-1]
        if isinstance(value, list):
            fail('native tensor entries must be scalars')
        return value
    return [{p: index4(t,p) for p in itertools.product((1,2), repeat=4)} for t in arr]


def zero(expr) -> bool:
    return sp.simplify(expr) == 0


def check_actual_native_arrays(native_dir: Path, stem: str):
    """Check all 32 native array entries, not just a Wolfram-recorded boolean."""
    base = native_dir / f'{stem}_nativeCG'
    path = base.with_name(base.name + '_native_cross_cg_arrays.inputform.txt')
    record_path = base.with_name(base.name + '_native_cross_cg_check.json')
    if not path.is_file() or not record_path.is_file():
        fail(f'missing actual native CG arrays/JSON for {stem} at {native_dir}')
    record = json.loads(record_path.read_text(encoding='utf-8'))
    if record.get('format') != 't3-bc-native-cross-cg-exact-component-check-v1' or record.get('model') != MODELS[stem]:
        fail(f'invalid native CG metadata for {stem}')
    t1,t2 = parse_native_arrays(path.read_text(encoding='utf-8'))
    p0 = (1,1,1,1)
    a1 = t1[(1,1,2,2)]
    b1 = t1[p0]-a1
    a2 = t2[(1,1,2,2)]
    b2 = t2[p0]-a2
    if not zero(a1-1) or not zero(b1+1):
        fail(stem + ': CrossInv1 has incorrect coefficients')
    if not zero(a2-b2) or zero(a2):
        fail(stem + ': CrossInv2 must have a nonzero common coefficient')
    for i,j,k,l in itertools.product((1,2),repeat=4):
        p = (i,j,k,l)
        A = int(i==j)*int(k==l)
        B = int(i==k)*int(j==l)
        if not zero(t1[p]-(A-B)):
            fail(stem + ': CrossInv1 component mismatch at ' + str(p))
        if not zero(t2[p]-a2*(A+B)):
            fail(stem + ': CrossInv2 component mismatch at ' + str(p))
    for key,expected in (
        ('CrossInv1_coeff_A_inputform',a1),('CrossInv1_coeff_B_inputform',b1),
        ('CrossInv2_coeff_A_inputform',a2),('CrossInv2_coeff_B_inputform',b2),
    ):
        # Parse symbolic coeff strings safely (e.g. Sqrt[2]/2).
        actual = numeric_ast(Parser(str(record.get(key,''))).parse())
        if not zero(actual-expected):
            fail(stem + ': coefficient disagrees with native array: '+key)
    required = ('native_cross_conditions_met','CrossInv1_in_singlet_span',
                'CrossInv2_in_singlet_span','CrossInv1_exact_required_A_minus_B',
                'CrossInv2_required_equal_A_B_coefficients','CrossInv1_nonzero',
                'CrossInv2_nonzero')
    if any(record.get(k) is not True for k in required):
        fail(stem + ': Wolfram source CG validation did not pass')
    return {
        'native_tensor_file': str(path),
        'native_tensor_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'cross_inv1_coefficients_A_B': [str(sp.simplify(a1)),str(sp.simplify(b1))],
        'cross_inv2_common_coefficient': str(sp.simplify(a2)),
        'native_entries_independently_verified': 32,
        'cross_inv1_equals_A_minus_B': True,
        'cross_inv2_equals_gamma_times_A_plus_B': True,
    }


def inspect_y_evidence(z: zipfile.ZipFile, stem: str):
    p = (f'output/general_t3_mass_eigen_matchete/eft1_scalar_Y_only_canonical_probe/'
         f'{stem}_final_y_only_canonical_match')
    data=json.loads(z.read(p+'.json'))
    recon=json.loads(z.read(p+'_reconstruction_audit.json'))
    herm=json.loads(z.read(p+'_hermiticity_diagnostics.json'))
    if data.get('model') != MODELS[stem]:
        fail('Y matching model mismatch: '+stem)
    if recon.get('independent_matchete_check_passed') is not True or recon.get('independent_hermiticity') is not True:
        fail(stem+': canonical candidate did not pass native Hermiticity check')
    if herm.get('heavy_EFT1_original_hermiticity') is not False or herm.get('light_EFT1_hermiticity') is not False:
        fail(stem+': recorded original Hermiticity failure not found')
    if data.get('tree_matching_completed') is not True or data.get('one_loop_matching_completed') is not True:
        fail(stem+': missing tree/one-loop Y matching')
    if data.get('raw_scalar_llhh_term_count') != 132:
        fail(stem+': unexpected raw Y LLHH inventory')
    if data.get('hard_pole_counterterms_combined') is not False or data.get('eft1_running_performed') is not False:
        fail(stem+': archived renormalization status changed; re-audit before conclusions')
    return {
        'original_matchete_hermiticity_check': False,
        'independent_canonical_matchete_hermiticity_check': True,
        'original_canonical_check_passing': bool(herm.get('heavy_EFT1_original_canonical')),
        'canonical_reconstruction_preserves_saved_majorana_terms': True,
        'raw_Y_LLHH_terms': 132,
        'Y_tree_and_one_loop_executed': True,
        'Y_UV_EFT_counterterms_combined': False,
        'EFT1_running_performed': False,
        'renormalized_Y_C5_available': False,
        'Y_input_json_member': p+'.json',
    }


def run(args):
    ev=Evidence(args.residual_archive)
    summary=json.loads(ev.text('summary.json'))
    if summary.get('format') != 't3-bc-complete-eft1-residual-audit-v1':
        fail('residual ZIP summary format changed')
    with zipfile.ZipFile(args.matching_archive) as matching:
        results=[]
        for stem,model in MODELS.items():
            observed=next((m for m in summary['models'] if m['model']==model),None)
            if observed is None: fail('missing '+model+' residual summary')
            pre=f'{stem}_full_eft1_residual_'
            for k,f in [('saved_eomeft','saved_eomeft_full.inputform.txt'),
                        ('reconstructed_candidate','reconstructed_candidate_full.inputform.txt'),
                        ('symbolic_residual','full_symbolic_residual.inputform.txt')]:
                if ev.sha(pre+f)!=observed['full_expression_sha256'][k]:
                    fail(stem+': source residual hash mismatch: '+k)
            manifest=json.loads(ev.text(stem+'_full_eft1_residual_full_residual_manifest.json'))
            if manifest.get('contains_majorana_residual') is not False:
                fail(stem+': reconstructed candidate does not preserve the Majorana sector')
            sectors,solutions=audit_model(ev,stem)
            if len(sectors)!=25 or sum(q['raw_terms'] for q in sectors)!=137:
                fail(stem+': unexpected sector inventory')
            other=[v for k,v in solutions.items() if k not in ('T3MENativeCrossCGInv1','T3MENativeCrossCGInv2')]
            if len(other)!=7 or any(v['matrix_rank']!=0 or v['augmented_rank']!=0 for v in other):
                fail(stem+': failed noncross Bose symmetry cancellation')
            a=solutions['T3MENativeCrossCGInv1']
            b=solutions['T3MENativeCrossCGInv2']
            if a['required_coefficients']!=['1','-1'] or a['status']!='UNIQUE' or b['status']!='NONUNIQUE':
                fail(stem+': mixed-scalar tensor conditions inconsistent')
            native=check_actual_native_arrays(Path(args.native_dir),stem)
            y=inspect_y_evidence(matching,stem)
            results.append({
                'model':model,
                'source_residual_sha256':observed['full_expression_sha256']['symbolic_residual'],
                'full_residual_terms':137,
                'structural_sectors':25,
                'component_residual_zero':True,
                'original_vs_canonical_component_equivalence_established':True,
                'original_hermiticity_by_component_equivalence':True,
                'native_matchete_original_hermiticity_check_still_false':True,
                'representation_discrepancy_interpretation':'Consistent with native CG/index simplifier not identifying Bose-equivalent expressions; source mechanism not independently isolated',
                'native':native,
                'y_evidence':y,
                'physical_Y_C5':None,
                'physical_Y_matching_certified':False,
            })
        if results[0]['source_residual_sha256']!=results[1]['source_residual_sha256']:
            fail('B/C residual hashes differ')
    output={
        'format':FORMAT,
        'status':'EFT1ComponentHermiticityEstablishedYRenormalizationOutstanding',
        'scope':'All 25 renormalizable residual sectors plus exact inherited Majorana operators',
        'assumptions':['Matchete source fields interpreted as commuting complex scalar components',
                       'SU(2) CG tensors as exported from the original DefineCG registration',
                       'Identical-boson exchange symmetry and standard complex conjugation',
                       'Canonical reconstruction built using the same physical scalar basis'],
        'not_a_matchete_CheckLagrangian_pass_for_saved_eomeft':True,
        'source_replacement_performed':False,
        'renormalized_Y_amplitude_subtraction_performed':False,
        'MSbar_counterterm_subtraction_performed':False,
        'EFT1_running_performed':False,
        'physical_C5':None,
        'models':results,
    }
    p=Path(args.output)
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    for m in results:
        print(m['model']+': 25/25 scalar residual sectors cancel; original component Hermiticity=PROVEN; Matchete syntactic check=FALSE; Y finite C5=NOT CERTIFIED')
    print('PASS: native arrays, B/C component equivalence, original-vs-canonical Hermiticity handoff; report='+str(p))


def self_test():
    import tempfile
    zero_tensor={(i,j,k,l): int(i==j)*int(k==l)-int(i==k)*int(j==l)
        for i,j,k,l in itertools.product((1,2),repeat=4)}
    def arr(t): return '{'+','.join(arr(t[k]) for k in range(2))+'}' if isinstance(t,list) else str(t)
    def make(t):
        return [[[[t[(i,j,k,l)] for l in (1,2)] for k in (1,2)]
                 for j in (1,2)] for i in (1,2)]
    t2={(i,j,k,l): 2*(int(i==j)*int(k==l)+int(i==k)*int(j==l))
        for i,j,k,l in itertools.product((1,2),repeat=4)}
    text=arr([make(zero_tensor),make(t2)])
    x,y=parse_native_arrays(text)
    assert all(zero(x[p]-zero_tensor[p]) and zero(y[p]-t2[p]) for p in zero_tensor)
    for expr,expected in [('Sqrt[2]^2',sp.Integer(2)),('Rational[3,2]',sp.Rational(3,2)),('Complex[0,1]',sp.I)]:
        assert zero(numeric_ast(Parser(expr).parse())-expected)
    # Independent rejection test: a forged true JSON flag cannot hide a
    # corrupted native array. Only the 32 component entries are authoritative.
    with tempfile.TemporaryDirectory() as root:
        nd = Path(root)
        array_path = nd / 'B_m1_nativeCG_native_cross_cg_arrays.inputform.txt'
        meta_path = nd / 'B_m1_nativeCG_native_cross_cg_check.json'
        array_path.write_text(text, encoding='utf-8')
        rec = {
          'format':'t3-bc-native-cross-cg-exact-component-check-v1',
          'model':'B_-1', 'native_cross_conditions_met':True,
          'CrossInv1_in_singlet_span':True,
          'CrossInv2_in_singlet_span':True,
          'CrossInv1_exact_required_A_minus_B':True,
          'CrossInv2_required_equal_A_B_coefficients':True,
          'CrossInv1_nonzero':True,'CrossInv2_nonzero':True,
          'CrossInv1_coeff_A_inputform':'1',
          'CrossInv1_coeff_B_inputform':'-1',
          'CrossInv2_coeff_A_inputform':'2',
          'CrossInv2_coeff_B_inputform':'2'
        }
        meta_path.write_text(json.dumps(rec),encoding='utf-8')
        assert check_actual_native_arrays(nd,'B_m1')['native_entries_independently_verified']==32
        # Modify one exact tensor component while keeping the JSON success flag.
        bad = [make(zero_tensor),make(t2)]
        bad[0][0][0][0][0] = 99
        array_path.write_text(arr(bad), encoding='utf-8')
        try:
            check_actual_native_arrays(nd,'B_m1')
        except ValueError:
            pass
        else:
            raise AssertionError('corrupted native tensor was accepted')
    print('PASS: safe parser, exact 32-component checks, native JSON spoof rejection')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--residual-archive',type=Path)
    ap.add_argument('--native-dir',type=Path)
    ap.add_argument('--matching-archive',type=Path)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    try:
        if args.self_test:self_test();return
        for k in ('residual_archive','native_dir','matching_archive','output'):
            if getattr(args,k) is None:ap.error('--'+k.replace('_','-')+' required')
        run(args)
    except (ValueError,KeyError,IndexError,FileNotFoundError,zipfile.BadZipFile) as e:
        print('FAIL (Y matching uncertified): '+str(e),file=sys.stderr)
        sys.exit(1)


if __name__=='__main__':
    main()
