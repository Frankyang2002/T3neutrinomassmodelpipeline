#!/usr/bin/env python3
"""B(-1)/C(-1) conditional Y scalar-bubble -> Warsaw/Weinberg interface.

This is NOT complete physical EFT threshold matching. It composes *existing*
independent v5/v6 scalar-bubble subtraction, v8 flavour and v11.2 live
Matchete operator identity. It deliberately does not compute missing full
EFT counterterms, inherited hard corrections, running, or X matching.

At strict dimension 5 and muY=MY, the audited formal scalar-bubble sector is
A_PL = [1/(16*pi**2*MF)] sum_{a<=b} k_ab conj(y_a)_i conj(y_b)_j.
The installed Warsaw operator has no explicit 1/2 and v11.2 shows the
native PL operator equals it (positive sign); consequently the conventional
1/2-normalised symmetric coefficient is C5=A_PL+A_PL^T. Do NOT sum PR again.

Usage:
  python scripts/BuildT3BCYFiniteBubbleInterface.experimental.py --self-test
  python scripts/BuildT3BCYFiniteBubbleInterface.experimental.py --output output/audits/t3_bc_y_conditional_finite_interface_v1/summary.json
  python scripts/BuildT3BCYFiniteBubbleInterface.experimental.py --numeric-input path/to/inputs.json --output output/.../numeric.json
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
DEFAULT_AUDIT = Path('output/audits')
QUARTICS = {'lambdaT3', 'lambdaT3ExtraQ020020I1', 'lambdaT3ExtraQ020200I1'}

class MatchingGateError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise MatchingGateError(message)


def load_json(path: Path, fmt: str | None = None) -> dict:
    require(path.is_file(), 'Required local evidence missing: ' + str(path))
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(data, dict), 'Invalid JSON object: ' + str(path))
    if fmt is not None:
        require(data.get('format') == fmt,
                'Wrong evidence schema at ' + str(path) + ': ' + str(data.get('format')))
    return data


def file_sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def get_v8_helper(root: Path):
    """Reuse project-owned strict ledger parser, rather than reinterpret it."""
    script = root / 'scripts' / 'AnalyzeT3BCYFlavorSymmetry.experimental.py'
    require(script.is_file(), 'Required existing v8 script missing: ' + str(script))
    spec = importlib.util.spec_from_file_location('t3bc_v8_conditional_interface', script)
    require(spec is not None and spec.loader is not None, 'Cannot load v8 parser')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def polynomial_coefficients(ledger: dict, chirality: str):
    """Directly retain native PL/PR symbol conjugation and ordered y1y2 slot.

    Returns three SymPy polynomials with hbar/MF FACTORED OUT. The source has
    no y2y1 record. The second term is not silently doubled.
    """
    import sympy as sp
    require(chirality in {'PL','PR'}, 'Unknown chirality')
    out = defaultdict(lambda: sp.Integer(0))
    for mon in ledger['monomials']:
        if mon['chirality'] != chirality:
            continue
        pair = (mon['yukawa_i'], mon['yukawa_j'])
        require(pair in {('y1','y1'), ('y1','y2'), ('y2','y2')},
                'Unexpected source Yukawa order')
        require(mon['quartic'] in QUARTICS, 'Unexpected quartic source')
        value = sp.Rational(mon['coefficient_numerator'],mon['coefficient_denominator'])
        value *= sp.sqrt(3)**mon['sqrt3_exponent']
        value *= sp.Symbol(('bar_' if mon['quartic_conjugated'] else '')+mon['quartic'])
        for u in mon['mixing']:
            value *= sp.Symbol(('bar_' if u['conjugated'] else '')+u['name'])**u['exponent']
        out[pair] += value
    require(all(sp.expand(out[p]) != 0 for p in [('y1','y1'),('y1','y2'),('y2','y2')]),
            'Missing one of three Yukawa structures')
    return tuple(sp.expand(out[p]) for p in [('y1','y1'),('y1','y2'),('y2','y2')])


def polynomial_conjugate(expr):
    """Exact formal involution on native barred/unbarred coupling symbols."""
    import sympy as sp
    repl = {}
    for symbol in expr.free_symbols:
        nm = symbol.name
        repl[symbol] = sp.Symbol(nm[4:] if nm.startswith('bar_') else 'bar_'+nm)
    return sp.expand(expr.xreplace(repl))


def check_symmetrization(A):
    """No Hermitian conjugate here: flavour symmetry uses TRANSPOSE."""
    C5 = A + A.T
    require(C5.shape[0] == C5.shape[1], 'Expected square ordered coefficient')
    require((C5 - C5.T).is_zero_matrix, 'C5 is not symmetric')
    return C5


def gates_for_model(stem: str, audit_dir: Path, v5: dict, v6: dict,
                    v8: dict, v11: dict, helper):
    import sympy as sp
    model = STEMS[stem]
    v5item = v5.get('model_results', {}).get(model)
    v6item = v6.get('models', {}).get(model)
    v8item = v8.get('models', {}).get(stem)
    require(all(isinstance(x,dict) for x in (v5item, v6item, v8item)),
            model+': incomplete v5/v6/v8 evidence')
    v5flags = ('tree_WQ_kernel_equals_native_Y_pole',
               'independent_bubble_MSbar_CT_cancels_native_Y_pole',
               'Y_formal_finite_equals_XY_kernel_at_muY',
               'UV_IR_separated_for_YY_XY_scalar_bubbles')
    v6flags = ('exact_LLXX_component_equality','exact_HHXX_component_equality',
               'exact_XX_loop_UV_component_equality',
               'EFT1_minus_EFT2_scalar_bubble_CT_equals_minus_Y',
               'EFT1_EFT2_scalar_bubble_CT_difference_cancels_native_Y_pole',
               'XX_finite_MSbar_log_cancelled_between_EFTs',
               'previous_v5_CT_matches_v6_CT_difference')
    for nm in v5flags:
        require(v5item.get(nm) is True, model+': missing/failed v5 gate '+nm)
    for nm in v6flags:
        require(v6item.get(nm) is True, model+': missing/failed v6 gate '+nm)
    require(v5item.get('nonzero_component_monomials', {}).get('YY') == 18 and
            v5item.get('nonzero_component_monomials', {}).get('XY') == 32,
            model+': YY/XY support changed')
    require(v5item.get('physical_matching_certified') is False and
            v6item.get('physical_matching_certified') is False,
            model+': unexpected changed certification flags; inspect source')
    finite_file = audit_dir/'t3_bc_y_renormalization_ledger_v1'/f'{stem}_y_llhh_constant_neutral_flavour.json'
    ledger = load_json(finite_file)
    require(ledger.get('model') == model, model+': source ledger model mismatch')
    require(v8item.get('ledgers', {}).get('constant', {}).get('sha256') == file_sha(finite_file),
            model+': v8-validated source ledger hash changed')
    parsed = helper.parse_ledger(ledger, stem, 'constant')
    require(parsed['conjugacy'] is True, model+': source H.c. mismatch')
    require(v11['native'][stem]['Warsaw_chirality_and_H_orientations_verified'] is True,
            model+': v11 PL/PR orientation missing')
    pl = polynomial_coefficients(ledger, 'PL')
    pr = polynomial_coefficients(ledger, 'PR')
    for pl_expr, pr_expr in zip(pl,pr):
        require(sp.expand(polynomial_conjugate(pr_expr)-pl_expr)==0,
                model+': exact PL/PR polynomial mismatch')
    # Cross-check directly against the v8 independent PR polynomial map.
    for pr_ours, pr_v8 in zip(pr,parsed['coefficients']):
        require(sp.expand(pr_ours-pr_v8)==0, model+': differing flavour parser result')
    # The source terms are ORDERED in Yukawa space: A12 exists, A21=0.
    a11,a12,a22 = pl
    A = sp.Matrix([[a11,a12],[0,a22]])
    C5 = check_symmetrization(A)
    return {
        'model': model,
        'source_ledger': str(finite_file), 'source_sha256': file_sha(finite_file),
        'status': 'CONDITIONAL_RESTRICTED_SCALAR_BUBBLE_ONLY',
        'PR_PL_exact_source_conjugacy': True,
        'native_PL_equals_Warsaw_positive_sign': True,
        'scalar_bubble_YY_XY_MSbar_and_XX_overlap_verified': True,
        'native_mass_and_loop_factor': 'hbar/MF, hbar=1/(16*pi**2)',
        'Warsaw_A_ordered_2x2_coefficient_without_hbar_over_MF':
            [[str(A[i,j]) for j in range(2)] for i in range(2)],
        'conditional_conventional_C5_2x2_coefficient_without_hbar_over_MF':
            [[str(C5[i,j]) for j in range(2)] for i in range(2)],
        'Yukawa_matrix_orientation': 'L=[conj(y1),conj(y2)] (3x2), C5=(hbar/MF)*L@C5_2x2@L.T',
        'includes_Y_bubble_hard_finite_dimension5_at_muY': True,
        'soft_power_corrections_included': False,
        'full_operator_counterterms_certified': False,
        'inherited_F_one_loop_and_EFT_running_included': False,
        'physical_C5': None,
        'physical_C5_certified': False,
    }


def validate_all(root: Path, audit_dir: Path):
    v5path=audit_dir/'t3_bc_y_scalar_bubble_msbar_v1'/'summary.json'
    v6path=audit_dir/'t3_bc_y_eft2_light_xx_overlap_v1'/'summary.json'
    v8path=audit_dir/'t3_bc_y_flavour_symmetry_v1'/'summary.json'
    v11path=audit_dir/'t3_bc_y_warsaw_direct_identity_v11'/'summary_live.json'
    v5=load_json(v5path,'t3bc-y-scalar-bubble-msbar-counterterm-audit-v1')
    v6=load_json(v6path,'t3bc-y-e1-e2-light-xx-overlap-report-v1')
    v8=load_json(v8path,'t3-bc-y-ordered-flavour-factorization-v1')
    v11=load_json(v11path,'t3-bc-warsaw-spinor-orientation-v11')
    require(v11.get('absolute_installed_Matchete_operator_factor_certified') is True and
            v11.get('live_Matchete_operator_identity',{}).get('direct_equivalence_status')=='ExactPlusEquality' and
            v11.get('native_spinor_relative_sign')=='+1',
            'Live Matchete v11.2 positive-sign operator equivalence missing')
    require(v11.get('physical_C5_certified') is False and
            v5.get('physical_C5_Y') is None and v6.get('physical_C5_Y') is None,
            'Pre-existing reports changed physical certification assumptions')
    helper = get_v8_helper(root)
    models={stem:gates_for_model(stem,audit_dir,v5,v6,v8,v11,helper) for stem in STEMS}
    sources={str(path):file_sha(path) for path in (v5path,v6path,v8path,v11path)}
    return {'format':'t3bc-y-restricted-conditional-finite-warsaw-interface-v1',
            'scope':'strict dim-5 scalar-quartic bubble Y threshold at muY=MY; native PL mapped exactly to Warsaw Q',
            'input_audit_sha256':sources,
            'models':models,
            'normalized_Weinberg_convention':'L5=(1/2) C5_ij Q_ij + h.c.; C5=A+A.T for positive native/Warsaw identity',
            'candidate_is_only_one_loop_scalar_bubble_sector':True,
            'not_included':[
                'complete EFT1/EFT2 field and operator counterterms outside scalar-quartic bubble subset',
                'fermion-hard one-loop boundary and intermediate running',
                'Y-inherited LLXX/LLHX pieces propagated to X and X threshold',
                'dimension-7 and higher mixed-XY soft-region terms',
                'low-scale SM running and neutrino pole-mass matching'],
            'full_physical_matching_certified':False,
            'physical_C5':None,
            'physical_neutrino_mass_matrix':None}


def parse_complex(value, where: str):
    """Explicit [real,imag] pairs: no eval/implicit imaginary sign conventions."""
    require(isinstance(value,list) and len(value)==2 and all(type(z) in (int,float) for z in value),
            'Complex input must be [real,imag] at '+where)
    z=complex(value[0],value[1])
    import math
    require(math.isfinite(z.real) and math.isfinite(z.imag), 'Nonfinite complex input: '+where)
    return z


def numeric_one(stem: str, input_record: dict, ledger: dict, source_params: dict):
    """Callable conditional 3x3 matrix; only if Yukawas and quartics supplied."""
    import numpy as np
    # Executing `python scripts/name.py` does not reliably put the repository
    # root on sys.path; the existing Numerical package lives at that root.
    project_root = Path(__file__).resolve().parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0,str(project_root))
    from Numerical.experimental.T3BCNativeYLLHHFlavorKernel import evaluate_spurion_kernel
    model=STEMS[stem]
    require(input_record.get('model')==model, model+': numeric input model mismatch')
    MF=input_record.get('MF_GeV')
    require(type(MF) in (int,float) and 0<MF<1e15, 'MF_GeV must be positive real')
    ux=source_params['unitary_U']
    U=np.asarray(ux['real'],dtype=float)+1j*np.asarray(ux['imag'],dtype=float)
    require(U.shape==(2,2),'Unexpected benchmark mixing U')
    y=[]
    for nm in ('y1','y2'):
        vals=input_record.get(nm)
        require(isinstance(vals,list) and len(vals)==3, nm+': require three flavour components')
        y.append(np.array([parse_complex(z,nm) for z in vals],dtype=complex))
    lam=input_record.get('quartics',{})
    require(isinstance(lam,dict) and set(lam)==QUARTICS,
            'Exactly the three source quartics must be supplied')
    quartics={nm:parse_complex(v,nm) for nm,v in lam.items()}
    hbar=1.0/(16.0*np.pi**2)
    ms=evaluate_spurion_kernel(ledger,U=U,quartics=quartics,
                              y1=y[0],y2=y[1],MF=MF,hbar=hbar)
    PL,PR=ms['PL'],ms['PR']
    scale=max(1.0,float(np.max(np.abs(PL))),float(np.max(np.abs(PR))))
    require(np.max(np.abs(PR-PL.conjugate()))<=1e-12*scale,
            'PR and PL numerical conjugacy not satisfied')
    C5=PL+PL.T
    require(np.max(np.abs(C5-C5.T))<=1e-12*scale, 'C5 symmetry check failed')
    def pairs(arr):
        return [[[float(v.real),float(v.imag)] for v in row] for row in arr]
    return {'model':model, 'input_MF_GeV':MF,
            'normalization':'1/(16*pi^2*MF); C5=A_PL+A_PL.T',
            'A_PL_ordered_native_and_Warsaw_GeV_inverse':pairs(PL),
            'C5_Y_restricted_scalar_bubble_only_GeV_inverse':pairs(C5),
            'complex_symmetric':True,'exact_source_Hc_consistency_checked':True,
            'is_physical_C5':False,'full_physical_matching_certified':False,
            'physical_neutrino_mass_matrix':None}


def self_test():
    import sympy as sp
    xx,xy,yy=sp.symbols('a b c')
    A=sp.Matrix([[xx,xy],[0,yy]])
    C=check_symmetrization(A)
    require(C==sp.Matrix([[2*xx,xy],[xy,2*yy]]),'Wrong symmetrization/1/2 factor')
    q=sp.Symbol('lambdaT3');bar=sp.Symbol('bar_lambdaT3')
    require(polynomial_conjugate(q+2*bar)==bar+2*q,'Wrong native bar mapping')
    require(polynomial_conjugate(polynomial_conjugate(q+bar))==q+bar,
            'Non-involutive conjugation')
    try:
        parse_complex([1.0, 'i'],'invalid')
    except MatchingGateError: pass
    else: raise MatchingGateError('Bad complex input was accepted')
    try:
        require(False, 'tamper')
    except MatchingGateError: pass
    else: raise MatchingGateError('Gate tamper rejection failed')
    print('PASS: ordered Warsaw A -> C5=A+A.T; mixed Yukawa not double-counted')
    print('PASS: explicit PL/PR conjugation, native-sign and numeric-input rejection')
    print('PASS: no complete physical matching or neutrino mass claimed')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--root',type=Path,default=Path('.'))
    ap.add_argument('--audit-dir',type=Path,default=DEFAULT_AUDIT)
    ap.add_argument('--numeric-input',type=Path,
                    help='Optional JSON with model, MF_GeV, three-component y1/y2 and three complex quartics')
    ap.add_argument('--require-physical',action='store_true',
                    help='Fail by design: incomplete EFT counterterms/running/X matching')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    try:
        if args.self_test: self_test()
        if args.require_physical:
            raise MatchingGateError('Physical C5 is NOT certified: complete EFT1/EFT2 counterterms, running and X matching missing')
        if args.output:
            root=args.root.resolve()
            audit=args.audit_dir if args.audit_dir.is_absolute() else root/args.audit_dir
            report=validate_all(root,audit)
            if args.numeric_input:
                input_path=args.numeric_input if args.numeric_input.is_absolute() else root/args.numeric_input
                provided=load_json(input_path)
                stem=next((s for s,n in STEMS.items() if n==provided.get('model')),None)
                require(stem is not None,'Unsupported numeric model')
                source_params=load_json(root/'output/general_t3_mass_eigen_matchete'/f'{stem}_final_parameters.json')
                require(source_params.get('model')==STEMS[stem], 'Benchmark source parameter mismatch')
                ledger=load_json(audit/'t3_bc_y_renormalization_ledger_v1'/f'{stem}_y_llhh_constant_neutral_flavour.json')
                report['optional_numeric_candidate']=numeric_one(stem,provided,ledger,source_params)
                report['optional_numeric_input_sha256']=file_sha(input_path)
            out=args.output if args.output.is_absolute() else root/args.output
            out.parent.mkdir(parents=True,exist_ok=True)
            out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
            for stem,model in STEMS.items():
                print(f'{model}: v5/v6/v8/v11.2 gates PASS; normalized conditional scalar-bubble C5 constructed')
            print('FULL PHYSICAL C5: NOT CERTIFIED')
            print('Saved:',out)
        elif not args.self_test:
            ap.error('--output and/or --self-test required')
        return 0
    except (MatchingGateError,KeyError,ValueError,TypeError,AttributeError,ImportError,OSError) as exc:
        print('FAIL (not physical C5): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    raise SystemExit(main())
