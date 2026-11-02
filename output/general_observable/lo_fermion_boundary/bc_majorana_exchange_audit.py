"""Source-bound, conditional signed Majorana p=0 test for T3 B_-1/C_-1.

Project path: output/general_observable/lo_fermion_boundary/bc_majorana_exchange_audit.py

The exact Grassmann UV saddle uses DECLARED epsilon_12=+1 and mass -M/2 F^T epsilon F.
The native Matchete sources are independently checked term-by-term; this is NOT
an independently computed native B/C Feynman-rule four-point amplitude or a
Matchete GammaCC -> master Weyl phase certificate. No production files changed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import re
from pathlib import Path
import sys

import build_boundary as bound
from signed_zero_momentum_check import G, bil, integrate_out

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RELEASE = HERE / 'releases' / 'source_021dcc1e8fe8ba0f'
REPORT_FORMAT = 't3-bc-m1-conditional-majorana-p0-native-source-audit-v1'

# Fields in the SOURCE, not an assumed universal operator decomposition.
# (index, y1/y2 pair, U spurions (unique names), scalar eigenfields,
#  explicit leading sign, denominator, extra numerator sqrt(3), CG family)
# Denominator represents only the literal native expression. In particular,
# C_-1 #23/24 contain sqrt(3)/(2*M), NOT a bare 1/(2*M).
ROWS = {
    'B_-1': [
        (11,'11','11 12','XY', 1,'M',False,'Y'),
        (12,'11','11','XX', 1,'2M',False,'Y'),
        (13,'11','12','YY', 1,'2M',False,'Y'),
        (14,'12','11 21','XX', 1,'M',False,'Y'),
        (15,'12','12 21','XY', 1,'M',False,'Y'),
        (16,'12','12 22','YY', 1,'M',False,'Y'),
        (17,'22','21 22','XY', 1,'M',False,'Y'),
        (18,'12','11 22','XY', 1,'M',False,'Y'),
        (19,'22','21','XX', 1,'2M',False,'Y'),
        (20,'22','22','YY', 1,'2M',False,'Y'),
    ],
    'C_-1': [
        (13,'11','11','XX', 1,'2sqrt3M',False,'Q'),
        (14,'12','11 21','XX',-1,'sqrt3M',False,'Q'),
        (15,'22','21','XX', 1,'2sqrt3M',False,'Q'),
        (16,'11','12','YY', 1,'2sqrt3M',False,'Q'),
        (17,'12','12 22','YY',-1,'sqrt3M',False,'Q'),
        (18,'22','22','YY', 1,'2sqrt3M',False,'Q'),
        (19,'11','11 12','XY', 1,'sqrt3M',False,'Q'),
        (20,'12','12 21','XY',-1,'sqrt3M',False,'Q'),
        (21,'12','11 22','XY',-1,'sqrt3M',False,'Q'),
        (22,'22','21 22','XY', 1,'sqrt3M',False,'Q'),
        (23,'12','12 21','XY', 1,'2M',True,'Bridge'),
        (24,'12','11 22','XY',-1,'2M',True,'Bridge'),
    ]
}
DENOMS = {
    'M':'/Coupling[MF, {}, 0]',
    '2M':'/(2*Coupling[MF, {}, 0])',
    'sqrt3M':'/(Sqrt[3]*Coupling[MF, {}, 0])',
    '2sqrt3M':'/(2*Sqrt[3]*Coupling[MF, {}, 0])',
}
CG_MARKER = {
    'Y': ('CG[Bar[T3MENativeY',),
    'Q': ('CG[Bar[T3ExtraCGQ020020I1]',),
    'Bridge': ('CG[T3MEBridgeBB', 'CG[Bar[eps[SU2L]]'),
}


def exact_majorana_saddle(mass=Q(11), weights=(Q(2), Q(-3), Q(5), Q(-7))):
    """Exact LEFT-Grassmann Euler saddle with 4 independent source currents.

    Each current is a two-component external Weyl spinor; its complex-scalar
    factor and SU(2) CG are held OUTSIDE this spinor-only calculation.
    """
    mass = Q(mass)
    weights = tuple(map(Q, weights))
    if mass <= 0 or len(weights) != 4 or not all(weights):
        raise ValueError('Need nonzero mass and four nonzero independent weights')
    external = [[G.var(f'{j}a'), G.var(f'{j}b')] for j in 'abcd']
    heavy = [G.var('ha'), G.var('hb')]
    uv = -mass * bil(heavy, heavy) / 2
    for w, field in zip(weights, external):
        uv = uv + w * bil(heavy, field)
    effective, solution = integrate_out(uv, ('ha', 'hb'))
    diag = sum((weights[i]**2 * bil(external[i],external[i])/(2*mass)
                for i in range(4)), G())
    cross = sum((weights[i]*weights[j]*bil(external[i],external[j])/mass
                 for i in range(4) for j in range(i+1,4)), G())
    if effective != diag + cross:
        raise AssertionError('Independent Grassmann saddle differs from full signed contact')
    if any((uv.left_derivative(x).substitute(solution) != G()) for x in ('ha','hb')):
        raise AssertionError('Unsolved heavy Grassmann equation')
    if bil(external[0],external[1]) != bil(external[1],external[0]):
        raise AssertionError('Weyl bilinear exchange is not symmetric')
    # Verify each of 4 diagonal and 6 unordered mixed monomial coefficients.
    for i in range(4):
        key = tuple(sorted((f'{"abcd"[i]}a',f'{"abcd"[i]}b')))
        if effective.c.get(key) != weights[i]**2/mass:
            raise AssertionError('Diagonal coefficient wrong (bil(aa)=2*a1*a2)')
    for i in range(4):
        for j in range(i+1,4):
            key = tuple(sorted((f'{"abcd"[i]}a',f'{"abcd"[j]}b')))
            if effective.c.get(key) != weights[i]*weights[j]/mass:
                raise AssertionError(f'Mixed coefficient wrong: {i},{j}')
    return {
        'zero_momentum_saddle_verified': True,
        'explicit_spinor_convention': 'epsilon12=+1; Lmass=-(M/2)h^T epsilon h; Lint=sum y_a h^T epsilon l_a',
        'mass': str(mass), 'weights': list(map(str, weights)),
        'independent_diagonal_channels': 4,
        'independent_cross_channels': 6,
        'diagonal_lagrangian_factor': '1/(2M)',
        'cross_lagrangian_factor': '1/M',
        'identical_spinor_bilinear_yields_two_components': True,
        'native_GammaCC_to_master_Weyl_phase_certified': False,
    }


def verify_native_term(term, model, row):
    idx, ys, spur, scalars, sign, denom, numer_sqrt, cg = row
    if term.get('source_index') != idx or term.get('chirality') != 'PL':
        raise ValueError(f'{model} #{idx}: wrong index/chirality')
    raw = term.get('source_term_inputform')
    if not isinstance(raw,str) or not raw:
        raise ValueError(f'{model} #{idx}: source text missing')
    if bound.chirality(raw) != 'PL' or 'GammaCC, Proj[-1]' not in raw:
        raise ValueError(f'{model} #{idx}: wrong native spinor structure')
    if [x['name'] for x in term.get('factors',())] != [f'y{n}' for n in ys]:
        raise ValueError(f'{model} #{idx}: Yukawa factors changed')
    # Regex returns full T3MEUxx; compare ordered unique spurions. A square
    # contains one coupling symbol with power 2, not two text instances.
    expected_spurs = [f'T3MEU{n}' for n in spur.split()]
    native_spurs = re.findall(r'Bar\[Coupling\[(T3MEU\d\d),',raw)
    if native_spurs != expected_spurs or term.get('spurion_names') != expected_spurs:
        raise ValueError(f'{model} #{idx}: scalar mixing spurions changed')
    got_scalars = re.findall(r'Bar\[Field\[(NewScalar[XY]),\s*Scalar,',raw)
    want_scalars = [f'NewScalar{c}' for c in scalars]
    if sorted(got_scalars) != sorted(want_scalars):
        raise ValueError(f'{model} #{idx}: X/Y mass-eigenstate fields changed')
    if len(got_scalars) != 2 or 'NewScalar1' in raw or 'NewScalar2' in raw:
        raise ValueError(f'{model} #{idx}: stale original scalars or wrong leg count')
    if (raw.lstrip().startswith('-')) != (sign < 0):
        raise ValueError(f'{model} #{idx}: explicit native sign changed')
    has_numer = bool(re.match(r'^\s*-?\s*\(Sqrt\[3\]\*',raw))
    if has_numer != numer_sqrt or not raw.rstrip().endswith(DENOMS[denom]):
        raise ValueError(f'{model} #{idx}: denominator/sqrt(3) prefactor changed')
    if not all(x in raw for x in CG_MARKER[cg]):
        raise ValueError(f'{model} #{idx}: expected CG family absent')
    if cg == 'Bridge' and 'T3ExtraCGQ020020I1' in raw:
        raise ValueError(f'{model} #{idx}: bridge CG conflated with quartic CG')
    if cg == 'Q' and 'T3MEBridgeBB' in raw:
        raise ValueError(f'{model} #{idx}: quartic CG conflated with bridge')
    if not all(f'Bar[Coupling[y{n},' in raw for n in set(ys)):
        raise ValueError(f'{model} #{idx}: holomorphic barred Yukawa coupling absent')
    return {
        'source_index':idx,'yukawa_channels':[f'y{n}' for n in ys],
        'mixing_spurions':expected_spurs,'scalar_eigenfields':want_scalars,
        'explicit_sign':sign,'native_denominator':denom,
        'numerator_sqrt3':numer_sqrt,'cg_family':cg,
        'source_expression_sha256':hashlib.sha256(raw.encode('utf8')).hexdigest(),
        'spinor_factor_vs_su2_cg_factor_not_conflated':True,
    }


def source_audit(root=ROOT, release=RELEASE):
    root, release = Path(root).resolve(),Path(release).resolve()
    base = root/'output'/'general_observable'
    ordered_root=base/'ordered_wilson_tensor_v4'
    poly_root=base/'general_wilson_polynomial_batch_v2'
    audit_root=base/'general_wilson_normalization_audit_v1'
    summary = bound.read(release/'summary.json')
    if summary.get('n_models_verified') != 16:
        raise ValueError('Source-resolved 16-model release is missing/incomplete')
    models={}
    for model,rows in ROWS.items():
        saved=bound.read(release/bound.tag(model)/'boundary.json')
        built,_=bound.compile_model(model,root,ordered_root,poly_root,audit_root)
        if saved['checksums_sha256'] != built['checksums_sha256']:
            raise ValueError(f'{model}: source chain differs from prior release')
        ordered=bound.read(ordered_root/bound.tag(model)/'ordered_tensor.json')
        native_path=bound.source_file(model,root,ordered_root,ordered)
        native=bound.read(native_path)
        if bound.digest(native_path) != saved['checksums_sha256']['native_source']:
            raise ValueError(f'{model}: native source digest mismatch')
        if native.get('model') != model or native.get('format') != 't3-mixed-factorised-LLSS-source-v1':
            raise ValueError(f'{model}: wrong native source type')
        terms=native.get('operator_terms',[])
        pl={t['source_index']:t for t in terms if bound.chirality(t['source_term_inputform'])=='PL'}
        pr=[t for t in terms if bound.chirality(t['source_term_inputform'])=='PR']
        if len(pl)!=len(rows) or len(pr)!=len(rows) or set(pl)!={r[0] for r in rows}:
            raise ValueError(f'{model}: unexpected PL/PR source inventory')
        if len(native.get('heavy_fermion_masses_gev',[])) != 3 or not all(
           isinstance(m,(int,float)) and m>0 for m in native['heavy_fermion_masses_gev']):
            raise ValueError(f'{model}: invalid diagonal Majorana mass eigenvalues')
        items=[verify_native_term(pl[r[0]],model,r) for r in rows]
        models[model]={
            'native_source_sha256':bound.digest(native_path),
            'release_boundary_sha256':bound.digest(release/bound.tag(model)/'boundary.json'),
            'heavy_mass_eigenvalues_gev':native['heavy_fermion_masses_gev'],
            'pl_source_terms_checked':len(items),'pr_terms_excluded':len(pr),
            'native_source_terms':items,
            'native_cg_tensor_contraction_independently_evaluated':False,
            'connected_native_mixed_eigenfield_feynman_rule_amplitude_evaluated':False,
        }
    return {'release_summary_sha256':bound.digest(release/'summary.json'),'models':models}


def run(args):
    grassmann=exact_majorana_saddle()
    if args.offline:
        print('PASS: exact signed 4-current Majorana Grassmann p=0 saddle (4 diagonal + 6 cross)')
        print('CONDITIONAL ONLY: source/CG/Matchete phase not checked in offline mode')
        return 0
    source=source_audit(args.repo_root,args.release or RELEASE)
    result={
        'format':REPORT_FORMAT,
        'status':'ConditionalMajoranaExchangeAndAllNativePLSourceSignaturesVerified_AmplitudeCGAndMasterPhaseOpen',
        'grassmann_benchmark':grassmann,
        'native_source_audit':source,
        'four_point_native_BC_feynman_rule_exchange_computed':False,
        'native_mixed_scalar_su2_cg_contraction_computed':False,
        'matchete_gamma_cc_to_master_weyl_phase_certified':False,
        'physical_16_model_boundary_certified':False,
        'previous_release_modified':False,
        'scope_note': ('The 1/2 vs 1 factors are independently derived for a generic four-current Majorana source. '
                       'The 22 actual Matchete source structures are authenticated against the release; '
                       'especially C_-1 #23/24 have a separate sqrt(3)/2 bridge CG structure, '
                       'so native CG tensor products and full Feynman amplitudes remain open.'),
    }
    payload=json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+'\n'
    digest=hashlib.sha256(payload.encode('utf8')).hexdigest()[:16]
    print('PASS: B_-1 (10 PL + 10 PR) and C_-1 (12 PL + 12 PR) pinned source signatures')
    print('PASS: exact conditional Majorana 4-current saddle, 4 diagonal / 6 mixed spinor channels')
    print('OPEN: independent mixed eigenfield Feynman-rule amplitudes, SU(2) CG products, Weyl phase')
    if args.check_only:
        print('CHECK ONLY: no output written')
        return 0
    dest=Path(args.output_dir or HERE/'normalization_probes').resolve()
    if dest == Path(args.release or RELEASE).resolve() or Path(args.release or RELEASE).resolve() in dest.parents:
        raise ValueError('Refusing to write inside the certified release')
    dest.mkdir(parents=True,exist_ok=True)
    output=dest/f'bc_majorana_exchange_{digest}.json'
    if output.exists():
        if output.read_text(encoding='utf8') != payload:
            raise ValueError('Deterministic report exists but contains different content')
        print('UNCHANGED: existing identical report:',output)
    else:
        with output.open('x',encoding='utf8') as fp:fp.write(payload)
        print('CREATED:',output)
    return 0


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--offline',action='store_true')
    ap.add_argument('--check-only',action='store_true')
    ap.add_argument('--repo-root',type=Path,default=ROOT)
    ap.add_argument('--release',type=Path,default=None)
    ap.add_argument('--output-dir',type=Path,default=None)
    args=ap.parse_args(argv)
    try:return run(args)
    except (AssertionError,KeyError,ValueError,TypeError,OSError) as exc:
        print(f'BLOCKED: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2


if __name__=='__main__':raise SystemExit(main())
