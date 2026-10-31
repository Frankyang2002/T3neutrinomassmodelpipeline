#!/usr/bin/env python3
"""Workstream C B(-1)/C(-1) F/Y/X fixed-order audit and conditional interfaces.

Source anchored; deliberately does not construct renormalized physical C5.

Batch operations:
  (i) validate certified v11.2/v12/v13 source evidence;
 (ii) check Y/LLXX and X/LLHH raw pole-log inventories against saved SHA256;
(iii) confirm the missing X inheritance and full CT/running gates are explicit;
 (iv) provide separate numeric F-formal and Y-restricted 3x3 C5 tensors on
      explicit user couplings, never their unqualified physical sum.

No Matchete re-run; no production sources, UV inputs or RG equations edited.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
RAW_EXPECTATIONS = {
    ('Y', 'LLXX'): {'B_m1': (666, 240, 240, 186), 'C_m1': (582, 198, 198, 186)},
    ('X', 'LLHH'): {'B_m1': (36, 18, 18, 0), 'C_m1': (36, 18, 18, 0)},
}
MASS_LOGS = {
    'Y': 'Log[μbar2/Coupling[MSY, {}, 0]^2]',
    'X': 'Log[μbar2/Coupling[MSX, {}, 0]^2]',
}
AUDIT = Path('output/audits')


class BatchError(ValueError):
    pass


def check(condition, message):
    if not condition:
        raise BatchError(message)


def get_json(path: Path, fmt=None):
    check(path.is_file(), f'Missing required local evidence: {path}')
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    check(isinstance(data, dict), f'Expected JSON object: {path}')
    if fmt is not None:
        check(data.get('format') == fmt, f'Wrong evidence schema: {path}; {data.get("format")}')
    return data


def hash_file(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def import_file(path: Path, name: str):
    check(path.is_file(), f'Missing prerequisite interface: {path}')
    spec = importlib.util.spec_from_file_location(name, path)
    check(spec is not None and spec.loader is not None, f'Cannot load: {path}')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def import_local_project(root: Path):
    r = str(root.resolve())
    if r in sys.path:
        sys.path.remove(r)
    sys.path.insert(0, r)
    from Numerical.experimental import T3BCPoleLogConsistency as pole
    from Numerical.experimental import T3BCNativeYLLHHFlavorKernel as kernel
    for obj, rel in ((pole, 'Numerical/experimental/T3BCPoleLogConsistency.py'),
                     (kernel, 'Numerical/experimental/T3BCNativeYLLHHFlavorKernel.py')):
        check(Path(obj.__file__).resolve() == (root / rel).resolve(),
              f'Unexpected imported parser source: {obj.__file__}')
    return pole, kernel


def unpack_sectors(raw: str, stage: str, pole_module):
    check(stage in MASS_LOGS, 'Unsupported heavy scalar')
    source_terms = pole_module.terms(raw)
    sectors = {'pole': [], 'log': [], 'constant': []}
    lgname = MASS_LOGS[stage]
    for t in source_terms:
        e = t.count('ϵ')
        l = t.count('Log[')
        check((e, l) in {(0, 0), (0, 1), (1, 0)},
              f'{stage}: mixed or unsupported pole/log structure')
        if e:
            check(t.count('ϵ*') == 1, f'{stage}: unsupported epsilon denominator')
            sectors['pole'].append(t.replace('ϵ*', '', 1))
        elif l:
            check(t.count(lgname) == 1, f'{stage}: unexpected mass log')
            sectors['log'].append(t.replace(lgname, '1', 1))
        else:
            sectors['constant'].append(t)
    return len(source_terms), sectors


def certify_raw_sector(path, *, stage, sector, stem, expected_hash, pole_module):
    check(path.is_file(), f'Missing local raw matching input: {path}')
    observed_hash = hash_file(path)
    check(observed_hash == expected_hash,
          f'{stem}/{stage}/{sector}: source SHA256 differs from audited finite extraction')
    total, sec = unpack_sectors(path.read_text(encoding='utf-8-sig'), stage, pole_module)
    n = (total, len(sec['pole']), len(sec['log']), len(sec['constant']))
    expected = RAW_EXPECTATIONS[(stage, sector)][stem]
    check(n == expected, f'{stem}/{stage}/{sector}: expected {expected}, observed {n}')
    table = {}
    p = pole_module.polynomial(tuple(sec['pole']), table=table)
    l = pole_module.polynomial(tuple(sec['log']), table=table)
    check(p == l, f'{stem}/{stage}/{sector}: exact native pole-log mismatch')
    return {
        'source': str(path), 'source_sha256': observed_hash,
        'raw_terms': n[0], 'pole_terms': n[1], 'log_terms': n[2],
        'scale_independent_constant_terms': n[3],
        'nonzero_opaque_pole_monomials': len(p),
        'exact_opaque_pole_equals_mass_log': True,
        'UV_IR_pole_assignment_certified': False,
        'full_EFT_counterterm_subtraction_certified': False,
        'formal_constant_is_physical_matching': False,
    }


def evidence_for_model(root, stem, extraction, physical_gate, pole_module):
    model = STEMS[stem]
    rows = extraction.get('sector_results', [])
    check(isinstance(rows, list), 'Malformed scalar finite extraction report')
    out = {}
    for stage, sector in RAW_EXPECTATIONS:
        match = [x for x in rows if x.get('model') == model and
                 x.get('threshold') == stage and x.get('sector') == sector]
        check(len(match) == 1, f'{model}: missing or duplicate {stage}/{sector} extraction')
        row = match[0]
        src = root / row['source_file']
        record = certify_raw_sector(src, stage=stage, sector=sector, stem=stem,
                                    expected_hash=row['source_sha256'], pole_module=pole_module)
        check(record['raw_terms'] == row['terms_total'] and
              record['pole_terms'] == row['pole_terms'] and
              record['log_terms'] == row['log_terms'] and
              record['scale_independent_constant_terms'] == row['constant_terms'],
              f'{model}: source differs from reported finite sector inventory')
        out[f'{stage}_{sector}'] = record
    gate = physical_gate.get('models', {}).get(model)
    check(isinstance(gate, dict), f'{model}: physical Y evidence gate missing')
    check(gate.get('source_algebraic_equivalence_proven') is False and
          gate.get('physical_matching_certified') is False,
          f'{model}: changed/corrupted physical Y evidence gate status')
    # Check that the X native source omitted the Y one-loop EFT2 inheritance.
    # Source code of RunT3BCXFromSavedY explicitly sets this False; it must
    # not be inferred from the zero standalone X constant.
    x_y_tree = root / 'output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3' / \
        f'{stem}_final_x_from_saved_y_match_x_followup_y_tree_eomeft.inputform.txt'
    x_y_full = root / 'output/general_t3_mass_eigen_matchete/eft1_scalar_xy_fresh_kernel_v3' / \
        f'{stem}_final_x_from_saved_y_match_x_followup_y_full_eomeft.inputform.txt'
    check(x_y_tree.is_file() and x_y_full.is_file(),
          f'{model}: missing saved Y tree/full EFT2 continuation')
    htree, hfull = hash_file(x_y_tree), hash_file(x_y_full)
    check(htree != hfull, f'{model}: Y tree and full EFT2 continuation unexpectedly equal')
    out['Y_EFT2_inheritance'] = {
        'Y_tree_EFT2_source': str(x_y_tree), 'Y_full_EFT2_source': str(x_y_full),
        'Y_tree_sha256': htree, 'Y_full_sha256': hfull,
        'Y_full_contains_corrections_beyond_tree': True,
        'X_native_match_was_tree_only': True,
        'Y_one_loop_EFT2_in_X_matching': False,
        'requires_consistent_fixed_order_inheritance': True,
        'formal_Y_LLXX_at_one_loop_is_not_automatically_X_one_loop_matching': True,
    }
    return out


def numeric_formal_f(root, stem, record, v12, v13, kernel):
    """Use the exact local native F term parser, not a fabricated coefficient."""
    import numpy as np
    model = STEMS[stem]
    check(record.get('model') == model, 'F numeric input model mismatch')
    mF = record.get('MF_GeV')
    check(type(mF) in (float, int) and np.isfinite(mF) and 0 < mF < 1e15,
          'MF_GeV must be a finite positive number')
    params = get_json(root/'output/general_t3_mass_eigen_matchete'/f'{stem}_final_parameters.json')
    check(params.get('model') == model, 'Mass-eigen benchmark model mismatch')
    source_path = root/'output/general_t3_mass_eigen_matchete'/'eft1_fermion_first_oneloop'/\
        f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
    sec = v13.split_sectors(source_path.read_text(encoding='utf-8-sig'))
    led, _ = v13.make_ledger(sec['constant'], model)
    U0 = params['unitary_U']
    U = np.asarray(U0['real'], dtype=float) + 1j*np.asarray(U0['imag'], dtype=float)
    check(U.shape == (2, 2) and np.linalg.norm(U.conj().T@U-np.eye(2)) < 1e-9,
          'Mass-eigen U mixing is not unitary')
    q = record.get('quartics')
    check(isinstance(q, dict) and set(q) == v12.QUARTICS,
          'Must supply exactly three source quartics')
    y1 = record.get('y1')
    y2 = record.get('y2')
    check(isinstance(y1, list) and isinstance(y2, list) and len(y1) == 3 and len(y2) == 3,
          'Must provide three complex Yukawas per vector')
    y = [np.array([v12.parse_complex(z, f'y{k}[{i}]') for i,z in enumerate(vals)],
                  dtype=complex) for k,vals in ((1,y1),(2,y2))]
    quartics = {key: v12.parse_complex(q[key], key) for key in v12.QUARTICS}
    matrices = kernel.evaluate_spurion_kernel(
        led, U=U, quartics=quartics, y1=y[0], y2=y[1],
        MF=mF, hbar=1.0/(16.0*np.pi**2))
    PL, PR = matrices['PL'], matrices['PR']
    scale = max(1., float(np.max(np.abs(PL))),float(np.max(np.abs(PR))))
    check(np.max(np.abs(PL.conj()-PR)) < 1e-12*scale,
          'Native formal F PR/PL coefficient conjugacy failed')
    C = PL + PL.T
    check(np.max(np.abs(C-C.T)) < 1e-12*scale, 'F formal C5 tensor is not symmetric')
    return {
      'model':model, 'status':'FORMAL_F_HARD_ONE_LOOP_CONSTANT_ONLY',
      'common_factor':'hbar/MF, hbar=1/(16*pi**2)',
      'A_F_PL_ordered_Warsaw_GeV_inverse': pairs(PL),
      'C5_F_FORMAL_GeV_inverse': pairs(C),
      'PR_PL_complex_conjugacy_verified':True,
      'physical_matching_certified':False,
    }


def pairs(arr):
    return [[[float(z.real),float(z.imag)] for z in row] for row in arr]


def numeric_for_model(root, stem, record, v12, v13, kernel):
    # Return individual matrices, deliberately NEVER their numerical sum.
    model = STEMS[stem]
    params = get_json(root/'output/general_t3_mass_eigen_matchete'/f'{stem}_final_parameters.json')
    yledger = get_json(root/AUDIT/'t3_bc_y_renormalization_ledger_v1'/
                       f'{stem}_y_llhh_constant_neutral_flavour.json')
    ynumeric = v12.numeric_one(stem, record, yledger, params)
    fnumeric = numeric_formal_f(root, stem, record, v12, v13, kernel)
    return {
       'model': model,
       'F_boundary_FORMAL': fnumeric,
       'Y_threshold_RESTRICTED_BUBBLE': ynumeric,
       'F_plus_Y_physical_sum_constructed': False,
       'X_threshold_correction_included': False,
       'SM_running_included': False,
       'physical_C5': None,
    }


def run_batch(root: Path, output: Path, numeric_input: Path | None, numeric_manifest: Path | None):
    pole, kernel = import_local_project(root)
    v12 = import_file(root/'scripts/BuildT3BCYFiniteBubbleInterface.experimental.py', 't3bc_v14_v12')
    v13 = import_file(root/'scripts/BuildT3BCFHardBoundaryInterface.experimental.py', 't3bc_v14_v13')
    audits = root / AUDIT
    y = get_json(audits/'t3_bc_y_conditional_finite_interface_v1/summary.json',
                 't3bc-y-restricted-conditional-finite-warsaw-interface-v1')
    f = get_json(audits/'t3_bc_f_hard_boundary_v1/summary.json',
                 't3-bc-F-hard-formal-boundary-plus-Y-conditional-ledger-v1')
    w = get_json(audits/'t3_bc_y_warsaw_direct_identity_v11/summary_live.json',
                 't3-bc-warsaw-spinor-orientation-v11')
    extraction_path = audits/'t3_bc_scalar_finite_term_extraction_v1/summary.json'
    physical_gate_path = audits/'t3_bc_physical_y_matching_gate_v1/summary.json'
    extraction = get_json(extraction_path, 't3-bc-native-local-matchete-pole-log-split-v1')
    physical_gate = get_json(physical_gate_path, 't3-bc-physical-y-evidence-gate-v1')
    check(w.get('live_Matchete_operator_identity', {}).get('direct_equivalence_status') ==
          'ExactPlusEquality' and w.get('absolute_installed_Matchete_operator_factor_certified') is True,
          'Native positive-sign Warsaw operator identity not certified')
    check(y.get('full_physical_matching_certified') is False and
          f.get('full_physical_matching_certified') is False,
          'Unexpected certification status in input reports')
    check(set(f.get('models',{})) == set(STEMS) and set(y.get('models',{})) == set(STEMS),
          'Missing B/C model reports')
    result = {}
    for stem, model in STEMS.items():
        fr = f['models'][stem]
        yr = y['models'][stem]
        check(fr['model'] == model and yr['model'] == model,
              f'{model}: incorrect v12/v13 model identities')
        check(fr.get('pole_equals_MF_log_coefficient') is True and
              fr.get('F_hard_finite_physical_matching_certified') is False,
              f'{model}: F hard coefficient provenance changed')
        check(yr.get('native_PL_equals_Warsaw_positive_sign') is True and
              yr.get('full_operator_counterterms_certified') is False,
              f'{model}: Y restricted coefficient provenance changed')
        source = root/fr['source_path']
        check(source.is_file() and hash_file(source) == fr['source_sha256'],
              f'{model}: original F raw source changed')
        sub = evidence_for_model(root, stem, extraction, physical_gate, pole)
        params = get_json(root/'output/general_t3_mass_eigen_matchete'/f'{stem}_final_parameters.json')
        mass = params.get('heavy_masses_gev')
        check(isinstance(mass,list) and len(mass)==2 and all(type(x) in (float,int) and x>0 for x in mass),
              f'{model}: missing scalar threshold masses')
        check(mass[0] < mass[1], f'{model}: expected X lighter than Y')
        result[stem] = {
            'model':model,
            'thresholds_GeV':{'Y':mass[1],'X':mass[0]},
            'F_boundary_status':'FORMAL_CONSTANT_NOT_RENORMALIZED',
            'Y_boundary_status':'MSBAR_RESTRICTED_SCALAR_BUBBLE_ONLY',
            'X_boundary_status':'RAW_POLE_LOG_NO_STANDALONE_CONSTANT',
            'F_boundary_input_sha256':fr['source_sha256'],
            'Y_boundary_input_sha256':yr['source_sha256'],
            'raw_scalar_sector_evidence':sub,
            'F_renormalized_finite_missing':True,
            'Y_full_EFT1_minus_EFT2_missing':True,
            'Y_LLXX_inherited_one_loop_into_X_missing':True,
            'X_full_finite_matched_subtraction_missing':True,
            'EFT_running_missing':True,
            'physical_C5':None,
        }
        print(f'{model}: Y/LLXX {sub["Y_LLXX"]["pole_terms"]} pole + '
              f'{sub["Y_LLXX"]["log_terms"]} log + '
              f'{sub["Y_LLXX"]["scale_independent_constant_terms"]} formal finite; '
              f'X/LLHH 18 pole + 18 log + 0 constant; exact pole/log PASS')
    report = {
        'format':'t3-bc-sequential-matching-fixed-order-batch-v14',
        'status':'F_Y_XEvidenceAssembled_CompletePhysicalMatchingOutstanding',
        'provenance_sha256':{
            str(p.relative_to(root)):hash_file(p) for p in (
                audits/'t3_bc_y_conditional_finite_interface_v1/summary.json',
                audits/'t3_bc_f_hard_boundary_v1/summary.json',
                audits/'t3_bc_y_warsaw_direct_identity_v11/summary_live.json',
                extraction_path, physical_gate_path)},
        'models':result,
        'renormalized_matching_definition':
            'Delta C5 at each boundary requires [Gamma^(1)_high - Gamma^(1)_low + delta_CT] finite at same scale, plus finite EOM/field/parameter shifts; formal pole deletion is insufficient',
        'one_loop_sequential_matching_completed':False,
        'never_add_F_and_Y_as_physical_Wilson_without_subtraction':True,
        'X_zero_raw_constant_does_not_imply_zero_physical_threshold':True,
        'physical_C5':None,
        'physical_neutrino_mass_matrix':None,
        'full_physical_matching_certified':False,
    }
    if numeric_input is not None:
        supplied = get_json(numeric_input)
        check(set(supplied) >= {'model', 'MF_GeV', 'y1', 'y2', 'quartics'},
              'Numeric input missing model, MF_GeV, y1/y2 or quartics')
        stem = next((s for s,m in STEMS.items() if m==supplied['model']),None)
        check(stem is not None, 'Numeric input model must be B_-1 or C_-1')
        report['numeric_components'] = {stem: numeric_for_model(root, stem, supplied, v12, v13, kernel)}
        report['numeric_input_sha256'] = hash_file(numeric_input)
    if numeric_manifest is not None:
        manifest = get_json(numeric_manifest)
        check(set(manifest) == set(STEMS.values()),
              'Numeric manifest requires exactly B_-1 and C_-1 records')
        report['numeric_components'] = {}
        for stem, model in STEMS.items():
            item=manifest[model]
            check(isinstance(item, dict) and item.get('model') == model,
                  f'{model}: numeric manifest identity mismatch')
            report['numeric_components'][stem] = numeric_for_model(
                root, stem, item, v12, v13, kernel)
        report['numeric_manifest_sha256'] = hash_file(numeric_manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: v11.2 + v12 + v13 cross-boundary provenance and raw X/Y pole-log')
    print('F and Y kept SEPARATE; X Y-loop inheritance incomplete; physical C5 NOT CERTIFIED')
    print('Saved: '+str(output.resolve()))


def self_test():
    # The test does not need local project parser; use a strict toy parser.
    class Toy:
        @staticmethod
        def terms(s):
            return tuple(s.split(' + '))
    raw='a/(ϵ*m) + a*Log[μbar2/Coupling[MSX, {}, 0]^2]/m + b/m'
    n, rows=unpack_sectors(raw, 'X', Toy)
    check(n==3 and len(rows['pole'])==1 and len(rows['log'])==1 and
          len(rows['constant'])==1 and 'ϵ' not in rows['pole'][0] and
          'Log[' not in rows['log'][0], 'Synthetic scalar pole/log splitting failed')
    for bad in ('a/(ϵ*m)*Log[μbar2/Coupling[MSX, {}, 0]^2]',
                'a*Log[μbar2/Coupling[MSY, {}, 0]^2]'):
        try:
            unpack_sectors(bad, 'X', Toy)
        except BatchError:
            pass
        else:
            raise BatchError('Tampered raw sector was accepted')
    check(RAW_EXPECTATIONS[('Y','LLXX')]['B_m1']==(666,240,240,186),
          'B Y LLXX inventory regression')
    check(RAW_EXPECTATIONS[('X','LLHH')]['C_m1']==(36,18,18,0),
          'C X LLHH inventory regression')
    import numpy as np
    a=np.array([[1+1j,2],[0,3j]],dtype=complex)
    C=a+a.T
    check(np.allclose(C,C.T) and C[0,1]==2 and C[1,0]==2,
          'Ordered Warsaw symmetrization regression')
    print('PASS: strict raw pole/log splitter and wrong-scale/mixed-pole rejection')
    print('PASS: B/C Y LLXX and X LLHH expected inventories; symmetric Warsaw normalization')
    print('PASS: no physical F+Y sum, no inferred zero X physical coefficient')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--numeric-input',type=Path,
                    help='Optional v12-compatible explicit complex couplings JSON, one B/C model')
    ap.add_argument('--numeric-manifest',type=Path,
                    help='Optional B_-1 and C_-1 JSON keyed records, numeric matrices for both models')
    ap.add_argument('--output',type=Path,default=Path('output/audits/t3_bc_sequential_batch_v14/summary.json'))
    ap.add_argument('--require-physical',action='store_true')
    args=ap.parse_args()
    try:
        if args.self_test:
            self_test()
            return 0
        if args.require_physical:
            raise BatchError('Incomplete EFT1/EFT2 UV-IR subtraction, F finite matching, Y inheritance, running and X finite matching; physical C5 cannot be certified')
        check(not(args.numeric_input is not None and args.numeric_manifest is not None),
              'Choose either --numeric-input or --numeric-manifest')
        root=args.root.resolve()
        out=args.output if args.output.is_absolute() else root/args.output
        path=args.numeric_input
        if path is not None and not path.is_absolute():
            path=root/path
        manifest=args.numeric_manifest
        if manifest is not None and not manifest.is_absolute():
            manifest=root/manifest
        run_batch(root,out,path,manifest)
        return 0
    except (BatchError,ValueError,TypeError,KeyError,IndexError,OSError,ImportError) as exc:
        print('FAIL (no physical C5): '+str(exc), file=sys.stderr)
        return 1

if __name__ == '__main__':
    sys.exit(main())
