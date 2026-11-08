"""A(-4) scalar-loop proof + isolated native RGBeta dummy-index intervention.

Run ONCE from the local repository root, after extracting the batch ZIP:
  python scripts/RunAminus4FreshIndexCausalBatch.py \
    --output-root output/a_minus4_next_batch/fresh_index_causal

The module Numerical.running.Aminus4QuarticHessianIdentities derives all 55
one-loop scalar-loop monomials from 1/2 Tr(V''^2), with exact SymPy algebra;
no native RGBeta beta expressions are used in that derivation.

The Wolfram probe starts 4 fresh RGBeta processes and *only in memory*
re-registers either/both adjoint quartics with Module-local dummy indices.
No Matchete, installed RGBeta, original RGBeta source or production ODE is edited.
The resulting 55 native polynomials are compared to the exact Hessian result.
A change is source-native experimental evidence, not a proof of the causal
mechanism inside the RGBeta projector without kernel-level inspection.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
import shutil
import subprocess
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from Numerical.running import Aminus4EFT1Audit as audit
from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.running.Aminus4QuarticHessianIdentities import (
    COORDINATES, beta_monomials_exact, evaluated_beta,
    group_identities, monomial_json, restricted_coefficients,
)
from scripts.VerifyAminus4GlobalQuarticProjection import unit_point

PROJECTED_RAW = Path('output/a_minus4_eft1_batch2/A_m4_projected_raw.json')
SOURCE_AUDIT = Path('output/a_minus4_eft1_batch2/A_m4_projection_audit.json')
DECLARED_SPINOR = Path('output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json')
MODEL = Path('RGE/running/rgbeta/T3RGBetaModel.wl')
WOLFRAM_SCRIPT = Path('scripts/RunAminus4RGBetaFreshIndexProbe.wl')
MODEL_SHA256 = 'a58c79f18b2abdf9d476588552dfa52e2f131f322094e5b72d1198162aa20abc'
MODES = ('baseline', 'triplet', 'portal', 'both')
QUARTICS = frozenset({'lambdaH','lambdaS1','lambdaS2','lambdaS2Adj',
                       'lambdaH1','lambdaH2','lambdaH2Adj','lambda12','lambdaT3'})


def require(condition: bool, error: str) -> None:
    if not condition:
        raise ValueError(error)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing source {path}')
    obj = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj, dict), f'Expected JSON object: {path}')
    return obj


def file_sha(path: Path) -> str:
    require(path.is_file(), f'Missing required file: {path}')
    return sha256(path.read_bytes()).hexdigest()


def support_point(names: tuple[str, ...]) -> dict:
    return {name: int(name in names) for name in COORDINATES}


def expected_reference(point: dict, exact: dict) -> np.ndarray:
    predicted = evaluated_beta(exact, support_point(tuple(point)))
    return np.array([float(predicted[name]) for name in COORDINATES])


def payload_for_native(raw: dict, native: dict) -> dict:
    require(native.get('model') == 'A_-4' and native.get('stage') == 'EFT1',
            'Wrong native model/stage')
    require(native.get('source_model_SHA256') == MODEL_SHA256,
            'Unrecognized native model source')
    require(native.get('status') ==
            'NATIVE_QUARTIC_BETAS_EXPORTED_NOT_PHYSICAL_CERTIFICATION',
            'Native beta calculation failed or changed status')
    beta = native.get('report_betas')
    require(isinstance(beta, dict) and set(beta) == QUARTICS and
            all(isinstance(v, str) and v.strip() for v in beta.values()),
            'Missing native quartic beta')
    result = deepcopy(raw)
    result['metadata'].update({'IntegratedField':'F','ActiveBSMFields':['S1','S2']})
    result['report_betas'].update(beta)
    eft.validate_beta_payload(result)
    return result


def inspect_polynomial(payload: dict, exact: dict, *, original: dict | None = None) -> dict:
    """Independently compare 55 native scalar-only probes to exact Hessian.

    The shared real-scalar 4th-derivative audit is used ONLY as an additional
    cross-check that the declared 12-real-scalar potential basis agrees with
    the independently constructed symbolic polynomial.
    """
    single = [(name,) for name in COORDINATES]
    pairs = list(combinations(COORDINATES, 2))
    all_probes = single + pairs
    measured: dict[tuple[str, ...], np.ndarray] = {}
    failed = []
    max_ref_err = 0.
    max_base_diff = 0.
    max_closure = 0.
    rows = []
    for active in all_probes:
        point = unit_point(active)
        comp = audit.scalar_beta_comparison(point, payload)
        reference = np.array([comp['coupling_beta_comparison'][name]['reference_scalar_loop']
                              for name in COORDINATES])
        native = np.array([comp['coupling_beta_comparison'][name]['native_rgbeta_scalar_loop']
                           for name in COORDINATES])
        expected = expected_reference(dict.fromkeys(active), exact)
        ref_error = float(np.max(np.abs(reference - expected)))
        require(ref_error < 2e-8, f'Independent Hessian disagrees with real-scalar tensor at {active}: {ref_error}')
        max_ref_err = max(max_ref_err, ref_error)
        max_closure = max(max_closure,
            float(comp['scalar_loop_projection_closure_relative']))
        deviation = native - expected
        measured[active] = deviation
        if float(np.max(np.abs(deviation))) > 2e-8:
            failed.append('+'.join(active))
        if original is not None:
            prev = audit.scalar_beta_comparison(point, original)
            prev_native = np.array([prev['coupling_beta_comparison'][name]['native_rgbeta_scalar_loop']
                                    for name in COORDINATES])
            max_base_diff = max(max_base_diff, float(np.max(np.abs(prev_native-native))))
        rows.append({'active':list(active),
                     'native_minus_exact_hessian':{
                         key: float(value) for key, value in zip(COORDINATES,deviation)
                         if abs(value) > 2e-8}})

    # Polarisation: isolate all 55 independent quadratic monomial differences.
    found = {}
    singles = {k: measured[(k,)] for k in COORDINATES}
    for key in COORDINATES:
        found[(key,key)] = singles[key]
    for a,b in pairs:
        found[(a,b)] = measured[(a,b)] - singles[a] - singles[b]
    support=[]
    for (a,b),vector in found.items():
        diff={key:float(v) for key,v in zip(COORDINATES,vector)
              if abs(v)>2e-8}
        if diff:support.append({'couplings':[a,b],'beta_differences':diff})
    return {
        'native_vs_exact_hessian_failing_probes':len(failed),
        'failing_probe_names':failed,
        'independent_discrepant_monomials':support,
        'discrepant_monomial_count':len(support),
        'max_reference_hessian_minus_fourth_derivative_abs':max_ref_err,
        'max_reference_projection_closure':max_closure,
        'max_difference_from_archived_raw_global_beta':max_base_diff if original is not None else None,
        'all_55_probe_details':rows,
    }


def signed_weyl_convention_check(declared: dict) -> dict:
    require(declared.get('model') == 'A_-4' and
            declared.get('native_matchete_to_declared_Weyl_factor') is None and
            declared.get('physical_certification') is False,
            'Native signed spinor evidence changed; source re-audit required')
    plus=declared.get('declared_rho_plus_one_independent_coordinates',[])
    minus=declared.get('rho_minus_one_fixed_native_y1_alternative_coordinates',[])
    require(len(plus)==len(minus)==36,'Changed native/master 36 coordinate support')
    def decode(rows):
        ans={}
        for row in rows:
            key=tuple(row[n] for n in ('i','j','a','b'))
            require(key not in ans, 'Duplicate Wilson coordinate')
            ans[key]=complex(row['value_re_gev_inv'],row['value_im_gev_inv'])
        return ans
    a,b=decode(plus),decode(minus)
    require(set(a)==set(b),'Opposite sign changes master coordinate support')
    scale=max(abs(v) for v in a.values())
    relative=max(abs(a[k]+b[k]) for k in a)/scale
    require(scale>0 and relative<1e-12,'Failed signed Weyl-branch algebra')
    return {
      'source_two_fixed_native_y1_branches_have_opposite_tensor':True,
      'max_relative_opposite_sign_residual':relative,
      'rho_plus_one_is_declared_not_measured':True,
      'native_4x4_GammaCC_signed_factor_measured':False,
      'linear_Wilson_running_respects_global_C_to_minus_C':True,
      'statement_of_last_field_is_conditional_on_linearity_of_beta_C_and_C_conjugate':True,
      'Takagi_masses_invariant_under_global_mnu_to_minus_mnu':True,
      'does_not_certify_relative_phase_of_unmatched_hard_loop_amplitudes':True,
      'physical_C_at_muF':None,
    }


def source_checks(root: Path) -> dict:
    src = root/MODEL
    got = file_sha(src)
    require(got == MODEL_SHA256,
            f'Native T3RGBetaModel.wl differs from audited source: {got}')
    text = src.read_text(encoding='utf-8-sig')
    require('T3RGBetaPortalAdjInvariant[3]' in text and
            'lambdaH2Adj,' in text and 'lambdaS2Adj,' in text and
            'tGen[SU2L @ S2, A,' in text,
            'Relevant fixed adjoint-index source declarations changed')
    require((root/WOLFRAM_SCRIPT).is_file(), 'Native Wolfram probe not installed')
    return {'model_path':str(MODEL),'model_sha256':got,
            'original_fixed_adjoint_dummy_present':True,
            'no_source_mutation_or_global_rule_override':True}


def run_wolfram(root: Path, output_root: Path, mode: str, timeout: int) -> dict:
    exe=shutil.which('wolframscript')
    require(exe is not None,'wolframscript not installed or not on PATH')
    outfile=output_root/'native'/f'{mode}.json'
    require(not outfile.exists() and not Path(str(outfile)+'.partial').exists(),
            f'Refusing native report overwrite: {outfile}')
    outfile.parent.mkdir(parents=True, exist_ok=True)
    command=[exe, '-file',str(root/WOLFRAM_SCRIPT),mode,str(outfile)]
    try:
        process=subprocess.run(command,cwd=root,text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=timeout, check=False)
    except subprocess.TimeoutExpired as e:
        raise RuntimeError(f'{mode}: native RGBeta timeout after {timeout}s') from e
    (output_root/'native'/f'{mode}.log').write_text(
        'Command: '+subprocess.list2cmdline(command)+'\n'+
        f'Exit code: {process.returncode}\nSTDOUT:\n{process.stdout}\n'+
        f'STDERR:\n{process.stderr}\n', encoding='utf8')
    if process.returncode:
        raise RuntimeError(f'{mode}: WolframScript exit {process.returncode}; '
                           f'see native/{mode}.log')
    data=read_json(outfile)
    require(data.get('mode')==mode, 'Native mode differs from requested experiment')
    return data


def evaluate(root: Path, output_root: Path, *, modes: tuple[str, ...], timeout: int) -> dict:
    require(all(m in MODES for m in modes) and 'baseline' in modes and len(set(modes))==len(modes),
            'Native modes must be unique and include baseline')
    root=root.resolve(); output_root=output_root.resolve()
    output_root.mkdir(parents=True,exist_ok=True)
    report_file=output_root/'causal_comparison.json'
    require(not report_file.exists(),'Refusing to overwrite completed native comparison')
    hessian=beta_monomials_exact()
    group=group_identities()
    require(group['fundamental_Fierz_exact'] and group['spin_one_adjoint_square_exact'],
            'Exact SU2 group identities failed')
    source=source_checks(root)
    raw=read_json(root/PROJECTED_RAW)
    audited=read_json(root/SOURCE_AUDIT)
    require(audited.get('new_global_beta_export_sha256')==file_sha(root/PROJECTED_RAW),
            'Native global source beta differs from archived audit')
    declared=read_json(root/DECLARED_SPINOR)
    signed=signed_weyl_convention_check(declared)
    # Save a fully independent exact Hessian derivation regardless of Wolfram.
    proof={'format':'a-minus4-exact-component-Hessian-quartic-proof-v1',
        'model':'A_-4',
        'loop_identity': '16*pi^2*dV4/dln(mu)|scalar = Tr[(d2 V4 / dphi_a dphi_b)^2]/2',
        'potential_convention': 'V=(lambdaH/2)(Hdag H)^2 + (lambdaS1/2)|S1|^4 + (lambdaS2/2)(S2dag S2)^2 + (lambdaS2Adj/2)sum_A(S2dag T_A S2)^2 + lambdaH1(Hdag H)|S1|^2 + lambdaH2(Hdag H)(S2dag S2) + lambdaH2Adj sum_A(Hdag t_A H)(S2dag T_A S2) + lambda12|S1|^2(S2dag S2) + lambdaT3(HH)_3 S1 S2dag + h.c.',
        'independent_12_real_scalar_components':12,
        'independent_real_quartic_coordinates':list(COORDINATES),
        'one_loop_quadratic_monomials':len(hessian),
        'full_12_component_polynomial_closure_all_monomials':True,
        'SU2_group_identites':group,
        'all_exact_coefficients':monomial_json(hessian),
        'five_disputed_reference_channels':monomial_json(restricted_coefficients(hessian)),
        'physical_mapped_to_Matchete':False,
        'native_RGBeta_used_to_derive':False,
    }
    (output_root/'exact_hessian_reference.json').write_text(
        json.dumps(proof,indent=2)+'\n',encoding='utf8')
    report={'format':'a-minus4-fresh-index-native-causal-batch-v1',
       'model':'A_-4','EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars -> SMEFT',
       'source_verification':source,'fresh_index_modes_requested':list(modes),
       'exact_hessian_proof':'exact_hessian_reference.json',
       'spinor_signed_branch_convention':signed,
       'native_modes':{},
       'native_origin_of_discrepancy_proven':False,
       'matchete_to_RGBeta_invariant_map_certified':False,
       'physical_Wilson_boundary_certified':False,
       'production_source_modified':False,
       'status':'EXACT_HESSIAN_PROVEN_NATIVE_EXPERIMENT_PENDING',
       'errors':[],
    }
    baseline=deepcopy(raw)
    baseline['metadata'].update({'IntegratedField':'F','ActiveBSMFields':['S1','S2']})
    eft.validate_beta_payload(baseline)
    # First re-check the archived original RGBeta coefficients with independent Hessian.
    old=inspect_polynomial(baseline,hessian)
    require(old['native_vs_exact_hessian_failing_probes']==34 and
            old['discrepant_monomial_count']==5,
            'Archived RGBeta discrepancy no longer matches 34/5 evidence')
    report['archived_global_beta_vs_exact_hessian']={
        k:v for k,v in old.items() if k!='all_55_probe_details'
    }
    try:
        for mode in modes:
            print(f'[native RGBeta] {mode}: global quartic projection',flush=True)
            native=run_wolfram(root,output_root,mode,timeout)
            candidate=payload_for_native(raw,native)
            comparison=inspect_polynomial(candidate,hessian,
                    original=baseline if mode=='baseline' else None)
            report['native_modes'][mode]={
                'fresh_triplet_adjoint_dummy':native['fresh_triplet_adjoint_dummy'],
                'fresh_portal_adjoint_dummy':native['fresh_portal_adjoint_dummy'],
                'source_sha256':native['source_model_SHA256'],
                'beta_output':f'native/{mode}.json',
                'analysis':{k:v for k,v in comparison.items()
                            if k!='all_55_probe_details'},
            }
            print(f'  failed probes: {comparison["native_vs_exact_hessian_failing_probes"]}/55; '
                  f'discrepant monomials: {comparison["discrepant_monomial_count"]}',flush=True)
            if mode=='baseline':
                maxdifference=comparison['max_difference_from_archived_raw_global_beta']
                require(maxdifference<2e-8,
                        f'Installed baseline RGBeta differs from archived source '
                        f'(max {maxdifference:.3e}); intervention uninterpretable')
        report['status']='NATIVE_FRESH_INDEX_INTERVENTIONS_EVALUATED'
    except (Exception) as exc:
        report['status']='NATIVE_CAUSAL_EXPERIMENT_BLOCKED'
        report['errors'].append(str(exc))
    (report_file).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',
                            encoding='utf8')
    return report


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root',type=Path,
        default=Path('output/a_minus4_next_batch/fresh_index_causal'))
    p.add_argument('--modes',nargs='+',choices=MODES,default=list(MODES))
    p.add_argument('--timeout-seconds',type=int,default=900)
    args=p.parse_args(argv)
    if args.timeout_seconds<30:
        p.error('Require at least 30s for an independent native kernel')
    results=evaluate(ROOT,args.output_root,modes=tuple(args.modes),
                     timeout=args.timeout_seconds)
    print('Independent exact Hessian: all 55 scalar-only monomials derived')
    print('SU(2) generator identities: exact')
    print('Archived source beta disagreements: 34/55 = five monomials')
    for mode,data in results['native_modes'].items():
        print(f'{mode}: {data["analysis"]["native_vs_exact_hessian_failing_probes"]}/55 discrepancies')
    print('Native experiment status:',results['status'])
    print('Matchete scalar map/physical Wilson boundary: NOT CERTIFIED')
    print('Output:',args.output_root.resolve()/'causal_comparison.json')
    return 0 if results['status']=='NATIVE_FRESH_INDEX_INTERVENTIONS_EVALUATED' else 1


if __name__=='__main__':
    raise SystemExit(main())
