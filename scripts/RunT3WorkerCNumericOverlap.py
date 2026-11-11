"""Read-only numerical D(-2) restricted-Yukawa/finite scalar overlap diagnostic.

Reuses the repository's existing, provisional Wolfram numeric evaluator against
isolated *copies* of previously restricted WXF candidates. No Match[] is run.
The hard-only matrix (not assembled total) is compared with the optional
canonical Worker-C scalar bubble input. No difference is a certified additive
nonbubble correction until the Matchete/real-tensor normalization map is proved.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from uuid import uuid4

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = ROOT / 'output/general_sequential_pilot/D_m2_stage1/data/general_d_m2_scalar_threshold/profile_tree'
DEFAULT_CONFIG = ROOT / 'benchmarks/general_d_m2_rough_benchmark.json'
DEFAULT_EVALUATOR = ROOT / 'Lagrangian/RunGeneralDMinus2C5RoughBenchmark.wl'
CANDIDATE = 'profile_c5_flavour_candidate_total_provisional.wxf'
AUDIT_INPUTS = ('profile_c5_flavour_candidate_result.json', 'profile_c5_normalization_reference_result.json')
COMPONENTS = {'hard': 'restricted_hard_finite_candidate.wxf',
              'total': 'restricted_total_finite_candidate.wxf'}

# The legacy benchmark demands three entries from SingularValueList, which
# may return only the nonzero singular values for a rank-deficient NF=1 C5.
# In the isolated diagnostic copy, form ALL three singular values from M^dag M.
# Eigenvalues of a 3x3 matrix preserve the zero mode; Clip only removes
# negative eigenvalues from numerical roundoff. The physical matching
# candidate itself is never modified by this evaluator change.
LEGACY_MASS_LINE = 'massesEV = Sort[N[SingularValueList[mnuEV], 12]];'
DIAGNOSTIC_MASS_LINE = (
    'massesEV = Sort[N[Sqrt[Clip[Re[Eigenvalues['
    'ConjugateTranspose[mnuEV].mnuEV]], {0, Infinity}]], 12]];'
)


def isolated_numeric_evaluator(source: Path, destination: Path) -> Path:
    original = source.read_text(encoding='utf-8-sig')
    if original.count(LEGACY_MASS_LINE) != 1:
        raise ValueError('Numeric evaluator source changed: expected unique legacy Takagi line')
    if 'mnuEV = N[-(vev^2/2) matrixC5 10^9, 15];' not in original:
        raise ValueError('Numeric evaluator mass convention is no longer recognized')
    rendered = original.replace(LEGACY_MASS_LINE, DIAGNOSTIC_MASS_LINE)
    destination.write_text(
        '(* Worker C diagnostic copy: preserve zero Takagi mode for NF=1.\n'
        '   Source evaluator and physical C5 expressions remain unchanged. *)\n'
        + rendered, encoding='utf-8')
    return destination


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def preflight(restricted: Path, profile: Path, config: Path, evaluator: Path, output: Path) -> dict:
    restricted = restricted.resolve()
    output = output.resolve()
    if output.exists():
        raise ValueError(f'Refusing to overwrite: {output}')
    if output == restricted or restricted in output.parents or output == profile.resolve() or profile.resolve() in output.parents:
        raise ValueError('Output must not be inside preserved source directories')
    required = [restricted / 'result.json', config, evaluator, *(profile / n for n in AUDIT_INPUTS),
                *(restricted / n for n in COMPONENTS.values())]
    if evaluator.is_file():
        source=evaluator.read_text(encoding='utf-8-sig')
        if source.count(LEGACY_MASS_LINE)!=1 or 'mnuEV = N[-(vev^2/2) matrixC5 10^9, 15];' not in source:
            raise ValueError('Numeric evaluator source changed; zero-mode diagnostic cannot be applied safely')
    absent = [str(p) for p in required if not p.is_file()]
    if absent:
        raise ValueError(f'Missing required files: {absent}')
    restriction = json.loads((restricted / 'result.json').read_text(encoding='utf-8-sig'))
    if (restriction.get('status') != 'RestrictedNF1PoleFreeCandidateExported'
        or restriction.get('restricted_total_pole_zero_exact_NF1') is not True
        or restriction.get('safe_to_add_to_worker_c_joint') is not False):
        raise ValueError('Restricted saved-expression report has not passed its exact NF=1 pole gate')
    params = json.loads(config.read_text(encoding='utf-8-sig'))
    nf = params.get('NF')
    yuk = params.get('Yukawas')
    if nf != 1 or not isinstance(yuk, dict) or not all(k in yuk for k in ('y1','y2','y2Conj')):
        raise ValueError('Only explicitly supplied NF=1 y1/y2/y2Conj benchmark matrices supported')
    for key in ('y1','y2','y2Conj'):
        v = np.asarray(yuk[key])
        if v.shape != (3,1) or not np.isfinite(v).all():
            raise ValueError(f'Invalid numeric 3x1 Yukawa matrix {key}')
    if params.get('PoleHandling', 'Reject') not in ('Reject', 'FinitePartDiagnostic'):
        raise ValueError('Unknown dimensional-pole handling flag')
    if not all(isinstance(params.get(k), (int,float)) and np.isfinite(params[k]) for k in
               ('MF_GeV','MS1_GeV','MS2_GeV','AbsKappaS_GeV2','MuScalar_GeV')):
        raise ValueError('Missing finite masses or scale')
    if min(params['MF_GeV'],params['MS1_GeV'],params['MS2_GeV'],params['MuScalar_GeV']) <= 0:
        raise ValueError('Masses/scale must be positive')
    if params['AbsKappaS_GeV2'] < 0 or params['MS2_GeV']**2 <= 2*params['AbsKappaS_GeV2']:
        raise ValueError('Physical real-S2 mass squared is nonpositive')
    # Never modify original benchmark. An explicitly recorded restricted copy sets
    # the *independent* y2Conj Yukawa to exactly zero for the saved restricted EFT.
    restricted_params = json.loads(json.dumps(params))
    independent_source_nonzero = bool(np.any(np.asarray(yuk['y2Conj']) != 0))
    restricted_params['Yukawas']['y2Conj'] = [[0.0] for _ in range(3)]
    restricted_params['PoleHandling'] = 'FinitePartDiagnostic'  # hard-alone retains a pole; no renormalization claim
    return {'source_report': restriction, 'benchmark': params,
            'restricted_benchmark': restricted_params,
            'independent_y2Conj_zeroed_in_isolated_benchmark': independent_source_nonzero}


def load_complex_matrix(real_path: Path, imag_path: Path) -> np.ndarray:
    def read(path: Path) -> np.ndarray:
        with path.open('r', encoding='utf-8-sig', newline='') as f:
            rows = [[float(x) for x in row] for row in csv.reader(f)]
        a = np.asarray(rows, dtype=float)
        if a.shape != (3,3) or not np.isfinite(a).all():
            raise ValueError(f'Invalid 3x3 finite matrix in {path}')
        return a
    out = read(real_path) + 1j * read(imag_path)
    resid = np.linalg.norm(out-out.T) / max(np.linalg.norm(out),1e-100)
    if resid > 1e-8:
        raise ValueError(f'C5 matrix is not symmetric: relative residual {resid}')
    return out


def compare_bubble(meta_path: Path, hard: np.ndarray, params: dict) -> dict:
    from Numerical.joint_scalar_matching.adapter import load_input
    from Numerical.joint_scalar_matching.core import match_c5

    inp = load_input(meta_path)
    if inp.model != 'D_-2' or inp.is_fixture:
        raise ValueError('Comparison requires non-fixture D_-2 real-tensor input')
    mu = float(params['MuScalar_GeV'])
    if not np.isclose(inp.mu_s_gev, mu, rtol=1e-12, atol=0):
        raise ValueError('Matchete benchmark and tensor input have different mu_S')
    if not np.isclose(inp.mu_f_gev, float(params['MF_GeV']),rtol=1e-12,atol=0):
        raise ValueError('Matchete benchmark and tensor input have different fermion masses / boundary')
    ms1,ms2,k = (float(params[n]) for n in ('MS1_GeV','MS2_GeV','AbsKappaS_GeV2'))
    expected = np.sort(np.array([ms1**2]*6 + [ms2**2-2*k,ms2**2+2*k]))
    actual = np.sort(np.linalg.eigvalsh(inp.mass_squared))
    if not np.allclose(actual,expected,rtol=1e-9,atol=1e-6):
        raise ValueError('Real-heavy-scalar mass spectrum does not match D_-2 benchmark')
    result=match_c5(inp.c_llss,inp.lambda_hhss,inp.mass_squared,mu,
                    group_names=['S1']*6+['S2']*2,strict_projection=True)
    bubble=result['c5_scalar_bubble']
    residual=hard-bubble
    nh=float(np.linalg.norm(hard)); nb=float(np.linalg.norm(bubble)); nr=float(np.linalg.norm(residual))
    return {
        'status':'ConditionalNumericalOverlapOnly',
        'tensor_input':str(meta_path.resolve()),
        'tensor_input_sha256':digest(meta_path),
        'tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
        'same_model_muF_muS_and_scalar_spectrum':True,
        'quartic_and_LLSS_Matchete_to_real_basis_normalization_certified':False,
        'same_physical_couplings_certified':False,
        'matchete_hard_finite_c5':complex_grid(hard),
        'worker_c_bubble_c5':complex_grid(bubble),
        'conditional_difference_not_for_addition':complex_grid(residual),
        'hard_frobenius_norm_gev_inv':nh,
        'bubble_frobenius_norm_gev_inv':nb,
        'difference_frobenius_norm_gev_inv':nr,
        'relative_difference_to_max_norm':nr/max(nh,nb,1e-100),
        'projection_residual':result['projection_residual'],
        'independent_nonbubble_remainder_certified':False,
        'safe_to_add_to_worker_c_joint':False,
    }


def complex_grid(a: np.ndarray) -> dict:
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def run(*, restricted: Path, profile: Path, config: Path, evaluator: Path,
        output: Path, bubble_input: Path | None = None,
        wolframscript: str = 'wolframscript', dry_run: bool = False) -> Path:
    inputs = preflight(restricted,profile,config,evaluator,output)
    if dry_run:
        return output
    output.mkdir(parents=True, exist_ok=False)
    params = inputs['restricted_benchmark']
    benchmark_file = output/'restricted_benchmark.json'
    benchmark_file.write_text(json.dumps(params,indent=2)+'\n', encoding='utf-8')
    evaluator_copy=isolated_numeric_evaluator(evaluator,output/'worker_c_benchmark_evaluator.wl')
    report = {
        'status':'InProgress', 'model':'D_-2','NF':1,
        'source_restricted_report':str((restricted/'result.json').resolve()),
        'source_sha256':{'restricted_result':digest(restricted/'result.json'),
                         'source_benchmark':digest(config),
                         'source_evaluator':digest(evaluator),
                         'isolated_evaluator':digest(evaluator_copy),
                         **{n:digest(restricted/f) for n,f in COMPONENTS.items()}},
        'restricted_yukawa_config_sha256':digest(benchmark_file),
        'independent_y2Conj_zeroed_in_isolated_benchmark':inputs['independent_y2Conj_zeroed_in_isolated_benchmark'],
        'reused_existing_numeric_evaluator_no_source_change':True,
        'isolated_takagi_zero_mode_workaround':True,
        'isolated_evaluator_path':str(evaluator_copy),
        'no_new_Matchete_matching':True,
        'renormalization_certified':False,
        'operator_and_CG_normalization_certified':False,
        'nonbubble_remainder_certified':False,
        'safe_to_add_to_worker_c_joint':False,
        'components':{},
        'pole_bookkeeping':'Hard component alone has a nonzero pole; full restricted total has exact NF=1 pole cancellation in source. Both are finite-part diagnostics, not renormalized hard Wilsons.'
    }
    report_path=output/'result.json'
    try:
        for label, name in COMPONENTS.items():
            # Original evaluator writes ALL files to its profileDir. Stage ONLY
            # the minimum source artifacts under the *fresh* audit directory.
            stage=output/label/'data/general_d_m2_scalar_threshold/profile_tree'
            stage.mkdir(parents=True,exist_ok=False)
            shutil.copy2(restricted/name,stage/CANDIDATE)
            for audit in AUDIT_INPUTS:
                shutil.copy2(profile/audit,stage/audit)
            source_root=output/label
            command=[wolframscript,'-file',str(evaluator_copy),str(source_root),str(benchmark_file)]
            log=output/f'{label}_wolfram.log'
            with log.open('w',encoding='utf-8') as fh:
                proc=subprocess.run(command,cwd=ROOT,stdout=fh,stderr=subprocess.STDOUT,
                                    check=False)
            if proc.returncode:
                raise RuntimeError(f'{label} Wolfram evaluator exited {proc.returncode}; see {log}')
            result_file=stage/'profile_c5_rough_benchmark_result.json'
            if not result_file.is_file():
                raise RuntimeError(f'{label} Wolfram exited 0 but result is missing; see {log}')
            entry=json.loads(result_file.read_text(encoding='utf-8-sig'))
            if (entry.get('Status') != 'ProvisionalFinitePartBenchmarkEvaluated'
                or entry.get('PoleHandling') != 'FinitePartDiagnostic'
                or entry.get('FinitePartOnly') is not True
                or entry.get('NoNewMatchPerformed') is not True
                or entry.get('NumericalBenchmarkCertifiedPhysical') is not False
                or entry.get('RenormalizationCertified') is not False
                or not isinstance(entry.get('C5PoleResidueMaxAbs_GeVinv'), (int,float))
                or not np.isfinite(entry['C5PoleResidueMaxAbs_GeVinv'])
                or entry['C5PoleResidueMaxAbs_GeVinv'] < 0):
                raise RuntimeError(f'{label} Wolfram result failed finite-part diagnostic gate; see {result_file}')
            real=stage/'profile_c5_rough_benchmark_C5_real_GeVinv.csv'
            imag=stage/'profile_c5_rough_benchmark_C5_imag_GeVinv.csv'
            matrix=load_complex_matrix(real,imag)
            report['components'][label]={
                'wolfram_status':entry['Status'],'wolfram_report':str(result_file),
                'real_csv':str(real),'imag_csv':str(imag),
                'finite_c5_gev_inv':complex_grid(matrix),
                'norm_gev_inv':float(np.linalg.norm(matrix)),
                'pole_residue_max_abs_gev_inv':float(entry['C5PoleResidueMaxAbs_GeVinv']),
                'pole_reported':bool(entry['DimRegPolePresent']),
                'finite_part_only':True, 'renormalization_certified':False}
        hard=load_complex_matrix(Path(report['components']['hard']['real_csv']),Path(report['components']['hard']['imag_csv']))
        total=load_complex_matrix(Path(report['components']['total']['real_csv']),Path(report['components']['total']['imag_csv']))
        hard_pole=report['components']['hard']['pole_residue_max_abs_gev_inv']
        total_pole=report['components']['total']['pole_residue_max_abs_gev_inv']
        # The complete restricted NF=1 expression has an independently proven
        # exact symbolic zero residue; floating point evaluation can leave noise.
        # This tolerance is diagnostic only, never a renormalization proof.
        noise_bound=max(1e-24,1e-10*hard_pole)
        report['total_numeric_pole_consistent_with_exact_symbolic_zero']=(total_pole <= noise_bound)
        report['total_numeric_pole_noise_bound_gev_inv']=noise_bound
        if total_pole > noise_bound:
            raise RuntimeError(f'Total numeric pole {total_pole} exceeds diagnostic noise bound {noise_bound}; no valid total finite-part comparison')
        report['total_minus_hard_not_nonbubble_correction']=complex_grid(total-hard)
        report['total_contains_inherited_and_running']=True
        if bubble_input:
            report['overlap']=compare_bubble(bubble_input,hard,params)
        else:
            report['overlap']={'status':'NotAttemptedNoCanonicalTensorInput',
                               'reason':'A D_-2 Worker C canonical NPZ+JSON tensor input at the exact same benchmark couplings is required',
                               'safe_to_add_to_worker_c_joint':False}
        report['status']='RestrictedHardAndTotalFinitePartsEvaluatedDiagnosticOnly'
    except (OSError,ValueError,RuntimeError) as exc:
        report['status']='Blocked'
        report['error']=str(exc)
        report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        raise
    report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report_path


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restricted-dir',type=Path,required=True)
    parser.add_argument('--profile',type=Path,default=DEFAULT_PROFILE)
    parser.add_argument('--benchmark',type=Path,default=DEFAULT_CONFIG)
    parser.add_argument('--evaluator',type=Path,default=DEFAULT_EVALUATOR)
    parser.add_argument('--bubble-input',type=Path,default=None,
                        help='Optional canonical non-fixture D_-2 metadata JSON with sibling NPZ')
    parser.add_argument('--output-dir',type=Path,default=None)
    parser.add_argument('--preflight',action='store_true')
    args=parser.parse_args()
    output=args.output_dir or (ROOT/'output/audits'/('worker_c_numeric_overlap_D_m2_'+
        datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid4().hex[:8]))
    try:
        result=run(restricted=args.restricted_dir,profile=args.profile,config=args.benchmark,
                   evaluator=args.evaluator,output=output,bubble_input=args.bubble_input,
                   dry_run=args.preflight)
    except (OSError,ValueError,RuntimeError) as exc:
        print('BLOCKED:',exc,file=sys.stderr)
        return 2
    if args.preflight:
        print('PASS: isolated Matchete evaluator inputs valid; no numeric evaluation attempted')
    else:
        print('PASS: restricted hard and total FINITE-PART numerical diagnostics:',result)
        print('Independent nonbubble correction safe to add: False')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
