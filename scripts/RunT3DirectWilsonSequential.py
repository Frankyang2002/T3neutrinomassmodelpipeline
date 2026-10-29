"""Candidate F-first direct Weinberg Wilson source with the *running* EFT1 scalar sector.

For models with a saved hierarchical MSbar hard threshold and a certified
EFT1 direct Weinberg beta, this batch computes

    K(muF) = sum_r[K_hard(MF_r) + beta_direct,r(lambdaF)*ln(muF/MF_r)],
    K(muS) = K(muF) + integral[ln muF, ln muS] beta_direct(t) dt,
    K(mu_low) = SM+WeinbergRG(K(muS)).

The direct source is lifted to full flavour using the repository's symmetric
Yukawa convention; the *evolved* fermion-threshold Yukawas/MF are held fixed
below MF. The full RGBeta general scalar EFT1 state is re-integrated from
its SAVED fermion-threshold endpoint, with endpoint parity verification.
There is NO new UV integration and NO retuning of historical benchmarks.

IMPORTANT: this is *NOT* complete general-scalar sequential Wilson matching:
other LLSS operators, homogeneous Wilson mixing, bilinear/extra-quartic
sources, finite scalar-threshold corrections, and complete matching scheme
certification are omitted. Four models without saved hierarchical hard
thresholds remain marked as explicit legacy controls, not invented matches.

    python -m scripts.RunT3DirectWilsonSequential --all --run-id continuous_v1
    python -m scripts.RunT3DirectWilsonSequential --all --run-id continuous_v1 --resume
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import math
from pathlib import Path
import re
from typing import Any

import numpy as np
import sympy as sp
from scipy.integrate import simpson

from scripts import RunT3FullScalarContinuous as continuous
from scripts import RunT3FullScalarNeutrino as neutrino

ROOT = continuous.ROOT
MODELS = continuous.MODELS
FORMAT = 't3-direct-source-sequential-wilson-candidate-v1'
STATUS = 'CandidateHierarchicalHardPlusTimeResolvedDirectWilsonTransport'
SKIP = 'LegacyOnlyNoSavedHierarchicalHardBoundary'
DEFAULT_OUTPUT_ID = 'direct_wilson_sequential_v1'
HBAR = 1.0 / (16.0 * math.pi**2)


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('T3 direct sequential Wilson: ' + message)


def sha(path: Path) -> str:
    require(Path(path).is_file(), 'missing input '+str(path))
    return sha256(Path(path).read_bytes()).hexdigest()


def model_sources(model: str) -> tuple[Path, Path]:
    require(model in MODELS, 'unknown T3 model '+model)
    group, alpha = model.split('_')
    token = ('m' if alpha.startswith('-') else 'p') + alpha[1:]
    base = ROOT/f'output/full/hypercharge/T3_{group}_alpha_{token}/data'
    return base/'final_weinberg_coefficient.json', base/'eft1_wilson_rge.json'


def pair_matrix(weights: np.ndarray, y1: np.ndarray, y2: np.ndarray) -> np.ndarray:
    """Match project's Majorana orientation: (y1bar y2bar + y2bar y1bar)."""
    w = np.asarray(weights, dtype=complex)
    a, b = np.asarray(y1, dtype=complex), np.asarray(y2, dtype=complex)
    require(a.shape == (3, 3) and b.shape == (3, 3) and w.shape == (3,),
            'invalid heavy-flavour kernel input')
    require(np.isfinite(w).all() and np.isfinite(a).all() and np.isfinite(b).all(),
            'nonfinite heavy flavour input')
    x, z = np.conj(a), np.conj(b)
    return neutrino.validated_c5((x*w)@z.T+(z*w)@x.T)


@dataclass(frozen=True)
class DirectWilsonSource:
    hard: Any
    derivative: Any
    hard_text: str
    beta_text: str
    running_text: str

    def matrix(self, masses: np.ndarray, y1: np.ndarray, y2: np.ndarray,
               lam: complex, *, derivative: bool=False) -> np.ndarray:
        fn = self.derivative if derivative else self.hard
        masses = np.asarray(masses, dtype=float)
        require(masses.shape==(3,) and np.isfinite(masses).all() and np.min(masses)>0,
                'bad fermion threshold masses')
        values=np.asarray([complex(fn(float(m), complex(lam), HBAR)) for m in masses],
                          dtype=complex)
        return pair_matrix(values,y1,y2)


def compile_direct_source(matching: Path, wilson: Path) -> DirectWilsonSource:
    """Source-derived hard term and d/dln(mu) of its direct EFT1 log.

    Different flavour 1/2 conventions are a common source of accidental
    factors of two. We differentiate the authoritative *ordered* running
    kernel and compare its coefficient to the *independent* one-generation
    Wilson RGE, then use the repository's Majorana symmetrisation verbatim.
    """
    from RGE.matching.FinalWeinbergAdapter import load_hierarchical_majorana_c5
    hard_json = neutrino.load(matching)
    wilson_json = neutrino.load(wilson)
    require(hard_json.get('status')=='Success' and
            hard_json.get('scheme')=='fixed_order_one_loop_MSbar_hierarchical' and
            hard_json.get('combined',{}).get('ready_for_physical_majorana_numerics') is True,
            'non-certified saved hierarchical hard coefficient')
    check = wilson_json.get('weinberg_subspace_validation',{})
    require(wilson_json.get('status')=='Success' and wilson_json.get('loop_order')==1
            and wilson_json.get('fixed_order_input')=='tree-level EFT1 Wilson coefficients only'
            and check.get('matches_weinberg_subspace') is True
            and check.get('failure_count')==0,
            'one-generation direct Wilson source not certified')
    symbolic=load_hierarchical_majorana_c5(matching,n_lepton=1,n_heavy=1)
    hard=sp.simplify(symbolic['hard_kernel'])
    running=sp.simplify(symbolic['running_kernel'])
    MF, MS, lam, hbar = sp.symbols('MF MS lambdaT3 hbar')
    derivative=sp.simplify(MS*sp.diff(running,MS))
    free={MF,lam,hbar}
    require(hard.free_symbols <= free and derivative.free_symbols <= free,
            'hard or direct source depends on unhandled scalar masses/couplings')
    require(sp.simplify(running-derivative*sp.log(MS/MF))==0,
            'direct source is not a pure fixed-order log')
    require(sp.simplify(derivative-lam*derivative.subs(lam,1))==0 and
            sp.simplify(hard-lam*hard.subs(lam,1))==0,
            'C5 hard/direct source is not linear in lambdaT3')
    y1,y2=sp.symbols('y1 y2')
    rge_expr=sp.sympify(check['beta_kappa_16pi2'],
                 locals={'sqrt':sp.sqrt,'conjugate':sp.conjugate})
    expected=sp.simplify(rge_expr/(lam*sp.conjugate(y1)*sp.conjugate(y2)/MF))
    actual=sp.simplify(derivative/(hbar*lam/MF))
    require(not expected.free_symbols and sp.simplify(expected-actual)==0,
            'Wilson source beta differs from hierarchical running-log derivative')
    variables=(MF,lam,hbar)
    hfn=sp.lambdify(variables,hard,modules='numpy',cse=True)
    sfn=sp.lambdify(variables,derivative,modules='numpy',cse=True)
    return DirectWilsonSource(hfn,sfn,str(hard),str(expected),str(running))


def heavy_inputs(uv: Any) -> tuple[np.ndarray,np.ndarray,np.ndarray,float]:
    from Numerical.running.FinalC5TrajectoryAdapter import _heavy_mass_basis_data
    x=_heavy_mass_basis_data(uv,offdiagonal_atol=1e-7,imaginary_atol=1e-9,
                 takagi_residual_rtol=1e-8,takagi_residual_atol=1e-7)
    require(len(x.masses_gev)==3 and np.min(x.masses_gev)>0,
            'heavy spectrum invalid at fermion threshold')
    return (np.asarray(x.masses_gev,dtype=float),
            np.asarray(x.y1,dtype=complex),np.asarray(x.y2,dtype=complex),
            float(max(x.takagi_residual,x.basis_residual)))


def grouped_hard_boundary(source: DirectWilsonSource, masses: np.ndarray,
                          y1: np.ndarray, y2: np.ndarray, lam: complex,
                          *, mu_f: float) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Relocate per-eigenmass MSbar hard matching to the grouped F scale.

    The saved analytic direct kernel logs ln(MS/MF_r), whereas the numerical
    EFT1 trajectory begins at mu_f. To avoid dropping or double counting
    ln(mu_f/MF_r), put that term in the grouped F boundary before evolving.
    """
    require(np.isfinite(mu_f) and mu_f>0, 'invalid grouped F threshold')
    mf=np.asarray(masses,dtype=float)
    require(mf.shape==(3,) and np.isfinite(mf).all() and np.min(mf)>0,
            'invalid heavy eigenmasses')
    hard=source.matrix(mf,y1,y2,lam)
    weights=np.asarray([complex(source.derivative(float(m),complex(lam),HBAR)) *
                        math.log(mu_f/float(m)) for m in mf],dtype=complex)
    relocation=pair_matrix(weights,y1,y2)
    grouped=neutrino.validated_c5(hard+relocation)
    return hard,relocation,grouped


def historical_frozen_parity(source: DirectWilsonSource, match: Path,
                             uv: Any, eft: Any, masses: np.ndarray,
                             y1: np.ndarray, y2: np.ndarray, lam: complex) -> float:
    """Numerically reproduce original saved frozen-coupling hard+log kernel."""
    from dataclasses import replace
    from scripts.FitT3BPlus1GeneralThreeStudies import MatchedKernel
    frozen=replace(eft,lambdaT3=complex(lam)).validated()
    legacy=neutrino.validated_c5(MatchedKernel(match).evaluate(uv,frozen))
    _hard,_relocation,grouped=grouped_hard_boundary(source,masses,y1,y2,lam,
                                                      mu_f=continuous.s34.START)
    predicted=grouped+math.log(continuous.s34.STOP/continuous.s34.START)*source.matrix(masses,y1,y2,lam,derivative=True)
    error=float(np.linalg.norm(predicted-legacy)/max(np.linalg.norm(legacy),1e-300))
    require(np.isfinite(error) and error<3e-8,
            'hierarchical hard+running parity failed: '+str(error))
    return error


def scalar_eft1_trajectory(model: str, record: dict, sidecar: dict, *, points: int=65) -> dict:
    """Continue the SAVED fermion-threshold state with original full scalar EFT1 betas.

    The continuous report stores only the new quartics at interior samples,
    not lambdaT3. Reintegrating EFT1 (NOT UV) is required to obtain the
    source term lambdaT3(t) without making up interpolation assumptions.
    """
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState, pack_general_scalar_state,
        unpack_general_scalar_state)
    from scripts import IntermediateEFTRunning as s34
    require(model in MODELS and 16<=points<=200, 'invalid model/point count')
    require(sidecar.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated'
            and sidecar.get('model')==model and sidecar.get('format')==continuous.qeft.FORMAT,
            'saved intermediate EFT1 endpoint incompatible')
    registry=continuous.load_registry(continuous.REGISTRY)
    base=s34.decode_intermediate_state(sidecar['initial_baseline'])
    extras={k:neutrino.complex_pair(v,k)
            for k,v in sidecar['initial_new_quartics'].items()}
    full=GeneralT3ScalarRGState(base,extras).validated(registry)
    y0,layout=pack_general_scalar_state(full,registry)
    kind=record['quadratic_kind']
    if kind!='none':
        mass=neutrino.complex_pair(sidecar['mu2_initial_gev2'],'initial mu2')
        y0=np.concatenate((y0,[mass.real,mass.imag]))
    beta_path=continuous.BETA_DIR/f'{continuous.model_tag(model)}_eft1_quadratic_live_rgbeta.json'
    betas=continuous.beta_inventory(continuous.read(beta_path),model=model,stage='EFT1',
             representation=base.representation,layout=layout,registry=registry,kind=kind)
    result=continuous.integrate_segment(y0,mu_start=s34.START,mu_end=s34.STOP,
            layout=layout,registry=registry,stage='EFT1',kind=kind,betas=betas,
            points=points)
    # Verify independently saved endpoint *before* using the re-run trajectory.
    final=unpack_general_scalar_state(result.y[:layout.size,-1],layout,registry,
                                     uv=False,mu_gev=s34.STOP)
    saved=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    chk={}
    for key in ('mS1Sq','mS2Sq','lambdaT3'):
        a=complex(getattr(final.baseline,key)); b=complex(getattr(saved,key))
        err=abs(a-b)/max(abs(b),1e-12)
        require(err<2e-6, model+': saved scalar endpoint parity failed for '+key)
        chk[key]=float(err)
    saved_extras=sidecar['final_new_quartics']
    require(set(saved_extras)==set(final.extra), 'extra quartic inventory parity failed')
    for name,val in final.extra.items():
        target=neutrino.complex_pair(saved_extras[name],name)
        require(abs(val-target)<=2e-6*max(abs(target),1e-10),
                model+': saved extra quartic endpoint parity failed: '+name)
    if kind!='none':
        target=neutrino.complex_pair(sidecar['mu2_final_gev2'],'final mu2')
        current=complex(result.y[layout.size,-1],result.y[layout.size+1,-1])
        require(abs(current-target)<=2e-6*max(abs(target),1e-10),
                model+': quadratic beta endpoint parity failed')
    lams=np.asarray([unpack_general_scalar_state(result.y[:layout.size,i],
             layout,registry,uv=False,mu_gev=math.exp(float(t))).baseline.lambdaT3
             for i,t in enumerate(result.t)],dtype=complex)
    require(lams.shape==(points,) and np.isfinite(lams).all(), 'nonfinite EFT1 lambdaT3 trajectory')
    integral=complex(simpson(lams,x=result.t))
    frozen=complex(lams[0])*(math.log(s34.STOP)-math.log(s34.START))
    return {'lambda_at_F':complex(lams[0]), 'lambda_at_S':complex(lams[-1]),
            'integral_lambda_dlnmu':integral,
            'frozen_lambda_integral_dlnmu':frozen,
            'trajectory_samples':points,'eft1_nfev':int(result.nfev),
            'endpoint_parity':chk, 'beta_path':beta_path}


def calculate(model: str, *,source_root: Path,output_root: Path,
              points: int=65,resume: bool=False) -> dict:
    require(model in MODELS,'unknown model')
    tag=continuous.model_tag(model)
    source_file=source_root/tag/'continuous_trajectory.json'
    record=neutrino.verify_continuous(model,source_file)
    sidecar_file=source_root/tag/record['saved_candidate_eft1_endpoint']
    sidecar=neutrino.load(sidecar_file)
    matching,wilson=model_sources(model)
    target=output_root/f'{tag}_direct_wilson.json'
    hashes={'continuous':sha(source_file),'endpoint':sha(sidecar_file),
            'runner':sha(Path(__file__).resolve()),
            'continuous_runner':sha(ROOT/'scripts/RunT3FullScalarContinuous.py'),
            'neutrino_runner':sha(ROOT/'scripts/RunT3FullScalarNeutrino.py')}
    if matching.is_file() and wilson.is_file():
        hashes['hierarchical_matching']=sha(matching)
        hashes['direct_wilson_beta']=sha(wilson)
    if target.is_file():
        existing=neutrino.load(target)
        require(resume and existing.get('format')==FORMAT and
                existing.get('input_sha256')==hashes and
                existing.get('saved_points')==points,
                'existing output or changed inputs: use --resume or new --output-id')
        return existing
    if not matching.is_file() or not wilson.is_file():
        result={'format':FORMAT,'status':SKIP,'model':model,
            'input_sha256':hashes,'saved_points':points,
            'reason':'Missing saved hierarchical MSbar hard threshold or EFT1 direct Weinberg beta; legacy one-scale C5 is not a valid substitute',
            'hierarchical_source_available':matching.is_file(),
            'wilson_beta_available':wilson.is_file(),
            'full_general_scalar_matching_certified':False,
            'full_sequential_wilson_matching_performed':False}
        neutrino.atomic_json(target,result)
        return result
    source=compile_direct_source(matching,wilson)
    uv=neutrino.decode_uv(record['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    mf,y1,y2,residual=heavy_inputs(uv)
    running=scalar_eft1_trajectory(model,record,sidecar,points=points)
    lam_f=running['lambda_at_F']
    hard,relocation,grouped=grouped_hard_boundary(source,mf,y1,y2,lam_f,
                                                      mu_f=continuous.s34.START)
    unit=source.matrix(mf,y1,y2,1+0j,derivative=True)
    direct=neutrino.validated_c5(unit*running['integral_lambda_dlnmu'])
    c5=neutrino.validated_c5(grouped+direct)
    frozen_direct=unit*running['frozen_lambda_integral_dlnmu']
    frozen_c5=neutrino.validated_c5(grouped+frozen_direct)
    parity=historical_frozen_parity(source,matching,uv,eft,mf,y1,y2,lam_f)
    low=neutrino.low_energy(model,c5,eft)
    low_frozen=neutrino.low_energy(model,frozen_c5,eft)
    delta=np.linalg.norm(c5-frozen_c5)/max(np.linalg.norm(frozen_c5),1e-300)
    threshold_hierarchy={
        'heavy_masses_above_fermion_threshold':bool(np.min(mf)>continuous.s34.START),
        'scalar_masses_above_scalar_endpoint':bool(min(eft.mS1Sq,eft.mS2Sq)>continuous.s34.STOP**2),
    }
    out={
        'format':FORMAT,'status':STATUS,'model':model,
        'input_sha256':hashes,'saved_points':points,
        'matching_source':matching.relative_to(ROOT).as_posix(),
        'wilson_source':wilson.relative_to(ROOT).as_posix(),
        'hard_matching_scheme':'saved MSbar per-eigenmass hard relocated to grouped fermion threshold',
        'direct_eft1_source_scheme':'one-loop tree-LLSS source, full-flavour lift; frozen UV Yukawas/MF at fermion threshold',
        'source_beta_16pi2_without_couplings':source.beta_text,
        'hard_kernel':source.hard_text,'running_kernel':source.running_text,
        'hard_c5_at_each_heavy_mass_gev_inverse':neutrino.matrix_payload(hard),
        'grouped_hard_scale_relocation_c5_gev_inverse':neutrino.matrix_payload(relocation),
        'c5_at_grouped_fermion_threshold_gev_inverse':neutrino.matrix_payload(grouped),
        'direct_wilson_c5_added_between_thresholds_gev_inverse':neutrino.matrix_payload(direct),
        'c5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(c5),
        'frozen_coefficient_c5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(frozen_c5),
        'lambdaT3_at_F':continuous.cjson(running['lambda_at_F']),
        'lambdaT3_at_S':continuous.cjson(running['lambda_at_S']),
        'lambdaT3_log_integral':continuous.cjson(running['integral_lambda_dlnmu']),
        'relative_C5_change_vs_frozen_lambda':float(delta),
        'frozen_hierarchical_kernel_parity_error':parity,
        'heavy_mass_basis_residual':residual,
        'eft1_endpoint_parity':running['endpoint_parity'],
        'eft1_source_transport_ode_evaluations':running['eft1_nfev'],
        'eft1_source_transport_samples':running['trajectory_samples'],
        'threshold_hierarchy_diagnostic':threshold_hierarchy,
        'uv_integration_repeated':False,
        'eft1_full_scalar_integration_repeated_for_wilson_source':True,
        'fermion_threshold_hard_C5_matched_from_saved_source':True,
        'full_fermion_threshold_LLSS_wilson_tensors_evolved':False,
        'scalar_threshold_full_one_loop_matching_performed':False,
        'additional_scalar_quartic_wilson_sources_included':False,
        'additional_scalar_bilinear_wilson_sources_included':False,
        'homogeneous_wilson_rge_mixing_included':False,
        'full_general_scalar_matching_certified':False,
        'full_sequential_wilson_matching_performed':False,
        'matching_compatible_with_one_scale_massdependent_six_model_table_certified':False,
        'note':('Source-derived, first sequential *direct* Wilson contribution. The EFT1 scalar background is fully RG-evolved, '
                'but LLSS tensor/matrix mixing, self bilinears and extra-quartic Wilson sources, '
                'finite scalar matching and cross-scheme harmonization are omitted. '
                'Do NOT identify this with the previous mass-dependent one-scale six-model C5 result.'),
        **low,
        'frozen_takagi_masses_ev_descending':low_frozen['takagi_masses_ev_descending'],
    }
    neutrino.atomic_json(target,out)
    return out


def summarize(source_root: Path, output_root: Path, selected: tuple[str,...]) -> dict:
    rows=[]
    for model in selected:
        path=output_root/f'{continuous.model_tag(model)}_direct_wilson.json'
        report=neutrino.load(path)
        require(report.get('format')==FORMAT and report.get('model')==model and
                report.get('full_sequential_wilson_matching_performed') is False,
                model+': incompatible report')
        masses=report.get('takagi_masses_ev_descending')
        row={'model':model,'status':report['status'],
             'direct_wilson_source_evolved':report['status']==STATUS,
             'max_mnu_ev':masses[0] if masses else None,
             'middle_mnu_ev':masses[1] if masses else None,
             'min_mnu_ev':masses[2] if masses else None,
             'relative_C5_change_vs_frozen_lambda':report.get('relative_C5_change_vs_frozen_lambda'),
             'frozen_hierarchical_kernel_parity_error':report.get('frozen_hierarchical_kernel_parity_error'),
             **report.get('oscillation_observables_normal_ordering',{})}
        rows.append(row)
    success=sum(x['status']==STATUS for x in rows)
    payload={'format':FORMAT,'status':'DirectWilsonSourceTransportSummary',
             'model_count':len(rows),'direct_source_models':success,
             'missing_hierarchical_source_models':len(rows)-success,
             'all_16_models_sequentially_matched':False,
             'full_sequential_matching_performed':False,
             'physical_full_general_scalar_matching_certified':False,
             'rows':rows}
    neutrino.atomic_json(output_root/'summary.json',payload)
    cols=['model','status','direct_wilson_source_evolved','max_mnu_ev','middle_mnu_ev',
          'min_mnu_ev','relative_C5_change_vs_frozen_lambda',
          'frozen_hierarchical_kernel_parity_error',
          'delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23']
    output_root.mkdir(parents=True,exist_ok=True)
    temp=output_root/'summary.csv.tmp'
    with temp.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
    temp.replace(output_root/'summary.csv')
    return payload


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    selector=parser.add_mutually_exclusive_group(required=True)
    selector.add_argument('--all',action='store_true')
    selector.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True,type=continuous.run_id)
    parser.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    parser.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    parser.add_argument('--samples',type=int,default=65)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--combine-only',action='store_true')
    args=parser.parse_args()
    continuous.run_id(args.output_id)
    require(args.output_id not in ('neutrino_v1','self_bilinear_v1',
              'self_anomalous_v1','self_anomalous_massdep_v1'),
            'output-id overlaps historical results')
    require(16<=args.samples<=200 and args.samples%2==1,
            '--samples must be an odd integer 17..199 for composite Simpson quadrature')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing continuous run '+str(folder))
    out=folder/args.output_id
    selected=MODELS if args.all else (args.model,)
    require(not args.combine_only or args.all, '--combine-only requires --all')
    failed=[]
    if not args.combine_only:
        for model in selected:
            try:
                report=calculate(model,source_root=folder,output_root=out,
                                 points=args.samples,resume=args.resume)
                if report['status']==STATUS:
                    print(f'{model}: direct hierarchical Wilson source OK; '
                          f'frozen-parity={report["frozen_hierarchical_kernel_parity_error"]:.3e}; '
                          f'delta(C5)/frozen={report["relative_C5_change_vs_frozen_lambda"]:.6g}; '
                          f'masses={report["takagi_masses_ev_descending"]}',flush=True)
                else:
                    print(f'{model}: {report["status"]} (legacy one-scale control only)',flush=True)
            except (ValueError,KeyError,TypeError,OSError,RuntimeError,OverflowError,
                    FloatingPointError,np.linalg.LinAlgError) as exc:
                failed.append(model)
                print(f'{model}: BLOCKED: {type(exc).__name__}: {exc}',flush=True)
    if failed:
        print('Blocked models:',', '.join(failed),flush=True)
        return 2
    if args.all:
        summary=summarize(folder,out,selected)
        print('Summary:',out/'summary.csv',
              'transported=',summary['direct_source_models'],
              'no_hierarchical_source=',summary['missing_hierarchical_source_models'],flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
