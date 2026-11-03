#!/usr/bin/env python3
"""B1 read-only, fail-closed 16-benchmark source-spurion phase transport.

This script is a conditional *source-action* numerator/phase audit, not a
physical T3 loop calculation. It reads the previously executed B1 witnesses;
all output goes to a new directory which must not exist.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.lo_certification.physical_normalisation.b1_phase_transport import (  # noqa: E402
    PhaseAuditError, require, exact_invariant_witness, sign_orbit_witness,
    complex_matrix, complex_scalar, source_flavour_numerator,
    test_flavour_covariance, relative,
)
from Numerical.lo_certification.physical_normalisation.b1_matchete_weyl import model_specs  # noqa: E402
from Numerical.lo_certification.physical_normalisation.b1_native_action_contract import verify_record  # noqa: E402
from scripts.RunPhysicalB1NativeCoefficientAudit import (  # noqa: E402
    check_all16, check_cc, load, sha, MODELS,
)
from scripts.RunPhysicalB1WeylAudit import guard_source, git_revision  # noqa: E402

STATUS = 'All16SourceSpurionFlavourAndPhaseCovarianceVerified_PhysicalDictionaryPending'
ACTION_STATUS = 'ThreeNativeActionsCoefficientAndChiralityVerified_AbsoluteCCAndB2Pending'
BENCHMARK = Path('Numerical/benchmarks/fitted_reference')
FROZEN = Path('output/general_observable/lo_direct_one_loop/provisional_inputs')


def source_contract(root: Path, action_path: Path, prior_path: Path,
                    cc_path: Path, native_dir: Path) -> dict:
    inventory = guard_source(root)
    action = load(action_path)
    require(action.get('format')=='t3-b1-native-source-factor-audit-v2' and
            action.get('status')==ACTION_STATUS and
            action.get('native_actions_structurally_verified_count')==3 and
            action.get('physical_source_certified') is False and
            action.get('safe_to_integrate_into_production') is False,
            'Previous native coefficient audit not accepted or overclaims physics')
    require(action.get('previous_conditional_audit_sha256')==sha(prior_path),
            'Prior 16-model evidence differs from native coefficient audit provenance')
    previous = load(prior_path)
    summary=check_all16(previous)
    require(action.get('previous_conditional_16')==summary,
            'Conditional residuals/model inventory inconsistent with signed native report')
    a_inventory=action.get('guarded_source_blobs')
    require(isinstance(a_inventory,dict) and set(a_inventory)==set(inventory),
            'Native report source inventory differs from current guard')
    for path, rec in inventory.items():
        require(rec.get('github_main_blob_match') is True and
                a_inventory[path].get('git_blob_sha1')==rec.get('git_blob_sha1'),
                'Local source drift since native action audit: '+path)
    witnesses=action.get('native_action_witnesses')
    require(isinstance(witnesses,dict) and set(witnesses)==set(MODELS),
            'Native action witnesses incomplete')
    native_hashes={}
    for model in MODELS:
        p=native_dir/(model+'.json')
        actual=verify_record(load(p),model)
        record=witnesses[model]
        require(record.get('sha256')==sha(p) and
                record.get('couplings_and_chiral_projectors')==actual['couplings_and_chiral_projectors'],
                'Native action replaced or structural evidence drift: '+model)
        native_hashes[model]=sha(p)
    cc=load(cc_path)
    abstract=check_cc(cc)
    require(action.get('native_cc_abstract_probe')==abstract,
            'Native abstract GammaCC witness changed or was not included in prior audit')
    require(action.get('chosen_basis_algebra',{}).get('absolute_native_GammaCC_phase_measured') is False,
            'Previous basis calculation unexpectedly claims concrete Matchete phase')
    return {'current_guarded_source_blobs':{k:v['git_blob_sha1'] for k,v in inventory.items()},
            'native_action_audit_sha256':sha(action_path),
            'all16_conditional_audit_sha256':sha(prior_path),
            'native_abstract_cc_sha256':sha(cc_path),
            'native_actions_sha256':native_hashes,
            'previous_conditional_max_residuals': {
                k:summary[k] for k in ('max_mass_relative_residual',
                                        'max_Y_relative_residual',
                                        'max_Z_relative_residual')},
            'previous_conditional_count':summary['conditional_models']}


def benchmark_filename(dims: tuple[int,int,int,int]) -> str:
    d1,d2,df,alpha=dims
    suffix=f'{"m" if alpha<0 else "p"}{abs(alpha)}'
    return f'T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{suffix}.json'


def run_benchmarks(root: Path, prior: dict) -> list[dict]:
    records=prior['models']
    byname={rec['model']:rec for rec in records}
    require(len(byname)==16,'Duplicate model in conditional audit')
    output=[]
    for i,(model,dims) in enumerate(model_specs()):
        p=root/BENCHMARK/benchmark_filename(dims)
        frozen=root/FROZEN/(model+'.json')
        source=byname[model]
        require(sha(p)==source.get('benchmark_sha256'),
                model+': benchmark changed since conditional B1 run')
        require(sha(frozen)==source.get('frozen_input_sha256'),
                model+': frozen component input changed since conditional B1 run')
        bench=load(p)
        rep=bench['representation']
        require(tuple(int(rep[key]) for key in ('d_s1','d_s2','d_f','alpha'))==dims,
                model+': benchmark representation mismatch')
        o=bench['base_state']['ordinary']
        y1=complex_matrix(o['y1'],'y1')
        y2=complex_matrix(o['y2'],'y2')
        mf=complex_matrix(o['MF'],'MF')
        lam=complex_scalar(o['lambdaT3'],'lambdaT3')
        majorana=dims[2]%2==1 and dims[3]==-1
        if majorana:
            require(relative(mf,mf.T)<1e-11,model+': Majorana benchmark MF not symmetric')
        p0=source_flavour_numerator(y1,y2,mf,lam)
        require(float(abs(p0).max())>0.,model+': zero spurion numerator; sign control inconclusive')
        # Paired signs preserve the source spurion; a one-vertex sign does not.
        paired=source_flavour_numerator(-y1,-y2,mf,lam)
        single=source_flavour_numerator(-y1,y2,mf,lam)
        masspair=source_flavour_numerator(-y1,y2,-mf,lam)
        scalpair=source_flavour_numerator(-y1,y2,mf,-lam)
        require(relative(paired,p0)<1e-12 and relative(masspair,p0)<1e-12 and
                relative(scalpair,p0)<1e-12 and relative(single,-p0)<1e-12,
                model+': sign orbit inconsistency')
        cov=test_flavour_covariance(y1,y2,mf,lam,majorana,seed=301+i)
        output.append({'model':model,'fermion_type':'Majorana' if majorana else 'Dirac',
                       'source_benchmark_sha256':sha(p),'frozen_candidate_sha256':sha(frozen),
                       'source_spurion_numerator_frobenius':cov['source_spurion_numerator_frobenius'],
                       'U3_flavour_covariance_residual':cov['U3_flavour_covariance_residual'],
                       'symmetric_nu_pair_covariance_residual':cov['symmetric_nu_pair_covariance_residual'],
                       'double_yukawa_sign_rephase_verified':True,
                       'single_yukawa_sign_flips_source_spurion':True,
                       'mass_and_yukawa_coherent_sign_verified':True,
                       'scalar_mixing_and_yukawa_coherent_sign_verified':True,
                       'physical_amplitude_computed':False,
                       'physical_source_certified':False})
    require(len(output)==16,'Not all sixteen restricted T3 benchmarks observed')
    return output


def build(root: Path, action_path: Path, prior_path: Path, cc_path: Path,
          native_dir: Path, b2_path: Path | None=None) -> dict:
    provenance=source_contract(root,action_path,prior_path,cc_path,native_dir)
    previous=load(prior_path)
    dirac=exact_invariant_witness(False)
    majorana=exact_invariant_witness(True)
    orbit=sign_orbit_witness()
    benchmark_results=run_benchmarks(root,previous)
    b2={'supplied':False,'relative_SU2_physical_phase_measured':False}
    if b2_path is not None:
        rec=load(b2_path)
        require(rec.get('physical_source_certified') is False and
                rec.get('physical_crosspackage_complete') is False,
                'B2 evidence unexpectedly claims completed physical source map')
        b2={'supplied':True,'sha256':sha(b2_path),
            'reported_status':rec.get('status'),
            'relative_SU2_physical_phase_measured':False,
            'not_merged_into_B1_physical_certificate':True}
    return {
      'format':'t3-b1-phase-transport-source-spurion-v1',
      'status':STATUS,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
      'local_git_revision':git_revision(root), 'provenance':provenance,
      'exact_source_monomial_phase_witnesses':{'Dirac':dirac,'Majorana':majorana},
      'sign_rephasing_orbit':orbit,
      'models':benchmark_results,
      'source_flavour_benchmarks_checked_count':len(benchmark_results),
      'maximum_numerator_covariance_residual':max(x['U3_flavour_covariance_residual'] for x in benchmark_results),
      'maximum_symmetrised_covariance_residual':max(x['symmetric_nu_pair_covariance_residual'] for x in benchmark_results),
      'B2_evidence':b2,
      'absolute_native_GammaCC_phase_measured':False,
      'physical_crosspackage_SU2_dictionary_verified':False,
      'non_diagonal_MF_native_source_verified':False,
      'full_one_loop_neutrino_mass_computed_here':False,
      'physically_normalized_loop_amplitude_verified':False,
      'relative_physical_factors_verified':False,
      'physical_B1_scope_complete':False,
      'physical_source_certified':False,
      'physical_source_certified_models':0,
      'safe_to_integrate_into_production':False,
      'running_or_EFT_performed':False,
      'interpretation':[
          'The source rephasing invariant is a conditional source-coupling numerator, NOT a computed neutrino mass.',
          'For non-diagonal 3x3 heavy masses the U(3) result is an exact algebraic extension, not a native source validation.',
          'Global neutrino phase produces expected nu-nu congruence; fixed external neutrino basis gives an internal invariant.',
          'A single Yukawa sign change is not an internal rephasing of the four source coefficients if the other three remain fixed.',
          'Concrete installed GammaCC sign and RGBeta CG phase still require source-linked physical measurement.'
      ],
    }


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--action-audit',type=Path,required=True)
    p.add_argument('--conditional-audit',type=Path,required=True)
    p.add_argument('--cc-json',type=Path,required=True)
    p.add_argument('--native-dir',type=Path,required=True)
    p.add_argument('--b2-ward-json',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args()
    root=args.project_root.resolve()
    def resolve(path:Path|None)->Path|None:
        return None if path is None else (path if path.is_absolute() else root/path)
    out=resolve(args.output_dir)
    if out.exists():
        print('BLOCKED: refusing to overwrite '+str(out),file=sys.stderr);return 3
    try:
        report=build(root,resolve(args.action_audit),resolve(args.conditional_audit),
                     resolve(args.cc_json),resolve(args.native_dir),resolve(args.b2_ward_json))
    except (PhaseAuditError,ValueError,TypeError,KeyError,OSError) as exc:
        print('BLOCKED: '+str(exc),file=sys.stderr);return 2
    out.mkdir(parents=True,exist_ok=False)
    (out/'phase_transport_all16.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    (out/'HANDOFF.md').write_text(
        '# B1 source-spurion phase transport (conditional)\n\n'
        '- 16/16 benchmark coupling arrays passed independent U(3) heavy-flavour, U(3) external-neutrino, and scalar U(1) covariance.\n'
        f'- Maximum relative flavour numerator covariance residual: `{report["maximum_numerator_covariance_residual"]:.3e}`.\n'
        '- Unique internal invariant exponent: `(1,1,-1,-1)` for `y1,y2,MF,lambdaT3`.\n'
        '- Both Yukawa signs can flip together; flipping one alone while keeping the other three coefficients fixed changes the source invariant.\n'
        '- The benchmark source numerator is not a one-loop neutrino mass; no loop factor, scalar diagonalisation, or RGBeta matching is certified.\n'
        '- Native abstract GammaCC orientation and B2 physical SU(2) dictionary still unresolved. Physical certifications **0/16**.\n'
        '- No production sources, benchmarks, frozen inputs, GitHub, EFT, or RGE altered.\n',
        encoding='utf-8')
    print(STATUS+'; source_benchmarks=16/16; physical_source_certified=FALSE')
    print('max covariance residual:',report['maximum_numerator_covariance_residual'])
    print('Evidence:',out/'phase_transport_all16.json')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
