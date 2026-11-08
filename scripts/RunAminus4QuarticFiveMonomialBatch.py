"""Recompute A(-4) RGBeta/reference disagreement; produce NONPHYSICAL comparison inputs.

Re-evaluates the existing 55 independent scalar-quartic directions, polarises
all their native-minus-reference differences, and checks the exact five-
monomial identity. Two copies of the native global RGBeta 17-beta payload
are emitted: a byte-equivalent-compatible baseline and a DIAGNOSTIC candidate
whose scalar-only polynomial is made equal to the chosen 12-real-field
fourth-derivative reference. It is NOT a correction to the source RGBeta
package, and NOT a physical Matchete coupling or scalar-basis matching.

Run from the project root:
  python scripts/RunAminus4QuarticFiveMonomialBatch.py \
    --output-root output/a_minus4_next_batch/quartic_five_monomial
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from math import pi, log
from pathlib import Path
import sys
from typing import Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.running import Aminus4EFT1Audit as audit
from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.running.Aminus4QuarticDiscrepancy import (
    COORDINATES, DELTA_INPUTFORM, assert_expected_five, coordinates_of,
    discrepancy_vector, eval_polynomial, extract_monomials, monomials_to_json,
)
from scripts.VerifyAminus4GlobalQuarticProjection import compatibility_payload

RAW = Path('output/a_minus4_eft1_batch2/A_m4_projected_raw.json')
OLD_AUDIT = Path('output/a_minus4_eft1_batch2/A_m4_projection_audit.json')
BASIS = Path('output/a_minus4_eft1_batch2/A_m4_eft1_tensor_basis.json')
NATIVE_SPINOR = Path('output/general_observable/a_minus4_boundary_batch1/native_spinor_action.json')
MASTER_SPINOR = Path('output/general_observable/a_minus4_boundary_batch1/master_basis_candidate.json')
DECLARED_SPINOR = Path('output/general_observable/a_minus4_boundary_batch1/declared_convention_boundary.json')


def require(condition, reason: str) -> None:
    if not condition: raise ValueError(reason)


def load(root: Path, relative: Path):
    file = root/relative
    require(file.is_file(), f'Missing project source: {file}')
    raw = file.read_bytes()
    data = json.loads(raw.decode('utf-8-sig'))
    require(isinstance(data, dict), f'Expected JSON object: {file}')
    return data, sha256(raw).hexdigest()


def point(names: tuple[str, ...]) -> dict:
    from scripts.VerifyAminus4GlobalQuarticProjection import unit_point
    return unit_point(names)


def validate_spinor_status(root: Path) -> dict:
    native, nh = load(root, NATIVE_SPINOR)
    master, mh = load(root, MASTER_SPINOR)
    declared, dh = load(root, DECLARED_SPINOR)
    require(all(x.get('model')=='A_-4' for x in (native,master,declared)),
            'Source spinor witnesses describe different models')
    require(native.get('native_concrete_4x4_cc_matrix') is None and
            native.get('native_4x4_CC_to_Weyl_epsilon_ratio') is None and
            native.get('physical_certification') is False,
            'Native signed GammaCC status changed; direct re-audit required')
    require(master.get('physical_C_ijab_at_MF') is None and
            master.get('physical_certification') is False,
            'Master Wilson physical status changed; direct re-audit required')
    require(declared.get('physical_C_ijab_at_MF') is None and
            declared.get('native_matchete_to_declared_Weyl_factor') is None and
            declared.get('physical_certification') is False,
            'Declared Weyl source status changed; direct re-audit required')
    a=declared.get('declared_rho_plus_one_independent_coordinates', [])
    b=declared.get('rho_minus_one_fixed_native_y1_alternative_coordinates', [])
    require(len(a)==len(b)==36, 'Declared Weyl 36-coordinate inventory changed')
    def map_values(rows):
        ans={}
        for row in rows:
            key=tuple(row[k] for k in ('i','j','a','b'))
            require(key not in ans, 'Duplicated spinor candidate component')
            ans[key]=complex(row['value_re_gev_inv'],row['value_im_gev_inv'])
        return ans
    ap,bm=map_values(a),map_values(b)
    require(set(ap)==set(bm), 'Signed Weyl candidate supports differ')
    err=max(abs(ap[k]+bm[k]) for k in ap)
    size=max(abs(v) for v in ap.values())
    require(size>0 and err/size<1e-12, 'Signed Weyl candidate sign algebra changed')
    numeric=declared.get('numerical_checks',{}).get(
        'synthetic_full_three_family_UV_to_EFT_Grassmann_rephasings',{})
    require(numeric.get('gamma_cc_branch_physically_discriminated') is False and
            float(numeric.get('max_equivalent_rho_y1_branch_relative_residual',1))<1e-12,
            'UV rephasing evidence changed')
    return {
        'source_sha256':{'native_action':nh,'master_candidate':mh,'declared':dh},
        'provisional_master_support':len(ap),
        'two_fixed_native_y1_signed_branches_are_opposites_rel_max':float(err/size),
        'opposite_sign_with_simultaneously_flipped_native_y1_is_convention_equivalent':True,
        'native_signed_GammaCC_or_connected_amplitude_measured':False,
        'Matchete_to_master_Weyl_normalisation_physically_certified':False,
        'physical_C_ijab_at_MF':None,
        'required_for_physical_sign':'Installed-Matchete signed spinor matrix/3-point rule in a fixed Weyl representation or independent connected UV four-point amplitude with the same source couplings',
    }


def quadratic_impact(monomials: Mapping, basis: Mapping) -> dict:
    states=basis.get('running_couplings',[])
    scales=np.asarray(basis.get('sampled_scales_gev',[]),float)
    require(len(states)>=3 and len(states)==len(scales) and
            np.all(np.diff(scales)<0) and np.all(scales>0),
            'Expected saved descending A(-4) numerical running grid')
    delta=np.asarray([eval_polynomial(coordinates_of(s),monomials) for s in states])
    integral=np.zeros(10)
    for k in range(len(scales)-1):
        dt=log(scales[k+1]/scales[k])
        integral+=dt*(delta[k+1]+delta[k])/(2*16*pi*pi)
    return {
        'purpose':'One-insertion fixed-background comparison only, not rerun ODE or physical threshold',
        'scale_high_gev':float(scales[0]),'scale_low_gev':float(scales[-1]),
        'samples':len(scales),
        'definition':'integral of (native minus reference scalar-only beta)/(16*pi^2) dln(mu), evaluated on archived native reference trajectory with no feedback',
        'integrated_native_minus_reference_fixed_background':{
            n:float(z) for n,z in zip(COORDINATES,integral)},
        'boundary_delta_beta_16pi2':{n:float(z) for n,z in zip(COORDINATES,delta[0])},
        'endpoint_delta_beta_16pi2':{n:float(z) for n,z in zip(COORDINATES,delta[-1])},
        'not_certified_physics':True,
    }


def compare_random_points(payload: dict, monomials: Mapping, *, seed=20261108,
                          n=5, atol=2e-8) -> float:
    from scripts.VerifyAminus4GlobalQuarticProjection import unit_point
    rng=np.random.default_rng(seed)
    worst=0.0
    for _ in range(n):
        values=rng.uniform(-.5,.5,size=len(COORDINATES))
        p=unit_point(())
        for k,v in zip(COORDINATES,values):
            if k=='lambdaT3.real':p['lambdaT3']+=v
            elif k=='lambdaT3.imag':p['lambdaT3']+=1j*v
            else:p[k]=float(v)
        row=audit.scalar_beta_comparison(p,payload)
        err=np.max(np.abs(discrepancy_vector(row)-eval_polynomial(values,monomials)))
        worst=max(worst,float(err))
        require(err<atol,f'Five-monomial identity failed at random quartic probe: {err:.3e}')
    return worst


def build_diagnostic_beta(baseline: dict) -> dict:
    candidate=deepcopy(baseline)
    for name,delta in DELTA_INPUTFORM.items():
        candidate['report_betas'][name]=f"({baseline['report_betas'][name]}) - ({delta})"
    candidate['metadata'].update({
        'Aminus4ScalarOnlyReferenceDiagnostic':'native scalar-only 1-loop quartic polynomial minus five verified discrepancies',
        'Aminus4PhysicalQuarticBetaCertified':False,
        'Aminus4MatcheteInvariantMappingCertified':False,
        'Aminus4GammaCCToWeylNormalizationCertified':False,
        'Aminus4DoNotTreatAsProductionOrPhysical':True,
    })
    eft.validate_beta_payload(candidate)
    return candidate


def validate_reference_candidate(candidate: dict) -> dict:
    failures=[]
    max_closure=0.
    max_coeff_error=0.
    from itertools import combinations
    probes=[(k,) for k in COORDINATES]+list(combinations(COORDINATES,2))
    for names in probes:
        row=audit.scalar_beta_comparison(point(names),candidate)
        max_closure=max(max_closure,float(row['scalar_loop_projection_closure_relative']))
        max_coeff_error=max(max_coeff_error,float(np.max(abs(discrepancy_vector(row)))))
        if not row['reference_and_native_scalar_loop_match']:
            failures.append('+'.join(names))
    require(not failures and max_closure<1e-9,
            f'Diagnostic candidate fails independent real-scalar quadratic comparison: {failures}')
    return {'probe_count':len(probes),'remaining_disagreements':failures,
            'maximum_abs_scalar_beta_difference':max_coeff_error,
            'maximum_reference_projection_closure_relative':max_closure}


def evaluate(root: Path, output_root: Path) -> dict:
    output_root=output_root.resolve()
    paths={
        'report':output_root/'five_monomial_report.json',
        'baseline':output_root/'native_global_projection_compatible.json',
        'candidate':output_root/'scalar_only_reference_DIAGNOSTIC_beta.json',
    }
    for key,path in paths.items():
        require(not path.exists(), f'Refusing to overwrite existing {key}: {path}')
    raw, rh=load(root,RAW)
    old, ah=load(root,OLD_AUDIT)
    basis, bh=load(root,BASIS)
    require(old.get('new_global_beta_export_sha256')==rh and
            old.get('model')==basis.get('model')=='A_-4',
            'Global raw beta is not the audited source; investigate changed files')
    require(old.get('scalar',{}).get('scalar_quadratic_probe_count')==55,
            'Archived RGBeta quadratic audit inventory changed')
    require(len(old['scalar']['failures'])==34,
            'Archived RGBeta disagreement count changed; audit directly')
    compatible=compatibility_payload(raw)
    eft.validate_beta_payload(compatible)

    def probe(names):return audit.scalar_beta_comparison(point(names),compatible)
    extraction=extract_monomials(probe)
    require(extraction['probe_count']==55 and extraction['failing_probes']==34 and
            extraction['max_scalar_reference_closure_relative']<1e-9 and
            extraction['max_probed_reconstruction_abs']<2e-8,
            'Source RGBeta quadratic decomposition is not the known A(-4) one')
    coeff_err=assert_expected_five(extraction['monomials'])
    random_error=compare_random_points(compatible,extraction['monomials'])

    candidate=build_diagnostic_beta(compatible)
    candidate_check=validate_reference_candidate(candidate)
    # At generic complex phases, the Bar[lambdaT3]*lambdaT3 replacement is
    # independently tested by the actual InputForm evaluator, not assumed.
    candidate_random_error=compare_random_points(candidate,{},n=5)
    spinor=validate_spinor_status(root)
    impact=quadratic_impact(extraction['monomials'],basis)
    report={
        'format':'t3-a-minus4-five-monomial-quartic-reference-DIAGNOSTIC-v1',
        'model':'A_-4','status':'DIAGNOSTIC_SCALAR_ONLY_POLYNOMIAL_RECONCILED_PHYSICAL_BASIS_OPEN',
        'input_sha256':{'global_rgbeta_projection':rh,'archived_audit':ah,
                        'archived_running_basis':bh},
        'scalar_only':{
            'native_minus_real_scalar_reference':True,
            'independent_unit_and_pair_probe_count':55,
            'failed_unit_pair_probes_before':extraction['failing_probes'],
            'unaffected_passing_probes_before':55-extraction['failing_probes'],
            'actual_independent_bad_monomial_count':len(extraction['monomials']),
            'polynomial_terms':monomials_to_json(extraction['monomials']),
            'max_observed_vs_expected_exact_coeff_error':coeff_err,
            'max_55_point_reconstruction_abs':extraction['max_probed_reconstruction_abs'],
            'max_five_random_points_reconstruction_abs':random_error,
            'candidate_55_probes_against_reference':candidate_check,
            'candidate_five_random_points_vs_reference_max_abs':candidate_random_error,
            'diagnostic_correction_inputform':DELTA_INPUTFORM,
            'interpretation':'These discrepancies can arise from a native scalar tensor/group-index convention issue; no package bug or physical basis correspondence is asserted',
        },
        'fixed_background_impact':impact,
        'fermion_spinor_certification':spinor,
        'produced':{name:str(p) for name,p in paths.items()},
        'physical_quartic_beta_certified':False,
        'physical_EFT1_Wilson_boundary_certified':False,
        'production_beta_files_modified':False,
        'EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars together -> SMEFT',
    }
    output_root.mkdir(parents=True, exist_ok=True)
    for name, obj in [('baseline',compatible),('candidate',candidate),('report',report)]:
        paths[name].write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--output-root',type=Path,
                   default=Path('output/a_minus4_next_batch/quartic_five_monomial'))
    args=p.parse_args(argv)
    result=evaluate(args.root.resolve(),args.output_root)
    scalar=result['scalar_only']
    print('Native/reference scalar-only probe disagreements:',
          scalar['failed_unit_pair_probes_before'],'/55')
    print('Independent discrepant quadratic monomials:',
          scalar['actual_independent_bad_monomial_count'])
    print('Diagnostic candidate remaining scalar-only disagreements:',
          len(scalar['candidate_55_probes_against_reference']['remaining_disagreements']))
    print('Native physical quartic map: NOT CERTIFIED')
    print('Matchete GammaCC/Weyl factor: NOT CERTIFIED')
    print('Audit:',result['produced']['report'])
    return 0


if __name__=='__main__':
    raise SystemExit(main())
