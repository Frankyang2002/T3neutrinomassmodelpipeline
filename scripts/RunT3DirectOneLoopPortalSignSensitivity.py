"""All-16, source-guarded, *conditional* Higgs-portal sign sensitivity.

Project path: scripts/RunT3DirectOneLoopPortalSignSensitivity.py

Frozen lambdaH{1,2} are treated as the baseline RGBeta (MBL) couplings.
RGBeta's model-builder quartic Lagrangian uses a minus sign, whereas the
Matchete portal builder uses +1/2 PlusHc. With the all-32 native identity
span already verified, the candidate signed dictionary is
lambdaH{j}Inv{k}_Matchete = -c_{jk} lambdaH{j}_RGBeta.

The existing direct-one-loop adapter uses A_ii = mSiSq - lambdaHi*v**2/2.
For this conditional dictionary the shift is +lambdaHi*v**2/2, so the
variant is A_variant = A_baseline + lambdaHi*v**2 * Identity for each Si.

THIS IS A SENSITIVITY STUDY, NOT A CERTIFIED PHYSICAL CONVERSION:
fermionic phases, Yukawa normalization, lambdaT3 sign, and general
symmetry-allowed vertices remain open. Only the portal diagonal is changed.
No EFT code, frozen benchmark, prior UV input, or prior UV result is written.
"""
from __future__ import annotations

import argparse
import copy
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.direct_one_loop.core import encoded, evaluate  # noqa: E402
from scripts.RunT3DirectOneLoop import models, CLASSES  # noqa: E402

OUT_FORMAT = 't3-direct-one-loop-portal-sign-sensitivity-v1'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(path: Path) -> dict:
    require(path.is_file(), f'Missing input: {path}')
    value = json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(value, dict), f'Malformed JSON: {path}')
    return value


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_contract(root: Path) -> dict:
    """Fail closed if the relevant Matchete or RGBeta field words change.

    RGBeta L_quartic=-lambda*invariant is a MBL convention described in
    Thomsen, EPJC 81 (2021) 408, Secs. 2 and 4; this repository check
    guards the specific tensor and field-word declarations only.
    """
    mat = root/'Lagrangian/interactions/ScalarPotential.wl'
    rgb = root/'RGE/running/rgbeta/T3RGBetaModel.wl'
    m = mat.read_text(encoding='utf-8-sig')
    r = rgb.read_text(encoding='utf-8-sig')
    first = m.find('BuildHiggsPortalInvariantCandidates[')
    last = m.find('BuildCrossScalarInvariantCandidates[', first)
    require(first >= 0 and last > first, 'Matchete portal function absent')
    fn = m[first:last]
    for needle in ('1/2 PlusHc[', 'coupling[] fieldProduct',
                   '{True, False, True, False}',
                   '"lambdaH" <> ToString[which]'):
        require(needle in fn, f'Matchete portal source contract changed: {needle}')
    for j in (1, 2):
        pat = (r'AddQuartic\s*\[\s*lambdaH'+str(j)+
               r'\s*,\s*\{\s*Bar\s*@\s*H\s*,\s*H\s*,\s*Bar\s*@\s*S'+str(j)+
               r'\s*,\s*S'+str(j)+r'\s*\}')
        require(bool(re.search(pat, r)), f'RGBeta lambdaH{j} field-word contract changed')
    return {str(mat.relative_to(root)): digest(mat),
            str(rgb.relative_to(root)): digest(rgb)}


def complex_matrix(obj: dict) -> np.ndarray:
    require(isinstance(obj, dict) and set(obj)=={'real', 'imag'},
            'Expected encoded complex matrix')
    return np.asarray(obj['real'], float)+1j*np.asarray(obj['imag'], float)


def make_variant(source: dict, frozen: dict, identity_row: dict, *, tolerance: float=1e-11) -> dict:
    """Check baseline diagonal convention, then change only its portal sign."""
    rep = frozen['representation']
    d1, d2 = int(rep['d_s1']), int(rep['d_s2'])
    expected_model = identity_row['model']
    require(source.get('model') == expected_model, 'Direct input model mismatch')
    require(len(identity_row.get('scalars', [])) == 2, 'Missing native portal sectors')
    for j, dim in ((1,d1),(2,d2)):
        row = next((z for z in identity_row['scalars'] if z['scalar']==f'S{j}'), None)
        require(row is not None and row['dimension']==dim, 'Native portal dimension mismatch')
        require(row['full_tensor_identity_in_real_native_span'] is True,
                'Native four-index identity did not pass')
        require(float(row['max_abs_four_index_residual']) <= 2e-10,
                'Nonzero native identity residual')
        require(row['native_invariant_count']==len(row['native_coupling_coefficients_for_unit_identity']),
                'Native identity coefficients incomplete')
    ordinary = frozen['base_state']['ordinary']
    vev = float(frozen['numerical']['vev_gev'])
    require(vev > 0 and np.isfinite(vev), 'Invalid benchmark Higgs VEV')
    lam1 = float(ordinary['lambdaH1']); lam2 = float(ordinary['lambdaH2'])
    require(all(np.isfinite(x) for x in (lam1,lam2)), 'Nonfinite portal coupling')
    original = complex_matrix(source['scalar_A_gev2'])
    require(original.shape==(d1+d2,d1+d2), 'Wrong scalar mass matrix shape')
    require(np.linalg.norm(original-original.conj().T) < 1e-10*max(np.linalg.norm(original),1),
            'Scalar A is not Hermitian')
    # The original adapter's ordinary scalar mass diagonal is exactly
    # mSiSq-lambdaHi*v^2/2, not merely an inferred prescription.
    baseline_diag = np.array([float(ordinary['mS1Sq'])-lam1*vev**2/2]*d1 +
                             [float(ordinary['mS2Sq'])-lam2*vev**2/2]*d2)
    mismatch = max(abs(original.diagonal()-baseline_diag))
    require(mismatch < tolerance*max(1,float(np.max(abs(baseline_diag)))),
            'Existing scalar A does not match adapter universal-portal baseline')
    diag_shift = np.array([lam1*vev**2]*d1+[lam2*vev**2]*d2)
    result = copy.deepcopy(source)
    updated = original+np.diag(diag_shift)
    require(np.allclose(updated-np.diag(updated.diagonal()),
                        original-np.diag(original.diagonal()), atol=0,rtol=0),
            'Unexpected offdiagonal scalar modification')
    result['scalar_A_gev2'] = encoded(updated)
    prov = result.setdefault('provenance', {})
    require(prov.get('validated_uv_component_conventions') is False,
            'Do not downgrade a previously certified source input')
    prov['portal_sign_sensitivity_only'] = True
    prov['portal_lagrangian_map_conditional'] = (
        'lambdaH{j}Inv{k}_Matchete = -c_{jk} lambdaH{j}_RGBeta; '
        'other UV interaction sign/phase mappings not established')
    prov['portal_diagonal_shift_relative_to_original_gev2'] = diag_shift.tolist()
    prov['native_identity_coefficients'] = {
        row['scalar']:row['native_coupling_coefficients_for_unit_identity']
        for row in identity_row['scalars']}
    return result


def run(root: Path, inp: Path, baseline: Path, identity: Path, out: Path) -> dict:
    source_hashes = source_contract(root)
    report = read_json(identity)
    require(report.get('format')=='t3-direct-uv-portal-offshell-identity-span-v1'
            and report.get('models_checked')==16
            and report.get('full_tensor_identities_representable')==32
            and report.get('full_tensor_identities_blocked')==0,
            'All-32 source native portal identity audit must pass')
    table = {v['model']:v for v in report['models']}
    require(len(table)==16, 'Duplicate or missing source identity model')
    prior = read_json(baseline/'all_16_status.json')
    prior_rows = {r['model']:r for r in prior['models']}
    require(len(prior_rows)==16, 'Prior source mass inventory incomplete')
    rows=[]
    variants=[]
    for short, name in models():
        fpath = root/'Numerical/benchmarks/fitted_reference'/f'{name}.json'
        ipath = inp/f'{short}.json'
        frozen = read_json(fpath)
        source = read_json(ipath)
        rep=frozen['representation']
        expected=(*CLASSES[short[0]][:3],int(short[2:]))
        require(tuple(rep[k] for k in ('d_s1','d_s2','d_f','alpha'))==expected,
                f'Benchmark metadata mismatch: {short}')
        require(source['provenance']['benchmark_sha256']==digest(fpath),
                f'Frozen input source SHA mismatch: {short}')
        require(table[short]['frozen_sha256']==digest(fpath),
                f'Native portal evidence benchmark SHA mismatch: {short}')
        old=prior_rows[short]
        require(old['status']=='uv_mass_evaluated_uncertified_source_mapping',
                f'Existing mass input is not provisional: {short}')
        oldfile = baseline/old['output_file']
        olddata = read_json(oldfile)
        require(olddata['model']==short, 'Prior physical output model mismatch')
        variant=make_variant(source,frozen,table[short])
        candidate=evaluate(variant,retain_individual_terms=False)
        require(candidate['diagnostics']['uv_finite'] is True,
                f'Variant UV pole failed: {short}')
        masses0=np.asarray(olddata['total']['neutrino_masses_ev'],float)
        masses1=np.asarray(candidate['total']['neutrino_masses_ev'],float)
        require(masses0.shape==(3,) and masses1.shape==(3,),
                f'Malformed mass spectrum: {short}')
        mat0=complex_matrix(olddata['total']['complex_symmetric_matrix_gev'])
        mat1=complex_matrix(candidate['total']['complex_symmetric_matrix_gev'])
        frob0=float(np.linalg.norm(mat0)); frob1=float(np.linalg.norm(mat1))
        require(frob0>0, f'Zero baseline mass: {short}')
        row={'model':short,
             'baseline_masses_ev':masses0.tolist(),
             'portal_signed_variant_masses_ev':masses1.tolist(),
             'relative_matrix_change_frobenius':float(np.linalg.norm(mat1-mat0)/frob0),
             'relative_m3_change':float((masses1[-1]-masses0[-1])/masses0[-1]),
             'baseline_matrix_frobenius_gev':frob0,
             'variant_matrix_frobenius_gev':frob1,
             'variant_relative_uv_pole':candidate['diagnostics']['relative_uv_pole_coefficient'],
             'physical_source_mapping_certified':False}
        rows.append(row)
        variants.append((short,variant,candidate))
        print(f'{short:6s} m3: {masses0[-1]:.8g} -> {masses1[-1]:.8g} eV; '
              f'delta_m3={100*row["relative_m3_change"]:.5g}%; '
              f'delta_matrix={100*row["relative_matrix_change_frobenius"]:.5g}%',flush=True)
    require(len(rows)==16, 'Did not evaluate all 16 models')
    # Write ONLY after every model has passed all checks.
    out.mkdir(parents=True,exist_ok=True)
    for short,variant,candidate in variants:
        (out/f'{short}_portal_signed_input.json').write_text(
            json.dumps(variant,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        (out/f'{short}_portal_signed_one_loop.json').write_text(
            json.dumps(candidate,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary={'format':OUT_FORMAT,
             'status':'16_of_16_portal_sign_sensitivity_evaluated_physical_mapping_pending',
             'source_hashes':source_hashes,
             'native_identity_audit_sha256':digest(identity),
             'assumption':'RGBeta MBL Lquartic=-lambda*(Hdag H)(Sdag S); '
                          'Matchete portal Lquartic=+1/2 PlusHc[lambdaInvk T_k]; '
                          'lambda_Matchete=-c*lambda_RGBeta under native full-tensor span',
             'only_change':'scalar_A diagonal -lambdaH*v^2/2 -> +lambdaH*v^2/2',
             'no_eft_matching_rge_threshold_performed':True,
             'physical_source_mapping_certified':False,
             'models':rows}
    (out/'all_16_portal_sign_sensitivity.json').write_text(
        json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (out/'all_16_portal_sign_sensitivity.csv').open('w',newline='',encoding='utf8') as fp:
        wr=csv.writer(fp)
        wr.writerow(['model','baseline_m1_ev','baseline_m2_ev','baseline_m3_ev',
                     'variant_m1_ev','variant_m2_ev','variant_m3_ev',
                     'relative_m3_change','relative_matrix_change_frobenius',
                     'variant_relative_uv_pole'])
        for row in rows:
            wr.writerow([row['model'],*row['baseline_masses_ev'],
                         *row['portal_signed_variant_masses_ev'],
                         row['relative_m3_change'],
                         row['relative_matrix_change_frobenius'],
                         row['variant_relative_uv_pole']])
    print(f'Portal sign sensitivity evaluated: {len(rows)}/16 (physical mapping still uncertified)')
    print('Output:',out/'all_16_portal_sign_sensitivity.csv')
    return summary


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=ROOT)
    p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--baseline-dir',type=Path,required=True)
    p.add_argument('--identity-audit',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args()
    try:
        run(a.project_root.resolve(),a.input_dir.resolve(),a.baseline_dir.resolve(),
            a.identity_audit.resolve(),a.output_dir.resolve())
    except (KeyError,ValueError,TypeError,IndexError,OSError,RuntimeError,ArithmeticError) as e:
        p.exit(2,'BLOCKED: '+str(e)+'\n')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
