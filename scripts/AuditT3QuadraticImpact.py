"""Batch 3: audit scalar quadratic spectra and leading-order omission limits.

READ-ONLY with respect to the existing scalar/RGE/matching calculations.
Consumes saved candidate 16-model quadratic EFT1 reports, checks their reported
quadratic spectra against independent closed-form mass eigenvalues, and
compares mu2 != 0 against mu2=0 at IDENTICAL diagonal masses.

The comparison is NOT a new RGE run, C5 matching, or a neutrino-mass prediction.
A small scalar eigenvalue shift is not evidence of a small C5 shift.

Usage from repository root:
  python -m scripts.AuditT3QuadraticImpact
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path('docs/general_scalar_trial_operator_manifest.json')
BATCH2 = Path('output/scalar_relevance/leading_c5_vertex_audit.json')
SOURCE = Path('output/general_t3_quadratic_eft1')
DEST = Path('output/scalar_relevance')
MODELS = ('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3','C_+1','C_-1','C_-3',
          'D_+0','D_+2','D_-2','E_+0','E_+2','E_-2','E_-4')
REPS = {'A':(1,3), 'B':(2,2), 'C':(2,2),'D':(3,1),'E':(3,3)}
KINDS = {'self_s1': 'muT3GeneralS1Sq', 'self_s2': 'muT3GeneralS2Sq',
         'mixed': 'muT3GeneralS1S2', 'none': None}


class AuditError(ValueError):
    """Unverified or contradictory saved physics inputs."""


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise AuditError(f'Missing input: {path}')
    data = json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(data, dict):
        raise AuditError(f'Expected a JSON object: {path}')
    return data


def finite(value, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise AuditError(f'{label}: expected numeric value')
    v = float(value)
    if not math.isfinite(v):
        raise AuditError(f'{label}: nonfinite numeric value')
    return v


def complex_value(value, label: str) -> complex:
    if not isinstance(value, dict) or set(value) != {'real','imag'}:
        raise AuditError(f'{label}: requires a real/imag complex representation')
    return complex(finite(value['real'], label+'.real'),
                   finite(value['imag'], label+'.imag'))


def kind_from_manifest(record: dict, model: str) -> str:
    structures = record.get('quadratic_structures_not_enabled')
    count = record.get('quadratic_complex_couplings_not_enabled')
    if not isinstance(structures,list) or len(structures)>1 or count!=len(structures):
        raise AuditError(f'{model}: invalid quadratic manifest cardinality')
    if not structures:
        return 'none'
    structure = structures[0]
    counts = structure.get('counts')
    if not isinstance(counts,list) or len(counts)!=6 or sum(counts)!=2 or \
            not all(isinstance(n,int) and n>=0 for n in counts) or \
            structure.get('invariants')!=1:
        raise AuditError(f'{model}: invalid quadratic field counts')
    possibilities = {
       (0,0,0,2,0,0):'self_s1',
       (0,0,0,0,0,2):'self_s2',
       (0,0,0,1,0,1):'mixed',
    }
    kind = possibilities.get(tuple(counts))
    if kind is None:
        raise AuditError(f'{model}: unexpected quadratic field structure {counts}')
    dim1,dim2 = REPS[model[0]]
    alpha = int(model.split('_')[1])
    valid = ((kind=='self_s1' and alpha==0 and dim1%2==1) or
             (kind=='self_s2' and alpha==-2 and dim2%2==1) or
             (kind=='mixed' and alpha==-1 and dim1==dim2==2))
    if not valid:
        raise AuditError(f'{model}: quadratic operator inconsistent with representation')
    return kind


def analytic_spectrum(kind: str, m1: float, m2: float, mu: complex,
                      dim1: int, dim2: int) -> list[float]:
    """Gauge-unbroken spherical-basis spectrum, with NO hidden 1/2."""
    if kind not in KINDS:
        raise AuditError('Unknown quadratic kind')
    m1=finite(m1,'mS1Sq');m2=finite(m2,'mS2Sq')
    if not all(math.isfinite(x) for x in (mu.real,mu.imag)):
        raise AuditError('Nonfinite mu2')
    if kind=='none':
        if mu != 0:
            raise AuditError('No allowed quadratic operator but mu2 nonzero')
        return sorted([m1]*(2*dim1)+[m2]*(2*dim2))
    if kind=='self_s1':
        return sorted([m1-2*abs(mu)]*dim1+[m1+2*abs(mu)]*dim1+
                      [m2]*(2*dim2))
    if kind=='self_s2':
        return sorted([m1]*(2*dim1)+[m2-2*abs(mu)]*dim2+
                      [m2+2*abs(mu)]*dim2)
    mid=(m1+m2)/2
    rad=math.hypot((m1-m2)/2,abs(mu))
    return sorted([mid-rad]*(2*dim1)+[mid+rad]*(2*dim1))


def close_spectrum(a, b, label: str) -> None:
    if not isinstance(a,list) or len(a)!=len(b):
        raise AuditError(f'{label}: incorrect eigenvalue multiplicity')
    sorted_a=sorted(finite(x,label) for x in a)
    for x,y in zip(sorted_a,b):
        if not math.isclose(x,y,rel_tol=3e-10,abs_tol=2e-6):
            raise AuditError(f'{label}: mass-spectrum mismatch {x} vs {y}')


def stage_diagnostics(record: dict, model: str, kind: str, stage: str) -> dict:
    if stage not in ('initial','final'):
        raise AuditError('Invalid stage')
    base = record.get(stage+'_baseline')
    spectrum=record.get(stage+'_scalar_mass_spectrum')
    if not isinstance(base,dict) or not isinstance(spectrum,dict):
        raise AuditError(f'{model}/{stage}: missing state or spectrum')
    m1=finite(base.get('mS1Sq'),model+'/'+stage+'/mS1Sq')
    m2=finite(base.get('mS2Sq'),model+'/'+stage+'/mS2Sq')
    mu=complex_value(record.get('mu2_'+stage+'_gev2'),model+'/'+stage+'/mu2')
    if spectrum.get('quadratic_kind')!=kind or spectrum.get('model')!=model or \
            spectrum.get('analytic_spectrum_verified') is not True or \
            spectrum.get('included_in_threshold_matching') is not False or \
            spectrum.get('matchete_mass_basis_verified') is not False or \
            spectrum.get('rgbeta_mass_normalisation_verified') is not False:
        raise AuditError(f'{model}/{stage}: invalid candidate spectrum provenance')
    if abs(mu-complex_value(spectrum.get('extra_mass_coefficient_gev2'), 'spectrum mu2'))>1e-10:
        raise AuditError(f'{model}/{stage}: inconsistent mu2 and spectrum')
    dims=REPS[model[0]]
    if spectrum.get('dimensions')!=list(dims):
        raise AuditError(f'{model}/{stage}: incorrect representation dimensions')
    computed=analytic_spectrum(kind,m1,m2,mu,*dims)
    zero=analytic_spectrum(kind,m1,m2,0j,*dims)
    close_spectrum(spectrum.get('mass_squared_eigenvalues_gev2'),computed,
                   model+'/'+stage+'/numeric_spectrum')
    close_spectrum(spectrum.get('analytic_mass_squared_eigenvalues_gev2'),computed,
                   model+'/'+stage+'/analytic_spectrum')
    tach=bool(any(x<0 for x in computed))
    if spectrum.get('tachyonic') is not tach:
        raise AuditError(f'{model}/{stage}: tachyon status contradicts spectrum')
    if not math.isclose(finite(spectrum.get('minimum_mass_squared_gev2'),'spectrum min'),
                        computed[0],rel_tol=3e-10,abs_tol=2e-6):
        raise AuditError(f'{model}/{stage}: minimum eigenvalue inconsistent')
    normalise=lambda v: max(abs(v),1e-20)
    max_shift_sq=max(abs(a-b)/normalise(b) for a,b in zip(computed,zero))
    # Masses are real only for non-tachyonic states.
    max_shift_mass=(max(abs(math.sqrt(a)-math.sqrt(b))/math.sqrt(b)
                     for a,b in zip(computed,zero))
                     if all(a>0 for a in computed+zero) else None)
    min_diag=min(abs(m1),abs(m2))
    fraction=(abs(mu)/min_diag if min_diag>0 else None)
    mixing_sin2theta=None
    if kind=='mixed':
        denominator=math.hypot(m1-m2,2*abs(mu))
        mixing_sin2theta=(2*abs(mu)/denominator if denominator else None)
    return {'mu2_abs_gev2':abs(mu),
            'mu2_complex_gev2':{'real':mu.real,'imag':mu.imag},
            'mS1Sq_gev2':m1,'mS2Sq_gev2':m2,
            'mu2_over_min_abs_diagonal_m2':fraction,
            'sin_2theta_mixed_doublets':mixing_sin2theta,
            'zero_mu_mass_squared_eigenvalues_gev2':zero,
            'nonzero_mu_mass_squared_eigenvalues_gev2':computed,
            'max_relative_mass_squared_shift_at_fixed_diagonals':max_shift_sq,
            'max_relative_mass_shift_at_fixed_diagonals':max_shift_mass,
            'all_scalar_mass_squared_positive':all(v>0 for v in computed),
            'matching_neutrino_mass_shift_evaluated':False}


def audit_model(model: str, data: dict, manifest: dict, fraction: float) -> dict:
    if model not in MODELS:
        raise AuditError('Invalid ordinary T3 model')
    kind=kind_from_manifest(manifest['models'][model],model)
    if (data.get('format')!='general-t3-candidate-quadratic-eft1-v1' or
            data.get('status')!='CandidateQuadraticAndFullQuarticEFT1Integrated' or
            data.get('model')!=model or data.get('quadratic_kind')!=kind or
            data.get('mu2_rgbeta_symbol') != KINDS[kind] or
            data.get('mu2_beta_included') is not (kind!='none') or
            data.get('complete_one_loop_quartic_class_projected') is not True or
            data.get('physical_matchete_rgbeta_mapping_certified') is not False or
            data.get('scalar_mass_threshold_matching_performed') is not False or
            data.get('wilson_coefficient_matching_performed') is not False or
            data.get('physical_neutrino_mass_prediction') is not False):
        raise AuditError(f'{model}: expected uncertified candidate EFT1 report')
    if (data.get('quadratic_operators_evolved')!=int(kind!='none') or
            not math.isclose(finite(data.get('mu2_fraction_initial'),'mu2_fraction'),
                             fraction if kind!='none' else 0.0,rel_tol=1e-12,abs_tol=1e-12)):
        raise AuditError(f'{model}: mismatched quadratic input scenario')
    initial=stage_diagnostics(data,model,kind,'initial')
    final=stage_diagnostics(data,model,kind,'final')
    source_fraction=initial['mu2_over_min_abs_diagonal_m2']
    if kind!='none':
        if not math.isclose(source_fraction,fraction,rel_tol=5e-11,abs_tol=1e-12):
            raise AuditError(f'{model}: initial quadratic boundary is inconsistent')
    elif (initial['mu2_abs_gev2'] or final['mu2_abs_gev2']):
        raise AuditError(f'{model}: forbidden nonzero quadratic coefficient')
    return {
        'model':model,'quadratic_kind':kind,
        'symbol':KINDS[kind],
        'can_be_discarded_by_loop_order_alone':kind=='none',
        'direct_effect_on_scalar_propagators_at_tree_level':kind!='none',
        'potential_leading_one_loop_C5_effect_if_nonzero':kind!='none',
        'proven_nonzero_C5_contribution':False,
        'initial':initial,'final':final,
        'interpretation':(
           'No additional gauge-invariant holomorphic scalar bilinear for this model.'
           if kind=='none' else
           'Allowed tree-level scalar mass splitting; whether it modifies C5 requires physical matching.'
           if kind.startswith('self') else
           'Allowed tree-level doublet mixing; mass shifts can be small even when Yukawa basis rotation is not.'
        ),
        'omission_condition':(
          'No additional quadratic coefficient exists in the audited field basis.'
          if kind=='none' else
          'Set this parameter to zero at a declared input scale, and check radiative regeneration; or match nonzero mixing explicitly.'
        )
    }


def run(root: Path=ROOT, output_dir: Path | None=None) -> dict:
    root=root.resolve()
    output_dir=output_dir or root/DEST
    manifest_path=root/MANIFEST
    batch2_path=root/BATCH2
    folder=root/SOURCE
    summary_path=folder/'all_models_quadratic_eft1_summary.json'
    manifest=read_json(manifest_path)
    batch2=read_json(batch2_path)
    total=read_json(summary_path)
    if (manifest.get('model_count')!=16 or set(manifest.get('models',{}))!=set(MODELS)):
        raise AuditError('Expected all 16 ordinary T3 operator manifests')
    kinds={m:kind_from_manifest(manifest['models'][m],m) for m in MODELS}
    if sum(k!='none' for k in kinds.values())!=8:
        raise AuditError('Expected exactly eight additional quadratic slots')
    counts=batch2.get('counts',{})
    if (batch2.get('status')!='WickNecessaryConditionsReproduceSavedC5Inventory' or
            counts.get('models_checked')!=16 or counts.get('additional_quartics_total')!=50 or
            counts.get('new_direct_wick_eligible')!=4 or counts.get('extra_quadratics_not_screened')!=8):
        raise AuditError('Batch 2 prerequisite is missing or inconsistent')
    if (total.get('status')!='Candidate16ModelQuadraticEFT1Integrated' or
            total.get('models_checked')!=16 or total.get('new_quartics_evolved')!=38 or
            total.get('reused_quartics_evolved')!=12 or total.get('quadratic_model_slots_evolved')!=8 or
            total.get('physical_mapping_certified') is not False or
            total.get('matching_extended') is not False or
            set(total.get('models',{}))!=set(MODELS)):
        raise AuditError('Uncertified/missing 16-model candidate quadratic summary')
    fraction=finite(total.get('mu2_fraction'),'summary mu2_fraction')
    if not 0<=fraction<=0.45:
        raise AuditError('Unphysical/out-of-range input quadratic fraction')
    results=[]
    hashes={str(MANIFEST):sha256(manifest_path.read_bytes()).hexdigest(),
            str(BATCH2):sha256(batch2_path.read_bytes()).hexdigest(),
            str(SOURCE/'all_models_quadratic_eft1_summary.json'):sha256(summary_path.read_bytes()).hexdigest()}
    for model in MODELS:
        filename=model.replace('+','p').replace('-','m').replace('_','_')
        path=folder/f'{filename}_eft1_quadratic_running.json'
        record=read_json(path)
        row=audit_model(model,record,manifest,fraction)
        record_summary=total['models'][model]
        if (record_summary.get('quadratic_kind')!=row['quadratic_kind'] or
                record_summary.get('tachyonic') is not (not row['final']['all_scalar_mass_squared_positive'])):
            raise AuditError(f'{model}: saved cross-file summary does not agree')
        results.append(row)
        hashes[str(SOURCE/path.name)]=sha256(path.read_bytes()).hexdigest()
    if len(results)!=16 or sum(r['quadratic_kind']!='none' for r in results)!=8:
        raise AuditError('Quadratic report coverage mismatch')
    impact={
      'format':'t3-scalar-quadratic-propagator-impact-audit-v1',
      'status':'CandidateQuadraticMassSpectraValidatedNotC5',
      'scope':'gauge-unbroken scalar mass spectra at identical saved diagonal m^2 (not C5)',
      'models_checked':16,'additional_quadratic_slots':8,
      'self_bilinears':6,'mixed_bilinears':2,
      'initial_mu2_fraction':fraction,
      'physical_matchete_rgbeta_mapping_certified':False,
      'general_neutrino_C5_matching_completed':False,
      'zero_mu_RGE_trajectory_recomputed':False,
      'omission_by_loop_counting_for_any_nonzero_quadratic':False,
      'no_tachyons_in_reported_candidate_spectra':all(
          r[s]['all_scalar_mass_squared_positive'] for r in results for s in ('initial','final')),
      'source_sha256':hashes,
      'models':results,
      'limitations':[
         'Nonzero dimension-two mass insertions carry no additional loop factor; the eight slots cannot be generically discarded at leading one-loop order.',
         'Mass spectra were compared at identical diagonal masses; no independent mu2=0 RG evolution was run.',
         'Even a tiny mass-eigenvalue change does not bound the change in C5 (Yukawa rotations and cancellations can be linear in mu2).',
         'Spherical SU2 mass conventions are not certified against physical Matchete/RGBeta components.',
         'Electroweak-symmetry breaking, finite threshold matching, general Weinberg C5 and full vacuum stability are not calculated.',
         'Whether a zero quadratic boundary stays zero under the actual RGEs is not checked.',
      ]}
    output_dir.mkdir(parents=True,exist_ok=True)
    (output_dir/'quadratic_impact_audit.json').write_text(json.dumps(impact,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    table_path=output_dir/'quadratic_impact_audit.csv'
    columns=['model','quadratic_kind','mu2_initial_abs_gev2','mu2_final_abs_gev2',
             'mu2_final_over_min_diagonal','max_final_mass_squared_shift_fraction',
             'max_final_mass_shift_fraction','final_sin_2theta','final_spectrum_positive',
             'extra_mass_term_may_enter_one_loop_C5']
    with table_path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader()
        for r in results:
            a,b=r['initial'],r['final'];writer.writerow({
                'model':r['model'],'quadratic_kind':r['quadratic_kind'],
                'mu2_initial_abs_gev2':a['mu2_abs_gev2'],
                'mu2_final_abs_gev2':b['mu2_abs_gev2'],
                'mu2_final_over_min_diagonal':b['mu2_over_min_abs_diagonal_m2'],
                'max_final_mass_squared_shift_fraction':b['max_relative_mass_squared_shift_at_fixed_diagonals'],
                'max_final_mass_shift_fraction':b['max_relative_mass_shift_at_fixed_diagonals'],
                'final_sin_2theta':b['sin_2theta_mixed_doublets'],
                'final_spectrum_positive':b['all_scalar_mass_squared_positive'],
                'extra_mass_term_may_enter_one_loop_C5':r['potential_leading_one_loop_C5_effect_if_nonzero']})
    doc=['# T3 scalar relevance - Batch 3: quadratic mass impact','',
         '- Models: **16**; allowed additional bilinears: **8** (6 self, 2 mixed).',
         f'- Candidate initial mass coefficient: **{fraction:.6g} x min(mS1Sq,mS2Sq)** for allowed cases.',
         f'- All candidate spectra have positive mass squared: **{impact["no_tachyons_in_reported_candidate_spectra"]}**.',
         '- Independent zero-mu2 RG trajectory computed: **no**.',
         '- Matchete/RGBeta physical normalization certified: **no**.',
         '- General C5 with nonzero mixing matched: **no**.','',
         '| Model | Kind | Final |mu2| [GeV^2] | Max relative m^2 shift | sin(2theta) |',
         '|---|---|---:|---:|---:|']
    for r in results:
        f=r['final'];s=f['sin_2theta_mixed_doublets'];doc.append(
             f'| {r["model"]} | {r["quadratic_kind"]} | {f["mu2_abs_gev2"]:.6g} | '
             f'{f["max_relative_mass_squared_shift_at_fixed_diagonals"]:.3g} | '
             f'{"-" if s is None else format(s,".3g")} |')
    doc+=['','## Interpretation','',
           'A nonzero quadratic is a tree-level insertion in a scalar propagator; it can enter an otherwise one-loop neutrino-mass diagram at one-loop order. Its omission is NOT justified by perturbative loop counting.',
           'A smaller eigenvalue shift is not a neutrino-mass error estimate: in the mixed B/C models the eigenvector/Yukawa rotation is linear in mu2/(m1^2-m2^2) while eigenvalue shifts begin quadratically away from degeneracy.',
           'Setting an allowed quadratic coefficient to zero is a benchmark boundary assumption. Establish RG generation and nonzero-mixing C5 before claiming a quantitative neutrino-mass error bound.','',
           '## What was verified','',
           'The saved candidate quadratic spectra were checked against independent closed-form eigenvalues and their zero-mu2 counterparts at fixed diagonal masses. This does NOT re-run native Matchete/RGBeta or certify the physical component basis.', '']
    (output_dir/'quadratic_impact_summary.md').write_text('\n'.join(doc),encoding='utf-8')
    return impact


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    args=p.parse_args()
    result=run(args.root)
    print(f'Batch 3 candidate quadratic audit: {result["models_checked"]} models, '
          f'{result["additional_quadratic_slots"]} allowed bilinears; '
          'no physical C5 conclusion')
    print('Output:',(args.root/DEST/'quadratic_impact_summary.md').resolve())


if __name__=='__main__':
    main()
