"""Batch 6: consolidate five saved audits into an explicitly restricted T3 potential.

Read-only over existing scientific reports; does not run beta functions, C5
matching, or alter the production Lagrangian. Requires all five completed
Batch 1--5 outputs from the same local repository.

    python -m scripts.FinalizeT3LeadingScalarPotential
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('output/scalar_relevance')
INPUTS = {
    'batch1': OUT/'operator_classification.json',
    'batch2': OUT/'leading_c5_vertex_audit.json',
    'batch3': OUT/'quadratic_impact_audit.json',
    'batch4': OUT/'quadratic_rg_zero_boundary_audit.json',
    'batch5': OUT/'reduced_scalar_rg_order_audit.json',
}
MODELS = ('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1',
          'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')
DIRECT = frozenset({'B_-1','C_-1'})
EXPECTED_2H = {'A_-4':0,'A_-2':2,'A_+0':2,'B_-3':1,'B_-1':4,'B_+1':1,
               'C_-3':1,'C_-1':4,'C_+1':1,'D_-2':2,'D_+0':2,'D_+2':0,
               'E_-4':0,'E_-2':2,'E_+0':2,'E_+2':0}

class FinalAuditError(ValueError):
    pass


def require(ok: bool, why: str) -> None:
    if not ok:
        raise FinalAuditError('Batch 6: ' + why)


def load_reports(root: Path) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    reports: dict[str, dict[str, Any]] = {}
    hashes: dict[str, str] = {}
    for name, relative in INPUTS.items():
        path = root/relative
        require(path.is_file(), 'missing report '+str(relative))
        blob = path.read_bytes()
        try:
            obj = json.loads(blob.decode('utf-8-sig'))
        except (ValueError, UnicodeError) as e:
            raise FinalAuditError('Batch 6: invalid JSON '+str(relative)) from e
        require(isinstance(obj,dict), 'report must be a JSON object: '+name)
        reports[name] = obj
        hashes[str(relative).replace('\\','/')] = sha256(blob).hexdigest()
    return reports, hashes


def validate_headlines(d: Mapping[str, dict]) -> None:
    a,b,c,e,f=(d[k] for k in ('batch1','batch2','batch3','batch4','batch5'))
    summ=a.get('summary',{})
    require(summ.get('strict_complete') is True and summ.get('models_checked')==16
            and summ.get('extra_quartic_slots')==50 and summ.get('extra_quadratic_slots')==8
            and summ.get('observed_additional_quartic_C5_slots')==4,
            'Batch 1 strict model/operator inventory inconsistent')
    k=b.get('counts',{})
    require(b.get('status')=='WickNecessaryConditionsReproduceSavedC5Inventory'
            and b.get('source_provenance_verified_on_local_checkout') is True
            and b.get('new_matching_calculation_performed') is False
            and all(k.get(n)==v for n,v in {
                'models_checked':16,'additional_quartics_total':50,
                'new_two_H_vertices':24,'new_direct_wick_eligible':4,
                'saved_C5_observed_new':4,
                'excluded_H_Hdag_by_higgs_orientation':10,
                'excluded_Hdag2_by_scalar_orientation':10,
                'new_no_H_vertex':26,'extra_quadratics_not_screened':8}.items()),
            'Batch 2 Wick-screen evidence or provenance inconsistent')
    require(c.get('status')=='CandidateQuadraticMassSpectraValidatedNotC5'
            and c.get('models_checked')==16 and c.get('additional_quadratic_slots')==8
            and c.get('self_bilinears')==6 and c.get('mixed_bilinears')==2
            and c.get('no_tachyons_in_reported_candidate_spectra') is True
            and c.get('general_neutrino_C5_matching_completed') is False
            and c.get('physical_matchete_rgbeta_mapping_certified') is False,
            'Batch 3 candidate quadratic interpretation inconsistent')
    require(e.get('model_count')==16 and e.get('quadratic_slots')==8
            and e.get('uv_and_eft1_beta_exports_checked') is True
            and e.get('source_vanishes_if_all_additional_quartics_zero_at_one_loop') is True
            and e.get('not_a_zero_mu2_trajectory') is True
            and e.get('not_a_neutrino_mass_error_estimate') is True,
            'Batch 4 zero-boundary RG evidence inconsistent')
    require(f.get('status')=='CandidateReducedScalarOneLoopRGGenerationOrderAudit'
            and f.get('models_checked')==16 and f.get('quartic_slots_checked')==50
            and f.get('quadratic_slots_checked')==8
            and f.get('uv_immediately_sourced_quartic_slots')==6
            and f.get('uv_quadratic_first_generated_at_order_2_count')==0,
            'Batch 5 RG source headline inconsistent')


def build(reports: Mapping[str, dict], hashes: Mapping[str,str]) -> dict:
    validate_headlines(reports)
    b1,b2,b3,b4,b5 = (reports[x] for x in ('batch1','batch2','batch3','batch4','batch5'))
    p2=b2.get('per_model',{})
    p3={r.get('model'):r for r in b3.get('models',[])}
    p4=b4.get('models',{})
    p5=b5.get('model_results',{})
    require(all(set(z)==set(MODELS) for z in (p2,p3,p4,p5)),
            'reports do not describe the same 16 models')
    ops=b1.get('operators')
    require(isinstance(ops,list) and len(ops)==186,
            'Batch 1 raw operator inventory missing (expected 186 rows)')
    tally=Counter(x.get('kind') for x in ops)
    require(tally=={'baseline_topology':16,'baseline_quadratic':32,
                    'baseline_quartic_family':80,'additional_quartic':50,
                    'additional_quadratic':8}, 'source operator kind counts mismatch')
    # The model-by-model output is intentionally about operator RELEVANCE,
    # not a claim to calculate a general-parameter neutrino mass.
    model_rows=[]
    for model in MODELS:
        quartics=[o for o in ops if o['model']==model and o['kind']=='additional_quartic']
        masses=[o for o in ops if o['model']==model and o['kind']=='additional_quadratic']
        observed=[o for o in quartics if o.get('observed_in_saved_C5') is True]
        per=p2[model]
        q3=p3[model]
        q4=p4[model]
        g=p5[model]
        uv=g['stages']['uv']
        eft=g['stages']['eft1']
        require(per.get('new_two_H_vertices')==EXPECTED_2H[model]
                and per.get('wick_eligible')==(2 if model in DIRECT else 0)
                and per.get('observed_in_saved_C5')==len(observed),
                model+': Wick and saved C5 model counts conflict')
        require(len(quartics)==g.get('extra_quartics')
                and len(masses)==int(g.get('extra_quadratic')),
                model+': quartic/mass registry conflict')
        require(q3.get('quadratic_kind') in ('none','self_s1','self_s2','mixed')
                and (q3['quadratic_kind']=='none')==(len(masses)==0),
                model+': quadratic kind / allowed bilinear conflict')
        require(q4.get('extra_quadratic_exists')==(len(masses)==1),
                model+': quadratic RG inventory conflict')
        direct=uv.get('direct_C5_quartics',[])
        direct_orders=uv.get('direct_C5_quartic_orders',{})
        require(len(direct)==(2 if model in DIRECT else 0)
                and set(direct)==set(direct_orders)
                and all(direct_orders[q]==1 for q in direct),
                model+': direct quartic generation does not match model C5 evidence')
        require(len(uv.get('immediately_generated_quartics',[]))==(3 if model in DIRECT else 0)
                and len(eft.get('immediately_generated_quartics',[]))==0,
                model+': source quartic generations differ from Batch 5')
        require((len(observed)==2)==(model in DIRECT),
                model+': baseline matching of B/C leading quartics inconsistent')
        require(uv.get('quadratic_first_generation_order') is None
                and eft.get('quadratic_first_generation_order') is None,
                model+': a quadratic can be generated in the candidate one-loop chain')
        require(q3.get('proven_nonzero_C5_contribution') is False,
                model+': candidate spectrum falsely claims C5 verification')
        if model in DIRECT:
            require(all(q['higgs_legs']==2 for q in observed),
                    model+': direct observed quartic has no two Higgs legs')
        model_rows.append({
            'model':model,
            'extra_quartic_slots':len(quartics),
            'extra_two_H_quartics':per['new_two_H_vertices'],
            'extra_direct_C5_quartics_if_nonzero':len(observed),
            'extra_direct_C5_couplings_Matchete':[x['coupling'] for x in observed],
            'extra_quartics_not_direct_at_selected_order':len(quartics)-len(observed),
            'extra_bilinear_kind':q3['quadratic_kind'],
            'extra_bilinear_exists':bool(masses),
            'uv_quartics_sourced_at_order_one':len(uv['immediately_generated_quartics']),
            'uv_direct_C5_quartic_rg_order':{q:direct_orders[q] for q in sorted(direct)},
            'candidate_uv_bilinear_rg_order':None,
            'candidate_eft1_bilinear_rg_order':None,
            'conclusion':(
                'Additional quartics or bilinear may affect leading C5 if nonzero at the UV boundary; zero boundary essential.'
                if model in DIRECT else
                'No newly enumerated scalar operators in this catalogue.'
                if not quartics and not masses else
                'Additional quadratic may affect leading C5 if nonzero; zero boundary essential.'
                if masses else
                'Additional quartics absent from the specified direct one-vertex C5 graph; retain for scalar-sector/RG questions.'
            ),
        })
    require(sum(r['extra_quartic_slots'] for r in model_rows)==50
            and sum(r['extra_direct_C5_quartics_if_nonzero'] for r in model_rows)==4
            and sum(r['uv_quartics_sourced_at_order_one'] for r in model_rows)==6
            and sum(r['extra_bilinear_exists'] for r in model_rows)==8,
            'final cross-model totals do not match the physics audits')
    return {
        'status':'ConditionalLeadingOneLoopT3ScalarSectorAudit',
        'scope':'16 ordinary T3 models; gauge-unbroken leading dimension-five, fixed-order one-loop C5, diagonal BSM scalar propagators at tree level',
        'source_sha256':dict(hashes),
        'input_audits':list(INPUTS),
        'conditions':{
            'additional_quartics_zero_at_uv_boundary':True,
            'additional_quadratics_zero_at_uv_boundary':True,
            'zero_extra_couplings_are_boundary_assumptions_not_symmetry_prohibitions':True,
            'additional_Higgs_vev_insertions_and_higher_dimension_operators_excluded':True,
            'full_physical_matchete_rgbeta_mapping_certified':False,
            'general_two_threshold_matching_certified':False,
            'general_C5_neutrino_mass_error_computed':False,
        },
        'quartic_result':{
            'extra_quartic_slots':50,
            'not_direct_at_selected_order':46,
            'direct_extra_if_tree_nonzero':4,
            'immediately_sourced_in_candidate_UV_one_loop_betas':6,
            'direct_extra_first_generated_in_candidate_UV_running':4,
            'candidate_EFT1_immediately_sourced_at_zero_extra_boundary':0,
        },
        'quadratic_result':{
            'allowed_complex_bilinear_slots':8,
            'unsuppressed_if_tree_nonzero':8,
            'first_generated_order_in_separate_candidate_one_loop_UV_and_EFT1_beta_chains':None,
            'no_finite_iterated_one_loop_RGE_source_found':True,
        },
        'model_results':model_rows,
        'mathematical_conclusion':(
            'Under an all-additional-couplings-zero UV boundary, the direct leading one-loop '
            'C5 is supplied by the ordinary T3 quartic and original Yukawas in the '
            'specified unbroken, diagonal-propagator setup. The four additional '
            'B/C direct quartics can be radiatively generated at one RGBeta loop, '
            'so inserting them into a one-loop amplitude is formally of second '
            'loop order. No extra bilinear is reached through the saved iterated '
            'one-loop candidate RG chains; this is not a higher-loop or '
            'threshold-matching non-generation theorem.'
        ),
        'limitations':[
            'Do not claim a fully general symmetry-allowed potential: additional terms exist and have been set zero at a specified UV boundary.',
            'Nonzero additional B/C direct quartics contribute at leading one-loop order, not at two loops.',
            'Nonzero quadratic bilinears can change tree-level scalar propagators, mixing and leading one-loop neutrino-mass amplitudes.',
            'The 46 quartics not direct in the selected topology can affect RG evolution, scalar spectra, vacuum stability and other amplitudes.',
            'The candidate UV and EFT1 RG generation analyses are separate; no newly integrated zero-boundary RG trajectory, finite threshold matching or general physical C5 correction was evaluated.',
            'Finite threshold generation, genuine two-loop beta functions, EWSB v^2/M^2 corrections, and large logs are outside the demonstrated scope.',
            'The Matchete/RGBeta physical scalar field and coupling conventions remain uncertified.',
        ],
    }


def render_markdown(report: dict) -> str:
    r=report['quartic_result']
    lines=[
        '# T3 leading-order scalar potential: conditional justification (Batch 6)', '',
        '## Claim and assumptions', '',
        'We work in the electroweak-unbroken theory, in the leading dimension-five, fixed-order one-loop matching approximation, with diagonal tree-level scalar propagators. All symmetry-allowed **additional** quartic and quadratic coefficients are set to zero at the specified UV boundary scale. This is a boundary condition, **not** a new symmetry.', '',
        'The ordinary scalar potential and original T3 topology quartic remain present. Ordinary self/portal quartics are not discarded from the UV scalar potential simply because they do not enter the selected single-vertex Weinberg graph.', '',
        'Schematically, with invariant contractions understood:', '',
        '`V_reduced = V_SM + m1^2 S1dag S1 + m2^2 S2dag S2 + V_baseline_hermitian_quartics + [lambdaT3 (HH S1 S2dag)_1 + h.c.]`', '',
        'The full symmetry-allowed potential also permits the catalogued `V_additional_quartics + V_additional_bilinears`; their UV coefficients are set zero *only in this benchmark approximation*.', '',
        '## Verified operator classification', '',
        f"- {r['extra_quartic_slots']} additional quartic model-slots: {r['not_direct_at_selected_order']} do not enter the selected direct one-quartic-vertex matching graph, and {r['direct_extra_if_tree_nonzero']} can enter (two in each of B(-1), C(-1)).",
        f"- {r['immediately_sourced_in_candidate_UV_one_loop_betas']} UV quartic slots have nonzero one-loop source when all additional scalar couplings are zero: 3 per B(-1)/C(-1). Four are direct-C5 slots; the other two are four-BSM-scalar slots.",
        '- No additional quadratic mass coefficient is reached in the finite iterated one-loop beta-source chains of either separately audited candidate UV or EFT1 system. This does not establish non-generation under genuine two-loop RGEs or threshold matching.',
        '- All eight allowed bilinears could affect leading one-loop masses if specified nonzero as independent tree-level input.', '',
        '## Model-by-model classifications', '',
        '| Model | Extra quartics | Direct extra C5 if nonzero | UV order-1 quartics | Extra bilinear |',
        '|---|---:|---:|---:|---|',
    ]
    for row in report['model_results']:
        lines.append(f"| {row['model']} | {row['extra_quartic_slots']} | {row['extra_direct_C5_quartics_if_nonzero']} | {row['uv_quartics_sourced_at_order_one']} | {row['extra_bilinear_kind']} |")
    lines += ['', '## Perturbative accounting', '',
        'With `epsilon = 1/(16 pi^2)` and `L = log(muF/muUV)`, a direct extra quartic satisfying a zero UV boundary but having a nonzero one-loop beta source obeys `delta lambda_extra ~ epsilon L beta_extra^(1)`. Its insertion into the standard one-loop C5 graph is formally `delta C5 ~ epsilon^2 L`, multiplied by Yukawa, tensor and heavy-mass factors. This is a **formal order**, not a numerical neutrino-mass error bound.', '',
        'The candidate one-loop UV beta functions generate the two direct extra quartics in each B(-1)/C(-1), and also one four-BSM-scalar quartic per model. The two extra quartics that directly source the B/C mixed bilinear are **not reached** along this zero-extra iterated one-loop source chain. The bilinear is therefore not regenerated in this specific one-loop candidate RG truncation, but must not be assumed absent in threshold effects or the full quantum theory.', '',
        'For a 100:1 RG scale interval, `|L|/(16 pi^2) = 0.02916`; the coefficient and possible cancellations determine the actual correction. Two successive such intervals are not equivalent to a certified continuous physical RG+matching computation.', '',
        '## What can be claimed', '',
        '**Supported:** A restricted leading-dimension-five, fixed-order one-loop scalar-sector benchmark with all additional symmetry-allowed couplings set to zero at the UV boundary. The 46 additional quartic slots excluded by the direct graph do not give independent single-vertex contributions in that setup; the four other slots would be leading if nonzero, but are radiatively generated by one-loop RG flow and thus enter a strict one-loop fixed-order C5 only beyond that order.', '',
        '**Not supported:** That the restricted scalar potential is the full gauge/Z2-allowed potential; that no omitted interactions are produced at higher order or thresholds; that general mixed-mass one-loop C5 is computed; that the physical neutrino-mass error is below any chosen percentage.', '',
        '## Remaining checks before numerical precision claims', '',
        '1. Choose and record the exact UV boundary scale and zero boundary conditions used by the benchmark pipeline.',
        '2. Verify the Matchete/RGBeta coupling conventions for any RG-improved numerical comparison; the present RG audit is in an uncertified candidate basis.',
        '3. If quantitative precision is desired, integrate a zero-extra UV-to-EFT1 running trajectory and compare matched C5 with the restricted result (including threshold effects); no such result is claimed here.',
        '4. Check scalar vacuum viability for the chosen retained self/portal quartics; positive mass-squared eigenvalues alone are insufficient.', '',
        '## Provenance and limitations', '',
        'This report synthesizes five saved local audits. No new symbolic matching, beta generation, RG integration, mass-basis threshold matching, or neutrino-mass calculation was performed.', '',
    ]
    lines += ['- '+s for s in report['limitations']]
    lines+=['','Inputs are SHA-256 fingerprinted in `final_leading_scalar_potential.json`.','']
    return '\n'.join(lines)


def write(report: dict, destination: Path) -> tuple[Path,Path]:
    destination.mkdir(parents=True,exist_ok=True)
    j=destination/'final_leading_scalar_potential.json'
    md=destination/'final_leading_scalar_potential.md'
    j.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    md.write_text(render_markdown(report),encoding='utf-8')
    return j,md


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output-dir',type=Path,default=OUT)
    args=parser.parse_args()
    root=args.root.resolve()
    dest=args.output_dir if args.output_dir.is_absolute() else root/args.output_dir
    try:
        reports, hashes=load_reports(root)
        final=build(reports, hashes)
        files=write(final,dest)
    except (FinalAuditError,KeyError,TypeError,ValueError,IndexError) as exc:
        parser.exit(1, f'Batch 6 FAILED: {exc}\n')
    print('Batch 6 reduced potential synthesis: 16 models, 50 quartics, 8 bilinears; no general physical C5 claim')
    print('Input audit files:',len(hashes),'SHA-256 fingerprinted')
    for f in files:print('Output:',f)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
