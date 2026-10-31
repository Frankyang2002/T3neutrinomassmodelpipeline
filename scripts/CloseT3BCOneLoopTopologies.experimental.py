#!/usr/bin/env python3
"""v21: close two-trivalent-vertex branch of B/C(-1) LLHH one-loop topology audit.

Source-verified conditional claim: archived B/C renormalizable UV source with the
T3 two lepton Yukawas, unbroken SU(2)L x U(1)Y, Z2, zero extra mixed quadratic,
and no nonrenormalizable or EWSB-induced cubic insertions has only the three
one-quartic-vertex one-loop LLHH triangles at order y_i y_j lambda.

Does NOT prove full sequential renormalized C5, running, matching scheme, or
validity for nonzero independent quadratic or changed interaction inventory.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
STEMS={'B_m1':'B_-1','C_m1':'C_-1'}
AUD=Path('output/audits')
MODEL_DIR=Path('output/general_t3_mass_eigen_matchete')
SOURCE_MANIFEST=Path('output/scalar_relevance/operator_classification.json')
SOURCE_UV=Path('Lagrangian/interactions/T3Topology.wl')
SOURCE_MODEL=Path('Lagrangian/model/T3ModelCatalog.wl')
SOURCE_BUILDER=Path('Lagrangian/model/LagrangianBuilder.wl')
SOURCE_SCALAR=Path('Lagrangian/interactions/GeneralScalarOperators.wl')
V20=AUD/'t3_bc_one_loop_topology_completeness_v20/summary.json'
V19=AUD/'t3_bc_hard_soft_subtraction_v19/summary.json'

class TopologyError(ValueError): pass

def need(ok,why):
    if not ok:raise TopologyError(why)

def source(root,path):
    p=root/path
    need(p.is_file(),'Missing local evidence: '+str(p))
    return p.read_text(encoding='utf-8-sig')

def js(root,path):
    v=json.loads(source(root,path));need(isinstance(v,dict),'Expected JSON object: '+str(path));return v

def checksum(root,path):
    return hashlib.sha256((root/path).read_bytes()).hexdigest()

def has_cubic_su2_singlet(dims):
    """Exact SU(2) angular-momentum addition for three irreps, using 2j integers."""
    need(len(dims)==3 and all(type(x) is int and x>=1 for x in dims),'Bad SU2 dimensions')
    a,b,c=(x-1 for x in dims)
    return abs(a-b)<=c<=a+b and ((a+b+c)%2==0)

def hypercharge_singlet(charges_twice):
    need(len(charges_twice)==3 and all(type(x) is int for x in charges_twice),'Bad Y2 charges')
    return sum(charges_twice)==0

def theoretical_graphs(d_f):
    """Enumerate all nonnegative cubic/ quartic renormalizable vertex budgets.

    2 original BSM Yukawa vertices have three legs each; exactly 4
    external LLHH legs and one loop give sum(v-2)=4. The Yukawas consume 2.
    With no other higher-valence UV vertices, remaining 2 is [2] or [1,1].
    """
    need(d_f in (1,3),'Only B(-1) singlet/C(-1) triplet F')
    # In units of twice hypercharge, S1 and Bar[S2] both carry -1;
    # NewScalarX,NewScalarY both transform as SU(2) doublets, Y=-1/2.
    needed_higgs=(-1,-1)  # Hdag + Hdag, the PL holomorphic orientation
    counts={'F_dimension':d_f,'S1_dimension':2,'S2_dimension':2,
            'hypercharge_twice':{'H':1,'Hdag':-1,'S1':-1,
                                 'BarS2':-1,'F':0},
            'graph_valence_budget':{'external_legs':4,'loops':1,
                                    'sum_valence_minus_two':4,
                                    'two_T3_Yukawa_vertices_contribution':2,
                                    'remaining':2,
                                    'partitions':[[4],[3,3]]}}
    # HSS and HFF candidates: list all conjugate charge orientations of
    # two identical-representation scalar/fermion endpoints, then check
    # both SU2 tensor parity and twice-hypercharge conservation.
    forbidden=[]
    for kind,d,y2 in (('H_S_S',2,-1),('H_F_F',d_f,0)):
        for nbar in (0,1,2):
            scalar_charges=(y2 if nbar<2 else -y2,
                            y2 if nbar==0 else -y2)
            for H in (-1,1):
                charges=(H,*scalar_charges)
                su2=has_cubic_su2_singlet((2,d,d))
                hyp=hypercharge_singlet(charges)
                need(not(su2 and hyp),'Unexpected renormalizable H BSM BSM trilinear')
                forbidden.append({'class':kind,'external_H_twice_Y':H,
                   'internal_pair_twice_Y':list(scalar_charges),
                   'SU2_singlet':su2,'U1Y_neutral':hyp,
                   'allowed':su2 and hyp})
    # H_SMfermion_BSMfermion is gauge-allowed for F=1 or 3, but is
    # excluded by the model's imposed Z2 interaction selection; the actual
    # archived UV interaction inventory is checked separately below.
    # Gauge H Hdag V gives one H and one Hdag; cannot provide Hdag Hdag.
    # SM Yukawa H f f contains no BSM scalar legs and cannot saturate the
    # two T3 scalar ends at fixed four-vertex budget.
    return {'graph_counting':counts['graph_valence_budget'],
       'SU2_HSS_and_HFF_tensor_search':forbidden,
       'Z2_H_SMfermion_BSMfermion_vertex_allowed':False,
       'exclusion_of_gauge_allowed_H_l_F_relies_on_model_Z2_and_archived_UV_inventory':True,
       'unbroken_phase_Higgs_odd_scalar_cubic_allowed':False,
       'gauge_Higgs_current_H_Hdag_V_compatible_with_Hdag_Hdag':False,
       'two_SM_Yukawa_Hff_vertices_can_absorb_BSM_scalar_ends':False,
       'additional_two_trivalent_graph_contributes_to_LLHH':False,
       'proof':'Two required T3 Yukawas supply the external leptons and an F line; both scalar ends have twice-hypercharge -1. Without an anomalous scalar propagator they require a vertex supplying two conjugate scalar ends. Any two-trivalent graph must attach external holomorphic Higgses through H-containing cubics, but HSS/HFF are forbidden by SU2/Y, H-SM-BSM fermion is excluded by the defining Z2 and checked absent in archived UV interactions (not forbidden by gauge symmetry), gauge Hdag-H current has mixed Higgs orientation and SM Yukawa cubics contain no BSM scalar ends.',
       'scope':'unbroken renormalizable T3 source with two original Yukawas, no extra quadratic or new-fermion Yukawas'}

def require_source_contract(root, v20):
    prior=v20.get('source_contract_sha256')
    need(isinstance(prior,dict) and len(prior)==5,'v20 UV source contract missing')
    for path,expect in prior.items():
        p=Path(path)
        need(not p.is_absolute() and '..' not in p.parts,'Invalid source path')
        need(checksum(root,p)==expect,'UV source changed since v20: '+path)
    topology=source(root,SOURCE_UV)
    code=source(root,SOURCE_MODEL)
    builder=source(root,SOURCE_BUILDER)
    scalar=source(root,SOURCE_SCALAR)
    for token in ('useCConj = which === 1;','scalarBar = which === 2;',
                  'T3GeneralExtraYukawaCount[model_Association]',
                  'model["Scalar1", "Y"] === 0',
                  'model["Scalar2", "Y"] === 0',
                  'T3ExperimentalGeneralFermionsQ[]',
                  'lambdaT3[] H[i] H[j] s1 s2'):
        need(token in topology,'Original Yukawa or extra-Yukawa source changed: '+token)
    for token in ('"B" -> {2, 2, 1}', '"C" -> {2, 2, 3}',
                  'alpha/2', '(alpha + 2)/2','(alpha + 1)/2'):
        need(token in code,'B/C SU2/hypercharge contract changed: '+token)
    need('LInt = Total[Lookup[assembled, "Expression", {}]]' in builder,
         'UV interaction assembly changed')
    need('T3ExtendedBuildExtraQuarticCandidates' in scalar and
         'T3ExtendedBuildCandidates' in scalar,
         'Extended scalar source changed')
    return {'source_files_sha256':dict(prior),
            'extra_conjugate_scalar_Yukawas_for_B_C_alpha_m1':0,
            'T3_Yukawa_source_fields':['S1','BarS2'],
            'SU2_dimensions':{'H':2,'S1':2,'S2':2,'B_F':1,'C_F':3},
            'twice_hypercharges':{'H':1,'S1':-1,'BarS2':-1,'F':0}}

def validate_manifest(root):
    j=js(root,SOURCE_MANIFEST)
    need(j.get('summary',{}).get('strict_complete') is True and
         j['summary'].get('models_checked')==16,'Original operator inventory not strictly complete')
    rows={}
    for model in STEMS.values():
        a=[x for x in j['operators'] if x.get('model')==model]
        triples={x['coupling']:x for x in a if x['kind'] in ('baseline_topology','additional_quartic')}
        need(len(triples)==8,model+': quartic inventory changed')
        expected={'lambdaT3':(2,0,1,0,0,1),
          'lambdaT3ExtraQ020200I1':(0,2,0,2,0,0),
          'lambdaT3ExtraQ020020I1':(0,2,0,0,2,0)}
        for q,spec in expected.items():
            need(q in triples and tuple(triples[q]['counts_H_Hdag_S1_S1dag_S2_S2dag'])==spec,
                 model+': triangle vertex signature changed: '+q)
        mass=[x for x in a if x.get('kind')=='additional_quadratic']
        need(len(mass)==1 and mass[0]['coupling']=='mu2_extra_000101_I1' and
             tuple(mass[0]['counts_H_Hdag_S1_S1dag_S2_S2dag'])==(0,0,0,1,0,1),
             model+': missing/new independent mixed quadratic')
        rows[model]={'required_HHSS_vertices':list(expected),
           'extra_mixed_quadratic':mass[0]['coupling'],
           'extra_quartic_count':7,'original_baseline_topology_count':1}
    return rows,checksum(root,SOURCE_MANIFEST)

def validate_v20(v20):
    need(v20.get('format')=='t3-bc-one-loop-llhh-topology-completeness-v20' and
         v20.get('status')=='RestrictedY2LambdaOneLoopGraphInventoryAudited' and
         v20.get('all_one_loop_SM_BSM_diagrams_exhaustively_generated') is False and
         v20.get('physical_C5_certified') is False,'v20 gate changed')
    result={}
    for stem,model in STEMS.items():
        row=v20.get('models',{}).get(stem)
        need(isinstance(row,dict) and row.get('model')==model and
             row.get('direct_quartic_wick_sector_exhaustive') is True and
             row.get('restricted_one_loop_hard_soft_identity_v19') is True,
             model+': incomplete v20 source gate')
        stages=row.get('F_Y_X_matching_sector_scan',{})
        for s,n in (('F',204),('Y',132),('X',36)):
            need(stages.get(s+'_full',{}).get('selected_llhh_GammaCC_terms')==n,
                 model+': v20 matched-EFT sector census changed: '+s)
            need(stages.get(s+'_raw',{}).get('selected_llhh_GammaCC_terms')==n,
                 model+': v20 raw sector census changed: '+s)
        graphs=row.get('one_loop_triangle_graphs',[])
        need(len(graphs)==3 and {g.get('channel') for g in graphs}=={'y1_y1','y1_y2','y2_y2'},
             model+': incomplete Wick channels')
        result[stem]={'full_matching_EFT_sector_counts':{'F':204,'Y':132,'X':36},
                      'triangle_channels':[g['channel'] for g in graphs]}
    return result

def verify_matched_source_hashes(root,v20):
    """Do not accept a stale v20 census after a local source change."""
    verified={}
    for stem,model in STEMS.items():
        row=v20['models'][stem]
        saved=row.get('source_matching_sha256',{})
        need(set(saved)=={'F_raw','F_full','Y_raw','Y_full','X_raw','X_full'},
             model+': v20 source hash contract incomplete')
        base_f=f'{stem}_final_eft1_fermion_first_oneloop'
        base_y=f'{stem}_final_y_only_canonical_match'
        base_x=f'{stem}_final_x_from_saved_y_match_x_followup'
        d_f=MODEL_DIR/'eft1_fermion_first_oneloop'
        d_y=MODEL_DIR/'eft1_scalar_Y_only_canonical_probe'
        d_x=MODEL_DIR/'eft1_scalar_xy_fresh_kernel_v3'
        paths={
          'F_raw':d_f/f'{base_f}_raw_hard_llhh_sector.inputform.txt',
          'F_full':d_f/f'{base_f}_eft1_one_loop_correction.inputform.txt',
          'Y_raw':d_y/f'{base_y}_raw_scalar_llhh.inputform.txt',
          'Y_full':d_y/f'{base_y}_scalar_one_loop_correction.inputform.txt',
          'X_raw':d_x/f'{base_x}_x_raw_llhh.inputform.txt',
          'X_full':d_x/f'{base_x}_x_loop_difference.inputform.txt',
        }
        for name,path in paths.items():
            need((root/path).is_file(),'Missing current matched source: '+str(path))
            need(checksum(root,path)==saved[name],
                 model+': prior v20 match census stale for '+name)
        verified[stem]=dict(saved)
    return verified

def run(root,output):
    v20=js(root,V20);prev=validate_v20(v20)
    source_hashes=verify_matched_source_hashes(root,v20)
    proven=require_source_contract(root,v20)
    rows,man_hash=validate_manifest(root)
    wick=js(root,Path('output/scalar_relevance/leading_c5_vertex_audit.json'))
    need(wick.get('status')=='WickNecessaryConditionsReproduceSavedC5Inventory' and
         wick.get('batch1_report_sha256')==man_hash and
         wick.get('source_provenance_verified_on_local_checkout') is True,
         'Prior source-verified complete Wick inventory no longer matches')
    zero=js(root,Path('output/scalar_relevance/final_leading_scalar_potential.json'))
    need(zero.get('conditions',{}).get('additional_quadratics_zero_at_uv_boundary') is True,
         'Source zero-extra-quadratic benchmark assumption changed')
    v19=js(root,V19)
    need(v19.get('format')=='t3-bc-one-loop-hard-soft-msbar-restricted-v19' and
         v19.get('full_T3_renormalized_matching_certified') is False,
         'v19 limited hard-soft subtraction gate changed')
    results={}
    for stem,model in STEMS.items():
        archive=js(root,MODEL_DIR/f'{stem}_final_matchete.json')
        need(archive.get('model')==model and archive.get('mass_basis_ordinary_check_passed') is True and
             archive.get('heavy_mass_basis_detailed_check_passed') is True and
             archive.get('all_17_quartics_present') is True and
             archive.get('both_yukawa_interactions_present') is True and
             archive.get('source_extra_quartics')==7,
             model+': archived UV source incomplete/changed')
        interactions=archive.get('source_interactions')
        need(isinstance(interactions,list) and len(interactions)==18 and
             archive.get('source_new_interaction_count')==18,
             model+': unexpected UV source interaction count')
        names=[x['name'] for x in interactions]
        cls=[x['class'] for x in interactions]
        need(len(names)==len(set(names)) and names.count('Yukawa1_GenericCG')==1 and
             names.count('Yukawa2_GenericCG')==1 and
             cls.count('AdditionalQuartic')==7 and cls.count('Baseline')==11 and
             not any('Yukawa' in n and n not in ('Yukawa1_GenericCG','Yukawa2_GenericCG') for n in names) and
             not any('Quadratic' in n or 'ExtraM' in n for n in names),
             model+': additional Yukawa or mixed quadratic present in archived UV source')
        extra=archive.get('native_additional_quartic_tensors',{}).get('Certificates',{})
        need(isinstance(extra,dict) and len(extra)==7 and
             all(v.get('ExactTensorProjectionVerified') is True for v in extra.values()),
             model+': incomplete rotated quartic SU2 tensors')
        d_f=1 if model=='B_-1' else 3
        theory=theoretical_graphs(d_f)
        results[stem]={'model':model,'UV_archive_sha256':checksum(root,MODEL_DIR/f'{stem}_final_matchete.json'),
           'independent_quadratic_coupling_in_archived_UV_interactions':False,
           'extra_quadratic_permitted_by_symmetry_but_absent_from_archived_UV':True,
           'additional_fermion_Yukawa_in_archived_UV':False,
           'one_loop_two_trivalent_graphs_excluded_for_archived_UV':True,
           'one_loop_one_quartic_channels_exhaustive_for_archived_UV':True,
           'renormalizable_one_loop_topology_exhaustive_at_order_y2_lambda_within_source':True,
           'one_loop_topology_exhaustive_for_general_nonzero_extra_mass':False,
           'prior_matched_EFT_check':prev[stem],
           'scalar_operator_inventory':rows[model],
           'trivalent_graph_proof':theory,
           'physical_C5_certified':False}
        print(f'{model}: trivalent graph class excluded by SU2/U1/Z2 and field-flow; archived extra quadratic absent')
    record={'format':'t3-bc-one-loop-renormalizable-topology-closure-v21',
      'status':'ConditionalArchivedUVLeadingLLHHGraphClassesClosed',
      'source_verified':True,
      'UV_operator_manifest_sha256':man_hash,
      'source_matched_F_Y_X_sha256_from_v20_verified_unchanged':source_hashes,
      'original_UV_source':proven,
      'models':results,
      'scope':'unbroken-gauge-phase, renormalizable B(-1)/C(-1) archived 2-Yukawa, 7-extra-quartic UV interaction set with no independent extra mixed quadratic; leading O(y^2 lambda)/(16 pi^2 MF) Weinberg topology only',
      'one_loop_diagram_topology_complete_for_archived_UV_under_scope':True,
      'not_complete_for_general_T3_Lagrangian_with_nonzero_extra_quadratic':True,
      'field_EOM_finite_scheme_and_running_certified':False,
      'full_sequential_physical_matching_certified':False,
      'physical_C5':None,
      'remaining_work':['Independent finite UV->EFT1->EFT2 matching prescription/equivalence for complete operators and scheme',
                        'Consistent RG running and boundary values of Yukawas/quartics',
                        'Nonzero independent mixed quadratic if generalizing beyond archived source',
                        'Potential additional higher-order terms or EWSB insertions outside leading unbroken dimension-five scope']}
    target=output if output.is_absolute() else root/output
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: both renormalizable one-loop topology classes closed for archived B/C UV source')
    print('FULL PHYSICAL C5: NOT CERTIFIED')
    print('Saved: '+str(target))

def self_test():
    need(has_cubic_su2_singlet((2,2,2)) is False,'Odd-half-integer triple failed')
    need(has_cubic_su2_singlet((2,1,1)) is False,'H singlet-fermion pair permitted')
    need(has_cubic_su2_singlet((2,3,3)) is False,'H triplet-fermion pair permitted')
    need(has_cubic_su2_singlet((2,2,1)) is True,'Known 2x2x1 invariant was rejected')
    need(has_cubic_su2_singlet((3,2,2)) is True,'Known 3x2x2 invariant was rejected')
    need(hypercharge_singlet((1,-1,0)) is True and
         hypercharge_singlet((-1,-1,-1)) is False,'Y2 charge audit broken')
    for f in (1,3):
        result=theoretical_graphs(f)
        need(len(result['SU2_HSS_and_HFF_tensor_search'])==12 and
             not any(x['allowed'] for x in result['SU2_HSS_and_HFF_tensor_search']),
             'Trivalent group-theory scan failed')
        need(result['additional_two_trivalent_graph_contributes_to_LLHH'] is False,
             'Unknown trivalent branch')
    # Check fail-closed on an incorrect or incomplete prior match census.
    gate={'format':'t3-bc-one-loop-llhh-topology-completeness-v20',
          'status':'RestrictedY2LambdaOneLoopGraphInventoryAudited',
          'all_one_loop_SM_BSM_diagrams_exhaustively_generated':False,
          'physical_C5_certified':False,'models':{}}
    for stem,model in STEMS.items():
        gate['models'][stem]={'model':model,'direct_quartic_wick_sector_exhaustive':True,
            'restricted_one_loop_hard_soft_identity_v19':True,
            'F_Y_X_matching_sector_scan':{f'{s}_{typ}':{'selected_llhh_GammaCC_terms':n}
                for s,n in (('F',204),('Y',132),('X',36)) for typ in ('full','raw')},
            'one_loop_triangle_graphs':[{'channel':c} for c in ('y1_y1','y1_y2','y2_y2')]}
    validate_v20(gate)
    gate['models']['B_m1']['F_Y_X_matching_sector_scan']['X_raw']['selected_llhh_GammaCC_terms']=35
    try:validate_v20(gate)
    except TopologyError:pass
    else:raise TopologyError('Incomplete v20 matched source was accepted')
    print('PASS: exact SU2 triple products, U1Y, gauge H/Hdag and Z2 two-trivalent exclusions')
    print('PASS: both F-dimensions and v20 provenance/tampering gates')
    print('PASS: general mixed quadratic explicitly remains a scope limitation')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--output',type=Path,default=AUD/'t3_bc_topology_closure_v21/summary.json')
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--require-physical',action='store_true')
    args=ap.parse_args()
    try:
        need(not args.require_physical,'Topology closure does not certify fully physical C5')
        if args.self_test:self_test()
        else:run(args.root.resolve(),args.output)
        return 0
    except (TopologyError,ValueError,OSError,TypeError,KeyError,AttributeError,IndexError) as exc:
        print('FAIL (physical C5 not certified): '+str(exc),file=sys.stderr)
        return 2

if __name__=='__main__':raise SystemExit(main())
