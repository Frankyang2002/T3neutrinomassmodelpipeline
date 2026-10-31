#!/usr/bin/env python3
"""B/C(-1) one-loop LLHH topology-completeness audit, source-anchored.

Proves exhaustiveness ONLY in the defined y_i y_j lambda single-quartic
triangle sector with diagonal free scalar propagation and no extra quadratic
insertions. Reads the complete saved F/Y/X one-loop EFT differences and checks
that their LLHH GammaCC sectors are exhausted by the three predicted quartics.

Does NOT certify an exhaustive SM+BSM amplitude or a physical C5.
No Mathematica evaluation, GitHub writes, source replacement, or RG evolution.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
STEMS = {'B_m1': 'B_-1', 'C_m1': 'C_-1'}
EXPECTED = {
    'y1_y1': 'lambdaT3ExtraQ020200I1',
    'y1_y2': 'lambdaT3',
    'y2_y2': 'lambdaT3ExtraQ020020I1',
}
AUD = Path('output/audits')
MATCHING = Path('output/general_t3_mass_eigen_matchete')
SOURCE_OPS = Path('output/scalar_relevance/operator_classification.json')

class CompletenessError(RuntimeError):
    pass

def need(cond, message):
    if not cond:
        raise CompletenessError(message)

def read_text(root, path):
    p = root / path
    need(p.is_file(), 'Missing local source: ' + str(p))
    return p.read_text(encoding='utf-8-sig')

def read_json(root, path):
    return json.loads(read_text(root,path))

def hashfile(root, path):
    return sha256((root/path).read_bytes()).hexdigest()

def import_existing(root, filename, name):
    path = root / filename
    need(path.is_file(), 'Missing existing project module: ' + str(path))
    spec = importlib.util.spec_from_file_location(name, path)
    need(spec is not None and spec.loader is not None, 'Cannot import ' + str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def conjugate(counts):
    return tuple(counts[i] for i in (1,0,3,2,5,4))

def direct_triangle_channel(counts):
    """Two external Hdag, two scalar legs conjugate to y1:S1 / y2:S2dag.

    Checks BOTH Hermitian orientations and rejects ambiguous vertices.
    This is a necessary Wick contraction condition, not a CG proof.
    """
    need(len(counts)==6 and all(isinstance(x,int) and x>=0 for x in counts)
         and sum(counts)==4, 'Bad quartic species counts')
    targets={(0,2,0,2,0,0):'y1_y1',
             (0,2,0,1,1,0):'y1_y2',
             (0,2,0,0,2,0):'y2_y2'}
    possibilities=[(name,targets[orientation]) for name,orientation in
                   (('native',tuple(counts)),('h.c.',conjugate(counts)))
                   if orientation in targets]
    need(len(possibilities)<=1,'Ambiguous Hermitian orientations')
    return possibilities[0] if possibilities else None

def classify_model(operators, model):
    selected = [o for o in operators if o.get('model')==model]
    need(selected,'No operator manifest for '+model)
    kinds=Counter(o['kind'] for o in selected)
    need(kinds=={'baseline_topology':1,'baseline_quadratic':2,
                 'baseline_quartic_family':5,'additional_quartic':7,
                 'additional_quadratic':1},
         model+': unexpected UV bosonic operator inventory '+str(kinds))
    diagrams=[]; screened=[]; mass=[]
    for o in selected:
        kind=o['kind']; name=o['coupling']
        # The project's actual 16-model operator manifest intentionally uses
        # JSON null for baseline masses and baseline quartic *families*.  Only
        # vertices with explicit species multiplicities can be Wick-counted.
        # Requiring counts for those three single-vertex classes avoids ever
        # interpreting a missing physical quartic/mixed-bilinear signature as
        # a screened zero contribution.
        raw_counts=o.get('counts_H_Hdag_S1_S1dag_S2_S2dag')
        if kind in ('baseline_topology','additional_quartic','additional_quadratic'):
            need(isinstance(raw_counts,list) and len(raw_counts)==6 and
                 all(type(x) is int and x>=0 for x in raw_counts),
                 model+': invalid explicit boson species counts for '+name)
            counts=tuple(raw_counts)
        else:
            need(raw_counts is None,
                 model+': baseline family changed from expected null inventory: '+name)
            counts=None
        if kind in ('baseline_topology','additional_quartic'):
            if sum(counts[:2])==2:
                candidate=direct_triangle_channel(counts)
                if candidate is not None:
                    orient,channel=candidate
                    need(EXPECTED[channel]==name,model+': new or missing quartic for '+channel+': '+name)
                    need((kind=='baseline_topology' and o.get('classification')=='leading_one_loop_topology')
                         or (kind=='additional_quartic' and o.get('observed_in_saved_C5') is True),
                         model+': required quartic not seen in archived C5: '+name)
                    diagrams.append({'channel':channel,'vertex':name,
                                     'vertex_orientation':orient,
                                     'external':['l','l','Hdag','Hdag'],
                                     'internal_lines':['F','scalar','scalar'],
                                     'V':3,'I':3,'loops':1,
                                     'quartic_species_counts':counts})
                else:
                    screened.append({'vertex':name,'reason':'Higgs-or-scalar Wick orientation',
                                     'kind':kind})
            else:
                screened.append({'vertex':name,'reason':'No pair of Higgs legs',
                                 'kind':kind})
        elif kind=='additional_quadratic':
            mass.append({'vertex':name,'counts':counts,'status':'unexcluded tree-level mass insertion'})
        elif kind=='baseline_quartic_family':
            screened.append({'vertex':name,'reason':'portal or scalar quartic not direct holomorphic single-vertex triangle',
                             'kind':kind})
    need(len(diagrams)==3 and {d['channel'] for d in diagrams}==set(EXPECTED),
         model+': three exclusive triangle channels not reproduced')
    need(len(mass)==1 and mass[0]['counts']==(0,0,0,1,0,1),
         model+': unaccounted additional scalar quadratic')
    for d in diagrams:
        need(d['I']-d['V']+1==d['loops']==1,'Incorrect topological loop count')
    return {'model':model,'UV_operator_family_counts':dict(kinds),
            'one_loop_triangle_graphs':sorted(diagrams,key=lambda a:a['channel']),
            'direct_quartic_wick_sector_exhaustive':True,
            'not_direct_in_this_graph':screened,
            'additional_tree_quadratic_not_excluded':mass,
            'source_quadratic_requires_zero_boundary_or_resummed_verified_mass_basis':True}

# This parses saved Matchete's expanded expressions, not source-code Mathematica.
# The complete F and Y loop-differences are scanned to detect missing sectors.
# No Fierz or other operator-basis transformation is claimed.
H = re.compile(r'Field\[H, Scalar')
L = re.compile(r'Field\[l, Fermion')
COUPLING = re.compile(r'Coupling\[(lambdaT3(?:ExtraQ[A-Za-z0-9]+)?), \{\}, 0\]')
Y = re.compile(r'Coupling\[(y[12]), \{')

def scan_llhh(raw, *, term_parser, stage, expect_sector):
    ts=term_parser(raw)
    candidates=[]; other=[]
    for t in ts:
        nh=len(H.findall(t)); nl=len(L.findall(t))
        has_new=any('Field['+x in t for x in ('NewScalarX','NewScalarY','NewFermion'))
        if nh==2 and nl==2 and not has_new:
            if 'GammaCC' not in t:
                other.append(t[:160])
            else:
                candidates.append(t)
    by_kind=Counter()
    quartics=Counter(); yukawas=Counter()
    for t in candidates:
        qs=set(COUPLING.findall(t))
        ys=Y.findall(t)
        need(len(qs)==1 and next(iter(qs)) in EXPECTED.values(),
             stage+': missing/unknown quartic in LLHH; '+str(qs))
        need(len(ys)==2 and all(y in ('y1','y2') for y in ys),
             stage+': unexpected LLHH Yukawa structure '+str(ys))
        need(t.count('GammaCC')==1 and (('Proj[-1]' in t) != ('Proj[1]' in t)),
             stage+': ambiguous Majorana bilinear')
        quartics[next(iter(qs))]+=1
        yukawas['_'.join(sorted(ys))]+=1
        c='pole' if 'ϵ' in t else 'log' if 'Log[' in t else 'constant'
        need(not ('ϵ' in t and 'Log[' in t),stage+': mixed pole/log')
        by_kind[c]+=1
    need(len(candidates)==expect_sector,
         stage+f': {len(candidates)} selected LLHH GammaCC terms; expected {expect_sector}')
    need(not other,stage+': new LLHH non-GammaCC structures require inspection')
    need(set(quartics)==set(EXPECTED.values()),stage+': quartic channel coverage incomplete')
    need(set(yukawas)=={'y1_y1','y1_y2','y2_y2'},stage+': Yukawa coverage incomplete')
    return {'source_additive_terms':len(ts),'selected_llhh_GammaCC_terms':len(candidates),
            'pole_log_constant_counts':dict(by_kind),
            'quartic_term_counts':dict(sorted(quartics.items())),
            'Yukawa_channel_term_counts':dict(sorted(yukawas.items())),
            'unclassified_LLHH_without_GammaCC':len(other),
            'operator_basis_Fierz_completeness_certified':False}

def source_paths(stem):
    fdir=MATCHING/'eft1_fermion_first_oneloop'
    ydir=MATCHING/'eft1_scalar_Y_only_canonical_probe'
    xdir=MATCHING/'eft1_scalar_xy_fresh_kernel_v3'
    fstem=f'{stem}_final_eft1_fermion_first_oneloop'
    ystem=f'{stem}_final_y_only_canonical_match'
    xstem=f'{stem}_final_x_from_saved_y_match_x_followup'
    return {
        'F_full': fdir/f'{fstem}_eft1_one_loop_correction.inputform.txt',
        'F_raw': fdir/f'{fstem}_raw_hard_llhh_sector.inputform.txt',
        'F_report': fdir/f'{fstem}.json',
        'Y_full': ydir/f'{ystem}_scalar_one_loop_correction.inputform.txt',
        'Y_raw': ydir/f'{ystem}_raw_scalar_llhh.inputform.txt',
        'Y_report': ydir/f'{ystem}.json',
        'X_full': xdir/f'{xstem}_x_loop_difference.inputform.txt',
        'X_raw': xdir/f'{xstem}_x_raw_llhh.inputform.txt',
        'X_report': xdir/f'{xstem}.json',
    }

def verify_uv_code(root):
    paths=['Lagrangian/model/T3ModelCatalog.wl',
           'Lagrangian/model/LagrangianBuilder.wl',
           'Lagrangian/interactions/T3Topology.wl',
           'Lagrangian/interactions/GeneralScalarOperators.wl',
           'Lagrangian/RunGeneralT3FermionFirstEFT1OneLoop.wl']
    content={p:read_text(root,Path(p)) for p in paths}
    checks={
        paths[0]:['"Scalar1" -> alpha/2','"Scalar2" -> (alpha + 2)/2',
                  '"Fermion" -> (alpha + 1)/2'],
        paths[1]:['extraYExpected = T3GeneralExtraYukawaCount[model]',
                  'LInt = Total[Lookup[assembled, "Expression", {}]]'],
        paths[2]:['useCConj = which === 1;', 'scalarBar = which === 2;',
                  'T3GeneralExtraYukawaCount','BuildT3ExtraYukawaCandidates',
                  'lambdaT3[] H[i] H[j] s1 s2'],
        paths[3]:['T3ExtendedBuildExtraQuarticCandidates','T3ExtendedBuildCandidates',
                  'T3ExtendedExtraQuadraticStructures'],
        paths[4]:['oneLoopMatchingResult=CheckAbort[RunT3Matching[lag,5,1]',
                  'StringContainsQ[s,"GammaCC"]','hardCandidates=Select[loopTerms'],
    }
    for p,fragments in checks.items():
        for frag in fragments:
            need(frag in content[p],'UV source contract changed: '+p+' / '+frag)
    return {p:hashfile(root,Path(p)) for p in paths}

def run(root,output):
    old=import_existing(root,Path('scripts/AuditT3OneLoopScalarVertices.py'),'t3_v20_prior_wick_audit')
    parser=import_existing(root,Path('Numerical/experimental/T3BCPoleLogConsistency.py'),'t3_v20_terms')
    base=read_json(root,SOURCE_OPS)
    old.validate_batch1(base)
    old.verify_source_hashes(base,root)
    uvhash=verify_uv_code(root)
    old_wick=read_json(root,Path('output/scalar_relevance/leading_c5_vertex_audit.json'))
    need(old_wick['status']=='WickNecessaryConditionsReproduceSavedC5Inventory',
         'Prior Wick classification changed')
    need(old_wick.get('batch1_report_sha256')==hashfile(root,SOURCE_OPS),
         'Old vertex audit does not correspond to current operator manifest')
    v19=read_json(root,AUD/'t3_bc_hard_soft_subtraction_v19/summary.json')
    need(v19.get('format')=='t3-bc-one-loop-hard-soft-msbar-restricted-v19' and
         v19.get('full_T3_renormalized_matching_certified') is False,
         'v19 source/scope changed')
    warsaw=read_json(root,AUD/'t3_bc_y_warsaw_direct_identity_v11/native_identity.json')
    need(warsaw.get('probe_revision')=='11.2' and
         warsaw.get('direct_equivalence_status')=='ExactPlusEquality' and
         warsaw.get('positive_sign_equality') is True and
         warsaw.get('negative_sign_equality') is False,
         'Matchete v11.2 operator-basis normalization no longer verified')
    retained=read_json(root,Path('output/scalar_relevance/final_leading_scalar_potential.json'))
    need(retained.get('conditions',{}).get('additional_quadratics_zero_at_uv_boundary') is True,
         'Missing benchmark zero-extra-quadratic boundary assumption')
    result={'format':'t3-bc-one-loop-llhh-topology-completeness-v20',
            'status':'RestrictedY2LambdaOneLoopGraphInventoryAudited',
            'scope':'unbroken gauge phase, two original T3 Yukawas and one HHSS quartic, LNV Majorana F, diagonal tree-level scalar propagators',
            'source_contract_sha256':uvhash,
            'preexisting_Wick_audit_reused_and_source_hashes_verified':True,
            'new_graph_enumeration_performed':True,
            'renormalizable_one_loop_topological_valence_budget': {
               'external_legs': 4, 'I_minus_V_plus_one': 1,
               'sum_over_vertices_valence_minus_two': 4,
               'two_required_BSM_Yukawas_contribute': 2,
               'remaining_valence_budget': 2,
               'combinatorial_classes': [
                  {'extra_vertices':'one quartic','status':'enumerated_all_HHSS_in_operator_manifest'},
                  {'extra_vertices':'two trivalent','status':'SM_and_BSM_single_H_attachments_not_exhaustively_certified'},
               ],
               'quadratic_mass_insertions':'arbitrary number at same graph-loop order; require propagator resummation'
            },
            'full_matchete_loop_source_files_scanned':True,
            'three_triangle_channels_exhaustive_within_specified_vertex_truncation':True,
            'all_one_loop_SM_BSM_diagrams_exhaustively_generated':False,
            'two_single_H_trivalent_vertex_graphs_exhaustively_excluded':False,
            'extra_scalar_quadratic_nonzero_allowed_generally':True,
            'assumed_zero_extra_quadratic_at_this_benchmark_boundary':True,
            'no_extra_original_Yukawas_BC_alpha_minus_one':True,
            'one_loop_field_counterterm_times_C5_tree_necessary':False,
            'physical_C5_certified':False,'physical_C5':None,
            'models':{},'remaining_certification':[
             'The remaining two-trivalent-vertex one-loop graph class needs an explicit SM+BSM charge/chirality connectivity exclusion, including SM charged-lepton Yukawa and gauge interactions',
             'Independent UV-to-EFT coefficient comparison for all allowed tensor and EOM structures, not just raw GammaCC LLHH selection',
             'Nonzero additional quadratic mass mixing must be incorporated in the propagator/eigenbasis before completeness can be extended',
             'Renormalized EFT1/EFT2 running and finite scheme/field/parameter maps not certified',
            ]}
    for stem,model in STEMS.items():
        data=classify_model(base['operators'],model)
        path=source_paths(stem)
        meta={k:read_json(root,p) for k,p in path.items() if k.endswith('_report')}
        need(meta['F_report'].get('source_yukawa_count')==2 and
             meta['F_report'].get('source_quartic_count')==17 and
             meta['F_report'].get('matching_loop_order')==1,
             model+': F matching source census changed')
        need(meta['F_report'].get('raw_hard_llhh_candidate_term_count')==204 and
             meta['Y_report'].get('raw_scalar_llhh_term_count')==132 and
             meta['X_report'].get('x_tree_matching_completed') is True,
             model+': F/Y/X source inventory changed')
        need(meta['X_report'].get('y_loop_inherited_in_x_matching') is False,
             model+': X loop source is not the known tree-only follow-up')
        expected={'F':204,'Y':132,'X':36}
        audits={}; hashes={}
        for stage in ('F','Y','X'):
            for suffix in ('raw','full'):
                p=path[stage+'_'+suffix]
                source=read_text(root,p)
                summary=scan_llhh(source,term_parser=parser.terms,stage=model+'/'+stage+'/'+suffix,
                                  expect_sector=expected[stage])
                audits[stage+'_'+suffix]=summary
                hashes[stage+'_'+suffix]=hashfile(root,p)
            need(audits[stage+'_raw']['quartic_term_counts']==audits[stage+'_full']['quartic_term_counts'],
                 model+': mismatch in quartic coverage of full vs exported '+stage+' matching')
            need(audits[stage+'_raw']['Yukawa_channel_term_counts']==audits[stage+'_full']['Yukawa_channel_term_counts'],
                 model+': mismatch in Yukawa coverage of full vs exported '+stage+' matching')
        need(audits['F_full']['source_additive_terms']==meta['F_report']['loop_difference_term_count'],
             model+': full F matched EFT one-loop term inventory changed')
        need(audits['Y_full']['source_additive_terms']==meta['Y_report']['scalar_loop_difference_term_count'],
             model+': full Y matched EFT one-loop term inventory changed')
        exact_counts={'F':{'pole':68,'log':68,'constant':68},
                      'Y':{'pole':50,'log':50,'constant':32},
                      'X':{'pole':18,'log':18}}
        for stage,ct in exact_counts.items():
            for suffix in ('full','raw'):
                observed=audits[stage+'_'+suffix]['pole_log_constant_counts']
                need(all(observed.get(k,0)==val for k,val in ct.items()) and
                     sum(observed.values())==sum(ct.values()),
                     model+': raw/full pole-log-constant partition changed for '+stage+'/'+suffix)
        need(v19['models'][stem].get('symbolic_f_hard_plus_eft_soft_poles_zero_in_each_yukawa_structure') is True,
             model+': v19 pole+subtraction evidence failed')
        need(data['direct_quartic_wick_sector_exhaustive'],model+': graph basis incomplete')
        data.update({'F_Y_X_matching_sector_scan':audits,
                     'source_matching_sha256':hashes,
                     'source_verified_three_channels_are_present_at_every_boundary':True,
                     'restricted_one_loop_hard_soft_identity_v19':True,
                     'complete_BC_leading_one_loop_C5_certified':False})
        result['models'][stem]=data
        print(f'{model}: three vertex-level triangle channels PASS; full F/Y/X matched-EFT LLHH source extraction PASS')
        print(f'{model}: F={audits["F_full"]["selected_llhh_GammaCC_terms"]}, '
              f'Y={audits["Y_full"]["selected_llhh_GammaCC_terms"]}, '
              f'X={audits["X_full"]["selected_llhh_GammaCC_terms"]} additive terms')
    out=output if output.is_absolute() else root/output
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: y_i*y_j*lambda single-quartic topology completeness within declared B/C source sector')
    print('FULL UV/EFT one-loop physical C5: NOT CERTIFIED')
    print('Saved: '+str(out))

def self_test():
    for counts,name in (((2,0,1,0,0,1),'y1_y2'),
                        ((0,2,0,2,0,0),'y1_y1'),
                        ((0,2,0,0,2,0),'y2_y2')):
        r=direct_triangle_channel(counts)
        need(r is not None and r[1]==name,'Failed graph orientation '+name)
    for counts in ((1,1,1,0,0,1),(0,2,1,0,0,1),(0,0,0,2,2,0)):
        need(direct_triangle_channel(counts) is None,'Wrongly retained forbidden quartic')
    def parts(s):
        return tuple(s.split(' + '))
    w=('hbar*Coupling[lambdaT3, {}, 0]*Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*'
       'Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]*'
       'Field[H, Scalar, {a}, {}]*Field[H, Scalar, {b}, {}]*'
       'NCM[Field[l, Fermion, {x}, {}], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {y}, {}]]')
    # Other two channels are also present in source, not just the mixed term.
    synthetic=' + '.join([w,w.replace('lambdaT3,','lambdaT3ExtraQ020200I1,').replace('Coupling[y2,','Coupling[y1,'),
                              w.replace('lambdaT3,','lambdaT3ExtraQ020020I1,').replace('Coupling[y1,','Coupling[y2,')])
    r=scan_llhh(synthetic,term_parser=parts,stage='self-test',expect_sector=3)
    need(r['selected_llhh_GammaCC_terms']==3,'Synthetic scan failed')
    for x in ('GammaCC','lambdaT3,'):
        try:
            scan_llhh(synthetic.replace(x,'BAD',1),term_parser=parts,stage='tamper',expect_sector=3)
        except CompletenessError:
            continue
        raise CompletenessError('Tampered source not rejected: '+x)
    # Regression for the actual GitHub operator_classification.json schema:
    # both baseline quadratic masses and all five baseline quartic families
    # have counts=null; individual physical quartics and the mixed quadratic
    # always supply six integer species counts.
    def entry(kind,name,counts=None,classification=None,observed=None):
        record={'model':'B_-1','kind':kind,'coupling':name,
                'counts_H_Hdag_S1_S1dag_S2_S2dag':counts}
        if classification is not None:record['classification']=classification
        if observed is not None:record['observed_in_saved_C5']=observed
        return record
    mock=[entry('baseline_topology','lambdaT3',[2,0,1,0,0,1],
                'leading_one_loop_topology')]
    mock += [entry('baseline_quadratic',n) for n in ('MS1','MS2')]
    mock += [entry('baseline_quartic_family',n) for n in
             ('lambdaS1Inv*','lambdaS2Inv*','lambdaH1Inv*',
              'lambdaH2Inv*','lambda12Inv*')]
    mock += [entry('additional_quartic',n,c,'observed_in_saved_one_loop_C5',True)
             for n,c in (('lambdaT3ExtraQ020020I1',[0,2,0,0,2,0]),
                         ('lambdaT3ExtraQ020200I1',[0,2,0,2,0,0]))]
    mock += [entry('additional_quartic',n,c,'not_direct',False) for n,c in (
        ('lambdaT3ExtraQ000112I1',[0,0,0,1,1,2]),
        ('lambdaT3ExtraQ000202I1',[0,0,0,2,0,2]),
        ('lambdaT3ExtraQ001201I1',[0,0,1,2,0,1]),
        ('lambdaT3ExtraQ110101I1',[1,1,0,1,0,1]),
        ('lambdaT3ExtraQ110101I2',[1,1,0,1,0,1]))]
    mock.append(entry('additional_quadratic','mu2_extra_000101_I1',
                      [0,0,0,1,0,1]))
    model_audit=classify_model(mock,'B_-1')
    need(len(model_audit['one_loop_triangle_graphs'])==3 and
         len(model_audit['additional_tree_quadratic_not_excluded'])==1,
         'Actual-manifest null baseline schema regression')
    for wrong in ('additional_quartic','baseline_topology','additional_quadratic'):
        altered=[dict(v) for v in mock]
        victim=next(v for v in altered if v['kind']==wrong)
        victim['counts_H_Hdag_S1_S1dag_S2_S2dag']=None
        try: classify_model(altered,'B_-1')
        except CompletenessError:pass
        else:raise CompletenessError('Required counts-null mutation accepted: '+wrong)
    print('PASS: actual manifest null baseline fields; required quartic and quadratic signatures remain strict')
    print('PASS: independent quartic Wick orientations and three graph channels')
    print('PASS: full matched-source LLHH selection, missing quartic and bilinear tamper rejected')
    need(2*(3-2)+(4-2)==4,'Wrong one-loop renormalizable valence budget')
    need(2*(3-2)+2*(3-2)==4,'Two-trivalent branch omitted')
    print('PASS: renormalizable one-loop graph valence partitions (one quartic OR two trivalent)')
    print('PASS: no complete physical C5 certification')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,default=AUD/'t3_bc_one_loop_topology_completeness_v20/summary.json')
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--require-physical',action='store_true')
    a=p.parse_args()
    try:
        need(not a.require_physical,'Full physical C5 cannot be certified by topology counting')
        if a.self_test:self_test()
        else:run(a.root.resolve(),a.output)
        return 0
    except (CompletenessError,ValueError,TypeError,OSError,KeyError,AttributeError) as exc:
        print('FAIL (no physical C5): '+str(exc),file=sys.stderr)
        return 2

if __name__=='__main__':
    raise SystemExit(main())
