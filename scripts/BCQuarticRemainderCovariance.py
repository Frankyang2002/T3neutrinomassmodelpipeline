"""Stage 24 read-only B/C SU(2) audit of four additionally reused quartics.

Source: fresh Stage22 Matchete/RGBeta native exports and Stage23 provenance.
Tests canonical U / conjugate(U) action on four fields, independent
proposed del/eps substitutions, and exact native Matchete tensor spans.
No RGBeta/Matchete production source is modified or certified.

python -m scripts.BCQuarticRemainderCovariance
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import product
import json
from pathlib import Path

import sympy as sp

from scripts.RGBetaExplicitBasis import parse_fields, eval_component, SOURCE_ORDER
from scripts.RGBetaAlignedGeometry import reindexed_doublet_tensor
from scripts.RGBetaTensorAlgebra import solve_basis, TensorCertificateError, exact
from scripts.QuarticPrefactorAudit import FOLDER, run as run_stage23, compact

MODELS = ('B_-1', 'C_-1')
# A *source fingerprint*: refuse unreviewed modifications or stale exports.
OLD = {
 'lambdaS1barS2S2bar2': '(del[SU2L[fund],t3i1,t3i4]*del[SU2L[fund],t3i2,t3i3])/2+(del[SU2L[fund],t3i1,t3i3]*del[SU2L[fund],t3i2,t3i4])/2',
 'lambdaS1bar2S2bar2': '(del[SU2L[fund],t3i1,t3i4]*del[SU2L[fund],t3i2,t3i3])/2+(del[SU2L[fund],t3i1,t3i3]*del[SU2L[fund],t3i2,t3i4])/2',
 'lambdaS1S1bar2S2bar': '(del[SU2L[fund],t3i1,t3i3]*del[SU2L[fund],t3i2,t3i4])/2+(del[SU2L[fund],t3i1,t3i2]*del[SU2L[fund],t3i3,t3i4])/2',
 'lambdaHHdagS1barS1bar': '(del[SU2L[fund],t3i1,t3i4]*del[SU2L[fund],t3i2,t3i3])/2+(del[SU2L[fund],t3i1,t3i3]*del[SU2L[fund],t3i2,t3i4])/2',
}
# Tensors proposed for correction; functions are tested against the native
# Matchete component basis, NOT inferred to be registered in current RGBeta.
PROPOSED = {
 'lambdaS1barS2S2bar2': '(eps[SU2L[fund],t3i1,t3i4]*del[SU2L[fund],t3i2,t3i3]+eps[SU2L[fund],t3i1,t3i3]*del[SU2L[fund],t3i2,t3i4])/2',
 'lambdaS1bar2S2bar2': '(eps[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]+eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4])/2',
 'lambdaS1S1bar2S2bar': '(del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]+del[SU2L[fund],t3i1,t3i2]*eps[SU2L[fund],t3i3,t3i4])/2',
 'lambdaHHdagS1barS1bar': '(eps[SU2L[fund],t3i1,t3i4]*eps[SU2L[fund],t3i2,t3i3]+eps[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4])/2',
}
FIELDS = {
 'lambdaS1barS2S2bar2': ['S1dag','S2','S2dag','S2dag'],
 'lambdaS1bar2S2bar2': ['S1dag','S1dag','S2dag','S2dag'],
 'lambdaS1S1bar2S2bar': ['S1','S1dag','S1dag','S2dag'],
 'lambdaHHdagS1barS1bar': ['Hdag','Hdag','S1dag','S1dag'],
}
GENS = (sp.Matrix([[0,1],[1,0]])/2,
        sp.Matrix([[0,-sp.I],[sp.I,0]])/2,
        sp.diag(sp.Rational(1,2), -sp.Rational(1,2)))
INDEXES = tuple(product(range(2), repeat=4))


def fail(msg):
    raise ValueError('Stage24 fail-closed: '+msg)


def components(expr):
    return {i:sp.simplify(eval_component(expr,i)) for i in INDEXES}


def covariance(expr, fields):
    if len(fields) != 4 or any(x not in SOURCE_ORDER for x in fields):
        fail('Unsupported field representation for B/C SU2 covariance')
    tensor = components(expr)
    failures = []
    witnesses = []
    for g in GENS:
        gs = tuple(-sp.conjugate(g) if f.endswith('dag') else g for f in fields)
        residual = {}
        for i in INDEXES:
            val = sp.simplify(sum(gs[axis][j,i[axis]] * tensor[i[:axis]+(j,)+i[axis+1:]]
                                  for axis in range(4) for j in range(2)))
            if val != 0: residual[i] = val
        failures.append(len(residual))
        witnesses.append({'indices_0_based':list(next(iter(residual))),
                          'exact_residual':str(next(iter(residual.values())))} if residual else None)
    return {'su2_covariant': failures == [0,0,0],
            'nonzero_residuals_T1_T2_T3':failures,
            'first_counterexamples_T1_T2_T3':witnesses}


def group_blocks(mat_model, inv_model):
    m=defaultdict(list); r=defaultdict(list)
    for t in mat_model['tensors']:
        m[''.join(map(str,t['counts']))].append(t)
    for t in inv_model['invariants']:
        species=parse_fields(t['native_fields_inputform'])
        r[''.join(str(species.count(k)) for k in SOURCE_ORDER)].append(t)
    return m,r


def global_doublet_scan(m,r):
    # A *restricted* but physically coherent test: each complex doublet and its
    # conjugate are reindexed by the same real epsilon matrix J. No independent
    # field-specific, block-specific, or arbitrary unitary transformations.
    scores=[]
    for bits in product((0,1), repeat=3):
        selected=tuple(field for i,field in enumerate(SOURCE_ORDER) if bits[i//2])
        aligned=[]
        for block in sorted(m):
            if any(t['dimensions'] != [2,2,2,2] for t in m[block]):
                fail('Unexpected non-doublet in B/C global scan')
            cols=[]
            for v in r[block]:
                v=dict(v)
                if v['symbol']=='lambdaHHdagS1barS2barCross':
                    v['native_invariant_expression_inputform']='del[SU2L[fund],t3i1,t3i3]*eps[SU2L[fund],t3i2,t3i4]'
                cols.append(reindexed_doublet_tensor(v,m[block][0],selected))
            try:
                solve_basis(m[block],cols)
                aligned.append(block)
            except (TensorCertificateError,ValueError):
                pass
        scores.append({'rotated_species':list(selected),'aligned_blocks':aligned})
    return {'choices':len(scores),'maximum_aligned_blocks':max(len(x['aligned_blocks']) for x in scores),
            'common_global_basis_exists':any(len(x['aligned_blocks'])==len(m) for x in scores),
            'witnesses':scores}


def calculate(mat, inv, stage23):
    if stage23.get('totals',{}).get('slots') != 50 or stage23.get('totals',{}).get('blocks') != 46 or stage23.get('signed_physical_coupling_map_certified') is not False:
        fail('Expected Stage23 uncertified 50/46 evidence')
    results={}
    for model in MODELS:
        m,r=group_blocks(mat['models'][model],inv['models'][model])
        rr={v['symbol']:v for records in r.values() for v in records}
        if len(m) != 6 or sum(map(len,r.values())) != 7 or len(rr) != 7:
            fail(model+': native six-block/seven-slot inventory changed')
        checks=[]
        for symbol in OLD:
            rec=rr.get(symbol)
            if not rec or parse_fields(rec['native_fields_inputform']) != FIELDS[symbol]:
                fail(model+'/'+symbol+': field order or registration changed')
            if compact(rec['native_invariant_expression_inputform']) != OLD[symbol]:
                fail(model+'/'+symbol+': stale or unexpectedly changed source expression')
            old=covariance(OLD[symbol],FIELDS[symbol]); new=covariance(PROPOSED[symbol],FIELDS[symbol])
            if old['su2_covariant'] or not new['su2_covariant']:
                fail(model+'/'+symbol+': unexpected SU2 covariance audit')
            key=''.join(str(FIELDS[symbol].count(f)) for f in SOURCE_ORDER)
            mats=m[key]
            if len(mats) != 1:
                fail(model+'/'+key+': unexpected Matchete basis count')
            candidate={**{k:mats[0][k] for k in ('counts','dimensions')},
                       'species':FIELDS[symbol],'symbol':symbol,
                       'components':[str(components(PROPOSED[symbol])[i]) for i in INDEXES]}
            geo=solve_basis(mats,[candidate])
            if geo['r_in_m'] != [['1']]:
                fail(model+'/'+symbol+': proposed replacement changes expected Matchete normalisation')
            twisted=reindexed_doublet_tensor(rec,mats[0],('S1','S1dag'))
            if [exact(x) for x in candidate['components']] != [exact(x) for x in twisted['components']]:
                fail(model+'/'+symbol+': proposed tensor differs from prior J reindex witness')
            checks.append({'symbol':symbol,'block':key,'native':old,'proposed':new,
                           'proposed_invariant_expression':PROPOSED[symbol],
                           'exact_matchete_geometry_r_in_m':geo['r_in_m'],
                           'physical_coupling_certified':False})
        scan=global_doublet_scan(m,r)
        if scan['common_global_basis_exists'] or scan['maximum_aligned_blocks'] != 5:
            fail(model+': global coherent I/J obstruction changed')
        results[model]={'quartics':checks,'global_basis_scan':scan}
    return {'format':'stage24-bc-reused-quartic-su2-covariance-v1',
            'status':'AdditionalReusedInvariantCovarianceFailuresAndConditionalRepairCandidates',
            'source_phase':'after 110101 fix, before four additional native fixes',
            'field_convention':'fundamental U; conjugate fields conjugate(U); epsilon[0,1]=+1',
            'models':results,
            'totals':{'models':2,'old_non_covariant_slots':8,'proposed_covariant_slots':8,
                      'proposed_exact_matchete_1d_matches':8,
                      'restricted_global_I_or_J_basis_exists':False},
            'global_su2_component_basis_verified':False,
            'live_wolfram_candidate_covariance_checked':False,
            'physical_coupling_mapping_certified':False,'production_modified':False}


def run(base=FOLDER, output=None):
    stage23=run_stage23(base)
    names=('matchete_tensors.json','rgbeta_native_invariants.json')
    evidence={f:json.loads((base/f).read_text(encoding='utf-8-sig')) for f in names}
    res=calculate(evidence[names[0]],evidence[names[1]],stage23)
    res['input_sha256']={f:hashlib.sha256((base/f).read_bytes()).hexdigest() for f in names}
    if output:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(res,indent=2)+'\n',encoding='utf-8')
    return res


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-dir',type=Path,default=FOLDER)
    ap.add_argument('--output',type=Path,default=FOLDER/'stage24_bc_reused_covariance.json')
    args=ap.parse_args()
    try:res=run(args.input_dir,args.output)
    except (KeyError,IndexError,ValueError,TypeError,OSError,TensorCertificateError) as ex:ap.exit(2,str(ex)+'\n')
    print('Stage 24 B/C native SU2 failures:',res['totals']['old_non_covariant_slots'],'of 8')
    print('Proposed invariant SU2 + exact Matchete matches:',res['totals']['proposed_covariant_slots'],'of 8')
    print('Global coherent I/J reindexing:',res['totals']['restricted_global_I_or_J_basis_exists'])
    print('Physical coupling mapping: UNCERTIFIED')
    print('Report:',args.output)

if __name__=='__main__':main()
