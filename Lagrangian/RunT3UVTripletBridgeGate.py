#!/usr/bin/env python3
"""Workstream A: live Matchete triplet intertwiner + independent RGBeta coordinate gate.

Read 16 already-verified live Matchete CG reports; derive exact triplet representation,
metric and basis intertwiner from six HH-triplet tensors (A/D). Independently ask
installed RGBeta for numerical delS2 and tGen indices. Never manufacture RGBeta
component values or promote a conditional quartic map into physical certification.

python Lagrangian/RunT3UVTripletBridgeGate.py --self-test
python Lagrangian/RunT3UVTripletBridgeGate.py --all
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from Lagrangian import RunT3UVNativeBridge as NB

CG_ROOT = ROOT / 'Reports/output/workstream_a/native_cg_bridge'
OUT = ROOT / 'Reports/output/workstream_a/triplet_coordinate_gate'
SCRIPT = ROOT / 'Lagrangian/RunT3RGBetaTripletCoordinateProbe.wl'
STATUS_BLOCKED = 'LiveMatcheteTripletRepresentationVerifiedRGBetaCoordinatesUnavailable'
STATUS_CANDIDATE = 'CrossPackageTripletCoordinateCandidateVerifiedNotPhysicalDictionary'
MODEL_KEYS = ('A_m4','A_m2','A_p0','D_m2','D_p0','D_p2')
E_KEYS = ('E_m4','E_m2','E_p0','E_p2')
I = sp.I
SQRT2 = sp.sqrt(2)
ZERO3 = sp.zeros(3)
E = sp.Matrix([[1, 0, 0], [0, 1/SQRT2, 0],
               [0, 1/SQRT2, 0], [0, 0, 1]])
JREF = sp.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
TF = (
    sp.Matrix([[0,1],[1,0]])/2,
    sp.Matrix([[0,-I],[I,0]])/2,
    sp.diag(1,-1)/2,
)
D = tuple(sp.kronecker_product(t.T,sp.eye(2)) +
          sp.kronecker_product(sp.eye(2),t.T) for t in TF)


def need(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError('Triplet coordinate gate: '+message)


def zero(expr: Any) -> bool:
    if isinstance(expr, sp.MatrixBase):
        return all(sp.simplify(v) == 0 for v in expr)
    return sp.simplify(expr) == 0


def number(s: str) -> sp.Expr:
    return NB.exact_number(s)


def source_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def mat_json(mat: sp.MatrixBase) -> list[list[str]]:
    return [[str(sp.simplify(mat[i,j])) for j in range(mat.cols)] for i in range(mat.rows)]


def verify_triplet(q: sp.MatrixBase) -> dict[str,Any]:
    """Triplet in the dual contraction of two ordinary (unbarred) H doublets.

    Q[(i,j),a] = CG[H_i,H_j,(S or Sdag)_a]. Literal field index order is
    supplied by the saved live Matchete topology tensor, not a guessed CG.
    D is the transposed fundamental generator acting on two covariant H slots.
    """
    need(q.shape == (4,3), '2x2x3 HH-triplet tensor required')
    need(any(not zero(v) for v in q), 'zero topology tensor')
    need(all(zero(q[1,k]-q[2,k]) for k in range(3)),
         'topology tensor is not HH-symmetric')
    gram = sp.simplify(q.H*q)
    g = sp.simplify(gram[0,0])
    need(g.is_positive is True and zero(gram-g*sp.eye(3)),
         'native HH-triplet Gram matrix is not positive scalar identity')
    f = q/sp.sqrt(g)
    # E is a declared independent orthonormal symmetric-doublet reference.
    u = sp.simplify(E.H*f)
    need(zero(u.H*u-sp.eye(3)) and zero(E*u-f),
         'native triplet subspace not related to symmetric pair by a unitary')
    reps = tuple(sp.simplify(-f.H*d*f) for d in D)
    for d,r in zip(D,reps):
        need(zero(d*f+f*r), 'native topology fails exact SU2 generator covariance')
        need(zero(r.H-r) and zero(r.trace()),
             'triplet generators not Hermitian/traceless')
    for a,b,c in ((0,1,2),(1,2,0),(2,0,1)):
        need(zero(reps[a]*reps[b]-reps[b]*reps[a]-I*reps[c]),
             'triplet generator Lie algebra failed')
    need(zero(sum((r*r for r in reps), ZERO3)-2*sp.eye(3)),
         'triplet quadratic Casimir != 2')
    for a in range(3):
        for b in range(3):
            need(zero(sp.trace(reps[a]*reps[b])-(2 if a==b else 0)),
                 'triplet Dynkin/trace normalization != 2 delta')
    metric = sp.simplify(u.T*JREF*u)
    need(zero(metric.T-metric) and zero(metric.H*metric-sp.eye(3)),
         'derived conjugate-intertwiner metric not symmetric unitary')
    for r in reps:
        need(zero(r.T*metric+metric*r),
             'triplet invariant bilinear metric does not intertwine conjugate')
    scalar_identity = zero(metric-metric[0,0]*sp.eye(3)) and not zero(metric[0,0])
    return {'native_topology_norm_squared':str(g),
            'reference_to_matchete_triplet_unitary':mat_json(u),
            'derived_cartesian_triplet_generators':[mat_json(r) for r in reps],
            'derived_triplet_symmetric_bilinear_metric':mat_json(metric),
            'bilinear_metric_scalar_identity_in_this_basis':bool(scalar_identity),
            'fundamental_generator_convention':'sigma_a/2; H covariant-index contraction',
            'unitarity_su2_commutator_casimir_metric_exact':True,
            'matchete_triplet_generator_independently_exported':False}


def q_from_live(report: dict[str,Any], label: str) -> sp.Matrix:
    _, tensors = NB.check_model(label,report)
    # Only A/D topology has exactly one triplet index: H_i H_j T_a.
    _,_,(d1,d2,_) = NB.model_specs(label)
    need(sorted((d1,d2))==[1,3],label+': not a single-triplet topology')
    values = tensors['TopologyHHSS']
    need(len(values)==12,label+': topology tensor should have 12 components')
    return sp.Matrix(4,3,list(values))


def e_two_triplet_covariance(report:dict[str,Any], label:str,
                             rd:tuple[sp.Matrix,...], ra:tuple[sp.Matrix,...])->dict[str,Any]:
    """Check E topology under native A/D-derived S1 and S2dag triplet action.

    This is an exact consistency probe across representation registrations,
    not an independent native RGBeta generator export. An obstruction is
    reported, not silently repaired by choosing a new field-basis phase.
    """
    _, tensors=NB.check_model(label,report)
    flat=tensors['TopologyHHSS']
    need(len(flat)==36, label+': E topology expected 2x2x3x3 components')
    ce=sp.Matrix(4,9,list(flat))
    rr=tuple(sp.kronecker_product(rd[k],sp.eye(3)) +
             sp.kronecker_product(sp.eye(3),ra[k]) for k in range(3))
    residuals=[sp.simplify(D[k]*ce + ce*rr[k]) for k in range(3)]
    nonzero=[sum(not zero(v) for v in r) for r in residuals]
    return {'model':label, 'components_checked':36,
            'common_matchete_A_D_triplet_generators_su2_covariant':all(x==0 for x in nonzero),
            'nonzero_residual_entries_per_generator':nonzero,
            'conjugate_triplet_field_assignment':'S1 unbarred / S2 barred',
            'cross_package_rgbeta_tested':False}


def read_json(path: Path) -> dict[str,Any]:
    return NB.read(path)


def validate_live_source(native_root: Path=CG_ROOT) -> dict[str,Any]:
    old = read_json(native_root/'native_cg_bridge.json')
    need(old.get('status')==NB.STATUS and old.get('native_registered_CG_count')==48 and
         old.get('cross_model_projective_CG_comparisons_passed')==36 and
         old.get('matchete_rgbeta_component_basis_certified') is False,
         'preceding 16-model live Matchete CG audit not passing')
    samples={}
    results={}
    digests={}
    for label in MODEL_KEYS:
        p=native_root/'native'/f'{label}.json'
        document=read_json(p)
        q=q_from_live(document,label)
        samples[label]=q
        results[label]=verify_triplet(q)
        digests[label]=source_sha256(p)
    # All three hypercharges within each family must implement the same SU(2)
    # representation up to the projective CG normalization (prior audit proof).
    for f in ('A','D'):
        labels=[m for m in MODEL_KEYS if m.startswith(f+'_')]
        ref=samples[labels[0]]
        for label in labels[1:]:
            need(NB.projectively_equal(tuple(ref),tuple(samples[label])),
                 label+': changed SU2 HH-triplet tensor across hypercharges')
    ad=tuple(sp.Matrix([[number(v) for v in row] for row in mat])
             for mat in results['D_m2']['derived_cartesian_triplet_generators'])
    aa=tuple(sp.Matrix([[number(v) for v in row] for row in mat])
             for mat in results['A_m2']['derived_cartesian_triplet_generators'])
    e_rows={}
    for label in E_KEYS:
        p=native_root/'native'/f'{label}.json'
        e_rows[label]=e_two_triplet_covariance(read_json(p),label,ad,aa)
        digests[label]=source_sha256(p)
    return {'six_models':results,
            'four_E_model_cross_family_covariance_probes':e_rows,
            'E_model_cross_family_covariance_pass_count':sum(
                bool(x['common_matchete_A_D_triplet_generators_su2_covariant'])
                for x in e_rows.values()),
            'ten_native_report_sha256':digests,
            'previous_native_cg_summary_sha256':source_sha256(native_root/'native_cg_bridge.json'),
            'all_six_exact_triplet_witnesses_passed':True}


def probe_value(record: dict[str,Any], shape: list[int], name: str) -> list[sp.Expr]|None:
    need(isinstance(record,dict) and record.get('expected_shape')==shape and
         record.get('shape')==shape, name+': malformed native shape')
    exposed=record.get('exact_numeric_components_exposed')
    need(type(exposed) is bool,name+': native availability is not boolean')
    arr=record.get('components_inputform')
    need(isinstance(arr,list),name+': malformed native components')
    if not exposed:
        need(arr==[],name+': non-numeric native probe cannot supply components')
        return None
    need(len(arr)==int(sp.prod(shape)),name+': incomplete native numeric array')
    return [number(x) for x in arr]


def compare_rgbeta(probe:dict[str,Any], matchete:dict[str,Any])->dict[str,Any]:
    need(probe.get('format')=='t3-workstream-a-rgbeta-triplet-coordinate-probe-v1' and
         probe.get('source')=='installed RGBeta runtime after T3RGBetaBuild[1,3,2,-2]' and
         probe.get('model')=='A_-2' and
         probe.get('integer_index_probe_only') is True and
         probe.get('production_modified') is False and
         all(probe.get(k) is False for k in ('cross_package_component_basis_certified',
                  'barred_field_map_certified','signed_physical_coupling_dictionary_certified')),
         'unexpected RGBeta native probe provenance or certification')
    vals=probe_value(probe['native_triplet_CG_delS2'],[3,2,2],'delS2')
    genvals=probe_value(probe['native_triplet_generators_tGen'],[3,3,3],'tGen')
    epsvals=probe_value(probe['native_fundamental_epsilon'],[2,2],'eps')
    both=vals is not None and genvals is not None
    need(probe.get('numeric_triplet_CG_and_generators_both_exposed') is both,
         'RGBeta numeric availability flag inconsistent')
    if not both:
        return {'status':'RGBetaNativeNumericalTripletCoordinatesUnavailable',
                'native_rgbeta_delS2_component_array_available':vals is not None,
                'native_rgbeta_tGen_component_array_available':genvals is not None,
                'native_rgbeta_epsilon_component_array_available':epsvals is not None,
                'cross_package_triplet_generator_identity_verified':False,
                'native_RGBeta_coordinates_inserted_by_hand':False}
    # A native numerical array is not automatically a valid SU2 representation:
    # assert its exact algebra and compatibility with the native delS2 tensor.
    rg=[sp.Matrix(3,3,genvals[9*k:9*(k+1)]) for k in range(3)]
    for k,r in enumerate(rg):
        need(zero(r.H-r) and zero(r.trace()),f'native RGBeta tGen[{k+1}] not Hermitian/traceless')
    for a,b,c in ((0,1,2),(1,2,0),(2,0,1)):
        need(zero(rg[a]*rg[b]-rg[b]*rg[a]-I*rg[c]),
             'native RGBeta tGen does not have standard oriented su2 algebra')
    need(zero(sum((x*x for x in rg),ZERO3)-2*sp.eye(3)),
         'native RGBeta tGen triplet Casimir != 2')
    # RGBeta delS2 uses the first slot for S2; matrix ordering k,i,j.
    q=sp.Matrix(4,3,lambda ij,k:vals[k*4+ij])
    rgbeta=verify_triplet(q)
    ma=matchete['six_models']['A_m2']
    um=sp.Matrix([[number(x) for x in row] for row in ma['reference_to_matchete_triplet_unitary']])
    ur=sp.Matrix([[number(x) for x in row] for row in rgbeta['reference_to_matchete_triplet_unitary']])
    w=sp.simplify(ur.H*um)
    need(zero(w.H*w-sp.eye(3)),'Matchete/RGBeta triplet basis change not unitary')
    rref=[sp.Matrix([[number(v) for v in row] for row in mat])
          for mat in rgbeta['derived_cartesian_triplet_generators']]
    # Literal integer RGBeta indices may not be physical index coordinates.
    # Accept both expected dual-generator orientations but do not promote
    # to a physical coupling dictionary without live action agreement.
    direct=all(zero(rg[k]-rref[k]) for k in range(3))
    conjugated=all(zero(rg[k]+rref[k].T) for k in range(3))
    need(direct or conjugated,
         'native RGBeta tGen does not match delS2 generator under direct/conjugate orientation')
    return {'status':'NativeRGBetaCGAndGeneratorsInternallyConsistent',
            'native_rgbeta_delS2_component_array_available':True,
            'native_rgbeta_tGen_component_array_available':True,
            'native_rgbeta_epsilon_component_array_available':epsvals is not None,
            'native_rgbeta_generator_orientation':('direct' if direct else 'minus_transpose'),
            'matchete_to_rgbeta_triplet_unitary_candidate':mat_json(w),
            'matchete_rgbeta_unitary_intertwiner_exact':True,
            'native_RGBeta_coordinates_inserted_by_hand':False,
            'physical_signed_coupling_map_certified':False}


def self_test() -> None:
    # Synthetic EXACT 2x2 -> 3 standard SU2 embedding with complex unitary.
    unitary=sp.Matrix([[1,0,0],[0,0,I],[0,1,0]])
    fake=sp.simplify(E*unitary)
    result=verify_triplet(fake)
    need(result['unitarity_su2_commutator_casimir_metric_exact'],
         'synthetic exact triplet proof failed')
    bad=sp.Matrix(fake)
    bad[0,0]=2
    try:
        verify_triplet(bad)
    except ValueError:
        pass
    else:
        raise AssertionError('non-unitary/anisotropic native CG falsely accepted')
    probe={'format':'t3-workstream-a-rgbeta-triplet-coordinate-probe-v1',
           'source':'installed RGBeta runtime after T3RGBetaBuild[1,3,2,-2]',
           'model':'A_-2','integer_index_probe_only':True,'production_modified':False,
           'cross_package_component_basis_certified':False,
           'barred_field_map_certified':False,
           'signed_physical_coupling_dictionary_certified':False,
           'native_triplet_CG_delS2':{'shape':[3,2,2],'expected_shape':[3,2,2],
              'exact_numeric_components_exposed':False,'components_inputform':[]},
           'native_triplet_generators_tGen':{'shape':[3,3,3],'expected_shape':[3,3,3],
              'exact_numeric_components_exposed':False,'components_inputform':[]},
           'native_fundamental_epsilon':{'shape':[2,2],'expected_shape':[2,2],
              'exact_numeric_components_exposed':False,'components_inputform':[]},
           'numeric_triplet_CG_and_generators_both_exposed':False}
    checked=compare_rgbeta(probe,{'six_models':{'A_m2':result}})
    need(checked['status']=='RGBetaNativeNumericalTripletCoordinatesUnavailable',
         'missing native components were incorrectly accepted')
    # Synthetic fully numeric RGBeta reference data with direct generator match.
    q0=sp.Matrix(E)
    cg0=[str(q0[ij,k]) for k in range(3) for ij in range(4)]
    gr0=[-E.H*d*E for d in D]
    probe['native_triplet_CG_delS2']['exact_numeric_components_exposed']=True
    probe['native_triplet_CG_delS2']['components_inputform']=cg0
    probe['native_triplet_generators_tGen']['exact_numeric_components_exposed']=True
    probe['native_triplet_generators_tGen']['components_inputform']=[str(m[i,j]) for m in gr0 for i in range(3) for j in range(3)]
    probe['numeric_triplet_CG_and_generators_both_exposed']=True
    checked=compare_rgbeta(probe,{'six_models':{'A_m2':result}})
    need(checked['status']=='NativeRGBetaCGAndGeneratorsInternallyConsistent',
         'synthetic known-good native component bridge failed')
    probe['native_triplet_generators_tGen']['components_inputform'][0]='7'
    try:
        compare_rgbeta(probe,{'six_models':{'A_m2':result}})
    except ValueError:
        pass
    else:
        raise AssertionError('tampered native RGBeta generator falsely accepted')
    print('PASS: exact SU2 triplet generator, Casimir, conjugate metric and basis intertwiner')
    print('PASS: native RGBeta missing-data gate; synthetic valid/invalid numerical native probes')
    print('NOTE: synthetic tests are not evidence that installed RGBeta exports numerical components')


def run(args:argparse.Namespace) -> None:
    actual=validate_live_source()
    print('Matchete live HH-triplet: 6/6 exact SU2 representation and conjugation-metric witnesses PASS',flush=True)
    out=Path(args.output_dir).resolve()
    out.mkdir(parents=True,exist_ok=True)
    js=out/'rgbeta_native_probe.json'
    log=out/'rgbeta_native_probe.log'
    if args.use_existing:
        need(js.is_file(),'RGBeta native probe does not exist: '+str(js))
    else:
        need(not js.exists() and not log.exists(),
             'Native probe already exists; use --use-existing or choose new --output-dir')
        cmd=[args.wolframscript,'-file',str(SCRIPT),str(js)]
        env=os.environ.copy()
        print('Probing installed RGBeta native delS2 and tGen integer-index components',flush=True)
        try:
            res=subprocess.run(cmd,cwd=ROOT,env=env,text=True,stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,timeout=args.timeout_seconds,check=False)
        except subprocess.TimeoutExpired as exc:
            log.write_text('Timeout: '+str(exc),encoding='utf-8')
            raise ValueError('RGBeta probe timed out: '+str(log)) from exc
        log.write_text(res.stdout,encoding='utf-8')
        need(res.returncode==0,'RGBeta Wolfram probe failed (exit '+str(res.returncode)+'): '+str(log))
    probe=read_json(js)
    checked=compare_rgbeta(probe,actual)
    is_candidate=checked['status']=='NativeRGBetaCGAndGeneratorsInternallyConsistent'
    status=STATUS_CANDIDATE if is_candidate else STATUS_BLOCKED
    all_results={'format':'t3-workstream-a-triplet-coordinate-gate-v1',
                 'status':status, 'matchete_native_baseline_CGs_in_previous_stage':48,
                 'matchete_HH_triplet_models_exactly_checked':6,
                 'E_two_triplet_topology_covariance_probes':4,
                 'matchete_triplet_S2_generators_independently_native_exported':False,
                 'matchete_triplet_witness':actual,
                 'native_RGBeta_coordinate_probe':checked,
                 'native_RGBeta_probe_sha256':source_sha256(js),
                 'rgbeta_native_numeric_triplet_basis_available':is_candidate,
                 'matchete_rgbeta_physical_component_identification_certified':False,
                 'signed_physical_coupling_dictionary_certified':False,
                 'full_UV_physical_normalization_certified':False,
                 'explicit_Z2_registration_certified':False,
                 'matching_performed':False,'production_modified':False}
    (out/'triplet_coordinate_gate.json').write_text(json.dumps(all_results,indent=2)+'\n',encoding='utf-8')
    lines=['# Workstream A — native triplet coordinate gate', '',
           '**Status:** `'+status+'`', '',
           '- Native Matchete HH-triplet CG representation: **6/6**, exact SU(2) tests.',
           '- Previous baseline native Matchete CGs: **48/48**.',
           '- E-family cross-representation covariance probes passed: **'+str(actual['E_model_cross_family_covariance_pass_count'])+'/4** (diagnostic).',
           '- RGBeta `delS2` exact numerical components exposed: **'+str(checked['native_rgbeta_delS2_component_array_available'])+'**.',
           '- RGBeta `tGen` exact numerical components exposed: **'+str(checked['native_rgbeta_tGen_component_array_available'])+'**.',
           '- Installed RGBeta is probed directly; **no numerical component assignments injected**.', '',
           '## Verified/conditional distinction', '',
           'The Matchete triplet representation, SU(2) Casimir, generator algebra, conjugate',
           'bilinear metric and exact unitary map to the declared symmetric-doublet reference',
           'are established from the **live registered** A/D topology Clebsch–Gordan tensors.',
           'This fixes a Matchete-side reference convention, **not an RGBeta physical basis**.',
           'Literal integer indices are only a diagnostic probe of RGBeta. A failure to',
           'resolve numeric components does not invalidate its abstract group algebra.',
           'Even if the diagnostic exports arrays, the candidate is not a certified',
           'cross-package full Lagrangian identity, signed quartic/Yukawa dictionary,',
           'or Matchete matching. The source Stage23 maps remain conditional.', '',
           '| Model | Native CG norm squared | Metric scalar identity | SU2 exact |',
           '|---|---|---|---|']
    for model, row in actual['six_models'].items():
        lines.append(f'| {model} | `{row["native_topology_norm_squared"]}` | '
                     f'{row["bilinear_metric_scalar_identity_in_this_basis"]} | PASS |')
    lines += ['', '## E-family cross-representation diagnostic', '',
              'This separately checks HH S1 S2dag against the A-derived S2dag and',
              'D-derived S1 triplet generators. Any nonzero residual is retained',
              'rather than fixed with an assumed phase or basis rotation.', '',
              '| Model | Exact generator covariance | Nonzero residual counts |',
              '|---|---|---|']
    for model,row in actual['four_E_model_cross_family_covariance_probes'].items():
        lines.append(f'| {model} | {row["common_matchete_A_D_triplet_generators_su2_covariant"]} | '+
                     f'`{row["nonzero_residual_entries_per_generator"]}` |')
    lines += ['', 'No production files, RGE beta functions or EFT outputs were modified.', '']
    (out/'triplet_coordinate_gate.md').write_text('\n'.join(lines),encoding='utf-8')
    print('RESULT:',status,flush=True)
    print('JSON:',out/'triplet_coordinate_gate.json')
    print('REPORT:',out/'triplet_coordinate_gate.md')


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    action=p.add_mutually_exclusive_group(required=True)
    action.add_argument('--self-test',action='store_true')
    action.add_argument('--all',action='store_true')
    p.add_argument('--use-existing',action='store_true',help='Validate a saved probe without running Wolfram')
    p.add_argument('--output-dir',default=str(OUT))
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--timeout-seconds',type=int,default=900)
    args=p.parse_args()
    try:
        if args.self_test:self_test()
        else:run(args)
        return 0
    except (ValueError, FileNotFoundError, TypeError, KeyError, subprocess.SubprocessError) as exc:
        print('FAILED:',exc,file=sys.stderr)
        return 1

if __name__=='__main__':
    sys.exit(main())
