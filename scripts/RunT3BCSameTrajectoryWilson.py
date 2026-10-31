"""Same-five-mass-trajectory B/C UV -> EFT1 Wilson -> grouped SM C5 candidate.

One EFT1 only; no Y->X EFT. All EFT1 background nodes are taken from the
saved Workstream-B five-mass trajectory. No historical continuous_v1 scalars.
This is a conditional operator evolution, NOT renormalized finite matching.
"""
from __future__ import annotations
import argparse
from dataclasses import replace
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
MODELS=('B_-1','C_-1')
FORMAT='t3-bc-same-five-mass-eft1-wilson-grouped-scalar-v1'
STATUS='SameTrajectoryEFT1WilsonIntegratedFiniteMatchingPending'
FIVE=ROOT/'output/studies/t3_full_stack_overnight_v1/bc_full_mass'
AUDIT=ROOT/'output/studies/t3_bc_fermion_threshold_batch_v1/summary.json'
OUT=ROOT/'output/studies/t3_bc_same_trajectory_wilson_v1'

def require(test,msg):
    if not bool(test):raise ValueError('Same-trajectory B/C Wilson: '+msg)

def read(p):
    require(p.is_file(),'missing input '+str(p))
    o=json.loads(p.read_text(encoding='utf-8-sig'))
    require(isinstance(o,dict),'invalid JSON object '+str(p))
    return o

def digest(p):
    require(p.is_file(),'missing provenance input '+str(p))
    return sha256(p.read_bytes()).hexdigest()

def check_five(five,model):
    from scripts import RunT3BCCoupledUVToEFT1 as work
    require(five.get('format')==work.FORMAT and five.get('model')==model and
            five.get('status')=='CandidateCoupledUVAndEFT1FiveRealMassesIntegrated' and
            five.get('solver_success') is True,'invalid five-mass trajectory')
    require(five.get('physical_matchete_rgbeta_mapping_certified') is False and
            five.get('physical_neutrino_mass_prediction') is False and
            five.get('scalar_mass_eigenstate_threshold_matching_performed') is False,
            'five-mass physical status changed; re-audit required')
    require(five.get('uv_range_gev')==[1e7,1e5] and
            five.get('eft1_range_gev')==[1e5,1e3], 'unexpected UV/EFT1 scales')
    require(five.get('threshold',{}).get('wilson_coefficient_matching_performed') is False,
            'original five-mass fermion matching status changed')
    for raw,expected in five['input_sha256'].items():
        p=Path(raw)
        if not p.is_absolute():p=ROOT/p
        require(digest(p)==expected,'five-mass upstream source changed '+str(p))
    nodes=five['eft1_trajectory']
    require(len(nodes)>=9,'too few EFT1 nodes for Wilson solver')
    mu=np.array([n['mu_gev'] for n in nodes],float)
    require(np.all(np.isfinite(mu)) and np.all(mu>0) and np.all(np.diff(mu)<0),
            'EFT1 scales invalid')
    require(abs(mu[0]/1e5-1)<1e-10 and abs(mu[-1]/1e3-1)<1e-10,
            'EFT1 end points changed')
    for key in ('mHSq_gev2','mS1Sq_gev2','mS2Sq_gev2','muT3GeneralS1S2_gev2'):
        require(five['uv_threshold'][key]==five['eft1_threshold'][key],
                'F-threshold scalar mass discontinuity '+key)
    return mu

def check_audit(audit):
    require(audit.get('status')=='ConditionalUVToEFT1SourceBoundaryAudited',
            'Batch1 source audit not passed')
    require(audit.get('physical_C5_certified') is not True,'unexpected physical certificate')

def decode_states(five,model):
    from scripts import IntermediateEFTRunning as s34
    from Numerical.experimental.GeneralT3ScalarRGState import GeneralT3ScalarRGState,load_registry
    from scripts import RunT3FullScalarNeutrino as neutrino
    registry=load_registry(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    states=[]
    for node in five['eft1_trajectory']:
        raw=dict(node['state'],mu_gev=float(s34.START))
        base=replace(s34.decode_intermediate_state(raw),mu_gev=float(node['mu_gev'])).validated()
        extra={key:neutrino.complex_pair(val,key) for key,val in node['additional_quartics'].items()}
        states.append(GeneralT3ScalarRGState(base,extra).validated(registry))
    return states

def candidate_summary(model,mu,grouped,c5,residual,ode,states,heavy_res,
                      mass,low,source_hashes,rtol):
    from scripts import RunT3FullScalarNeutrino as neutrino
    require(np.isfinite(c5).all() and np.linalg.norm(c5)>0,'empty/nonfinite C5')
    require(np.linalg.norm(c5-c5.T)/np.linalg.norm(c5)<1e-9,'C5 asymmetric')
    require(math.isfinite(residual) and residual<1e-4,'large Weinberg projection residual')
    return dict(format=FORMAT,status=STATUS,model=model,
        source_trajectory='WorkstreamBFiveRealMassesWithRunningUVHeavyFermions',
        source_provenance_sha256=source_hashes,
        UV_to_EFT1_to_SM_one_intermediate_EFT=True,
        old_continuous_v1_scalar_trajectory_used=False,
        uv_fermion_threshold_gev=float(mu[0]),grouped_scalar_threshold_gev=float(mu[-1]),
        source_eft1_nodes=len(mu),eft1_wilson_rtol=float(rtol),
        eft1_wilson_nfev=int(ode.nfev),eft1_wilson_solver_status=ode.message,
        heavy_basis_residual=float(heavy_res),
        fermion_grouped_hard_C5_candidate=neutrino.matrix_payload(grouped),
        grouped_scalar_threshold_C5_candidate=neutrino.matrix_payload(c5),
        scalar_threshold_weinberg_projection_residual=float(residual),
        masses_ev_descending=[float(x) for x in mass],
        low_energy_sm_proxy=low,
        source_matching_boundary='historical grouped hard+log convention, not independently renormalized',
        seven_term_wilson_master_used=True,
        finite_F_matching_counterterms_computed=False,
        finite_scalar_threshold_matching_computed=False,
        full_three_family_wilson_anomalous_dimension=False,
        physical_C5_certified=False,physical_neutrino_spectrum_certified=False)

def compute(model, five_path, audit_path, *,rtol=5e-7):
    from scripts import RunT3FullScalarNeutrino as neutrino
    from scripts import RunT3DirectWilsonSequential as direct
    from scripts import RunT3BCExtraWilsonSequential as selected
    from scripts import RunT3BCCoupledWilsonSequential as coupled
    from Numerical.experimental import T3BCCoupledWilsonTensor as tensor
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    require(model in MODELS,'unsupported model')
    require(1e-10<rtol<1e-3,'invalid Wilson tolerance')
    five=read(five_path);mu=check_five(five,model)
    audit=read(audit_path);check_audit(audit)
    uv=neutrino.decode_uv(five['uv_threshold']['state'])
    states=decode_states(five,model)
    # Verify exact saved scalar nodes, not an unrelated historical background.
    for state,node in zip(states,five['eft1_trajectory']):
        for field in ('mS1Sq','mS2Sq'):
            target=float(node[field+'_gev2'])
            require(abs(float(getattr(state.baseline,field))-target)<1e-8*max(abs(target),1.),
                    'decoded scalar state disagrees with five-mass trajectory: '+field)
    paths=coupled.sources(model)
    # Existing native LLSS tree, grouped F hard LLHH and Wilson master seeds.
    for label in ('seed','quartic_seed','old_eft1_rgbeta','matching','wilson_rge','stage22','stage29'):
        require(label in paths,'missing project source mapping '+label)
        require(paths[label].is_file(),'missing source '+label+' '+str(paths[label]))
    source=direct.compile_direct_source(paths['matching'],paths['wilson_rge'])
    mf,y1,y2,heavy_res=direct.heavy_inputs(uv)
    context=build_scalar_only_rge_context(paths['seed'],paths['quartic_seed'],paths['old_eft1_rgbeta'])
    static=tensor.build_static_kinetics(context)
    quartic=tensor.QuarticEvaluator.build(context.quartic,
              neutrino.load(paths['stage22']),neutrino.load(paths['stage29']),model)
    template=coupled.template_at_context(context)
    ids=static['indices']
    tree=tensor.independent_tensor(context.coefficient,tensor.flavor_pairs(mf,y1,y2))
    initial=tree.copy()
    hard,relocated,grouped=direct.grouped_hard_boundary(
        source,mf,y1,y2,complex(states[0].baseline.lambdaT3),mu_f=float(mu[0]))
    coupled.embed_weinberg(initial,template,ids,grouped)
    boundary_projected,boundary_res=tensor.weinberg_projection(initial,template,ids)
    require(boundary_res<1e-9 and np.allclose(tensor.flavor_to_matrix(boundary_projected),grouped,
        rtol=1e-10,atol=1e-20),'fermion boundary embedding failed')
    kernels=[tensor.kernels(state,static,quartic.evaluate(state)) for state in states]
    evolved,ode=tensor.integrate_coupled(initial,mu,kernels,rtol=rtol)
    projected,residual=tensor.weinberg_projection(evolved[...,-1],template,ids)
    c5=neutrino.validated_c5(tensor.flavor_to_matrix(projected))
    # Existing SM one-loop Weinberg evolution is started at the SAME EFT1 endpoint.
    low=neutrino.low_energy(model,c5,states[-1].baseline)
    mass=np.asarray(low['takagi_masses_ev_descending'],float)
    require(mass.shape==(3,) and np.isfinite(mass).all() and np.all(mass>=0),
            'SM mass spectrum invalid')
    hashes={'five_mass':digest(five_path),'fermion_source_audit':digest(audit_path),
            **{key:digest(path) for key,path in paths.items()}}
    return candidate_summary(model,mu,grouped,c5,residual,ode,states,heavy_res,
                             mass,low,hashes,rtol)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--five-mass-dir',type=Path,default=FIVE)
    ap.add_argument('--audit',type=Path,default=AUDIT)
    ap.add_argument('--output-dir',type=Path,default=OUT)
    ap.add_argument('--rtol',type=float,default=5e-7)
    ns=ap.parse_args()
    ns.output_dir.mkdir(parents=True,exist_ok=True)
    rows=[];failures=[]
    for model in MODELS:
        tag=model.replace('-','m')
        try:
            path=ns.five_mass_dir/f'{tag}_uv_to_eft1_coupled_five_mass_v1.json'
            result=compute(model,path,ns.audit,rtol=ns.rtol)
            target=ns.output_dir/f'{tag}_same_trajectory_wilson.json'
            target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            rows.append({'model':model,'masses_ev_descending':result['masses_ev_descending'],
               'weinberg_projection_residual':result['scalar_threshold_weinberg_projection_residual'],
               'nfev':result['eft1_wilson_nfev']})
            print(model,'SAME five-mass EFT1 Wilson candidate (eV)',result['masses_ev_descending'],
                  'projection',result['scalar_threshold_weinberg_projection_residual'],flush=True)
        except (OSError,ValueError,KeyError,TypeError,RuntimeError,ArithmeticError,AssertionError) as err:
            traceback.print_exc();failures.append({'model':model,'error':str(err)})
    status=STATUS if len(rows)==2 and not failures else 'BlockedNoPhysicalC5'
    (ns.output_dir/'summary.json').write_text(json.dumps({'format':FORMAT,'status':status,
       'models':rows,'failures':failures,'one_intermediate_EFT_only':True,
       'same_trajectory_Wilson_completed':status==STATUS,
       'finite_threshold_matching_computed':False,'physical_C5_certified':False},
       indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('SUMMARY',status,'physical C5=NOT CERTIFIED',ns.output_dir/'summary.json')
    return 0 if status==STATUS else 1

if __name__=='__main__':raise SystemExit(main())
