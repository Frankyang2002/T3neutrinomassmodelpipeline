"""Large same-five-mass B/C running + grouped X,Y finite-bubble sensitivity.

Uses exactly UV -> EFT1 -> SM, with one grouped scalar matching step.
Baseline and finite-bubble variants share the SAME saved UV/EFT1 trajectory,
fermion boundary, EFT1 Wilson evolution, and SM starting fields.
Finite scalar bubble is only a RENORMALIZATION-LOG-CHECKED SUBSET.
Not a complete physical finite matching calculation.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import traceback
import numpy as np

from scripts import RunT3BCSameTrajectoryWilson as previous
from scripts import RunT3BCCoupledWilsonSequential as coupled
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import RunT3DirectWilsonSequential as direct
from Numerical.experimental import T3BCCoupledWilsonTensor as tensor
from Numerical.experimental import T3BCGroupedFiniteScalarBubble as bubble
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from Numerical.experimental import T3BCScalarYThresholdProjection as projection
from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context

ROOT=Path(__file__).resolve().parents[1]
MODELS=('B_-1','C_-1')
FORMAT='t3-bc-single-eft-grouped-finite-bubble-megabatch-v1'
STATUS='SameTrajectoryWilsonAndGroupedScalarBubbleComputedPhysicalMatchingPending'
DEFAULT_OUT=ROOT/'output/studies/t3_bc_matching_megabatch_v1'
FIVE=previous.FIVE
AUDIT=previous.AUDIT


def require(ok,reason):
    if not bool(ok):raise ValueError('B/C matching mega batch: '+reason)


def hashfile(path):
    path=Path(path)
    require(path.is_file(),'source file missing '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def array_c5(v):return np.asarray(v,dtype=complex)


def norm_rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))


def endpoint_mass_basis(model,node):
    z=node['muT3GeneralS1S2_gev2']
    require(isinstance(z,dict) and set(z)=={'real','imag'},'mixed scalar parameter format changed')
    bs=MixedDoubletMassBasis.build(model,float(node['mS1Sq_gev2']),
         float(node['mS2Sq_gev2']),complex(z['real'],z['imag']))
    require(not bs.tachyonic and not bs.degenerate,'scalar endpoint tachyonic or degenerate')
    return bs


def exact_sources(model,five,audit_path):
    path=FIVE/(model.replace('-','m')+'_uv_to_eft1_coupled_five_mass_v1.json')
    record=previous.read(path)
    mu=previous.check_five(record,model)
    previous.check_audit(previous.read(audit_path))
    return record,mu,path


def evaluate(model,*,five_dir=FIVE,audit_path=AUDIT,rtol=5e-7,
             secondary_rtol=2e-7,compute_convergence=True):
    require(model in MODELS,'unsupported model')
    require(1e-10<secondary_rtol<rtol<1e-3,'invalid two-tolerance choices')
    five_path=Path(five_dir)/(model.replace('-','m')+'_uv_to_eft1_coupled_five_mass_v1.json')
    five=previous.read(five_path)
    mu=previous.check_five(five,model)
    previous.check_audit(previous.read(audit_path))
    uv=neutrino.decode_uv(five['uv_threshold']['state'])
    states=previous.decode_states(five,model)
    pathmap=coupled.sources(model)
    for name in ('seed','quartic_seed','old_eft1_rgbeta','matching','wilson_rge','stage22','stage29'):
        require(pathmap[name].is_file(),'required source missing '+name)
    source=direct.compile_direct_source(pathmap['matching'],pathmap['wilson_rge'])
    mf,y1,y2,heavyres=direct.heavy_inputs(uv)
    context=build_scalar_only_rge_context(pathmap['seed'],pathmap['quartic_seed'],pathmap['old_eft1_rgbeta'])
    static=tensor.build_static_kinetics(context)
    evaluator=tensor.QuarticEvaluator.build(context.quartic,
                    neutrino.load(pathmap['stage22']),neutrino.load(pathmap['stage29']),model)
    template=coupled.template_at_context(context)
    ids=static['indices']
    tree=tensor.independent_tensor(context.coefficient,tensor.flavor_pairs(mf,y1,y2))
    hard,relocation,grouped=direct.grouped_hard_boundary(source,mf,y1,y2,
                          complex(states[0].baseline.lambdaT3),mu_f=float(mu[0]))
    initial=tree.copy()
    coupled.embed_weinberg(initial,template,ids,grouped)
    start,res0=tensor.weinberg_projection(initial,template,ids)
    require(res0<1e-9 and np.allclose(tensor.flavor_to_matrix(start),grouped,
                                     rtol=1e-10,atol=1e-20),'F hard embedding failed')
    # Every kernel uses a node of the saved five-real-mass EFT1 background.
    quartics=[evaluator.evaluate(st) for st in states]
    kernels=[tensor.kernels(st,static,q) for st,q in zip(states,quartics)]
    result,ode=tensor.integrate_coupled(initial,mu,kernels,rtol=rtol)
    Cend=result[...,-1]
    p0,r0=tensor.weinberg_projection(Cend,template,ids)
    c5=neutrino.validated_c5(tensor.flavor_to_matrix(p0))
    basis=endpoint_mass_basis(model,five['eft1_final'])
    T=projection.mixing_real_matrix(basis.U)
    maperr=projection.check_field_identity(basis.U,T)
    gaugeerr=projection.check_gauge_covariance(T)
    mass_C=projection.rotate_wilson(Cend,T)
    # Rotate all 4 quartic legs. Coefficients are native real-scalar V=Q/4!.
    mass_q=np.einsum('abcd,ap,bq,cr,ds->pqrs',quartics[-1],T,T,T,T,optimize=True)
    require(np.linalg.norm(mass_q.imag)<1e-8*max(np.linalg.norm(mass_q),1e-25),
            'rotated quartic not a real scalar potential')
    mass_q=np.real(mass_q)
    originalHH=mass_C[:,:,:4,:4,:]
    require(np.allclose(originalHH,Cend[:,:,:4,:4,:],rtol=1e-12,atol=1e-25),
            'Higgs block changed under massive field mixing')
    # Evaluate common-threshold bubble at the saved EFT1 endpoint, not the
    # historical separately-decoupled Y/X EFT2 chain.
    m2=tuple(float(a) for a in basis.eigenvalues_sq)
    mu_s=float(mu[-1])
    delta,diag=bubble.group_bubble(mass_C,mass_q,m2,mu_s)
    logcheck=bubble.rge_log_check(mass_C,mass_q,m2,mu_s)
    delta5,res_delta=bubble.project_delta(delta,template,ids)
    changed=neutrino.validated_c5(c5+delta5)
    # All SM evolutions are re-run with the SAME EFT1 endpoint SM fields.
    baseline=neutrino.low_energy(model,c5,states[-1].baseline)
    updated=neutrino.low_energy(model,changed,states[-1].baseline)
    frozen_scale={}
    for factor in (0.5,1.,2.):
        shifted,d=bubble.group_bubble(mass_C,mass_q,m2,mu_s*factor)
        shift5,_=bubble.project_delta(shifted,template,ids)
        cm=neutrino.validated_c5(c5+shift5)
        low=neutrino.low_energy(model,cm,states[-1].baseline)
        frozen_scale[str(factor)]={'relative_C5_shift':norm_rel(cm,c5),
             'mass_eV_descending':low['takagi_masses_ev_descending'],
             'matching_scale_gev_for_loop_log_only':mu_s*factor}
    # Local RGE sensitivity of the actual tensor to native quartic sources.
    qb=evaluator.evaluate(states[-1],include_extra=False)
    qa=evaluator.evaluate(states[-1],include_extra=True)
    fullbeta=tensor.rhs_master(Cend,kernels[-1])
    deltaq=np.einsum('abcd,ijcdf->ijabf',qa-qb,Cend,optimize=True)/tensor.LOOP
    localquartic={'full_beta_norm':float(np.linalg.norm(fullbeta/tensor.LOOP)),
       'extra_quartic_beta_norm':float(np.linalg.norm(deltaq)),
       'extra_fraction_of_full_beta':float(np.linalg.norm(deltaq)/max(np.linalg.norm(fullbeta/tensor.LOOP),1e-30))}
    convergence=None
    if compute_convergence:
        fine,ode2=tensor.integrate_coupled(initial,mu,kernels,rtol=secondary_rtol)
        pp,pr=tensor.weinberg_projection(fine[...,-1],template,ids)
        refine=neutrino.validated_c5(tensor.flavor_to_matrix(pp))
        convergence={'wilson_rtol_primary':rtol,'wilson_rtol_secondary':secondary_rtol,
          'relative_C5_difference':norm_rel(refine,c5),
          'second_projection_residual':float(pr),'second_solver_nfev':int(ode2.nfev)}
    provenance={'five_mass_trajectory':hashfile(five_path),'batch1_audit':hashfile(audit_path),
      'same_trajectory_runner':hashfile(Path(previous.__file__)),
      'mega_runner':hashfile(Path(__file__)),
      'bubble_kernel':hashfile(Path(bubble.__file__)),
      **{k:hashfile(v) for k,v in pathmap.items()}}
    return {'format':FORMAT,'status':STATUS,'model':model,
      'one_intermediate_EFT_only':True,'same_five_mass_eft1_scalar_nodes':len(mu),
      'running_eft1_interval_gev':[float(mu[0]),mu_s],
      'grouped_scalar_matching_scale_gev':mu_s,
      'heavy_masses_at_fermion_threshold_gev':[float(x) for x in mf],
      'scalar_eigenmasses_at_matching_gev':[math.sqrt(x) for x in m2],
      'heavy_basis_residual':float(heavyres),'mass_basis_field_identity':maperr,
      'mass_basis_gauge_covariance':gaugeerr,
      'fermion_hard_boundary_frobenius_norm':float(np.linalg.norm(grouped)),
      'weinberg_projection_residual':float(r0),
      'bubble_weinberg_projection_residual':float(res_delta),
      'baseline_C5_at_muS':neutrino.matrix_payload(c5),
      'bubble_delta_C5_at_muS':neutrino.matrix_payload(delta5),
      'bubble_adjusted_C5_at_muS':neutrino.matrix_payload(changed),
      'baseline_masses_ev_descending':baseline['takagi_masses_ev_descending'],
      'bubble_adjusted_masses_ev_descending':updated['takagi_masses_ev_descending'],
      'relative_bubble_delta_C5':float(np.linalg.norm(delta5)/max(np.linalg.norm(c5),1e-30)),
      'bubble_sector_diagnostics':diag,'bubble_RGE_log_consistency':logcheck,
      'matching_log_scale_variation_frozen_background':frozen_scale,
      'local_native_extra_quartic_beta_sensitivity':localquartic,
      'primary_eft1_wilson_nfev':int(ode.nfev),
      'convergence':convergence,'input_sha256':provenance,
      'full_finite_F_matching_calculated':False,
      'full_finite_scalar_matching_calculated':False,
      'independent_grouped_scalar_bubble_finite_normalization_certified':False,
      'full_three_family_wilson_anomalous_dimension':False,
      'physical_C5_certified':False,'physical_neutrino_spectrum_certified':False,
      'notes':('Only the one-quartic scalar-pair bubble finite sector is included. '
        'Matching scale scan varies the bubble logarithm at FROZEN endpoint couplings; '
        'it is not complete physical matching-scale variation. No EFT2.')}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--five-mass-dir',type=Path,default=FIVE)
    ap.add_argument('--audit',type=Path,default=AUDIT)
    ap.add_argument('--output-dir',type=Path,default=DEFAULT_OUT)
    ap.add_argument('--rtol',type=float,default=5e-7)
    ap.add_argument('--secondary-rtol',type=float,default=2e-7)
    ap.add_argument('--skip-convergence',action='store_true')
    ns=ap.parse_args()
    ns.output_dir.mkdir(parents=True,exist_ok=True)
    rows=[];failures=[]
    for model in MODELS:
        try:
            record=evaluate(model,five_dir=ns.five_mass_dir,audit_path=ns.audit,
                    rtol=ns.rtol,secondary_rtol=ns.secondary_rtol,
                    compute_convergence=not ns.skip_convergence)
            target=ns.output_dir/(model.replace('-','m')+'_grouped_finite_bubble.json')
            require(not target.exists(),'refusing to overwrite '+str(target))
            temp=target.with_suffix('.json.tmp')
            temp.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            temp.replace(target)
            rows.append({'model':model,'baseline_masses_ev_descending':record['baseline_masses_ev_descending'],
               'bubble_adjusted_masses_ev_descending':record['bubble_adjusted_masses_ev_descending'],
               'relative_bubble_delta_C5':record['relative_bubble_delta_C5'],
               'relative_log_check_residual':record['bubble_RGE_log_consistency']['relative_log_derivative_residual'],
               'wilson_relative_convergence':None if record['convergence'] is None else record['convergence']['relative_C5_difference']})
            print(model,'same-trajectory baseline',record['baseline_masses_ev_descending'],
                  'bubble candidate',record['bubble_adjusted_masses_ev_descending'],
                  'relative delta C5',record['relative_bubble_delta_C5'],flush=True)
        except Exception as e:
            traceback.print_exc();failures.append({'model':model,'exception':type(e).__name__,'error':str(e)})
    summary={'format':FORMAT,'status':STATUS if not failures and len(rows)==2 else 'BlockedReviewRequired',
       'models':rows,'failures':failures,'one_intermediate_EFT_only':True,
       'computed_finite_matching_scope':'single quartic scalar-pair bubble only',
       'physical_C5_certified':False}
    path=ns.output_dir/'summary.json'
    path.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print('MEGA',summary['status'],'physical C5=NOT CERTIFIED',path,flush=True)
    return 0 if not failures else 1

if __name__=='__main__':raise SystemExit(main())
