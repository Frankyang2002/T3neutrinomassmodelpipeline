#!/usr/bin/env python3
"""One genuinely complex direct-UV T3 one-loop convention check for A_+0,B_-1,C_-1.

Inputs are actual saved installed-Matchete native action/CG witnesses and local
benchmark Yukawas. Not the provisional inputs or core.py. The Feynman-parameter
integral is cross-checked against a separate analytic B0 divided difference.

This does NOT identify a numerical Matchete GammaCC matrix/phase; it tests the
conjugation pattern in a DECLARED two-component Weyl convention. Outputs remain
physical_source_certified=False, irrespective of numerical agreement.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
import mpmath as mp
import numpy as np

REL_BASE = Path('output/general_observable/physical_normalisation')
MODEL_FILES = {
 'A_p0':'Numerical/benchmarks/fitted_reference/T3_dS1_1_dS2_3_dF_2_alpha_p0.json',
 'B_m1':'Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_1_alpha_m1.json',
 'C_m1':'Numerical/benchmarks/fitted_reference/T3_dS1_2_dS2_2_dF_3_alpha_m1.json',
}
WITNESS_FORMATS={
 'A_p0':'t3-b1-a0-live-native-scalar-quadratic-v1',
 'B_m1':'t3-b1-bm1-live-native-scalar-quadratic-v2',
 'C_m1':'t3-b1-cm1-live-native-scalar-quadratic-v1',
}
SPINOR_FORMAT='t3-b1-native-spinor-convention-v1'
SPINOR_STATUS='LiveNativeB1SpinorOperatorsAudited_AbsoluteGammaCCPhaseUnresolved'


def require(cond, message):
 if not cond: raise ValueError(message)


def read_json(p: Path) -> dict:
 require(p.is_file(),f'Missing input: {p}')
 result=json.loads(p.read_text(encoding='utf-8-sig'))
 require(isinstance(result,dict),f'Non-object JSON: {p}')
 return result


def fingerprint(p: Path)->str:
 return sha256(p.read_bytes()).hexdigest()


def check_sources(root:Path, row:dict, label:str)->None:
 source=row['source_sha256']
 require(isinstance(source,dict) and len(source)>=5, f'{label}: source hashes absent')
 for rel,hash_expected in source.items():
  p=root/rel
  require(p.is_file() and fingerprint(p).lower()==hash_expected.lower(),
          f'{label}: source changed or absent after native witness: {rel}')


def native_array(row:dict, name:str, shape:tuple[int,...])->np.ndarray:
 obj=row['native_tensors'][name]
 require(tuple(obj['shape'])==shape, f'{name}: native shape mismatch')
 real=np.asarray(obj['real_flat'],float)
 imag=np.asarray(obj['imag_flat'],float)
 require(real.size==math.prod(shape) and imag.size==real.size,f'{name}: invalid payload')
 z=(real+1j*imag).reshape(shape)
 require(np.isfinite(z).all(),f'{name}: nonfinite tensor')
 return z


def check_witness(root:Path, label:str, witness:dict)->tuple[complex,complex,complex]:
 require(witness.get('format')==WITNESS_FORMATS[label],f'{label}: wrong native witness format')
 require(witness.get('model')=={'A_p0':'A_+0','B_m1':'B_-1','C_m1':'C_-1'}[label],
         f'{label}: wrong model')
 require(witness.get('physical_source_certified') is False,f'{label}: overstated witness')
 checks=witness.get('native_component_checks' if label=='A_p0' else 'native_checks',{})
 require(len(checks)>6 and all(x is True for x in checks.values()),
         f'{label}: failed native component witness')
 check_sources(root,witness,label)
 if label=='A_p0':
  y1=native_array(witness,'T3Y1CG',(2,2))
  y2=native_array(witness,'T3Y2CG',(2,2,3))
  mix=native_array(witness,'T3MixCG',(2,2,3))
  # u_- fixed by the independent T3 generator witness; no provisional CGs.
  uminus=np.array([1j,0,1],complex)/np.sqrt(2)
  c1=complex(y1[0,1]);c2=complex(np.dot(y2[0,1,:],uminus.conj()))
  cm=complex(np.dot(mix[1,1,:],uminus.conj()))
  require(abs(c1-1)<2e-13 and abs(c2+1j*np.sqrt(2)*3**(-.25))<2e-13
          and abs(cm-1j*np.sqrt(2)*3**(-.25))<2e-13,
          'A_+0: native neutral CGs conflict with chosen weight convention')
 elif label=='B_m1':
  y1=native_array(witness,'T3Y1CG',(2,2));y2=native_array(witness,'T3Y2CG',(2,2))
  mix=native_array(witness,'T3MixCG',(2,2,2,2))
  c1=complex(y1[0,0]);c2=complex(y2[0,1]);cm=complex(mix[1,1,0,1])
  require(abs(c1-1)<2e-13 and abs(c2-1)<2e-13 and abs(cm+2/np.sqrt(3))<2e-13,
          'B_-1: native neutral CGs changed')
 else:
  y1=native_array(witness,'T3Y1CG',(2,3,2))
  y2=native_array(witness,'T3Y2CG',(2,3,2))
  mix=native_array(witness,'T3MixCG',(2,2,2,2))
  c1=complex(y1[0,1,0]);c2=complex(y2[0,1,1]);cm=complex(mix[1,1,0,1])
  require(abs(c1-3**(-.25))<2e-13 and abs(c2+3**(-.25))<2e-13
          and abs(cm+2/np.sqrt(3))<2e-13,'C_-1: native neutral CGs changed')
  require(np.count_nonzero(abs(y1[0,:,0])>1e-12)==1 and
          np.count_nonzero(abs(y2[0,:,1])>1e-12)==1,
          'C_-1: fermion multiplicity incorrectly inferred')
 return c1,c2,cm


def ordinary_matrix(obj)->np.ndarray:
 if isinstance(obj,dict):
  require(set(obj)=={'real','imag'},'Bad complex source encoding')
  out=np.asarray(obj['real'],float)+1j*np.asarray(obj['imag'],float)
 else: out=np.asarray(obj,complex)
 require(np.isfinite(out).all(),'Nonfinite source')
 return out


def scalar_source(obj)->complex:
 z=ordinary_matrix(obj)
 require(z.ndim==0,'Expected scalar source coupling')
 return complex(z.item())


def mpz(z:complex)->mp.mpc:
 z=complex(z)
 return mp.mpc(repr(float(z.real)),repr(float(z.imag)))


def mp_feynman(M2, lo, hi):
 """Scalar divided difference from x integral, not from production B0."""
 d=hi-lo
 return -mp.quad(lambda x: mp.log1p((1-x)*d/(x*M2+(1-x)*lo))/d,[0,1])


def mp_analytic(M2, lo, hi):
 """Independent exact zero-momentum B0 finite difference at mu^2=M^2."""
 def B0(s): return 1+s*mp.log(s/M2)/(M2-s)
 return (B0(hi)-B0(lo))/(hi-lo)


def complex_data(bench:dict):
 o=bench['base_state']['ordinary']
 y1=ordinary_matrix(o['y1']); y2=ordinary_matrix(o['y2'])
 mf=ordinary_matrix(o['MF'])
 require(y1.shape==y2.shape==mf.shape==(3,3),'Expected 3x3 source flavour')
 require(np.linalg.norm(mf-np.diag(np.diag(mf)))<1e-12 and
         np.max(abs(np.diag(mf).imag))<1e-12 and (np.diag(mf).real>0).all(),
         'Only source benchmark real diagonal masses supported')
 for k in ('lambdaH1Adj','lambdaH2Adj','lambda12Adj',
           'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
           'lambdaHHdagS1barS2barCross'):
  if k in o: require(abs(scalar_source(o[k]))<1e-15,f'Extra quadratic unsupported: {k}')
 return y1,y2,np.diag(mf).real,scalar_source(o['lambdaT3'])


def evaluate(bench:dict, cg:tuple[complex,complex,complex], y1:np.ndarray,
             y2:np.ndarray, lam:complex, portal='benchmark-identity'):
 mp.mp.dps=65
 o=bench['base_state']['ordinary']; v=mp.mpf(str(bench['numerical']['vev_gev']))
 m1=mp.mpf(str(o['mS1Sq']));m2=mp.mpf(str(o['mS2Sq']))
 if portal=='benchmark-identity':
  m1-=mp.mpf(str(o['lambdaH1']))*v*v/2
  m2-=mp.mpf(str(o['lambdaH2']))*v*v/2
 else: require(portal=='zero','Unexpected portal assignment')
 c1,c2,cm=cg
 A21=-mpz(lam)*v*v*mpz(cm)/2
 spread=mp.sqrt(((m1-m2)/2)**2+abs(A21)**2)
 lo=(m1+m2)/2-spread;hi=(m1+m2)/2+spread
 require(lo>0 and hi>lo,'Nonpositive/degenerate scalar spectrum')
 masses=ordinary_matrix(o['MF']).diagonal().real
 g1=np.empty((3,3),dtype=complex);g2=np.empty_like(g1)
 for a in range(3):
  for r in range(3):
   g1[a,r]=-np.conj(y1[a,r]*c1)
   g2[a,r]=-np.conj(y2[a,r]*c2)
 def mass_matrix(Ds):
  result=mp.matrix(3)
  for r in range(3):
   pref=mp.mpf(str(masses[r]))*A21*Ds[r]/(16*mp.pi**2)
   for a in range(3):
    for b in range(3):
     result[a,b]+=pref*(mpz(g1[a,r])*mpz(g2[b,r])+mpz(g2[a,r])*mpz(g1[b,r]))
  return np.array([[complex(result[a,b])*1e9 for b in range(3)] for a in range(3)])
 feynman=[];analytic=[]
 for M in masses:
  M2=mp.mpf(str(M))**2
  feynman.append(mp_feynman(M2,lo,hi))
  analytic.append(mp_analytic(M2,lo,hi))
 x=mass_matrix(feynman);y=mass_matrix(analytic)
 residue=float(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-40))
 require(residue<1e-12, f'Independent Feynman integral versus analytic B0 mismatch: {residue}')
 require(float(np.linalg.norm(x-x.T))<1e-12*max(float(np.linalg.norm(x)),1e-40),
         'Majorana matrix not symmetric')
 return x, {'method_residual':residue,'scalar_squares_gev2':[float(lo),float(hi)],
            'native_mixing_A21_gev2':encode_scalar(A21),
            'takagi_ev_descending':np.linalg.svd(x,compute_uv=False).tolist()}


def encode_scalar(z):
 return {'real':float(z.real),'imag':float(z.imag)}


def encode_matrix(z:np.ndarray):
 return {'real':z.real.tolist(),'imag':z.imag.tolist()}


def cp_invariant(m:np.ndarray)->float:
 H=m@m.conj().T
 numer=np.imag(H[0,1]*H[1,2]*H[2,0]);denom=(np.trace(H).real/3)**3
 return float(numer/denom)


def relative(a:np.ndarray,b:np.ndarray)->float:
 return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-40))


def check_spinor_record(root:Path, row:dict):
 require(row.get('format')==SPINOR_FORMAT and row.get('status')==SPINOR_STATUS,
         'Live Matchete spinor witness missing/incorrect')
 require(row.get('physical_source_certified') is False,'Improper native certification')
 check_sources(root,row,'shared-spinor')
 cases=row['models']
 require(set(cases)=={'A_p0','B_m1','C_m1'},'Missing native model in spinor audit')
 for label,model in cases.items():
  require(model.get('native_action_passed') is True and
          len(model['checks'])>=6 and all(v is True for v in model['checks'].values()),
          f'{label}: native operator check failed')
  require(model['mass']['coefficients']==(['-1'] if label=='A_p0' else ['-1/2']),
          f'{label}: unexpected native mass normalisation')
  for k in ('y1_original','y2_original','y1_hc','y2_hc'):
   require(model[k]['coefficients']==['1'],f'{label}: unexpected Yukawa factor {k}')
 return {'gammaCC_concrete_phase_measured':False,
         'native_relative_spinor_sign_from_concrete_components_measured':False,
         'action_operator_signs_and_HC_checked':True}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--repo-root',type=Path,default=Path('.'))
 parser.add_argument('--spinor-native',type=Path,required=True)
 parser.add_argument('--a0-native',type=Path,required=True)
 parser.add_argument('--bm1-native',type=Path,required=True)
 parser.add_argument('--cm1-native',type=Path,required=True)
 parser.add_argument('--output-dir',type=Path,required=True)
 args=parser.parse_args()
 root=args.repo_root.resolve();out=args.output_dir.resolve()
 require(str(out).startswith(str((root/REL_BASE).resolve())+str(Path('/')) ) or
         out.parent==(root/REL_BASE).resolve(),
         'Output directory must be inside output/general_observable/physical_normalisation')
 require(not out.exists(),f'Refusing to overwrite existing output: {out}')
 inputs={'spinor':args.spinor_native,'A_p0':args.a0_native,'B_m1':args.bm1_native,'C_m1':args.cm1_native}
 for name,path in inputs.items():
  require(path.resolve().is_relative_to(root),f'Witness outside local project: {name}')
 native={name:read_json(path) for name,path in inputs.items()}
 shared=check_spinor_record(root,native['spinor'])
 results={}
 for label in ('A_p0','B_m1','C_m1'):
  cg=check_witness(root,label,native[label])
  bp=root/MODEL_FILES[label]; bench=read_json(bp)
  y1,y2,_,lam=complex_data(bench)
  base, bmeta=evaluate(bench,cg,y1,y2,lam)
  # A genuinely complex, non-factorisable benchmark deformation. These are
  # explicitly CHOSEN physical inputs, not run or fitted parameters.
  yy1=y1.copy();yy2=y2.copy()
  yy1[0,1]*=np.exp(0.37j); yy1[1,2]*=np.exp(-0.19j)
  yy2[2,0]*=np.exp(0.41j); yy2[0,2]*=np.exp(-0.23j)
  ll=lam*np.exp(0.29j)
  complex_result,cmeta=evaluate(bench,cg,yy1,yy2,ll)
  # A global source-coupling phase test must use conjugates of BOTH Yukawas
  # and the UNconjugated scalar mixing lambda, in the stated Weyl convention.
  th1,th2,thlam=0.21,-0.17,0.43
  rotated,_=evaluate(bench,cg,yy1*np.exp(1j*th1),yy2*np.exp(1j*th2),
                     ll*np.exp(1j*thlam))
  prediction=complex_result*np.exp(1j*(thlam-th1-th2))
  covres=relative(rotated,prediction)
  require(covres<2e-12,f'{label}: incorrect complex-source spurion phases: {covres}')
  wrong_y1,_=evaluate(bench,cg,yy1.conj(),yy2,ll)
  wrong_ys_residual=relative(complex_result,wrong_y1)
  require(wrong_ys_residual>1e-4,
          f'{label}: deliberately wrong Yukawa-conjugation control did not distinguish cases')
  results[label]={
   'native_CG':{'y1':encode_scalar(cg[0]),'y2':encode_scalar(cg[1]),'mix':encode_scalar(cg[2])},
   'unmodified_benchmark_identity_portal':{'matrix_ev':encode_matrix(base),**bmeta},
   'complex_nondiagonal_phases_identity_portal':{'matrix_ev':encode_matrix(complex_result),**cmeta,
      'cp_odd_matrix_invariant':cp_invariant(complex_result)},
   'complex_phase_input':{'y1_entrywise':{'(0,1)':0.37,'(1,2)':-0.19},
      'y2_entrywise':{'(2,0)':0.41,'(0,2)':-0.23},'lambdaT3_global':0.29},
   'global_source_phase_transport':{'y1':th1,'y2':th2,'lambdaT3':thlam,
      'expected_matrix_phase_radians':thlam-th1-th2,'relative_residual':covres},
   'wrong_y1_conjugation_negative_control_residual':wrong_ys_residual,
   'native_fermion_masses':{'single_symbol_MF':True,
      'benchmark_diagonal_MF_generalisation_is_conditional':True},
   'benchmark_identity_portal_assignment_is_conditional':True,
   'explicit_numerical_Matchete_GammaCC_phase_verified':False,
   'physical_source_certified':False
  }
 report={'format':'t3-b1-independent-complex-uv-phase-stress-v1',
  'status':'ComplexSourcePhaseCovarianceAndIndependentLoopCrosscheck_PendingNativeSpinorPhase',
  'input_files_sha256':{k:fingerprint(v) for k,v in inputs.items()},
  'benchmark_sha256':{k:fingerprint(root/v) for k,v in MODEL_FILES.items()},
  'native_spinor_contract':shared,'models':results,
  'no_provisional_amplitude_or_engine_used':True,
  'installed_GammaCC_absolute_phase_measured':False,
  'benchmark_portal_and_mass_mapping_conditional':True,
  'physical_source_certified':False,'physics_certified_models':0,
  'matching_or_running_performed':False}
 out.mkdir(parents=True,exist_ok=False)
 target=out/'B1_complex_spinor_convention_witness.json'
 target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
 for name,v in results.items():
  print(f'{name}: complex mnu (eV):\n'+str(np.round(np.asarray(v['complex_nondiagonal_phases_identity_portal']['matrix_ev']['real'])+1j*np.asarray(v['complex_nondiagonal_phases_identity_portal']['matrix_ev']['imag']),10)))
  print(f"  masses (eV): {v['complex_nondiagonal_phases_identity_portal']['takagi_ev_descending']}")
  print(f"  phase covariance residual: {v['global_source_phase_transport']['relative_residual']:.3e}")
 print('Evidence saved:',target)
 print('Absolute installed GammaCC numerical phase: NOT ESTABLISHED; physical source-certified: 0/16')
 return 0


if __name__=='__main__':
 try: raise SystemExit(main())
 except (ValueError,KeyError,TypeError,ArithmeticError,IndexError) as exc:
  print(f'BLOCKED: {exc}',file=sys.stderr)
  raise SystemExit(2)
