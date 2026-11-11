"""Two-stage coupled, conjugate-linear Wilson IVP numerical invariance experiment.

Independent synthetic laboratory system, NOT the production T3 beta functions.
It tests IVP uniqueness under common real scaling across a stage transition.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
try:
    from scripts.worker3_normalisation.audit import obs,synthetic_matrix,encode
except ModuleNotFoundError:
    from audit import obs,synthetic_matrix,encode


def pack(m):return np.r_[m.real.flatten(),m.imag.flatten()]
def unpack(y):return y[:9].reshape(3,3)+1j*y[9:].reshape(3,3)

def integration(initial,scale=1.0,phase=1.0+0j,channels=(1.0,1.0)):
    """Stage I: C'=a C + b C* + k D ; D'=u D + v D* + w C.
    Stage II: C'=c C + d C*; matching C+= C- + 0.3 D-.
    All coefficients synthetic, all rates real except pure imaginary mixing.
    """
    c,d=initial
    c=scale*phase*channels[0]*c
    d=scale*phase*channels[1]*d
    def rhs(t,y):
        x=unpack(y[:18]);z=unpack(y[18:])
        dx=(-0.12+0.03j)*x + (0.04-0.01j)*x.conj()+0.2*z
        dz=0.05*x + (-0.17+0.04j)*z+(0.02+0.03j)*z.conj()
        return np.r_[pack(dx),pack(dz)]
    result=solve_ivp(rhs,(0,1.7),np.r_[pack(c),pack(d)],rtol=1e-11,atol=1e-12,dense_output=True)
    if not result.success:raise RuntimeError(result.message)
    cs=unpack(result.y[:18,-1]);ds=unpack(result.y[18:,-1]);matched=cs+0.3*ds
    def rhs2(t,y):
        x=unpack(y)
        return pack(-0.21*x + (0.07+0.02j)*x.conj())
    result2=solve_ivp(rhs2,(1.7,3),pack(matched),rtol=1e-11,atol=1e-12)
    if not result2.success:raise RuntimeError(result2.message)
    return [c,cs,matched,unpack(result2.y[:,-1])]


def run(out,seed=8301,n=16):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cases=[]
    for i in range(n):
        x=(synthetic_matrix(seed+i),synthetic_matrix(seed+100+i))
        reference=integration(x)
        real=integration(x,scale=-2.4)
        complex_rescale=integration(x,phase=0.6+0.9j)
        unequal=integration(x,channels=(1.8,0.75))
        errors=[];angle_error=[];ratio_error=[];stage_errors=[]
        for a,b in zip(reference,real):
            errors.append(float(np.linalg.norm(b+2.4*a)))
            oa,ob=obs(a),obs(b)
            angle_error.append(float(np.max(np.abs(np.array(oa['sin2_angles'])-ob['sin2_angles']))))
            ratio_error.append(float(np.max(np.abs(np.array(oa['ratios'])-ob['ratios']))))
        for j in range(1,len(reference)):
            p=np.linalg.norm(reference[j])/np.linalg.norm(reference[j-1])
            q=np.linalg.norm(real[j])/np.linalg.norm(real[j-1])
            stage_errors.append(float(abs(p-q)))
        complex_obstruction=float(np.linalg.norm(complex_rescale[-1]-(0.6+0.9j)*reference[-1]))
        different_channel_flavour=float(np.max(np.abs(np.asarray(obs(unequal[-1])['sin2_angles'])-obs(reference[-1])['sin2_angles'])))
        assert max(errors)<1e-8
        assert max(angle_error)<1e-8
        assert max(ratio_error)<1e-8
        assert max(stage_errors)<1e-12
        assert complex_obstruction>1e-6
        cases.append({'id':i,'common_real_state_error_max':max(errors),
         'common_real_angles_error_max':max(angle_error),
         'common_real_ratios_error_max':max(ratio_error),
         'stage_norm_ratio_error_max':max(stage_errors),
         'complex_rescaling_state_obstruction':complex_obstruction,
         'unequal_channel_angle_shift':different_channel_flavour,
         'source':'independent synthetic two-stage real-linear coupled IVP'})
    output={'format':'worker3-synthetic-coupled-ivp-v1','count':n,'pass':True,
      'physical_t3_ivp':False,'production_ode_rerun':False,'cases':cases}
    (out/'coupled_ivp_regression.json').write_text(json.dumps(output,indent=2))
    return output

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--output',default='output/audits/worker3_normalisation')
    z=a.parse_args();result=run(z.output)
    print(json.dumps({'cases':result['count'],'pass':result['pass'],'physical_t3_ivp':False}))
