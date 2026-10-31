"""Field-phase covariance diagnostics in the *RGBeta* left-Weyl basis.

The heavy chiral fields F and Fc have charges (+1,-1) under a vectorlike
flavour-independent rephasing, leaving MF fixed. Since AddYukawa has
Chirality->Right, the named coupling symbols are phase-conjugate to the
left-Weyl vertex coefficients. For zero scalar/lepton phase the expected
coupling transformation is y1,y1Conj -> exp(+iθ)*coupling and
 y2,y2Conj -> exp(-iθ)*coupling.

This is NOT an absolute Matchete-to-RGBeta dictionary: it checks only an
internal U(1) reparameterisation covariance of RGBeta's exported RGEs.
"""
from __future__ import annotations

from typing import Any, Mapping
import numpy as np

from Numerical.experimental.T3GeneralYukawaCoupledUV import DynamicLayout, evaluate_betas, norm, require


def vectorlike_phase_transform(values:Mapping[str,Any], theta:float) -> dict:
    require(bool(np.isfinite(theta)),'Nonfinite rephase angle')
    plus=np.exp(1j*theta)
    minus=np.conjugate(plus)
    transformed={k:(v.copy() if isinstance(v,np.ndarray) else v)
                 for k,v in values.items()}
    for key in ('y1','y1Conj'):
        if key in transformed:transformed[key] *= plus
    for key in ('y2','y2Conj'):
        if key in transformed:transformed[key] *= minus
    return transformed


def test_rephase_covariance(expressions:Mapping[str,str], state:Mapping[str,Any], *,
                            theta:float=0.43, evaluator=evaluate_betas,
                            tol:float=5e-7) ->dict:
    require(set(expressions)==set(state),'Incompatible beta/state keys')
    layout=DynamicLayout(expressions.keys())
    # These are identical group elements acting on input and output.
    b0=evaluator(expressions,state)
    rotated=vectorlike_phase_transform(state,theta)
    b1=evaluator(expressions,rotated)
    expected=vectorlike_phase_transform(b0,theta)
    layout.pack(b0);layout.pack(b1);layout.pack(expected)
    residual={k:norm(b1[k]-expected[k]) / max(norm(b0[k]),norm(b1[k]),1e-10)
              for k in expressions}
    bad={k:v for k,v in residual.items() if v>tol}
    return {
        'status':'RGBetaVectorlikeFieldRephaseCovariant' if not bad else 'RGBetaVectorlikeFieldRephaseDiagnosticFailed',
        'theta_radians':theta,'tolerance':tol,
        'max_relative_residual':max(residual.values(),default=0.),
        'relative_residual_by_coupling':residual,
        'mismatched_couplings':bad,
        'absolute_matchete_rgbeta_coupling_dictionary_certified':False,
        'note':'Independent RGBeta F/Fc vectorlike rephase; does not certify four-to-two component spinor sign, SU2 triplet orientation, or physical threshold matching.',
    }


def self_test():
    raw={'y1':np.eye(3,dtype=complex),'y2':np.eye(3,dtype=complex),
         'MF':np.eye(3,dtype=complex)*1e5,'gY':.37,'lambdaT3':.1+0.01j,
         'y1Conj':np.eye(3,dtype=complex)*.1}
    ex={k:'unused' for k in raw}
    def covariant(expressions,state):
        # U(1) equivariant cubic for each complex Yukawa, identity for others.
        return {k:v*np.linalg.norm(v)**2 if k in ('y1','y2','y1Conj') else
                np.asarray(v)*1.2 if isinstance(v,np.ndarray) else v*1.2
                for k,v in state.items()}
    result=test_rephase_covariance(ex,raw,evaluator=covariant)
    require(result['status']=='RGBetaVectorlikeFieldRephaseCovariant','Expected phase covariance')
    def broken(expressions,state):
        r=covariant(expressions,state)
        r['y1']=r['y1']+state['y2']
        return r
    result=test_rephase_covariance(ex,raw,evaluator=broken)
    require(result['status']=='RGBetaVectorlikeFieldRephaseDiagnosticFailed','Broken phase covariance undetected')
    print('PASS: RGBeta F/Fc mass-preserving rephase, Yukawa ±θ covariance, broken-beta negative control')


if __name__=='__main__':self_test()
