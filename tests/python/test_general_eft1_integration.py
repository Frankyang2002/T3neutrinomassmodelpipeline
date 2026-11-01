import numpy as np
import pytest
from Numerical.core.State import SMNumericalState,T3Representation,T3UVState
from Numerical.core.StateVector import pack_uv_state
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state,unpack_intermediate_scalar_state
from Numerical.running.IntermediateScalarRGBetaEvaluator import evaluate_intermediate_scalar_rgbeta_payload

def sample():
    r=T3Representation(1,3,2,-2);zeros=np.zeros((3,3),complex)
    sm=SMNumericalState(.36,.65,1.17,.13,zeros,zeros,zeros)
    return T3UVState(mu_gev=1e6,representation=r,sm=sm,y1=zeros,y2=zeros,MF=np.eye(3)*1e5,
        mS1Sq=1e6,mS2Sq=2e6,lambdaS1=.1,lambdaS2=.2,lambdaH1=.03,lambdaH2=.04,
        lambda12=.05,lambdaT3=.01+.02j,lambdaH2Adj=.01,lambdaS2Adj=.02,
        general_couplings={'lambdaT3GeneralQ000004I1':.12+.3j,
          'muT3GeneralS2Sq':500+100j})

def test_projection_preserves_bilinear_and_quartic():
    uv=sample().validated();e=build_intermediate_scalar_boundary(uv)
    assert e.general_couplings['muT3GeneralS2Sq']==500+100j
    v,layout=pack_intermediate_scalar_state(e);rev=unpack_intermediate_scalar_state(v,layout)
    assert rev.general_couplings==e.general_couplings
    assert len(v)>0 and 'y1' not in layout.names

def test_stage_and_complete_coverage():
    e=build_intermediate_scalar_boundary(sample());_,layout=pack_intermediate_scalar_state(e)
    r=e.representation
    meta={'SharedScalar':False,'dS1':r.d_s1,'dS2':r.d_s2,'dF':r.d_f,'alpha':r.alpha,
          'IntegratedField':'F','ActiveBSMFields':['S1','S2'],'GeneralScalarStage':'EFT1',
          'ReportBetaConvention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g'}
    payload={'status':'SuccessCandidateMassBeta','format':'general-t3-quadratic-rgbeta-experimental-v1',
             'metadata':meta,'report_betas':{n:'0' for n in layout.names}}
    out=evaluate_intermediate_scalar_rgbeta_payload(payload,e)
    assert out['muT3GeneralS2Sq']==0
    payload['report_betas'].pop('muT3GeneralS2Sq')
    with pytest.raises(ValueError,match='missing'):
        evaluate_intermediate_scalar_rgbeta_payload(payload,e)
    payload['metadata']['GeneralScalarStage']='UV'
    with pytest.raises(ValueError,match='stage'):
        evaluate_intermediate_scalar_rgbeta_payload(payload,e)
