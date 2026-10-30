from scripts.BuildT3FrozenLambdaOriginal import collect_scans
from Numerical.plotting.BuildFittedModelEnvelopes import OBS

def test_discovery_two_groups():
    inputs=[]
    fits=[{'model':'A_+0','case':'small_uv_frozen_5e-9','chi2':0.0}]
    for d in ('T1','T2','T3'):
        for g in (0.9,1.0,1.1):
            inputs.append(dict(model='A_+0',case='small_uv_frozen_5e-9',direction=d,gamma=g,**{k:float(g) for k in OBS}))
    models,skipped=collect_scans(inputs,fits)
    assert len(models)==1 and not skipped
    assert models[0]['chi2']==0.0
