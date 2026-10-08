from Numerical.diagnostics.ContractT3LoopComponents import candidates,analyze,pairing_map


def _case():
    def vertex(indices):
        return {'components':[{'indices':list(k),'cg':v} for k,v in indices]}
    return {'model':{'dS1':2,'dS2':2,'dF':1,'alpha':-1},'vertices':{
        'T3Y1CG':vertex([((1,1,1),'1')]),
        'T3Y2CG':vertex([((1,1,2),'1')]),
        'T3MixCG':vertex([((2,2,1,2),'sqrt(3)/2')])}}


def test_candidate_weights_preserved_and_no_multiplicity():
    out=candidates(_case())
    assert out['channel_count']==1
    assert out['complex_sum']=='sqrt(3)/2'
    assert out['loop_multiplicity'] is None
    assert out['compensation_factor'] is None


def test_fermion_epsilon_has_antisymmetric_sign():
    assert pairing_map(2,'epsilon')[(1,2)]==1
    assert pairing_map(2,'epsilon')[(2,1)]==-1


def test_candidate_analysis_explicitly_not_physical():
    out=analyze({'cases':[_case()]})
    assert out['count']==1
    assert out['models'][0]['verified_physical_multiplicity'] is None
