from validation.ValidateIntermediateWilsonCouplings import validate, symbols_in_matrix

def test_unknown_symbol_fails_closed():
    v={'anomalous_dimension_matrix':{'gamma_5_5':'2*λ+x'}}
    result=validate(v,{'ScalarQuarticTerms':[]},{'report_betas':{'lambdaH':'0'}})
    assert result['status']=='NeedsMappingReview'
    assert set(result['unknown_symbol_names'])=={'λ','x'}

def test_preserves_raw_matchete_symbol_and_evidence():
    v={'anomalous_dimension_matrix':{'gamma_5_5':'2*λ'}}
    quartic={'ScalarQuarticTerms':[{'Name':'SM','TermInputForm':'Coupling[λ, {}, 0]*Field[H, Scalar, {}]^4'}]}
    result=validate(v,quartic,{'report_betas':{'lambdaH':'0'}})
    assert result['candidate_mapping']=={'λ':'lambdaH'}
    assert result['status']=='CandidateMappingNeedsNormalizationCheck'
