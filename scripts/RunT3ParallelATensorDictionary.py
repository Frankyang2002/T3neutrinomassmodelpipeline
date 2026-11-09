"""Run Batch1 independent-polynomial Matchete/RGBeta matrix diagnostics.

Only accept FULL native components *and* explicit common real field axes,
normalization, signed potential convention and H.c. terms. The old RGBeta
DeclarationMetadataOnly export is rejected before any numeric solving.
"""
from __future__ import annotations
from pathlib import Path
import argparse,json


def validated_records(path:Path,model:str,kind:str):
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if data.get('format')=='stage22-native-rgbeta-declarations-v1':
        raise ValueError('Blocked_RGBetaDeclarationsAreNotComponentTensors')
    if data.get('format')!='t3-parallel-a-full-native-polynomials-v1':
        raise ValueError('Blocked_WrongComponentTensorFormat')
    if data.get('scalar_frame')!='RGBeta_reference_complex_highest_weight':
        raise ValueError('Blocked_ComponentFrameNotProvedEqual')
    if data.get('normalizations_explicit_and_native_verified') is not True:
        raise ValueError('Blocked_QuarticNativePrefactorAndPotentialSignUnverified')
    if model not in data.get('models',{}):
        raise ValueError('Blocked_MissingModelTensorComponents: '+model)
    records=data['models'][model].get('records',[])
    if not records:raise ValueError('Blocked_Zero'+kind+'PolynomialRecords')
    for record in records:
        if not all(k in record for k in ('species','dimensions','components',
                                        'potential_prefactor','phase','add_hc','source_is_lagrangian')):
            raise ValueError('Blocked_PolynomialRecordIncomplete: '+record.get('symbol','?'))
    return records


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',required=True)
    parser.add_argument('--matchete-components',type=Path,required=True)
    parser.add_argument('--rgbeta-components',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--samples',type=int,default=256)
    a=parser.parse_args(argv)
    if a.output.exists():raise FileExistsError(a.output)
    m=validated_records(a.matchete_components,a.model,'Matchete')
    r=validated_records(a.rgbeta_components,a.model,'RGBeta')
    from Numerical.parallel_t3.wilson.dictionary import group_fit
    result=group_fit(m,r,samples=a.samples)
    result.update({'model':a.model,
                   'matrix_maps':'RGBeta coupling coefficients = map_rgbeta_from_matchete @ Matchete coefficients',
                   'matchete_native_file':str(a.matchete_components),
                   'rgbeta_native_file':str(a.rgbeta_components),
                   'spinor_signed_factor_certified':False,
                   'mass_Yukawa_normalizations_certified':False,
                   'physical_certificate':False,
                   'antiHermitianPortal_Aminus4':
                      'keep I inside PlusHc; do not replace ScalarPotential.wl'})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(a.model,result['status'],'rank_M',result['rank_matchete'],
          'rank_R',result['rank_rgbeta'],'residual',result['holdout_relative_residual'])
    return int(result['status']!='Candidate_ExactPolynomialSubspace')


if __name__=='__main__':raise SystemExit(main())
