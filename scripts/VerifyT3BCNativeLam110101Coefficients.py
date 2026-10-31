"""Validate live RGBeta Lam 110101 component evidence, never promote physical matching.

Run from project root after the Wolfram probe:
python -m scripts.VerifyT3BCNativeLam110101Coefficients --input output/audits/t3_bc_native_lam/native_lam_coeff_v1.json
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

FORMAT = 't3-bc-native-rgbeta-lam-110101-component-coefficients-v1'
MODELS = ('B_-1', 'C_-1')
SYMBOLS = ('lambdaT3GeneralQ110101I1', 'lambdaHHdagS1barS2barCross')
STATUS = 'LiveLam110101SourceComponentCoefficientsPassed'


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError('Live RGBeta Lam source coefficient audit: ' + message)


def verify(document: dict[str, Any]) -> dict[str, Any]:
    require(isinstance(document, dict) and document.get('format') == FORMAT,
            'wrong source evidence format')
    require(document.get('status') == 'LiveLam110101SourceCoefficientComponentAudit',
            'wrong native stage status')
    require(document.get('source') == 'installed RGBeta`PackageScope`Lam',
            'not independent runtime Lam')
    for key in ('action_level_matchete_independently_evaluated',
                'physical_coupling_mapping_certified', 'physical_matching_certified'):
        require(document.get(key) is False, 'unexpected global certification: '+key)
    models = document.get('models')
    require(isinstance(models, dict) and set(models) == set(MODELS),
            'missing/extra B/C cases')
    rows = []
    for model in MODELS:
        row = models[model]
        require(isinstance(row,dict) and row.get('model')==model,model+': wrong model')
        require(row.get('status') == STATUS, model+': native source coefficient check failed')
        require(row.get('live_native_lam_called') is True and
                row.get('symbol') == 'RGBeta`PackageScope`Lam',
                model+': native package function not called')
        require(row.get('source_symbols') == list(SYMBOLS),
                model+': selected 110101 quartics changed')
        require(row.get('registered_complex_quartic_normalization') == '6*Matrix[...]' and
                row.get('unique_arrangements_per_coupling') == 24,
                model+': source normalization/arrangement discrepancy')
        require(row.get('expected_lam_over_registered_invariant') == '1/4',
                model+': wrong expected native coefficient')
        require(row.get('coefficient_extraction') ==
                'symbolic Coefficient on evaluated native Lam',
                model+': did not use native evaluated coupling coefficients')
        forms=row.get('native_evaluated_coupling_forms')
        require(isinstance(forms,list) and len(forms)==4,
                model+': four direct/conjugate coefficient forms required')
        expected_pairs={(orientation,symbol) for orientation in ('direct','conjugate')
                        for symbol in SYMBOLS}
        require({(x.get('orientation'),x.get('coupling')) for x in forms
                 if isinstance(x,dict)} == expected_pairs,
                model+': incomplete extracted form inventory')
        require(all(isinstance(x,dict) and isinstance(x.get('evaluated_coupling_form'),str)
                    and x['evaluated_coupling_form'] for x in forms),
                model+': missing native coefficient form')
        require(row.get('component_tests') == 64 and
                row.get('expected_component_tests') == 64 and
                row.get('failed_component_count') == 0 and
                row.get('first_mismatches') == [],
                model+': incomplete or failed component comparison')
        positions = row.get('source_scalar_field_positions')
        require(isinstance(positions,list) and len(positions)==2 and
                all(isinstance(r,list) and len(r)==4 and
                    all(type(v)==int and v>0 for v in r) for r in positions),
                model+': native scalar positions missing')
        residual=row.get('maximum_absolute_component_residual')
        require(isinstance(residual,(int,float)) and 0 <= residual < 1e-12,
                model+': coefficient residual too large')
        for key in ('matchete_runtime_evaluated','physical_coupling_mapping_certified',
                    'physical_matching_certified'):
            require(row.get(key) is False,model+': unexpected physical certification '+key)
        rows.append({'model':model,'status':STATUS,
                     'component_tests':64,'max_residual':float(residual)})
    return {
        'format':'t3-bc-native-lam-110101-component-verification-v1',
        'status':'NativeLam110101ComponentCoefficientsVerifiedNotPhysicalMapping',
        'results':rows,
        'total_components_checked':128,
        'physical_coupling_mapping_certified':False,
        'matchete_runtime_evaluated':False,
        'physical_matching_certified':False,
    }


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=args.input
    if not source.is_file():
        raise FileNotFoundError(source)
    payload=json.loads(source.read_text(encoding='utf-8-sig'))
    result=verify(payload)
    result['input_file']=source.as_posix()
    result['input_sha256']=sha256(source.read_bytes()).hexdigest()
    if args.output:
        if args.output.exists():
            raise FileExistsError('Refusing to overwrite '+str(args.output))
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for row in result['results']:
        print(f"{row['model']}: live Lam component coefficients {row['component_tests']}/64; "
              f"max residual={row['max_residual']:.3e}")
    print('Native component coefficients: PASSED; full physical matching: NOT CERTIFIED')
    if args.output:
        print('Summary:',args.output)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
