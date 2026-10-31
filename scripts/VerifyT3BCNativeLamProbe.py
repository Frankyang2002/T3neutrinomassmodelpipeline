"""Fail-closed inspection of a *live* B/C RGBeta Lam diagnostic JSON.

Checks what the local Wolfram runtime actually returned. This does NOT
certify a potential-level mapping, Matchete equivalence, or one-loop matching.

python -X utf8 -m scripts.VerifyT3BCNativeLamProbe \
    --input output/audits/t3_bc_native_lam/native_lam_v1.json \
    --output output/audits/t3_bc_native_lam/native_lam_v1_verified.json
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path

EXPECTED = ('lambdaT3GeneralQ110101I1', 'lambdaHHdagS1barS2barCross')
MODELS = ('B_-1', 'C_-1')
FORMAT = 't3-bc-native-rgbeta-lam-diagnostic-v1'


def verify(doc: dict) -> dict:
    if doc.get('format') != FORMAT or doc.get('source') != 'installed RGBeta runtime Lam':
        raise ValueError('Not a native installed-RGBeta Lam diagnostic')
    if doc.get('physical_coupling_mapping_certified') is not False or doc.get('physical_matching_certified') is not False:
        raise ValueError('Unexpected physical certification claim')
    models = doc.get('models')
    if not isinstance(models, dict) or set(models) != set(MODELS):
        raise ValueError('Expected B_-1 and C_-1')
    rows = []
    for model in MODELS:
        m = models[model]
        if m.get('model') != model or m.get('actual_installed_rgbeta_lam_called') is not True:
            raise ValueError(f'{model}: missing native runtime evaluation')
        if m.get('native_registered_symbols') != list(EXPECTED):
            raise ValueError(f'{model}: unexpected quartic inventory')
        if m.get('lam_action_potential_reconstructed') is not False or m.get('physical_coupling_mapping_certified') is not False or m.get('matchete_runtime_evaluated') is not False:
            raise ValueError(f'{model}: diagnostic contains an invalid certification claim')
        if not isinstance(m.get('lam_expression_sha256'), str) or len(m['lam_expression_sha256']) != 64:
            raise ValueError(f'{model}: missing full expression hash')
        factors = m.get('native_coupling_factors')
        if not isinstance(factors, list) or len(factors) != 2 or any('6*Matrix['+name+']' not in s for name, s in zip(EXPECTED, factors)):
            raise ValueError(f'{model}: registered native factors changed')
        if m.get('native_self_conjugate') != [False,False]:
            raise ValueError(f'{model}: conjugation metadata changed')
        if m.get('status') != 'NativeLamExpressionAndSparseTensorCaptured' or m.get('sparse_array_count', 0) < 1 or not all(m.get('lam_contains_couplings', {}).get(name) is True for name in EXPECTED):
            raise ValueError(f'{model}: Lam was evaluated but tensor evidence incomplete; inspect output excerpt before changing code')
        rows.append({'model':model,'native_lam_evaluated':True,
                     'sparse_array_count':m['sparse_array_count'],
                     'lam_expression_sha256':m['lam_expression_sha256']})
    return {
        'status':'NativeRGBetaLamExpressionEvidencePresent',
        'scope':'B/C 110101 structural runtime evaluation ONLY',
        'results':rows,
        'native_potential_coefficient_verified':False,
        'physical_coupling_mapping_certified':False,
        'matchete_runtime_equivalence_certified':False,
        'one_loop_matching_certified':False,
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(argv)
    if not a.input.is_file():
        raise FileNotFoundError(a.input)
    if a.output.exists() or a.input.resolve() == a.output.resolve():
        raise FileExistsError('Output exists or equals input; use another path')
    data = json.loads(a.input.read_text(encoding='utf-8-sig'))
    result = verify(data)
    result['source_sha256'] = sha256(a.input.read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for row in result['results']:
        print(f"{row['model']}: live Lam captured; sparse arrays={row['sparse_array_count']}")
    print('Potential prefactor/physical coupling normalization: STILL UNVERIFIED')
    print('Summary:', a.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
