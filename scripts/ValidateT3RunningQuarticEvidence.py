"""Validate full Cartesian quartic response evidence against a model registry.

Usage: python -m scripts.ValidateT3RunningQuarticEvidence \
    --model A_-4 --responses output/.../A_-4_cartesian_responses.json \
    --state output/.../A_-4_state.json --registry RGE/general/generated/general_rgbeta_quartic_registry.json \
    --baseline-symbols lambdaH,lambdaS1,lambdaS2,lambda12,lambdaH1,lambdaH2,lambdaT3

The --state JSON is a flat {coupling_name: value} mapping; values can be
numbers or {real, imag}. This command does not evolve RGEs or infer tensors.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from Numerical.running.T3RunningQuarticTensor import (
    RunningQuarticBasis, expected_scalar_symbols
)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',required=True)
    p.add_argument('--responses',required=True,type=Path)
    p.add_argument('--state',required=True,type=Path)
    p.add_argument('--registry',required=True,type=Path)
    p.add_argument('--baseline-symbols',required=True,
                   help='Comma-separated complete baseline quartic coupling names')
    p.add_argument('--out',type=Path,default=None)
    a=p.parse_args(argv)
    registry=json.loads(a.registry.read_text(encoding='utf-8-sig'))
    baseline={x.strip() for x in a.baseline_symbols.split(',') if x.strip()}
    expected=expected_scalar_symbols(registry,a.model,baseline)
    data=json.loads(a.responses.read_text(encoding='utf-8-sig'))
    if data.get('model') != a.model:
        raise ValueError('Response model mismatch')
    basis=RunningQuarticBasis.from_payload(data,expected)
    state=json.loads(a.state.read_text(encoding='utf-8-sig'))
    tensor=basis.evaluate(state)
    result={'status':'ResponseReconstructionPassedNotPhysicalCertification',
            'model':a.model,'n_real_scalars':int(tensor.shape[0]),
            'n_couplings':len(basis.symbols),
            'n_directions':len(basis.directions),
            'max_abs_lambda':float(abs(tensor).max()),
            'source':basis.source,
            'physical_normalisation_certified':False}
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
