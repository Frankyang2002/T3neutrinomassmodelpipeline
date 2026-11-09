"""Audit 16-model native portal CG phases and conditional Matchete/RGBeta map.

Runs after BOTH v5 native Wolfram exporters; refuses stale source hashes,
missing A(-4) I phases, rank-deficient potentials, and output overwrites.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from Numerical.parallel_t3.wilson.portal_tensor_audit import audit_from_paths
from Numerical.parallel_t3.wilson.matchete_rgbeta_dictionary import DictionaryEvidenceError


def main(argv: list[str] | None = None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matchete',type=Path,required=True)
    p.add_argument('--rgbeta',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--project-root',type=Path,default=Path.cwd())
    p.add_argument('--samples',type=int,default=5)
    a=p.parse_args(argv)
    try:
        if a.output.exists():
            raise DictionaryEvidenceError(f'Refusing to overwrite {a.output}')
        result=audit_from_paths(a.matchete,a.rgbeta,a.project_root,samples=a.samples)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(f"{result['portal_coupling_slots']}/58 source portal operators, "
              f"{result['model_count']}/16 models, "
              f"{result['native_anti_Hermitian_CGs_given_I_inside_PlusHc']} "
              f"anti-Hermitian phases; residual="
              f"{result['max_relative_tensor_and_potential_residual']:.3e}")
        print(f"{result['status']} -> {a.output}")
        print('NOT a physically certified RGBeta generator orientation or scalar dictionary.')
        return 0
    except (DictionaryEvidenceError,OSError,ValueError,KeyError,TypeError,
            json.JSONDecodeError) as exc:
        print(f'Portal v5 BLOCKED: {exc}',file=sys.stderr)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
