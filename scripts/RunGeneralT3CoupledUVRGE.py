"""Stage 13: coupled B/C(-1) UV RGBeta RGE and tensor identity audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Numerical.experimental.GeneralT3CoupledUV import run_coupled_uv
from Numerical.experimental.GeneralT3TensorMap import build_tensor_mapping


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, default=Path('docs/general_c5_toy_benchmark.json'))
    ap.add_argument('--matchete-tensors', type=Path, default=Path('output/general_t3_trial/matchete_active_quartic_tensors.json'))
    ap.add_argument('--rgbeta-model', type=Path, default=Path('RGE/running/rgbeta/T3RGBetaModel.wl'))
    ap.add_argument('--source-proof', type=Path, default=Path('output/general_rgbeta_source_proof.json'))
    ap.add_argument('--rgbeta-root', type=Path, default=Path('output'))
    ap.add_argument('--output', type=Path, default=Path('output/general_t3_stage13'))
    ap.add_argument('--tensor-only', action='store_true')
    ns = ap.parse_args()
    ns.output.mkdir(parents=True, exist_ok=True)
    errors = []
    try:
        mapping = build_tensor_mapping(ns.matchete_tensors, ns.rgbeta_model)
        (ns.output / 'doublet_tensor_map.json').write_text(json.dumps(mapping, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        print('Doublet SU(2) tensor identities:', mapping['status'])
        for v in mapping['mappings']:
            print(' ', v['model'], v['matchete_coupling'], 'RGBeta/Matchete tensor factor =', v['rgbeta_over_matchete_tensor_factor'])
    except (ValueError, OSError, KeyError) as e:
        errors.append({'stage':'tensor_map', 'error':str(e)})
        print('BLOCKED tensor map:', e, file=sys.stderr)

    if not ns.tensor_only:
        pars = json.loads(ns.config.read_text(encoding='utf-8-sig'))
        for model in ('B_m1', 'C_m1'):
            klass = model[0]
            beta = ns.rgbeta_root / f'T3_{klass}_alpha_m1' / 'data' / 'uv_rgbeta_rge.json'
            if not beta.is_file():
                beta = ns.rgbeta_root / 'general_t3_trial' / model / 'baseline' / 'data' / 'uv_rgbeta_rge.json'
            try:
                result = run_coupled_uv(model, pars, beta, ns.source_proof)
                path = ns.output / f'{model}_coupled_uv.json'
                path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
                print(model, 'coupled UV evolution PASSED, nfev=', result['solver_nfev'])
                for row in result['active_couplings']:
                    print(' ',row['coupling'],row['change_coupled'])
                print('Report:',path)
            except (ValueError, KeyError, TypeError, OSError, RuntimeError) as e:
                errors.append({'stage':model+'_coupled_uv', 'source':str(beta), 'error':str(e)})
                print('BLOCKED',model,'coupled UV evolution:',e,file=sys.stderr)
    summary={'status':'Success' if not errors else 'ReviewRequired','errors':errors,'full_c5_matching_rge_verified':False}
    (ns.output/'stage13_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return 2 if errors else 0


if __name__ == '__main__':
    raise SystemExit(main())
