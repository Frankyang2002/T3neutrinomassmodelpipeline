"""Read-only v23 B/C six-of-seven source-axis order audit.

Compare the *named* seven historical Matchete Stage22 source directions in the
saved v19 conditional map with v18's independently exported native CG order.
Reorder both the 6x7 map columns and its 7-vector kernel, without changing
original witness files or treating the conditional map as physical.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np

from Numerical.parallel_t3.wilson.bc_mass_eigen_native_bridge_v21 import (
    BCBridgeBlocked, BC, V19_FORMAT, V19_STATUS, V20_FORMAT, V20_STATUS,
    read_json, digest, require, verify_v19_input_sha256,
    special_source_permutation, validate_special_matrix,
)
from Numerical.parallel_t3.wilson.complex_scalar_bridge_v19 import native_document

BASE = Path('output/parallel_t3/wilson')


def preflight(root: Path) -> dict:
    root = Path(root).resolve()
    top_path = root/BASE/'matchete_topology_native_v18.json'
    cg_path = root/BASE/'matchete_bc_extra_native_v18.json'
    v19_path = root/BASE/'complex_scalar_bridge_conditional_v19.json'
    v20_path = root/BASE/'full_scalar_census_conditional_v20.json'
    top, native, v19, v20 = map(read_json, (top_path, cg_path, v19_path, v20_path))
    native_document(top,root,'TOPOLOGY')
    native_document(native,root,'BCEXTRA')
    require(v19.get('format') == V19_FORMAT and v19.get('status') == V19_STATUS and
            v19.get('baseline_EFT1_complex_slots') == 28 and
            v19.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'v23: source v19 conditional map has wrong identity/status')
    hashes = verify_v19_input_sha256(root,v19,(top_path,cg_path))
    require(v20.get('format') == V20_FORMAT and v20.get('status') == V20_STATUS and
            v20.get('total_EFT1_scalar_quartic_registration_slots') == 178 and
            v20.get('input_sha256',{}).get('complex_bridge_v19') == digest(v19_path) and
            v20.get('certification',{}).get('physical_Wilson_evolution_authorized') is False,
            'v23: v20 source hash or 178-slot census has changed')
    rows = {}
    for model in BC:
        require(model in native.get('models',{}) and model in top.get('models',{}),
                'v23: missing native B/C source '+model)
        records=native['models'][model].get('records',[])
        require(len(records)==7 and all(isinstance(rec,dict) for rec in records),
                'v23: missing seven native Matchete source records')
        source=[rec.get('symbol') for rec in records]
        old, perm=special_source_permutation(v19,model,source)
        B,N=validate_special_matrix(v19,model,source)
        rows[model]={
            'v19_historical_source_order':old,
            'v18_native_source_order':source,
            'v19_column_indices_in_native_order_zero_based':perm,
            'reordered':perm!=list(range(7)),
            'baseline_rank_six':int(np.linalg.matrix_rank(B,tol=1e-10)),
            'kernel_norm':float(np.linalg.norm(N)),
            'kernel_residual':float(np.linalg.norm(B@N)),
        }
    return {'status':'SourceNamedSevenColumnReindexingChecked_NoPhysicalMatching',
            'v19_input_sha256_verified':hashes,
            'models':rows,'physical_matching_authorized':False}


def main(argv=None) -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=Path('.'))
    ns=ap.parse_args(argv)
    try:
        report=preflight(ns.project_root)
        print('v23 source-order check: both B/C historical-v19 and native-v18 coupling sets match')
        for model,row in report['models'].items():
            print('  {}: native index -> historical v19 column (zero-based) {} ; '
                  'rank {}; ||B N||={:.3e}'.format(
                      model,row['v19_column_indices_in_native_order_zero_based'],
                      row['baseline_rank_six'],row['kernel_residual']))
        print('Source-only matrix column reordering; no physical matching or Wilson evolution')
        return 0
    except (BCBridgeBlocked,ValueError,KeyError,TypeError,OSError) as exc:
        print('v23 source-order BLOCKED:',exc,file=sys.stderr)
        return 7

if __name__=='__main__':
    raise SystemExit(main())
