#!/usr/bin/env python3
"""Audit the completed A(-4) three-F/joint-two-S direct hard/soft channel.

No Wolfram/Matchete invocation. Reuses the successfully saved source outputs
and the prior full-general Batch 1 running boundary. Writes ONLY a new
strictly provisional diagnostic JSON, NEVER a physical C5 input document.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Numerical.parallel_t3.complete.a4_direct_hard_soft import (
    DirectBridgeError, audit_direct_sources, digest_file,
    formal_numeric_direct_remainder,
)

A4 = ROOT/'output/parallel_t3/full_general_final_one_batch/thesis_a4_fresh_f_v3/A_-4'
B1 = ROOT/'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/A_-4'
DEFAULT_OUT = ROOT/'output/parallel_t3/a4_direct_hard_soft_bridge_v1/A_-4'


def _json(path: Path) -> dict:
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise DirectBridgeError('Not a JSON object: '+str(path))
    return data


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise DirectBridgeError(message)


def run(a4: Path, batch1: Path) -> dict:
    status_path = a4/'native_two_threshold_status.json'
    audit_path = a4/'A_joint_scalar_lagrangian_audit.json'
    match_path = a4/'A_joint_scalar_match_diagnostics.json'
    hard_path = a4/'data/eft1_after_F_one_loop_inputform.txt'
    soft_path = a4/'A_scalar_soft_only_unrenormalized.txt'
    w_path = a4/'A_EFT1_Wilson_definitions.json'
    f_path = batch1/'fermion_matching_input.json'
    s_path = batch1/'scalar_matching_input.json'
    all_paths = {
        'fresh_native_status': status_path,
        'fresh_lagrangian_audit': audit_path,
        'fresh_joint_match_audit': match_path,
        'fresh_full_F_one_loop_source': hard_path,
        'fresh_scalar_soft_Weinberg_source': soft_path,
        'fresh_three_flavour_W_definition': w_path,
        'batch1_fermion_endpoint': f_path,
        'batch1_scalar_endpoint': s_path,
    }
    for p in all_paths.values():
        _require(p.is_file(), 'Required source missing: '+str(p))
    status = _json(status_path)
    lag = _json(audit_path)
    match = _json(match_path)
    w = _json(w_path)
    _require(status.get('status') ==
             'A_PhysicalScalarSoftOnly_HardInheritanceAndWilsonRGPending',
             'Fresh A4 joint scalar-soft run did not succeed')
    _require(status.get('model') == 'A_-4' and
             status.get('heavy_fermion_flavours') == 3 and
             status.get('fermion_groups') == [['F']] and
             status.get('scalar_groups') == [['S1','S2']] and
             status.get('first_threshold_reused') is False and
             status.get('A_joint_scalar_match_performed') is True and
             status.get('A_physical_scalar_soft_Weinberg_terms') == 10 and
             status.get('exploratory_source_theory_used') is False,
             'Full three-F/fresh-F/joint-two-S source conditions changed')
    _require(lag.get('CheckLagrangian_original') == 'True' and
             lag.get('accepted_canonicalisation') == 'ExactSavedTree' and
             lag.get('projection_used_for_physical_match') is False,
             'A4 saved native EFT1 did not validate exactly')
    _require(match.get('matching_at_one_loop_succeeded') is True and
             match.get('source_theory') == 'ExactSavedTree' and
             match.get('exploratory_non_equivalent_source') is False,
             'Native scalar matching was not source-exact')
    _require(w.get('all_three_fermions_integrated_together') is True and
             w.get('W12[p,q]') == 'Sum[y1[p,r]*y2[q,r]/MF[r],{r,1,3}]',
             'F boundary three-flavour Wilson definition changed')
    symbolic = audit_direct_sources(
        hard_path.read_text(encoding='utf-8-sig'),
        soft_path.read_text(encoding='utf-8-sig'))
    numeric = formal_numeric_direct_remainder(_json(f_path), _json(s_path))
    return {
        'format':'t3-a4-direct-hard-soft-source-bridge-v1',
        'status':'SourceVerifiedFormalDirectChannel_PolePairedNotFullC5',
        'model':'A_-4',
        'fresh_first_threshold_recomputed':True,
        'matching_engine':'Matchete; existing exact saved outputs only; no re-execution',
        'source_files':{key:{'project_path':str(path.relative_to(ROOT))
                             if path.is_relative_to(ROOT) else str(path),
                             'sha256':digest_file(path)}
                        for key, path in all_paths.items()},
        'source_audit':asdict(symbolic),
        'numerical_formal_channel':numeric,
        'full_F_finite_MSbar_Wilson_boundary_computed':False,
        'all_F_one_loop_operator_mixing_inherited':False,
        'one_loop_EFT1_Wilson_RG_performed':False,
        'full_joint_scalar_hard_plus_soft_matching_performed':False,
        'physical_C5_computed':False,
        'next_physics_blocker':(
            'Prove complete operator and SU2/Spinor basis mapping, compute local '
            'MSbar counterterms and EFT1 Wilson mixing, run Wilson RG, and '
            'assemble full finite F-hard + S-soft at common scalar threshold.'),
    }


def main(argv: list[str] | None=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--a4-dir', type=Path, default=A4)
    p.add_argument('--batch1-a4-dir', type=Path, default=B1)
    p.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = p.parse_args(argv)
    dest = args.output_dir/'A4_direct_hard_soft_formal_audit.json'
    if dest.exists():
        p.error('Refusing to overwrite existing result: '+str(dest))
    try:
        report = run(args.a4_dir.resolve(), args.batch1_a4_dir.resolve())
    except (OSError, ValueError, KeyError, TypeError, DirectBridgeError) as exc:
        p.exit(2, 'BLOCKED: '+str(exc)+'\n')
    dest.parent.mkdir(parents=True, exist_ok=True)
    temporary = dest.with_suffix('.json.partial')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True,
                                    allow_nan=False)+'\n', encoding='utf-8')
    temporary.replace(dest)
    norm = report['numerical_formal_channel']['ordered_formal_finite_frobenius_GeVinv']
    print('A(-4) direct hard/soft Laurent source audit: PASSED')
    print('A(-4) ordered formal finite remainder norm [GeV^-1]:',f'{norm:.9e}')
    print('Full EFT matching/renormalised physical C5: NOT COMPUTED')
    print('Output:',dest)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
