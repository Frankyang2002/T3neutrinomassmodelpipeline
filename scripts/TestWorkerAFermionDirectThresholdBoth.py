"""Regression tests for source-audited A(-4)/D(-2) direct hard/soft diagnostics."""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import sys
import zipfile

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from WorkerAFermionDirectThresholdBoth import (
    SCHEMA, S, LOOP, demand, evaluate_both, evaluate_d, parse_d_soft,
    restricted_hard, real_sympy, assert_restricted_yukawa_inventory,
)
from WorkerAFermionEFT1SoftA4 import scalar_bubble_finite
from WorkerASoftMixingKernel import extract


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def run(sources_zip: Path, soft_zip: Path) -> None:
    both = evaluate_both(sources_zip, soft_zip)
    check(both['schema']==SCHEMA and both['production_boundary_permitted'] is False,
          'Accidental production boundary status')
    a,d = both['results']['A_-4'],both['results']['D_-2']
    check(abs(a['hard_direct_pole_frobenius_GeV_inv']-2.445231389e-16)<2.e-25,
          'A4 baseline hard pole changed')
    check(abs(a['paired_direct_finite_frobenius_GeV_inv']-2.045493967e-15)<2.e-24,
          'A4 baseline paired diagnostic changed')
    check(d['hard_llhh_restriction']['direct_llhh_total']==12 and
          d['hard_llhh_restriction']['direct_llhh_restricted']==6 and
          d['hard_llhh_restriction']['direct_llhh_excluded_y2Conj']==6,
          'D2 restricted hard operator count changed')
    check(abs(d['uv_direct_pole_frobenius_GeV_inv']-2.550267561e-16)<2.e-25,
          'D2 hard pole changed')
    check(abs(d['paired_direct_finite_frobenius_GeV_inv']-2.103739257e-15)<2.e-24,
          'D2 paired direct diagnostic changed')
    check(abs(d['kernel_2x2_hbar1'][0][1]+0.00432064056169)<1.e-11,
          'D2 nonzero mixed kernel changed')
    check(abs(d['candidate_eft1_soft_finite_mixed_frobenius_GeV_inv']-9.220920798e-20)<2.e-28,
          'D2 mixed contribution changed')
    check(d['direct_ordered_pole_residual_GeV_inv']<1e-30 and
          d['explicit_scale_derivative_residual_GeV_inv']<1e-30,
          'D2 pole or scale cancellation changed')
    check(d['restricted_WExtra_is_zero'] is True and
          d['candidate_eft1_soft_finite_mixed_frobenius_GeV_inv']>0,
          'WExtra-zero did not preserve mixed quartic effect')
    check(d['finite_one_loop_fermion_boundary_certified'] is False and
          d['production_boundary_permitted'] is False,
          'D2 output falsely certified')
    check(math.isclose(scalar_bubble_finite(1.e6,1.e6,1.e10),
                       math.log(1.e10/1.e6),rel_tol=1e-12),
          'A4 degenerate scalar bubble limit failed')
    # Independent source verification: derivative residues and source symmetry.
    with zipfile.ZipFile(soft_zip) as z:
        raw,_=extract(z,'D_-2')
    _,residue,finite,fd,fm=parse_d_soft(raw)
    vals={S['A']:math.sqrt(d['scalar_masses_squared_GeV2_at_muF'][0]),
          S['B']:math.sqrt(d['scalar_masses_squared_GeV2_at_muF'][1]),
          S['C']:math.sqrt(d['scalar_masses_squared_GeV2_at_muF'][2]),
          S['Q']:1.e10}
    fd2=real_sympy(fd,{**vals,S['Q']:vals[S['Q']]*math.e**2})
    fm2=real_sympy(fm,{**vals,S['Q']:vals[S['Q']]*math.e**2})
    check(abs(fd2-d['kernel_2x2_hbar1'][0][0]+4/math.sqrt(3)) < 2e-9,
          'D2 wrong logarithmic scale evolution')
    check(abs(fm2-d['kernel_2x2_hbar1'][0][1]) < 2e-10,
          'D2 mixed kernel changed with scale')
    # Negative test: tampered source must not pass SHA audit.
    with zipfile.ZipFile(sources_zip) as z:
        root='output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1/D_-2'
        ledger=json.loads(z.read(root+'/one_loop_operator_ledger.json'))
        raw=z.read(root+'/native_matching/data/eft1_after_F_one_loop_inputform.txt')
    try:
        restricted_hard(ledger,raw+b' ')
    except ValueError as e:
        check('hash' in str(e).lower(), 'Incorrect tamper error')
    else:
        raise AssertionError('Tampered Matchete source accepted')
    # Negative test: altering restricted triplet inventory must fail.
    tampered=copy.deepcopy(ledger)
    tampered['terms']=[t for t in tampered['terms'] if t.get('index') != 1013]
    try:
        restricted_hard(tampered,raw)
    except ValueError:
        pass
    else:
        raise AssertionError('Incomplete direct triplet accepted')
    # Negative test: nonzero independent conjugate Yukawa must be rejected.
    with zipfile.ZipFile(sources_zip) as z:
        root='output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1/D_-2'
        benchmark=json.loads(z.read(root+'/fermion_matching_input.json'))
    y=benchmark['native_uv_rgbeta_couplings_at_muF']['y2Conj']
    check(all(x['real']==x['imag']==0 for row in y for x in row),
          'Benchmark does not have restricted independent conjugate Yukawa')
    assert_restricted_yukawa_inventory(benchmark['native_uv_rgbeta_couplings_at_muF'])
    mutant=copy.deepcopy(benchmark['native_uv_rgbeta_couplings_at_muF'])
    mutant['y2Conj'][0][0]['real']=1e-3
    try:
        assert_restricted_yukawa_inventory(mutant)
    except ValueError as e:
        check('y2Conj' in str(e), 'Incorrect conjugate-Yukawa gate error')
    else:
        raise AssertionError('Nonzero independent y2Conj was accepted')
    print('PASS A4 baseline; D2 restricted 12->6 hard terms; D2 benchmark;')
    print('PASS exact symbolic soft kernel; mixed quartic WExtra=0; pole pairing;')
    print('PASS logarithmic scale derivative; source hash tamper; term removal;')
    print('PASS scalar mass limit; no physical Wilson boundary exported.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources-zip',type=Path,required=True)
    p.add_argument('--soft-zip',type=Path,required=True)
    a=p.parse_args()
    run(a.sources_zip,a.soft_zip)


if __name__=='__main__':main()
