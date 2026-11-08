"""Export the full sampled A(-4) EFT1 running tensors (RGBeta basis only).

Usage: python scripts/ExportAminus4EFT1Couplings.py --initial-file path.json
       --beta-file path.json --mu-s 1000 --output-dir path --samples 7

Input JSON is a native T3IntermediateScalarState-compatible dictionary:
{ "mu_gev":100000, "sm":{"gY":...,"g2":...,"g3":...,
  "lambdaH":...,"yu":3x3,"yd":3x3,"ye":3x3},
  "mS1Sq":...,"mS2Sq":...,"lambdaS1":...,"lambdaS2":...,
  "lambdaH1":...,"lambdaH2":...,"lambda12":...,
  "lambdaT3":{"real":...,"imag":...},
  "lambdaH2Adj":...,"lambdaS2Adj":...}
Complex SM 3x3 matrices may be {"real": matrix, "imag": matrix}.
Also accepts the raw A(-4) report
  output/general_observable/lo_certification_couplings/ws2_all16_v2/
  A_m4/running_couplings_muF_muS.json
using its muF_renormalizable_state as input, WITHOUT accepting any stored
source/status flag as a physical certification. The stored muS entries are
not used to replace solving the equations.
All values MUST be supplied at the fermion threshold (not blindly reused
from the UV fitting scale). No couplings or threshold masses are invented.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from Numerical.running import Aminus4EFT1Couplings as eft


def _matrix(v):return eft._complex_array(v,(3,3),'SM Yukawa')


def make_state(doc):
    # Also accept the repository's raw, 17-coupling boundary-data file;
    # numerical values are INPUTS, never physical-basis certificates.
    if 'muF_renormalizable_state' in doc:
        if doc.get('model')!='A_-4':
            raise eft.EFT1BasisError('Raw running data does not belong to A(-4)')
        raw=doc['muF_renormalizable_state']
        if set(raw)!=eft.BETA_NAMES:
            raise eft.EFT1BasisError('Raw muF input does not have 17 EFT1 beta couplings')
        sm_names={'gY','g2','g3','lambdaH','yu','yd','ye'}
        doc={'mu_gev':doc['MSbar_mu_F_gev'],
             'sm':{k:raw[k] for k in sm_names},
             **{k:raw[k] for k in set(raw)-sm_names}}
    from Numerical.core.State import SMNumericalState,T3Representation
    from Numerical.running.IntermediateScalarState import T3IntermediateScalarState
    expected={
        'mu_gev','sm','mS1Sq','mS2Sq','lambdaS1','lambdaS2',
        'lambdaH1','lambdaH2','lambda12','lambdaT3','lambdaH2Adj','lambdaS2Adj'}
    if set(doc)!=expected:
        raise eft.EFT1BasisError(f'Initial state keys must be {sorted(expected)}')
    sm=doc['sm']
    if set(sm)!={'gY','g2','g3','lambdaH','yu','yd','ye'}:
        raise eft.EFT1BasisError('Incomplete SM 3-generation source at mu_F')
    state=T3IntermediateScalarState(
        mu_gev=float(doc['mu_gev']),
        representation=T3Representation(1,3,2,-4,False),
        sm=SMNumericalState(gY=float(sm['gY']),g2=float(sm['g2']),
            g3=float(sm['g3']),lambdaH=float(sm['lambdaH']),
            yu=_matrix(sm['yu']),yd=_matrix(sm['yd']),ye=_matrix(sm['ye'])),
        mS1Sq=float(doc['mS1Sq']),mS2Sq=float(doc['mS2Sq']),
        lambdaS1=float(doc['lambdaS1']),lambdaS2=float(doc['lambdaS2']),
        lambdaH1=float(doc['lambdaH1']),lambdaH2=float(doc['lambdaH2']),
        lambda12=float(doc['lambda12']),
        lambdaT3=eft._complex_scalar(doc['lambdaT3'],'lambdaT3'),
        lambdaH2Adj=float(doc['lambdaH2Adj']),
        lambdaS2Adj=float(doc['lambdaS2Adj']))
    state=state.validated()
    eft.state_point(state)
    return state


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--initial-file',type=Path,required=True)
    ap.add_argument('--beta-file',type=Path,required=True)
    ap.add_argument('--mu-s',type=float,required=True)
    ap.add_argument('--samples',type=int,default=7)
    ap.add_argument('--output-dir',type=Path,
                    default=ROOT/'output/a_minus4_eft1_batch2')
    args=ap.parse_args(argv)
    initial=json.loads(args.initial_file.read_text(encoding='utf-8-sig'))
    state=make_state(initial)
    if not 0<args.mu_s<state.mu_gev:
        raise eft.EFT1BasisError('Require 0 < mu_S < mu_F')
    output=eft.integrate_eft1(ROOT,state,args.beta_file,args.mu_s,
                              n_samples=args.samples)
    info=eft.write_sampled_npz(output,args.output_dir,source_root=ROOT,
                               source_beta=args.beta_file)
    print('A(-4) RGBeta EFT1 tensor coordinates:',info['sampled_scales_gev'])
    print('Wrote:',args.output_dir/'A_m4_eft1_tensors.npz')
    print('Wrote:',args.output_dir/'A_m4_eft1_tensor_basis.json')
    print('Matchete physical basis map certified: NO')
    return 0

if __name__=='__main__':raise SystemExit(main())
