"""Compare A(-4) ORIGINAL and five-monomial REFERENCE-DIAGNOSTIC EFT1 running.

Both trajectories start from exactly the same saved A(-4) mu_F numerical state.
This uses the repository's *existing* 17-coupling RG/ODE engine twice and
exports full 45-Weyl / 12-real-scalar tensors for downstream diagnostics.
Neither trajectory is physically certified or a Matchete scalar dictionary.
The alternative beta differs ONLY in the five scalar-only polynomial terms
isolated by RunAminus4QuarticFiveMonomialBatch.py.

Run from repository root, after the five-monomial script:
 python scripts/RunAminus4ConditionalQuarticODE.py \
    --input-root output/a_minus4_next_batch/quartic_five_monomial \
    --output-root output/a_minus4_next_batch/conditional_quartic_ode
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from Numerical.running import Aminus4EFT1Couplings as eft
from Numerical.running.Aminus4QuarticDiscrepancy import COORDINATES, coordinates_of
from scripts.ExportAminus4EFT1Couplings import make_state

BASIS = Path('output/a_minus4_eft1_batch2/A_m4_eft1_tensor_basis.json')
ORIGINAL_RAW = Path('output/a_minus4_eft1_batch2/A_m4_projected_raw.json')


def load(path:Path):
    if not path.is_file(): raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path:Path)->str:
    return sha256(path.read_bytes()).hexdigest()


def evaluate(root:Path,input_root:Path,output_root:Path,*,samples:int=7,
             rtol:float=1e-8,atol:float=1e-11)->dict:
    root=root.resolve(); input_root=input_root.resolve(); output_root=output_root.resolve()
    report=load(input_root/'five_monomial_report.json')
    if report.get('format')!='t3-a-minus4-five-monomial-quartic-reference-DIAGNOSTIC-v1' or \
       report.get('physical_quartic_beta_certified') is not False or \
       report.get('physical_EFT1_Wilson_boundary_certified') is not False:
        raise ValueError('Missing safely labelled five-monomial diagnostic audit')
    if report['input_sha256']['global_rgbeta_projection']!=sha(root/ORIGINAL_RAW):
        raise ValueError('Native global projected RGBeta source has changed')
    baseline=input_root/'native_global_projection_compatible.json'
    candidate=input_root/'scalar_only_reference_DIAGNOSTIC_beta.json'
    base=load(baseline); diag=load(candidate)
    eft.validate_beta_payload(base);eft.validate_beta_payload(diag)
    if diag['metadata'].get('Aminus4DoNotTreatAsProductionOrPhysical') is not True or \
       diag['metadata'].get('Aminus4PhysicalQuarticBetaCertified') is not False:
        raise ValueError('Candidate missing nonphysical diagnostic safety labels')
    if any(base['report_betas'][k]!=diag['report_betas'][k]
           for k in eft.BETA_NAMES - {'lambdaH','lambdaS2','lambdaS2Adj','lambdaH2','lambdaH2Adj'}):
        raise ValueError('Unexpected changes outside the five isolated scalar-beta channels')
    saved=load(root/BASIS)
    if saved.get('model')!='A_-4' or len(saved.get('running_couplings',[]))<3:
        raise ValueError('Saved A(-4) boundary is incomplete')
    muF=float(saved['sampled_scales_gev'][0]); muS=float(saved['sampled_scales_gev'][-1])
    if not (0<muS<muF):raise ValueError('Invalid A(-4) threshold order')
    seed={'model':'A_-4','MSbar_mu_F_gev':muF,
          'muF_renormalizable_state':saved['running_couplings'][0]}
    initial=make_state(seed)
    native=eft.integrate_eft1(root,initial,baseline,muS,n_samples=samples,
                               rtol=rtol,atol=atol)
    reference=eft.integrate_eft1(root,initial,candidate,muS,n_samples=samples,
                                  rtol=rtol,atol=atol)
    mu_n=np.asarray(native.mu_gev,float);mu_r=np.asarray(reference.mu_gev,float)
    if len(mu_n)!=samples or not np.allclose(mu_n,mu_r,rtol=0,atol=1e-9):
        raise RuntimeError('Compared EFT1 ODE grids disagree')
    for folder in (output_root,output_root/'native',output_root/'reference_diagnostic'):
        if folder.exists():raise FileExistsError('Refusing to overwrite conditional ODE output: '+str(folder))
    output_root.mkdir(parents=True,exist_ok=False)
    native_info=eft.write_sampled_npz(native,output_root/'native',source_root=root,
                                      source_beta=baseline)
    reference_info=eft.write_sampled_npz(reference,output_root/'reference_diagnostic',
                                         source_root=root,source_beta=candidate)
    rows=[]; max_abs={n:0. for n in COORDINATES}
    for k,mu in enumerate(mu_n):
        a=coordinates_of(native_info['running_couplings'][k])
        b=coordinates_of(reference_info['running_couplings'][k])
        change=b-a
        line={'mu_gev':float(mu)}
        for n,value in zip(COORDINATES,change):
            line['reference_minus_native_'+n]=float(value)
            max_abs[n]=max(max_abs[n],abs(float(value)))
        rows.append(line)
    with (output_root/'coupling_differences.csv').open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    endpoints={n:{'native':float(x),'reference_diagnostic':float(y),
                  'reference_minus_native':float(y-x)}
               for n,x,y in zip(COORDINATES,
                    coordinates_of(native_info['running_couplings'][-1]),
                    coordinates_of(reference_info['running_couplings'][-1]))}
    summary={
        'format':'t3-a-minus4-native-vs-reference-scalar-only-CONDITIONAL-ode-v1',
        'status':'TWO_CONDITIONAL_ODE_TRAJECTORIES_NUMERICALLY_COMPUTED_NOT_PHYSICAL',
        'model':'A_-4','muF_gev':muF,'muS_gev':muS,
        'n_samples':len(mu_n),
        'native_beta_sha256':sha(baseline),'reference_diagnostic_beta_sha256':sha(candidate),
        'reference_diagnostic_differs_from_native_only_in_five_scalar_betas':True,
        'endpoint_quartic_comparison':endpoints,
        'max_absolute_along_saved_grid_reference_minus_native':max_abs,
        'native_tensor_output':str(output_root/'native'/'A_m4_eft1_tensors.npz'),
        'reference_tensor_output':str(output_root/'reference_diagnostic'/'A_m4_eft1_tensors.npz'),
        'csv':str(output_root/'coupling_differences.csv'),
        'physically_certified_scalar_running':False,
        'physically_certified_Wilson_boundary':False,
        'physical_C5_generated':False,
        'no_source_or_production_beta_modified':True,
        'EFT_order':'UV -> integrate out ALL fermions -> EFT1(SM+S1+S2) -> integrate out ALL scalars together -> SMEFT',
    }
    (output_root/'conditional_ode_summary.json').write_text(
        json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8')
    return summary


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,default=ROOT)
    ap.add_argument('--input-root',type=Path,
                    default=Path('output/a_minus4_next_batch/quartic_five_monomial'))
    ap.add_argument('--output-root',type=Path,
                    default=Path('output/a_minus4_next_batch/conditional_quartic_ode'))
    ap.add_argument('--samples',type=int,default=7)
    ap.add_argument('--rtol',type=float,default=1e-8)
    ap.add_argument('--atol',type=float,default=1e-11)
    args=ap.parse_args(argv)
    report=evaluate(args.root,args.input_root,args.output_root,
                    samples=args.samples,rtol=args.rtol,atol=args.atol)
    print('Compared original and diagnostic reference A(-4) EFT1 running,',report['n_samples'],'samples')
    print('Endpoint reference-minus-original lambdaH2Adj:',
          report['endpoint_quartic_comparison']['lambdaH2Adj']['reference_minus_native'])
    print('Physical scalar/Wilson certification: NO')
    print('Outputs:', args.output_root.resolve())
    return 0


if __name__=='__main__':
    raise SystemExit(main())
