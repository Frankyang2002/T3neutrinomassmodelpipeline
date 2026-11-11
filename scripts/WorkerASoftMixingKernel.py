"""Source-exact numerical scalar-soft mixing kernel for A_-4 / D_-2.

Reads original Matchete joint-scalar-threshold expressions and the saved
RGBeta scalar endpoint from worker_a_sources_*.zip. This IS NOT the one-loop
fermion-threshold UV-EFT1 matching, and its matrices MUST NOT be used as the
fermion matching boundary. The scalar mass-to-Matchete and quartic dictionary
are provisional, explicitly recorded in output.

Only stdlib + numpy + sympy required. Never changes production outputs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import zipfile
import numpy as np
import sympy as sp

NAMES = ('L', 'LE', 'W', 'WE', 'A', 'B', 'C', 'Q', 'eps')
S = {k: sp.Symbol(k) for k in NAMES}
REPLACEMENTS = {
    'Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]]': 'LE',
    'Bar[Coupling[T3EFT1WExtra, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]': 'WE',
    'Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]': 'W',
    'Coupling[lambdaT3, {}, 0]': 'L',
    'Coupling[MS1, {}, 0]': 'A',
    'Coupling[MS2, {}, 0]': 'B',
    'Coupling[T32MinusMass, {}, 0]': 'B',
    'Coupling[T32PlusMass, {}, 0]': 'C',
    'μbar2': 'Q', 'ϵ': 'eps', 'hbar': '1', 'Sqrt[3]': 'sqrt(3)',
}


def source_parse(data: bytes) -> sp.Expr:
    text = data.decode('utf-8-sig').strip()
    if 'T3EFT1W12' not in text:
        raise ValueError('Not a native W12 scalar-soft Matchete export')
    for old, new in REPLACEMENTS.items():
        text = text.replace(old, new)
    text = text.replace('Log[', 'log(').replace('[', '(').replace(']', ')').replace('^', '**')
    if not re.fullmatch(r'[A-Za-z0-9_+*/().\s-]+', text):
        raise ValueError('Unparsed source token remains')
    return sp.sympify(text, locals={**S, 'sqrt':sp.sqrt, 'log':sp.log})


def numeric(expr: sp.Expr, values: dict) -> float:
    x = complex(sp.N(expr.subs(values), 16))
    if abs(x.imag) > 1e-10 * max(1., abs(x.real)) or not np.isfinite([x.real, x.imag]).all():
        raise ValueError('Nonreal or nonfinite evaluated coefficient')
    return x.real


def real_coupling(x, label: str) -> float:
    if isinstance(x, dict) and set(x) >= {'real', 'imag'}:
        z = complex(x['real'], x['imag'])
    else:
        z = complex(x)
    if not np.isfinite([z.real,z.imag]).all() or abs(z.imag)>1e-9*max(1., abs(z.real)):
        raise ValueError(label + ': complex coupling outside checked real benchmark')
    return float(z.real)


def read_inputs(sources: zipfile.ZipFile, model: str) -> tuple[dict, str, bytes]:
    suffix = '/'+model+'/scalar_matching_input.json'
    options = [n for n in sources.namelist() if n.endswith(suffix) and '/full_general_batch1/' in n]
    if len(options)!=1: raise ValueError(f'Expected exactly one {model} scalar endpoint, got {options}')
    path=options[0]; payload=sources.read(path)
    data=json.loads(payload)
    if data.get('model')!=model or data.get('scale_GeV')!=1000.0 or data.get('both_heavy_scalar_multiplets_to_integrate_out_together') is not True:
        raise ValueError('Scalar endpoint metadata mismatch')
    return data,path,payload


def extract(soft: zipfile.ZipFile, model: str) -> tuple[bytes,str]:
    prefix = 'A_scalar_' if model=='A_-4' else 'D_scalar_'
    names=sorted(n for n in soft.namelist() if f'/{model}/' in n and n.endswith(prefix+'soft_only_unrenormalized.txt'))
    if len(names)!=(1 if model=='A_-4' else 2):raise ValueError(f'Unexpected {model} source count: {names}')
    payloads=[soft.read(n) for n in names]
    if any(p != payloads[0] for p in payloads):raise ValueError('Duplicate D_-2 soft sources disagree')
    return payloads[0], names[0]


def build(model: str, soft: zipfile.ZipFile, sources: zipfile.ZipFile) -> dict:
    raw, rawpath=extract(soft,model)
    inp, inputpath,inputbytes=read_inputs(sources,model)
    expr=source_parse(raw)
    eps=S['eps']; W=S['W'];WE=S['WE'];L=S['L'];LE=S['LE']
    residue=sp.simplify(sp.limit(eps*expr,eps,0))
    finite=sp.cancel(expr-residue/eps)
    finite=sp.simplify(sp.limit(finite,eps,0))
    expected=-2/sp.sqrt(3)*(L*W if model=='A_-4' else L*W+LE*WE)
    if sp.simplify(residue-expected)!=0:raise ValueError(model+': source pole structure changed')
    fdiag=sp.factor(sp.diff(sp.diff(finite,L),W))
    fmix=sp.S(0) if model=='A_-4' else sp.factor(sp.diff(sp.diff(finite,LE),W))
    reconstructed=fdiag*L*W if model=='A_-4' else fdiag*(L*W+LE*WE)+fmix*(LE*W+L*WE)
    if sp.simplify(finite-reconstructed)!=0:raise ValueError('Finite source decomposition did not reconstruct')
    # Q*d/dQ checks the soft logarithmic scale dependence against the UV pole.
    # The mixed operator kernel is scale-independent in the native source.
    scale_diag = sp.simplify(S['Q']*sp.diff(fdiag,S['Q']))
    scale_mix = sp.simplify(S['Q']*sp.diff(fmix,S['Q']))
    if sp.simplify(scale_diag+2/sp.sqrt(3)) != 0 or sp.simplify(scale_mix) != 0:
        raise ValueError('Unexpected logarithmic scale derivative / pole relation')
    if model=='D_-2':
        if sp.simplify(sp.diff(sp.diff(finite,LE),WE)-fdiag)!=0:raise ValueError('D diagonal asymmetry')
        if sp.simplify(sp.diff(sp.diff(finite,L),WE)-fmix)!=0:raise ValueError('D offdiagonal asymmetry')
    sq=inp['scalar_quadratic_terms_GeV2']
    m1=float(sq['mS1Sq']);m2=float(sq['mS2Sq'])
    if m1<=0 or m2<=0: raise ValueError('Nonpositive scalar mass squared')
    values={S['A']:math.sqrt(m1),S['Q']:1.e6}
    if model=='A_-4':
        values[S['B']]=math.sqrt(m2)
        mass_sq=[m1,m2]
    else:
        # Native holomorphic mass convention: m^2_{2,+/-}=mS2Sq +/- 2|mu2|.
        mu=real_coupling(inp['native_eft1_rgbeta_couplings_at_muS']['muT3GeneralS2Sq'],'muT3GeneralS2Sq')
        minussq=m2-2*abs(mu); plussq=m2+2*abs(mu)
        if minussq<=0:raise ValueError('Tachyonic split scalar')
        values[S['B']]=math.sqrt(minussq);values[S['C']]=math.sqrt(plussq)
        mass_sq=[m1,minussq,plussq]
    diagonal=numeric(fdiag,values)
    mixed=numeric(fmix,values)
    base=real_coupling(inp['native_eft1_rgbeta_couplings_at_muS']['lambdaT3'],'lambdaT3')
    extra=0.0 if model=='A_-4' else real_coupling(inp['general_quartics']['lambdaT3GeneralQ020101I1'],'lambdaT3GeneralQ020101I1')
    # At source level finite = W*(diag*L+mixed*LE) + WE*(mixed*L+diag*LE).
    # This is the coefficient row of the formal scalar-soft functional in W,WE.
    row=[diagonal*base+mixed*extra, mixed*base+diagonal*extra]
    coefficients=np.asarray([[diagonal,mixed],[mixed,diagonal]],dtype=float)
    if not np.isfinite(coefficients).all() or not np.isfinite(row).all():raise ValueError('Bad numerical mixing kernel')
    # Compare direct expression evaluated at unit W and unit WE against matrix contraction.
    for w,we in [(1.,0.),(0.,1.),(.37,-.29)]:
        lhs=numeric(finite,{**values,L:base,LE:extra,W:w,WE:we})
        rhs=w*row[0]+we*row[1]
        if not np.isclose(lhs,rhs,rtol=2e-9,atol=1e-12):raise ValueError('Direct/source matrix reconstruction failed')
    return {
        'model':model, 'status':'ScalarThresholdRawFiniteMixingKernel_SourceVerified_NotFermionMatching',
        'soft_source_path':rawpath, 'soft_source_sha256':hashlib.sha256(raw).hexdigest(),
        'scalar_endpoint_path':inputpath,'scalar_endpoint_sha256':hashlib.sha256(inputbytes).hexdigest(),
        'scale_gev':1000.0,'formal_mu_bar_gev':1000.0,
        'scalar_mass_squared_gev2':mass_sq, 'quartics':{'lambdaT3':base,'lambdaT3ExtraQ020101I1':extra},
        'source_pole_residue':sp.sstr(residue),
        'source_finite_diagonal_kernel':sp.sstr(fdiag),
        'source_finite_mixed_kernel':sp.sstr(fmix),
        'kernel_matrix_hbar1':coefficients.tolist(),
        'finite_row_coefficients_for_conjugated_W_WExtra_hbar1':row,
        'finite_coefficient_for_WExtra_zero_hbar1':row[0],
        'mixed_quartic_induced_contribution_with_WExtra_zero_hbar1':mixed*extra,
        'source_expression_reconstruction_verified':True,
        'd_diagonal_d_log_mubar2':sp.sstr(scale_diag),
        'd_mixed_d_log_mubar2':sp.sstr(scale_mix),
        'source_scale_derivative_pole_relation_verified':True,
        'physical_scalar_mass_mapping_certified':False,
        'matchete_rgbeta_quartic_dictionary_certified':False,
        'source_hbar_to_physical_loop_normalization_certified':False,
        'fermion_threshold_eft1_soft_subtraction_done':False,
        'fermion_threshold_counterterms_done':False,
        'gammaCC_weyl_normalization_certified':False,
        'production_boundary_permitted':False,
        'note':'Numerical scalar-threshold kernel diagnostic only. Do not insert into EFT1 at mu_F.',
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--soft-zip',type=Path,required=True)
    ap.add_argument('--sources-zip',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    with zipfile.ZipFile(args.soft_zip) as soft,zipfile.ZipFile(args.sources_zip) as src:
        results=[build(m,soft,src) for m in ('A_-4','D_-2')]
    out={'status':'NumericalRawScalarSoftKernelAudited_NotRenormalisedFermionMatching',
         'soft_archive_sha256':hashlib.sha256(args.soft_zip.read_bytes()).hexdigest(),
         'source_archive_sha256':hashlib.sha256(args.sources_zip.read_bytes()).hexdigest(),
         'results':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf8')
    for r in results:
        print(f"{r['model']}: f_diag={r['kernel_matrix_hbar1'][0][0]:.12g}, f_mix={r['kernel_matrix_hbar1'][0][1]:.12g}, W coefficient={r['finite_coefficient_for_WExtra_zero_hbar1']:.12g}")
    print('Wrote',args.output.resolve())

if __name__=='__main__':main()
