"""Source-checked A(-4) and D(-2) direct fermion-scale hard/soft diagnostics.

This deliberately does NOT implement or export a renormalised one-loop EFT1
Wilson boundary. The scalar-soft sources originate at the *scalar* threshold;
reevaluating their topology at mu_F is a candidate EFT1 diagram, not proof of
UV-minus-EFT1 subtraction, full counterterms, or GammaCC -> Weyl conversion.

Use --sources-zip and --soft-zip from the previously collected local outputs.
Only the independent UV Yukawas y1,y2 are allowed. The Matchete source complex
conjugates Bar[y1], Bar[y2] are NOT independent y1Conj/y2Conj couplings.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
import zipfile

import numpy as np
import sympy as sp

# Resolve helper modules next to this script; no editable global repository state.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from WorkerAFermionEFT1SoftA4 import (  # noqa: E402
    audit_hard_direct, complex_scalar, matrix, native_encode,
    read_member, scalar_bubble_finite, evaluate as evaluate_a4,
)
from WorkerASoftMixingKernel import (  # noqa: E402
    S, source_parse, extract,
)

MODELS = ('A_-4', 'D_-2')
RAW_BASE = 'output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1'
INPUT_BASE = 'output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1'
SCHEMA = 'worker-a-direct-threshold-both-v1-uncertified'
LOOP = 1 / (16 * math.pi**2)


def demand(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def decode_complex_matrix(value, label: str) -> np.ndarray:
    a = matrix(value)
    demand(a.shape == (3, 3) and np.isfinite(a).all(), label + ' must be finite 3x3')
    return a


def nearly_real(value, label: str) -> float:
    z = complex_scalar(value)
    demand(np.isfinite([z.real, z.imag]).all() and
           abs(z.imag) < 1e-9 * max(1.0, abs(z.real)),
           label + ' has unsupported phase')
    return float(z.real)


def real_sympy(expr: sp.Expr, numbers: dict[sp.Symbol, float]) -> float:
    remaining = expr.free_symbols - set(numbers)
    demand(not remaining, f'Unsubstituted symbolic variables: {remaining}')
    val = complex(sp.N(expr.subs(numbers), 45))
    demand(np.isfinite([val.real, val.imag]).all() and
           abs(val.imag) < 1e-10 * max(1., abs(val.real)),
           'Nonreal/nonfinite kernel')
    return val.real


def restricted_hard(ledger: dict, raw: bytes) -> dict:
    """Audit 6 allowed D direct terms and exact triplets of each chirality."""
    demand(ledger.get('source_sha256') == sha256(raw).hexdigest(),
           'Native Matchete source hash mismatch')
    demand(ledger.get('exact_reconstruction_by_expansion') is True,
           'Ledger not reconstructed exactly')
    direct = [t for t in ledger['terms']
              if t.get('sector') == 'CandidateLeptonNumberViolatingLLHH']
    demand(len(direct) == 12, f'Expected 12 D direct LLHH, got {len(direct)}')
    kept = [t for t in direct if 'Coupling[y2Conj,' not in t['raw_inputform']]
    discarded = [t for t in direct if 'Coupling[y2Conj,' in t['raw_inputform']]
    demand(len(kept) == len(discarded) == 6,
           'Unexpected restricted D direct LLHH inventory')
    # Reuse A's strict six-term source signature / exact Laurent triplets.
    audit = audit_hard_direct({**ledger, 'terms': kept}, raw)
    return {'source_term_count': len(ledger['terms']),
            'direct_llhh_total': len(direct),
            'direct_llhh_restricted': len(kept),
            'direct_llhh_excluded_y2Conj': len(discarded),
            'direct_hard_source_indices': audit,
            'excluded_source_indices': [t['index'] for t in discarded]}


def assert_restricted_yukawa_inventory(params: dict) -> None:
    """Reject independent y1Conj/y2Conj while allowing Bar[y1], Bar[y2]."""
    demand('y1Conj' not in params,
           'Independent y1Conj UV coupling is prohibited')
    if 'y2Conj' in params:
        y2conj = decode_complex_matrix(params['y2Conj'], 'y2Conj')
        demand(np.count_nonzero(y2conj) == 0,
               'Independent y2Conj is nonzero; restricted model does not apply')


def parse_d_soft(data: bytes) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Derive 2x2 operator kernel from native D(-2) scalar-soft expression."""
    expr = source_parse(data)
    eps, W, WE, L, LE, Q = (S[s] for s in ('eps', 'W', 'WE', 'L', 'LE', 'Q'))
    residue = sp.simplify(sp.limit(eps * expr, eps, 0))
    expected = -2 * (L*W + LE*WE) / sp.sqrt(3)
    demand(sp.simplify(residue - expected) == 0,
           'D source soft residue does not reconstruct')
    finite = sp.simplify(sp.limit(sp.cancel(expr - residue/eps), eps, 0))
    fd = sp.factor(sp.diff(sp.diff(finite, L), W))
    fm = sp.factor(sp.diff(sp.diff(finite, LE), W))
    reconstructed = fd*(L*W + LE*WE) + fm*(LE*W + L*WE)
    demand(sp.simplify(finite - reconstructed) == 0,
           'D mixed soft finite expression does not reconstruct')
    demand(sp.simplify(Q * sp.diff(fd, Q) + 2/sp.sqrt(3)) == 0,
           'D diagonal logarithmic derivative not tied to pole')
    demand(sp.simplify(Q * sp.diff(fm, Q)) == 0,
           'D mixed logarithmic derivative changed')
    return expr, residue, finite, fd, fm


def evaluate_d(sources_zip: Path, soft_zip: Path, mu_f: float) -> dict:
    demand(mu_f > 0 and np.isfinite(mu_f), 'Invalid fermion matching scale')
    root = RAW_BASE + '/D_-2'
    root_in = INPUT_BASE + '/D_-2'
    with zipfile.ZipFile(sources_zip) as z:
        source = read_member(z, root + '/native_matching/data/eft1_after_F_one_loop_inputform.txt')
        ledger_data = read_member(z, root + '/one_loop_operator_ledger.json')
        input_data = read_member(z, root_in + '/fermion_matching_input.json')
    with zipfile.ZipFile(soft_zip) as z:
        soft, soft_member = extract(z, 'D_-2')  # Checks identical duplicate exports.
    f = json.loads(input_data)
    ledger = json.loads(ledger_data)
    hard_source = restricted_hard(ledger, source)
    demand(f.get('model') == 'D_-2' and f.get('scale_GeV') == mu_f,
           'D benchmark or scale mismatch')
    demand(f.get('all_heavy_fermions_to_integrate_out_together') is True and
           f.get('both_heavy_scalar_multiplets_remain_light') is True,
           'D EFT1 must retain both scalars and integrate all fermions')
    params = f['native_uv_rgbeta_couplings_at_muF']
    assert_restricted_yukawa_inventory(params)
    mf = decode_complex_matrix(params['MF'], 'MF')
    y1 = decode_complex_matrix(params['y1'], 'y1')
    y2 = decode_complex_matrix(params['y2'], 'y2')
    offdiag = float(np.linalg.norm(mf - np.diag(np.diag(mf))) / np.linalg.norm(mf))
    demand(offdiag < 1e-8, 'General heavy chiral mass-basis rotation not certified')
    eig = np.diag(mf)
    demand(np.min(eig.real) > 0 and
           np.max(abs(eig.imag)) <= 1e-8 * max(1., np.max(abs(eig.real))),
           'Heavy eigenmasses must be real positive')
    masses = eig.real
    a = nearly_real(params['mS1Sq'], 'mS1Sq')
    m2 = nearly_real(params['mS2Sq'], 'mS2Sq')
    mu2 = nearly_real(params['muT3GeneralS2Sq'], 'muT3GeneralS2Sq')
    b, c = m2 - 2*abs(mu2), m2 + 2*abs(mu2)
    demand(min(a,b,c)>0, 'Nonpositive scalar mass squared')
    lam = nearly_real(params['lambdaT3'], 'lambdaT3')
    extra = nearly_real(params['lambdaT3GeneralQ020101I1'], 'extra quartic')
    full, pole, finite, fd, fm = parse_d_soft(soft)
    vals = {S['A']: math.sqrt(a), S['B']: math.sqrt(b),
            S['C']: math.sqrt(c), S['Q']: mu_f**2}
    diag = real_sympy(fd, vals)
    mix = real_sympy(fm, vals)
    # The native scalar source has two formal Wilson coordinates W12, WExtra.
    # The restricted independent Yukawa setup sets WExtra = 0 for this benchmark;
    # f_mix * conjugate(lambdaExtra) * W12 nevertheless remains.
    row_w = diag*lam + mix*extra
    row_extra = mix*lam + diag*extra
    demand(np.isfinite([diag,mix,row_w,row_extra]).all(), 'Bad numerical kernel')
    for w, we in [(1.,0.), (0.,1.), (.37,-.29)]:
        lhs = real_sympy(finite, {**vals,S['L']:lam,S['LE']:extra,
                                  S['W']:w,S['WE']:we})
        demand(np.isclose(lhs, row_w*w+row_extra*we, rtol=5e-10, atol=1e-11),
               'Direct Matchete expression disagrees with numerical 2x2 kernel')
    wbar = (y1.conj()/masses[np.newaxis, :]) @ y2.conj().T
    # Physical loop factor is assumed, not independently GammaCC-normalised.
    wbar_extra = np.zeros_like(wbar)
    prefactor = 2/(math.sqrt(3)*16*math.pi**2)
    hard_pole = prefactor*lam*wbar
    # Candidate scalar-soft Laurent source, using same ordered native Wbar.
    soft_pole = -prefactor*(lam*wbar + extra*wbar_extra)
    soft_fin = LOOP*(row_w*wbar + row_extra*wbar_extra)
    soft_fin_diag = LOOP*diag*lam*wbar
    soft_fin_mix = LOOP*mix*extra*wbar
    hard_logs = 1+np.log(mu_f**2/masses**2)
    hard_fin = prefactor*lam*((y1.conj()*hard_logs[np.newaxis,:]/masses[np.newaxis,:]) @ y2.conj().T)
    paired = hard_fin + soft_fin
    pole_residual = float(np.linalg.norm(hard_pole+soft_pole))
    demand(pole_residual <= 1e-13*max(float(np.linalg.norm(hard_pole)),1e-40),
           'Restricted D direct pole does not pair in ordered basis')
    # At fixed masses and couplings, fd's derivative is -2/sqrt(3),
    # cancelling derivative of the hard log +2/sqrt(3).
    dlogq_hard = prefactor*lam*wbar
    dlogq_soft = LOOP*(-2/math.sqrt(3))*lam*wbar
    scale_residual = float(np.linalg.norm(dlogq_hard+dlogq_soft))
    demand(scale_residual <= 1e-13*max(float(np.linalg.norm(dlogq_hard)),1e-40),
           'Explicit log-scale cancellation check failed')
    return {
        'schema':SCHEMA, 'model':'D_-2',
        'status':'D2_RESTRICTED_DIRECT_HARD_SOFT_FORMAL_DIAGNOSTIC_NOT_MATCHED',
        'mu_F_GeV':mu_f,
        'integrated_fields':['F'], 'retained_scalar_multiplets':['S1','S2'],
        'original_source_sha256':{
            'uv_one_loop':sha256(source).hexdigest(),
            'operator_ledger':sha256(ledger_data).hexdigest(),
            'fermion_benchmark':sha256(input_data).hexdigest(),
            'later_scalar_soft':sha256(soft).hexdigest()},
        'later_scalar_soft_source_path':soft_member,
        'hard_llhh_restriction':hard_source,
        'all_three_fermion_masses_GeV':masses.tolist(),
        'fermion_mass_offdiagonal_relative_norm':offdiag,
        'scalar_masses_squared_GeV2_at_muF':[a,b,c],
        'scalar_bilinear_GeV2':mu2,
        'quartic_lambdaT3':lam,
        'quartic_lambdaExtraQ020101I1':extra,
        'soft_pole_reduced_sympy':sp.sstr(pole),
        'soft_diagonal_finite_kernel_sympy':sp.sstr(fd),
        'soft_mixed_finite_kernel_sympy':sp.sstr(fm),
        'kernel_2x2_hbar1':[[diag,mix],[mix,diag]],
        'restricted_W12_coefficient_hbar1':row_w,
        'extra_quartic_W12_mixed_coefficient_hbar1':mix*extra,
        'restricted_WExtra_is_zero':True,
        'ordered_Wbar12_GeV_inv':native_encode(wbar),
        'ordered_WbarExtra_GeV_inv':native_encode(wbar_extra),
        'uv_direct_pole_GeV_inv':native_encode(hard_pole),
        'candidate_eft1_soft_pole_GeV_inv':native_encode(soft_pole),
        'uv_direct_finite_GeV_inv':native_encode(hard_fin),
        'candidate_eft1_soft_finite_GeV_inv':native_encode(soft_fin),
        'candidate_eft1_soft_finite_diagonal_GeV_inv':native_encode(soft_fin_diag),
        'candidate_eft1_soft_finite_mixed_GeV_inv':native_encode(soft_fin_mix),
        'paired_direct_finite_GeV_inv':native_encode(paired),
        'uv_direct_pole_frobenius_GeV_inv':float(np.linalg.norm(hard_pole)),
        'candidate_eft1_soft_finite_frobenius_GeV_inv':float(np.linalg.norm(soft_fin)),
        'candidate_eft1_soft_finite_mixed_frobenius_GeV_inv':float(np.linalg.norm(soft_fin_mix)),
        'paired_direct_finite_frobenius_GeV_inv':float(np.linalg.norm(paired)),
        'direct_ordered_pole_residual_GeV_inv':pole_residual,
        'explicit_scale_derivative_residual_GeV_inv':scale_residual,
        'symbolic_source_reconstruction_exact':True,
        'hbar_identification_assumed':'1/(16*pi*pi)',
        'scalar_mass_basis_mapping_certified':False,
        'full_eft1_soft_amplitude_certified':False,
        'counterterms_verified':False,
        'full_operator_mixing_verified':False,
        'gammaCC_to_weyl_normalization_verified':False,
        'finite_one_loop_fermion_boundary_certified':False,
        'production_boundary_permitted':False,
    }


def evaluate_both(sources_zip: Path, soft_zip: Path, mu_f: float = 1e5) -> dict:
    a = evaluate_a4(sources_zip,soft_zip,mu_f)
    d = evaluate_d(sources_zip,soft_zip,mu_f)
    return {'schema':SCHEMA,
            'status':'SOURCE_CHECKED_DIRECT_CHANNELS_ONLY_NO_FINITE_WILSON_MATCHING',
            'input_archives_sha256':{'sources':sha256(sources_zip.read_bytes()).hexdigest(),
                                     'soft':sha256(soft_zip.read_bytes()).hexdigest()},
            'matching_scale_GeV':mu_f,
            'results':{'A_-4':a,'D_-2':d},
            'production_boundary_permitted':False}


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--sources-zip',type=Path,required=True)
    ap.add_argument('--soft-zip',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--mu-f',type=float,default=1.e5)
    args=ap.parse_args()
    result=evaluate_both(args.sources_zip,args.soft_zip,args.mu_f)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    a=result['results']['A_-4'];d=result['results']['D_-2']
    print('A_-4: direct hard pole {:.9e}, paired diagnostic {:.9e} GeV^-1'.format(
        a['hard_direct_pole_frobenius_GeV_inv'],a['paired_direct_finite_frobenius_GeV_inv']))
    print('D_-2: restricted hard pole {:.9e}, paired diagnostic {:.9e} GeV^-1'.format(
        d['uv_direct_pole_frobenius_GeV_inv'],d['paired_direct_finite_frobenius_GeV_inv']))
    print('D_-2: mixed kernel {:.12g}, induced W12 contribution norm {:.9e} GeV^-1'.format(
        d['kernel_2x2_hbar1'][0][1],d['candidate_eft1_soft_finite_mixed_frobenius_GeV_inv']))
    print('No renormalised Wilson coefficient was produced.')
    print('Wrote',args.output)


if __name__=='__main__':main()
