"""A(-4) direct LLHH: source-verified EFT1 scalar-bubble candidate at mu_F.

No production Wilson coefficient is returned.  All calculations are in the
ordered Matchete GammaCC coefficient convention, before Weyl/CG conversion.
Only UV F-hard direct LLHH and the tree-W12 scalar-loop topology are included.
No general EFT1 one-loop operator mixing or MSbar counterterms are supplied.

Usage:
  python scripts/WorkerAFermionEFT1SoftA4.py --sources-zip ... --soft-zip ... --output ...
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path, PurePosixPath
import re
import zipfile

import numpy as np
import sympy as sp

MODEL = 'A_-4'
RAW_ROOT = ('output/parallel_t3/full_general_batch2/'
            'thesis_full_general_batch2_raw_v1/A_-4')
BATCH1_ROOT = ('output/parallel_t3/full_general_batch1/'
               'thesis_full_general_batch1_v1/A_-4')
SOFT_MEMBER = ('output/parallel_t3/full_general_final_one_batch/'
               'thesis_a4_fresh_f_v3/A_-4/A_scalar_soft_only_unrenormalized.txt')
LEDGER = RAW_ROOT + '/one_loop_operator_ledger.json'
RAW = RAW_ROOT + '/native_matching/data/eft1_after_F_one_loop_inputform.txt'
F_INPUT = BATCH1_ROOT + '/fermion_matching_input.json'
SCHEMA = 'worker-a-a4-eFT1-soft-candidate-v1'


def demand(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def read_member(z: zipfile.ZipFile, name: str) -> bytes:
    names = z.namelist()
    demand(name in names, 'Missing exact source: ' + name)
    demand(not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts,
           'Unsafe source path')
    return z.read(name)


def complex_scalar(x):
    if isinstance(x, dict) and set(x) == {'real', 'imag'}:
        return complex(x['real'], x['imag'])
    return complex(x)


def matrix(x):
    return np.array([[complex_scalar(c) for c in row] for row in x], complex)


def native_encode(x):
    x = np.asarray(x, complex)
    return {'real': x.real.tolist(), 'imag': x.imag.tolist()}


def scalar_bubble_finite(a: float, b: float, mu_sq: float) -> float:
    """1+[a log(mu²/a)-b log(mu²/b)]/(a-b), stable near a=b."""
    demand(min(a, b, mu_sq) > 0 and np.isfinite([a, b, mu_sq]).all(),
           'Nonpositive or nonfinite squared mass / scale')
    t = (a - b) / b
    if abs(t) < 1e-5:
        # a ln(a/b)/(a-b) = (1+t) ln(1+t)/t.
        mass_diff = 1 + t * (.5 + t * (-1/6 + t * (1/12 + t * (-1/20 + t/30))))
    else:
        mass_diff = (1 + t) * math.log1p(t) / t
    return 1 + math.log(mu_sq / b) - mass_diff


def audit_soft_source(raw: str) -> dict:
    """Exact Laurent extraction of the native scalar-soft InputForm bubble."""
    demand('T3EFT1W12' in raw and 'lambdaT3' in raw, 'Wrong soft operator')
    replacement = re.compile(
        r'Bar\[Coupling\[T3EFT1W12, \{Index\[[^,\]]+, Flavor\], '
        r'Index\[[^,\]]+, Flavor\]\}, 0\]\]')
    s, n = replacement.subn('W', raw.strip().lstrip('\ufeff'))
    demand(n >= 2, 'Unexpected soft Wilson flavour syntax')
    for old, new in (
        ('Coupling[lambdaT3, {}, 0]', 'L'),
        ('Coupling[MS1, {}, 0]', 'M1'),
        ('Coupling[MS2, {}, 0]', 'M2'),
        ('μbar2', 'Q'), ('ϵ', 'eps'), ('ε', 'eps'), ('hbar', 'HB'),
    ):
        s = s.replace(old, new)
    demand(not re.search(r'(Coupling\[|Bar\[|Field\[|CG\[|Index\[)', s),
           'Unrecognised operator in soft source')
    s = s.replace('Sqrt[3]', 'sqrt(3)').replace('Log[', 'log(')
    s = s.replace('[', '(').replace(']', ')').replace('^', '**')
    demand(re.fullmatch(r'[A-Za-z0-9_+*/().\s-]+', s) is not None,
           'Unsafe soft expression characters')
    W, L, M1, M2, Q, eps, HB = sp.symbols('W L M1 M2 Q eps HB', nonzero=True)
    symbols = {str(v): v for v in (W, L, M1, M2, Q, eps, HB)}
    expression = sp.cancel(sp.sympify(s, locals={**symbols, 'sqrt': sp.sqrt, 'log': sp.log}).subs(HB, 1))
    residue = sp.simplify(sp.limit(eps * expression, eps, 0))
    finite = sp.simplify(sp.limit(expression - residue / eps, eps, 0))
    wanted_residue = -2 * L * W / sp.sqrt(3)
    wanted_finite = wanted_residue * (
        1 + (M1**2 * sp.log(Q / M1**2) - M2**2 * sp.log(Q / M2**2)) /
        (M1**2 - M2**2))
    demand(sp.simplify(residue - wanted_residue) == 0,
           'Scalar bubble pole does not match source')
    demand(sp.simplify(finite - wanted_finite) == 0,
           'Scalar bubble finite expression does not match source')
    return {'soft_pole_reduced': sp.sstr(residue),
            'soft_finite_reduced': sp.sstr(finite),
            'symbolic_reconstruction_exact': True}


def audit_hard_direct(ledger: dict, raw: bytes) -> dict:
    demand(ledger.get('source_sha256') == sha256(raw).hexdigest(),
           'Native raw source SHA mismatch')
    demand(ledger.get('exact_reconstruction_by_expansion') is True,
           'Raw ledger not fully reconstructed')
    rows = [t for t in ledger['terms']
            if t['sector'] == 'CandidateLeptonNumberViolatingLLHH']
    demand(len(rows) == 6, 'Expected six direct LLHH structures')
    out = {}
    canonical = {}
    for x in rows:
        term = x['raw_inputform']
        demand(term.startswith('(2*hbar*') and 'GammaCC' in term and
               'Coupling[lambdaT3' in term and 'Sqrt[3]' in term and
               'y1Conj' not in term and 'y2Conj' not in term,
               'Direct LLHH source changed')
        demand('NewScalar' not in term and term.count('Field[l, Fermion') == 2,
               'Not a direct LLHH operator')
        chi = 'HH' if 'Proj[-1]' in term else 'HdagHdag'
        demand(('Proj[1]' not in term) if chi == 'HH' else ('Proj[-1]' not in term),
               'Ambiguous chirality')
        kind = 'pole' if 'ϵ' in term else ('log' if 'Log[' in term else 'constant')
        demand(kind not in out.get(chi, {}), 'Duplicate Laurent structure')
        # Remove only the source's dimension-regulator or scalar mass log factor;
        # the remaining ordered GammaCC operator must be exactly identical.
        base = term.replace('ϵ*', '').replace('ε*', '')
        if kind == 'log':
            base, count = re.subn(
                r'\*Log\[μbar2/Coupling\[MF, \{Index\[[^,\]]+, NFlavor\]\}, 0\]\^2\]',
                '', base)
            demand(count == 1, 'Unrecognised native hard log')
        canonical.setdefault(chi, {})[kind] = base
        out.setdefault(chi, {})[kind] = x['index']
    demand(set(out) == {'HH', 'HdagHdag'}, 'Missing direct conjugate')
    demand(all(set(v) == {'pole', 'log', 'constant'} for v in out.values()),
           'Incomplete hard Laurent triplet')
    demand(all(len(set(v.values())) == 1 for v in canonical.values()),
           'Hard Laurent triplets contain different ordered operators')
    return out


def evaluate(sources_zip: Path, soft_zip: Path, mu_f: float) -> dict:
    with zipfile.ZipFile(sources_zip) as z:
        source = read_member(z, RAW)
        ledger_data = read_member(z, LEDGER)
        f_data = read_member(z, F_INPUT)
    with zipfile.ZipFile(soft_zip) as z:
        soft_data = read_member(z, SOFT_MEMBER)
    ledger = json.loads(ledger_data)
    f = json.loads(f_data)
    hard_triplets = audit_hard_direct(ledger, source)
    soft_symbolic = audit_soft_source(soft_data.decode('utf-8-sig'))
    demand(f['model'] == MODEL and f['scale_GeV'] == mu_f,
           'Benchmark must be A(-4) at the common fermion scale')
    demand(f['all_heavy_fermions_to_integrate_out_together'] is True and
           f['both_heavy_scalar_multiplets_remain_light'] is True,
           'Incorrect first-threshold field inventory')
    params = f['native_uv_rgbeta_couplings_at_muF']
    demand('y1Conj' not in params and 'y2Conj' not in params,
           'Independent conjugate Yukawa present')
    mf = matrix(params['MF'])
    y1 = matrix(params['y1'])
    y2 = matrix(params['y2'])
    demand(mf.shape == y1.shape == y2.shape == (3, 3),
           'Expected three heavy generations and three lepton flavours')
    demand(all(np.isfinite(t).all() for t in (mf, y1, y2)), 'Nonfinite UV matrices')
    relative_offdiag = float(np.linalg.norm(mf - np.diag(np.diag(mf))) /
                             np.linalg.norm(mf))
    demand(relative_offdiag < 1e-8, 'Non-diagonal mass basis not supported')
    masses = np.diag(mf)
    demand(np.max(abs(masses.imag)) < 1e-8 * max(1., np.max(abs(masses.real)))
           and np.min(masses.real) > 0, 'Heavy masses need a verified rotation')
    masses = masses.real
    lam = complex_scalar(params['lambdaT3'])
    a, b = float(params['mS1Sq']), float(params['mS2Sq'])
    demand(min(a, b) > 0, 'Tachyonic scalar input at fermion threshold')
    mu2 = mu_f**2
    bubble = scalar_bubble_finite(a, b, mu2)
    # Ordered native Wbar = sum_r conj(y1_pr)*conj(y2_qr)/M_r.
    # Assumes the verified near-diagonal heavy basis, not a general SVD mapping.
    wbar = (y1.conj() / masses[np.newaxis, :]) @ y2.conj().T
    normalization = 2 / (math.sqrt(3) * 16 * math.pi**2)
    # The explicit hbar-to-1/(16pi²) identification is a convention choice.
    hard_pole = normalization * lam * wbar
    soft_pole = -hard_pole
    soft_finite = soft_pole * bubble
    # Triplet log piece is source-verified, but only a direct-channel hard term.
    logs = 1 + np.log(mu2 / masses**2)
    hard_finite = normalization * lam * (
        (y1.conj() * logs[np.newaxis, :] / masses[np.newaxis, :]) @ y2.conj().T)
    # This sum is a diagnostic channel combination, NOT subtraction-certified.
    paired_finite = hard_finite + soft_finite
    cancellation = float(np.linalg.norm(hard_pole + soft_pole))
    demand(cancellation <= 1e-13 * max(1e-40, np.linalg.norm(hard_pole)),
           'Direct reduced pole pairing failed')
    return {
        'schema': SCHEMA,
        'model': MODEL,
        'mu_F_GeV': mu_f,
        'status': 'A4_DIRECT_EFT1_SCALAR_BUBBLE_CANDIDATE_NOT_RENORMALISED_MATCHING',
        'integrated_fields_at_muF': ['F'],
        'retained_scalar_multiplets': ['S1', 'S2'],
        'source_sha256': {
            'uv_raw_one_loop': sha256(source).hexdigest(),
            'uv_operator_ledger': sha256(ledger_data).hexdigest(),
            'fermion_benchmark': sha256(f_data).hexdigest(),
            'scalar_soft_native_expression': sha256(soft_data).hexdigest(),
        },
        'direct_hard_source_indices': hard_triplets,
        'soft_symbolic_source': soft_symbolic,
        'scalar_mass_squared_GeV2_at_muF': [a, b],
        'heavy_masses_GeV_at_muF': masses.tolist(),
        'heavy_offdiagonal_relative_norm': relative_offdiag,
        'scalar_bubble_finite_dimensionless': bubble,
        'wbar_ordered_GeV_inv': native_encode(wbar),
        'hard_direct_pole_GeV_inv': native_encode(hard_pole),
        'soft_candidate_pole_GeV_inv': native_encode(soft_pole),
        'hard_direct_finite_GeV_inv': native_encode(hard_finite),
        'soft_candidate_finite_GeV_inv': native_encode(soft_finite),
        'paired_direct_finite_diagnostic_GeV_inv': native_encode(paired_finite),
        'hard_direct_pole_frobenius_GeV_inv': float(np.linalg.norm(hard_pole)),
        'soft_candidate_finite_frobenius_GeV_inv': float(np.linalg.norm(soft_finite)),
        'paired_direct_finite_frobenius_GeV_inv': float(np.linalg.norm(paired_finite)),
        'direct_ordered_pole_sum_frobenius_GeV_inv': cancellation,
        'hbar_identified_as_one_over_16pi_squared': True,
        'counterterms_verified': False,
        'full_eft1_operator_basis_subtracted': False,
        'uv_hard_minus_eft1_soft_verified': False,
        'operator_GammaCC_to_Weyl_verified': False,
        'fermion_threshold_finite_wilson_array_available': False,
        'source_limitation': (
            'The only original scalar-soft formula supplied is a later joint-scalar '
            'threshold expression. It is used as a same-topology EFT1 bubble '
            'candidate at muF, not proof of the UV-minus-EFT1 matching prescription.'),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources-zip', type=Path, required=True)
    p.add_argument('--soft-zip', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mu-f', type=float, default=100000.0)
    args = p.parse_args()
    result = evaluate(args.sources_zip, args.soft_zip, args.mu_f)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print('A_-4 direct hard pole: {:.9e} GeV^-1'.format(result['hard_direct_pole_frobenius_GeV_inv']))
    print('A_-4 EFT1 scalar bubble: {:.12g}'.format(result['scalar_bubble_finite_dimensionless']))
    print('A_-4 candidate soft finite: {:.9e} GeV^-1'.format(result['soft_candidate_finite_frobenius_GeV_inv']))
    print('A_-4 paired direct diagnostic: {:.9e} GeV^-1'.format(result['paired_direct_finite_frobenius_GeV_inv']))
    print('Not a renormalised fermion-threshold Wilson boundary')
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
