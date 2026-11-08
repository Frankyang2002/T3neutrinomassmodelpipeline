"""Read-only independent one-loop comparison, gaugeless/full complex tensors.

Run from repository root (no physical input needed):
    python scripts/CheckBatch34IndependentMaster.py \
      --output output/batch34_independent/comparison.json

Exit status 2 is INTENTIONAL: source-v1 and candidate seven-term beta
fail the independent one-loop reference, and neither may be certified.
Pass --diagnostic to return 0 for exploration/report generation only.
"""
from __future__ import annotations

import argparse
import json
from itertools import permutations
from pathlib import Path
import sys

import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.GeneralEFT1IndependentOneLoopReference import (
    equation22_gaugeless, sum_reference, CHANNELS)
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings, master_beta, total_beta)
from RGE.running.intermediate.GeneralEFT1CrossedWardCandidate import candidate_beta

FORMAT='t3-batch34-independent-one-loop-audit-v1'


def fixture(seed=311,nf=3,ns=4):
    """Ungauged synthetic Weyl/real-scalar theory: unrestricted complex phases."""
    rng=np.random.default_rng(seed)
    def rand(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    y=rand((nf,nf,ns));y=.2*(y+y.swapaxes(0,1))*.5
    q=rng.normal(size=(ns,)*4)
    lam=.05*sum(q.transpose(p) for p in permutations(range(4)))/24
    w=rand((nf,nf,ns,ns))
    c=.3*(w+w.swapaxes(0,1)+w.swapaxes(2,3)+w.swapaxes(0,1).swapaxes(2,3))/4
    return c,Couplings(y,lam,())


def diagnostics(seed=311) -> dict:
    c,inputs=fixture(seed)
    independent=equation22_gaugeless(c,inputs.yukawa,inputs.quartic)
    reference=sum_reference(independent)
    source4=master_beta(c,inputs,corrected_crossed=True)
    source2=master_beta(c,inputs,corrected_crossed=False)
    candidate=candidate_beta(c,inputs)
    without_crossed=total_beta(source4)-source4['crossed_yukawa']
    def residual(x):
        return float(np.max(np.abs(x-reference)))
    scale=float(np.max(np.abs(reference)))
    per_term={key:float(np.max(np.abs(independent[key]-source4[key])))
              for key in CHANNELS}
    result={
        'format':FORMAT,
        'status':'BLOCKED_CONFLICTING_ONE_LOOP_FORMULAE',
        'source_git_revision':'12128de602ac6bf109dbad50dbac8446ed629434',
        'independent_reference':'Fonseca et al JHEP07(2026)094 Eq.(22) arXiv:2512.15866v1',
        'compared_reference':'Aebischer et al arXiv:2512.16890v1 Eq.(4.85)',
        'reference_operator_normalization':'L5=(1/4)a_ijab psi_i psi_j phi_a phi_b+hc',
        'engine_operator_normalization':'L5=C_ijab psi_i psi_j phi_a phi_b+hc',
        'mapping':'a_ijab=4*C_ijab; linear/real-linear beta scales by the same factor',
        'scope':'g=0, other dimension-five coefficients=0; three Weyl, four real scalars; ALL complex indices',
        'synthetic':True,
        'random_seed':seed,
        'independent_total_norm':float(np.linalg.norm(reference)),
        'independent_total_max_absolute':scale,
        'per_non_crossed_term_absolute_discrepancy':per_term,
        'non_crossed_max_absolute_discrepancy':residual(without_crossed),
        'source_v1_four_exchange_max_absolute_discrepancy':residual(total_beta(source4)),
        'source_repository_two_exchange_max_absolute_discrepancy':residual(total_beta(source2)),
        'Cstar_candidate_max_absolute_discrepancy':residual(total_beta(candidate)),
        'source_four_exchange_crossed_norm':float(np.linalg.norm(source4['crossed_yukawa'])),
        'source_two_exchange_crossed_norm':float(np.linalg.norm(source2['crossed_yukawa'])),
        'Cstar_candidate_crossed_norm':float(np.linalg.norm(candidate['crossed_yukawa'])),
        'relative_without_crossed':residual(without_crossed)/scale,
        'relative_source_v1_four_exchange':residual(total_beta(source4))/scale,
        'relative_Cstar_candidate':residual(total_beta(candidate))/scale,
        'operator_basis_or_loop_formula_physically_certified':False,
        'note':('Eq.(22) reproduces six terms exactly in the gaugeless limit. '
                'No independent crossed Yukawa structure exists in that reference. '
                'Ward covariance is necessary but cannot settle a discrepant loop coefficient. '
                'Do not alter source-v1 or promote the candidate until this conflict is resolved.'),
    }
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--diagnostic',action='store_true',help='Return zero despite conflicting formulae')
    parser.add_argument('--seed',type=int,default=311)
    args=parser.parse_args(argv)
    report=diagnostics(args.seed)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x',encoding='utf8') as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write('\n')
    print(report['status'])
    for key in ('non_crossed_max_absolute_discrepancy',
                'source_v1_four_exchange_max_absolute_discrepancy',
                'source_repository_two_exchange_max_absolute_discrepancy',
                'Cstar_candidate_max_absolute_discrepancy'):
        print(f'{key}: {report[key]:.12g}')
    # Failure by design unless --diagnostic; this is a block on physical promotion.
    return 0 if args.diagnostic else 2

if __name__=='__main__':
    raise SystemExit(main())
