"""Read-only A(-4) mixed S1-S2 Clebsch-Gordan structure audit."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical


def ward(c, f, t):
    return (np.einsum('ki,kjab->ijab', f, c, optimize=True)
            + np.einsum('kj,ikab->ijab', f, c, optimize=True)
            + np.einsum('ca,ijcb->ijab', t, c, optimize=True)
            + np.einsum('cb,ijac->ijab', t, c, optimize=True))


def norm(x):
    return float(np.linalg.norm(x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--capture', type=Path,
                    default=Path('output/parallel_t3/all16_ward_capture01'))
    args = ap.parse_args()
    root = Path.cwd()
    for d in sorted(args.capture.glob('fit16_*_ward_diagnostic.json')):
        token = d.name.removesuffix('_ward_diagnostic.json')
        case = load_canonical(root,
            root / 'output/parallel_t3/rough_inputs' / token / 'a/A_m4',
            'A_-4', c5_boundary='require', assume_formal_identity=True)
        b = case.basis
        if b.dims != (1, 3, 2):
            raise ValueError('Only A model (S1 singlet, S2 triplet) supported')
        lidx = [i for pair in b.leptons for i in pair]
        hidx = list(range(4, b.ns))
        c = case.initial[np.ix_(lidx, lidx, hidx, hidx)]
        generators = [(f[np.ix_(lidx,lidx)], t[np.ix_(hidx,hidx)])
                      for f,t in b.generators()['SU2']]
        def casimir(x):
            return sum((ward(ward(x,f,t),f,t) for f,t in generators),
                       np.zeros_like(x))
        q = casimir(c)
        invariant = c - (2/3)*q + casimir(q)/12
        # Fixed normalized weak/scalar pattern for each lepton-flavour pair.
        # The physical mixed operator has one SU2 singlet and one triplet.
        # Thus its invariant SU2 tensor space is one-dimensional (up to
        # scalar R/I complex structures and the flavour coefficient).
        details = []
        for p in range(3):
            for r in range(p,3):
                ix = np.ix_(range(2*p,2*p+2),range(2*r,2*r+2),
                            range(0,2),range(2,8))
                x = c[ix]
                sing = invariant[ix]
                bad = x-sing
                details.append({
                    'flavours':[p+1,r+1],
                    'mixed_norm':norm(x),
                    'singlet_norm':norm(sing),
                    'non_singlet_norm':norm(bad),
                    'singlet_squared_fraction':norm(sing)**2/max(norm(x)**2,1e-200),
                    'singlet_overlap_complex':[
                        float(np.vdot(sing,x).real/max(norm(sing)*norm(x),1e-200)),
                        float(np.vdot(sing,x).imag/max(norm(sing)*norm(x),1e-200))],
                    'raw_weak_block_norms':[[norm(x[i,j]) for j in range(2)] for i in range(2)],
                    'singlet_weak_block_norms':[[norm(sing[i,j]) for j in range(2)] for i in range(2)],
                    'non_singlet_weak_block_norms':[[norm(bad[i,j]) for j in range(2)] for i in range(2)],
                })
        # Characterise how much contamination resides in each Casimir spin.
        # P0=(Q-2)(Q-6)/12; P1=-Q(Q-6)/8; P2=Q(Q-2)/24.
        spin1 = (6*q-casimir(q))/8
        spin2 = (casimir(q)-2*q)/24
        fractions = {'J0':norm(invariant)**2/max(norm(c)**2,1e-200),
                     'J1':norm(spin1)**2/max(norm(c)**2,1e-200),
                     'J2':norm(spin2)**2/max(norm(c)**2,1e-200)}
        report = {
            'trial':token, 'status':'READ_ONLY_NOT_A_REPAIR',
            'SU2_total_spin_squared_fractions':fractions,
            'spin_sector_reconstruction_error':norm(invariant+spin1+spin2-c)/max(norm(c),1e-100),
            'singlet_ward_relative':[norm(ward(invariant,f,t))/max(norm(invariant),1e-100)
                                     for f,t in generators],
            'flavour_pair_mixed_CG_comparison':details,
            'interpretation':'Measures the actual mixed LL S1 S2 tensor against the unique SU2 singlet subspace. Does not determine the correct UV spinor normalization, complex-field conjugations, or original source contraction.'}
        out = args.capture / f'{token}_mixed_CG_comparison.json'
        out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('\n',token,'spin fractions:',fractions,
              'reconstruction:',report['spin_sector_reconstruction_error'])
        for item in details:
            print('flavours',item['flavours'], 'norm',item['mixed_norm'],
                  'singlet fraction',item['singlet_squared_fraction'])
        print('report:',out)

if __name__ == '__main__':
    main()
