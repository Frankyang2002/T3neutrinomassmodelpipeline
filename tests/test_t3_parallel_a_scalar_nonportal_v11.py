"""v11 basis-independent rank, source complex action and negative controls."""
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from Numerical.parallel_t3.wilson.scalar_nonportal_rank_v11 import (
    V11Blocked,audit,source_invariant_values,su2_generators,
    special_invariant_tensor,source_complex_action,constructed_complex_lam,
    triplet_self_values,rank_report,relative,
)
from Numerical.parallel_t3.wilson.native_quartic_projectors_v10 import SPECIAL

class TestSU2Nonportal(unittest.TestCase):
    def test_16_models_178_registration_14_10_4_12(self):
        r=audit(samples=40)
        self.assertEqual(r['source_registration_slots_total'],178)
        self.assertEqual(r['rank_sector_counts'],{
            'triplet_self_sectors':14,'mixed_adjoint_sectors':10,
            'mixed_crossed_triplet_sectors':4,'special_complex_slots':12})
        self.assertLess(r['max_invariant_SU2_covariance_error'],1e-10)
        self.assertLess(r['max_special_constructed_action_error'],1e-10)

    def test_fundamental_and_real_cartesian_triplet_indices(self):
        for d,index in ((2,.5),(3,2)):
            g=su2_generators(d)
            self.assertLess(relative(sum(x@x for x in g),{2:.75,3:2}[d]*np.eye(d)),1.e-12)
            self.assertAlmostEqual(np.trace(g[0]@g[0]).real,index)
            self.assertLess(np.linalg.norm(g[0]-g[0].conj().T),1.e-12)
        self.assertLess(np.linalg.norm(su2_generators(3).real),1.e-12)

    def test_triplet_cross_polynomial_real_positive(self):
        s=np.array([1.+1j,0.2-1j,2.+.5j]);t=np.array([.4-1j,1.+.2j,.6-1j])
        v=source_invariant_values(3,3,s,t)
        self.assertAlmostEqual(v['lambda12Cross'],abs(np.dot(s,t))**2+abs(np.vdot(s,t))**2)
        self.assertGreater(v['lambda12Cross'],0.)

    def test_triplet_self_representation_identity(self):
        rng=np.random.default_rng(90)
        for _ in range(30):
            s=rng.normal(size=3)+1j*rng.normal(size=3)
            values=triplet_self_values(s)
            self.assertAlmostEqual(2*values['self_adj'],
                np.vdot(s,s).real**2-abs(np.dot(s,s))**2,places=9)

    def test_cross_rank_negative_duplicate_invariant(self):
        x=np.linspace(.1,2,64)
        with self.assertRaisesRegex(V11Blocked,'rank deficient'):
            rank_report(np.array([x,x]).T.tolist(),['norm','duplicate'],'synthetic')

    def test_6_special_source_tensors_and_conjugate_actions(self):
        rng=np.random.default_rng(22)
        for name in SPECIAL:
            fields,t=special_invariant_tensor(name)
            self.assertEqual(t.shape,(2,2,2,2))
            n=4*len({f[4:-1] if f.startswith('Bar[') else f for f in fields})
            x=rng.normal(size=n)
            c=.37+.59j
            a=source_complex_action(t,fields,c,x)
            l=constructed_complex_lam(t,fields,c)
            self.assertLess(abs(a-np.einsum('ijkl,i,j,k,l->',l,x,x,x,x)/24),1.e-10)

    def test_missing_complex_conjugate_or_six_prefactor_rejected_by_action(self):
        fields,t=special_invariant_tensor('lambdaHHdagS1barS1bar')
        x=np.array([.8,.4,.2,-.4, 1.,-.3,-.7,.3])
        c=.9+.6j
        want=source_complex_action(t,fields,c,x)
        good=constructed_complex_lam(t,fields,c)
        got=np.einsum('ijkl,i,j,k,l->',good,x,x,x,x)/24
        self.assertLess(abs(want-got),1.e-12)
        self.assertGreater(abs(want-got/6),1.e-5)

    def test_unsupported_scalar_dimensions_rejected(self):
        with self.assertRaisesRegex(V11Blocked,'Unsupported generator dimension'):
            su2_generators(4)

    def test_source_hash_guard_rejects_missing_current_model(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(V11Blocked,'model source missing'):
                audit(samples=32,root=Path(td))

if __name__=='__main__':unittest.main()
