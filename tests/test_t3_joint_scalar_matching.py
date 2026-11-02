"""No Wolfram dependency. Synthetic algebra tests; NOT UV/SMEFT amplitude certification."""
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from Numerical.joint_scalar_matching.adapter import (
    ALL_MODELS, CONVENTION, MIXED, SELF_SPLIT, load_input, save_input,
    normalize_model, scalar_groups,
)
from Numerical.joint_scalar_matching.fixtures import make_fixture
from Numerical.joint_scalar_matching.core import (
    b0bar, LOOP, EigenBasis, match_c5, project_weinberg,
    scalar_bubble_matching, weinberg_tensor,
)


class JointScalarTests(unittest.TestCase):
    def test_all_16_fixture_models_and_complex_flavour(self):
        self.assertEqual(len(ALL_MODELS),16)
        self.assertEqual(len(SELF_SPLIT),6)
        self.assertEqual(len(MIXED),2)
        self.assertEqual(len(set(ALL_MODELS)-SELF_SPLIT-MIXED),8)
        for model in ALL_MODELS:
            with self.subTest(model=model):
                inp=make_fixture(model)
                r=match_c5(inp.c_llss,inp.lambda_hhss,inp.mass_squared,
                           inp.mu_s_gev,c5_eft1=inp.c5_eft1,
                           reference_scale=inp.reference_scale_gev,
                           group_names=scalar_groups(model))
                self.assertLess(r["projection_residual"],1e-11)
                self.assertEqual(r["c5_total"].shape,(3,3))
                self.assertTrue(np.allclose(r["c5_total"],r["c5_total"].T,atol=1e-25))
                self.assertGreater(abs(r["c5_total"].imag).max(),1e-16)
                self.assertEqual(len(r["pair_data"]),
                                 inp.mass_squared.shape[0]*(inp.mass_squared.shape[0]+1)//2)
                self.assertTrue(np.allclose(r["c5_scalar_bubble"],
                                sum(r["pair_c5"].values()),rtol=1e-10,atol=1e-25))
                self.assertTrue(np.allclose(r["c5_total"],
                                r["c5_eft1"]+r["c5_scalar_bubble"],rtol=1e-11,atol=1e-28))
                self.assertTrue(np.allclose(r["c5_scalar_bubble"],
                                r["c5_reference_finite"]+r["c5_explicit_mu_log"],
                                rtol=1e-11,atol=1e-28))

    def test_loop_degeneracies_massless_limit_not_required(self):
        m2=1000.0**2
        self.assertAlmostEqual(b0bar(m2,m2,1000.),0.,places=13)
        self.assertAlmostEqual(b0bar(m2,m2,500.),math.log(0.25),places=13)
        for diff in (1e-2,1e-5,1e-8,1e-12,1e-14):
            x=m2*(1+diff); y=m2*(1-diff)
            val=b0bar(x,y,1000.)
            self.assertTrue(math.isfinite(val))
            self.assertAlmostEqual(val, diff*diff/6 + diff**4/20 + diff**6/42,delta=4e-14)
        a,b=1e6,1.6e6
        self.assertAlmostEqual(b0bar(a,b,1300.)-b0bar(a,b,700.),
                               2*math.log(1300/700),places=12)

    def test_vanishing_couplings_and_symmetric_cancellation(self):
        p=make_fixture('A_-4')
        zero=np.zeros_like(p.c_llss)
        a=match_c5(zero,p.lambda_hhss,p.mass_squared,p.mu_s_gev)
        self.assertTrue(np.array_equal(a["c5_total"],np.zeros((3,3))))
        b=match_c5(p.c_llss,np.zeros_like(p.lambda_hhss),p.mass_squared,p.mu_s_gev)
        self.assertTrue(np.array_equal(b["c5_total"],np.zeros((3,3))))
        # Equal-mass diagonal scalar pairs at mu_S have B0bar=0, even with C lambda !=0.
        n=p.mass_squared.shape[0]
        m2=np.eye(n)*p.mu_s_gev**2
        c=match_c5(p.c_llss,p.lambda_hhss,m2,p.mu_s_gev)
        self.assertLess(np.linalg.norm(c["c5_total"]),1e-30)
        # Opposing independently nonzero eigenstate channels cancel.
        C=np.zeros((3,3,2,2,2,2),complex)
        q=np.zeros((4,4,2,2))
        C[0,0,0,0,0,0]=C[0,0,0,0,1,1]=2e-10
        q[2,2,0,0]=1.0;q[2,2,1,1]=-1.0
        r=scalar_bubble_matching(C,q,np.eye(2)*(1250**2),1000.)
        self.assertLess(np.linalg.norm(r["amplitude"]),1e-28)
        self.assertGreater(np.linalg.norm(r["pair_amplitudes"]["e00_e00"]),1e-16)

    def test_mixed_basis_rotation_equivalence(self):
        for model in ('B_-1','C_-1'):
            inp=make_fixture(model)
            r=match_c5(inp.c_llss,inp.lambda_hhss,inp.mass_squared,inp.mu_s_gev)
            # Arbitrary real orthogonal change of original scalar coordinates.
            n=inp.mass_squared.shape[0]
            from numpy.linalg import qr
            O=qr(np.random.default_rng(47).normal(size=(n,n)))[0]
            newC=np.einsum('fgijab,ka,lb->fgijkl',inp.c_llss,O,O,optimize=True)
            newQ=np.einsum('ABab,ka,lb->ABkl',inp.lambda_hhss,O,O,optimize=True)
            newM=O@inp.mass_squared@O.T
            s=match_c5(newC,newQ,newM,inp.mu_s_gev)
            self.assertTrue(np.allclose(r["c5_scalar_bubble"],s["c5_scalar_bubble"],
                                        rtol=1e-10,atol=1e-24))

    def test_self_bilinear_split_real_scalar_contribution(self):
        # Nondegenerate real components require each separate eigenmass in B0bar.
        C=np.zeros((3,3,2,2,2,2),complex)
        Q=np.zeros((4,4,2,2))
        C[0,0,0,0,0,0]=1e-10*(1+.3j)
        C[0,0,0,0,1,1]=2e-10*(1+.3j)
        Q[2,2,0,0]=.2; Q[2,2,1,1]=.1
        masses2=np.diag([800.**2,1300.**2])
        r=scalar_bubble_matching(C,Q,masses2,1000.)
        expectation=-(.2*1e-10*b0bar(masses2[0,0],masses2[0,0],1000.)
                      +.1*2e-10*b0bar(masses2[1,1],masses2[1,1],1000.))*(1+.3j)/(2*LOOP)
        self.assertAlmostEqual(abs(r["amplitude"][0,0,0,0,2,2]-expectation),0.,delta=1e-26)
        # Partial, non-gauge-invariant toy amplitude must not be mislabelled C5.
        with self.assertRaises(ValueError):
            project_weinberg(r["amplitude"],strict=True)

    def test_rge_log_identity_matches_heavy_beta_not_double_counted(self):
        inp=make_fixture('E_+2')
        mu=inp.mu_s_gev
        h=1e-4
        plus=scalar_bubble_matching(inp.c_llss,inp.lambda_hhss,inp.mass_squared,mu*math.exp(h))
        minus=scalar_bubble_matching(inp.c_llss,inp.lambda_hhss,inp.mass_squared,mu*math.exp(-h))
        numerical=(plus["amplitude"]-minus["amplitude"])/(2*h)
        base=scalar_bubble_matching(inp.c_llss,inp.lambda_hhss,inp.mass_squared,mu)
        # Compare exact weighted tensor after rotation; no EFT1 running is fed back.
        U=base["mass_rotation_U"]
        from Numerical.joint_scalar_matching.core import rotate_to_eigenbasis
        CM,QM=rotate_to_eigenbasis(inp.c_llss,inp.lambda_hhss,U)
        analytic=-np.einsum('fgijab,ABab->fgijAB',CM,QM,optimize=True)/LOOP
        self.assertTrue(np.allclose(numerical,analytic,rtol=1e-9,atol=1e-26))

    def test_separately_labelled_quartic_channels_and_all_pairs(self):
        from scripts.RunT3JointScalarMatching import compute
        p=make_fixture("B_-1")
        result=compute(p)
        self.assertEqual(len(result["independent_quartic_channel_contributions"]),6)
        self.assertTrue(result["quartic_channels_explicitly_supplied"])
        def unpack(item):
            return np.asarray(item["real"])+1j*np.asarray(item["imag"])
        pieces=[unpack(row["C5_gev_inverse"])
                for row in result["independent_quartic_channel_contributions"]]
        total=unpack(result["contributions_gev_inverse"]["scalar_bubble_all_eigenstate_pairs"])
        self.assertTrue(np.allclose(np.sum(pieces,axis=0),total,rtol=1e-11,atol=1e-26))
        self.assertFalse(result["physical_full_one_loop_certified"])
        self.assertEqual(result["matching"],"both scalars together at mu_S")

    def test_strict_adapter_and_roundtrip(self):
        with tempfile.TemporaryDirectory() as td:
            path=save_input(Path(td)/"model.json",make_fixture('A_+0'))
            inp=load_input(path)
            self.assertEqual(inp.model,'A_+0')
            self.assertEqual(inp.scalar_mass_origin,'synthetic_algebra_only')
            self.assertEqual(inp.c_llss.shape,(3,3,2,2,8,8))
            self.assertEqual(inp.quartic_channels.shape,(6,4,4,8,8))
            self.assertTrue(np.allclose(np.sum(inp.quartic_channels,axis=0),inp.lambda_hhss))
            obj=__import__('json').loads(path.read_text(encoding='utf-8'))
            obj['tensor_convention']='UNKNOWN'
            path.write_text(__import__('json').dumps(obj),encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'tensor convention'):
                load_input(path)
        self.assertEqual(normalize_model('T3_dS1_2_dS2_2_dF_3_alpha_m1'),'C_-1')
        self.assertEqual(normalize_model('A_p0'),'A_+0')
        self.assertEqual(CONVENTION,'T3_REAL_LLSS_HHSS_WEINBERG_v1')


if __name__=='__main__':
    unittest.main()
