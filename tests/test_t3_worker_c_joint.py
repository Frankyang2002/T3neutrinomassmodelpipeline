"""Run from repository root: python -m unittest tests.test_t3_worker_c_joint -v"""
import unittest
import numpy as np
from Numerical.joint_scalar_matching.worker_c_joint import match_joint_scalars, CONVENTION


def data():
    n = 8  # A class: S1 singlet + S2 triplet, both complex
    meta = dict(mu_s_gev=1000.0, tensor_convention=CONVENTION,
                renormalization_scheme="MSbar", provenance={"source":"synthetic regression"})
    wilsons = {**meta, "C5_EFT1": np.diag([1e-12, 2e-12, 3e-12]).astype(complex),
               "C_LLSS": np.zeros((3,3,2,2,n,n), complex)}
    scalars = {**meta, "mass_squared_gev2": np.diag(np.linspace(900**2, 1500**2,n)),
               "lambda_HHSS": np.zeros((4,4,n,n)),
               "scalar_groups": ["S1"]*2 + ["S2"]*6,
               "scalar_mass_origin": "unbroken_electroweak"}
    return wilsons, scalars


class TestWorkerCJoint(unittest.TestCase):
    def test_inherited_once(self):
        w,s = data()
        result = match_joint_scalars(w,s,1000.)
        np.testing.assert_allclose(result["c5_total"],w["C5_EFT1"],rtol=0,atol=0)
        np.testing.assert_array_equal(result["c5_finite_scalar"],np.zeros((3,3)))
        self.assertFalse(result["full_one_loop_scalar_matching_completed"])
        self.assertEqual(len(result["contribution_breakdown"]["eigenstate_pairs"]),36)

    def test_channel_decomposition(self):
        w,s=data()
        s["lambda_HHSS_channels"]=np.zeros((2,4,4,8,8))
        s["quartic_channel_names"]=["portal_a","portal_b"]
        result=match_joint_scalars(w,s,1000.)
        self.assertEqual(set(result["contribution_breakdown"]["quartic_channels"]),{"portal_a","portal_b"})

    def test_wrong_scale_rejected(self):
        w,s=data()
        with self.assertRaisesRegex(ValueError,"evaluated at"):
            match_joint_scalars(w,s,999.)

    def test_missing_scalar_sector_rejected(self):
        w,s=data()
        s["scalar_groups"]=["S1"]*8
        with self.assertRaisesRegex(ValueError,"both S1 and S2"):
            match_joint_scalars(w,s,1000.)

    def test_inconsistent_quartic_channels_rejected(self):
        w,s=data()
        s["lambda_HHSS_channels"]=np.ones((1,4,4,8,8))
        s["quartic_channel_names"]=["wrong"]
        with self.assertRaisesRegex(ValueError,"channel sum"):
            match_joint_scalars(w,s,1000.)


if __name__ == "__main__":
    unittest.main()
