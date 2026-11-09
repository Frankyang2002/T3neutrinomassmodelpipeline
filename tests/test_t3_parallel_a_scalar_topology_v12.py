"""v12 lambdaT3 source-CG identities, all16 actions and negative controls."""
from __future__ import annotations
from pathlib import Path
import tempfile,unittest
import numpy as np
from Numerical.parallel_t3.wilson.scalar_topology_action_v12 import (
 V12Blocked,audit,doublet_triplet_cg,triple_invariant,topology_tensor,
 topology_source_complex_action,topology_constructed_real_lam,relative,
)

class TestTopologyCG(unittest.TestCase):
 def test_symmetric_cg_isometry_and_triple_epsilon(self):
    cg=doublet_triplet_cg()
    self.assertEqual(cg.shape,(3,2,2))
    self.assertLess(relative(cg.reshape(3,4)@cg.reshape(3,4).conj().T,
                             np.eye(3)),1e-13)
    tri=triple_invariant(cg)
    self.assertAlmostEqual(tri[0,1,2].imag,-1/np.sqrt(2),places=12)
    self.assertAlmostEqual(tri[1,0,2].imag,1/np.sqrt(2),places=12)

 def test_source_topologies_have_correct_axes(self):
    for x in ((1,3),(3,1),(2,2),(3,3)):
        t=topology_tensor(*x)
        self.assertEqual(t.shape,(2,2,*x))
        self.assertGreater(np.linalg.norm(t),.1)

 def test_cg_wrong_scale_rejected(self):
    with self.assertRaisesRegex(V12Blocked,'three-triplet'):
        triple_invariant(doublet_triplet_cg()*2)

 def test_lambdaT3_16_actions_and_su2_covariance(self):
    r=audit(samples=2)
    self.assertEqual(r['lambdaT3_complex_quartic_slots_constructed'],16)
    self.assertLess(r['max_constructed_complex_action_residual'],1.e-11)
    self.assertLess(max(r['SU2_declared_frame_CG_checks'].values()),1.e-12)
    self.assertFalse(r['certification']['physical_Wilson_evolution_authorized'])

 def test_missing_source_plus_conjugate_or_six_rejected(self):
    t=topology_tensor(2,2)
    x=np.array([.5,.2,-.3,.7,.6,-.1,.9,.4,.2,-.2,.1,.5])
    c=.9+.4j
    expect=topology_source_complex_action(t,c,x)
    lam=topology_constructed_real_lam(t,c)
    got=np.einsum('ijkl,i,j,k,l->',lam,x,x,x,x)/24
    self.assertLess(abs(expect-got),1.e-12)
    self.assertGreater(abs(expect-got/6),1.e-7)

 def test_invalid_topology_dimension_rejected(self):
    with self.assertRaisesRegex(V12Blocked,'unsupported'):
        topology_tensor(1,1)

 def test_provenance_without_source_blocked(self):
    with tempfile.TemporaryDirectory() as td:
        with self.assertRaisesRegex(V12Blocked,'topology source missing'):
            audit(root=Path(td),samples=1)

if __name__=='__main__':unittest.main()
