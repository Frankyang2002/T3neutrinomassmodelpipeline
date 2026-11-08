"""Fast algebra/representation contracts for reduced-potential all-16 trial."""
from __future__ import annotations
import unittest
import numpy as np
from Numerical.running.All16ConditionalFullWilsonTrial import (
    SCALAR_CONVENTION_FACTOR,LEPTON_POSITIONS,original_scalar_generators,
    rotated_scalar_generators,higgs_quartic,hhss_full_quartic,
    embedded_boundary,generic_leptons,midpoint_transport,
    native_triplet_su2_generators, source_ward_residual,
    quartic_ward_residual, DEFAULT_RELEASE, DEFAULT_RUNNING, ROOT,
)
from Numerical.joint_scalar_matching.adapter import ALL_MODELS,CLASS_DIMS

class All16ReducedTrialTests(unittest.TestCase):
    def test_all16_inventory(self):
        self.assertEqual(len(ALL_MODELS),16)
        self.assertEqual(len(set(ALL_MODELS)),16)
        self.assertEqual(set(CLASS_DIMS),set('ABCDE'))

    def test_45_weyl_lepton_positions(self):
        self.assertEqual(LEPTON_POSITIONS.shape,(3,2))
        self.assertEqual(len(set(LEPTON_POSITIONS.flat)),6)
        self.assertTrue(np.all(LEPTON_POSITIONS>=0))
        self.assertTrue(np.all(LEPTON_POSITIONS<45))

    def test_half_dictionary_is_explicit(self):
        self.assertEqual(SCALAR_CONVENTION_FACTOR,0.5)

    def test_generic_SU2_and_U1_generator_shapes(self):
        for model in ALL_MODELS:
            n=4+2*sum(CLASS_DIMS[model[0]][:2])
            g=original_scalar_generators(model)
            self.assertEqual(len(g['SU2']),3)
            self.assertEqual(len(g['U1']),1)
            self.assertEqual(len(g['SU3']),8)
            for grp in g.values():
                for x in grp:
                    self.assertEqual(x.shape,(n,n))
                    np.testing.assert_allclose(x,x.conj().T,atol=1e-13,rtol=0)

    def test_real_orthogonal_rotation_preserves_generator_hermiticity(self):
        rotation=np.eye(8)
        rotation[0,0]=0;rotation[0,1]=1
        rotation[1,0]=-1;rotation[1,1]=0
        g=rotated_scalar_generators('B_-1',rotation)
        for group in g.values():
            for t in group:
                np.testing.assert_allclose(t,t.conj().T,atol=1e-13,rtol=0)

    def test_higgs_quartic_exact_columns(self):
        for n in (12,16):
            q=higgs_quartic(n)
            self.assertEqual(q.shape,(n,)*4)
            self.assertEqual(q[0,0,0,0],3)
            self.assertEqual(q[0,0,1,1],1)
            self.assertEqual(q[0,1,0,1],1)
            self.assertEqual(q[4,4,4,4],0)

    def test_HHSS_tensor_has_all_permutations_and_no_other_channels(self):
        q=np.zeros((4,4,8,8));q[0,0,1,2]=q[0,0,2,1]=0.5
        full=hhss_full_quartic(q)
        for inds in [(0,0,5,6),(0,5,0,6),(5,0,6,0),(6,5,0,0)]:
            self.assertAlmostEqual(full[inds],0.5)
        self.assertEqual(full[0,0,0,0],0)
        self.assertEqual(full[5,5,5,5],0)

    def test_full_master_embedding_preserves_both_symmetries(self):
        source=np.zeros((3,3,2,2,8,8),complex)
        source[0,1,0,1,0,3]=source[1,0,1,0,3,0]=2+1j
        source[1,0,1,0,0,3]=source[0,1,0,1,3,0]=2+1j
        full=embedded_boundary(source,np.eye(8))
        self.assertEqual(full.shape,(45,45,12,12))
        np.testing.assert_allclose(full,full.swapaxes(0,1),rtol=0,atol=0)
        np.testing.assert_allclose(full,full.swapaxes(2,3),rtol=0,atol=0)
        heavy=generic_leptons(full,True)
        self.assertEqual(heavy[0,1,0,1,0,3],1+0.5j)

    def test_matchete_cartesian_triplet_hermitian_su2_algebra(self):
        gens=native_triplet_su2_generators()
        self.assertEqual(len(gens),3)
        for t in gens:
            np.testing.assert_allclose(t,t.conj().T,atol=1e-14,rtol=0)
            np.testing.assert_allclose(t,-t.T,atol=1e-14,rtol=0)
        np.testing.assert_allclose(gens[0]@gens[1]-gens[1]@gens[0],
                                   1j*gens[2],atol=1e-13,rtol=0)

    def test_actual_A_D_E_source_component_ward_witness(self):
        # Actual A(-4), D(-2), E(-4) p=q=1 source-component witnesses
        # read from the model source release, scaled by the common amplitude.
        # First row: i,j,heavy-a,heavy-b,real,imag, all 1-based local indices.
        witnesses={
            'A_-4': ((1,3), [
                (1,1,1,3,1,0),(1,1,1,4,0,1),(1,1,1,7,0,1),
                (1,1,1,8,-1,0),(1,1,2,3,0,-1),(1,1,2,4,1,0),
                (1,1,2,7,1,0),(1,1,2,8,0,1),
                (1,2,1,5,-1,0),(1,2,1,6,0,-1),
                (1,2,2,5,0,1),(1,2,2,6,-1,0),
                (2,2,1,3,-1,0),(2,2,1,4,0,-1),
                (2,2,1,7,0,1),(2,2,1,8,-1,0),
                (2,2,2,3,0,1),(2,2,2,4,-1,0),
                (2,2,2,7,1,0),(2,2,2,8,0,1),
            ]),
            'D_-2': ((3,1), [
                (1,1,1,7,1,0),(1,1,2,7,0,-1),
                (1,1,5,7,0,1),(1,1,6,7,1,0),
                (1,2,3,7,-1,0),(1,2,4,7,0,1),
                (2,2,1,7,-1,0),(2,2,2,7,0,1),
                (2,2,5,7,0,1),(2,2,6,7,1,0),
            ]),
            'E_-4': ((3,3), [
                (1,1,1,9,1,0),(1,1,1,10,0,1),
                (1,1,2,9,0,-1),(1,1,2,10,1,0),
                (1,1,3,7,-1,0),(1,1,3,8,0,-1),
                (1,1,3,11,0,-1),(1,1,3,12,1,0),
                (1,1,4,7,0,1),(1,1,4,8,-1,0),
                (1,1,4,11,-1,0),(1,1,4,12,0,-1),
                (1,1,5,9,0,1),(1,1,5,10,-1,0),
                (1,1,6,9,1,0),(1,1,6,10,0,1),
                (1,2,1,11,0,-1),(1,2,1,12,1,0),
                (1,2,2,11,-1,0),(1,2,2,12,0,-1),
                (1,2,5,7,0,1),(1,2,5,8,-1,0),
                (1,2,6,7,1,0),(1,2,6,8,0,1),
                (2,2,1,9,1,0),(2,2,1,10,0,1),
                (2,2,2,9,0,-1),(2,2,2,10,1,0),
                (2,2,3,7,-1,0),(2,2,3,8,0,-1),
                (2,2,3,11,0,1),(2,2,3,12,-1,0),
                (2,2,4,7,0,1),(2,2,4,8,-1,0),
                (2,2,4,11,1,0),(2,2,4,12,0,1),
                (2,2,5,9,0,-1),(2,2,5,10,1,0),
                (2,2,6,9,-1,0),(2,2,6,10,0,-1),
            ]),
        }
        from Numerical.running import Aminus4EFT1Couplings as eft
        for model,(dims,rows) in witnesses.items():
            nh=2*sum(dims)
            src=np.zeros((3,3,2,2,nh,nh),complex)
            for i,j,a,b,re,im in rows:
                for w,x in {(i-1,j-1),(j-1,i-1)}:
                    for c,d in {(a-1,b-1),(b-1,a-1)}:
                        src[0,0,w,x,c,d]=re+1j*im
            native=original_scalar_generators(model)
            audit=source_ward_residual(src,native)
            self.assertLess(max(audit.values()),5e-13,(model,audit))
            # The previous incorrect spherical generators should fail badly.
            wrong={k:list(v) for k,v in native.items()}
            off=0
            for dim in dims:
                if dim==3:
                    for index,t in enumerate(eft.su2_generators(3)):
                        wrong['SU2'][index]=wrong['SU2'][index].copy()
                        wrong['SU2'][index][4+off:4+off+6,4+off:4+off+6]=eft.real_scalar_generator(t)
                off+=2*dim
            wrong_audit=source_ward_residual(src,wrong)
            self.assertGreater(wrong_audit['SU2'],1.,model)

    def test_native_triplet_source_ward_covariant_after_orthogonal_rotation(self):
        from Numerical.running import Aminus4EFT1Couplings as eft
        # Construct one nonzero component in a known SU(2)-invariant tensor:
        # C(LL,S1,S2) with the actual A(-4) coefficients is tested above.
        # Here the generators must transform by R^T T R, not R T R^T.
        rotation=np.eye(8)
        rotation[2,2]=0;rotation[2,3]=1
        rotation[3,3]=0;rotation[3,2]=-1
        g=rotated_scalar_generators('A_-4',rotation)
        old=original_scalar_generators('A_-4')
        full=np.eye(12)
        full[4:,4:]=rotation
        for name in ('SU2','U1'):
            for a,b in zip(g[name],old[name]):
                np.testing.assert_allclose(a,full.T@b@full,atol=1e-13,rtol=0)

    def test_zero_beta_midpoint_keeps_full_C(self):
        import unittest.mock as mock
        c=np.zeros((45,45,12,12),complex)
        c[5,5,4,4]=1j
        y=np.zeros((45,45,12),complex)
        f={k:.1 for k in ('gY','g2','g3')}
        sf={'SU3':(np.zeros((45,45),complex),)*8,
            'SU2':(np.zeros((45,45),complex),)*3,
            'U1':(np.zeros((45,45),complex),)}
        ss={'SU3':(np.zeros((12,12),complex),)*8,
            'SU2':(np.zeros((12,12),complex),)*3,
            'U1':(np.zeros((12,12),complex),)}
        with mock.patch('Numerical.running.All16ConditionalFullWilsonTrial.six_term_at',
                        side_effect=lambda x,*args:np.zeros_like(x)):
            out,info=midpoint_transport(c,y,y,f,f,np.zeros((4,4,8,8)),
                                        0.,0.,sf,ss,1e5,1e3,3)
        np.testing.assert_array_equal(out,c)
        self.assertEqual(info['beta_rhs_evaluations'],6)

if __name__=='__main__':unittest.main()
