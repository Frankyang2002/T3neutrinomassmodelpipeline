import unittest
import numpy as np
from types import SimpleNamespace

from scripts.RunT3WorkerBBatchReconcile import (
    component_diagnostic, species_blocks, direction_geometry, hybrid_tensor,
    source_direction_columns, reconcile,
)


def h(a,b,i,j,value=1.):
    x=np.zeros((4,4,8,8))
    for A,B in {(a,b),(b,a)}:
        for I,J in {(i,j),(j,i)}:x[A,B,I,J]=value
    return x


class BatchWorkerBTests(unittest.TestCase):
    def test_source_hhss_overrides_only_six_orbits(self):
        from Numerical.running.T3CoupledQuarticContext import _put_hhss
        x=np.zeros((12,)*4); x[5,5,5,5]=17.;x[0,0,0,0]=11.
        x=_put_hhss(x,h(0,0,0,0,-1))
        final,record=hybrid_tensor(x,h(1,2,1,2,3))
        self.assertEqual(final[5,5,5,5],17.)
        self.assertEqual(final[0,0,0,0],11.)
        self.assertEqual(final[1,2,5,6],3.)
        self.assertTrue(record['non_hhss_unchanged_exactly'])
        self.assertEqual(record['full_bose_symmetry_relative'],0.)

    def test_subspace_detects_pure_sign(self):
        a=h(0,0,0,0);b=h(1,1,3,3)
        d=direction_geometry([('x','real'),('y','real')],np.stack((-a,b)),np.stack((a,b)))
        self.assertEqual(d['source_column_rank'],2)
        self.assertLess(d['candidate_direction_subspace_relative_residual'],1e-12)
        self.assertAlmostEqual(d['best_linear_direction_map_source_to_candidate'][0][0],-1.)

    def test_subspace_detects_missing_tensor_shape(self):
        a=h(0,0,0,0);b=h(1,1,3,3)
        d=direction_geometry([('x','real')],np.stack((b,)),np.stack((a,)))
        self.assertAlmostEqual(d['candidate_direction_subspace_relative_residual'],1.)

    def test_cross_direction_basis_mix(self):
        a=h(0,0,0,0);b=h(1,1,3,3)
        d=direction_geometry([('x','real'),('y','real')],np.stack((a+b,a-b)),np.stack((a,b)))
        self.assertLess(d['candidate_direction_subspace_relative_residual'],1e-12)
        self.assertAlmostEqual(d['best_linear_direction_map_source_to_candidate'][0][0],1.)
        self.assertAlmostEqual(d['best_linear_direction_map_source_to_candidate'][1][0],1.)

    def test_complex_source_axes(self):
        state={'lambdaT3':{'real':1.,'imag':2.}}
        def f(x):
            z=complex(x['lambdaT3']['real'],x['lambdaT3']['imag'])
            return z.real*h(0,0,0,0)+z.imag*h(1,1,3,3)
        base,cols=source_direction_columns(state,f,[('lambdaT3','real'),('lambdaT3','imag')])
        self.assertTrue(np.allclose(base,f(state)))
        self.assertTrue(np.allclose(cols[0],h(0,0,0,0)))
        self.assertTrue(np.allclose(cols[1],h(1,1,3,3)))

    def test_component_and_species(self):
        c=h(0,0,0,0,4.);s=h(0,0,0,0,1.)
        record=component_diagnostic(c,s,model='A_-4',top=2)
        self.assertAlmostEqual(record['maximum_absolute_difference'],3.)
        blocks=species_blocks(c,s,model='A_-4')
        self.assertAlmostEqual(blocks['S1S1']['relative_difference'],3.)
        self.assertEqual(blocks['S2S2']['source_norm'],0.)

    def test_d2_species_six_two(self):
        c=h(0,0,5,5,2.)+h(0,0,6,6,3.)
        s=h(0,0,5,5,1.)+h(0,0,6,6,1.)
        blocks=species_blocks(c,s,model='D_-2')
        self.assertAlmostEqual(blocks['S1S1']['best_scale_candidate_over_source'],2.)
        self.assertAlmostEqual(blocks['S2S2']['best_scale_candidate_over_source'],3.)
        diag=component_diagnostic(c,s,model='D_-2',top=1)
        self.assertEqual(diag['largest_components'][0]['heavy_species_pair'],'S2S2')

    def test_reject_non_symmetric_hybrid(self):
        from Numerical.running.T3CoupledQuarticContext import _put_hhss
        q=np.zeros((12,)*4); q[0,1,2,3]=1.
        with self.assertRaisesRegex(ValueError,'symmetry'):
            hybrid_tensor(q,h(0,0,0,0))

    def test_end_to_end_synthetic_source_hybrid(self):
        from Numerical.running.T3CoupledQuarticContext import _put_hhss
        a=h(0,0,0,0);b=h(1,1,3,3)
        def source(st):return float(st['lambdaA'])*a+float(st['lambdaB'])*b
        state={'lambdaA':0.2,'lambdaB':0.3}
        full=_put_hhss(np.zeros((12,)*4),-source(state))
        full[5,5,5,5]=3.
        cols=[]
        for direction in [a,b]:
            cols.append(_put_hhss(np.zeros((12,)*4),-direction))
        ctx=SimpleNamespace(model='A_-4',directions=(('lambdaA','real'),('lambdaB','real')),
                            columns=np.stack(cols),full_reference=lambda _:full)
        arrays,diag=reconcile('A_-4',ctx,state,source)
        self.assertTrue(np.allclose(arrays['hybrid_full_tensor'][:4,:4,4:,4:],source(state)))
        self.assertEqual(arrays['hybrid_full_tensor'][5,5,5,5],3.)
        self.assertLess(diag['source_affine_interpolation_relative'],1e-12)
        self.assertEqual(diag['hybrid_verification']['full_bose_symmetry_relative'],0.)


if __name__=='__main__':unittest.main()
