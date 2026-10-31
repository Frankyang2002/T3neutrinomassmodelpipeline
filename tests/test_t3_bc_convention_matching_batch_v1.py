import unittest
import numpy as np
from Numerical.experimental import T3BCConventionMatchingBatch as m

class ConventionTest(unittest.TestCase):
    def test_spectrum(self):
        masses=m.spectrum(np.diag([1e-14,2e-14,3e-14]),246.)
        self.assertEqual(len(masses),3)
        self.assertTrue(np.all(np.diff(masses)<=0))

    def test_recompute_and_negative_control(self):
        from types import SimpleNamespace
        V=np.diag([np.exp(.4j),np.exp(.7j),np.exp(1.1j)])
        y1=np.ones((3,3),complex)@V.conj()
        ys=np.ones((3,3),complex)@V
        boundary={'model':'B_-1','mu_gev':1e5,'mf':np.array([5000.,6000.,7000.]),
          'y1':y1,'y2':np.eye(3, dtype=complex),
          'y1_same_V_alternative':ys,'heavy_basis':V,
          'scalar_basis':SimpleNamespace(eigenvalues_sq=np.array([100.,400.]),U=np.eye(2)),
          'quartics':{'lambdaT3':1.,'lambdaT3ExtraQ020020I1':0.,'lambdaT3ExtraQ020200I1':0.}}
        def coeff(mf,mx,my,U,q):return 1e-12,2e-12,3e-12
        from Numerical.experimental.T3BCRunningOneScaleTriangle import assemble_ordered_c5,packed
        c5,_=assemble_ordered_c5(boundary['mf'],y1,boundary['y2'],10.,20.,np.eye(2),boundary['quartics'],coeff)
        archived={'candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL':packed(c5)}
        row=m.evaluate(boundary,coeff,246.,archived)
        self.assertLess(row['recomputed_production_vs_archived_C5_relative'],1e-12)
        self.assertGreater(row['relative_C5_same_U_minus_production'],.01)
        self.assertFalse(row['physical_C5_certified'])

    def test_real_takagi_rotation_is_valid_zero_sensitivity(self):
        from types import SimpleNamespace
        from Numerical.experimental.T3BCRunningOneScaleTriangle import assemble_ordered_c5,packed
        V=np.eye(3,dtype=complex)
        y1=np.array([[1.,2.,0.],[0.,1.,3.],[1.,0.,2.]],complex)
        y2=np.eye(3,dtype=complex)
        q={'lambdaT3':1.,'lambdaT3ExtraQ020020I1':0.,
           'lambdaT3ExtraQ020200I1':0.}
        basis=SimpleNamespace(eigenvalues_sq=np.array([100.,400.]),U=np.eye(2))
        boundary={'model':'B_-1','mu_gev':1e5,'mf':np.array([5000.,6000.,7000.]),
          'y1':y1,'y2':y2,'y1_same_V_alternative':y1.copy(),
          'heavy_basis':V,'scalar_basis':basis,'quartics':q}
        def coeff(mf,mx,my,U,quartics):return 1e-12,2e-12,3e-12
        c5,_=assemble_ordered_c5(boundary['mf'],y1,y2,10.,20.,np.eye(2),q,coeff)
        archived={'candidate_ordered_source_C5_3x3_gev_inverse_NOT_PHYSICAL':packed(c5)}
        row=m.evaluate(boundary,coeff,246.,archived)
        self.assertFalse(row['rotation_prescriptions_numerically_distinguishable'])
        self.assertEqual(row['relative_C5_same_U_minus_production'],0.0)
        self.assertEqual(row['relative_neutrino_spectrum_same_U_minus_production'],0.0)
        self.assertFalse(row['physical_C5_certified'])

    def test_gate_rejects_physical_claim(self):
        native={'status':'NativeAbstractChiralChecksPassedPhysicalDictionaryPending',
          'models':{k:{'passed':True,'checks':{'a':True}} for k in ('B_m1','C_m1')},
          'rgbeta_cross_package_dictionary_certified':False}
        action={'status':'MajoranaWeylActionCovariancePassedNativeDictionaryPending',
                'physical_C5_certified':False}
        source={'physical_matchete_rgbeta_mapping_certified':False}
        archived={'physical_C5_certified':False,'full_fermion_X_Y_sequential_RG_transport':False}
        self.assertEqual(m.check_gate(source,archived,native,action)['weyl_action_check'],'Passed')
        with self.assertRaises(ValueError):
            m.check_gate(source,{**archived,'physical_C5_certified':True},native,action)

if __name__=='__main__':unittest.main()
