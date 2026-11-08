"""Targeted A(-4) native-tree and complex 3x3 flavour reconstruction tests."""
from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
import numpy as np

from Numerical.lo_certification.spinor import a_minus4_boundary_verify as a4

CG = ('SparseArray[Automatic, {2, 2, 3}, 0, {1, '
      '{{0, 3, 6}, {{1, 1}, {1, 3}, {2, 2}, {1, 2}, {2, 1}, {2, 3}}}, '
      '{3^(-1/4), I/3^(1/4), -3^(-1/4), -3^(-1/4), -3^(-1/4), I/3^(1/4)}}]')
PL = ('(Bar[Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]]*'
      'Bar[Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]]*'
      'Bar[Field[NewScalar1, Scalar, {}, {}]]*'
      'CG[T3MixCG, {Bar[Index[d$$1, SU2L[fund]]], Bar[Index[d$$2, SU2L[fund]]], '
      'Index[d$$1, T3BSMd3]}]*'
      'Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]*'
      'NCM[Transp[Field[l, Fermion, {Index[d$$2, SU2L[fund]], '
      'Index[d$$2, Flavor]}, {}]], DiracProduct[GammaCC, Proj[-1]], '
      'Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]])'
      '/Coupling[MF, {Index[d$$1, NFlavor]}, 0]')
PR = ('(Bar[Field[NewScalar2, Scalar, {Index[d$$1, T3BSMd3]}, {}]]*'
      'CG[Bar[T3MixCG], {Index[d$$1, SU2L[fund]], Index[d$$2, SU2L[fund]], '
      'Index[d$$1, T3BSMd3]}]*'
      'Coupling[y1, {Index[d$$1, Flavor], Index[d$$1, NFlavor]}, 0]*'
      'Coupling[y2, {Index[d$$2, Flavor], Index[d$$1, NFlavor]}, 0]*'
      'Field[NewScalar1, Scalar, {}, {}]*'
      'NCM[Bar[Field[l, Fermion, {Index[d$$1, SU2L[fund]], Index[d$$1, Flavor]}, {}]], '
      'DiracProduct[GammaCC, Proj[1]], '
      'Transp[Bar[Field[l, Fermion, {Index[d$$2, SU2L[fund]], '
      'Index[d$$2, Flavor]}, {}]]]])'
      '/Coupling[MF, {Index[d$$1, NFlavor]}, 0]')
REG = [{'Name':'T3MixCG',
        'RepsInputForm':'{Bar[SU2L[fund]], Bar[SU2L[fund]], T3BSMd3}',
        'TensorInputForm':CG}]


def make_fixture(root: Path):
    masses = np.array([1200.,2000.,4500.])
    y1 = np.array([[1+2j, 2-1j, 3+0.2j],
                   [-3j, .3+.1j, -1+.5j], [3-2j, 4+.2j, 2j]])
    y2 = np.array([[.7+2j, 1j, 2+1j],[-2+.3j, 3j, 4], [1-2j, .5, 1+3j]])
    w=a4.dense_flavour(y1,y2,masses)
    cg,_=a4.parse_native_cg(REG)
    original=a4.raw_components(w,cg)
    independent,_=a4.symmetrise(original)
    seed={'Level':1,'IntegratedFields':['F'],'ActiveHeavyFields':['S1','S2'],
          'TreeWilsonTermCount':2,
          'TreeWilsonTerms':[{'Index':1,'Chirality':'PR','TermInputForm':PR,
                              'Scalar1Count':1,'Scalar2Count':1,'CGNames':['T3MixCG']},
                             {'Index':2,'Chirality':'PL','TermInputForm':PL,
                              'Scalar1Count':1,'Scalar2Count':1,'CGNames':['T3MixCG']}],
          'CGRegistry':REG}
    source={'model':'A_-4','mu_F_gev':1000.,'heavy_masses_at_F_gev':masses.tolist(),
            'operator_terms':[{'source_index':2,
                               'source_term_inputform':PL,
                               'flavour_matrix_axes_external_indices':['d$$1','d$$2'],
                               'factors':[{'name':'y1','conjugated':True,'external_index':'d$$1'},
                                          {'name':'y2','conjugated':True,'external_index':'d$$2'}],
                               'ordered_flavour_matrix_gev_inv':{'real':w.real.tolist(),'imag':w.imag.tolist()}}],
            'source_CG_registry':REG}
    components={}
    for (i,q,j,p,a,b),z in original.items():
        mat=components.setdefault((i,j,a,b),np.zeros((3,3),complex))
        mat[p-1,q-1] += z
    ordered={'model':'A_-4','operators':[{'source_index':2,
             'flavour_matrix_axes_external_indices':['d$$1','d$$2'],
             'components':[{'lepton_su2':[i,j],'scalar_real':[a,b],
                            'ordered_flavour_coefficient_gev_inv':
                                {'real':m.real.tolist(),'imag':m.imag.tolist()}}
                            for (i,j,a,b),m in sorted(components.items())]}]}
    poly={'model':'A_-4','independent_symmetric_coordinates':[]}
    for (i,j,a,b),z in sorted(independent.items()):
        def invert(k): return (k-1)%2+1, (k-1)//2+1
        k1,p=invert(i);k2,q=invert(j)
        poly['independent_symmetric_coordinates'].append({
            'coordinate':[k1,p,k2,q,a,b],
            'coefficient_real':float(z.real),'coefficient_imag':float(z.imag),
            'orbit_size':len(a4.orbit(i,j,a,b))})
    items = {a4.RAW/'eft1_after_F_wilson_seed.json':seed,
             a4.SOURCE:source,a4.ORDERED:ordered,a4.POLY:poly}
    for path,obj in items.items():
        dest=root/path;dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(obj),encoding='utf-8')
    # Native full tree must have the original matched terms, not only a report.
    dest=root/a4.RAW/'eft1_after_F_tree_inputform.txt'
    dest.write_text(PR+' + '+PL,encoding='utf-8')
    return w,independent,source,ordered,poly


class TestAMinus4SourceBoundary(unittest.TestCase):
    def test_cg_exact_complex_phases_and_symmetric_triplet(self):
        g,proof=a4.parse_native_cg(REG)
        v=3**(-.25)
        self.assertAlmostEqual(g[0,0,0].real,v)
        self.assertAlmostEqual(g[0,0,2].imag,v)
        self.assertAlmostEqual(g[0,1,1].real,-v)
        self.assertAlmostEqual(g[1,1,0].real,-v)
        self.assertEqual(proof['native_triplet_rank'],3)
        invalid=json.loads(json.dumps(REG)); invalid[0]['TensorInputForm']=CG.replace('I/3^(1/4)','-I/3^(1/4)',1)
        with self.assertRaises(a4.BoundaryBlocked): a4.parse_native_cg(invalid)

    def test_general_complex_three_by_three_and_flavour_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            w,independent,*_ = make_fixture(root)
            self.assertGreater(abs(w[0,1]),0)
            self.assertGreater(abs(w[0,1]-w[1,0]),1e-6)
            obj=a4.verify_from_sources(root,check_running=False)
            self.assertFalse(obj['physical_certification'])
            self.assertIsNone(obj['physical_C_ijab_at_MF'])
            self.assertIsNone(obj['convention']['native_to_Weyl_factor'])
            self.assertEqual(obj['basis_contract']['full_tensor_shape'],[45,45,12,12])
            self.assertEqual(len(obj['independent_conditional_tensor_entries']),len(independent))
            self.assertEqual(obj['independent_component_checks']['native_ordered_component_count'],24)
            self.assertEqual(obj['native_PR_formal_conjugate_check']['PL_PR_scalar_CG_flavour_conjugate_components_checked'],216)
            self.assertFalse(obj['native_PR_formal_conjugate_check']['PR_native_spinor_HC_phase_certified'])
            cg,_=a4.parse_native_cg(REG)
            raw=a4.raw_components(w,cg)
            self.assertAlmostEqual(abs(raw[(1,2,1,1,5,7)]-w[0,1]*cg[0,0,0]/2),0.)
            self.assertEqual(obj['polynomial_reconstruction_max_abs_residual'],0.)

    def test_native_source_change_and_missing_conjugate_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            make_fixture(root)
            path=root/a4.RAW/'eft1_after_F_wilson_seed.json'
            seed=json.loads(path.read_text())
            seed['TreeWilsonTerms'][0]['TermInputForm']=seed['TreeWilsonTerms'][0]['TermInputForm'].replace('Proj[1]','Proj[-1]')
            path.write_text(json.dumps(seed))
            with self.assertRaises(a4.BoundaryBlocked):
                a4.verify_from_sources(root,check_running=False)

    def test_explicit_weyl_C_versus_C_inverse_sign(self):
        # Fixed +--- Weyl gamma convention, epsilon12=+1.
        zero=np.zeros((2,2),complex)
        unit=np.eye(2,dtype=complex)
        sigma2=np.array([[0,-1j],[1j,0]])
        gamma0=np.block([[zero,unit],[unit,zero]])
        gamma2=np.block([[zero,sigma2],[-sigma2,zero]])
        gamma5=np.block([[-unit,zero],[zero,unit]])
        pl=(np.eye(4)-gamma5)/2
        eps=1j*sigma2
        c=1j*gamma2@gamma0
        a=np.array([1+2j, 3-1j])
        b=np.array([-2+.4j, 1+5j])
        a4=np.r_[a, [0j,0j]]
        b4=np.r_[b, [0j,0j]]
        left=a.T@eps@b
        self.assertAlmostEqual(abs(a4.T@c@pl@b4-left),0.)
        self.assertAlmostEqual(abs(a4.T@np.linalg.inv(c)@pl@b4+left),0.)
        self.assertGreater(abs(left),0.1)
        # The Matchete input expression contains only the SYMBOL GammaCC;
        # neither of these two matrix choices is fixed by the raw tree export.

    def test_unmeasured_native_phase_cannot_be_inferred(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            _,independent,*_=make_fixture(root)
            one=dict(independent)
            minus={k:-v for k,v in one.items()}
            self.assertEqual(a4.reconstruct(minus),
                             {k:-v for k,v in a4.reconstruct(one).items()})
            self.assertNotEqual(one,minus)
            # Both candidate signs pass the same formal polynomial identity.
            self.assertFalse(a4.verify_from_sources(root,check_running=False)['physical_certification'])

    def test_bad_ordered_triplet_phase_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            make_fixture(root)
            path=root/a4.ORDERED
            obj=json.loads(path.read_text())
            obj['operators'][0]['components'][0]['ordered_flavour_coefficient_gev_inv']['imag'][0][1] += 0.01
            path.write_text(json.dumps(obj))
            with self.assertRaises(a4.BoundaryBlocked):
                a4.verify_from_sources(root,check_running=False)


if __name__=='__main__':
    unittest.main()
