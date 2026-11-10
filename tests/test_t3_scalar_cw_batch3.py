"""One-loop CW matching formula and native source inventory regression tests.

These are pure arithmetic/source-format checks; full Windows 16-model data
are exercised by scripts/RunT3ScalarCWBatch3.py.
"""
import importlib.util
import json
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

CW=module('Numerical/parallel_t3/threshold/scalar_cw_batch3.py','CWBatch3')
NATIVE=module('Numerical/parallel_t3/threshold/native_source_inventory_batch3.py','NativeBatch3')


class TestScalarCW(unittest.TestCase):
    def test_divided_difference_degenerate_and_non_degenerate(self):
        for a,b in ((1.e6,1.e6),(1.e6,1.00000001e6),(1.5e6,2.7e6),
                    (4e4,2e7)):
            got=CW.divided_fprime(a,b,1e6)
            expected=2*np.log(a/1e6) if a==b else (CW.fprime(a,1e6)-CW.fprime(b,1e6))/(a-b)
            self.assertAlmostEqual(got,expected,delta=5e-7)

    def test_single_real_mode_closed_form_and_mass(self):
        # Use two heavy components, only the first coupled. This is equivalent
        # to one real field with M2(h)=M2+lambda*h^2/2 in CW potential.
        n=2
        masses=np.diag([2.1e6,3.5e6])
        q=np.zeros((4+n,)*4)
        lam=0.14
        for i in range(4):
            for a,b in ((i,0),(0,i)):
                for perm in ((a,b,4,4),(a,4,b,4),(a,4,4,b),
                             (4,a,b,4),(4,a,4,b),(4,4,a,b)):
                    q[perm]=lam if a==b==i else 0.
        result=CW.one_loop_higgs_potential(masses,q,1000.)
        analytic_mass=CW.fprime(masses[0,0],1e6)*lam/(64*np.pi**2)
        self.assertAlmostEqual(result['higgs_mass2_shift_gev2'][0,0],analytic_mass,places=9)
        analytic_quartic=3*lam**2*np.log(masses[0,0]/1e6)/(32*np.pi**2)
        self.assertAlmostEqual(result['higgs_quartic_shift_tensor'][0,0,0,0],analytic_quartic,places=12)

    def test_singlet_and_channel_pair_interference(self):
        n=4
        masses=np.diag([1.7e6,1.7e6,2.3e6,2.3e6])
        arr=np.zeros((2,4+n,4+n,4+n,4+n))
        for k,g in enumerate((.10,.045)):
            for i in range(4):
                for a in range(n):
                    ids=(i,i,a+4,a+4)
                    from itertools import permutations
                    for perm in set(permutations(ids)):
                        arr[k][perm]=g
        full=arr.sum(axis=0)
        r=CW.one_loop_higgs_potential(masses,full,1000.,channels=arr,names=('portalA','portalB'))
        self.assertLess(r['quartic_singlet_residual'],1e-12)
        self.assertLess(r['quadratic_singlet_residual'],1e-12)
        self.assertAlmostEqual(r['quartic_projected_channel_pairs'].sum(),r['delta_lambdaH_source_convention'],places=12)
        self.assertAlmostEqual(r['quartic_projected_channel_pairs'][0,1],r['quartic_projected_channel_pairs'][1,0],places=13)
        self.assertAlmostEqual(r['delta_mH_squared_source_convention_gev2'],r['higgs_mass2_shift_gev2'][1,1],places=9)

    def test_nonphysical_quartic_rejected(self):
        q=np.zeros((6,)*4)
        q[0,1,4,5]=1
        with self.assertRaises(ValueError):
            CW.one_loop_higgs_potential(np.eye(2)*1e6,q,1000.)

    def test_native_laurent_inventory_not_injected(self):
        kernel='+ (2*hbar*Coupling[lambdaT3, {}, 0]*Coupling[MF, {}, 0]*GammaCC)'
        d=sha256(kernel.encode('utf8')).hexdigest()
        def rec(ch,i):
            return {'chirality':ch,'coefficient_factor':'1 + 1/ϵ + Log[μbar2/MF^2]',
                    'kernel_inputform':kernel,'kernel_sha256':d,
                    'pole_residue_equals_constant_kernel':True,
                    'log_scale_coefficient_equals_pole_residue':True,
                    'source_indices':{'constant':i,'pole':i+1,'scale_log':i+2},
                    'mass_log_inputform':'Log[μbar2/MF^2]'}
        doc={'format':'t3-f-first-source-resolved-pole-log-bridge-v1','model':'A_-4',
             'EFT1_soft_subtraction_performed':False,
             'counterterm_subtraction_performed':False,
             'physical_MSbar_hard_coefficient':None,
             'exact_kernel_groups':2,'PL_groups':1,'PR_groups':1,
             'records':[rec('PL',1),rec('PR',4)]}
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/NATIVE.ROOT_REL/'A_m4'
            p.mkdir(parents=True)
            (p/'trios.json').write_text(json.dumps(doc),encoding='utf-8')
            out=NATIVE.source_inventory(Path(t),'A_-4','A_m4',{'lambdaT3'})
            self.assertEqual(out['num_native_PL_kernels'],1)
            self.assertEqual(out['quartics_with_conditional_dictionary_record'],['lambdaT3'])
            self.assertFalse(out['raw_native_terms_added_to_C5'])

if __name__=='__main__':unittest.main()
