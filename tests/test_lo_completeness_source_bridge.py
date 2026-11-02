"""Algebraic tests of the diagonal-benchmark-only source bridge."""
from __future__ import annotations
from dataclasses import dataclass
import unittest
import numpy as np
from Numerical.lo_certification.completeness.source_bridge import (
    diagonal_fermion_benchmark,split_formal_by_flavour,prepare_conditional_input,SourceBridgeError
)


def encoded_matrix(m):
    return [[{'real':float(m[i,j].real),'imag':float(m[i,j].imag)} for j in range(3)] for i in range(3)]


def running(y1=None,y2=None,mf=None):
    return {'status':'NumericalRenormalisableRunningOnly','model':'A_-4',
            'scales_gev':{'F':1e5,'S':1e3},
            'uv_at_fermion_threshold':{
               'y1':encoded_matrix(np.diag([.1,.2,.3]) if y1 is None else y1),
               'y2':encoded_matrix(np.diag([.4,.5,.6]) if y2 is None else y2),
               'MF':encoded_matrix(np.diag([1e5,1.1e5,1.2e5]) if mf is None else mf)}}


@dataclass
class MockSource:
    model:str
    raw:np.ndarray
    by_source:np.ndarray
    labels:tuple
    source_indices:tuple
    mu_f_gev:float=1e5
    hash:str='sourcehash'


def good_source():
    n=8
    z=np.zeros((3,3,2,2,n,n),complex)
    for f in range(3):
        for i in range(2):
            z[f,f,i,i,0,4]=z[f,f,i,i,4,0]=(1+1j*f)*1e-9
    sources=np.stack([z*.1,z*.9])
    return MockSource('A_-4',z,sources,tuple(['S1']*2+['S2']*6),(1,2))


class SourceBridgeTests(unittest.TestCase):
    def test_verified_three_flavour_decomposition(self):
        s=good_source()
        d=split_formal_by_flavour(s,running(),'runhash')
        self.assertEqual(d.c_by_fermion.shape,(3,3,3,2,2,8,8))
        self.assertLess(np.max(abs(d.c_total-d.c_by_fermion.sum(0))),1e-25)
        for k in range(3):
            for j in range(3):
                if j!=k:self.assertEqual(np.max(abs(d.c_by_fermion[k,j])),0)
        self.assertEqual(d.fermion_masses_gev.tolist(),[1e5,1.1e5,1.2e5])
        self.assertEqual(d.c_total[2,2,1,1,0,4],(1+2j)*1e-9)

    def test_reject_mixed_heavy_flavour(self):
        mf=np.diag([1e5,1.1e5,1.2e5]).astype(complex)
        mf[0,1]=mf[1,0]=50
        with self.assertRaisesRegex(SourceBridgeError,'non-diagonal'):
            diagonal_fermion_benchmark(running(mf=mf),1e5)
        y=np.diag([.1,.2,.3]).astype(complex);y[0,2]=.00001
        with self.assertRaisesRegex(SourceBridgeError,'non-diagonal'):
            diagonal_fermion_benchmark(running(y1=y),1e5)

    def test_reject_off_diagonal_formal_source(self):
        s=good_source()
        s.raw[0,1,0,0,0,4]=1e-9
        s.by_source[0,0,1,0,0,0,4]=1e-9
        with self.assertRaisesRegex(SourceBridgeError,'Off-diagonal'):
            split_formal_by_flavour(s,running(),'hash')

    def test_reject_wrong_scales_and_mass_phases(self):
        r=running();r['scales_gev']['F']=120000
        with self.assertRaisesRegex(SourceBridgeError,'muF'):
            diagonal_fermion_benchmark(r,1e5)
        mf=np.diag([1e5,1.1e5,1.2e5]).astype(complex);mf[0,0]*=1j
        with self.assertRaisesRegex(SourceBridgeError,'Complex'):
            diagonal_fermion_benchmark(running(mf=mf),1e5)

    def test_missing_scalar_bundle_is_not_inferred(self):
        d=split_formal_by_flavour(good_source(),running(),'hash')
        with self.assertRaisesRegex(SourceBridgeError,'Explicit real-component'):
            prepare_conditional_input(d,{'model':'A_-4'},{},'scalarhash')

    def test_quartic_input_shape_guard(self):
        d=split_formal_by_flavour(good_source(),running(),'hash')
        metadata={'model':'A_-4','tensor_convention':'T3_REAL_LLSS_HHSS_WEINBERG_v1',
          'scheme':'MSbar','source_scalar_real_labels':list(d.labels),
          'complete_scalar_HHSS_channels_supplied':True,
          'provenance':{'source':'unit test'},'mu_s_gev':1e3,
          'quartic_channel_names':['a']}
        with self.assertRaisesRegex(SourceBridgeError,'Explicit mass matrix'):
            prepare_conditional_input(d,metadata,{},'hash')

if __name__=='__main__': unittest.main()
