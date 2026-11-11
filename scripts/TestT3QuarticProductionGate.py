"""Synthetic evidence tests: successful binding and all critical failure gates."""
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from Numerical.running.T3QuarticProductionGate import audit_gate
from Numerical.running.T3RunningQuarticTensor import QuarticEvidenceError
import numpy as np


def payload():
    n = 12
    a = np.zeros((n,)*4)
    b = np.zeros((n,)*4)
    a[4,4,4,4] = 24.
    return dict(format='t3-running-real-quartic-response-v1',
                model='A_-4',
                convention='V4=lambda_abcd phi_a phi_b phi_c phi_d/4!',
                basis='fixed-original-Cartesian',source='synthetic-test-only',
                constant=np.zeros((n,)*4).tolist(),
                responses=[
                    {'symbol':'lambdaS1','component':'real','tensor':a.tolist()},
                    {'symbol':'lambdaH1','component':'real','tensor':b.tolist()}])


def save(p, data):
    p.write_text(json.dumps(data), encoding='utf8')
    return hashlib.sha256(p.read_bytes()).hexdigest()


class TestProductionGate(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d=Path(self.tmp.name)
        self.registry=d/'reg.json'
        self.stage=d/'stage.json'
        self.response=d/'response.json'
        self.binding=d/'binding.json'
        self.model='A_-4'
        self.reg={'models':{'A_-4':{'new_rgbeta_symbols':[],
                                  'preexisting_rgbeta_symbols':[],
                                  'total_additional_quartics':0}}}
        self.cert={'normalisation_certified':True,
                   'certification_status':'ExactTensorAndExplicitPrefactorAlgebraPassed',
                   'certification':{'models':{'A_-4':{
                       'status':'ExactAlgebraPassed','slots':0,'blocks':[]}}}}
        self.resp=payload()
        self._write()

    def _write(self):
        h1=save(self.registry,self.reg)
        h2=save(self.stage,self.cert)
        h3=save(self.response,self.resp)
        self.bind={'format':'t3-quartic-response-certificate-binding-v1',
                   'model':'A_-4','status':'IndependentlyVerifiedResponseMapping',
                   'registry_sha256':h1,'stage22_sha256':h2,
                   'response_sha256':h3,'checks':{
                       'field_component_alignment':True,
                       'potential_prefactors':True,
                       'rgbeta_to_real_cartesian':True,
                       'coupling_direction_completeness':True}}
        save(self.binding,self.bind)

    def check(self):
        return audit_gate(self.registry,self.stage,self.response,self.binding,
                          self.model,{'lambdaS1','lambdaH1'})

    def test_bound_synthetic_evidence(self):
        self.assertEqual(self.check()['coupling_count'],2)

    def test_uncertified_stage22(self):
        self.cert['normalisation_certified']=False
        self._write()
        with self.assertRaises(QuarticEvidenceError): self.check()

    def test_response_tampering(self):
        self.response.write_bytes(self.response.read_bytes()+b' ')
        with self.assertRaises(QuarticEvidenceError): self.check()

    def test_registry_tampering(self):
        self.registry.write_bytes(self.registry.read_bytes()+b' ')
        with self.assertRaises(QuarticEvidenceError): self.check()

    def test_binding_false_claim(self):
        self.bind['checks']['rgbeta_to_real_cartesian']=False
        save(self.binding,self.bind)
        with self.assertRaises(QuarticEvidenceError): self.check()

    def test_missing_quartic(self):
        self.reg['models']['A_-4']['new_rgbeta_symbols']=[{'name':'lambdaExtra'}]
        self.reg['models']['A_-4']['total_additional_quartics']=1
        self.cert['certification']['models']['A_-4']['slots']=1
        self._write()
        with self.assertRaises(QuarticEvidenceError): self.check()


if __name__ == '__main__': unittest.main()
