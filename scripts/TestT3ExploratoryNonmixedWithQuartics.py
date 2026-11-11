"""Small isolated tests for Batch 10 nonmixed source and tensor adapter."""
from pathlib import Path
from unittest.mock import patch
import tempfile
import unittest
import json
import numpy as np

from scripts import RunT3ExploratoryNonmixedWithQuartics as batch


class FakeResponse:
    model = 'A_+0'
    symbols = ('lambdaX',)
    source = 'synthetic'
    mapping_assumption = 'test'
    def __init__(self):
        self.provider_value = None
    def evaluate(self, state):
        self.provider_value = float(state['lambdaX'])
        return np.full((12,)*4, self.provider_value)
    def validate_hhss(self, state, expected, *, rtol=1e-9):
        return 0.


class Batch10Tests(unittest.TestCase):
    def test_hhss_wrapper_delegates_to_inner(self):
        fake = FakeResponse()
        slot={'value':None}
        wrapper=batch.HHSSOverlay(fake,slot)
        self.assertAlmostEqual(wrapper.evaluate({'lambdaX':2.})[0,0,0,0],2.)
        hhss=np.zeros((4,4,8,8))
        self.assertEqual(wrapper.validate_hhss({'lambdaX':1.},hhss),0.)
        self.assertIs(slot['value'],hhss)
        self.assertEqual(wrapper.symbols,('lambdaX',))

    def test_complete_all16_candidate_inventory_required(self):
        with tempfile.TemporaryDirectory() as tmp:
            reg=Path(tmp)/'registry.json'
            candidate=Path(tmp)/'candidate.json'
            reg.write_text(json.dumps({'models':{'A_-4':{}}}),encoding='utf8')
            candidate.write_text(json.dumps({
                'format':'t3-matchete-candidate-cartesian-quartics-v1',
                'models':{'A_-4':{},'B_-1':{}}
            }),encoding='utf8')
            with patch.object(batch,'REGISTRY',reg), patch.object(batch,'CANDIDATES',candidate), \
                 patch.object(batch,'load_registry',return_value={'models':{'A_-4':{}}}):
                with self.assertRaises(ValueError):
                    batch.require_sources()

    def test_zero_extra_boundary_explicit(self):
        state=object()
        with patch.object(batch,'GeneralT3ScalarRGState') as state_cls:
            state_cls.return_value.validated.return_value='validated'
            registry={'models':{'A_+0':{'new_rgbeta_symbols':[
                {'name':'lambdaNew1'}, {'name':'lambdaNew2'}]}}}
            self.assertEqual(batch.extra_boundary('A_+0',state,registry),'validated')
            self.assertEqual(state_cls.call_args.args[1],
                             {'lambdaNew1':0j,'lambdaNew2':0j})

    def test_missing_only_quadratic_is_explicitly_frozen(self):
        with tempfile.TemporaryDirectory() as tmp:
            native=Path(tmp)/'native.json'; effective=Path(tmp)/'effective.json'
            native.write_text(json.dumps({'status':'Success',
                'report_betas':{'lambdaH':'lambdaH'}}),encoding='utf8')
            fake_layout=type('Layout',(),{'names':('lambdaH','muT3GeneralS1Sq')})()
            with patch.object(batch,'pack_general_scalar_state',return_value=(None,fake_layout)), \
                 patch.object(batch.full_beta,'checked_payload',return_value=({'lambdaH':'lambdaH','muT3GeneralS1Sq':'0'},fake_layout.names)) as checked:
                beta,names,frozen=batch.checked_eft1_betas_with_frozen_quadratics(
                    'A_+0',object(),{},native,effective)
            self.assertEqual(frozen,('muT3GeneralS1Sq',))
            self.assertEqual(beta['muT3GeneralS1Sq'],'0')
            self.assertEqual(json.loads(native.read_text())['report_betas'],{'lambdaH':'lambdaH'})
            self.assertEqual(json.loads(effective.read_text())['report_betas']['muT3GeneralS1Sq'],'0')
            self.assertEqual(checked.call_args.args[3],effective)

    def test_missing_quartic_cannot_be_frozen(self):
        with tempfile.TemporaryDirectory() as tmp:
            native=Path(tmp)/'native.json'
            native.write_text(json.dumps({'status':'Success','report_betas':{'lambdaH':'1'}}))
            layout=type('Layout',(),{'names':('lambdaH','lambdaT3GeneralQ123')})()
            with patch.object(batch,'pack_general_scalar_state',return_value=(None,layout)):
                with self.assertRaisesRegex(ValueError,'unexpected missing/extra'):
                    batch.checked_eft1_betas_with_frozen_quadratics('A_+0',object(),{},native,Path(tmp)/'out.json')

    def test_unexpected_native_beta_excess_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            native=Path(tmp)/'native.json'
            native.write_text(json.dumps({'status':'Success','report_betas':{'lambdaH':'1','EXCESS':'1'}}))
            layout=type('Layout',(),{'names':('lambdaH',)})()
            with patch.object(batch,'pack_general_scalar_state',return_value=(None,layout)):
                with self.assertRaisesRegex(ValueError,'unexpected missing/extra'):
                    batch.checked_eft1_betas_with_frozen_quadratics('A_+0',object(),{},native,Path(tmp)/'out.json')


if __name__=='__main__':
    unittest.main()
