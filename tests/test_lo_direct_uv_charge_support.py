"""Regression: source T3 loop channels differ from all Yukawa-only channels."""
import json
from pathlib import Path
import unittest
import numpy as np
from Numerical.direct_one_loop.core import encoded,evaluate
from Numerical.lo_certification.direct_uv.certification import all_models,charge_bases,validate_components
from Numerical.lo_certification.direct_uv.channels import partition_exact
from Numerical.lo_certification.direct_uv.runner import verify_source_charge_support
from test_lo_direct_uv_certification import synthetic_all16,_matrix_from_encoded
ROOT=Path(__file__).resolve().parents[1]

def reports(charges):
    return ({'internal_charge_channels':list(charges)},
            {'channels':[{'charge':q,'nonzero':True,'complex_weight':'conditional'} for q in charges],
             'total_group_coefficient':'conditional','verified_physical_loop_multiplicity':None})

def extra(data,short,expected):
    base=charge_bases(short[0],int(short[2:]));qs=base['scalar'];qf=base['weyl'];d1=base['dimensions'][0]
    added=[]
    for index,field,lo,hi,sign in [(0,'conjugate',0,d1,1),(1,'field',d1,len(qs),-1)]:
        v=data['vertices'][index];arr=_matrix_from_encoded(v[field])
        for s in range(lo,hi):
            q=float(qs[s]);f=np.where(np.isclose(qf,sign*q))[0]
            if q not in expected and len(f):
                arr[0,int(f[0]),s]=1e-6;added.append(q)
        v[field]=encoded(arr)
    return set(added)

class NativeChargeRegression(unittest.TestCase):
    def test_all16_keep_yukawa_only_channels(self):
        p=ROOT/'output/general_observable/lo_certification_direct_uv/source_coverage_all16.json'
        entries={e['model']:e for e in json.loads(p.read_text())['rows']}
        self.assertEqual(len(entries),16)
        with_extras=0
        for short,_ in all_models():
            with self.subTest(model=short):
                expected=set(map(float,entries[short]['common_internal_electric_charge_channels']))
                data=synthetic_all16(short);addition=extra(data,short,expected)
                kind,alpha=short[0],int(short[2:]);struct=validate_components(data,kind,alpha)
                result=verify_source_charge_support(struct,data,kind,alpha,*reports(expected))
                self.assertTrue(result['source_expected_channel_coverage_checked'])
                self.assertEqual(set(result['both_yukawa_charge_support']),expected)
                self.assertEqual(set(result['yukawa_only_charge_channels_not_in_T3_contraction']),addition)
                if addition:
                    with_extras+=1
                    loop=evaluate(data)
                    parts=partition_exact(data,kind,alpha,loop)
                    self.assertLess(parts['reconstruction_relative_error'],1e-5)
        self.assertGreaterEqual(with_extras,9)

    def test_mismatched_source_audits_rejected(self):
        data=synthetic_all16('A_-4');s=validate_components(data,'A',-4)
        q,c=reports({-2.});c['channels'][0]['charge']=-1.
        with self.assertRaisesRegex(ValueError,'audits disagree'):
            verify_source_charge_support(s,data,'A',-4,q,c)

    def test_missing_required_yukawa_topology_channel_rejected(self):
        data=synthetic_all16('A_-4');v=data['vertices'][0]['conjugate']
        data['vertices'][0]['conjugate']=encoded(np.zeros_like(_matrix_from_encoded(v)))
        s=validate_components(data,'A',-4)
        with self.assertRaisesRegex(ValueError,'Missing source T3 topology charge'):
            verify_source_charge_support(s,data,'A',-4,*reports({-2.}))

    def test_missing_audit_is_explicitly_not_certified(self):
        data=synthetic_all16('A_-4');s=validate_components(data,'A',-4)
        info=verify_source_charge_support(s,data,'A',-4,None,None)
        self.assertFalse(info['source_expected_channel_coverage_checked'])
        self.assertFalse(info['native_CG_candidate_channel_coverage_checked'])

    def test_full_runner_with_16_native_charge_audit_fixtures(self):
        """Exercise runner's source-CG branch, missed by the older unit tests."""
        import hashlib
        import tempfile
        from contextlib import redirect_stdout
        import io
        from Numerical.lo_certification.direct_uv.runner import run
        from Numerical.lo_certification.direct_uv.certification import CLASSES
        coverage=ROOT/'output/general_observable/lo_certification_direct_uv/source_coverage_all16.json'
        entries={r['model']:r for r in json.loads(coverage.read_text())['rows']}
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            inpdir=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
            auditdir=root/'output/general_observable/lo_direct_one_loop/source_audits'
            inpdir.mkdir(parents=True);auditdir.mkdir(parents=True)
            ch=[];cg=[]
            for short,key in all_models():
                kind,alpha=short[0],int(short[2:]);d1,d2,df,_=CLASSES[kind]
                rep={'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':alpha}
                benchmark={'representation':rep,'base_state':{'mu_gev':1e7,'ordinary':{}},
                           'scales':{'mu_fermion_threshold_gev':1e5,
                                     'mu_scalar_threshold_gev':1e3}}
                bp=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bp.parent.mkdir(parents=True,exist_ok=True)
                bp.write_text(json.dumps(benchmark))
                data=synthetic_all16(short)
                expected=set(map(float,entries[short]['common_internal_electric_charge_channels']))
                extra(data,short,expected)
                data['provenance']={'benchmark_sha256':hashlib.sha256(bp.read_bytes()).hexdigest(),
                                    'validated_uv_component_conventions':False}
                (inpdir/(short+'.json')).write_text(json.dumps(data))
                model={'dS1':d1,'dS2':d2,'dF':df,'alpha':alpha}
                qe,ce=reports(expected)
                ch.append({'model':model,**qe})
                cg.append({'model':model,**ce})
            (auditdir/'charge_channels.json').write_text(json.dumps({'models':ch}))
            (auditdir/'charge_resolved_loop_weights.json').write_text(json.dumps({'models':cg}))
            out=root/'output/general_observable/lo_certification_direct_uv'
            with redirect_stdout(io.StringIO()):summary=run(root,out)
            self.assertEqual(summary['counts']['exact_uv_provisional'],16)
            self.assertEqual(summary['counts']['blocked'],0)
            self.assertEqual(summary['counts']['physical_source_certified'],0)
            self.assertEqual(len(list((out/'models').glob('*.json'))),16)

if __name__=='__main__':unittest.main()
