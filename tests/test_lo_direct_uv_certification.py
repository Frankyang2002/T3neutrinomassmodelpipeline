"""Tests use mathematical fixtures; no source certification is inferred."""
import json
import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.direct_one_loop.core import evaluate, encoded
from Numerical.lo_certification.direct_uv.certification import (
    all_models, charge_bases, witness_quadratic, validate_components,source_evidence,CLASSES)
from Numerical.lo_certification.direct_uv.channels import partition_exact


def _matrix_from_encoded(a):
    return np.asarray(a["real"])+1j*np.asarray(a["imag"])


def toy_B():
    # B alpha=-1 is the neutral Majorana singlet with two scalar doublets.
    m=np.diag([100.,110.,120.]).astype(complex)
    a=np.diag([300**2,320**2,330**2,400**2]).astype(complex)
    a[0,3]=a[3,0]=25.0
    y=np.zeros((3,3,4),complex)
    z=np.zeros_like(y)
    z[:,0,0]=[.03+.01j,.02,-.01j]
    y[:,0,3]=[.01-.02j,.03j,.04]
    return {'model':'B_-1','fermion_mass_gev':encoded(m),
            'scalar_A_gev2':encoded(a),
            'scalar_B_gev2':encoded(np.zeros_like(a)),
            'renormalization_scale_gev':1000.,
            'vertices':[{'name':'y1_HC','conjugate':encoded(z)},
                        {'name':'y2_HC','field':encoded(y)}],
            'provenance':{'validated_uv_component_conventions':False}}


def synthetic_all16(short):
    """Pure mathematical input, NOT actual source-CG data or fitted benchmark."""
    kind,alpha=short[0],int(short[2:])
    base=charge_bases(kind,alpha)
    d1,d2,df=base['dimensions']; qs=base['scalar'];qf=base['weyl']
    nf=len(qf);ns=len(qs);majorana=base['self_conjugate']
    mm=np.zeros((nf,nf),complex)
    if majorana:
        for generation in range(3):
            for a in range(df):
                b=df-1-a
                mm[generation*df+a,generation*df+b]=(100000+10000*generation)
    else:
        nc=3*df
        for generation in range(3):
            for a in range(df):
                mm[generation*df+a,nc+generation*df+a]=(100000+10000*generation)
                mm[nc+generation*df+a,generation*df+a]=(100000+10000*generation)
    A=np.diag(np.array([900.**2+800*s for s in range(d1)]+[1020.**2+900*s for s in range(d2)])).astype(complex)
    allowed=sorted(set(qs[:d1])&set(qs[d1:]))
    y1=np.zeros((3,nf,ns),complex)
    y2=np.zeros_like(y1)
    for q in allowed:
        s1=next(s for s in range(d1) if qs[s]==q)
        s2=next(s for s in range(d1,ns) if qs[s]==q)
        A[s1,s2]=A[s2,s1]=3200.
        f1=next(i for i in range(nf) if qf[i]==q)
        f2=next(i for i in range(nf) if qf[i]==-q)
        for a in range(3):
            y1[a,f1,s1]=(a+1)*.001
            y2[a,f2,s2]=(3-a)*.001
    return {'model':short,'fermion_mass_gev':encoded(mm),'scalar_A_gev2':encoded(A),
            'scalar_B_gev2':encoded(np.zeros_like(A)),'renormalization_scale_gev':1e7,
            'vertices':[{'name':'y1_HC','conjugate':encoded(y1)},
                        {'name':'y2_HC','field':encoded(y2)}]}


class SourceAndChargeChecks(unittest.TestCase):
    def test_all_sixteen(self):
        cases=list(all_models())
        self.assertEqual(len(cases),16)
        self.assertEqual(len({x for x,_ in cases}),16)
        self.assertEqual(len({x for _,x in cases}),16)
        for name,_ in cases:
            c=charge_bases(name[0],int(name[2:]))
            self.assertEqual(len(c['scalar']),sum(c['dimensions'][:2]))
            self.assertEqual(len(c['weyl']),3*c['dimensions'][2]*(1 if c['self_conjugate'] else 2))

    def test_neutral_fermion_majorana_vs_charged_dirac(self):
        self.assertTrue(charge_bases('B',-1)['self_conjugate'])
        self.assertTrue(charge_bases('C',-1)['self_conjugate'])
        self.assertFalse(charge_bases('C',1)['self_conjugate'])
        self.assertEqual(charge_bases('B',-1)['weyl'].tolist(),[0]*3)

    def test_complex_bilinear_hessian(self):
        A=np.array([[5.,.2+.4j],[.2-.4j,10.]],complex)
        B=np.array([[.3+1j, .2-.3j],[.2-.3j,2+.4j]])
        r=witness_quadratic(A,B)
        self.assertLess(r['quadratic_action_max_relative_residual'],1e-13)
        self.assertLess(r['real_hessian_symmetry_residual'],1e-13)

    def test_charge_and_actions(self):
        d=toy_B();r=validate_components(d,'B',-1)
        self.assertTrue(all(r['checks'].values()))
        self.assertEqual(r['observed_scalar_charge_channels'],[0.])
        self.assertEqual(len(r['vertex_dictionary']),2)

    def test_wrong_scalar_vertex_charge_is_rejected(self):
        d=toy_B();bad=np.asarray(d['vertices'][0]['conjugate']['real'])
        bad[0,0,1]=.5
        d['vertices'][0]['conjugate']['real']=bad.tolist()
        with self.assertRaisesRegex(ValueError,'charge'):
            validate_components(d,'B',-1)

    def test_wrong_fermion_mass_charge_is_rejected(self):
        # A, alpha=-2, Dirac charge assignments: [0,+1, 0,-1] each x3.
        c=charge_bases('A',-2)
        nf=len(c['weyl']);ns=len(c['scalar'])
        d={'fermion_mass_gev':encoded(np.eye(nf,dtype=complex)),
           'scalar_A_gev2':encoded(np.eye(ns,dtype=complex)),
           'scalar_B_gev2':encoded(np.zeros((ns,ns),complex)),
           'vertices':[{'name':'dummy','field':encoded(np.zeros((3,nf,ns),complex))}]}
        with self.assertRaisesRegex(ValueError,'charge'):
            validate_components(d,'A',-2)

    def test_source_never_automatically_certified(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);bp=root/'bench.json';bp.write_text('{}')
            inp=toy_B();inp['provenance']['benchmark_sha256']='not-the-source-hash'
            benchmark={'representation':{'d_s1':2,'d_s2':2,'d_f':1,'alpha':-1},
                       'base_state':{'mu_gev':1000,'ordinary':{}}}
            ev=source_evidence(root,'B_-1','x',inp,benchmark,None,None,{'benchmark':bp})
            self.assertFalse(ev['physically_certified'])
            self.assertFalse(ev['checks']['source_benchmark_matches'])
            self.assertGreater(len(ev['blocking_physical_conventions']),2)


class LoopDecompositionChecks(unittest.TestCase):
    def test_full_matrix_and_charge_decomposition(self):
        d=toy_B();r=evaluate(d)
        p=partition_exact(d,'B',-1,r)
        self.assertEqual(p['active_charge_channels'],['0|0'])
        self.assertLess(p['reconstruction_relative_error'],1e-7)
        self.assertTrue(r['diagnostics']['uv_finite'])
        self.assertLess(r['total']['takagi_residual'],1e-9)
        self.assertEqual(len(r['individual_terms']),72)  # 3 vertex pairs x 3 fermions x 8 real scalars

    def test_all16_synthetic_charge_hessian_and_exact_channel_reconstruction(self):
        models=list(all_models())
        tested=0
        for short,_key in models:
            with self.subTest(model=short):
                d=synthetic_all16(short)
                structural=validate_components(d,short[0],int(short[2:]))
                self.assertTrue(all(structural['checks'].values()))
                direct=evaluate(d)
                self.assertTrue(direct['diagnostics']['uv_finite'])
                channels=partition_exact(d,short[0],int(short[2:]),direct)
                self.assertLess(channels['reconstruction_relative_error'],1e-5)
                self.assertEqual(len(direct['total']['neutrino_masses_ev']),3)
                tested+=1
        self.assertEqual(tested,16)

    def test_full_all16_runner_in_separate_synthetic_fixture_repo(self):
        # Integration of output contracts only. NOT source-derived neutrino masses.
        import hashlib
        from Numerical.lo_certification.direct_uv.runner import run
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            inputs=root/'output/general_observable/lo_direct_one_loop/provisional_inputs'
            inputs.mkdir(parents=True)
            for short,key in all_models():
                bpath=root/'Numerical/benchmarks/fitted_reference'/(key+'.json')
                bpath.parent.mkdir(parents=True,exist_ok=True)
                c=charge_bases(short[0],int(short[2:]))
                d1,d2,df=c['dimensions']
                benchmark={'representation':{'d_s1':d1,'d_s2':d2,'d_f':df,'alpha':int(short[2:])},
                    'base_state':{'mu_gev':1e7,'ordinary':{}},
                    'scales':{'mu_fermion_threshold_gev':1e5,
                              'mu_scalar_threshold_gev':1e3}}
                bpath.write_text(json.dumps(benchmark))
                data=synthetic_all16(short)
                data['provenance']={'benchmark_sha256':hashlib.sha256(bpath.read_bytes()).hexdigest(),
                                    'validated_uv_component_conventions':False}
                (inputs/(short+'.json')).write_text(json.dumps(data))
            output=root/'output/general_observable/lo_certification_direct_uv'
            # Keep test output quiet, not misleading supervisor-facing source results.
            from contextlib import redirect_stdout
            import io
            with redirect_stdout(io.StringIO()):result=run(root,output)
            self.assertEqual(result['counts']['catalogue'],16)
            self.assertEqual(result['counts']['exact_uv_provisional'],16)
            self.assertEqual(result['counts']['physical_source_certified'],0)
            self.assertEqual(result['counts']['blocked'],0)
            self.assertEqual(len(list((output/'models').glob('*.json'))),16)
            self.assertTrue((output/'all_16_summary.csv').exists())
            with self.assertRaises(FileExistsError):run(root,output)

    def test_complex_ma_like_one_scalar_exact_B0_and_phase(self):
        from Numerical.direct_one_loop.core import LOOP,b0_finite
        M=100.;m2=300.**2;delta=700.;h=.03+.02j;mu=600.
        A=np.array([[m2]],complex);B=np.array([[delta]],complex)
        field=np.zeros((3,1,1),complex);field[0,0,0]=h
        d={'model':'single_scalar_math_check','fermion_mass_gev':encoded([[M]]),
           'scalar_A_gev2':encoded(A),'scalar_B_gev2':encoded(B),
           'vertices':[{'name':'h','field':encoded(field)}],
           'renormalization_scale_gev':mu}
        out=evaluate(d)
        observed=_matrix_from_encoded(out['total']['complex_symmetric_matrix_gev'])[0,0]
        expected=.5*LOOP*h*h*M*(b0_finite(M*M,m2+delta,mu*mu)
                                   -b0_finite(M*M,m2-delta,mu*mu))
        self.assertAlmostEqual(abs(observed-expected),0,delta=max(abs(expected)*1e-8,1e-16))

    def test_same_point_eft_comparison_guard_and_discrepancy(self):
        from Numerical.direct_one_loop.comparison import compare_same_point
        d=toy_B();d['provenance']['benchmark_sha256']='a'*64
        out=evaluate(d)
        eft={'physical_dimension':'mass_matrix_gev','scheme':'MSbar',
             'benchmark_sha256':'a'*64,'renormalization_scale_gev':1000.,
             'same_renormalized_uv_parameters':True,
             'matrix_gev':out['total']['complex_symmetric_matrix_gev']}
        result=compare_same_point(out,eft)
        self.assertEqual(result['status'],'same_scale_comparison_completed')
        self.assertLess(result['relative_mass_matrix_frobenius_difference'],1e-12)
        other=dict(eft,renormalization_scale_gev=100.)
        self.assertEqual(compare_same_point(out,other)['status'],'not_compared_different_scales')
        with self.assertRaisesRegex(ValueError,'Different'):
            compare_same_point(out,dict(eft,benchmark_sha256='b'*64))

    def test_complex_hessian_mass_and_zero_mixing(self):
        d=toy_B()
        d['scalar_A_gev2']['real'][0][3]=0.
        d['scalar_A_gev2']['real'][3][0]=0.
        r=evaluate(d)
        self.assertLess(np.linalg.norm(r['total']['complex_symmetric_matrix_gev']['real']),1e-12)

    def test_tiny_vertex_retained_even_if_suppressed(self):
        d=toy_B();row={'name':'tiny_suppressed',
            'conjugate':encoded(1e-6*np.array(
                np.asarray(d['vertices'][0]['conjugate']['real'])+1j*np.asarray(d['vertices'][0]['conjugate']['imag'])))}
        d['vertices'].append(row)
        r=evaluate(d)
        p=partition_exact(d,'B',-1,r)
        self.assertIn('tiny_suppressed x tiny_suppressed',p['per_vertex_pair'])
        self.assertTrue(p['per_vertex_pair']['tiny_suppressed x tiny_suppressed']['numerically_suppressed'])
        self.assertEqual(len(r['vertex_pair_matrices_gev']),6)

    def test_serialization_of_full_report(self):
        d=toy_B();r=evaluate(d);p=partition_exact(d,'B',-1,r)
        self.assertTrue(json.loads(json.dumps(p,allow_nan=False)))

if __name__=='__main__':unittest.main()
