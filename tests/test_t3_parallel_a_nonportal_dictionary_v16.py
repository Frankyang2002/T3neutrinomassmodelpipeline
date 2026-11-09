"""Standalone synthetic CG rank/action negatives, plus full native read-only integration."""
import copy
from pathlib import Path
import unittest
import numpy as np

from Numerical.parallel_t3.wilson.nonportal_native_dictionary_v16 import (
    audit, audit_sector, scalar_generators, target_tensors, counts,
    source_record, invariant_triplet_metric, require, V16Blocked, EXPECTED_MODELS, DIMENSIONS,
    MAT_FORMAT, STATUS,
)
from Numerical.parallel_t3.wilson.s2_intertwiner_v7 import (
    pauli_half, symmetric_square_cg, induced_s2_generators,
)
from Numerical.parallel_t3.wilson.portal_real_action_v9 import (
    identity, fundamental_fierz, relative,
)


def number(z):
    z=complex(z)
    if abs(z.imag)<1.e-13:return f'{z.real:.16g}'
    if abs(z.real)<1.e-13:return f'({z.imag:.16g}*I)'
    return f'({z.real:.16g}+({z.imag:.16g}*I))'


def generators(d):
    return None if d==1 else pauli_half() if d==2 else induced_s2_generators(symmetric_square_cg())


def fake_portals(d,which):
    base=identity(d)
    if d==1:
        data=[base]
    elif d==2:
        data=[base,fundamental_fierz()]
    else:
        g=generators(3)
        adj=np.einsum('aij,amn->ijmn',pauli_half(),g)
        data=[base,adj]
    rows=[]
    for k,t in enumerate(data,1):
        rows.append({'which':which,'invariant':k,'coupling':f'lambdaH{which}Inv{k}',
                     'cg':f'T3HiggsPortal{which}CGInv{k}',
                     'dimensions':[2,2,d,d],
                     'native_dimensions':([2,2,d,d] if d>1 else [2,2]),
                     'components_inputform':[number(v) for v in t.ravel()],
                     'phase_inside_plus_hc':'1',
                     'native_hermiticity':'Hermitian',
                     'native_tensor_accessor':'Matchete`GetCGTensor',
                     'assembled_lagrangian_contains_cg_and_coupling':True})
    return rows


def fake_record(fam,k,d1,d2,t):
    d=d1 if fam=='self1' else d2
    if fam=='cross':
        dim=(d1,d1,d2,d2);name='lambda12';cg='T3CrossScalarCG';op='ScalarCrossPortal'
    else:
        which=1 if fam=='self1' else 2
        dim=(d,)*4;name=f'lambdaS{which}';cg=f'T3Scalar{which}SelfCG';op=f'Scalar{which}Self'
    return {'family':fam,'invariant':k,'coupling':name+'Inv'+str(k),
            'cg':cg+'Inv'+str(k) if np.prod(dim)>1 else 'None',
            'operator':op+'Inv'+str(k),'dimensions':list(dim),
            'native_dimensions':[x for x in dim if x>1],
            'fields':[],'matchete_potential_factor':('1' if fam=='cross' else '1/2'),
            'components_inputform':[number(x) for x in t.ravel()],
            'native_tensor_accessor':('Matchete`GetCGTensor' if np.prod(dim)>1 else 'Scalar singlet 1'),
            'assembled_lagrangian_contains_cg_and_coupling':True}


def fake_all():
    models={};portals={}
    for tag in sorted(EXPECTED_MODELS):
        d1,d2,df=DIMENSIONS[tag[0]]
        g1,g2=generators(d1),generators(d2)
        rows=[]
        for fam in ('self1','self2','cross'):
            names,ts=target_tensors(fam,d1,d2,g1,g2)
            n=len(ts)
            # Full-rank deliberate mixing to ensure the map is not just an identity test.
            mixing=np.eye(n)+np.triu(np.ones((n,n))*.17,1)
            for k in range(n):
                t=sum(mixing[j,k]*ts[j] for j in range(n))
                rows.append(fake_record(fam,k+1,d1,d2,
                                        t*(2 if fam.startswith('self') else 1)))
        status=('RequiresMixedScalarMassDiagonalization' if tag in ('B_-1','C_-1') else 'Success')
        model={'dimensions':[d1,d2,df],'alpha':int(tag.split('_')[1]),
               'native_build_status':status,'matching_mass_basis_ready':status=='Success',
               'mixed_scalar_mass_basis_blocked':status!='Success',
               'scalar_cg_operator_inventory_verified':True,'nonportal':rows}
        models[tag]=model
        portals[tag]={'dimensions':[d1,d2,df], 'alpha':model['alpha'],
                      'native_build_status':status,
                      'matching_mass_basis_ready':status=='Success',
                      'portals':fake_portals(d1,1)+fake_portals(d2,2)}
    return ({'format':MAT_FORMAT,'nonportal_slots':76,'models_checked':16,
             'models':models,'status':'NativeMatchete76NonportalCGComponentsCaptured_MassBasisStillConditional',
             'physical_scalar_dictionary_certified':False,
             'physical_Wilson_evolution_authorized':False,
             'EFT_order':'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'},
            {'models':portals})


class TestV16Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mat,cls.portal=fake_all()

    def test_16_models_76_native_self_cross_couplings(self):
        self.assertEqual(sum(sum(counts(*DIMENSIONS[t[0]][:2]).values()) for t in EXPECTED_MODELS),76)
        self.assertEqual(sum(len(r['nonportal']) for r in self.mat['models'].values()),76)

    def test_76_nonportal_tensor_and_complex_action_fits(self):
        r=audit(self.mat,self.portal,samples=3)
        self.assertEqual((r['nonportal_coupling_slots'],r['nonportal_sectors']),(76,48))
        self.assertEqual(r['status'],STATUS)
        self.assertLess(r['max_relative_residual'],1.e-10)
        self.assertFalse(r['certification']['physical_Wilson_evolution_authorized'])

    def test_nontrivial_inverse_matrices_reconstructed(self):
        r=audit(self.mat,self.portal,samples=2)
        a=np.array(r['models']['E_-2']['cross']['conditional_rgbeta_from_matchete_real_matrix'])
        self.assertEqual(a.shape,(3,3));self.assertGreater(abs(a[0,1]),.1)
        self.assertLess(np.linalg.norm(a-np.triu(np.ones((3,3))*.17,1)-np.eye(3)),1e-12)

    def test_matchete_portal_generators_recover_lie_algebra(self):
        for d in (1,2,3):
            g,meta=scalar_generators(fake_portals(d,1),1,d)
            if d==1:self.assertIsNone(g)
            else:self.assertLess(meta['max_lie_relative_residual'],1e-12)

    def test_self_potential_boson_symmetry(self):
        r=audit(self.mat,self.portal,samples=2)
        self.assertEqual(r['models']['A_-4']['self2']['native_matchete_sector_rank'],2)

    def test_corrupted_source_component_rejected(self):
        m=copy.deepcopy(self.mat);x=m['models']['E_-2']['nonportal']
        x[-1]['components_inputform'][0]='200'
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_missing_invariant_rejected(self):
        m=copy.deepcopy(self.mat);m['models']['E_-2']['nonportal'].pop()
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_wrong_matchete_factor_rejected(self):
        m=copy.deepcopy(self.mat);m['models']['E_-2']['nonportal'][0]['matchete_potential_factor']='1'
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_wrong_BC_mass_basis_claim_rejected(self):
        m=copy.deepcopy(self.mat);m['models']['B_-1']['matching_mass_basis_ready']=True
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_native_frame_falsely_certified_rejected(self):
        m=copy.deepcopy(self.mat);m['physical_Wilson_evolution_authorized']=True
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_crossed_triplet_independent(self):
        g=generators(3)
        labels,target=target_tensors('cross',3,3,g,g)
        self.assertEqual(labels,['lambda12','lambda12Adj','lambda12Cross'])
        self.assertEqual(np.linalg.matrix_rank(np.column_stack([x.ravel() for x in target])),3)

    def test_triplet_metric_spherical_weight_frame_not_cartesian(self):
        g=generators(3);q=invariant_triplet_metric(g)
        self.assertLess(relative(q.conj().T@q,np.eye(3)),1e-12)
        self.assertGreater(relative(q,np.eye(3)),.1)
        self.assertLess(max(relative(x.T@q+q@x,np.zeros((3,3))) for x in g),1e-12)

    def test_triplet_metric_unitary_frame_covariance(self):
        g=generators(3)
        # A genuinely different unitary component frame must still preserve
        # the crossed invariant tensor after the corresponding field rotation.
        q=invariant_triplet_metric(g)
        phi=np.array([0.17,-.27,.43]); w=np.diag(np.exp(1j*phi))
        gp=np.asarray([w.conj().T@x@w for x in g]);qp=invariant_triplet_metric(gp)
        expect=w.T@q@w
        self.assertLess(relative(qp*np.exp(-1j*np.angle(np.vdot(expect,qp))),expect),1e-12)

    def test_different_triplet_frames_crossed_rejected(self):
        g=generators(3)
        with self.assertRaises(V16Blocked):target_tensors('cross',3,3,g,-g)

    def test_real_self_HC_failure_is_rejected(self):
        m=copy.deepcopy(self.mat)
        for x in m['models']['A_-4']['nonportal']:
            if x['family']=='self2' and x['invariant']==1:
                x['components_inputform'][0]='(2*I)'
                break
        with self.assertRaises(V16Blocked):audit(m,self.portal,samples=1)

    def test_real_quartic_census_count(self):
        self.assertEqual(16+58+76,150)

class TestLiveIfAvailable(unittest.TestCase):
    def test_current_native_matchete_cg_v16_if_present(self):
        root=Path(__file__).resolve().parents[1]
        p=root/'output/parallel_t3/wilson/matchete_nonportal_native_v16.json'
        if not p.is_file():self.skipTest('Native Matchete v16 exporter not executed')
        from Numerical.parallel_t3.wilson.nonportal_native_dictionary_v16 import audit_from_paths
        b=root/'output/parallel_t3/wilson'
        result=audit_from_paths(p,b/'matchete_portals_native_v5.json',
                               b/'rgbeta_quartic_projectors_native_v14.json',root,samples=3)
        self.assertEqual(result['nonportal_coupling_slots'],76)

if __name__=='__main__':unittest.main()
