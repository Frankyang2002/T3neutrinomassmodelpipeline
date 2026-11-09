"""Basic Workstream A algebra, Ward, numerics and provenance gates.
Run from repo root: python -m unittest discover -s tests -p test_t3_parallel_a.py -v
"""
from __future__ import annotations
import hashlib,json,math,tempfile,unittest
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.basis import ALL_MODELS,build_basis
from Numerical.parallel_t3.wilson.coordinates import (
    embed_c5,project_c5,symmetrize,symmetry_residual,embed_formal_llss)
from Numerical.parallel_t3.wilson.adm import beta,total_beta,gauge_ward_residual
from Numerical.parallel_t3.wilson.workflow import (
    synthetic_case,run_wilson,write_results,run_from_files)
from Numerical.parallel_t3.wilson.interchange import (
    _source_check,_base_check,BlockedInput,BASIS_ID,sha256)
from Numerical.parallel_t3.wilson.transforms import (
    wilson_to_canonical,quartic_to_canonical,scalar_mass_to_canonical)


class TestT3ParallelA(unittest.TestCase):
    def test_all16_basis(self):
        self.assertEqual(len(ALL_MODELS),16)
        for model in ALL_MODELS:
            b=build_basis(model)
            self.assertEqual(b.nf,45)
            self.assertEqual(b.ns,4+2*(b.dims[0]+b.dims[1]))
            self.assertEqual(len(set(b.fermion_labels)),45)
            self.assertEqual(len(set(b.scalar_labels)),b.ns)
            g=b.generators()
            for name,n in (('SU3',8),('SU2',3),('U1',1)):
                self.assertEqual(len(g[name]),n)
                for tf,ts in g[name]:
                    self.assertLess(np.max(abs(tf-tf.conj().T)),1e-12)
                    self.assertLess(np.max(abs(ts-ts.conj().T)),1e-12)

    def test_all16_synthetic_invariant_sources_and_portals(self):
        from Numerical.parallel_t3.wilson.synthetic_invariants import (
            synthetic_llss,synthetic_portal)
        for model in ALL_MODELS:
            b=build_basis(model)
            mat=np.diag([1+1j,2,3-1j])*1e-10
            c=synthetic_llss(b,mat)
            lam=synthetic_portal(b,.02)
            self.assertGreater(np.linalg.norm(c[:,:,4:,4:]),0)
            self.assertGreater(np.linalg.norm(lam[:4,:4,4:,4:]),0)
            self.assertLess(symmetry_residual(c),1e-21)
            for group in b.generators().values():
                for tf,ts in group:
                    self.assertLess(gauge_ward_residual(c,tf,ts),1e-11)
                    v=(np.einsum('ca,cbde->abde',ts,lam,optimize=True)
                       +np.einsum('cb,acde->abde',ts,lam,optimize=True)
                       +np.einsum('cd,abce->abde',ts,lam,optimize=True)
                       +np.einsum('ce,abdc->abde',ts,lam,optimize=True))
                    self.assertLess(np.max(abs(v)),1e-10)

    def test_weinberg_full_embedding_and_ward(self):
        b=build_basis('A_-4')
        m=np.array([[1+2j,3-1j,4],[3-1j,-2j,5j],[4,5j,-1]],complex)*1e-12
        c=embed_c5(m,b)
        found, resid=project_c5(c,b)
        self.assertTrue(np.allclose(found,m,rtol=1e-12,atol=1e-25))
        self.assertLess(resid,1e-25)
        self.assertEqual(symmetry_residual(c),0)
        for sector in b.generators().values():
            for tf,ts in sector:
                self.assertLess(gauge_ward_residual(c,tf,ts),1e-12)

    def test_formal_source_transform_and_required_factor(self):
        b=build_basis('A_-4')
        raw=np.zeros((3,3,2,2,b.ns-4,b.ns-4),complex)
        raw[0,0,0,0,0,0]=3+2j
        c=embed_formal_llss(raw,b,np.eye(b.ns-4),-2.)
        self.assertEqual(symmetry_residual(c),0)
        self.assertAlmostEqual(c[b.leptons[0][0],b.leptons[0][0],4,4],-6-4j)
        with self.assertRaises(ValueError):
            embed_formal_llss(raw,b,np.ones((b.ns-4,b.ns-4)),1.)

    def test_explicit_basis_maps(self):
        U=np.array([[0,1],[1,0]],complex)
        R=np.array([[0,1],[1,0]],float)
        c=np.zeros((2,2,2,2),complex);c[0,0,0,0]=1+2j
        d=wilson_to_canonical(c,U,R)
        self.assertEqual(d[1,1,1,1],1+2j)
        m=np.diag([2.,5.]);out=scalar_mass_to_canonical(m,R)
        self.assertTrue(np.allclose(out,np.diag([5.,2.])))
        q=np.zeros((2,2,2,2));q[0,0,0,0]=4
        self.assertEqual(quartic_to_canonical(q,R)[1,1,1,1],4)

    def test_reallinear_and_not_complex_linear(self):
        rng=np.random.default_rng(21);nf=3;ns=2
        y=rng.normal(size=(nf,nf,ns))+1j*rng.normal(size=(nf,nf,ns))
        y=(y+y.swapaxes(0,1))/2
        a=rng.normal(size=(nf,nf,ns,ns))+1j*rng.normal(size=(nf,nf,ns,ns))
        a=symmetrize(a)
        b=symmetrize(rng.normal(size=(nf,nf,ns,ns))+1j*rng.normal(size=(nf,nf,ns,ns)))
        lam=np.zeros((ns,)*4)
        evaluate=lambda x:total_beta(beta(x,y,lam,variant='fonseca_eq22_six_term',prefer_repo=False)[0])
        lhs=evaluate(a+2*b);rhs=evaluate(a)+2*evaluate(b)
        self.assertTrue(np.allclose(lhs,rhs,rtol=1e-12,atol=1e-12))
        self.assertGreater(np.linalg.norm(evaluate(1j*a)-1j*evaluate(a)),1e-5)

    def test_unitary_basis_covariance_rejects_unstarred_crossed(self):
        rng=np.random.default_rng(7);nf,ns=3,2
        Y=rng.normal(size=(nf,nf,ns))+1j*rng.normal(size=(nf,nf,ns))
        Y=(Y+Y.swapaxes(0,1))/2
        C=symmetrize(rng.normal(size=(nf,nf,ns,ns))+1j*rng.normal(size=(nf,nf,ns,ns)))
        U,_=np.linalg.qr(rng.normal(size=(nf,nf))+1j*rng.normal(size=(nf,nf)))
        T=rng.normal(size=(nf,nf))+1j*rng.normal(size=(nf,nf))
        T=(T+T.conj().T)/2
        theta=np.array([[0,-1j],[1j,0]])
        rotate=lambda x:np.einsum('pi,qj,pqab->ijab',U,U,x,optimize=True)
        Y2=np.einsum('pi,qj,pqa->ija',U,U,Y,optimize=True)
        g=((.4,((T,theta),)),)
        grot=((.4,((U.conj().T@T@U,theta),)),)
        lam=np.zeros((ns,)*4)
        errors={}
        for variant in ('fonseca_eq22_six_term','ward_covariant_candidate','source_v1_unstarred'):
            b=total_beta(beta(C,Y,lam,g,variant=variant,prefer_repo=False)[0])
            b2=total_beta(beta(rotate(C),Y2,lam,grot,variant=variant,prefer_repo=False)[0])
            errors[variant]=np.linalg.norm(b2-rotate(b))/np.linalg.norm(rotate(b))
        self.assertLess(errors['fonseca_eq22_six_term'],1e-12)
        self.assertLess(errors['ward_covariant_candidate'],1e-12)
        self.assertGreater(errors['source_v1_unstarred'],1e-2)

    def test_ward_on_candidate_full_beta(self):
        basis,c,bg,_=synthetic_case('A_-4')
        gs=((0.61,bg.generators['SU2']),(0.36,bg.generators['U1']))
        channels,_=beta(c,np.zeros_like(bg.y[0]),bg.quartic[0],gs,
                        variant='fonseca_eq22_six_term',prefer_repo=False)
        bt=total_beta(channels)
        self.assertLess(symmetry_residual(bt),1e-20)
        for _,g in gs:
            for tf,ts in g:
                self.assertLess(gauge_ward_residual(bt,tf,ts),1e-10)

    def test_Aminus4_finite_full_coupled_C5_and_reverse(self):
        basis,full,bg,_=synthetic_case('A_-4')
        self.assertGreater(np.linalg.norm(project_c5(full,basis)[0]),0)
        heavy=full-embed_c5(project_c5(full,basis)[0],basis)
        # C_LLSS is an exact SU2xU1 invariant and mixes into C5 via quartic.
        zero_hard=run_wilson(basis,heavy,bg,1e5,1e3,steps=2)
        self.assertGreater(np.linalg.norm(zero_hard['C5_final']),1e-15)
        hard=run_wilson(basis,full,bg,1e5,1e3,steps=2)
        self.assertFalse(hard['metadata']['physical_certified'])
        self.assertGreater(np.linalg.norm(hard['C5_initial']),0)
        self.assertLess(hard['metadata']['fermion_pair_symmetry_residual'],1e-21)
        just_hard=embed_c5(project_c5(full,basis)[0],basis)
        term,_=beta(just_hard,bg.y[0],bg.quartic[0],variant='fonseca_eq22_six_term',prefer_repo=False)
        # Weinberg can source heavy LLSS through the same quartic tensor.
        self.assertGreater(np.linalg.norm(total_beta(term)[:,:,4:,4:]),0)

    def test_provenance_and_basis_fail_closed(self):
        b=build_basis('A_-4')
        with self.assertRaisesRegex(BlockedInput,'BasisIdentifier'):
            _base_check({'model':b.model,'basis_id':'Matchete',
                         'fermion_labels':list(b.fermion_labels),
                         'scalar_labels':list(b.scalar_labels)},b)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);src=p/'evidence.txt';src.write_text('x')
            artifact=p/'b.npz';artifact.write_bytes(b'bytes')
            meta={'npz_sha256':hashlib.sha256(b'bytes').hexdigest(),
                  'source_sha256':{'evidence.txt':'0'*64}}
            with self.assertRaisesRegex(BlockedInput,'SourceHashMismatch'):
                _source_check(meta,p,artifact)

    def test_real_B_to_A_contract_files_conditional(self):
        # Generate actual disk-backed B-style inputs with corrected native-source
        # provenance; numbers deliberately remain a synthetic formal fixture.
        basis,full,bg,p=synthetic_case('A_-4')
        hard=project_c5(full,basis)[0]
        heavy=full-embed_c5(hard,basis)
        raw=np.zeros((3,3,2,2,basis.ns-4,basis.ns-4),complex)
        for pi in range(3):
            for qi in range(3):
                for i in range(2):
                    for j in range(2):
                        fi=basis.leptons[pi][i]; fj=basis.leptons[qi][j]
                        raw[pi,qi,i,j]=heavy[fi,fj,4:,4:]
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            portal=root/'Lagrangian/interactions/ScalarPotential.wl'
            portal.parent.mkdir(parents=True)
            portal.write_text('T3HiggsPortalHermitianPhase phase coupling[] fieldProduct phase === I PlusHc')
            native=root/'output/general_observable/fermion_first_hermitian_portal_v1/A_m4'/'native_evidence.json'
            native.parent.mkdir(parents=True)
            native.write_text('{"toy_normalization_only": true}')
            source={'Lagrangian/interactions/ScalarPotential.wl':sha256(portal),
                    native.relative_to(root).as_posix():sha256(native)}
            shared={'model':'A_-4','basis_id':BASIS_ID,
                    'fermion_labels':list(basis.fermion_labels),
                    'scalar_labels':list(basis.scalar_labels),
                    'full_SM_three_generation_Weyl_basis':True,
                    'MF_gev':p['MF_gev'],'MS_gev':p['MS_gev'],
                    'source_sha256':source}
            a=root/'b_background.npz'
            np.savez_compressed(a,mu_gev=np.exp(bg.t),Y=bg.y,lambda4=bg.quartic,
                                g2=bg.gauge_couplings[:,0],gY=bg.gauge_couplings[:,1],
                                g3=bg.gauge_couplings[:,2])
            am=root/'b_background.json'
            am.write_text(json.dumps({**shared,'npz_sha256':sha256(a)}))
            b=root/'b_formal.npz'
            np.savez_compressed(b,C_LLSS_formal_muF=raw,
                                heavy_source_from_canonical=np.eye(basis.ns-4),
                                C5_hard=hard)
            bm=root/'b_formal.json'
            mb={**shared,'npz_sha256':sha256(b),'operator_basis':'psi2phi2_no_factorial',
                'native_source_path':native.parent.relative_to(root).as_posix(),
                'corrected_Aminus4_portal_evidence':True,
                'gamma_cc_to_weyl_signed_factor':1.0}
            bm.write_text(json.dumps(mb))
            out=root/'a_result'
            status=run_from_files(root,'A_-4',a,am,b,bm,out,steps=1)
            self.assertEqual(status['status'],'Candidate_Qualitative')
            self.assertTrue((out/'wilson_at_MS.npz').is_file())
            self.assertEqual(status['provenance']['input_boundary']['origin'],
                             'formal_LLSS_EXPLICIT_ORTHOGONAL_INPUT')
            self.assertGreater(status['C5_at_MF_norm'],0)
            # A stale or missing corrected-native source must not be accepted.
            mb['native_source_path']='output/eft1_stageB_saved_a_m4/initial_ode'
            bm.write_text(json.dumps(mb))
            with self.assertRaisesRegex(BlockedInput,'StaleAminus4FormalBoundary'):
                run_from_files(root,'A_-4',a,am,b,bm,root/'rejected',steps=1)

    def test_save_metadata_and_C5_files(self):
        b,c,bg,p=synthetic_case('A_-4')
        r=run_wilson(b,c,bg,p['MF_gev'],p['MS_gev'],steps=1)
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'A_m4';manifest=write_results(r,b,path,{'synthetic':True})
            self.assertEqual(manifest['status'],'Candidate_Qualitative')
            self.assertFalse(manifest['physical_certified'])
            with np.load(path/'wilson_at_MS.npz',allow_pickle=False) as f:
                self.assertEqual(f['C_full'].shape,(45,45,12,12))
                self.assertEqual(f['C5_EFT1_MS'].shape,(3,3))
            with self.assertRaisesRegex(BlockedInput,'OutputDirectoryNonempty'):
                write_results(r,b,path,{})


if __name__=='__main__':unittest.main()
