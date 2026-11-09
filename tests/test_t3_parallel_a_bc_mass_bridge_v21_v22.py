"""v21/v22 synthetic invariance tests, negative controls, optional live integration.

The random tensors are test fixtures, NOT native Matchete component evidence.
"""
from __future__ import annotations
import copy
import tempfile
import unittest
from pathlib import Path
from unittest import mock
import json
import numpy as np
from Numerical.experimental.GeneralT3MixedScalarMassBasis import EPS, MixedDoubletMassBasis
from Numerical.parallel_t3.wilson import bc_mass_eigen_native_bridge_v21 as v21
from Numerical.parallel_t3.wilson import bc_mass_bridge_gate_v22 as v22

BASE=Path(__file__).resolve().parents[1]/'output/parallel_t3/wilson'
ROOT=Path(__file__).resolve().parents[1]


def mass(model='B_-1'):
    return MixedDoubletMassBasis.build(model,2.4e6,3.1e6,80000+42000j)


def synthetic_native():
    rng=np.random.default_rng(141)
    names=[f'lambdaT3ExtraQ{i}I1' for i in range(7)]
    extra=[]
    for i,name in enumerate(names):
        # Real integer CG fixture: avoids source-parser assumptions about floats.
        a=rng.integers(-2,3,size=(2,2,2,2))
        if not np.any(a):a[0,0,0,0]=1
        extra.append({'symbol':name,'fields':['S1dag','S1dag','S2','S2dag'],
                      'dimensions':[2]*4,
                      'assembled_lagrangian_contains_operator':True,
                      'components_inputform':[str(int(x)) for x in a.reshape(-1)]})
    topology_arr=rng.integers(-2,3,size=(2,2,2,2))
    topology_arr=(topology_arr+topology_arr.swapaxes(0,1))/2
    # Stick to integer values for restricted Mathematica source-parser.
    topology_arr=(2*topology_arr).astype(int)
    topology={'symbol':'lambdaT3','fields':['H','H','S1','Bar[S2]'],
              'dimensions':[2]*4,'assembled_lagrangian_contains_operator':True,
              'components_inputform':[str(int(x)) for x in topology_arr.reshape(-1)]}
    top={'models':{k:{'records':[topology]} for k in ('B_-1','C_-1')}}
    bc={'models':{k:{'records':extra} for k in ('B_-1','C_-1')}}
    B=np.eye(6,7)
    m={'BC_baseline_six_of_seven_conditional_source_maps':{}}
    for k in ('B_-1','C_-1'):
        m['BC_baseline_six_of_seven_conditional_source_maps'][k]={
            'matchete_seven_source_symbols':names,
            'baseline_rgbeta_six_symbols':list(v21.SPECIAL),
            'optional_extra_rgbeta_symbol':'lambdaT3GeneralQ110101I1',
            'baseline_map_rank':6,'optional_map_rank':7,
            'baseline_kernel_dimension_in_Matchete_seven_space':1,
            'source_monomial_not_native_projector_normalization':True,
            'baseline_rgbeta_from_Matchete_seven_complex_matrix':{'real':B.tolist(),'imag':(B*0).tolist()},
            'kernel_direction_normalized_complex':{'real':[0,0,0,0,0,0,1],'imag':[0]*7},
        }
    return top,bc,m

class TestV24HistoricalNullVectorSerialization(unittest.TestCase):
    """Exact producer shape of v19's NumPy SVD null vector."""

    def test_native_v19_svd_producer_serializes_flat(self):
        for model in ('B_-1','C_-1'):
            _,bc,witness=synthetic_native()
            d=witness['BC_baseline_six_of_seven_conditional_source_maps'][model]
            B=np.asarray(d['baseline_rgbeta_from_Matchete_seven_complex_matrix']['real'])
            _,_,vh=np.linalg.svd(B,full_matrices=True)
            null=vh[-1].conj();null/=np.linalg.norm(null)
            d['kernel_direction_normalized_complex']={
                'real':null.real.tolist(),'imag':null.imag.tolist()}
            source=[r['symbol'] for r in bc['models'][model]['records']]
            bn,n=v21.validate_special_matrix(witness,model,source)
            self.assertEqual(n.shape,(7,))
            self.assertLess(np.linalg.norm(bn@n),1e-12)

    def test_nontrivial_bc_source_matrix_flat_svd_null_after_native_reorder(self):
        # The historical B/C 110101 two-dimensional sector is not diagonal.
        # Its optional seventh direction must survive the named permutation.
        historical=[
            'lambdaT3ExtraQ110101I1', 'lambdaT3ExtraQ110101I2',
            'lambdaT3ExtraQ000112I1', 'lambdaT3ExtraQ000202I1',
            'lambdaT3ExtraQ001201I1', 'lambdaT3ExtraQ020020I1',
            'lambdaT3ExtraQ020200I1',
        ]
        g=np.sqrt(3.)
        A=np.zeros((7,7),dtype=complex)
        A[0,1]=-2/g
        A[1,2]=A[2,3]=A[3,4]=A[4,5]=A[5,6]=1
        A[6,0]=1;A[6,1]=1/g
        rgb=[
            'lambdaT3GeneralQ110101I1','lambdaS1barS2S2bar2',
            'lambdaS1bar2S2bar2','lambdaS1S1bar2S2bar',
            'lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
            'lambdaHHdagS1barS2barCross',
        ]
        B=A[[rgb.index(name) for name in v21.SPECIAL]]
        _,_,vh=np.linalg.svd(B,full_matrices=True)
        null=vh[-1].conj()
        null/=np.linalg.norm(null)
        _,bc,witness=synthetic_native()
        for model in ('B_-1','C_-1'):
            d=witness['BC_baseline_six_of_seven_conditional_source_maps'][model]
            d['matchete_seven_source_symbols']=historical
            d['baseline_rgbeta_from_Matchete_seven_complex_matrix']={
                'real':B.real.tolist(),'imag':B.imag.tolist()}
            d['kernel_direction_normalized_complex']={
                'real':null.real.tolist(),'imag':null.imag.tolist()}
            # v18's verified-structure order: the 110101 pair is last.
            native=historical[2:]+historical[:2]
            gotB,gotN=v21.validate_special_matrix(witness,model,native)
            np.testing.assert_allclose(gotB,B[:,[2,3,4,5,6,0,1]])
            np.testing.assert_allclose(gotN,null[[2,3,4,5,6,0,1]])
            self.assertLess(np.linalg.norm(gotB@gotN),1e-12)
            self.assertEqual(np.linalg.matrix_rank(gotB),6)

    def test_previous_test_only_single_row_remains_accepted(self):
        _,bc,witness=synthetic_native()
        model='B_-1'
        d=witness['BC_baseline_six_of_seven_conditional_source_maps'][model]
        d['kernel_direction_normalized_complex']={'real':[[0,0,0,0,0,0,1]],'imag':[[0]*7]}
        source=[r['symbol'] for r in bc['models'][model]['records']]
        b,n=v21.validate_special_matrix(witness,model,source)
        self.assertEqual(n.shape,(7,))
        self.assertLess(np.linalg.norm(b@n),1e-12)

    def test_malformed_vector_shapes_nonfinite_and_broadcast_rejected(self):
        _,bc,witness=synthetic_native();model='C_-1'
        source=[r['symbol'] for r in bc['models'][model]['records']]
        wrong=[
            {'real':[0,0,0,0,0,0,1],'imag':[[0]*7]},
            {'real':[[0],[0],[0],[0],[0],[0],[1]],'imag':[[0],[0],[0],[0],[0],[0],[0]]},
            {'real':[[[0,0,0,0,0,0,1]]],'imag':[[[0]*7]]},
            {'real':[0]*6,'imag':[0]*6},
            {'real':[0,0,0,0,0,0,float('nan')],'imag':[0]*7},
            {'real':[0,0,0,0,0,0,float('inf')],'imag':[0]*7},
            {'real':[0,0,0,0,0,0,'bad'],'imag':[0]*7},
        ]
        for bad in wrong:
            with self.subTest(bad=repr(bad)[:50]):
                obj=copy.deepcopy(witness)
                obj['BC_baseline_six_of_seven_conditional_source_maps'][model][
                    'kernel_direction_normalized_complex']=bad
                with self.assertRaises(v21.BCBridgeBlocked):
                    v21.validate_special_matrix(obj,model,source)

    def test_flat_null_tampering_detected(self):
        _,bc,witness=synthetic_native();model='B_-1'
        source=[r['symbol'] for r in bc['models'][model]['records']]
        witness['BC_baseline_six_of_seven_conditional_source_maps'][model][
            'kernel_direction_normalized_complex']['real']=[1,0,0,0,0,0,0]
        with self.assertRaisesRegex(v21.BCBridgeBlocked,'rank defective'):
            v21.validate_special_matrix(witness,model,source)

    def test_matrix_real_imag_shape_mismatch_rejected(self):
        _,bc,witness=synthetic_native()
        d=witness['BC_baseline_six_of_seven_conditional_source_maps']['B_-1']
        d['baseline_rgbeta_from_Matchete_seven_complex_matrix']['imag']=[0]*7
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.validate_special_matrix(witness,'B_-1',
                [r['symbol'] for r in bc['models']['B_-1']['records']])


class TestV21Math(unittest.TestCase):
    def test_complex_mass_eigenvalues_and_quadratic_transport_both(self):
        for model in ('B_-1','C_-1'):
            basis=mass(model)
            info=v21.gauge_covariance_recheck(basis)
            self.assertLess(info['su2_covariance_max_absolute'],1e-12)
            self.assertEqual(len(info['quadratic_masses_gev']),2)
            self.assertLess(info['quadratic_tensor_transport_max_absolute_gev2'],1e-7)

    def test_epsilon_pseudoreality_and_inverse(self):
        basis=mass()
        h=np.array([1+2j,-3+0.2j])
        x=np.array([.2+.8j,2-.9j]);y=np.array([-1+.9j,3j])
        s1,s2=basis.original_fields(x,y)
        xx,yy=basis.eigen_fields(s1,s2)
        np.testing.assert_allclose(xx,x,atol=1e-12)
        np.testing.assert_allclose(yy,y,atol=1e-12)
        self.assertFalse(np.allclose(s2,basis.U[1,0]*x+basis.U[1,1]*y))

    def test_all_four_doublet_species_rotated_with_H_holomorphic(self):
        basis=mass();rng=np.random.default_rng(114)
        for names in [('H','H','S1','S2dag'),('S1dag','S1dag','S2','S2dag'),
                      ('S1','S1dag','S2','S2dag'),('Hdag','Hdag','S2','S2'),
                      ('H','Hdag','S1dag','S2dag')]:
            t=rng.normal(size=(2,)*4)+1j*rng.normal(size=(2,)*4)
            rot=v21.rotate_tensor(t,names,basis.U)
            self.assertLess(v21.action_residual(t,names,rot,basis,trials=7),1e-12)

    def test_wrong_sign_in_epsilon_fails(self):
        basis=mass();rng=np.random.default_rng(410)
        t=rng.normal(size=(2,)*4)+1j*rng.normal(size=(2,)*4)
        fields=('H','H','S1','S2dag')
        correct=v21.rotate_tensor(t,fields,basis.U)
        wrong={k:(-z if k[-1] in ('X','Y') else z) for k,z in correct.items()}
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.action_residual(t,fields,wrong,basis,seed=144,trials=4)

    def test_tensor_shape_and_unknown_species_block(self):
        basis=mass()
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.rotate_tensor(np.zeros((2,2,2,3)),('H','H','S1','S2dag'),basis.U)
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.convert_fields(['H','H','S1','F'])

    def test_topology_native_source_record_parses(self):
        top,bc,m=synthetic_native()
        sym,fields,t=v21.tensor_record(top['models']['B_-1']['records'][0])
        self.assertEqual(sym,'lambdaT3')
        self.assertEqual(fields,('H','H','S1','S2dag'))
        np.testing.assert_allclose(t,t.swapaxes(0,1))

    def test_complete_synthetic_8_complex_native_CGs_each(self):
        top,bc,v19=synthetic_native()
        for model in ('B_-1','C_-1'):
            x=v21.audit_model(model,mass(model),top,bc,v19,seed=777,trials=8)
            self.assertEqual(x['source_native_complex_CG_count'],8)
            self.assertEqual(x['baseline_six_from_seven_rank'],6)
            self.assertTrue(x['optional_seventh_source_direction_action_nonzero'])
            self.assertLess(x['max_source_native_action_relative_residual'],2e-11)

    def test_source_seventh_subspace_rank_guard(self):
        top,bc,v19=synthetic_native()
        src=[x['symbol'] for x in bc['models']['B_-1']['records']]
        b,n=v21.validate_special_matrix(v19,'B_-1',src)
        self.assertEqual(np.linalg.matrix_rank(b),6)
        self.assertLess(np.linalg.norm(b@n),1e-12)
        v19['BC_baseline_six_of_seven_conditional_source_maps']['B_-1']['baseline_map_rank']=7
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.validate_special_matrix(v19,'B_-1',src)

    def test_missing_complex_conjugate_rejected(self):
        top,bc,v19=synthetic_native()
        orig=v21.sample_vectors
        def fake(basis,rng):
            old,new=orig(basis,rng)
            # Wrong S2 complex-field conjugation in direct action.
            old=dict(old);old['S2dag']=old['S2']
            return old,new
        with mock.patch.object(v21,'sample_vectors',side_effect=fake):
            with self.assertRaises(v21.BCBridgeBlocked):
                v21.audit_model('B_-1',mass(),top,bc,v19,seed=333,trials=6)

    def test_false_matchete_native_component_fields_rejected(self):
        top,bc,v19=synthetic_native()
        rec=copy.deepcopy(top['models']['B_-1']['records'][0]);rec['assembled_lagrangian_contains_operator']=False
        with self.assertRaises(v21.BCBridgeBlocked):v21.tensor_record(rec)

    def test_special_map_wrong_null_direction_rejected(self):
        top,bc,v19=synthetic_native()
        x=v19['BC_baseline_six_of_seven_conditional_source_maps']['B_-1']
        x['kernel_direction_normalized_complex']['real']=[1,0,0,0,0,0,0]
        with self.assertRaises(v21.BCBridgeBlocked):
            v21.audit_model('B_-1',mass(),top,bc,v19)

    def test_saved_mass_basis_stale_hash_fails_closed(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'source.txt').write_text('new',encoding='utf-8')
            witness=root/'m.json'
            witness.write_text(json.dumps({
                'format':v21.MASS_FORMAT,'status':v21.MASS_STATUS,
                'model':'B_-1','endpoint':'final',
                'physical_mass_matching_performed':False,
                'MatcheteRGBetaPhysicalCouplingMapVerified':False,
                'MatcheteHeavyMassBasisCheckOnRedefinedLagrangian':False,
                'quadratic_hessian_agreement':True,'source_quartic_count':17,
                'max_source_tensor_transport_relative_error':0.,
                'source_sha256':{'source.txt':'0'*64}}),encoding='utf-8')
            with self.assertRaises(v21.BCBridgeBlocked):
                v21.read_mass_basis(root,'B_-1',witness)

    def test_source_json_nan_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'x.json';p.write_text('{"x":NaN}',encoding='utf-8')
            with self.assertRaises(v21.BCBridgeBlocked):v21.read_json(p)


class TestV22Gate(unittest.TestCase):
    def test_gate_source_guard_physical_false(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            v20p=root/'v20.json';v21p=root/'v21.json'
            v20p.write_text(json.dumps({'format':v22.V20_FORMAT,'status':v22.V20_STATUS,
                'total_EFT1_scalar_quartic_registration_slots':178,
                'native_RGBeta_real_quartics':150,'native_RGBeta_complex_quartics':28,
                'certification':{'physical_Wilson_evolution_authorized':False}}))
            v21p.write_text(json.dumps({'format':v22.V21_FORMAT,'status':v22.V21_STATUS,
                'models_checked':2,'native_Matchete_complex_CG_sources_rotated':16,
                'native_topology_CGs_rotated':2,
                'native_BC_optional_extra_CGs_rotated':14,
                'baseline_RGBeta_special_complex_quartic_slots_conditionally_mapped':12,
                'max_complex_quartic_action_transport_relative_residual':0,
                'certification':{'physical_Wilson_evolution_authorized':True}}))
            with self.assertRaises(v21.BCBridgeBlocked):v22.audit(root,v20p,v21p)

    def test_source_inventory_invalid_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);p=root/'a.json';q=root/'b.json'
            p.write_text(json.dumps({'format':v22.V20_FORMAT,'status':v22.V20_STATUS,
                                     'total_EFT1_scalar_quartic_registration_slots':177}))
            q.write_text('{}')
            with self.assertRaises(v21.BCBridgeBlocked):v22.audit(root,p,q)

    def test_all16_physical_mode_remains_off(self):
        self.assertIn('physical_one_loop_dim5_ADM',v22.BLOCKERS)
        self.assertIn('B_C_mixed_scalar_HardMatching',v22.BLOCKERS)
        self.assertEqual(len(v22.BLOCKERS),8)

    @unittest.skipUnless((BASE/'full_scalar_census_conditional_v20.json').exists()
                         and (ROOT/'output/general_t3_mass_basis/B_m1_final_mass_basis.json').exists()
                         and (ROOT/'output/general_t3_mass_basis/C_m1_final_mass_basis.json').exists(),
                         'requires full current local repository and v20/native B/C evidence')
    def test_current_local_repository_sources_end_to_end(self):
        result=v21.audit(ROOT,
           topology_path=BASE/'matchete_topology_native_v18.json',
           bc_path=BASE/'matchete_bc_extra_native_v18.json',
           v19_path=BASE/'complex_scalar_bridge_conditional_v19.json',
           v20_path=BASE/'full_scalar_census_conditional_v20.json',
           mass_paths={tag:ROOT/f'output/general_t3_mass_basis/{tag[0]}_m1_final_mass_basis.json'
                       for tag in ('B_-1','C_-1')},trials=6)
        self.assertLess(result['max_complex_quartic_action_transport_relative_residual'],2e-11)
        self.assertFalse(result['certification']['physical_Wilson_evolution_authorized'])


class TestV21WindowsProvenanceRegression(unittest.TestCase):
    """Reproduces saved v19 on Windows: str(relative path) uses backslashes."""

    def fixture(self, root: Path, *, windows: bool = True):
        from hashlib import sha256
        base = Path('output/parallel_t3/wilson')
        names = (
            'matchete_topology_native_v18.json',
            'matchete_bc_extra_native_v18.json',
            'rgbeta_quartic_projectors_native_v14.json',
            'rgbeta_complex_conjugate_native_v15.json',
            'matchete_native_tensors.json',
            'rgbeta_abstract_components_conditional_v2.json',
            'additional_quartic_dictionary_conditional_v3.json',
        )
        stamps = {}
        for i, name in enumerate(names):
            relative = (base / name).as_posix()
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(f'fixture {i} {name}'.encode('utf-8'))
            key = relative.replace('/', chr(92)) if windows else relative
            stamps[key] = sha256(path.read_bytes()).hexdigest()
        return {'input_sha256': stamps}, tuple(root / base / name for name in names)

    def test_windows_v19_source_keys_verify_all_seven_unmodified_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root, windows=True)
            result = v21.verify_v19_input_sha256(root, witness, paths[:2])
            self.assertEqual(result['source_files_sha256_rechecked'], 7)
            self.assertEqual(result['legacy_windows_separator_keys_accepted'], 7)
            self.assertTrue(result['provenance_hashes_unchanged'])
            # The saved v19 file itself is NOT rewritten to POSIX paths.
            self.assertTrue(all(chr(92) in k for k in witness['input_sha256']))

    def test_posix_v19_source_keys_still_supported(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root, windows=False)
            result = v21.verify_v19_input_sha256(root, witness, paths[:2])
            self.assertEqual(result['source_files_sha256_rechecked'], 7)
            self.assertEqual(result['legacy_windows_separator_keys_accepted'], 0)

    def test_genuine_changed_source_rejected_not_excused_by_separator(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root)
            # This is a real byte change, not merely a spelling difference.
            paths[0].write_bytes(paths[0].read_bytes() + b' changed')
            with self.assertRaisesRegex(v21.BCBridgeBlocked, 'SHA256 actually differs'):
                v21.verify_v19_input_sha256(root, witness, paths[:2])

    def test_missing_source_or_digest_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root)
            paths[-1].unlink()
            with self.assertRaisesRegex(v21.BCBridgeBlocked, 'input missing'):
                v21.verify_v19_input_sha256(root, witness, paths[:2])
            witness, paths = self.fixture(root)
            witness['input_sha256'].pop(next(iter(witness['input_sha256'])))
            with self.assertRaisesRegex(v21.BCBridgeBlocked, 'exactly seven'):
                v21.verify_v19_input_sha256(root, witness, paths[:2])

    def test_ambiguous_slash_alias_traversal_and_absolute_path_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root)
            example = next(iter(witness['input_sha256']))
            cases = [
                r'..\output\parallel_t3\wilson\evil.json',
                r'C:\output\parallel_t3\wilson\evil.json',
            ]
            # Both Windows and POSIX spellings of the same source are not
            # allowed to create two contradictory provenance entries.
            second = list(witness['input_sha256'])[1]
            del witness['input_sha256'][second]
            witness['input_sha256'][example.replace(chr(92), '/')] = witness['input_sha256'][example]
            with self.assertRaisesRegex(v21.BCBridgeBlocked, 'aliases'):
                v21.verify_v19_input_sha256(root, witness, paths[:2])
            for bad in cases:
                doc, paths = self.fixture(root)
                del doc['input_sha256'][next(iter(doc['input_sha256']))]
                doc['input_sha256'][bad] = '0' * 64
                with self.assertRaises(v21.BCBridgeBlocked):
                    v21.verify_v19_input_sha256(root, doc, paths[:2])

    def test_readonly_preflight_reuses_windows_sources_without_rewriting(self):
        from scripts import RunT3ParallelABCMassProvenanceV22 as pre
        from hashlib import sha256
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root)
            v19_path = root/'output/parallel_t3/wilson/complex_scalar_bridge_conditional_v19.json'
            witness.update({'format':v21.V19_FORMAT, 'status':v21.V19_STATUS,
                            'baseline_EFT1_complex_slots':28,
                            'certification':{'physical_Wilson_evolution_authorized':False}})
            v19_path.write_text(json.dumps(witness),encoding='utf-8')
            v20_path = root/'output/parallel_t3/wilson/full_scalar_census_conditional_v20.json'
            v20_path.write_text(json.dumps({
                'format':v21.V20_FORMAT,'status':v21.V20_STATUS,
                'total_EFT1_scalar_quartic_registration_slots':178,
                'input_sha256':{'complex_bridge_v19':sha256(v19_path.read_bytes()).hexdigest()},
                'certification':{'physical_Wilson_evolution_authorized':False},
            }),encoding='utf-8')
            sources_before = {p:sha256(p.read_bytes()).hexdigest() for p in paths}
            def fake_read(repo,tag,p):
                return mass(tag), {'source_sha256':{'existing_source':'0'*64}}
            with mock.patch.object(pre,'read_mass_basis',side_effect=fake_read):
                report = pre.preflight(root)
            self.assertEqual(report['v19']['source_files_sha256_rechecked'],7)
            self.assertEqual(report['v19']['legacy_windows_separator_keys_accepted'],7)
            self.assertTrue(report['v20_v19_sha256_chain_verified'])
            self.assertFalse(report['physical_matching_authorized'])
            self.assertEqual(sources_before,{p:sha256(p.read_bytes()).hexdigest() for p in paths})

    def test_full_v21_v22_synthetic_bridge_with_windows_v19_hashes(self):
        """Exercise the real audit entrypoints beyond the originally failing line."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            witness, paths = self.fixture(root, windows=True)
            top, bc, map19 = synthetic_native()
            # Simulate the actual v18/v19 discrepancy at the complete
            # source-hash-checked audit entrypoint, not just audit_model.
            native_order=[2,3,4,5,6,0,1]
            for tag in ('B_-1','C_-1'):
                records=bc['models'][tag]['records']
                bc['models'][tag]['records']=[records[i] for i in native_order]
            top_path, bc_path = paths[:2]
            top_path.write_text(json.dumps(top), encoding='utf-8')
            bc_path.write_text(json.dumps(bc), encoding='utf-8')
            from hashlib import sha256
            for p in (top_path, bc_path):
                witness['input_sha256'][p.relative_to(root).as_posix().replace('/',chr(92))] = sha256(p.read_bytes()).hexdigest()
            map19.update(witness)
            map19.update({'format':v21.V19_FORMAT, 'status':v21.V19_STATUS,
                          'baseline_EFT1_complex_slots':28,
                          'certification':{'physical_Wilson_evolution_authorized':False}})
            bridge = root/'output/parallel_t3/wilson/complex_scalar_bridge_conditional_v19.json'
            bridge.write_text(json.dumps(map19),encoding='utf-8')
            v20_path = root/'output/parallel_t3/wilson/full_scalar_census_conditional_v20.json'
            doc20 = {'format':v21.V20_FORMAT, 'status':v21.V20_STATUS,
                     'total_EFT1_scalar_quartic_registration_slots':178,
                     'native_RGBeta_real_quartics':150,
                     'native_RGBeta_complex_quartics':28,
                     'certification':{'physical_Wilson_evolution_authorized':False,
                                      'BC_minus1_mixed_scalar_mass_basis_certified':False},
                     'input_sha256':{'complex_bridge_v19':sha256(bridge.read_bytes()).hexdigest()}}
            v20_path.write_text(json.dumps(doc20),encoding='utf-8')
            masspaths = {}
            for model in ('B_-1','C_-1'):
                p = root/f'output/general_t3_mass_basis/{model[0]}_m1_final_mass_basis.json'
                p.parent.mkdir(parents=True,exist_ok=True)
                p.write_text('{}',encoding='utf-8')
                masspaths[model] = p
            def saved_mass(repo, model, file):
                return mass(model), {'synthetic_fixture': True}
            with mock.patch.object(v21,'native_document',return_value={}), \
                 mock.patch.object(v21,'read_mass_basis',side_effect=saved_mass):
                result=v21.audit(root,topology_path=top_path,bc_path=bc_path,
                                 mass_paths=masspaths,v19_path=bridge,
                                 v20_path=v20_path,trials=6)
            self.assertEqual(result['v19_input_provenance']['source_files_sha256_rechecked'],7)
            self.assertEqual(result['v19_input_provenance']['legacy_windows_separator_keys_accepted'],7)
            self.assertEqual(result['native_Matchete_complex_CG_sources_rotated'],16)
            for model in ('B_-1','C_-1'):
                self.assertEqual(result['models'][model]['v19_column_indices_in_native_v18_order'],native_order)
            self.assertLess(result['max_complex_quartic_action_transport_relative_residual'],2e-11)
            v21_path = root/'output/parallel_t3/wilson/bc_native_mass_eigen_complex_conditional_v21.json'
            v21_path.write_text(json.dumps(result),encoding='utf-8')
            gate = v22.audit(root,v20_path,v21_path)
            self.assertEqual(gate['full_scalar_quartic_slots'],178)
            self.assertFalse(gate['certification']['physical_Wilson_evolution_authorized'])


if __name__=='__main__':unittest.main()


class TestHistoricalV19NativeV18SourceOrderRegression(unittest.TestCase):
    """v19 saved historical Stage22 order; v18 enumerated structure-sorted CGs.

    The physics map is a linear map on NAMED source couplings, not a map on
    list positions. Test both matrix and its null direction after reordering.
    """

    @staticmethod
    def stage22_native_reorder():
        top, bc, old = synthetic_native()
        # Match the observed historical/v18 difference: older Stage22 starts
        # with 110101I1,I2, while the v18 verified-structure order puts them last.
        permutation = [2, 3, 4, 5, 6, 0, 1]
        for model in ('B_-1', 'C_-1'):
            recs = bc['models'][model]['records']
            bc['models'][model]['records'] = [recs[j] for j in permutation]
        return top, bc, old, permutation

    def test_exact_named_column_permutation_of_matrix_and_kernel(self):
        top, bc, witness, permutation = self.stage22_native_reorder()
        model = 'B_-1'
        native = [x['symbol'] for x in bc['models'][model]['records']]
        old_B = np.asarray(witness['BC_baseline_six_of_seven_conditional_source_maps'][model][
            'baseline_rgbeta_from_Matchete_seven_complex_matrix']['real'], float)
        old_N = np.asarray(witness['BC_baseline_six_of_seven_conditional_source_maps'][model][
            'kernel_direction_normalized_complex']['real'],float).reshape(7)
        B, N = v21.validate_special_matrix(witness, model, native)
        np.testing.assert_array_equal(B, old_B[:, permutation])
        np.testing.assert_array_equal(N, old_N[permutation])
        old_q = np.arange(1.,8.)+0.07j*np.arange(7)
        native_q = old_q[permutation]
        np.testing.assert_allclose(B@native_q,old_B@old_q)
        np.testing.assert_allclose(np.linalg.norm(B@N),0.,atol=1e-13)

    def test_complete_both_models_with_reordered_v18_records(self):
        top, bc, witness, permutation = self.stage22_native_reorder()
        for model in ('B_-1', 'C_-1'):
            out = v21.audit_model(model, mass(model), top, bc, witness, trials=7)
            self.assertEqual(out['v19_column_indices_in_native_v18_order'],permutation)
            self.assertTrue(out['source_basis_reordered_by_exact_coupling_symbol'])
            self.assertTrue(out['optional_seventh_source_direction_action_nonzero'])
            self.assertEqual(out['baseline_six_from_seven_rank'],6)
            self.assertLess(out['max_source_native_action_relative_residual'],2e-11)
            self.assertLess(out['optional_seventh_source_direction_worst_holdout_residual'],2e-11)

    def test_unknown_duplicate_or_missing_coupling_not_permitted(self):
        top, bc, witness, _ = self.stage22_native_reorder()
        model='B_-1'
        native=[r['symbol'] for r in bc['models'][model]['records']]
        with self.assertRaisesRegex(v21.BCBridgeBlocked,'inventories differ'):
            v21.validate_special_matrix(witness,model,native[:-1]+['lambdaT3ExtraQ_FAKEI1'])
        with self.assertRaisesRegex(v21.BCBridgeBlocked,'inventories differ'):
            v21.validate_special_matrix(witness,model,native[:-1]+[native[0]])
        with self.assertRaisesRegex(v21.BCBridgeBlocked,'inventories differ'):
            v21.validate_special_matrix(witness,model,native[:-1])

    def test_reordering_does_not_mask_tampered_native_matrix(self):
        top,bc,witness,_=self.stage22_native_reorder()
        model='B_-1'
        native=[r['symbol'] for r in bc['models'][model]['records']]
        d=witness['BC_baseline_six_of_seven_conditional_source_maps'][model]
        d['kernel_direction_normalized_complex']['real']=[1,0,0,0,0,0,0]
        with self.assertRaisesRegex(v21.BCBridgeBlocked,'rank defective'):
            v21.validate_special_matrix(witness,model,native)

    def test_actual_b_minus_one_historical_coupling_names(self):
        _,bc,witness,_=self.stage22_native_reorder()
        model='B_-1'
        old = witness['BC_baseline_six_of_seven_conditional_source_maps'][model]
        actual_old = [
            'lambdaT3ExtraQ110101I1', 'lambdaT3ExtraQ110101I2',
            'lambdaT3ExtraQ000112I1', 'lambdaT3ExtraQ000202I1',
            'lambdaT3ExtraQ001201I1', 'lambdaT3ExtraQ020020I1',
            'lambdaT3ExtraQ020200I1',
        ]
        # Substitute real historical symbol inventory while retaining test CGs.
        test_old=old['matchete_seven_source_symbols']
        remap=dict(zip(test_old,actual_old))
        old['matchete_seven_source_symbols']=actual_old
        for rec in bc['models'][model]['records']:
            rec['symbol']=remap[rec['symbol']]
        new=[r['symbol'] for r in bc['models'][model]['records']]
        _, indices=v21.special_source_permutation(witness,model,new)
        self.assertEqual(indices,[2,3,4,5,6,0,1])



class TestV23NativeStage22SourceAxisPreflight(unittest.TestCase):
    def test_source_preflight_validates_reordered_basis_and_hash_chain(self):
        from scripts import RunT3ParallelABCMassSubspaceOrderV23 as pre
        from hashlib import sha256
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            hashes, files=TestV21WindowsProvenanceRegression().fixture(root)
            top,bc,mapping,_=TestHistoricalV19NativeV18SourceOrderRegression.stage22_native_reorder()
            top_path,bc_path=files[:2]
            top_path.write_text(json.dumps(top),encoding='utf-8')
            bc_path.write_text(json.dumps(bc),encoding='utf-8')
            for f in (top_path,bc_path):
                key=f.relative_to(root).as_posix().replace('/',chr(92))
                hashes['input_sha256'][key]=sha256(f.read_bytes()).hexdigest()
            mapping.update(hashes)
            mapping.update({'format':v21.V19_FORMAT,'status':v21.V19_STATUS,
                            'baseline_EFT1_complex_slots':28,
                            'certification':{'physical_Wilson_evolution_authorized':False}})
            v19_path=root/'output/parallel_t3/wilson/complex_scalar_bridge_conditional_v19.json'
            v19_path.write_text(json.dumps(mapping),encoding='utf-8')
            v20_path=root/'output/parallel_t3/wilson/full_scalar_census_conditional_v20.json'
            v20_path.write_text(json.dumps({'format':v21.V20_FORMAT,'status':v21.V20_STATUS,
                'total_EFT1_scalar_quartic_registration_slots':178,
                'input_sha256':{'complex_bridge_v19':sha256(v19_path.read_bytes()).hexdigest()},
                'certification':{'physical_Wilson_evolution_authorized':False}}),encoding='utf-8')
            with mock.patch.object(pre,'native_document',return_value={}):
                report=pre.preflight(root)
            self.assertFalse(report['physical_matching_authorized'])
            for row in report['models'].values():
                self.assertEqual(row['v19_column_indices_in_native_order_zero_based'],
                                 [2,3,4,5,6,0,1])
                self.assertLess(row['kernel_residual'],1e-10)
            # Really stale v20 hash still blocks; reordering cannot excuse it.
            v19_path.write_text(v19_path.read_text()+" ",encoding='utf-8')
            with mock.patch.object(pre,'native_document',return_value={}):
                with self.assertRaises(v21.BCBridgeBlocked):pre.preflight(root)

    def test_runner_performs_source_order_check_and_preserves_existing_data(self):
        source=(ROOT/'scripts/RunT3ParallelABCMassEigenBatchV21V22.ps1').read_text()
        self.assertIn('RunT3ParallelABCMassSubspaceOrderV23',source)
        self.assertIn('RunT3ParallelABCMassProvenanceV22',source)
        self.assertIn('RunT3ParallelABCMassEigenNativeV21',source)
        self.assertNotIn('Remove-Item',source)
        self.assertNotIn('Set-Content',source)
