"""All-model offline validation: exercises the release parser using synthetic data.

These are NOT numerical re-runs of the user's saved Matchete/CG coefficients.
Actual runs are performed by scripts/RunLOSpinorCertification.py on local input.
"""
import copy
from pathlib import Path
import tempfile
import unittest
from Numerical.lo_certification.spinor.boundary import (
    MODELS, EXPECTED, BoundaryError, SOURCE_FORMAT, canonical, certify_source_release,
    extract_rows, native_evidence_probe, orbit, tag, validate_model, weave_source_handoff,
)


def row(key=(1, 1, 1, 1, 5, 5), real=2., imag=3.):
    return {'coordinate': list(key), 'coefficient_real': real,
            'coefficient_imag': imag, 'orbit_size': len(orbit(key))}


def mock_source(model='A_-4', *, split=False, diagonal=False, phases=False):
    dims = ({'NewScalar1': 2, 'NewScalar2': 2} if not split else
            {'NewScalar1': 1, 'T3Scalar2Minus': 3, 'T3Scalar2Plus': 3})
    offsets = ({'NewScalar1': 5, 'NewScalar2': 9} if not split else
               {'NewScalar1': 5, 'T3Scalar2Minus': 7, 'T3Scalar2Plus': 10})
    key = (1, 1, 2, 3, 5, 5) if diagonal else (1, 1, 2, 3, 5, 9)
    values = [1+2j, 3-1j] if phases else [2+3j]
    return {
        'format': SOURCE_FORMAT, 'model': model,
        'status': 'NumericalLLSSBySourceVerified_SpinorNormalizationOpen',
        'mu_F_gev': 100000, 'units': 'GeV^-1',
        'field_content': {'heavy_fermions_integrated_out': 'all at mu_F',
                          'heavy_scalar_multiplets': 'both retained in EFT1'},
        'index_order': ['lepton_su2_i', 'lepton_flavour_p', 'lepton_su2_j',
                        'lepton_flavour_q', 'real_scalar_a', 'real_scalar_b'],
        'index_base': 1, 'real_scalar_dimension': 12,
        'scalar_field_complex_dimensions': dims,
        'scalar_field_offsets_one_based': offsets,
        'completeness': {
            'F_tree_LLSS_included': True,
            'distinct_PL_sources_retained': True,
            'PR_Hermitian_conjugate_not_double_counted': True,
            'one_loop_F_hard_matching_included': False,
            'complete_physical_leading_one_loop_neutrino_mass_certified': False,
        },
        'spinor_normalization': {'certified': False,
            'gamma_cc_to_master_weyl_factor': None,
            'gamma_cc_to_master_weyl_sign': None},
        'formal_polynomial': 'sum C[i,p,j,q,a,b] psi[i,p] psi[j,q] phi[a] phi[b]',
        'sources': [{
            'matchete_source_index': i+1,
            'native_chirality': 'PL',
            'flavour_permutation_lepton_order': [1, 0],
            'cg_registry_used': ['RealNativeCG'],
            'inputform_operator': 'Coupling[y1, {}]*NCM[Transp[Field[l, Fermion, {}, {}]], DiracProduct[GammaCC, Proj[-1]], Field[l, Fermion, {}, {}]]',
            'zero_at_this_benchmark': False,
            'independent_symmetric_coordinates': [row(key, x.real, x.imag)]
        } for i, x in enumerate(values)],
        'summed_independent_symmetric_coordinates': [row(key, sum(x.real for x in values),
                                                            sum(x.imag for x in values))],
    }


class TensorTests(unittest.TestCase):
    def test_orbit_sizes_and_fermion_scalar_exchange(self):
        self.assertEqual(len(orbit((1, 1, 1, 1, 5, 5))), 1)
        self.assertEqual(len(orbit((1, 1, 1, 1, 5, 9))), 2)
        self.assertEqual(len(orbit((1, 1, 2, 2, 5, 5))), 2)
        self.assertEqual(len(orbit((1, 1, 2, 2, 5, 9))), 4)
        self.assertEqual(canonical((2, 2, 1, 1, 9, 5)), (1, 1, 2, 2, 5, 9))

    def test_all_16_model_parsers_synthetic(self):
        self.assertEqual(len(MODELS), 16)
        self.assertEqual(set(MODELS), set(EXPECTED))
        for model in MODELS:
            with self.subTest(model=model):
                data = mock_source(model, split=('A_-2' == model))
                result = validate_model(data, model, strict_counts=False)
                self.assertEqual(result['verified_source_count_PL'], 1)
                self.assertEqual(result['summed_coordinate_rows'], 1)
                self.assertEqual(result['max_relative_source_sum_residual'], 0)
                h = weave_source_handoff(data, result, 'a'*64)
                self.assertIsNone(h['physical_master_weyl_boundary'])
                self.assertFalse(h['downstream_usage']['safe_for_physical_EFT1_or_scalar_matching'])

    def test_complex_source_sum_cancellation_and_phases(self):
        d = mock_source('B_-1', phases=True)
        r = validate_model(d, 'B_-1', strict_counts=False)
        self.assertEqual(r['verified_source_count_PL'], 2)
        self.assertEqual(r['summed_coordinate_rows'], 1)
        self.assertEqual(r['source_terms'][1]['independent_coordinates'][0]['coefficient_imag'], -1)

    def test_native_PR_never_double_counted(self):
        d = mock_source()
        d['sources'][0]['inputform_operator'] = d['sources'][0]['inputform_operator'].replace('Proj[-1]', 'Proj[1]')
        with self.assertRaises(BoundaryError):
            validate_model(d, 'A_-4', strict_counts=False)

    def test_bad_orbit_rejected(self):
        d = mock_source()
        d['sources'][0]['independent_symmetric_coordinates'][0]['orbit_size'] += 1
        with self.assertRaises(BoundaryError):
            validate_model(d, 'A_-4', strict_counts=False)

    def test_double_counting_rejected(self):
        d = mock_source()
        d['summed_independent_symmetric_coordinates'][0]['coefficient_real'] *= 2
        with self.assertRaises(BoundaryError):
            validate_model(d, 'A_-4', strict_counts=False)

    def test_requires_two_original_multiplets(self):
        d = mock_source()
        d['field_content']['heavy_scalar_multiplets'] = 'one retained in EFT1'
        with self.assertRaises(BoundaryError):
            validate_model(d, 'A_-4', strict_counts=False)

    def test_cannot_silently_set_phase(self):
        d = mock_source()
        d['spinor_normalization']['gamma_cc_to_master_weyl_factor'] = 1
        with self.assertRaises(BoundaryError):
            validate_model(d, 'A_-4', strict_counts=False)

    def test_split_plus_minus_real_basis(self):
        r = validate_model(mock_source('A_-2', split=True), 'A_-2', strict_counts=False)
        self.assertIn('T3Scalar2Minus', r['scalar_multiplets_active_in_EFT1'])
        self.assertIn('T3Scalar2Plus', r['scalar_multiplets_active_in_EFT1'])

    def test_zero_source_retained(self):
        d = mock_source('C_-1')
        empty = copy.deepcopy(d['sources'][0]); empty['matchete_source_index'] = 19
        empty['zero_at_this_benchmark'] = True; empty['independent_symmetric_coordinates'] = []
        d['sources'].append(empty)
        result = validate_model(d, 'C_-1', strict_counts=False)
        self.assertEqual(result['verified_source_count_PL'], 2)
        self.assertEqual(result['source_coordinate_rows'], 1)
        self.assertTrue(result['source_terms'][1]['zero_at_benchmark'])

    def test_all_16_batch_handoff_synthetic_and_idempotent(self):
        import json
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            release = root/'output/general_observable/lo_fermion_boundary/releases/source_fake'
            release.mkdir(parents=True)
            inventory = []
            for model in MODELS:
                sample = mock_source(model)
                path = release/tag(model)/'boundary.json'
                path.parent.mkdir(parents=True)
                path.write_text(json.dumps(sample), encoding='utf-8')
                inventory.append({
                    'model': model, 'path': f'{tag(model)}/boundary.json',
                    'PL_sources': 1, 'PR_conjugates_excluded': EXPECTED[model][1],
                    'source_resolved_nonzero_coordinates': 1,
                    'summed_nonzero_coordinates': 1,
                })
            summary = {
                'format': SOURCE_FORMAT, 'n_models_verified': 16,
                'absolute_spinor_normalization_certified': False,
                'complete_leading_one_loop_mass_boundary_certified': False,
                'source_snapshot_sha256': 'artificial-synthetic-fixture',
                'models': inventory,
            }
            (release/'summary.json').write_text(json.dumps(summary), encoding='utf-8')
            destination = root/'output/general_observable/lo_certification_spinor'
            first = certify_source_release(root, release, destination,
                                           strict_references=False, strict_counts=False)
            second = certify_source_release(root, release, destination,
                                            strict_references=False, strict_counts=False)
            self.assertEqual(first, second)
            self.assertEqual(first['n_models_structurally_verified'], 16)
            self.assertEqual(first['n_models_physically_normalized'], 0)
            self.assertEqual(first['n_models_ready_for_physical_EFT1_running'], 0)
            outdir = destination/('source_' + first['source_release_summary_sha256'][:16])
            self.assertEqual(len(list(outdir.glob('*/boundary.json'))), 16)
            for item in first['models']:
                with self.subTest(model=item['model']):
                    data = json.loads((outdir/tag(item['model'])/'boundary.json').read_text())
                    self.assertEqual(len(data['source_resolved_unscaled']), 1)
                    self.assertIsNone(data['normalization']['matchete_GammaCC_to_master_Weyl_complex_factor'])
                    self.assertIsNone(data['physical_master_weyl_boundary'])

    def test_evidence_located_is_not_phase_certificate(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'evidence.txt'
            path.write_text('GammaCC::usage = "charge conjugation matrix";\nPackageExport["GammaCC"]', encoding='utf-8')
            r = native_evidence_probe(path)
            self.assertFalse(r['GammaCC_Phase_certified'])
            self.assertIsNone(r['GammaCC_to_master_weyl_factor'])


if __name__ == '__main__':
    unittest.main()
