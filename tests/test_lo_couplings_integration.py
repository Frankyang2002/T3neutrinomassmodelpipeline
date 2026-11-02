"""All-16 real-repository gates (skipped only if repository sources are absent)."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from Numerical.lo_certification.couplings.batch import (
    model_paths, check_source_versions, read_json, observed_running_states,
    calculate_candidate_quartics, nonzero_sources, SOURCE_BLOBS,
    local_source_recheck, RECHECKABLE_PATHS)

ROOT=Path(__file__).resolve().parents[1]
RELEASE=ROOT/'output/general_observable/lo_fermion_boundary/releases/source_021dcc1e8fe8ba0f'
CACHE=ROOT/'output/general_observable/lo_joint_scalar_matching/source_quartic_v4/native_quartic_channels'
REGISTRY=ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json'


@unittest.skipUnless((RELEASE/'summary.json').is_file() and REGISTRY.is_file(),
                     'Real local repository/release not mounted: algebra-only tests still run')
class TestRealRepositoryAll16(unittest.TestCase):
    def test_all_16_matchete_release_and_native_CG_cache(self):
        from Numerical.joint_scalar_matching.adapter import ALL_MODELS
        from Numerical.joint_scalar_matching.fermion_boundary_bridge import (
            load_release_manifest, load_source)
        from Numerical.joint_scalar_matching.numerical_benchmark_v5 import (
            read_v4_cached_channels)
        from Numerical.joint_scalar_matching.source_quartic_v4 import POTENTIAL_CONVENTION
        self.assertEqual(len(ALL_MODELS),16)
        refs=load_release_manifest(RELEASE)
        registry=read_json(REGISTRY)
        self.assertEqual(set(refs),set(ALL_MODELS))
        for model in ALL_MODELS:
            with self.subTest(model=model):
                p=model_paths(ROOT,model,RELEASE,CACHE)
                src=load_source(p['boundary'],expected_model=model)
                ch=read_v4_cached_channels(p['quartic_cache'],p['quartic_seed'],src)
                self.assertEqual(len(ch.source_basis_labels),len(src.labels))
                self.assertFalse(src.normalization_certified)
                self.assertGreater(len(ch.by_name),0)
                self.assertTrue((p['running']).is_file())
                doc=observed_running_states(p['running'],model)
                self.assertEqual(doc['scales_gev']['F'],src.mu_f_gev)
                self.assertGreater(doc['scales_gev']['F'],doc['scales_gev']['S'])
                extra=nonzero_sources(registry['models'][model],set(ch.by_name))
                self.assertTrue(all(x['represented_as_certified_native_HHSS'] is False for x in extra))

    def test_local_source_drift_is_reaudited_before_output(self):
        # Unlike a byte-for-byte GitHub reproducibility check, this tests the
        # locally installed sources. A changed parser must reproduce the native
        # seed-to-CG tensors for every one of the sixteen models.
        from Numerical.joint_scalar_matching.adapter import ALL_MODELS
        observed=check_source_versions(ROOT)
        self.assertEqual(set(observed),set(SOURCE_BLOBS))
        drift={rel for rel,row in observed.items() if not row['matches_inspected_main']}
        self.assertFalse(drift-RECHECKABLE_PATHS,
                         'Other source files changed: '+str(sorted(drift-RECHECKABLE_PATHS)))
        registry=read_json(REGISTRY)
        result=local_source_recheck(ROOT,tuple(ALL_MODELS),observed,RELEASE,CACHE,registry)
        self.assertEqual(set(result['drift']),drift)
        self.assertFalse(result['physical_source_certified'])
        if 'RGE/running/intermediate/ScalarOnlyTensorAdapters.py' in drift:
            self.assertEqual(set(result['parser_reparsed_models']),set(ALL_MODELS))
        if 'RGE/general/generated/general_rgbeta_quartic_registry.json' in drift:
            self.assertEqual(set(result['registry_models_checked']),set(ALL_MODELS))


if __name__=='__main__':
    unittest.main()
