"""Regression: all-16 runner must not overwrite its scalar-background directory.

In the original runner `old` held a Path before the loop, but the first
model's result was assigned to `old`, an ExtractedScalar. Every subsequent
model then failed with `TypeError: unsupported operand type(s) for /`.
The mocked all-16 test exercises the real main() loop without external data.
"""
from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np

from Numerical.lo_certification.completeness.inventory import ALL_MODELS
from Numerical.lo_certification.completeness.scalar_bridge import ExtractedScalar
from scripts import RunLOPortalCompletionResponse as runner


EXTRA_COUNTS = {
    'A_-4': 0, 'A_-2': 2, 'A_+0': 2,
    'B_-3': 1, 'B_-1': 4, 'B_+1': 1,
    'C_-3': 1, 'C_-1': 4, 'C_+1': 1,
    'D_-2': 2, 'D_+0': 2, 'D_+2': 0,
    'E_-4': 0, 'E_-2': 2, 'E_+0': 2, 'E_+2': 0,
}


class PortalRunnerLoopRegressionTests(unittest.TestCase):
    def test_check_only_all16_keeps_input_path_after_first_model(self):
        self.assertEqual(set(EXTRA_COUNTS), set(ALL_MODELS))
        self.assertEqual(sum(EXTRA_COUNTS.values()), 24)
        baseline_counts = {m: (2 if m.startswith('E_') else 1 if
                               m.startswith(('A_', 'D_')) else 0) for m in ALL_MODELS}
        self.assertEqual(sum(baseline_counts.values()), 14)

        fake_source_pkg = types.ModuleType('Numerical.joint_scalar_matching')
        fake_source_pkg.__path__ = []
        fake_source = types.ModuleType('Numerical.joint_scalar_matching.fermion_boundary_bridge')
        fake_source.load_release_manifest = lambda path: {m: Path('unused') for m in ALL_MODELS}
        observed_paths: list[Path] = []

        def fake_model_input(model, release, running_root, old_root, v3, v4, allow_special):
            # This fails with exactly the reported TypeError if old_root is an
            # ExtractedScalar produced by an earlier loop iteration.
            actual = old_root / 'scalar_backgrounds' / (runner.slug(model) + '.json')
            observed_paths.append(actual)
            scalar = ExtractedScalar(
                model=model, mu_f_gev=1e5, mu_s_gev=1e3, labels=(),
                lambda_channels_source=np.empty((0,)),
                lambda_source=np.empty((0,)), mass2_source=np.empty((0,)),
                channel_names=(), metadata_sha256='source', arrays_sha256='arrays',
                provenance={})
            return (None, scalar, None, None, (), np.empty((0,)), [], {
                'model': model,
                'baseline_completion': {'added_baseline_directions': baseline_counts[model]},
                'extra_named_inventory': {'extra_named_symbols': EXTRA_COUNTS[model]},
                'physical_C5_certified': False,
            })

        with patch.dict(sys.modules, {
             'Numerical.joint_scalar_matching': fake_source_pkg,
             'Numerical.joint_scalar_matching.fermion_boundary_bridge': fake_source,
        }), patch.object(runner, 'source_inventory', return_value={'source':'mock'}), \
             patch.object(runner, 'model_input', side_effect=fake_model_input), \
             redirect_stdout(StringIO()) as stdout:
            code = runner.main(['--check-only'])

        self.assertEqual(code, 0, stdout.getvalue())
        self.assertEqual(len(observed_paths), 16)
        root = runner.resolved(runner.OLD)
        for model, path in zip(ALL_MODELS, observed_paths):
            self.assertEqual(path, root / 'scalar_backgrounds' / (runner.slug(model) + '.json'))
        self.assertIn('"validated_input_models": 16', stdout.getvalue())
        self.assertIn('"blocked_source_models": 0', stdout.getvalue())
        self.assertIn('"baseline_directions_verified": 14', stdout.getvalue())
        self.assertIn('"extra_named_symbols_verified": 24', stdout.getvalue())


if __name__ == '__main__':
    unittest.main()
