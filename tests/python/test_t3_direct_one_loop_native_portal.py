"""Pure regression tests of independent native-portal projection and source guards."""
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
import sympy as sp

from scripts.AuditT3DirectUVPortalBasis import (
    project_model, projected_neutral_portal, source_contract,
)


class DirectUVNativePortalTests(unittest.TestCase):
    def test_higgs_neutral_portal_projects_diagonal_charge_masses(self):
        a = np.zeros((2, 2, 3, 3), dtype=complex)
        a[1, 1] = np.diag([0.75, 1.0, 1.25])
        self.assertEqual(projected_neutral_portal(a, 3), [0.75, 1.0, 1.25])
        self.assertNotEqual(projected_neutral_portal(a, 3), [1., 1., 1.])

    def test_reject_charge_offdiagonal_and_apply_plus_hc_to_complex_diagonal(self):
        arr = np.zeros((2, 2, 2, 2), dtype=complex)
        arr[1, 1, 0, 1] = 0.1
        with self.assertRaisesRegex(ValueError, "off-diagonal"):
            projected_neutral_portal(arr, 2)
        arr[1, 1, 0, 1] = 0
        # Imaginary native CG components cancel in source 1/2 PlusHc.
        arr[1, 1, 0, 0] = 0.1j
        self.assertEqual(projected_neutral_portal(arr, 2), [0.0, 0.0])

    def test_source_contract_fails_closed_when_prefactor_modified(self):
        with TemporaryDirectory() as dirname:
            root = Path(dirname)
            base = root / "Lagrangian/interactions"
            base.mkdir(parents=True)
            portal = (
                'BuildHiggsPortalInvariantCandidates[which_] := Module[{},\n'
                '"T3HiggsPortal" <> ToString[which] <> "CG";\n'
                '"lambdaH" <> ToString[which]; {2, 2, d, d};\n'
                '{True, False, True, False}; {1, 2};\n'
                '"Expression" -> 1/2 PlusHc[coupling[] fieldProduct];\n'
                'BuildCrossScalarInvariantCandidates[which_] := 1;\n'
            )
            (base / "ScalarPotential.wl").write_text(portal, encoding="utf8")
            (base / "T3Topology.wl").write_text(
                'lambdaT3[] H[i] H[j] s1 s2 CG[T3MixCG, labels]', encoding="utf8")
            self.assertEqual(len(source_contract(root)), 2)
            (base / "ScalarPotential.wl").write_text(
                portal.replace('1/2 PlusHc[', '2 PlusHc['), encoding="utf8")
            with self.assertRaisesRegex(ValueError, "source changed"):
                source_contract(root)

    def test_metadata_rejects_wrong_representation_before_parsing_native_tensors(self):
        obj = {"status": "Success", "format": "t3-native-direct-uv-portal-cg-v1",
               "metadata": {"dS1": 2, "dS2": 2, "dF": 1, "alpha": 1},
               "CGRegistry": []}
        with self.assertRaisesRegex(ValueError, "metadata mismatch"):
            project_model(obj, "B_-1", (2, 2, 1, -1),
                          {1: sp.eye(1), 2: sp.eye(2), 3: sp.eye(3)},
                          {"lambdaH1": .02, "lambdaH2": .03}, 246.22,
                          Path("unread-synthetic-export.json"))


if __name__ == '__main__':
    unittest.main()
