"""Checks the standalone direct-one-loop CG inventory guard and generators."""
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from Numerical.direct_one_loop.PrepareSourceChargeAudits import (
    EXPECTED_CG_NAMES,
    _write_new_or_identical,
    generator_candidate,
    verify_source_inventory,
)
from Numerical.diagnostics.AnalyseMatcheteGenerators import analyse


class DirectOneLoopSourceAuditTests(unittest.TestCase):
    def test_test_fixture_generators_have_correct_su2_algebra(self):
        obj = analyse(generator_candidate())
        self.assertEqual(len(obj["representations"]), 2)
        self.assertEqual([e["status"] for e in obj["representations"]],
                         ["generator_algebra_verified"] * 2)
        self.assertFalse(generator_candidate()["generator_basis_is_live_export"])

    def test_inventory_rejects_wrong_alpha(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            bench = root / "Numerical/benchmarks/fitted_reference"
            cgroot = root / "source_cg_all16"
            bench.mkdir(parents=True)
            # A wrong-alpha export is rejected before any physics calculation.
            for i in range(16):
                name = f"T3_benchmark_{i:02d}"
                (bench / f"{name}.json").write_text(json.dumps({"representation": {
                    "d_s1": 1, "d_s2": 3, "d_f": 2, "alpha": i,
                }}), encoding="utf-8")
                path = cgroot / name
                path.mkdir(parents=True)
                (path / "uv_cg_registry.json").write_text(json.dumps({
                    "status": "Success",
                    "metadata": {"dS1": 1, "dS2": 3, "dF": 2, "alpha": i+1},
                    "CGRegistry": [{"Name": k} for k in sorted(EXPECTED_CG_NAMES)],
                }), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "metadata mismatch"):
                verify_source_inventory(root, cgroot, {"cases": [], "errors": []})

    def test_existing_audit_never_overwritten(self):
        with TemporaryDirectory() as temp:
            p = Path(temp) / "audit.json"
            _write_new_or_identical(p, {"value": 1})
            _write_new_or_identical(p, {"value": 1})
            with self.assertRaises(FileExistsError):
                _write_new_or_identical(p, {"value": 2})
            self.assertEqual(json.loads(p.read_text()), {"value": 1})


if __name__ == "__main__":
    unittest.main()
