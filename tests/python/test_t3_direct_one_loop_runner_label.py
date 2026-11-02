"""Regression: running provisional UV inputs must never print source certified."""
from contextlib import redirect_stdout
from io import StringIO
from unittest import TestCase
from unittest.mock import patch

from scripts import RunT3DirectOneLoop as runner


class DirectUVRunnerLabelTest(TestCase):
    def test_provisional_summary_does_not_claim_certification(self):
        report = {"evaluated": 16, "models": [
            {"status": "uv_mass_evaluated_uncertified_source_mapping"} for _ in range(16)
        ]}
        stream = StringIO()
        with patch.object(runner, "run", return_value=report), \
             patch("sys.argv", ["RunT3DirectOneLoop.py", "--allow-uncertified"]), \
             redirect_stdout(stream):
            code = runner.main()
        self.assertEqual(code, 0)
        self.assertIn("Direct UV numerical results: 16/16", stream.getvalue())
        self.assertIn("physically source-certified: 0/16", stream.getvalue())
        self.assertNotIn("Source-certified direct UV results: 16/16", stream.getvalue())
