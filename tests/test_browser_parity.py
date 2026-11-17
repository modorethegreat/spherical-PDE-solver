import json
import shutil
import subprocess
import unittest
from pathlib import Path

from mode_benchmark import ModeCase, benchmark


@unittest.skipUnless(shutil.which("node"), "Node is required to check browser/Python parity")
class BrowserParityTests(unittest.TestCase):
    def test_static_browser_uses_the_same_numerical_model_and_export_schema(self):
        cases = [(4, 2, 40, 0.5), (12, 12, 320, 0), (1, 0, 80, 4), (2, 1, 20, 4)]
        module = Path(__file__).resolve().parents[1] / "docs" / "benchmark.js"
        script = """
import fs from 'node:fs';
const {modeBenchmark} = await import('data:text/javascript,' + encodeURIComponent(fs.readFileSync(process.argv[1], 'utf8')));
console.log(JSON.stringify(JSON.parse(fs.readFileSync(0, 'utf8')).map(args => modeBenchmark(...args))));
"""
        actual = json.loads(subprocess.check_output(
            ["node", "--input-type=module", "-e", script, str(module)],
            input=json.dumps(cases), text=True,
        ))
        for got, case in zip(actual, cases, strict=True):
            with self.subTest(case=case):
                self.check_values(got, benchmark(ModeCase(*case)))

    def check_values(self, actual, expected):
        if isinstance(expected, dict):
            self.assertEqual(actual.keys(), expected.keys())
            for key in expected:
                self.check_values(actual[key], expected[key])
        elif isinstance(expected, list):
            for x, y in zip(actual, expected, strict=True):
                self.check_values(x, y)
        elif isinstance(expected, float):
            self.assertAlmostEqual(actual, expected, delta=2e-12)
        else:
            self.assertEqual(actual, expected)
