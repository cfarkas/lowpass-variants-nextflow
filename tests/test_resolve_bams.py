#!/usr/bin/env python3
"""Unit tests for deterministic BAM/sample resolution."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
RESOLVER = PACKAGE_DIR / "bin" / "resolve_bams.py"


class ResolveBamsTests(unittest.TestCase):
    def run_resolver(self, input_path: Path, output_path: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(RESOLVER),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def test_rejects_two_bams_with_the_same_logical_sample_id(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            tmp = Path(raw_tmp)
            (tmp / "SAMPLE.bam").write_bytes(b"bam-a")
            (tmp / "SAMPLE.sorted.bam").write_bytes(b"bam-b")
            output = tmp / "selected.tsv"

            result = self.run_resolver(tmp, output)

            self.assertNotEqual(0, result.returncode)
            self.assertIn("same logical sample ID", result.stderr)
            self.assertIn("SAMPLE", result.stderr)
            self.assertFalse(output.exists())

    def test_emits_one_row_per_unique_sample(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            tmp = Path(raw_tmp)
            first = tmp / "SAMPLE_A.bam"
            second = tmp / "SAMPLE_B.bqsr.bam"
            first.write_bytes(b"bam-a")
            second.write_bytes(b"bam-b")
            output = tmp / "selected.tsv"

            result = self.run_resolver(tmp, output)

            self.assertEqual(0, result.returncode, result.stderr)
            rows = output.read_text(encoding="utf-8").splitlines()
            self.assertEqual("sample\tbam", rows[0])
            self.assertEqual(3, len(rows))
            self.assertIn(f"SAMPLE_A\t{first.resolve()}", rows)
            self.assertIn(f"SAMPLE_B\t{second.resolve()}", rows)


if __name__ == "__main__":
    unittest.main(verbosity=2)
