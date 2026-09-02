#!/usr/bin/env python3
"""Behavioral unit test for fresh versus FFPE C>T filtering."""

from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


sys.dont_write_bytecode = True

PACKAGE_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = PACKAGE_DIR / "bin" / "build_final_per_sample_vcf.py"
SPEC = importlib.util.spec_from_file_location("build_final_per_sample_vcf", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"could not load {MODULE_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def arguments(sample_mode: str) -> SimpleNamespace:
    return SimpleNamespace(
        sample_mode=sample_mode,
        min_dp=10,
        min_alt_reads=4,
        min_af=0.05,
        vote_threshold=1,
        ffpe_ct_vaf_keep=0.35,
        ffpe_ct_min_alt=8,
    )


def low_vaf_ct_call() -> dict[str, object]:
    return {
        "caller": "MUTECT2",
        "filter": "PASS",
        "caller_pass": True,
        "status": "completed",
        "dp": 100,
        "alt_reads": 20,
        "af_observed": 0.20,
        "af_calculated": 0.20,
    }


class SampleModeTests(unittest.TestCase):
    def test_fresh_bypasses_ffpe_ct_review_while_ffpe_applies_it(self) -> None:
        calls = {"MUTECT2": low_vaf_ct_call()}

        fresh = MODULE.decide_variant("C", "T", calls, arguments("fresh"))
        ffpe = MODULE.decide_variant("C", "T", calls, arguments("ffpe"))

        self.assertEqual("HIGH_CONFIDENCE", fresh["decision"])
        self.assertNotIn("ffpe_risk", fresh)
        self.assertNotIn("ffpe_low_support", fresh)
        self.assertNotIn("FFPE_CtoT_GtoA_low_support", fresh["drop_reason"])

        self.assertEqual("MODERATE_REVIEW", ffpe["decision"])
        self.assertEqual("yes", ffpe["ffpe_low_support"])
        self.assertIn("FFPE_CtoT_GtoA_low_support", ffpe["drop_reason"])

    def test_fresh_cli_omits_ffpe_tsv_columns_and_vcf_info(self) -> None:
        with tempfile.TemporaryDirectory(prefix="lowpass-sample-mode-") as tmp_name:
            tmp = Path(tmp_name)
            input_vcf = tmp / "input.vcf"
            ref_fai = tmp / "tiny.fa.fai"
            input_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=4>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr1\t2\t.\tC\tT\t60\tPASS\tDP=100;AF=0.20\tGT:DP:AD:AF\t0/1:100:80,20:0.20\n",
                encoding="utf-8",
            )
            ref_fai.write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")

            fresh_dir = tmp / "fresh"
            fresh_vcf = fresh_dir / "SAMPLE.final.vcf"
            base_command = [
                sys.executable,
                str(MODULE_PATH),
                "--sample",
                "SAMPLE",
                "--mutect-vcf",
                str(input_vcf),
                "--mutect-status",
                "completed",
                "--ref-fai",
                str(ref_fai),
                "--min-dp",
                "10",
                "--min-alt-reads",
                "4",
                "--min-af",
                "0.05",
                "--vote-threshold",
                "1",
            ]
            subprocess.run(
                base_command
                + [
                    "--sample-mode",
                    "fresh",
                    "--outdir",
                    str(fresh_dir),
                    "--out-vcf",
                    str(fresh_vcf),
                ],
                check=True,
            )

            fresh_tsv = fresh_dir / "SAMPLE.final_variants.all_candidates.tsv"
            with fresh_tsv.open(newline="", encoding="utf-8") as handle:
                fresh_header = next(csv.reader(handle, delimiter="\t"))
            fresh_vcf_text = fresh_vcf.read_text(encoding="utf-8")
            self.assertFalse(any(column.lower().startswith("ffpe") for column in fresh_header))
            self.assertNotIn("is_CT_or_GA", fresh_header)
            self.assertNotIn("FFPE_RISK", fresh_vcf_text)
            self.assertNotIn("FFPE_LOW_SUPPORT", fresh_vcf_text)

            ffpe_dir = tmp / "ffpe"
            ffpe_vcf = ffpe_dir / "SAMPLE.final.vcf"
            subprocess.run(
                base_command
                + [
                    "--sample-mode",
                    "ffpe",
                    "--ffpe-ct-vaf-keep",
                    "0.35",
                    "--ffpe-ct-min-alt",
                    "8",
                    "--outdir",
                    str(ffpe_dir),
                    "--out-vcf",
                    str(ffpe_vcf),
                ],
                check=True,
            )

            ffpe_tsv = ffpe_dir / "SAMPLE.final_variants.all_candidates.tsv"
            with ffpe_tsv.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                ffpe_header = reader.fieldnames or []
                ffpe_row = next(reader)
            ffpe_vcf_text = ffpe_vcf.read_text(encoding="utf-8")
            self.assertIn("ffpe_risk", ffpe_header)
            self.assertIn("ffpe_low_support", ffpe_header)
            self.assertEqual("MODERATE_REVIEW", ffpe_row["decision"])
            self.assertIn("FFPE_RISK", ffpe_vcf_text)
            self.assertIn("FFPE_LOW_SUPPORT", ffpe_vcf_text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
