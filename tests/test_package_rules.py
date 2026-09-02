#!/usr/bin/env python3
"""Static release rules for the integrated low-pass variants package."""

from __future__ import annotations

import re
import subprocess
import sys
import unittest
from pathlib import Path


sys.dont_write_bytecode = True

PACKAGE_DIR = Path(__file__).resolve().parents[1]
MAIN = PACKAGE_DIR / "main.nf"
CONFIG = PACKAGE_DIR / "nextflow.config"
BUILDER = PACKAGE_DIR / "bin" / "build_final_per_sample_vcf.py"
COMMON = PACKAGE_DIR / "bin" / "nf_common.sh"
FILES_MANIFEST = PACKAGE_DIR / "FILES.txt"


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def without_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"(?m)^\s*//.*$", "", text)


class PackageRulesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.main = read(MAIN)
        cls.main_code = without_comments(cls.main)
        cls.config = read(CONFIG)
        cls.builder = read(BUILDER)
        cls.common = read(COMMON)

    def test_ffpe_and_fresh_modes_are_declared_and_xor_validated(self) -> None:
        self.assertRegex(self.config, r"(?m)^\s*ffpe\s*=\s*false\b")
        self.assertRegex(self.config, r"(?m)^\s*fresh\s*=\s*false\b")
        self.assertRegex(self.main_code, r"truthyParam\(\s*params\.ffpe\s*\)")
        self.assertRegex(self.main_code, r"truthyParam\(\s*params\.fresh\s*\)")

        error_messages = " ".join(
            re.findall(r"\berror\s+[\"']([^\"']+)[\"']", self.main_code)
        )
        self.assertIn("--ffpe", error_messages)
        self.assertIn("--fresh", error_messages)
        self.assertRegex(
            error_messages.lower(),
            r"exactly\s+one|mutually\s+exclusive|choose\s+one",
            "workflow must reject both-mode and neither-mode invocations",
        )
        self.assertIn('add_argument("--sample-mode"', self.builder)

    def test_bqsr_is_real_and_explicitly_skippable(self) -> None:
        self.assertRegex(self.config, r"(?m)^\s*skip_bqsr\s*=\s*false\b")
        self.assertRegex(self.config, r"(?m)^\s*known_sites\s*=")
        self.assertIn("params.skip_bqsr", self.main_code)
        self.assertIn("params.known_sites", self.main_code)
        self.assertRegex(self.main_code, r"\bgatk\s+BaseRecalibrator\b")
        self.assertRegex(self.main_code, r"\bgatk\s+ApplyBQSR\b")
        self.assertRegex(self.main_code, r"--known-sites\b")

    def test_ffperase_process_is_called_only_from_an_ffpe_branch(self) -> None:
        process_names = re.findall(
            r"(?m)^\s*process\s+([A-Z0-9_]*FFPERASE[A-Z0-9_]*)\s*\{",
            self.main_code,
        )
        self.assertTrue(process_names, "no FFPErase process declaration found")

        lines = self.main_code.splitlines()
        branch_snippets = []
        for index, line in enumerate(lines):
            if re.search(r"\bif\b.*ffpe", line, flags=re.IGNORECASE):
                branch_snippets.append("\n".join(lines[index : index + 400]))

        self.assertTrue(branch_snippets, "no executable FFPE-mode conditional found")
        invoked = any(
            re.search(rf"(?m)^\s*{re.escape(name)}\s*\(", snippet)
            for name in process_names
            for snippet in branch_snippets
        )
        self.assertTrue(invoked, "FFPErase process is not invoked inside an FFPE branch")

    def test_normalization_uses_the_reference(self) -> None:
        self.assertRegex(self.common, r"\bbcftools\s+norm\b")
        self.assertRegex(self.common, r"(?m)^\s*-f\s+[\"']?\$ref")

    def test_literal_project_files_and_manifest_entries_exist(self) -> None:
        references = set(
            re.findall(r"\$\{projectDir\}/((?:assets|bin|envs)/[A-Za-z0-9_./-]+)", self.main)
        )
        self.assertTrue(references, "no project-local runtime references discovered")
        missing = sorted(path for path in references if not (PACKAGE_DIR / path).exists())
        self.assertEqual([], missing, f"missing projectDir references: {missing}")

        self.assertIn("envs/lowpass_variants.yml", references)
        self.assertNotIn("envs/ffpe_lowpass.yml", references)

        manifest_entries = []
        for raw_line in read(FILES_MANIFEST).splitlines():
            entry = raw_line.strip()
            if entry and not entry.startswith("#"):
                manifest_entries.append(entry)
        missing_manifest = sorted(
            entry for entry in manifest_entries if not (PACKAGE_DIR / entry).exists()
        )
        self.assertEqual([], missing_manifest, f"FILES.txt contains missing paths: {missing_manifest}")

    def test_reserved_terminology_is_absent(self) -> None:
        reserved_terms = ("con" + "tract", "mini" + "map2")
        git_probe = subprocess.run(
            ["git", "-C", str(PACKAGE_DIR), "rev-parse", "--is-inside-work-tree"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if git_probe.returncode == 0:
            manifest_entries = subprocess.run(
                ["git", "-C", str(PACKAGE_DIR), "ls-files", "-z"],
                check=True,
                stdout=subprocess.PIPE,
            ).stdout.decode("utf-8").split("\0")
        else:
            manifest_entries = [
                line.strip()
                for line in read(FILES_MANIFEST).splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            ]
        matches = []
        for entry in filter(None, manifest_entries):
            path = PACKAGE_DIR / entry
            if any(term in entry.lower() for term in reserved_terms):
                matches.append(entry)
                continue
            try:
                content = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            if any(term in content.lower() for term in reserved_terms):
                matches.append(entry)
        self.assertEqual([], matches, "reserved terminology found in: " + ", ".join(matches))

    def test_container_release_files_and_docker_profile_exist(self) -> None:
        required = [
            Path("Dockerfile"),
            Path(".dockerignore"),
            Path("docker/environment.yml"),
            Path(".github/workflows/container.yml"),
        ]
        missing = sorted(str(path) for path in required if not (PACKAGE_DIR / path).is_file())
        self.assertEqual([], missing, "missing container release files: " + ", ".join(missing))

        dockerfile = read(PACKAGE_DIR / "Dockerfile")
        workflow = read(PACKAGE_DIR / ".github/workflows/container.yml")
        from_lines = [line.strip() for line in dockerfile.splitlines() if line.startswith("FROM ")]
        self.assertTrue(from_lines, "Dockerfile has no base image")
        for line in from_lines:
            self.assertRegex(line, r"@sha256:[0-9a-f]{64}(?:\s|$)")

        self.assertRegex(self.config, r"(?m)^\s*docker\s*\{\s*$")
        self.assertRegex(self.config, r"(?m)^\s*docker\.enabled\s*=\s*true\s*$")
        self.assertIn("container = params.pipeline_container", self.config)
        self.assertRegex(
            self.config,
            r"(?m)^\s*pipeline_container\s*=\s*['\"]ghcr\.io/",
        )
        self.assertIn("packages: write", workflow)
        self.assertIn("${{ secrets.GITHUB_TOKEN }}", workflow)
        self.assertIn("linux/amd64", workflow)
        self.assertIn("org.opencontainers.image.source", dockerfile)

    def test_no_generated_runtime_artifacts_are_release_candidates(self) -> None:
        git_probe = subprocess.run(
            ["git", "-C", str(PACKAGE_DIR), "rev-parse", "--is-inside-work-tree"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if git_probe.returncode == 0:
            tracked = subprocess.run(
                ["git", "-C", str(PACKAGE_DIR), "ls-files", "-z"],
                check=True,
                stdout=subprocess.PIPE,
            ).stdout.decode("utf-8").split("\0")
            candidates = [Path(path) for path in tracked if path]
        else:
            # Before git init, FILES.txt is the explicit release-candidate set.
            # Local ignored Nextflow work/cache files must not make this test
            # depend on the developer's runtime state.
            candidates = [
                Path(line.strip())
                for line in read(FILES_MANIFEST).splitlines()
                if line.strip() and not line.lstrip().startswith("#")
            ]

        forbidden_parts = {".nextflow", ".conda", "work", "__pycache__", ".pytest_cache"}

        def stale(path: Path) -> bool:
            if forbidden_parts.intersection(path.parts):
                return True
            name = path.name
            return (
                name.startswith(".nextflow.log")
                or name.endswith((".pyc", ".pyo"))
                or name in {"trace.txt", "timeline.html", "report.html", "dag.html"}
            )

        stale_candidates = sorted(str(path) for path in candidates if stale(path))
        self.assertEqual(
            [],
            stale_candidates,
            "generated runtime artifacts must be removed or remain untracked: "
            + ", ".join(stale_candidates),
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
