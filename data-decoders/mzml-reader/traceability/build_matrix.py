"""Build the requirements -> tests -> results traceability matrix from JUnit XML.

Standard library only. Inputs:

* ``requirements.tsv`` (this folder): one row per requirement (id, area,
  requirement, source);
* one or more JUnit XML files written by ``pytest --junitxml=...`` in this
  folder's test suite. ``conftest.py`` copies each test's
  ``@pytest.mark.req(...)`` ids into the XML as ``requirement`` properties, and
  ``pytest.ini`` selects the ``xunit1`` format that keeps them.

Output: a Markdown matrix with, for every requirement, the tests that verify it
and their results; requirements without tests; tests without a requirement.
The header records when and from what it was generated, so a committed copy is
a dated snapshot of one real run, never a hand-written table.

Usage (from data-decoders/mzml-reader/)::

    python -m pytest -q --junitxml=junit.xml
    python traceability/build_matrix.py junit.xml --out docs/traceability-matrix.md

Exit status: 0 when the matrix was written (even if tests failed: the matrix
shows it), 2 for unreadable input or a marker naming an unknown requirement.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

HERE = Path(__file__).resolve().parent
REQUIREMENTS = HERE / "requirements.tsv"

__all__ = ["Requirement", "TestResult", "Matrix", "load_requirements", "load_junit", "build", "render", "main"]


class MatrixError(ValueError):
    """Unreadable input, or a test marker naming an unknown requirement."""


@dataclass
class Requirement:
    id: str
    area: str
    text: str
    source: str


@dataclass
class TestResult:
    node: str                  # e.g. test_mzml.py::test_binary_round_trip[32-False]
    outcome: str               # passed, failed, error or skipped
    requirements: List[str]


@dataclass
class Matrix:
    requirements: List[Requirement]
    tests: List[TestResult]
    by_requirement: Dict[str, List[TestResult]] = field(default_factory=dict)
    untraced_tests: List[TestResult] = field(default_factory=list)
    run_timestamp: Optional[str] = None

    def status(self, req_id: str) -> str:
        tests = self.by_requirement.get(req_id, [])
        if not tests:
            return "not tested"
        if any(t.outcome in ("failed", "error") for t in tests):
            return "FAILED"
        if all(t.outcome == "skipped" for t in tests):
            return "skipped"
        return "verified"


def load_requirements(path: Path = REQUIREMENTS) -> List[Requirement]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    reqs, seen = [], set()
    for row in rows:
        rid = (row.get("id") or "").strip()
        if not rid or rid in seen:
            raise MatrixError(f"{path}: missing or duplicate requirement id {rid!r}")
        seen.add(rid)
        reqs.append(Requirement(rid, row["area"].strip(), row["requirement"].strip(), row["source"].strip()))
    return reqs


def _node(classname: str, name: str) -> str:
    """'test_mzml' + 'test_x[1]' -> 'test_mzml.py::test_x[1]'."""
    parts = classname.split(".")
    for i in range(len(parts), 0, -1):              # module path, then optional class names
        if parts[i - 1].startswith("test"):
            return "/".join(parts[:i]) + ".py::" + "::".join(parts[i:] + [name])
    return classname.replace(".", "/") + "::" + name


def load_junit(paths: Sequence[Path]):
    """All test cases from the JUnit files, and the earliest run timestamp."""
    tests: List[TestResult] = []
    stamps = []
    for path in paths:
        try:
            root = ET.parse(str(path)).getroot()
        except (ET.ParseError, OSError) as exc:
            raise MatrixError(f"{path}: not a readable JUnit XML file ({exc})") from None
        suites = [root] if root.tag == "testsuite" else root.iter("testsuite")
        for suite in suites:
            if suite.get("timestamp"):
                stamps.append(suite.get("timestamp"))
            for case in suite.iter("testcase"):
                outcome = "passed"
                for tag in ("failure", "error", "skipped"):
                    if case.find(tag) is not None:
                        outcome = {"failure": "failed", "error": "error", "skipped": "skipped"}[tag]
                        break
                reqs = [p.get("value", "") for p in case.iter("property") if p.get("name") == "requirement"]
                tests.append(TestResult(_node(case.get("classname", ""), case.get("name", "")), outcome, reqs))
    return tests, (min(stamps) if stamps else None)


def build(requirements: List[Requirement], tests: List[TestResult], run_timestamp: Optional[str] = None) -> Matrix:
    known = {r.id for r in requirements}
    unknown = sorted({rid for t in tests for rid in t.requirements if rid not in known})
    if unknown:
        raise MatrixError(f"test markers name unknown requirement ids: {', '.join(unknown)}")
    matrix = Matrix(requirements, sorted(tests, key=lambda t: t.node), run_timestamp=run_timestamp)
    for r in requirements:
        matrix.by_requirement[r.id] = []
    for t in matrix.tests:
        if not t.requirements:
            matrix.untraced_tests.append(t)
        for rid in dict.fromkeys(t.requirements):        # each test once per requirement
            matrix.by_requirement[rid].append(t)
    return matrix


def _md(text: str) -> str:
    return text.replace("|", "\\|").replace("<", "&lt;").replace(">", "&gt;")


def render(matrix: Matrix, generated: str, source: str, python: str) -> str:
    counts = {k: sum(1 for t in matrix.tests if t.outcome == k) for k in ("passed", "failed", "error", "skipped")}
    statuses = {r.id: matrix.status(r.id) for r in matrix.requirements}
    tally = {k: sum(1 for s in statuses.values() if s == k) for k in ("verified", "FAILED", "skipped", "not tested")}
    lines = [
        "# Traceability matrix: mzML tools",
        "",
        "Generated by `traceability/build_matrix.py` from the JUnit XML of a real test run. Do not edit by hand: regenerate it.",
        "",
        "| | |",
        "|---|---|",
        f"| Generated (UTC) | {generated} |",
        f"| Test run started | {matrix.run_timestamp or 'not recorded'} |",
        f"| Source | {_md(source)} |",
        f"| Python | {python} |",
        f"| Tests | {len(matrix.tests)}: {counts['passed']} passed, {counts['failed']} failed, "
        f"{counts['error']} errors, {counts['skipped']} skipped |",
        f"| Requirements | {len(matrix.requirements)}: {tally['verified']} verified, {tally['FAILED']} with failures, "
        f"{tally['skipped']} only skipped, {tally['not tested']} not tested |",
        "",
        "A requirement is **verified** when at least one test is linked to it and all its linked tests passed in this run. "
        "Requirements come from `traceability/requirements.tsv`; links come from each test's `@pytest.mark.req(...)` marker.",
        "",
        "## Summary",
        "",
        "| ID | Area | Requirement | Tests | Passed | Failed / error | Skipped | Status |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in matrix.requirements:
        tests = matrix.by_requirement[r.id]
        p = sum(t.outcome == "passed" for t in tests)
        f = sum(t.outcome in ("failed", "error") for t in tests)
        s = sum(t.outcome == "skipped" for t in tests)
        status = statuses[r.id]
        shown = f"**{status}**" if status in ("FAILED", "not tested") else status
        lines.append(f"| {r.id} | {r.area} | {_md(r.text)} | {len(tests)} | {p} | {f} | {s} | {shown} |")
    lines += ["", "## Requirements without tests", ""]
    missing = [r for r in matrix.requirements if not matrix.by_requirement[r.id]]
    lines += [f"- **{r.id}**: {_md(r.text)} (source: {_md(r.source)})" for r in missing] or ["None."]
    lines += ["", "## Tests without a requirement marker", ""]
    lines += [f"- `{t.node}` ({t.outcome})" for t in matrix.untraced_tests] or ["None."]
    lines += ["", "## Tests per requirement", ""]
    for r in matrix.requirements:
        tests = matrix.by_requirement[r.id]
        if not tests:
            continue
        lines += [f"### {r.id}: {_md(r.text)}", "", f"Source: {_md(r.source)}", ""]
        lines += [f"- `{t.node}`: {t.outcome}" for t in tests]
        lines.append("")
    return "\n".join(lines).rstrip("\n") + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the traceability matrix from JUnit XML.")
    ap.add_argument("junit", nargs="+", type=Path, help="JUnit XML files from pytest --junitxml")
    ap.add_argument("--requirements", type=Path, default=REQUIREMENTS)
    ap.add_argument("--out", type=Path, help="write the Markdown here (default: print it)")
    ap.add_argument("--source", help="what was tested, e.g. a commit id (default: $GITHUB_SHA, else 'not recorded')")
    ap.add_argument("--python", help="Python version of the test run (default: this interpreter's)")
    ap.add_argument("--generated", help="generation time to record (default: now, UTC)")
    args = ap.parse_args(argv)
    try:
        reqs = load_requirements(args.requirements)
        tests, stamp = load_junit(args.junit)
        matrix = build(reqs, tests, stamp)
    except MatrixError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    source = args.source or (f"commit {os.environ['GITHUB_SHA']}" if os.environ.get("GITHUB_SHA") else "not recorded")
    text = render(matrix, args.generated or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), source,
                  args.python or sys.version.split()[0])
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_bytes(text.encode("utf-8"))
        print(f"wrote {args.out}: {len(matrix.tests)} tests, {len(reqs)} requirements")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
