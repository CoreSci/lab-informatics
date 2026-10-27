"""Unit tests for traceability/build_matrix.py.

Run from this folder: ``python -m pytest -q``
"""
import csv
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "traceability"))
import build_matrix as bm  # noqa: E402

JUNIT = """<?xml version="1.0" encoding="utf-8"?>
<testsuites><testsuite name="pytest" tests="4" timestamp="2026-01-02T03:04:05">
  <testcase classname="test_demo" name="test_ok">
    <properties><property name="requirement" value="A-01"/></properties></testcase>
  <testcase classname="test_demo" name="test_bad[1]">
    <properties><property name="requirement" value="A-02"/><property name="requirement" value="A-01"/></properties>
    <failure message="boom">trace</failure></testcase>
  <testcase classname="test_demo" name="test_skip">
    <properties><property name="requirement" value="A-03"/></properties><skipped message="later"/></testcase>
  <testcase classname="sub.test_other" name="test_untagged"/>
</testsuite></testsuites>
"""


@pytest.fixture
def files(tmp_path):
    req = tmp_path / "req.tsv"
    with open(req, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "area", "requirement", "source"])
        for i in range(1, 5):
            w.writerow([f"A-0{i}", "demo", f"requirement {i} with a | pipe", "here"])
    junit = tmp_path / "junit.xml"
    junit.write_text(JUNIT, encoding="utf-8")
    return req, junit


@pytest.mark.req("TR-01")
def test_matrix_links_results_to_requirements(files):
    req, junit = files
    tests, stamp = bm.load_junit([junit])
    m = bm.build(bm.load_requirements(req), tests, stamp)
    assert stamp == "2026-01-02T03:04:05"
    assert [m.status(f"A-0{i}") for i in range(1, 5)] == ["FAILED", "FAILED", "skipped", "not tested"]
    assert [t.node for t in m.by_requirement["A-01"]] == ["test_demo.py::test_bad[1]", "test_demo.py::test_ok"]
    assert [t.node for t in m.untraced_tests] == ["sub/test_other.py::test_untagged"]


@pytest.mark.req("TR-01")
def test_rendered_matrix_lists_gaps_and_is_reproducible(files, tmp_path):
    req, junit = files
    out = tmp_path / "m.md"
    args = [str(junit), "--requirements", str(req), "--out", str(out), "--source", "commit abc123",
            "--python", "3.12.0", "--generated", "2026-01-02T04:00:00Z"]
    assert bm.main(args) == 0
    text = out.read_text(encoding="utf-8")
    assert "| Tests | 4: 2 passed, 1 failed, 0 errors, 1 skipped |" in text
    assert "| Requirements | 4: 0 verified, 2 with failures, 1 only skipped, 1 not tested |" in text
    assert "- **A-04**: requirement 4 with a \\| pipe (source: here)" in text
    assert "- `sub/test_other.py::test_untagged` (passed)" in text
    assert "| Source | commit abc123 |" in text
    first = out.read_bytes()
    assert bm.main(args) == 0 and out.read_bytes() == first


@pytest.mark.req("TR-01")
def test_unknown_requirement_ids_are_rejected(files, tmp_path, capsys):
    req, junit = files
    junit.write_text(JUNIT.replace('value="A-03"', 'value="Z-99"'), encoding="utf-8")
    assert bm.main([str(junit), "--requirements", str(req)]) == 2
    assert "Z-99" in capsys.readouterr().err


@pytest.mark.req("TR-01")
def test_unreadable_junit_is_an_error(files, tmp_path):
    req, _ = files
    bad = tmp_path / "bad.xml"
    bad.write_text("<testsuite>", encoding="utf-8")
    assert bm.main([str(bad), "--requirements", str(req)]) == 2


@pytest.mark.req("TR-01")
def test_shipped_requirements_are_well_formed():
    reqs = bm.load_requirements()
    assert len(reqs) == len({r.id for r in reqs}) >= 30
    assert all(r.area and r.text and r.source for r in reqs)
