"""Unit tests for jcamp_dx.py.

Run from this folder: ``python -m pytest -q``
"""
import math
from pathlib import Path

import numpy as np
import pytest

import jcamp_dx as jdx

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "sample-data"


def minimal(table_lines, **overrides):
    """A small valid XYDATA file (3 points, x = 100, 101, 102) with optional header overrides."""
    hdr = {"TITLE": "t", "JCAMP-DX": "4.24", "DATA TYPE": "INFRARED SPECTRUM", "ORIGIN": "test",
           "OWNER": "public domain", "XUNITS": "1/CM", "YUNITS": "ABSORBANCE", "XFACTOR": "1",
           "YFACTOR": "1", "FIRSTX": "100", "LASTX": "102", "NPOINTS": "3", "FIRSTY": "10"}
    hdr.update(overrides)
    lines = [f"##{k}={v}" for k, v in hdr.items() if v is not None]
    return "\n".join(lines + ["##XYDATA=(X++(Y..Y))", *table_lines, "##END="]) + "\n"


# ------------------------------------------------------------------ ASDF decoding
@pytest.mark.parametrize("line, expected", [
    ("100 10 20 30", [10, 20, 30]),                  # AFFN, spaces
    ("100,10,20,30", [10, 20, 30]),                  # AFFN, commas
    ("100 5-3+2", [5, -3, 2]),                       # AFFN, signs as separators
    ("100 1.5E+01 2.5", [15.0, 2.5]),                # AFFN exponent (explicit sign)
    ("100@A1B2", [0, 11, 22]),                       # SQZ
    ("100a1b2", [-11, -22]),                         # SQZ negatives
    ("100A0JJ%", [10, 11, 12, 12]),                  # DIF (incl. zero difference)
    ("100A0JT", [10, 11, 12]),                       # DUP repeats the DIF (T = 2 in total)
    ("100A0U", [10, 10, 10]),                        # DUP repeats an absolute value
    ("100A0Js", [10] + list(range(11, 20))),         # DUP 9 ('s')
    ("100A0JS2", [10] + list(range(11, 23))),        # multi-digit DUP (12)
    ("100 5 ? 7", [5, math.nan, 7]),                 # missing value
])
def test_decode_line(line, expected):
    x, values, _, _ = jdx.decode_line(line)
    assert x == 100
    np.testing.assert_allclose(values, expected, equal_nan=True)


def test_decode_reports_dif_state():
    assert jdx.decode_line("100A0JJ")[3] is True
    assert jdx.decode_line("100A0JJA2")[3] is False


def test_decode_rejects_bad_characters():
    with pytest.raises(jdx.JcampError, match="invalid character"):
        jdx.decode_line("100 A0 #5")


def test_label_normalisation():
    assert jdx.normalize_label("XY DATA") == jdx.normalize_label("xy_data") == "XYDATA"
    assert jdx.normalize_label("JCAMP-DX") == "JCAMPDX"


# ------------------------------------------------------------------ parsing / validation
def test_minimal_file_is_valid():
    sp = jdx.parse(minimal(["100 10 20 30"]))
    assert sp.ok, sp.issues
    np.testing.assert_allclose(sp.x, [100, 101, 102])
    np.testing.assert_allclose(sp.y, [10, 20, 30])


def test_factors_are_applied():
    sp = jdx.parse(minimal(["50 10 20 30"], XFACTOR="2", YFACTOR="0.5", FIRSTY="5"))
    assert sp.ok, sp.issues
    np.testing.assert_allclose(sp.y, [5, 10, 15])


def test_comments_and_spaced_labels():
    text = minimal(["100 10 20 30 $$ trailing comment"]).replace("##DATA TYPE", "##Data_Type")
    sp = jdx.parse(text)
    assert sp.ok, sp.issues
    assert sp.header("data type") == "INFRARED SPECTRUM"


def test_dif_lines_with_y_check():
    # second line repeats the last ordinate (12) as its Y-check; x of that line is point 3
    sp = jdx.parse(minimal(["100A0JJ", "102A2"]))
    assert sp.ok, sp.issues
    np.testing.assert_allclose(sp.y, [10, 11, 12])


def test_y_check_failure_is_an_error():
    sp = jdx.parse(minimal(["100A0JJ", "102A5"]))
    assert any("Y-check failed" in i.message for i in sp.errors)


def test_x_check_failure_is_an_error():
    sp = jdx.parse(minimal(["100 10", "105 20 30"]))
    assert any("X-check failed" in i.message for i in sp.errors)


def test_npoints_mismatch_is_an_error():
    sp = jdx.parse(minimal(["100 10 20 30"], NPOINTS="4", LASTX="103"))
    assert any("NPOINTS=4" in i.message for i in sp.errors)


def test_firsty_mismatch_is_an_error():
    sp = jdx.parse(minimal(["100 10 20 30"], FIRSTY="99"))
    assert any("FIRSTY" in i.message for i in sp.errors)


def test_missing_required_label_is_an_error():
    sp = jdx.parse(minimal(["100 10 20 30"], YUNITS=None))
    assert any("##YUNITS" in i.message for i in sp.errors)


def test_missing_end_and_owner_are_warnings():
    text = minimal(["100 10 20 30"], OWNER=None).replace("##END=\n", "")
    sp = jdx.parse(text)
    assert sp.ok
    messages = [i.message for i in sp.issues if i.level == "warning"]
    assert any("##END" in m for m in messages) and any("##OWNER" in m for m in messages)


def test_strict_mode_raises():
    with pytest.raises(jdx.JcampError, match="X-check"):
        jdx.parse(minimal(["100 10", "105 20 30"]), strict=True)


def test_xypoints_table():
    text = minimal([], NPOINTS="3", FIRSTX=None, LASTX=None, FIRSTY=None).replace(
        "##XYDATA=(X++(Y..Y))", "##XYPOINTS=(XY..XY)\n1.0, 5; 2.5, 6\n4.0, 7")
    sp = jdx.parse(text)
    assert sp.ok, sp.issues
    np.testing.assert_allclose(sp.x, [1.0, 2.5, 4.0])
    np.testing.assert_allclose(sp.y, [5, 6, 7])


@pytest.mark.parametrize("label, message", [("##BLOCKS=2", "BLOCKS"), ("##NTUPLES=NMR SPECTRUM", "NTUPLES")])
def test_unsupported_structures(label, message):
    with pytest.raises(jdx.JcampError, match=message):
        jdx.parse(label + "\n" + minimal(["100 10 20 30"]))


def test_unsupported_variable_list():
    with pytest.raises(jdx.JcampError, match="variable list"):
        jdx.parse(minimal(["100 10 20 30"]).replace("(X++(Y..Y))", "(R++(I..I))"))


# ------------------------------------------------------------------ writer round trips
@pytest.mark.parametrize("form", ["AFFN", "DIFDUP"])
@pytest.mark.parametrize("decreasing", [False, True])
def test_round_trip(form, decreasing):
    rng = np.random.default_rng(1)
    x = np.linspace(400, 4000, 901)
    if decreasing:
        x = x[::-1]
    y = np.exp(-0.5 * ((x - 1700) / 40) ** 2) + rng.normal(0, 0.01, x.size)
    y[100:130] = y[100]                       # a flat run exercises DUP
    text = jdx.write(x, y, title="round trip", form=form)
    assert max(len(l) for l in text.splitlines()) <= 80
    sp = jdx.parse(text, strict=True)
    yfactor = float(sp.header("YFACTOR"))
    np.testing.assert_allclose(sp.x, x, rtol=0, atol=1e-6)
    assert np.max(np.abs(sp.y - y)) <= yfactor / 2 + 1e-12


def test_difdup_is_smaller_than_affn():
    x = np.linspace(0, 100, 1001)
    y = np.sin(x / 5)
    assert len(jdx.write(x, y, title="c", form="DIFDUP")) < len(jdx.write(x, y, title="c", form="AFFN"))


def test_writer_rejects_uneven_x():
    with pytest.raises(ValueError, match="evenly spaced"):
        jdx.write(np.array([0, 1, 3.0]), np.array([1, 2, 3.0]), title="bad")


# ------------------------------------------------------------------ sample data + CLI
@pytest.mark.parametrize("name", ["synthetic_ir_absorbance.jdx", "synthetic_uv_vis.jdx", "synthetic_ir_peak_table.jdx"])
def test_sample_files_are_valid(name):
    sp = jdx.read(SAMPLES / name, strict=True)
    assert len(sp.x) == int(sp.header("NPOINTS"))


def test_invalid_example_reports_all_defects():
    sp = jdx.read(SAMPLES / "invalid_example.jdx")
    text = " ".join(i.message for i in sp.issues)
    assert not sp.ok
    for expected in ("Y-check failed", "X-check failed", "NPOINTS=1800", "##OWNER"):
        assert expected in text, expected
    # one wrong NPOINTS must not cascade into an X-check error on every line
    assert sum("X-check failed" in i.message for i in sp.errors) == 1
    np.testing.assert_allclose(sp.x[[0, -1]], [4000, 400])


def test_cli_exit_codes(tmp_path, capsys):
    good = SAMPLES / "synthetic_uv_vis.jdx"
    assert jdx.main([str(good), "--csv", str(tmp_path / "uv.csv"), "--plot", str(tmp_path / "uv.png")]) == 0
    assert (tmp_path / "uv.csv").read_text().count("\n") == 602       # header + 601 points
    assert (tmp_path / "uv.png").stat().st_size > 1000
    assert jdx.main([str(SAMPLES / "invalid_example.jdx")]) == 1
    assert "validation : " in capsys.readouterr().out
