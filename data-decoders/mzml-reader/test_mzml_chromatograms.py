"""Unit tests for mzml_chromatograms.py and generate_example_plots.py.

Run from this folder: ``python -m pytest -q``
"""
import csv
import io
import math
import xml.etree.ElementTree as ET
from array import array
from pathlib import Path

import pytest

import generate_example_plots as plots
import mzml
import mzml_chromatograms as chrom
import mzml_writer as writer

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "sample-data"
CENTROID = SAMPLES / "synthetic_caffeine_lc_ms_centroid.mzML"
PROFILE = SAMPLES / "synthetic_caffeine_profile.mzML"
INDEXED = SAMPLES / "synthetic_caffeine_ms2_indexed.mzML"
CAFFEINE_MH = 195.0877


@pytest.fixture(scope="module")
def run():
    return mzml.read(CENTROID)


@pytest.fixture(scope="module")
def result():
    return chrom.extract(CENTROID, [(CAFFEINE_MH, 10, "ppm"), (500.0, 0.01, "Da")])


def doc(spectra_args, chromatograms=(), unit="minute") -> io.BytesIO:
    """An in-memory mzML file from (index, time, mz, intensity, ms_level) tuples, written by mzml_writer."""
    spectra = []
    for i, t, mz, inten, level in spectra_args:
        pre = [mzml.Precursor(selected_ion_mz=CAFFEINE_MH, activation=["collision-induced dissociation"])] if level > 1 else []
        spectra.append(mzml.Spectrum(id=f"scan={i + 1}", index=i, ms_level=level, centroid=True, scan_start_time=t,
                                     scan_start_time_unit=unit if isinstance(unit, str) else unit[i],
                                     mz=array("d", mz), intensity=array("d", inten), total_ion_current=sum(inten),
                                     precursors=pre, intensity_unit="number of detector counts"))
    return io.BytesIO(writer.to_bytes(spectra, chromatograms, zlib_compressed=False))


def tic_chromatogram(time, values, unit="minute"):
    return mzml.Chromatogram(id="TIC", index=0, kind="total ion current chromatogram", time=array("d", time),
                             time_unit=unit, intensity=array("d", values))


# ------------------------------------------------------------------ windows
def test_mz_window_ppm_and_da():
    lo, hi = chrom.mz_window(200.0, 10, "ppm")
    assert (lo, hi) == pytest.approx((199.998, 200.002))
    assert chrom.mz_window(200.0, 0.002, "Da") == pytest.approx((lo, hi))
    assert chrom.mz_window(200.0, 0.002, "da") == chrom.mz_window(200.0, 0.002, "Da")


@pytest.mark.parametrize("tol, unit", [(-1, "ppm"), (5, "mmu")])
def test_mz_window_rejects_bad_input(tol, unit):
    with pytest.raises(ValueError):
        chrom.mz_window(200.0, tol, unit)


# ------------------------------------------------------------------ TIC / BPC / XIC on the synthetic run
def test_tic_is_the_sum_of_each_ms1_spectrum(run, result):
    assert len(result.tic.time) == 60 and result.ms1_spectra == 60 and result.skipped_spectra == 0
    assert list(result.tic.time) == [s.scan_start_time for s in run.spectra]
    assert list(result.tic.intensity) == [math.fsum(s.intensity) for s in run.spectra]
    assert result.tic.spectrum_ids == [s.id for s in run.spectra]
    assert result.tic.time_unit == "minute" and result.tic.intensity_unit == "number of detector counts"


def test_tic_peaks_at_the_elution_apex(result):
    t, _ = result.tic.apex()
    assert t == pytest.approx(2.80, abs=0.02)


def test_bpc_is_the_largest_peak_with_its_mz(run, result):
    for s, value, mz in zip(run.spectra, result.bpc.intensity, result.bpc.mz):
        i = max(range(len(s.intensity)), key=s.intensity.__getitem__)
        assert (value, mz) == (s.intensity[i], s.mz[i])
    t, _ = result.bpc.apex()
    i = list(result.bpc.time).index(t)
    assert t == pytest.approx(2.80, abs=0.02)
    assert result.bpc.mz[i] == pytest.approx(CAFFEINE_MH, abs=0.002)


def test_xic_of_caffeine_peaks_at_2_80(result):
    x = result.xics[0]
    assert x.kind == "XIC" and x.window == pytest.approx(chrom.mz_window(CAFFEINE_MH, 10, "ppm"))
    t, peak = x.apex()
    assert t == pytest.approx(2.80, abs=0.02)
    assert peak > 1e6
    # the elution profile is a Gaussian (sigma 0.10 min): far from the apex the XIC falls to background
    assert x.intensity[0] < peak * 1e-3 and x.intensity[-1] < peak * 1e-3


def test_xic_away_from_the_compound_is_flat(result):
    assert set(result.xics[1].intensity) == {0.0}
    assert len(result.xics[1].intensity) == 60


def test_ppm_and_da_tolerances_agree_at_the_same_width(run):
    da = CAFFEINE_MH * 10 / 1e6
    by_ppm = chrom.extract(CENTROID, [(CAFFEINE_MH, 10, "ppm"), (CAFFEINE_MH, da, "Da")])
    assert list(by_ppm.xics[0].intensity) == list(by_ppm.xics[1].intensity)
    assert by_ppm.xics[0].window == pytest.approx(by_ppm.xics[1].window)


def test_xic_sums_every_peak_inside_the_window(run):
    lo, hi = chrom.mz_window(196.0, 1.5, "Da")             # holds M, M+1 and M+2 of caffeine
    x = chrom.xic(CENTROID, 196.0, 1.5, "Da")
    expected = [math.fsum(v for m, v in zip(s.mz, s.intensity) if lo <= m <= hi) for s in run.spectra]
    assert list(x.intensity) == expected


def test_convenience_functions_match_extract(result):
    assert chrom.tic(CENTROID) == result.tic
    assert chrom.bpc(CENTROID) == result.bpc
    assert chrom.xic(CENTROID, CAFFEINE_MH) == result.xics[0]


# ------------------------------------------------------------------ edge cases
def test_empty_window_gives_zeros():
    x = chrom.xic(CENTROID, CAFFEINE_MH + 0.5, 0.0, "Da")     # zero width, between peaks
    assert len(x.intensity) == 60 and set(x.intensity) == {0.0}


def test_ms2_spectra_are_excluded():
    result = chrom.extract(INDEXED, [(138.0662, 20, "ppm")])     # a product ion: only in MS2 spectra
    assert result.ms1_spectra == 3 and result.skipped_spectra == 3
    assert result.tic.spectrum_ids == ["scan=1", "scan=3", "scan=5"]
    assert set(result.xics[0].intensity) == {0.0}


def test_units_are_kept_as_stored():
    t = chrom.tic(doc([(0, 30.0, [100.0], [5.0], 1), (1, 31.5, [100.0], [7.0], 1)], unit="second"))
    assert t.time_unit == "second" and list(t.time) == [30.0, 31.5] and list(t.intensity) == [5.0, 7.0]


def test_mixed_time_units_are_an_error():
    source = doc([(0, 1.0, [100.0], [5.0], 1), (1, 70.0, [100.0], [7.0], 1)], unit=["minute", "second"])
    with pytest.raises(mzml.MzMLError, match="mixed units"):
        chrom.tic(source)


def test_empty_spectrum_and_no_ms1():
    result = chrom.extract(doc([(0, 1.0, [], [], 1)]))
    assert list(result.tic.intensity) == [0.0] and math.isnan(result.bpc.mz[0])
    only_ms2 = chrom.extract(doc([(0, 1.0, [100.0], [5.0], 2)]))
    assert only_ms2.ms1_spectra == 0 and len(only_ms2.tic.time) == 0 and only_ms2.tic.apex() is None


# ------------------------------------------------------------------ TIC cross-check
def test_tic_check_agrees_on_the_centroid_sample(result):
    sources = [c.source for c in result.checks]
    assert sources == ["declared TIC chromatogram 'TIC'", "spectrum total ion current"]
    assert all(c.ok and c.compared == 60 for c in result.checks)


def test_tic_check_reports_a_declared_tic_that_includes_ms2():
    declared, per_spectrum = chrom.extract(INDEXED).checks
    assert not declared.ok and declared.compared == 3 and declared.mismatches == 0
    assert "3 declared point(s) have no MS1 spectrum" in declared.notes[0]
    assert per_spectrum.ok and per_spectrum.compared == 6        # MS2 spectra declare their own TIC


def test_tic_check_reports_differences():
    tic = tic_chromatogram([1.0, 2.0], [5.0, 99.0])
    result = chrom.extract(doc([(0, 1.0, [100.0], [5.0], 1), (1, 2.0, [100.0], [7.0], 1)], [tic]))
    declared = result.checks[0]
    assert (declared.compared, declared.mismatches) == (2, 1)
    assert declared.max_rel_diff == pytest.approx(92 / 99)
    per_spectrum = result.checks[1]                          # doc() declares the true per-spectrum TIC
    assert per_spectrum.ok


def test_tic_check_wrong_spectrum_value():
    xml = doc([(0, 1.0, [100.0], [5.0], 1)]).getvalue().replace(b'name="total ion current" value="5.0"',
                                                                b'name="total ion current" value="6.0"')
    (check,) = chrom.extract(io.BytesIO(xml)).checks
    assert check.mismatches == 1 and not check.ok


def test_tic_check_does_not_convert_units():
    tic = tic_chromatogram([60.0], [5.0], unit="second")
    declared = chrom.extract(doc([(0, 1.0, [100.0], [5.0], 1)], [tic])).checks[0]
    assert declared.compared == 0 and "time units differ" in declared.notes[0]


# ------------------------------------------------------------------ CSV
def test_csv_output(result, tmp_path):
    path = tmp_path / "out.csv"
    chrom.write_csv(result, path)
    rows = list(csv.reader(io.StringIO(path.read_text(encoding="utf-8"))))
    assert rows[0] == ["spectrum id", "time (minute)", "TIC", "BPC", "BPC m/z",
                       "XIC m/z 195.0877 +/- 10 ppm", "XIC m/z 500.0000 +/- 0.01 Da"]
    assert len(rows) == 61
    first = rows[1]
    assert first[0] == "scan=1" and float(first[1]) == result.tic.time[0]
    assert float(first[2]) == result.tic.intensity[0] and float(first[4]) == result.bpc.mz[0]
    assert b"\r\n" not in path.read_bytes()


def test_csv_leaves_missing_bpc_mz_empty():
    text = chrom.csv_text(chrom.extract(doc([(0, 1.0, [], [], 1)])))
    assert text.splitlines()[1].split(",")[4] == ""


# ------------------------------------------------------------------ SVG
def test_svg_is_well_formed_and_labelled(result):
    svg = chrom.svg_plot(result.xics[0])
    root = ET.fromstring(svg)
    texts = " ".join("".join(t.itertext()) for t in root.iter("{http://www.w3.org/2000/svg}text"))
    assert "Time (minute)" in texts and "Intensity, number of detector counts" in texts
    assert "± 10 ppm" in texts
    (line,) = root.iter("{http://www.w3.org/2000/svg}polyline")
    assert len(line.get("points").split()) == 60


def test_svg_is_deterministic(result):
    assert chrom.svg_plot(result.tic) == chrom.svg_plot(chrom.extract(CENTROID).tic)


def test_svg_of_an_empty_trace():
    svg = chrom.svg_plot(chrom.extract(doc([(0, 1.0, [100.0], [5.0], 2)])).tic)
    root = ET.fromstring(svg)
    assert not list(root.iter("{http://www.w3.org/2000/svg}polyline"))


def test_example_plots_are_reproducible(tmp_path):
    for fresh in plots.generate(tmp_path):
        shipped = (SAMPLES / "output" / fresh.name).read_bytes()
        assert fresh.read_bytes() == shipped, fresh.name
        assert len(shipped) < 10_000


def test_ticks_are_round_numbers():
    assert chrom._ticks(0, 2.03) == [0, 0.5, 1, 1.5, 2]
    assert chrom._ticks(1.6, 3.96) == [2, 2.5, 3, 3.5]
    assert chrom._ticks(0, 0) == [0, 0.2, 0.4, 0.6, 0.8, 1]


# ------------------------------------------------------------------ trace -> chromatogram
def test_trace_to_chromatogram(result):
    c = result.bpc.to_chromatogram()
    assert (c.id, c.kind, c.time_unit) == ("BPC", "basepeak chromatogram", "minute")
    assert result.xics[0].to_chromatogram("XIC195", 2).kind == "selected ion current chromatogram"


# ------------------------------------------------------------------ CLI
def test_cli(tmp_path, capsys):
    csv_path, svg_path = tmp_path / "c.csv", tmp_path / "x.svg"
    assert chrom.main([str(CENTROID), "--xic", "195.0877,10,ppm", "--csv", str(csv_path), "--svg", str(svg_path)]) == 0
    out = capsys.readouterr().out
    assert "MS1 spectra : 60" in out and "XIC m/z 195.0877 +/- 10 ppm: max" in out and "TIC check   : OK" in out
    assert csv_path.exists() and "XIC m/z 195.0877" in svg_path.read_text(encoding="utf-8")
    assert chrom.main([str(INDEXED)]) == 0
    assert "TIC check   : DIFFERENCE" in capsys.readouterr().out
    assert chrom.main([str(tmp_path / "missing.mzML")]) == 1


def test_cli_rejects_bad_xic():
    with pytest.raises(SystemExit):
        chrom.main([str(CENTROID), "--xic", "195.0877,10,mmu"])
    with pytest.raises(SystemExit):
        chrom.main([str(CENTROID), "--plot", "xic", "--svg", "x.svg"])
