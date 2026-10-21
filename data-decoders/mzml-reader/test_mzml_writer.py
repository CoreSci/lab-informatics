"""Unit tests for mzml_writer.py.

Run from this folder: ``python -m pytest -q``
"""
import hashlib
import io
import re
from array import array
from pathlib import Path

import pytest

import mzml
import mzml_chromatograms as chrom
import mzml_validator as mzv
import mzml_writer as writer

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "sample-data"
SAMPLE_NAMES = ["synthetic_caffeine_lc_ms_centroid.mzML", "synthetic_caffeine_profile.mzML",
                "synthetic_caffeine_ms2_indexed.mzML"]
MZ = [100.0, 150.123456789, 195.08765, 399.999999]
INTENSITY = [0.0, 12.5, 1.0e6 / 3, 3.25]


def spectrum(sid="s1", ms_level=1, mz=MZ, intensity=INTENSITY, **kw):
    fields = dict(id=sid, index=0, ms_level=ms_level, centroid=True, scan_start_time=1.5,
                  scan_start_time_unit="minute", mz=array("d", mz), intensity=array("d", intensity),
                  total_ion_current=sum(intensity), intensity_unit="number of detector counts")
    fields.update(kw)
    return mzml.Spectrum(**fields)


def ms2(sid="s2", ref="s1"):
    pre = mzml.Precursor(spectrum_ref=ref, selected_ion_mz=195.08765, charge=1, intensity=1234.5,
                         isolation_target_mz=195.08765, isolation_lower_offset=0.5, isolation_upper_offset=0.5,
                         activation=["collision-induced dissociation"], collision_energy=20.0)
    return spectrum(sid, ms_level=2, mz=[56.05, 138.07], intensity=[5.0, 20.0], total_ion_current=25.0,
                    precursors=[pre], scan_start_time=1.51, index=1)


def chromatogram(cid="TIC", time=(1.5, 1.51), values=(1.0, 2.0)):
    return mzml.Chromatogram(id=cid, index=0, kind="total ion current chromatogram", time=array("d", time),
                             time_unit="minute", intensity=array("d", values), intensity_unit="number of detector counts")


def write_and_read(tmp_path, spectra, chromatograms=(), **options):
    path = tmp_path / "w.mzML"
    writer.write(path, spectra, chromatograms, **options)
    return path, mzml.read(path)


def as_float32(values):
    return list(array("d", array("f", values)))


# ------------------------------------------------------------------ round trips
@pytest.mark.parametrize("bits", [32, 64])
@pytest.mark.parametrize("zlib_compressed", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
def test_round_trip_every_encoding(tmp_path, bits, zlib_compressed, indexed):
    spectra, chroms = [spectrum(), ms2()], [chromatogram()]
    path, run = write_and_read(tmp_path, spectra, chroms, bits=bits, zlib_compressed=zlib_compressed, indexed=indexed)
    assert run.indexed is indexed
    expect = (lambda v: list(v)) if bits == 64 else as_float32
    for got, sent in zip(run.spectra, spectra):
        assert list(got.mz) == expect(sent.mz) and list(got.intensity) == expect(sent.intensity)
        assert (got.id, got.ms_level, got.centroid, got.scan_start_time, got.scan_start_time_unit) == \
            (sent.id, sent.ms_level, sent.centroid, sent.scan_start_time, sent.scan_start_time_unit)
        assert got.total_ion_current == sent.total_ion_current and got.precursors == sent.precursors
    assert list(run.chromatograms[0].time) == expect(chroms[0].time)
    if bits == 64:
        assert run.spectra == spectra and run.chromatograms == chroms
    assert mzv.validate(path, strict=True).findings == []


@pytest.mark.parametrize("name", SAMPLE_NAMES)
@pytest.mark.parametrize("indexed", [False, True])
def test_samples_round_trip_and_validate(tmp_path, name, indexed):
    source = mzml.read(SAMPLES / name)
    path, back = write_and_read(tmp_path, source.spectra, source.chromatograms, indexed=indexed)
    assert back.spectra == source.spectra and back.chromatograms == source.chromatograms
    assert mzv.validate(path, strict=True).findings == []


def test_profile_and_no_time():
    data = writer.to_bytes([spectrum(centroid=False, scan_start_time=None, scan_start_time_unit=None)])
    (sp,) = mzml.read(io.BytesIO(data)).spectra
    assert sp.centroid is False and sp.scan_start_time is None


def test_empty_spectrum_and_empty_file(tmp_path):
    path, run = write_and_read(tmp_path, [spectrum(mz=[], intensity=[], total_ion_current=0.0)], indexed=True)
    assert len(run.spectra[0].mz) == 0 and mzv.validate(path, strict=True).findings == []
    path, run = write_and_read(tmp_path, [], indexed=True)
    assert run.spectra == [] and mzv.validate(path, strict=True).findings == []


def test_traces_can_be_written_as_chromatograms(tmp_path):
    result = chrom.extract(SAMPLES / SAMPLE_NAMES[0], [(195.0877, 10, "ppm")])
    chroms = [result.tic.to_chromatogram("TIC", 0), result.bpc.to_chromatogram("BPC", 1),
              result.xics[0].to_chromatogram("XIC", 2)]
    path, run = write_and_read(tmp_path, [], chroms, indexed=True)
    assert [c.kind for c in run.chromatograms] == ["total ion current chromatogram", "basepeak chromatogram",
                                                   "selected ion current chromatogram"]
    assert list(run.chromatograms[2].intensity) == list(result.xics[0].intensity)
    assert mzv.validate(path, strict=True).findings == []


def test_ids_are_escaped(tmp_path):
    path, run = write_and_read(tmp_path, [spectrum('a"<&>b')], indexed=True)
    assert run.spectra[0].id == 'a"<&>b' and mzv.validate(path, strict=True).findings == []


def test_writes_to_a_file_object():
    buf = io.BytesIO()
    writer.write(buf, [spectrum()])
    assert buf.getvalue() == writer.to_bytes([spectrum()])


def test_output_is_deterministic_and_lf_only():
    a = writer.to_bytes([spectrum(), ms2()], [chromatogram()], indexed=True, zlib_compressed=False)
    assert a == writer.to_bytes([spectrum(), ms2()], [chromatogram()], indexed=True, zlib_compressed=False)
    assert b"\r" not in a


def test_metadata_is_generic():
    text = writer.to_bytes([spectrum()]).decode("utf-8")
    assert 'name="instrument model" value=""' in text and 'value="mzml_writer.py"' in text
    assert "<componentList" not in text


# ------------------------------------------------------------------ indexed output
def test_indexed_offsets_and_checksum():
    data = writer.to_bytes([spectrum(), ms2()], [chromatogram()], indexed=True)
    offsets = re.findall(rb'<offset idRef="([^"]+)">(\d+)</offset>', data)
    assert [ref for ref, _ in offsets] == [b"s1", b"s2", b"TIC"]
    for ref, pos in offsets:
        assert re.match(rb'<(spectrum|chromatogram) index="\d+" id="' + ref + b'"', data[int(pos):])
    index_list = int(re.search(rb"<indexListOffset>(\d+)<", data).group(1))
    assert data[index_list:].startswith(b'<indexList count="2">')
    end = data.index(b"<fileChecksum>") + len(b"<fileChecksum>")
    assert re.search(rb"<fileChecksum>(\w+)<", data).group(1).decode() == hashlib.sha1(data[:end]).hexdigest()


# ------------------------------------------------------------------ input that would make an invalid file
@pytest.mark.parametrize("bad, message", [
    (dict(ms_level=None), "ms_level"),
    (dict(centroid=None), "centroid"),
    (dict(ms_level=2), "needs at least one precursor"),
    (dict(intensity=[1.0, 2.0]), "4 m/z values but 2 intensities"),
    (dict(mz=[100.0, 99.0, 101.0, 102.0]), "m/z array decreases"),
    (dict(intensity=[1.0, -2.0, 3.0, 4.0]), "negative intensities"),
    (dict(mz=[100.0, float("nan"), 101.0, 102.0]), "NaN"),
    (dict(scan_start_time_unit=None), "needs its unit"),
    (dict(scan_start_time_unit="hour"), "time unit 'hour'"),
    (dict(intensity_unit="volts"), "intensity unit 'volts'"),
])
def test_invalid_spectrum_is_refused(bad, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        writer.to_bytes([spectrum(**bad)])


def test_invalid_input_is_refused():
    with pytest.raises(ValueError, match="duplicate spectrum id"):
        writer.to_bytes([spectrum("x"), spectrum("x")])
    with pytest.raises(ValueError, match="non-empty string id"):
        writer.to_bytes([spectrum("")])
    with pytest.raises(ValueError, match="time array decreases"):
        writer.to_bytes([], [chromatogram(time=(2.0, 1.0))])
    with pytest.raises(ValueError, match="chromatogram type"):
        writer.to_bytes([], [mzml.Chromatogram("c", 0, "made-up chromatogram", array("d"), None, array("d"))])
    with pytest.raises(ValueError, match="activation"):
        sp = ms2()
        sp.precursors[0].activation = ["made-up activation"]
        writer.to_bytes([spectrum(), sp])
    with pytest.raises(ValueError, match="bits"):
        writer.to_bytes([spectrum()], bits=16)
