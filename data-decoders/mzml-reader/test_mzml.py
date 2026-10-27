"""Unit tests for mzml.py and generate_synthetic_mzml.py.

Run from this folder: ``python -m pytest -q``
"""
import io
import random
import tracemalloc
from array import array
from pathlib import Path

import pytest

import generate_synthetic_mzml as gen
import mzml
import psi_ms_cv as cv

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "sample-data"
CENTROID = SAMPLES / "synthetic_caffeine_lc_ms_centroid.mzML"
PROFILE = SAMPLES / "synthetic_caffeine_profile.mzML"
INDEXED = SAMPLES / "synthetic_caffeine_ms2_indexed.mzML"

MZ = [100.0, 150.25, 195.08765, 399.999999]
INTENSITY = [0.0, 12.5, 1.0e6, 3.25]


# ------------------------------------------------------------------ in-memory documents
def param(term, unit=None):
    return gen.cv_param(term, unit=unit)


def bda(values, bits=64, compression=cv.ZLIB_COMPRESSION, array_type=cv.MZ_ARRAY, *,
        text=None, array_length=None, encoded_length=None, params=None):
    """One <binaryDataArray>; ``params`` replaces the default cvParams when given."""
    if text is None:
        text = mzml.encode_binary(values, bits, compression == cv.ZLIB_COMPRESSION)
    if params is None:
        params = [param(cv.FLOAT_32 if bits == 32 else cv.FLOAT_64), param(compression), param(array_type)]
    attrs = f' encodedLength="{len(text) if encoded_length is None else encoded_length}"'
    if array_length is not None:
        attrs += f' arrayLength="{array_length}"'
    return f"<binaryDataArray{attrs}>{''.join(params)}<binary>{text}</binary></binaryDataArray>"


def spectrum(arrays, n=len(MZ), index=0, params=""):
    return (f'<spectrum index="{index}" id="scan={index + 1}" defaultArrayLength="{n}">{params}'
            f'<binaryDataArrayList count="{len(arrays)}">{"".join(arrays)}</binaryDataArrayList></spectrum>')


def document(*spectra, chromatograms="", groups=""):
    chrom = f'<chromatogramList count="1" defaultDataProcessingRef="DP1">{chromatograms}</chromatogramList>' \
        if chromatograms else ""
    return io.BytesIO((
        '<?xml version="1.0" encoding="utf-8"?>'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">'
        f'{groups}<run id="r"><spectrumList count="{len(spectra)}" defaultDataProcessingRef="DP1">'
        f'{"".join(spectra)}</spectrumList>{chrom}</run></mzML>').encode("ascii"))


def one(doc):
    return next(mzml.iter_spectra(doc))


def as_float32(values):
    return list(array("d", array("f", values)))


# ------------------------------------------------------------------ binary decoding round trips
@pytest.mark.parametrize("bits", [32, 64])
@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.req("RD-03")
def test_binary_round_trip(bits, compressed):
    text = mzml.encode_binary(MZ, bits, compressed)
    decoded = mzml.decode_binary(text, bits, compressed, expected_length=len(MZ))
    assert list(decoded) == (MZ if bits == 64 else as_float32(MZ))


@pytest.mark.parametrize("bits", [32, 64])
@pytest.mark.parametrize("compression", [cv.NO_COMPRESSION, cv.ZLIB_COMPRESSION])
@pytest.mark.req("RD-03")
def test_spectrum_round_trip(bits, compression):
    sp = one(document(spectrum([bda(MZ, bits, compression, cv.MZ_ARRAY),
                                bda(INTENSITY, bits, compression, cv.INTENSITY_ARRAY)])))
    expect = (lambda v: v) if bits == 64 else as_float32
    assert list(sp.mz) == expect(MZ)
    assert list(sp.intensity) == expect(INTENSITY)


@pytest.mark.req("RD-03")
def test_mixed_precision_arrays():
    sp = one(document(spectrum([bda(MZ, 64, array_type=cv.MZ_ARRAY), bda(INTENSITY, 32, array_type=cv.INTENSITY_ARRAY)])))
    assert list(sp.mz) == MZ and list(sp.intensity) == as_float32(INTENSITY)


@pytest.mark.req("RD-03")
def test_base64_whitespace_is_ignored():
    text = mzml.encode_binary(MZ)
    wrapped = "\n".join(text[i:i + 16] for i in range(0, len(text), 16))
    assert list(mzml.decode_binary(wrapped)) == MZ


@pytest.mark.req("RD-03")
def test_empty_spectrum():
    sp = one(document(spectrum([bda([], array_type=cv.MZ_ARRAY), bda([], array_type=cv.INTENSITY_ARRAY)], n=0)))
    assert len(sp.mz) == 0 and len(sp.intensity) == 0


@pytest.mark.req("RD-03")
def test_array_length_overrides_default():
    arrays = [bda(MZ[:2], array_type=cv.MZ_ARRAY, array_length=2),
              bda(INTENSITY[:2], array_type=cv.INTENSITY_ARRAY, array_length=2)]
    assert list(one(document(spectrum(arrays))).mz) == MZ[:2]


# ------------------------------------------------------------------ error paths
def assert_error(doc, *fragments):
    with pytest.raises(mzml.MzMLError) as info:
        list(mzml.iterparse(doc))
    for frag in fragments:
        assert frag in str(info.value), str(info.value)
    return str(info.value)


def with_mz_params(*params):
    return document(spectrum([bda(MZ, params=list(params)), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))


@pytest.mark.parametrize("acc, name", [
    ("MS:1002312", "MS-Numpress linear prediction compression"),
    ("MS:1002746", "MS-Numpress linear prediction compression followed by zlib compression"),
    ("MS:1003785", "MS-Numpress short logged float compression followed by zstd compression"),
])
@pytest.mark.req("RD-08")
def test_numpress_is_out_of_scope(acc, name):
    p = f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value=""/>'
    assert_error(with_mz_params(param(cv.FLOAT_64), p, param(cv.MZ_ARRAY)), "out of scope", "MS-Numpress", acc)


@pytest.mark.parametrize("acc, name", [("MS:1003780", "zstd compression"),
                                       ("MS:9999999", "some future compression")])
@pytest.mark.req("RD-08")
def test_other_compression_is_out_of_scope(acc, name):
    p = f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value=""/>'
    assert_error(with_mz_params(param(cv.FLOAT_64), p, param(cv.MZ_ARRAY)), "out of scope", name)


@pytest.mark.req("RD-08")
def test_missing_compression_term():
    assert_error(with_mz_params(param(cv.FLOAT_64), param(cv.MZ_ARRAY)), "missing compression term")


@pytest.mark.req("RD-08")
def test_missing_array_type():
    msg = assert_error(with_mz_params(param(cv.FLOAT_64), param(cv.ZLIB_COMPRESSION)), "missing array type")
    assert "scan=1" in msg                                  # errors say where they happened


@pytest.mark.req("RD-08")
def test_unknown_array_type():
    p = '<cvParam cvRef="MS" accession="MS:1000516" name="charge array" value=""/>'
    assert_error(with_mz_params(param(cv.FLOAT_64), param(cv.ZLIB_COMPRESSION), p),
                 "unsupported array type", "charge array")


@pytest.mark.req("RD-08")
def test_time_array_not_allowed_in_a_spectrum():
    assert_error(document(spectrum([bda(MZ, array_type=cv.TIME_ARRAY), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)])),
                 "time array", "unexpected")


@pytest.mark.req("RD-08")
def test_duplicate_array():
    assert_error(document(spectrum([bda(MZ), bda(MZ)])), "appears twice")


@pytest.mark.req("RD-08")
def test_missing_intensity_array():
    assert_error(document(spectrum([bda(MZ)])), "no intensity array")


@pytest.mark.req("RD-08")
def test_missing_data_type():
    assert_error(with_mz_params(param(cv.ZLIB_COMPRESSION), param(cv.MZ_ARRAY)), "missing binary data type")


@pytest.mark.req("RD-08")
def test_unsupported_data_type():
    p = '<cvParam cvRef="MS" accession="MS:1000519" name="32-bit integer" value=""/>'
    assert_error(with_mz_params(p, param(cv.ZLIB_COMPRESSION), param(cv.MZ_ARRAY)),
                 "unsupported binary data type", "32-bit integer")


@pytest.mark.req("RD-09")
def test_length_mismatch_default_array_length():
    doc = document(spectrum([bda(MZ), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)], n=5))
    assert_error(doc, "decoded 4 values", "declared array length is 5", "m/z array")


@pytest.mark.req("RD-09")
def test_length_mismatch_array_length_attribute():
    doc = document(spectrum([bda(MZ, array_length=3), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))
    assert_error(doc, "decoded 4 values", "declared array length is 3")


@pytest.mark.req("RD-09")
def test_encoded_length_mismatch():
    doc = document(spectrum([bda(MZ, encoded_length=7), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))
    assert_error(doc, "encodedLength is 7")


@pytest.mark.req("RD-09")
def test_partial_value_bytes():
    text = mzml.encode_binary(MZ, 64, False)
    import base64
    short = base64.b64encode(base64.b64decode(text)[:-3]).decode()
    doc = document(spectrum([bda(MZ, compression=cv.NO_COMPRESSION, text=short),
                             bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))
    assert_error(doc, "not a whole number of 64-bit values")


@pytest.mark.req("RD-09")
def test_invalid_base64():
    doc = document(spectrum([bda(MZ, text="@@@@"), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))
    assert_error(doc, "invalid base64")


@pytest.mark.req("RD-09")
def test_corrupt_zlib():
    text = mzml.encode_binary(MZ, 64, False)                 # declared as zlib, but isn't
    doc = document(spectrum([bda(MZ, compression=cv.ZLIB_COMPRESSION, text=text),
                             bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)]))
    assert_error(doc, "zlib data could not be decompressed")


@pytest.mark.parametrize("xml", [
    b"",                                                     # empty file
    b'<mzML xmlns="http://psi.hupo.org/ms/mzml"><run id="r"><spectrumList count="1">',   # truncated
    b'<mzML xmlns="http://psi.hupo.org/ms/mzml"><run></spectrumList></mzML>',            # mismatched tag
])
@pytest.mark.req("RD-10")
def test_malformed_xml(xml):
    assert_error(io.BytesIO(xml), "malformed XML")


@pytest.mark.req("RD-10")
def test_not_mzml():
    assert_error(io.BytesIO(b"<html><body/></html>"), "not an mzML document", "<html>")


@pytest.mark.req("RD-10")
def test_missing_index_attribute():
    xml = document(spectrum([bda(MZ), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)])).getvalue()
    assert_error(io.BytesIO(xml.replace(b'index="0" ', b"")), "missing required attribute 'index'")


# ------------------------------------------------------------------ metadata
@pytest.mark.req("RD-04", "RD-07")
def test_spectrum_metadata_and_param_groups():
    groups = ('<referenceableParamGroupList count="1"><referenceableParamGroup id="common">'
              f'{param(cv.MS1_SPECTRUM)}{param(cv.PROFILE_SPECTRUM)}'
              '</referenceableParamGroup></referenceableParamGroupList>')
    params = ('<referenceableParamGroupRef ref="common"/>' + gen.cv_param(cv.MS_LEVEL, 1) +
              gen.cv_param(cv.TOTAL_ION_CURRENT, "1012515.75") +
              '<scanList count="1"><scan>' + gen.cv_param(cv.SCAN_START_TIME, "42.5", cv.UO_SECOND) +
              "</scan></scanList>")
    sp = one(document(spectrum([bda(MZ), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)], params=params),
                      groups=groups))
    assert (sp.id, sp.index, sp.ms_level, sp.centroid) == ("scan=1", 0, 1, False)
    assert (sp.scan_start_time, sp.scan_start_time_unit) == (42.5, "second")
    assert sp.total_ion_current == 1012515.75
    assert sp.precursors == []


@pytest.mark.req("RD-07")
def test_unknown_param_group():
    params = '<referenceableParamGroupRef ref="nope"/>'
    assert_error(document(spectrum([bda(MZ), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)], params=params)),
                 "unknown referenceableParamGroup 'nope'")


# ------------------------------------------------------------------ streaming
class CountingReader(io.RawIOBase):
    """A binary file wrapper that counts the bytes handed to the parser."""

    def __init__(self, path):
        self.fh, self.read_bytes = open(path, "rb"), 0

    def readable(self):
        return True

    def readinto(self, buf):
        data = self.fh.read(min(len(buf), 16384))
        buf[:len(data)] = data
        self.read_bytes += len(data)
        return len(data)

    def close(self):
        self.fh.close()
        super().close()


def big_file(path: Path, n: int) -> Path:
    rng = random.Random(1)
    spectra = []
    for i in range(n):
        mz = sorted(round(rng.uniform(100, 1000), 4) for _ in range(200))
        inten = [round(rng.uniform(0, 1e5), 1) for _ in mz]
        spectra.append(gen.spectrum_xml(i, i / 60, mz, inten, ms_level=1, centroid=True, window=(100.0, 1000.0),
                                        mz_bits=64, int_bits=32, zlib_compressed=True))
    path.write_bytes(gen.plain_mzml("big", spectra, []))
    return path


@pytest.mark.req("RD-02")
def test_first_spectrum_arrives_before_the_file_is_read(tmp_path):
    path = big_file(tmp_path / "big.mzML", 400)
    reader = CountingReader(path)
    try:
        first = next(mzml.iter_spectra(reader))
        assert first.index == 0
        assert reader.read_bytes < path.stat().st_size / 10
    finally:
        reader.close()


def peak_memory(path: Path) -> int:
    tracemalloc.start()
    try:
        count = sum(1 for _ in mzml.iter_spectra(path))
        return count, tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


@pytest.mark.req("RD-02")
def test_memory_stays_flat(tmp_path):
    small_n, small_peak = peak_memory(big_file(tmp_path / "small.mzML", 100))
    large_n, large_peak = peak_memory(big_file(tmp_path / "large.mzML", 1500))
    assert (small_n, large_n) == (100, 1500)
    assert large_peak < 2 * small_peak                      # 15x the spectra, not 15x the memory


@pytest.mark.req("RD-02")
def test_items_arrive_before_a_later_error():
    good = spectrum([bda(MZ), bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)])
    xml = document(good).getvalue().replace(b"</spectrumList></run></mzML>", b"<spectrum")
    it = mzml.iter_spectra(io.BytesIO(xml))
    assert next(it).id == "scan=1"
    with pytest.raises(mzml.MzMLError, match="malformed XML"):
        next(it)


# ------------------------------------------------------------------ sample files: plain and indexed
@pytest.mark.req("RD-04", "DG-02")
def test_centroid_sample():
    run = mzml.read(CENTROID)
    assert not run.indexed and run.run_id == "synthetic_caffeine_lc_ms_centroid"
    assert len(run.spectra) == 60
    assert all(s.ms_level == 1 and s.centroid is True and s.scan_start_time_unit == "minute" for s in run.spectra)
    assert [s.index for s in run.spectra] == list(range(60))
    for s in run.spectra:
        assert len(s.mz) == len(s.intensity) > 0 and list(s.mz) == sorted(s.mz)
        assert s.total_ion_current == pytest.approx(sum(s.intensity), abs=1)
    apex = max(run.spectra, key=lambda s: s.total_ion_current)
    assert apex.scan_start_time == pytest.approx(2.80, abs=0.05)
    # isotope ratio at the apex: M+1 / M of [M+H]+ is about 10.3 %
    peaks = dict(zip((round(m, 2) for m in apex.mz), apex.intensity))
    assert 0.09 < peaks[196.09] / peaks[195.09] < 0.115


@pytest.mark.req("RD-03", "RD-04")
def test_profile_sample():
    run = mzml.read(PROFILE)
    assert len(run.spectra) == 5 and not run.chromatograms
    s = run.spectra[0]
    assert s.centroid is False and len(s.mz) == 1201
    assert list(s.mz) == as_float32([round(193.5 + 0.005 * k, 4) for k in range(1201)])   # 32-bit storage
    top = max(range(len(s.mz)), key=s.intensity.__getitem__)
    assert s.mz[top] == pytest.approx(195.0877, abs=0.003)


@pytest.mark.req("RD-01", "RD-05")
def test_indexed_sample_and_precursors():
    run = mzml.read(INDEXED)
    assert run.indexed
    assert [s.ms_level for s in run.spectra] == [1, 2, 1, 2, 1, 2]
    for ms1, ms2 in zip(run.spectra[::2], run.spectra[1::2]):
        assert ms1.precursors == []
        (pre,) = ms2.precursors
        assert pre.spectrum_ref == ms1.id
        assert pre.selected_ion_mz == pytest.approx(195.08765, abs=1e-4)
        assert pre.isolation_target_mz == pre.selected_ion_mz
        assert (pre.isolation_lower_offset, pre.isolation_upper_offset) == (0.5, 0.5)
        assert pre.charge == 1 and pre.collision_energy == 20.0 and pre.intensity > 0
        assert pre.activation == ["collision-induced dissociation"]
        assert any(abs(m - 138.0662) < 0.001 for m in ms2.mz)


@pytest.mark.req("DG-01")
def test_indexed_sample_offsets_and_checksum():
    """The index the generator writes is valid (the reader ignores it in this unit)."""
    import hashlib
    import re
    data = INDEXED.read_bytes()
    offsets = re.findall(rb'<offset idRef="([^"]+)">(\d+)</offset>', data)
    assert len(offsets) == 7                                # 6 spectra + 1 chromatogram
    for ref, pos in offsets:
        assert re.match(rb'<(spectrum|chromatogram) index="\d+" id="' + re.escape(ref) + b'"', data[int(pos):])
    index_list = int(re.search(rb"<indexListOffset>(\d+)<", data).group(1))
    assert data[index_list:].startswith(b"<indexList ")
    end = data.index(b"<fileChecksum>") + len(b"<fileChecksum>")
    assert re.search(rb"<fileChecksum>(\w+)<", data).group(1).decode() == hashlib.sha1(data[:end]).hexdigest()


@pytest.mark.req("RD-01")
def test_plain_and_indexed_wrappers_give_the_same_spectra(tmp_path):
    ions = gen.caffeine_ions()
    spectra = [gen.spectrum_xml(0, 1.0, [m for m, _ in ions], [r for _, r in ions], ms_level=1, centroid=True,
                                window=(100.0, 400.0), mz_bits=64, int_bits=64, zlib_compressed=True)]
    plain = mzml.read(io.BytesIO(gen.plain_mzml("x", spectra, [])))
    indexed = mzml.read(io.BytesIO(gen.indexed_mzml("x", spectra, [])))
    assert (plain.indexed, indexed.indexed) == (False, True)
    assert plain.spectra == indexed.spectra


# ------------------------------------------------------------------ chromatograms
@pytest.mark.req("RD-06")
def test_tic_chromatogram():
    run = mzml.read(CENTROID)
    (tic,) = run.chromatograms
    assert (tic.id, tic.index, tic.kind, tic.time_unit) == ("TIC", 0, "total ion current chromatogram", "minute")
    assert list(tic.time) == [s.scan_start_time for s in run.spectra]
    assert list(tic.intensity) == pytest.approx([sum(s.intensity) for s in run.spectra], abs=1)


@pytest.mark.req("RD-06")
def test_iter_chromatograms_and_iterparse_order():
    assert [c.id for c in mzml.iter_chromatograms(INDEXED)] == ["TIC"]
    kinds = [type(x).__name__ for x in mzml.iterparse(INDEXED)]
    assert kinds == ["Spectrum"] * 6 + ["Chromatogram"]


@pytest.mark.req("RD-06", "RD-08")
def test_chromatogram_needs_a_time_array():
    chrom = ('<chromatogram index="0" id="TIC" defaultArrayLength="4">' + param(cv.TIC_CHROMATOGRAM) +
             f'<binaryDataArrayList count="2">{bda(MZ)}{bda(INTENSITY, array_type=cv.INTENSITY_ARRAY)}'
             "</binaryDataArrayList></chromatogram>")
    assert_error(document(chromatograms=chrom), "chromatogram 'TIC'", "m/z array", "unexpected")


# ------------------------------------------------------------------ generator
@pytest.mark.req("DG-01")
def test_generator_is_deterministic(tmp_path):
    a = gen.generate(tmp_path / "a")
    b = gen.generate(tmp_path / "b")
    for pa, pb in zip(a, b):
        assert pa.read_bytes() == pb.read_bytes(), pa.name


@pytest.mark.req("DG-01")
def test_shipped_samples_match_the_generator(tmp_path):
    for fresh in gen.generate(tmp_path):
        shipped = SAMPLES / fresh.name
        if b"zlib compression" in fresh.read_bytes():
            # zlib output may differ between zlib builds; the decoded content must not
            assert mzml.read(fresh) == mzml.read(shipped), fresh.name
        else:
            assert fresh.read_bytes() == shipped.read_bytes(), fresh.name


@pytest.mark.req("DG-01")
def test_other_seed_gives_other_data(tmp_path):
    a = gen.generate(tmp_path / "a", seed=1)
    b = gen.generate(tmp_path / "b", seed=2)
    assert all(pa.read_bytes() != pb.read_bytes() for pa, pb in zip(a, b))


@pytest.mark.req("DG-03")
def test_samples_are_small():
    for path in (CENTROID, PROFILE, INDEXED):
        assert path.stat().st_size < 500_000


@pytest.mark.req("DG-02")
def test_caffeine_isotope_pattern():
    (m0, r0), (m1, r1), (m2, r2) = gen.caffeine_ions()
    assert m0 == pytest.approx(195.08765, abs=5e-5)          # [M+H]+ of C8H10N4O2
    assert r0 == 1.0
    assert m1 - m0 == pytest.approx(1.0030, abs=0.002)
    assert 0.100 < r1 < 0.106                                  # mostly 8 x 13C
    assert 0.007 < r2 < 0.011


@pytest.mark.req("DG-03")
def test_sample_metadata_is_generic():
    for path in (CENTROID, PROFILE, INDEXED):
        text = path.read_text(encoding="ascii")
        assert f'accession="{cv.INSTRUMENT_MODEL[0]}" name="instrument model" value=""' in text
        assert text.count("<cvParam") == text.count('cvRef="MS"') + text.count('cvRef="UO"')
