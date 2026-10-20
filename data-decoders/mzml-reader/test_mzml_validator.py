"""Unit tests for mzml_validator.py and generate_invalid_mzml.py.

Run from this folder: ``python -m pytest -q``

Each check code has a negative case (the valid base document yields no
findings) and a positive case (one targeted edit yields exactly that code).
"""
import re
from pathlib import Path

import pytest

import generate_invalid_mzml as bad
import generate_synthetic_mzml as gen
import mzml
import mzml_validator as mzv
import psi_ms_cv as cv

HERE = Path(__file__).resolve().parent
SAMPLES = HERE / "sample-data"
INVALID = SAMPLES / "invalid"


# ------------------------------------------------------------------ a small valid base document
def base_spectra():
    ions = gen.caffeine_ions()
    mz, inten = [m for m, _ in ions], [1000.0 * r for _, r in ions]
    ms1 = gen.spectrum_xml(0, 1.0, mz, inten, ms_level=1, centroid=True, window=(100.0, 400.0),
                           mz_bits=64, int_bits=32, zlib_compressed=True)
    ms2 = gen.spectrum_xml(1, 1.01, [56.05, 138.07], [5.0, 20.0], ms_level=2, centroid=True,
                           window=(50.0, 200.0), mz_bits=64, int_bits=64, zlib_compressed=False,
                           precursor={"ref": "scan=1", "mz": mz[0], "intensity": inten[0]})
    tic = gen.tic_chromatogram_xml([1.0, 1.01], [sum(inten), 25.0], 64)
    return [ms1, ms2], [tic]


def base_plain() -> bytes:
    spectra, chroms = base_spectra()
    return gen.plain_mzml("base", spectra, chroms, has_msn=True)


def base_indexed() -> bytes:
    spectra, chroms = base_spectra()
    return gen.indexed_mzml("base", spectra, chroms, has_msn=True)


def check(tmp_path, data: bytes, strict=False, name="t.mzML"):
    path = tmp_path / name
    path.write_bytes(data)
    return mzv.validate(path, strict=strict)


def edit(data: bytes, old: bytes, new: bytes, indexed=False) -> bytes:
    """Replace the first occurrence. In an indexed file, keep the index and checksum valid."""
    assert old in data, old
    data = data.replace(old, new, 1)
    if indexed:
        data = reindex(data)
    return data


def reindex(data: bytes) -> bytes:
    """Recompute offsets, indexListOffset and checksum after an edit that shifts bytes."""
    out = data
    for ref in re.findall(rb'<offset idRef="([^"]+)">', data):
        pos = re.search(rb'<(spectrum|chromatogram) index="\d+" id="' + re.escape(ref) + b'"', out).start()
        out = re.sub(rb'(<offset idRef="' + re.escape(ref) + rb'">)\d+', rb"\g<1>" + str(pos).encode(), out, count=1)
    out = re.sub(rb"<indexListOffset>\d+", b"<indexListOffset>" + str(out.index(b"<indexList ")).encode(), out)
    return bad.fix_checksum(out)


def only(report, code):
    assert report.codes == [code], [str(f) for f in report.findings]
    return report.findings[0]


def array_cv(term):
    return f'accession="{term[0]}" name="{term[1]}"'.encode()


# ------------------------------------------------------------------ negative cases
def test_valid_base_documents_have_no_findings(tmp_path):
    assert check(tmp_path, base_plain()).findings == []
    report = check(tmp_path, base_indexed())
    assert report.findings == [] and report.indexed and (report.spectra, report.chromatograms) == (2, 1)


@pytest.mark.parametrize("name", ["synthetic_caffeine_lc_ms_centroid.mzML", "synthetic_caffeine_profile.mzML",
                                  "synthetic_caffeine_ms2_indexed.mzML"])
def test_valid_samples_have_no_findings(name):
    report = mzv.validate(SAMPLES / name, strict=True)
    assert report.findings == [] and report.ok


# ------------------------------------------------------------------ one positive case per code
def test_mzv001_malformed_xml(tmp_path):
    f = only(check(tmp_path, base_plain()[:-200]), "MZV001")
    assert "malformed XML" in f.message


def test_mzv001_not_mzml(tmp_path):
    assert "<html>" in only(check(tmp_path, b"<html/>"), "MZV001").message


@pytest.mark.parametrize("section", mzv.REQUIRED_SECTIONS[:-1])
def test_mzv002_missing_section(tmp_path, section):
    data = base_plain()
    s = section.encode()
    start = data.index(b"<" + s)
    end = data.index(b"</" + s + b">") + len(s) + 3
    report = check(tmp_path, data[:start] + data[end:])
    assert "MZV002" in report.codes
    assert any(f"<{section}>" in f.message for f in report.findings if f.code == "MZV002")


def test_mzv003_undeclared_cv(tmp_path):
    data = edit(base_plain(), b'<cv id="UO"', b'<cv id="XX"')
    f = only(check(tmp_path, data), "MZV003")
    assert "cvRef 'UO'" in f.message and "spectrum 'scan=1'" in f.message


def test_mzv004_param_group(tmp_path):
    data = edit(base_plain(), b'<spectrum index="1" id="scan=2" defaultArrayLength="2">',
                b'<spectrum index="1" id="scan=2" defaultArrayLength="2"><referenceableParamGroupRef ref="nope"/>')
    f = only(check(tmp_path, data), "MZV004")
    assert "'nope'" in f.message and "spectrum 'scan=2'" in f.message


def test_mzv004_resolved_group_is_fine(tmp_path):
    data = edit(base_plain(), b"<softwareList",
                b'<referenceableParamGroupList count="1"><referenceableParamGroup id="g">'
                + gen.cv_param(cv.POSITIVE_SCAN).encode() + b"</referenceableParamGroup></referenceableParamGroupList><softwareList")
    data = edit(data, b'id="scan=2" defaultArrayLength="2">', b'id="scan=2" defaultArrayLength="2"><referenceableParamGroupRef ref="g"/>')
    assert check(tmp_path, data).findings == []


def test_mzv005_instrument_configuration(tmp_path):
    data = edit(base_plain(), b'defaultInstrumentConfigurationRef="IC1"', b'defaultInstrumentConfigurationRef="IC9"')
    assert "'IC9'" in only(check(tmp_path, data), "MZV005").message
    data = edit(base_plain(), b"<scan>", b'<scan instrumentConfigurationRef="IC2">')
    assert "spectrum 'scan=1'" in only(check(tmp_path, data), "MZV005").message


def test_mzv006_data_processing(tmp_path):
    data = edit(base_plain(), b'<chromatogramList count="1" defaultDataProcessingRef="DP1"',
                b'<chromatogramList count="1" defaultDataProcessingRef="DP7"')
    assert "'DP7'" in only(check(tmp_path, data), "MZV006").message
    data = edit(base_plain(), b'id="scan=2" defaultArrayLength="2"', b'id="scan=2" defaultArrayLength="2" dataProcessingRef="DP3"')
    assert "spectrum 'scan=2'" in only(check(tmp_path, data), "MZV006").message


def test_mzv010_list_counts(tmp_path):
    only(check(tmp_path, edit(base_plain(), b'<spectrumList count="2"', b'<spectrumList count="3"')), "MZV010")
    f = only(check(tmp_path, edit(base_plain(), b'<chromatogramList count="1"', b'<chromatogramList count="0"')), "MZV010")
    assert "chromatogramList" in f.message


def test_mzv011_index_not_consecutive(tmp_path):
    data = edit(base_plain(), b'<spectrum index="1" id="scan=2"', b'<spectrum index="2" id="scan=2"')
    f = only(check(tmp_path, data), "MZV011")
    assert "spectrum 'scan=2'" in f.message and "expected 1" in f.message


def test_mzv012_duplicate_id(tmp_path):
    data = edit(base_plain(), b'<spectrum index="1" id="scan=2"', b'<spectrum index="1" id="scan=1"')
    assert "spectrum 'scan=1'" in only(check(tmp_path, data), "MZV012").message


def test_mzv020_length_mismatch(tmp_path):
    data = edit(base_plain(), b'id="scan=2" defaultArrayLength="2"', b'id="scan=2" defaultArrayLength="3"')
    report = check(tmp_path, data)
    assert report.codes == ["MZV020", "MZV020"]                 # both arrays of scan=2
    assert all("spectrum 'scan=2'" in f.message for f in report.findings)


def test_mzv021_undecodable(tmp_path):
    numpress = b'accession="MS:1002312" name="MS-Numpress linear prediction compression"'
    data = edit(base_plain(), array_cv(cv.NO_COMPRESSION), numpress)
    f = only(check(tmp_path, data), "MZV021")
    assert "spectrum 'scan=2'" in f.message and "out of scope" in f.message


def test_mzv021_bad_base64(tmp_path):
    data = base_plain()
    i = data.index(b"<binary>", data.index(b'id="scan=2"')) + len(b"<binary>")
    data = data[:i] + b"!" + data[i + 1:]
    assert "invalid base64" in only(check(tmp_path, data), "MZV021").message


def test_mzv022_encoded_length(tmp_path):
    data = base_plain()
    m = re.search(rb'<binaryDataArray encodedLength="(\d+)"', data)
    data = data[:m.start(1)] + str(int(m.group(1)) + 4).encode() + data[m.end(1):]
    assert "encodedLength" in only(check(tmp_path, data), "MZV022").message


def test_mzv023_missing_array(tmp_path):
    data = base_plain()
    start = data.index(b'id="scan=2"')
    second = data.index(b"<binaryDataArray ", data.index(b"<binaryDataArray ", start) + 1)
    end = data.index(b"</binaryDataArray>", second) + len(b"</binaryDataArray>")
    f = only(check(tmp_path, data[:second] + data[end:]), "MZV023")
    assert "spectrum 'scan=2'" in f.message and "intensity array" in f.message


def test_mzv024_mz_and_intensity_lengths_differ(tmp_path):
    # scan=2 holds 2 points; drop one intensity value and declare it with arrayLength
    data = base_plain()
    start = data.index(b'id="scan=2"')
    second = data.index(b"<binaryDataArray ", data.index(b"<binaryDataArray ", start) + 1)
    b0 = data.index(b"<binary>", second) + len(b"<binary>")
    b1 = data.index(b"</binary>", b0)
    text = mzml.encode_binary([5.0], 64, False).encode()
    data = data[:b0] + text + data[b1:]
    data = data[:second] + b'<binaryDataArray arrayLength="1" encodedLength="%d">' % len(text) \
        + data[data.index(b">", second) + 1:]
    f = only(check(tmp_path, data), "MZV024")
    assert "2 m/z values but 1 intensities" in f.message


def test_mzv030_no_representation(tmp_path):
    data = edit(base_plain(), gen.cv_param(cv.CENTROID_SPECTRUM).encode(), b"")
    f = only(check(tmp_path, data), "MZV030")
    assert f.level == "warning" and "spectrum 'scan=1'" in f.message


def test_mzv031_mz_decreasing(tmp_path):
    report = mzv.validate(INVALID / "invalid_mz_order.mzML")
    assert "point 11" in only(report, "MZV031").message


def test_mzv032_negative_intensity(tmp_path):
    data = base_plain()
    start = data.index(b'id="scan=2"')
    second = data.index(b"<binaryDataArray ", data.index(b"<binaryDataArray ", start) + 1)
    b0 = data.index(b"<binary>", second) + len(b"<binary>")
    b1 = data.index(b"</binary>", b0)
    data = data[:b0] + mzml.encode_binary([-5.0, 20.0], 64, False).encode() + data[b1:]
    f = only(check(tmp_path, data), "MZV032")
    assert f.level == "warning" and "1 negative" in f.message


def test_mzv033_no_ms_level(tmp_path):
    data = edit(base_plain(), gen.cv_param(cv.MS_LEVEL, 1).encode(), b"")
    assert "spectrum 'scan=1'" in only(check(tmp_path, data), "MZV033").message


def test_mzv034_msn_without_precursor(tmp_path):
    data = base_plain()
    data = data[:data.index(b"<precursorList")] + data[data.index(b"</precursorList>") + len(b"</precursorList>"):]
    assert "spectrum 'scan=2'" in only(check(tmp_path, data), "MZV034").message


def test_mzv040_chromatogram_lengths(tmp_path):
    chrom_start = base_plain().index(b"<chromatogram ")
    data = base_plain()
    second = data.index(b"<binaryDataArray ", data.index(b"<binaryDataArray ", chrom_start) + 1)
    b0 = data.index(b"<binary>", second) + len(b"<binary>")
    b1 = data.index(b"</binary>", b0)
    text = mzml.encode_binary([1.0], 64, True).encode()
    data = data[:b0] + text + data[b1:]
    data = data[:second] + b'<binaryDataArray arrayLength="1" encodedLength="%d">' % len(text) \
        + data[data.index(b">", second) + 1:]
    assert "chromatogram 'TIC'" in only(check(tmp_path, data), "MZV040").message


def test_mzv041_chromatogram_time_decreasing(tmp_path):
    spectra, _ = base_spectra()
    data = gen.plain_mzml("t", spectra, [gen.tic_chromatogram_xml([1.01, 1.0], [1.0, 2.0], 64)], has_msn=True)
    assert "chromatogram 'TIC'" in only(check(tmp_path, data), "MZV041").message


def test_mzv050_bad_offset(tmp_path):
    report = mzv.validate(INVALID / "invalid_index_offset.mzML")
    assert "spectrum 'scan=3'" in only(report, "MZV050").message


def test_mzv050_offset_points_at_other_id(tmp_path):
    data = base_indexed()
    pos2 = re.search(rb'<offset idRef="scan=2">(\d+)<', data).group(1)
    data = re.sub(rb'(<offset idRef="scan=1">)\d+', rb"\g<1>" + pos2, data)
    data = bad.fix_checksum(data)
    assert "spectrum 'scan=1'" in only(check(tmp_path, data), "MZV050").message


def test_mzv051_index_list_offset(tmp_path):
    data = base_indexed()
    m = re.search(rb"<indexListOffset>(\d+)<", data)
    data = bad.fix_checksum(data[:m.start(1)] + str(int(m.group(1)) - 1).encode() + data[m.end(1):])
    only(check(tmp_path, data), "MZV051")


def test_mzv052_checksum(tmp_path):
    report = mzv.validate(INVALID / "invalid_checksum.mzML")
    assert "fileChecksum" in only(report, "MZV052").message


def test_mzv053_not_in_index(tmp_path):
    data = base_indexed()
    data = re.sub(rb'\s*<index name="chromatogram">.*?</index>', b"", data, flags=re.S)
    data = edit(data, b'<indexList count="2">', b'<indexList count="1">')
    data = bad.fix_checksum(data)
    f = only(check(tmp_path, data), "MZV053")
    assert f.level == "warning" and "'TIC'" in f.message


def test_edit_in_indexed_file_keeps_index_valid(tmp_path):
    """The test helper itself: an edit that shifts bytes, then reindex, leaves only the intended finding."""
    data = edit(base_indexed(), b'<spectrumList count="2"', b'<spectrumList count="12"', indexed=True)
    only(check(tmp_path, data), "MZV010")


def test_every_code_is_tested():
    tested = set(re.findall(r"def test_(mzv\d{3})_", Path(__file__).read_text(encoding="utf-8")))
    assert {c.lower() for c in mzv.CHECKS} <= tested


# ------------------------------------------------------------------ invalid samples
@pytest.mark.parametrize("name, code", [(n, c) for n, (_, c) in bad.INVALID.items()])
def test_invalid_sample_yields_exactly_its_finding(name, code):
    only(mzv.validate(INVALID / name), code)


def test_invalid_samples_are_reproducible(tmp_path):
    for fresh in bad.generate(tmp_path):
        shipped = INVALID / fresh.name
        if b"zlib compression" in fresh.read_bytes():
            # zlib output may differ between zlib builds; the defect found must not
            assert mzv.validate(fresh).codes == mzv.validate(shipped).codes
        else:
            assert fresh.read_bytes() == shipped.read_bytes(), fresh.name


def test_invalid_samples_are_otherwise_readable():
    """Apart from the indexed defects, the reader still opens the corrupted files where it can."""
    for name in ("invalid_spectrum_count.mzML", "invalid_duplicate_id.mzML", "invalid_mz_order.mzML",
                 "invalid_index_offset.mzML", "invalid_checksum.mzML"):
        assert mzml.read(INVALID / name).spectra
    with pytest.raises(mzml.MzMLError, match="arrayLength|declared array length"):
        mzml.read(INVALID / "invalid_array_length.mzML")


# ------------------------------------------------------------------ strict mode and CLI
def test_strict_turns_warnings_into_errors(tmp_path):
    data = edit(base_plain(), gen.cv_param(cv.CENTROID_SPECTRUM).encode(), b"")
    assert check(tmp_path, data).ok
    report = check(tmp_path, data, strict=True)
    assert not report.ok and report.findings[0].level == "error"


def test_cli_exit_codes(tmp_path, capsys):
    good = [str(SAMPLES / n) for n in ("synthetic_caffeine_profile.mzML", "synthetic_caffeine_ms2_indexed.mzML")]
    assert mzv.main(good) == 0
    assert "validation : OK" in capsys.readouterr().out
    assert mzv.main([str(INVALID / "invalid_checksum.mzML")]) == 1
    assert "ERROR: MZV052:" in capsys.readouterr().out
    assert mzv.main([str(tmp_path / "missing.mzML")]) == 1


def test_cli_strict(tmp_path, capsys):
    path = tmp_path / "warn.mzML"
    path.write_bytes(edit(base_plain(), gen.cv_param(cv.CENTROID_SPECTRUM).encode(), b""))
    assert mzv.main([str(path)]) == 0
    assert "1 warning(s)" in capsys.readouterr().out
    assert mzv.main([str(path), "--strict"]) == 1
    assert "ERROR: MZV030:" in capsys.readouterr().out
