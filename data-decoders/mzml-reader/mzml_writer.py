"""Minimal mzML 1.1 writer, standard library only.

Writes :class:`mzml.Spectrum` and :class:`mzml.Chromatogram` objects (as the
reader returns them, or built by hand) to a plain ``<mzML>`` file or an
``<indexedmzML>`` file with byte offsets, ``indexListOffset`` and a SHA-1
``fileChecksum``. Arrays are 64-bit floats by default (32-bit on request),
with zlib or no compression.

The writer checks its input first and raises :class:`ValueError` for anything
that would make an invalid or ambiguous file: duplicate ids, a missing MS level
or centroid/profile flag, an MS2 spectrum without a precursor, arrays of
unequal length, decreasing m/z or time values, negative intensities, and units
or CV terms it has no accession for. So every file it writes reads back to the
same values and passes ``mzml_validator.py --strict`` with no findings.

File metadata uses generic PSI-MS terms only (*instrument model*, *custom
unreleased software tool*, *data processing action*); it describes no real
instrument.

Not supported: MS-Numpress and other compression schemes, integer and string
arrays, and any instrument-specific metadata.
"""
from __future__ import annotations

import hashlib
import io
import math
from pathlib import Path
from typing import BinaryIO, Dict, Iterable, List, Optional, Sequence, Tuple, Union
from xml.sax.saxutils import escape

import mzml
import psi_ms_cv as cv

__all__ = ["write", "to_bytes"]

MS_NS = "http://psi.hupo.org/ms/mzml"
CV_LIST = [
    ("MS", "Proteomics Standards Initiative Mass Spectrometry Ontology", cv.PSI_MS_VERSION,
     "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"),
    ("UO", "Unit Ontology", None,
     "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"),
]

_TIME_UNITS = {name: (acc, name) for acc, name in cv.TIME_UNITS.items()}
_INTENSITY_UNITS = {cv.DETECTOR_COUNTS[1]: cv.DETECTOR_COUNTS}
_CHROMATOGRAM_KINDS = {t[1]: t for t in (cv.TIC_CHROMATOGRAM, cv.BASEPEAK_CHROMATOGRAM,
                                         cv.SELECTED_ION_CURRENT_CHROMATOGRAM)}
_ACTIVATIONS = {t[1]: t for t in (cv.COLLISION_INDUCED_DISSOCIATION, cv.BEAM_TYPE_COLLISION_INDUCED_DISSOCIATION)}


# ------------------------------------------------------------------ input checks
def _unit(table: Dict[str, tuple], name: Optional[str], what: str, where: str) -> Optional[tuple]:
    if name is None:
        return None
    if name not in table:
        raise ValueError(f"{where}: no CV accession known for {what} {name!r} (known: {', '.join(sorted(table))})")
    return table[name]


def _check_finite(values, what: str, where: str) -> None:
    if any(not math.isfinite(v) for v in values):
        raise ValueError(f"{where}: {what} contains NaN or infinite values")


def _non_decreasing(values, what: str, where: str) -> None:
    for i in range(1, len(values)):
        if values[i] < values[i - 1]:
            raise ValueError(f"{where}: {what} decreases at point {i}")


def _check_spectrum(sp: mzml.Spectrum) -> None:
    where = f"spectrum '{sp.id}'"
    if not isinstance(sp.ms_level, int) or sp.ms_level < 1:
        raise ValueError(f"{where}: ms_level must be a positive integer, not {sp.ms_level!r}")
    if sp.centroid is None:
        raise ValueError(f"{where}: centroid must be True (centroid) or False (profile)")
    if sp.ms_level >= 2 and not sp.precursors:
        raise ValueError(f"{where}: an MS{sp.ms_level} spectrum needs at least one precursor")
    if len(sp.mz) != len(sp.intensity):
        raise ValueError(f"{where}: {len(sp.mz)} m/z values but {len(sp.intensity)} intensities")
    _check_finite(sp.mz, "the m/z array", where)
    _check_finite(sp.intensity, "the intensity array", where)
    _non_decreasing(sp.mz, "the m/z array", where)
    if any(v < 0 for v in sp.intensity):
        raise ValueError(f"{where}: negative intensities")
    if sp.scan_start_time is not None:
        _check_finite([sp.scan_start_time], "the scan start time", where)
        if sp.scan_start_time_unit is None:
            raise ValueError(f"{where}: a scan start time needs its unit (minute, second or millisecond)")
        _unit(_TIME_UNITS, sp.scan_start_time_unit, "time unit", where)
    _unit(_INTENSITY_UNITS, sp.intensity_unit, "intensity unit", where)
    for p in sp.precursors:
        for name in p.activation:
            _unit(_ACTIVATIONS, name, "activation", where)


def _check_chromatogram(ch: mzml.Chromatogram) -> None:
    where = f"chromatogram '{ch.id}'"
    if len(ch.time) != len(ch.intensity):
        raise ValueError(f"{where}: {len(ch.time)} time values but {len(ch.intensity)} intensities")
    _check_finite(ch.time, "the time array", where)
    _check_finite(ch.intensity, "the intensity array", where)
    _non_decreasing(ch.time, "the time array", where)
    _unit(_TIME_UNITS, ch.time_unit, "time unit", where)
    _unit(_INTENSITY_UNITS, ch.intensity_unit, "intensity unit", where)
    _unit(_CHROMATOGRAM_KINDS, ch.kind, "chromatogram type", where)


def _check_ids(items: Sequence, kind: str) -> None:
    seen = set()
    for it in items:
        if not isinstance(it.id, str) or not it.id:
            raise ValueError(f"every {kind} needs a non-empty string id, got {it.id!r}")
        if it.id in seen:
            raise ValueError(f"duplicate {kind} id {it.id!r}")
        seen.add(it.id)


# ------------------------------------------------------------------ XML text
def _attr(value) -> str:
    return escape(str(value), {'"': "&quot;"})


def _num(v) -> str:
    """Exact text for a number: repr() round-trips every float."""
    return repr(float(v)) if not isinstance(v, int) else str(v)


def _cv(term: tuple, value="", unit: Optional[tuple] = None) -> str:
    s = (f'<cvParam cvRef="{term[0].split(":")[0]}" accession="{term[0]}" name="{_attr(term[1])}" '
         f'value="{_attr(value)}"')
    if unit:
        s += f' unitCvRef="{unit[0].split(":")[0]}" unitAccession="{unit[0]}" unitName="{_attr(unit[1])}"'
    return s + "/>"


class _Sink:
    """Writes text as bytes, tracking the byte position and a running SHA-1."""

    def __init__(self, fh: BinaryIO):
        self.fh, self.pos, self.sha1 = fh, 0, hashlib.sha1()

    def write(self, text: str) -> None:
        data = text.encode("utf-8")
        self.fh.write(data)
        self.sha1.update(data)
        self.pos += len(data)


def _array_xml(values, array_type: tuple, unit: Optional[tuple], bits: int, zlib_compressed: bool, pad: str) -> List[str]:
    text = mzml.encode_binary(values, bits, zlib_compressed)
    return [f'{pad}<binaryDataArray encodedLength="{len(text)}">',
            f"{pad}  " + _cv(cv.FLOAT_32 if bits == 32 else cv.FLOAT_64),
            f"{pad}  " + _cv(cv.ZLIB_COMPRESSION if zlib_compressed else cv.NO_COMPRESSION),
            f"{pad}  " + _cv(array_type, unit=unit),
            f"{pad}  <binary>{text}</binary>",
            f"{pad}</binaryDataArray>"]


def _precursor_xml(p: mzml.Precursor, pad: str) -> List[str]:
    ref = f' spectrumRef="{_attr(p.spectrum_ref)}"' if p.spectrum_ref else ""
    out = [f"{pad}<precursor{ref}>"]
    window = [(cv.ISOLATION_WINDOW_TARGET_MZ, p.isolation_target_mz), (cv.ISOLATION_WINDOW_LOWER_OFFSET, p.isolation_lower_offset),
              (cv.ISOLATION_WINDOW_UPPER_OFFSET, p.isolation_upper_offset)]
    if any(v is not None for _, v in window):
        out.append(f"{pad}  <isolationWindow>")
        out += [f"{pad}    " + _cv(t, _num(v), cv.MZ_UNIT) for t, v in window if v is not None]
        out.append(f"{pad}  </isolationWindow>")
    ion = []
    if p.selected_ion_mz is not None:
        ion.append(_cv(cv.SELECTED_ION_MZ, _num(p.selected_ion_mz), cv.MZ_UNIT))
    if p.charge is not None:
        ion.append(_cv(cv.CHARGE_STATE, int(p.charge)))
    if p.intensity is not None:
        ion.append(_cv(cv.PEAK_INTENSITY, _num(p.intensity), cv.DETECTOR_COUNTS))
    if ion:
        out += [f'{pad}  <selectedIonList count="1">', f"{pad}    <selectedIon>"]
        out += [f"{pad}      {x}" for x in ion]
        out += [f"{pad}    </selectedIon>", f"{pad}  </selectedIonList>"]
    out.append(f"{pad}  <activation>")
    out += [f"{pad}    " + _cv(_ACTIVATIONS[name]) for name in p.activation]
    if p.collision_energy is not None:
        out.append(f"{pad}    " + _cv(cv.COLLISION_ENERGY, _num(p.collision_energy), cv.UO_ELECTRONVOLT))
    out += [f"{pad}  </activation>", f"{pad}</precursor>"]
    return out


def _spectrum_xml(sp: mzml.Spectrum, index: int, bits: int, zlib_compressed: bool) -> Tuple[str, List[str]]:
    pad = "      "
    head = f'{pad}<spectrum index="{index}" id="{_attr(sp.id)}" defaultArrayLength="{len(sp.mz)}">'
    body = [f"{pad}  " + _cv(cv.MS_LEVEL, sp.ms_level),
            f"{pad}  " + _cv(cv.MS1_SPECTRUM if sp.ms_level == 1 else cv.MSN_SPECTRUM),
            f"{pad}  " + _cv(cv.CENTROID_SPECTRUM if sp.centroid else cv.PROFILE_SPECTRUM)]
    if sp.total_ion_current is not None:
        body.append(f"{pad}  " + _cv(cv.TOTAL_ION_CURRENT, _num(sp.total_ion_current)))
    if sp.scan_start_time is not None:
        body += [f'{pad}  <scanList count="1">', f"{pad}    " + _cv(cv.NO_COMBINATION), f"{pad}    <scan>",
                 f"{pad}      " + _cv(cv.SCAN_START_TIME, _num(sp.scan_start_time), _TIME_UNITS[sp.scan_start_time_unit]),
                 f"{pad}    </scan>", f"{pad}  </scanList>"]
    if sp.precursors:
        body.append(f'{pad}  <precursorList count="{len(sp.precursors)}">')
        for p in sp.precursors:
            body += _precursor_xml(p, pad + "    ")
        body.append(f"{pad}  </precursorList>")
    body.append(f'{pad}  <binaryDataArrayList count="2">')
    body += _array_xml(sp.mz, cv.MZ_ARRAY, cv.MZ_UNIT, bits, zlib_compressed, pad + "    ")
    body += _array_xml(sp.intensity, cv.INTENSITY_ARRAY, _INTENSITY_UNITS.get(sp.intensity_unit),
                       bits, zlib_compressed, pad + "    ")
    body += [f"{pad}  </binaryDataArrayList>", f"{pad}</spectrum>"]
    return head, body


def _chromatogram_xml(ch: mzml.Chromatogram, index: int, bits: int, zlib_compressed: bool) -> Tuple[str, List[str]]:
    pad = "      "
    head = f'{pad}<chromatogram index="{index}" id="{_attr(ch.id)}" defaultArrayLength="{len(ch.time)}">'
    body = [f"{pad}  " + _cv(_CHROMATOGRAM_KINDS[ch.kind])] if ch.kind else []
    body.append(f'{pad}  <binaryDataArrayList count="2">')
    body += _array_xml(ch.time, cv.TIME_ARRAY, _TIME_UNITS.get(ch.time_unit) if ch.time_unit else None,
                       bits, zlib_compressed, pad + "    ")
    body += _array_xml(ch.intensity, cv.INTENSITY_ARRAY, _INTENSITY_UNITS.get(ch.intensity_unit),
                       bits, zlib_compressed, pad + "    ")
    body += [f"{pad}  </binaryDataArrayList>", f"{pad}</chromatogram>"]
    return head, body


def _header(run_id: str, spectra: List[mzml.Spectrum]) -> List[str]:
    content = []
    if any(s.ms_level == 1 for s in spectra):
        content.append("      " + _cv(cv.MS1_SPECTRUM))
    if any(s.ms_level >= 2 for s in spectra):
        content.append("      " + _cv(cv.MSN_SPECTRUM))
    cvs = [f'    <cv id="{cid}" fullName="{full}"' + (f' version="{v}"' if v else "") + f' URI="{uri}"/>'
           for cid, full, v, uri in CV_LIST]
    return [f'  <cvList count="{len(cvs)}">', *cvs, "  </cvList>",
            "  <fileDescription>", "    <fileContent>", *content, "    </fileContent>", "  </fileDescription>",
            '  <softwareList count="1">', '    <software id="mzml_writer" version="1.0">',
            "      " + _cv(cv.CUSTOM_SOFTWARE, "mzml_writer.py"), "    </software>", "  </softwareList>",
            '  <instrumentConfigurationList count="1">', '    <instrumentConfiguration id="IC1">',
            "      " + _cv(cv.INSTRUMENT_MODEL), "    </instrumentConfiguration>", "  </instrumentConfigurationList>",
            '  <dataProcessingList count="1">', '    <dataProcessing id="DP1">',
            '      <processingMethod order="1" softwareRef="mzml_writer">',
            "        " + _cv(cv.DATA_PROCESSING_ACTION, "written by mzml_writer.py"),
            "      </processingMethod>", "    </dataProcessing>", "  </dataProcessingList>",
            f'  <run id="{_attr(run_id)}" defaultInstrumentConfigurationRef="IC1">']


# ------------------------------------------------------------------ public API
def write(target: Union[str, Path, BinaryIO], spectra: Iterable[mzml.Spectrum],
          chromatograms: Iterable[mzml.Chromatogram] = (), *, indexed: bool = False, bits: int = 64,
          zlib_compressed: bool = True, run_id: str = "run") -> None:
    """Write spectra (and chromatograms) as mzML 1.1 to a path or binary file.

    ``indexed=True`` wraps the file in ``<indexedmzML>`` with byte offsets for
    every spectrum and chromatogram, ``indexListOffset`` and a SHA-1 checksum.
    Index values are renumbered 0, 1, 2 ... in the order given. Raises
    :class:`ValueError` on input that would give an invalid file.
    """
    spectra, chromatograms = list(spectra), list(chromatograms)
    if bits not in (32, 64):
        raise ValueError(f"bits must be 32 or 64, not {bits!r}")
    _check_ids(spectra, "spectrum")
    _check_ids(chromatograms, "chromatogram")
    for sp in spectra:
        _check_spectrum(sp)
    for ch in chromatograms:
        _check_chromatogram(ch)
    if isinstance(target, (str, Path)):
        with open(target, "wb") as fh:
            _write(_Sink(fh), spectra, chromatograms, indexed, bits, zlib_compressed, run_id)
    else:
        _write(_Sink(target), spectra, chromatograms, indexed, bits, zlib_compressed, run_id)


def to_bytes(spectra: Iterable[mzml.Spectrum], chromatograms: Iterable[mzml.Chromatogram] = (), **options) -> bytes:
    """The same as :func:`write`, returning the file content."""
    buf = io.BytesIO()
    write(buf, spectra, chromatograms, **options)
    return buf.getvalue()


def _write(out: _Sink, spectra, chromatograms, indexed, bits, zlib_compressed, run_id) -> None:
    offsets: Dict[str, List[Tuple[str, int]]] = {"spectrum": [], "chromatogram": []}
    out.write('<?xml version="1.0" encoding="utf-8"?>\n')
    if indexed:
        out.write(f'<indexedmzML xmlns="{MS_NS}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation="{MS_NS} http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">\n')
        out.write(f'<mzML id="{_attr(run_id)}" version="1.1.0">\n')
    else:
        out.write(f'<mzML xmlns="{MS_NS}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                  f'xsi:schemaLocation="{MS_NS} http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" '
                  f'id="{_attr(run_id)}" version="1.1.0">\n')
    for line in _header(run_id, spectra):
        out.write(line + "\n")

    def items(kind, objs, build):
        for i, obj in enumerate(objs):
            head, body = build(obj, i, bits, zlib_compressed)
            out.write(head[: len(head) - len(head.lstrip())])
            offsets[kind].append((obj.id, out.pos))
            out.write(head.lstrip() + "\n")
            out.write("\n".join(body) + "\n")

    out.write(f'    <spectrumList count="{len(spectra)}" defaultDataProcessingRef="DP1">\n')
    items("spectrum", spectra, _spectrum_xml)
    out.write("    </spectrumList>\n")
    if chromatograms:
        out.write(f'    <chromatogramList count="{len(chromatograms)}" defaultDataProcessingRef="DP1">\n')
        items("chromatogram", chromatograms, _chromatogram_xml)
        out.write("    </chromatogramList>\n")
    out.write("  </run>\n</mzML>\n")
    if not indexed:
        return
    kinds = [k for k in ("spectrum", "chromatogram") if offsets[k]] or ["spectrum"]
    out.write("  ")
    index_list_offset = out.pos
    out.write(f'<indexList count="{len(kinds)}">\n')
    for kind in kinds:
        out.write(f'    <index name="{kind}">\n')
        for ref, pos in offsets[kind]:
            out.write(f'      <offset idRef="{_attr(ref)}">{pos}</offset>\n')
        out.write("    </index>\n")
    out.write(f"  </indexList>\n  <indexListOffset>{index_list_offset}</indexListOffset>\n  <fileChecksum>")
    out.write(f"{out.sha1.hexdigest()}</fileChecksum>\n</indexedmzML>\n")
