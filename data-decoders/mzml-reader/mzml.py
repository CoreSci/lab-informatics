"""Streaming mzML reader, standard library only.

mzML is the open HUPO-PSI XML format for mass spectrometry data. Each spectrum
carries its metadata as controlled-vocabulary terms (``cvParam``) and its
m/z and intensity values as base64-encoded binary arrays. This module
implements the part of the mzML 1.1 specification needed to read spectra:

* plain ``<mzML>`` files and the ``<indexedmzML>`` wrapper (the index itself
  is ignored; spectra are read in document order);
* binary arrays: base64, no compression or zlib, 32-bit or 64-bit
  little-endian floats, with declared lengths checked;
* per spectrum: id, index, MS level, centroid or profile, scan start time and
  its unit, m/z and intensity arrays (with the intensity unit), total ion
  current and precursor information (MS2);
* chromatograms: id, index, type, time and intensity arrays and their units;
* ``referenceableParamGroup`` references, expanded where they are used.

Parsing is streaming: :func:`iterparse` yields one spectrum or chromatogram at
a time and frees each element after use, so memory stays flat on large files.

Not supported (rejected with a clear :class:`MzMLError`): MS-Numpress and any
other compression except zlib, data types other than 32/64-bit floats, and
array types other than m/z, intensity and time.

Command line::

    python mzml.py FILE [FILE ...]

prints a summary of each file. It exits with status 1 when any file can't be read.
"""
from __future__ import annotations

import argparse
import base64
import binascii
import sys
import xml.etree.ElementTree as ET
import zlib
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO, Dict, Iterator, List, Optional, Union

import psi_ms_cv as cv

__all__ = ["MzMLError", "Precursor", "Spectrum", "Chromatogram", "Run",
           "decode_binary", "encode_binary", "iterparse", "iter_spectra", "iter_chromatograms", "read"]


class MzMLError(ValueError):
    """Unreadable or unsupported mzML content."""


@dataclass
class Precursor:
    """Precursor of an MSn spectrum: the isolation window, the selected ion and the activation."""

    spectrum_ref: Optional[str] = None
    selected_ion_mz: Optional[float] = None
    charge: Optional[int] = None
    intensity: Optional[float] = None
    isolation_target_mz: Optional[float] = None
    isolation_lower_offset: Optional[float] = None
    isolation_upper_offset: Optional[float] = None
    activation: List[str] = field(default_factory=list)
    collision_energy: Optional[float] = None


@dataclass
class Spectrum:
    """One spectrum. ``centroid`` is True (centroid), False (profile) or None (not stated)."""

    id: str
    index: int
    ms_level: Optional[int]
    centroid: Optional[bool]
    scan_start_time: Optional[float]
    scan_start_time_unit: Optional[str]
    mz: array
    intensity: array
    total_ion_current: Optional[float] = None
    precursors: List[Precursor] = field(default_factory=list)
    intensity_unit: Optional[str] = None


@dataclass
class Chromatogram:
    """One chromatogram. ``kind`` is its CV type, e.g. 'total ion current chromatogram'."""

    id: str
    index: int
    kind: Optional[str]
    time: array
    time_unit: Optional[str]
    intensity: array
    intensity_unit: Optional[str] = None


@dataclass
class Run:
    """Everything in a file, as returned by :func:`read`."""

    spectra: List[Spectrum]
    chromatograms: List[Chromatogram]
    indexed: bool
    run_id: Optional[str]


Source = Union[str, Path, BinaryIO]


# ------------------------------------------------------------------ binary arrays
def _bits_typecode(bits: int) -> str:
    if bits == 32:
        return "f"
    if bits == 64:
        return "d"
    raise MzMLError(f"unsupported float size: {bits}-bit (only 32-bit and 64-bit floats are handled)")


def decode_binary(text: str, bits: int = 64, zlib_compressed: bool = False,
                  expected_length: Optional[int] = None) -> array:
    """Decode one mzML binary array (base64, optional zlib, little-endian floats).

    Returns an ``array('d')``. Raises :class:`MzMLError` on invalid base64,
    zlib data that won't decompress, a byte count that isn't a whole number of
    values, or a value count different from ``expected_length``.
    """
    typecode = _bits_typecode(bits)
    try:
        raw = base64.b64decode("".join(text.split()), validate=True)
    except (binascii.Error, ValueError) as exc:
        raise MzMLError(f"invalid base64 data ({exc})") from None
    if zlib_compressed:
        try:
            raw = zlib.decompress(raw)
        except zlib.error as exc:
            raise MzMLError(f"zlib data could not be decompressed ({exc})") from None
    size = bits // 8
    if len(raw) % size:
        raise MzMLError(f"{len(raw)} bytes is not a whole number of {bits}-bit values")
    values = array(typecode)
    values.frombytes(raw)
    if sys.byteorder == "big":
        values.byteswap()                      # mzML binary data is always little-endian
    if expected_length is not None and len(values) != expected_length:
        raise MzMLError(f"decoded {len(values)} values, but the declared array length is {expected_length}")
    return values if typecode == "d" else array("d", values)


def encode_binary(values, bits: int = 64, zlib_compressed: bool = False) -> str:
    """Encode floats as an mzML binary array: little-endian, optional zlib, base64.

    The counterpart of :func:`decode_binary`, used by the sample generator and the tests.
    """
    packed = array(_bits_typecode(bits), values)
    if sys.byteorder == "big":
        packed.byteswap()
    raw = packed.tobytes()
    if zlib_compressed:
        raw = zlib.compress(raw)
    return base64.b64encode(raw).decode("ascii")


# ------------------------------------------------------------------ XML helpers
def _local(tag: str) -> str:
    """Tag name without its namespace: '{http://psi.hupo.org/ms/mzml}spectrum' -> 'spectrum'."""
    return tag.rsplit("}", 1)[-1]


def _children(elem: ET.Element, name: str) -> List[ET.Element]:
    return [c for c in elem if _local(c.tag) == name]


def _child(elem: ET.Element, name: str) -> Optional[ET.Element]:
    for c in elem:
        if _local(c.tag) == name:
            return c
    return None


def _params(elem: Optional[ET.Element], groups: Dict[str, list]) -> List[dict]:
    """The element's cvParams, with referenced param groups expanded in place."""
    if elem is None:
        return []
    out: List[dict] = []
    for c in elem:
        tag = _local(c.tag)
        if tag == "cvParam":
            out.append(dict(c.attrib))
        elif tag == "referenceableParamGroupRef":
            ref = c.get("ref")
            if ref not in groups:
                raise MzMLError(f"reference to unknown referenceableParamGroup '{ref}'")
            out.extend(groups[ref])
    return out


def _find(params: List[dict], term: tuple) -> Optional[dict]:
    for p in params:
        if p.get("accession") == term[0]:
            return p
    return None


def _float(params: List[dict], term: tuple) -> Optional[float]:
    p = _find(params, term)
    if p is None or p.get("value") in (None, ""):
        return None
    try:
        return float(p["value"])
    except ValueError:
        raise MzMLError(f"'{term[1]}' has a non-numeric value {p['value']!r}") from None


def _int_attr(elem: ET.Element, name: str, where: str) -> int:
    value = elem.get(name)
    if value is None:
        raise MzMLError(f"{where}: missing required attribute '{name}'")
    try:
        return int(value)
    except ValueError:
        raise MzMLError(f"{where}: attribute '{name}' is not an integer: {value!r}") from None


def _term_label(p: dict) -> str:
    return f"'{p.get('name', '?')}' ({p.get('accession', '?')})"


# ------------------------------------------------------------------ binaryDataArray
_ARRAY_TYPES = {cv.MZ_ARRAY[0]: cv.MZ_ARRAY[1], cv.INTENSITY_ARRAY[0]: cv.INTENSITY_ARRAY[1],
                cv.TIME_ARRAY[0]: cv.TIME_ARRAY[1]}
_FLOAT_BITS = {cv.FLOAT_32[0]: 32, cv.FLOAT_64[0]: 64}
_COMPRESSIONS = {cv.NO_COMPRESSION[0]: False, cv.ZLIB_COMPRESSION[0]: True}


def _array_spec(params: List[dict], where: str):
    """(array name, bits, zlib?, unit) from a binaryDataArray's cvParams."""
    name = bits = compressed = None
    unit = None
    for p in params:
        acc, pname = p.get("accession", ""), p.get("name", "")
        if acc in _ARRAY_TYPES:
            if name is not None:
                raise MzMLError(f"{where}: more than one array type")
            name, unit = _ARRAY_TYPES[acc], p.get("unitName")
        elif pname.endswith(" array"):
            # Every array-type term in the PSI-MS CV is named '... array'
            # (children of 'binary data array', MS:1000513).
            raise MzMLError(f"{where}: unsupported array type {_term_label(p)}. "
                            "This reader handles m/z, intensity and time arrays only")
        elif acc in _FLOAT_BITS:
            bits = _FLOAT_BITS[acc]
        elif pname.endswith(("-bit float", "-bit integer", "ASCII string")):
            raise MzMLError(f"{where}: unsupported binary data type {_term_label(p)}. "
                            "Only 32-bit and 64-bit floats are handled")
        elif acc in _COMPRESSIONS:
            compressed = _COMPRESSIONS[acc]
        elif acc in cv.NUMPRESS_ACCESSIONS or "numpress" in pname.lower():
            raise MzMLError(f"{where}: {_term_label(p)} is out of scope for this reader. "
                            "Only zlib and no compression are handled; MS-Numpress is not")
        elif "compression" in pname or pname.endswith("encoding"):
            # Other children of 'binary data compression type' (MS:1000572).
            raise MzMLError(f"{where}: compression {_term_label(p)} is out of scope for this reader. "
                            "Only zlib and no compression are handled")
    if name is None:
        raise MzMLError(f"{where}: missing array type (expected m/z, intensity or time array)")
    if bits is None:
        raise MzMLError(f"{where}: missing binary data type (32-bit or 64-bit float)")
    if compressed is None:
        raise MzMLError(f"{where}: missing compression term (zlib or no compression)")
    return name, bits, compressed, unit


def _arrays(elem: ET.Element, default_length: int, allowed: tuple, groups, where: str) -> Dict[str, tuple]:
    """Decode every binaryDataArray of a spectrum or chromatogram: {name: (values, unit)}."""
    out: Dict[str, tuple] = {}
    bdal = _child(elem, "binaryDataArrayList")
    for i, bda in enumerate(_children(bdal, "binaryDataArray") if bdal is not None else []):
        here = f"{where}, binary array {i + 1}"
        name, bits, compressed, unit = _array_spec(_params(bda, groups), here)
        here = f"{where}, {name}"
        if name not in allowed:
            raise MzMLError(f"{here}: unexpected here (expected {' and '.join(allowed)})")
        if name in out:
            raise MzMLError(f"{here}: appears twice")
        length = _int_attr(bda, "arrayLength", here) if bda.get("arrayLength") is not None else default_length
        binary = _child(bda, "binary")
        text = "".join((binary.text or "").split()) if binary is not None else ""
        declared = bda.get("encodedLength")
        if declared is not None and declared.isdigit() and int(declared) != len(text):
            raise MzMLError(f"{here}: encodedLength is {declared}, but the base64 text has {len(text)} characters")
        try:
            values = decode_binary(text, bits, compressed, expected_length=length)
        except MzMLError as exc:
            raise MzMLError(f"{here}: {exc}") from None
        out[name] = (values, unit)
    for name in allowed:
        if name not in out:
            if default_length:
                raise MzMLError(f"{where}: no {name}")
            out[name] = (array("d"), None)
    return out


# ------------------------------------------------------------------ spectra and chromatograms
def _precursor(elem: ET.Element, groups) -> Precursor:
    pre = Precursor(spectrum_ref=elem.get("spectrumRef"))
    window = _params(_child(elem, "isolationWindow"), groups)
    pre.isolation_target_mz = _float(window, cv.ISOLATION_WINDOW_TARGET_MZ)
    pre.isolation_lower_offset = _float(window, cv.ISOLATION_WINDOW_LOWER_OFFSET)
    pre.isolation_upper_offset = _float(window, cv.ISOLATION_WINDOW_UPPER_OFFSET)
    sil = _child(elem, "selectedIonList")
    ions = _children(sil, "selectedIon") if sil is not None else []
    if ions:
        ion = _params(ions[0], groups)
        pre.selected_ion_mz = _float(ion, cv.SELECTED_ION_MZ)
        charge = _float(ion, cv.CHARGE_STATE)
        pre.charge = int(charge) if charge is not None else None
        pre.intensity = _float(ion, cv.PEAK_INTENSITY)
    activation = _params(_child(elem, "activation"), groups)
    pre.collision_energy = _float(activation, cv.COLLISION_ENERGY)
    pre.activation = [p.get("name", "") for p in activation if p.get("accession") != cv.COLLISION_ENERGY[0]]
    return pre


def _spectrum(elem: ET.Element, groups) -> Spectrum:
    sid = elem.get("id") or "?"
    index = _int_attr(elem, "index", f"spectrum '{sid}'")
    where = f"spectrum '{sid}' (index {index})"
    params = _params(elem, groups)
    level = _float(params, cv.MS_LEVEL)
    centroid = (True if _find(params, cv.CENTROID_SPECTRUM)
                else False if _find(params, cv.PROFILE_SPECTRUM) else None)
    rt = rt_unit = None
    scan_list = _child(elem, "scanList")
    scans = _children(scan_list, "scan") if scan_list is not None else []
    if scans:
        p = _find(_params(scans[0], groups), cv.SCAN_START_TIME)
        if p is not None:
            rt = _float([p], cv.SCAN_START_TIME)
            rt_unit = cv.TIME_UNITS.get(p.get("unitAccession", ""), p.get("unitName"))
    plist = _child(elem, "precursorList")
    precursors = [_precursor(p, groups) for p in _children(plist, "precursor")] if plist is not None else []
    arrays = _arrays(elem, _int_attr(elem, "defaultArrayLength", where),
                     (cv.MZ_ARRAY[1], cv.INTENSITY_ARRAY[1]), groups, where)
    return Spectrum(id=sid, index=index, ms_level=int(level) if level is not None else None,
                    centroid=centroid, scan_start_time=rt, scan_start_time_unit=rt_unit,
                    mz=arrays[cv.MZ_ARRAY[1]][0], intensity=arrays[cv.INTENSITY_ARRAY[1]][0],
                    total_ion_current=_float(params, cv.TOTAL_ION_CURRENT), precursors=precursors,
                    intensity_unit=arrays[cv.INTENSITY_ARRAY[1]][1])


def _chromatogram(elem: ET.Element, groups) -> Chromatogram:
    cid = elem.get("id") or "?"
    index = _int_attr(elem, "index", f"chromatogram '{cid}'")
    where = f"chromatogram '{cid}' (index {index})"
    kind = next((p.get("name") for p in _params(elem, groups)
                 if p.get("name", "").endswith("chromatogram")), None)
    arrays = _arrays(elem, _int_attr(elem, "defaultArrayLength", where),
                     (cv.TIME_ARRAY[1], cv.INTENSITY_ARRAY[1]), groups, where)
    time, unit = arrays[cv.TIME_ARRAY[1]]
    return Chromatogram(id=cid, index=index, kind=kind, time=time, time_unit=unit,
                        intensity=arrays[cv.INTENSITY_ARRAY[1]][0],
                        intensity_unit=arrays[cv.INTENSITY_ARRAY[1]][1])


# ------------------------------------------------------------------ streaming parser
def _stream(source: Source, info: dict) -> Iterator[Union[Spectrum, Chromatogram]]:
    if isinstance(source, (str, Path)):
        with open(source, "rb") as fh:
            yield from _stream(fh, info)
        return
    groups: Dict[str, list] = {}
    stack: List[ET.Element] = []
    try:
        for event, elem in ET.iterparse(source, events=("start", "end")):
            tag = _local(elem.tag)
            if event == "start":
                if not stack:
                    if tag not in ("mzML", "indexedmzML"):
                        raise MzMLError(f"not an mzML document: the root element is <{tag}>")
                    info["indexed"] = tag == "indexedmzML"
                elif tag == "run":
                    info["run_id"] = elem.get("id")
                stack.append(elem)
                continue
            stack.pop()
            if tag == "referenceableParamGroup":
                groups[elem.get("id", "")] = _params(elem, groups)
                continue
            if tag == "spectrum":
                item = _spectrum(elem, groups)
            elif tag == "chromatogram":
                item = _chromatogram(elem, groups)
            elif tag == "indexList":           # the index isn't used in this reader
                item = None
            else:
                continue
            # Free the element and detach it from its parent, so the tree never grows.
            elem.clear()
            if stack:
                stack[-1].remove(elem)
            if item is not None:
                yield item
    except ET.ParseError as exc:
        raise MzMLError(f"malformed XML: {exc}") from None


def iterparse(source: Source) -> Iterator[Union[Spectrum, Chromatogram]]:
    """Yield each spectrum and chromatogram in document order, one at a time.

    ``source`` is a path or a binary file object. Raises :class:`MzMLError`.
    """
    return _stream(source, {})


def iter_spectra(source: Source) -> Iterator[Spectrum]:
    """Yield each spectrum, one at a time (chromatograms are skipped)."""
    return (x for x in _stream(source, {}) if isinstance(x, Spectrum))


def iter_chromatograms(source: Source) -> Iterator[Chromatogram]:
    """Yield each chromatogram, one at a time (spectra are parsed and dropped)."""
    return (x for x in _stream(source, {}) if isinstance(x, Chromatogram))


def read(source: Source) -> Run:
    """Read a whole file into memory. Use :func:`iterparse` for large files."""
    info: dict = {}
    spectra: List[Spectrum] = []
    chroms: List[Chromatogram] = []
    for item in _stream(source, info):
        (spectra if isinstance(item, Spectrum) else chroms).append(item)
    return Run(spectra=spectra, chromatograms=chroms,
               indexed=info.get("indexed", False), run_id=info.get("run_id"))


# ------------------------------------------------------------------ command line
def summarize(path: Path) -> str:
    """A short, streaming summary of one file."""
    info: dict = {}
    levels: Dict[int, int] = {}
    modes = set()
    n = n_chrom = 0
    t_min = t_max = mz_min = mz_max = None
    t_unit = None
    chrom_lines = []
    for item in _stream(path, info):
        if isinstance(item, Chromatogram):
            n_chrom += 1
            chrom_lines.append(f"{item.id} ({item.kind or 'unknown type'}, {len(item.time)} points)")
            continue
        n += 1
        levels[item.ms_level or 0] = levels.get(item.ms_level or 0, 0) + 1
        modes.add({True: "centroid", False: "profile", None: "not stated"}[item.centroid])
        if item.scan_start_time is not None:
            t_unit = item.scan_start_time_unit
            t_min = item.scan_start_time if t_min is None else min(t_min, item.scan_start_time)
            t_max = item.scan_start_time if t_max is None else max(t_max, item.scan_start_time)
        if len(item.mz):
            mz_min = min(item.mz) if mz_min is None else min(mz_min, min(item.mz))
            mz_max = max(item.mz) if mz_max is None else max(mz_max, max(item.mz))
    by_level = ", ".join(f"MS{k}: {v}" for k, v in sorted(levels.items()))
    lines = [str(path),
             f"  format        : {'indexed mzML' if info.get('indexed') else 'mzML'} (run '{info.get('run_id')}')",
             f"  spectra       : {n} ({by_level}), {', '.join(sorted(modes)) or '-'}"]
    if t_min is not None:
        lines.append(f"  scan time     : {t_min:g} to {t_max:g} {t_unit or ''}".rstrip())
    if mz_min is not None:
        lines.append(f"  m/z           : {mz_min:.4f} to {mz_max:.4f}")
    lines.append(f"  chromatograms : {n_chrom}" + (f": {'; '.join(chrom_lines)}" if chrom_lines else ""))
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Summarize mzML files (streaming, standard library only).")
    ap.add_argument("files", nargs="+", type=Path)
    args = ap.parse_args(argv)
    status = 0
    for path in args.files:
        try:
            print(summarize(path))
        except (MzMLError, OSError) as exc:
            print(f"{path}\n  ERROR: {exc}")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
