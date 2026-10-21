"""Structural and integrity validator for mzML files, standard library only.

Checks what a reader relies on, in one streaming pass plus a few seeks for
indexed files: required sections, cross-references (CV, param groups,
instrument configurations, data processing), list counts, spectrum indices and
ids, decoded array lengths, spectrum and chromatogram content rules, and the
``indexedmzML`` byte offsets and SHA-1 checksum.

Each finding has a stable code (``MZV001`` ...), a level (``error`` or
``warning``) and a message naming the spectrum or chromatogram concerned.
See ``CHECKS`` for the full list.

Not covered: XSD schema validation, the full PSI semantic mapping rules (the
official PSI semantic validator does that), and MS-Numpress data, which can't
be decoded here and is reported as undecodable.

Command line::

    python mzml_validator.py FILE [FILE ...] [--strict]

prints a summary and every finding for each file. It exits with status 1 when
any file has errors. ``--strict`` turns warnings into errors.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union

import mzml
import psi_ms_cv as cv

__all__ = ["CHECKS", "Finding", "Report", "validate", "main"]

# code: (level, what it checks)
CHECKS = {
    "MZV001": ("error", "the file is well-formed XML with an <mzML> or <indexedmzML> root"),
    "MZV002": ("error", "required top-level elements are present"),
    "MZV003": ("error", "every cvRef / unitCvRef is declared in cvList"),
    "MZV004": ("error", "every referenceableParamGroupRef resolves"),
    "MZV005": ("error", "every instrumentConfigurationRef resolves"),
    "MZV006": ("error", "every dataProcessingRef resolves"),
    "MZV010": ("error", "spectrumList / chromatogramList count matches the number of elements"),
    "MZV011": ("error", "index values are consecutive from 0"),
    "MZV012": ("error", "spectrum and chromatogram ids are unique"),
    "MZV020": ("error", "each array's decoded length equals defaultArrayLength / arrayLength"),
    "MZV021": ("error", "each binary array can be decoded (type, data type, compression, base64, zlib)"),
    "MZV022": ("error", "encodedLength matches the base64 text"),
    "MZV023": ("error", "required arrays are present (m/z + intensity, or time + intensity)"),
    "MZV024": ("error", "m/z and intensity arrays have the same length"),
    "MZV030": ("warning", "a spectrum states centroid or profile"),
    "MZV031": ("error", "the m/z array is non-decreasing"),
    "MZV032": ("warning", "intensities are non-negative"),
    "MZV033": ("error", "a spectrum states its MS level"),
    "MZV034": ("error", "an MS2 (or higher) spectrum has a precursor"),
    "MZV040": ("error", "a chromatogram's time and intensity arrays have the same length"),
    "MZV041": ("error", "a chromatogram's time array is non-decreasing"),
    "MZV050": ("error", "each index offset points at the start of the element with that id"),
    "MZV051": ("error", "indexListOffset points at <indexList>"),
    "MZV052": ("error", "the SHA-1 fileChecksum matches the file"),
    "MZV053": ("warning", "every spectrum and chromatogram is in the index"),
}

REQUIRED_SECTIONS = ("cvList", "fileDescription", "softwareList", "instrumentConfigurationList",
                     "dataProcessingList", "run")


@dataclass
class Finding:
    """A validation finding. ``level`` is ``"error"`` or ``"warning"``."""

    code: str
    level: str
    message: str

    def __str__(self) -> str:
        return f"{self.level.upper()}: {self.code}: {self.message}"


@dataclass
class Report:
    """Validation result for one file."""

    path: str
    findings: List[Finding] = field(default_factory=list)
    indexed: bool = False
    spectra: int = 0
    chromatograms: int = 0

    @property
    def errors(self) -> List[Finding]:
        return [f for f in self.findings if f.level == "error"]

    @property
    def warnings(self) -> List[Finding]:
        return [f for f in self.findings if f.level == "warning"]

    @property
    def ok(self) -> bool:
        """True when validation found no errors (warnings are allowed)."""
        return not self.errors

    @property
    def codes(self) -> List[str]:
        return [f.code for f in self.findings]


class _Validator:
    def __init__(self, path: Path):
        self.path = path
        self.report = Report(str(path))
        self.cvs: set = set()
        self.cv_used: Dict[str, str] = {}             # cvRef -> first place it's used
        self.groups: Dict[str, list] = {}
        self.group_used: Dict[str, str] = {}
        self.instruments: set = set()
        self.instrument_used: Dict[str, str] = {}
        self.processing: set = set()
        self.processing_used: Dict[str, str] = {}
        self.sections: set = set()
        self.ids: Dict[str, set] = {"spectrum": set(), "chromatogram": set()}
        self.next_index = {"spectrum": 0, "chromatogram": 0}
        self.counted = {"spectrum": 0, "chromatogram": 0}
        self.declared: Dict[str, Optional[str]] = {}
        self.offsets: List[tuple] = []                 # (index name, idRef, offset text)
        self.index_list_offset: Optional[str] = None
        self.checksum: Optional[str] = None
        self.where = "file header"

    def add(self, code: str, message: str) -> None:
        self.report.findings.append(Finding(code, CHECKS[code][0], message))

    # ---------------------------------------------------------------- streaming pass
    def run(self) -> Report:
        stack: List[ET.Element] = []
        index_name = None
        try:
            with open(self.path, "rb") as fh:
                for event, elem in ET.iterparse(fh, events=("start", "end")):
                    tag = mzml._local(elem.tag)
                    if event == "start":
                        if not stack:
                            if tag not in ("mzML", "indexedmzML"):
                                self.add("MZV001", f"not an mzML document: the root element is <{tag}>")
                                return self.report
                            self.report.indexed = tag == "indexedmzML"
                        elif mzml._local(stack[-1].tag) == "mzML":
                            self.sections.add(tag)
                        if tag == "index":
                            index_name = elem.get("name")
                        self._on_start(tag, elem)
                        stack.append(elem)
                        continue
                    stack.pop()
                    if self._on_end(tag, elem, index_name):
                        elem.clear()
                        if stack:
                            stack[-1].remove(elem)
        except ET.ParseError as exc:
            self.add("MZV001", f"malformed XML: {exc}")
            return self.report
        self._check_references()
        if self.report.indexed:
            self._check_index()
        return self.report

    def _on_start(self, tag: str, elem: ET.Element) -> None:
        a = elem.attrib
        if tag in ("spectrum", "chromatogram"):
            self.where = f"{tag} '{a.get('id', '?')}'"
        elif tag in ("spectrumList", "chromatogramList"):
            self.declared[tag] = a.get("count")
            self._use(self.processing_used, a.get("defaultDataProcessingRef"), f"<{tag}>")
        elif tag == "run":
            self._use(self.instrument_used, a.get("defaultInstrumentConfigurationRef"), "<run>")
        elif tag == "scan":
            self._use(self.instrument_used, a.get("instrumentConfigurationRef"), self.where)
        elif tag == "cvParam":
            self._use(self.cv_used, a.get("cvRef"), self.where)
            self._use(self.cv_used, a.get("unitCvRef"), self.where)
        elif tag == "referenceableParamGroupRef":
            self._use(self.group_used, a.get("ref"), self.where)
        elif tag == "cv":
            self.cvs.add(a.get("id"))
        elif tag == "instrumentConfiguration":
            self.instruments.add(a.get("id"))
        elif tag == "dataProcessing":
            self.processing.add(a.get("id"))
        if tag in ("spectrum", "chromatogram", "binaryDataArray"):
            self._use(self.processing_used, a.get("dataProcessingRef"), self.where)

    def _on_end(self, tag: str, elem: ET.Element, index_name: Optional[str]) -> bool:
        """Handle a finished element. Returns True when it can be freed."""
        if tag == "referenceableParamGroup":
            self.groups[elem.get("id", "")] = self._params(elem)
        elif tag in ("spectrumList", "chromatogramList"):
            kind = tag[:-4]
            declared, real = self.declared.get(tag), self.counted[kind]
            if declared is None or not declared.isdigit() or int(declared) != real:
                plural = "spectra" if kind == "spectrum" else "chromatograms"
                self.add("MZV010", f"<{tag}> count is {declared!r}, but it holds {real} {plural}")
        elif tag == "spectrum":
            self._spectrum(elem)
            self.where = "run"
            return True
        elif tag == "chromatogram":
            self._chromatogram(elem)
            self.where = "run"
            return True
        elif tag == "offset":
            self.offsets.append((index_name, elem.get("idRef"), (elem.text or "").strip()))
        elif tag == "indexListOffset":
            self.index_list_offset = (elem.text or "").strip()
        elif tag == "fileChecksum":
            self.checksum = (elem.text or "").strip().lower()
        return False

    @staticmethod
    def _use(table: Dict[str, str], ref: Optional[str], where: str) -> None:
        if ref is not None and ref not in table:
            table[ref] = where

    def _params(self, elem: Optional[ET.Element]) -> List[dict]:
        """cvParams with param groups expanded; unknown groups are reported separately (MZV004)."""
        out: List[dict] = []
        for c in (elem if elem is not None else []):
            t = mzml._local(c.tag)
            if t == "cvParam":
                out.append(dict(c.attrib))
            elif t == "referenceableParamGroupRef":
                out.extend(self.groups.get(c.get("ref"), []))
        return out

    # ---------------------------------------------------------------- items
    def _identity(self, kind: str, elem: ET.Element) -> str:
        """Index and id checks shared by spectra and chromatograms. Returns the label for messages."""
        self.counted[kind] += 1
        self.report.spectra += kind == "spectrum"
        self.report.chromatograms += kind == "chromatogram"
        item_id = elem.get("id")
        label = f"{kind} '{item_id}'"
        if item_id is None:
            label = f"{kind} #{self.counted[kind]} (no id)"
        elif item_id in self.ids[kind]:
            self.add("MZV012", f"{label}: duplicate id")
        else:
            self.ids[kind].add(item_id)
        expected, index = self.next_index[kind], elem.get("index")
        if index != str(expected):
            self.add("MZV011", f"{label}: index is {index!r}, expected {expected}")
        self.next_index[kind] += 1
        return label

    def _arrays(self, elem: ET.Element, label: str, allowed: tuple) -> Dict[str, list]:
        """Decode every binary array, reporting MZV020-MZV023. Returns {array name: values}."""
        out: Dict[str, list] = {}
        default = elem.get("defaultArrayLength")
        if default is None or not default.isdigit():
            self.add("MZV020", f"{label}: defaultArrayLength is {default!r}, not a whole number")
            default = None
        bdal = mzml._child(elem, "binaryDataArrayList")
        for i, bda in enumerate(mzml._children(bdal, "binaryDataArray") if bdal is not None else []):
            try:
                name, bits, compressed, _ = mzml._array_spec(self._params(bda), f"{label}, binary array {i + 1}")
            except mzml.MzMLError as exc:
                self.add("MZV021", str(exc))
                continue
            here = f"{label}, {name}"
            if name not in allowed or name in out:
                self.add("MZV021", f"{here}: " + ("appears twice" if name in out else "not expected here"))
                continue
            binary = mzml._child(bda, "binary")
            text = "".join((binary.text or "").split()) if binary is not None else ""
            enc = bda.get("encodedLength")
            if enc is not None and (not enc.isdigit() or int(enc) != len(text)):
                self.add("MZV022", f"{here}: encodedLength is {enc!r}, but the base64 text has {len(text)} characters")
            try:
                values = mzml.decode_binary(text, bits, compressed)
            except mzml.MzMLError as exc:
                self.add("MZV021", f"{here}: {exc}")
                continue
            declared = bda.get("arrayLength", default)
            if declared is not None and (not declared.isdigit() or int(declared) != len(values)):
                which = "arrayLength" if bda.get("arrayLength") is not None else "defaultArrayLength"
                self.add("MZV020", f"{here}: decoded {len(values)} values, but {which} is {declared}")
            out[name] = values
        for name in allowed:
            if name not in out and not any(f.code == "MZV021" and f.message.startswith(label) for f in self.report.findings):
                self.add("MZV023", f"{label}: no {name}")
        return out

    def _spectrum(self, elem: ET.Element) -> None:
        label = self._identity("spectrum", elem)
        params = self._params(elem)
        accs = {p.get("accession") for p in params}
        level = next((p.get("value") for p in params if p.get("accession") == cv.MS_LEVEL[0]), None)
        if level is None or not level.strip().isdigit():
            self.add("MZV033", f"{label}: no MS level" if level is None else f"{label}: MS level is {level!r}")
        elif int(level) >= 2:
            plist = mzml._child(elem, "precursorList")
            if plist is None or not mzml._children(plist, "precursor"):
                self.add("MZV034", f"{label}: MS level {level} but no precursor")
        if not accs & {cv.CENTROID_SPECTRUM[0], cv.PROFILE_SPECTRUM[0]}:
            self.add("MZV030", f"{label}: neither centroid spectrum nor profile spectrum is stated")
        arrays = self._arrays(elem, label, (cv.MZ_ARRAY[1], cv.INTENSITY_ARRAY[1]))
        mz, inten = arrays.get(cv.MZ_ARRAY[1]), arrays.get(cv.INTENSITY_ARRAY[1])
        if mz is not None and inten is not None and len(mz) != len(inten):
            self.add("MZV024", f"{label}: {len(mz)} m/z values but {len(inten)} intensities")
        if mz is not None:
            i = _first_decrease(mz)
            if i is not None:
                self.add("MZV031", f"{label}: m/z decreases at point {i} ({mz[i]:.5f} after {mz[i - 1]:.5f})")
        if inten is not None:
            negative = sum(1 for v in inten if v < 0)
            if negative:
                self.add("MZV032", f"{label}: {negative} negative intensity value(s)")

    def _chromatogram(self, elem: ET.Element) -> None:
        label = self._identity("chromatogram", elem)
        arrays = self._arrays(elem, label, (cv.TIME_ARRAY[1], cv.INTENSITY_ARRAY[1]))
        time, inten = arrays.get(cv.TIME_ARRAY[1]), arrays.get(cv.INTENSITY_ARRAY[1])
        if time is not None and inten is not None and len(time) != len(inten):
            self.add("MZV040", f"{label}: {len(time)} time values but {len(inten)} intensities")
        if time is not None:
            i = _first_decrease(time)
            if i is not None:
                self.add("MZV041", f"{label}: time decreases at point {i} ({time[i]:g} after {time[i - 1]:g})")

    # ---------------------------------------------------------------- file-level checks
    def _check_references(self) -> None:
        for name in REQUIRED_SECTIONS:
            if name not in self.sections:
                self.add("MZV002", f"missing required element <{name}>")
        for ref, where in self.cv_used.items():
            if ref not in self.cvs:
                self.add("MZV003", f"{where}: cvRef '{ref}' is not declared in <cvList>")
        for ref, where in self.group_used.items():
            if ref not in self.groups:
                self.add("MZV004", f"{where}: referenceableParamGroupRef '{ref}' does not resolve")
        for ref, where in self.instrument_used.items():
            if ref not in self.instruments:
                self.add("MZV005", f"{where}: instrumentConfigurationRef '{ref}' does not resolve")
        for ref, where in self.processing_used.items():
            if ref not in self.processing:
                self.add("MZV006", f"{where}: dataProcessingRef '{ref}' does not resolve")

    def _check_index(self) -> None:
        size = self.path.stat().st_size
        indexed = {"spectrum": set(), "chromatogram": set()}
        with open(self.path, "rb") as fh:
            for name, ref, text in self.offsets:
                label = f"{name} '{ref}'"
                if name in indexed:
                    indexed[name].add(ref)
                if not text.isdigit() or int(text) >= size:
                    self.add("MZV050", f"{label}: index offset {text!r} is not a position in the file")
                    continue
                fh.seek(int(text))
                head = fh.read(512)
                m = re.match(rb'<(spectrum|chromatogram)\s[^>]*?\bid="([^"]*)"', head)
                # the id attribute is raw XML text: undo entity escapes (&quot; &amp; ...) before comparing
                if not m or m.group(1).decode() != name or html.unescape(m.group(2).decode("utf-8")) != ref:
                    self.add("MZV050", f"{label}: index offset {text} does not point at <{name} ... id=\"{ref}\">")
            if self.index_list_offset is not None:
                ok = self.index_list_offset.isdigit()
                if ok:
                    fh.seek(int(self.index_list_offset))
                    ok = fh.read(11) == b"<indexList "
                if not ok:
                    self.add("MZV051", f"indexListOffset {self.index_list_offset!r} does not point at <indexList>")
            if self.checksum is not None:
                actual = _sha1_through_checksum_tag(fh, size)
                if actual is None or actual != self.checksum:
                    self.add("MZV052", f"fileChecksum is {self.checksum}, but the file's SHA-1 is {actual}")
        for kind in ("spectrum", "chromatogram"):
            missing = sorted(self.ids[kind] - indexed[kind])
            if missing:
                shown = ", ".join(f"'{m}'" for m in missing[:5]) + (" ..." if len(missing) > 5 else "")
                self.add("MZV053", f"{len(missing)} {kind}(s) not in the index: {shown}")


def _first_decrease(values) -> Optional[int]:
    for i in range(1, len(values)):
        if values[i] < values[i - 1]:
            return i
    return None


def _sha1_through_checksum_tag(fh, size: int) -> Optional[str]:
    """SHA-1 from the start of the file to the end of the '<fileChecksum>' opening tag."""
    tag = b"<fileChecksum>"
    tail = min(size, 65536)
    fh.seek(size - tail)
    pos = fh.read(tail).rfind(tag)
    if pos < 0:
        return None
    end = size - tail + pos + len(tag)
    digest = hashlib.sha1()
    fh.seek(0)
    remaining = end
    while remaining:
        chunk = fh.read(min(remaining, 1 << 20))
        if not chunk:
            break
        digest.update(chunk)
        remaining -= len(chunk)
    return digest.hexdigest()


def validate(path: Union[str, Path], strict: bool = False) -> Report:
    """Validate one mzML file. With ``strict=True``, warnings are reported as errors."""
    report = _Validator(Path(path)).run()
    if strict:
        for f in report.findings:
            f.level = "error"
    return report


def _summary(report: Report) -> str:
    kind = "indexed mzML" if report.indexed else "mzML"
    lines = [report.path,
             f"  format     : {kind}, {report.spectra} spectra, {report.chromatograms} chromatogram(s)",
             "  validation : " + ("OK" if report.ok else f"{len(report.errors)} error(s)")
             + (f", {len(report.warnings)} warning(s)" if report.warnings else "")]
    lines += [f"    {f}" for f in report.findings]
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Validate mzML files (structure, references, arrays, index).")
    ap.add_argument("files", nargs="+", help="mzML files (.mzML)")
    ap.add_argument("--strict", action="store_true", help="treat warnings as errors")
    args = ap.parse_args(argv)
    status = 0
    for path in args.files:
        try:
            report = validate(path, strict=args.strict)
        except OSError as exc:
            print(f"{path}\n  ERROR: {exc}")
            status = 1
            continue
        print(_summary(report))
        status |= 0 if report.ok else 1
    return status


if __name__ == "__main__":
    sys.exit(main())
