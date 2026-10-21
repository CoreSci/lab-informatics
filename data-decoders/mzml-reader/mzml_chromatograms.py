"""Chromatograms from mzML spectra: TIC, BPC and XIC, a cross-check, CSV and SVG output.

Standard library only. Spectra are streamed with :mod:`mzml`, one at a time,
and every trace is built in a single pass over the file:

* **TIC** (total ion current): the sum of all intensities in each MS1 spectrum;
* **BPC** (base peak chromatogram): the largest intensity in each MS1 spectrum,
  with the m/z of that peak;
* **XIC** (extracted ion chromatogram): the sum of the intensities inside an
  m/z window, given as a centre +/- a tolerance in Da or in ppm (bounds
  included).

Only MS1 spectra contribute (MS2 and higher, and spectra without an MS level,
are skipped and counted). Times and their unit are kept exactly as stored.

The computed TIC is compared with what the file itself declares, when present:
its total ion current chromatogram and each spectrum's *total ion current*
value. Differences are reported; neither side is trusted silently.

:func:`svg_plot` draws one trace as a small, dependency-free SVG line plot.

Command line::

    python mzml_chromatograms.py FILE [--xic MZ[,TOL[,ppm|Da]]] ... [--csv OUT.csv]
                                      [--svg OUT.svg [--plot tic|bpc|xic]] [--rel-tol 1e-5]

prints a summary and the TIC cross-check. It exits with status 1 when the file
can't be read.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import io
import math
import sys
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union
from xml.sax.saxutils import escape

import mzml
import psi_ms_cv as cv

__all__ = ["Trace", "TicCheck", "Extraction", "mz_window", "extract", "tic", "bpc", "xic",
           "csv_text", "write_csv", "svg_plot", "main"]

_KINDS = {"TIC": cv.TIC_CHROMATOGRAM[1], "BPC": cv.BASEPEAK_CHROMATOGRAM[1],
          "XIC": cv.SELECTED_ION_CURRENT_CHROMATOGRAM[1]}


@dataclass
class Trace:
    """One chromatogram computed from MS1 spectra.

    ``kind`` is ``"TIC"``, ``"BPC"`` or ``"XIC"``. ``mz`` holds the base peak m/z
    of each point (BPC only; NaN for an empty spectrum). ``window`` is the XIC's
    (low, high) m/z bounds.
    """

    kind: str
    label: str
    time: array
    time_unit: Optional[str]
    intensity: array
    intensity_unit: Optional[str]
    spectrum_ids: List[str]
    mz: Optional[array] = None
    window: Optional[Tuple[float, float]] = None

    def apex(self) -> Optional[Tuple[float, float]]:
        """(time, intensity) of the highest point, or None for an empty trace."""
        if not len(self.intensity):
            return None
        i = max(range(len(self.intensity)), key=self.intensity.__getitem__)
        return self.time[i], self.intensity[i]

    def to_chromatogram(self, chrom_id: Optional[str] = None, index: int = 0) -> mzml.Chromatogram:
        """This trace as an :class:`mzml.Chromatogram`, e.g. to write it with :mod:`mzml_writer`."""
        return mzml.Chromatogram(id=chrom_id or self.kind, index=index, kind=_KINDS[self.kind],
                                 time=array("d", self.time), time_unit=self.time_unit,
                                 intensity=array("d", self.intensity), intensity_unit=self.intensity_unit)


@dataclass
class TicCheck:
    """Computed TIC compared with a value the file declares."""

    source: str
    compared: int = 0
    mismatches: int = 0
    max_rel_diff: float = 0.0
    notes: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.mismatches and not self.notes

    def __str__(self) -> str:
        head = (f"{self.source}: {self.compared} point(s) compared, "
                f"{self.mismatches} differ (max relative difference {self.max_rel_diff:.2g})")
        return "; ".join([head] + self.notes)


@dataclass
class Extraction:
    """Everything :func:`extract` computes from one file."""

    tic: Trace
    bpc: Trace
    xics: List[Trace]
    checks: List[TicCheck]
    ms1_spectra: int
    skipped_spectra: int


# ------------------------------------------------------------------ windows
def mz_window(center: float, tolerance: float, unit: str = "ppm") -> Tuple[float, float]:
    """(low, high) m/z bounds for ``center`` +/- ``tolerance`` in ``"ppm"`` or ``"Da"``."""
    if tolerance < 0:
        raise ValueError(f"tolerance must not be negative: {tolerance}")
    u = unit.lower()
    if u == "ppm":
        delta = center * tolerance / 1e6
    elif u == "da":
        delta = tolerance
    else:
        raise ValueError(f"unknown tolerance unit {unit!r} (use 'ppm' or 'Da')")
    return center - delta, center + delta


def _xic_label(center: float, tolerance: float, unit: str) -> str:
    return f"XIC m/z {center:.4f} +/- {tolerance:g} {'ppm' if unit.lower() == 'ppm' else 'Da'}"


# ------------------------------------------------------------------ extraction
XicSpec = Tuple[float, float, str]


def extract(source, xics: Iterable[XicSpec] = (), rel_tol: float = 1e-5) -> Extraction:
    """Compute TIC, BPC and one XIC per ``(center, tolerance, unit)`` in a single streaming pass.

    ``rel_tol`` is the relative difference above which a declared TIC value
    counts as different from the computed one.
    """
    specs = [(c, t, u, mz_window(c, t, u)) for c, t, u in xics]
    times, tics, bpcs, bpc_mz, ids = array("d"), array("d"), array("d"), array("d"), []
    xic_values = [array("d") for _ in specs]
    time_unit = intensity_unit = None
    skipped = 0
    declared: List[mzml.Chromatogram] = []
    spectrum_check = TicCheck("spectrum total ion current")
    for item in mzml.iterparse(source):
        if isinstance(item, mzml.Chromatogram):
            if item.kind == cv.TIC_CHROMATOGRAM[1]:
                declared.append(item)
            continue
        total = math.fsum(item.intensity)
        if item.total_ion_current is not None:
            _compare(spectrum_check, item.total_ion_current, total, rel_tol)
        if item.ms_level != 1:
            skipped += 1
            continue
        if item.scan_start_time is None:
            raise mzml.MzMLError(f"spectrum '{item.id}': no scan start time, so it can't be placed on a chromatogram")
        if not ids:
            time_unit, intensity_unit = item.scan_start_time_unit, item.intensity_unit
        elif item.scan_start_time_unit != time_unit:
            raise mzml.MzMLError(f"spectrum '{item.id}': scan start time in {item.scan_start_time_unit!r}, "
                                 f"earlier spectra in {time_unit!r}; mixed units are not converted")
        ids.append(item.id)
        times.append(item.scan_start_time)
        tics.append(total)
        if len(item.intensity):
            i = max(range(len(item.intensity)), key=item.intensity.__getitem__)
            bpcs.append(item.intensity[i])
            bpc_mz.append(item.mz[i])
        else:
            bpcs.append(0.0)
            bpc_mz.append(math.nan)
        for values, (_, _, _, (lo, hi)) in zip(xic_values, specs):
            values.append(math.fsum(v for m, v in zip(item.mz, item.intensity) if lo <= m <= hi))

    def trace(kind, label, values, **extra):
        return Trace(kind=kind, label=label, time=array("d", times), time_unit=time_unit, intensity=values,
                     intensity_unit=intensity_unit, spectrum_ids=list(ids), **extra)

    result = Extraction(
        tic=trace("TIC", "TIC (MS1)", tics),
        bpc=trace("BPC", "BPC (MS1)", bpcs, mz=bpc_mz),
        xics=[trace("XIC", _xic_label(c, t, u), v, window=w) for (c, t, u, w), v in zip(specs, xic_values)],
        checks=[], ms1_spectra=len(ids), skipped_spectra=skipped)
    for chrom in declared:
        result.checks.append(_check_declared(chrom, result.tic, rel_tol))
    if spectrum_check.compared:
        result.checks.append(spectrum_check)
    return result


def _compare(check: TicCheck, declared: float, computed: float, rel_tol: float) -> None:
    check.compared += 1
    scale = max(abs(declared), abs(computed))
    diff = abs(declared - computed) / scale if scale else 0.0
    check.max_rel_diff = max(check.max_rel_diff, diff)
    if diff > rel_tol:
        check.mismatches += 1


def _check_declared(chrom: mzml.Chromatogram, computed: Trace, rel_tol: float) -> TicCheck:
    check = TicCheck(f"declared TIC chromatogram '{chrom.id}'")
    if chrom.time_unit != computed.time_unit:
        check.notes.append(f"time units differ (declared {chrom.time_unit!r}, spectra {computed.time_unit!r}); "
                           "values not compared")
        return check
    order = sorted(range(len(computed.time)), key=computed.time.__getitem__)
    sorted_times = [computed.time[i] for i in order]
    matched, extra = set(), 0
    for t, value in zip(chrom.time, chrom.intensity):
        k = bisect.bisect_left(sorted_times, t)
        near = [j for j in (k - 1, k) if 0 <= j < len(sorted_times) and abs(sorted_times[j] - t) <= 1e-9 * max(1.0, abs(t))]
        if not near:
            extra += 1
            continue
        i = order[near[0]]
        matched.add(i)
        _compare(check, value, computed.intensity[i], rel_tol)
    if extra:
        check.notes.append(f"{extra} declared point(s) have no MS1 spectrum at that time "
                           "(the declared TIC may include MS2 or higher spectra)")
    missing = len(computed.time) - len(matched)
    if missing:
        check.notes.append(f"{missing} MS1 spectrum time(s) missing from the declared TIC")
    return check


def tic(source) -> Trace:
    """TIC of the MS1 spectra."""
    return extract(source).tic


def bpc(source) -> Trace:
    """Base peak chromatogram of the MS1 spectra (``.mz`` holds the base peak m/z)."""
    return extract(source).bpc


def xic(source, center: float, tolerance: float = 10.0, unit: str = "ppm") -> Trace:
    """Extracted ion chromatogram of ``center`` +/- ``tolerance`` (``"ppm"`` or ``"Da"``)."""
    return extract(source, [(center, tolerance, unit)]).xics[0]


# ------------------------------------------------------------------ CSV
def _cell(v: float) -> str:
    return "" if math.isnan(v) else repr(v)


def csv_text(result: Extraction) -> str:
    """All traces as CSV: one row per MS1 spectrum."""
    unit = result.tic.time_unit or "unknown unit"
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    w.writerow(["spectrum id", f"time ({unit})", "TIC", "BPC", "BPC m/z"]
               + [x.label for x in result.xics])
    for i, sid in enumerate(result.tic.spectrum_ids):
        w.writerow([sid, repr(result.tic.time[i]), _cell(result.tic.intensity[i]), _cell(result.bpc.intensity[i]),
                    _cell(result.bpc.mz[i])] + [_cell(x.intensity[i]) for x in result.xics])
    return out.getvalue()


def write_csv(result: Extraction, path: Union[str, Path]) -> None:
    Path(path).write_bytes(csv_text(result).encode("utf-8"))


# ------------------------------------------------------------------ SVG
def _ticks(lo: float, hi: float, n: int = 5) -> List[float]:
    """Round tick values (steps of 1, 2, 2.5 or 5 x 10^k) covering [lo, hi]."""
    if hi <= lo:
        hi = lo + 1.0
    raw = (hi - lo) / n
    mag = 10.0 ** math.floor(math.log10(raw))
    step = next(m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw)
    first = math.ceil(lo / step - 1e-9) * step
    count = int(math.floor((hi - first) / step + 1e-9)) + 1
    return [round(first + k * step, 12) for k in range(count)]


def _num(v: float) -> str:
    return f"{v:.2f}"


def _label(v: float) -> str:
    return f"{round(v, 9):g}"


def svg_plot(trace: Trace, width: int = 640, height: int = 360, title: Optional[str] = None) -> str:
    """A line plot of one trace as SVG text: labelled axes with units, no dependencies.

    The output depends only on the trace values, so the same trace always
    gives the same bytes.
    """
    left, right, top, bottom = 78, 20, 36, 52
    pw, ph = width - left - right, height - top - bottom
    t, y = list(trace.time), list(trace.intensity)
    t0, t1 = (min(t), max(t)) if t else (0.0, 1.0)
    if t1 <= t0:
        t1 = t0 + 1.0
    ymax = max(y) if y and max(y) > 0 else 1.0
    exp = int(math.floor(math.log10(ymax))) if ymax >= 1e4 else 0
    scale = 10.0 ** exp
    xt = [v for v in _ticks(t0, t1) if t0 - 1e-9 <= v <= t1 + 1e-9]
    yt = _ticks(0.0, ymax / scale)
    ytop = max(yt[-1], ymax / scale)

    def sx(v):
        return left + (v - t0) / (t1 - t0) * pw

    def sy(v):
        return top + ph - v / scale / ytop * ph

    unit = trace.time_unit or "unknown unit"
    yunit = f", {trace.intensity_unit}" if trace.intensity_unit else ""
    ylab = escape(f"Intensity{yunit}")
    if exp:
        ylab += f' (×10<tspan baseline-shift="super" font-size="9">{exp}</tspan>)'
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="11">',
        f'<rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="{left + pw / 2:.1f}" y="22" text-anchor="middle" font-size="13">{escape(title or trace.label).replace("+/-", "±")}</text>',
    ]
    for v in yt:
        yy = _num(sy(v * scale))
        lines.append(f'<line x1="{left}" y1="{yy}" x2="{left + pw}" y2="{yy}" stroke="#e5e5e5"/>')
        lines.append(f'<text x="{left - 6}" y="{yy}" text-anchor="end" dominant-baseline="middle">{_label(v)}</text>')
    for v in xt:
        xx = _num(sx(v))
        lines.append(f'<line x1="{xx}" y1="{top + ph}" x2="{xx}" y2="{top + ph + 4}" stroke="#333333"/>')
        lines.append(f'<text x="{xx}" y="{top + ph + 16}" text-anchor="middle">{_label(v)}</text>')
    lines += [
        f'<line x1="{left}" y1="{top + ph}" x2="{left + pw}" y2="{top + ph}" stroke="#333333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + ph}" stroke="#333333"/>',
        f'<text x="{left + pw / 2:.1f}" y="{height - 12}" text-anchor="middle">Time ({escape(unit)})</text>',
        f'<text transform="translate(16 {top + ph / 2:.1f}) rotate(-90)" text-anchor="middle">{ylab}</text>',
    ]
    if t:
        pts = " ".join(f"{_num(sx(a))},{_num(sy(b))}" for a, b in zip(t, y))
        lines.append(f'<polyline fill="none" stroke="#1f5fa8" stroke-width="1.5" points="{pts}"/>')
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------ command line
def _parse_xic(text: str) -> XicSpec:
    parts = [p.strip() for p in text.split(",")]
    try:
        center = float(parts[0])
        tol = float(parts[1]) if len(parts) > 1 else 10.0
    except ValueError:
        raise argparse.ArgumentTypeError(f"--xic expects MZ[,TOL[,ppm|Da]], got {text!r}") from None
    unit = parts[2] if len(parts) > 2 else "ppm"
    if unit.lower() not in ("ppm", "da") or len(parts) > 3:
        raise argparse.ArgumentTypeError(f"--xic expects MZ[,TOL[,ppm|Da]], got {text!r}")
    return center, tol, unit


def _describe(trace: Trace) -> str:
    peak = trace.apex()
    if peak is None:
        return "empty"
    if peak[1] == 0:
        return "all zero"
    unit = trace.time_unit or ""
    text = f"max {peak[1]:.6g} at {peak[0]:g} {unit}".rstrip()
    if trace.kind == "BPC":
        i = list(trace.time).index(peak[0])
        text += f" (m/z {trace.mz[i]:.4f})"
    if trace.window:
        text += f"  [window {trace.window[0]:.4f} to {trace.window[1]:.4f}]"
    return text


def summary(path, result: Extraction) -> str:
    t = result.tic
    span = f", {t.time[0]:g} to {t.time[-1]:g} {t.time_unit or ''}".rstrip() if len(t.time) else ""
    lines = [str(path),
             f"  MS1 spectra : {result.ms1_spectra}{span} ({result.skipped_spectra} other spectra skipped)",
             f"  TIC         : {_describe(result.tic)}",
             f"  BPC         : {_describe(result.bpc)}"]
    lines += [f"  {x.label}: {_describe(x)}" for x in result.xics]
    if not result.checks:
        lines.append("  TIC check   : nothing declared to compare with")
    for c in result.checks:
        lines.append(f"  TIC check   : {'OK' if c.ok else 'DIFFERENCE'}: {c}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="TIC, BPC and XIC from the MS1 spectra of an mzML file.")
    ap.add_argument("file", type=Path)
    ap.add_argument("--xic", action="append", type=_parse_xic, default=[], metavar="MZ[,TOL[,ppm|Da]]",
                    help="extracted ion chromatogram; tolerance defaults to 10 ppm (repeatable)")
    ap.add_argument("--csv", type=Path, help="write all traces to this CSV file")
    ap.add_argument("--svg", type=Path, help="write one trace as an SVG line plot")
    ap.add_argument("--plot", choices=("tic", "bpc", "xic"), help="trace for --svg (default: the first XIC, else the TIC)")
    ap.add_argument("--rel-tol", type=float, default=1e-5, help="relative tolerance of the TIC check (default 1e-5)")
    args = ap.parse_args(argv)
    if args.plot == "xic" and not args.xic:
        ap.error("--plot xic needs --xic")
    try:
        result = extract(args.file, args.xic, rel_tol=args.rel_tol)
    except (mzml.MzMLError, OSError) as exc:
        print(f"{args.file}\n  ERROR: {exc}")
        return 1
    print(summary(args.file, result))
    if args.csv:
        write_csv(result, args.csv)
        print(f"  wrote {args.csv}")
    if args.svg:
        choice = args.plot or ("xic" if result.xics else "tic")
        trace = {"tic": result.tic, "bpc": result.bpc}.get(choice) or result.xics[0]
        args.svg.write_bytes(svg_plot(trace).encode("utf-8"))
        print(f"  wrote {args.svg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
