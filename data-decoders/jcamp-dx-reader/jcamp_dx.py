"""JCAMP-DX reader, validator and writer for single-block 1-D spectra.

JCAMP-DX is the open IUPAC exchange format for spectra: a plain-text file of
``##LABEL=value`` records followed by a data table. This module implements the
parts of the published specification (McDonald & Wilks, *Applied Spectroscopy*
42(1), 1988; JCAMP-DX 4.24 conventions) needed for 1-D spectra such as IR,
UV-Vis and Raman:

* data tables ``##XYDATA=(X++(Y..Y))``, ``##XYPOINTS=(XY..XY)`` and
  ``##PEAK TABLE=(XY..XY)`` / ``(XYW..XYW)``;
* ordinate encodings AFFN (free-format numbers) and the compressed ASDF forms
  SQZ, DIF and DUP, including DIF "Y-checks" between lines;
* validation: required labels, X-checks (line abscissa vs. FIRSTX/LASTX/NPOINTS),
  Y-checks, point counts, FIRSTY, and ``##END``;
* a writer (AFFN or DIF/DUP compression) used to produce test data.

Not supported (rejected with a clear :class:`JcampError`): compound/linked
files (``##BLOCKS``) and ``##NTUPLES`` (e.g. multi-dimensional NMR).

Command line::

    python jcamp_dx.py FILE [FILE ...] [--csv OUT.csv] [--plot OUT.png] [--strict]

prints a summary and every validation issue for each file. It exits with
status 1 when any file has errors.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["JcampError", "Issue", "Spectrum", "normalize_label", "decode_line", "parse", "read", "write"]


class JcampError(ValueError):
    """Unreadable or unsupported JCAMP-DX content."""


@dataclass
class Issue:
    """A validation finding. ``level`` is ``"error"`` or ``"warning"``."""

    level: str
    message: str
    line: int | None = None

    def __str__(self) -> str:
        where = f"line {self.line}: " if self.line else ""
        return f"{self.level.upper()}: {where}{self.message}"


@dataclass
class Spectrum:
    """A parsed spectrum: labelled header values, x/y arrays and validation issues."""

    headers: dict[str, str]
    x: np.ndarray
    y: np.ndarray
    table: str                      # "XYDATA", "XYPOINTS" or "PEAKTABLE"
    issues: list[Issue] = field(default_factory=list)

    def header(self, label: str, default: str | None = None) -> str | None:
        """Header value by label, compared the JCAMP way (case, spaces, -, /, _ ignored)."""
        return self.headers.get(normalize_label(label), default)

    @property
    def title(self) -> str:
        return self.header("TITLE", "") or ""

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "error"]

    @property
    def ok(self) -> bool:
        """True when validation found no errors (warnings are allowed)."""
        return not self.errors


def normalize_label(label: str) -> str:
    """Canonical label form: JCAMP-DX treats ``##XY DATA``, ``##xydata`` and ``##XY_DATA`` alike."""
    return re.sub(r"[\s\-/_]", "", label).upper()


# ----------------------------------------------------------------------------- ASDF decoding
# Pseudo-digits replace the first digit of a number and encode its sign/meaning.
_SQZ = {"@": 0, **{c: i + 1 for i, c in enumerate("ABCDEFGHI")}, **{c: -(i + 1) for i, c in enumerate("abcdefghi")}}
_DIF = {"%": 0, **{c: i + 1 for i, c in enumerate("JKLMNOPQR")}, **{c: -(i + 1) for i, c in enumerate("jklmnopqr")}}
_DUP = {**{c: i + 1 for i, c in enumerate("STUVWXYZ")}, "s": 9}

# Exponents are only recognised with an explicit sign (1.5E+03) because a bare
# 'E'/'e' after a digit is a valid SQZ pseudo-digit (+5 / -5).
_TOKEN = re.compile(r"""
    (?P<affn>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]\d+)?)
  | (?P<sqz>[@A-Ia-i]\d*(?:\.\d*)?)
  | (?P<dif>[%J-Rj-r]\d*(?:\.\d*)?)
  | (?P<dup>[S-Zs]\d*)
  | (?P<missing>\?)
  | (?P<sep>[\s,;]+)
""", re.VERBOSE)


def _pseudo_value(text: str, table: dict[str, int]) -> float:
    """'A23' -> 123, 'b5' -> -25, '@' -> 0 (SQZ); the same rule applies to DIF letters."""
    d = table[text[0]]
    magnitude = float(str(abs(d)) + text[1:]) if text[1:] else float(abs(d))
    return -magnitude if d < 0 else magnitude


def _tokens(line: str) -> list[tuple[str, str]]:
    out, pos = [], 0
    for m in _TOKEN.finditer(line):
        if m.start() != pos:
            raise JcampError(f"invalid character {line[pos]!r} at column {pos + 1}")
        pos = m.end()
        if m.lastgroup != "sep":
            out.append((m.lastgroup, m.group()))
    if pos != len(line):
        raise JcampError(f"invalid character {line[pos]!r} at column {pos + 1}")
    return out


def decode_line(line: str, carry: float | None = None) -> tuple[float, list[float], str | None, bool]:
    """Decode one ``(X++(Y..Y))`` data line.

    Returns ``(x, ordinates, first_kind, ends_in_dif)``. ``first_kind`` is ``"abs"``
    or ``"dif"`` for the first ordinate token (``None`` if the line has no ordinates).
    ``carry`` is the last ordinate of the previous line, used when a line starts
    with a DIF token.
    """
    toks = _tokens(line.strip())
    if not toks:
        raise JcampError("empty data line")
    kind, text = toks[0]
    if kind == "affn":
        x = float(text)
    elif kind == "sqz":
        x = _pseudo_value(text, _SQZ)
    else:
        raise JcampError(f"data line must start with an abscissa value, found {text!r}")

    values: list[float] = []
    last = carry
    prev_kind: str | None = None      # "abs" or "dif": what a DUP repeats
    prev_delta = 0.0
    first_kind: str | None = None
    for kind, text in toks[1:]:
        if kind in ("affn", "sqz", "missing"):
            v = float(text) if kind == "affn" else (math.nan if kind == "missing" else _pseudo_value(text, _SQZ))
            values.append(v)
            last, prev_kind = v, "abs"
        elif kind == "dif":
            if last is None:
                raise JcampError("DIF value without a preceding ordinate")
            prev_delta = _pseudo_value(text, _DIF)
            last = last + prev_delta
            values.append(last)
            prev_kind = "dif"
        elif kind == "dup":
            if prev_kind is None:
                raise JcampError("DUP count without a preceding value")
            count = int(str(_DUP[text[0]]) + text[1:])
            for _ in range(count - 1):          # the count includes the original value
                if prev_kind == "dif":
                    last = last + prev_delta
                values.append(last)
            continue
        first_kind = first_kind or prev_kind
    return x, values, first_kind, prev_kind == "dif"


# ----------------------------------------------------------------------------- parsing
_XYDATA_VARS = "(X++(Y..Y))"
_POINT_VARS = {"(XY..XY)": 2, "(XYW..XYW)": 3, "(XYM..XYM)": 3}


def _number(headers: dict[str, str], label: str, issues: list[Issue], required: bool = True) -> float | None:
    raw = headers.get(label)
    if raw is None:
        if required:
            issues.append(Issue("error", f"missing required label ##{label}"))
        return None
    try:
        return float(raw.split()[0])
    except (ValueError, IndexError):
        issues.append(Issue("error", f"##{label} is not a number: {raw!r}"))
        return None


def _parse_xydata(lines: list[tuple[int, str]], headers: dict[str, str], issues: list[Issue]) -> tuple[np.ndarray, np.ndarray]:
    for label in ("XUNITS", "YUNITS"):
        if label not in headers:
            issues.append(Issue("error", f"missing required label ##{label}"))
    firstx = _number(headers, "FIRSTX", issues)
    lastx = _number(headers, "LASTX", issues)
    npoints = _number(headers, "NPOINTS", issues)
    xfactor = _number(headers, "XFACTOR", issues)
    yfactor = _number(headers, "YFACTOR", issues)
    if None in (firstx, lastx, npoints):
        raise JcampError("cannot build the abscissa without FIRSTX, LASTX and NPOINTS")
    npoints = int(npoints)
    xfactor = xfactor if xfactor is not None else 1.0
    yfactor = yfactor if yfactor is not None else 1.0
    ys: list[float] = []
    line_starts: list[tuple[int, float, int]] = []   # (line number, abscissa, index of its first ordinate)
    prev_ended_dif = False
    last_y: float | None = None
    for lineno, text in lines:
        if not text.strip():
            continue
        try:
            x_line, vals, first_kind, ends_dif = decode_line(text, carry=last_y)
        except JcampError as exc:
            issues.append(Issue("error", str(exc), lineno))
            continue
        index = len(ys)
        if prev_ended_dif:
            if first_kind == "abs":
                check = vals.pop(0)                      # DIF Y-check: repeats the previous ordinate
                if last_y is not None and not math.isclose(check, last_y, rel_tol=0, abs_tol=1e-9 * max(1.0, abs(last_y))):
                    issues.append(Issue("error", f"Y-check failed: line starts with {check:g}, previous line ended with {last_y:g}", lineno))
                index = len(ys) - 1
            else:
                issues.append(Issue("warning", "DIF line without a Y-check value", lineno))
        line_starts.append((lineno, x_line * xfactor, index))
        ys.extend(vals)
        if vals:
            last_y = vals[-1]
        prev_ended_dif = ends_dif

    # X-checks run after decoding: if NPOINTS disagrees with the data, the spacing is
    # derived from the decoded count so that one wrong header doesn't cascade into an
    # X-check error on every line.
    count = len(ys)
    if count != npoints:
        issues.append(Issue("error", f"NPOINTS={npoints} but the table holds {count} ordinates "
                                     "(X-checks use the decoded count)"))
    n_grid = count if count > 1 else npoints
    dx = (lastx - firstx) / (n_grid - 1) if n_grid > 1 else 0.0
    for lineno, x_real, index in line_starts:
        expected = firstx + index * dx
        tol = max(abs(dx) * 0.5, 1e-9 * max(1.0, abs(expected)))
        if abs(x_real - expected) > tol:
            issues.append(Issue("error", f"X-check failed: line abscissa {x_real:g}, expected {expected:g} for point {index + 1}", lineno))
    issues.sort(key=lambda i: (i.line is None, i.line or 0))
    y = np.asarray(ys, dtype=float) * yfactor
    x = firstx + np.arange(len(y)) * dx
    firsty = _number(headers, "FIRSTY", issues, required=False)
    if firsty is None and "FIRSTY" not in headers:
        issues.append(Issue("warning", "##FIRSTY is missing (recommended as a check on YFACTOR)"))
    elif firsty is not None and len(y) and not math.isclose(y[0], firsty, abs_tol=abs(yfactor) + 1e-12):
        issues.append(Issue("error", f"FIRSTY={firsty:g} but the first ordinate decodes to {y[0]:g}"))
    return x, y


def _parse_points(lines: list[tuple[int, str]], headers: dict[str, str], table: str, issues: list[Issue]) -> tuple[np.ndarray, np.ndarray]:
    group = _POINT_VARS[headers[table].replace(" ", "").upper()]
    numbers: list[float] = []
    for lineno, text in lines:
        for tok in (t for t in re.split(r"[\s,;]+", text.strip()) if t):
            try:
                numbers.append(math.nan if tok == "?" else float(tok))
            except ValueError:
                issues.append(Issue("error", f"not a number in {table}: {tok!r}", lineno))
    if len(numbers) % group:
        issues.append(Issue("error", f"{table} has {len(numbers)} values, not a multiple of {group}"))
        numbers = numbers[: len(numbers) - len(numbers) % group]
    arr = np.asarray(numbers, dtype=float).reshape(-1, group) if numbers else np.empty((0, group))
    xfactor = _number(headers, "XFACTOR", issues, required=False) or 1.0
    yfactor = _number(headers, "YFACTOR", issues, required=False) or 1.0
    npoints = _number(headers, "NPOINTS", issues, required=False)
    if npoints is not None and int(npoints) != len(arr):
        issues.append(Issue("error", f"NPOINTS={int(npoints)} but {table} holds {len(arr)} points"))
    return arr[:, 0] * xfactor, arr[:, 1] * yfactor


def parse(text: str, strict: bool = False) -> Spectrum:
    """Parse JCAMP-DX text. With ``strict=True``, raise :class:`JcampError` on any validation error."""
    issues: list[Issue] = []
    headers: dict[str, str] = {}
    data_lines: list[tuple[int, str]] = []
    table: str | None = None
    current: str | None = None
    seen_end = False
    first_label: str | None = None

    for lineno, raw in enumerate(text.splitlines(), 1):
        cut = raw.find("$$")                     # $$ starts a comment anywhere on a line
        line = raw if cut < 0 else raw[:cut]
        if seen_end:
            if line.strip():
                if normalize_label(line[2:].partition("=")[0]) == "TITLE":
                    raise JcampError("multiple blocks after ##END are not supported (single-block files only)")
                issues.append(Issue("warning", "content after ##END ignored", lineno))
                break
            continue
        if line.startswith("##"):
            label, _, value = line[2:].partition("=")
            key = normalize_label(label)
            first_label = first_label or key
            if key == "END":
                seen_end = True
                current = None
                continue
            if key == "BLOCKS":
                raise JcampError("compound (##BLOCKS / LINK) files are not supported")
            if key == "NTUPLES":
                raise JcampError("##NTUPLES files (e.g. NMR) are not supported")
            if key in ("XYDATA", "XYPOINTS", "PEAKTABLE"):
                if table is not None:
                    issues.append(Issue("error", f"second data table ##{key} ignored (one table per block)", lineno))
                    current = None
                    continue
                table, current = key, "__data__"
                headers[key] = value.strip()
                continue
            if key in headers:
                issues.append(Issue("warning", f"duplicate label ##{key}; the last value is used", lineno))
            headers[key] = value.strip()
            current = key
        elif current == "__data__":
            data_lines.append((lineno, line))
        elif current and line.strip():
            headers[current] = (headers[current] + "\n" + line.strip()).strip()
        elif line.strip():
            issues.append(Issue("warning", "text outside any labelled record ignored", lineno))

    # ---- header validation
    if not seen_end:
        issues.append(Issue("warning", "missing ##END= (file may be truncated)"))
    if first_label != "TITLE":
        issues.append(Issue("warning", "##TITLE should be the first label"))
    for label in ("TITLE", "JCAMPDX", "DATATYPE"):
        if label not in headers:
            issues.append(Issue("error", f"missing required label ##{'JCAMP-DX' if label == 'JCAMPDX' else label}"))
    for label in ("ORIGIN", "OWNER"):
        if label not in headers:
            issues.append(Issue("warning", f"missing ##{label} (required since JCAMP-DX 4.24)"))
    version = headers.get("JCAMPDX", "")
    if version and not re.match(r"\s*[45]\.", version):
        issues.append(Issue("warning", f"untested JCAMP-DX version {version!r} (4.x/5.x expected)"))
    if table is None:
        raise JcampError("no data table found (##XYDATA, ##XYPOINTS or ##PEAK TABLE)")

    variables = headers[table].replace(" ", "").upper()
    if table == "XYDATA":
        if variables != _XYDATA_VARS:
            raise JcampError(f"unsupported ##XYDATA variable list {headers[table]!r}; only {_XYDATA_VARS} is supported")
        x, y = _parse_xydata(data_lines, headers, issues)
    else:
        if variables not in _POINT_VARS:
            raise JcampError(f"unsupported ##{table} variable list {headers[table]!r}")
        x, y = _parse_points(data_lines, headers, table, issues)

    spectrum = Spectrum(headers=headers, x=x, y=y, table=table, issues=issues)
    if strict and spectrum.errors:
        raise JcampError("; ".join(str(e) for e in spectrum.errors))
    return spectrum


def read(path: str | Path, strict: bool = False) -> Spectrum:
    """Parse a JCAMP-DX file (ASCII/Latin-1 text)."""
    return parse(Path(path).read_text(encoding="latin-1"), strict=strict)


# ----------------------------------------------------------------------------- writing
def _pseudo(v: int, pos: str, neg: str, zero: str) -> str:
    if v == 0:
        return zero
    s = str(abs(v))
    return chr(ord(pos if v > 0 else neg) + int(s[0]) - 1) + s[1:]


def _sqz(v: int) -> str:
    return _pseudo(v, "A", "a", "@")


def _dif(v: int) -> str:
    return _pseudo(v, "J", "j", "%")


def _dup(n: int) -> str:
    s = str(n)
    return ("s" if s[0] == "9" else chr(ord("S") + int(s[0]) - 1)) + s[1:]


def _fmt(v: float) -> str:
    return f"{v:.10g}"


def write(x: np.ndarray, y: np.ndarray, *, title: str, data_type: str = "INFRARED SPECTRUM",
          xunits: str = "1/CM", yunits: str = "ABSORBANCE", form: str = "DIFDUP",
          yfactor: float | None = None, origin: str = "synthetic data", owner: str = "public domain",
          extra: dict[str, str] | None = None, width: int = 80) -> str:
    """Return JCAMP-DX text for evenly spaced ``x`` / ``y`` as ``##XYDATA=(X++(Y..Y))``.

    ``form`` is ``"AFFN"`` (plain integers) or ``"DIFDUP"`` (DIF + DUP compression with
    Y-checks). Ordinates are stored as integers ``round(y / yfactor)``. By default
    ``yfactor`` keeps about 6 significant digits.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 2:
        raise ValueError("x and y must be 1-D arrays of the same length (>= 2)")
    steps = np.diff(x)
    if not np.allclose(steps, steps[0], rtol=1e-6, atol=0):
        raise ValueError("XYDATA requires evenly spaced x values")
    if form not in ("AFFN", "DIFDUP"):
        raise ValueError("form must be 'AFFN' or 'DIFDUP'")
    if yfactor is None:
        peak = float(np.nanmax(np.abs(y))) or 1.0
        yfactor = float(f"{peak / 1_000_000:.3g}")
    yi = np.rint(y / yfactor).astype(np.int64)
    n = len(x)

    head = [f"##TITLE={title}", "##JCAMP-DX=4.24", f"##DATA TYPE={data_type}",
            f"##ORIGIN={origin}", f"##OWNER={owner}"]
    for k, v in (extra or {}).items():
        head.append(f"##{k}={v}")
    head += [f"##XUNITS={xunits}", f"##YUNITS={yunits}", "##XFACTOR=1", f"##YFACTOR={_fmt(yfactor)}",
             f"##FIRSTX={_fmt(x[0])}", f"##LASTX={_fmt(x[-1])}", f"##NPOINTS={n}",
             f"##FIRSTY={_fmt(yi[0] * yfactor)}", f"##MAXY={_fmt(yi.max() * yfactor)}",
             f"##MINY={_fmt(yi.min() * yfactor)}", f"##XYDATA={_XYDATA_VARS}"]

    body: list[str] = []
    if form == "AFFN":
        i = 0
        while i < n:
            line = _fmt(x[i])
            while i < n and len(line) + 1 + len(str(yi[i])) <= width:
                line += " " + str(yi[i])
                i += 1
            body.append(line)
    else:
        i = 0
        while True:
            line = f"{_fmt(x[i])} {_sqz(int(yi[i]))}"
            j = i
            while j + 1 < n:
                d = int(yi[j + 1] - yi[j])
                k = j + 1
                while k + 1 < n and int(yi[k + 1] - yi[k]) == d:
                    k += 1
                run = k - j
                token = _dif(d) + (_dup(run) if run > 1 else "")
                if len(line) + len(token) > width and j > i:   # always emit >= 1 token per line
                    break
                line += token
                j = k
            body.append(line)
            if j == n - 1:
                break
            i = j                     # next line repeats point j as its Y-check
        body.append(f"{_fmt(x[n - 1])} {_sqz(int(yi[n - 1]))}")   # final Y-check line
    return "\n".join(head + body + ["##END="]) + "\n"


# ----------------------------------------------------------------------------- CLI
def _summary(path: str, sp: Spectrum) -> str:
    xu, yu = sp.header("XUNITS", "?"), sp.header("YUNITS", "?")
    lines = [f"{path}",
             f"  title      : {sp.title}",
             f"  data type  : {sp.header('DATA TYPE', '?')}  (table ##{sp.table}, {len(sp.x)} points)"]
    if len(sp.x):
        lines += [f"  x range    : {np.nanmin(sp.x):g} .. {np.nanmax(sp.x):g} {xu}",
                  f"  y range    : {np.nanmin(sp.y):g} .. {np.nanmax(sp.y):g} {yu}"]
    lines.append("  validation : " + ("OK" if sp.ok else f"{len(sp.errors)} error(s)")
                 + (f", {len(sp.issues) - len(sp.errors)} warning(s)" if len(sp.issues) > len(sp.errors) else ""))
    lines += [f"    {issue}" for issue in sp.issues]
    return "\n".join(lines)


def _plot(sp: Spectrum, out: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.5))
    if sp.table == "XYDATA":
        ax.plot(sp.x, sp.y, lw=1)
    else:
        ax.vlines(sp.x, 0, sp.y)
    if "INFRARED" in (sp.header("DATA TYPE", "") or "").upper():
        ax.invert_xaxis()             # IR convention: wavenumber decreasing to the right
    ax.set_title(sp.title or "JCAMP-DX spectrum")
    ax.set_xlabel(sp.header("XUNITS", "x") or "x")
    ax.set_ylabel(sp.header("YUNITS", "y") or "y")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Read and validate JCAMP-DX spectra.")
    ap.add_argument("files", nargs="+", help="JCAMP-DX files (.jdx / .dx)")
    ap.add_argument("--csv", help="write x,y of the (single) input file to this CSV")
    ap.add_argument("--plot", help="save a plot of the (single) input file to this PNG")
    ap.add_argument("--strict", action="store_true", help="treat validation errors as fatal")
    args = ap.parse_args(argv)
    if (args.csv or args.plot) and len(args.files) != 1:
        ap.error("--csv/--plot need exactly one input file")

    status = 0
    for path in args.files:
        try:
            sp = read(path, strict=args.strict)
        except (JcampError, OSError) as exc:
            print(f"{path}\n  ERROR: {exc}")
            status = 1
            continue
        print(_summary(path, sp))
        status |= 0 if sp.ok else 1
        if args.csv:
            np.savetxt(args.csv, np.column_stack([sp.x, sp.y]), delimiter=",",
                       header=f"{sp.header('XUNITS', 'x')},{sp.header('YUNITS', 'y')}", comments="")
            print(f"  wrote {args.csv}")
        if args.plot:
            _plot(sp, args.plot)
            print(f"  wrote {args.plot}")
    return status


if __name__ == "__main__":
    sys.exit(main())
