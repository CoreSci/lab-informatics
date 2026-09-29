"""Generate the synthetic JCAMP-DX sample files in ``sample-data/``.

Every spectrum is synthetic: a sum of Gaussian bands at arbitrary positions,
plus a baseline and seeded noise. None is modelled on a real compound or
instrument. Re-running the script reproduces the files byte-for-byte.

Files written:
  * ``synthetic_ir_absorbance.jdx``: IR, 4000 -> 400 cm-1 (decreasing x), DIF/DUP compressed
  * ``synthetic_uv_vis.jdx``: UV-Vis, 200 -> 800 nm, AFFN (plain integers)
  * ``synthetic_ir_peak_table.jdx``: the IR band list as ``##PEAK TABLE=(XY..XY)``
  * ``invalid_example.jdx``: the IR file with deliberate defects, to demo the validator

Run: ``python generate_synthetic_spectra.py``
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

import jcamp_dx

OUT = Path(__file__).resolve().parent / "sample-data"

# (centre, height, width) of synthetic absorption bands
IR_BANDS = [(3350, 0.35, 120), (2930, 0.55, 25), (2860, 0.40, 20), (1715, 0.90, 18),
            (1460, 0.30, 15), (1375, 0.22, 12), (1100, 0.60, 30), (720, 0.15, 10)]
UV_BANDS = [(262, 0.85, 14), (410, 0.35, 35), (560, 0.12, 40)]


def bands(x: np.ndarray, spec: list[tuple[float, float, float]]) -> np.ndarray:
    """Sum of Gaussian bands evaluated at x."""
    return sum(h * np.exp(-0.5 * ((x - c) / w) ** 2) for c, h, w in spec)


def ir_spectrum(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(4000, 400, 1801)                       # 2 cm-1 steps, decreasing
    baseline = 0.02 + 0.01 * (4000 - x) / 3600
    return x, baseline + bands(x, IR_BANDS) + rng.normal(0, 0.002, x.size)


def uv_spectrum(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    x = np.linspace(200, 800, 601)                         # 1 nm steps
    return x, bands(x, UV_BANDS) + rng.normal(0, 0.003, x.size)


def peak_table_text() -> str:
    lines = ["##TITLE=Synthetic IR spectrum: band list", "##JCAMP-DX=4.24", "##DATA TYPE=INFRARED PEAK TABLE",
             "##ORIGIN=synthetic data", "##OWNER=public domain", "##XUNITS=1/CM", "##YUNITS=ABSORBANCE",
             "##XFACTOR=1", "##YFACTOR=1", f"##NPOINTS={len(IR_BANDS)}", "##PEAK TABLE=(XY..XY)"]
    lines += [f"{c},{h + 0.02:.3f}" for c, h, _ in IR_BANDS]
    return "\n".join(lines + ["##END="]) + "\n"


def corrupt(text: str) -> str:
    """Introduce three defects the validator should report."""
    lines = text.splitlines()
    lines = [("##NPOINTS=1800" if l.startswith("##NPOINTS=") else l) for l in lines]   # wrong point count
    data_start = next(i for i, l in enumerate(lines) if l.startswith("##XYDATA")) + 1
    # break the Y-check of the 3rd data line: change its first ordinate (the SQZ check value)
    x_tok, rest = lines[data_start + 2].split(" ", 1)
    lines[data_start + 2] = f"{x_tok} {'B' if rest[0] != 'B' else 'C'}{rest[1:]}"
    # shift the abscissa of the 5th data line (X-check)
    x_tok, rest = lines[data_start + 4].split(" ", 1)
    lines[data_start + 4] = f"{float(x_tok) + 50:g} {rest}"
    lines = [l for l in lines if not l.startswith("##OWNER")]                            # missing label
    return "\n".join(lines) + "\n"


def main() -> None:
    rng = np.random.default_rng(2024)
    OUT.mkdir(exist_ok=True)
    x, y = ir_spectrum(rng)
    ir = jcamp_dx.write(x, y, title="Synthetic IR spectrum (absorbance)", data_type="INFRARED SPECTRUM",
                        xunits="1/CM", yunits="ABSORBANCE", form="DIFDUP",
                        extra={"SAMPLING PROCEDURE": "simulated", "RESOLUTION": "2"})
    x2, y2 = uv_spectrum(rng)
    uv = jcamp_dx.write(x2, y2, title="Synthetic UV-Vis spectrum", data_type="UV/VIS SPECTRUM",
                        xunits="NANOMETERS", yunits="ABSORBANCE", form="AFFN")
    files = {"synthetic_ir_absorbance.jdx": ir, "synthetic_uv_vis.jdx": uv,
             "synthetic_ir_peak_table.jdx": peak_table_text(), "invalid_example.jdx": corrupt(ir)}
    for name, text in files.items():
        with open(OUT / name, "w", encoding="ascii", newline="\n") as f:
            f.write(text)
        print(f"wrote {OUT / name} ({len(text.splitlines())} lines)")


if __name__ == "__main__":
    main()
