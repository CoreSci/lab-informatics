"""Generate the example chromatogram plots in ``sample-data/output/``.

Both plots are computed from the synthetic LC-MS sample written by
``generate_synthetic_mzml.py``, and drawn by ``mzml_chromatograms.svg_plot``
(standard library only). The output depends only on the decoded values, so
re-running the script reproduces the files byte for byte.

Files written:
  * ``synthetic_caffeine_tic.svg``: TIC of the 60 MS1 spectra
  * ``synthetic_caffeine_xic_195.svg``: XIC of caffeine [M+H]+, m/z 195.0877 +/- 10 ppm

Run: ``python generate_example_plots.py``
"""
from __future__ import annotations

from pathlib import Path
from typing import List

import mzml_chromatograms as chrom

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "sample-data" / "synthetic_caffeine_lc_ms_centroid.mzML"
OUT = HERE / "sample-data" / "output"
CAFFEINE_MH = 195.0877
XIC_TOLERANCE_PPM = 10.0


def generate(out_dir: Path = OUT, source: Path = SOURCE) -> List[Path]:
    """Write both SVG files to ``out_dir``."""
    result = chrom.extract(source, [(CAFFEINE_MH, XIC_TOLERANCE_PPM, "ppm")])
    out_dir.mkdir(parents=True, exist_ok=True)
    plots = {
        "synthetic_caffeine_tic.svg": chrom.svg_plot(result.tic, title="Synthetic LC-MS run: TIC (MS1)"),
        "synthetic_caffeine_xic_195.svg": chrom.svg_plot(result.xics[0], title=f"Synthetic caffeine [M+H]+: {result.xics[0].label}"),
    }
    paths = []
    for name, svg in plots.items():
        path = out_dir / name
        path.write_bytes(svg.encode("utf-8"))           # bytes: LF line endings on every OS
        paths.append(path)
    return paths


if __name__ == "__main__":
    for p in generate():
        print(f"wrote output/{p.name} ({p.stat().st_size:,} bytes)")
