"""Generate the synthetic mzML sample files in ``sample-data/``.

Every file is synthetic. The model compound is caffeine (C8H10N4O2), observed
as [M+H]+ at m/z 195.0877. Its isotope pattern is computed below from natural
isotopic abundances. Elution profile, intensities, background ions and noise
are invented and seeded, so re-running the script reproduces the files byte
for byte (with the same zlib build: zlib output can differ between zlib
implementations, while the decoded values do not).

Instrument metadata uses generic PSI-MS terms only ("instrument model",
electrospray ionization, quadrupole, electron multiplier): it describes no
real instrument.

Files written:
  * ``synthetic_caffeine_lc_ms_centroid.mzML``: 60 centroid MS1 spectra of a
    simulated LC-MS run with a Gaussian elution peak, 64-bit floats, zlib, plus
    a total ion current chromatogram
  * ``synthetic_caffeine_profile.mzML``: 5 profile-mode MS1 spectra, 32-bit
    floats, no compression
  * ``synthetic_caffeine_ms2_indexed.mzML``: indexedmzML wrapper (valid offsets
    and SHA-1 checksum), 3 MS1 + 3 MS2 spectra, 64-bit m/z and 32-bit
    intensities, zlib, plus a total ion current chromatogram

Run: ``python generate_synthetic_mzml.py``
"""
from __future__ import annotations

import hashlib
import math
import random
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import mzml
import psi_ms_cv as cv

OUT = Path(__file__).resolve().parent / "sample-data"
SEED = 20261019

# Isotope masses (u) and natural abundances (mole fraction): NIST "Atomic Weights
# and Isotopic Compositions" (representative terrestrial compositions).
ISOTOPES: Dict[str, List[Tuple[float, float]]] = {
    "C": [(12.0, 0.9893), (13.00335483507, 0.0107)],
    "H": [(1.00782503223, 0.999885), (2.01410177812, 0.000115)],
    "N": [(14.00307400443, 0.99636), (15.00010889888, 0.00364)],
    "O": [(15.99491461957, 0.99757), (16.99913175650, 0.00038), (17.99915961286, 0.00205)],
}
ELECTRON = 0.000548579909
PROTON = ISOTOPES["H"][0][0] - ELECTRON

CAFFEINE = {"C": 8, "H": 10, "N": 4, "O": 2}
# Illustrative product ions for the MS2 spectra (cation formulas); relative
# intensities are invented.
FRAGMENTS = [({"C": 6, "H": 8, "N": 3, "O": 1}, 1.00), ({"C": 5, "H": 8, "N": 3}, 0.35),
             ({"C": 4, "H": 7, "N": 2}, 0.12), ({"C": 3, "H": 6, "N": 1}, 0.08)]

MS_NS = "http://psi.hupo.org/ms/mzml"
CV_LIST = [
    ("MS", "Proteomics Standards Initiative Mass Spectrometry Ontology", cv.PSI_MS_VERSION,
     "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"),
    ("UO", "Unit Ontology", None,
     "https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"),
]


# ------------------------------------------------------------------ chemistry
def monoisotopic_mass(formula: Dict[str, int]) -> float:
    return sum(ISOTOPES[el][0][0] * n for el, n in formula.items())


def isotope_pattern(formula: Dict[str, int], min_rel: float = 1e-3) -> List[Tuple[float, float]]:
    """Neutral isotope pattern: [(mass, abundance relative to the monoisotopic peak)].

    Convolves every atom's isotope distribution, then merges the fine structure
    by nominal mass shift (abundance-weighted mean mass), as a mass spectrometer
    of moderate resolution would see it.
    """
    dist = [(0.0, 1.0)]
    for el, n in formula.items():
        for _ in range(n):
            nxt: Dict[float, float] = {}
            for m, p in dist:
                for im, ip in ISOTOPES[el]:
                    key = round(m + im, 9)
                    nxt[key] = nxt.get(key, 0.0) + p * ip
            dist = [(m, p) for m, p in nxt.items() if p > 1e-12]
    mono = monoisotopic_mass(formula)
    groups: Dict[int, List[float]] = {}
    for m, p in dist:
        g = groups.setdefault(round(m - mono), [0.0, 0.0])
        g[0] += m * p
        g[1] += p
    base = groups[0][1]
    return [(g[0] / g[1], g[1] / base) for _, g in sorted(groups.items()) if g[1] / base >= min_rel]


def protonated(pattern: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """[M+H]+ m/z values for a neutral isotope pattern."""
    return [(m + PROTON, r) for m, r in pattern]


def gaussian(x: float, centre: float, sigma: float) -> float:
    return math.exp(-0.5 * ((x - centre) / sigma) ** 2)


# ------------------------------------------------------------------ mzML text
def _cv_ref(acc: str) -> str:
    return acc.split(":")[0]


def cv_param(term: Tuple[str, str], value="", unit: Tuple[str, str] = None) -> str:
    s = f'<cvParam cvRef="{_cv_ref(term[0])}" accession="{term[0]}" name="{term[1]}" value="{value}"'
    if unit:
        s += f' unitCvRef="{_cv_ref(unit[0])}" unitAccession="{unit[0]}" unitName="{unit[1]}"'
    return s + "/>"


def binary_array(values: Sequence[float], bits: int, zlib_compressed: bool,
                 array_type: Tuple[str, str], unit: Tuple[str, str]) -> List[str]:
    text = mzml.encode_binary(values, bits, zlib_compressed)
    return [f'<binaryDataArray encodedLength="{len(text)}">',
            "  " + cv_param(cv.FLOAT_32 if bits == 32 else cv.FLOAT_64),
            "  " + cv_param(cv.ZLIB_COMPRESSION if zlib_compressed else cv.NO_COMPRESSION),
            "  " + cv_param(array_type, unit=unit),
            f"  <binary>{text}</binary>",
            "</binaryDataArray>"]


def indent(lines: List[str], n: int) -> List[str]:
    return [" " * n + line for line in lines]


def spectrum_xml(index: int, rt_min: float, mz: List[float], inten: List[float], *, ms_level: int,
                 centroid: bool, window: Tuple[float, float], mz_bits: int, int_bits: int,
                 zlib_compressed: bool, precursor: dict = None) -> List[str]:
    base = max(range(len(inten)), key=inten.__getitem__)
    lines = [f'<spectrum index="{index}" id="scan={index + 1}" defaultArrayLength="{len(mz)}">',
             "  " + cv_param(cv.MS_LEVEL, ms_level),
             "  " + cv_param(cv.MS1_SPECTRUM if ms_level == 1 else cv.MSN_SPECTRUM),
             "  " + cv_param(cv.CENTROID_SPECTRUM if centroid else cv.PROFILE_SPECTRUM),
             "  " + cv_param(cv.POSITIVE_SCAN),
             "  " + cv_param(cv.BASE_PEAK_MZ, f"{mz[base]:.5f}", cv.MZ_UNIT),
             "  " + cv_param(cv.BASE_PEAK_INTENSITY, f"{inten[base]:.1f}", cv.DETECTOR_COUNTS),
             "  " + cv_param(cv.TOTAL_ION_CURRENT, f"{sum(inten):.1f}"),
             "  " + cv_param(cv.LOWEST_OBSERVED_MZ, f"{mz[0]:.5f}", cv.MZ_UNIT),
             "  " + cv_param(cv.HIGHEST_OBSERVED_MZ, f"{mz[-1]:.5f}", cv.MZ_UNIT),
             '  <scanList count="1">',
             "    " + cv_param(cv.NO_COMBINATION),
             "    <scan>",
             "      " + cv_param(cv.SCAN_START_TIME, f"{rt_min:.4f}", cv.UO_MINUTE),
             '      <scanWindowList count="1">',
             "        <scanWindow>",
             "          " + cv_param(cv.SCAN_WINDOW_LOWER_LIMIT, f"{window[0]:.1f}", cv.MZ_UNIT),
             "          " + cv_param(cv.SCAN_WINDOW_UPPER_LIMIT, f"{window[1]:.1f}", cv.MZ_UNIT),
             "        </scanWindow>",
             "      </scanWindowList>",
             "    </scan>",
             "  </scanList>"]
    if precursor:
        p = precursor
        lines += ['  <precursorList count="1">',
                  f'    <precursor spectrumRef="{p["ref"]}">',
                  "      <isolationWindow>",
                  "        " + cv_param(cv.ISOLATION_WINDOW_TARGET_MZ, f"{p['mz']:.5f}", cv.MZ_UNIT),
                  "        " + cv_param(cv.ISOLATION_WINDOW_LOWER_OFFSET, "0.5", cv.MZ_UNIT),
                  "        " + cv_param(cv.ISOLATION_WINDOW_UPPER_OFFSET, "0.5", cv.MZ_UNIT),
                  "      </isolationWindow>",
                  '      <selectedIonList count="1">',
                  "        <selectedIon>",
                  "          " + cv_param(cv.SELECTED_ION_MZ, f"{p['mz']:.5f}", cv.MZ_UNIT),
                  "          " + cv_param(cv.CHARGE_STATE, 1),
                  "          " + cv_param(cv.PEAK_INTENSITY, f"{p['intensity']:.1f}", cv.DETECTOR_COUNTS),
                  "        </selectedIon>",
                  "      </selectedIonList>",
                  "      <activation>",
                  "        " + cv_param(cv.COLLISION_INDUCED_DISSOCIATION),
                  "        " + cv_param(cv.COLLISION_ENERGY, "20.0", cv.UO_ELECTRONVOLT),
                  "      </activation>",
                  "    </precursor>",
                  "  </precursorList>"]
    lines += ['  <binaryDataArrayList count="2">',
              *indent(binary_array(mz, mz_bits, zlib_compressed, cv.MZ_ARRAY, cv.MZ_UNIT), 4),
              *indent(binary_array(inten, int_bits, zlib_compressed, cv.INTENSITY_ARRAY, cv.DETECTOR_COUNTS), 4),
              "  </binaryDataArrayList>",
              "</spectrum>"]
    return lines


def tic_chromatogram_xml(times: List[float], tic: List[float], int_bits: int) -> List[str]:
    return [f'<chromatogram index="0" id="TIC" defaultArrayLength="{len(times)}">',
            "  " + cv_param(cv.TIC_CHROMATOGRAM),
            '  <binaryDataArrayList count="2">',
            *indent(binary_array(times, 64, True, cv.TIME_ARRAY, cv.UO_MINUTE), 4),
            *indent(binary_array(tic, int_bits, True, cv.INTENSITY_ARRAY, cv.DETECTOR_COUNTS), 4),
            "  </binaryDataArrayList>",
            "</chromatogram>"]


def mzml_lines(run_id: str, spectra: List[List[str]], chromatograms: List[List[str]],
               has_msn: bool, root_attrs: str) -> List[str]:
    """A complete <mzML> element (generic metadata) around the given spectra and chromatograms."""
    cvs = []
    for cid, full, version, uri in CV_LIST:
        v = f' version="{version}"' if version else ""
        cvs.append(f'    <cv id="{cid}" fullName="{full}"{v} URI="{uri}"/>')
    content = ["      " + cv_param(cv.MS1_SPECTRUM)] + (["      " + cv_param(cv.MSN_SPECTRUM)] if has_msn else [])
    lines = [f"<mzML {root_attrs}>",
             f'  <cvList count="{len(cvs)}">', *cvs, "  </cvList>",
             "  <fileDescription>", "    <fileContent>", *content, "    </fileContent>", "  </fileDescription>",
             '  <softwareList count="1">',
             '    <software id="synthetic_generator" version="1.0">',
             "      " + cv_param(cv.CUSTOM_SOFTWARE, "generate_synthetic_mzml.py"),
             "    </software>",
             "  </softwareList>",
             '  <instrumentConfigurationList count="1">',
             '    <instrumentConfiguration id="IC1">',
             "      " + cv_param(cv.INSTRUMENT_MODEL),
             '      <componentList count="3">',
             '        <source order="1">', "          " + cv_param(cv.ELECTROSPRAY_IONIZATION), "        </source>",
             '        <analyzer order="2">', "          " + cv_param(cv.QUADRUPOLE), "        </analyzer>",
             '        <detector order="3">', "          " + cv_param(cv.ELECTRON_MULTIPLIER), "        </detector>",
             "      </componentList>",
             "    </instrumentConfiguration>",
             "  </instrumentConfigurationList>",
             '  <dataProcessingList count="1">',
             '    <dataProcessing id="DP1">',
             '      <processingMethod order="1" softwareRef="synthetic_generator">',
             "        " + cv_param(cv.DATA_PROCESSING_ACTION, "synthetic data generation"),
             "      </processingMethod>",
             "    </dataProcessing>",
             "  </dataProcessingList>",
             f'  <run id="{run_id}" defaultInstrumentConfigurationRef="IC1">',
             f'    <spectrumList count="{len(spectra)}" defaultDataProcessingRef="DP1">']
    for s in spectra:
        lines += indent(s, 6)
    lines.append("    </spectrumList>")
    if chromatograms:
        lines.append(f'    <chromatogramList count="{len(chromatograms)}" defaultDataProcessingRef="DP1">')
        for c in chromatograms:
            lines += indent(c, 6)
        lines.append("    </chromatogramList>")
    return lines + ["  </run>", "</mzML>"]


XML_DECL = '<?xml version="1.0" encoding="utf-8"?>'
XSI = 'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"'


def plain_mzml(run_id: str, spectra, chromatograms, has_msn=False) -> bytes:
    attrs = (f'xmlns="{MS_NS}" {XSI} xsi:schemaLocation="{MS_NS} '
             f'http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" id="{run_id}" version="1.1.0"')
    lines = [XML_DECL] + mzml_lines(run_id, spectra, chromatograms, has_msn, attrs)
    return ("\n".join(lines) + "\n").encode("ascii")


def indexed_mzml(run_id: str, spectra, chromatograms, has_msn=False) -> bytes:
    """An <indexedmzML> file with byte offsets for every spectrum and chromatogram and a SHA-1 checksum."""
    head = [XML_DECL,
            f'<indexedmzML xmlns="{MS_NS}" {XSI} xsi:schemaLocation="{MS_NS} '
            f'http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">']
    body = indent(mzml_lines(run_id, spectra, chromatograms, has_msn, f'id="{run_id}" version="1.1.0"'), 2)
    out = bytearray()
    offsets: Dict[str, List[Tuple[str, int]]] = {"spectrum": [], "chromatogram": []}
    for line in head + body:
        stripped = line.lstrip()
        for kind in offsets:
            if stripped.startswith(f"<{kind} "):
                ref = stripped.split(' id="', 1)[1].split('"', 1)[0]
                offsets[kind].append((ref, len(out) + len(line) - len(stripped)))
        out += (line + "\n").encode("ascii")
    index_offset = len(out) + 2                          # the '<' of '  <indexList'
    kinds = [k for k in offsets if offsets[k]]
    tail = [f'  <indexList count="{len(kinds)}">']
    for kind in kinds:
        tail.append(f'    <index name="{kind}">')
        tail += [f'      <offset idRef="{ref}">{pos}</offset>' for ref, pos in offsets[kind]]
        tail.append("    </index>")
    tail += ["  </indexList>", f"  <indexListOffset>{index_offset}</indexListOffset>", "  <fileChecksum>"]
    out += "\n".join(tail).encode("ascii")
    checksum = hashlib.sha1(bytes(out)).hexdigest()   # from the start of the file to the end of <fileChecksum>
    out += f"{checksum}</fileChecksum>\n</indexedmzML>\n".encode("ascii")
    return bytes(out)


# ------------------------------------------------------------------ the three runs
def caffeine_ions() -> List[Tuple[float, float]]:
    return protonated(isotope_pattern(CAFFEINE))


def background_ions(rng: random.Random, n: int, lo: float, hi: float) -> List[Tuple[float, float]]:
    """Seeded background ions: (m/z, typical intensity)."""
    return sorted((round(rng.uniform(lo, hi), 4), rng.uniform(500, 5000)) for _ in range(n))


def centroid_peaks(rng: random.Random, ions, scale: float, background, threshold: float = 100.0):
    """Centroids of the compound (scaled) and the background, with seeded m/z and intensity noise."""
    peaks = []
    for m, rel in ions:
        inten = scale * rel * rng.uniform(0.95, 1.05)
        if inten >= threshold:
            peaks.append((m, inten))
    for m, typical in background:
        peaks.append((m, typical * rng.uniform(0.8, 1.2)))
    peaks = [(round(m * (1 + rng.gauss(0, 2e-6)), 5), round(i, 1)) for m, i in peaks]
    peaks.sort()
    return [m for m, _ in peaks], [i for _, i in peaks]


def lc_ms_centroid(rng: random.Random) -> bytes:
    """60 MS1 centroid spectra, 1.60 to 3.96 min, caffeine eluting at 2.80 min."""
    ions, background = caffeine_ions(), background_ions(rng, 25, 100.0, 400.0)
    spectra, times, tics = [], [], []
    for i in range(60):
        rt = round(1.60 + 0.04 * i, 4)
        mz, inten = centroid_peaks(rng, ions, 2.0e6 * gaussian(rt, 2.80, 0.10), background)
        spectra.append(spectrum_xml(i, rt, mz, inten, ms_level=1, centroid=True, window=(100.0, 400.0),
                                    mz_bits=64, int_bits=64, zlib_compressed=True))
        times.append(rt)
        tics.append(round(sum(inten), 1))
    return plain_mzml("synthetic_caffeine_lc_ms_centroid", spectra, [tic_chromatogram_xml(times, tics, 64)])


def profile_spectra(rng: random.Random) -> bytes:
    """5 profile-mode MS1 spectra of the [M+H]+ region, 193.5 to 199.5 in 0.005 steps."""
    ions = caffeine_ions()
    sigma = 0.02 / 2.3548                                 # FWHM 0.02
    grid = [round(193.5 + 0.005 * k, 4) for k in range(1201)]
    spectra = []
    for i in range(5):
        scale = 8.0e5 * (1 + 0.05 * i)
        inten = [round(sum(scale * r * gaussian(x, m, sigma) for m, r in ions) + rng.uniform(0, 50), 1)
                 for x in grid]
        spectra.append(spectrum_xml(i, round(0.10 * (i + 1), 4), grid, inten, ms_level=1, centroid=False,
                                    window=(193.5, 199.5), mz_bits=32, int_bits=32, zlib_compressed=False))
    return plain_mzml("synthetic_caffeine_profile", spectra, [])


def ms2_indexed(rng: random.Random) -> bytes:
    """3 cycles of MS1 + MS2 of [M+H]+, as an indexedmzML file."""
    ions, background = caffeine_ions(), background_ions(rng, 10, 100.0, 400.0)
    precursor_mz = ions[0][0]
    products = [(monoisotopic_mass(f) - ELECTRON, rel) for f, rel in FRAGMENTS] + [(precursor_mz, 0.15)]
    spectra, times, tics = [], [], []
    for cycle in range(3):
        rt = round(2.70 + 0.05 * cycle, 4)
        scale = 2.0e6 * gaussian(rt, 2.75, 0.10)
        mz, inten = centroid_peaks(rng, ions, scale, background)
        ms1 = len(spectra)
        spectra.append(spectrum_xml(ms1, rt, mz, inten, ms_level=1, centroid=True, window=(100.0, 400.0),
                                    mz_bits=64, int_bits=32, zlib_compressed=True))
        times.append(rt)
        tics.append(round(sum(inten), 1))
        rt2 = round(rt + 0.01, 4)
        mz2, inten2 = centroid_peaks(rng, sorted(products), 0.3 * scale, [])
        spectra.append(spectrum_xml(ms1 + 1, rt2, mz2, inten2, ms_level=2, centroid=True,
                                    window=(50.0, 200.0), mz_bits=64, int_bits=32, zlib_compressed=True,
                                    precursor={"ref": f"scan={ms1 + 1}", "mz": precursor_mz,
                                               "intensity": inten[mz.index(min(mz, key=lambda m: abs(m - precursor_mz)))]}))
        times.append(rt2)
        tics.append(round(sum(inten2), 1))
    return indexed_mzml("synthetic_caffeine_ms2_indexed", spectra, [tic_chromatogram_xml(times, tics, 32)],
                        has_msn=True)


FILES = {
    "synthetic_caffeine_lc_ms_centroid.mzML": lc_ms_centroid,
    "synthetic_caffeine_profile.mzML": profile_spectra,
    "synthetic_caffeine_ms2_indexed.mzML": ms2_indexed,
}


def generate(out_dir: Path = OUT, seed: int = SEED) -> List[Path]:
    """Write every sample file to ``out_dir``. Each file gets its own seeded generator."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for k, (name, build) in enumerate(FILES.items()):
        path = out_dir / name
        path.write_bytes(build(random.Random(seed + k)))   # bytes: LF line endings on every OS
        paths.append(path)
    return paths


if __name__ == "__main__":
    for p in generate():
        print(f"wrote {p.name} ({p.stat().st_size:,} bytes)")
