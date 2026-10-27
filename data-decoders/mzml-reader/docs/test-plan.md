# Test plan: mzML tools

| | |
|---|---|
| Version | 1.0 |
| Date | 2026-10-23 |
| Covers | `data-decoders/mzml-reader/`: `mzml.py` (reader), `mzml_validator.py`, `mzml_chromatograms.py`, `mzml_writer.py`, and the sample generators |
| Results | [`traceability-matrix.md`](traceability-matrix.md) (dated snapshot), and the JUnit reports and matrix of every CI run (artifacts of the `tests` workflow) |

## 1. Purpose and scope

This plan says what is tested, why, how, and what is **not** tested. The tools are educational and work on synthetic data. The plan is about making sure the tools do what their README claims, and that they fail loudly instead of returning wrong numbers.

**In scope:** reading plain and indexed mzML; decoding binary arrays; the validator's checks (MZV001–MZV053); TIC, BPC and XIC extraction with CSV and SVG output; the writer; determinism of the sample generators.

**Out of scope:** files from real acquisitions; MS-Numpress and other compression schemes; integer and string arrays; the mzML semantic mapping rules; any other file format. These are outside the tools' own scope, as their README states.

## 2. Requirements

The requirements are listed, each with its source, in [`../traceability/requirements.tsv`](../traceability/requirements.tsv). 37 of them are claims already made in the README; one (TR-01) covers the traceability tooling described in this plan. The 38 ids are grouped by area:

- `RD` reader
- `VL` validator
- `CH` chromatograms
- `WR` writer
- `DG` sample data
- `TR` the traceability tooling itself

## 3. Approach: risk-based

Effort goes where a defect would be silent, and so most harmful. The highest risk is a wrong number returned without any error. A crash or a clear error message is the lesser risk.

| Risk | Impact | How it is tested | Residual risk |
|---|---|---|---|
| Wrong decoded values (byte order, precision, compression) | High: silent | Round trip of every encoding combination (32/64-bit × none/zlib), in both the reader and the writer; sample files checked against values computed independently in the tests | Real files with encodings not in the samples; big-endian hosts (the byte-swap path is not exercised) |
| Bad input accepted as if valid | High: silent | One negative test per reader error (unsupported compression and data types, length mismatches, invalid base64/zlib, malformed XML); each validator code has a failing case and a passing case | Defect types not anticipated in the tests |
| Validator false positives or false negatives | Medium | Each code fires alone on a targeted edit; valid documents and samples give zero findings; writer output must pass `--strict`. *Example:* the writer tests found a false MZV050 for ids with XML-escaped characters, fixed in its own commit | Combinations of several defects in one file |
| Wrong chromatogram values | High: silent | TIC, BPC and XIC compared with sums computed directly from the spectra; the caffeine XIC must peak at the 2.80 min apex; a window away from the compound must stay at zero; ppm and Da widths must agree | Profile-mode data and overlapping isotope windows only partly covered |
| Unit confusion (minutes vs seconds) | Medium | Units returned as stored; mixed units refused; the TIC cross-check refuses to compare differing units | None known |
| Memory growth on large files | Medium | Peak traced memory for 1,500 spectra must stay below twice that for 100 spectra; the first spectrum must arrive before 10 % of the file is read | Absolute speed is not measured |
| Corrupt index or checksum in written files | Medium | Offsets, `indexListOffset` and SHA-1 recomputed independently in the tests | — |
| Non-reproducible sample data | Low | Generators run twice must give byte-identical files; shipped zlib files are compared by decoded values, the others byte for byte | zlib byte-identity across zlib implementations is not tested |
| Behaviour differs between Python versions | Medium | CI runs the full suite on every supported version (see 5) | Operating systems other than Linux in CI |

## 4. Test levels and techniques

- **Unit tests** on small in-memory documents, one per error case or check, built by the tests or by the writer.
- **Sample-based tests** on the three synthetic files, with expected values computed independently (sums, isotope ratio, elution apex).
- **Negative tests** for every documented error and every validator code.
- **Round-trip tests:** writer → reader → validator.
- **Determinism tests** for the generators and the SVG plots.
- **Command-line tests:** exit codes and `--strict` for the validator; summary, CSV and SVG for the chromatogram tool.

## 5. Environment

- **CI:** GitHub Actions workflow `tests`, `ubuntu-24.04`, Python 3.9, 3.10, 3.11, 3.12, 3.13, 3.14 and 3.15. Dependencies pinned in `.github/requirements-ci.txt` (the mzML tools need only pytest).
- **Before the workflow was added:** the same suites passed in Docker `python:3.x-slim` images for all seven versions on 2026-10-23.
- **During development:** the mzML suite also ran on Windows with Python 3.7.4.

## 6. Entry and exit criteria

- **Entry:** the change is complete, and the suite runs locally on at least one supported version.
- **Exit:** every test passes on every Python version in the CI matrix; no requirement has status FAILED in the generated matrix. Requirements with status "not tested" are acceptable only when listed as known gaps in section 8.

## 7. Traceability

1. Each test carries `@pytest.mark.req("ID", ...)` naming the requirements it verifies.
2. `conftest.py` copies those ids into the JUnit XML (`pytest.ini` selects the `xunit1` format, which keeps per-test properties).
3. `traceability/build_matrix.py` reads the JUnit XML and `requirements.tsv`, and writes the matrix:
   - every requirement with its tests and their results;
   - requirements without tests;
   - tests without a requirement.

   It stops with an error if a marker names an unknown id.
4. The CI workflow generates the matrix for every Python version, as an artifact.
5. The committed [`traceability-matrix.md`](traceability-matrix.md) is one such generated snapshot. Its header records when, from what and on which Python version it was generated. It is never edited by hand.

## 8. Known gaps

- **No real data:** never run on files from real acquisitions. All test data is synthetic or built in the tests.
- **RD-11:** the reader's command-line summary (`python mzml.py FILE`) has no automated test.
- **Schema check:** the XSD check of the sample files and writer output was run when they were made, with a third-party library; it is not part of the suite.
- **Platforms:** big-endian platforms and operating systems other than Linux (in CI) are not tested.
- **Performance:** checked only as relative memory growth; no timing requirements.
- **zlib builds:** byte-identity of zlib-compressed samples across zlib implementations is not tested.

## 9. Change history

| Version | Date | Change |
|---|---|---|
| 1.0 | 2026-10-23 | First version: requirements list, markers on every test, generated matrix, CI artifacts |
