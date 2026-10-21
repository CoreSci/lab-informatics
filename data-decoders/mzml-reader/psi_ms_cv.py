"""PSI-MS and Unit Ontology terms used by the mzML reader and the sample generator.

Every accession below was checked against the PSI-MS controlled vocabulary
(``psi-ms.obo``, data-version 4.2.6, https://github.com/HUPO-PSI/psi-ms-CV):
the name matches and the term is not obsolete. The UO accessions are the units
that ``psi-ms.obo`` itself declares with ``has_units`` relationships.

Each constant is an ``(accession, name)`` pair, so code and generated files
always carry both.
"""
from __future__ import annotations

PSI_MS_VERSION = "4.2.6"

# ------------------------------------------------------------------ binary data arrays
MZ_ARRAY = ("MS:1000514", "m/z array")
INTENSITY_ARRAY = ("MS:1000515", "intensity array")
TIME_ARRAY = ("MS:1000595", "time array")

FLOAT_32 = ("MS:1000521", "32-bit float")
FLOAT_64 = ("MS:1000523", "64-bit float")

ZLIB_COMPRESSION = ("MS:1000574", "zlib compression")
NO_COMPRESSION = ("MS:1000576", "no compression")

# Parent terms, used to recognise array-type and compression terms this reader
# does not support (any other child of these is rejected with a clear error).
BINARY_DATA_ARRAY = ("MS:1000513", "binary data array")
BINARY_DATA_COMPRESSION_TYPE = ("MS:1000572", "binary data compression type")

# MS-Numpress encodings (alone or followed by zlib / zstd). Out of scope here.
NUMPRESS_ACCESSIONS = frozenset({
    "MS:1002312", "MS:1002313", "MS:1002314",   # linear prediction, positive integer, short logged float
    "MS:1002746", "MS:1002747", "MS:1002748",   # ... followed by zlib compression
    "MS:1003783", "MS:1003784", "MS:1003785",   # ... followed by zstd compression
})

# ------------------------------------------------------------------ spectrum description
MS_LEVEL = ("MS:1000511", "ms level")
MS1_SPECTRUM = ("MS:1000579", "MS1 spectrum")
MSN_SPECTRUM = ("MS:1000580", "MSn spectrum")
CENTROID_SPECTRUM = ("MS:1000127", "centroid spectrum")
PROFILE_SPECTRUM = ("MS:1000128", "profile spectrum")
POSITIVE_SCAN = ("MS:1000130", "positive scan")
TOTAL_ION_CURRENT = ("MS:1000285", "total ion current")
BASE_PEAK_MZ = ("MS:1000504", "base peak m/z")
BASE_PEAK_INTENSITY = ("MS:1000505", "base peak intensity")
LOWEST_OBSERVED_MZ = ("MS:1000528", "lowest observed m/z")
HIGHEST_OBSERVED_MZ = ("MS:1000527", "highest observed m/z")

SCAN_START_TIME = ("MS:1000016", "scan start time")
NO_COMBINATION = ("MS:1000795", "no combination")
SCAN_WINDOW_LOWER_LIMIT = ("MS:1000501", "scan window lower limit")
SCAN_WINDOW_UPPER_LIMIT = ("MS:1000500", "scan window upper limit")

# ------------------------------------------------------------------ precursors (MS2)
ISOLATION_WINDOW_TARGET_MZ = ("MS:1000827", "isolation window target m/z")
ISOLATION_WINDOW_LOWER_OFFSET = ("MS:1000828", "isolation window lower offset")
ISOLATION_WINDOW_UPPER_OFFSET = ("MS:1000829", "isolation window upper offset")
SELECTED_ION_MZ = ("MS:1000744", "selected ion m/z")
CHARGE_STATE = ("MS:1000041", "charge state")
PEAK_INTENSITY = ("MS:1000042", "peak intensity")
COLLISION_INDUCED_DISSOCIATION = ("MS:1000133", "collision-induced dissociation")
BEAM_TYPE_COLLISION_INDUCED_DISSOCIATION = ("MS:1000422", "beam-type collision-induced dissociation")
COLLISION_ENERGY = ("MS:1000045", "collision energy")

# ------------------------------------------------------------------ chromatograms
TIC_CHROMATOGRAM = ("MS:1000235", "total ion current chromatogram")
BASEPEAK_CHROMATOGRAM = ("MS:1000628", "basepeak chromatogram")
SELECTED_ION_CURRENT_CHROMATOGRAM = ("MS:1000627", "selected ion current chromatogram")

# ------------------------------------------------------------------ generic file metadata
# Generic terms only: no instrument model, manufacturer or acquisition software.
INSTRUMENT_MODEL = ("MS:1000031", "instrument model")
ELECTROSPRAY_IONIZATION = ("MS:1000073", "electrospray ionization")
QUADRUPOLE = ("MS:1000081", "quadrupole")
ELECTRON_MULTIPLIER = ("MS:1000253", "electron multiplier")
CUSTOM_SOFTWARE = ("MS:1000799", "custom unreleased software tool")
DATA_PROCESSING_ACTION = ("MS:1000543", "data processing action")

# ------------------------------------------------------------------ units
MZ_UNIT = ("MS:1000040", "m/z")
DETECTOR_COUNTS = ("MS:1000131", "number of detector counts")
UO_SECOND = ("UO:0000010", "second")
UO_MILLISECOND = ("UO:0000028", "millisecond")
UO_MINUTE = ("UO:0000031", "minute")
UO_ELECTRONVOLT = ("UO:0000266", "electronvolt")

TIME_UNITS = {acc: name for acc, name in (UO_SECOND, UO_MILLISECOND, UO_MINUTE)}
