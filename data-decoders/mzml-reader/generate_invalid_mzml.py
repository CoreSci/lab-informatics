"""Generate the deliberately invalid mzML samples in ``sample-data/invalid/``.

Each file is a valid synthetic file from ``generate_synthetic_mzml.py``,
rebuilt with the same seed, then corrupted deterministically with exactly one
defect for the validator to find. The valid samples are not touched.

Files written (defect -> expected finding):
  * ``invalid_spectrum_count.mzML``: spectrumList count says 6, the list holds 5 (MZV010)
  * ``invalid_array_length.mzML``: one intensity array declares arrayLength 1200 but holds 1201 values (MZV020)
  * ``invalid_duplicate_id.mzML``: the second spectrum reuses the first spectrum's id (MZV012)
  * ``invalid_mz_order.mzML``: two neighbouring m/z values swapped in the first spectrum (MZV031)
  * ``invalid_index_offset.mzML``: one index offset is 1 byte too far; checksum recomputed (MZV050)
  * ``invalid_checksum.mzML``: the SHA-1 fileChecksum doesn't match the file (MZV052)

The first four start from the profile sample (no compression, so they are
byte-identical on every platform). The last two start from the indexed sample.

Run: ``python generate_invalid_mzml.py``
"""
from __future__ import annotations

import hashlib
import random
import re
from pathlib import Path
from typing import Dict, List

import generate_synthetic_mzml as gen
import mzml

OUT = gen.OUT / "invalid"


def valid(name: str, seed: int = gen.SEED) -> bytes:
    """Rebuild one valid sample in memory, exactly as generate_synthetic_mzml.generate() writes it."""
    k = list(gen.FILES).index(name)
    return gen.FILES[name](random.Random(seed + k))


def replace_once(data: bytes, old: bytes, new: bytes) -> bytes:
    if data.count(old) < 1:
        raise ValueError(f"pattern not found: {old!r}")
    return data.replace(old, new, 1)


def fix_checksum(data: bytes) -> bytes:
    """Recompute an indexedmzML SHA-1 checksum after an edit, so only the intended defect remains."""
    tag = b"<fileChecksum>"
    end = data.rindex(tag) + len(tag)
    return data[:end] + hashlib.sha1(data[:end]).hexdigest().encode() + data[end + 40:]


def spectrum_count(seed: int) -> bytes:
    return replace_once(valid("synthetic_caffeine_profile.mzML", seed),
                        b'<spectrumList count="5"', b'<spectrumList count="6"')


def array_length(seed: int) -> bytes:
    data = valid("synthetic_caffeine_profile.mzML", seed)
    first = data.index(b'<spectrum index="0"')
    second_array = data.index(b"<binaryDataArray ", data.index(b"<binaryDataArray ", first) + 1)
    return data[:second_array] + b'<binaryDataArray arrayLength="1200" ' + data[second_array + len(b"<binaryDataArray "):]


def duplicate_id(seed: int) -> bytes:
    return replace_once(valid("synthetic_caffeine_profile.mzML", seed),
                        b'<spectrum index="1" id="scan=2"', b'<spectrum index="1" id="scan=1"')


def mz_order(seed: int) -> bytes:
    """Swap m/z points 10 and 11 of the first spectrum (32-bit, no compression)."""
    data = valid("synthetic_caffeine_profile.mzML", seed)
    start = data.index(b"<binary>", data.index(b'<spectrum index="0"')) + len(b"<binary>")
    end = data.index(b"</binary>", start)
    values = list(mzml.decode_binary(data[start:end].decode("ascii"), 32, False))
    values[10], values[11] = values[11], values[10]
    text = mzml.encode_binary(values, 32, False).encode("ascii")
    assert len(text) == end - start                          # same byte count, so encodedLength still holds
    return data[:start] + text + data[end:]


def index_offset(seed: int) -> bytes:
    data = valid("synthetic_caffeine_ms2_indexed.mzML", seed)
    m = re.search(rb'<offset idRef="scan=3">(\d+)</offset>', data)
    data = data[:m.start(1)] + str(int(m.group(1)) + 1).encode() + data[m.end(1):]
    return fix_checksum(data)


def checksum(seed: int) -> bytes:
    data = valid("synthetic_caffeine_ms2_indexed.mzML", seed)
    end = data.rindex(b"<fileChecksum>") + len(b"<fileChecksum>")
    digest = data[end:end + 40]
    wrong = hashlib.sha1(digest).hexdigest().encode()      # a well-formed SHA-1 of something else
    return data[:end] + wrong + data[end + 40:]


# file name -> (builder, expected finding code)
INVALID: Dict[str, tuple] = {
    "invalid_spectrum_count.mzML": (spectrum_count, "MZV010"),
    "invalid_array_length.mzML": (array_length, "MZV020"),
    "invalid_duplicate_id.mzML": (duplicate_id, "MZV012"),
    "invalid_mz_order.mzML": (mz_order, "MZV031"),
    "invalid_index_offset.mzML": (index_offset, "MZV050"),
    "invalid_checksum.mzML": (checksum, "MZV052"),
}


def generate(out_dir: Path = OUT, seed: int = gen.SEED) -> List[Path]:
    """Write every invalid sample to ``out_dir``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, (build, _) in INVALID.items():
        path = out_dir / name
        path.write_bytes(build(seed))                       # bytes: LF line endings on every OS
        paths.append(path)
    return paths


if __name__ == "__main__":
    for p in generate():
        print(f"wrote invalid/{p.name} ({p.stat().st_size:,} bytes)")
