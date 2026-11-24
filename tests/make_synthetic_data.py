"""Regenerate the tiny synthetic inputs in tests/data/.

The files are random sequence, not biological data. They exist only so that
`snakemake --dry-run` can resolve rule inputs. Output is deterministic
(fixed seed, gzip mtime=0).

Usage: python3 tests/make_synthetic_data.py [output_dir]
"""

import gzip
import pathlib
import random
import sys

READ_LENGTH = 100
FRAGMENT_LENGTH = 300
N_PAIRS = 25
REFERENCE_LENGTH = 2000


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main(out_dir: pathlib.Path) -> None:
    rng = random.Random(0)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference = "".join(rng.choice("ACGT") for _ in range(REFERENCE_LENGTH))
    wrapped = "\n".join(reference[i : i + 70] for i in range(0, len(reference), 70))
    (out_dir / "reference.fna").write_text(
        ">synthetic_reference_1 synthetic test sequence, not a real genome\n"
        + wrapped
        + "\n"
    )

    forward, reverse = [], []
    quality = "I" * READ_LENGTH
    for i in range(N_PAIRS):
        start = rng.randint(0, REFERENCE_LENGTH - FRAGMENT_LENGTH)
        fragment = reference[start : start + FRAGMENT_LENGTH]
        forward.append(f"@synthetic_read_{i}/1\n{fragment[:READ_LENGTH]}\n+\n{quality}\n")
        reverse.append(
            f"@synthetic_read_{i}/2\n{reverse_complement(fragment[-READ_LENGTH:])}\n+\n{quality}\n"
        )

    for name, records in (("sample_1.fastq.gz", forward), ("sample_2.fastq.gz", reverse)):
        with open(out_dir / name, "wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
                handle.write("".join(records).encode())


if __name__ == "__main__":
    default = pathlib.Path(__file__).resolve().parent / "data"
    main(pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else default)
