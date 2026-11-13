#!/usr/bin/env python3

from pathlib import Path
import sys

ROOT = Path("final_genotypes")


def read_raw(path):
    with open(path, "r") as f:
        lines = [line.rstrip("\n") for line in f]

    if not lines:
        raise RuntimeError(f"Empty file: {path}")

    header = lines[0].split()
    rows = [line.split() for line in lines[1:]]

    if len(header) < 6:
        raise RuntimeError(f"Invalid RAW header (<6 columns): {path}")

    for i, row in enumerate(rows, start=2):
        if len(row) != len(header):
            raise RuntimeError(
                f"Column number mismatch in {path}, line {i}: "
                f"expected {len(header)}, got {len(row)}"
            )

    return header, rows


def check_sample_columns(chip_rows, extra_rows, label):
    if len(chip_rows) != len(extra_rows):
        raise RuntimeError(
            f"{label}: sample number mismatch: "
            f"chip={len(chip_rows)}, extra={len(extra_rows)}"
        )

    for i, (a, b) in enumerate(zip(chip_rows, extra_rows), start=1):
        if a[:6] != b[:6]:
            raise RuntimeError(
                f"{label}: sample metadata/order mismatch at sample row {i}\n"
                f"chip  : {' '.join(a[:6])}\n"
                f"extra : {' '.join(b[:6])}"
            )


def make_combined(chip_path, extra_path, out_path):
    chip_header, chip_rows = read_raw(chip_path)

    # No QTN: combined is exactly the chip file.
    if extra_path is None:
        out_path.write_text(Path(chip_path).read_text())
        n_chip = len(chip_header) - 6
        return n_chip, 0, n_chip

    extra_header, extra_rows = read_raw(extra_path)

    check_sample_columns(chip_rows, extra_rows, str(out_path))

    chip_markers = chip_header[6:]
    extra_markers = extra_header[6:]

    chip_set = set(chip_markers)
    extra_set = set(extra_markers)

    duplicated = sorted(chip_set & extra_set)

    if duplicated:
        raise RuntimeError(
            f"Duplicate marker IDs between chip and extra in {out_path}:\n"
            + "\n".join(duplicated[:20])
        )

    combined_header = chip_header[:6] + chip_markers + extra_markers

    combined_rows = []
    for chip_row, extra_row in zip(chip_rows, extra_rows):
        combined_rows.append(
            chip_row[:6] +
            chip_row[6:] +
            extra_row[6:]
        )

    with open(out_path, "w") as f:
        f.write(" ".join(combined_header) + "\n")
        for row in combined_rows:
            f.write(" ".join(row) + "\n")

    return len(chip_markers), len(extra_markers), len(combined_header) - 6


def process_trait(outer, trait):
    d = ROOT / f"outer_{outer}" / trait

    chip_train = d / "chip_train.raw"
    chip_test = d / "chip_test.raw"

    extra_train = d / "extra_train.raw"
    extra_test = d / "extra_test.raw"

    combined_train = d / "combined_train.raw"
    combined_test = d / "combined_test.raw"

    if not chip_train.exists() or not chip_test.exists():
        print(f"WARNING: {trait} outer_{outer}: CHIP missing, skip")
        return

    if extra_train.exists() != extra_test.exists():
        raise RuntimeError(
            f"{trait} outer_{outer}: only one of extra_train/extra_test exists"
        )

    has_extra = extra_train.exists()

    if combined_train.exists() or combined_test.exists():
        raise RuntimeError(
            f"{trait} outer_{outer}: combined file already exists. "
            f"Remove it first if you want to regenerate."
        )

    if has_extra:
        n_chip_train, n_extra_train, n_comb_train = make_combined(
            chip_train, extra_train, combined_train
        )

        n_chip_test, n_extra_test, n_comb_test = make_combined(
            chip_test, extra_test, combined_test
        )

        if n_chip_train != n_chip_test:
            raise RuntimeError(
                f"{trait} outer_{outer}: chip marker number "
                f"train/test mismatch"
            )

        if n_extra_train != n_extra_test:
            raise RuntimeError(
                f"{trait} outer_{outer}: extra marker number "
                f"train/test mismatch"
            )

        if n_comb_train != n_comb_test:
            raise RuntimeError(
                f"{trait} outer_{outer}: combined marker number "
                f"train/test mismatch"
            )

        print(
            f"{trait} outer_{outer}: "
            f"CHIP={n_chip_train}, "
            f"QTN={n_extra_train}, "
            f"COMBINED={n_comb_train}"
        )

    else:
        n_chip_train, n_extra_train, n_comb_train = make_combined(
            chip_train, None, combined_train
        )

        n_chip_test, n_extra_test, n_comb_test = make_combined(
            chip_test, None, combined_test
        )

        print(
            f"{trait} outer_{outer}: "
            f"CHIP={n_chip_train}, "
            f"QTN=0, "
            f"COMBINED={n_comb_train} "
            f"(same as CHIP)"
        )


for outer in [f"{i:02d}" for i in range(1, 11)]:
    for trait in ["BB", "H6"]:
        process_trait(outer, trait)

print()
print("===== ALL COMBINED FILES CREATED SUCCESSFULLY =====")
