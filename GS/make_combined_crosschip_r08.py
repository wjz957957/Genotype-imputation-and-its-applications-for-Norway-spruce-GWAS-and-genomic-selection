#!/usr/bin/env python3
from pathlib import Path
from itertools import zip_longest
import shutil
import sys

ROOT = Path("07_final_genotypes_crosschip_r2_0p8")
SENTINEL = object()


def read_header(path):
    with open(path, "r") as f:
        line = f.readline()
    if not line:
        raise RuntimeError(f"Empty file: {path}")
    header = line.split()
    if len(header) < 7:
        raise RuntimeError(f"Invalid RAW header (<7 columns): {path}")
    if len(set(header[6:])) != len(header[6:]):
        raise RuntimeError(f"Duplicated marker names in header: {path}")
    return header


def check_header_compatibility(chip_train, chip_test, extra_train=None, extra_test=None):
    if chip_train[:6] != chip_test[:6]:
        raise RuntimeError("CHIP train/test metadata column names differ")
    if chip_train[6:] != chip_test[6:]:
        raise RuntimeError("CHIP train/test marker names/order differ")

    if extra_train is not None:
        if extra_train[:6] != extra_test[:6]:
            raise RuntimeError("EXTRA train/test metadata column names differ")
        if extra_train[6:] != extra_test[6:]:
            raise RuntimeError("EXTRA train/test marker names/order differ")
        if chip_train[:6] != extra_train[:6]:
            raise RuntimeError("CHIP and EXTRA metadata column names differ")

        overlap = set(chip_train[6:]).intersection(extra_train[6:])
        if overlap:
            examples = "\n".join(sorted(overlap)[:20])
            raise RuntimeError(
                f"CHIP/EXTRA marker overlap ({len(overlap)} markers):\n{examples}"
            )


def combine_streaming(chip_path, extra_path, out_path):
    chip_header = read_header(chip_path)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    if extra_path is None:
        shutil.copyfile(chip_path, tmp_path)
        tmp_path.replace(out_path)
        return len(chip_header) - 6, 0, len(chip_header) - 6

    extra_header = read_header(extra_path)
    if chip_header[:6] != extra_header[:6]:
        raise RuntimeError(f"Metadata column names differ: {chip_path} vs {extra_path}")

    chip_n = len(chip_header) - 6
    extra_n = len(extra_header) - 6

    with open(chip_path, "r") as cf, open(extra_path, "r") as ef, open(tmp_path, "w") as out:
        # Headers have already been checked by the caller.
        cf.readline()
        ef.readline()
        out.write(" ".join(chip_header[:6] + chip_header[6:] + extra_header[6:]) + "\n")

        row_num = 1
        for c_line, e_line in zip_longest(cf, ef, fillvalue=SENTINEL):
            row_num += 1
            if c_line is SENTINEL or e_line is SENTINEL:
                raise RuntimeError(
                    f"Row count mismatch while combining {chip_path} and {extra_path}"
                )

            c_row = c_line.split()
            e_row = e_line.split()

            if len(c_row) != len(chip_header):
                raise RuntimeError(
                    f"Column mismatch in {chip_path}, line {row_num}: "
                    f"expected {len(chip_header)}, got {len(c_row)}"
                )
            if len(e_row) != len(extra_header):
                raise RuntimeError(
                    f"Column mismatch in {extra_path}, line {row_num}: "
                    f"expected {len(extra_header)}, got {len(e_row)}"
                )
            if c_row[:6] != e_row[:6]:
                raise RuntimeError(
                    f"Sample metadata/order mismatch at data row {row_num - 1}:\n"
                    f"CHIP : {' '.join(c_row[:6])}\n"
                    f"EXTRA: {' '.join(e_row[:6])}"
                )

            out.write(" ".join(c_row[:6] + c_row[6:] + e_row[6:]) + "\n")

    tmp_path.replace(out_path)
    return chip_n, extra_n, chip_n + extra_n


def process_trait(outer, trait):
    d = ROOT / f"outer_{outer}" / trait

    chip_train = d / "chip_train.raw"
    chip_test = d / "chip_test.raw"
    extra_train = d / "extra_train.raw"
    extra_test = d / "extra_test.raw"
    combined_train = d / "combined_train.raw"
    combined_test = d / "combined_test.raw"

    if not chip_train.exists() or not chip_test.exists():
        raise RuntimeError(f"Missing CHIP train/test file in {d}")

    has_extra_train = extra_train.exists()
    has_extra_test = extra_test.exists()
    if has_extra_train != has_extra_test:
        raise RuntimeError(
            f"{trait} outer_{outer}: only one of extra_train.raw/extra_test.raw exists"
        )

    chip_train_header = read_header(chip_train)
    chip_test_header = read_header(chip_test)

    if has_extra_train:
        extra_train_header = read_header(extra_train)
        extra_test_header = read_header(extra_test)
        check_header_compatibility(
            chip_train_header, chip_test_header,
            extra_train_header, extra_test_header
        )
    else:
        check_header_compatibility(chip_train_header, chip_test_header)

    # Do not overwrite pre-existing combined data silently.
    if combined_train.exists() or combined_test.exists():
        raise RuntimeError(
            f"Combined file already exists for {trait} outer_{outer}. "
            f"Inspect it first; remove both combined files manually only if regeneration is intended."
        )

    if has_extra_train:
        n_chip_train, n_extra_train, n_comb_train = combine_streaming(
            chip_train, extra_train, combined_train
        )
        n_chip_test, n_extra_test, n_comb_test = combine_streaming(
            chip_test, extra_test, combined_test
        )

        if (n_chip_train, n_extra_train, n_comb_train) != (
            n_chip_test, n_extra_test, n_comb_test
        ):
            raise RuntimeError(f"Train/test marker counts differ for {trait} outer_{outer}")

        print(
            f"{trait} outer_{outer}: CHIP={n_chip_train}, "
            f"QTN={n_extra_train}, COMBINED={n_comb_train}"
        )
    else:
        n_chip_train, _, n_comb_train = combine_streaming(
            chip_train, None, combined_train
        )
        n_chip_test, _, n_comb_test = combine_streaming(
            chip_test, None, combined_test
        )
        if n_chip_train != n_chip_test or n_comb_train != n_comb_test:
            raise RuntimeError(f"Train/test marker counts differ for {trait} outer_{outer}")
        print(
            f"{trait} outer_{outer}: CHIP={n_chip_train}, QTN=0, "
            f"COMBINED={n_comb_train} (identical to CHIP)"
        )

    # Final strict check: the generated train/test marker headers must match.
    final_train_header = read_header(combined_train)
    final_test_header = read_header(combined_test)
    if final_train_header != final_test_header:
        raise RuntimeError(f"Combined train/test marker headers differ: {trait} outer_{outer}")


def main():
    if not ROOT.is_dir():
        sys.exit(f"ERROR: input directory not found: {ROOT}")

    print(f"Input/output root: {ROOT}")
    for outer in [f"{i:02d}" for i in range(1, 11)]:
        for trait in ("BB", "H6"):
            process_trait(outer, trait)

    print("\n===== ALL COMBINED FILES CREATED SUCCESSFULLY =====")


if __name__ == "__main__":
    main()
