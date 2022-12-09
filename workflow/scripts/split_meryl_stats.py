#!/usr/bin/env python

import argparse as argp
import pathlib as pl
import re as re

import pandas as pd


HISTOGRAM_HEADER = [
    "frequency",
    "num_distinct_kmers",
    "cum_fraction_distinct",
    "cum_fraction_total",
    "presence_in_dataset",
]

HISTOGRAM_VALUE_TYPES = [int, int, float, float, float]


def parse_command_line():

    parser = argp.ArgumentParser(add_help=True)
    parser.add_argument(
        "--meryl-stats",
        "-s",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="meryl_stats",
        help="Output file of meryl statistics subcommand.",
    )

    parser.add_argument(
        "--out-stats",
        "-os",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_stats",
        help="Path to output file for statistics.",
    )

    parser.add_argument(
        "--out-hist",
        "-oh",
        type=lambda x: pl.Path(x).resolve(),
        dest="out_hist",
        help="Path to output file for histogram.",
    )

    args = parser.parse_args()

    return args


def parse_histogram_line(line):
    """
    Parse a histogram line and return
    typed info record as dictionary.

    Args:
        line (str): histogram line from meryl stats output

    Returns:
        bool|None, dict|None: keep record, data record
    """

    if line.startswith("---"):
        info_good = None
        record = None
    else:
        info_good = True
        values = [entry.strip() for entry in line.split()]
        assert len(values) == len(HISTOGRAM_HEADER), f"Value length mismatch: {values}"

        record = dict(
            (h, t(v))
            for h, t, v in zip(HISTOGRAM_HEADER, HISTOGRAM_VALUE_TYPES, values)
        )

    return info_good, record


def parse_statistics_line(line):
    """_summary_

    Args:
        line (str): statistics line from meryl stats output

    Raises:
        ValueError: raised if line cannot be parsed

    Returns:
        bool|None, tuple|None: use record, (k-mer count statistic, value)
    """

    if line.startswith("frequency"):
        info_good = False
        info = None
        assert "(1e-6)" in line, f"Counts not scaled to million: {line}"
    elif any(x in line for x in ["fraction", "cumulative"]):
        info_good = None
        info = None
    else:
        info_good = True
        mobj = re.search("[0-9]+", line)
        assert mobj is not None, f"Cannot extract number from line: {line}"
        extracted_number = int(mobj.group(0))
        if line.startswith("Number of"):
            info = "kmer_size", extracted_number
        elif line.startswith("unique"):
            info = "unique_kmers", extracted_number
        elif line.startswith("distinct"):
            info = "distinct_kmers", extracted_number
        elif line.startswith("present"):
            info = "present_kmers", extracted_number
        elif line.startswith("missing"):
            info = "missing_kmers", extracted_number
        else:
            raise ValueError(f"Cannot process statistics line: {line}")
    return info_good, info


def create_output_path(file_path):
    file_path.parent.mkdir(exist_ok=True, parents=True)
    return file_path


def main():

    args = parse_command_line()

    # collect records for count statistics
    # and k-mer histogram
    statistics = []
    histogram = []

    # meryl statistic output file starts with
    # count statistics, followed by histogram
    line_parser = parse_statistics_line
    records = statistics

    with open(args.meryl_stats, "r") as stats_out:
        for line in stats_out:
            if not line.strip():
                continue
            info_good, info = line_parser(line.strip())
            if info_good is None:
                continue
            if not info_good:
                line_parser = parse_histogram_line
                records = histogram
                continue
            records.append(info)

    df_stats = pd.DataFrame.from_records(statistics, columns=["statistic", "value"])
    df_stats.set_index("statistic", inplace=True)

    df_hist = pd.DataFrame.from_records(histogram)

    # some sanity checks that could reveal missed
    # data entries in input
    stats_unique = df_stats.loc["unique_kmers"]["value"]
    stats_distinct = df_stats.loc["distinct_kmers"]["value"]
    stats_total = df_stats.loc["present_kmers"]["value"]

    histogram_unique = df_hist.loc[
        df_hist["frequency"] == 1, "num_distinct_kmers"
    ].values[0]
    histogram_distinct = df_hist["num_distinct_kmers"].sum()
    histogram_total = (df_hist["num_distinct_kmers"] * df_hist["frequency"]).sum()

    same_unique = stats_unique == histogram_unique
    same_distinct = stats_distinct == histogram_distinct
    same_total = stats_total == histogram_total

    assert same_unique and same_distinct and same_total, (
        f"Not all counts identical:\n"
        f"{stats_unique} - {histogram_unique}\n"
        f"{stats_distinct} - {histogram_distinct}\n"
        f"{stats_total} - {histogram_total}\n"
    )

    df_stats.to_csv(
        create_output_path(args.out_stats), sep="\t", index=True, header=True
    )
    df_hist.to_csv(
        create_output_path(args.out_hist), sep="\t", index=False, header=True
    )

    return


if __name__ == "__main__":
    main()
