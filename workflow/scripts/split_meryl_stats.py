#!/usr/bin/env python

import argparse as argp
import pathlib as pl
import re as re
import sys

import numpy as np
import pandas as pd
import scipy.signal as scisig


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
        "--forced-kmerfreq-threshold",
        "-f",
        type=int,
        default=-1,
        dest="forced_threshold",
        help="Record forced threshold if set externally. "
             "If a forced threshold is set, failing to "
             "determine a threshold using the data will "
             "not raise an exception. Default: -1 (none set)"
    )

    parser.add_argument(
        "--no-threshold-warning",
        action="store_true",
        default=False,
        dest="no_threshold_warning",
        help="If set, if no reliable k-mer threshold can be identified AND "
             "no threshold is forced, just print a warning instead of "
             "raising an exception (see previous option). Default: False"
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


def prettify_db_name(db_name):

    ext_to_split = [
        ".txt",
        ".tsv",
        ".stat",
        ".stats",
        ".statistic",
        ".statistics",
        ".meryl",
        ".meryl-stat",
        ".meryl-statistic",
        ".meryl-stats",
        ".meryl-statistics",
        ".meryl_stat",
        ".meryl_statistic",
        ".meryl_stats",
        ".meryl_statistics",
    ]

    prettified_name = db_name
    for ext in ext_to_split:
        if not prettified_name.endswith(ext):
            continue
        prettified_name = prettified_name.rsplit(".", 1)[0]
    if not prettified_name.strip():
        # some bogus name for testing?
        prettified_name = db_name
    return prettified_name


def determine_reliable_threshold(histogram):
    """Determine threshold value after which the slope
    over 'num_distinct_kmers' turns positive, i.e.,
    the leftmost index before the inflection point
    of the distribution marking the transition
    to more trustworthy/reliable k-mers.
    (cf. Merqury, Methods, section "K-mer completeness")

    Args:
        histogram (pandas.DataFrame):
    """

    dist_peaks, _ = scisig.find_peaks(
        histogram["num_distinct_kmers"].values,
        height=(None, None),
        distance=5
    )
    first_peak_index = dist_peaks[0]
    assert first_peak_index > 0
    # NB: histogram must be default indexed 0, 1, 2, ...
    first_peak_freq = int(histogram.loc[first_peak_index, "frequency"])

    hist_subset = histogram["num_distinct_kmers"].values[:first_peak_index]
    hist_gradient = np.gradient(hist_subset)

    inflection_point_index = None
    inflection_point_freq = None
    if any(hist_gradient <= 0):
        inflection_point_index = hist_gradient[hist_gradient <= 0].argmax()
        inflection_point_freq = int(histogram.loc[inflection_point_index, "frequency"])

    # NB: inflection point freq. can still be None,
    # e.g., for prefiltered data. That needs
    # to be checked in the calling scope.
    return inflection_point_freq, first_peak_freq


def main():

    args = parse_command_line()
    db_name = prettify_db_name(
        args.meryl_stats.name
    )
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

    df_hist = pd.DataFrame.from_records(histogram)

    min_kmer_frequency = df_hist["frequency"].min()
    assert min_kmer_frequency > 0
    statistics.append(("kmer_frequency_min", min_kmer_frequency))
    max_kmer_frequency = df_hist["frequency"].max()
    statistics.append(("kmer_frequency_max", max_kmer_frequency))

    reliable_kmer_threshold, first_peak_freq = determine_reliable_threshold(df_hist)
    statistics.append(("dist_first_peak_kmer_freq", first_peak_freq))
    threshold_is_forced = args.forced_threshold > -1
    # if this is a k-mer set that has already been filtered,
    # computing the boundary towards more reliable k-mers
    # will likely fail and that should be acceptable
    kmerset_is_filtered = min_kmer_frequency != 1

    if reliable_kmer_threshold is None:
        if threshold_is_forced:
            statistics.append(("kmer_reliable_greater", -1))
            statistics.append(("kmer_forced_threshold", args.forced_threshold))
        elif args.no_threshold_warning or kmerset_is_filtered:
            # no threshold id'ed, no threshold forced
            warn_msg = f"\nWarning: no 'reliable k-mer' threshold could be identified for {db_name}\n"
            sys.stderr.write(warn_msg)
            statistics.append(("kmer_reliable_greater", -1))
            statistics.append(("kmer_forced_threshold", args.forced_threshold))
        else:
            # no threshold id'ed, no threshold forced, and not just warning => raise
            err_msg = f"ERROR: no 'reliable k-mer' threshold could be identified for {db_name}"
            raise ValueError(err_msg)
    else:
        statistics.append(("kmer_reliable_greater", reliable_kmer_threshold))
        statistics.append(("kmer_forced_threshold", args.forced_threshold))

    df_stats = pd.DataFrame.from_records(statistics, columns=["statistic", "value"])
    df_stats["db_name"] = db_name
    df_stats.set_index(["db_name", "statistic"], inplace=True)

    # some sanity checks that could reveal missed
    # data entries in input
    stats_unique = df_stats.loc[(db_name, "unique_kmers")]["value"]
    stats_distinct = df_stats.loc[(db_name, "distinct_kmers")]["value"]
    stats_total = df_stats.loc[(db_name, "present_kmers")]["value"]

    try:
        histogram_unique = df_hist.loc[
            df_hist["frequency"] == 1, "num_distinct_kmers"
        ].values[0]
    except IndexError:
        assert stats_unique == 0
        histogram_unique = 0
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
