#!/usr/bin/env python3

import argparse as argp
import pathlib as pl

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--histogram", "-hist",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="histogram",
        help="Path to clean meryl histogram file.",
        required=True
    )

    parser.add_argument(
        "--statistics", "-stats",
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest="statistics",
        help="Path to histogram statistics file.",
        required=True
    )

    parser.add_argument(
        "--out-pdf", "-pdf",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="plot_pdf",
        required=False,
        help="Output path to histogram plot (PDF format)."
    )

    parser.add_argument(
        "--num-bins", "-b",
        type=int,
        default=100,
        dest="num_bins",
        help="Number of bins of the full histogram to plot. Default: 100"
    )

    args = parser.parse_args()

    return args


def compute_range_limits(histogram, num_bins):

    lowest_freq = int(histogram.at[0, "frequency"])
    highest_freq = int(histogram.loc[histogram.index[-1], "frequency"])

    bin_indices = np.arange(0, num_bins, dtype=int)
    bin_positions = bin_indices + 1

    assert highest_freq - lowest_freq > num_bins

    bin_labels = np.arange(lowest_freq, lowest_freq + num_bins, dtype=int)
    last_label = bin_labels[-1]
    bin_labels = list(map(str, bin_labels[:-1])) + [f"{last_label}+"]

    xlimit = (bin_indices.min()-1, bin_indices.max()+2)

    assert bin_indices.size == bin_positions.size == len(bin_labels)

    return bin_indices, bin_positions, bin_labels, xlimit


def extract_annotations(statistics):

    db_name = str(statistics.loc[0, "db_name"])

    select_stat = statistics["statistic"] == "dist_first_peak_kmer_freq"
    freq_peak = statistics.loc[select_stat, "value"].iloc[0]

    select_stat = statistics["statistic"] == "kmer_reliable_greater"
    freq_reliable = statistics.loc[select_stat, "value"].iloc[0]

    select_stat = statistics["statistic"] == "kmer_forced_threshold"
    freq_forced = statistics.loc[select_stat, "value"].iloc[0]

    return db_name, freq_peak, freq_reliable, freq_forced


def create_histogram_plot(statistics, histogram, num_bins):

    bin_indices, bar_positions, bar_labels, xlimit = compute_range_limits(
        histogram, num_bins
    )

    db_name, freq_peak, freq_reliable, freq_forced = extract_annotations(
        statistics
    )

    reduced_hist = np.zeros(num_bins, dtype=float)
    reduced_hist[:-1] = histogram.loc[bin_indices[:-1], "num_distinct_kmers"].values
    select_rest = histogram.index[bin_indices[-1]:]
    reduced_hist[-1] = histogram.loc[select_rest, "num_distinct_kmers"].sum()
    reduced_hist = np.log10(reduced_hist)

    fig, ax = plt.subplots(figsize=(10,8))

    ax.bar(
        bar_positions,
        reduced_hist,
        color="lightgrey",
        width=0.75
    )
    ax.set_title(db_name, fontsize=14)
    ax.set_xlim(xlimit)
    ax.set_xticks(bar_positions)
    # sparsify labels a bit for proper plotting
    sparse_x_labels = []
    for bl in bar_labels:
        try:
            if int(bl) % 5 == 0:
                sparse_x_labels.append(bl)
            else:
                sparse_x_labels.append("")
        except ValueError:  # last one
            sparse_x_labels.append(bl)
    ax.set_xticklabels(sparse_x_labels, fontsize=10)
    ax.set_xlabel("k-mer frequency", fontsize=12)
    ax.set_ylabel("count (log10)", fontsize=12)

    # note here: depending on the actual k-mer freqs. in the
    # dataset, there is not necessarily a 1-to-1 correspondence
    # between bar position and the frequency value at that position
    # (although this is the most likely case).
    # Hence, build a lookup table for convenience - exclude last position
    # because it refers to the rest of the dataset (label is, e.g., '100+')
    freq_pos_lut = dict((int(freq), pos) for freq, pos in zip(bar_labels[:-1], bar_positions[:-1]))

    peak_pos = freq_pos_lut[freq_peak]
    ax.axvline(
        peak_pos, 0, 0.99,
        color="black",
        ls="dashed",
        label=f"1st peak / freq. = {freq_peak}"
    )

    if freq_reliable > -1:
        inflection_point = freq_pos_lut[freq_reliable]
        # the reliable threshold is set for 'greater than'
        # filtering, hence make the inflection point visually
        # more appealing by setting it between this position
        # and the next (gradient already positive)
        ax.axvline(
            inflection_point+0.5, 0, 0.99,
            color="blue",
            ls="dashed",
            label=f"Inferred threshold\nkeep k-mer freq. > {freq_reliable}"
        )

    if freq_forced > -1:
        forced_pos = freq_pos_lut[freq_forced]
        ax.axvline(
            forced_pos+0.5, 0, 0.99,
            color="red",
            ls="dotted",
            label=f"Forced threshold\nkeep k-mer freq. > {freq_forced}"
        )

    ax.legend(loc="best", fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return fig, ax




def main():

    args = parse_command_line()

    stats = pd.read_csv(args.statistics, sep="\t", header=0)
    hist = pd.read_csv(args.histogram, sep="\t", header=0)

    fig, ax = create_histogram_plot(stats, hist, args.num_bins)

    if args.plot_pdf:
        args.plot_pdf.parent.mkdir(exist_ok=True, parents=True)
        plt.savefig(args.plot_pdf, transparent=None, facecolor="w", format="pdf")



if __name__ == "__main__":
    main()
