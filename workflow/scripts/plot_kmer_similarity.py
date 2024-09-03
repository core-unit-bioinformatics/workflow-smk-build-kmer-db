#!/usr/bin/env python3

import argparse as argp

import pandas as pd
import pathlib as pl
import matplotlib as mpl
import numpy as np
import seaborn as sns

mpl.use('Agg')
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--input", "-i",
        type=lambda x: pl.Path(x).resolve(strict=True),
        nargs="+",
        dest="input",
        help="Meryl stats dump files or single folder containing such files."
    )

    parser.add_argument(
        "--output", "-o",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="output",
        default="jaccard_similarity.pdf",
        help="Path to output PDF. Default: jaccard_similarity.pdf"
    )

    args = parser.parse_args()

    return args


def collect_input_files(folder):

    if not folder.is_dir():
        raise ValueError("Path is not a folder")

    input_files = []
    for tsv_file in folder.glob("*.meryl-stats.tsv"):
        if "1-and-2" in tsv_file.name or "1-or-2" in tsv_file.name:
            input_files.append(tsv_file)

    if not input_files:
        raise FileNotFoundError(f"No meryl stats tsv files collected from path: {folder}")

    return sorted(input_files)


def main():

    args = parse_command_line()

    if len(args.input) == 1:
        # must be folder
        input_files = collect_input_files(args.input[0])
    else:
        input_files = args.input

    if len(input_files) < 3:
        raise RuntimeError(f"Less than 3 input files provided - aborting cluster map creation: {input_files}")

    shared = dict()
    union = dict()
    samples = set()

    db_parameters = set()

    for stats_file in input_files:
        df = pd.read_csv(stats_file, sep="\t", header=0)
        db_name, db_op = df["db_name"].iloc[0].rsplit(".", 1)
        db1, db2 = db_name.split("_vs_")
        sample1, kms1, hpc1 = db1.rsplit(".", 2)
        sample2, kms2, hpc2 = db2.rsplit(".", 2)
        # collect k-mer size (kms) and
        # homopolymer compression parameters
        # to assert that they are identical
        # over all samples
        db_parameters.update({kms1, kms2, hpc1, hpc2})
        distinct_count = df.loc[df["statistic"] == "distinct_kmers", "value"].values[0]
        assert sample1 < sample2
        samples.add(sample1)
        samples.add(sample2)
        if db_op == "1-and-2":
            shared[(sample1, sample2)] = distinct_count
        elif db_op == "1-or-2":
            union[(sample1, sample2)] = distinct_count
        else:
            raise ValueError(f"Unknown DB op: {db_op}")

    if len(db_parameters) != 2:
        raise RuntimeError(f"Incoherent database parameters: {db_parameters}")

    num_samples = len(samples)

    jaccard = pd.DataFrame(
        np.zeros((num_samples, num_samples), dtype=float),
        index=sorted(samples),
        columns=sorted(samples)
    )

    for (s1, s2) in shared.keys():
        n_union = union[(s1, s2)]
        n_shared = shared[(s1, s2)]
        jacc = round(n_shared / n_union, 3)
        jaccard.loc[s1, s2] = jacc
        jaccard.loc[s2, s1] = jacc
        jaccard.loc[s1, s1] = 1.
        jaccard.loc[s2, s2] = 1.

    cm = sns.clustermap(jaccard, cmap="vlag", **{"annot": True})
    cm.figure.suptitle(f"Jaccard similarity / k={kms1[1:]}", fontsize=14, y=1.02,x=0.5)

    args.output.parent.mkdir(exist_ok=True, parents=True)
    cm.savefig(args.output)

    return 0


if __name__ == "__main__":
    main()
