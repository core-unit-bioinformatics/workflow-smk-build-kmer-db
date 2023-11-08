import pathlib


def load_reliable_kmer_threshold(tsv_file):

    if isinstance(tsv_file, list):
        # can happen if expand() is used to
        # dynamically derive rule input
        assert len(tsv_file) == 1
        tsv_file = tsv_file[0]

    # case always force threshold, no matter what's
    # in the stats output file
    if FORCE_KMERFREQ_THRESHOLD_ALWAYS:
        assert FORCE_KMERFREQ_THRESHOLD > 0, (
            "Always enforcing a k-mer frequency threshold "
            "requires a positive integer to be specified "
            "as threshold in the config. You specified: "
            f"FORCE_KMERFREQ_THRESHOLD = {FORCE_KMERFREQ_THRESHOLD}"
        )
        kmer_t = FORCE_KMERFREQ_THRESHOLD
    elif not pathlib.Path(tsv_file).is_file():
        # file does not exist yet, e.g., during
        # a Snakemake dry run. Return dummy
        # value of zero
        kmer_t = 0
    else:
        # not always enforcing and the stats file exists
        # implies we either load the data-derived threshold
        # or the enforced threshold
        data_derived_t = -1
        enforced_t = -1
        assert pathlib.Path(tsv_file).name.endswith(".tsv")
        with open(tsv_file, "r") as table:
            for line in table:
                if not line.strip() or line.startswith("#"):
                    continue
                db_name, statistic, value = line.strip().split("\t")
                if statistic == "kmer_reliable_greater":
                    data_derived_t = int(value)
                if statistic == "kmer_forced_threshold":
                    enforced_t = int(value)
        if data_derived_t != -1:
            kmer_t = data_derived_t
        elif enforced_t != -1:
            kmer_t = enforced_t
        else:
            kmer_t = FORCE_KMERFREQ_THRESHOLD
        if kmer_t < 0:
            raise ValueError(
                f"Cannot extract k-mer reliability threshold from file: {tsv_file}\n"
                f"Data-derived threshold: {data_derived_t}\n"
                f"Enforced threshold loaded from file: {enforced_t}\n"
                f"Enforced threshold loaded from workflow config: {FORCE_KMERFREQ_THRESHOLD}\n"
            )
    return kmer_t
