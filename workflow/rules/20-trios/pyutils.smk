import pathlib


def load_reliable_kmer_threshold(tsv_file):

    if isinstance(tsv_file, list):
        # can happen if expand() is used to
        # dynamically derive rule input
        assert len(tsv_file) == 1
        tsv_file = tsv_file[0]

    if FORCE_KMERFREQ_THRESHOLD > -1:
        kmer_t = FORCE_KMERFREQ_THRESHOLD
    elif not pathlib.Path(tsv_file).is_file():
        # file does not exist yet, e.g., during
        # a Snakemake dry run. Return dummy
        # value of zero
        kmer_t = 0
    else:
        kmer_t = None
        assert pathlib.Path(tsv_file).name.endswith(".tsv")
        with open(tsv_file, "r") as table:
            for line in table:
                if not line.strip() or line.startswith("#"):
                    continue
                db_name, statistic, value = line.strip().split("\t")
                if statistic == "kmer_reliable_greater":
                    kmer_t = int(value)
                    break
        if kmer_t is None:
            raise ValueError(
                f"Cannot extract k-mer reliability threshold from file: {tsv_file}"
            )
    return kmer_t
