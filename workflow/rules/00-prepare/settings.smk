
# CONFIG SETTINGS

PROCESS_SINGLETONS = config.get("process_singletons", True)
PROCESS_TRIOS = config.get("process_trios", True)

COMPARE_PAIRWISE_BY_FILE = config.get("compare_pairwise_by_file", False)
COMPARE_PAIRWISE_BY_SAMPLE = config.get("compare_pairwise_by_sample", False)

DISCARD_KMER_DATABASES = config.get("discard_kmer_databases", False)

RUN_MERYL = config.get("run_meryl", True)

FORCE_KMERFREQ_THRESHOLD = config.get("force_kmerfreq_threshold", -1)
assert isinstance(FORCE_KMERFREQ_THRESHOLD, int)

# TOOL SETTINGS

#======================
# MERYL k-mer databases
#======================

# produce homopolymer-compressed
# output yes/no/both
MERYL_COMPRESS_CONFIG_OPTION = config["meryl_compress"]
assert isinstance(MERYL_COMPRESS_CONFIG_OPTION, int)

MAP_MERYL_COMPRESS_CONFIG_OPTION = {
    1: ["is-hpc"],
    0: ["is-hpc", "no-hpc"],
    -1: ["no-hpc"]
}

# this is then simply a list of wildcard
# replacement values to be applied during
# expand() operations by Snakemake
MERYL_COMPRESS_WILDCARD_VALUES = \
    MAP_MERYL_COMPRESS_CONFIG_OPTION[MERYL_COMPRESS_CONFIG_OPTION]

MERYL_KMER_VALUES = config["meryl_kmer_size"]
assert all(isinstance(value, int) for value in MERYL_KMER_VALUES)
