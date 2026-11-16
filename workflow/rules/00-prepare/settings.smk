
# CONFIG SETTINGS

# input data file paths read from fofn
# files may come as a relative path;
# prepend DATA_ROOT if path does not
# resolve
DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

PROCESS_SINGLETONS = config.get("process_singletons", True)
PROCESS_TRIOS = config.get("process_trios", True)

COMPARE_PAIRWISE_BY_FILE = config.get("compare_pairwise_by_file", False)
# 2024-09-02
# temporary debug setting
# comparisons by file are currently not support
if COMPARE_PAIRWISE_BY_FILE:
    logerr("k-mer database comparisons by file are currently not supported - deactivating switch")
    COMPARE_PAIRWISE_BY_FILE = False

COMPARE_PAIRWISE_BY_SAMPLE = config.get("compare_pairwise_by_sample", False)

DOWNSAMPLE_COMPARISONS = config.get("downsample_comparisons", 0.)
assert isinstance(DOWNSAMPLE_COMPARISONS, float)
assert 0. <= DOWNSAMPLE_COMPARISONS < 1

DISCARD_KMER_DATABASES = config.get("discard_kmer_databases", False)

FORCE_KMERFREQ_THRESHOLD = config.get("force_kmerfreq_threshold", -1)
assert isinstance(FORCE_KMERFREQ_THRESHOLD, int)
FORCE_KMERFREQ_THRESHOLD_ALWAYS = config.get("force_kmerfreq_threshold_always", False)
assert isinstance(FORCE_KMERFREQ_THRESHOLD_ALWAYS, bool)

# TOOL SETTINGS

#======================
# MERYL k-mer databases
#======================

RUN_MERYL = config.get("run_meryl", False)

# produce homopolymer-compressed
# output yes/no/both
MERYL_COMPRESS_CONFIG_OPTION = config.get("meryl_compress", -1)
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

MERYL_KMER_VALUES = config.get("meryl_kmer_size", [21])
assert all(isinstance(value, int) for value in MERYL_KMER_VALUES)


#======================
# YAK k-mer databases
#======================

RUN_YAK = config.get("run_yak", False)
assert isinstance(RUN_YAK, bool)

# hifiasm recommends k=31 (yak default)
YAK_KMER_VALUES = config.get("yak_kmer_size", [31])
if isinstance(YAK_KMER_VALUES, int):
    YAK_KMER_VALUES = [YAK_KMER_VALUES]
assert all(isinstance(value, int) for value in YAK_KMER_VALUES)

# size of the Bloom filter in log2(bits);
# 0 deactivates the Bloom filter (singleton
# k-mers are then kept)
YAK_BLOOM_FILTER_BITS = config.get("yak_bloom_filter_bits", 37)
assert isinstance(YAK_BLOOM_FILTER_BITS, int)
assert 0 <= YAK_BLOOM_FILTER_BITS < 64

