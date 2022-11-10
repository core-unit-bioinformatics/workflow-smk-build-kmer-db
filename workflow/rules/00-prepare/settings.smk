
# CONFIG SETTINGS

PROCESS_SINGLETONS = config.get("process_singletons", True)
PROCESS_TRIOS = config.get("process_trios", True)

RUN_MERYL = config.get("run_meryl", True)

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
