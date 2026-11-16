"""
Build yak k-mer databases per input sample.

Source: https://github.com/chhylp123/hifiasm

# Trio binning assembly (requiring https://github.com/lh3/yak)
yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
yak count -b37 -t16 -o mat.yak <(cat mat_1.fq.gz mat_2.fq.gz) <(cat mat_1.fq.gz mat_2.fq.gz)

Notes on yak count:
- with a Bloom filter (-b > 0), yak performs two passes over the
  input (pass 1: populate Bloom filter; pass 2: count k-mers seen
  at least twice, i.e. singletons are discarded). If the input is
  a regular file, yak re-opens it for the second pass. If the
  input is a stream (process substitution), the same stream has
  to be given twice.
- yak has no operation to merge/union databases. Hence, samples
  with several input files are counted in one go by concatenating
  (decompressed) input files instead of counting per part and
  merging afterwards (as done for meryl).
- yak does not support homopolymer-compressed counting; there is
  no {hpc} wildcard for yak databases.
"""

#=====================
# YAK k-mer databases
# (module-local settings; candidates for 00-prepare/settings.smk)
#=====================

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


def get_yak_sample_input_files(wildcards):
    """Collect all input files for a sample,
    irrespective of the sample being a single-
    or multi-input sample. Files are returned
    in a stable (sorted) order.
    """
    sample_files = sorted(
        str(input_file) for (sample, part), input_file
        in MAP_SAMPLE_TO_INPUT_FILE.items()
        if sample == wildcards.sample
    )
    if not sample_files:
        raise ValueError(f"No input files for sample: {wildcards.sample}")
    return sample_files


def make_yak_count_input_spec(input_files, bloom_bits):
    """Build the input part of the 'yak count'
    command line. A single (regular) file can be
    re-read by yak for the second (Bloom filter)
    pass; several files are streamed via process
    substitution, which requires repeating the
    stream if the Bloom filter is active.
    'gzip -dcf' handles compressed and plain
    files alike, avoiding broken streams when
    plain and compressed inputs are mixed.
    """
    if len(input_files) == 1:
        return f"{input_files[0]}"
    file_list = " ".join(input_files)
    stream = f"<(gzip -dcf {file_list})"
    if bloom_bits > 0:
        return f"{stream} {stream}"
    return stream


_YAK_RESULT_DB = DIR_RES.joinpath("databases", "singletons", "{sample}.k{size_k}.yak")


rule yak_build_kmer_db:
    """
    yak databases are already a compact binary format,
    hence they are written directly to the results folder
    (no compression step as for meryl). If the databases
    shall be discarded, the output is marked as temp (it is
    still needed to produce the k-mer histogram) and is not
    registered in the results manifest (the checksum files
    of the manifest would otherwise require the database).
    """
    input:
        seq = get_yak_sample_input_files
    output:
        db = temp(_YAK_RESULT_DB) if DISCARD_KMER_DATABASES else _YAK_RESULT_DB
    log:
        DIR_LOG.joinpath("10-count", "yak", "10-build", "{sample}.k{size_k}.yak.log")
    benchmark:
        DIR_RSRC.joinpath("10-count", "yak", "10-build", "{sample}.k{size_k}.yak.rsrc")
    wildcard_constraints:
        sample = "(" + "|".join(map(re.escape, map(str, SAMPLES))) + ")",
        size_k = "[0-9]+"
    threads: CPU_MEDIUM
    resources:
        # -b37 => Bloom filter of 2^37 bits = 16 GiB,
        # plus the k-mer hash table
        mem_mb = lambda wildcards, attempt: (
            max(0, 2 ** YAK_BLOOM_FILTER_BITS // 8 // 1024 ** 2) + 32768
        ) * attempt,
        time_hrs = lambda wildcards, attempt: 4 * attempt,
    conda:
        DIR_ENVS.joinpath("yak.yaml")
    params:
        bloom = YAK_BLOOM_FILTER_BITS,
        input_spec = lambda wildcards, input: make_yak_count_input_spec(
            list(input.seq), YAK_BLOOM_FILTER_BITS
        ),
        acc_in = lambda wildcards, input: register_input(input.seq, allow_non_existing=True),
        acc_out = lambda wildcards, output: (
            None if DISCARD_KMER_DATABASES else register_result(output.db)
        ),
    shell:
        "yak count -k {wildcards.size_k} -b {params.bloom} -t {threads} "
            "-o {output.db} {params.input_spec} &> {log}"


rule yak_dump_kmer_histogram:
    """
    'yak inspect' on a single database prints the
    k-mer frequency histogram (lines prefixed with 'HS')
    """
    input:
        db = _YAK_RESULT_DB
    output:
        hist = DIR_RES.joinpath(
            "statistics", "singletons", "{sample}.k{size_k}.yak-hist.tsv"
        )
    log:
        DIR_LOG.joinpath("10-count", "yak", "20-hist", "{sample}.k{size_k}.yak-hist.log")
    wildcard_constraints:
        size_k = "[0-9]+"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("yak.yaml")
    params:
        acc_out = lambda wildcards, output: register_result(output.hist),
    shell:
        "yak inspect {input.db} > {output.hist} 2> {log}"


_SINGLETONS_RESULT_YAK_DB = expand(
    _YAK_RESULT_DB,
    sample=SAMPLES,
    size_k=YAK_KMER_VALUES,
)
if DISCARD_KMER_DATABASES:
    _SINGLETONS_RESULT_YAK_DB = []


rule yak_run_singletons:
    input:
        dbs = _SINGLETONS_RESULT_YAK_DB,
        hists = expand(
            rules.yak_dump_kmer_histogram.output.hist,
            sample=SAMPLES,
            size_k=YAK_KMER_VALUES,
        ),


YAK_WORKFLOW_OUTPUT = []
if RUN_YAK and PROCESS_SINGLETONS:
    YAK_WORKFLOW_OUTPUT.extend(rules.yak_run_singletons.input)
