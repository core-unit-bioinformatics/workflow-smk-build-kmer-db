
rule meryl_dump_part_file_dbs:
    input:
        seq = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILE[(wildcards.sample, wildcards.part)]
    output:
        db = temp(directory(
            DIR_PROC.joinpath("10-singles", "meryl", "parts", "{sample}.{part}.k{size_k}.{hpc}.meryl")
        ))
    log:
        DIR_LOG.joinpath("10-singles", "meryl", "parts", "{sample}.{part}.k{size_k}.{hpc}.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("10-singles", "meryl", "parts", "{sample}.{part}.k{size_k}.{hpc}.meryl.rsrc")
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES_MULTI_INPUT
    threads: CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_gb = lambda wildcards, attempt: 32 * attempt,
        time_hrs = lambda wildcards, attempt: 4 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    params:
        compress = lambda wildcards: "compress" if wildcards.hpc == "is-hpc" else "",
        acc_in=lambda wildcards, input: register_input(input.seq, allow_non_existing=True),
    shell:
        "meryl count k={wildcards.size_k} memory={resources.mem_gb} threads={threads} "
            "{params.compress} {input.seq} output {output} &> {log}"


rule meryl_union_part_file_dbs:
    input:
        dbs = lambda wildcards:
            expand(
                DIR_PROC.joinpath("10-singles", "meryl", "parts", "{{sample}}.{part}.k{{size_k}}.{{hpc}}.meryl"),
                part=MAP_SAMPLE_TO_PART_IDS[wildcards.sample]
            )
    output:
        db = temp(directory(
            DIR_PROC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl")
        ))
    log:
        DIR_LOG.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl.rsrc")
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES_MULTI_INPUT
    threads: CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_gb = lambda wildcards, attempt: 32 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl union-sum {input.dbs} output {output.db} &> {log}"


rule meryl_dump_kmers:
    input:
        seq = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILE[(wildcards.sample, None)]
    output:
        db = temp(directory(
            DIR_PROC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl")
        ))
    log:
        DIR_LOG.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl.rsrc")
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES_SINGLE_INPUT
    threads: CPU_MEDIUM
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_gb = lambda wildcards, attempt: 32 * attempt,
        time_hrs = lambda wildcards, attempt: 4 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    params:
        compress = lambda wildcards: "compress" if wildcards.hpc == "is-hpc" else "",
        acc_in=lambda wildcards, input: register_input(input.seq, allow_non_existing=True),
    shell:
        "meryl count k={wildcards.size_k} memory={resources.mem_gb} threads={threads} "
            "{params.compress} {input.seq} output {output.db} &> {log}"


_SINGLES_RESULT_MERYL_DB = expand(
    DIR_RES.joinpath("databases", "singles", "{sample}.k{size_k}.{hpc}.meryl.tar.gz"),
    sample=SAMPLES,
    size_k=MERYL_KMER_VALUES,
    hpc=MERYL_COMPRESS_WILDCARD_VALUES
),
if DISCARD_KMER_DATABASES:
    _SINGLES_RESULT_MERYL_DB = []


rule meryl_run_singles:
    input:
        dbs = _SINGLES_RESULT_MERYL_DB,
        stats = expand(
            DIR_RES.joinpath("statistics", "singles", "{sample}.k{size_k}.{hpc}.meryl-stats.tsv"),
            sample=SAMPLES,
            size_k=MERYL_KMER_VALUES,
            hpc=MERYL_COMPRESS_WILDCARD_VALUES,
        ),
