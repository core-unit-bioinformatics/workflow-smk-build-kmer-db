
rule meryl_dump_part_file_dbs:
    input:
        seq = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILE[(wildcards.sample, wildcards.part)]
    output:
        db = directory(
            DIR_PROC.joinpath("10-singles", "meryl", "parts", "{sample}.{part}.k{size_k}.{hpc}.meryl")
        )
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
        "../../envs/meryl.yaml"
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
        db = directory(
            DIR_PROC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl")
        )
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
        "../../envs/meryl.yaml"
    shell:
        "meryl union {input.dbs} output {output.db} &> {log}"


rule meryl_dump_kmers:
    input:
        seq = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILE[(wildcards.sample, None)]
    output:
        db = directory(
            DIR_PROC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl")
        )
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
        "../../envs/meryl.yaml"
    params:
        compress = lambda wildcards: "compress" if wildcards.hpc == "is-hpc" else "",
        acc_in=lambda wildcards, input: register_input(input.seq, allow_non_existing=True),
    shell:
        "meryl count k={wildcards.size_k} memory={resources.mem_gb} threads={threads} "
            "{params.compress} {input.seq} output {output.db} &> {log}"


rule meryl_run_singles:
    input:
        meryl_single_dbs = expand(
            DIR_PROC.joinpath("10-singles", "meryl", "{sample}.k{size_k}.{hpc}.meryl"),
            sample=SAMPLES,
            size_k=config["meryl_kmer_size"],
            hpc=MERYL_COMPRESS_WILDCARD_VALUES
        )
