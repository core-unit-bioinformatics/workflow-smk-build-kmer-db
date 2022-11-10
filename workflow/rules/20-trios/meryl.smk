
rule meryl_create_inherited_kmer_dbs:
    input:
        child = DIR_PROC.joinpath("10-singles", "meryl", "{child}.k{size_k}.{hpc}.meryl"),
        parent = lambda wildcards: expand(
            DIR_PROC.joinpath("10-singles", "meryl", "{parent}.k{{size_k}}.{{hpc}}.meryl"),
            parent=MAP_TRIOS[wildcards.child]["mother"] \
                if wildcards.hap == "maternal" \
                else MAP_TRIOS[wildcards.child]["father"]
        )
    output:
        db = directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.geq{min_freq}.meryl"
            )
        )
    log:
        DIR_LOG.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.geq{min_freq}.meryl.log"
            )
    benchmark:
        DIR_RSRC.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.geq{min_freq}.meryl.rsrc"
            )
    wildcard_constraints:
        hap = "(maternal|paternal)"
    conda:
        "../../envs/meryl.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    params:
        min_kfreq = lambda wildcards: int(wildcards.min_freq) - 1  # meryl only knows "greater-than"
    shell:
        "meryl greater-than {params.min_kfreq} "
            "[ intersect-min {input.child} {input.parent} ] output {output.db} "
            "&> {log}"


rule meryl_create_specific_and_shared_kmer_dbs:
    input:
        mat = DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.maternal.k{size_k}.{hpc}.geq{min_freq}.meryl"
            ),
        pat = DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.paternal.k{size_k}.{hpc}.geq{min_freq}.meryl"
            )
    output:
        mat_only = directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.maternal-only.k{size_k}.{hpc}.geq{min_freq}.meryl"
            )
        ),
        pat_only = directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.paternal-only.k{size_k}.{hpc}.geq{min_freq}.meryl"
            )
        ),
        shared = directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.parental-shared.k{size_k}.{hpc}.geq{min_freq}.meryl"
            )
        )
    log:
        DIR_LOG.joinpath(
            "20-trios", "meryl",
            "{child}.parents.k{size_k}.{hpc}.geq{min_freq}.meryl.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-trios", "meryl",
            "{child}.parents.k{size_k}.{hpc}.geq{min_freq}.meryl.rsrc"
        )
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    conda:
        "../../envs/meryl.yaml"
    params:
        tmp_shared = lambda wildcards, output: output.shared.replace(".meryl", ".tmp-db.meryl")
    shell:
        "meryl difference {input.mat} {input.pat} output {output.mat_only} &> {log}"
            " && "
        "meryl difference {input.pat} {input.mat} output {output.pat_only} &> {log}"
            " && "
        "meryl union-sum {input.mat} {input.pat} output {params.tmp_shared} &> {log}"
            " && "
        "meryl difference {params.tmp_shared} {output.mat_only} {output.pat_only} output {output.shared} &> {log}"
            " && "
        "rm -rf {params.tmp_shared}"


rule meryl_run_trios:
    input:
        meryl_trio_dbs = expand(
            DIR_PROC.joinpath("20-trios", "meryl", "{sample}.parental-shared.k{size_k}.{hpc}.geq{min_kfreq}.meryl"),
            sample=CHILDREN,
            size_k=config["meryl_kmer_size"],
            hpc=MERYL_COMPRESS_WILDCARD_VALUES,
            min_kfreq=config["meryl_min_kfreq"]
        )
