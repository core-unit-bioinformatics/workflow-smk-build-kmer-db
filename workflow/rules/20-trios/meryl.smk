
rule meryl_create_inherited_kmer_dbs:
    input:
        child = DIR_PROC.joinpath("10-singles", "meryl", "{child}.k{size_k}.{hpc}.meryl"),
        child_stats = DIR_RES.joinpath("statistics", "singles", "{child}.k{size_k}.{hpc}.meryl-stats.tsv"),
        parent = lambda wildcards: expand(
            DIR_PROC.joinpath("10-singles", "meryl", "{parent}.k{{size_k}}.{{hpc}}.meryl"),
            parent=MAP_TRIOS[wildcards.child]["mother"] \
                if wildcards.hap == "maternal" \
                else MAP_TRIOS[wildcards.child]["father"]
        ),
        parent_stats = lambda wildcards: expand(
            DIR_RES.joinpath("statistics", "singles", "{parent}.k{{size_k}}.{{hpc}}.meryl-stats.tsv"),
            parent=MAP_TRIOS[wildcards.child]["mother"] \
                if wildcards.hap == "maternal" \
                else MAP_TRIOS[wildcards.child]["father"]
        ),
    output:
        db = temp(directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.meryl"
            )
        ))
    log:
        DIR_LOG.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.meryl.log"
            )
    benchmark:
        DIR_RSRC.joinpath(
                "20-trios", "meryl",
                "{child}.{hap}.k{size_k}.{hpc}.meryl.rsrc"
            )
    wildcard_constraints:
        hap = "(maternal|paternal)"
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    params:
        min_kfreq_child = lambda wildcards, input: load_reliable_kmer_threshold(input.child_stats),
        min_kfreq_parent = lambda wildcards, input: load_reliable_kmer_threshold(input.parent_stats),
    shell:
        "meryl intersect-min "
            "[ greater-than {params.min_kfreq_child} {input.child} ] "
            "[ greater-than {params.min_kfreq_parent} {input.parent} ] "
            " output {output.db} "
            "&> {log}"


rule meryl_create_specific_and_shared_kmer_dbs:
    input:
        mat = DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.maternal.k{size_k}.{hpc}.meryl"
            ),
        pat = DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.paternal.k{size_k}.{hpc}.meryl"
            )
    output:
        mat_only = temp(directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.maternal-only.k{size_k}.{hpc}.meryl"
            )
        )),
        pat_only = temp(directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.paternal-only.k{size_k}.{hpc}.meryl"
            )
        )),
        shared = temp(directory(
            DIR_PROC.joinpath(
                "20-trios", "meryl",
                "{child}.parental-shared.k{size_k}.{hpc}.meryl"
            )
        ))
    log:
        DIR_LOG.joinpath(
            "20-trios", "meryl",
            "{child}.parents.k{size_k}.{hpc}.meryl.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-trios", "meryl",
            "{child}.parents.k{size_k}.{hpc}.meryl.rsrc"
        )
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
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


_TRIOS_RESULT_MERYL_DB = expand(
    DIR_RES.joinpath("databases", "trios", "{sample}.{kmer_set}.k{size_k}.{hpc}.meryl.tar.gz"),
    sample=CHILDREN,
    kmer_set=["maternal-only", "paternal-only", "parental-shared"],
    size_k=config["meryl_kmer_size"],
    hpc=MERYL_COMPRESS_WILDCARD_VALUES,
)
if DISCARD_KMER_DATABASES:
    _TRIOS_RESULT_MERYL_DB = []


rule meryl_run_trios:
    input:
        dbs = [],
        stats = expand(
            DIR_RES.joinpath("statistics", "trios", "{sample}.{kmer_set}.k{size_k}.{hpc}.meryl-stats.tsv"),
            sample=CHILDREN,
            kmer_set=["maternal-only", "paternal-only", "parental-shared"],
            size_k=config["meryl_kmer_size"],
            hpc=MERYL_COMPRESS_WILDCARD_VALUES,
        ),
