
rule meryl_trio_disjoin_parent_kmer_dbs:
    input:
        mat_db = lambda wildcards: expand(
            DIR_PROC.joinpath("10-count", "meryl", "10-build", "{parent}.k{{size_k}}.{{hpc}}.meryl"),
            parent=MAP_TRIOS[wildcards.child]["mother"]
        ),
        pat_db = lambda wildcards: expand(
            DIR_PROC.joinpath("10-count", "meryl", "10-build", "{parent}.k{{size_k}}.{{hpc}}.meryl"),
            parent=MAP_TRIOS[wildcards.child]["father"]
        ),
    output:
        mat_only = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "10-disjoin", "{child}.k{size_k}.{hpc}.maternal-only.meryl"))),
        pat_only = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "10-disjoin", "{child}.k{size_k}.{hpc}.paternal-only.meryl"))),
        shared = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "10-disjoin", "{child}.k{size_k}.{hpc}.parental-shared.meryl"))),
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    params:
        tmp_shared = lambda wildcards, output: output.shared.replace(".meryl", ".tmp-db.meryl")
    shell:
        "meryl difference {input.mat_db} {input.pat_db} output {output.mat_only} &> {log}"
            " && "
        "meryl difference {input.pat_db} {input.mat_db} output {output.pat_only} &> {log}"
            " && "
        "meryl union-sum {input.mat_db} {input.pat_db} output {params.tmp_shared} &> {log}"
            " && "
        "meryl difference {params.tmp_shared} {output.mat_only} {output.pat_only} output {output.shared} &> {log}"
            " && "
        "rm -rf {params.tmp_shared}"


rule meryl_trio_create_inherited_kmer_dbs:
    """Following the Merqury reference implementation,
    the 'intersect' operation here sets the k-mer count
    to the first database (= the child)
    """
    input:
        child_db = DIR_PROC.joinpath("10-count", "meryl", "10-build", "{child}.k{size_k}.{hpc}.meryl"),
        mat_only = rules.meryl_trio_disjoin_parent_kmer_dbs.output.mat_only,
        pat_only = rules.meryl_trio_disjoin_parent_kmer_dbs.output.pat_only,
        shared = rules.meryl_trio_disjoin_parent_kmer_dbs.output.shared,
    output:
        mat_inherit = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "20-inherit", "{child}.k{size_k}.{hpc}.maternal-inherit.meryl"))),
        pat_inherit = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "20-inherit", "{child}.k{size_k}.{hpc}.paternal-inherit.meryl"))),
        shared_inherit = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "20-inherit", "{child}.k{size_k}.{hpc}.parental-inherit.meryl"))),
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl intersect {input.child_db} {input.mat_only} output {output.mat_inherit} &> {log}"
            " && "
        "meryl intersect {input.child_db} {input.pat_only} output {output.pat_inherit} &> {log}"
            " && "
        "meryl intersect {input.child_db} {input.shared} output {output.shared_inherit} &> {log}"


rule meryl_trio_create_hapmer_dbs:
    input:
        stats = DIR_RES.joinpath(
            "statistics", "trios", "{child}.k{size_k}.{hpc}.{kmer_set}-inherit.meryl-stats.tsv"
        ),
        inherit = DIR_PROC.joinpath("20-trios", "meryl", "20-inherit", "{child}.k{size_k}.{hpc}.{kmer_set}-inherit.meryl"),
    output:
        hapmer = temp(directory(DIR_PROC.joinpath("20-trios", "meryl", "30-filter", "{child}.k{size_k}.{hpc}.{kmer_set}-hapmer.meryl"))),
    wildcard_constraints:
        kmer_set="(maternal|paternal|parental)"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    params:
        kmer_t = lambda wildcards, input: load_reliable_kmer_threshold(input.stats)
    shell:
        "meryl greater-than {params.kmer_t} {input.inherit} output {output.hapmer}"


_TRIOS_RESULT_MERYL_DB = expand(
    DIR_RES.joinpath("databases", "trios", "{sample}.k{size_k}.{hpc}.{kmer_set}.meryl.tar.gz"),
    sample=CHILDREN,
    kmer_set=["maternal-hapmer", "paternal-hapmer", "parental-hapmer"],
    size_k=config["meryl_kmer_size"],
    hpc=MERYL_COMPRESS_WILDCARD_VALUES,
)
if DISCARD_KMER_DATABASES:
    _TRIOS_RESULT_MERYL_DB = []


rule meryl_run_trios:
    input:
        dbs = _TRIOS_RESULT_MERYL_DB,
        stats = expand(
            DIR_RES.joinpath("statistics", "trios", "{sample}.k{size_k}.{hpc}.{kmer_set}.meryl-stats.tsv"),
            sample=CHILDREN,
            kmer_set=["maternal-inherit", "paternal-inherit", "parental-inherit"],
            size_k=config["meryl_kmer_size"],
            hpc=MERYL_COMPRESS_WILDCARD_VALUES,
        ),
        plots = expand(
            DIR_RES.joinpath("plots", "trios", "{sample}.k{size_k}.{hpc}.{kmer_set}.meryl-thresholds.pdf"),
            sample=CHILDREN,
            kmer_set=["maternal-inherit", "paternal-inherit", "parental-inherit"],
            size_k=config["meryl_kmer_size"],
            hpc=MERYL_COMPRESS_WILDCARD_VALUES,

        )
