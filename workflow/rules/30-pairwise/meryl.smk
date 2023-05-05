

rule meryl_pairwise_union:
    input:
        db1 = lambda wildcards: infer_meryl_single_database_path(wildcards.db1),
        db2 = lambda wildcards: infer_meryl_single_database_path(wildcards.db2),
    output:
        db = temp(directory(
            DIR_PROC.joinpath("30-pairwise", "meryl", "union-sum",
                "{db1}_vs_{db2}.1-or-2.meryl")
        ))
    log:
        DIR_LOG.joinpath("30-pairwise", "meryl", "union-sum", "{db1}_vs_{db2}.1-or-2.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("30-pairwise", "meryl", "union-sum", "{db1}_vs_{db2}.1-or-2.meryl.rsrc")
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl union-sum {input.db1} {input.db2} output {output.db} &> {log}"


rule meryl_pairwise_intersect:
    input:
        db1 = lambda wildcards: infer_meryl_single_database_path(wildcards.db1),
        db2 = lambda wildcards: infer_meryl_single_database_path(wildcards.db2),
    output:
        db = temp(directory(
            DIR_PROC.joinpath("30-pairwise", "meryl", "intersect-min",
                "{db1}_vs_{db2}.1-and-2.meryl")
        ))
    log:
        DIR_LOG.joinpath("30-pairwise", "meryl", "intersect-min", "{db1}_vs_{db2}.1-and-2.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("30-pairwise", "meryl", "intersect-min", "{db1}_vs_{db2}.1-and-2.meryl.rsrc")
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl intersect-min {input.db1} {input.db2} output {output.db} &> {log}"


rule meryl_pairwise_difference:
    input:
        db1 = lambda wildcards: infer_meryl_single_database_path(wildcards.db1),
        db2 = lambda wildcards: infer_meryl_single_database_path(wildcards.db2),
    output:
        db12 = temp(directory(
            DIR_PROC.joinpath("30-pairwise", "meryl", "difference",
                "{db1}_vs_{db2}.1-not-2.meryl")
        )),
        db21 = temp(directory(
            DIR_PROC.joinpath("30-pairwise", "meryl", "difference",
                "{db1}_vs_{db2}.2-not-1.meryl")
        ))
    log:
        DIR_LOG.joinpath("30-pairwise", "meryl", "difference", "{db1}_vs_{db2}.N-not-M.meryl.log")
    benchmark:
        DIR_RSRC.joinpath("30-pairwise", "meryl", "difference", "{db1}_vs_{db2}.N-not-M.meryl.rsrc")
    resources:
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl difference {input.db1} {input.db2} output {output.db12} &> {log}"
            " && "
        "meryl difference {input.db2} {input.db1} output {output.db21} &>> {log}"


##########################
## Prepare output pull
##########################


_PAIRWISE_BY_FILE_RESULT_MERYL_DB = expand(
    DIR_RES.joinpath(
        "databases", "pairwise", "{db_pair}.{operation}.meryl.tar.gz"
    ),
    db_pair=pair_all_inputs(
        True,
        MERYL_KMER_VALUES,
        MERYL_COMPRESS_WILDCARD_VALUES
    ),
    operation=["1-or-2", "1-and-2", "1-not-2", "2-not-1"]
)

_PAIRWISE_BY_SAMPLE_RESULT_MERYL_DB = expand(
    DIR_RES.joinpath(
        "databases", "pairwise", "{db_pair}.{operation}.meryl.tar.gz"
    ),
    db_pair=pair_all_inputs(
        False,
        MERYL_KMER_VALUES,
        MERYL_COMPRESS_WILDCARD_VALUES
    ),
    operation=["1-or-2", "1-and-2", "1-not-2", "2-not-1"]
)
if DISCARD_KMER_DATABASES:
    _PAIRWISE_BY_FILE_RESULT_MERYL_DB = []
    _PAIRWISE_BY_SAMPLE_RESULT_MERYL_DB = []


rule meryl_run_pairwise_by_file:
    """
    db_pair(True ...) -> pair by file
    """
    input:
        dbs = _PAIRWISE_BY_FILE_RESULT_MERYL_DB,
        stats = expand(
            DIR_RES.joinpath(
                "statistics", "pairwise", "{db_pair}.{operation}.meryl-stats.tsv"
            ),
            db_pair=pair_all_inputs(
                True,
                MERYL_KMER_VALUES,
                MERYL_COMPRESS_WILDCARD_VALUES
            ),
            operation=["1-or-2", "1-and-2", "1-not-2", "2-not-1"]
        )


rule meryl_run_pairwise_by_sample:
    """
    db_pair(False ...) -> pair by sample
    """
    input:
        dbs = _PAIRWISE_BY_SAMPLE_RESULT_MERYL_DB,
        stats = expand(
            DIR_RES.joinpath(
                "statistics", "pairwise", "{db_pair}.{operation}.meryl-stats.tsv"
            ),
            db_pair=pair_all_inputs(
                False,
                MERYL_KMER_VALUES,
                MERYL_COMPRESS_WILDCARD_VALUES
            ),
            operation=["1-or-2", "1-and-2", "1-not-2", "2-not-1"]
        )
