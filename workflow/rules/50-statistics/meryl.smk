
rule meryl_dump_db_statistics:
    input:
        db = lambda wildcards: infer_any_meryl_database_path(
            wildcards.setting, wildcards.db_name
        )
    output:
        stats_dump = temp(
            DIR_RES.joinpath(
                "statistics", "{setting}", "{db_name}.meryl-stats.txt"
            )
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "{setting}", "{db_name}.meryl-stats.txt")
    wildcard_constraints:
        setting="(singles|trios|pairwise)"
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * max(0, attempt - 1),
    conda:
        DIR_ENVS.joinpath("meryl.yaml")
    shell:
        "meryl statistics {input.db} > {output.stats_dump}"


rule meryl_split_statistics_dump:
    input:
        stats_dump = DIR_RES.joinpath(
            "statistics", "{setting}", "{db_name}.meryl-stats.txt"
        )
    output:
        stats = DIR_RES.joinpath(
            "statistics", "{setting}", "{db_name}.meryl-stats.tsv"
        )
        hist = DIR_RES.joinpath(
            "statistics", "{setting}", "{db_name}.meryl-hist.tsv.gz"
        )
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        time_hrs=lambda wildcards, attempt: 1 * attempt,
    conda:
        DIR_ENVS.joinpath("pyscript.yaml")
    params:
        script=find_script("split_meryl_stats", extension="py"),
        forced_thres=FORCE_KMERFREQ_THRESHOLD,
        threshold_warning=lambda wildcards: (
            " --no-threshold-warning " if wildcards.setting != "singles"
            else " "
        ),
        acc_out=lambda wildcards, output: register_result(output),
    shell:
        "{params.script} --meryl-stats {input.stats_dump} "
            "--forced-kmerfreq-threshold {params.forced_thres} "
            "{params.threshold_warning}"
            "--out-stats {output.stats} --out-hist {output.hist}"
