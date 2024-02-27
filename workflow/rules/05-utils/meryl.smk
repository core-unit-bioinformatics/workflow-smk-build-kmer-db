"""
This module sits at the beginning of the pipeline
to provide base functionality for all subsequent
(meryl) operations.
"""


rule meryl_dump_db_statistics:
    input:
        db = DIR_PROC.joinpath("{filepath}", "{db_name}.meryl")
    output:
        stats_dump = temp(
            DIR_PROC.joinpath("{filepath}", "{db_name}.meryl-stats.txt")
        )
    benchmark:
        DIR_RSRC.joinpath("{filepath}", "{db_name}.meryl-stats.rsrc")
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
        stats_dump = lambda wildcards: infer_meryl_data_path(
            wildcards.setting, wildcards.db_name, "stats-dump"
        )
    output:
        stats = DIR_RES.joinpath(
            "statistics", "{setting}", "{db_name}.meryl-stats.tsv"
        ),
        hist = DIR_RES.joinpath(
            "statistics", "{setting}", "{db_name}.meryl-hist.tsv.gz"
        )
    wildcard_constraints:
        setting="(singletons|trios|pairwise)"
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
            " --no-threshold-warning " if wildcards.setting != "singletons"
            else " "
        ),
        acc_out=lambda wildcards, output: register_result(output),
    shell:
        "{params.script} --meryl-stats {input.stats_dump} "
            "--forced-kmerfreq-threshold {params.forced_thres} "
            "{params.threshold_warning} "
            "--out-stats {output.stats} --out-hist {output.hist}"


localrules: meryl_plot_thresholds
rule meryl_plot_thresholds:
    input:
        stats = rules.meryl_split_statistics_dump.output.stats,
        hist = rules.meryl_split_statistics_dump.output.hist
    output:
        pdf = DIR_RES.joinpath(
            "plots", "{setting}", "{db_name}.meryl-thresholds.pdf"
        )
    conda:
        DIR_ENVS.joinpath("pyscript.yaml")
    params:
        acc_out=lambda wildcards, output: register_result(output),
    shell:
        "exit 1"
