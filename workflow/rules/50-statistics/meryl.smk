
rule meryl_dump_db_statistics:
    input:
        db = DIR_PROC.joinpath("{scenario}", "meryl", "{sample_db}.meryl")
    output:
        stats_dump = DIR_PROC.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-stats.txt"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-stats.rsrc"
        )
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * max(0, attempt - 1),
    conda:
        "../../envs/meryl.yaml"
    shell:
        "meryl statistics {input.db} > {output.stats_dump}"


rule meryl_split_statistics_dump:
    input:
        stats_dump = DIR_PROC.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-stats.txt"
        )
    output:
        stats = DIR_RES.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-stats.tsv"
        ),
        hist = DIR_RES.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-hist.tsv.gz"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "50-statistics", "{scenario}", "meryl", "{sample_db}.meryl-stats.rsrc"
        )
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        time_hrs=lambda wildcards, attempt: 1 * attempt,
    conda:
        "../../envs/pyscript.yaml"
    params:
        script=find_script("split_meryl_stats", extension="py"),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "{params.script} --meryl-stats {input.stats_dump} "
            "--out-stats {output.stats} --out-hist {output.hist}"
