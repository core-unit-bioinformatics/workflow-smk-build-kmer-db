
rule compress_meryl_database:
    """
    Keeping the database in compressed
    form saves about ~50% disk space.
    """
    input:
        db = lambda wildcards: infer_meryl_data_path(
            wildcards.setting, wildcards.db_name, "db", True
        )
    output:
        tar = DIR_RES.joinpath(
                "databases", "{setting}", "{db_name}.meryl.tar.gz"
            )
    benchmark:
        DIR_RSRC.joinpath("databases", "{setting}", "{db_name}.meryl-compress.rsrc")
    wildcard_constraints:
        setting="(singletons|trios|pairwise)"
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: 1 * max(0, attempt - 1),
    params:
        source_dir=lambda wildcards, input: pathlib.Path(input.db).parent,
        source_name=lambda wildcards, input: pathlib.Path(input.db).name,
        acc_out=lambda wildcards, output: register_result(output.tar),
    shell:
        "tar czf {output} -C {params.source_dir} {params.source_name}"
