
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


rule plot_jaccard_similarity:
    """TODO
    This will break if the workflow is executed with more
    than one k-mer size (or hpc setting) because different
    parameter combinations cannot be plotted in a single
    clustermap. One workaround would be to adapt the
    plotting script to internally sort all files by
    their parameterization and then produce one plot
    per parameter setting.
    """
    input:
        stats = expand(
            DIR_RES.joinpath(
                "statistics", "pairwise", "{db_pair}.{db_op}.meryl-stats.tsv"
            ),
            db_pair=pair_all_inputs(
                False,
                MERYL_KMER_VALUES,
                MERYL_COMPRESS_WILDCARD_VALUES
            ),
            db_op=["1-or-2", "1-and-2"]
            allow_missing=True
        )
    output:
        pdf = expand(
            DIR_RES.joinpath(
                "plots", "pairwise", "SAMPLES.k{size_k}.{hpc}.jaccard-similarity.pdf"
            ),
            size_k=MERYL_KMER_VALUES,
            hpc=MERYL_COMPRESS_WILDCARD_VALUES
        )
    conda:
        DIR_ENVS.joinpath("pyplot.yaml")
    params:
        script=find_script("plot_kmer_similarity")
    shell:
        "{params.script} --input {input} --output {output}"
