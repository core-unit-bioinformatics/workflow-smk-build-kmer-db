import itertools


def pair_all_inputs(by_file, kmer_sizes, hpc_setting):
    """
    hpc_setting: count in homopolymer-compressed space
    """

    _this_fun = "30-pairwise::pyutils.smk::pair_all_inputs"

    if by_file:
        # pair all input files / created databases
        dbs_to_pair = list(MAP_SAMPLE_TO_INPUT_FILE.keys())
    else:
        dbs_to_pair = [(sample, None) for sample in SAMPLES]
    assert dbs_to_pair

    dbs_to_pair = sorted(dbs_to_pair)
    db_pair_specs = []
    for (s1, p1), (s2, p2) in itertools.combinations(dbs_to_pair, 2):
        for size_k in kmer_sizes:
            for hpc in hpc_setting:
                if p1 is None and p2 is None:
                    db1 = f"{s1}.k{size_k}.{hpc}"
                    db2 = f"{s2}.k{size_k}.{hpc}"
                elif p1 is None:
                    db1 = f"{s1}.k{size_k}.{hpc}"
                    db2 = f"{s2}.{p2}.k{size_k}.{hpc}"
                elif p2 is None:
                    db1 = f"{s1}.{p1}.k{size_k}.{hpc}"
                    db2 = f"{s2}.k{size_k}.{hpc}"
                else:
                    db1 = f"{s1}.{p1}.k{size_k}.{hpc}"
                    db2 = f"{s2}.{p2}.k{size_k}.{hpc}"
                pair_spec = f"{db1}_vs_{db2}"
                db_pair_specs.append(pair_spec)
    assert db_pair_specs

    if DOWNSAMPLE_COMPARISONS > 0.:
        import random
        random.seed(None)
        select_subset = int(len(db_pair_specs) * DOWNSAMPLE_COMPARISONS)
        if 0 <= select_subset < 1:
            if VERBOSE:
                warn_msg = (
                    f"WARNING: {_this_fun}\n"
                    f"Downsampling comparisons (total: {len(db_pair_specs)}) "
                    f"to fraction {DOWNSAMPLE_COMPARISONS} resulted in 0-length "
                    "subset. Increasing subset size to 1\n"
                )
                logerr(warn_msg)
            select_subset = 1
        elif select_subset >= 1:
            if VERBOSE:
                info_msg = (
                    f"INFO: {_this_fun}\n"
                    f"Downsampling comparisons from a total of {len(db_pair_specs)} "
                    f"to a subset of size {select_subset} to be processed in this run.\n"
                )
                logout(info_msg)
        else:
            raise ValueError(f"ERROR {_this_fun}: illegal subset size: {select_subset}")
        db_pair_specs = random.sample(db_pair_specs, select_subset)

    return db_pair_specs


def infer_meryl_single_database_path(database_name):
    """ Pairwise computations may require
    working with individual input files
    (e.g., SMRT or Flow cells), hence the required
    database might be stored under 'parts' before
    the per-sample merging is performed.
    """

    assert not database_name.endswith("meryl")
    name_components = database_name.split(".")
    assert name_components[-1] != "meryl"
    if len(name_components) == 4:
        # DB name has part (hash) identifier
        db_path = DIR_PROC.joinpath(
            "10-singles", "meryl", "parts", f"{database_name}.meryl"
        )
    elif len(name_components) == 3:
        # single input for sample
        db_path = DIR_PROC.joinpath(
            "10-singles", "meryl", f"{database_name}.meryl"
        )
    else:
        errmsg = (
            "30-pairwise::pyutils.smk::infer_database_path\n"
            f"Cannot infer database path: {database_name}"
        )
        logerr(errmsg)
        raise ValueError(f"Cannot process database name: {database_name}")

    return db_path
