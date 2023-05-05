import itertools


def pair_all_inputs(by_file, kmer_sizes, hpc_setting):
    """
    hpc_setting: count in homopolymer-compressed space
    """

    if by_file:
        # pair all input files / created databases
        dbs_to_pair = list(MAP_SAMPLE_TO_INPUT_FILE.keys())
    else:
        dbs_to_pair = [(sample, None) for sample in SAMPLES]
    assert dbs_to_pair

    dbs_to_pair = sorted(dbs_to_pair)
    db_pair_specs = []
    for (s1, p1), (s2, p2) in itertools.combinations(dbs_to_pair):
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
    if len(name_components) == 5:
        # DB name has part (hash) identifier
        db_path = DIR_PROC.joinpath(
            "10-singles", "meryl", "parts", f"{database_name}.meryl"
        )
    elif len(name_components) == 4:
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
