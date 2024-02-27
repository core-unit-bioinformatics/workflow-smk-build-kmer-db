

def infer_meryl_data_path(setting, db_name, which_output=None):

    assert which_output in [
        "db", "database", "stats-dump", "stats-table", "stats-hist"
        ] or which_output is None
    assert setting in ["singletons", "trios", "pairwise"]

    if setting == "singletons":
        proc_base = DIR_PROC.joinpath("10-count", "meryl", "20-filter")
    elif setting == "trios":
        proc_base = DIR_PROC.joinpath("20-trios", "meryl", "30-filter")
    elif setting == "pairwise":
        raise NotImplementedError()
    else:
        raise RuntimeError(f"Should not happen: {setting}")

    select_path = None
    if which_output in ["db", "database"] or which_output is None:
        select_path = get_meryl_database(proc_base, db_name)
    if which_output in ["stats-dump"]:
        select_path = get_meryl_statistics_dump(proc_base, db_name)
    if which_output in ["stats-table"]:
        select_path = DIR_RES.joinpath(
            "statistics", setting, db_name + ".meryl-stats.tsv"
        )
    if select_path is None:
        raise ValueError(f"Cannot find meryl data: {setting} / {db_name} / {which_output}")

    return select_path


def get_meryl_database(top_level, db_name):

    all_dbs = list(top_level.glob("*.meryl"))
    select_db = [db for db in all_dbs if db_name in db.name]
    if len(select_db) > 1:
        raise ValueError(f"Ambiguous meryl DB: {top_level} / {db_name}")
    elif len(select_db) < 1:
        # may not be created yet
        select_db = top_level.joinpath(db_name + ".meryl")
    else:
        select_db = select_db[0]
        assert select_db.is_dir()
    return select_db


def get_meryl_statistics_dump(top_level, db_name):

    db_path = get_meryl_database(top_level, db_name)
    stats_dump = db_path.with_suffix(".meryl-stats.txt")
    return stats_dump
