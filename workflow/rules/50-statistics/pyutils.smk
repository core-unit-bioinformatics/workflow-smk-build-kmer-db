import pathlib


def infer_any_meryl_database_path(setting, db_name):
    """
    """
    db_source_path = None
    if setting == "pairwise":
        pairwise_indicators = {
            "1-or-2": pathlib.Path("30-pairwise", "meryl", "union-sum"),
            "1-and-2": pathlib.Path("30-pairwise", "meryl", "intersect-min"),
            "1-not-2": pathlib.Path("30-pairwise", "meryl", "difference"),
            "2-not-1": pathlib.Path("30-pairwise", "meryl", "difference"),
        }
        for pwi, pw_path in pairwise_indicators.items():
            if pwi not in db_name:
                continue
            db_source_path = DIR_PROC.joinpath(
                pw_path, f"{db_name}.meryl"
            )
            break
    elif setting == "trios":
        db_source_path = DIR_PROC.joinpath(
            "20-trios", "meryl", f"{db_name}.meryl"
        )
    else:
        assert setting == "singles"
        # TODO this should be refactored into one
        # assume single
        # reuse utility function from module
        # 30-pairwise::pyutils.smk
        db_source_path = infer_meryl_single_database_path(db_name)

    if db_source_path is None:
        err_msg = (
            "ERROR: cannot infer meryl database source path for name:\n"
            f"{db_name} /// {setting}"
        )
        raise ValueError(err_msg)

    return db_source_path
