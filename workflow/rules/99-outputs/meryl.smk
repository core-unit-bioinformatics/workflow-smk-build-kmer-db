
MERYL_WORKFLOW_OUTPUT = []

if RUN_MERYL:

    if PROCESS_SINGLETONS:
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_singles.input.meryl_single_db_stats
        )

    if CHILDREN and PROCESS_TRIOS:
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_trios.input.meryl_trio_db_stats
        )
