
MERYL_WORKFLOW_OUTPUT = []

if RUN_MERYL:

    if PROCESS_SINGLETONS:
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_singles.input.dbs
        )
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_singles.input.stats
        )

    if CHILDREN and PROCESS_TRIOS:
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_trios.input.dbs
        )
        MERYL_WORKFLOW_OUTPUT.extend(
            rules.meryl_run_trios.input.stats
        )

    # PAIRWISE mode needs refactoring
    # if COMPARE_PAIRWISE_BY_FILE:
    #     MERYL_WORKFLOW_OUTPUT.extend(
    #         rules.meryl_run_pairwise_by_file.input.dbs
    #     )
    #     MERYL_WORKFLOW_OUTPUT.extend(
    #         rules.meryl_run_pairwise_by_file.input.stats
    #     )

    # if COMPARE_PAIRWISE_BY_SAMPLE:
    #     MERYL_WORKFLOW_OUTPUT.extend(
    #         rules.meryl_run_pairwise_by_sample.input.dbs
    #     )
    #     MERYL_WORKFLOW_OUTPUT.extend(
    #         rules.meryl_run_pairwise_by_sample.input.stats
    #     )
