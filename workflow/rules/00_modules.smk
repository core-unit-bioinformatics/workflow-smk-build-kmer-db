"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"

include: "05-utils/pyutils.smk"
include: "05-utils/meryl.smk"

include: "10-count/meryl.smk"

include: "20-trios/pyutils.smk"
include: "20-trios/meryl.smk"

# PAIRWISE mode needs refactoring
#include: "30-pairwise/pyutils.smk"
#include: "30-pairwise/meryl.smk"

include: "60-finalize/meryl.smk"

include: "99-outputs/meryl.smk"
