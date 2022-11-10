"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/settings.smk"
include: "00-prepare/sample_table.smk"
include: "10-singles/meryl.smk"
include: "20-trios/meryl.smk"

include: "99-outputs/meryl.smk"

include: "99-outputs/aggregate.smk"
