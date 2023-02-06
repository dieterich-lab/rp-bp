# Controls for webapp

# for read filtering counts
STACK_CTS_ORDER = [
    "length_count",
    "unique_count",
    "genome_count",
    "without_rrna_count",
    "without_adapters_count",
    "raw_data_count",
]

STACK_CTS_NAME = [
    "Usable",
    "Non-periodic",
    "Multimappers",
    "No alignment",
    "Ribosomal",
    "Poor quality",
]

# these are now in the right order
FUNNEL_CTS_NAME = [
    "Raw reads",
    "Adapters/low quality removed",
    "rRNA-free reads",
    "Aligned to genome",
    "Uniquely aligned",
    "Periodic reads (usable)",
]
