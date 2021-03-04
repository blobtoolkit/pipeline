def reads_by_prefix(config):
    """Return read meta by prefix"""
    reads = {}
    if "reads" not in config:
        return {}

    for strategy in ("paired", "single"):
        if not strategy in config["reads"] or not config["reads"][strategy]:
            continue
        reads.update({entry[0]: entry for entry in config["reads"][strategy]})

    return reads


def minimap_tuning(config, prefix):
    """Set minimap2 mapping parameter."""
    reads = reads_by_prefix(config)
    tunings = {"ILLUMINA": "sr", "OXFORD_NANOPORE": "map-ont", "PACBIO_SMRT": "map-pb"}
    return tunings[reads[prefix][1]]


def read_files(config, prefix):
    """Set minimap2 mapping parameter."""
    reads = reads_by_prefix(config)
    return reads[prefix][3].split(";")


def seqtk_sample_input(config, prefix):
    """Generate seqtk command to subsamplereads if required."""
    meta = reads_by_prefix(config)[prefix]
    filenames = meta[3].split(";")
    ratio = 1
    if "coverage" in config["reads"] and "max" in config["reads"]["coverage"]:
        base_count = meta[2]
        if isinstance(base_count, int):
            ratio = (
                config["assembly"]["span"]
                * config["reads"]["coverage"]["max"]
                / base_count
            )
    if ratio <= 0.8:
        command = " ".join(
            [
                "<(seqtk sample -s 100 %s %.2f)" % (filename, ratio)
                for filename in filenames
            ]
        )
    else:
        command = " ".join(filenames)
    return command
