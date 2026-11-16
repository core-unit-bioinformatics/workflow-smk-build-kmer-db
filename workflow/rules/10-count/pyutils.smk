
def get_yak_sample_input_files(wildcards):
    """Collect all input files for a sample,
    irrespective of the sample being a single-
    or multi-input sample. Files are returned
    in a stable (sorted) order.
    """
    sample_files = sorted(
        str(input_file) for (sample, part), input_file
        in MAP_SAMPLE_TO_INPUT_FILE.items()
        if sample == wildcards.sample
    )
    if not sample_files:
        raise ValueError(f"No input files for sample: {wildcards.sample}")
    return sample_files


def make_yak_count_input_spec(input_files, bloom_bits):
    """Build the input part of the 'yak count'
    command line. A single (regular) file can be
    re-read by yak for the second (Bloom filter)
    pass; several files are streamed via process
    substitution, which requires repeating the
    stream if the Bloom filter is active.
    'gzip -dcf' handles compressed and plain
    files alike, avoiding broken streams when
    plain and compressed inputs are mixed.
    """
    if len(input_files) == 1:
        return f"{input_files[0]}"
    file_list = " ".join(input_files)
    stream = f"<(gzip -dcf {file_list})"
    if bloom_bits > 0:
        return f"{stream} {stream}"
    return stream
