process PARSE_SAMPLE_METADATA {
    tag "$input_metadata"

    input:
    path input_metadata
    val base_bam_dir

    output:
    stdout
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], row[1], row[2]) }

    script:
    """
    ./modules/local/sample_metadata/parse_metadata.py ${input_metadata} ${base_bam_dir}
    """
}