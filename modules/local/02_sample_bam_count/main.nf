process PRINT_SAMPLE_BAM_COUNT {
    input:
    tuple val(sample_id), val(condition), val(bam_paths)

    output:
    file "sample_bam_count.txt"

    script:
    """
    echo "Sample: ${sample_id}, Condition: ${condition}, BAM files: ${bam_paths.size()}" > sample_bam_count.txt
    """
}
