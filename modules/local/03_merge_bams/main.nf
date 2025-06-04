process MERGE_BAMS {
    
    input:
    tuple val(sample_id), val(condition), val(bam_paths)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_merged.bam")

    script:
    """
    samtools merge -f ${sample_id}_merged.bam ${bam_paths.join(' ')}
    """
}
