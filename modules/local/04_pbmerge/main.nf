process PBMERGE_BAMS {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(condition), path(bam_files)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_pbmerged.flnc.bam"), path("${sample_id}_pbmerged.flnc.bam.pbi"), emit: pbmerge_output

    publishDir "data/${sample_id}", mode: 'copy'
    
    script:
    """
    echo "[PBMERGE] Merging files:"
    ls ${bam_files}

    pbmerge -o ${sample_id}_pbmerged.flnc.bam ${bam_files}
    """
}