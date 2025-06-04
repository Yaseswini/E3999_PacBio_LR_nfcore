process PBMM2_ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(condition), path(clustered_bam), val(reference_genome_fasta)

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}_minimap2_aligned.bam"),
          path("${sample_id}_minimap2_aligned.bam.bai"),
          emit: aligned_bam

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    echo "Aligning clustered reads for sample ${sample_id} (${condition})..."

    pbmm2 align \\
        ${reference_genome_fasta} \\
        ${clustered_bam} \\
        ${sample_id}_minimap2_aligned.bam \\
        --preset ISOSEQ \\
        --sort \\
        -j ${task.cpus} \\
        --log-level INFO

    """
}



