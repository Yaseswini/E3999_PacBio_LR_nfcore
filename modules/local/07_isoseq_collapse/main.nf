process ISOSEQ_COLLAPSE {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(condition), path(aligned_bam), path(bam_for_analysis), path(pbi_file), path(output_dir)

    output:
    tuple val(sample_id), val(condition),
          path("${sample_id}_isoseq_collapsed.abundance.txt"),
          path("${sample_id}_isoseq_collapsed.fasta"),
          path("${sample_id}_isoseq_collapsed.flnc_count.txt"),
          path("${sample_id}_isoseq_collapsed.gff"),
          path("${sample_id}_isoseq_collapsed.group.txt"),
          path("${sample_id}_isoseq_collapsed.read_stat.txt"),
          path("${sample_id}_isoseq_collapsed.report.json"),
          emit: collapsed

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    echo "Running IsoSeq collapse for ${sample_id} (${condition})..."
    
    isoseq collapse --do-not-collapse-extra-5exons \\
      ${aligned_bam} \\
      ${bam_for_analysis} \\
      ${sample_id}_isoseq_collapsed.gff

    """
}
