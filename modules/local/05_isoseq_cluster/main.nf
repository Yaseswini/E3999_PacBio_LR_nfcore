process ISOSEQ_CLUSTER {
    tag "$sample_id"

    input:
    tuple val(sample_id), val(condition), path(bam_file)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_isoseq_clustered.bam"), path("${sample_id}_isoseq_clustered.bam.pbi") , emit: clustered_bam

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    echo "Clustering using IsoSeq for sample ${sample_id} (${condition}) started..."

    # Run IsoSeq clustering
    isoseq cluster2 ${bam_file} ${sample_id}_isoseq_clustered.bam --use-qvs

    """
}

