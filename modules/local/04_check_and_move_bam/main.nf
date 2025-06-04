process CHECK_AND_MOVE_BAM {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(condition), path(bam_file), path(pbi_file)
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    mkdir -p data/${sample_id}
    cp ${bam_file} data/${sample_id}/$bam_file
    cp ${pbi_file} data/${sample_id}/$pbi_file
    """
}