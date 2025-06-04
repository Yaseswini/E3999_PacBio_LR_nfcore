process SIX_FRAME_TRANSLATION {

    tag "$sample_id"

    input:
    tuple val(sample_id), 
          val(condition), 
          path(filtered_classification_file), 
          path(ensg_gene_file), 
          path(filtered_fasta_file), 
          path(output_dir)

    output:
    path("*"), emit: sixframe

    publishDir "${output_dir}", mode: 'copy'

    script:
    """
    echo "Running six-frame translation for ${sample_id}..."

    python3 ${moduleDir}/04_six_frame_translation.py \\
        --iso_annot ${filtered_classification_file} \\
        --ensg_gene ${ensg_gene_file} \\
        --sample_fasta ${filtered_fasta_file} \\
        --output_fasta ${sample_id}_6frame_database_gene_grouped.fasta \\
        2> ${sample_id}_6frame_database_gene_grouped.error
    """
}

