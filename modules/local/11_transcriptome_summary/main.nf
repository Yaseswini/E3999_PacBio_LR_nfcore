process TRANSCRIPTOME_SUMMARY {

    tag "$sample_id"

    input:
    tuple val(sample_id), path(sqanti_classification), path(ensg_to_gene), path(enst_to_isoname),  path(output_dir)

    output:
    path("*"), emit: summary

    publishDir "${output_dir}", mode: 'copy'

    script:
    """
    echo "Running transcriptome summary for ${sample_id}..."

    python3 ${moduleDir}/04_transcriptome_summary_gene_table_only.py \\
        --sq_out ${sqanti_classification} \\
        --ensg_to_gene ${ensg_to_gene} \\
        --enst_to_isoname ${enst_to_isoname} \\
        --odir ${output_dir}
    """
}
