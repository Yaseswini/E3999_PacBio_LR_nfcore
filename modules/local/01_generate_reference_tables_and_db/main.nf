process GENERATE_REFERENCE_TABLES_AND_DB {
    tag "Generating reference tables and clustering database"

    input:
    tuple val(reference_gtf), val(reference_fasta), val(reference_dir)

    output:
    path("*.tsv"), emit: reference_outputs
    path("*.txt")
    path("*.fasta")

    publishDir "${reference_dir}", mode: 'copy', copyOptions: 'rsync'

    script:
    """
    echo "Reference directory: ${reference_dir}"
    mkdir -p ${reference_dir}
    echo "Reference directory created at ${reference_dir}"

    python3 ${moduleDir}/01_prepare_reference_tables.py \\
        --gtf ${reference_gtf} \\
        --fa ${reference_fasta} \\
        --ensg_gene ensg_gene.tsv \\
        --enst_isoname enst_isoname.tsv \\
        --gene_ensp gene_ensp.tsv \\
        --gene_isoname gene_isoname.tsv \\
        --isoname_lens isoname_lens.tsv \\
        --gene_lens gene_lens.tsv \\
        --protein_coding_genes protein_coding_genes.txt

    python3 ${moduleDir}/02_make_gencode_database.py \\
        --gencode_fasta ${reference_fasta} \\
        --protein_coding_genes protein_coding_genes.txt \\
        --output_fasta gencode_clusters.fasta \\
        --output_cluster gencode_isoname_clusters.tsv
    """

}
