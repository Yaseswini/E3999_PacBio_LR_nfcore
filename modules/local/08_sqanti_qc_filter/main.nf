process SQANTI_QC_AND_FILTER {
    container 'containers/E3999_PacBio_LR_YN_BuildContainer_sqanti_cupcake.sif'

    tag "SQANTI QC and filter for ${output_prefix}"

    input:
    tuple val(collapsed_gff),
          val(reference_gtf),
          val(reference_genome_fasta),
          val(parsed_reference_dir),
          val(collapsed_abundance),
          val(output_prefix),
          val(output_dir)

    output:
    path("*"), emit: sqanti_outputs
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    echo "Activating conda env and running SQANTI3 for ${output_prefix}"
    
    . /opt/conda/etc/profile.d/conda.sh
    conda activate SQANTI3.env

    ## 1. SQANTI3 QC
    /opt/conda/envs/SQANTI3.env/bin/python /opt/SQANTI3/sqanti3_qc.py \\
        ${collapsed_gff} \\
        ${reference_gtf} \\
        ${reference_genome_fasta} \\
        --skipORF \\
        -o ${output_prefix} \\
        -d . \\
        --fl_count ${collapsed_abundance}

     ## 2. Filtering SQANTI outputs
    /opt/conda/envs/SQANTI3.env/bin/python ${moduleDir}/03_filter_sqanti.py \
        --sqanti_classification ${output_prefix}_classification.txt \
        --sqanti_corrected_fasta ${output_prefix}_corrected.fasta \
        --sqanti_corrected_gtf ${output_prefix}_corrected.gtf \
        --protein_coding_genes ${parsed_reference_dir}/protein_coding_genes.txt \
        --ensg_gene ${parsed_reference_dir}/ensg_gene.tsv \
        --filter_protein_coding yes \
        --filter_intra_polyA yes \
        --filter_template_switching yes \
        --percent_A_downstream_threshold 95 \
        --structural_categories_level strict \
        --minimum_illumina_coverage 3

    ## 3. Collapsing filtered isoforms
    /opt/conda/envs/SQANTI3.env/bin/python ${moduleDir}/03_collapse_isoforms.py \
        --name ${output_prefix} \
        --sqanti_gtf filtered_${output_prefix}_corrected.gtf \
        --sqanti_fasta filtered_${output_prefix}_corrected.fasta

    ## 4. Collapsing classification
    /opt/conda/envs/SQANTI3.env/bin/python ${moduleDir}/03_collapse_classification.py \
        --name ${output_prefix} \
        --collapsed_fasta ${output_prefix}_corrected.5degfilter.fasta \
        --classification filtered_${output_prefix}_classification.tsv 

    """
}
