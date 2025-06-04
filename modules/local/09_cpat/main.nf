process CPAT_PREDICTION {

    container 'containers/E3999_PacBio_LR_YN_Container_CPAT.sif'
    tag "$sample_id"

    input:
    tuple val(sample_id), val(condition), path(filtered_fasta_file), path(output_dir)

    output:
    path("*"), emit: cpat

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "*"

    script:
    """
    echo "Running CPAT for ${sample_id}..."

    cpat -x ${moduleDir}/Human_Hexamer.tsv -d ${moduleDir}/Human_logitModel.RData -g ${filtered_fasta_file} --min-orf=50 \\
        --top-orf=50 \\
        -o ${sample_id}_cpat \\
        2> ${sample_id}_cpat.error
    """
}
