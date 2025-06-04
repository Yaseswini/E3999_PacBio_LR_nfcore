nextflow.enable.dsl=2

include { GENERATE_REFERENCE_TABLES_AND_DB } from './modules/local/01_generate_reference_tables_and_db/main.nf';
include { PRINT_SAMPLE_BAM_COUNT } from './modules/local/02_sample_bam_count/main.nf';
include { MERGE_BAMS }             from './modules/local/03_merge_bams/main.nf';
include { PBMERGE_BAMS }           from './modules/local/04_pbmerge/main.nf';
include { CHECK_AND_MOVE_BAM }     from './modules/local/04_check_and_move_bam/main.nf';
include { ISOSEQ_CLUSTER }         from './modules/local/05_isoseq_cluster/main.nf';
include { PBMM2_ALIGN }            from './modules/local/06_align/main.nf';
include { ISOSEQ_COLLAPSE }        from './modules/local/07_isoseq_collapse/main.nf';
include { SQANTI_QC_AND_FILTER }   from './modules/local/08_sqanti_qc_filter/main.nf';
include { CPAT_PREDICTION }        from './modules/local/09_cpat/main.nf';
include { SIX_FRAME_TRANSLATION }  from './modules/local/10_six_frame_translation/main.nf';
include { TRANSCRIPTOME_SUMMARY }  from './modules/local/11_transcriptome_summary/main.nf';


workflow {

    // ---- REFERENCE GENOME INPUTS ---- //
    // This chunk of code 

    def reference_gtf   = file("${params.reference_genomedir}/gencode.v42.annotation.gtf")
    def reference_genome_fasta = file("${params.reference_genomedir}/GRCh38.p13.genome.fa")
    def reference_pc_fasta = file("${params.reference_genomedir}/gencode.v42.pc_transcripts.fa")
    def parsed_reference_dir   = file("${workflow.projectDir}/reference_genome/reference_db_for_analysis")

    // Generate reference tables and clustering database
    GENERATE_REFERENCE_TABLES_AND_DB(
        tuple(reference_gtf, reference_pc_fasta, parsed_reference_dir)
    )
    .reference_outputs
    .view { ref_output -> 
        println "Reference files created in: ${ref_output}" 
    }
    .set { reference_outputs }

    // --- END OF REFERENCE GENOME INPUTS --- // 


    // --- PARSED METADATA FILE --- // 
    Channel
        .fromPath(params.input_metadata)
        .flatMap { file ->
            file.text
                .split('\n')
                .drop(1)
                .findAll { it.trim() }
                .collect { line ->
                    def (bam, sample, condition) = line.trim().split(',').collect { it.trim() }
                    tuple(sample, condition, "${params.base_bam_dir}/${bam}")
                }
        }
        .map { sample, condition, bam_path ->
            tuple(sample, condition, file(bam_path))  // Wrap the path in the file() function
        }
        .filter { sample, condition, bam_file ->
            if (!bam_file.exists()) {
                println "WARN: Missing file: ${bam_file}"
                return false
            }
            return true
        }
        .groupTuple()
        .set { grouped_bams }

    PRINT_SAMPLE_BAM_COUNT(grouped_bams)
        .set { sample_bam_count_files }

    // Process merged bams
    grouped_bams
        .filter { sample, condition, bams -> bams.size() > 1 }
        .map { sample, condition, bams -> tuple(sample, condition, bams) }
        .set { to_merge }

    MERGE_BAMS(to_merge)
        .set { merged_bams }

    // Process original bams (only one per sample)
    grouped_bams
        .filter { sample, condition, bams -> bams.size() == 1 }
        .map { sample, condition, bams -> tuple(sample, condition, bams[0]) }
        .set { original_bams }

    // Combine merged and original
    merged_bams.mix(original_bams)
        .set { _all_bams }


    // *** START : USE PBMERGE OR NOT *** //
    // Function : This code chunk will merge all the bam files across all conditions in the sample metadata //

           _all_bams.view { sample_id, condition, bam_file ->
                println "[DEBUG] _all_bams entry - Sample: $sample_id, Condition: $condition, BAM: $bam_file"
            }
                if (params.use_pbmerge) {
                    _all_bams
                        .map { sample_id, condition, bam_file -> bam_file }
                        .collect()
                        .map { bam_files ->
                            tuple("pbmerged_bams", "all_conditions", bam_files)
                        }
                        .set { merged_pbmerge_input }

                    merged_pbmerge_input.view { sample_id, condition, bam_files ->
                        println "[DEBUG] Merging ${bam_files.size()} BAM files for ${sample_id}"
                        bam_files.each { println "  File: $it" }
                    }

                    PBMERGE_BAMS(merged_pbmerge_input)
                        .map { sample_id, condition, merged_bam_file, pbi_file ->
                            tuple(sample_id, condition, merged_bam_file, pbi_file)
                        }
                        .set { all_bams }

                } else {
                    _all_bams
                        .map { sample_id, condition, bam_file ->
                            def pbi_file = file("${bam_file}.pbi")
                            tuple(sample_id, condition, bam_file, pbi_file)
                        }
                        .set { all_bams }
                }

                all_bams.view { sample_id, condition, bam_file, pbi_file ->
                    println "[FINAL OUTPUT] Sample: ${sample_id}, Condition: ${condition}"
                    println "  BAM: ${bam_file}"
                    println "  PBI: ${pbi_file}"
                }

    // *** END : USE PBMERGE OR NOT *** // 



    // *** START : CHECK AND MOVE BAMS TO THEIR SAMPLE-SPECIFIC FOLDERS *** //
        
        CHECK_AND_MOVE_BAM(all_bams)

    // *** END : CHECK AND MOVE BAMS TO THEIR SAMPLE-SPECIFIC FOLDERS *** //


    // *** START : CLUSTER FLNC READS USING ISOSEQ *** //
    // // 
            ISOSEQ_CLUSTER(
                all_bams.map { sample_id, condition, bam_file , pbi_file ->
                    tuple(sample_id, condition, bam_file)
                }
            )
            .view { sample_id, condition, clustered_bam, clustered_bai ->  // ✅ fixed arity
                println "Clustering using IsoSeq for sample ${sample_id} (${condition}) completed. Output: ${clustered_bam}"
            }
            .set { isoseq_bam_files }

            isoseq_bam_files
                .map { sample_id, condition, clustered_bam , clustered_bai ->
                    tuple(sample_id, condition, clustered_bam, reference_genome_fasta)
                }
                .set { pbmm2_input }
                
    // *** END : CLUSTER FLNC READS USING ISOSEQ *** // 

     // *** START : ALIGN FLNC READS USING PBMM2 *** //
    // //
        isoseq_bam_files
            .map { sample_id, condition, clustered_bam , clustered_bai ->
                tuple(sample_id, condition, clustered_bam, reference_genome_fasta )
            }
            .set { pbmm2_input }

        PBMM2_ALIGN(pbmm2_input)
            .view { sample_id, condition, bam, bai ->
                println "Alignment done for ${sample_id} (${condition}): ${bam}"
            }
            .set { aligned_bams }

    // *** END : ALIGN FLNC READS USING PBMM2 *** // 

    // *** START : COLLAPSE FLNC READS USING ISOSEQ *** //
    // //
            aligned_bams
                .map { sample_id, condition, aligned_bam, bai ->
                    tuple([sample_id, condition], aligned_bam)
                }
                .join(
                    all_bams.map { sample_id, condition, bam_for_analysis, pbi_file ->
                        tuple([sample_id, condition], bam_for_analysis, pbi_file)
                    }
                )
                .map { key, aligned_bam, bam_for_analysis, pbi_file ->
                    def (sample_id, condition) = key
                    def output_dir = file("data/${sample_id}")
                    tuple(sample_id, condition, aligned_bam, bam_for_analysis, pbi_file, output_dir)
                }
                .set { isoseq_collapse_input }

            ISOSEQ_COLLAPSE(isoseq_collapse_input)
                .view { sample_id, condition, abundance, fasta, flnc_count, gff, group, read_stat, report_json ->
                    println "IsoSeq collapse done for ${sample_id} (${condition})"
                    println "GFF:         ${gff}"
                    println "FASTA:       ${fasta}"
                    println "Abundance:   ${abundance}"
                    println "Group:       ${group}"
                    println "FLNC Count:  ${flnc_count}"
                    println "Read Stat:   ${read_stat}"
                    println "Report JSON: ${report_json}"
                }
                .set { collapsed_gff_files }

    // *** END : COLLAPSE FLNC READS USING ISOSEQ ***
    

    // *** START : PERFROM QC AND FILTER READS USING SQANTI *** //
    // Prepare inputs for SQANTI QC and filtering //
            collapsed_gff_files
                .map { sample_id, condition,
                       collapsed_abundance,
                       collapsed_fasta,
                       collapsed_flnc_count,
                       collapsed_gff,
                       collapsed_group,
                       collapsed_read_stat,
                       collapsed_report ->

                    def output_dir     = file("data/${sample_id}")
                    def output_prefix  = "${sample_id}_sqanti"

                    tuple(
                        collapsed_gff,
                        reference_gtf,
                        reference_genome_fasta,
                        parsed_reference_dir,
                        collapsed_abundance,
                        output_prefix,
                        output_dir
                    )
                }
                .set { sqanti_qc_inputs }

            SQANTI_QC_AND_FILTER(sqanti_qc_inputs)
                .view { row -> println "SQANTI QC and filter completed for sample ${row[5]}. Outputs in: ${row[6]}" }
                .set { sqanti_qc_outputs }

    // *** END : PERFROM QC AND FILTER READS USING SQANTI *** //

    
    // *** START : CPAT  *** // 
    // Function : //
            sqanti_qc_outputs
                .map { files ->
                    def fasta = files.find { it.name ==~ /filtered_.+?_sqanti_corrected\.fasta/ }
                    def fname = fasta.getName()
                    def matcher = fname =~ /filtered_(.+?)_sqanti_corrected\.fasta/
                    def sample_id = matcher ? matcher[0][1] : null
                    def output_dir = file("data/${sample_id}")
                    def condition = null  // set later if needed
                    tuple(sample_id, condition, fasta, output_dir)
                }
                .set { cpat_inputs }

            CPAT_PREDICTION(cpat_inputs).set { cpat_output }

    // *** END : CPAT *** // 
    // CPAT is a bioinformatics tool to predict an RNA’s coding probability based on the RNA sequence characteristics //  
    // To achieve this goal, CPAT calculates scores of sequence-based features from a set of known protein-coding genes and background set of non-coding genes //


    // -- SIX FRAME TRANSLATION --- //
    // -- TRANSCRIPTOME SUMMARY --- //
    // -- ORF CALLING -- //
    // -- REFINE ORF DATABASE -- //
    // -- MAKE CDS GTF -- //
    // -- RENAME CDS TO EXON --//
    // -- SQANTI PROTEIN -- //
    // -- 5P UTR -- //
    // -- PROTEIN CLASSIFICATION -- //


}