#!/usr/bin/env nextflow

/*
========================================================================================
                         uct-cbio/bacterial_transcriptomics
========================================================================================
 Based on the nf-core Analysis Pipeline.
 #### Homepage / Documentation
 LINK TO GIT REPOSITORY
----------------------------------------------------------------------------------------


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    =========================================

    Bacterial transcriptomics pipeline

    Developed by the bioinformatics support team at the University of Cape Town

    =========================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow main.nf --reads sample_sheet.csv --genome <path to fasta> -with-docker <docker image>

    or

    nextflow main.nf --reads sample_sheet.csv --genome <path to fasta> -with-singularity <singularity image>


    Mandatory arguments:
        --reads                       The sample sheet containing the paths to the fastq files, as well as sample names.
        --genome                      The reference genome to be used in fasta format. Also acts as an outgroup.
        --gff                         Path to GFF3 file OR (see next arg)
        --gtf                         Path to GTF file
        --transcripts                 Will try to automate transcripts extraction based on gff / gtf
        -profile                      Hardware config to use. local / uct_hex

    Optional arguments:
        --minQuality                  The minimum quality to be passed to vcf-tools for filtering variants.
        --vcf_qual_cutoff             Soon to be removed
        --aligner                     Currently bwa-mem, star
        --quantification              Currently Salmon
        --srst_min_gene_cov           Minimum coverage for srst2 (default 90)
        --srst_max_gene_divergence    Maximum %divergence cutoff for gene reporting (default 10)


    Other arguments:
        --SRAdir                      The directory where reads downloaded from the SRA will be stored
        --vf_db                       Whether to look for virulence factors
        --outdir                      The output directory where the results will be saved
        --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.




    """.stripIndent()
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name             = false
params.project          = false
params.email            = false
params.plaintext_email  = false

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}


// Elvis syntax
// Reference index path configuration
// Define these here - after the profiles are loaded with the iGenomes paths
//params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
//params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
//params.gff = params.genome ? params.genomes[ params.genome ].gff ?: false : false
//params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
//params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false


Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
       .into{ch_where_trim_galore; ch_where_star; ch_where_hisat2; ch_where_hisat2_sort}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}



// Preset trimming options
if (params.pico){
    clip_r1 = 3
    clip_r2 = 0
    three_prime_clip_r1 = 0
    three_prime_clip_r2 = 3
    forward_stranded = true
    reverse_stranded = false
    unstranded = false
}


//Validate inputs
if ( params.genome == false ) {
    exit 1, "Must set a reference genome fasta file (--genome)"
}

if ( params.transcripts == false ) {
    exit 1, "Must set a reference transcripts fasta file (--transcripts)"
}

if ( params.reads == false ) {
    exit 1, "Must set the path to the sample file (--reads) in csv format"
}

// SNPeff needs a gff, all else gtf
if( params.gtf ){
    Channel
        .fromPath(params.gtf)
        .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
        .into { gtfFile }
} else if( params.gff ){
    Channel
        .fromPath(params.gff)
        .ifEmpty { exit 1, "GFF annotation file not found: ${params.gff}" }
        .into { gffFile }
} else {
    exit 1, "No GTF or GFF3 annotation specified!"
}




// Has the run name been specified by the user?
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}




log.info """\
Bacterial transcriptomics pipeline v0.1
================================
genome      : $params.genome
reads       : $params.reads
transcripts : $params.transcripts
Output      : $params.outdir
SRA dir     : $params.SRAdir
"""

/*
 *  Parse the input parameters
 */

genome_file             = file(params.genome)
transcripts_file        = file(params.transcripts)
sample_sheet            = file(params.reads)
reads_ch                = Channel.fromFilePairs(params.reads)
threads                 = 4
aligner                 = params.aligner
SRAdir                  = params.SRAdir
file_ext                = 'int'
quantification          = params.quantification

// Skip some QC not directly related to transcriptomics, move to QC pipeline
params.skip_qc          = false
params.skip_rseqc       = false
params.skip_preseq      = true
params.skip_multiqc     = false
params.subsampFilesizeThreshold = 10000000000

get_software_versions   = false


// Read clipping and strandedness
clip_r1                 = params.clip_r1
clip_r2                 = params.clip_r2
three_prime_clip_r1     = params.three_prime_clip_r1
three_prime_clip_r2     = params.three_prime_clip_r2
forward_stranded        = params.forward_stranded
reverse_stranded        = params.reverse_stranded
unstranded              = params.unstranded


// SRST and MLST parameters
srst_min_gene_cov           = params.srst_min_gene_cov
srst_max_gene_divergence    = params.srst_max_gene_divergence


// From https://pubmlst.org/data/dbases.xml             <----------------------- This needs a tweak to be generalised
mlst_species_srst2 = "mycobacteria"
mlst_definitions_srst2 = "mycobacteria"
mlst_seperator_srst2 = "_"


// Header log info
log.info nfcoreHeader()
def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'



/*
 *  ------------------------------------- PREPROCESSING -------------------------------------
 *
 * Convert GFF3 to GTF
 */

if(params.gff){
  process convertGFFtoGTF {
      tag "$gff"

      input:
      file gff from gffFile

      output:
      file "${gff.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts
      file "${gff.baseName}.gff3" into snpeff_gff

      script:
      """
      gffread $gff -T -o ${gff.baseName}.gtf
      """
  }
} else {
  process convertGTFtoGFF {

  input:
  file gtf from gtfFile

  output:

  file "${gtf.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12, gtf_star, gtf_dupradar, gtf_featureCounts
  file "${gtf.baseName}.gff" into snpeff_gff

  script:
  """
  gffread $gtf -o ${gtf.baseName}.gff
  """

  }

}

/*
 * Build BED12 file
 */
if(!params.bed12){
    process makeBED12 {
        tag "$gtf"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file gtf from gtf_makeBED12

        output:
        file "${gtf.baseName}.bed" into bed_rseqc, bed_genebody_coverage

        script: // This script is bundled with the pipeline, in nfcore/rnaseq/bin/
        """
        gtf2bed $gtf > ${gtf.baseName}.bed
        """
    }
}




/*
 * Parse software version numbers -- The scrape_software_versions.py needs updating, REMOVING FOR NOW
 */

if(get_software_versions){
    process '0A_get_software_versions' {

        output:
        file 'software_versions_mqc.yaml' into software_versions_yaml

        script:
        """
        echo $workflow.manifest.version &> v_ngi_rnaseq.txt
        echo $workflow.nextflow.version &> v_nextflow.txt
        fastqc --version &> v_fastqc.txt                        # Not working, works in Docker
        cutadapt --version &> v_cutadapt.txt                    # Working
        trim_galore --version &> v_trim_galore.txt              # Working
        #bwa &> v_bwa.txt                                        # Working, not parsing
        #preseq &> v_preseq.txt                                  # Not working libgsl.so.0: cannot open shared object file also in docker
        read_duplication.py --version &> v_rseqc.txt            # Working
        echo \$(bamCoverage --version 2>&1) > v_deeptools.txt       # unknown
        picard MarkDuplicates --version &> v_markduplicates.txt  || true    # Not working, not in docker either
        samtools --version &> v_samtools.txt                    # Working
        multiqc --version &> v_multiqc.txt                      # Working
        #scrape_software_versions.py &> software_versions_mqc.yaml   # unknown
        echo "this" &> software_versions_mqc.yaml
        """
    }
}

/*
 * ------------------------------------- ANALYSIS PART 1: Data preparation -------------------------------------
 *
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process '1A_prepare_genome_samtools' { 
  tag "$genome.baseName"
  
  input: 
      file genome from genome_file 
 
  output: 
      file "${genome}.fai" into genome_index_ch  
  
  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  picard -XX:ParallelGCThreads=5 -Xmx16G -Xms16G CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
  """
}

/*
 * Process 1C: Create a FASTA genome sequence dictionary for BWA
 */

process '1C_prepare_genome_bwa' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome}.amb" into genome_bwa_amb
      file "${genome}.ann" into genome_bwa_ann
      file "${genome}.bwt" into genome_bwa_bwt
      file "${genome}.pac" into genome_bwa_pac
      file "${genome}.sa" into genome_bwa_sa

  script:
  """
  bwa index $genome 
  """
}

/*
 * Process 1D: Prepare and download samples as per sample sheet
 */

process '1D_prepare_samples' {

  publishDir "$params.SRAdir", mode: "link"

  input:
      file samples from sample_sheet
  output:
      file "sample_sheet_new.csv" into newSampleSheet
      file "sample_sheet_new.csv" into newSampleSheetFastQC
      file "*.fastq" optional true into SRA_new_reads
  script:
  """
  echo $params.SRAdir > out.txt
  process_samples.py -i $samples -f $params.SRAdir
  """
}


newSampleSheet
  .splitCsv(header:true)
  .map { row-> tuple(row.number, file(row.R1), file(row.R2), row.isolate) }
  .set { newSampleChannel }

newSampleSheetFastQC
  .splitCsv(header:true)
  .map { row-> tuple(row.number, file(row.R1), file(row.R2), row.isolate) }
  .set { newSampleChannelFastQC }




/*
 * Process 1E: FastQC
 */

 process '1E_fastqc' {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set number, file(R1), file(R2), isolate from newSampleChannelFastQC

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    #  MiSeq file naming convention (samplename_S1_L001_[R1]_001)

    mv $R1 sample_${isolate}_R1_001.fq.gz
    mv $R2 sample_${isolate}_R2_001.fq.gz

    fastqc -q sample_${isolate}_R1_001.fq.gz sample_${isolate}_R2_001.fq.gz
    """
}



/*
 * Process 1F: Trim Galore!
 */

process '1F_trim_galore' {
    label 'low_memory'
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set number, file(R1), file(R2), isolate from newSampleChannel

    output:
    file "*_R1_001.fq.gz" into forwardTrimmedAlignment
    file "*_R2_001.fq.gz" into reverseTrimmedAlignment
    file "*_R1_001.fq.gz" into forwardTrimmedQuant
    file "*_R2_001.fq.gz" into reverseTrimmedQuant
    file "*_R1_001.fq.gz" into forward_trimmed_reads_for_srst2
    file "*_R2_001.fq.gz" into reverse_trimmed_reads_for_srst2

    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    val "$number" into sampleNumber_srst2
    val "$number" into sampleNumber
    val "$isolate" into sampleIsolate
    val "$isolate" into sampleIsolateQuant
    val "$isolate" into sampleIsolate_srst2
    set number, file("*_R1_001.fq.gz"), file("*_R2_001.fq.gz") into vf_read_pairs

    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $R1 $R2
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $R1 $R2

        #  MiSeq file naming convention (samplename_S1_L001_[R1]_001)

        rename 's/fastq.gz/fq.gz/' *.fastq.gz

        rename 's/_val_1/_R1_001/' *.fq.gz
        rename 's/_val_2/_R2_001/' *.fq.gz
        """
    }
}


/*
 * ------------------------------------- ANALYSIS PART 2: Alignment -------------------------------------
 *
 * Process 2B: Align reads to the reference genome  ------------------------------------- EDIT TO INCLUDE OTHER ALIGNERS
 *
 */

process '2A_read_mapping' {
  input:
    file forwardTrimmedAlignment
    file reverseTrimmedAlignment
    val sampleNumber
    val sampleIsolate
    file genome from genome_file
    file genome_bwa_amb
    file genome_bwa_ann
    file genome_bwa_bwt
    file genome_bwa_pac
    file genome_bwa_sa
  output:
    file "sample_${sampleIsolate}.sorted.bam" into bamfiles
    file "sample_${sampleIsolate}.sorted.bai" into bamindexfiles
    file "sample_${sampleIsolate}.sorted.bam" into bam_rseqc
    file "sample_${sampleIsolate}.sorted.bai" into bamindexfiles_rseqc
    file "sample_${sampleIsolate}.sorted.bam" into bam_preseq
    file "sample_${sampleIsolate}.sorted.bam" into bam_forSubsamp
    file "sample_${sampleIsolate}.sorted.bam" into bam_skipSubsamp
    file "sample_${sampleIsolate}.sorted.bam" into bam_featurecounts
  script:
  if( aligner == 'bwa-mem' )
    """
    bwa mem $genome $forwardTrimmedAlignment $reverseTrimmedAlignment | samtools sort -O BAM -o sample_${sampleIsolate}.sorted.bam
    samtools index sample_${sampleIsolate}.sorted.bam sample_${sampleIsolate}.sorted.bai
    """

  else
    error "Invalid aligner: ${aligner}"

}



/*
 * Process 2B: RSeQC analysis -- Appears working
 */

process '2B_rseqc' {
    label 'high_memory'
    tag "${bam_rseqc.baseName - '.sorted'}"
    publishDir "${params.outdir}/rseqc" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.curves.pdf") > 0)       "geneBodyCoverage/$filename"
            else if (filename.indexOf("geneBodyCoverage.r") > 0)                "geneBodyCoverage/rscripts/$filename"
            else if (filename.indexOf("geneBodyCoverage.txt") > 0)              "geneBodyCoverage/data/$filename"
            else if (filename.indexOf("log.txt") > -1) false
            else filename
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_rseqc
    file index from bamindexfiles_rseqc
    file bed12 from bed_rseqc.collect()

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    infer_experiment.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.infer_experiment.txt
    bam_stat.py -i $bam_rseqc 2> ${bam_rseqc.baseName}.bam_stat.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    read_distribution.py -i $bam_rseqc -r $bed12 > ${bam_rseqc.baseName}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc.baseName}.read_duplication

    geneBody_coverage.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    mv log.txt ${bam_rseqc.baseName}.rseqc.log.txt

    # Not applicable for bacteria
    #junction_annotation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12
    #junction_saturation.py -i $bam_rseqc -o ${bam_rseqc.baseName}.rseqc -r $bed12 2> ${bam_rseqc.baseName}.junction_annotation_log.txt

    """
}




/*
 * Process 2F: Mark duplicate reads  --- Edit file naming
 */

process '2F_mark_duplicates' {
  label 'high_memory'
  publishDir "${params.outdir}/picard", mode: "copy"

  input:
    file sample_bam from bamfiles
  output:
    set file("${sample_bam.baseName}.dedup.bam"), file("${sample_bam.baseName}.dedup.bam.bai") into dedup_bamfiles
    set file("${sample_bam.baseName}.dedup.bam"), file("${sample_bam.baseName}.dedup.bam.bai") into dupradar_bamfiles
    file "${sample_bam.baseName}.txt" into picard_results
  script:
    """
    picard MarkDuplicates INPUT=$sample_bam OUTPUT=${sample_bam.baseName}.dedup.bam METRICS_FILE=${sample_bam.baseName}.txt ASSUME_SORTED=true REMOVE_DUPLICATES=false
    samtools index ${sample_bam.baseName}.dedup.bam
    """
}



/*
 * Process 2G: dupradar
 */

process '2G_dupradar' {
    label 'low_memory'
    tag "${bamfile.baseName}"
    publishDir "${params.outdir}/dupradar", mode: "link", overwrite: true

    input:
      set file(bamfile), file(bamindex) from dupradar_bamfiles
      file gtf from gtf_dupradar
    output:
      file "*.{pdf,txt}" into dupradar_results

    script:

    // This script was bundled with the nfcore/rnaseq pipeline in nfcore/rnaseq/bin/
    def dupradar_direction = 0
    if (forward_stranded && !unstranded) {
        dupradar_direction = 1
    } else if (reverse_stranded && !unstranded){
        dupradar_direction = 2
    }
    def paired = params.singleEnd ? 'single' :  'paired'

    """
    dupRadar.r $bamfile $gtf $dupradar_direction $paired ${task.cpus}
    """
}


/*
 * ------------------------------------ ANALYSIS PART 3: Virulence and DB analysis ------------------------------------
 *
 * Process 3A: srst2 (run per sample)
 * https://github.com/kviljoen/uct-srst2/blob/master/main.nf
 */

process '3A_srst2' {
    tag { "srst2.${sampleNumber_srst2}" }
    publishDir "${params.outdir}/srst2_mlst", mode: "copy"
    label 'high_memory'

    input:
    file forward_trimmed_reads_for_srst2
    file reverse_trimmed_reads_for_srst2
    val sampleNumber_srst2
    val sampleIsolate_srst2
    val srst_min_gene_cov
    val srst_max_gene_divergence
    val mlst_species_srst2
    val mlst_definitions_srst2
    val mlst_seperator_srst2

    output:
    file("${sampleIsolate_srst2}_srst2__mlst*")

    script:
    geneDB = params.gene_db ? "--gene_db $gene_db" : ''
    mlstDB = params.mlst_db ? "--mlst_db $mlst_db" : ''
    mlstdef = params.mlst_db ? "--mlst_definitions $mlst_definitions" : ''
    mlstdelim = params.mlst_db ? "--mlst_delimiter $params.mlst_delimiter" : ''
    mlstfasta = mlst_species_srst2.replace(" ", "_")

    """
    # /samtools-0.1.18/
    export SRST2_SAMTOOLS="/samtools-0.1.18/samtools"
    getmlst.py --species "${mlst_species_srst2}"
    srst2 --output ${sampleIsolate_srst2}_srst2 --input_pe $forward_trimmed_reads_for_srst2 $reverse_trimmed_reads_for_srst2 --mlst_db ${mlstfasta}.fasta --mlst_definitions profiles_csv --mlst_delimiter '_' --min_coverage $srst_min_gene_cov --max_divergence $srst_max_gene_divergence
    """
}

/*
 * Process 3B: AMR Resistance
 */

if( params.amr_db ) {
    /*
     * Build resistance database index with Bowtie2
     */
    process '3B_BuildAMRIndex' {
        tag { "${amr_db.baseName}" }

        input:
            file amr_db

            output:
            file 'amr.index*' into amr_index

            """
            bowtie2-build $amr_db amr.index --threads ${threads}
            """
    }

    /*
         * Align reads to resistance database with Bowtie2
         */
    process '3B_AMRAlignment' {
            publishDir "${params.outdir}/Alignment", pattern: "*.bam"

            tag { dataset_id }

            input:
            set dataset_id, file(forward), file(reverse) from amr_read_pairs
            file index from amr_index.first()

            output:
            set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files
            set dataset_id, file("${dataset_id}_amr_alignment.bam") into amr_bam_files

            """
            bowtie2 -p ${threads} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
            samtools view -bS ${dataset_id}_amr_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_amr_alignment.bam
            """
    }

    process '3B_AMRResistome' {
            publishDir "${params.outdir}/Resistome"

            tag { dataset_id }

            input:
            file amr_db
            set dataset_id, file(amr_sam) from amr_sam_files

            output:
            set dataset_id, file("${dataset_id}_amr_gene_resistome.tsv") into amr_gene_level

            """
            csa -ref_fp ${vf_db} -sam_fp ${vf_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_amr_gene_resistome.tsv"
            """
    }
}

/*
 * Process 3C: Virulence factor analysis
 */

if( params.vf_db ) {
    /*
    * Build resistance database index with Bowtie2
    */
    process '3C_BuildVFIndex' {
        tag { "Building index" }

        input:

        output:
        file 'vf.index*' into vf_index
        file 'VFDB_setB_nt.fa' into vf_fa

        """
        wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
        gunzip VFDB_setB_nt.fas.gz
        mv VFDB_setB_nt.fas VFDB_setB_nt.fa
        sed -i 's/(/_/g' VFDB_setB_nt.fa
        sed -i 's/)/_/g' VFDB_setB_nt.fa
        bowtie2-build VFDB_setB_nt.fa vf.index
        """
    }
    /*
         * Align reads to virulence factor database with Bowtie2
         */
    process '3C_VFAlignment' {
            publishDir "${params.outdir}/Alignment", pattern: "*.bam"

            tag { dataset_id }

            input:
            set dataset_id, file(forward), file(reverse) from vf_read_pairs
            file index from vf_index.first()
            file vf_fasta from vf_fa

            output:
            set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files
            set dataset_id, file("${dataset_id}_vf_alignment.bam") into vf_bam_files

            """
            bowtie2 -p ${threads} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
            samtools view -bS ${dataset_id}_vf_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_vf_alignment.bam
            """
    }

    process '3C_VFResistome' {
            publishDir "${params.outdir}/Resistome"

            label 'high_memory'
            tag { dataset_id }

            input:
            file vf_db from vf_fa
            set dataset_id, file(vf_bam) from vf_bam_files

            output:
            set dataset_id, file("${dataset_id}_raw_wgs_metrics.txt") into vf_gene_level

            """
            picard CollectWgsMetrics I=$vf_bam O=${dataset_id}_raw_wgs_metrics.txt R=${vf_db} INCLUDE_BQ_HISTOGRAM=true
            """
    }
}

/*
 * Process 3D: Plasmid resistome analysis - CHECK WHERE THESE FILES ARE PULLED FROM
 */

if( params.plasmid_db ) {
    /*
         * Build plasmid index with Bowtie2
         */
    process '3D_BuildPlasmidIndex' {
        tag { "${plasmid_db.baseName}" }

        input:
            file plasmid_db

            output:
            file 'plasmid.index*' into plasmid_index

            """
            bowtie2-build $plasmid_db plasmid.index --threads ${threads}
            """
    }
    /*
         * Align reads to plasmid database with Bowtie2
         */
    process '3D_PlasmidAlignment' {
            publishDir "${params.outdir}/Alignment", pattern: "*.bam"

            tag { dataset_id }

            input:
            set dataset_id, file(forward), file(reverse) from plasmid_read_pairs
            file index from plasmid_index.first()

            output:
            set dataset_id, file("${dataset_id}_plasmid_alignment.sam") into plasmid_sam_files
            set dataset_id, file("${dataset_id}_plasmid_alignment.bam") into plasmid_bam_files

            """
            bowtie2 -p ${threads} -x plasmid.index -1 $forward -2 $reverse -S ${dataset_id}_plasmid_alignment.sam
            samtools view -bS ${dataset_id}_plasmid_alignment.sam | samtools sort -@ ${threads} -o ${dataset_id}_plasmid_alignment.bam
            """
    }

    process '3D_PlasmidResistome' {
            publishDir "${params.outdir}/Resistome"

            tag { dataset_id }

            input:
            file plasmid_db
            set dataset_id, file(plasmid_sam) from plasmid_sam_files

            output:
            set dataset_id, file("${dataset_id}_plasmid_gene_resistome.tsv") into plasmid_gene_level

            """
            csa -ref_fp ${plasmid_db} -sam_fp ${plasmid_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_plasmid_gene_resistome.tsv"
            """
    }
}




/*
 * ------------------------------------ ANALYSIS PART 4: Quantification ------------------------------------
 *
 * Process 4A: Quantification of the reads
 *
 */


process '4A_quantify_reads' {

  publishDir "${params.outdir}/salmon", mode: "link", overwrite: true

  input:
    file genome from genome_file
    file transcripts from transcripts_file
    file forwardTrimmedQuant
    file reverseTrimmedQuant
    val sampleIsolateQuant
    file gtf from gtf_featureCounts
  output:
    file "*" into salmon_results

  script:
  if( quantification == 'salmon' )
    """
    salmon index -t $transcripts -i transcripts_index -k 31
    salmon quant \\
        --geneMap $gtf \\
        --threads $task.cpus \\
        -l A \\
        -i transcripts_index \\
        $genome \\
        --validateMappings \\
        -1 $forwardTrimmedQuant \\
        -2 $reverseTrimmedQuant \\
        -o $sampleIsolateQuant
    """
  else if( quantification == 'other-option' )
    """
    other-option
    """
  else
    error "Invalid quantification method: ${quantification}"
}


/*
 *  END OF PART 4
 *********/




/*
 *  END OF PART 5
 *********/


/*
 * ------------------------------------ ANALYSIS PART 6: MultiQC ------------------------------------
 *
 *
 */

process '6A_multiqc' {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    when:
    !params.skip_multiqc

    input:
    file multiqc_config from ch_multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    //file ('software_versions/*') from software_versions_yaml
    file ('picard/*') from picard_results

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}



def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/rnaseq v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
