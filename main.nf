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
mlst_species_srst2 = "Mycobacteria spp."
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

  file "${gtf.baseName}.gtf" into gtf_makeSTARindex, gtf_makeBED12, gtf_star
  file "${gtf.baseName}.gff" into snpeff_gff
  file "${gtf.baseName}.gtf" into gtf_featureCounts
  file "${gtf.baseName}.gtf" into gtf_dupradar
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
      set file(genome), file("${genome.baseName}.dict") into genome_quant_ch

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
    label 'high_memory'
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
        set val(number), file(R1), file(R2), isolate from newSampleChannel

    output:
        set val(number), file("*_R1_001.fq.gz"), file("*_R2_001.fq.gz") into QuantInput

        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

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

        #rename 's/_val_1/_R1_001/' *.fq.gz
        #rename 's/_val_2/_R2_001/' *.fq.gz

        rename 's/_trimmed/.trimmed/' *.fq.gz

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
  label 'high_memory'
  publishDir "${params.outdir}/salmon", mode: "link", overwrite: true

  input:
    set val(number), file(R1_reads), file(R2_reads) from QuantInput
    set genome_fasta, genome_dict from genome_quant_ch
    file transcripts from transcripts_file
  output:
    file "*" into salmon_results

  script:
  if( quantification == 'salmon' )
    """
    salmon index -t $transcripts -i transcripts_index -k 31
    salmon quant \\
        --geneMap $genome_fasta \\
        --threads $task.cpus \\
        -l A \\
        -i transcripts_index \\
        $genome_fasta \\
        -1 $R1_reads \\
        -2 $R2_reads \\
        -o $number
    """
  else if( quantification == 'other-option' )
    """
    other-option
    """
  else
    error "Invalid quantification method: ${quantification}"
}




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
    //file ('rseqc/*') from rseqc_results.collect().ifEmpty([])
    //file ('dupradar/*') from dupradar_results.collect().ifEmpty([])
    //file ('software_versions/*') from software_versions_yaml
    //file ('picard/*') from picard_results

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
