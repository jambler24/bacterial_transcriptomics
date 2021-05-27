# Bacterial RNA seq pipeline (Salmon)

This pipeline can be used to generate the salmon quant files required by deseq2 for differential expression analysis.

A example shell script run template is provided named "run_template.sh"

The current version has been thinned out to only cover read trimming, QC reports, and quantification with salmon.

This is due to compatibility issues, other tools may be added back in future. 

## Basic usage: 
The typical command for running the pipeline is as follows:

    nextflow run jambler24/bac_pangenome --reads sample_sheet.csv --genome refgenome.fa -profile ilifu
    Mandatory arguments:
      --reads                       Path to sample sheet
      --genome                      Path to reference genome against which the reads will be aligned (in fasta format) for use in QC steps.
      --gtf                         Path to the GTF formatted annotation file. Salmon does not work with mony of the gff formats.
      --transcripts                 Path to the transcripts fasta file. 
      -profile                      Hardware config to use. Currently profile available for ilifu and UCT's HPC 'uct_hex' - create your own if necessary
      
      
    Other arguments:
      --outdir                      The output directory where the results will be saved
      --SRAdir                      The directory where reads downloaded from the SRA will be stored
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name      


## Sample file
To allow for both local reads and reads from the [SRA](https://www.ncbi.nlm.nih.gov/sra) to be used, the pipeline has the 
ability to pull reads from the SRA based on the accession number (eg, [SRR5989977](https://www.ncbi.nlm.nih.gov/sra/SRX3145707[accn])). 

The 'number' column must contain a unique value. 

number | origin | replicate | isolate | R1 | R2
------------ | ------------- | ------------- | ------------- | ------------- | -------------
1 | genomic | 1 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
2 | genomic | 2 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
3 | genomic | 3 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
4 | genomic | 1 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
5 | genomic | 2 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
6 | genomic | 3 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
7 | genomic | 1 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
8 | genomic | 2 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
9 | genomic | 3 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
10 | genomic | 1 | H37Rv | SRR5989977 | 


In the above example, samples 1-9 are locally stored where sample 10 is a control sample from the SRA. 
Including the accession number in the R1 column will result in the reads from the SRA to be downloaded and used in the analysis. 
This must be exported to a csv file, with a comma ',' separating the columns:

    number,origin,replicate,isolate,R1,R2
    1,genomic,1,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    2,genomic,2,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    ...
    10,genomic,1,H37Rv,SRR5989977
