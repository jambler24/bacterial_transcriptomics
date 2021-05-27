refGenomeFolder='/path/to/reference/genome/fasta/file/'
refGenomeFasta='ref_genome.fa'
refGenomeGTF='ref_annotation.gtf'
refGenomeSTAR='star/'
refGenomeTranscripts='ref_genome_transcripts.fa'

/cbio/soft/nextflow/nextflow run jambler24/bacterial_transcriptomics \
        -latest -resume \
        --reads '/path/to/sample_sheet.csv' \
        --quantification salmon \
        --genome $refGenomeFolder$refGenomeFasta \
        --gtf $refGenomeFolder$refGenomeGTF \
        --transcripts $refGenomeFolder$refGenomeTranscripts \
        --SRAdir /path/to/SRA/download/folder/ \
        -profile ilifu \
        --outdir /path/to/output/directory/