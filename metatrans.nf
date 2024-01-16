#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.readsPath = '/mnt/gwu/metatrans_nf/samples/*.R{1,2}.fastq.gz'
params.trimmomaticPath = '/home/gwu/bin//trimmomatic.jar'
params.adaptersPath = '/home/gwu/bin/trimmomatic/adapters/NexteraPE-PE.fa'
params.sourmashdbPath = '/mnt/gwu/metatrans_nf/reference_db.sbt.json'

// Define processes
process trimReads {
    container 'quay.io/biocontainers/trimmomatic'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*_trimmed.fastq.gz')

    script:
    """
    java -jar ${params.trimmomaticPath} PE -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${sample}_1_trimmed.fastq.gz ${sample}_1_unpaired.fastq.gz \\
        ${sample}_2_trimmed.fastq.gz ${sample}_2_unpaired.fastq.gz \\
        ILLUMINACLIP:${params.adaptersPath}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process removeRrna {
    container 'quay.io/biocontainers/ribodetector'
    
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*_rRNA_removed.fastq.gz')

    script:
    """
    ribodetector_cpu -t 8 -l 80 -i ${reads[0]} -o ${sample}_1_rRNA_removed.fastq.gz
    ribodetector_cpu -t 8 -l 80 -i ${reads[1]} -o ${sample}_2_rRNA_removed.fastq.gz
    """
}

process pairReads {
    container 'quay.io/biocontainers/fastp'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path('*_paired.fastq.gz')

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample}_1_paired.fastq.gz -O  ${sample}_2_paired.fastq.gz
    """
}

process assembleTranscripts {
    container 'quay.io/biocontainers/metaspades'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_assembled.fasta")

    script:
    """
    rnaspades.py -1 ${reads[0]} -2 ${reads[1]} -o ${sample}_spades
    mv ${sample}_spades/transcripts.fasta ${sample}_assembled.fasta
    """
}

process annotateGenes {
    container 'quay.io/biocontainers/prokka'
    publishDir '.', mode: 'copy', pattern: '*.faa'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.faa")

    script:
//    def sample = transcripts.simpleName - '_assembled.fasta'
    """
    /home/gwu/prokka/bin/prokka --metagenome ${reads} --centre 37 --compliant --outdir ${sample}_prokka_out
    mv ./${sample}_prokka_out/PROKKA*.faa ${sample}.faa
    """
}

process functionalAnnotation {
    container 'quay.io/biocontainers/diamond'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_diamond.out")

    script:
    """
    /home/gwu/diamond blastp -d /mnt/gwu/bac_test -q ${reads} -o ${sample}_diamond.out
    """
}

process identifySpecies {
    container 'quay.io/biocontainers/sourmash'
    publishDir '.', mode: 'copy', pattern: '*_species.txt'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_species.txt")

    script:
    """
    cat ${reads[0]} ${reads[1]} > ${sample}_combined.fastq.gz
    sourmash compute -k 31 --scaled 1000 ${sample}_combined.fastq.gz
    sourmash search ${sample}_combined.fastq.gz.sig ${params.sourmashdbPath} --output ${sample}_species.txt
    """
}

// Define workflow
workflow {
    read_pairs = Channel.fromFilePairs(params.readsPath)

    trimmed_reads = read_pairs
        | trimReads

    rna_removed_reads = trimmed_reads
        | removeRrna

    paired_reads = rna_removed_reads
        | pairReads

    assembled_transcripts = paired_reads
        | assembleTranscripts

    gene_annotations = assembled_transcripts
        | annotateGenes

    functional_annotations = gene_annotations
        | functionalAnnotation

    species_identification = rna_removed_reads
        | identifySpecies
}
