process FASTQ2FASTA {

    maxForks 8
    tag "$id"
    label 'process_low'
    
    publishDir "trimmed_fasta", mode: "symlink", overwrite: true


    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("*fasta")       , emit: trimmed_fasta
             


    script:


    """

        gunzip ${reads[0]}

        fastq2fasta.pl ${reads[0]} > ${id}.fasta


    """


}