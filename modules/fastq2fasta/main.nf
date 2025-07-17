process FASTQ2FASTA {

    maxForks 8
    tag "$id"
    label 'process_low'
    
    publishDir "trimmed_fasta", mode: "symlink", overwrite: true


    input:
        tuple val(id), path(reads), val(config)
    
    output:
        tuple val(id), path("*fasta"), val(config)     , emit: trimmed_fasta
             


    script:


    """

        gunzip ${reads}

        fastq2fasta.pl *trimmed.fq > ${id}.fasta


    """


}