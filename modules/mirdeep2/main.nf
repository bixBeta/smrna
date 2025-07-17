process FASTQ2FASTA {

    maxForks 8
    tag "$id"
    label 'process_high'
    
    publishDir "trimmed_fasta", mode: "symlink", overwrite: true


    input:
        tuple val(id), path(reads), val(config)
    
    output:
        tuple val(id), path("*fasta"), val(config)     , emit: trimmed_fasta
             


    script:


    """

        mapper.pl $CONFIG -d -c -m -s ${id}.collapsed.fa


    """


}