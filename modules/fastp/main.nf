runmode = params.instrument

process FASTP {
    maxForks 8
    tag "$id"
    label 'process_high'
    
    publishDir "trimmed_fastqs", mode: "symlink"        , overwrite: true


    input:
        tuple val(id), path(reads), val(config)
    
    output:
        tuple val(id), path("*trimmed.fq.gz"), path("*.fasta"), val(config)         , emit: trimmed_fqs
             
        
    script:

    if ( runmode == "nova" ){
        
        """
        fastp \
        -z 4 -w 16 \
        --length_required 10 --qualified_quality_phred 20 \
        --trim_poly_g \
        -i ${reads} \
        -o ${id}_trimmed.fq.gz \
        -h ${id}.fastp.html \
        -j ${id}.fastp.json
    
        gunzip ${id}_trimmed.fq.gz
        fastq2fasta.pl ${id}_trimmed.fq > ${id}.fasta
        gzip ${id}_trimmed.fq

        """

    }

    else if ( runmode == "hiseq" ){

        """
            fastp \
            -z 4 -w 16 \
            --length_required 10 --qualified_quality_phred 20 \
            -i ${reads} \
            -o ${id}_trimmed.fq.gz \
            -h ${id}.fastp.html \
            -j ${id}.fastp.json
        

        gunzip ${id}_trimmed.fq.gz
        fastq2fasta.pl ${id}_trimmed.fq > ${id}.fasta
        ${id}_trimmed.fq
        
        """


    } else {

        error "Runmode ${runmode} is not supported"
        exit 0
    }



}

