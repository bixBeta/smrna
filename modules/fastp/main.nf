runmode = params.instrument

process FASTPM {
    maxForks 8
    tag "$id"
    label 'process_high'
    
    publishDir "trimmed_fastqs", mode: "symlink", overwrite: true


    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("*gz")       , emit: trimmed_fqs
             
        
    script:

    if ( runmode == "nova" ){
        
        """
        fastp \
        -z 4 -w 16 \
        --length_required 10 --qualified_quality_phred 20 \
        --trim_poly_g \
        -i ${reads[0]} \
        -o ${id}_val_1.fq.gz \
        -h ${id}.fastp.html \
        -j ${id}.fastp.json
    
        """

    }

    else if ( runmode == "hiseq" ){

        """
            fastp \
            -z 4 -w 16 \
            --length_required 10 --qualified_quality_phred 20 \
            -i ${reads[0]} \
            -o ${id}_val_1.fq.gz \
            -h ${id}.fastp.html \
            -j ${id}.fastp.json
        
        """


    } else {

        error "Runmode ${runmode} is not supported"
        exit 0
    }



}

