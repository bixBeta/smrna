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
        path("${id}.fastp.json")                                                    , emit: json
        path("${id}.fastp.html")                                                    , emit: html
        path("${id}_fastp_mqc_versions.yml")                                        , emit: versions


    script:

    // fastp only records its version in the JSON from v0.22 onward, so capture
    // it here instead of relying on MultiQC finding it there.
    def version_cmd = """
        FASTP_VERSION=\$(fastp --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1) || true
        [ -n "\$FASTP_VERSION" ] || FASTP_VERSION="unknown"
        printf 'fastp: ["%s"]\\n' "\$FASTP_VERSION" > ${id}_fastp_mqc_versions.yml
    """

    if ( runmode == "nova" ){

        """
        fastp \
        -z 4 -w 16 \
        --adapter_sequence ${params.adapter} \
        --length_required ${params.min_len} --qualified_quality_phred 20 \
        --trim_poly_g \
        -i ${reads} \
        -o ${id}_trimmed.fq.gz \
        -h ${id}.fastp.html \
        -j ${id}.fastp.json

        gunzip ${id}_trimmed.fq.gz
        fastq2fasta.pl ${id}_trimmed.fq > ${id}.fasta
        gzip ${id}_trimmed.fq
        ${version_cmd}
        """

    }

    else if ( runmode == "hiseq" ){

        """
            fastp \
            -z 4 -w 16 \
            --adapter_sequence ${params.adapter} \
            --length_required ${params.min_len} --qualified_quality_phred 20 \
            -i ${reads} \
            -o ${id}_trimmed.fq.gz \
            -h ${id}.fastp.html \
            -j ${id}.fastp.json


        gunzip ${id}_trimmed.fq.gz
        fastq2fasta.pl ${id}_trimmed.fq > ${id}.fasta
        gzip ${id}_trimmed.fq
        ${version_cmd}
        """


    } else {

        error "Runmode ${runmode} is not supported"
        exit 0
    }



}


// Alternative to FASTP for --fastp false: mapper.pl needs fasta, so the
// conversion still has to happen even when no trimming is wanted.
process FASTQ2FASTA {
    maxForks 8
    tag "$id"
    label 'process_high'

    publishDir "raw_fastas", mode: "symlink"            , overwrite: true


    input:
        tuple val(id), path(reads), val(config)

    output:
        tuple val(id), path("${id}_untrimmed.fq.gz"), path("*.fasta"), val(config)  , emit: trimmed_fqs


    script:

    """
        gunzip -c ${reads} > ${id}_untrimmed.fq
        fastq2fasta.pl ${id}_untrimmed.fq > ${id}.fasta
        gzip ${id}_untrimmed.fq
    """
}
