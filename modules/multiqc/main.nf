process SMRNA_MQC_TABLES {

    tag "$pin"
    label 'process_mqc'

    publishDir "multiqc/custom_content", mode: "copy", overwrite: true


    input:
        val(pin)
        path(table)
        path(sheet)

    output:
        path("mqc/*_mqc.yaml")      , emit: mqc_files


    script:

    """
        mkdir -p mqc

        smrna_mqc_tables.py ${table} \\
            --sample-sheet ${sheet} \\
            --outdir mqc \\
            --prefix smrna
    """
}


process MULTIQC {

    tag "$pin"
    label 'process_mqc'

    publishDir "multiqc", mode: "copy", overwrite: true


    input:
        val(pin)
        path('fastp/*')
        path('custom_content/*')
        path(mqc_versions)
        path(mqc_config)
        path(logo)

    output:
        path("${pin}_multiqc_report.html")      , emit: report
        path("*_data")                          , emit: data


    script:

    """
        multiqc -f \\
            -c ${mqc_config} \\
            --cl-config "custom_logo: ${logo}" \\
            -n ${pin}_multiqc_report.html \\
            .
    """
}
