process SMRNA_MQC_TABLES {

    tag "$pin"
    label 'process_mqc'

    publishDir "multiqc/custom_content", mode: "copy", overwrite: true


    input:
        val(pin)
        path(table)
        path(sheet)
        path('fastp/*')
        val(mirbase)

    output:
        path("mqc/*_mqc.yaml")      , emit: mqc_files


    script:

    // Set when QUANT ran. The miRBase panel and metric are then always emitted,
    // so a run that mapped nothing shows 0% instead of dropping the section.
    def mirbase_arg = mirbase ? "--mirbase" : ""

    """
        mkdir -p mqc

        smrna_mqc_tables.py ${table} \\
            --sample-sheet ${sheet} \\
            --outdir mqc \\
            --fastp-dir fastp \\
            --prefix smrna ${mirbase_arg}
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
        val(mqcgenome)

    output:
        path("${pin}_multiqc_report.html")      , emit: report
        path("*_data")                          , emit: data


    script:

    // Search the staged inputs by name rather than ".". Scanning the whole work
    // directory also picks up anything else that lands in it, so a section that
    // is no longer generated can reappear from a stale file.
    """
        export MQC_GENOME=${mqcgenome}

        multiqc -f \\
            -c ${mqc_config} \\
            --cl-config "custom_logo: ${logo}" \\
            -n ${pin}_multiqc_report.html \\
            fastp custom_content ${mqc_versions}
    """
}
