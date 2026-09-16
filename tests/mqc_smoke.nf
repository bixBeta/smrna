// Runs the real SMRNA_MQC_TABLES and MULTIQC processes against fixtures, so CI
// checks what actually gets written and published rather than only parsing the
// DAG. Everything upstream (fastp, mapper.pl, quantifier.pl) is replaced by the
// fixture table, which is what those steps would have produced.

nextflow.enable.dsl=2

include { SMRNA_MQC_TABLES; MULTIQC } from '../modules/multiqc'

workflow {

    SMRNA_MQC_TABLES(
        channel.value('CI'),
        channel.value(file("$projectDir/fixtures/table.txt")),
        channel.value(file("$projectDir/fixtures/sample-sheet.csv")),
        channel.fromPath("$projectDir/fixtures/fastp/*.json").collect(),
        channel.value(file("$projectDir/../bin/smrna_mqc_tables.py")),
        true
    )

    MULTIQC(
        channel.value('CI'),
        channel.fromPath("$projectDir/fixtures/fastp/*.json").collect(),
        SMRNA_MQC_TABLES.out.mqc_files.collect(),
        channel.value(file("$projectDir/fixtures/software_versions_mqc.yml")),
        channel.value(file("$projectDir/../assets/multiqc_config.yaml")),
        channel.value(file("$projectDir/../img/trex-extended-logo.png")),
        'hsa'
    )
}
