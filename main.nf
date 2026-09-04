nextflow.enable.dsl=2


// Params:

params.help                     = false

params.sheet                    = "sample-sheet.csv"
params.id                       = "TREX_ID"
params.genome                   = null
params.instrument               = "nova"

// Opt-in, matching the convention in the bixBeta nextflow repo. Without it the
// fastqs are only converted to fasta for mapper.pl, with no trimming.
params.fastp                    = false

// NEBNext Small RNA 3' SR Adaptor. The kit uses fixed adapters with no
// randomised ends or UMI, so nothing needs trimming off the read termini.
params.adapter                  = "AGATCGGAAGAGCACACGTCT"
params.min_len                  = 10

if( params.help ) {

log.info """
s  m  R  N  A  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/smrna -r ${workflow.revision ?: 'main'} < args ... >

Args:
    * --id             : TREx Project ID 
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -------------------------------------------
        Sample Sheet Example:    
        label   fastq1          config
        SS1     SS1_R1.fastq.gz 101
        SS2     SS2_R1.fastq.gz 102 
        .
        .
        . etc.
        -------------------------------------------
    * --fastp           : Invokes fastp trimming module < default: false >
                         Strongly recommended for smRNA: inserts are ~22 nt, so reads run
                         through into adapter and will not map to miRBase untrimmed.
    * --genome          : Invokes Quant + specifies reference genome; available options < hsa, mmu, cel > 
    * --instrument      : Use 'nova' for 2 channel chemistry, else use 'hiseq'
    * --adapter         : 3' adapter to trim < default: AGATCGGAAGAGCACACGTCT, NEBNext Small RNA 3' SR Adaptor >
    * --min_len         : Minimum read length kept after trimming < default: 10 >


"""

    exit 0
}


// Channels:

ch_pin      = channel.value(params.id)
ch_genome   = channel.value(params.genome)
ch_sheet    = channel.fromPath(params.sheet)
ch_meta     = ch_sheet
                | splitCsv( header:true )
                | map { row -> [row.label, file(row.fastq1), row.config]}
                | view

// Re-read as a value channel: ch_sheet is consumed by splitCsv above, and
// SMRNA_MQC_TABLES needs the sheet again to map config codes back to labels.
ch_sheet_f  = channel.value(file(params.sheet))
ch_mqc_conf = channel.value(file("$projectDir/assets/multiqc_config.yaml"))
ch_mqc_logo = channel.value(file("$projectDir/img/trex-extended-logo.png"))

// Shown in the MultiQC report header; NA when the run had no --genome.
ch_mqc_genome = channel.value(params.genome ?: 'NA')

// Import Modules:

include { FASTP              } from './modules/fastp'
include { FASTQ2FASTA        } from './modules/fastp'
include { MAPPER             } from './modules/mirdeep2'
include { QUANT              } from './modules/mirdeep2'
include { SMRNA_MQC_TABLES   } from './modules/multiqc'
include { MULTIQC            } from './modules/multiqc'
include { DUMP_VERSIONS      } from './modules/versions'



process WCONFIG {

    tag "writing_config.txt"
    publishDir "mirdeep2", mode: "symlink"        

    input:
    val lines

    output:
    path "config.txt"         , emit:config_file

    script:
    """
    echo -e '${lines.join("\\n")}' >> config.txt
    """
}


workflow {

    if( params.fastp ){

        FASTP(ch_meta)

        ch_fastp_out  = FASTP.out.trimmed_fqs
        ch_fastp_json = FASTP.out.json
        ch_trim_vers  = FASTP.out.versions

    } else {

        log.warn "Running without --fastp: reads are not adapter trimmed. smRNA " +
                 "inserts are shorter than the read, so untrimmed reads carry adapter " +
                 "and will not map to miRBase."

        FASTQ2FASTA(ch_meta)

        ch_fastp_out  = FASTQ2FASTA.out.trimmed_fqs
        ch_fastp_json = channel.empty()
        ch_trim_vers  = channel.empty()

    }

    ch_config = ch_fastp_out
                    // .collect()
                    .map{id, fq, fa, config -> "$fa\t$config" }
                    .flatMap { it.split(/,\s*/) }   // split on comma + optional space
                    .map { it.trim() }
                    .collect()
                    .view()

    WCONFIG(ch_config)

    MAPPER(ch_pin, WCONFIG.out.config_file)

    // Collect software versions. FASTP runs per sample so its channel needs
    // .first(); MAPPER and QUANT run once each and are already value channels.
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(ch_trim_vers.first())
    ch_versions = ch_versions.mix(MAPPER.out.versions)

    // QUANT fills in the miRBaseMatch column; without it MAPPER's table has the
    // flag hard-coded to 0, so the miRBase panel is simply omitted downstream.
    if( params.genome != null ){

        QUANT(ch_pin, ch_genome, MAPPER.out.collapsed_out)

        ch_awk_table = QUANT.out.awk_quant_out
        ch_versions  = ch_versions.mix(QUANT.out.versions)

    } else {

        ch_awk_table = MAPPER.out.awk_out

    }

    DUMP_VERSIONS(ch_versions.collect())

    SMRNA_MQC_TABLES(ch_pin, ch_awk_table, ch_sheet_f, params.genome != null)

    MULTIQC(
        ch_pin,
        ch_fastp_json.collect().ifEmpty([]),
        SMRNA_MQC_TABLES.out.mqc_files.collect(),
        DUMP_VERSIONS.out.mqc_yml,
        ch_mqc_conf,
        ch_mqc_logo,
        ch_mqc_genome
    )
}