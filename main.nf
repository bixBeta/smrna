nextflow.enable.dsl=2


// Params:

params.help                     = false

params.sheet                    = "sample-sheet.csv"
params.id                       = "TREX_ID"
params.genome                   = null
params.instrument               = "nova"
params.genome                   = null

if( params.help ) {

log.info """
R  N  A  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/smrna -r g2 < args ... >

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
    * --fastp           : Invokes fastp trimming module.
    * --genome          : Invokes Quant + specifies reference genome; available options < hsa, mmu, cel > 
    * --instrument      : Use 'nova' for 2 channel chemistry, else use 'hiseq'  


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

ch_genome   = channel.value(params.genome)

// Import Modules:

include { FASTP              } from './modules/fastp'
include { MAPPER             } from './modules/mirdeep2'
include { QUANT              } from './modules/mirdeep2'



process writeChannelToFile {

    publishDir "trimmed_fastqs", mode: "symlink"        

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

    FASTP(ch_meta)
        .set {ch_fastp_out}
        

    // FASTQ2FASTA(ch_fastp_out)
    //     .set {ch_fasta_config}

    ch_fastp_out 
        // | view   
        | flatten
        // | view


    ch_config = ch_fastp_out
                    // .collect()
                    .map{id, fq, fa, config -> "$fa\t$config" }
                    .flatMap { it.split(/,\s*/) }   // split on comma + optional space
                    .map { it.trim() }
                    .collect()
                    .view()

    writeChannelToFile(ch_config)

    MAPPER(ch_pin, writeChannelToFile.out.config_file)

    if( params.genome != null ){

        QUANT(ch_pin, ch_genome, MAPPER.out.collapsed_out)

    } 
}