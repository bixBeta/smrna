# smrna

Nextflow pipeline for smRNA seq runs.

## Usage

```
s  m  R  N  A  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/smrna -r main < args ... >

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
```

Outputs land in `multiqc/`, with the report as `<id>_multiqc_report.html` and the
generated custom-content sections under `multiqc/custom_content/`.
