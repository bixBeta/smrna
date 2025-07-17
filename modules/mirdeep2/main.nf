process MAPPER {

    maxForks 8
    tag "$id"
    label 'process_high'
    
    publishDir "mirdeep2", mode: "symlink", overwrite: true


    input:
        val(pin)
        path(config)
    
    output:
        path(*)     , emit: mapper_outs
             


    script:


    """

        mapper.pl $config -d -c -m -s ${pin}.collapsed.fa


    """


}