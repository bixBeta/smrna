process MAPPER{

    maxForks 1
    tag "$pin"
    label 'process_mirdeep2'
    
    publishDir "mirdeep2", mode: "symlink", overwrite: true


    input:
        val(pin)
        path(config)
    
    output:
        path("mirmap_firstbase_readlengthcounts.txt")     , emit: awk_out
        path("${pin}.collapsed.fa.gz")                    , emit: collapsed_out
             


    script:


    """

        mapper.pl $config -d -c -m -s ${pin}.collapsed.fa

        gzip ${pin}.collapsed.fa

        awk 'BEGIN {OFS="\t"; m=0} {if(FNR%2==1) {s=substr(\$1,2,3); c=substr(\$1,match(\$1,"_x")+2,15);} else { l=length(\$1); f=substr(\$1,1,1); a[s,l,f,m]++; b[s,l,f,m]=b[s,l,f,m]+c;s=0;c=0;l=0;f=0; m=0}} END {print "library\treadlength\tbase1\tmiRBaseMatch\t#distinctReads\t#reads"; for (var in a) {split(var,q,SUBSEP); print q[1], q[2], q[3], q[4], a[var], b[var]} }' <(zcat ${pin}.collapsed.fa.gz) > mirmap_firstbase_readlengthcounts.txt 
        

    """
}


process QUANT{

    maxForks 1
    tag "$genome"
    label 'process_mirdeep2'
    
    publishDir "mirdeep2/expression_analyses/", mode: "symlink", overwrite: true


    input:
        val(pin)
        val(genome)
        path(collapsed)
    
    output:
        path("quant_mirmap_firstbase_readlengthcounts.txt")     , emit: awk_quant_out
        path("*.arf.gz")                                        , emit: arf_out
        path("*.mrd.gz")                                        , emit: mrd_out
        path("*.csv")                                           , emit: counts_out
        path("*.html")                                          , emit: html_out


    script:


    """
        gunzip -f ${collapsed}

        quantifier.pl -p /workdir/genomes/smRNA/hairpin.fa \\
            -m /workdir/genomes/smRNA/mature.fa \\
            -t ${genome} -y ${pin} -r ${pin}.collapsed.fa -W -d

        
        cd expression_analyses/expression_analyses_${pin}
        mv miRBase.mrd ${pin}_miRBase.mrd
        mv *.mrd *arf ../../
        cd ../../
        
        gzip *.arf *.fa *.mrd

        awk 'BEGIN {OFS="\t"; m=0} \\
            FNR==NR {d[\$1]=1; next} { if(FNR%2==1) {s=substr(\$1,2,3); \\
            c=substr(\$1,match(\$1,"_x")+2,15); \\
            if(substr(\$1,2,25) in d) m=1} \\
            else { l=length(\$1); f=substr(\$1,1,1); \\
            a[s,l,f,m]++; b[s,l,f,m]=b[s,l,f,m]+c;s=0;c=0;l=0;f=0; m=0}} \\
            END {print "library\treadlength\tbase1\tmiRBaseMatch\t#distinctReads\t#reads"; \\
            for (var in a) {split(var,q,SUBSEP); print q[1], q[2], q[3], q[4], a[var], b[var]} }' <(zcat *miRBase.mrd.gz) <(zcat *collapsed.fa.gz) > quant_mirmap_firstbase_readlengthcounts.txt

    """
}