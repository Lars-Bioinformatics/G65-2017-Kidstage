shell("mkdir -p "+output_dir+"raw/varscan2cna/")

rule VarScan2_CNA_flagstat_normal:
    input:
        # tumor=input_dir+"{tumor}"+bamsuffix,
        normal=input_dir+"{normal}"+bamsuffix
    output:
        # flagstat_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.flagstat",
        flagstat_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.flagstat"
    params:
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem_mb=20000,
        cpus=1
    shell:
        """
        mkdir -p {params.outdir}

        # Samtools flagstat
        conda run --name varscan2cna \
        samtools flagstat {input.normal} > {output.flagstat_normal}
        """
        # samtools flagstat {input.tumor} > {output.flagstat_tumor}

rule VarScan2_CNA_flagstat_tumor:
    input:
        tumor=input_dir+"{tumor}"+bamsuffix,
        # normal=input_dir+"{normal}"+bamsuffix
    output:
        flagstat_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.flagstat",
        # flagstat_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.flagstat"
    params:
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem_mb=20000,
        cpus=1
    shell:
        """
        mkdir -p {params.outdir}

        # Samtools flagstat
        conda run --name varscan2cna \
        samtools flagstat {input.tumor} > {output.flagstat_tumor}
        """
        # samtools flagstat {input.normal} > {output.flagstat_normal}


rule VarScan2_CNA_mpileup_normal:
    input:
        # tumor=input_dir+"{tumor}"+bamsuffix,
        normal=input_dir+"{normal}"+bamsuffix
    output:
        # pileup_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.mpileup",
        pileup_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.mpileup"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem_mb=20000,
        cpus=2
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;

        # Samtools mpileup
        samtools mpileup \
        -R \
        -B \
        -q 10 \
        -Q 30 \
        -f {ref} \
        {input.normal} > {output.pileup_normal}
        """

rule VarScan2_CNA_mpileup_tumor:
    input:
        tumor=input_dir+"{tumor}"+bamsuffix,
        # normal=input_dir+"{normal}"+bamsuffix
    output:
        pileup_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.mpileup",
        # pileup_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.mpileup"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem_mb=20000,
        cpus=2
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;

        # Samtools mpileup
        samtools mpileup \
        -R \
        -B \
        -q 10 \
        -Q 30 \
        -f {ref} \
        {input.tumor} > {output.pileup_tumor}
        """

rule VarScan2_CNA_copynumber:
    input:
        flagstat_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.flagstat",
        flagstat_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.flagstat",
        pileup_tumor=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}.mpileup",
        pileup_normal=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{normal}.mpileup"
    output:
        copynumber=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copynumber"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/",
        copynumber=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # Varscan copynumber
        # must calculate data ratio from flagstat output
        # also must move to output dir to run this because varscan doesn't parse the output name
        tnum=$(grep -m 1 mapped {input.flagstat_tumor} | cut -f1 -d' ')
        cnum=$(grep -m 1 mapped {input.flagstat_normal} | cut -f1 -d' ')
        dratio=$(echo "scale=2;$cnum/$tnum" | bc)

        varscan copynumber {input.pileup_normal} \
        {input.pileup_tumor} \
        {params.copynumber} \
        --data-ratio $dratio
        """

rule VarScan2_CNA_copycaller_filter:
    input:
        copynumber=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copynumber"
    output:
        copynumber_cov=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copynumber.cov"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # From the output, filter any segments for which the tumor coverage is less than 20
        # and the control coverage is less than 20
        awk -v x=20 '$6 >= x' {input.copynumber} | \
        awk -v x=20 '$5 >= x' > {output.copynumber_cov}
        """

rule VarScan2_CNA_copyCaller_first:
    input:
        copynumber_cov=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copynumber.cov"
    output:
        copycalled1=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled",
        copycalled1_homdel=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled_homdel"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # Varscan copycaller
        varscan copyCaller {input.copynumber_cov} --output-file {output.copycalled1} --output-homdel-file {output.copycalled1_homdel}
        """

rule VarScan2_CNA_copyCaller_second:
    input:
        copynumber_cov=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copynumber.cov",
        copycalled1=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled",
        copycalled1_homdel=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled_homdel"
    output:
        copycalled2=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter",
        copycalled2_homdel=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled_homdel.recenter"
    params:
        env="varscan2cna",
        scriptdir="scripts/varscan2-master/",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # Calculate recenter amount
        delta=$({params.scriptdir}meanLogRatioByChromosome.py {input.copycalled1})

        # Rerun copycaller
        cmp=$(awk -v delta=$delta 'END{{if (delta < -0.2) {{print "lt"}} else {{if (delta > 0.2) {{print "gt"}} else {{print "eq"}}}}}}' < /dev/null)
        if [[ "$cmp" == "lt" ]]; then
            rd=$(echo $delta | sed 's/-//')
            varscan copyCaller {input.copynumber_cov} --output-file {output.copycalled2} --output-homdel-file {output.copycalled2_homdel} --recenter-down $rd
        elif [[ "$cmp" == "gt" ]]; then
            varscan copyCaller {input.copynumber_cov} --output-file {output.copycalled2} --output-homdel-file {output.copycalled2_homdel} --recenter-up $delta
        else
            cp {input.copycalled1} {output.copycalled2}
            cp {input.copycalled1_homdel} {output.copycalled2_homdel}
        fi
        """

rule VarScan2_CNA_separateArms:
    input:
        copycalled2=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter",
        copycalled2_homdel=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled_homdel.recenter"
    output:
        separatearms=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter.sep"
    params:
        env="varscan2cna",
        scriptdir="scripts/varscan2-master/",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # add p and q to chromosome arms
        {params.scriptdir}separateArms.py {input.copycalled2} {centromeres} > {output.separatearms}
        """

rule VarScan2_CNA_DNAcopy:
    input:
        separatearms=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter.sep"
    output:
        cirbinseg=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter.sep.SD.2.5.dnacopy.out"
    params:
        env="varscan2cna",
        scriptdir="scripts/varscan2-master/",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # Circular binary segmentation
        Rscript {params.scriptdir}basicDNAcopy.R {input.separatearms} 2.5 >/dev/null
        """

rule VarScan2_CNA_final:
    input:
        cirbinseg=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter.sep.SD.2.5.dnacopy.out"
    output:
        final=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv"
    params:
        env="varscan2cna",
        outdir=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/"
    resources:
        mem=240000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        
        # remove the arms and print to stdout
        sed 's/\.[pq]	/	/' {input.cirbinseg} | \
            sed "s/^sample/{wildcards.tumor}_vs_{wildcards.normal}/" > {output.final}
        """

rule VarScan2_CNA_convBED:
    input:
        separatearms=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv"
    output:
        separatearms=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed"
    resources:
        mem=12000
    threads:
        1
    shell:
        """
        awk 'OFS=\"\t\" {{
            if(NR!=1) {{
                if($6 > 0.25) {{
                    print $2,$3,$4,"DUP",$6
                }} else if($6 < -0.25) {{
                    print $2,$3,$4,"DEL",$6
                }} else {{
                    print $2,$3,$4,"NEU",$6
                }}
            }} else {{
                print "#chr","start","end","type","log2"
            }}
        }}' {input} > {output}
        """