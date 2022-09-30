#########################################################
####                    MRDetect                     ####
#########################################################

script_dir_mrdetectcna=os.getcwd()+"/code/MRDetectCNA/"

MRDetectCNA_outdir=output_dir+"raw/varscan2cna/mrdetectcna/max150/"

rule MRDetectCNA_DepthOfCoverage:
    input:
        plasma=plasma_dir+"{plasma}"+bamsuffix,
        bed=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed"
    output:
        dummy=MRDetectCNA_outdir+"{tumor}_vs_{normal}/{plasma}/MRDetectCNA_DepthOfCoverage.complete"
    params:
        output_dir=MRDetectCNA_outdir+"{tumor}_vs_{normal}/",
        windowsize=500
    resources:
        mem_mb=256000
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate mrdetect \
        && echo $CONDA_PREFIX;

        sample_dir={params.output_dir}/{wildcards.plasma}
        if [ ! -d $sample_dir ]; then
            mkdir -p $sample_dir
        fi

        interval_dir={params.output_dir}{wildcards.plasma}/tmp/intervals
        segments_out_dir={params.output_dir}{wildcards.plasma}/tmp/segments

        if [ ! -d $interval_dir ]; then
            mkdir -p $interval_dir
        fi

        if [ ! -d $segments_out_dir ]; then
            mkdir -p $segments_out_dir
        fi

        # Run DOC for each CNV interval - running in parallel for each interval
        doit() {{
            if [[ ! ${{1:0:1}} == '#' ]]; then
                v1=$(echo $1 | awk -F" " '{{print $1}}')
                v2=$(echo $1 | awk -F" " '{{print $2}}')
                v3=$(echo $1 | awk -F" " '{{print $3}}')
                v4=$(echo $1 | awk -F" " '{{print $4}}')
                DICT={ref}
                DICT=${{DICT%.fasta}}.dict
                cat $DICT > {params.output_dir}{wildcards.plasma}/tmp/intervals/${{v1}}_${{v2}}_${{v3}}.${{v4}}.cnv.interval_list; printf "$v1\t$v2\t$v3\t+\twindow\n" >> {params.output_dir}{wildcards.plasma}/tmp/intervals/${{v1}}_${{v2}}_${{v3}}.${{v4}}.cnv.interval_list
                if [ ! -f "{params.output_dir}{wildcards.plasma}/tmp/segments/${{v1}}_${{v2}}_${{v3}}.${{v4}}.{wildcards.plasma}.DepthOfCoverage.sample_summary" ]; then
                    java \
                    -jar {script_dir_mrdetectcna}gakt3.4/GenomeAnalysisTK.jar \
                    -T DepthOfCoverage \
                    -dels \
                    -rf BadCigar \
                    -ct 100 \
                    -L {params.output_dir}{wildcards.plasma}/tmp/intervals/${{v1}}_${{v2}}_${{v3}}.${{v4}}.cnv.interval_list \
                    --interval_merging OVERLAPPING_ONLY \
                    -R {ref} \
                    -I {input.plasma} \
                    -o {params.output_dir}{wildcards.plasma}/tmp/segments/${{v1}}_${{v2}}_${{v3}}.${{v4}}.{wildcards.plasma}.DepthOfCoverage \
                    &>{params.output_dir}{wildcards.plasma}/tmp/segments/${{v1}}_${{v2}}_${{v3}}.${{v4}}.{wildcards.plasma}.log
                fi
                rm -f {params.output_dir}{wildcards.plasma}/tmp/intervals/${{v1}}_${{v2}}_${{v3}}.${{v4}}.cnv.interval_list
            fi
        }}
        export -f doit
        parallel -a {input.bed} -j {threads} doit

        # Generate dummy output
        echo "success" > {output.dummy}
        """


rule MRDetectCNA_PrintMedian_CNV:
    input:
        dummy=MRDetectCNA_outdir+"{tumor}_vs_{normal}/{plasma}/MRDetectCNA_DepthOfCoverage.complete",
        bed=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed"
    output:
        cnv=MRDetectCNA_outdir+"{tumor}_vs_{normal}/{plasma}/{plasma}.win500.delsdups.perwindow.median"
    params:
        output_dir=MRDetectCNA_outdir+"{tumor}_vs_{normal}/",
        windowsize=500
    resources:
        mem_mb=12000
    threads: 
        1
    shell:
        """
        # Generate median counts for each window - depends on the output of the previous step
        python {script_dir_mrdetectcna}printmedian.py -i {params.output_dir}{wildcards.plasma}/tmp/segments -o {params.output_dir}{wildcards.plasma} -s {wildcards.plasma} -f {input.bed} -w {params.windowsize} -t CNV -z
        """


rule MRDetectCNA_PrintMedian_NEU:
    input:
        dummy=MRDetectCNA_outdir+"{tumor}_vs_{normal}/{plasma}/MRDetectCNA_DepthOfCoverage.complete",
        bed=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed"
    output:
        MRDetectCNA_outdir+"{tumor}_vs_{normal}/{plasma}/{plasma}.win500.neu.perwindow.median"
    params:
        output_dir=MRDetectCNA_outdir+"{tumor}_vs_{normal}/",
        windowsize=500
    resources:
        mem_mb=12000
    threads: 
        1
    shell:
        """
        # Generate median counts for each window - depends on the output of the previous step
        python {script_dir_mrdetectcna}printmedian.py -i {params.output_dir}{wildcards.plasma}/tmp/segments -o {params.output_dir}{wildcards.plasma} -s {wildcards.plasma} -f {input.bed} -w {params.windowsize} -t NEU -z
        """


# rule MRDetectCNA_CreatePairs:


#rule MRDetectCNA_PairRunner:
#    input:
#        plasma=ALL_PLASMA,
#        bed=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.copyCalled.recenter.sep.bed"
#    output:
#        end=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/MRDetectCNA_SampleRunner_complete"
#    params:
#        output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/mrdetect/"
#    threads:
#        N_ALL_PLASMA
#    shell:
#        """
#        set +eu \
#        && PS1=dummy \
#        && . $(conda info --base)/etc/profile.d/conda.sh \
#        && conda activate mrdetect \
#        && echo $CONDA_PREFIX;
#
#        mkdir -p {params}
#
#        for plasma in {input.plasma}
#        do
#            tmp=${{plasma##*/}}
#            mkdir -p {params}${{tmp%.*}}/
#            scripts/MRDetect-master/MRDetectCNA/MRDetectCNA_sample_runner.sh $plasma {ref} 10000 {input.bed} {params} ${{tmp%.*}} all &
#        done
#        echo "SUCCESS" > {output.end}
#        """
