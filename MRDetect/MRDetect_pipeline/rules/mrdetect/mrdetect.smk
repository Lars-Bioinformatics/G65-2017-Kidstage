#########################################################
####                    MRDetect                     ####
#########################################################

#ALL_PLASMA=[glob.glob("/mnt/d/master/data/plasmanotfiltered/*.bam")]
#ALL_PLASMA=[glob.glob("/mnt/d/master/data/plasmafilteredmaxinsertsize/*.bam")]
# ALL_PLASMA=[glob.glob("/mnt/d/master/data/plasmafilteredminmaxinsertsize/*.bam")]
# N_ALL_PLASMA=len(ALL_PLASMA)
# print(ALL_PLASMA)

rule copyConsensusVariants:
    input:
        vcf=output_dir+"variants_ensemble/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf",
    output:
        vcf=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf",
    params:
        out_dir=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}"
    shell:
        """
        mkdir -p {params.out_dir}
        cp {input} {output}
        cp scripts/MRDetect-master/MRDetectSNV/* {params}
        """

rule MRDetectSNV:
    input:
        vcf=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        plasma = input_dir+"{plasma}"+bamsuffix
    output:
        # end=output_dir+"variants_final/{ensemble}/{tumor}_vs_{plasma}_vs_{normal}/MRDetectSNV_complete"
        plasma = output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/{plasma}"+bamsuffix+"_VS_Consensus.sSNV.vcf_RESULT.csv"
    params:
        output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/"
    shell:
        """
        cd {params}
        conda run --name mrdetect \
        bash MRDetectSNV_runner.sh {input.plasma} {input.vcf}
        """

rule MRDetectSNV_plasmaControls:
    input:
        plasma = expand(output_dir+"variants_final/{{ensemble}}/{{tumor}}_vs_{{normal}}/{plasma}"+bamsuffix+"_VS_Consensus.sSNV.vcf_RESULT.csv", plasma=ALL_PLASMA)
    output:
        end=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/MRDetectSNV_complete"
        # end=output_dir+"variants_final/{ensemble}/{tumor}_vs_{plasma}_vs_{normal}/MRDetectSNV_complete"
    shell:
        """
        echo "SUCCESS" > {output.end}
        """

# checkpoint MRDetectSNV:
#     input:
#         plasma=ALL_PLASMA,
#         vcf=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf"
#     output:
#         end=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/MRDetectSNV_complete"
#     params:
#         output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/"
#     threads:
#         N_ALL_PLASMA
#     shell:
#         """
#         conda run --name mrdetect \
#         DIR=$PWD
#         cd {params}
#         for plasma in {input.plasma}
#         do
#             scripts/MRDetect-master/MRDetectSNV/MRDetectSNV_runner.sh $plasma {input.vcf} &
#         done
#         echo "SUCCESS" > {output.end}
#         cd $DIR
#         """

# checkpoint MRDetectSNV:
#     input:
#         plasma=ALL_PLASMA,
#         vcf=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/Consensus.sSNV.vcf"
#     output:
#         end=output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/MRDetectSNV_complete"
#     params:
#         output_dir+"variants_final/{ensemble}/{tumor}_vs_{normal}/"
#     threads:
#         N_ALL_PLASMA
#     shell:
#         """
#         set +eu \
#         && PS1=dummy \
#         && . $(conda info --base)/etc/profile.d/conda.sh \
#         && conda activate mrdetect \
#         && echo $CONDA_PREFIX;
#         cp scripts/MRDetect-master/MRDetectSNV/* {params}
#         DIR=$PWD
#         cd {params}
#         for plasma in {input.plasma}
#         do
#             {params}MRDetectSNV_runner.sh $plasma {input.vcf} &
#         done
#         echo "SUCCESS" > {output.end}
#         cd $DIR
#         """

#rule MRDetectCNA_SampleRunner:
#    input:
#        plasma=ALL_PLASMA,
#        bed=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/{tumor}_vs_{normal}.cnv.bed"
#    output:
#        end=output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/maxinsertsize150/MRDetectCNA_SampleRunner_complete"
#    params:
#        output_dir+"raw/varscan2cna/{tumor}_vs_{normal}/maxinsertsize150/mrdetect/"
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
#            scripts/MRDetect-master/MRDetectCNA/MRDetectCNA_sample_runner.sh $plasma {ref} 500 {input.bed} {params} ${{tmp%.*}} all &
#        done
#        echo "SUCCESS" > {output.end}
#        """



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