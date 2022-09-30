#########################################################
####                   SomaticSeq                    ####
#########################################################

# Initial operations
ensemble_output_dir = output_dir+"variants_ensemble/"
somaticseq_output_dir = ensemble_output_dir+"somaticseq/"
somaticseqmincallers1_output_dir = ensemble_output_dir+"somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[0])+"/"
somaticseqmincallers2_output_dir = ensemble_output_dir+"somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[1])+"/"
somaticseqmincallers3_output_dir = ensemble_output_dir+"somaticseq_mincallers"+str(SOMATICSEQ_MINCALLERS[2])+"/"
mrdetectensemble_output_dir = ensemble_output_dir+"mrdetect/"
mutect2_output_dir = ensemble_output_dir+"mutect2/"
varscan2_output_dir = ensemble_output_dir+"varscan2/"


onstart:
    shell("mkdir -p " + ensemble_output_dir)
    shell("mkdir -p " + somaticseq_output_dir)
    shell("mkdir -p " + somaticseq_output_dir)
    shell("mkdir -p " + somaticseq_output_dir)
    shell("mkdir -p " + mrdetectensemble_output_dir)
    shell("mkdir -p " + mutect2_output_dir)
    shell("mkdir -p " + varscan2_output_dir)

rule SomaticSeq:
    input:
        normal = input_dir + "{normal}" + bamsuffix,
        tumor = input_dir + "{tumor}" + bamsuffix,
        lofreq_snvs = output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else [],
        lofreq_indels = output_variants_filtered+"{tumor}_vs_{normal}/lofreq_indels.vcf" if config["variantcaller"]["lofreq"]["include"] else [],
        muse_snvs = output_variants_filtered+"{tumor}_vs_{normal}/muse_snvs.vcf" if config["variantcaller"]["muse"]["include"] else [],
        mutect2_vcf = output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else [],
        scalpel_indels = output_variants_filtered+"{tumor}_vs_{normal}/scalpel_indels.vcf" if config["variantcaller"]["scalpel"]["include"] else [],
        somaticsniper_snvs = output_variants_filtered+"{tumor}_vs_{normal}/somaticsniper_snvs.vcf" if config["variantcaller"]["somaticsniper"]["include"] else [],
        strelka_snvs = output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else [],
        strelka_indels = output_variants_filtered+"{tumor}_vs_{normal}/strelka_indels.vcf" if config["variantcaller"]["strelka"]["include"] else [],
        vardict_snvs = output_variants_filtered+"{tumor}_vs_{normal}/vardict_bothcombined.vcf" if config["variantcaller"]["vardict"]["include"] else [],
        varscan_snvs = output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else [],
        varscan_indels = output_variants_filtered+"{tumor}_vs_{normal}/varscan_indels.vcf" if config["variantcaller"]["varscan"]["include"] else [],
    output:
        cons_snvs = somaticseq_output_dir + "{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        cons_indels = somaticseq_output_dir + "{tumor}_vs_{normal}/Consensus.sINDEL.vcf",
        ense_snvs = somaticseq_output_dir + "{tumor}_vs_{normal}/Ensemble.sSNV.tsv",
        ense_indels = somaticseq_output_dir + "{tumor}_vs_{normal}/Ensemble.sINDEL.tsv"
    params:
        outdir = somaticseq_output_dir + "{tumor}_vs_{normal}/",
        env="somaticseq",
        #inclusion = "-include {medexome}" if config["ngsdatatype"] == "exome" else "",
        lofreq_snvs = "-lofreqsnv "+output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else "",
        lofreq_indels = "-lofreqindel "+output_variants_filtered+"{tumor}_vs_{normal}/lofreq_indels.vcf" if config["variantcaller"]["lofreq"]["include"] else "",
        muse_snvs = "-muse "+output_variants_filtered+"{tumor}_vs_{normal}/muse_snvs.vcf" if config["variantcaller"]["muse"]["include"] else "",
        mutect2_vcf = "-mutect2 "+output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else "",
        scalpel_indels = "-scalpel "+output_variants_filtered+"{tumor}_vs_{normal}/scalpel_indels.vcf" if config["variantcaller"]["scalpel"]["include"] else "",
        somaticsniper_snvs = "-sniper "+output_variants_filtered+"{tumor}_vs_{normal}/somaticsniper_snvs.vcf" if config["variantcaller"]["somaticsniper"]["include"] else "",
        strelka_snvs = "-strelkasnv "+output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else "",
        strelka_indels = "-strelkaindel "+output_variants_filtered+"{tumor}_vs_{normal}/strelka_indels.vcf" if config["variantcaller"]["strelka"]["include"] else "",
        vardict_snvs = "-vardict "+output_variants_filtered+"{tumor}_vs_{normal}/vardict_bothcombined.vcf" if config["variantcaller"]["vardict"]["include"] else "",
        varscan_snvs = "-varscansnv "+output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else "",
        varscan_indels = "-varscanindel "+output_variants_filtered+"{tumor}_vs_{normal}/varscan_indels.vcf" if config["variantcaller"]["varscan"]["include"] else "",
    resources:
        mem_mb = 350000,
        cpus = threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        mkdir -p {params.outdir}
        somaticseq_parallel.py \
        -outdir {params.outdir} \
        -ref {ref} \
        -include {medexome} \
        -minMQ 10 \
        -minBQ 20 \
        -dbsnp {dbsnp} \
        -nt {threads} \
        paired \
        -tbam {input.tumor} \
        -nbam {input.normal} \
        {params.lofreq_snvs} \
        {params.lofreq_indels} \
        {params.muse_snvs} \
        {params.mutect2_vcf} \
        {params.scalpel_indels} \
        {params.somaticsniper_snvs} \
        {params.strelka_snvs} \
        {params.strelka_indels} \
        {params.vardict_snvs} \
        {params.varscan_snvs} \
        {params.varscan_indels}
        """

rule SomaticSeq_mincallers:
    input:
        normal = input_dir + "{normal}" + bamsuffix,
        tumor = input_dir + "{tumor}" + bamsuffix,
        lofreq_snvs = output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else [],
        lofreq_indels = output_variants_filtered+"{tumor}_vs_{normal}/lofreq_indels.vcf" if config["variantcaller"]["lofreq"]["include"] else [],
        muse_snvs = output_variants_filtered+"{tumor}_vs_{normal}/muse_snvs.vcf" if config["variantcaller"]["muse"]["include"] else [],
        mutect2_vcf = output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else [],
        scalpel_indels = output_variants_filtered+"{tumor}_vs_{normal}/scalpel_indels.vcf" if config["variantcaller"]["scalpel"]["include"] else [],
        somaticsniper_snvs = output_variants_filtered+"{tumor}_vs_{normal}/somaticsniper_snvs.vcf" if config["variantcaller"]["somaticsniper"]["include"] else [],
        strelka_snvs = output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else [],
        strelka_indels = output_variants_filtered+"{tumor}_vs_{normal}/strelka_indels.vcf" if config["variantcaller"]["strelka"]["include"] else [],
        vardict_snvs = output_variants_filtered+"{tumor}_vs_{normal}/vardict_bothcombined.vcf" if config["variantcaller"]["vardict"]["include"] else [],
        varscan_snvs = output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else [],
        varscan_indels = output_variants_filtered+"{tumor}_vs_{normal}/varscan_indels.vcf" if config["variantcaller"]["varscan"]["include"] else [],
    output:
        cons_snvs = ensemble_output_dir+"somaticseq_mincallers{SomaticSeq_mincallers}/{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        cons_indels = ensemble_output_dir+"somaticseq_mincallers{SomaticSeq_mincallers}/{tumor}_vs_{normal}/Consensus.sINDEL.vcf",
        ense_snvs = ensemble_output_dir+"somaticseq_mincallers{SomaticSeq_mincallers}/{tumor}_vs_{normal}/Ensemble.sSNV.tsv",
        ense_indels = ensemble_output_dir+"somaticseq_mincallers{SomaticSeq_mincallers}/{tumor}_vs_{normal}/Ensemble.sINDEL.tsv"
    params:
        outdir = ensemble_output_dir+"somaticseq_mincallers{SomaticSeq_mincallers}/{tumor}_vs_{normal}/",
        env="somaticseq",
        #inclusion = "-include {medexome}" if config["ngsdatatype"] == "exome" else "",
        lofreq_snvs = "-lofreqsnv "+output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else "",
        lofreq_indels = "-lofreqindel "+output_variants_filtered+"{tumor}_vs_{normal}/lofreq_indels.vcf" if config["variantcaller"]["lofreq"]["include"] else "",
        muse_snvs = "-muse "+output_variants_filtered+"{tumor}_vs_{normal}/muse_snvs.vcf" if config["variantcaller"]["muse"]["include"] else "",
        mutect2_vcf = "-mutect2 "+output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else "",
        scalpel_indels = "-scalpel "+output_variants_filtered+"{tumor}_vs_{normal}/scalpel_indels.vcf" if config["variantcaller"]["scalpel"]["include"] else "",
        somaticsniper_snvs = "-sniper "+output_variants_filtered+"{tumor}_vs_{normal}/somaticsniper_snvs.vcf" if config["variantcaller"]["somaticsniper"]["include"] else "",
        strelka_snvs = "-strelkasnv "+output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else "",
        strelka_indels = "-strelkaindel "+output_variants_filtered+"{tumor}_vs_{normal}/strelka_indels.vcf" if config["variantcaller"]["strelka"]["include"] else "",
        vardict_snvs = "-vardict "+output_variants_filtered+"{tumor}_vs_{normal}/vardict_bothcombined.vcf" if config["variantcaller"]["vardict"]["include"] else "",
        varscan_snvs = "-varscansnv "+output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else "",
        varscan_indels = "-varscanindel "+output_variants_filtered+"{tumor}_vs_{normal}/varscan_indels.vcf" if config["variantcaller"]["varscan"]["include"] else "",
    resources:
        mem_mb = 350000,
        cpus = threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        mkdir -p {params.outdir}
        somaticseq_parallel.py \
        -outdir {params.outdir} \
        -ref {ref} \
        -include {medexome} \
        -mincaller {wildcards.SomaticSeq_mincallers} \
        -minMQ 10 \
        -minBQ 20 \
        -dbsnp {dbsnp} \
        -nt {threads} \
        paired \
        -tbam {input.tumor} \
        -nbam {input.normal} \
        {params.lofreq_snvs} \
        {params.lofreq_indels} \
        {params.muse_snvs} \
        {params.mutect2_vcf} \
        {params.scalpel_indels} \
        {params.somaticsniper_snvs} \
        {params.strelka_snvs} \
        {params.strelka_indels} \
        {params.vardict_snvs} \
        {params.varscan_snvs} \
        {params.varscan_indels}
        """

MRDetectENSEMBLE_mincallers=0
for variantcaller in VARIANTCALLERS:
    if variantcaller == "lofreq":
        MRDetectENSEMBLE_mincallers=MRDetectENSEMBLE_mincallers+1
    elif variantcaller == "mutect":
        MRDetectENSEMBLE_mincallers=MRDetectENSEMBLE_mincallers+1
    elif variantcaller == "strelka":
        MRDetectENSEMBLE_mincallers=MRDetectENSEMBLE_mincallers+1

rule MRDetectENSEMBLE:
    input:
        normal = input_dir + "{normal}" + bamsuffix,
        tumor = input_dir + "{tumor}" + bamsuffix,
        lofreq_snvs = output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else [],
        mutect2_vcf = output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else [],
        strelka_snvs = output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else []
    output:
        cons_snvs = mrdetectensemble_output_dir + "{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        cons_indels = mrdetectensemble_output_dir + "{tumor}_vs_{normal}/Consensus.sINDEL.vcf",
        ense_snvs = mrdetectensemble_output_dir + "{tumor}_vs_{normal}/Ensemble.sSNV.tsv",
        ense_indels = mrdetectensemble_output_dir + "{tumor}_vs_{normal}/Ensemble.sINDEL.tsv"
    params:
        outdir = mrdetectensemble_output_dir + "{tumor}_vs_{normal}/",
        env="somaticseq",
        #inclusion = "-include {medexome}" if config["ngsdatatype"] == "exome" else "",
        lofreq_snvs = "-lofreqsnv "+output_variants_filtered+"{tumor}_vs_{normal}/lofreq_snvs.vcf" if config["variantcaller"]["lofreq"]["include"] else "",
        mutect2_vcf = "-mutect2 "+output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else "",
        strelka_snvs = "-strelkasnv "+output_variants_filtered+"{tumor}_vs_{normal}/strelka_snvs.vcf" if config["variantcaller"]["strelka"]["include"] else ""
    resources:
        mem_mb = 350000,
        cpus = threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        mkdir -p {params.outdir}
        somaticseq_parallel.py \
        -outdir {params.outdir} \
        -ref {ref} \
        -include {medexome} \
        -mincaller {MRDetectENSEMBLE_mincallers} \
        -minMQ 10 \
        -minBQ 20 \
        -dbsnp {dbsnp} \
        -nt {threads} \
        paired \
        -tbam {input.tumor} \
        -nbam {input.normal} \
        {params.lofreq_snvs} \
        {params.mutect2_vcf} \
        {params.strelka_snvs} \
        """


rule Mutect2_final:
    input:
        normal = input_dir + "{normal}" + bamsuffix,
        tumor = input_dir + "{tumor}" + bamsuffix,
        mutect2_vcf = output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else [],
    output:
        cons_snvs = mutect2_output_dir + "{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        cons_indels = mutect2_output_dir + "{tumor}_vs_{normal}/Consensus.sINDEL.vcf",
        ense_snvs = mutect2_output_dir + "{tumor}_vs_{normal}/Ensemble.sSNV.tsv",
        ense_indels = mutect2_output_dir + "{tumor}_vs_{normal}/Ensemble.sINDEL.tsv"
    params:
        outdir = mutect2_output_dir + "{tumor}_vs_{normal}/",
        env="somaticseq",
        #inclusion = "-include {medexome}" if config["ngsdatatype"] == "exome" else "",
        mutect2_vcf = "-mutect2 "+output_variants_filtered+"{tumor}_vs_{normal}/mutect_bothcombined.vcf" if config["variantcaller"]["mutect"]["include"] else ""
    resources:
        mem_mb = 350000,
        cpus = threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        mkdir -p {params.outdir}
        somaticseq_parallel.py \
        -outdir {params.outdir} \
        -ref {ref} \
        -include {medexome} \
        -mincaller 1 \
        -minMQ 10 \
        -minBQ 20 \
        -dbsnp {dbsnp} \
        -nt {threads} \
        paired \
        -tbam {input.tumor} \
        -nbam {input.normal} \
        {params.mutect2_vcf}
        """


rule VarScan2_final:
    input:
        normal = input_dir + "{normal}" + bamsuffix,
        tumor = input_dir + "{tumor}" + bamsuffix,
        varscan_snvs = output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else [],
    output:
        cons_snvs = varscan2_output_dir + "{tumor}_vs_{normal}/Consensus.sSNV.vcf",
        cons_indels = varscan2_output_dir + "{tumor}_vs_{normal}/Consensus.sINDEL.vcf",
        ense_snvs = varscan2_output_dir + "{tumor}_vs_{normal}/Ensemble.sSNV.tsv",
        ense_indels = varscan2_output_dir + "{tumor}_vs_{normal}/Ensemble.sINDEL.tsv"
    params:
        outdir = varscan2_output_dir + "{tumor}_vs_{normal}/",
        env="somaticseq",
        #inclusion = "-include {medexome}" if config["ngsdatatype"] == "exome" else "",
        varscan_snvs = "-varscansnv "+output_variants_filtered+"{tumor}_vs_{normal}/varscan_snvs.vcf" if config["variantcaller"]["varscan"]["include"] else ""
    resources:
        mem_mb = 350000,
        cpus = threadsparallel
    threads:
        threadsparallel
    shell:
        """
        set +eu \
        && PS1=dummy \
        && . $(conda info --base)/etc/profile.d/conda.sh \
        && conda activate {params.env} \
        && echo $CONDA_PREFIX;
        mkdir -p {params.outdir}
        somaticseq_parallel.py \
        -outdir {params.outdir} \
        -ref {ref} \
        -include {medexome} \
        -mincaller 1 \
        -minMQ 10 \
        -minBQ 20 \
        -dbsnp {dbsnp} \
        -nt {threads} \
        paired \
        -tbam {input.tumor} \
        -nbam {input.normal} \
        {params.varscan_snvs} 
        """