outdir_plots = output_dir + "plots/"
outdir_mrdetectplots = outdir_plots + "MRDetect-Z-score-barplots/"

shell("mkdir -p " + outdir_plots)
shell("mkdir -p " + outdir_mrdetectplots)

def input_snv_barplots(wildcards):
    return sorted(glob.glob(output_dir+"variants_final/" + wildcards.ensemble + "/*/"))

rule snv_barplots:
    input:
        input_snv_barplots
    output:
        csv = outdir_mrdetectplots + "snv_{ensemble}.csv",
        pdf = outdir_mrdetectplots + "snv_{ensemble}.pdf"
    params:
        outdir_mrdetectplots
    resources:
        mem_mb = 12000
    threads:
        threadsparallel
    shell:
        """
        conda run --name r Rscript code/MRDetectSNV-Z-scores-barplots.R {params} {input}
        """ 
        