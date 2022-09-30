outdir_plots = output_dir + "plots/"
outdir_meanVAFplots = outdir_plots + "meanVAF-Z-score-barplots/"
outdir_BinaryVAFplots = outdir_plots + "BinaryVAF-Z-score-barplots/"

shell("mkdir -p " + outdir_plots)
shell("mkdir -p " + outdir_meanVAFplots)

rule MeanVAF_snv_barplots:
    input:
        glob.glob(output_dir+"meanvaf/{ensemble}/*")
    output:
        csv1 = outdir_meanVAFplots + "snv_{ensemble}.csv",
        pdf1 = outdir_meanVAFplots + "snv_{ensemble}.pdf",
        csv2 = outdir_meanVAFplots + "qual_snv_{ensemble}.csv",
        pdf2 = outdir_meanVAFplots + "qual_snv_{ensemble}.pdf"
    params:
        indir = output_dir+"meanvaf/{ensemble}",
        outdir = outdir_meanVAFplots
    resources:
        mem_mb = 12000
    threads:
        threadsparallel
    shell:
        """
        DIRS=$(ls -d {params.indir}/*)
        conda run --name r Rscript code/meanVAF-Z-scores-barplots.R {params.outdir} $DIRS
        """ 


shell("mkdir -p " + outdir_plots)
shell("mkdir -p " + outdir_meanVAFplots)

rule BinaryVAF_snv_barplots:
    input:
        glob.glob(output_dir+"meanvaf/maxinsertsize150/{ensemble}/*")
    output:
        csv1 = outdir_BinaryVAFplots + "snv_{ensemble}.csv",
        pdf1 = outdir_BinaryVAFplots + "snv_{ensemble}.pdf",
        csv2 = outdir_BinaryVAFplots + "qual_snv_{ensemble}.csv",
        pdf2 = outdir_BinaryVAFplots + "qual_snv_{ensemble}.pdf"
    params:
        indir = output_dir+"meanvaf/maxinsertsize150/{ensemble}",
        outdir = outdir_BinaryVAFplots
    resources:
        mem_mb = 12000
    threads:
        threadsparallel
    shell:
        """
        DIRS=$(ls -d {params.indir}/*)
        conda run --name r Rscript code/BinaryVAF-Z-scores-barplots.R {params.outdir} $DIRS
        """ 