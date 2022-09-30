rule filter_variants_triallelic:
    input:
        output_variants_notfiltered+"{tumor}_vs_{normal}/{variantcalleroutput}"
    output:
        output_variants_filtered+"{tumor}_vs_{normal}/{variantcalleroutput}"
    resources:
        mem_mb=12000
    threads:
        threadsparallel
    shell:
        """
        conda run -n bcftools bcftools view -m2 -M2 -v snps --apply-filters PASS --threads {threads} -Ov -o {output} {input}
        """