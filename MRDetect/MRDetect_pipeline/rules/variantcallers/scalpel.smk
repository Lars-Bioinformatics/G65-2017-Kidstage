#########################################################
####                     Scalpel                     ####
#########################################################

# Initial operations
scalpel_output_dir = output_raw+"scalpel/"
onstart:
    shell("mkdir -p "+scalpel_output_dir)

rule Scalpel_discovery:
    input:
        normal=input_dir+"{normal}"+bamsuffix,
        tumor=input_dir+"{tumor}"+bamsuffix
    output:
        db=scalpel_output_dir+"{tumor}_vs_{normal}/twopass/somatic.db.pag",
        vcf=scalpel_output_dir+"{tumor}_vs_{normal}/somatic.indel.vcf"
    params:
        outdir=scalpel_output_dir+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=240000,
        cpus=1
    threads: 
        1
    shell:
        """
        mkdir {params.outdir} -p
        conda run --name scalpel \
        scalpel-discovery \
        --somatic \
        --normal {input.normal} \
        --tumor {input.tumor} \
        --bed {medexome} \
        --ref {ref} \
        --dir {params.outdir} \
        --numprocs {threads} \
        --two-pass
        """

rule Scalpel_export:
    input:
        db=scalpel_output_dir+"{tumor}_vs_{normal}/twopass/somatic.db.pag"
    output:
        vcf=scalpel_output_dir+"{tumor}_vs_{normal}/somatic_filtered.indel.vcf"
    resources:
        mem_mb=100000,
        cpus=1
    threads: 
        1
    shell:
        """
        conda run --name scalpel \
        scalpel-export \
        --somatic \
        --db {input.db}
        --bed {medexome} \
        --ref {ref} \
        --max-ins-size 100 \
        --max-del-size 100 \
        --min-coverage-normal 20 \
        --min-phred-fisher 20 > {output.vcf}
        """

checkpoint Scalpel_copy:
    input:
        indels=scalpel_output_dir+"{tumor}_vs_{normal}/somatic_filtered.indel.vcf"
    output:
        indels=output_variants_notfiltered+"{tumor}_vs_{normal}/scalpel_indels.vcf"
    params:
        folder=output_variants_notfiltered+"{tumor}_vs_{normal}/"
    resources:
        mem_mb=12000,
        cpus=1
    threads:
        1
    shell:
        """
        mkdir -p {params.folder}
        cp {input.indels} {output.indels}
        """