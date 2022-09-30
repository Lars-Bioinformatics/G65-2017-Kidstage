#!/usr/bin/python

import sys
import csv
import subprocess

snp_database     = sys.argv[1] # /scratcha/nrlab/resources/1000genomes/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
cosmic_database  = sys.argv[2] # /scratcha/nrlab/resources/COSMIC/v82/CosmicCodingMuts.vcf.gz
mutation_file    = sys.argv[3]

# Variables for trinucleotide context
samtools="/home/nrlab/tools/samtools/current/samtools"
reference_genome="/scratcha/nrlab/resources/ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta"

# The mutation file needs the following column headings
# CHROM    POS    REF    ALT

mutation_fileObj = open(mutation_file, 'r')
mutation_fileReader = csv.DictReader(mutation_fileObj, delimiter="\t")

print "\t".join(['CHROM','POS', 'REF', 'ALT',
                     'DP', 'DP4', 'REF_F', 'ALT_F', 'REF_R', 'ALT_R', 'MQSB', 'SLX', 'BARCODE',
                     'COSMIC_MUTATIONS', 'COSMIC_SNP', '1KG_AF', 'TRINUCLEOTIDE'])

for row in mutation_fileReader:
    chromosome = row['CHROM'].upper()
    position   = row['POS']
    reference  = row['REF']
    mutation   = row['ALT']

    # ensembl format chromosome for cosmic vcf
    chromosome_short  = chromosome.replace("CHR", "")
    mutation_position = chromosome_short+":"+position+"-"+position

    # assume tabix is available on the $PATH
    tabix_snp_command     = " ".join(["tabix",snp_database,mutation_position])
    tabix_snp_command_out = subprocess.check_output(tabix_snp_command, shell=True).rstrip()

    k1g_af = 0
    if tabix_snp_command_out:
        for vcf_line in tabix_snp_command_out.split("\n"):
            vcf_values    = vcf_line.split("\t")
            k1g_chromsome = vcf_values[0]
            k1g_position  = vcf_values[1]
            k1g_reference = vcf_values[3]
            k1g_alt       = vcf_values[4]
            k1g_filter    = vcf_values[6]
            k1g_info      = vcf_values[7]

            ## for boolean field use 1 as the value
            infoDict = {}
            for pair in k1g_info.split(';'):
                try:
                    key,value = pair.split('=')
                    infoDict[key] = value
                except ValueError:
                    infoDict[pair] = '1'

            k1g_vt = infoDict['VT']

            # sanity check values match
            if (chromosome_short == k1g_chromsome and
                position   == k1g_position and
                reference  == k1g_reference and
                mutation   == k1g_alt and
                k1g_filter == "PASS" and
                k1g_vt     == "SNP"):

                k1g_af = infoDict['AF']
                break


    tabix_cosmic_command     = " ".join(["tabix",cosmic_database,mutation_position])
    tabix_cosmic_command_out = subprocess.check_output(tabix_cosmic_command, shell=True).rstrip()

    cosmic_sampleN = 0
    cosmic_snp = 0
    if tabix_cosmic_command_out:
        for cosmic_line in tabix_cosmic_command_out.split("\n"):
            cosmic_values    = cosmic_line.split("\t")
            comic_chromsome  = cosmic_values[0]
            cosmic_position  = cosmic_values[1]
            cosmic_reference = cosmic_values[3]
            cosmic_alt       = cosmic_values[4]
            cosmic_filter    = cosmic_values[6]
            cosmic_info      = cosmic_values[7]
            ## for boolean field use 1 as the value
            infoDict = {}
            for pair in cosmic_info.split(';'):
                try:
                    key,value = pair.split('=')
                    infoDict[key] = value
                except ValueError:
                    infoDict[pair] = '1'
            cosmic_sampleN = infoDict['CNT']

            try:
                cosmic_snp = infoDict['SNP']
            except KeyError:
                cosmic_snp = 0

            # sanity check values match
            if (chromosome == comic_chromsome and
                position   == cosmic_position and
                reference  == cosmic_reference and
                mutation   == cosmic_alt and
                cosmic_filter == "PASS"):

                cosmic_sampleN = infoDict['CNT']
                break

    trinucleotide_position = chromosome[0:3].lower() + chromosome[3:] + ":" + str(int(position) - 1) + "-" + str(int(position) + 1)
    trinucleotide_command     = " ".join([samtools, "faidx", reference_genome, trinucleotide_position,  "| tail -n -1"])
    trinucleotide_command_out = subprocess.check_output(trinucleotide_command, shell=True).rstrip()


    print "\t".join([row['CHROM'], row['POS'], row['REF'], row['ALT'],
                         row['DP'], row['DP4'],
                         row['REF_F'], row['ALT_F'], row['REF_R'], row['ALT_R'],
                         row['MQSB'], row['SLX'], row['BARCODE'],
                         str(cosmic_sampleN), str(cosmic_snp), str(k1g_af), trinucleotide_command_out.upper()])
