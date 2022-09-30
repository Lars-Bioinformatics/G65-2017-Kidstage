import sys, re

print(snakemake.input[0])
print(snakemake.input[1])
print(snakemake.params[0])
print(len(snakemake.input))


# input_dir = "/work/G65-2017-Kidstage/Connor/mutect2_somatic_variants/"
# bed = input_dir+"regions_bam_readcount/G65-T33B-4315_nimblegen-medexome_HYLKFDSXX_snvs_pass_fixed.bed"
# snvs_readcount = input_dir+"bam_readcount_output_q20_b20/G65-T33B-4315_nimblegen-medexome_HYLKFDSXX/snvs/G65-T33B-4315_nimblegen-medexome_HYLKFDSXX_snvs_readcounts_q20_b20.txt"

# Output file destination
outfile = snakemake.output[0]
# outfile = input_dir+"bam_readcount_output_q20_b20/G65-T33B-4315_nimblegen-medexome_HYLKFDSXX/G65-T33B-4315_bam-readcount_VAF.txt"

# Mutation type - snv, insertion, deletion
type = snakemake.params[0]

# Define header for output
header = ["sample","tumor","tissue","type","chr","pos","ref","alt","RD","AD","VAF"]

# Get ALT alleles from bed file
bed = snakemake.input[0]
with open(bed,"r") as f:
    alt_alleles = [line.split()[4] for line in f]

# print(alt_alleles)

def get_sample_name(filename):
    sample = filename.rsplit("/",1)[1].split("_",1)[0]
    return(sample)

def get_tissue(sample, tumor):
    if (sample == tumor):
        tissue = "tumor"
    # Check if sample number is the same i.e. it's plasma belonging to the tumor 
    elif (re.sub("\D","",sample.split("-")[1]) == re.sub("\D","",tumor.split("-")[1])):
        tissue = "plasma"
    else:
        tissue = "plasma-pool"
    return(tissue)

# Extract readcount from bam-readcount line, already split on whitespace
def get_read_depth(base, readcount_split, type="snvs"):
    base = base[0]
    if type == "insertions":
        # Expects exactly 11 fields - however, 
        # not sure if multiple insertions at give position could exists
        if len(readcount_split) > 10:
            index = 10
        else:
            return(-1)
    elif base == "A":
        index = 5
    elif base == "C":
        index = 6
    elif base == "G":
        index = 7
    elif base == "T":
        index = 8
    DP = readcount_split[index].split(":",2)[1]
    return(DP)

with open(outfile, "w") as output:
    # write header to output
    output.write("\t".join(header)+"\n")
    # get tumor name
    tumor = get_sample_name(bed)
    for i in range(1,len(snakemake.input)):
        snvs_readcount = snakemake.input[i]
        # get sample name
        sample = get_sample_name(snvs_readcount)
        # Get tissue type
        tissue = get_tissue(sample, tumor)
        with open(snvs_readcount,"r") as f:
            for linenum, line in enumerate(f):
                readcount = line.split()
                # create output with tumor, sample, chr, pos
                res = [sample, tumor, tissue, type] + readcount[0:2]
                # add ref and alt to output
                ref = readcount[2]
                alt = alt_alleles[linenum]
                # get RD and AD
                RD = get_read_depth(ref, readcount, "snvs")
                AD = get_read_depth(alt, readcount, type)
                # Compute VAF
                if AD == -1: 
                    continue
                elif AD == '0':
                    VAF = 0
                else:
                    VAF = int(AD) / (int(RD)+int(AD))
                # Append to output
                res += [ref,alt,RD,AD,str(VAF)]
                print('\t'.join(res))
                output.write("\t".join(res)+"\n")



