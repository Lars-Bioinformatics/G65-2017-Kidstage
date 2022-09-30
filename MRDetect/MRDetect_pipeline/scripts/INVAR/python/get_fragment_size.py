#! /usr/bin python

# Purpose: 	Split a BAM file at a SNV. Ensure that the BAM is indexed
#			Write the reads containing alternate alleles at SNVs to new BAM
#			Write the reads containing reference alleles to 2nd new BAM
# Author: Dineika Chandrananda (dineika.chandrananda@cruk.cam.ac.uk)
#			with help from Pete Hickey at John Hopkins (phickey4@jhu.edu)
# Date: 25/04/2017

# To run: python2.7 splitByAllele.py <BAM_FILE> <OUTPUT_PREFIX> chrZ:POSITION:REF:ALT

# Notes: Not compatible with the print function of python 3

import argparse
import pysam
import os # OS package to interact with the filesystem
import sys # allows system commands

def can_create_file(folder_path):
    try:
        tempfile.TemporaryFile(dir=folder_path)
        return True
    except OSError:
        return False

# required variables
parser = argparse.ArgumentParser(description="splitByAllele.py - "
    "Extracts reads from a BAM file that cover a specified SNV and "
    "writes the reference and alternate allele containing reads to separate BAM files.")
parser.add_argument("input_bam", help="Input BAM file. Has to be indexed.")
parser.add_argument("SNV_list", help="List of SNV positions (.txt), reference and alternate "
    "allele of interest in chr:position:ref:alt format, e.g. chr21:11106932:A:G. "
    "Coordinates are 1-based.")
parser.add_argument("output", help="Name of your output file")


# optional variables
parser.add_argument("--print_only_mut", help="Set to print only mut read sizes", default="TRUE")
parser.add_argument("--pair_distance", type=int, default=10000, help="The distance in "
    "basepairs to search up and downstream from the specified SNV/deletion for "
    "the pair of overlapping reads (default is 10000)")
parser.add_argument("--max_depth", type=int, default=1000000, help="Maximum number "
    "of reads to process at the specified SNV position (default is 1000000)")
args = parser.parse_args()

if args.max_depth < 10000:
    print("Specified max_depth is too low - changing to 10000")
    args.max_depth = 10000

# open the bam or die
try:
    samfile = pysam.AlignmentFile(args.input_bam, "rb")
except IOError:
    sys.exit("BAM file %s doesn't exist. Please provide one" % args.input_bam)

# open output or die
try:
    output_file = open(args.output, "w")
except IOError:
    sys.exit("Cannot write %s" % args.output)

# start loop to pileup at each of the SNV sites based on SNV_list
# this is to get the ref/alt status, and the fragment length (TLEN)
nucleotides = ['A', 'C', 'G', 'T']

with open(args.SNV_list, "rU") as SNVs:
    for SNV in SNVs:

        # do check on the format of the SNV_list
        try:
            chrom, pos, ref, alt = SNV.strip().split(":")
        except ValueError:
            print("SNV specified '%s' not in chr:position:ref:alt format" % SNV)

        # set case
        SNV_ref = ref.upper()
        SNV_alt = alt.upper()

        # check nucleotide validity
        if ref not in nucleotides:
            print("Reference allele %s is not A, C, G or T" % ref)
        if alt not in nucleotides:
            print("Alternate allele %s is not A, C, G or T" % alt)

        # Check validity of pos - is it an integer?
        try:
            pos = int(pos)
        except ValueError:
            print("Position '%s' is not valid" % pos)

        # get the read
        for pileup in samfile.pileup(chrom, pos-1, pos, max_depth=args.max_depth, stepper="nofilter"):
            if pileup.reference_pos == pos-1: # filter for position of interest
                print "Processing %s reads covering SNV position %s:%s in %s" % (
                    len(pileup.pileups), chrom, pos, args.input_bam)
                for read in pileup.pileups:
                    if read.query_position is not None:
                        SNV_base = read.alignment.query_sequence[read.query_position]
                        mutant_status = "."

                        if SNV_base == ref:
                            mutant_status = 0
                        elif SNV_base == alt:
                            mutant_status = 1

                        #print SNV_base
                        #print SNV.strip()
                        #print mutant_status
                        #print abs(read.alignment.template_length)

                        # write to a dataframe
                        output_file.write("%s,%s,%s,%d\n" % (SNV_base, str(SNV.strip()), mutant_status, abs(read.alignment.template_length)))


samfile.close()
output_file.close()
