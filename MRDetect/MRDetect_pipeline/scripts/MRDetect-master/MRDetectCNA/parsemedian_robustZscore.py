#!/usr/bin/python

########################################
# Minita Shah (mshah@nygenome.org)
# New York Genome Center
########################################

import os
import sys
import argparse
import re
import pandas as pd
import numpy as np

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))

def normalize(col, median, mad):
    return float((col - median)/(1.4826*mad))

def determine_segment_direction(col):
    if col == 'DEL':
        return -1
    elif col=='DUP':
        return 1

################### Main ###################
parser = ArgumentParser(prog='parsemedian_robustZscore.py', description='', epilog='')
parser.add_argument('-i', '--inp-dir', help='Input directory.', required=True)
parser.add_argument('-t', '--tumor', help='Tumor sample.', required=True)
parser.add_argument('-n', '--normal', help='Normal sample.', required=True)
parser.add_argument('-in', '--interval-name', help='Interval name.', required=True)
parser.add_argument('-o', '--output-dir', help='Output Directory', required=True)
parser.add_argument('-w', '--wsize', help='Segment window size.', required=True)
parser.add_argument('-mt', '--med-thresh', help='Windows above this median threshold in the normal will be removed.', required=True)
parser.add_argument('-tc', '--tumor-cvg', help='Tumor sample avg. coverage.', required=True)
parser.add_argument('-nc', '--normal-cvg', help='Normal sample avg. coverage.', required=True)

args = parser.parse_args()

INPDIR = args.inp_dir
OUTDIR = args.output_dir
WSIZE = int(args.wsize)
tumor = args.tumor
normal = args.normal
interval_name = args.interval_name
med_thresh = float(args.med_thresh)
tumor_avg_cvg = float(args.tumor_cvg)
normal_avg_cvg = float(args.normal_cvg)

#for each plasma-normal pair, aggregate score across segments
out_summary = os.path.join(OUTDIR, tumor+"."+normal+".win"+str(WSIZE)+".robustZscore.delsdups.summary.txt")
f_out = open(out_summary, 'w')
f_out.write("Tumor\tNormal\tWindowSize(bp)\tFinalScore\tMeanScore\tDelFinalScore\tDelMeanScore\tAmpFinalScore\tAmpMeanScore\n")

outfile = os.path.join(OUTDIR, tumor+"."+normal+".win"+str(WSIZE)+".robustZscore.delsdups.perwindow.mediandiff")

DOC_T = os.path.join(INPDIR, tumor, interval_name, "segments", tumor+".win"+str(WSIZE)+".delsdups.perwindow.median")
DOC_N = os.path.join(INPDIR, normal, interval_name, "segments", normal+".win"+str(WSIZE)+".delsdups.perwindow.median")

print "DOC_T=" + DOC_T
print "DOC_N=" + DOC_N

df_tumor = pd.read_table(DOC_T, sep="\t", header=0)
df_normal = pd.read_table(DOC_N, sep="\t", header=0)
df_tn = pd.concat([df_tumor, df_normal['Median']], axis=1)
df_tn.columns = ['Chunk_start', 'Orig_seg_type', 'Median_tumor', 'Median_normal'] #Rename column headers

df_tn['Median_tumor'] = df_tn['Median_tumor']/tumor_avg_cvg
df_tn['Median_normal'] = df_tn['Median_normal']/normal_avg_cvg

#Subset of rows where median tumor and median normal not 0
df_tn = df_tn[(df_tn['Median_tumor'] != 0) & (df_tn['Median_normal'] != 0)]

#CNV med, mad
tumor_median_new = df_tn['Median_tumor'].median()
normal_median_new = df_tn['Median_normal'].median()
tumor_mad_new = (abs(df_tn['Median_tumor'] - tumor_median_new)).median()
normal_mad_new = (abs(df_tn['Median_normal'] - normal_median_new)).median()

print "Tumor_med_new = " + str(tumor_median_new) + "\t" + "Normal_med_new = " + str(normal_median_new) + "\t" + "Tumor_mad_new = " + str(tumor_mad_new) + "\t" + "Normal_mad_new = " + str(normal_mad_new)

print "Input windows = " + str(df_tumor.count()[0])
print "After removing 0 bins = " + str(df_tn.count()[0])

df_tn['Norm_tumor'] = df_tn['Median_tumor'].apply(normalize, args=(tumor_median_new, tumor_mad_new))
df_tn['Norm_normal'] = df_tn['Median_normal'].apply(normalize, args=(normal_median_new, normal_mad_new))

df_tn = df_tn[(df_tn['Norm_normal'] <= med_thresh) & (df_tn['Norm_normal'] >= -med_thresh)]

print "Windows after high coverage filtering: \t" + str(df_tn.count()[0])

df_tn['Seg_direction'] = df_tn['Orig_seg_type'].apply(determine_segment_direction)

df_tn['Signal_diff'] = (df_tn['Norm_tumor'] - df_tn['Norm_normal'])
sum_of_signal = ((df_tn['Norm_tumor'] - df_tn['Norm_normal']) * df_tn['Seg_direction']).sum()
chunk_cnt = df_tn.count()[0]
df_rows = df_tn[df_tn['Seg_direction'] == 1]
sum_of_signal_dup = ((df_rows['Norm_tumor'] - df_rows['Norm_normal']) * df_rows['Seg_direction']).sum()
chunk_cnt_dup = df_rows.count()[0]
df_rows = df_tn[df_tn['Seg_direction'] == -1]
sum_of_signal_del = ((df_rows['Norm_tumor'] - df_rows['Norm_normal']) * df_rows['Seg_direction']).sum()
chunk_cnt_del = df_rows.count()[0]

df_tn.to_csv(outfile, sep='\t', index=False)

f_out.write(tumor + "\t" + normal + "\t" + str(WSIZE) + "\t" + str(sum_of_signal) + "\t" + str(sum_of_signal/chunk_cnt) + "\t" + str(sum_of_signal_del) + "\t" + str(sum_of_signal_del/chunk_cnt_del) + "\t" + str(sum_of_signal_dup) + "\t" + str(sum_of_signal_dup/chunk_cnt_dup) +"\n")
f_out.close()
