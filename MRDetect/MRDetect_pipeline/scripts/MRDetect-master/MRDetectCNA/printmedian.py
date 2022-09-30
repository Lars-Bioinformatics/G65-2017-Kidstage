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

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
	self.print_help(sys.stderr)
	self.exit(2, '\nERROR: %s\n\n' % (message))

def read_seg_file(seg_file, seg_type, seg_size):
	#chr	start	end	type	info
	#1	10025	2246016	DUP	0.490794488822515
    seg = dict()
    f = open(seg_file)
    for line in f:
        if line.startswith("#"):
            continue
        else:
            toks=line.strip().split("\t")
            if seg_size:
                if seg_type == "CNV":
                    if ((toks[3] == "DUP" or toks[3] == "DEL")):
                        key = toks[0]+"_"+toks[1]+"_"+toks[2]
                        seg[key] = toks[3]
                elif seg_type == "NEU":
                    if ((toks[3] == "NEU")):
                        key = toks[0]+"_"+toks[1]+"_"+toks[2]
                        seg[key] = toks[3]
            else:
                if seg_type == "CNV":
                    if ((int(toks[2]) - int(toks[1]) + 1 > 5000000) and (toks[3] == "DUP" or toks[3] == "DEL")):
                        key = toks[0]+"_"+toks[1]+"_"+toks[2]
                        seg[key] = toks[3] 
                elif seg_type == "NEU":
                    if ((int(toks[2]) - int(toks[1]) + 1 > 5000000) and (toks[3] == "NEU")):
                        key = toks[0]+"_"+toks[1]+"_"+toks[2]
                        seg[key] = toks[3]
    f.close()
    return seg

################### Main ###################
parser = ArgumentParser(prog='printmedian.py', description='', epilog='')
parser.add_argument('-i', '--inp-dir', help='Input directory with DOC files.', required=True)
parser.add_argument('-o', '--out-dir', help='Output directory.', required=True)
parser.add_argument('-s', '--sample-name', help='Sample name (tumor/normal).', required=True)
parser.add_argument('-f', '--seg-file', help='CNV segment file.', required=True)
parser.add_argument('-w', '--wsize', help='Segment window size.', required=True)
parser.add_argument('-t', '--seg-type', help='Segment type (CNV/NEU).', required=True)
parser.add_argument('-z', help='Set if you want to include all segs including <5Mb.', action='store_true')

args = parser.parse_args()

inp_dir = args.inp_dir
out_dir = args.out_dir
SEGFILE = args.seg_file
WSIZE = int(args.wsize)
seg_type = args.seg_type
sample = args.sample_name
seg_size = args.z
seg_file = read_seg_file(SEGFILE, seg_type, seg_size)

if seg_type == "CNV":
    outfile = os.path.join(out_dir, sample+".win"+str(WSIZE)+".delsdups.perwindow.median")
elif seg_type == "NEU":
    outfile = os.path.join(out_dir, sample+".win"+str(WSIZE)+".neu.perwindow.median")
f_win = open(outfile, 'w')
f_win.write("Chunk_start\tOrig_seg_type\tMedian\n")

for s in seg_file:
    DOC_S = os.path.join(inp_dir, s+"."+seg_file[s]+"."+sample+".DepthOfCoverage")
    df_sample = pd.read_table(DOC_S, sep="\t", header=0, chunksize=WSIZE)

    for chunk_s in df_sample:
        # Generate median counts for each window
        median = chunk_s['Total_Depth'].median()
        f_win.write(str(chunk_s['Locus'].iloc[0]) + "\t" + str(seg_file[s]) + "\t" + str(median) + "\n")

f_win.close()
