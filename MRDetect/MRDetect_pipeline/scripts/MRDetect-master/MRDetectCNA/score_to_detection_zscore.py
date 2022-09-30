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
from collections import defaultdict
import numpy as np

########## Classes and Functions ##########

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
	self.print_help(sys.stderr)
	self.exit(2, "\nERROR: %s\n\n" % (message))

################### Main ###################
parser = ArgumentParser(prog="score_to_detection_zscore.py", description="", epilog="")
parser.add_argument("-s", "--summary-dir", help="Path to input summary directory.", required=True)
parser.add_argument("-t", "--tn-pairs", help="Patient specific tumor-normal pairs.", required=True)
parser.add_argument("-o", "--out-file-name", help="Output file name.", required=True)

args = parser.parse_args()
summary_dir = args.summary_dir
tn_pairs = args.tn_pairs
out_file_name = args.out_file_name

control_stats = defaultdict(dict)
patient_stats = defaultdict(dict)

tf = open(tn_pairs)
for tf_line in tf:
    if tf_line.startswith("#"):
        continue
    else:
        tf_toks = tf_line.strip().split("\t")
        pon_name = tf_toks[1]
        inp_file = os.path.join(summary_dir, tf_toks[0]+"."+tf_toks[1]+".win500.robustZscore.delsdups.summary.txt")
        f = open(inp_file)
        for line in f:
            if line.startswith("Tumor"):
                continue
            else:
                toks = line.strip().split("\t")
                if tf_toks[-1] == "control":
                    control_stats['total'][toks[0]] = float(toks[4])
                    control_stats['del'][toks[0]] = float(toks[6])
                    control_stats['amp'][toks[0]] = float(toks[8])
                else:
                    patient_stats[toks[0]]['total'] = float(toks[4]) 
                    patient_stats[toks[0]]['del'] = float(toks[6])
                    patient_stats[toks[0]]['amp'] = float(toks[8])
        f.close()
tf.close()

mean_noise_total = np.mean(control_stats['total'].values())
mean_noise_del = np.mean(control_stats['del'].values())
mean_noise_amp = np.mean(control_stats['amp'].values())
std_noise_total = np.std(control_stats['total'].values(), ddof=1)
std_noise_del = np.std(control_stats['del'].values(), ddof=1)
std_noise_amp = np.std(control_stats['amp'].values(), ddof=1)

out_file = os.path.join(summary_dir, out_file_name)
f_out = open(out_file, 'w')
f_out.write("Tumor\tNormal\tTotalScore\tDelMeanScore\tAmpMeanScore\tTotalZScore\tDelZScore\tAmpZScore\n")

for i in patient_stats.keys():
    f_out.write(i + "\t" + pon_name + "\t" + str(patient_stats[i]['total']) + "\t" + str(patient_stats[i]['del']) + "\t" + str(patient_stats[i]['amp']) + "\t" + str((patient_stats[i]['total']-mean_noise_total)/std_noise_total) + "\t" + str((patient_stats[i]['del']-mean_noise_del)/std_noise_del) + "\t" + str((patient_stats[i]['amp']-mean_noise_amp)/std_noise_amp) + "\n")

for i in sorted(control_stats['total'].keys()):
    f_out.write(i + "\t" + pon_name + "\t" + str(control_stats['total'][i]) + "\t" + str(control_stats['del'][i]) + "\t" + str(control_stats['amp'][i]) + "\t" + str((control_stats['total'][i]-mean_noise_total)/std_noise_total) + "\t" + str((control_stats['del'][i]-mean_noise_del)/std_noise_del) + "\t" + str((control_stats['amp'][i]-mean_noise_amp)/std_noise_amp) + "\n")

f_out.close()
