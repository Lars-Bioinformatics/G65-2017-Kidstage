#!/usr/bin/env python
import sys
import numpy as np
from pickle import load

def load_model(filename):
    '''
    Load in Sklearn trained SVM
    '''
    with open(filename, 'r') as fd:
        svm, threshold, mean, std = load(fd)
        return svm, threshold, mean, std
    raise Exception("Unable to open pickle. {}".format(filename))

def transform_x(Xs, mean, std):
    ''' 
    Transform SVM features
    '''
    Xs[:, :3] = (Xs[:, :3] - mean)/std #Standardize VBQ, MRBQ, and MAPQ by trained SVM's means and stds
    Xs[:, 3] = (Xs[:, 3] - 75.5) / 151. #Normalize PIR to be on interval [-0.5, 0.5]
    return Xs 

def pir_transform(pir, orient):
    ''' 
    Convert reverse reads to the proper position in read
    '''
    return abs(float(pir) - (int(orient) * 151))

def prepare_read(read, mean, std):
    '''
    Split read into features and prepare feature vector
    '''
    chrm, pos, ref, alt, readid, mrbq, vbq, mapq, pir, frag, orient, template_pos = read.strip().split("\t")
    features = np.asarray([[float(vbq), float(mrbq), int(mapq), pir_transform(pir, orient), int(template_pos)]])
    return transform_x(features, mean, std)


def predict(svm, X_read):
    '''
    Predict a reads call status with its given feature vector
    '''
    score = svm.decision_function(X_read)[0]
    return score


def get_reads(filestream):
    '''
    Generator that returns reads converted into read-level feature vectors
    '''
    for line in filestream:
        yield line

def main():
    '''
    Arguments: 
        --detections - A TSV output from BamReadFilt
        --pickle-name - A trained SVM model pickle
    Output:
        A stdout of each line in the TSV with its SVM score and call
    '''
    try:
        import argparse
        parser= argparse.ArgumentParser()
        parser.add_argument('--detections', dest='detections', action='store',
                                    required=True,
                                    help='Detections from BamReadFilt using the set of detections for the bam of interest')
        parser.add_argument('--pickle-name', dest='picklename', action='store', required=True,
                                    help='Name of the pickle with the SVM parameters')

        parser.add_argument('--output_file', dest='output_file', action='store', required=True, help='output file with svm scores appended')
        flags = parser.parse_args()

    except ImportError:
        flags = None
        return -1

    svm, threshold, mean, std = load_model(flags.picklename)

    open_output_file=open(flags.output_file,"w")

    with open(flags.detections, 'r') as filestream: 
        for read in get_reads(filestream):
            read_feature_vector = prepare_read(read, mean, std)
            prediction = predict(svm, read_feature_vector)
            call = int(prediction > threshold)
            open_output_file.write(read.strip() + "\t{}\t{}".format(prediction, call)+"\n")
    open_output_file.close()


    print "Done."
if __name__ == "__main__":
    main()
