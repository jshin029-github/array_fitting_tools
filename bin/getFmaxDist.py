#!/usr/bin/env python
#
# Sarah Denny
# July 2015

##### IMPORT #####
import numpy as np
import pandas as pd
import sys
import os
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from fittinglibs import fileio, initfits, plotting

NPoints = 8
### MAIN ###

################ Parse input parameters ################

#set up command line argument parser
parser = argparse.ArgumentParser(description='get fmax dist')
group = parser.add_argument_group(description="Returns a csv file with all cluster normalized fluorescence used for fitting annotated by variant number")

group.add_argument('-fp', '--fmax', 
                   help='File containing Fmax distribution (fmaxdist.p)')

group = parser.add_argument_group(description='Optional additional arguments.')
group.add_argument('-out', '--out_dir', default='./result.csv.gz',
                   help='output filename. default is current directory result.csv.gz')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.fmax:
        fmax = fileio.loadFile(args.fmax)
        fmaxparams = fmax.params
        fout = open(args.out_dir, "w")
        for k in fmaxparams.keys():
            fout.write("%s\t%f\n"%(fmaxparams[k].name, fmaxparams[k].value))
