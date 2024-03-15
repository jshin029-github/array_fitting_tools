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
parser = argparse.ArgumentParser(description='bootstrap fits')
group = parser.add_argument_group(description="Returns a csv file with all cluster normalized fluorescence used for fitting annotated by variant number")

group.add_argument('-vp', '--variantparams', 
                   help='file containing the variant params (variantParams.p)')

group.add_argument('-n', '--num_points', default='8' ,
                   help='Number of concentration points in binding curve')

group = parser.add_argument_group(description='Optional additional arguments.')
group.add_argument('-out', '--out_dir', default='./result.csv.gz',
                   help='output filename. default is current directory result.csv.gz')

if __name__ == '__main__':
    args = parser.parse_args()
    
    if args.variantparams:
        NPoints = int(args.num_points)
        print "Using N = ", str(NPoints)
        # load only the variant params if given
        variantParams = fileio.loadFile(args.variantparams)
        variants = variantParams.variants
        F = {}
        Variants = []
        for i in range(NPoints):
            F[i] = []
        clusterIDs = []
        for variant in variants:
            print "Processing variant %d / %d : "%(variant, variants[-1])
            FL = variantParams.get_ys(variant)
            Indices = FL.index.values.tolist()
            clusterIDs += Indices
            for i in range(NPoints):
                F[i] += FL[str(i)].tolist()
            Variants += [variant] * FL.shape[0]

        # Create Data Frame containing all info:
        DF = pd.DataFrame()
        DF['clusterID'] = clusterIDs
        DF['variant'] = Variants
        for i in range(NPoints):
            DF[str(i)] = F[i]
        print DF.head(10)
        DF.to_csv(args.out_dir, sep="\t", index=False, compression="gzip")
