#!/usr/bin/env python
""" Extracts data from merged, CPseries file for red and green and normalizes.

Sarah Denny

"""

import os
import numpy as np
import pandas as pd
import argparse
import sys
import itertools
import scipy.stats as st
import logging
# from fittinglibs import (plotting, fileio, processing)
from fittinglibs import (fileio, processing)


##### PARSE ARGUMENTS #####
parser = argparse.ArgumentParser(description='normalize fluorescence series in green channel '
                                 '(associated with binding signal) by the signal in the red channel '
                                 '(associated with total transcribed RNA)')
processing.add_common_args(parser.add_argument_group('common arguments'))

group = parser.add_argument_group('optional arguments for normalize series')
group.add_argument('--no_bounds', action="store_true",
                   help='By default, all cluster signal is bounded to prevent '
                   'dividing by zero in cases where signal is low. Flag to prevent this.')
group.add_argument('--bounds', nargs=2, metavar='N N',
                   help='use these lower and upper bounds if provided. ',
                   type=float)

# John Edits 8/6/2021
group.add_argument('--first_condition', action="store",
                   help='By default, each condition has its green fluorescence divided by red fluorescence. This option makes it such that only the first red fluorescence is used',
                   type = int)
# End John Edits


##### FUNCITONS #####
def boundFluorescence(signalF, plot=False, bounds=None, figDir = None):
    # take i.e. all cluster signal and bound it
    signalFull = signalF.copy()
    for i in range(signalFull.shape[1]):
        signal = signalFull.iloc[:,i]
        print "For set: ", i
        # check if at least one element of signal is not nan
        if np.isfinite(signal).sum() > 0:
            if bounds is None:
                lowerbound = np.percentile(signal.dropna(), 1)
                upperbound = signal.median() + 5*signal.std()
                logging.info('Found bounds: %4.3f, %4.3f'%(lowerbound, upperbound))
            else:
                lowerbound = bounds[0]
                upperbound = bounds[1]
                logging.info('Using given bounds: %4.3f, %4.3f'%(lowerbound, upperbound))

            # if plot and figDir:
            #     plotting.plotBoundFluorescence(signal, [lowerbound, upperbound])
            #     plotting.savefig(os.path.join(figDir, 'bound_all_cluster_signal_%d.pdf'%(i)))

            signal.loc[signal < lowerbound] = lowerbound
            signal.loc[signal > upperbound] = upperbound
            signalFull.iloc[:,i] = signal

        else:
            #if they are all nan, set to 1 for division
            signal.loc[:] = 1
    return signalFull

##### MAIN #####
if __name__=="__main__":


    args = parser.parse_args()
    processing.update_logger(logging, args.log)

    ##### DEFINE OUTPUT FILE AND DIRECTORY #####
    # find out file and fig directory
    if args.out_file is None:
        args.out_file = fileio.stripExtension(args.binding_series) + '_normalized.CPseries.gz'

    figDirectory = os.path.join(os.path.dirname(args.out_file), fileio.returnFigDirectory())
    if not os.path.exists(figDirectory):
        os.mkdir(figDirectory)
    #print "1: ", args.ref_fluor_series
    #print "2: ", args.binding_series

    # use only first columns of allCluster signal file
    logging.info("Loading all RNA signal")
    #allClusterSignal = fileio.loadFile(args.ref_fluor_series).iloc[:, 0]
    allClusterSignal = fileio.loadFile(args.ref_fluor_series)

    # laod whole binding Series
    logging.info("Loading series...")
    bindingSeries = fileio.loadFile(args.binding_series)

    # DEBUGGING:
    #print "DEBUGGING INFORMATION: "
    #print allClusterSignal
    #print bindingSeries

    # John Edits 8/6/2021

    if args.first_condition | args.first_condition == 0:
        first_condition = allClusterSignal.iloc[:,args.first_condition]

        for col in allClusterSignal.columns:
            allClusterSignal[col] = first_condition

    # End John Edits

    # make normalized binding series
    logging.info("Normalizing...")
    if not args.no_bounds:
        allClusterSignal = boundFluorescence(allClusterSignal, plot=True, bounds=args.bounds, figDir = figDirectory)

    #DEBUGGING:
   # print "POST filtering / masking of Red channel"


    #bindingSeriesNorm = np.divide(bindingSeries, np.vstack(allClusterSignal))
    bindingSeriesNorm = np.divide(bindingSeries, allClusterSignal)

    # save
    logging.info( "Saving...")
    bindingSeriesNorm.to_csv(args.out_file, sep='\t', compression='gzip')
    #bindingSeriesNorm.to_csv(args.out_file, sep='\t', compression='gzip', na_rep='nan')
    #plotting.savefig(os.path.join(figDirectory, 'bound_all_cluster_signal.pdf'))

