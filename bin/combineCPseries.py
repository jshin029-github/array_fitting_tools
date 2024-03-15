## Usage: python combineCPseries.py [N tiles] [N concs] [CPseries Dirs] [Out Dir]
import sys, os
import numpy as np
import pandas as pd

Ntiles = int(sys.argv[1])

TileFilenames = []
for i in range(1, Ntiles+1):
    TileFilenames.append("ALL_tile%03d_Bottom_filtered.CPseries"%(i))
print "Tiled Filenames:"
for n in TileFilenames:
    print "\t%s"%(n)

Nconcs = int(sys.argv[2])

CPseriesDirs = []
for j in range(Nconcs):
    CPseriesDirs.append(sys.argv[3+j])
print "CPseries Dir Filenames:"
for n in CPseriesDirs:
    print "\t%s"%(n)

print "\n\n"

outDir = sys.argv[-1]
outfilenames = []

for TileFilename in TileFilenames:
    outFilename = os.path.join(outDir, TileFilename)
    counter = 0
    FinalDF = None
    filenotfound = False
    for CPseriesDir in CPseriesDirs:
        print "Reading in CPseries file: %s/%s"%(CPseriesDir, TileFilename) 
        Filename = os.path.join(CPseriesDir, TileFilename)
        if not os.path.exists(Filename):
            print "Did not find the above file. Skipping ahead."
            #FinalDF[counter] = pd.Series()
	    filenotfound = True
            continue
        df = pd.read_csv(Filename, sep = "\t", header = 0)
        df.columns = ['1', '2', '3', '4', '5', '6', '7', '8' ]
	df = df.loc[df['2'] == 'anyRNA']
        print "\t No of clusters found  = %d"%df.shape[0]
        if counter == 0:
            FinalDF = pd.DataFrame(df.iloc[:, [0, -1]])
            FinalDF.columns = ['clusterID', 0]

        else:
            FinalDF[counter] = df.iloc[:, -1]
        
        counter += 1
    if not filenotfound :
        outfilenames.append(outFilename)
        FinalDF.to_csv(outFilename, sep = '\t', header=True, index=False, na_rep="nan")

         
