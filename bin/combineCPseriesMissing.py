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

# Concat all data frames vertically first:
DFseries = []
cpcount = 0

print " ***** Merging Vertically ***** "

for CPseriesDir in CPseriesDirs:
    print "Reading tiles from : %s" % (CPseriesDir)
    cpcount += 1
    # For each conc. point => read all available tile CPseries files and stack vertically
    firstPass = True
    DF  = None
    for TileFilename in TileFilenames:
	# If file doens't exist - continue on to next file
	if not os.path.exists(os.path.join(CPseriesDir, TileFilename)):
	   print  "	File not found : %s " % (TileFilename)
   	   print "	Skipping processing steps"
	   
	else:
	    # Otherwise, files exists => Read and format it;
	    print "		Reading in file: %s "%(TileFilename)
	    d = pd.read_csv(os.path.join(CPseriesDir, TileFilename), sep='\t', header=None)
	    d = d.loc[d[1] == 'anyRNA']
	    d = d.iloc[:, [0, -1]]
	    d.columns = ['clusterID', cpcount-1]
	    # File exists and is the first file in line
	    if  firstPass: ## Fist tile => read into DF
	        firstPass = False
	        DF = pd.DataFrame(d)
	    else:
	        DF = pd.concat([DF, d])
    DF.set_index('clusterID')
    print "		Final nrows = %d\n\n"%(DF.shape[0])
    print DF.head(2)
    print " -------------------------------------------- "
    DFseries.append(DF)

# Merge vertically
print " ***** Merging Horizontally ***** "
DFCombined = pd.DataFrame(DFseries[0])

for i in range(1, len(DFseries)):
    print "Merging conc point %d"%(i)
    DFCombined = pd.merge(DFCombined, DFseries[i], how='outer', on='clusterID')

DFCombined.to_csv(outDir, sep='\t', header=True, index=False)


