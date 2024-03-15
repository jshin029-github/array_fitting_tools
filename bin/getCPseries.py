from cpfiletools import generate_CPseries_files
import sys

## usage: python getCPseries.py [CPseq] [series.txt] [CPseries] [tileno] 

cpSeqFilename = sys.argv[1]
seriesInfo = sys.argv[2]
CPseriesFilename = sys.argv[3]
tile = int(sys.argv[4])
allRNA = "" # Not sure what this is supposed to be

fin = open(seriesInfo)
bindingSeries = []
for f in fin:
    bindingSeries.append(f)

print bindingSeries

generate_CPseries_files(cpSeqFilename, allRNA, bindingSeries, CPseriesFilename, tile)
    


