import numpy as np
import os

# find all csv files in the current directory
allfiles = os.listdir(os.getcwd())
infiles = list()
for i in allfiles:
    if '.csv' in str.lower(i):
        infiles.append(i)

for cf in infiles:
    print 'opening --> ' + cf
    indat = open(cf,'r').readlines()
    innames = indat.pop(0)
    innames = innames.strip().split(',')
    ofp = open(cf[:-4] + '.gslib','w')
    ofp.write('HUCs\n%d\n' %(len(innames)))
    for cn in innames:
        ofp.write('%s\n' %(cn))
    for line in indat:
        for i in line.strip().split(','):
            ofp.write('%s ' %(i))
        ofp.write('\n')
    ofp.close()