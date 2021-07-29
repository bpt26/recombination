#!/usr/bin/env python2
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import numpy
from numpy import random
import gzip
import math
import re
import matplotlib

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg

import numpy as np
import scipy.stats as stats
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg
import seaborn as sns

"""
Get only unique recombinant nodes, pull russ p-val from them, rank them, then subject them to FDR, then make that cool plot

- need to get total number of nodes with starting parsimony
- then multiply that number by the proportion of nulls that we see with that improvement
- then the number in the final data that we actually see
^ this is our final table! then we can easily report everything
"""

##########################
##### MAIN FUNCTIONS #####
##########################

def getFDR():
    parsimonyToTotalN2 = {}
    with open('pars_leafcounts.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0] == 'parsimony':
                if not int(splitLine[0]) in parsimonyToTotalN2:
                    parsimonyToTotalN2[int(splitLine[0])] = 0
                if int(splitLine[1]) >= 2:
                    parsimonyToTotalN2[int(splitLine[0])] += int(splitLine[2])

    russOrigParsToTotal = {}
    russNull = {}
    with open('russ_null.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            if len(splitLine) == 1 and splitLine[0].isdigit():
                myOrigPars = int(splitLine[0])
                russNull[myOrigPars] = {}
                russOrigParsToTotal[myOrigPars] = 0
            elif len(splitLine) == 2:
                (russNull[myOrigPars])[int(splitLine[0])] = int(splitLine[1])
                russOrigParsToTotal[myOrigPars] += int(splitLine[1])

    parsImpToProportion = {}
    parsImpToDiscovery = {}
    for n in sorted(russOrigParsToTotal):
        for i in sorted(russNull[n]):
            if (russNull[n])[i] == 0:
                parsImpToProportion[str(n)+':'+str(i)] = 1.0/(russOrigParsToTotal[n])
            else:
                parsImpToProportion[str(n)+':'+str(i)] = (russNull[n])[i]/(russOrigParsToTotal[n])
            parsImpToDiscovery[str(n)+':'+str(i)] = {}

    with open('finalRecombNodesSet.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not str(splitLine[10])+':'+str(int(splitLine[10])-int(splitLine[11])) in parsImpToDiscovery:
                    (parsImpToDiscovery[str(splitLine[10])+':'+str(int(splitLine[10])-int(splitLine[11]))]) = {}
                (parsImpToDiscovery[str(splitLine[10])+':'+str(int(splitLine[10])-int(splitLine[11]))])[splitLine[0]] = True

    myOutString = ''
    for ni in sorted(parsImpToDiscovery):
        temp = toInt(ni.split(':'))
        n, i = temp[0], temp[1]
        if len(parsImpToDiscovery[ni]) > 0:
            if n in sorted(russOrigParsToTotal):
                if not ni in parsImpToProportion:
                    parsImpToProportion[ni] = 1.0/(russOrigParsToTotal[n])
                if n in parsimonyToTotalN2:
                    if parsImpToProportion[ni] < 0.05:
                        myOutString += joiner([n,i,parsimonyToTotalN2[n],parsImpToProportion[ni],(float(parsimonyToTotalN2[n])*parsImpToProportion[ni]),len(parsImpToDiscovery[ni])])+'\n'
                else:
                    if parsImpToProportion[ni] < 0.05:
                        myOutString += joiner([n,i,0,parsImpToProportion[ni],(0.0*parsImpToProportion[ni]),len(parsImpToDiscovery[ni])])+'\n'
    open('fdr_table_russ.txt','w').write(myOutString)

    for k in parsImpToDiscovery:
        print(k, parsImpToDiscovery[k])

##########################
#### HELPER FUNCTIONS ####
##########################

def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))


def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
    return(myReturn)

def joinerU(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '_'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ','.join(newList)

#########################
##### FUNCTION CALL #####
#########################

def main():
    getFDR()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit







