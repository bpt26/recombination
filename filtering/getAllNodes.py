#!/usr/bin/env python3
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

##########################
##### MAIN FUNCTIONS #####
##########################

def getAllNodes():
    myNodes = {}
    needParent = {}
    with open('combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                myNodes[int(splitLine[0])] = True
                myNodes[int(splitLine[3])] = True
                myNodes[int(splitLine[6])] = True
                if splitLine[4] == 'y':
                    needParent[int(splitLine[3])] = True
                if splitLine[7] == 'y':
                    needParent[int(splitLine[6])] = True

    alreadyGotParent = {}
    lineCounter = 0
    with open('sample_paths.txt') as f:
        for line in f:
            lineCounter += 1
            splitLine = (line.strip()).split('\t')
            for k in (needParent.keys()):
                if not k in alreadyGotParent:
                    if '('+str(k)+')' in splitLine[1]:
                        myCheck = (splitLine[1].split('('+str(k)+')')[0]).split()
                        myPlace = 0
                        while int(k) not in alreadyGotParent:
                            myPlace -= 1
                            if '(' in myCheck[myPlace]:
                                alreadyGotParent[int(k)] = int(myCheck[myPlace][1:-1])
                                myNodes[int(myCheck[myPlace][1:-1])] = True
            if lineCounter % 25000 == 0:
                print(lineCounter)

    myOutString = 'node\tparent\n'
    for k in alreadyGotParent:
        myOutString += str(k)+'\t'+str(alreadyGotParent[k])+'\n'
    open('nodeToParent.txt','w').write(myOutString)

    myOutString = ''
    for k in sorted(myNodes.keys()):
        myOutString += str(k)+'\n'
    open('allRelevantNodes.txt','w').write(myOutString)

##########################
#### HELPER FUNCTIONS ####
##########################

def getPos(myInds, intLineNumToPos):
    myReturn = []
    for k in myInds:
        myReturn.append(intLineNumToPos[k])
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return('\t'.join(newList))

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return(','.join(newList))

#########################
##### FUNCTION CALL #####
#########################

def main():
    getAllNodes()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
