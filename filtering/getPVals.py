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


def getPVals():
    nodeToDesc = {}
    with open('catDescendants.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                nodeToDesc[int(splitLine[0])] = splitLine[1]

    russNull = {}
    russOrigParsToTotal = {}
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

    robNull = {}
    robOrigParsToTotal = {}
    with open('rob_null.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            print(splitLine)
            if len(splitLine) == 1 and splitLine[0].isdigit():
                myOrigPars = int(splitLine[0])
                robNull[myOrigPars] = {}
                robOrigParsToTotal[myOrigPars] = 0
            elif len(splitLine) == 2:
                (robNull[myOrigPars])[int(splitLine[0])] = int(splitLine[1])
                robOrigParsToTotal[myOrigPars] += int(splitLine[1])

    print(robOrigParsToTotal)
    print(robNull)

    myOutString = '#recomb_node_id\tbreakpoint-1_interval\tbreakpoint-2_interval\tdonor_node_id\tdonor_is_sibling\tdonor_parsimony\t'
    myOutString += 'acceptor_node_id\tacceptor_is_sibling\tacceptor_parsimony\toriginal_parsimony\tmin_starting_parsimony\trecomb_parsimony'
    myOutString += '\trob_pval\truss_pval\tdescendants\n'
    with open('combinedCatOnlyBest.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if int(splitLine[-2]) > 0 and (int(splitLine[-2])-int(splitLine[-1])) >= 3:

                if int(splitLine[-2]) not in robNull:
                    myRobNull = 'NA'
                else:
                    myTotal = robOrigParsToTotal[int(splitLine[-2])]
                    myImprovement = int(splitLine[-2])-int(splitLine[-1])
                    for k in sorted(robNull[(int(splitLine[-2]))].keys()):
                        if k < myImprovement:
                            myTotal -= robNull[(int(splitLine[-2]))][k]
                    if myTotal == 0:
                        myRobNull = '0/'+str(robOrigParsToTotal[int(splitLine[-2])])
                    else:
                        myRobNull = float(myTotal)/float(robOrigParsToTotal[int(splitLine[-2])])


                if int(splitLine[-2]) not in russNull:
                    myRussNull = 'NA'
                else:
                    myTotal = russOrigParsToTotal[int(splitLine[-2])]
                    myImprovement = int(splitLine[-2])-int(splitLine[-1])
                    for k in sorted(russNull[(int(splitLine[-2]))].keys()):
                        if k < myImprovement:
                            myTotal -= russNull[(int(splitLine[-2]))][k]
                    if myTotal == 0:
                        myRussNull = '0/'+str(russOrigParsToTotal[int(splitLine[-2])])
                    else:
                        myRussNull = float(myTotal)/float(russOrigParsToTotal[int(splitLine[-2])])
                splitLine.append(myRobNull)
                splitLine.append(myRussNull)
                splitLine.append(nodeToDesc[int(splitLine[0])])
                myOutString += joiner(splitLine)+'\n'
    print("Writing with Pvals")
    open('combinedCatOnlyBestWithPVals.txt','w').write(myOutString)


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
    getPVals()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit




