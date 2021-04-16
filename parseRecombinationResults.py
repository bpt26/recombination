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


##########################
###### COMMAND LINE ######
##########################

class CommandLine(object):
    """Handles the input arguments from the command line. Manages 
    the argument parser."""

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line input using argparse.
        '''
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-f", "--fasta", help="Fasta file containing sequences for recombinant genomes. [REQUIRED]",default='')
        self.parser.add_argument("-l", "--log", help="Log file containing recombinant genomes. Format: recombinantSample sample1 sample2 bp1 (bp2)... [REQUIRED]",default='')
        self.parser.add_argument('-d', '--descendants', help="Descendants file output by findRecombination. (Default = 'descendants.tsv').",default='descendants.tsv')
        self.parser.add_argument('-r', '--recombination', help="Recombination file output by findRecombination. (Default = 'recombination.tsv').",default='recombination.tsv')
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()        

##########################
##### MAIN FUNCTIONS #####
##########################

def parseDescRec(myF, myL, myD, myR):
    mySamples = {}
    with open(myF) as f:
        for line in f:
            l = line.strip()
            if l.startswith('>RECOMB'):
                mySamples[l[1:]] = True

    sampleToBreakPoints = {}
    with open(myL) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if splitLine[0]+'_X0' in mySamples:
                sampleToBreakPoints[splitLine[0]+'_X0'] = int(splitLine[3])

    sampleToRecNodeIds = {}
    with open(myD) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                mySamples = splitLine[1].split(',')
                for k in mySamples:
                    if k in sampleToBreakPoints:
                        if not k in sampleToRecNodeIds:
                            sampleToRecNodeIds[k] = []
                        sampleToRecNodeIds[k].append(int(splitLine[0]))

    nodeIdToBreakPoints = {}
    with open(myR) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in nodeIdToBreakPoints:
                    nodeIdToBreakPoints[int(splitLine[0])] = []
                if 'GENOME_SIZE)' in splitLine[2]:
                    splitLine[2] = splitLine[2].replace('GENOME_SIZE)', '29903)')
                tempList = []
                for k in (splitLine[1][1:-1]).split(','):
                    tempList.append(int(k))
                for k in (splitLine[2][1:-1]).split(','):
                    tempList.append(int(k))
                nodeIdToBreakPoints[int(splitLine[0])].append(tempList)

    sampleToCorrectBreakPoints = {}
    for s in sampleToBreakPoints:
        sampleToCorrectBreakPoints[s] = False
        bp = sampleToBreakPoints[s]
        for n in sampleToRecNodeIds[s]:
            for l in nodeIdToBreakPoints[n]:
                if (bp > l[0] and bp < l[1]) or (bp > l[2] and bp < l[3]):
                    sampleToCorrectBreakPoints[s] = True

    tf = [0,0]
    for k in sampleToCorrectBreakPoints:
        if sampleToCorrectBreakPoints[k] == True:
            tf[0] += 1
        else:
            tf[1] += 1
            print(k)
    print(tf)


def parseDescRec2(myF, myL, myD, myR):
    mySamples = {}
    with open(myF) as f:
        for line in f:
            l = line.strip()
            if l.startswith('>RECOMB'):
                mySamples[l[1:]] = True

    sampleToBreakPoints = {}
    with open(myL) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if splitLine[0]+'_X0' in mySamples:
                sampleToBreakPoints[splitLine[0]+'_X0'] = [int(splitLine[3]), int(splitLine[4])]

    sampleToRecNodeIds = {}
    with open(myD) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                mySamples = splitLine[1].split(',')
                for k in mySamples:
                    if k in sampleToBreakPoints:
                        if not k in sampleToRecNodeIds:
                            sampleToRecNodeIds[k] = []
                        sampleToRecNodeIds[k].append(int(splitLine[0]))

    nodeIdToBreakPoints = {}
    with open(myR) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in nodeIdToBreakPoints:
                    nodeIdToBreakPoints[int(splitLine[0])] = []
                if 'GENOME_SIZE)' in splitLine[2]:
                    splitLine[2] = splitLine[2].replace('GENOME_SIZE)', '29903)')
                tempList = []
                for k in (splitLine[1][1:-1]).split(','):
                    tempList.append(int(k))
                for k in (splitLine[2][1:-1]).split(','):
                    tempList.append(int(k))
                nodeIdToBreakPoints[int(splitLine[0])].append(tempList)

    sampleToCorrectBreakPoints = {}
    for s in sampleToBreakPoints:
        sampleToCorrectBreakPoints[s] = False
        bp1 = sampleToBreakPoints[s][0]
        bp2 = sampleToBreakPoints[s][1]
        for n in sampleToRecNodeIds[s]:
            for l in nodeIdToBreakPoints[n]:
                if (bp1 > l[0] and bp1 < l[1] and bp2 > l[2] and bp2 < l[3]):
                    sampleToCorrectBreakPoints[s] = True

    tf = [0,0]
    for k in sampleToCorrectBreakPoints:
        if sampleToCorrectBreakPoints[k] == True:
            tf[0] += 1
        else:
            tf[1] += 1
            print(k)
    print(tf)

##########################
#### HELPER FUNCTIONS ####
##########################

def replaceSymbols(myEntry):
    myEntry = myEntry.replace('|', '_')
    myEntry = myEntry.replace('/', '_')
    return(myEntry)

def toInt(entry):
    myReturn = []
    for k in entry:
        myReturn.append(int(k))
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def joinerU(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '_'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append('"'+str(k)+'"')
    return ','.join(newList)

#########################
##### FUNCTION CALL #####
#########################

def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    # Necessary files:
    if myCommandLine.args.fasta:
        myF = myCommandLine.args.fasta
    if myCommandLine.args.log:
        myL = myCommandLine.args.log
    if myCommandLine.args.descendants:
        myD = myCommandLine.args.descendants
    if myCommandLine.args.recombination:
        myR = myCommandLine.args.recombination

    with open(myL) as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myLen = len(splitLine)
    if myLen == 4:
        parseDescRec(myF, myL, myD, myR)
    else:
        parseDescRec2(myF, myL, myD, myR)

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

























