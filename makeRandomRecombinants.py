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
import argparse

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
        self.parser.add_argument("-b", "--breakpoints", help="Number of breakpoints that each recombinant sample will have. Must be 1 or 2 (Default = 1).", default=1, type=int)
        self.parser.add_argument("-s", "--samples", help="Number of recombinant samples to create (Default = 100).", default=100, type=int)
        self.parser.add_argument("-c", "--copies", help="Number of identical copies to make for each recombinant sample (Default = 10).", default=10, type=int)
        self.parser.add_argument("-d", "--differences", help="Minimum mutational distance for acceptor/donor samples (Default = 10).", default=10, type=int)
        self.parser.add_argument("-f", "--fasta", help="Fasta file containing sequences for acceptor/donor samples. [REQUIRED]", default='')
        self.parser.add_argument("-r", "--ref", help="Fasta file containing reference genome for use in creating VCF. (Default = 'wuhan.ref.fa').", default='wuhan.ref.fa')
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()     
        if self.args.breakpoints < 1 or self.args.breakpoints > 2:
            sys.stderr.write("Please retry with either 1 or 2 breakpoints.\n")
            sys.exit(1)

##########################
##### MAIN FUNCTIONS #####
##########################

def makeExamples(myS, myB, myC, myD, myF, myR):
    posToRef = {}
    with open(myR) as f:
        for line in f:
            l = line.strip().upper()
            if not l.startswith('>'):
                myReference = l
                for i in range(0,len(l)):
                    l[i] != 'N'
                    posToRef[i] = l[i]
            else:
                myRefName = l[1:]

    sampleToSeq = {}
    with open(myF) as f:
        for line in f:
            l = line.strip()
            if l.startswith('>'):
                mySample = l[1:]
            else:
                sampleToSeq[mySample] = l

    recSampleToLog = {}
    recSampleToSeq = {}
    recSampleToDiffBetweenBps = {}
    while len(recSampleToSeq) < myS:
        samples = numpy.random.choice(list(sampleToSeq.keys()), size=2, replace=False)
        mySampleName = 'RECOMB_'+str(myB)+'_'+str(r)+'_'+str(samples[0])+'_'+str(samples[1])
        s1 = samples[0]
        s2 = samples[1]
        if myB == 1:
            bps = numpy.random.choice(sorted(list(posToRef.keys()))[5000:-5000],size=1, replace=False)
            bp1 = bps[0]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:]
            myDiff = minLen(getDiff(sampleToSeq[s1][:bp1], sampleToSeq[s2][:bp1], 0), getDiff(sampleToSeq[s1][bp1:], sampleToSeq[s2][bp1:], bp1))
        elif myB == 2:
            bps = sorted(numpy.random.choice(sorted(list(posToRef.keys()))[5000:-5000],size=2, replace=False))
            if bps[1]-bps[0] <= 1000:
                bps = sorted(numpy.random.choice(sorted(list(posToRef.keys()))[5000:-5000],size=2, replace=False))
            bp1 = bps[0]
            bp2 = bps[1]
            mySeq = sampleToSeq[s1][:bp1]+sampleToSeq[s2][bp1:bp2]+sampleToSeq[s1][bp2:]
            myDiff =  minLen(getDiff(sampleToSeq[s1][bp1:bp2], sampleToSeq[s2][bp1:bp2], bp1), getDiff(sampleToSeq[s1][:bp1]+sampleToSeq[s1][bp2:], sampleToSeq[s2][:bp1]+sampleToSeq[s2][bp2:], 0))
        if myDiff >= myD:
            recSampleToLog[mySampleName] = [samples, bps]
            recSampleToSeq[mySampleName] = mySeq
            recSampleToDiffBetweenBps[mySampleName] = myDiff

    myOutMSA = '>'+myRefName+'\n'+myReference+'\n'
    myOutLog = ''
    myOutDiff = ''

    for s in recSampleToSeq:
        for x in range(0,myC):
            myOutMSA += '>'+s+'_X'+str(x)+'\n'+recSampleToSeq[s]+'\n'
        myOutLog += s+'\t'+doubleJoiner(recSampleToLog[s])+'\n'
        myOutDiff += s+'\t'+joinerC(recSampleToDiffBetweenBps[s])+'\n'
    open('recombination_'+str(myB)+'.msa.fa','w').write(myOutMSA)
    open('recombination_'+str(myB)+'.log','w').write(myOutLog)
    open('recombination_'+str(myB)+'.differences.txt','w').write(myOutDiff)

##########################
#### HELPER FUNCTIONS ####
##########################

def addMut(seq, pos):
    print(pos, len(seq))
    myReturn = []
    for i in range(0,len(seq)):
        if i != int(pos):
            myReturn.append(seq[i])
        else:
            if seq[i] == 'A':
                mySub = numpy.random.choice(['C','G','T'], size=1)[0]
            elif seq[i] == 'C':
                mySub = numpy.random.choice(['A','G','T'], size=1)[0]
            elif seq[i] == 'G':
                mySub = numpy.random.choice(['C','A','T'], size=1)[0]
            elif seq[i] == 'T':
                mySub = numpy.random.choice(['A','G','C'], size=1)[0]
            else:
                mySub = numpy.random.choice(['A','G','C','T'], size=1)[0]
            myMut = seq[i]+str(pos)+mySub
            myReturn.append(mySub)
    return(''.join(myReturn), myMut)

def getDiff(s1, s2, add):
    myReturn = []
    for i in range(0,len(s1)):
        if s1[i] != s2[i]:
            myReturn.append(add+i)
    return(myReturn)

def minLen(l1, l2):
    if len(l1) <= len(l2):
        return(l1)
    else:
        return(l2)

def doubleJoiner(myList):
    myReturn = []
    for k in myList:
        myReturn.append(joiner(k))
    return(joiner(myReturn))

def replaceSymbols(myEntry):
    myEntry = myEntry.replace('|', '_')
    myEntry = myEntry.replace('/', '_')
    return(myEntry)

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
        newList.append(str(k))
    return(','.join(newList))

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
    if myCommandLine.args.samples:
        myS = myCommandLine.args.samples
    if myCommandLine.args.breakpoints:
        myB = myCommandLine.args.breakpoints
    if myCommandLine.args.copies:
        myC = myCommandLine.args.copies
    if myCommandLine.args.differences:
        myD = myCommandLine.args.differences
    if myCommandLine.args.fasta:
        myF = myCommandLine.args.fasta
    if myCommandLine.args.ref:
        myR = myCommandLine.args.ref

    makeExamples(myS, myB, myC, myF, myR)

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit























if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit






















