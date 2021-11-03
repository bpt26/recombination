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

"""
- for each recombinant node, get the sites where it matches one parent but not the other
- however, if a recombinant node does not match a parent, and that parent was placed as a sibling,
check to see if the grandparent on that side matches the recombinant node. if it does, act as if it matches the parent.
"""

##########################
##### MAIN FUNCTIONS #####
##########################

def getABABA():
    recombToParents = {}
    recombToEndRow = {}
    recombToParentSib = {}
    with open('combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in recombToParents:
                    recombToParents[int(splitLine[0])] = []
                    recombToEndRow[int(splitLine[0])] = []
                recombToParents[int(splitLine[0])].append([int(splitLine[3]),int(splitLine[6])])
                recombToEndRow[int(splitLine[0])].append(splitLine[10:-1])
                if splitLine[4] == 'y':
                    if not int(splitLine[0]) in recombToParentSib:
                        recombToParentSib[(int(splitLine[0]))] = {}
                    recombToParentSib[(int(splitLine[0]))][int(splitLine[3])] = True
                if splitLine[7] == 'y':
                    if not int(splitLine[0]) in recombToParentSib:
                        recombToParentSib[(int(splitLine[0]))] = {}
                    recombToParentSib[(int(splitLine[0]))][int(splitLine[6])] = True

    parentToGrand = {}
    with open('nodeToParent.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0] == 'node':
                parentToGrand[int(splitLine[0])] = int(splitLine[1])

    nodeToIndex = {}
    indexToNode = {}
    myRecombIndices = {}
    recombToInformativeSites = {}
    recombToInformativeSeq = {}
    recombToSiteChanges = {}
    with open('allRelevantNodes.vcf') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if splitLine[0] == '#CHROM':
                for i in range(9,len(splitLine)):
                    if 'node_' in splitLine[i]:
                        splitLine[i] = splitLine[i].split('_')[-1]
                    indexToNode[i] = int(splitLine[i])
                    nodeToIndex[int(splitLine[i])] = i
                    if int(splitLine[i]) in recombToParents:
                        myRecombIndices[i] = []
                        recombToInformativeSites[i] = []
                        recombToInformativeSeq[i] = []
                        recombToSiteChanges[i] = []
                for i in myRecombIndices:
                    for p in recombToParents[int(indexToNode[i])]:
                        myRecombIndices[i].append([nodeToIndex[p[0]], nodeToIndex[p[1]]])
                        recombToInformativeSeq[i].append([])
                        recombToInformativeSites[i].append([])
                        recombToSiteChanges[i].append([])

            elif splitLine[0] == 'NC_045512v2':
                myGTs = [splitLine[3]+str(splitLine[1])+splitLine[3]]+splitLine[2].split(',')
                for i in myRecombIndices.keys():
                    for p in range(0,len(myRecombIndices[i])):
                        parentInds = (myRecombIndices[i])[p]
                        parentInd0 = parentInds[0]
                        parentInd1 = parentInds[1]

                        """
                        IF IT IS PLACED AS A SIBLING, INCLUDE ALL MUTATIONS THAT ARE ON THE BRANCH LEADING TO THE PARENT NODE THAT MATCH THE RECOMBINANT,
                        NO MATTER WHAT THE VCF SAYS, regardless of if it's a match or not
                        """

                        if int(splitLine[i]) != (int(splitLine[parentInds[0]])): # if recomb does not match parent0, get grandparent and check for match:
                            if indexToNode[i] in recombToParentSib and indexToNode[(int(parentInds[0]))] in recombToParentSib[indexToNode[i]]: # if non-matching parent is sib, get grandparent
                                grandParentInd = nodeToIndex[parentToGrand[indexToNode[parentInds[0]]]]
                                if int(splitLine[i]) == (int(splitLine[grandParentInd])): # if it matches recombinant, use grandparent
                                    parentInd0 = nodeToIndex[parentToGrand[indexToNode[parentInds[0]]]] # else, keep parentInd0 as default
                                    ((recombToSiteChanges[i])[p]).append(str(indexToNode[int(parentInds[0])])+':'+myGTs[int(splitLine[i])])


                        if int(splitLine[i]) != (int(splitLine[parentInds[1]])):
                            if indexToNode[i] in recombToParentSib and indexToNode[(int(parentInds[1]))] in recombToParentSib[indexToNode[i]]: # if non-matching parent is sib, get grandparent
                                grandParentInd = nodeToIndex[parentToGrand[indexToNode[parentInds[1]]]]
                                if int(splitLine[i]) == (int(splitLine[grandParentInd])): # if it matches recombinant, use grandparent
                                    parentInd1 = nodeToIndex[parentToGrand[indexToNode[parentInds[1]]]] # else, keep parentInd0 as default
                                    #print(i, p, recombToSiteChanges[i], myGTs, int(splitLine[i]))
                                    ((recombToSiteChanges[i])[p]).append(str(indexToNode[int(parentInds[1])])+':'+myGTs[int(splitLine[i])])

                        if int(splitLine[i]) == int(splitLine[parentInd0]) and int(splitLine[i]) != int(splitLine[parentInd1]):
                            ((recombToInformativeSeq[i])[p]).append('A')
                            ((recombToInformativeSites[i])[p]).append(int(splitLine[1]))
                        if int(splitLine[i]) == int(splitLine[parentInd1]) and int(splitLine[i]) != int(splitLine[parentInd0]):
                            ((recombToInformativeSeq[i])[p]).append('B')
                            ((recombToInformativeSites[i])[p]).append(int(splitLine[1]))
    myOutInfSeq = ''
    myOutInfSites = ''
    myOutSiteChanges = ''
    for recInd in recombToInformativeSeq:
        myRecombNode = indexToNode[recInd]
        for i in range(0,len(recombToParents[myRecombNode])):
            myParents = recombToParents[myRecombNode][i]
            myOutInfSites += str(myRecombNode)+'\t'+joiner(myParents)+'\t'+joiner(recombToEndRow[myRecombNode][i])+'\t'+joinerC(recombToInformativeSites[recInd][i])+'\n'
            myOutInfSeq += str(myRecombNode)+'\t'+joiner(myParents)+'\t'+joiner(recombToEndRow[myRecombNode][i])+'\t'+''.join(recombToInformativeSeq[recInd][i])+'\n'
            myOutSiteChanges += str(myRecombNode)+'\t'+joiner(myParents)+'\t'+joiner(recombToEndRow[myRecombNode][i])+'\t'+joinerC(recombToSiteChanges[recInd][i])+'\n'
    open('allRelevantNodesInfSites.txt','w').write(myOutInfSites)
    open('allRelevantNodesInfSeq.txt','w').write(myOutInfSeq)
    open('allRelevantNodesSiteChanges.txt','w').write(myOutSiteChanges)



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
    getABABA()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit


