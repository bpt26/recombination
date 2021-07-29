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

def makeMNK():
    myOutString = ''
    with open('allRelevantNodesInfSeq.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            seq = splitLine[-1]
            if seq.startswith('A'):
                myOutString += joiner(splitLine)+'\t'+joiner([seq.count('A'),seq.count('B'),getK(seq,'A','B')])+'\n'
            else:
                myOutString += joiner(splitLine)+'\t'+joiner([seq.count('B'),seq.count('A'),getK(seq,'B','A')])+'\n'
    open('allRelevantNodesMNK.txt','w').write(myOutString)

def removeDups():
    alreadyDone = {}
    myOutString = ''
    with open('allRelevantNodesMNK.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not str(splitLine[-3])+'_'+str(splitLine[-2])+'_'+str(splitLine[-1]) in alreadyDone:
                myOutString += str(splitLine[-3])+' '+str(splitLine[-2])+' '+str(splitLine[-1])+'\n'
            alreadyDone[str(splitLine[-3])+'_'+str(splitLine[-2])+'_'+str(splitLine[-1])] = True
    open('mnk_no_dups.txt','w').write(myOutString)

def addPVals():
    keyToP = {}
    with open('mnk.log') as f:
        for line in f:
            splitLine = (line.strip()).split()
            if len(splitLine) > 1:
                if splitLine[0].startswith('Enter'):
                    myKey = '_'.join(splitLine[-3:])
                else:
                    keyToP[myKey] = float(splitLine[-1])

    myOutString = ''
    alreadyUsed = {}
    with open('allRelevantNodesMNK.txt') as f:
        for line in f:
            splitLine = (line.strip()).split()
            myKey = '_'.join(splitLine[-3:])
            if not myKey in keyToP:
                print(myKey)
            else:
                myOutString += joiner(splitLine)+'\t'+str(keyToP[myKey])+'\n'
    open('allRelevantNodesMNKPval.txt','w').write(myOutString)


def combinePValueFiles():
    recombToParents = {}
    with open('combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in recombToParents:
                    recombToParents[int(splitLine[0])] = {}
                recombToParents[int(splitLine[0])][str(splitLine[3])+'_'+str(splitLine[6])] = True

    recombTo3seqPval = {}
    recombToBestParents = {}
    recombToAB = {}
    recombToPrinted = {}
    with open('allRelevantNodesMNKPval.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if int(splitLine[0]) in recombToParents:
                    if (str(splitLine[1])+'_'+str(splitLine[2]) in recombToParents[int(splitLine[0])]):
                        if int(splitLine[0]) not in recombTo3seqPval or float(splitLine[11]) < recombTo3seqPval[int(splitLine[0])]:
                            recombTo3seqPval[int(splitLine[0])] = float(splitLine[11])
                            recombToBestParents[int(splitLine[0])] = str(splitLine[1])+'_'+str(splitLine[2])
                            recombToAB[int(splitLine[0])] = splitLine[7]
                            recombToPrinted[int(splitLine[0])] = False
                    #else:
                    #    print(splitLine[:3])

    #print(recombToParents)

    myOutString = ''
    myOutString2 = ''
    with open('combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if int(splitLine[0]) in recombTo3seqPval:
                    if recombToBestParents[int(splitLine[0])] == str(splitLine[3])+'_'+str(splitLine[6]):
                        myOutString += '\t'.join(splitLine[:-1])+'\t'+str(recombTo3seqPval[int(splitLine[0])])+'\t'+recombToAB[int(splitLine[0])]+'\t'+str(splitLine[-1])+'\n'
                        if splitLine[12].startswith('0/'):
                            splitLine[12] = '0.0'
                        if splitLine[13].startswith('0/'):
                            splitLine[13] = '0.0'
                        myOutString2 += '\t'.join(splitLine[:-1])+'\t'+str(recombTo3seqPval[int(splitLine[0])])+'\t'+recombToAB[int(splitLine[0])]+'\t'+str(splitLine[-1])+'\n'
                        recombToPrinted[int(splitLine[0])] = True
                #else:
                #    print(splitLine[0])
    open('combinedCatOnlyBestWithAll3PValsTiesBroken.txt','w').write(myOutString)
    open('combinedCatOnlyBestWithAll3PValsRealTiesBroken.txt','w').write(myOutString2)

    for k in recombToPrinted:
        if recombToPrinted[k] == False:
            print(k, recombToBestParents[k])

def addInfSites():
    trioToInfSites = {}
    with open('allRelevantNodesInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioToInfSites[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[7]

    trioToInfSeq = {}
    with open('allRelevantNodesInfSeq.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioToInfSeq[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[7]

    myOutString = ''
    with open('combinedCatOnlyBestWithPVals.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6]) in trioToInfSites:
                splitLine.append(trioToInfSites[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])])
                splitLine.append(trioToInfSeq[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])])
                myOutString += joiner(splitLine)+'\n'
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSites.txt','w').write(myOutString)


def checkClusters():
    trioTo3P = {}
    with open('allRelevantNodesMNKPval.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            trioTo3P[str(splitLine[0])+'_'+str(splitLine[1])+'_'+str(splitLine[2])] = splitLine[8:]

    myOutString = ''
    bp1 = {}
    bp2 = {}
    #myOKs = {28881:True,28882:True,28883:True,28280:True,28281:True,28282:True}
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSites.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myTrio = str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])
            myStart = int(splitLine[1].split(',')[0][1:])
            myEnd = int(splitLine[2].split(',')[1][:-1])
            mySeq = splitLine[16]
            mySites = toInt(splitLine[15].split(','))
            myA = []
            myB = []
            for i in range(0,len(mySites)):
                if mySeq[i] == 'A':
                    myA.append(mySites[i])
                elif mySeq[i] == 'B':
                    myB.append(mySites[i])
            #print(max(myA), min(myA), max(myB), min(myB))
            if max(myA)-min(myA) > 20 and max(myB)-min(myB) > 20:
                myOutString += line.strip()+'\t'+joiner(trioTo3P[myTrio])+'\n'
                if myStart == 0 or myEnd == 29903:
                    bp1[splitLine[0]] = True
                else:
                    bp2[splitLine[0]] = True
            else:
                if int(sys.argv[1]) == 1:
                    keep = False
                    for k in myOKs:
                        if k in myA or k in myB:
                            keep = True
                    if keep == True:
                        myOutString += line.strip()+'\t'+joiner(trioTo3P[myTrio])+'\n'
                        if myStart == 0 or myEnd == 29903:
                            bp1[splitLine[0]] = True
                        else:
                            bp2[splitLine[0]] = True
    print(len(bp1),len(bp2))
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt','w').write(myOutString)


def applyPval():
    myOutString = ''
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if float(splitLine[20]) <= 0.2:
                if splitLine[13].startswith('0/') or splitLine[13].startswith('NA') or float(splitLine[13]) < 0.05:
                    myOutString += joiner(splitLine)+'\n'
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3seqP02RussPval005.txt','w').write(myOutString)


def doNewTiebreakers():
    nodeToLeaves = {}
    with gzip.open('optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz') as f:
        for line in f:
            splitLine = (line.decode('utf8').strip()).split('\t')
            if splitLine[0].isdigit():
                nodeToLeaves[int(splitLine[0])] = int(splitLine[1])

    bp1 = {}
    bp2 = {}
    recombToBPs = {}
    recombToStringSize = {}
    recombToLines = {}
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClusters3seqP02RussPval005.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            if not splitLine[0].startswith('#'):
                if not int(splitLine[0]) in recombToBPs:
                    recombToBPs[int(splitLine[0])] = []
                    recombToStringSize[int(splitLine[0])] = []
                    recombToLines[int(splitLine[0])] = []
                tempPreStart = 0
                tempStart = 0
                tempMid = 0
                tempEnd = 0
                tempPostEnd = 0
                myInfSites = toInt(splitLine[15].split(','))
                myStart1 = int(splitLine[1].split(',')[0][1:])
                myStart2 = int(splitLine[1].split(',')[1][:-1])
                myEnd1 = int(splitLine[2].split(',')[1][:-1])
                myEnd2 = int(splitLine[2].split(',')[0][1:])
                for k in myInfSites:
                    if k <= myStart1:
                        tempPreStart += 1
                    elif k >= myStart1 and k <= myStart2:
                        tempStart += 1
                    elif k > myStart2 and k < myEnd1:
                        tempMid += 1
                    elif k >= myEnd1 and k <= myEnd2:
                        tempEnd += 1
                    elif k > myEnd2:
                        tempPostEnd += 1
                if tempPreStart == 0 or tempPostEnd == 0:
                    recombToBPs[int(splitLine[0])].append(1)
                    recombToStringSize[int(splitLine[0])].append(len(myInfSites))
                    recombToLines[int(splitLine[0])].append(splitLine)
                else:
                    recombToBPs[int(splitLine[0])].append(2)
                    recombToStringSize[int(splitLine[0])].append(len(myInfSites))
                    recombToLines[int(splitLine[0])].append(splitLine)

    myOutString = ''
    for k in recombToBPs:
        tempB = []
        tempS = []
        tempL = []
        for i in range(0,len(recombToBPs[k])):
            tempB.append(recombToBPs[k][i])
            tempS.append(recombToStringSize[k][i])
            tempL.append(recombToLines[k][i])

        if len(tempB) == 1: # if only one, just print it
            myOutString += joiner(tempL[0])+'\n'
        else:
            if tempB.count(1) == 1:
                myOutString += joiner(tempL[tempB.index(1)])+'\n' # if one best, print that

            else: # else, get all tied and go to next tiebreaker: smallest 3seq string
                if tempB.count(1) == 0:
                    newB = tempB
                    newS = tempS
                    newL = tempL
                elif tempB.count(1) > 1:
                    newB = []
                    newS = []
                    newL = []
                    for i in range(0,len(tempB)):
                        if tempB[i] == 1:
                            newB.append(tempB[i])
                            newS.append(tempS[i])
                            newL.append(tempL[i])

                minS = min(newS)
                if newS.count(minS) == 1:
                    myOutString += joiner(newL[newS.index(minS)])+'\n'
                elif newS.count(minS) > 1:
                    tempB = []
                    tempS = []
                    tempL = []
                    tempLeaves = []
                    tempP = []
                    tempPval = []
                    for i in range(0,len(newS)):
                        if newS[i] == minS:
                            tempB.append(newB[i])
                            tempS.append(newS[i])
                            tempL.append(newL[i])
                            tempLeaves.append(nodeToLeaves[int(newL[i][3])]+nodeToLeaves[int(newL[i][6])])
                            tempP.append(set([int(newL[i][3]),int(newL[i][6])]))
                            tempPval.append(float(newL[i][-1]))

                    minPval = min(tempPval)
                    if tempPval.count(minPval) == 1:
                        myOutString += joiner(tempL[tempPval.index(minPval)])+'\n'
                    elif tempPval.count(minPval) > 1:
                        if getNumUnique(tempP) == 1: # if all have the same parents, print biggest breakpoint interval
                            myOutString += joiner(getBiggestBreakpointInterval(tempL))+'\n'
                        else:
                            print(tempP)
                            newB = []
                            newS = []
                            newL = []
                            newLeaves = []
                            newP = []
                            newPval = []
                            for i in range(0,len(tempPval)):
                                if tempPval[i] == minPval:
                                    newB.append(tempB[i])
                                    newS.append(tempS[i])
                                    newL.append(tempL[i])
                                    newLeaves.append(nodeToLeaves[int(tempL[i][3])]+nodeToLeaves[int(tempL[i][6])])
                                    newP.append(set([int(tempL[i][3]),int(tempL[i][6])]))

                            minLeaves = min(newLeaves)
                            if newLeaves.count(minLeaves) == 1:
                                myOutString += joiner(newL[newLeaves.index(minLeaves)])+'\n'
                            elif newLeaves.count(minLeaves) > 1:
                                tempB = []
                                tempS = []
                                tempL = []
                                tempLeaves = []
                                tempP = []
                                tempPval = []
                                for i in range(0,len(newLeaves)):
                                    if newLeaves[i] == minLeaves:
                                        tempB.append(newB[i])
                                        tempS.append(newS[i])
                                        tempL.append(newL[i])
                                        tempLeaves.append(nodeToLeaves[int(newL[i][3])]+nodeToLeaves[int(newL[i][6])])
                                        tempP.append(set([int(newL[i][3]),int(newL[i][6])]))
                                myOutString += joiner(getBiggestBreakpointInterval(tempL))
    open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005.txt','w').write(myOutString)


def removeRedundantTrios():
    nodeToLeaves = {}
    with gzip.open('optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.leaves.txt.gz') as f:
        for line in f:
            splitLine = (line.decode('utf8').strip()).split('\t')
            if splitLine[0].isdigit():
                nodeToLeaves[int(splitLine[0])] = int(splitLine[1])

    myTrios = []
    trioToPVal = {}
    trioToSites = {}
    trioToLeaves = {}
    trioToLine = {}
    lc = 0
    with open('combinedCatOnlyBestWithPValsFinalReportWithInfSitesNoClustersNewTiebreak3seqP02RussPval005.txt') as f:
        for line in f:
            splitLine = (line.strip()).split('\t')
            myTrios.append([int(splitLine[0]),int(splitLine[3]),int(splitLine[6])])
            trioToLine[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = splitLine
            if splitLine[13].startswith('0/'):
                splitLine[13] = (1.0/float(splitLine[13][2:]))
            trioToPVal[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = float(splitLine[13])
            trioToSites[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = len(splitLine[16])
            trioToLeaves[str(splitLine[0])+'_'+str(splitLine[3])+'_'+str(splitLine[6])] = nodeToLeaves[int(splitLine[3])]+nodeToLeaves[int(splitLine[6])]
            lc += 1

    toRemove = {}
    for i in range(0,len(myTrios)):
        for j in range(i+1,len(myTrios)):
            if checkTwo(myTrios[i],myTrios[j]) == True: # if circular logic:
                if trioToPVal[joinerU(myTrios[i])] < trioToPVal[joinerU(myTrios[j])]: # remove case with lower pval
                    toRemove[joinerU(myTrios[j])] = True
                elif trioToPVal[joinerU(myTrios[i])] > trioToPVal[joinerU(myTrios[j])]:
                    toRemove[joinerU(myTrios[i])] = True
                elif trioToPVal[joinerU(myTrios[i])] == trioToPVal[joinerU(myTrios[j])]:

                    if trioToSites[joinerU(myTrios[i])] > trioToSites[joinerU(myTrios[j])]: # remove case with more informative sites
                        toRemove[joinerU(myTrios[i])] = True
                    elif trioToSites[joinerU(myTrios[i])] < trioToSites[joinerU(myTrios[j])]:
                        toRemove[joinerU(myTrios[j])] = True
                    elif trioToSites[joinerU(myTrios[i])] == trioToSites[joinerU(myTrios[j])]:

                        if trioToLeaves[joinerU(myTrios[i])] < trioToLeaves[joinerU(myTrios[j])]:
                             toRemove[joinerU(myTrios[i])] = True
                        if trioToLeaves[joinerU(myTrios[i])] > trioToLeaves[joinerU(myTrios[j])]:
                             toRemove[joinerU(myTrios[j])] = True
                        elif trioToLeaves[joinerU(myTrios[i])] == trioToLeaves[joinerU(myTrios[j])]:
                            print(myTrios[i], myTrios[j], trioToPVal[joinerU(myTrios[i])])

    myOutString = ''
    for t in trioToLine:
        if not t in toRemove:
            myOutString += joiner(trioToLine[t])+'\n'
    open('finalRecombNodesSet.txt','w').write(myOutString)



##########################
#### HELPER FUNCTIONS ####
##########################

def checkTwo(list1, list2):
    mySame = 0
    for k in list1:
        if k in list2:
            mySame += 1
    if mySame == 3:
        return(True)
    elif mySame == 2:
        if list1[0] not in list2 and list2[0] not in list1:
            return(False)
        else:
            return(True)
    else:
        return(False)

def getK(seq, a, b):
    myPath = []
    currentPlace = 0
    for k in seq:
        if k == a:
            currentPlace += 1
        else:
            currentPlace -= 1
        myPath.append(currentPlace)
    maxDesc = 0
    for i in range(1,len(myPath)):
        if max(myPath[:i])-myPath[i] > maxDesc:
            maxDesc = max(myPath[:i])-myPath[i]
    return(maxDesc)


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


def toInt(myList):
    myReturn = []
    for k in myList:
        myReturn.append(int(k))
    return(myReturn)

def getBiggestBreakpointInterval(myList):
    temp = []
    for k in myList:
        temp.append(numpy.sum(toInt(k[1][1:-1].split(',')))+numpy.sum(toInt(k[2][1:-1].split(','))))
    return(myList[temp.index(max(temp))])


def getNumUnique(myList):
    myReturn = 1
    for i in range(0,len(myList)):
        for j in range(i+1,len(myList)):
            if myList[i] != myList[j]:
                #print(myList[i], myList[j])
                myReturn += 1
    return(myReturn)

#########################
##### FUNCTION CALL #####
#########################

def main():
    makeMNK()
    removeDups()
    addPVals()
    combinePValueFiles()
    addInfSites()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit












