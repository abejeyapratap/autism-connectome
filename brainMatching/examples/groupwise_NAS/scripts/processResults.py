#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 15:11:23 2018

@author: yusuf
"""
import argparse
import numpy as np
from helpers import getOutlierIndices

# -r 5_combinedStrFuncMatchingScore/results/direct_none_edgesIncludeDiagZeroDiag_LinAss_0.res -o ./file.txt -s /home/yusuf/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -mt similarity -al group
#get parameters
parser = argparse.ArgumentParser(description='process raw matching data to produces similarity scores for subjects')
parser.add_argument('-r','--resultFile', help='file path to the results of the matching experiment', required=True)
parser.add_argument('-s','--subjectsList', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-o','--outputFile', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('-mt','--measureType', help='which measure to analyze: similarity score or matching accuracy', required=False,type=str,choices=['accuracy','similarity'],default='similarity')
parser.add_argument('-al','--analysisLevel', help='is this a subjectwise or a groupwise analysis?', required=True,type=str,choices=['subject','group'])
parser.add_argument('-rt','--relativeTo', help='calculate similarity relative to which group?', required=False,type=str,choices=['self','healthy'],default='healthy')
parser.add_argument('--noSymmetry', help='do not symmetrize distances', required=False,action='store_true',default=False)
#parser.add_argument('--discardOutliersFromHealthy', help='remove outliers from healthy control population in calculation of scores', required=False,action='store_true',default=False)


args = vars(parser.parse_args())
resultFilePath=args['resultFile']
subjectListPath=args['subjectsList']
outputFilePath=args['outputFile']
measureType=args['measureType']
analysisLevel=args['analysisLevel']
relativeTo=args['relativeTo']
noSymmetry=args['noSymmetry']
#discardOutliersFromHealthy=args['discardOutliersFromHealthy']

##############load subjet IDs for which we have a score#########################
with open(subjectListPath,"r") as f:
    subjectList =  f.read().splitlines()
numScans=len(subjectList)

patientIDs=[]
healthyIDs=[]
healthyOrder=[]
for order,ID in enumerate(subjectList):
    if (ID.split("_")[0][0]=='c'):
        healthyIDs.append(ID.split('_')[0])
        healthyOrder.append(order)
    elif (ID.split("_")[0][0]=='p') and (ID.split('_')[0] not in patientIDs):
        patientIDs.append(ID.split('_')[0])

patientLongitudinalOrder = np.full((len(patientIDs),3),-1,dtype=int)
### Get indices of healthy controls and patients(_s1,_s2,_s3)
for i in range(len(subjectList)):
    ID=subjectList[i].split("_")[0]
    if ID in patientIDs:
        timepoint=int(subjectList[i].split("_")[-1][-1])-1
        patientLongitudinalOrder[patientIDs.index(ID)][timepoint] = i

### Get indices of healthy controls and patients(_s1,_s2,_s3)
healthy=healthyOrder
patientS1=[]
patientS2=[]
patientS3=[]
for i in range(len(patientLongitudinalOrder)):
    if (patientLongitudinalOrder[i,0]!=-1):
        patientS1.append(patientLongitudinalOrder[i,0])
    if (patientLongitudinalOrder[i,1]!=-1):
        patientS2.append(patientLongitudinalOrder[i,1])
    if (patientLongitudinalOrder[i,2]!=-1):
        patientS3.append(patientLongitudinalOrder[i,2])


##############load results of the structure-function coupling experiment###################
with open(resultFilePath,"r") as resultFile:
    fileContent = resultFile.read().splitlines()
numNodes = int(fileContent[3].split('\t')[0])
numSubjects = int(fileContent[3].split('\t')[1])

###load similarity scores and the matchings
if(analysisLevel=="group"):
    if(measureType=='similarity'):
        # starting with 6th line, read numSubjects line for the pairwise similarity scores
        scores=np.array([np.fromstring(cont,dtype=float,sep='\t') for cont in fileContent[5:5+numSubjects]])
    elif(measureType=='accuracy'):
        #starting from 6+numSubjects line, read until the end of the file for matching nodes between subject pairs
        matchings=np.array([np.fromstring(cont,dtype=int,sep='\t') for cont in fileContent[6+numSubjects:-1]])
        scores=np.zeros((numSubjects,numSubjects),dtype=float)
        for i in range(len(matchings)):
            row=matchings[i][0] #first colulmn is the order number of the first subject
            col=matchings[i][1] #second colulmn is the order number of the second subject
            for j in range(numNodes):
                if(matchings[i][2+j]==j):
                    scores[row][col]+=1
        scores/=float(numNodes)
        scores*=100
    
    ### now, calculate average matching/similarity scores relative to healthy controls
    ### NOTE: we are discarding the matching of a healthy control subject to itself in 
    ###       calculation of each subjec't average score relative to healthy controls
#    print(healthyIDs)
#    print(healthy)
#    if(discardOutliersFromHealthy==True):
#        tempScores_avg=np.zeros(len(healthy))
#        for rowOrder in range(len(healthy)):
#            count=0
#            row=healthy[rowOrder]
#            for colOrder in range(len(healthy)):
#                col=healthy[colOrder]
#                if(row!=col):
#                    if(noSymmetry==True):
#                        tempScores_avg[row] += scores[row][col]
#                    else:
#                        tempScores_avg[row] += (scores[row][col]+scores[col][row])/2.0
#                    count+=1
#            tempScores_avg[row] /= float(count)
#        indices=getOutlierIndices(tempScores_avg)
#        scores=np.delete(scores,indices,0)
#        scores=np.delete(scores,indices,1)
#        for i in reversed(range(len(healthy))):
#            del healthy[i]
#            del healthyIDs[i]
#            
#    print(healthyIDs)
#    print(healthy)

    scores_avg = np.zeros(len(scores))
    if(relativeTo=="healthy"):
        for row in range(numSubjects):
            count=0
            for col in healthy:
                if(row!=col):
                    if(noSymmetry==True):
                        scores_avg[row] += scores[row][col]
                    else:
                        scores_avg[row] += (scores[row][col]+scores[col][row])/2.0
                    count+=1
            scores_avg[row] /= float(count)
    elif(relativeTo=="self"):
        groups=[healthyOrder,patientS1,patientS2,patientS3]
        for group in groups:
            for row in group:
                count = 0;
                for col in group:
                    if(row!=col):
                        if(noSymmetry==True):
                            scores_avg[row] += scores[row][col]
                        else:
                            scores_avg[row] += (scores[row][col]+scores[col][row])/2.0
                        count+=1
                scores_avg[row] /= float(count)
    scores=scores_avg
elif(analysisLevel=="subject"):
    if(measureType=='similarity'):
        #6th line of result file contains similarity scores between structure and function for each subject, 
        #i.e. similarity scores is a vector
        scores=np.fromstring(fileContent[5],dtype=float,sep='\t') 
    elif(measureType=='accuracy'):
        #8th line until the end of file, contains the matching nodes for str vs function for each subject
        #i.e., 8th line is the structure-function matching nodes of the first subject
        matchings=np.array([np.fromstring(cont,dtype=int,sep='\t') for cont in fileContent[7:-1]])
        
        scores=np.zeros(len(matchings),dtype=float)
        for i in range(len(matchings)):
            for j in range(numNodes):
                if(matchings[i][j]==j):
                    scores[i]+=1
        scores/=float(numNodes) 
        scores*=100

if(measureType=='similarity'):
    scoreName='dissimilarity'
elif(measureType=='accuracy'):
    scoreName='matching accuracy (%)'
    
outputFile=open(outputFilePath,'w')
outputFile.write("#numNodes,numSubjects\n%d\t%d\n" % (numNodes,numSubjects))
outputFile.write("#measureType\n"+measureType+"\n")
outputFile.write("#scoreName\n"+scoreName+"\n")
outputFile.write("#scores\n")
for i in range(len(scores)):
    outputFile.write("%0.4f\t" % scores[i])
outputFile.close()
