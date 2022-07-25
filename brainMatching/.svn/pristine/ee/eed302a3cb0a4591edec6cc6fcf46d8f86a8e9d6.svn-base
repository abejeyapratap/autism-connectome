#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 17:10:09 2017

@author: yusuf

sample: python pathCompare.py -res ../direct/positive/src/ -out ./
"""

import argparse
import os
import numpy as np
from helpers import fdr

########helper functions################
#this function reads the experiment information from the result file
#takes path to a result file @filePath as input
#and returns @numNodes and @numSubjects as output
def readExperimentSize(filePath):
    fileContent =  open(filePath,"r").read().splitlines()
    numNodes = int(fileContent[3].split('\t')[0])
    numSubjects = int(fileContent[3].split('\t')[1])
    return numNodes, numSubjects


#This function reads the average subjectwise matching matrices from a single file
#note that, to be able to use this function, you should have already averaged the subjectwise matching results using 
#calcAvgMatchingResults.py code
def readAverageMatchingsForAllSubjects(filePath,ignoreData):
    fileContent =  open(filePath,"r").read().splitlines()
    numNodes = int(fileContent[3].split('\t')[0])
    numSubjects = int(fileContent[3].split('\t')[1])
    rowOffset=4
    mat=np.empty((numSubjects,numNodes,numNodes),dtype=float)
    for subj in range(numSubjects):
        rowStart=rowOffset+1
        for row in range(numNodes):
            columns=fileContent[rowStart+row].split('\t')
            for col in range(len(columns)):
                mat[subj,row,col]=columns[col]
        rowOffset+=numNodes+1
    if(len(ignoreData)!=0):
        mat = np.delete(mat,ignoreData,axis=0)
    return mat

#This function reads the average permutation test matching matrices from a single file
#note that, to be able to use this function, you should have already averaged the permutationTest matching results using 
#calcAvgMatchingResults.py code
def readAverageMatchingsForAllExperiments(filePath):
    fileContent =  open(filePath,"r").read().splitlines()
    numNodes = int(fileContent[3].split('\t')[0])
    numExperiments = int(fileContent[3].split('\t')[1])
    rowOffset=4
    mat=np.empty((numExperiments,numNodes,numNodes),dtype=float)
    for subj in range(numExperiments):
        rowStart=rowOffset+1
        for row in range(numNodes):
            columns=fileContent[rowStart+row].split('\t')
            for col in range(len(columns)):
                mat[subj,row,col]=columns[col]
        rowOffset+=numNodes+1
   
    return mat


#given a matrix that contains matching matrices for several subjects (numSubjects x numNodes x numNodes),
#this function takes their average and returns a single structure function matching matrix
#which is the average across all subjects (numNodes x numNodes)
def averageMatchingMatrixAcrossSubjects(matchingMatrices):
    numNodes = len(matchingMatrices[0])
    numSubjects = len(matchingMatrices)
    avgMatchingMatrix = np.zeros(shape=(numNodes,numNodes),dtype=float)
    for i in range(numSubjects):
        avgMatchingMatrix += matchingMatrices[i]
    avgMatchingMatrix /= numSubjects
    return avgMatchingMatrix

#given @subjectwiseAverageMatchingMatrix of size (numNodes,numNodes) containing the average connectivity matrix for actual subjectwise matching
#and @permutationMatchingMatricesPerExperiment of size (numExperiments,numNodes,numNodes) containing subjectwise average connectivity matrix for shuffled matching
#this function calculates the pValue for each of the (numNodes,numNodes) cells of the subjectwise matrix comparing it to the numExperiments results obtained from permutation test
def calculatePValueForPermutationTest(subjectwiseAverageMatchingMatrix,permutationMatchingMatricesPerExperiment,systemSize,connectomeSize):
    numNodes = len(subjectwiseAverageMatchingMatrix)
    numExperiments = len(permutationMatchingMatricesPerExperiment)
    pValueMatrix = np.zeros(shape=(numNodes,numNodes),dtype=float)
    
    for i in range(numExperiments):
        for j in range(numNodes):
            for k in range(numNodes):
                #threshold = max(permutationMatchingMatricesPerExperiment[i][j].max(),float(systemSize[k])/float(len(systemSize))) # for row wise
                threshold = max(permutationMatchingMatricesPerExperiment[i][j][k],float(systemSize[k])/float(connectomeSize)) #for pairwise
                if threshold >= subjectwiseAverageMatchingMatrix[j][k]:
                    pValueMatrix[j][k] += 1
    pValueMatrix /= float(numExperiments)
    return pValueMatrix


#given a matching matrix and a mapping array that maps the nodes of one connectome 
#to the nodes of another (most likely smaller) connectome
#this function converts the matching matrix to the new connectome and returns new connectome matrix
def reduceMatchingMatrix(matchingMatrix,reductionList,systemSize,scaleToSize=True):
    size = len(matchingMatrix)
    newSize = len(systemSize)
    if size!=len(reductionList):
        exit
    reducedMatrix = np.zeros(shape=(newSize,newSize),dtype=float)
    for i in range(size):
        for j in range(size):
            reducedMatrix[reductionList[i]-1][reductionList[j]-1] += matchingMatrix[i][j]
    for i in range(newSize):
        for j in range(newSize):
            if(systemSize[i]==0 or systemSize[j]==0):#if one of the systems have no ROIs assigned to it, reduced matrix for the combination of these two systems should be zero as well
                reducedMatrix[i,j] = 0
                if(reducedMatrix[i,j]!=0):
                    print("There is something wrong with the reduction of matrices... Check reduceMatchingMatrix()\n")
            else:
                if(scaleToSize==True):
                    reducedMatrix[i,j] /= systemSize[i] # in row i, we have systemSize[i] nodes to get matched in total. In reduced matrix, we are interested in seeing what porition of these nodes gets assigned to which system. Thus something like sqrt(systemSize[i]*systemSize[j]) wouldn't give this information
    return reducedMatrix

#given a structure/function matching matrix, this function calculates the ratio of
#the summation of the diagonals and the summation of the whole matrix
def calculateIdentityMatchingRatio(matchingMatrix):
    diagSum = matchingMatrix.trace()
    allSum = matchingMatrix.sum()
    if allSum==0:
        matchingAccuracy = 0
    else:
        matchingAccuracy = diagSum*100/float(len(matchingMatrix[0])) #could have been divided by allSum as well
    return matchingAccuracy

######################main code####################
##### get command line parameters
parser = argparse.ArgumentParser(description='Gather values for a random variable in male/female populations')
parser.add_argument('-resS','--resultsSubj', help='File path that contains the average result matrix for subjectwise str-func matching', required=True,type=str)
parser.add_argument('-resP','--resultsPerm', help='File path that contains the average result matrix for permutation testing', required=True,type=str)
parser.add_argument('-out','--outputFolder', help='file path to save the distribution of values for the two populations', required=True,type=str)
parser.add_argument('--sign', help='matching is made btw structure and positive or negative function?', required=True,type=str,choices=['positive','negative'])
parser.add_argument('--yeoNetworkPath', help='File path that contains the mapping of full networks into Yeo 7 systems', required=True,type=str)
parser.add_argument('--numSystems', help='Number of systems to cluster nodes into', required=True,type=int)

args = vars(parser.parse_args())
subjectwiseResultsPath=args['resultsSubj']
permutationTestResultsPath=args['resultsPerm']
outputFolder=args['outputFolder']
sign=args['sign']
yeoNetworkPath=args['yeoNetworkPath']
numSystems=args['numSystems']

#if you would like to discard some of the subjects, enter their order number in here.
ignoreData=[]
    

#read number of nodes and subjects
numNodes, numSubjects = readExperimentSize(subjectwiseResultsPath) #these numbers needs to be the same for both experiments
numSubjects -= len(ignoreData) #reduce the size of the subjects wrt discarded ones from ignoreData list

#Load the mapping of full networks into Yeo 7 systems, and calculate sizes of systems
reductionList8 = np.array([int(x) for x in open(yeoNetworkPath,"r").read().splitlines()],dtype=int)

systemSizeFull=np.ones(numNodes);
systemSize8 = np.zeros(numSystems)
for i in range(numSystems):
    systemSize8[i] = np.sum(reductionList8==i+1)

############1)take average of subjectwise experiment
print("loading average subjectwise results "+subjectwiseResultsPath)
#subjectwiseMatchingMatricesPerSubject is of size (numSubjects x numNodes x numNodes)
#subjectwiseMatchingMatricesPerSubjectFull = averageMatchingMatricesPerSubjectAcrossSeveralExperimentRuns(subjectwiseFileList,numNodes, numSubjects,ignoreData)
subjectwiseMatchingMatricesPerSubjectFull = readAverageMatchingsForAllSubjects(subjectwiseResultsPath,ignoreData)

#subjectwiseAverageMatchingMatrix is of size (numNodes x numNodes)
subjectwiseAverageMatchingMatrixFull = averageMatchingMatrixAcrossSubjects(subjectwiseMatchingMatricesPerSubjectFull)
subjectwiseAverageMatchingMatrix8 = reduceMatchingMatrix(subjectwiseAverageMatchingMatrixFull,reductionList8,systemSize8)
#save average structure-function matching matrices for the experiment (we will use this for visualization of matching regions)
with open(outputFolder+"actualMatchingMatrixFull.txt",'w') as outputFile:
    np.savetxt(outputFile,subjectwiseAverageMatchingMatrixFull,fmt='%1.6f')
with open(outputFolder+"actualMatchingMatrix8.txt",'w') as outputFile:
    np.savetxt(outputFile,subjectwiseAverageMatchingMatrix8,fmt='%1.6f')

############2)take average of permutation test experiment
print("loading average permutation test results "+permutationTestResultsPath)
#permutationMatchingMatrices is of size (numExperiments x numNodes x numNodes)
#permutationMatchingMatricesPerExperimentFull = averageMatchingMatricesPerExperimentRunAcrossSubjects(permutationFileList,numNodes, numSubjects,ignoreData)
permutationMatchingMatricesPerExperimentFull = readAverageMatchingsForAllExperiments(permutationTestResultsPath)

#permutationAverageMatchingMatrixFull is of size (numNodes x numNodes)
permutationAverageMatchingMatrixFull = averageMatchingMatrixAcrossSubjects(permutationMatchingMatricesPerExperimentFull)
permutationAverageMatchingMatrix8 = reduceMatchingMatrix(permutationAverageMatchingMatrixFull,reductionList8,systemSize8)
#save average structure-function matching matrices for the experiment (we will use this for visualization of matching regions)
with open(outputFolder+"permutedMatchingMatrixFull.txt",'w') as outputFile:
    np.savetxt(outputFile,permutationAverageMatchingMatrixFull,fmt='%1.6f')
with open(outputFolder+"permutedMatchingMatrix8.txt",'w') as outputFile:
    np.savetxt(outputFile,permutationAverageMatchingMatrix8,fmt='%1.6f')

############3)calculate pValues by comparing the average of subjectwise test and the permutation test
#we have numExperiments x numSubjects x numNodes x numNodes data for both subjectwise and permutation test experiments
#for the subjectwise results: we will take the average of numExperiments x numSubjects results and obtain a numNodes x numNodes matrix
#  we will save this matrix for reference, and also use it in pValue calculation
#for the permutation test results: we will take the average of numSubjects results and obtain numExperiments x numNodes x numNodes matrix
#  we will use this matrix in pValue calculation
#  we will further take the average of numExperiments x numNodes x numNodes results and obtain numNodes x numNodes matrix to save it for reference
####thus, pValue function will take numNodes x numNodes matrix from subjectwise and numExperiments x numNodes x numNodes matrix from permutation test results.
print("calculating pValues ")
permutationMatchingMatricesPerExperiment8 = np.array([reduceMatchingMatrix(mat,reductionList8,systemSize8) for mat in permutationMatchingMatricesPerExperimentFull])

pValueFull = calculatePValueForPermutationTest(subjectwiseAverageMatchingMatrixFull,permutationMatchingMatricesPerExperimentFull,systemSizeFull,numNodes)
pValue8 = calculatePValueForPermutationTest(subjectwiseAverageMatchingMatrix8,permutationMatchingMatricesPerExperiment8,systemSize8,numNodes)

##do FDR correction for p values -- no need for them, as we are making permutation testing for handling multiple comparison correction
pValueFull_fdr=fdr(pValueFull)
pValue8_fdr=fdr(pValue8)

#save pValue matrices (we will use this while visualizing the statistically significantly matching regions of matching matrices)
with open(outputFolder+"pValuesFull.txt",'w') as outputFile:
    np.savetxt(outputFile,pValueFull,fmt='%1.6f')
with open(outputFolder+"pValues8.txt",'w') as outputFile:
    np.savetxt(outputFile,pValue8,fmt='%1.6f')
       
with open(outputFolder+"significantMatchingsFull.txt",'w') as outputFile:
    np.savetxt(outputFile,subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05),fmt='%1.6f')
with open(outputFolder+"significantMatchings8.txt",'w') as outputFile:
    np.savetxt(outputFile,subjectwiseAverageMatchingMatrix8*(pValue8<=0.05),fmt='%1.6f')
    
   
############4)calculate the overall matching accuracy for each subject by using the subjectwise experiment results
print("calculating matching accuracy for subjectwise matching experiment ")
subjectwiseMatchingMatricesPerSubject8=np.array([reduceMatchingMatrix(mat,reductionList8,systemSize8) for mat in subjectwiseMatchingMatricesPerSubjectFull])
    
identitiyMatchingRatiosFull = [calculateIdentityMatchingRatio(mat*(pValueFull<=0.05)) for mat in subjectwiseMatchingMatricesPerSubjectFull]
identitiyMatchingRatios8 = [calculateIdentityMatchingRatio(mat*(pValue8<=0.05)) for mat in subjectwiseMatchingMatricesPerSubject8]
with open(outputFolder+"subjectwiseAccuracyFull.txt",'w') as outputFile:
    np.savetxt(outputFile,identitiyMatchingRatiosFull,fmt='%1.6f')
with open(outputFolder+"subjectwiseAccuracy8.txt",'w') as outputFile:
    np.savetxt(outputFile,identitiyMatchingRatios8,fmt='%1.6f')    


###########5)save significant average matching accuracies per node
with open(outputFolder+"averageMatchingAccuracyPerNodeFull.txt",'w') as outputFile:
    outputFile.write('# Ratio of the matching accuracy per node\n')
    matchingRatioPerNode = [x for x in np.diag(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05))]
    np.savetxt(outputFile,matchingRatioPerNode,fmt='%1.5f',newline=" ")
with open(outputFolder+"averageMatchingAccuracyPerNode8.txt",'w') as outputFile:
    outputFile.write('# Ratio of the matching accuracy per node\n')
    matchingRatioPerNode = [x for x in np.diag(subjectwiseAverageMatchingMatrix8*(pValue8<=0.05))]
    np.savetxt(outputFile,matchingRatioPerNode,fmt='%1.5f',newline=" ")

############6)save the significant matching accuracy of each node for all subjects
print("save the significant matching accuracy of each node for all subjects")

if sign=="positive":
    for i in range(numNodes):
        with open(outputFolder+"Full/accuracyFull_node"+str(i)+".txt",'w') as outputFile:
            np.savetxt(outputFile,subjectwiseMatchingMatricesPerSubjectFull[:,i,i]*(pValueFull[i,i]<=0.05)*100,fmt='%1.3f')
    for i in range(numSystems):
        with open(outputFolder+"8/accuracy8_node"+str(i)+".txt",'w') as outputFile:
            np.savetxt(outputFile,subjectwiseMatchingMatricesPerSubject8[:,i,i]*(pValue8[i,i]<=0.05)*100,fmt='%1.3f')        
elif sign=="negative":
    for i in range(numSystems):
        with open(outputFolder+"8/accuracy8_node"+str(i)+".txt",'w') as outputFile:
            np.savetxt(outputFile,subjectwiseMatchingMatricesPerSubject8[:,i,6]*(pValue8[i,6]<=0.05)*100,fmt='%1.3f')
    
###########6)save overall statistics of average matching at connectome and systems level
##what is the average matching accuracy
ratio0_full=100*(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05)).trace()/float(len(subjectwiseAverageMatchingMatrixFull[0]))
ratio0_8=100*(subjectwiseAverageMatchingMatrix8*(pValue8<=0.05)).trace()/float(len(subjectwiseAverageMatchingMatrix8[0]))

##what percent of matchings were significant
ratio1_full=100*(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05)).sum()/float(len(subjectwiseAverageMatchingMatrixFull[0]))
ratio1_8=100*(subjectwiseAverageMatchingMatrix8*(pValue8<=0.05)).sum()/float(len(subjectwiseAverageMatchingMatrix8[0]))
##what percent of significant matchings were accurate
ratio2_full=100*(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05)).trace()/(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05)).sum()
ratio2_8=100*(subjectwiseAverageMatchingMatrix8*(pValue8<=0.05)).trace()/(subjectwiseAverageMatchingMatrix8*(pValue8<=0.05)).sum()
##what percent of significant connectome level matchings are within group
significantSubjectwiseAverageMatchingMatrix8 = reduceMatchingMatrix(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05),reductionList8,systemSize8,scaleToSize=False)
ratio3=100*significantSubjectwiseAverageMatchingMatrix8.trace()/(subjectwiseAverageMatchingMatrixFull*(pValueFull<=0.05)).sum()

with open(outputFolder+"overallStatistics.txt",'w') as outputFile:
    outputFile.write('# what is the average matching accuracy =>  full:%1.5f \t systems:%1.5f\n' %(ratio0_full,ratio0_8))
    outputFile.write('# what percent of matchings were significant =>  full:%1.5f \t systems:%1.5f\n' %(ratio1_full,ratio1_8))
    outputFile.write('# what percent of significant matchings were accurate =>  full:%1.5f \t systems:%1.5f\n' %(ratio2_full,ratio2_8))
    outputFile.write('# what percent of significant connectome level matchings are within group =>  %1.5f \n' % ratio3)

print("done ")

