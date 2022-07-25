#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 17:10:09 2017

@author: yusuf

sample: python calcAvgMatchingResults.py -res ../direct/positive/src/ -out ./
"""

import argparse
import os
import numpy as np

########helper functions################
#this function reads the experiment information from the result file
#takes path to a result file @filePath as input
#and returns @numNodes and @numSubjects as output
def readExperimentSize(filePath):
    fileContent =  open(filePath,"r").read().splitlines()
    numNodes = int(fileContent[3].split('\t')[0])
    numSubjects = int(fileContent[3].split('\t')[1])
    return numNodes, numSubjects

#this funciton reads the command that was entered to run the experiment on the cluster and returns it
def readCommandLine(filePath):
    fileContent =  open(filePath,"r").read().splitlines()
    commandLine = fileContent[1]
    return commandLine

#this function reads the expeirment results which has matching of 
#structural nodes to function nodes for each subject (a matrix of numSubjectsXnumNodes)
#and returns the content in an integer matrix
def readMatchingsOfAllSubjectsForSingleExperimentRun(filePath,ignoreData):
    fileContent =  open(filePath,"r").read().splitlines()
    for line in range(5): #discard the initial 5 lines where we keep the comments regarding the file content
        del fileContent[0]
    del fileContent[-1] #discard the last empty line
    for line in reversed(ignoreData):#discard the matching results for the subjects that do not meet our QA
        del fileContent[line]
    #at this point, fileContent keeps the matchings for each subject
    matchings=[np.fromstring(cont,dtype=int,sep='\t') for cont in fileContent]
    return matchings

#given the folder that contains results of the same experiment repeated several times, 
#where each result file contains the corresponding functional node for each structural node, for each subject
#this function calculates the average matching matrix for each subject, and returns the resulting matrix (numSubjects x numNodes x numNodes)
def averageMatchingMatricesPerSubjectAcrossSeveralExperimentRuns(fileList, numNodes, numSubjects,ignoreData):
    matchingMatrices = np.zeros(shape=(numSubjects,numNodes,numNodes),dtype=float)
    for fileName in fileList:
        matchings = readMatchingsOfAllSubjectsForSingleExperimentRun(fileName,ignoreData)
        for i in range(numSubjects):
            for j in range(numNodes):
                matchingMatrices[i][j][matchings[i][j]]+=1
    matchingMatrices /= float(len(fileList))
    return matchingMatrices

#given the folder that contains results of the same experiment repeated several times, 
#where each result file contains the corresponding functional node for each structural node, for each subject
#this function calculates the average matching matrix for each experiment, and returns the resulting matrix (numExperiments x numNodes x numNodes)
def averageMatchingMatricesPerExperimentRunAcrossSubjects(fileList, numNodes, numSubjects,ignoreData):
    numExperiments = len(fileList)
    matchingMatrices = np.zeros(shape=(numExperiments,numNodes,numNodes),dtype=float)
    for experimentId in range(numExperiments):
        matchings = readMatchingsOfAllSubjectsForSingleExperimentRun(fileList[experimentId],ignoreData)
        for i in range(numSubjects):
            for j in range(numNodes):
                matchingMatrices[experimentId][j][matchings[i][j]]+=1
    matchingMatrices /= float(numSubjects)
    return matchingMatrices


#given a file that contains the corresponding functional node for each structural node, for each subject
#this function calculates the average matching matrix for the experiment, and returns the resulting matrix (numNodes x numNodes)
def averageMatchingMatricesAcrossSubjectsForSingleExperimentRun(fileName, numNodes, numSubjects,ignoreData):
    avgMatchingMatrix = np.zeros(shape=(numNodes,numNodes),dtype=float)
    matchings = readMatchingsOfAllSubjectsForSingleExperimentRun(fileName,ignoreData)
    for i in range(numSubjects):
        for j in range(numNodes):
            avgMatchingMatrix[j][matchings[i][j]]+=1
    avgMatchingMatrix /= float(numSubjects)
    return avgMatchingMatrix

######################main code####################
##### get command line parameters
parser = argparse.ArgumentParser(description='calculate average matching matrices for subjectwise str/func matching and permutation test experiments')
parser.add_argument('-res','--resultsFolder', help='folder path that contains the result matrices', required=True)
parser.add_argument('-et','--experimentType', help='subjectwise and/or permutationTest', required=False,type=str,choices=['subjectwise','permutationTest','both'],default='both')

args = vars(parser.parse_args())
resultsFolder=args['resultsFolder']
experimentType=args['experimentType']


#paths for reading results of experiments
subjectwiseResultsPath=os.path.join(resultsFolder,"rawData/subjectwise/")
permutationTestResultsPath=os.path.join(resultsFolder,"rawData/permutationTest/")

#some of the subjects had disconnected structural connectomes. We ignore them in our analysis
ignoreData=[]

#load filenames for result files
subjectwiseFileList = [subjectwiseResultsPath + x for x in os.listdir(subjectwiseResultsPath)]
permutationFileList = [permutationTestResultsPath + x for x in os.listdir(permutationTestResultsPath)]
numNodes, numSubjects = readExperimentSize(subjectwiseFileList[0]) #these numbers needs to be the same for both experiments
numSubjects -= len(ignoreData) #reduce the size of the subjects wrt discarded ones from ignoreData list

############1)take average of subjectwise experiment
if(experimentType=="subjectwise" or experimentType=="both"):
    print("taking average of the subjectwise results "+resultsFolder)
    commandLine=readCommandLine(subjectwiseFileList[0])
    #subjectwiseMatchingMatricesPerSubject is of size (numSubjects x numNodes x numNodes)
    subjectwiseMatchingMatricesPerSubjectFull = averageMatchingMatricesPerSubjectAcrossSeveralExperimentRuns(subjectwiseFileList,numNodes, numSubjects,ignoreData)
    
    outputPath=os.path.join(resultsFolder,"rawData/averageSubjectwiseResults.txt")
    outputFile=open(outputPath,'w')
    outputFile.write("# Sample command line\n")
    outputFile.write(commandLine+"\n")
    outputFile.write("# <numNodes> <numSubjects>\n")
    outputFile.write(str(numNodes)+"\t"+str(numSubjects)+"\n")
    delimiter="\t"
    fmt="%.3f"
    for i in range(numSubjects):
        outputFile.write("# average matching matrix for subject "+str(i)+"\n")
        for row in subjectwiseMatchingMatricesPerSubjectFull[i]:
            line = delimiter.join("0" if value == 0 else fmt % value for value in row)
            outputFile.write(line + '\n')
    outputFile.close()


############2)take average of permutation test experiment
if(experimentType=="permutationTest" or experimentType=="both"):
    print("taking average of the permutation test results "+resultsFolder)
    commandLine=readCommandLine(permutationFileList[0])
    #permutationMatchingMatrices is of size (numExperiments x numNodes x numNodes)
    permutationMatchingMatricesPerExperimentFull = averageMatchingMatricesPerExperimentRunAcrossSubjects(permutationFileList,numNodes, numSubjects,ignoreData)
    
    numExperiments=len(permutationMatchingMatricesPerExperimentFull)
    
    outputPath=os.path.join(resultsFolder,"rawData/averagePermutationTestResults.txt")
    outputFile=open(outputPath,'w')
    outputFile.write("# Sample command line\n")
    outputFile.write(commandLine+"\n")
    outputFile.write("# <numNodes> <numExperiments>\n")
    outputFile.write(str(numNodes)+"\t"+str(numExperiments)+"\n")
    delimiter="\t"
    fmt="%.3f"
    for i in range(numExperiments):
        outputFile.write("# average permutation test matching matrix for experiment "+str(i)+"\n")
        for row in permutationMatchingMatricesPerExperimentFull[i]:
            line = delimiter.join("0" if value == 0 else fmt % value for value in row)
            outputFile.write(line + '\n')
    outputFile.close()
