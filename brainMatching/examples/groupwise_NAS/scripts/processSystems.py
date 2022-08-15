import argparse
import numpy as np
from helpers import getOutlierIndices

parser = argparse.ArgumentParser(description='process raw matching data to produces similarity scores for subjects')
parser.add_argument('-r','--resultFile', help='file path to the results of the matching experiment', required=True)
parser.add_argument('-s','--subjectsList', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-o','--outputFile', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('-mt','--measureType', help='which measure to analyze: similarity score or matching accuracy', required=False,type=str,choices=['accuracy','similarity'],default='similarity')
parser.add_argument('-al','--analysisLevel', help='is this a subjectwise or a groupwise analysis?', required=True,type=str,choices=['subject','group'])
parser.add_argument('-rt','--relativeTo', help='calculate similarity relative to which group?', required=False,type=str,choices=['self','healthy'],default='healthy')
parser.add_argument('--noSymmetry', help='do not symmetrize distances', required=False,action='store_true',default=False)


args = vars(parser.parse_args())
resultFilePath=args['resultFile']
subjectListPath=args['subjectsList']
outputFilePath=args['outputFile']
measureType=args['measureType']
analysisLevel=args['analysisLevel']
relativeTo=args['relativeTo']
noSymmetry=args['noSymmetry']

##############load subject IDs for which we have a score#########################

# load healthy & patient IDs from files
tdcPath = "../data/tdc_schaefer.txt"
asdPath = "../data/asd_schaefer.txt"
# tdcPath = "../data/tdc_schaefer_male.txt"
# asdPath = "../data/asd_schaefer_male.txt"
# tdcPath = "../data/tdc_desikan.txt"
# asdPath = "../data/asd_desikan.txt"
# tdcPath = "../data/tdc_desikan_male.txt"
# asdPath = "../data/asd_desikan_male.txt"

with open(tdcPath, "r") as f:
    healthyIDs = f.read().splitlines()

with open(asdPath, "r") as f:
    patientIDs = f.read().splitlines()

healthyOrder = list(range(len(healthyIDs)))
healthy = healthyOrder

# print(len(healthyIDs), len(patientIDs))
# print(healthyOrder[:])
# quit()

with open(subjectListPath,"r") as f:
    subjectList =  f.read().splitlines()
numScans=len(subjectList)


patientS1=[]
patientS2=[]
patientS3=[]

##############load results of the structure-function coupling experiment###################
with open(resultFilePath,"r") as resultFile:
    fileContent = resultFile.read().splitlines()
numNodes = int(fileContent[3].split('\t')[0])
numSubjects = int(fileContent[3].split('\t')[1])

# print(numNodes, numSubjects)
# quit()


###load similarity scores and the matchings
if(analysisLevel=="group"):
    if(measureType=='accuracy'):
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
    scores_avg = np.zeros(len(scores))
    if(relativeTo=="healthy"):
        for row in range(numSubjects):
            count=0
            for col in healthy:
                if(row!=col):
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


scoreName='matching accuracy (%)'
    
outputFile=open(outputFilePath,'w')
outputFile.write("#numNodes,numSubjects\n%d\t%d\n" % (numNodes,numSubjects))
outputFile.write("#measureType\n"+measureType+"\n")
outputFile.write("#scoreName\n"+scoreName+"\n")
outputFile.write("#scores\n")
for i in range(len(scores)):
    outputFile.write("%0.4f\t" % scores[i])
outputFile.close()
