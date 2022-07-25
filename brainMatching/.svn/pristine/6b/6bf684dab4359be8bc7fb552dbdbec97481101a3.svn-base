#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:51:28 2018

@author: yusuf
"""


import argparse
import numpy as np
from helpers import calculateZScore

# -r 2_a_strStrCoupling_Desikan86_det/results/direct_none_edgesIncludeDiagZeroDiag_LinAss_0.res -sl /home/yusuf/data/TBI/resources/qa_fullList_Deterministic_Desikan86.txt -tOut ./ -pOut ./ -mt accuracy -al subject --title someTitle --scoreType dist
 
#get parameters
parser = argparse.ArgumentParser(description='calculate correlation between the age and matching accuracies of the subjects')
parser.add_argument('-r','--resultFile', help='file path to the results of the matching experiment', required=True)
parser.add_argument('-sl','--subjectsList', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-sip','--subjectsInfoPath', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-st','--scoreType', help='zScore or the original similarity score (that is, l1-l2 dist etc.)', required=False,type=str,choices=['z','dist'],default='dist')
parser.add_argument('-tOut','--textOutputFolder', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('--noSymmetry', help='do not symmetrize distances', required=False,action='store_true',default=False)
parser.add_argument('-tp','--timePoints', help='patients having all time points or only with three time points?', required=False,type=str,choices=['three','any'],default='any')
parser.add_argument('-ds','--dropSubjects',nargs='*', help='list of subjects to from the analysis', required=False, type=str,default=[])

args = vars(parser.parse_args())
resultFile=args['resultFile']
subjectListPath=args['subjectsList']
subjectsInfoPath=args['subjectsInfoPath']
textOutputFolder=args['textOutputFolder']
scoreType=args['scoreType']
noSymmetry=args['noSymmetry']
timePoints=args['timePoints']
subjectsToDrop=args['dropSubjects']

##############load subjet IDs for which we have a score#########################
subjectList =  open(subjectListPath,"r").read().splitlines()
numSubjects = len(subjectList)
##############load results of the structure-function coupling experiment###################
scores=np.zeros(numSubjects)
if(resultFile!='none'):
    fileContent =  open(resultFile,"r").read().splitlines()
    numNodes = int(fileContent[1].split('\t')[0])
    numSubjects = int(fileContent[1].split('\t')[1])
    measureType = str(fileContent[3]) 
    scoreName = str(fileContent[5]) 
    for i in range(numSubjects):
        scores[i] = float(fileContent[7].split('\t')[i])

#######################load cognitive scores########################
###load subjects info
#subjectsInfoPath='/home/yusuf/data/TBI/resources/TBI_Data_20200325.csv'
subjectsInfo = np.genfromtxt(subjectsInfoPath, names=True, delimiter=',', dtype=None, missing_values='.', filling_values='-10000',encoding=None)
sbiaIDs = subjectsInfo['Subject'].tolist()

###remove the subjects from the score table for which we did not calculate score
removeList=[]
for i in reversed(range(len(sbiaIDs))):
    if (sbiaIDs[i] not in subjectList) or (sbiaIDs[i] in subjectsToDrop):
        removeList.append(i)
subjectsInfo = np.delete(subjectsInfo,removeList,axis=0)
sbiaIDs = subjectsInfo['Subject'].tolist()

removeList=[] #also remove the subjects which do not pass QA
for i in reversed(range(len(subjectList))):
    if (subjectList[i] not in sbiaIDs) or (subjectList[i] in subjectsToDrop) :
        removeList.append(i)
subjectList = np.delete(subjectList,removeList,axis=0)
scores = np.delete(scores,removeList,axis=0)

############## get subject lists
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

if timePoints=="three":
    rowsToDelete=[]
    for i in range(len(patientLongitudinalOrder)):
        if (patientLongitudinalOrder[i,0]==-1 or patientLongitudinalOrder[i,1]==-1 or patientLongitudinalOrder[i,2]==-1):
                rowsToDelete.append(i)
    patientLongitudinalOrder = np.delete(patientLongitudinalOrder,rowsToDelete,axis=0)
    for i in list(reversed(range(len(patientIDs)))):
        if(i in rowsToDelete):
            del patientIDs[i]

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

sampleGroups=[patientS1,patientS2,patientS3] #[healthy,patientS1,patientS2,patientS3]
sampleGroupNames=['3_months','6_months','12_months'] #['Healthy','3_months','6_months','12_months']

########### calculate z-score of matching scores if indicated
if(scoreType=="z"):
    scores = calculateZScore(scores,healthy)
    scoreName += " (z-score)"

### Now, get the cognitive scores for each subject
#executive = np.array(subjectsInfo[['Trails_B__TScore', 'DSB__Scaled', 'LNS__Scaled', 'COWA__Adjusted', 'CWIT_3__Scaled']].tolist())
#calculate an executive score based on the rank of individual scores
#executive = (executive.argsort(axis=0).argsort(axis=0) / float(numNodes)).sum(axis=1) 
executive = list(subjectsInfo['EF'])+[-1] #['Executive_Composite']  ### NOTE: we append '.' to the end of each score to use it as a replcaement for missing values while writing into file (i.e., index [-1] shold have '.' for all score types)
proSpeed = list(subjectsInfo['PSIT'])+[-1]#['PSI']
verbal = list(subjectsInfo['ReyT'])+[-1]#['Rey_Sum_I__V']
pta = list(subjectsInfo['PTA_Estimated'])+[-1]
gose = list(subjectsInfo['GOSE'])+[-1]
drs = list(subjectsInfo['DRS'])+[-1]
age = list(subjectsInfo['Age'])+[-1]
gender = list(subjectsInfo['Gender'])+['none']
daysSinceInjury = list(subjectsInfo['DaysSinceInjury'])+[-1]

scores=list(scores) + [-1]

# calculate baseline age at the first scan
age_baseline = np.zeros(len(age)-1)
unique_subjects=[]
min_age=[]
for order,ID in enumerate(subjectList):
    uniqueID = ID.split("_")[0]
    if(uniqueID not in unique_subjects):
        unique_subjects.append(uniqueID)
        min_age.append(age[order])
    else:
        if(min_age[-1] > age[order]):
            min_age[-1] = age[order]
for order,ID in enumerate(subjectList):
    uniqueID = ID.split("_")[0]
    age_baseline[order] = min_age[unique_subjects.index(uniqueID)]


### put all scores into a single array
cognitiveScores= [executive, proSpeed, verbal, pta, gose, drs]
cognitiveScoreNames = ['executive', 'proSpeed', 'verbal', 'pta', 'gose', 'drs']
numScores=len(cognitiveScoreNames)

### Save scores into a file

#for sc in range(numScores):
#    reportFile=open(textOutputFolder+cognitiveScoreNames[sc]+"_"+measureType+".txt",'w')
#    reportFile.write("subjectId\ttimePoint\tage\tgender\tdaysSinceInjury\tsimilarityScore\tcognitiveScore\n")
#    for i in range(len(healthy),len(subjectList)):
#        ID=subjectList[i].split("_")[0]
#        if(timePoints=="three"):
#            if(cognitiveScores[sc][i]>=0 and ID+"_s1" in subjectList and ID+"_s2" in subjectList and ID+"_s3" in subjectList):
#                reportFile.write("%s\t%s\t%d\t%s\t%d\t%f\t%f\n" % (subjectList[i].split("_")[0],subjectList[i].split("_")[1],age[i],gender[i],daysSinceInjury[i],scores[i], cognitiveScores[sc][i]))
#        else:
#            if(cognitiveScores[sc][i]>=0):
#                reportFile.write("%s\t%s\t%d\t%s\t%d\t%f\t%f\n" % (subjectList[i].split("_")[0],subjectList[i].split("_")[1],age[i],gender[i],daysSinceInjury[i],scores[i], cognitiveScores[sc][i]))
#    reportFile.close()    
    
reportFile=open(textOutputFolder+"allScores.txt",'w')
reportFile.write("subjectId\ttimePoint\tage\tageBaseline\tgender\tdaysSinceInjury\tsimilarityScore\texecutive\tproSpeed\tverbal\tpta\tgose\tdrs\n")
for i in range(len(healthy),len(subjectList)):
    ### Note: we look for scores having value greater than -5. This is becuase missing data will be replaced with -10000. The rest should be kept. In this dataset, the smallest valid data is greater than -5.
    if(timePoints=="three"):
        if(all(score[i]>=-5 for score in cognitiveScores) and ID+"_s1" in subjectList and ID+"_s2" in subjectList and ID+"_s3" in subjectList):
            reportFile.write("%s\t%s\t%d\t%d\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (subjectList[i].split("_")[0],subjectList[i].split("_")[1],age[i],age_baseline[i],gender[i],daysSinceInjury[i],scores[i], executive[i],proSpeed[i],verbal[i],pta[i],gose[i],drs[i]))
    else:
        if(all(score[i]>=-5 for score in cognitiveScores)):
            reportFile.write("%s\t%s\t%d\t%d\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (subjectList[i].split("_")[0],subjectList[i].split("_")[1],age[i],age_baseline[i],gender[i],daysSinceInjury[i],scores[i], executive[i],proSpeed[i],verbal[i],pta[i],gose[i],drs[i]))
reportFile.close()  

reportFile=open(textOutputFolder+"allScores_healthy.txt",'w')
reportFile.write("subjectId\tage\tgender\tsimilarity1\tEF1\tPSI1\tVL1\n")
for i in range(len(healthy)):
    reportFile.write("%s\t%d\t%s\t%f\t%f\t%f\t%f\n" % (subjectList[i].split("_")[0],age[i],gender[i],scores[i], executive[i],proSpeed[i],verbal[i]))
reportFile.close()  


reportFile=open(textOutputFolder+"allScores_subjectwise.txt",'w')
reportFile.write("subjectId\ttimePoints\tage1\tage2\tage3\tgender\tdsi1\tdsi2\tdsi3\tsimilarity1\tsimilarity2\tsimilarity3\tEF1\tEF2\tEF3\tPSI1\tPSI2\tPSI3\tVL1\tVL2\tVL3\tpta\tGOSE1\tGOSE2\tGOSE3\n")
for i in range(len(patientLongitudinalOrder)):
    ### keep the indices of timepoints for which this subject was kept in the experiment. Ex: If appeared only in first and third time points, then s1_s3...
    timePoints="s1"
    if(patientLongitudinalOrder[i,1]>=0):
        timePoints = timePoints+"_s2"
    if(patientLongitudinalOrder[i,2]>=0):
        timePoints = timePoints+"_s3"
        
    reportFile.write("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\n" 
                     % (subjectList[patientLongitudinalOrder[i,0]].split("_")[0],timePoints,
                        age[patientLongitudinalOrder[i,0]],age[patientLongitudinalOrder[i,1]],age[patientLongitudinalOrder[i,2]],
                        gender[patientLongitudinalOrder[i,0]],
                        daysSinceInjury[patientLongitudinalOrder[i,0]],daysSinceInjury[patientLongitudinalOrder[i,1]],daysSinceInjury[patientLongitudinalOrder[i,2]],
                        scores[patientLongitudinalOrder[i,0]],scores[patientLongitudinalOrder[i,1]],scores[patientLongitudinalOrder[i,2]], 
                        executive[patientLongitudinalOrder[i,0]],executive[patientLongitudinalOrder[i,1]],executive[patientLongitudinalOrder[i,2]],
                        proSpeed[patientLongitudinalOrder[i,0]],proSpeed[patientLongitudinalOrder[i,1]],proSpeed[patientLongitudinalOrder[i,2]],
                        verbal[patientLongitudinalOrder[i,0]],verbal[patientLongitudinalOrder[i,1]],verbal[patientLongitudinalOrder[i,2]],
                        pta[patientLongitudinalOrder[i,0]],
                        gose[patientLongitudinalOrder[i,0]],gose[patientLongitudinalOrder[i,1]],gose[patientLongitudinalOrder[i,2]]))
reportFile.close()   