#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:51:28 2018

@author: yusuf
"""


import argparse
import numpy as np
import scipy.stats as stt
import matplotlib.pylab as plt
from helpers import calculateZScore,fdr,calculateGroupDifference,drawHistogram2Dataset,drawBoxPlot,drawViolinPlot
from helpers import checkDifferenceOfVariance


# -r 2_b_strStrDistance_Desikan86_det/results/direct_none_subtract.res -o ./ -s /home/yusuf/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -mt similarity -st dist  -pe png -pt box -al group
# -r 1_a_strFuncCoupling_Desikan86_det/results/accuracy_direct_logScaleEdgesStructure_normalizeEdges_edgesIgnoreDiag_LinAss_0.res -o ./ -s /home/yusuf/data/TBI/resources/qa_commonList_Deterministic_Desikan86.txt -st dist  -pe png -pt violin
#get parameters
parser = argparse.ArgumentParser(description='calculate correlation between the age and matching accuracies of the subjects')
parser.add_argument('-r','--resultFile', help='file path to the results of the matching experiment', required=True)
parser.add_argument('-s','--subjectsList', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-o','--outputFolder', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('-st','--scoreType', help='zScore or the original similarity score (that is, l1-l2 dist etc.)', required=False,type=str,choices=['z','dist'],default='dist')
parser.add_argument('-t','--title', help='title for plot', required=False, type=str, default='') #You can escape white space in the title from command line with "\ " as in <two\ words>
parser.add_argument('-pe','--plotExtension', help='extension of the plot file', required=False,type=str,choices=['png','svg'],default='png')
parser.add_argument('-pt','--plotType', help='draw violin or box plot', required=False,type=str,choices=['violin','box'],default='violin')
parser.add_argument('-tp','--timePoints', help='patients having any time points or only with three time points?', required=False,type=str,choices=['three','any'],default='any')
parser.add_argument('--noHealthyPlot', help='do notplot healthy controls', required=False,action='store_true')

args = vars(parser.parse_args())
resultFile=args['resultFile']
subjectListPath=args['subjectsList']
outputFolder=args['outputFolder']
scoreType=args['scoreType']
title=args['title']
plotExtension=args['plotExtension']
plotType=args['plotType']
timePoints=args['timePoints']
noHealthyPlot=args['noHealthyPlot']

##############load subject IDs for which we have a score#########################

# load healthy & patient IDs from files
tdcPath = "../data/tdc.txt"
asdPath = "../data/asd.txt"
with open(tdcPath, "r") as f:
    healthyIDs = f.read().splitlines()

with open(asdPath, "r") as f:
    patientIDs = f.read().splitlines()

healthySize = len(healthyIDs)
patientSize = len(patientIDs)
healthyOrder = list(range(healthySize))
healthy = healthyOrder
patientsOrder = list(range(healthySize, healthySize+patientSize))
patients = patientsOrder
# print(len(healthy), len(patients))
# print(healthyIDs[healthySize-1])
# print(patientIDs[0], patientIDs[patientSize-1])
# quit()


""" subjectList =  open(subjectListPath,"r").read().splitlines()
numScans=len(subjectList)

patientIDs=[]
healthyIDs=[]
healthyOrder=[]
for order,ID in enumerate(subjectList):
    if (ID.split("_")[0][0]=='c'):
        healthyIDs.append(ID.split('_')[0])
        healthyOrder.append(order)
    elif (ID.split("_")[0][0]=='p') and (ID.split('_')[0] not in patientIDs):
        patientIDs.append(ID.split('_')[0]) """

patientLongitudinalOrder = np.full((len(patientIDs),3),-1,dtype=int)
### Get indices of healthy controls and patients(_s1,_s2,_s3)
""" for i in range(len(subjectList)):
    ID=subjectList[i].split("_")[0]
    if ID in patientIDs:
        timepoint=int(subjectList[i].split("_")[-1][-1])-1
        patientLongitudinalOrder[patientIDs.index(ID)][timepoint] = i """

""" if timePoints=="three":
    rowsToDelete=[]
    for i in range(len(patientLongitudinalOrder)):
        if (patientLongitudinalOrder[i,0]==-1 or patientLongitudinalOrder[i,1]==-1 or patientLongitudinalOrder[i,2]==-1):
                rowsToDelete.append(i)
    patientLongitudinalOrder = np.delete(patientLongitudinalOrder,rowsToDelete,axis=0)
    for i in list(reversed(range(len(patientIDs)))):
        if(i in rowsToDelete):
            del patientIDs[i] """
    

### Get indices of healthy controls and patients(_s1,_s2,_s3)
# healthy=healthyOrder
patientS1=[]
patientS2=[]
patientS3=[]
""" for i in range(len(patientLongitudinalOrder)):
    if (patientLongitudinalOrder[i,0]!=-1):
        patientS1.append(patientLongitudinalOrder[i,0])
    if (patientLongitudinalOrder[i,1]!=-1):
        patientS2.append(patientLongitudinalOrder[i,1])
    if (patientLongitudinalOrder[i,2]!=-1):
        patientS3.append(patientLongitudinalOrder[i,2]) """

if(timePoints=="three"):
    cutOff=6
else:
    cutOff=3


controlGroups = [healthy]
patientGroups = [patients]
controlGroupNames = ["TDC"]
patientGroupNames = ["ASD"]

""" controlGroups=[healthy,healthy,healthy,patientS1,patientS2,patientS1][:cutOff]
patientGroups=[patientS1,patientS2,patientS3,patientS2,patientS3,patientS3][:cutOff]
controlGroupNames=['Healthy','Healthy','Healthy','3_months','6_months','3_months'][:cutOff]
patientGroupNames=['3_months','6_months','12_months','6_months','12_months','12_months'][:cutOff] """

##############load results of the structure-function coupling experiment###################
fileContent =  open(resultFile,"r").read().splitlines()
numNodes = int(fileContent[1].split('\t')[0])
numSubjects = int(fileContent[1].split('\t')[1])
measureType = str(fileContent[3]) 
scoreName = str(fileContent[5]) 
scores = np.zeros(numSubjects)
# load matching accuracy scores
for i in range(numSubjects):
    scores[i] = float(fileContent[7].split('\t')[i])

# print(numNodes, numSubjects)
# print(len(scores))
# print(scores[:5], scores[len(scores)-1])
# quit()

if(scoreType=="z"):
    scores = calculateZScore(scores,healthy)
    scoreName += " (z-score)"


########################calculate group difference########################
colors=['aqua','darkorchid']
effectSize_parametric=np.zeros(len(patientGroups))
effectSize_nonparametric=np.zeros(len(patientGroups))
pValue_parametric=np.zeros(len(patientGroups))
pValue_nonparametric=np.zeros(len(patientGroups))

#calculate group difference btw healthy vs patients at three time points
numTimes = 1 # used to be 3!
for i in range(numTimes): 
    if len(scores[patientGroups[i]])==0:
        continue
    effectSize_parametric[i], pValue_parametric[i] = calculateGroupDifference(scores[controlGroups[i]],scores[patientGroups[i]],parametric=True,paired=False)
    effectSize_nonparametric[i], pValue_nonparametric[i] = calculateGroupDifference(scores[controlGroups[i]],scores[patientGroups[i]],parametric=False,paired=False)
    histogramPath=outputFolder+"histograms/"+timePoints+"_"+measureType+"_"+scoreType+"_"+str(patientGroupNames[i])+"_"+str(controlGroupNames[i])+"."+plotExtension
    drawHistogram2Dataset(scores[healthy],scores[patientGroups[i]],histogramPath,effectSize_parametric[i],pValue_parametric[i],'','',controlGroupNames[i],patientGroupNames[i],"",xLabel=measureType,yLabel="frequency",color1='dodgerblue',color2='magenta')


#calculate group difference among patients at three time points
if(timePoints=="three"):
    isPaired=True
    for i in range(3,6): 
        if len(scores[patientGroups[i]])==0:
            continue
        effectSize_parametric[i], pValue_parametric[i] = calculateGroupDifference(scores[controlGroups[i]],scores[patientGroups[i]],parametric=True,paired=isPaired)
        effectSize_nonparametric[i], pValue_nonparametric[i] = calculateGroupDifference(scores[controlGroups[i]],scores[patientGroups[i]],parametric=False,paired=isPaired)
        histogramPath=outputFolder+"histograms/"+timePoints+"_"+measureType+"_"+scoreType+"_"+str(patientGroupNames[i])+"_"+str(controlGroupNames[i])+".png"
        drawHistogram2Dataset(scores[controlGroups[i]],scores[patientGroups[i]],histogramPath,effectSize_parametric[i],pValue_parametric[i],'','',controlGroupNames[i],patientGroupNames[i],"",xLabel=measureType,yLabel="frequency",color1='dodgerblue',color2='magenta')


corrected_pValue_parametric=fdr(pValue_parametric)
corrected_pValue_nonparametric=fdr(pValue_nonparametric)
# corrected_pValue_parametric=pValue_parametric
# corrected_pValue_nonparametric=pValue_nonparametric

outputPath=outputFolder+timePoints+"_"+measureType+"_groupDifference.txt"

reportFile=open(outputPath,'w')
reportFile.write("============== Dataset statistics ==============\n")
reportFile.write("\tHealthy\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(healthy),np.median(scores[healthy]),scores[healthy].mean(),scores[healthy].var(),scores[healthy].std()))
reportFile.write("\tPatients\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(patients),np.median(scores[patients]),scores[patients].mean(),scores[patients].var(),scores[patients].std()))
""" reportFile.write("\t6 Months\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(patientS2),np.median(scores[patientS2]),scores[patientS2].mean(),scores[patientS2].var(),scores[patientS2].std()))
reportFile.write("\t12 Months\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(patientS3),np.median(scores[patientS3]),scores[patientS3].mean(),scores[patientS3].var(),scores[patientS3].std())) """
           

reportFile.write("============== F/Barlett/Levene Test for checking equalness of variance of patients relative to healthy controls ==============\n")
for discardOutliers in [True,False]:
    reportFile.write("Discard outliers:"+str(discardOutliers)+"\n")
    reportFile.write("\tHealthy vs Patients\t F-test:%f\tBartlett:%f\tLevene:%f \n" % (checkDifferenceOfVariance(scores[patients],scores[healthy],'F',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patients],'bartlett',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patients],'levene',discardOutliers)[1]))
    # reportFile.write("\tHealthy vs 6 months\t F-test:%f\tBartlett:%f\tLevene:%f \n" % (checkDifferenceOfVariance(scores[patientS2],scores[healthy],'F',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patientS2],'bartlett',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patientS2],'levene',discardOutliers)[1]))
    # reportFile.write("\tHealthy vs 12 months\t F-test:%f\tBartlett:%f\tLevene:%f \n" % (checkDifferenceOfVariance(scores[patientS3],scores[healthy],'F',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patientS3],'bartlett',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patientS3],'levene',discardOutliers)[1]))


reportFile.write("============== Group Difference: "+timePoints+" "+measureType+"  "+scoreType+" score==============\n")
reportFile.write("============== Parametric test (with Cohen's D) ==============\n")
for i in range(len(patientGroups)):
    if len(scores[patientGroups[i]])==0:
        continue
    if(corrected_pValue_parametric[i]<=0.05):
        reportFile.write("**")
    reportFile.write("\t"+str(controlGroupNames[i])+" vs "+str(patientGroupNames[i])+"\t:\t%0.2f (%f -- %f) \n" % (effectSize_parametric[i],pValue_parametric[i],corrected_pValue_parametric[i]))
reportFile.write("============== Non-parametric test (Wilcoxon Signed-rank (if paired) test or Mann-Whitney U test (if not paired) for repeated measures with non normal distribution) ==============\n")
for i in range(len(patientGroups)):
    if len(scores[patientGroups[i]])==0:
        continue
    if(corrected_pValue_nonparametric[i]<=0.05):
        reportFile.write("**")
    reportFile.write("\t"+str(controlGroupNames[i])+" vs "+str(patientGroupNames[i])+"\t:\t%0.2f (%f -- %f) \n" % (effectSize_nonparametric[i],pValue_nonparametric[i],corrected_pValue_nonparametric[i]))
reportFile.close()


##########draw boxplot or violin plot matching scores of systems###################
colors=['#351C4D', '#AB3E16','#849974','#2096BA','#F7DFD4','#F5AB99'] #nightfall, rust, fresh, shutter blue, macaron, tropical pink
if(noHealthyPlot==True):
    # numTimePoints=3
    numTimePoints = 1
    outputPath=outputFolder+timePoints+"_"+measureType+"_"+scoreType+"_noHealthyPlot."+plotExtension
else:
    # numTimePoints=4
    numTimePoints = 2
    outputPath=outputFolder+timePoints+"_"+measureType+"_"+scoreType+"."+plotExtension

# data=[scores[patientS1],scores[patientS2],scores[patientS3],scores[healthy]][:numTimePoints]
# dataLabels=['3 months','6 months','12 months','Healthy'][:numTimePoints]

data = [scores[patients], scores[healthy]][:numTimePoints]
dataLabels=['Patient','Healthy'][:numTimePoints]

minY=min([min(l) for l in data])
maxY=max([max(l) for l in data])
offset=(maxY-minY)/5.0
#generate enough empty space above and below boxes
yLim=[minY-offset/2.0,maxY+offset] ## use these lines to  make space specific to figure
# yLim=[67,102] ## use these lines to make the space constant (such as across different plots)
scoreName="network similarity (%)"
if plotType=="box":
    drawBoxPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=scoreName,colors=colors,rotation=0,plotScatter=True,yLim=yLim,middleLine='median') #since we use Mann-Whitney U test for group dofference, we should plot median line in boxplots
elif plotType=="violin":
    drawViolinPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=scoreName,colors=colors,rotation=0,yLim=yLim)
#def drawViolinPlot(data,dataLabels,title,outputPath,xLabel="",yLabel="matching accuracy (%)",colors=['#2096BA','#351C4D', '#AB3E16','#849974', '#F7DFD4','#F5AB99'],rotation=0,yLim=[]):
