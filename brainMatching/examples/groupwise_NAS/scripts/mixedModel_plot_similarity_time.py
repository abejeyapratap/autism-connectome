#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:51:28 2018

@author: yusuf
"""


import argparse
import numpy as np
import scipy.stats as stt
from helpers import fdr
import matplotlib.pylab as plt

   
# -r 2_a_strStrCoupling_Desikan86_det/results/direct_none_edgesIncludeDiagZeroDiag_LinAss_corrMixMod/data/drs_accuracy.txt -tOut ./ -pOut ./ --title someTitle
 
#get parameters
parser = argparse.ArgumentParser(description='calculate correlation between the age and matching accuracies of the subjects')
parser.add_argument('-m','--modelsSummary', help='file path to the results of the mixed models effect model', required=True)
parser.add_argument('-s','--scores', help='file path to the file containing cognitive and matching accuracy scores', required=True)
parser.add_argument('-pOut','--plotOutputFolder', help='file path to save the distribution of values for the two populations', required=True)
parser.add_argument('--plotExtension', help='extension of the plot file', required=False,type=str,choices=['png','svg'],default='png')
parser.add_argument('--coloringScheme', help='choose scatter plot dots coloring scheme', required=False,type=str,choices=['gender','group'],default='group')
parser.add_argument('--dimColorWithAge', help='dim color of scatter plot dots according to age', required=False,action='store_true')
parser.add_argument('--noSubjectName', help='do not write subject names on top of scatter plot dots', required=False,action='store_true')
parser.add_argument('--grayLines', help='plot gray lines across the same subjects dots', required=False,type=str,choices=['spaghetti','mixedModel','none'],default='none')
parser.add_argument('--lineCoef', help='coefficients of the gray lines to be plotted', required=False,type=str,default="")
parser.add_argument('--discardTimePoint', help='discard patients at a certain time point?', required=False,type=str,choices=['s1','s2','s3','none'],default='none')
parser.add_argument('--verbose', help='print debug comments', required=False,action='store_true')

args = vars(parser.parse_args())
modelsSummaryPath=args['modelsSummary']
scoresPath=args['scores']
plotOutputFolder=args['plotOutputFolder']
plotExtension=args['plotExtension']
coloringScheme=args['coloringScheme']
dimColorWithAge=args['dimColorWithAge']
noSubjectName=args['noSubjectName']
grayLines=args['grayLines']
lineCoefFile=args['lineCoef']
discardTimePoint=args['discardTimePoint']
verbose=args['verbose']

################## Read model summary for the cognitive scores: put all scores into a single array
modelsSummary =  open(modelsSummaryPath+discardTimePoint+"Discarded.txt","r").read().splitlines()

model=modelsSummary[1][2:]

if(verbose==True):
    valueNames=modelsSummary[5].split()
    print("%s %s %s %s %s %s %s %s\n" %(valueNames[0],valueNames[1],valueNames[2],valueNames[3],valueNames[4],valueNames[5],valueNames[6],valueNames[7],valueNames[8]))

values=[]
modelNames = []#['executive', 'proSpeed', 'verbal', 'pta', 'gose', 'drs']
for i in [5]: ### the model I'm interested in is recorded in 6th line of the summary file
    modelNames.append(modelsSummary[i].split()[0])
    
    values.append([float(x) for x in modelsSummary[i].split()[1:]])
    
    if(verbose==True):
        print("%s %f %f %f %f %f %f %f %f\n" %(modelNames[-1],values[-1][0],values[-1][1],values[-1][2],values[-1][3],values[-1][4],values[-1][5],values[-1][6],values[-1][7]))
    
###multiple comparison correction
pValues_LRT=[line[0] for line in values]
pValues_PBtest=[line[1] for line in values]
## enable lines below for FDR correction
#pValues_LRT=fdr(pValues_LRT)
#pValues_PBtest=fdr(pValues_PBtest)
for i in range(len(values)):
    values[i][0]=pValues_LRT[i]
    values[i][1]=pValues_PBtest[i]

##############load cognitive and similarity scores
scoresTable = np.genfromtxt(scoresPath, names=True, delimiter='\t', dtype=None,encoding=None)

###discard patients at a certain time point, if any is indicated to be discarded with --discardTimePoint flag
scoresTable = scoresTable[scoresTable['timePoint']!=discardTimePoint]

gender = [str(x) for x in scoresTable['gender'].tolist()]
subjectIDs = [str(x) for x in scoresTable['subjectId'].tolist()]
timePoints = [str(x) for x in scoresTable['timePoint'].tolist()]
similarity = scoresTable['similarityScore']
age = scoresTable['ageBaseline']
dsi = scoresTable['daysSinceInjury']
pta = scoresTable['pta']
genderBool = np.array([1 if x=="M" else 0 for x in scoresTable['gender'].tolist()])

numScores=len(similarity)
maxAge=float(max(age))

colors=['#351C4D', '#AB3E16','#849974','#2096BA','#F7DFD4','#F5AB99'] #nightfall, rust, fresh, shutter blue, macaron, tropical pink

        
########## do plotting
for i in range(len(modelNames)):
    yLabel='network similarity (%)'
    xLabel='days since injury'
    yData=similarity
    xData=dsi
    outputPath = plotOutputFolder+modelNames[i]+"_"+grayLines+"Lines_"+discardTimePoint+"Discarded."+plotExtension
    
    numSubjects=len(subjectIDs)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_ylabel(yLabel,fontsize=16)
    ax.set_xlabel(xLabel,fontsize=16)
    ax.text(0.05,0.92, '$Adj.R^{2} = $'+("{0:.3f}".format(values[i][7])) + '\n$p_{dsi} = $'+("{0:.3f}".format(values[i][2])), #'p = '+("{0:.3f}".format(values[i][0]))+'\np = '+("{0:.3f}".format(values[i][1])),
            horizontalalignment='left',verticalalignment='center',
            transform = ax.transAxes, color='grey',
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    plt.title(model)
    
    plt.ylim(67,102)
    
    #plot scatter plot dots
    yTextOffset=(np.nanmax(yData)-np.nanmin(yData))/45.0 ## move the label of each node up in the y axis by this much
    xTextOffset=(np.nanmax(xData)-np.nanmin(xData))/45.0 ## move the label of each node left in the x axis by this much
    #    plt.scatter(xData,yData,c='skyblue',edgecolors='none',alpha=0.7) # an alternative way of plotting points, where each node is colored single color
    for j in range(numScores):
        if(dimColorWithAge==True): # if age is to be taken into account in coloring, set alpha relative to age
            alpha=float(age[j])/maxAge
        else:#otherwise, set alpha to a constant transparency for all
            alpha=0.7
        if(coloringScheme=='gender'): #is coloring going to make males blue and females red? 
            if(gender[j]=="F"):
                plt.scatter(xData[j],yData[j],color='darkorchid',alpha=alpha,edgecolors='dimgray',linewidth=0.3)
            else:
                plt.scatter(xData[j],yData[j],color='dodgerblue',alpha=alpha,edgecolors='dimgray',linewidth=0.3)
        elif(coloringScheme=='group'): # or the coloring will consists of subjects in the same group having the same color
            plt.scatter(xData[j],yData[j],color=colors[int(timePoints[j][1])-1],alpha=alpha,edgecolors='dimgray',linewidth=0.3)
        
        if(noSubjectName==False): #write names of next to scatter plot dots
            ax.annotate(subjectIDs[j],(xData[j]-xTextOffset,yData[j]+yTextOffset),size=5,alpha=0.5) #ax.annotate(subjectIDs[j]+"_"+timePoints[j],(xData[j],yData[j]+yTextOffset),size=5,alpha=0.5)
    
    #plot fitted line
    intercept_fitted=values[i][3]
    slope_fitted_dsi=values[i][4]
    slope_fitted_age=values[i][5]
    slope_fitted_pta=values[i][6]
    
    regressedYData=intercept_fitted+ slope_fitted_dsi * xData + slope_fitted_age * np.average(age) + slope_fitted_pta * np.average(pta)
    plt.plot(xData,regressedYData,color='darkorchid',linewidth=2)
        

    # plot individual fitted lines per subject
    # first obtain the coordinates of nodes
    uniqueSubjects=[]
    similarityLongitudinal = []
    dsiLongitudinal = []
    simLine=[]
    dsiLine=[]
    for j in range(len(subjectIDs)):
        if(subjectIDs[j] not in uniqueSubjects):#if this is the first time we are encountering the subjectID in the subject list, generate a new line
            if(j!=0):
                similarityLongitudinal.append(simLine)
                dsiLongitudinal.append(dsiLine)
            uniqueSubjects.append(subjectIDs[j])
            simLine=[similarity[j]]
            dsiLine=[dsi[j]]
        else: #if I have seen this subject's ID before, than add this new timepoint to the line
            simLine.append(similarity[j])
            dsiLine.append(dsi[j])
        if(j==len(subjectIDs)-1):
            similarityLongitudinal.append(simLine)
            dsiLongitudinal.append(dsiLine)
    
    if(grayLines=="spaghetti"): #plot lines as a spaghettin plot
        for j in range(len(uniqueSubjects)):
            for k in range(len(similarityLongitudinal[j])-1):
                plt.plot(dsiLongitudinal[j][k:k+2],similarityLongitudinal[j][k:k+2],color='grey',alpha=0.2,linewidth=0.7)
    elif(grayLines=="mixedModel"): #or, plot a single line parallel to the main fitted line per subject
        line_parameterNames = open(lineCoefFile,"r").read().splitlines()[0].split(" ")
        lineCoefficients = np.genfromtxt(lineCoefFile, names=None, delimiter=' ', dtype=None,encoding=None,skip_header=1)
        
        for j in range(len(lineCoefficients)):
            subjectID = lineCoefficients[j][0]
            intercept_ = lineCoefficients[j][1]
            slope_dsi = lineCoefficients[j][2]
            slope_age = lineCoefficients[j][3]
            
            if scoresTable[scoresTable['subjectId']==subjectID]["gender"][0]=="M":
                genderSubject = 1
            else:
                genderSubject = 0
            ageSubject = scoresTable[scoresTable['subjectId']==subjectID]["ageBaseline"][0]
            
            ## to plot grey lines to span entire range of similarity scores
            dsiRange = np.array([min(xData),max(xData)]) #use this line to draw grey lines across the whole range
            # similarityRange= np.array(similarityLongitudinal[j]) #use this line to draw lines between the starting and ending points of the data per subject
            cogRegressed = list(intercept_ + slope_dsi*dsiRange + slope_age*ageSubject) 
            plt.plot(dsiRange,cogRegressed,color='grey',alpha=0.2,linewidth=0.7)
            

            if(verbose==True):
                print(cogRegressed)
            
#    plt.title("Correlation btw sim. and cog. scores")
    plt.savefig(outputPath, transparent=False, dpi=300, bbox_inches='tight')
    plt.close()