import os
import argparse
import numpy as np
import scipy.stats as stt
import matplotlib.pylab as plt
from helpers import calculateZScore,fdr,calculateGroupDifference,drawHistogram2Dataset,drawBoxPlot,drawViolinPlot
from helpers import checkDifferenceOfVariance

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
# resultFile=args['resultFile']
subjectListPath=args['subjectsList']
# outputFolder=args['outputFolder']
scoreType=args['scoreType']
title=args['title']
plotExtension=args['plotExtension']
plotType=args['plotType']
timePoints=args['timePoints']
noHealthyPlot=args['noHealthyPlot']


# Load healthy & patient IDs from files
tdcPath = "../data/tdc_schaefer.txt"
asdPath = "../data/asd_schaefer.txt"

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


controlGroups = [healthy]
patientGroups = [patients]
controlGroupNames = ["TDC"]
patientGroupNames = ["ASD"]


##############load results of average NNS processing ###################
systemMatchingPath = "../experiment/results/sys_level"
sysNames = ["visual", "somatomotor", "dorsal", "ventral", "limbic", "frontoparietal", "defaultmode", "subcortical"]
# sysNames = ["visual"]

# create output folders for each system
outputFolder = "../experiment/results/sys_level/plots"
for system in sysNames:
    sysPath = f"{outputFolder}/{system}"
    if not os.path.exists(sysPath):
        os.mkdir(sysPath)

    histPath = f"{outputFolder}/{system}/histogram"
    if not os.path.exists(histPath):
        os.mkdir(histPath)


colors=['aqua','darkorchid']
for ind, system in enumerate(sysNames):
    resultFile = f"{systemMatchingPath}/{system}.res"

    fileContent =  open(resultFile,"r").read().splitlines()
    numNodes = int(fileContent[1].split('\t')[0])
    numSubjects = int(fileContent[1].split('\t')[1])
    measureType = str(fileContent[3]) 
    scoreName = str(fileContent[5]) 
    scores = np.zeros(numSubjects)

    # load matching accuracy scores
    for i in range(numSubjects):
        scores[i] = float(fileContent[7].split('\t')[i])


    ########################calculate group difference########################
    effectSize_parametric=np.zeros(len(patientGroups))
    effectSize_nonparametric=np.zeros(len(patientGroups))
    pValue_parametric=np.zeros(len(patientGroups))
    pValue_nonparametric=np.zeros(len(patientGroups))


    # Calculate group difference b/w healthy vs patients
    effectSize_parametric[0], pValue_parametric[0] = calculateGroupDifference(scores[controlGroups[0]],scores[patientGroups[0]],parametric=True,paired=False)

    effectSize_nonparametric[0], pValue_nonparametric[0] = calculateGroupDifference(scores[controlGroups[0]],scores[patientGroups[0]],parametric=False,paired=False)

    outputPath = f"{outputFolder}/{sysNames[ind]}/"
    histogramPath=outputPath+"histogram/"+timePoints+"_"+measureType+"_"+scoreType+"_"+str(patientGroupNames[0])+"_"+str(controlGroupNames[0])+"."+plotExtension
    drawHistogram2Dataset(scores[healthy],scores[patientGroups[0]],histogramPath,effectSize_parametric[0],pValue_parametric[0],'','',controlGroupNames[0],patientGroupNames[0],"",xLabel=measureType,yLabel="frequency",color1='dodgerblue',color2='magenta')


    corrected_pValue_parametric=pValue_parametric
    corrected_pValue_nonparametric=pValue_nonparametric

    outputPath=outputPath+timePoints+"_"+measureType+"_groupDifference.txt"

    reportFile=open(outputPath,'w')
    reportFile.write("============== Dataset statistics ==============\n")
    reportFile.write("\tHealthy\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(healthy),np.median(scores[healthy]),scores[healthy].mean(),scores[healthy].var(),scores[healthy].std()))
    reportFile.write("\tPatients\t numSubjects:%d\tmedian:%0.3f\tmean:%0.3f\tvar:%0.3f\tstd:%0.3f \n" % (len(patients),np.median(scores[patients]),scores[patients].mean(),scores[patients].var(),scores[patients].std()))
           

    reportFile.write("============== F/Barlett/Levene Test for checking equalness of variance of patients relative to healthy controls ==============\n")
    for discardOutliers in [True,False]:
        reportFile.write("Discard outliers:"+str(discardOutliers)+"\n")
        reportFile.write("\tHealthy vs Patients\t F-test:%f\tBartlett:%f\tLevene:%f \n" % (checkDifferenceOfVariance(scores[patients],scores[healthy],'F',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patients],'bartlett',discardOutliers)[1],checkDifferenceOfVariance(scores[healthy],scores[patients],'levene',discardOutliers)[1]))


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
    outputPath = f"{outputFolder}/{sysNames[ind]}/"
    if(noHealthyPlot==True):
        numTimePoints = 1
        outputPath=outputPath+timePoints+"_"+measureType+"_"+scoreType+"_noHealthyPlot."+plotExtension
    else:
        numTimePoints = 2
        outputPath=outputPath+timePoints+"_"+measureType+"_"+scoreType+"."+plotExtension


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
        drawBoxPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=scoreName,colors=colors,rotation=0,plotScatter=True,yLim=yLim,middleLine='mean') #since we use Mann-Whitney U test for group dofference, we should plot median line in boxplots
    elif plotType=="violin":
        drawViolinPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=scoreName,colors=colors,rotation=0,yLim=yLim)