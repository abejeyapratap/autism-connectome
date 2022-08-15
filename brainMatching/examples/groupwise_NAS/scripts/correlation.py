import pandas as pd
import numpy as np
import scipy.stats as stt
from helpers import drawCorrelationPlot, drawBoxPlot

""" tdcPath = "../data/tdc_desikan.txt"
asdPath = "../data/asd_desikan.txt"

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
patientGroupNames = ["ASD"] """

### Load NAS results
""" # resultFile = "../experiment/results/NAS_rt_healthy.res"
resultFile = "../experiment/results_desikan_all/results_norm_desikan/NAS_rt_healthy.res"
fileContent =  open(resultFile,"r").read().splitlines()
numNodes = int(fileContent[1].split('\t')[0])
numSubjects = int(fileContent[1].split('\t')[1])
measureType = str(fileContent[3]) 
scoreName = str(fileContent[5]) 
scores = np.zeros(numSubjects)

# load matching accuracy scores
for i in range(numSubjects):
    scores[i] = float(fileContent[7].split('\t')[i]) """


path = "../data/ados/schaefer_all_ados.csv"
# path = "../data/ados/g6_schaefer_ados.csv"
patientsDf = pd.read_csv(path)
similarity, ados = patientsDf['NNS'], patientsDf['ados_css']

# Calculate correlation values
r_pearson, p_pearson = stt.pearsonr(similarity, ados)
r_spearman, p_spearman = stt.spearmanr(similarity, ados)

# save Scatter Plot
outputFolder = "../experiment/plots/correl/"
outputPath=outputFolder+"schaefer_all.png"

plot=drawCorrelationPlot(similarity, ados, r_pearson, p_pearson, "NNS","ADOS","", outputPath)
plot.close()

### BOXPLOT ### 
""" 
colors=['#351C4D', '#AB3E16','#849974','#2096BA','#F7DFD4','#F5AB99'] #nightfall, rust, fresh, shutter blue, macaron, tropical pink
# data = [scores[patients], scores[healthy]]
# dataLabels=['Patient','Healthy']
data = [scores[patients]]
dataLabels=['Patient']

minY=min([min(l) for l in data])
maxY=max([max(l) for l in data])
offset=(maxY-minY)/5.0

yLim=[minY-offset/2.0,maxY+offset] ## make space specific to figure
# yLim=[67,102] ## use these lines to make the space constant (such as across different plots)
scoreName="network similarity (%)"
title = ""
drawBoxPlot(data,dataLabels,title,outputPath,xLabel='',yLabel=scoreName,colors=colors,rotation=0,plotScatter=True,yLim=yLim,middleLine='median') #Mann-Whitney U test so plot median line in boxplots """