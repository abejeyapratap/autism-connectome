import pandas as pd
import numpy as np
import scipy.stats as stt
from helpers import drawCorrelationPlot, drawBoxPlot, calculateGroupDifference
import argparse
import os

parser = argparse.ArgumentParser(description='Calculate system-level correlations between NNS and ADOS/SCQ')
parser.add_argument('-sip','--subjectsInfoPath', help='path to the file that contains path to the connectomes', required=True)
parser.add_argument('-r','--resultFile', help='file path to the results of the matching experiment', required=True)
parser.add_argument('-o','--outputFile', help='file path to save the distribution of values for the two populations', required=True)

args = vars(parser.parse_args())
subjectsInfoPath=args['subjectsInfoPath']
resultFilePath=args['resultFile']
outputPath=args['outputFile']


# Load Tobacco data & processed NNS Scores
df = pd.read_csv(subjectsInfoPath)
filtered = df[['Subject', 'DX', 'ados_css', 'scq_total']]

# match patients with corresponding ADOS score
patients = []
fullP = []
ados = []
for i, val in enumerate(filtered['DX']):
  if val == 'ASD':
    patients.append(filtered['Subject'][i])
    fullP.append(f"{filtered['Subject'][i]}_DTI_Schaefer2018_200_7Networks_connmat_sift2")
    adosScore = filtered['ados_css'][i]
    if np.isnan(adosScore):
      ados.append(None)
    else:
      ados.append(adosScore)
""" print(sum([x is None for x in ados]))
print(len(patients), len(ados)) # 176 """

# Eliminate missing subjects
with open('../data/missing.txt', 'r') as f:
  missingSubs = f.read().splitlines()

# remove missing patients & ADOS
filteredP = []
filteredAdos = []

for i, p in enumerate(fullP):
  if p not in missingSubs:
    filteredP.append(patients[i])
    filteredAdos.append(ados[i])


sysNames = ["defaultmode", "frontoparietal", "visual", "somatomotor", "dorsal", "ventral", "limbic", "subcortical"]

# create output folders for each system
outputFolder = outputPath
for system in sysNames:
    sysPath = f"{outputFolder}/{system}"
    if not os.path.exists(sysPath):
        os.mkdir(sysPath)

systemMatchingPath = resultFilePath
for ind, system in enumerate(sysNames):
    resultFile = f"{systemMatchingPath}/{system}.res"

    fileContent =  open(resultFile,"r").read().splitlines()
    numSubjects = int(fileContent[1].split('\t')[1])
    scores = np.zeros(numSubjects)

    # load matching accuracy scores
    for i in range(numSubjects):
        scores[i] = float(fileContent[7].split('\t')[i])
    patientScores = scores[150:]


    # Create DataFrame with Patients, NNS, ADOS
    data = {'Subject':  filteredP, 'NNS': patientScores, 'ados_css': filteredAdos}
    newDf = pd.DataFrame(data)
    """ print(newDf.shape, len(patientScores))
    print(newDf[newDf['ados_css'].isnull()]) # check NaN ADOS vals
    # print(newDf.head()) """

    df = newDf[newDf['ados_css'].notna()]

    ### Split into sub-groups (3 or 5)
    # delimiting boundaries
    q1, q3 = np.percentile(df['NNS'], 25), np.percentile(df['NNS'], 75)
    iqr = q3-q1
    median = np.percentile(df['NNS'], 50)
    lowWhisker = q1 - (1.5*iqr)
    highWhisker = q3 + (1.5*iqr)
    # print(q1, q3, iqr, lowWhisker, highWhisker)

    
    g1 = df[df['NNS'] < lowWhisker] # Group 1 - all outliers
    # print(g1.shape, stt.pearsonr(g1['NNS'], g1['ados_css']))

    g2 = df[(df['NNS'] >= lowWhisker) & (df['NNS'] < q1)] # low whisker to 1st quartile
    g3 = df[df['NNS'] >= q1] # 1st quartile to rest

    ### other sub-groups    
    c4 = df[df['NNS'] >= lowWhisker] # Low whisker & up
    c5 = df[(df['NNS'] >= q3) & (df['NNS'] < highWhisker)] # Quartile 3 & up
    g4 = df[(df['NNS'] >= median) & (df['NNS'] < q3)] # b/w median & Q3
    g6 = df[df['NNS'] >= median] # median & up

    outFolder = f"{outputPath}/{system}"
    reportFile=open(f"{outFolder}/stats.txt",'w')
    reportFile.write(f"<<<{system}>>>\n")
    reportFile.write("============== Significant Correlations ==============\n")

    # Calculate and plot correlations for sub-groups
    groupNames = ["g1", "g2", "g3_quartile1", "c5_quartile3", "g6_median", "c4_lowwhisker", "g4_median_q3", "all_patients"]
    allGroups = [g1, g2, g3, c5, g6, c4, g4, df]
    for ind, patientsDf in enumerate(allGroups):
        nns, ados = patientsDf['NNS'], patientsDf['ados_css']

        # Calculate correlation values
        r_pearson, p_pearson = stt.pearsonr(nns, ados)
        r_spearman, p_spearman = stt.spearmanr(nns, ados)

        if p_pearson <= 0.05:
            reportFile.write(f"\t {groupNames[ind]}\t Pearson r-value:{r_pearson:.3f}\t p-value:{p_pearson:.5f}\n")
            outputFile = f"{outFolder}/{groupNames[ind]}.png"
            plot=drawCorrelationPlot(nns, ados, r_pearson, p_pearson, "NNS","ADOS","", outputFile)
            plot.close()
        elif p_spearman <= 0.05:
            reportFile.write(f"\t {groupNames[ind]}\t Spearman r-value:{r_spearman:.3f}\t p-value:{p_spearman:.5f}\n")
            outputFile = f"{outFolder}/spearman_{groupNames[ind]}.png"
            plot=drawCorrelationPlot(nns, ados, r_spearman, p_spearman, "NNS","ADOS","", outputFile)
            plot.close()


    reportFile.write("\n============== Group Difference between NNS of low ADOS vs high ADOS patients ==============\n")
    ### Group Difference between NNS of low ADOS vs high ADOS patients
    # see if NNS scores of ADOS <=5 is sig diff from ADOS >5
    lowAdos = df[df['ados_css'] <= 5]
    highAdos = df[df['ados_css'] > 5]

    paramDiff = calculateGroupDifference(lowAdos['NNS'], highAdos['NNS'], parametric=False, paired=False) # significant
    nonparamDiff = calculateGroupDifference(lowAdos['NNS'], highAdos['NNS'], parametric=True, paired=False)

    reportFile.write("============== Parametric test (with Cohen's D) ==============\n")
    if paramDiff[1] <= 0.05:
        reportFile.write("***")
    reportFile.write(f"\t Effect Size:{paramDiff[0]:.3f}\t p-value:{paramDiff[1]:.5f}\n")
    reportFile.write("============== Non-parametric test (Mann-Whitney U test (not paired) for repeated measures with non normal distribution) ==============\n")
    if nonparamDiff[1] <= 0.05:
        reportFile.write("***")
    reportFile.write(f"\t Effect Size:{nonparamDiff[0]:.3f}\t p-value:{nonparamDiff[1]:.5f}\n")

    dat = [lowAdos['NNS'].tolist(), highAdos['NNS'].tolist()]
    lab = ['<=5', '>5']
    minY=min([min(l) for l in dat])
    maxY=max([max(l) for l in dat])
    offset=(maxY-minY)/5.0
    yLim=[minY-offset/2.0,maxY+offset]
    colors=['#351C4D', '#AB3E16','#849974','#2096BA','#F7DFD4','#F5AB99']
    drawBoxPlot(dat,lab,"",f"{outFolder}/NNS_box.png",xLabel='',yLabel="NNS", colors=colors, rotation=0,plotScatter=True,yLim=yLim,middleLine='median')

    reportFile.close()