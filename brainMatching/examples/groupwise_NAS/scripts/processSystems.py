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


# load healthy & patient IDs from files
tdcPath = "../data/tdc_schaefer.txt"
asdPath = "../data/asd_schaefer.txt"

with open(tdcPath, "r") as f:
    healthyIDs = f.read().splitlines()

with open(asdPath, "r") as f:
    patientIDs = f.read().splitlines()

healthyOrder = list(range(len(healthyIDs)))
healthy = healthyOrder

# print(len(healthyIDs), len(patientIDs))
# print(healthyOrder[:])
# quit()

# Load system-level mappings
sysMapPath = "../data/yeo_7systems_schaefer.txt"
with open(sysMapPath, "r") as f:
    sysMaps = f.read().splitlines()

visualOne = []
somatomotorTwo = []
dorsalThree = []
ventralFour = []
limbicFive = []
frontoparietalSix = []
defaultmodeSeven = []
subcorticalEight = []
allSys = [visualOne, somatomotorTwo, dorsalThree, ventralFour, limbicFive, frontoparietalSix, defaultmodeSeven, subcorticalEight]
for index, mapping in enumerate(sysMaps):
    if mapping == '1':
        visualOne.append(index)
    elif mapping == '2':
        somatomotorTwo.append(index)
    elif mapping == '3':
        dorsalThree.append(index)
    elif mapping == '4':
        ventralFour.append(index)
    elif mapping == '5':
        limbicFive.append(index)
    elif mapping == '6':
        frontoparietalSix.append(index)
    elif mapping == '7':
        defaultmodeSeven.append(index)
    elif mapping == '8':
        subcorticalEight.append(index)
sysNames = ["visual", "somatomotor", "dorsal", "ventral", "limbic", "frontoparietal", "defaultmode", "subcortical"]

""" print(sum([len(mapping) for mapping in allSys]))
print(len(somatomotorTwo))
print(subcorticalEight[:]) """


##############load results of the structure-function coupling experiment###################
with open(resultFilePath,"r") as resultFile:
    fileContent = resultFile.read().splitlines()
numNodes = int(fileContent[3].split('\t')[0])
numSubjects = int(fileContent[3].split('\t')[1])

""" allSys = [frontoparietalSix]
print(len(allSys[0]))
print(allSys[0]) """

###load similarity scores and the matchings
matchings=np.array([np.fromstring(cont,dtype=int,sep='\t') for cont in fileContent[6+numSubjects:-1]])

# For each sub-system, calculate average NNS relative to healthy & save
for ind, system in enumerate(allSys):
    # create NxN pairwise-matching score matrix
    scores=np.zeros((numSubjects,numSubjects),dtype=float)
    for i in range(len(matchings)):
        row=matchings[i][0] #first column is the order number of the first subject
        col=matchings[i][1] #second column is the order number of the second subject

        for j in system:
            if(matchings[i][2+j]==j):
                scores[row][col]+=1
        """ print(scores[row][col])
        print(scores[row][col]/len(system) * 100)
        quit() """

    scores/=float(len(system))
    scores*=100

    ### calculate average matching/similarity scores relative to healthy controls
    ### NOTE: we are discarding the matching of a healthy control subject to itself in 
    ###       calculation of each subject average score relative to healthy controls
    scores_avg = np.zeros(len(scores))
    if(relativeTo=="healthy"):
        for row in range(numSubjects):
            count=0
            for col in healthy:
                if(row!=col):
                    scores_avg[row] += (scores[row][col]+scores[col][row])/2.0
                    count+=1
            scores_avg[row] /= float(count)

    scores=scores_avg


    scoreName='matching accuracy (%)'
    
    outputPath = "../experiment/results/sys_level"
    outputFile=open(f"{outputPath}/{sysNames[ind]}.res",'w')
    outputFile.write("#numNodes,numSubjects\n%d\t%d\n" % (numNodes,numSubjects))
    outputFile.write("#measureType\n"+measureType+"\n")
    outputFile.write("#scoreName\n"+scoreName+"\n")
    outputFile.write("#scores\n")
    for i in range(len(scores)):
        outputFile.write("%0.4f\t" % scores[i])
    outputFile.close()
