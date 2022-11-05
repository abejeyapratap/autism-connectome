import pandas as pd
import numpy as np
import scipy.stats as stt
import helpers as helper
import matplotlib.pyplot as plt
import matplotlib as mpl

dfPath = "../../data/cluster/connectome_clustered.csv"
df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

dfPath = "../../data/cluster/dmn_clustered.csv"
df1=pd.read_csv(dfPath)
df1 = df1.drop('Unnamed: 0', axis=1)

dfPath = "../../data/cluster/fpn_clustered.csv"
df2=pd.read_csv(dfPath)
df2 = df2.drop('Unnamed: 0', axis=1)

# connectome-level significant clusters
g0Con = df[df['cluster'] == 0]
g2Con = df[df['cluster'] == 2]

# system level sig clusters
g0Dmn = df1[df1['cluster'] == 0]
g2Dmn = df1[df1['cluster'] == 2]
g0Fpn = df2[df2['cluster'] == 0]
g2Fpn = df2[df2['cluster'] == 2]

allGroups = [g0Con, g2Con, g0Dmn, g2Dmn, g0Fpn, g2Fpn]
groupNames = ["cluster1_connectome", "cluster3_connectome", "cluster1_dmn", "cluster3_dmn", "cluster1_fpn", "cluster3_fpn"]
for ind, patientsDf in enumerate(allGroups):
    nns, severity = patientsDf['NNS'], patientsDf['ados']

    # Calculate correlation values
    r_pearson, p_pearson = stt.pearsonr(nns, severity)
    print(r_pearson, p_pearson)

    if p_pearson < 0.0005:
        text = 'r='+str("{0:.2f}".format(r_pearson))+'\np<5e-3'
    else:
        text = ""

    if groupNames[ind] == "cluster1":
        colors = ["dodgerblue", "#00b33c"]
    else:
        colors = ["dodgerblue", "#00b33c"]
    colors = ["#351C4D", "darkorchid"]

    cutoff = 0.05
    if p_pearson <= cutoff:
        outputFile = f"{groupNames[ind]}.svg"
        helper.drawCorrelationPlot(nns, severity, r_pearson, p_pearson, "NNS","ADOS","", outputFile, text=text, dotColor=colors[0], lineColor=colors[1])
    