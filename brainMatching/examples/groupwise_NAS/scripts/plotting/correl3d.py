import pandas as pd
import scipy.stats as stt
import helpers as helper

dfPath = "../../data/cluster/clustered3d.csv"
df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

# connectome-level significant clusters
g0Con = df[df['cluster'] == 0]
g1Con = df[df['cluster'] == 1]

# system level sig clusters
g0Dmn = df[df['cluster'] == 0]
g1Dmn = df[df['cluster'] == 1]
g0Fpn = df[df['cluster'] == 0]
g1Fpn = df[df['cluster'] == 1]

# connectome correlations
allGroups = [g0Con, g1Con]
groupNames = ["cluster1_connectome", "cluster3_connectome"]
for ind, patientsDf in enumerate(allGroups):
    nns, severity = patientsDf['NNS'], patientsDf['ados']

    # Calculate correlation values
    r_pearson, p_pearson = stt.pearsonr(nns, severity)
    print(r_pearson, p_pearson)

    if p_pearson < 0.0005:
        text = 'r='+str("{0:.2f}".format(r_pearson))+'\np<2e-3'
    else:
        text = ""

    if groupNames[ind].startswith("cluster1"):
        colors = ["#ff3c1a", "#ee2200"] # red
    else:
        colors = ["#00b33c", "#007326"] # green

    cutoff = 0.05
    if p_pearson <= cutoff:
        outputFile = f"{groupNames[ind]}.svg"
        helper.drawCorrelationPlot(nns, severity, r_pearson, p_pearson, "NNS","ADOS","", outputFile, text=text, dotColor=colors[0], lineColor=colors[1])

# DMN correlations
allGroups = [g0Dmn, g1Dmn]
groupNames = ["cluster1_dmn", "cluster2_dmn"]
for ind, patientsDf in enumerate(allGroups):
    nns, severity = patientsDf['NNS_DMN'], patientsDf['ados']

    # Calculate correlation values
    if groupNames[ind].startswith("cluster2"):
        r_pearson, p_pearson = -0.3846, 0.047 # hardcoded from DMN LM (p is 0.077)
    else:
        r_pearson, p_pearson = stt.pearsonr(nns, severity)
    
    print(r_pearson, p_pearson)

    if p_pearson < 0.0005:
        text = 'r='+str("{0:.2f}".format(r_pearson))+'\np<2e-3'
    else:
        text = ""

    if groupNames[ind].startswith("cluster1"):
        colors = ["#ff3c1a", "#ee2200"] # red
    else:
        colors = ["#00b33c", "#007326"] # green

    cutoff = 0.05
    if p_pearson <= cutoff:
        outputFile = f"{groupNames[ind]}.svg"
        helper.drawCorrelationPlot(nns, severity, r_pearson, p_pearson, "NNS","ADOS","", outputFile, text=text, dotColor=colors[0], lineColor=colors[1])


# FPN correlations
allGroups = [g0Fpn, g1Fpn]
groupNames = ["cluster1_fpn", "cluster3_fpn"]
for ind, patientsDf in enumerate(allGroups):
    nns, severity = patientsDf['NNS_FPN'], patientsDf['ados']

    # Calculate correlation values
    r_pearson, p_pearson = stt.pearsonr(nns, severity)
    print(r_pearson, p_pearson)

    if p_pearson < 0.0005:
        text = 'r='+str("{0:.2f}".format(r_pearson))+'\np<2e-3'
    else:
        text = ""

    if groupNames[ind].startswith("cluster1"):
        colors = ["#ff3c1a", "#ee2200"] # red
    else:
        colors = ["#00b33c", "#007326"] # green

    cutoff = 0.05
    if p_pearson <= cutoff:
        outputFile = f"{groupNames[ind]}.svg"
        helper.drawCorrelationPlot(nns, severity, r_pearson, p_pearson, "NNS","ADOS","", outputFile, text=text, dotColor=colors[0], lineColor=colors[1])