import numpy as np
import pandas as pd
import scipy.stats as stt
from statsmodels.formula.api import ols
from sklearn.preprocessing import MinMaxScaler

level = "fpn" # connectome, dmn, fpn
# dfPath = f"../data/cluster/{level}_clustered.csv"
dfPath = f"../data/cluster/subgroup_clustered.csv" # has IQ & clustered based on 3 NNS

df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

# print(df.head())

# Normalize 3 data columns for correct coefficients
df["Age"] = MinMaxScaler().fit_transform(np.array(df["Age"]).reshape(-1,1))
df["ados"] = MinMaxScaler().fit_transform(np.array(df["ados"]).reshape(-1,1))
df["NNS"] = MinMaxScaler().fit_transform(np.array(df["NNS"]).reshape(-1,1))
df["NNS_DMN"] = MinMaxScaler().fit_transform(np.array(df["NNS_DMN"]).reshape(-1,1))
df["NNS_FPN"] = MinMaxScaler().fit_transform(np.array(df["NNS_FPN"]).reshape(-1,1))

# print(df.head())

g0 = df[df['cluster'] == 0]
g1 = df[df['cluster'] == 1]
# g2 = df[df['cluster'] == 2]

""" allDf = [df, g0, g1, g2]
labels = ["All Patients", "G0 cluster", "G1 cluster", "G2 cluster"] """
allDf = [df, g0, g1]
labels = ["All Patients", "G0 cluster", "G1 cluster"]

# rString = "NNS ~ ados + Age" // Age doesn't seem to be related to NNS
middlePart = "_FPN " # blank, _DMN, _FPN
rString = "ados ~ NNS" + middlePart + "+ Age"

outFolder = "../experiment/results_main/lm/" + level
with open(f"{outFolder}/{level}_ols.txt", "w") as f:
    for i, df1 in enumerate(allDf):
        f.write(f"\t\t\t{labels[i]}\n")
        model = ols(rString, data=df1)
        results = model.fit()
        f.write(f'### Function call: ols("{rString}", data=df1) ###\n')
        f.write(results.summary().as_text())
        f.write("\n\n\n\n")