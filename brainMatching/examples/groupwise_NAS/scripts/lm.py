import numpy as np
import pandas as pd
import scipy.stats as stt
from statsmodels.formula.api import ols

dfPath = "../data/lm_clustered.csv"

df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

g0 = df[df['cluster'] == 0]
g1 = df[df['cluster'] == 1]
g2 = df[df['cluster'] == 2]

allDf = [df, g0, g1, g2]
labels = ["All Patients", "G0 cluster", "G1 cluster", "G2 cluster"]

""" model_all = ols("NNS ~ ados + Age", data=df)
results_all = model_all.fit()
f.write('### ols("NNS ~ ados + Age", data=df) ###\n')
f.write(results_all.summary())
f.write("\n\n") """

rString = "ados ~ NNS + Age"
rString = "NNS ~ ados + Age"

with open("ols_res.txt", "w") as f:
    for i, df1 in enumerate(allDf):
        f.write(f"\t\t\t{labels[i]}\n")
        model = ols(rString, data=df1)
        results = model.fit()
        f.write(f'### Function call: ols("{rString}", data=df1) ###\n')
        f.write(results.summary().as_text())
        f.write("\n\n\n\n")