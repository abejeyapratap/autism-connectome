import pandas as pd
import numpy as np
import helpers as helper
# import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# cluster plot
dfPath = "../../data/cluster/connectome_clustered.csv"
df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

u_labels = np.unique(df['cluster'])
labels = df['cluster']
patients = np.array(df["NNS"]).reshape(-1,1)
ados = np.array(df["ados"]).reshape(-1,1)

mpl.rcParams['font.sans-serif'] = "Times New Roman"
mpl.rcParams['font.family'] = "serif"

fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xlabel("Network Normality Score (%)",fontsize=12)
ax.set_ylabel("Autism Severity (ADOS)",fontsize=12)
ax.tick_params(axis='both',labelsize=10)

colors = ["limegreen", "deepskyblue", "tomato"]
colors = ["dodgerblue", "tomato", "#00b33c"]
for i in u_labels:
    plt.scatter(patients[labels==i, 0], ados[labels==i, 0], c=colors[i], alpha=0.9)

# plt.show()
outPath = "../../experiment/results_main/plots/cluster.svg"
plt.savefig(outPath, transparent=False, dpi=300)

# sns.lmplot(data=df, x="NNS", y="ados", hue="cluster", fit_reg=False)