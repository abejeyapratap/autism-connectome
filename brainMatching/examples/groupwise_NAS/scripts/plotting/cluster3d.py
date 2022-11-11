# Cluster all 3 NNS for ISBI

import pandas as pd
import numpy as np
import helpers as helper
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

# cluster plot
dfPath = "../../data/cluster/clustered3d.csv"
df=pd.read_csv(dfPath)
df = df.drop('Unnamed: 0', axis=1)

# set font stuff here
mpl.rcParams['font.sans-serif'] = "Times New Roman"
mpl.rcParams['font.family'] = "serif"
#mpl.rcParams['legend.fontsize'] = 15

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(projection='3d')

# makes background darker
ax.w_xaxis.set_pane_color((0.8, 0.8, 0.8, 0.35))
ax.w_yaxis.set_pane_color((0.8, 0.8, 0.8, 0.35))
ax.w_zaxis.set_pane_color((0.8, 0.8, 0.8, 0.35))

ax.azim = -195 # -85, -100, -200, 275, -220 (seems good)
ax.dist = 10
ax.elev = 25 # was 25

sVal = 45
mark = "o"
ax.scatter3D(df['NNS'][df['cluster']==0], df['NNS_DMN'][df['cluster']==0], df['NNS_FPN'][df['cluster']==0], s=sVal, c="#ff3c1a", marker=mark)
ax.scatter3D(df['NNS'][df['cluster']==1], df['NNS_DMN'][df['cluster']==1], df['NNS_FPN'][df['cluster']==1], s=sVal, c="#00b33c", marker=mark)

#ax.set_xlabel("NNS (%)")
ax.legend(["Red", "Green"])

ax.xaxis.set_tick_params(labelsize=10)
ax.yaxis.set_tick_params(labelsize=10)
ax.zaxis.set_tick_params(labelsize=10)

# plt.show()
plt.savefig("clustered3d.svg", transparent=False, dpi=600)