from operator import le
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#%%

def heatmap_process(data):
    df = pd.DataFrame({"value": [0.0]*64, "Outside Vehicles-Remaining Capacities": [""]*64, "Dynamic Events":[""]*64})
    df["value"] = np.mean(data, axis=1)

    with open("acronyms","r") as f:
        f.readline()
        for i in range(64):
            s = f.readline().strip().split(" ")
            OVRQ, dynamic = s[1], s[2]
            df.at[i, "Outside Vehicles-Remaining Capacities"] = OVRQ
            df.at[i, "Dynamic Events"] = dynamic
    data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "value")
    return data
    
    # plt.savefig(folder+'/usedTime.eps',bbox_inches = "tight", dpi=300,format="eps")
    # plt.show()



# %%
maxruns = 10

ratios = np.zeros((64, 10))
corre = np.zeros((64, 10))
for mapidx in range(maxruns):
    folder = "analysis/lor/egl-e2-A-{0}/".format(mapidx)
    for i in range(1, 65):
        filepath = folder+"instance{0}.txt".format(i)
        with open(filepath, "r") as f:
            for _ in range(8):
                f.readline()
            data = f.readline().split(",")
            dis_random, dis_optima = int(data[2]), int(data[4])
            ratios[i-1, mapidx] = (dis_optima-dis_random)/dis_random

            data = f.readline().split(",")
            data.pop(0)
            data.pop(-1)
            dis_list = list(map(lambda x:int(x), data))

            data = f.readline().split(",")
            data.pop(0)
            data.pop(-1)
            cost_list = list(map(lambda x:int(x), data))
            tmp = np.corrcoef(dis_list, cost_list)
            corre[i-1, mapidx] = tmp[0,1]
            print(tmp[0,1])
        # print(dis_random, dis_optima)
        # print(dis_list)
        # break
    # break

#%%
folder = "analysis/lor"
data = heatmap_process(ratios)
g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True, n_colors=2), center=0, linewidths=.1, fmt=".3f", annot=True, vmin=-0.2, vmax=0.2)
g.set(title="Ratio of Distances between Two Solutions")
plt.savefig(folder+'/dis_ratio.eps',bbox_inches = "tight", dpi=300,format="eps")

#%%
data = heatmap_process(corre)
g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True, reverse=True), center=0, linewidths=.1, fmt=".3f", annot=True, vmin=-1, vmax=1)
g.set(title="Correlation of Distances between Costs")
plt.savefig(folder+'/dis_correlation.eps',bbox_inches = "tight", dpi=300,format="eps")
# %%
returnp = np.zeros((64, 10))
for mapidx in range(maxruns):
    folder = "analysis/returnp/egl-e2-A-{0}/".format(mapidx)
    for i in range(1, 65):
        filepath = folder+"instance{0}.txt".format(i)
        with open(filepath, "r") as f:
            data = f.readline().split(",")
            if len(data) == 0:
                print(filepath)
            p = int(data[0])/1000
            returnp[i-1, mapidx] = p

df = pd.DataFrame({"value": [0.0]*64, "Outside Vehicles-Remaining Capacities": [""]*64, "Dynamic Events":[""]*64})
df["value"] = np.mean(returnp, axis=1)

with open("acronyms","r") as f:
    f.readline()
    for i in range(64):
        s = f.readline().strip().split(" ")
        OVRQ, dynamic = s[1], s[2]
        df.at[i, "Outside Vehicles-Remaining Capacities"] = OVRQ
        df.at[i, "Dynamic Events"] = dynamic
data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "value")
g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True, reverse=True), linewidths=.1, fmt=".3f", annot=True, vmin=0, vmax=0.3)
g.set(title="Return Probability of Global Optimum")
plt.savefig('analysis/lor/returnp.eps',bbox_inches = "tight", dpi=300,format="eps")
# %%
