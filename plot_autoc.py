
# %% import libarays
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.io.formats import style
import seaborn as sns

from lxml import etree as ET

sns.set_theme(style="whitegrid")

folder = "analysis/autoc"

# %% auto-correlation congestion
group_cgst = [0, 1, 2, 3, 4, 5, 6, 7, 8]
static = ["static", ]
OVRQ = ["few-low", "few-high", "many-low", "many-high"]
cgst = ["20%roads", "80%roads"]
ease = ["20%roads", "80%roads"]

# congestion factor
df = pd.DataFrame({"autoc": [], "step": [], "OV-RQ": [], "Congestion": []})
for i in group_cgst:
    filepath = folder+"/instance{0}.txt".format(i)
    with open(filepath, "r") as f:
        data = list(map(lambda x: float(x), f.readline().split(";")))
    corrlen = data.pop(-1)

    data = np.log(data)
    if i == 0:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": ["static"]*101, "Congestion": ["static"]*101}
    else:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": [OVRQ[(i-1)//2]]*101, "Congestion": [cgst[(i-1) % 2]]*101}
        print(OVRQ[(i-1)//2], cgst[(i-1) % 2], data[-1])

    dfTemp = pd.DataFrame(data=d)
    df = df.append(dfTemp)

    # break
ax = sns.lineplot(data=df, x="step", y="autoc", hue="OV-RQ",
                  style="Congestion", palette="tab10", linewidth=2.5)
ax.set(xlabel=r'tau $\tau$', ylabel='ln(Auto-Correlation)',
       title='Autocorrelation RQ-OV~Congestion')
plt.savefig('analysis/autoc1.eps', dpi=600, format="eps")


# %% auto-correlation roads ease
group_ease = [0, 9, 10, 11, 12, 13, 14, 15, 16]
static = ["static"]
OVRQ = ["few-low", "few-high", "many-low", "many-high"]
ease = ["20%roads", "80%roads"]
# congestion factor
df = pd.DataFrame({"autoc": [], "step": [], "OV-RQ": [], "Ease": []})
for i in group_ease:
    filepath = folder+"/instance{0}.txt".format(i)
    print(filepath)
    with open(filepath, "r") as f:
        data = list(map(lambda x: float(x), f.readline().split(";")))
    corrlen = data.pop(-1)

    data = np.log(data)
    if i == 0:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": ["static"]*101, "Ease": ["static"]*101}
    else:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": [OVRQ[(i-9)//2]]*101, "Ease": [ease[(i-9) % 2]]*101}
        print(OVRQ[(i-9)//2], ease[(i-9) % 2], data[-1])

    dfTemp = pd.DataFrame(data=d)
    df = df.append(dfTemp)

    # break
ax = sns.lineplot(data=df, x="step", y="autoc", hue="OV-RQ",
                  style="Ease", palette="tab10", linewidth=2.5)
ax.set(xlabel=r'tau $\tau$', ylabel='ln(Auto-Correlation)',
       title='Autocorrelation RQ-OV~Ease')
# legend = plt.legend()
# legend.get_frame().set_alpha(None)
# legend.get_frame().set_facecolor((1, 1, 1, 0.1))
plt.savefig('analysis/autoc2.eps', dpi=600, format="eps",transparent=True)

# %% auto-correlation new demand
group_new_demand = [0, 17, 18, 19, 20, 21, 22,
                    23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
static = ["static"]
OVRQ = ["few-low", "few-high", "many-low", "many-high"]
new_demand = ["few-low", "few-high", "many-low", "many-high"]  # num-value
# congestion factor
df = pd.DataFrame({"autoc": [], "step": [], "OV-RQ": [], "AddDemand": []})
for i in group_new_demand:
    filepath = folder+"/instance{0}.txt".format(i)
    print(filepath)
    with open(filepath, "r") as f:
        data = list(map(lambda x: float(x), f.readline().split(";")))
    corrlen = data.pop(-1)

    data = np.log(data)
    if i == 0:
        for k in range(4):
            d = {"autoc": data, "step": np.arange(
                101), "OV-RQ": [OVRQ[k]]*101, "AddDemand": ["static"]*101}
            dfTemp = pd.DataFrame(data=d)
            df = df.append(dfTemp)
    else:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": [OVRQ[(i-17)//4]]*101, "AddDemand": [new_demand[(i-17) % 4]]*101}
        print(OVRQ[(i-17)//4], new_demand[(i-17) % 4], data[-1])

        dfTemp = pd.DataFrame(data=d)
        df = df.append(dfTemp)

setting = {'color': ["blue", "orange", "green", "red"]}
g = sns.FacetGrid(df, col="OV-RQ")
g.map_dataframe(sns.lineplot, "step", "autoc", hue="AddDemand", style="AddDemand", linewidth=2.5)
g.add_legend(title="Add demand: Num-Value")
g.set_xlabels(r"tau $\tau$")
plt.savefig('analysis/autoc3.eps', dpi=600, format="eps")

# %%
group_new_task = [0, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64]
static = ["static"]
OVRQ = ["few-low", "few-high", "many-low", "many-high"]
new_task = ["few-low", "few-high", "many-low", "many-high"]
new_task_pos = ["near", "far"]
# congestion factor
df = pd.DataFrame({"autoc": [], "step": [], "OV-RQ": [], "AddTask": [], "Pos":[]})
for i in group_new_task:
    filepath = folder+"/instance{0}.txt".format(i)
    with open(filepath, "r") as f:
        data = list(map(lambda x: float(x), f.readline().split(";")))
    corrlen = data.pop(-1)

    data = np.log(data)
    if i == 0:
        for k in range(4):
            for kk in range(2):
                d = {"autoc": data, "step": np.arange(
                    101), "OV-RQ": [OVRQ[k]]*101, "AddTask": ["static"]*101, "Pos":[new_task_pos[kk]]*101,}
                dfTemp = pd.DataFrame(data=d)
                df = df.append(dfTemp)
    else:
        d = {"autoc": data, "step": np.arange(
            101), "OV-RQ": [OVRQ[(i-33)//8]]*101, "AddTask": [new_task[(i-33) % 4]]*101, "Pos": new_task_pos[(i-33)//4%2]}
        print(OVRQ[(i-33)//8], new_task[(i-33) % 4], new_task_pos[((i-33)//4)%2])

        dfTemp = pd.DataFrame(data=d)
        df = df.append(dfTemp)

setting = {'color': ["blue", "orange", "green", "red"]}
g = sns.FacetGrid(df, row="Pos", col="OV-RQ")
g.map_dataframe(sns.lineplot, "step", "autoc", hue="AddTask", style="AddTask", linewidth=2.5)
g.add_legend(title="Add Task: Num-Value")
g.set_xlabels(r"tau $\tau$")
plt.savefig('analysis/autoc4.eps', dpi=600, format="eps")

# %% correlation length


corrlens = [0]
task = [0]
for i in range(1, 65):
    tmp = np.zeros((10, ))
    tmp1 = np.zeros((10, ))
    for run in range(10):
        filepath = folder+"/egl-e2-A-{0}/instance{1}.txt".format(run, i)
        with open(filepath, "r") as f:
            data = list(map(lambda x: float(x), f.readline().split(";")))
            tmp[run] = data.pop(-1)
        

        xml_path = "instance/egl-e2-A-{0}/instance{1}.xml".format(run, i)
        tree = ET.ElementTree(file=xml_path)
        root = tree.getroot()
        task_num = int(root[1].attrib["num"])//2
        tmp1[run] = task_num


    corrlens.append(np.mean(tmp))
    task.append(np.mean(tmp1))



df = pd.DataFrame({"Correlation Length": [], "Outside Vehicles-Remaining Capacities": [], "Dynamic Events": [], "Number of Tasks":[]})
with open("acronyms","r") as f:
    f.readline()
    for i in range(1, 65):
        s = f.readline().strip().split(" ")
        
        # print(s)

        # xml_path = "instance/egl-e2-A-0/instance{0}.xml".format(i)
        # tree = ET.ElementTree(file=xml_path)
        # root = tree.getroot()
        # task_num = int(root[1].attrib["num"])//2
        # print("task num: ", task_num)

        task_num = task[i]

        d = {"Correlation Length": [corrlens[i]], "Outside Vehicles-Remaining Capacities": [s[1]], "Dynamic Events": [s[2]], "Number of Tasks":[task_num]}
        dfTemp = pd.DataFrame(data=d)
        df = df.append(dfTemp)


# data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "Correlation Length")
# ax = sns.heatmap(data, center=0, linewidths=.1, fmt=".1f", annot=True)
# ax.set(title='Heatmap of Correlation Length')

ax1 = sns.regplot(x = "Number of Tasks", y="Correlation Length", data=df, ci=0)
plt.savefig('analysis/autoc/corrlen_num.eps', bbox_inches = "tight", dpi=300, format="eps")

# plt.savefig('analysis/autoc/autoclength.eps', bbox_inches = "tight", dpi=300, format="eps")

# ax1 = sns.regplot(x="N", y="corrlength", data=df)
# ax1.set_xlabel("The Number of Tasks")
# ax1.set_ylabel("Relative Correlation Length")
# plt.savefig('analysis/N_corlen.png',bbox_inches = "tight", dpi=600)


# %%






# %%
