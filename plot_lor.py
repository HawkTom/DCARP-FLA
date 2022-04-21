# %% import libarays
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from lxml import etree as ET
# %% determine the folder

folder = "analysis/lorHy"
setting = {}
setting[0] = {"OVRQ": "static", "factor": "static", "factor-setting": "static","N":72}
with open("settings","r") as f:
    f.readline()
    for i in range(1, 65):
        s = f.readline().split(" ")
        ss = s[2].strip().split(":")
        setting[int(s[0])] = {"OVRQ": s[1], "factor": ss[0], "factor-setting": ss[1]}

        xml_path = "instance/egl-e2-A-0/instance{0}.xml".format(i)
        tree = ET.ElementTree(file=xml_path)
        root = tree.getroot()
        task_num = int(root[1].attrib["num"])//2
        setting[i]["N"] = task_num
        # print("task num: ", task_num)

        # break
print(setting)

bestCost = {}
addCost = {}
worseCost = {}
for i in range(65):
    with open("result/best/egl-e2-A-0/instance{0}.txt".format(i), "r") as f:
        bestCost[i] = int(f.readline())
        f.readline()
        addCost[i]=int(f.readline().split(":")[1])
        worseCost[i]=int(f.readline().split(":")[1])

pdf0 = pd.read_csv("{0}/instance0.csv".format(folder))
pdf0.pop("Unnamed: 0")
pdf0["dis"] = pdf0["dis"]/(setting[0]["N"]*2)
#%%

colors = ["#DB8458", "#57A669","#C15156","#8072B0"]
lstyles = ["--", ":"]

group_cgst = [1, 2, 3, 4, 5, 6, 7, 8]
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6,4))
g = axes
sns.kdeplot(data=pdf0, x="cost", bw_adjust=2, fill=False, color="blue", linewidth=2.5, ax=g)
for k in range(8):
    i = group_cgst[k]
    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    OVRQ = setting[i]['OVRQ']
    if i == 0:
        dynamic = "static"
        cost0 = pdf["cost"].mean()
    else:
        dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
    color = colors[k//2]
    lstyle = lstyles[k%2]
    print(i, OVRQ, color, lstyle)
    sns.kdeplot(data=pdf, x="cost", bw_adjust=2, color=color, linestyle=lstyle, linewidth=2.5, fill=False, ax=g, alpha=1)
    # break
# for l in g.get_lines():
#     print(l.get_data()[0].shape)

h = [Line2D([0], [0], color="blue", label='staitc'),
Line2D([0], [0], color="#DB8458", label='few-low'),
Line2D([0], [0], color="#57A669", label='few-high'),
Line2D([0], [0], color="#C15156", label='many-low'),
Line2D([0], [0], color="#8072B0", label='many-high')]
h1 = [Line2D([0], [0], label='static', color='k'),
Line2D([0], [0], label='20%roads', color='k', linestyle="--"),
Line2D([0], [0], label='20%roads', color='k', linestyle=":")]
legend1 = plt.legend(handles=h, title="OV-RQ", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
legend2 = plt.legend(handles=h1, title="Congest", bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0., fontsize=10)
axes.add_artist(legend2)
axes.add_artist(legend1)
g.set(title="Cost Distribution")
plt.savefig(folder+"/costPDFCgst", bbox_inches = "tight", dpi=600)
plt.show()

# %%
colors = ["#DB8458", "#57A669","#C15156","#8072B0"]
lstyles = ["--", ":"]

group_ease = list(range(9,17))
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(6,4))
g = axes
sns.kdeplot(data=pdf0, x="cost", bw_adjust=2, fill=False, color="blue", linewidth=2.5, ax=g)
for k in range(8):
    i = group_ease[k]
    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    OVRQ = setting[i]['OVRQ']
    if i == 0:
        dynamic = "static"
        cost0 = pdf["cost"].mean()
    else:
        dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
    color = colors[k//2]
    lstyle = lstyles[k%2]
    print(i, OVRQ, color, lstyle)
    sns.kdeplot(data=pdf, x="cost", bw_adjust=2, color=color, linestyle=lstyle, linewidth=2.5, fill=False, ax=g, alpha=1)
    # break

h = [Line2D([0], [0], color="blue", label='staitc'),
Line2D([0], [0], color="#DB8458", label='few-low'),
Line2D([0], [0], color="#57A669", label='few-high'),
Line2D([0], [0], color="#C15156", label='many-low'),
Line2D([0], [0], color="#8072B0", label='many-high')]
h1 = [Line2D([0], [0], label='static', color='k'),
Line2D([0], [0], label='20%roads', color='k', linestyle="--"),
Line2D([0], [0], label='20%roads', color='k', linestyle=":")]
legend1 = plt.legend(handles=h, title="OV-RQ", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=10)
legend2 = plt.legend(handles=h1, title="Ease", bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0., fontsize=10)
axes.add_artist(legend2)
axes.add_artist(legend1)
g.set(title="Cost Distribution")
plt.savefig(folder+"/costPDFEase", bbox_inches = "tight", dpi=600)
plt.show()

# %%

colors = ["#DB8458", "#57A669","#C15156","#8072B0"]
lstyles = ["solid", "dashed", "dotted", "dashdot"]

group_add_demand = list(range(17,33))
fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True, figsize=(15,4.5))
axIdx = -1
lstyleIdx = -1

prev_ovrq = None
for k in range(16):
    i = group_add_demand[k]
    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    OVRQ = setting[i]['OVRQ']
    if OVRQ != prev_ovrq:
        prev_ovrq = OVRQ
        axIdx += 1
        g = axes[axIdx]
        sns.kdeplot(data=pdf0, x="cost", bw_adjust=2, fill=False, color="blue", linewidth=2.5, ax=g)
        color = colors[axIdx]
        g.set_title("OV-RQ={0}".format(OVRQ))
        lstyleIdx = -1

    lstyleIdx += 1
    lstyle = lstyles[lstyleIdx]
    print(i, OVRQ, color, lstyle)
    sns.kdeplot(data=pdf, x="cost", bw_adjust=2, color=color, linestyle=lstyle, linewidth=2.5, fill=False, ax=g, alpha=1)
    # break

h = [Line2D([0], [0], color="blue", label='staitc'),
Line2D([0], [0], color="#DB8458", label='few-low'),
Line2D([0], [0], color="#57A669", label='few-high'),
Line2D([0], [0], color="#C15156", label='many-low'),
Line2D([0], [0], color="#8072B0", label='many-high')]

h1 = [Line2D([0], [0], label='static', color='k'),
Line2D([0], [0], label='few-low', color='k', linestyle="solid"),
Line2D([0], [0], label='few-high', color='k', linestyle="dashed"),
Line2D([0], [0], label='many-low', color='k', linestyle="dotted"),
Line2D([0], [0], label='many-high', color='k', linestyle="dashdot")]


legend1 = plt.legend(handles=h, title="OV-RQ", bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
legend2 = plt.legend(handles=h1, title="AddDemand", bbox_to_anchor=(1.05, 0.5), loc=2, borderaxespad=0., fontsize=10)
g.add_artist(legend2)
g.add_artist(legend1)
plt.suptitle('Cost Distribution',fontsize=15)
plt.savefig(folder+"/costPDFAddDemand", bbox_inches = "tight", dpi=600)
plt.show()

# %%
colors = ["#DB8458", "#57A669","#C15156","#8072B0"]
lstyles = ["solid", "dashed", "dotted", "dashdot"]

group_new_task = list(range(33,65))
fig, axes = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(15,7))

axIdx = -1
lstyleIdx = -1

for k in range(32):
    i = group_new_task[k]
    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    OVRQ = setting[i]['OVRQ']
    Pos = setting[i]["factor-setting"].split("-")[-1]
    if k%4 == 0:
        axIdx += 1
        g = axes[axIdx%2, axIdx//2]
        print(axIdx%2, axIdx//2)
        sns.kdeplot(data=pdf0, x="cost", bw_adjust=2, fill=False, color="blue", linewidth=2.5, ax=g)
        g.set_title("Pos={0} | OV-RQ={1}".format(Pos, OVRQ))
        color = colors[axIdx//2]
        lstyleIdx = -1

    lstyleIdx += 1
    lstyle = lstyles[lstyleIdx]
    print(i, OVRQ, color, lstyle)
    sns.kdeplot(data=pdf, x="cost", bw_adjust=2, color=color, linestyle=lstyle, linewidth=2.5, fill=False, ax=g, alpha=1)
    # break

h = [Line2D([0], [0], color="blue", label='staitc'),
Line2D([0], [0], color="#DB8458", label='few-low'),
Line2D([0], [0], color="#57A669", label='few-high'),
Line2D([0], [0], color="#C15156", label='many-low'),
Line2D([0], [0], color="#8072B0", label='many-high')]

h1 = [Line2D([0], [0], label='static', color='k'),
Line2D([0], [0], label='few-low', color='k', linestyle="solid"),
Line2D([0], [0], label='few-high', color='k', linestyle="dashed"),
Line2D([0], [0], label='many-low', color='k', linestyle="dotted"),
Line2D([0], [0], label='many-high', color='k', linestyle="dashdot")]


legend1 = plt.legend(handles=h, title="OV-RQ", bbox_to_anchor=(1.05, 1.85), loc=2, borderaxespad=0., fontsize=10)
legend2 = plt.legend(handles=h1, title="NewTasks", bbox_to_anchor=(1.05, 1.35), loc=2, borderaxespad=0., fontsize=10)
g.add_artist(legend2)
g.add_artist(legend1)

plt.suptitle('Cost Distribution',fontsize=15)
plt.savefig(folder+"/costPDFNewTasks", bbox_inches = "tight", dpi=600)
plt.show()


# %% plot distance to cost: scatter plot
pdf0 = pd.read_csv("{0}/instance0.csv".format(folder))
pdf0.pop("Unnamed: 0")
pdf0["dis"] = pdf0["dis"]/(setting[0]["N"]*2)

group_cgst = [1, 2, 3, 4, 5, 6, 7, 8]
color = ["darkred", "seagreen"]*4

fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True, figsize=(12,3))

prev_ovrq = None
axIdx = -1
titles = []
for k in range(8):
    i = group_cgst[k]

    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    pdf.pop("Unnamed: 0")
    pdf["dis"] = pdf["dis"]/(setting[i]["N"]*2)
    print(pdf["cost"].min(), pdf["cost"].max())
    OVRQ = setting[i]['OVRQ']
    if (OVRQ != prev_ovrq):
        axIdx += 1
        g = axes[axIdx]
        sns.kdeplot(data=pdf0, x="cost", y="dis", gridsize=20, cmap=sns.light_palette("blue", as_cmap=True), fill=False, ax=g)
        prev_ovrq = OVRQ
        titles.append(OVRQ)
    
    sns.kdeplot(data=pdf, x="cost", y="dis", gridsize=20, cmap=sns.light_palette(color[k], as_cmap=True), fill=False, ax=g, alpha=1)
    print(i)
    g.set(xlim=(0, 0.6))
    # g.set_aspect(1.0/g.get_data_ratio(), adjustable='box')

for j in range(4):
    axes[j].set_title("OV-RQ:"+titles[j])

h = [mpatches.Patch(color="blue", label='staitc'),
mpatches.Patch(color="darkred", label='20%roads'),
mpatches.Patch(color="seagreen", label='80%roads')]

plt.legend(handles=h, title="Roads Congest", bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
plt.savefig(folder+"/disCostCgst.png",bbox_inches = "tight", dpi=600)
plt.show()

# %%
group_ease= [9, 10, 11, 12, 13, 14, 15, 16]
color = ["darkred", "seagreen"]*4

fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True, figsize=(12,3))

prev_ovrq = None
axIdx = -1
titles = []
for k in range(8):
    i = group_ease[k]

    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    pdf.pop("Unnamed: 0")
    pdf["dis"] = pdf["dis"]/(setting[i]["N"]*2)
    print(pdf["cost"].min(), pdf["cost"].max())
    OVRQ = setting[i]['OVRQ']
    if (OVRQ != prev_ovrq):
        axIdx += 1
        g = axes[axIdx]
        sns.kdeplot(data=pdf0, x="cost", y="dis", gridsize=20, cmap=sns.light_palette("blue", as_cmap=True), fill=False, ax=g)
        prev_ovrq = OVRQ
        titles.append(OVRQ)
    
    sns.kdeplot(data=pdf, x="cost", y="dis", gridsize=20, cmap=sns.light_palette(color[k], as_cmap=True), fill=False, ax=g, alpha=1)
    print(i)
    g.set(xlim=(0, 0.6))
    # g.set_aspect(1.0/g.get_data_ratio(), adjustable='box')

for j in range(4):
    axes[j].set_title("OV-RQ:"+titles[j])

h = [mpatches.Patch(color="blue", label='staitc'),
mpatches.Patch(color="darkred", label='20%roads'),
mpatches.Patch(color="seagreen", label='80%roads')]

plt.legend(handles=h, title="Roads Ease", bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
plt.savefig(folder+"/disCostEase.png",bbox_inches = "tight", dpi=600)
plt.show()

# %%
group_add_demand= list(range(17, 33))
color = ["darkred", "seagreen", "black", "purple"]*4

fig, axes = plt.subplots(nrows=1, ncols=4, sharex=True, sharey=True, figsize=(12,3))

prev_ovrq = None
axIdx = -1
titles = []
for k in range(16):
    i = group_add_demand[k]

    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    pdf.pop("Unnamed: 0")
    pdf["dis"] = pdf["dis"]/(setting[i]["N"]*2)
    print(pdf["cost"].min(), pdf["cost"].max())
    OVRQ = setting[i]['OVRQ']
    if (OVRQ != prev_ovrq):
        axIdx += 1
        g = axes[axIdx]
        sns.kdeplot(data=pdf0, x="cost", y="dis", gridsize=20, cmap=sns.light_palette("blue", as_cmap=True), fill=False, ax=g)
        prev_ovrq = OVRQ
        titles.append(OVRQ)
    
    sns.kdeplot(data=pdf, x="cost", y="dis", gridsize=20, cmap=sns.light_palette(color[k], as_cmap=True), fill=False, ax=g, alpha=1)
    print(i)
    g.set(xlim=(0, 0.6))#ylim=(20,150)
    # g.set_aspect(1.0/g.get_data_ratio(), adjustable='box')

for j in range(4):
    axes[j].set_title("OV-RQ:"+titles[j])

h = [mpatches.Patch(color="blue", label='staitc'),
mpatches.Patch(color="darkred", label='few-low'),
mpatches.Patch(color="seagreen", label='few-high'),
mpatches.Patch(color="black", label='many-low'),
mpatches.Patch(color="purple", label='many-high')]

plt.legend(handles=h, title="AddDemand", bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
plt.savefig(folder+"/disCostAddDemand.png",bbox_inches = "tight", dpi=600)
plt.show()


# %%
group_new_task= list(range(33, 65))
color = ["darkred", "seagreen", "black", "purple"]*8

fig, axes = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(15,6))

axIdx = -1
titles = []
for k in range(32):
    i = group_new_task[k]

    pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
    pdf.pop("Unnamed: 0")
    pdf["dis"] = pdf["dis"]/(setting[i]["N"]*2)
    print(pdf["cost"].min(), pdf["cost"].max())
    OVRQ = setting[i]['OVRQ']
    dynamic = setting[i]["factor-setting"]
    
    if k % 4 == 0:
        axIdx += 1
        g = axes[axIdx%2, axIdx//2]
        print(axIdx%2, axIdx//2)
        sns.kdeplot(data=pdf0, x="cost", y="dis", gridsize=20, cmap=sns.light_palette("blue", as_cmap=True), fill=False, ax=g)
    
    sns.kdeplot(data=pdf, x="cost", y="dis", gridsize=20, cmap=sns.light_palette(color[k], as_cmap=True), fill=False, ax=g, alpha=1)
    print(i, OVRQ, dynamic)
    g.set(xlim=(0, 0.6)) #, ylim=(20,200)
    # g.set_aspect(1.0/g.get_data_ratio(), adjustable='box')

titles = ['few-low', 'few-high', 'many-low', 'many-high']
for j in range(4):
    axes[0, j].set_title("OV-RQ:"+titles[j])

h = [mpatches.Patch(color="blue", label='staitc'),
mpatches.Patch(color="darkred", label='few-low'),
mpatches.Patch(color="seagreen", label='few-high'),
mpatches.Patch(color="black", label='many-low'),
mpatches.Patch(color="purple", label='many-high')]

plt.legend(handles=h, title="NewTasks", bbox_to_anchor=(1.05, 1.45), loc=2, borderaxespad=0., fontsize=10)
plt.savefig(folder+"/disCostNewTasks.png", bbox_inches = "tight", dpi=600)
plt.show()


# %% FDC

filepath = folder+"/instance{0}.txt".format(0)
df = pd.DataFrame({"FDC": [], "OV-RQ": [], "factor": [], "Dynamic":[],  "TaskNum":[]})
dfTime = pd.DataFrame({"UsedTime": [], "OV-RQ": [],"factor": [], "Dynamic": [], "TaskNum":[]})
dfLoNum = pd.DataFrame({"LoNum": [], "OV-RQ": [],"factor": [], "Dynamic": [], "TaskNum":[]})
for i in range(65):
    filepath = folder+"/instance{0}.txt".format(i)
    OVRQ = setting[i]['OVRQ']
    factor = setting[i]['factor']
    dynamic = setting[i]['factor-setting']

    with open(filepath, "r") as f:
        hist = f.readline().strip().split(",")
        loNum = np.sum(list(map(lambda x:int(x), hist[1:-1])))

        hit = f.readline().strip().split(",")
        cost = f.readline().strip().split(",")
        ttime = f.readline().strip().split(",")
        avrtime = f.readline().strip().split(",")
        fdc = f.readline().strip().split(",")
        if i == 0:
            static_fdc= float(fdc[1])
            static_usedtime = float(avrtime[1])
            static_loNum = loNum
            continue


        xml_path = "instance/egl-e2-A-0/instance{0}.xml".format(i)
        tree = ET.ElementTree(file=xml_path)
        root = tree.getroot()
        task_num = int(root[1].attrib["num"])//2
        print("task num: ", task_num)

        d = {"FDC": float(fdc[1]), "OV-RQ": [OVRQ], "factor": [factor], "Dynamic": [dynamic], "TaskNum":[task_num]}
        dfTemp = pd.DataFrame(data=d)
        df = df.append(dfTemp)

        d = {"UsedTime": float(avrtime[1])-static_usedtime, "OV-RQ": [OVRQ],"factor": [factor],  "Dynamic": [dynamic], "TaskNum":[task_num]}
        dfTemp = pd.DataFrame(data=d)
        dfTime = dfTime.append(dfTemp)

        d = {"LoNum": loNum - static_loNum, "OV-RQ": [OVRQ], "factor": [factor], "Dynamic": [dynamic],"TaskNum":[task_num]}
        dfTemp = pd.DataFrame(data=d)
        dfLoNum = dfLoNum.append(dfTemp)
        # break

#%%
metric = "UsedTime"
inputdf = dfTime


fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True, sharex=True, figsize=(12,12),tight_layout=True)
g = axes[0, 0]
sns.regplot(x="TaskNum", y=metric, data=inputdf, scatter=False, ax = g)
sns.scatterplot(data=inputdf[inputdf["factor"] == "Congest"], x="TaskNum", y=metric,hue="Dynamic",ax=g, s=150)
g.legend(title="Congest",bbox_to_anchor=(0, 1), loc=2, borderaxespad=0., fontsize=10)
g.set_title("Dynamic Event: Congest")

g = axes[0, 1]
sns.regplot(x="TaskNum", y=metric, data=inputdf, scatter=False, ax = g)
sns.scatterplot(data=inputdf[inputdf["factor"] == "Ease"], x="TaskNum", y=metric,hue="Dynamic",ax=g, s=150)
g.legend(title="Ease",bbox_to_anchor=(0, 1), loc=2, borderaxespad=0., fontsize=10)
g.set_title("Dynamic Event: Ease")

g = axes[1, 0]
sns.regplot(x="TaskNum", y=metric, data=inputdf, scatter=False, ax = g)
sns.scatterplot(data=inputdf[inputdf["factor"] == "AddDemand"], x="TaskNum", y=metric,hue="Dynamic",ax=g, s=150)
g.legend(title="AddDemand",bbox_to_anchor=(0, 1), loc=2, borderaxespad=0., fontsize=10)
g.set_title("Dynamic Event: AddDemand")

g = axes[1, 1]
sns.regplot(x="TaskNum", y=metric, data=inputdf, scatter=False, ax = g)
sns.scatterplot(data=inputdf[inputdf["factor"] == "NewTask"], x="TaskNum", y=metric,hue="Dynamic",ax=g, s=150)
g.legend(title="NewTask",bbox_to_anchor=(0, 1), loc=2, borderaxespad=0., fontsize=10, ncol=2)
g.set_title("Dynamic Event: NewTask")
plt.savefig(folder+'/UsedTime_tasknum.png',bbox_inches = "tight", dpi=600)


ax1 = sns.regplot(x="TaskNum", y="LoNum", data=dfLoNum, scatter=False)
sns.scatterplot(data=dfLoNum, x="TaskNum",y="LoNum",ax=ax1,hue="OV-RQ", s=150)
plt.savefig(folder+'/LONum_tasknum.png',bbox_inches = "tight", dpi=600)

# data = df.pivot("Dynamic", "OV-RQ", "FDC")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True), linewidths=.5)
# g.set(title='Heatmap for Relative Fitness Distance Correlation')
# plt.savefig(folder+'/FDC.png',bbox_inches = "tight", dpi=600)
# plt.show()

# g.clear()
# data = dfTime.pivot("Dynamic", "OV-RQ", "UsedTime")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True),linewidths=.5)
# g.set(title='Heatmap for Relative UsedTime')
# plt.savefig(folder+'/UsedTime.png',bbox_inches = "tight", dpi=600)
# plt.show()

# g.clear()
# data = dfLoNum.pivot("Dynamic", "OV-RQ", "LoNum")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True),linewidths=.5)
# g.set(title='Heatmap for Relative Number of Local Optimal')
# plt.savefig(folder+'/LoNum.png',bbox_inches = "tight", dpi=600)
# plt.show()
# %%





# %% Temp codes
# static = pd.DataFrame({"mean":[], "skewness":[], "kurtosis": [], "OV-RQ": [], "Dynamic": []})
# group_cgst = [0, 1, 2, 3, 4, 5, 6, 7, 8]

# for i in group_cgst:
#     pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
#     OVRQ = setting[i]['OVRQ']
#     if i == 0:
#         dynamic = "static"
#         cost0 = pdf["cost"].mean()
#     else:
#         dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
#     dfTmp = pd.DataFrame({"mean":[pdf["cost"].mean()], "skewness":[pdf["cost"].skew()], "kurtosis": [pdf["cost"].kurtosis()], "OV-RQ": [OVRQ], "Dynamic": [dynamic]})
#     static = static.append(dfTmp, ignore_index=True)
#     # break

# g = sns.scatterplot(data=static, x="skewness", y="kurtosis", hue="OV-RQ", palette="deep", style="Dynamic", s = 200)
# g.set(ylim=(-1.0, 3.5), xlim=(-0.3, 1.5))
# plt.legend(bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
# filename = "/statisCgst.png"
# plt.savefig(folder+filename, bbox_inches = "tight", dpi=600)
# plt.show()
# static = pd.DataFrame({"mean":[], "skewness":[], "kurtosis": [], "OV-RQ": [], "Dynamic": []})
# group_ease= [0, 9, 10, 11, 12, 13, 14, 15, 16]
# for i in group_ease:
#     pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
#     OVRQ = setting[i]['OVRQ']
#     if i == 0:
#         dynamic = "static"
#         cost0 = pdf["cost"].mean()
#     else:
#         dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
#     dfTmp = pd.DataFrame({"mean":[pdf["cost"].mean()], "skewness":[pdf["cost"].skew()], "kurtosis": [pdf["cost"].kurtosis()], "OV-RQ": [OVRQ], "Dynamic": [dynamic]})
#     static = static.append(dfTmp, ignore_index=True)
#     # break

# g = sns.scatterplot(data=static, x="skewness", y="kurtosis", hue="OV-RQ", palette="deep", style="Dynamic", s = 200)
# g.set(ylim=(-1.0, 3.5), xlim=(-0.3, 1.5))
# plt.legend(bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)
# filename = "/statisEase.png"
# # filename = "/statisAddDemand.png"
# # filename = "/statisNewTask.png"
# plt.savefig(folder+filename, bbox_inches = "tight", dpi=600)
# plt.show()

# plt.clf()
# static = pd.DataFrame({"mean":[], "skewness":[], "kurtosis": [], "OV-RQ": [], "Dynamic": []})
# group_add_demand= [0] + list(range(17, 33))
# for i in group_add_demand:
#     pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
#     OVRQ = setting[i]['OVRQ']
#     if i == 0:
#         dynamic = "static"
#         cost0 = pdf["cost"].mean()
#     else:
#         dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
#     dfTmp = pd.DataFrame({"mean":[pdf["cost"].mean()], "skewness":[pdf["cost"].skew()], "kurtosis": [pdf["cost"].kurtosis()], "OV-RQ": [OVRQ], "Dynamic": [dynamic]})
#     static = static.append(dfTmp, ignore_index=True)
#     # break

# g = sns.scatterplot(data=static, x="skewness", y="kurtosis", hue="OV-RQ", palette="deep", style="Dynamic", s = 200)
# g.set(ylim=(-1.0, 3.5), xlim=(-0.3, 1.5))
# plt.legend(bbox_to_anchor=(1.05, 0.95), loc=2, borderaxespad=0., fontsize=10)

# filename = "/statisAddDemand.png"
# plt.savefig(folder+filename, bbox_inches = "tight", dpi=600)
# plt.show()


# static = pd.DataFrame({"mean":[], "skewness":[], "kurtosis": [], "OV-RQ": [], "Dynamic": []})
# group_new_task= [0] + list(range(33, 65))
# for i in group_new_task:
#     pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
#     OVRQ = setting[i]['OVRQ']
#     if i == 0:
#         dynamic = "static"
#         cost0 = pdf["cost"].mean()
#     else:
#         dynamic = "{0}:{1}".format(setting[i]['factor'], setting[i]['factor-setting'])
#     dfTmp = pd.DataFrame({"mean":[pdf["cost"].mean()], "skewness":[pdf["cost"].skew()], "kurtosis": [pdf["cost"].kurtosis()], "OV-RQ": [OVRQ], "Dynamic": [dynamic]})
#     static = static.append(dfTmp, ignore_index=True)
#     # break

# g = sns.scatterplot(data=static, x="skewness", y="kurtosis", hue="OV-RQ", palette="deep", style="Dynamic", s = 200)
# g.set(ylim=(-1.0, 3.5), xlim=(-0.3, 1.5))
# plt.legend(bbox_to_anchor=(1.05, 1.05), loc=2, borderaxespad=0., fontsize=10)
# filename = "/statisNewTask.png"
# plt.savefig(folder+filename, bbox_inches = "tight", dpi=600)
# plt.show()