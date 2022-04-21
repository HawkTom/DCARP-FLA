# %%
import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import mean
import pandas as pd
import seaborn as sns
import scipy.stats as scis
import scikit_posthocs as sp
from lxml import etree as ET
import Orange


#%%
def graph_ranks(avranks, names, cd=None, cdmethod=None, lowv=None, highv=None,
                width=6, textspace=1, reverse=False, filename=None, **kwargs):
    try:
        import matplotlib.pyplot as plt
        import math
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        plt.rcParams.update({
        "font.size":13,
        "text.usetex": True,
        "font.family": "monospace"
        })
    except ImportError:
        raise ImportError("Function graph_ranks requires matplotlib.")

    width = float(width)
    textspace = float(textspace)

    def nth(l, n):
        n = lloc(l, n)
        return [a[n] for a in l]

    def lloc(l, n):
        if n < 0:
            return len(l[0]) + n
        else:
            return n

    def mxrange(lr):
        if not len(lr):
            yield ()
        else:
            # it can work with single numbers
            index = lr[0]
            if isinstance(index, int):
                index = [index]
            for a in range(*index):
                for b in mxrange(lr[1:]):
                    yield tuple([a] + list(b))

    def print_figure(fig, *args, **kwargs):
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(bbox_inches="tight", dpi=600, *args, **kwargs)

    sums = avranks

    tempsort = sorted([(a, i) for i, a in enumerate(sums)], reverse=reverse)
    ssums = nth(tempsort, 0)
    sortidx = nth(tempsort, 1)
    nnames = [names[x] for x in sortidx]

    if lowv is None:
        lowv = min(1, int(math.floor(min(ssums))))
    if highv is None:
        highv = max(len(avranks), int(math.ceil(max(ssums))))

    cline = 0.4

    k = len(sums)

    lines = None

    linesblank = 0
    scalewidth = width - 2 * textspace

    def rankpos(rank):
        if not reverse:
            a = rank - lowv
        else:
            a = highv - rank
        return textspace + scalewidth / (highv - lowv) * a

    distanceh = 0.25

    if cd and cdmethod is None:
        # get pairs of non significant methods

        def get_lines(sums, hsd):
            # get all pairs
            lsums = len(sums)
            allpairs = [(i, j) for i, j in mxrange([[lsums], [lsums]]) if j > i]
            # remove not significant
            notSig = [(i, j) for i, j in allpairs
                      if abs(sums[i] - sums[j]) <= hsd]
            # keep only longest

            def no_longer(ij_tuple, notSig):
                i, j = ij_tuple
                for i1, j1 in notSig:
                    if (i1 <= i and j1 > j) or (i1 < i and j1 >= j):
                        return False
                return True

            longest = [(i, j) for i, j in notSig if no_longer((i, j), notSig)]

            return longest

        lines = get_lines(ssums, cd)
        linesblank = 0.2 + 0.2 + (len(lines) - 1) * 0.1

        # add scale
        distanceh = 0.25
        cline += distanceh

    # calculate height needed height of an image
    minnotsignificant = max(2 * 0.2, linesblank)
    height = cline + ((k + 1) / 2) * 0.2 + minnotsignificant

    fig = plt.figure(figsize=(width, height))
    fig.set_facecolor('white')
    ax = fig.add_axes([0, 0, 1, 1])  # reverse y axis
    ax.set_axis_off()

    hf = 1. / height  # height factor
    wf = 1. / width

    def hfl(l):
        return [a * hf for a in l]

    def wfl(l):
        return [a * wf for a in l]


    # Upper left corner is (0,0).
    ax.plot([0, 1], [0, 1], c="w")
    ax.set_xlim(0, 1)
    ax.set_ylim(1, 0)

    def line(l, color='k', **kwargs):
        """
        Input is a list of pairs of points.
        """
        ax.plot(wfl(nth(l, 0)), hfl(nth(l, 1)), color=color, **kwargs)

    def text(x, y, s, *args, **kwargs):
        ax.text(wf * x, hf * y, s, *args, **kwargs)

    line([(textspace, cline), (width - textspace, cline)], linewidth=0.7)

    bigtick = 0.1
    smalltick = 0.05

    tick = None
    for a in list(np.arange(lowv, highv, 0.5)) + [highv]:
        tick = smalltick
        if a == int(a):
            tick = bigtick
        line([(rankpos(a), cline - tick / 2),
              (rankpos(a), cline)],
             linewidth=0.7)

    for a in range(lowv, highv + 1):
        text(rankpos(a), cline - tick / 2 - 0.05, str(a),
             ha="center", va="bottom")

    k = len(ssums)

    for i in range(math.ceil(k / 2)):
        chei = cline + minnotsignificant + i * 0.2
        line([(rankpos(ssums[i]), cline),
              (rankpos(ssums[i]), chei),
              (textspace - 0.1, chei)],
             linewidth=0.7)
        text(textspace - 0.2, chei, nnames[i], ha="right", va="center")

    for i in range(math.ceil(k / 2), k):
        chei = cline + minnotsignificant + (k - i - 1) * 0.2
        line([(rankpos(ssums[i]), cline),
              (rankpos(ssums[i]), chei),
              (textspace + scalewidth + 0.1, chei)],
             linewidth=0.7)
        text(textspace + scalewidth + 0.2, chei, nnames[i],
             ha="left", va="center")

    if cd and cdmethod is None:
        # upper scale
        if not reverse:
            begin, end = rankpos(lowv), rankpos(lowv + cd)
        else:
            begin, end = rankpos(highv), rankpos(highv - cd)

        line([(begin, distanceh), (end, distanceh)], linewidth=0.7)
        line([(begin, distanceh + bigtick / 2),
              (begin, distanceh - bigtick / 2)],
             linewidth=0.7)
        line([(end, distanceh + bigtick / 2),
              (end, distanceh - bigtick / 2)],
             linewidth=0.7)
        text((begin + end) / 2, distanceh - 0.05, "CD",
             ha="center", va="bottom")

        # no-significance lines
        def draw_lines(lines, side=0.05, height=0.1):
            start = cline + 0.2
            for l, r in lines:
                line([(rankpos(ssums[l]) - side, start),
                      (rankpos(ssums[r]) + side, start)],
                     linewidth=2.5)
                start += height

        draw_lines(lines)

    elif cd:
        begin = rankpos(avranks[cdmethod] - cd)
        end = rankpos(avranks[cdmethod] + cd)
        line([(begin, cline), (end, cline)],
             linewidth=2.5)
        line([(begin, cline + bigtick / 2),
              (begin, cline - bigtick / 2)],
             linewidth=2.5)
        line([(end, cline + bigtick / 2),
              (end, cline - bigtick / 2)],
             linewidth=2.5)

    if filename:
        print_figure(fig, filename, **kwargs)

#%%
setting = {}
setting[0] = {"OVRQ": "static", "factor": "static", "factor-setting": "static","N":72}
with open("settings","r") as f:
    f.readline()
    for i in range(1, 65):
        s = f.readline().split(" ")
        ss = s[2].strip().split(":")
        setting[int(s[0])] = {"OVRQ": s[1], "factor": ss[0], "factor-setting": ss[1]}

        # xml_path = "instance/egl-e2-A-{0}/instance{1}.xml".format(mapidx, i)
        # tree = ET.ElementTree(file=xml_path)
        # root = tree.getroot()
        # task_num = int(root[1].attrib["num"])//2
        # setting[i]["N"] = task_num
        # print("task num: ", task_num)

settinglabels = {}
settinglabels[0] = "static"
with open("settingslabels","r") as f:
    f.readline()
    for i in range(1, 65):
        s = f.readline().strip().split(" ")
        settinglabels[int(s[0])] = s[1]

#%%
def plot_avranks(data, maxruns, isRank=True):

    def draw(df):
        lables = list(df.columns)
        g = sns.PairGrid(df, y_vars=lables[0],x_vars=lables[1:-1], height=3, aspect=1.2)
        g.map(sns.pointplot, color="xkcd:plum")
    
    def cal_statisitcs(x, group, event, metric):
        formula = []
        # print(x)
        if isRank:
            avranks = np.mean(scis.rankdata(x, axis=0), axis=1)
        else:
            avranks = np.mean(x, axis=1)
        print(avranks)
        for i in range(len(group)):
            formula.append("x[{0}, :]".format(i))
        output = eval("scis.friedmanchisquare(" + ", ".join(formula)+")")
        if output.pvalue < 0.05:
            names = []
            for k in group:
                names.append(settinglabels[k])

            if len(group) != 32:
                cd = Orange.evaluation.compute_CD(avranks, maxruns)
                # graph_ranks(avranks, names=names, cd=cd)
                # graph_ranks(avranks, names=names, cd=cd, filename="analysis/lor/{0}_{1}.png".format(metric, event))
            else:
                k = len(avranks)
                cd = 3.780 * (k * (k + 1) / (6.0 * maxruns)) ** 0.5
                # graph_ranks(avranks, names=names, cd=cd, width=10)
                # graph_ranks(avranks, names=names, cd=cd, width=10, filename="analysis/lor/{0}_{1}.png".format(metric, event))
            
            
            print("p:", output.pvalue, "cd:", cd)
        else:
            print("p:", output.pvalue)
        return avranks

    group_cgst = [1, 2, 3, 4, 5, 6, 7, 8]
    group_ease = list(range(9,17))
    group_add_demand = list(range(17,33))
    group_new_task = list(range(33,65))

    avranks1 = []
    tmp = cal_statisitcs(data[group_cgst, :], group_cgst, "cgst", "mean_cost")
    avranks1 += tmp.tolist()
    tmp = cal_statisitcs(data[group_ease, :], group_ease, "ease", "mean_cost")
    avranks1 += tmp.tolist()
    tmp = cal_statisitcs(data[group_add_demand, :], group_add_demand, "adddemand", "mean_cost")
    avranks1 += tmp.tolist()
    tmp = cal_statisitcs(data[group_new_task, :], group_new_task,"newtasks", "mean_cost")
    avranks1 += tmp.tolist()
    # print(avranks1, "\n")


    dfranksCgst = pd.DataFrame({"avrranks": [0.0]*8, "OV-RQ": [""]*8, "Dynamic": [""]*8, "TaskNum":[""]*8})
    for k in range(len(group_cgst)):
        i = group_cgst[k]
        OVRQ = setting[i]['OVRQ']
        factor = setting[i]['factor']
        dynamic = setting[i]['factor-setting']
        dfranksCgst.at[k, "avrranks"] = avranks1[i-1]
        dfranksCgst.at[k, "OV-RQ"] = OVRQ
        dfranksCgst.at[k, "Dynamic"] = dynamic

    draw(dfranksCgst)

    dfranksEase = pd.DataFrame({"avrranks": [0.0]*8, "OV-RQ": [""]*8, "Dynamic": [""]*8, "TaskNum":[""]*8})
    for k in range(len(group_ease)):
        i = group_ease[k]
        OVRQ = setting[i]['OVRQ']
        factor = setting[i]['factor']
        dynamic = setting[i]['factor-setting']
        dfranksEase.at[k, "avrranks"] = avranks1[i-1]
        dfranksEase.at[k, "OV-RQ"] = OVRQ
        dfranksEase.at[k, "Dynamic"] = dynamic
    draw(dfranksEase)

    dfranksAddDemand = pd.DataFrame({"avrranks": [0.0]*16, "OV-RQ": [""]*16, "Num": [""]*16,  "Value": [""]*16, "TaskNum":[""]*16})
    for k in range(len(group_add_demand)):
        i = group_add_demand[k]
        OVRQ = setting[i]['OVRQ']
        factor = setting[i]['factor']
        dynamic = setting[i]['factor-setting'].split("-")
        dfranksAddDemand.at[k, "avrranks"] = avranks1[i-1]
        dfranksAddDemand.at[k, "OV-RQ"] = OVRQ
        dfranksAddDemand.at[k, "Num"] = dynamic[0]
        dfranksAddDemand.at[k, "Value"] = dynamic[1]
    draw(dfranksAddDemand)


    dfranksNewTask = pd.DataFrame({"avrranks": [0.0]*32, "OV-RQ": [""]*32,  "Num": [""]*32,  "Value": [""]*32, "Pos":[""]*32, "TaskNum":[""]*32})
    for k in range(len(group_new_task)):
        i = group_new_task[k]
        OVRQ = setting[i]['OVRQ']
        factor = setting[i]['factor']
        dynamic = setting[i]['factor-setting'].split("-")
        dfranksNewTask.at[k, "avrranks"] = avranks1[i-1]
        dfranksNewTask.at[k, "OV-RQ"] = OVRQ
        dfranksNewTask.at[k, "Num"] = dynamic[0]
        dfranksNewTask.at[k, "Value"] = dynamic[1]
        dfranksNewTask.at[k, "Pos"] = dynamic[2]
    draw(dfranksNewTask)
print("function Ok.")    

#%%
maxruns = 10

pdf0 = pd.read_csv("analysis/lor/instance0.csv")
pdf0.pop("Unnamed: 0")
pdf0["dis"] = pdf0["dis"]/(setting[0]["N"]*2)
skewness0, kurtosis0, mean0 = pdf0["cost"].skew(), pdf0["cost"].kurt(), pdf0["cost"].mean()

skewness = np.zeros((65, maxruns))
kurtosis = np.zeros((65, maxruns))
means = np.zeros((65, maxruns))

for mapidx in range(maxruns):
    folder = "analysis/lor/egl-e2-A-{0}".format(mapidx)
    for i in range(1, 65):
        pdf = pd.read_csv("{0}/instance{1}.csv".format(folder, i))
        skewness_, kurtosis_, mean_ = pdf["cost"].skew(), pdf["cost"].kurt(), pdf["cost"].mean()
        skewness[i, mapidx] = skewness_
        kurtosis[i, mapidx] = kurtosis_
        means[i, mapidx] = mean_
    print(folder)

plot_avranks(means, maxruns, isRank=True)

# %% Heat map for three statistical measurements
df = pd.DataFrame({"mean": [], "skewness":[], "kurtosis":[], "Outside Vehicles-Remaining Capacities": [], "Dynamic Events": []})
with open("acronyms","r") as f:
    f.readline()
    means_avr = np.mean(means, axis=1)
    skewness_avr = np.mean(skewness, axis=1)
    kurtosis_avr = np.mean(kurtosis, axis=1)
    for i in range(1, 65):
        s = f.readline().strip().split(" ")
        
        # print(s)

        # xml_path = "instance/egl-e2-A-0/instance{0}.xml".format(i)
        # tree = ET.ElementTree(file=xml_path)
        # root = tree.getroot()
        # task_num = int(root[1].attrib["num"])//2
        # print("task num: ", task_num)
        a, b, c = means_avr[i], skewness_avr[i], kurtosis_avr[i]


        d = {"mean": [a], "skewness":[b], "kurtosis":[c], "Outside Vehicles-Remaining Capacities": [s[1]], "Dynamic Events": [s[2]]}
        dfTemp = pd.DataFrame(data=d)
        df = df.append(dfTemp)

data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "mean")
ax = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True), linewidths=.01, fmt=".2f", annot=True, vmin=0, vmax=0.6)
ax.set(title='Heatmap of Mean of Normalised Cost')
plt.savefig(folder+'/cost_mean.eps',bbox_inches = "tight", dpi=300,format="eps")

data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "skewness")
ax = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True), linewidths=.01, fmt=".2f", annot=True, vmin=-0.1, vmax=1.5)
ax.set(title='Heatmap of Skewness of Normalised Cost')
plt.savefig(folder+'/cost_skewness.eps',bbox_inches = "tight", dpi=300,format="eps")


# data = df.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "kurtosis")
# ax = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True), linewidths=.01, fmt=".2f", annot=True, vmin=-0.6, vmax=3.2)
# ax.set(title='Heatmap of Kurtosis of Normalised Cost')


# %%
maxruns = 10

folder = "analysis/lor"
filepath = folder+"/instance0.txt"
with open(filepath, "r") as f:
    hist = f.readline().strip().split(",")
    loNum = np.sum(list(map(lambda x:int(x), hist[1:-1])))

    hit = f.readline().strip().split(",")
    cost = f.readline().strip().split(",")
    ttime = f.readline().strip().split(",")
    avrtime = f.readline().strip().split(",")
    fdc0 = float(f.readline().strip().split(",")[1])

dfTime = pd.DataFrame({"UsedTime": [], "OV-RQ": [],"factor": [], "Dynamic": [], "TaskNum":[]})
dfLoNum = pd.DataFrame({"LoNum": [], "OV-RQ": [],"factor": [], "Dynamic": [], "TaskNum":[]})

fdcarray = np.zeros((65, maxruns))
usedtimearray = np.zeros((65, maxruns))
lonum = np.zeros((65, maxruns))

for mapidx in range(maxruns):
    for i in range(1, 65):
        filepath = folder+"/egl-e2-A-{0}/instance{1}.txt".format(mapidx, i)
        # OVRQ = setting[i]['OVRQ']
        # factor = setting[i]['factor']
        # dynamic = setting[i]['factor-setting']

        with open(filepath, "r") as f:
            hist = f.readline().strip().split(",")
            loNum = np.sum(list(map(lambda x:int(x), hist[1:-1])))

            hit = f.readline().strip().split(",")
            cost = f.readline().strip().split(",")
            ttime = f.readline().strip().split(",")
            avrtime = f.readline().strip().split(",")
            fdc = float(f.readline().strip().split(",")[1])
            

            xml_path = "instance/egl-e2-A-{0}/instance{1}.xml".format(mapidx, i)
            tree = ET.ElementTree(file=xml_path)
            root = tree.getroot()
            task_num = int(root[1].attrib["num"])//2
            # print(task_num)
            fdcarray[i, mapidx] = fdc

            usedtimearray[i, mapidx] = float(avrtime[1])

            lonum[i, mapidx] = loNum


dfTime = pd.DataFrame({"usedTime": [0.0]*64, "Outside Vehicles-Remaining Capacities": [""]*64, "factor": [""]*64, "Dynamic Events":[""]*64,  "TaskNum":[0]*64})
dfTime["usedTime"] = np.mean(usedtimearray[1:, :], axis=1)


dfLoNum = pd.DataFrame({"LoNum": [0.0]*64, "Outside Vehicles-Remaining Capacities": [""]*64, "factor": [""]*64, "Dynamic Events":[""]*64,  "TaskNum":[0]*64})
dfLoNum["LoNum"] = np.mean(lonum[1:, :], axis=1)

dffdc = pd.DataFrame({"FDC": [0.0]*64, "Outside Vehicles-Remaining Capacities": [""]*64, "factor": [""]*64, "Dynamic Events":[""]*64,  "TaskNum":[0]*64})
dffdc["FDC"] = np.mean(fdcarray[1:, :], axis=1)

with open("acronyms","r") as f:
    f.readline()
    for i in range(64):
        s = f.readline().strip().split(" ")
        OVRQ, dynamic = s[1], s[2]
        dfTime.at[i, "Outside Vehicles-Remaining Capacities"] = OVRQ
        dfTime.at[i, "Dynamic Events"] = dynamic
        dfLoNum.at[i, "Outside Vehicles-Remaining Capacities"] = OVRQ
        dfLoNum.at[i, "Dynamic Events"] = dynamic
        dffdc.at[i, "Outside Vehicles-Remaining Capacities"] = OVRQ
        dffdc.at[i, "Dynamic Events"] = dynamic


# for i in range(64):
#     OVRQ = setting[i+1]['OVRQ']
#     factor = setting[i+1]['factor']
#     dynamic = setting[i+1]['factor-setting']
#     dfTime.at[i, "OV-RQ"] = OVRQ
#     dfTime.at[i, "Dynamic"] = "{0}:{1}".format(factor, dynamic)
#     dfLoNum.at[i, "OV-RQ"] = OVRQ
#     dfLoNum.at[i, "Dynamic"] = "{0}:{1}".format(factor, dynamic)
#     dffdc.at[i, "OV-RQ"] = OVRQ
#     dffdc.at[i, "Dynamic"] = "{0}:{1}".format(factor, dynamic)

# data = dfTime.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "usedTime")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True,  reverse=True), linewidths=.1, fmt=".1f", annot=True, vmin=0, vmax=250)
# g.set(title='Heatmap for Time for One Time Local Search')
# plt.savefig(folder+'/usedTime.eps',bbox_inches = "tight", dpi=300,format="eps")
# # plt.show()

# data = dfLoNum.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "LoNum")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True, reverse=True),linewidths=.1, fmt=".0f", annot=True)
# g.set(title='Heatmap for Number of Local Optimal')
# plt.savefig(folder+'/LONum.eps',bbox_inches = "tight", dpi=300,format="eps")
# plt.show()

# data = dffdc.pivot("Dynamic Events", "Outside Vehicles-Remaining Capacities", "FDC")
# g = sns.heatmap(data, cmap=sns.cubehelix_palette(as_cmap=True), linewidths=.1, fmt=".2f", annot=True)
# g.set(title='Heatmap for Fitness Distance Correlation')
# plt.savefig(folder+'/FDC.eps',bbox_inches = "tight", dpi=300,format="eps")
# plt.show()

# plot_avranks(fdcarray, maxruns)

# %%
maxruns = 10
x = np.arange(10)*0.1+0.05

dis_ratio = np.zeros((65, maxruns))
pearcoef = np.zeros((65, maxruns))
for mapidx in range(maxruns):
    for i in range(1, 65):
        filepath = folder+"/egl-e2-A-{0}/instance{1}.txt".format(mapidx, i)
        with open(filepath, "r") as f:
            for k in range(6):
                f.readline()
            dis1 = f.readline().strip().split(",")
            dis2 = f.readline().strip().split(",")
        dis_ratio[i, mapidx] = int(dis1[2])/int(dis1[4])
        dis2.pop(0)
        dis2.pop(-1)
        dis2 = list(map(lambda x: int(x), dis2))
        pearcoef[i, mapidx] = scis.pearsonr(x, dis2)[0]

plot_avranks(dis_ratio, maxruns, isRank=True)
plot_avranks(pearcoef, maxruns)

# %% return probability
maxruns = 5
for mapidx in range(maxruns):
    for i in range(1, 65):
        filepath = "result/returnp/egl-e2-A-{0}/instance{1}.txt".format(mapidx, i)
        with open(filepath, "r") as f:
            data = f.readline().split(",")
            data.pop(-1)
            data.pop(-1)
            data = list(map(lambda x:int(x), data))
            if len(data) == 0:
                continue
            if data[0] < data[1]:
                print(filepath)
            # np.polyfit(np.arange(1, len(data)+1), np.log(data), 1)
            # sns.regplot(x=np.arange(len(data)), y=np.log(data))
        # print(data)
    break

# %%
2,3,4,7,12,20,23,24,33,35-49,52,53,54,56-64