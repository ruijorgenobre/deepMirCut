import sys, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
import statistics
import matplotlib.font_manager as fm
font_dirs = ['/raid1/home/bb/bellji/fonts', ]
font_files = fm.findSystemFonts(fontpaths=font_dirs)
font_list = fm.createFontList(font_files)
fm.fontManager.ttflist.extend(font_list)

DVMAX_GRAPH_INCLUDE_DOTS = True
DVCUTS_GRAPH_INCLUDE_DOTS = True
DELTA_GRAPH_INCLUDE_DOTS = True

def readCutSite(cutSite):
    if cutSite == '-':
        cutSite = ["-","-"]
    else:
        cutSite = [int(x)-1 for x in cutSite.split(',')]
    return cutSite

def read_cuts_file(dataSetFile):
    cuts = {}
    f = open(dataSetFile,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = row_data[0:12]
        drosha5p = readCutSite(drosha5p)
        dicer5p = readCutSite(dicer5p)
        dicer3p = readCutSite(dicer3p)
        drosha3p = readCutSite(drosha3p)
        cuts[id] = [drosha5p[0],dicer5p[0],dicer3p[0],drosha3p[0]]
    f.close()
    return cuts

def get_cuts_from_line(line):
    drosha5p = -1
    dicer5p = -1
    dicer3p = -1
    drosha3p = -1
    linechars = line.split(",")
    for i in range(0,len(linechars)):
        cut_positions = {}
        if linechars[i] == "DR5":
            drosha5p = i
        if linechars[i] == "DC5":
            dicer5p = i
        if linechars[i] == "DC3":
            dicer3p = i
        if linechars[i] == "DR3":
            drosha3p = i
    return drosha5p,dicer5p,dicer3p,drosha3p

def read_predictions_file(predictionsFile):
    predicted_cuts = {}
    f = open(predictionsFile,"r");
    for line in f.readlines():
        id,name,seq,fold,cuts,predictions = line.rstrip().split()
        predicted_cuts[id] = {}
        predicted_cuts[id]["DR5"],predicted_cuts[id]["DC5"],predicted_cuts[id]["DC3"],predicted_cuts[id]["DR3"] = get_cuts_from_line(predictions)
    f.close()
    return predicted_cuts

def test_check_classification_report(cuts,max_positions):
    validation_labels = []
    pred_labels = []
    for id in cuts:
        v_l = ['O' for _ in range(0,250)]
        dr5,dc5,dc3,dr3 = cuts[id]
        v_l[dr5] = "DR5"
        v_l[dc5] = "DC5"
        v_l[dc3] = "DC3"
        v_l[dr3] = "DR3"
        p_l = ['O' for _ in range(0,250)]
        for c in ["DR5","DC5","DC3","DR3"]:
            p_l[max_positions[id][c]] = c
        validation_labels.append(v_l)
        pred_labels.append(p_l)
    #print(validation_labels)
    #print(pred_labels)
    print(classification_report(validation_labels,pred_labels))


def read_predicted_values_file(predicted_values_file):
    COL_START = 2
    f = open(predicted_values_file,"r")
    headers = f.readline().rstrip().split()
    if headers[0][0] != "#":
        print("Error: first line in %s must begin with a header starting with #\n"%(predicted_values_file))
        exit()
    data = {}
    classes = headers[COL_START:]
    for line in f.readlines():
        row_data = line.rstrip().split()
        id = row_data[0]
        data[id] = {}
        class_data = row_data[COL_START:]
        for j,c in enumerate(classes):
            data[id][c] = [float(x) for x in class_data[j].split(",")]
    f.close()
    return data

def get_max_dvs(data):
    max_dvs = {}
    max_positions = {}
    for id in list(data):
        max_dvs[id] = {}
        max_positions[id] = {}
        for c in list(data[id]):
            if c != 'O':
                max_positions[id][c] = data[id][c].index(max(data[id][c]))
                max_dvs[id][c] = data[id][c][max_positions[id][c]]
    return max_dvs, max_positions

def get_cutsite_dvs(cuts,data):
    cutsite_dvs = {}
    for id in list(data):
        cutsite_dvs[id] = {}
        dr5,dc5,dc3,dr3 = cuts[id]
        cutsite_dvs[id]["DR5"] = data[id]["DR5"][dr5]
        cutsite_dvs[id]["DC5"] = data[id]["DC5"][dc5]
        cutsite_dvs[id]["DC3"] = data[id]["DC3"][dc3]
        cutsite_dvs[id]["DR3"] = data[id]["DR3"][dr3]
    return cutsite_dvs

def get_position_difference(cuts,max_positions):
    deltas = {}
    for id in list(cuts):
        deltas[id] = {}
        dr5,dc5,dc3,dr3 = cuts[id]
        deltas[id]["DR5"] = max_positions[id]["DR5"] - dr5
        deltas[id]["DC5"] = max_positions[id]["DC5"] - dc5
        deltas[id]["DC3"] = max_positions[id]["DC3"] - dc3
        deltas[id]["DR3"] = max_positions[id]["DR3"] - dr3
    return deltas;

def create_box_plot(data,labels,colors,dot_colors,yLabel,title,include_dots=True,output_file="boxplot.svg"):
    bp = plt.boxplot(data, labels = labels, zorder = 1, patch_artist=True, showfliers=False, medianprops=dict(color="black"))
    plt.rcParams["font.family"] = "Arial"
    plt.ylabel(yLabel, {'fontname':'Arial', 'size':'14'})
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
    if include_dots:
        for i in range(0,len(labels)):
            y = data[i]
            x = np.random.normal((i + 1),0.06,len(y))
            plt.plot(x, y, '.', markersize=2, zorder = 2, color = dot_colors[i], alpha=0.5)
    plt.xticks(fontsize=14);
    #plt.title(title)
    plt.savefig(output_file,bbox_inches='tight')
    plt.close()

def prepare_delta_boxplot_data(deltas,delta_range=[-9,9]):
    lb_key = "<" + str(delta_range[0])
    ub_key = ">" + str(delta_range[1])
    deltas_boxes = {}
    deltas_boxes["DR5"] = {}
    deltas_boxes["DC5"] = {}
    deltas_boxes["DC3"] = {}
    deltas_boxes["DR3"] = {}
    for cut in deltas_boxes:
        deltas_boxes[cut] = {}
        deltas_boxes[cut][lb_key] = 0
        deltas_boxes[cut][ub_key] = 0
        for i in range(delta_range[0],delta_range[1]+1):
            deltas_boxes[cut][i] = 0
    for id in deltas:
        for cut in deltas[id]:
            d = deltas[id][cut]
            if d < delta_range[0]:
                deltas_boxes[cut][lb_key] += 1
            elif d > delta_range[1]:
                deltas_boxes[cut][ub_key] += 1
            else:
                deltas_boxes[cut][d] += 1
    return deltas_boxes

def print_boxplot(delta_counts,delta_range=[-9,9],title="distance from annotated cutsite",output_file="dist.svg"):
    lb_key = "<" + str(delta_range[0])
    ub_key = ">" + str(delta_range[1])
    keys = [lb_key] + [i for i in range(delta_range[0],delta_range[1]+1)] + [ub_key]
    index = np.arange(len(keys))
    bar_width = 1
    rects = plt.bar(index + bar_width, [delta_counts[k] for k in keys], bar_width,
                    color='aqua',
                    edgecolor='black')
    plt.xticks(index + bar_width, keys)
    plt.ylabel('count', {'fontname':'Arial', 'size':'14'})
    plt.xlabel('delta', {'fontname':'Arial', 'size':'14'})
    plt.title(title)
    plt.savefig(output_file)
    plt.close()

def print_avg_distance(deltas,outputFile="stats.txt"):
    f = open(outputFile,"w")
    mean = {}
    median = {}
    mean['DR5'] = statistics.mean([abs(deltas[id]['DR5']) for id in list(deltas)])
    median['DR5'] = statistics.median([abs(deltas[id]['DR5']) for id in list(deltas)])
    mean['DC5'] = statistics.mean([abs(deltas[id]['DC5']) for id in list(deltas)])
    median['DC5'] = statistics.median([abs(deltas[id]['DC5']) for id in list(deltas)])
    mean['DC3'] = statistics.mean([abs(deltas[id]['DC3']) for id in list(deltas)])
    median['DC3'] = statistics.median([abs(deltas[id]['DC3']) for id in list(deltas)])
    mean['DR3'] = statistics.mean([abs(deltas[id]['DR3']) for id in list(deltas)])
    median['DR3'] = statistics.median([abs(deltas[id]['DR3']) for id in list(deltas)])
    f.write("Cutsite\tMean\tMedian\n")
    mean_sum = 0
    for c in ['DR5','DC5','DC3','DR3']:
        mean_sum += mean[c]
        f.write("%s\t%.03f\t%.03f\n"%(c,mean[c],median[c]))
    f.write("all\t%.03f"%(mean_sum/float(4)))
    f.close()

if __name__ == "__main__":
    if (len(sys.argv) < 4):
        print("USAGE: %s <dataset file> <decision values file> <output prefix>\n"%(sys.argv[0]))
        exit()
    cuts = read_cuts_file(sys.argv[1])
    data = read_predicted_values_file(sys.argv[2])
    output_prefix = sys.argv[3]
    max_dvs, max_positions = get_max_dvs(data)
    cutsite_dvs = get_cutsite_dvs(cuts,data)
    deltas = get_position_difference(cuts,max_positions)
    deltas_boxes = prepare_delta_boxplot_data(deltas)

    cutsites = ["DR5","DC5","DC3","DR3"]

    #print_avg_distance(deltas,output_prefix + "_stats.txt")
    cut_p_colors = {"DR5" : "#b2e3ff",
                    "DC5" : "#ffc798",
                    "DC3" : "#c7eac7",
                    "DR3" : "#fb9293"}
    cut_d_colors = {"DR5":"#1f77b4",
                    "DC5":"#ff7f0e",
                    "DC3":"#2ca02c",
                    "DR3":"#d62728"}
    label_aliases = {"DR5" : "Drosha\n5\' cutsite",
                     "DC5" : "Dicer\n5\' cutsite",
                     "DC3" : "Dicer\n3\' cutsite",
                     "DR3" : "Drosha\n3\' cutsite"}

    colors = [cut_p_colors[c] for c in cutsites]
    dot_colors = [cut_d_colors[c] for c in cutsites]
    labels = [label_aliases[c] for c in cutsites]

    max_dv_data = [np.asarray([max_dvs[id]['DR5'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DC5'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DC3'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DR3'] for id in list(max_dvs)])]
    create_box_plot(max_dv_data,labels,colors,dot_colors,yLabel="Decision Value at Predicted Cutsite\n(Ensemble Applied to Test Set)",title="Maximum Decision Values for Cleavage-Site Predictions",include_dots=DVMAX_GRAPH_INCLUDE_DOTS,output_file=output_prefix + "_max_dv_boxPlot.svg")

    cutsite_dv_data = [np.asarray([cutsite_dvs[id]['DR5'] for id in list(cutsite_dvs)]),
                       np.asarray([cutsite_dvs[id]['DC5'] for id in list(cutsite_dvs)]),
                       np.asarray([cutsite_dvs[id]['DC3'] for id in list(cutsite_dvs)]),
                       np.asarray([cutsite_dvs[id]['DR3'] for id in list(cutsite_dvs)])]
    create_box_plot(cutsite_dv_data,labels,colors,dot_colors,yLabel="DV",title="Decision Values at Cleavage-Sites",include_dots=DVCUTS_GRAPH_INCLUDE_DOTS,output_file=output_prefix + "_cutsite_dv_boxPlot.svg")

    delta_data = [np.asarray([deltas[id]['DR5'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DC5'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DC3'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DR3'] for id in list(deltas)])]
    create_box_plot(delta_data,labels,colors,dot_colors,yLabel="delta",title="Distances between Annotated and Predicted Cleavage-sites",include_dots=DELTA_GRAPH_INCLUDE_DOTS,output_file=output_prefix + "_deltas_boxPlot.svg")

    
    for cut in deltas_boxes:
        print(cut)
        print(deltas_boxes[cut])
        print_boxplot(deltas_boxes[cut],delta_range=[-9,9],title="Predicted vs Annotated Cleavage-Site Positions ("+cut+")", output_file=output_prefix + "_deltas_"+cut+".svg")

    #test_check_classification_report(cuts,max_positions)


