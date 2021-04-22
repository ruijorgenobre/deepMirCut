import sys, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
import statistics
import functools
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
    plt.ylabel(yLabel, {'fontname':'Arial', 'size':'14'})
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
    if include_dots:
        for i in range(0,len(labels)):
            y = data[i]
            x = np.random.normal((i + 1),0.06,len(y))
            plt.plot(x, y, '.', markersize=2, zorder = 2, color = dot_colors[i], alpha=0.5)
    plt.title(title)
    plt.savefig(output_file)
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

def prepare_precision_lineplot_data(deltas,total_examples,max_dist=9):
    precision = {}
    precision["DR5"] = [0 for i in range(0,max_dist+1)]
    precision["DC5"] = [0 for i in range(0,max_dist+1)]
    precision["DC3"] = [0 for i in range(0,max_dist+1)]
    precision["DR3"] = [0 for i in range(0,max_dist+1)]
    for id in deltas:
        for cut in deltas[id]:
            d = abs(deltas[id][cut])
            if d <= max_dist:
                precision[cut][d] += 1
    for cut in precision:
        for i in range(0,max_dist):
            precision[cut][i+1] += precision[cut][i]
    for cut in precision:
        for i in range(0,len(precision[cut])):
            precision[cut][i] /= total_examples
    return precision    

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

def createHistogramData_old(delta_counts,cutsite,keys,barplot_like=False):
    hist_gap = 1
    hist_out_width = 1
    hist_len = len(keys[1:-1])
    hist_xVals = []
    hist_tick_pos = []
    if barplot_like:
        hist_xVals = functools.reduce(lambda x,y :x+y,[(i+hist_out_width+hist_gap,i+1+hist_out_width+hist_gap) for i in range(len(keys[1:-1]))])
        hist_yVals = functools.reduce(lambda x,y :x+y,[(delta_counts[cutsite][k],delta_counts[cutsite][k]) for k in keys[1:-1]])
        hist_tick_pos = [(i+hist_out_width+hist_gap)/2 + (i+1+hist_out_width+hist_gap) / 2 for i in range(len(keys[1:-1]))]
    else:
        hist_xVals = [i + hist_out_width + hist_gap for i in range(0,len(keys[1:-1]))]
        hist_yVals = [delta_counts[cutsite][k] for k in keys[1:-1]]
        hist_tick_pos = [i + hist_out_width + hist_gap for i in range(0,len(keys[1:-1]))]
    lb_xVals = [0,hist_out_width]
    lb_yVals = [delta_counts[cutsite][keys[0]],delta_counts[cutsite][keys[0]]]
    lb_tick_pos = [lb_xVals[0]/2 + lb_xVals[1]/2]
    ub_xVals = [hist_xVals[-1] + hist_gap, hist_xVals[-1] + hist_out_width + hist_gap]
    ub_yVals = [delta_counts[cutsite][keys[-1]],delta_counts[cutsite][keys[-1]]]
    ub_tick_pos = [ub_xVals[0]/2 + ub_xVals[1]/2]
    return hist_xVals,hist_yVals,hist_tick_pos,lb_xVals,lb_yVals,lb_tick_pos,ub_xVals,ub_yVals,ub_tick_pos

def createHistogramData(delta_counts,keys,barplot_like=False):
    hist_gap = 1
    hist_out_width = 1
    hist_xVals = []
    hist_tick_pos = []
    hist_yVals = {}
    lb_yVals = {}
    ub_yVals = {}
    for c in delta_counts:
        if barplot_like:
            hist_xVals = functools.reduce(lambda x,y :x+y,[(i+hist_out_width+hist_gap,i+1+hist_out_width+hist_gap) for i in range(len(keys[1:-1]))])
            hist_tick_pos = [(i+hist_out_width+hist_gap)/2 + (i+1+hist_out_width+hist_gap) / 2 for i in range(len(keys[1:-1]))]
        else:
            hist_xVals = [i + hist_out_width + hist_gap for i in range(0,len(keys[1:-1]))]
            hist_tick_pos = [i + hist_out_width + hist_gap for i in range(0,len(keys[1:-1]))]
        lb_xVals = [0,hist_out_width]
        lb_tick_pos = [lb_xVals[0]/2 + lb_xVals[1]/2]
        ub_xVals = [hist_xVals[-1] + hist_gap, hist_xVals[-1] + hist_out_width + hist_gap]
        ub_tick_pos = [ub_xVals[0]/2 + ub_xVals[1]/2]
    for c in delta_counts:
        hist_yVals[c] = []
        if barplot_like:
            hist_yVals[c] = functools.reduce(lambda x,y :x+y,[(delta_counts[c][k],delta_counts[c][k]) for k in keys[1:-1]])
        else:
            hist_yVals[c] = [delta_counts[c][k] for k in keys[1:-1]]
        lb_yVals[c] = [delta_counts[c][keys[0]],delta_counts[c][keys[0]]]
        ub_yVals[c] = [delta_counts[c][keys[-1]],delta_counts[c][keys[-1]]]
    return hist_xVals,hist_yVals,hist_tick_pos,lb_xVals,lb_yVals,lb_tick_pos,ub_xVals,ub_yVals,ub_tick_pos


def normalize_deltaCounts(delta_counts,size):
    normalized_delta_counts = {}
    for c in delta_counts:
        normalized_delta_counts[c] = {}
        for k in delta_counts[c]:
            normalized_delta_counts[c][k] = delta_counts[c][k] / size
    return normalized_delta_counts

def get_vertical_lines_lengths(hist_xVals,hist_yVals,lb_xVals,lb_yVals,ub_xVals,ub_yVals):
    vertical_line_lengths = {}
    for i in range(0,len(lb_xVals)):
        xVal = lb_xVals[i]
        for c in lb_yVals:
            if xVal not in vertical_line_lengths:
                vertical_line_lengths[xVal] = lb_yVals[c][i]
            elif vertical_line_lengths[xVal] > lb_yVals[c][i]:
                vertical_line_lengths[xVal] = lb_yVals[c][i]    
    for i in range(0,len(hist_xVals)):
        xVal = hist_xVals[i]
        for c in hist_yVals:
            if xVal not in vertical_line_lengths:
                vertical_line_lengths[xVal] = hist_yVals[c][i]
            elif vertical_line_lengths[xVal] > hist_yVals[c][i]:
                vertical_line_lengths[xVal] = hist_yVals[c][i]
    for i in range(0,len(ub_xVals)):
        xVal = ub_xVals[i]
        for c in ub_yVals:
            if xVal not in vertical_line_lengths:
                vertical_line_lengths[xVal] = ub_yVals[c][i]
            elif vertical_line_lengths[xVal] > ub_yVals[c][i]:
                vertical_line_lengths[xVal] = ub_yVals[c][i]    
    return vertical_line_lengths


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
    delta_range=[-9,9]
    delta_counts = prepare_delta_boxplot_data(deltas,delta_range=delta_range)
    precision = prepare_precision_lineplot_data(deltas,len(cuts),max_dist=9)

    #print_avg_distance(deltas,output_prefix + "_stats.txt")

    labels = {"DR5": 'Drosha 5\' cutsite',
              "DC5": 'Dicer 5\' cutsite',
              "DC3": 'Dicer 3\' cutsite',
              "DR3": 'Drosha 3\' cutsite'}

    colors = {"DR5":"#1f77b4",
              "DC5":"#ff7f0e",
              "DC3":"#2ca02c",
              "DR3":"#d62728"}


    lb_key = "<" + str(delta_range[0])
    ub_key = ">" + str(delta_range[1])
    cutsites = ["DR5","DC5","DC3","DR3"]
    keys = [lb_key] + [i for i in range(delta_range[0],delta_range[1]+1)] + [ub_key]

    delta_counts = normalize_deltaCounts(delta_counts,size=len(cuts))
    #get_background_lines(delta_counts)

    plt.rcParams["font.family"] = "Arial"
    fig_size_mult = 1.35
    plt.rcParams["figure.figsize"][0] *= fig_size_mult
    plt.rcParams["figure.figsize"][1] *= fig_size_mult
    #plt.rcParams["figure.figsize"] = [8.5,8.5]

    #font_size = 14
    font_size = 22
    font_size_ticks = 16
    font_size_legend = 16

    alpha = 0.7

    hist_xVals,hist_yVals,hist_tick_pos,lb_xVals,lb_yVals,lb_tick_pos,ub_xVals,ub_yVals,ub_tick_pos = createHistogramData(delta_counts,keys,barplot_like=True)
    vertical_line_lengths = get_vertical_lines_lengths(hist_xVals,hist_yVals,lb_xVals,lb_yVals,ub_xVals,ub_yVals)
    #print(vertical_line_lengths)
    cutsite_linewidth = 2
    dot_linewidth = 1
    fig, ax = plt.subplots()
    ax.set_xlabel("Positional Shift of Cutsite Predictions\nRelative to Annotations",fontsize=font_size)
    ax.set_ylabel("Frequency of Predictions\n(Ensemble Applied to Test Set)",fontsize=font_size)
    for x in vertical_line_lengths:
        ax.plot([int(x),int(x)],[0,vertical_line_lengths[x]],linewidth=0.5,linestyle='-',color="black");
    for c in cutsites:
        ax.plot(lb_xVals,lb_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot([lb_xVals[-1],lb_xVals[-1]],[0,lb_yVals[c][1]],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot([lb_xVals[-1],hist_xVals[0]],[lb_yVals[c][-1],hist_yVals[c][0]],linewidth=dot_linewidth,linestyle=':',color="white")
        ax.plot([hist_xVals[0],hist_xVals[0]],[0,hist_yVals[c][-1]],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot(hist_xVals,hist_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot([hist_xVals[-1],hist_xVals[-1]],[0,hist_yVals[c][-1]],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot([hist_xVals[-1],ub_xVals[0]],[hist_yVals[c][-1],ub_yVals[c][-1]],linewidth=dot_linewidth,linestyle=':',color="white")
        ax.plot([ub_xVals[0],ub_xVals[0]],[0,ub_yVals[c][0]],linewidth=cutsite_linewidth,linestyle='-',color="white")
        ax.plot(ub_xVals,ub_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color="white")
    for c in cutsites:
        ax.plot(lb_xVals,lb_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot([lb_xVals[-1],lb_xVals[-1]],[0,lb_yVals[c][1]],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot([lb_xVals[-1],hist_xVals[0]],[lb_yVals[c][-1],hist_yVals[c][0]],linewidth=dot_linewidth,linestyle=':',color=colors[c],alpha=alpha)
        ax.plot([hist_xVals[0],hist_xVals[0]],[0,hist_yVals[c][-1]],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot(hist_xVals,hist_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot([hist_xVals[-1],hist_xVals[-1]],[0,hist_yVals[c][-1]],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot([hist_xVals[-1],ub_xVals[0]],[hist_yVals[c][-1],ub_yVals[c][-1]],linewidth=dot_linewidth,linestyle=':',color=colors[c],alpha=alpha)
        ax.plot([ub_xVals[0],ub_xVals[0]],[0,ub_yVals[c][0]],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
        ax.plot(ub_xVals,ub_yVals[c],linewidth=cutsite_linewidth,linestyle='-',color=colors[c],alpha=alpha)
    legend_lines = [Line2D([0], [0], color=colors[c], linewidth=cutsite_linewidth, linestyle='-', alpha=alpha) for c in cutsites]
    legend_labels = [labels[c] for c in cutsites]
    ax.legend(legend_lines, legend_labels,fontsize=font_size_legend)
    tickLocations = lb_tick_pos + hist_tick_pos + ub_tick_pos
    #keys[0] = "less\nthan %s" % (str(delta_range[0]));
    #keys[-1] = "more\nthan %s" % (str(delta_range[1]));
    ax.set_xticks(tickLocations)
    ax.set_xticklabels(["less\nthan %s" % (str(delta_range[0]))] + keys[1:-1] + ["more\nthan %s" % (str(delta_range[1]))])
    ax.set_ylim(ymin=0)
    ax.set_xlim(xmin=0,xmax=ub_xVals[1])
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    plt.savefig(output_prefix + "_hist.svg",bbox_inches='tight')
    plt.close()


    plt.xlabel("Positional Shift Error",fontsize=font_size)
    plt.ylabel("Frequency of Predictions\n(Ensemble Applied to Test Set)",fontsize=font_size)
    for c in cutsites:
        plt.plot(range(0,len(precision[c])),precision[c],cutsite_linewidth,linestyle="-",color=colors[c])
    legend_lines = [Line2D([0], [0], color=colors[c], linewidth=cutsite_linewidth, linestyle='-') for c in cutsites]
    legend_labels = [labels[c] for c in cutsites]
    plt.xticks(range(0,len(precision["DR5"])),range(0,len(precision["DR5"])))
    plt.ylim(ymin=0,ymax=1)
    plt.xlim(xmin=0)
    plt.legend(legend_lines, legend_labels,loc='lower right',fontsize=font_size_legend)
    plt.savefig(output_prefix + "_fractionIdentified.svg",bbox_inches='tight')
