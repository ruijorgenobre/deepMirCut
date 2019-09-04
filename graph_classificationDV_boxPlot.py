import sys, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import deepMirCut

DVMAX_GRAPH_INCLUDE_DOTS = False
DELTA_GRAPH_INCLUDE_DOTS = False


def read_cuts_file(dataSetFile):
    cuts = {}
    f = open(dataSetFile,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = row_data[0:12]
        drosha5p = deepMirCut.readCutSite(drosha5p)
        dicer5p = deepMirCut.readCutSite(dicer5p)
        dicer3p = deepMirCut.readCutSite(dicer3p)
        drosha3p = deepMirCut.readCutSite(drosha3p)
        cuts[id] = [drosha5p[0],dicer5p[0],dicer3p[0],drosha3p[0]]
    f.close()
    return cuts

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

def create_box_plot(data,labels,colors,dot_colors,yLabel,include_dots=True,output_file="boxplot.svg"):
    bp = plt.boxplot(data, labels = labels, zorder = 1, patch_artist=True, showfliers=False, medianprops=dict(color="black"))
    plt.ylabel(yLabel, {'fontname':'Arial', 'size':'16'})
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
    if include_dots:
        for i in range(0,len(labels)):
            y = data[i]
            x = np.random.normal((i + 1),0.06,len(y))
            plt.plot(x, y, '.', markersize=2, zorder = 2, color = dot_colors[i], alpha=0.5)
    plt.savefig(output_file)
    

if __name__ == "__main__":
    if (len(sys.argv) < 3):
        print("USAGE: %s <dataset file> <predicted values file> <output prefix>\n"%(sys.argv[0]))
        exit()
    print(sys.argv)
    cuts = read_cuts_file(sys.argv[1])
    data = read_predicted_values_file(sys.argv[2])
    output_prefix = sys.argv[3]
    max_dvs, max_positions = get_max_dvs(data)
    deltas = get_position_difference(cuts,max_positions)
    print(deltas)
    #print(deltas)

    colors = ['lightcyan','mistyrose','lavender','papayawhip']
    dot_colors = ['blue','red','darkviolet','orangered']
    labels = ['Drosha 5\' cutsite','Dicer 5\' cutsite','Dicer 3\' cutsite','Drosha 3\' cutsite']
    max_dv_data = [np.asarray([max_dvs[id]['DR5'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DC5'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DC3'] for id in list(max_dvs)]),
                   np.asarray([max_dvs[id]['DR3'] for id in list(max_dvs)])]

    create_box_plot(max_dv_data,labels,colors,dot_colors,yLabel="Maximum DV",include_dots=DVMAX_GRAPH_INCLUDE_DOTS,output_file=output_prefix + "_max_dv_boxPlot.svg")

    delta_data = [np.asarray([deltas[id]['DR5'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DC5'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DC3'] for id in list(deltas)]),
                  np.asarray([deltas[id]['DR3'] for id in list(deltas)])]

    create_box_plot(delta_data,labels,colors,dot_colors,yLabel="Difference in Predicted Position",include_dots=DELTA_GRAPH_INCLUDE_DOTS,output_file=output_prefix + "_deltas_boxPlot.svg")

