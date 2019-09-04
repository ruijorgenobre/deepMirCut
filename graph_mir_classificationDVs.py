
import sys, os
import matplotlib.pyplot as plt
import deepMirCut

def deduct_padding(seq):
    in_3p_padding = True
    seq_start = 0
    seq_stop = len(seq) - 1
    for i in range(0,len(seq)):
        if seq[i].upper() == 'P' and in_3p_padding == True:
            seq_start += 1
        elif in_3p_padding == True:
            in_3p_padding = False
        elif seq[i].upper() == 'P' and in_3p_padding == False:
            seq_stop = i - 1
            return seq_start,seq_stop
    return seq_start,seq_stop

def deduct_padding_from_dvs(dv_data,seq_start,seq_stop):
    new_dv_data = {}
    for c in dv_data:
        new_dv_data[c] = dv_data[c][seq_start:seq_stop+1]
    return new_dv_data

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

def read_decision_values_file(predicted_values_file):
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

def read_predicted_labels_file(predicted_labels_file):
    prediction_info = {}
    f = open(predicted_labels_file,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,mir_name,seq,fold,labels,predicted_labels = row_data[0:6]
        labels = [x for x in labels.split(",")]
        predicted_labels = [x for x in predicted_labels.split(",")]
        prediction_info[id] = [mir_name,seq,fold,labels,predicted_labels]
    f.close()
    return prediction_info

if __name__ == "__main__":
    if (len(sys.argv) < 6):
        print("USAGE: %s <dataset file> <decision values file> <predicted labels file> <mir id> <output prefix>\n"%(sys.argv[0]))
        exit()
    print(sys.argv)
    cuts = read_cuts_file(sys.argv[1])
    data = read_decision_values_file(sys.argv[2])
    prediction_info = read_predicted_labels_file(sys.argv[3])
    id = sys.argv[4]
    output_prefix = sys.argv[5]
    if id not in data:
        print("Error: %s not found in data"%(id))
        exit()

    mir_name,seq = prediction_info[id][0:2]
    seq_start,seq_stop = deduct_padding(seq)
    dv_data = deduct_padding_from_dvs(data[id],seq_start,seq_stop)
    new_cuts = [c - seq_start for c in cuts[id]]
    
    handles = []
    for xc in new_cuts:
        plt.vlines(x=xc,ymin=0,ymax=1,color='lightgray')
    for c in dv_data:
        if c != "O":
            h, = plt.plot(dv_data[c],label=c)
            handles.append(h)
    plt.legend(handles=handles)
    plt.xlabel("Position")
    plt.ylabel("Decision Value")
    plt.savefig(output_prefix + "_" + id + "_graph_classificationDVs.svg")
    plt.close()
