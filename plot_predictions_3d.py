import sys, os
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FixedLocator, FixedFormatter


def deduct_padding_from_dvs(dv_data,seq_len):
    new_dv_data = {}
    for c in dv_data:
        new_dv_data[c] = dv_data[c][0:seq_len]
    return new_dv_data

def get_seq_lengths(sequence_file):
    seq_lengths = {}
    f = open(sequence_file,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,name,seq = row_data[0:3]
        seq_lengths[id] = len(seq)
    f.close()
    return seq_lengths

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


if __name__ == "__main__":

    if (len(sys.argv) < 4):
        print("USAGE: %s <classification DVs file> <sequences file> <id> [output prefix]\n"%(sys.argv[0]))
        exit()
    data = read_decision_values_file(sys.argv[1])
    seq_lengths = get_seq_lengths(sys.argv[2])
    id = sys.argv[3]
    output_prefix = id
    if len(sys.argv) > 4:
        output_prefix = sys.argv[4]
    if id not in data:
        print("Error: %s not found in data"%(id))
        exit()

    colors = {"DR5":"#1f77b4",
              "DC5":"#ff7f0e",
              "DC3":"#2ca02c",
              "DR3":"#d62728"}
        
    dv_data = deduct_padding_from_dvs(data[id],seq_lengths[id])


    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_ylim([0,4])
    ax.set_zlim([0,1])
    ax.set_xlim([0,len(dv_data['O'])])

    ax.azim = -95
    count = 4

    #Plotting dv's for each cutsite
    for c in reversed([k for k in dv_data]):
        if c != "O":
            count -= 1
            #print(c,count)
            ys = [count for _ in range(len(dv_data[c])) ]
            #x values are a range from 0 to the number of decision values for that cutsite
            #y values are set as a single for each x value  {i.e. DR5 = [0,0,0,0,0,0...0], DC5 = [1,1,1,1,1,1...1]}
            #z values are the decision values
            h, = plt.plot(range(len(dv_data[c])),ys,dv_data[c],color=colors[c],label=c)
            #handles.append(h)
    ax.set_yticks(range(0,4))
    ax.set_yticklabels([c for c in dv_data if c != 'O'])
    ax.set_zlabel("Decision Value")
    plt.savefig(output_prefix + "_graph_3d.svg")
    plt.close()
