import sys, os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
font_dirs = ['/raid1/home/bb/bellji/fonts', ]
font_files = fm.findSystemFonts(fontpaths=font_dirs)
font_list = fm.createFontList(font_files)
fm.fontManager.ttflist.extend(font_list)

def readCutSite(cutSite):
    if cutSite == '-':
        cutSite = ["-","-"]
    else:
        cutSite = [int(x)-1 for x in cutSite.split(',')]
    return cutSite

def read_cuts_file(dataSetFile):
    cuts = {}
    mir_names = {}
    f = open(dataSetFile,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = row_data[0:12]
        drosha5p = readCutSite(drosha5p)
        dicer5p = readCutSite(dicer5p)
        dicer3p = readCutSite(dicer3p)
        drosha3p = readCutSite(drosha3p)
        cuts[id] = [drosha5p[0],dicer5p[0],dicer3p[0],drosha3p[0]]
        mir_names[id] = name
    f.close()
    return cuts,mir_names

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

def read_rnaSeq_file(rnaSeq_file,add_one=False):
    rnaSeq = {}
    heights = {}
    locations = {}
    f = open(rnaSeq_file,"r")
    for line in f.readlines():
        id,loc,rnaSeq_line = line.rstrip().split()
        if add_one:
            rnaSeq[id] = [float(x) + 1 for x in rnaSeq_line.split(",")]
        else:
            rnaSeq[id] = [float(x) for x in rnaSeq_line.split(",")]
        locations[id] = loc
        heights[id] = max(rnaSeq[id])
    f.close()
    return rnaSeq,locations,heights

def sort_by_readstack(cuts,rnaSeq,side):
    max_reads = {}
    for id in cuts:
        DR5_cut,DC5_cut,DC3_cut,DR3_cut = [c for c in cuts[id]]
        if side == "5p":
            max_reads[id] = max(rnaSeq[id][DR5_cut+1:DC5_cut+2])
        elif side == "3p":
            max_reads[id] = max(rnaSeq[id][DC3_cut+1:DR3_cut+2])
        else:
            print("Error: side not recognized:",side)
            exit()
    ids = [i for i,_ in sorted(max_reads.items(), key=lambda x: x[1], reverse = True)]
    return ids
 
        

if __name__ == "__main__":
    rnaSeq_file = sys.argv[1]
    cuts,mir_names = read_cuts_file(sys.argv[2])
    data = read_predicted_values_file(sys.argv[3])
    side = sys.argv[4]
    rnaSeq,locations,heights = read_rnaSeq_file(rnaSeq_file)
    max_dvs, max_positions = get_max_dvs(data)

    output_file = rnaSeq_file + ".svg"




    cutsite_labels = {"DR5": 'Drosha 5′-arm cut site',
              "DC5": 'Dicer 5′-arm cut site',
              "DC3": 'Dicer 3′-arm cut site',
              "DR3": 'Drosha 3′-arm cut site'}

    annot_cutsite_labels = {"DR5": 'Drosha 5′-arm annotated cut site (miRBase v22.1)',
                            "DC5": 'Dicer 5′-arm annotated cut site (miRBase v22.1)',
                            "DC3": 'Dicer 3′-arm annotated cut site (miRBase v22.1)',
                            "DR3": 'Drosha 3′-arm annotated cut site (miRBase v22.1)'}
    if side == '5p':
        annot_cutsite_labels["DR5"] = 'Drosha 5′-arm unannotated cut site (miRPreprocess.pl)'
        annot_cutsite_labels["DC5"] = 'Dicer 5′-arm unannotated cut site (miRPreprocess.pl)'
    elif side == '3p':
        annot_cutsite_labels["DC3"] = 'Dicer 3′-arm unannotated cut site (miRPreprocess.pl)'
        annot_cutsite_labels["DR3"] = 'Drosha 3′-arm unannotated cut site (miRPreprocess.pl)'
    else:
        print("Error: side not recognized:",side)
        exit()

    colors = {"DR5":"#1f77b4",
              "DC5":"#ff7f0e",
              "DC3":"#2ca02c",
              "DR3":"#d62728"}

    tick_fontsize = 22

    ids = sort_by_readstack(cuts,rnaSeq,side)

    rc={'font.family':'Arial'}
    plt.rcParams["figure.figsize"] = [22,38/22 * len(ids)]
    plt.rcParams["font.family"] = "Arial"

    scale_exp = len(str(int(heights[ids[0]]))) + 1
    yticks = [1] + [10**i for i in range(2,scale_exp+1,2)]
    fig, axs = plt.subplots(len(ids))
    for i in range(0,len(ids)):
        id = ids[i]
        name = mir_names[id]
        loc = locations[id]
        print(name,loc);
        x = [i for i in range(0,len(rnaSeq[id]))]
        y = rnaSeq[id]
        DR5_cut,DC5_cut,DC3_cut,DR3_cut = [c+1/2 for c in cuts[id]]
        bar_width = 1
        axs[i].bar(x,y,width=1.0,facecolor='cyan', edgecolor='cyan')
        axs[i].axvline(x = DR5_cut, linewidth=1, color = colors['DR5']) 
        axs[i].axvline(x = DC5_cut, linewidth=1, color = colors['DC5']) 
        axs[i].axvline(x = DC3_cut, linewidth=1, color = colors['DC3']) 
        axs[i].axvline(x = DR3_cut, linewidth=1, color = colors['DR3'])
        ax2 = axs[i].twinx()
        ax2.set_yticks([])
        ax2.set_ylim(bottom=0,top=1)
        for cs in ["DR5","DC5","DC3","DR3"]:
            ax2.plot(max_positions[id][cs] + 1/2,-0.115,"^", markersize=16,color = colors[cs],clip_on=False)
        axs[i].set_yscale('symlog')
        axs[i].set_ylim(top=10**scale_exp)
        axs[i].set_yticks(yticks)
        axs[i].text(0.01, 0.60, name, fontsize=36, transform=axs[i].transAxes)
        for tick in axs[i].xaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize) 
        for tick in axs[i].yaxis.get_major_ticks():
            tick.label.set_fontsize(tick_fontsize) 
    plt.tight_layout()
    plt.text(-0.02, 0.5, "Raw Read Count, symlog", fontsize=32, rotation=90, va='center', transform=fig.transFigure)
    cutsites = ['DR5','DC5','DC3','DR3']
    legend_lines = [Line2D([0], [0], marker='^', markersize=32, color='w', markerfacecolor=colors[c]) for c in cutsites]
    legend_cutsite_labels = [cutsite_labels[c] for c in cutsites]
    legend_annot_cutsite_labels = [annot_cutsite_labels[c] for c in cutsites]
    #legend_lines_2 = [Line2D([0,0], [0,1], linestyle='-', linewidth=1, color=colors[c]) for c in cutsites]
    legend_lines_2 = [Line2D([0,0], [0,1], marker='|', markersize=32, markeredgewidth=1.5, linestyle='None', color=colors[c]) for c in cutsites]
    legend1 = plt.legend(legend_lines, legend_cutsite_labels,loc=(0, -1.8 * len(ids)/11),mode = "expand",ncol=len(cutsites),fontsize=32,handletextpad=0.05,
                         title="Predicted Cut Sites:", title_fontsize=32)
    legend1._legend_box.align = "left"
    legend2 = plt.legend(legend_lines_2, legend_annot_cutsite_labels, loc=(0, -5.3 * len(ids)/11),mode = "expand",ncol=1,fontsize=32,handletextpad=0.05, 
                         title="Target Cut Sites:", title_fontsize=32)
    legend2._legend_box.align = "left"
    plt.gca().add_artist(legend1)
    #plt.gca().add_artist(legend2)
    plt.savefig(output_file,bbox_inches='tight')


