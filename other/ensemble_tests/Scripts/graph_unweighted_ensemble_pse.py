import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm
font_dirs = ['/raid1/home/bb/bellji/fonts', ]
font_files = fm.findSystemFonts(fontpaths=font_dirs)
font_list = fm.createFontList(font_files)
fm.fontManager.ttflist.extend(font_list)

def open_single_model_pse_file(pse_file):
    pseScores = {}
    for c in ["DR5","DC5","DC3","DR3","average"]:
        pseScores[c] = []
    f = open(pse_file,"r")
    for line in f.readlines():
        model,DR5_f1,DC5_f1,DC3_f1,DR3_f1,avg = line.rstrip().split()
        if model != "#model":
            pseScores["DR5"].append(float(DR5_f1))
            pseScores["DC5"].append(float(DC5_f1))
            pseScores["DC3"].append(float(DC3_f1))
            pseScores["DR3"].append(float(DR3_f1))
            pseScores["average"].append(float(avg))
    f.close()
    return(pseScores)

def open_ensemble_pse_file(pse_file):
    pseScores = {}
    for c in ["DR5","DC5","DC3","DR3","average"]:
        pseScores[c] = []
    num_models = []
    f = open(pse_file,"r")
    for line in f.readlines():
        members,DR5_f1,DC5_f1,DC3_f1,DR3_f1,avg = line.rstrip().split()
        if members != "#members":
            pseScores["DR5"].append(float(DR5_f1))
            pseScores["DC5"].append(float(DC5_f1))
            pseScores["DC3"].append(float(DC3_f1))
            pseScores["DR3"].append(float(DR3_f1))
            pseScores["average"].append(float(avg))
            num_models.append(int(members));
    f.close()
    return(pseScores,num_models)

if __name__ == "__main__":
    usage = "USAGE:\n%s <unweighted ensemble pse file> [single model pse file]\n"%sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        exit()

    add_single_model = False
    if len(sys.argv) > 2:
        pse_file_sm = sys.argv[2]
        add_single_model = True
    
    pse_file = sys.argv[1]
    pseScores,num_models = open_ensemble_pse_file(pse_file)
    print(pseScores)
    print(num_models)


    labels = {"DR5": 'Drosha 5\' cutsite',
              "DC5": 'Dicer 5\' cutsite',
              "DC3": 'Dicer 3\' cutsite',
              "DR3": 'Drosha 3\' cutsite',
              "average": "Average"}

    colors = {"DR5":"#1f77b4",
              "DC5":"#ff7f0e",
              "DC3":"#2ca02c",
              "DR3":"#d62728",
              "average":"black"}

    plt.rcParams["font.family"] = "Arial"
    fig_size_mult = 1.35
    plt.rcParams["figure.figsize"][0] *= fig_size_mult
    plt.rcParams["figure.figsize"][1] *= fig_size_mult
    #plt.rcParams["figure.figsize"] = [8.5,8.5]

    #plt.ylim(bottom=0.2,top=0.5)

    #font_size = 14
    #font_size_ticks = 12
    font_size = 22
    font_size_ticks = 16
    font_size_legend = 16


    fig, ax = plt.subplots()
    ax.set_xlabel("ensemble size",fontsize=font_size)
    ax.set_ylabel("Positional Shift Error\n(Ensembles Applied to Validation Set)",fontsize=font_size)

    linewidth=1
    
    cutsites = ["DR5","DC5","DC3","DR3","average"]
    for c in cutsites:
        ax.plot(num_models,pseScores[c],linewidth=linewidth,color=colors[c])
    if add_single_model == True:
        pseScores_single = open_single_model_pse_file(pse_file_sm)
        for c in cutsites:
            ax.scatter(num_models,pseScores_single[c],s=5,color=colors[c])    

    legend_lines = [Line2D([0], [0], color=colors[c], linewidth=linewidth, linestyle='-', alpha=0.5) for c in cutsites]
    legend_labels = [c for c in cutsites]
    ax.legend(legend_lines, legend_labels,fontsize=font_size_legend)
    ax.set_xticks(num_models)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    plt.savefig(pse_file + ".svg",bbox_inches='tight')
    plt.close()
