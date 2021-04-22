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

def open_single_model_pmf_file(pmf_file):
    pmfScores = {}
    for c in ["DR5","DC5","DC3","DR3","average"]:
        pmfScores[c] = []
    f = open(pmf_file,"r")
    for line in f.readlines():
        model,DR5_f1,DC5_f1,DC3_f1,DR3_f1,avg = line.rstrip().split()
        if model != "#model":
            pmfScores["DR5"].append(float(DR5_f1))
            pmfScores["DC5"].append(float(DC5_f1))
            pmfScores["DC3"].append(float(DC3_f1))
            pmfScores["DR3"].append(float(DR3_f1))
            pmfScores["average"].append(float(avg))
    f.close()
    return(pmfScores)

def open_ensemble_pmf_file(pmf_file):
    pmfScores = {}
    for c in ["DR5","DC5","DC3","DR3","average"]:
        pmfScores[c] = []
    num_models = []
    f = open(pmf_file,"r")
    for line in f.readlines():
        members,DR5_f1,DC5_f1,DC3_f1,DR3_f1,avg = line.rstrip().split()
        if members != "#members":
            pmfScores["DR5"].append(float(DR5_f1))
            pmfScores["DC5"].append(float(DC5_f1))
            pmfScores["DC3"].append(float(DC3_f1))
            pmfScores["DR3"].append(float(DR3_f1))
            pmfScores["average"].append(float(avg))
            num_models.append(int(members));
    f.close()
    return(pmfScores,num_models)

if __name__ == "__main__":
    usage = "USAGE:\n%s <unweighted ensemble pmf file> [single model pmf file]\n"%sys.argv[0]
    if len(sys.argv) < 2:
        print(usage)
        exit()

    add_single_model = False
    if len(sys.argv) > 2:
        pmf_file_sm = sys.argv[2]
        add_single_model = True
    
    pmf_file = sys.argv[1]
    pmfScores,num_models = open_ensemble_pmf_file(pmf_file)
    print(pmfScores)
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

    font_size = 22
    font_size_ticks = 16
    font_size_legend = 16


    fig, ax = plt.subplots()
    ax.set_ylim(bottom=0.1,top=0.5)
    ax.set_xlabel("ensemble size",fontsize=font_size)
    ax.set_ylabel("Perfect Match Fraction\n(Ensembles Applied to Validation Set)",fontsize=font_size)

    linewidth=1
    
    cutsites = ["DR5","DC5","DC3","DR3","average"]
    for c in cutsites:
        ax.plot(num_models,pmfScores[c],linewidth=linewidth,color=colors[c])
    if add_single_model == True:
        pmfScores_single = open_single_model_pmf_file(pmf_file_sm)
        for c in cutsites:
            ax.scatter(num_models,pmfScores_single[c],s=5,color=colors[c])    

    legend_lines = [Line2D([0], [0], color=colors[c], linewidth=linewidth, linestyle='-', alpha=0.5) for c in cutsites]
    legend_labels = [c for c in cutsites]
    ax.legend(legend_lines, legend_labels, fontsize=font_size_legend)
    ax.set_xticks(num_models)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(font_size_ticks) 
    plt.savefig(pmf_file + ".svg",bbox_inches='tight')
    plt.close()
