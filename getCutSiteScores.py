import deepMirCut
import sys, os
import numpy as np
import pandas as pd
from keras.models import load_model
from keras.models import Model, Input
from keras.layers import CuDNNLSTM, LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional
from keras_contrib.layers import CRF
from keras_contrib.utils import save_load_utils
from keras_contrib import losses
from keras_contrib import metrics
from deepMirCut_metrics import avg_proximity_metric
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
from scipy.special import logit
import getopt

def print_usage():
    print("USAGE: %s <testSet file> [-o OUTPUT_PREFIX] [-h]" %(sys.argv[0]));
    print("optional arguments:")
    print("-h, --help                  print usage");
    print("-o, --output_prefix         output prefix");
    print("-v, --verbose               verbose");
    print("-m, --model                 file containing weights for model");
    print("-d, --print_dvs             prints a list of decision values for label at each position");
    print("-i, --input_setting         set to one of the following");
    print("                            0 : USE_SEQ_ONLY")
    print("                            1 : USE_SEQ_FOLD")
    print("                            2 : USE_SEQ_BPRNA")
    print("")
    print("--max_seq_len               Maximum sequence length")
    print("--use_embedding_layer       Train with an embedding layer (True/False)")
    print("--use_crf_layer             Train with a crf layer (True/False)")
    print("--embedding_layer_output    Embedding layer output size")
    print("--embedding_dropout         Dropout for Embedding layer")
    print("--bi_lstm1_units            units for first bidirectional lstm")
    print("--bi_lstm2_units            units for second bidirectional lstm")
    print("--dense1_units              units for first dense layer")
    print("--dense2_units              units for second dense layer")

def get_tag_positions(labels):
    tag_pos = {}
    for i in range(0,len(labels)):
        if labels[i].upper() != 'O':
            tag_pos[labels[i]] = i
    return tag_pos

def convert_to_score(dv):
    min_dv = 0.00001
    max_dv = 0.99999
    if dv < min_dv:
        return logit(min_dv)
    elif dv > max_dv:
        return logit(max_dv)
    return logit(dv)

def print_cutsite_scores(validationSet,validation_labels,validation_pred,parameters):
    unsquash = False
    f = open(parameters["cutsite_scores"],"w")
    tags = [t for t in parameters["tags"] if t.upper() != 'O']
    f.write("#id\tname\t%s\n" % ("\t".join(tags)))
    classification_values = deepMirCut.get_classification_values(validation_pred,parameters)
    for i in range(0,len(validationSet)):
        id,name = validationSet[i][0:2]
        tag_pos = get_tag_positions(validation_labels[i])
        if unsquash:
            f.write("%s\t%s\t%s\n" %(id,name,"\t".join([str(convert_to_score(classification_values[t][i][tag_pos[t]])) for t in tags])))
        else:
            f.write("%s\t%s\t%s\n" %(id,name,"\t".join([str(classification_values[t][i][tag_pos[t]]) for t in tags])))
    
def load_input_output_file_parameters(parameters,opts={}):
    parameters["validation_file"] = opts['--validation_file'] if '--validation_file' in opts else None
    parameters["cutsite_scores"] = parameters["validation_file"] + "_cutsite_scores.txt"
    return parameters

if __name__ == "__main__":
    
    min_args = 2
    if (len(sys.argv) < min_args):
        print_usage()
        exit()
    try:
        validation_file = sys.argv[1]
        opts_array, args = getopt.getopt(sys.argv[min_args:], "hvdo:e:m:i:", ["help","verbose","print_dvs","output_prefix=","epochs=","model=","input_setting=","max_seq_len=",
                                                                               "use_embedding_layer=","use_crf_layer=","embedding_layer_output=","embedding_dropout=","bi_lstm1_units=","bi_lstm2_units=",
                                                                               "dense1_units=","dense2_units="])
        if len(args) > 0:
            print("%s is not an option\n"%(args[0]))
            print_usage()
            exit()
        longopts_map = {
            "-h" : "--help",
            "-d" : "--print_dvs",
            "-o" : "--output_prefix",
            "-e" : "--epochs",
            "-v" : "--verbose",
            "-m" : "--model",
            "-i" : "--input_setting",
        }
        opts = {longopts_map[o] : a for o,a in opts_array if o in longopts_map}
        opts.update({o : a for o,a in opts_array if o not in longopts_map})
        opts["--validation_file"] = validation_file
    except getopt.GetoptError as err:
        print(err)
        print_usage()
        exit()
    if '--help' in opts:
        print_usage()
        exit()
    
    parameters = deepMirCut.load_parameters(opts)
    parameters = load_input_output_file_parameters(parameters,opts)

    model = load_model(parameters["model"], custom_objects={'CRF': CRF, 'crf_loss': losses.crf_loss, 'crf_accuracy': metrics.crf_accuracy,'prox':avg_proximity_metric()})

    model.summary()
    validationSet = deepMirCut.readDataset(parameters["validation_file"],parameters)
    new_validationSet = deepMirCut.dropLongSequences(validationSet,parameters)
    X_vl,y_vl = deepMirCut.prepareData(new_validationSet,parameters)
    validation_labels = deepMirCut.pred2label(y_vl,parameters)
    validation_pred = model.predict(X_vl, verbose=parameters["verbose"])
    print_cutsite_scores(validationSet,validation_labels,validation_pred,parameters)
