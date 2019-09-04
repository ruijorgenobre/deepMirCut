import deepMirCut
import sys, os
import numpy as np
import pandas as pd
from keras_contrib.utils import save_load_utils
from keras.callbacks import EarlyStopping
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
import getopt

def print_usage():
    print("USAGE: %s <train file> [-e EPOCHS] [-o OUTPUT_PREFIX] [-h]" %(sys.argv[0]));
    print("optional arguments:")
    print("-h, --help                  print usage");
    print("-e, --epochs                number of epochs");
    print("-s, --validation_file       validation set file")
    print("-o, --output_prefix         output prefix");
    print("-v, --verbose               verbose");
    print("-m, --model                 file containing weights for model");
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

def load_input_output_file_parameters(parameters,opts={}):
    parameters["train_file"] = opts['--train_file']
    parameters["validation_file"] = opts['--validation_file'] if '--validation_file' in opts else None
    parameters["output_prefix"] = opts['--output_prefix'] if '--output_prefix' in opts else "results"
    parameters["predictions_file"] = parameters["output_prefix"] + "_validation_predictions.txt"
    parameters["results_file"] = parameters["output_prefix"] + "_classificationReport.txt"
    parameters["acc_fig"] = parameters["output_prefix"] + "_acc.svg"
    parameters["loss_fig"] = parameters["output_prefix"] + "_loss.svg"
    parameters["param_out_file"] = parameters["output_prefix"] + "_parameters.txt"
    return parameters

if __name__ == "__main__":
    
    min_args = 2
    if (len(sys.argv) < min_args):
        print_usage()
        exit()
    try:
        train_file = sys.argv[1]
        opts_array, args = getopt.getopt(sys.argv[min_args:], "hvo:e:s:m:i:p:", ["help","verbose","output_prefix=","epochs=","validation_file=","model=","input_setting=","max_seq_len=","patience=",
                                                                               "use_embedding_layer=","use_crf_layer=","embedding_layer_output=","embedding_dropout=","bi_lstm1_units=","bi_lstm2_units=",
                                                                               "dense1_units=","dense2_units="])
        if len(args) > 0:
            print("%s is not an option\n"%(args[0]))
            print_usage()
            exit()
        longopts_map = {
            "-h" : "--help",
            "-o" : "--output_prefix",
            "-e" : "--epochs",
            "-p" : "--patience",
            "-v" : "--verbose",
            "-s" : "--validation_file",
            "-m" : "--model",
            "-i" : "--input_setting",
        }
        opts = {longopts_map[o] : a for o,a in opts_array if o in longopts_map}
        opts.update({o : a for o,a in opts_array if o not in longopts_map})
        opts["--train_file"] = train_file
    except getopt.GetoptError as err:
        print(err)
        print_usage()
        exit()
    if '--help' in opts:
        print_usage()
        exit()
    
    parameters = deepMirCut.load_parameters(opts)
    parameters = load_input_output_file_parameters(parameters,opts)

    print("patience=",parameters["patience"])

    trainSet = deepMirCut.readDataset(parameters["train_file"],parameters)
    new_trainSet = deepMirCut.dropLongSequences(trainSet,parameters)
    X_tr,y_tr = deepMirCut.prepareData(new_trainSet,parameters)

    model = deepMirCut.createArchitecture(parameters)
    if parameters["validation_file"]:
        validationSet = deepMirCut.readDataset(parameters["validation_file"],parameters)
        new_validationSet = deepMirCut.dropLongSequences(validationSet,parameters)
        X_vl,y_vl = deepMirCut.prepareData(new_validationSet,parameters)
        validation_labels = deepMirCut.pred2label(y_vl,parameters)
        es = EarlyStopping(monitor='val_loss', min_delta=0, patience=parameters["patience"], verbose=parameters["verbose"], mode='min', restore_best_weights=True)
        history = model.fit(X_tr, np.array(y_tr),batch_size=128,epochs=parameters["epochs"],validation_data=(X_vl,np.array(y_vl)),verbose=parameters["verbose"],shuffle=True,callbacks=[es])
        validation_pred = model.predict(X_vl, verbose=parameters["verbose"])
        pred_labels = deepMirCut.pred2label(validation_pred,parameters)
        deepMirCut.print_validation_output_file(new_validationSet,X_vl,validation_labels,pred_labels,parameters["predictions_file"],parameters)
        deepMirCut.create_history_fig(history,parameters["acc_fig"],parameters["loss_fig"],parameters)
        deepMirCut.print_classification_report(validation_labels,pred_labels,outputFile=parameters["results_file"])
        print("F1-score: {:.1%}".format(f1_score(validation_labels, pred_labels)))
        print(classification_report(validation_labels, pred_labels))
    else:
        history = model.fit(X_tr, np.array(y_tr),batch_size=128,epochs=parameters["epochs"],verbose=parameters["verbose"],shuffle=True)        

    model.save(parameters["model"])
