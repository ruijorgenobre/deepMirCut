import sys, os
import numpy as np
import pandas as pd
import tensorflow as tf
from keras.preprocessing.sequence import pad_sequences
from keras.utils import to_categorical
from keras.models import Model, Input
from keras.layers import CuDNNLSTM, LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional, Softmax
from keras_contrib.layers import CRF
from keras_contrib import losses
from keras_contrib import metrics
from keras import backend as K
from keras import activations
from keras import optimizers
from deepMirCut_metrics import avg_proximity_metric
from sklearn.model_selection import train_test_split
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
import matplotlib.pyplot as plt


os.environ["CUDA_VISIBLE_DEVICES"]="0";

INPUT_SETTING_LIST = ["USE_SEQ_ONLY","USE_SEQ_FOLD","USE_SEQ_BPRNA"];

DEFAULTS = {"batch_size" : 64,
            "use_embedding_layer" : True,
            "use_crf_layer" : False,
            "optimizer" : 'adam',
            "learning_rate" : 4.14e-3,
            "epsilon" : 10**-9.56,
            "embedding_layer_output" : 32,
            "embedding_dropout" : 0.230,
            "bi_lstm1_units" : 32,
            "bi_lstm2_units" : 160,
            "bi_lstm3_units" : 0,
            "input_setting" : "USE_SEQ_BPRNA",
            "max_seq_len" : 250,
            "model" : "model.h5",
            "epochs" : 100,
            "verbose" : False,
            "patience" : 3,
            "label_by_max_position" : True
}

def str2bool(v):
    return v.lower() in ("yes", "y", "true", "t", "1")

def load_class_labels(parameters):
    tags = ['O','DR5','DC5','DC3','DR3']
    tag_idx = {x:y for y,x in enumerate(tags)}
    parameters["tags"] = tags
    parameters["tag_idx"] = tag_idx
    parameters["num_tags"] = len(tags)
    parameters["DR5P_LAB"] = tags[1]
    parameters["DC5P_LAB"] = tags[2]
    parameters["DC3P_LAB"] = tags[3]
    parameters["DR3P_LAB"] = tags[4]
    return parameters

def load_seq_tags(parameters,USE_EMBEDDING_LAYER=True,INPUT_SETTING="USE_SEQ_BPRNA"):
    nucs = ['A','C','G','U']
    if USE_EMBEDDING_LAYER:
        nucs = ['A','C','G','U','N']
    nuc_fld = ['(',')','.']
    bp_fld = ['(',')','H','M','X','I','B','E']
    nuc_list = ['PAD'] + [x for x in nucs]
    if INPUT_SETTING == "USE_SEQ_FOLD":
        nuc_list = ['PAD'] + [x + y for x in nucs for y in nuc_fld]
    elif INPUT_SETTING == "USE_SEQ_BPRNA":
        nuc_list = ['PAD'] + [x + y for x in nucs for y in bp_fld]
    one_hot_vector_len = len(nuc_list)
    nuc_idx = {x:y for y,x in enumerate(nuc_list)}
    parameters["one_hot_vector_len"] = one_hot_vector_len
    parameters["nuc_idx"] =  nuc_idx
    return parameters

def load_neural_network_parameters(parameters,opts={}):
    parameters["use_embedding_layer"] = str2bool(opts['--use_embedding_layer']) if '--use_embedding_layer' in opts else DEFAULTS["use_embedding_layer"]
    parameters["use_crf_layer"] = str2bool(opts['--use_crf_layer']) if '--use_crf_layer' in opts else DEFAULTS["use_crf_layer"]
    parameters["optimizer"] =  opts['--optimizer'] if '--optimizer' in opts else DEFAULTS["optimizer"]
    parameters["learning_rate"] =  float(opts['--learning_rate']) if '--learning_rate' in opts else DEFAULTS["learning_rate"]
    parameters["epsilon"] =  float(opts['--epsilon']) if '--epsilon' in opts else DEFAULTS["epsilon"]
    parameters["batch_size"] =  int(opts['--batch_size']) if '--batch_size' in opts else DEFAULTS["batch_size"]
    parameters["embedding_layer_output"] =  int(opts['--embedding_layer_output']) if '--embedding_layer_output' in opts else DEFAULTS["embedding_layer_output"]
    parameters["embedding_dropout"] = float(opts['--embedding_dropout']) if '--embedding_dropout' in opts else DEFAULTS["embedding_dropout"]
    parameters["bi_lstm1_units"] = int(opts['--bi_lstm1_units']) if '--bi_lstm1_units' in opts else DEFAULTS["bi_lstm1_units"]
    parameters["bi_lstm2_units"] = int(opts['--bi_lstm2_units']) if '--bi_lstm2_units' in opts else DEFAULTS["bi_lstm2_units"]
    parameters["bi_lstm3_units"] = int(opts['--bi_lstm3_units']) if '--bi_lstm3_units' in opts else DEFAULTS["bi_lstm3_units"]
    return parameters

def load_ensemble(parameters,ensemble_list):
    parameters["ensemble"] = []
    parameters["ensemble_weights"] = []
    f = open(ensemble_list,"r")
    ensemble_info = []
    last_num_entries = -1
    for line in f.readlines():
        line_entries = line.rstrip().split()
        if last_num_entries == -1:
            last_num_entries = len(line_entries)
        if len(line_entries) > last_num_entries:
            print("Error: the number of weights do not match for each line in %s.  Make sure each entry is seperated by a single tab\n" % ensemble_list)
            exit()
        if len(line_entries) > 6:
            print("Error: too many weights for line: %s\nMake sure each entry is seperated by a single tab\n" % str(line_entries))
            exit()
        if len(line_entries) <= 4 and len(line_entries) > 1:
            print("Error: not enought weights for line: %s\nMake sure each entry is seperated by a single tab\n" % str(line_entries))
            exit()
        if len(line_entries) == 1:
            print("warning: no weights found.  Equal weights will be applied accross all models.")
        ensemble_info.append(line_entries)
    for i in range(0,len(ensemble_info)):
        if len(ensemble_info[i]) == 1:
            parameters["ensemble"].append(str(ensemble_info[i][0]))
            parameters["ensemble_weights"].append([1/len(ensemble_info) for _ in range(0,5)])
        elif len(ensemble_info[i]) == 5:
            parameters["ensemble"].append(str(ensemble_info[i][0]))
            parameters["ensemble_weights"].append([1/len(ensemble_info)] + [float(ensemble_info[i][j]) for j in range(1,len(ensemble_info[i]))])
        else:
            parameters["ensemble"].append(str(ensemble_info[i][0]))
            parameters["ensemble_weights"].append([float(ensemble_info[i][j]) for j in range(1,len(ensemble_info[i]))])
    f.close()

def load_parameters(opts={}):
    parameters = {}
    input_settings = {str(x):setting for x,setting in enumerate(INPUT_SETTING_LIST)}
    input_settings.update({v:v for v in input_settings.values()})
    parameters["input_setting"] = input_settings[opts['--input_setting']] if '--input_setting' in opts else DEFAULTS["input_setting"]
    parameters["max_seq_len"] = int(opts['--max_seq_len']) if '--max_seq_len' in opts else DEFAULTS["max_seq_len"]
    #parameters["model"] = opts['--model'] if '--model' in opts else DEFAULTS["model"]
    if '--ensemble_list' in opts:
        load_ensemble(parameters,opts['--ensemble_list'])
    elif '--model' in opts:
        parameters["model"] = opts['--model']
    parameters["epochs"] = int(opts['--epochs']) if '--epochs' in opts else DEFAULTS["epochs"]
    parameters["verbose"] = True if '--verbose' in opts else DEFAULTS["verbose"]
    parameters["patience"] = int(opts['--patience']) if '--patience' in opts else DEFAULTS["patience"]
    parameters["label_by_max_position"] = int(opts['--label_by_max_position']) if '--label_by_max_position' in opts else DEFAULTS["label_by_max_position"]
    parameters = load_neural_network_parameters(parameters,opts)
    parameters = load_seq_tags(parameters,USE_EMBEDDING_LAYER=parameters["use_embedding_layer"],INPUT_SETTING=parameters["input_setting"])
    parameters = load_class_labels(parameters)
    return parameters

def readCutSite(cutSite):
    if cutSite == '-':
        cutSite = ["-","-"]
    else:
        cutSite = [int(x)-1 for x in cutSite.split(',')]
    return cutSite

def sanitizeSequence(seq):
    newSeq = ""
    for n in seq:
        n = n.upper()
        if n == 'T':
            n = 'U'
        if n != 'A' and n != 'C' and n != 'G' and n != 'U':
            n = 'N'
        newSeq += n
    return newSeq

def add_context(fold,bpRNA_context):
    new_fold = []
    for i in range(0,len(bpRNA_context)):
        if fold[i] != '(' and fold[i] != ')':
            new_fold.append(bpRNA_context[i])
        else:
            new_fold.append(fold[i])
    return "".join(new_fold)

def readDataset(dataSetFile,parameters):
    dataSet = []
    f = open(dataSetFile,"r")
    for line in f.readlines():
        row_data = line.rstrip().split()
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = row_data[0:12]
        seq = sanitizeSequence(seq)
        drosha5p = readCutSite(drosha5p)
        dicer5p = readCutSite(dicer5p)
        dicer3p = readCutSite(dicer3p)
        drosha3p = readCutSite(drosha3p)
        if hpStart != '-':
            hpStart = int(hpStart) - 1
        if hpStop != '-':
            hpStop = int(hpStop) - 1
        if parameters["input_setting"] == "USE_SEQ_ONLY":
            dataSet.append([id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,int(hpStart),int(hpStop),seq])
        elif parameters["input_setting"] == "USE_SEQ_FOLD":
            fold = row_data[12]
            dataSet.append([id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,int(hpStart),int(hpStop),seq,fold])            
        elif parameters["input_setting"] == "USE_SEQ_BPRNA":
            fold,bpRNA_context = row_data[12:14]
            fold = add_context(fold,bpRNA_context)
            dataSet.append([id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,int(hpStart),int(hpStop),seq,fold])
        else:
            print("Error: INPUT_SETTING parameter is set to an invalid setting.  Valid settings are (\"USE_SEQ_ONLY\",\"USE_SEQ_FOLD\",\"USE_SEQ_BPRNA\").")
            exit()
    f.close()
    return dataSet

def dropLongSequences(dataSet,parameters):
    newDataSet = []
    long_seq_count = 0
    if parameters["input_setting"] == "USE_SEQ_ONLY":
        for id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq in dataSet:
            if len(seq) <= parameters["max_seq_len"]:
                newDataSet.append([id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq])
            else:
                long_seq_count += 1
    else:
        for id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq,fold in dataSet:
            if len(seq) <= parameters["max_seq_len"]:
                newDataSet.append([id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq,fold])
            else:
                long_seq_count += 1
    if long_seq_count > 0:
        print("Warning: %d sequences were greater than %d and dropped" % (long_seq_count, parameters["max_seq_len"]))
    return newDataSet

def prepareSeqOneHot(seq,fold,parameters):
    nucs = ['A','C','G','U']
    w = [seq[i]+fold[i] for i in range(0,len(seq))]
    oneHotSeq = [[0 for i in range(0,parameters["one_hot_vector_len"])] for i in range(0,len(w))]
    if len(seq) != len(fold):
        print("Error: length of sequence and fold are not equal\n%s\n%s" %(seq,fold))
        exit()
    for i,nuc in enumerate(w):
        if nuc in parameters["nuc_idx"]:
            oneHotSeq[i][parameters["nuc_idx"][nuc]] = 1
        else:
            for nNuc in [parameters["nuc_idx"][x+fold[i]] for x in nucs]:
                oneHotSeq[i][nNuc] = 0.25
    return oneHotSeq

def prepareLabels(seq,drosha5p,dicer5p,dicer3p,drosha3p,parameters):
    labels = [parameters["tag_idx"]['O'] for i in range(0,len(seq))];
    labels[drosha5p[0]] = parameters["tag_idx"][parameters["DR5P_LAB"]]
    labels[dicer5p[0]] = parameters["tag_idx"][parameters["DC5P_LAB"]]
    labels[dicer3p[0]] = parameters["tag_idx"][parameters["DC3P_LAB"]]
    labels[drosha3p[0]] = parameters["tag_idx"][parameters["DR3P_LAB"]]
    return labels

def prepareData(dataSet,parameters):
    X = []
    Y = []
    for row_data in dataSet:
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = row_data[0:12]
        if parameters['input_setting'] != "USE_SEQ_ONLY":
            fold = row_data[12]
        if parameters["use_embedding_layer"]:
            if parameters['input_setting'] == "USE_SEQ_ONLY":
                seqVals = [parameters["nuc_idx"][seq[i]] for i in range(0,len(seq))]
            else:
                seqVals = [parameters["nuc_idx"][seq[i] + fold[i]] for i in range(0,len(seq))]
            X.append(seqVals)
        else:
            oneHotSeq = prepareSeqOneHot(seq,fold,parameters)
            X.append(oneHotSeq)
        labels = prepareLabels(seq,drosha5p,dicer5p,dicer3p,drosha3p,parameters)
        Y.append(labels)
    if parameters["use_embedding_layer"]:
        X = pad_sequences(maxlen=parameters["max_seq_len"], sequences=X, padding="post", value=parameters["nuc_idx"]['PAD'])
    else:
        xPadVal = [0 for i in range(0,parameters["one_hot_vector_len"])]
        xPadVal[parameters["nuc_idx"]['PAD']] = 1
        X = pad_sequences(maxlen=parameters["max_seq_len"], sequences=X, padding="post", value=xPadVal)
    Y = pad_sequences(maxlen=parameters["max_seq_len"], sequences=Y, padding="post", value=parameters["tag_idx"]['O'])
    Y = [to_categorical(i, num_classes=parameters["num_tags"]) for i in Y]
    return X,Y

def createArchitecture(parameters):

    optimizer = 0
    if parameters["optimizer"] == 'rmsprop':
        optimizer = optimizers.rmsprop(lr=parameters["learning_rate"],epsilon=parameters["epsilon"])
    elif parameters["optimizer"] == 'adam':
        optimizer = optimizers.adam(lr=parameters["learning_rate"],epsilon=parameters["epsilon"])
    elif parameters["optimizer"] == 'nadam':
        optimizer = optimizers.nadam(lr=parameters["learning_rate"],epsilon=parameters["epsilon"])
    elif parameters["optimizer"] == 'sgd':
        optimizer = optimizers.sgd(lr=parameters["learning_rate"])
    #else:
    #    optimizer = parameters["optimizer"]

    
    if parameters["use_embedding_layer"]:
        input = Input(shape=(parameters["max_seq_len"],))
        model = Embedding(input_dim=parameters["one_hot_vector_len"], output_dim=parameters["embedding_layer_output"],input_length=parameters["max_seq_len"])(input)
        if parameters["embedding_dropout"] > 0:
            model = Dropout(rate = parameters["embedding_dropout"])(model)
    else:
        input = Input(shape=(parameters["max_seq_len"],parameters["one_hot_vector_len"]))
        model = input
    if parameters["bi_lstm1_units"] > 0:
        model = Bidirectional(CuDNNLSTM(units=parameters["bi_lstm1_units"], return_sequences=True))(model)
    if parameters["bi_lstm2_units"] > 0:
        model = Bidirectional(CuDNNLSTM(units=parameters["bi_lstm2_units"], return_sequences=True))(model)
    if parameters["bi_lstm3_units"] > 0:
        model = Bidirectional(CuDNNLSTM(units=parameters["bi_lstm3_units"], return_sequences=True))(model)
    if parameters["use_crf_layer"]:
        crf = CRF(parameters["num_tags"],learn_mode="marginal")
        out = crf(model)  # output
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss=losses.crf_loss, metrics=[metrics.crf_accuracy,avg_proximity_metric()])
    else:
        out = TimeDistributed(Dense(parameters["num_tags"], activation="softmax"))(model)
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss="categorical_crossentropy", metrics=["accuracy",avg_proximity_metric()])
    model.summary()
    return model

def pred2label(pred,parameters):
    idx2tag = {i: w for w, i in parameters["tag_idx"].items()}
    out = []

    for pred_i in pred:
        out_i = []
        if parameters["label_by_max_position"]:
            #print(pred_i.shape[1])
            max_peaks = [-float('inf') for i in range(0,pred_i.shape[1])]
            max_peak_pos = [0 for i in range(0,pred_i.shape[1])]
            for i in range(0,pred_i.shape[0]):
                p = pred_i[i]
                for j in range(0,len(p)):
                    if p[j] > max_peaks[j]:
                        max_peaks[j] = p[j]
                        max_peak_pos[j] = i
                out_i.append("O")
            for i in range(1,len(max_peak_pos)):
                #print(i)
                out_i[max_peak_pos[i]] = idx2tag[i]
            out.append(out_i)            
        else:
            for p in pred_i:
                p_i = np.argmax(p)
                out_i.append(idx2tag[p_i].replace("PAD","O"))
            out.append(out_i)
    return out

def get_classification_values(pred,parameters):
    idx2tag = {i: w for w, i in parameters["tag_idx"].items()}
    out = {}
    for i in range(0,pred.shape[2]):
        out[idx2tag[i]] = []
        for _ in range(0,pred.shape[0]):
            out[idx2tag[i]].append([])
    for i in range(0,pred.shape[0]):
        for j in range(0,pred.shape[1]):
            for k in range(0,pred.shape[2]):
                out[idx2tag[k]][i].append(float(pred[i][j][k]))
    return out

def input2label(input,parameters):
    idx2tag = {i: w for w, i in parameters["nuc_idx"].items()}
    out = []
    for input_i in input:
        out_i = []
        for n in input_i:
            if parameters["use_embedding_layer"]:
                out_i.append(idx2tag[n])
            else:
                out_i.append(idx2tag[np.argmax(n)])
        out.append(out_i)
    return out

def create_history_fig(history,acc_outputFile,loss_outputFile,prox_outputFile,parameters):
    if parameters["use_crf_layer"]:
        handle_crf_acc, = plt.plot(history.history["crf_accuracy"],label="train_crf_Acc")
        handle_val_crf_acc, = plt.plot(history.history["val_crf_accuracy"],label="val_crf_Acc")
        plt.legend(handles=[handle_crf_acc,handle_val_crf_acc])
        plt.savefig(acc_outputFile)
        plt.close()
    else:
        handle_acc, = plt.plot(history.history["acc"],label="train_Acc")
        handle_val_acc, = plt.plot(history.history["val_acc"],label="val_Acc")
        plt.legend(handles=[handle_acc,handle_val_acc])
        plt.savefig(acc_outputFile)
        plt.close()

    handle_loss, = plt.plot(history.history["loss"],label="loss")
    handle_val_loss, = plt.plot(history.history["val_loss"],label="val_loss")
    plt.legend(handles=[handle_loss,handle_val_loss])
    plt.savefig(loss_outputFile)
    plt.close()

    handle_prox, = plt.plot(history.history["prox"],label="prox")
    handle_val_prox, = plt.plot(history.history["val_prox"],label="val_prox")
    plt.legend(handles=[handle_prox,handle_val_prox])
    plt.savefig(prox_outputFile)
    plt.close()


def print_classification_report(validation_labels,pred_labels,outputFile=None):
    print("F1-score: {:.1%}".format(f1_score(validation_labels, pred_labels)))
    print(classification_report(validation_labels, pred_labels))
    if outputFile:
        f = open(outputFile,"w");
        f.write(classification_report(validation_labels, pred_labels))
        f.write("\n\n%s\n"%("F1-score: {:.1%}".format(f1_score(validation_labels, pred_labels))))
        f.close()

def breakInputSeq(inputLabels,parameters):
    inputSequences = []
    inputFolds = []
    for i in range(0,len(inputLabels)):
        inputSeq = []
        inputFold = []
        for j in range(0,len(inputLabels[i])):
            if inputLabels[i][j] == 'PAD':
                inputSeq.append('P')
                if parameters["input_setting"] != "USE_SEQ_ONLY":
                    inputFold.append('P')
            else:
                inputSeq.append(inputLabels[i][j][0])
                if parameters["input_setting"] != "USE_SEQ_ONLY":
                    inputFold.append(inputLabels[i][j][1])
        inputSequences.append(inputSeq)
        if parameters["input_setting"] == "USE_SEQ_ONLY":
            inputFolds.append("NONE")
        else:
            inputFolds.append(inputFold)
    return inputSequences, inputFolds

def print_validation_output_file(validationSet,X_vl,validation_labels,pred_labels,outputFile,parameters):
    f = open(outputFile,"w")
    inputLabels = input2label(X_vl,parameters)
    inputSequences,inputFolds = breakInputSeq(inputLabels,parameters)
    for i in range(0,len(validationSet)):
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = validationSet[i][0:12]
        f.write("%s\n"%("\t".join([id,name,"".join(inputSequences[i]),"".join(inputFolds[i]),",".join(validation_labels[i]),",".join(pred_labels[i])])))
    f.close();

def print_classification_values_file(validationSet,pred,parameters,output_file="classification_values.txt"):
    f = open(output_file,"w")
    classification_values = get_classification_values(pred,parameters)
    tag_str = "\t".join(parameters["tags"])
    f.write("#id\tname\t%s\n"%(tag_str))
    for i in range(0,len(validationSet)):
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq = validationSet[i][0:12]
        tag_str_list = []
        for t in parameters["tags"]:
            tag_str_list.append(",".join([str(x) for x in classification_values[t][i]]))
        f.write("%s\t%s\t%s\n"%(id,name,"\t".join(tag_str_list)))
    f.close()
    
if __name__ == '__main__':
    #if (len(sys.argv) < 3):
    #    print("USAGE: %s <train set> <validation set>\n" % (sys.argv[0]))
    #    exit()
    parameters = load_parameters(opts={})
    
    print(parameters)
    #trainSet = readDataset(sys.argv[1],parameters)
    #validationSet = readDataset(sys.argv[2],parameters)

    #new_trainSet = dropLongSequences(trainSet,parameters)
    #new_validationSet = dropLongSequences(validationSet,parameters)
    
    #X_tr,y_tr = prepareData(new_trainSet,parameters)
    #X_vl,y_vl = prepareData(new_validationSet,parameters)
    #validation_labels = pred2label(y_vl,parameters)
    
    model = createArchitecture(parameters)
    #history = model.fit(X_tr, np.array(y_tr),batch_size=128,epochs=1,validation_data=(X_vl,np.array(y_vl)),verbose=1,shuffle=True)

    #validation_pred = model.predict(X_vl, verbose=1)    
    #pred_labels = pred2label(validation_pred,parameters)

    #print_validation_output_file(new_validationSet,X_vl,validation_labels,pred_labels,"validation_predictions.txt",parameters) 

    #create_history_fig(history,"acc.svg","loss.svg")

    #print_classification_report(validation_labels,pred_labels,outputFile="classificationReport.txt")
    
    #print("F1-score: {:.1%}".format(f1_score(validation_labels, pred_labels)))
    #print(classification_report(validation_labels, pred_labels))
