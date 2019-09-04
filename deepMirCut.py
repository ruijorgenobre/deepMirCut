import sys, os
import numpy as np
import pandas as pd
from keras.preprocessing.sequence import pad_sequences
from keras.utils import to_categorical
from keras.models import Model, Input
from keras.layers import CuDNNLSTM, LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional
from keras_contrib.layers import CRF
from keras_contrib import losses
from keras_contrib import metrics
from sklearn.model_selection import train_test_split
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
import matplotlib.pyplot as plt


os.environ["CUDA_VISIBLE_DEVICES"]="0";

INPUT_SETTING_LIST = ["USE_SEQ_ONLY","USE_SEQ_FOLD","USE_SEQ_BPRNA"];

DEFAULTS = {"use_embedding_layer" : True,
            "use_crf_layer" : False,
            "embedding_layer_output" : 64,
            "embedding_dropout" : 0.45,
            "bi_lstm1_units" : 64,
            "bi_lstm2_units" : 64,
            "dense1_units" : 64,
            "dense2_units" : 64,
            "input_setting" : "USE_SEQ_BPRNA",
            "max_seq_len" : 250,
            "model" : "model.h5",
            "epochs" : 20,
            "verbose" : False,
            "patience" : 1
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
    parameters["embedding_layer_output"] =  int(opts['--embedding_layer_output']) if '--embedding_layer_output' in opts else DEFAULTS["embedding_layer_output"]
    parameters["embedding_dropout"] = int(opts['--embedding_dropout']) if '--embedding_dropout' in opts else DEFAULTS["embedding_dropout"]
    parameters["bi_lstm1_units"] = int(opts['--bi_lstm1_units']) if '--bi_lstm1_units' in opts else DEFAULTS["bi_lstm1_units"]
    parameters["bi_lstm2_units"] = int(opts['--bi_lstm2_units']) if '--bi_lstm2_units' in opts else DEFAULTS["bi_lstm2_units"]
    parameters["dense1_units"] = int(opts['--dense1_units']) if '--dense1_units' in opts else DEFAULTS["dense1_units"]
    parameters["dense2_units"] = int(opts['--dense2_units']) if '--dense2_units' in opts else DEFAULTS["dense2_units"]
    return parameters

def load_parameters(opts={}):
    parameters = {}
    input_settings = {str(x):setting for x,setting in enumerate(INPUT_SETTING_LIST)}
    input_settings.update({v:v for v in input_settings.values()})
    parameters["input_setting"] = input_settings[opts['--input_setting']] if '--input_setting' in opts else DEFAULTS["input_setting"]
    parameters["max_seq_len"] = int(opts['--max_seq_len']) if '--max_seq_len' in opts else DEFAULTS["max_seq_len"]
    parameters["model"] = opts['--model'] if '--model' in opts else DEFAULTS["model"]
    parameters["epochs"] = int(opts['--epochs']) if '--epochs' in opts else DEFAULTS["epochs"]
    parameters["verbose"] = True if '--verbose' in opts else DEFAULTS["verbose"]
    parameters["patience"] = int(opts['--patience']) if '--patience' in opts else DEFAULTS["patience"]
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
            print("Error: INPUT_SETTING paramer is set to an invalid setting.  Valid settings are (\"USE_SEQ_ONLY\",\"USE_SEQ_FOLD\",\"USE_SEQ_BPRNA\").")
            exit()
    f.close()
    return dataSet

def dropLongSequences(dataSet,parameters):
    newDataSet = []
    long_seq_count = 0
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
    for id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq,fold in dataSet:
        if parameters["use_embedding_layer"]:
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
    if parameters["dense1_units"] > 0:
        model = TimeDistributed(Dense(parameters["dense1_units"], activation="relu"))(model)
    if parameters["bi_lstm2_units"] > 0:
        model = Bidirectional(CuDNNLSTM(units=parameters["bi_lstm2_units"], return_sequences=True))(model)
    if parameters["dense2_units"] > 0:
        model = TimeDistributed(Dense(parameters["dense2_units"], activation="relu"))(model)
    if parameters["use_crf_layer"]:
        crf = CRF(parameters["num_tags"])  # CRF layer
        out = crf(model)  # output
        model = Model(input, out)
        model.compile(optimizer="rmsprop", loss=losses.crf_loss, metrics=[metrics.crf_accuracy])
    else:
        out = TimeDistributed(Dense(parameters["num_tags"], activation="softmax"))(model)
        model = Model(input, out)
        model.compile(optimizer="rmsprop", loss="categorical_crossentropy", metrics=["accuracy"])
    model.summary()
    return model

def pred2label(pred,parameters):
    idx2tag = {i: w for w, i in parameters["tag_idx"].items()}
    out = []
    for pred_i in pred:
        out_i = []
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

def create_history_fig(history,acc_outputFile,loss_outputFile,parameters):
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


def print_classification_report(validation_labels,pred_labels,outputFile=None):
    print("F1-score: {:.1%}".format(f1_score(validation_labels, pred_labels)))
    print(classification_report(validation_labels, pred_labels))
    if outputFile:
        f = open(outputFile,"w");
        f.write(classification_report(validation_labels, pred_labels))
        f.close()

def breakInputSeq(inputLabels):
    inputSequences = []
    inputFolds = []
    for i in range(0,len(inputLabels)):
        inputSeq = []
        inputFold = []
        for j in range(0,len(inputLabels[i])):
            if inputLabels[i][j] == 'PAD':
                inputSeq.append('P')
                inputFold.append('P')
            else:
                inputSeq.append(inputLabels[i][j][0])
                inputFold.append(inputLabels[i][j][1])
        inputSequences.append(inputSeq)
        inputFolds.append(inputFold)
    return inputSequences, inputFolds

def print_validation_output_file(validationSet,X_vl,validation_labels,pred_labels,outputFile,parameters):
    f = open(outputFile,"w")
    inputLabels = input2label(X_vl,parameters)
    inputSequences,inputFolds = breakInputSeq(inputLabels)
    for i in range(0,len(validationSet)):
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq,fold = validationSet[i]
        f.write("%s\n"%("\t".join([id,name,"".join(inputSequences[i]),"".join(inputFolds[i]),",".join(validation_labels[i]),",".join(pred_labels[i])])))
    f.close();

def print_classification_values_file(validationSet,pred,parameters,output_file="classification_values.txt"):
    f = open(output_file,"w")
    classification_values = get_classification_values(pred,parameters)
    tag_str = "\t".join(parameters["tags"])
    f.write("#id\tname\t%s\n"%(tag_str))
    for i in range(0,len(validationSet)):
        id,name,mir_id,prod5p,prod3p,drosha5p,dicer5p,dicer3p,drosha3p,hpStart,hpStop,seq,fold = validationSet[i]
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
