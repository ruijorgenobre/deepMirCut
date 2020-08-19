import deepMirCut
from deepMirCut_metrics import avg_proximity_metric
import sys, os
import numpy as np
import pandas as pd
from keras import optimizers
from keras.models import load_model, model_from_json
from keras.models import Model, Input
from keras.layers import CuDNNLSTM, LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional
from keras_contrib.layers import CRF
from keras_contrib import losses
from keras_contrib import metrics
from keras_contrib.utils import save_load_utils
from keras.callbacks import EarlyStopping
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK, space_eval
import getopt
import pickle

os.environ["CUDA_VISIBLE_DEVICES"]="0";

opts = {}
dmc_parameters = {}
X_tr,y_tr,X_vl,y_vl = 0,0,0,0

space_keys = ['embedding_dropout','embedding_layer_output','bi_lstm1_units','bi_lstm2_units','learning_rate','epsilon_exp']

search_space = {'embedding_dropout' : hp.uniform('embedding_dropout',0, 0.5),
                'embedding_layer_output' : hp.choice('embedding_layer_output',[32,64,96,128,160]),
                'bi_lstm1_units' : hp.choice('bi_lstm1_units',[32,64,96,128,160]),
                'bi_lstm2_units' : hp.choice('bi_lstm2_units',[0,32,64,96,128,160]),
                'learning_rate' : hp.uniform('learning_rate',1e-5,0.1),
                'epsilon_exp' : hp.uniform('epsilon_exp',-10,-4)}


def print_to_pickle_file(trials,pickle_file="tuning_adam_seqOnly.pickle"):
    f = open(pickle_file, 'wb')
    pickle.dump(trials,f)

def get_trial_vals(search_space,trial):
    #returns values of space_keys in the same order
    trial_kv = space_eval(search_space,{k:v[0] for k,v in trial['misc']['vals'].items() if len(v) > 0})
    trial_vals = []
    for k in space_keys:
        if k in trial_kv:
            trial_vals.append(trial_kv[k])
        else:
            trial_vals.append("na")
    return trial_vals

def print_trials(search_space,trials,output_file="tuning_trials_seqOnly.txt"):
    trials_fd = open(output_file,"w");
    trials_fd.write("tid\t%s\toptimizer_loss\tloss\taccuracy\tproximity\tf1_score\ttrial_status\n" % ("\t".join(space_keys)))
    for trial in trials.trials:
        tid = trial['tid']
        trial_vals = get_trial_vals(search_space,trial)
        optimizer_loss = trial['result']['loss']
        trial_loss = trial['result']['real_loss']
        trial_accuracy = trial['result']['accuracy']
        trial_proximity = trial['result']['proximity']
        trial_status = trial['result']['status']
        trial_fScore = trial['result']['fScore']
        trials_fd.write("%d\t%s\t%g\t%g\t%g\t%g\t%g\t%s\n" % (tid,"\t".join([str(x) for x in trial_vals]),optimizer_loss,trial_loss,trial_accuracy,trial_proximity,trial_fScore,trial_status))
    trials_fd.close()


def print_best_trial(search_space,trials,output_file="tuning_best_trial_seqOnly.txt"):
    best_trial = trials.best_trial
    best_trial_tid = best_trial['tid']
    best_trial_vals = get_trial_vals(search_space,best_trial)
    best_optimizer_loss = best_trial['result']['loss']
    best_trial_loss = best_trial['result']['real_loss']
    best_trial_accuracy = best_trial['result']['accuracy']
    best_trial_proximity = best_trial['result']['proximity']
    best_trial_fScore = best_trial['result']['fScore']
    best_trial_status = best_trial['result']['status']
    trials_fd = open(output_file,"w");
    trials_fd.write("tid = %d\n" % best_trial_tid)
    for i in range(0,len(space_keys)):
        trials_fd.write("%s\t%s\n"%(space_keys[i],str(best_trial_vals[i])))
    trials_fd.write("optimizer_loss = %g\n"% (best_optimizer_loss))
    trials_fd.write("loss = %g\n" % (best_trial_loss))
    trials_fd.write("accuracy = %g\n" % (best_trial_accuracy))
    trials_fd.write("proximity = %g\n" % (best_trial_proximity))
    trials_fd.write("fScore = %g\n" % (best_trial_fScore))
    trials_fd.write("status = %s\n" % (best_trial_status))    
    trials_fd.close()

def hyperopt_train_test(params):

    epsilon = 10 ** params['epsilon_exp']
    optimizer = optimizers.adam(lr=params['learning_rate'],epsilon=epsilon)
    
    if dmc_parameters["use_embedding_layer"]:
        input = Input(shape=(dmc_parameters["max_seq_len"],))
        model = Embedding(input_dim=dmc_parameters["one_hot_vector_len"], output_dim=params['embedding_layer_output'],input_length=dmc_parameters["max_seq_len"])(input)
        model = Dropout(rate = params['embedding_dropout'])(model)
    else:
        input = Input(shape=(dmc_parameters["max_seq_len"],dmc_parameters["one_hot_vector_len"]))
        model = input
    if params['bi_lstm1_units'] > 0:
        model = Bidirectional(CuDNNLSTM(units=params['bi_lstm1_units'], return_sequences=True))(model)
    if params['bi_lstm2_units'] > 0:
        model = Bidirectional(CuDNNLSTM(units=params['bi_lstm2_units'], return_sequences=True))(model)
    if dmc_parameters["use_crf_layer"]:
        crf = CRF(dmc_parameters["num_tags"])  # CRF layer
        out = crf(model)  # output
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss=losses.crf_loss, metrics=[metrics.crf_accuracy,avg_proximity_metric()])
    else:
        out = TimeDistributed(Dense(dmc_parameters["num_tags"], activation="softmax"))(model)
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss="categorical_crossentropy", metrics=["accuracy",avg_proximity_metric()])
    model.summary()
    es = EarlyStopping(monitor='val_loss', min_delta=0, patience=dmc_parameters["patience"], verbose=False, mode='min', restore_best_weights=True)
    history = model.fit(X_tr, np.array(y_tr),batch_size=dmc_parameters['batch_size'],epochs=dmc_parameters["epochs"],validation_data=(X_vl,np.array(y_vl)),verbose=False,shuffle=True,callbacks=[es])
    loss, acc, prox = model.evaluate(x=X_vl, y=np.array(y_vl),batch_size=dmc_parameters['batch_size'], verbose=False)
    validation_labels = deepMirCut.pred2label(y_vl,dmc_parameters)
    validation_pred = model.predict(X_vl, verbose=False)
    pred_labels = deepMirCut.pred2label(validation_pred,dmc_parameters)
    fScore = f1_score(validation_labels, pred_labels)
    return loss, acc, prox, fScore

def f(params):
    loss, acc, prox,fScore = hyperopt_train_test(params)
    return {'loss': -fScore, 'real_loss': loss, 'accuracy':acc, 'proximity': prox, 'fScore': fScore , 'status': STATUS_OK}

def run_trials(pickle_file="tuning_adam_seqOnly.pickle",trials_step=25):

    max_evals = trials_step
    
    try:
        trials = pickle.load(open(pickle_file, "rb"))
        print("Loading trials from %s"%(pickle_file))
        max_evals = len(trials.trials) + trials_step
        print("Running trials %d to %d"%(len(trials.trials),max_evals))
    except:
        trials = Trials()
        print("Running trials 1 to %d"%(max_evals))
    
    best = fmin(f, search_space, algo=tpe.suggest, max_evals=max_evals, trials=trials)
    print(space_eval(search_space,best))

    print_to_pickle_file(trials,pickle_file=pickle_file)
    print_trials(search_space,trials,output_file="tuning_trials_adam_seqOnly.txt")
    print_best_trial(search_space,trials,output_file="tuning_best_trial_adam_seqOnly.txt")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("USAGE: %s <train file> <validation file> <evaluations>\n"%(sys.argv[0]))
        exit()
    
    new_evals = int(sys.argv[3])
    trial_step = 10000

    opts = {}
    opts['--verbose'] = True
    opts['--epochs'] = 20
    opts['--patience'] = 3
    opts['--use_embedding_layer'] = "True"
    opts['--use_crf_layer'] = "False"
    opts['--input_setting'] = '0' #USE_SEQ_ONLY
    dmc_parameters = deepMirCut.load_parameters(opts)
    dmc_parameters['train_file'] = sys.argv[1]
    dmc_parameters['validation_file'] = sys.argv[2]
    dmc_parameters['batch_size'] = 64

    print("input_setting=",dmc_parameters['input_setting'])

    trainSet = deepMirCut.readDataset(dmc_parameters["train_file"],dmc_parameters)
    new_trainSet = deepMirCut.dropLongSequences(trainSet,dmc_parameters)
    X_tr,y_tr = deepMirCut.prepareData(new_trainSet,dmc_parameters)
    validationSet = deepMirCut.readDataset(dmc_parameters["validation_file"],dmc_parameters)
    new_validationSet = deepMirCut.dropLongSequences(validationSet,dmc_parameters)
    X_vl,y_vl = deepMirCut.prepareData(new_validationSet,dmc_parameters)

    while (trial_step < new_evals):
        run_trials(pickle_file="tuning_adam_seqOnly.pickle",trials_step=trial_step)
        new_evals -= trial_step

    run_trials(pickle_file="tuning_adam_seqOnly.pickle",trials_step=new_evals)    
