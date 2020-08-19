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

opts = {}
dmc_parameters = {}
X_tr,y_tr,X_vl,y_vl = 0,0,0,0


space_keys = ['batch_size','embedding_dropout','embedding_layer_output','bi_lstm1_units','dense1_units','bi_lstm2_units','dense2_units','optimizer',
              'rms_learning_rate','rms_rho','rms_epsilon','adam_learning_rate','adam_beta_1','adam_beta_2','adam_epsilon','adam_amsgrad','nadam_learning_rate','nadam_beta_1','nadam_beta_2','nadam_epsilon','sgd_learning_rate',
              'sgd_momentum','sgd_nesterov']

search_space = {'batch_size': hp.choice('batch_size',[16,32,64,128]),
                'embedding_dropout' : hp.uniform('embedding_dropout',0, 0.5),
                'embedding_layer_output' : hp.choice('embedding_layer_output',[32,64,96,128,160]),
                'bi_lstm1_units' : hp.choice('bi_lstm1_units',[0,32,64,96,128,160]),
                'dense1_units' : hp.choice('dense1_units',[0,32,64,96,128,160]),
                'bi_lstm2_units' : hp.choice('bi_lstm2_units',[0,32,64,96,128,160]),
                'dense2_units' : hp.choice('dense2_units',[0,32,64,96,128,160]),
                'optimizer' : hp.choice('optimizer',[
                    {'type':'rmsprop',
                     'rms_learning_rate' : hp.uniform('rms_learning_rate',1e-5,0.1),
                     'rms_rho' : hp.uniform('rms_rho',0.6,1),
                     'rms_epsilon' : hp.uniform('rms_epsilon', 1e-10, 1e-4)
                    },
                    {'type':'adam',
                     'adam_learning_rate' : hp.uniform('adam_learning_rate',1e-5,0.1),
                     'adam_beta_1' : hp.uniform('adam_beta_1',0.7,1),
                     'adam_beta_2' : hp.uniform('adam_beta_2',0.7,1),
                     'adam_epsilon' : hp.uniform('adam_epsilon',1e-10, 1e-4),
                     'adam_amsgrad' : hp.choice('adam_amsgrad',[True,False]),
                    },
                    {'type':'nadam',
                     'nadam_learning_rate' : hp.uniform('nadam_learning_rate',1e-5,0.1),
                     'nadam_beta_1' : hp.uniform('nadam_beta_1',0.7,1),
                     'nadam_beta_2' : hp.uniform('nadam_beta_2',0.7,1),
                     'nadam_epsilon' : hp.uniform('nadam_epsilon',1e-10, 1e-4)
                    },
                    {'type':'sgd',
                     'sgd_learning_rate' : hp.uniform('sgd_learning_rate',1e-5,0.1),
                     'sgd_momentum' : hp.uniform('sgd_momentum',0.5,1),
                     'sgd_nesterov' : hp.choice('sgd_nesterov',[True,False]),
                    }])
                }



def print_to_pickle_file(trials,pickle_file="tuning.pickle"):
    f = open(pickle_file, 'wb')
    pickle.dump(trials,f)

def get_trial_vals(search_space,trial):
    #returns values of space_keys in the same order
    trial_kv = space_eval(search_space,{k:v[0] for k,v in trial['misc']['vals'].items() if len(v) > 0})
    trial_vals = []
    for k in space_keys:
        if k in trial_kv and k == 'optimizer':
            trial_vals.append(trial_kv[k]['type'])
        elif k in trial_kv['optimizer'].keys():
            trial_vals.append(trial_kv['optimizer'][k])
        elif k in trial_kv:
            trial_vals.append(trial_kv[k])
        else:
            trial_vals.append("na")
    return trial_vals


def print_trials(search_space,trials,output_file="tuning_trials.txt"):
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
        trials_fd.write("%d\t%s\t%g\t%g\t%g\t%g\t%d\t%s\n" % (tid,"\t".join([str(x) for x in trial_vals]),optimizer_loss,trial_loss,trial_accuracy,trial_proximity,trial_fScore,trial_status))
    trials_fd.close()


def print_best_trial(search_space,trials,output_file="tuning_best_trial.txt"):
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

    optimizer = 0
    if params['optimizer']['type'] == 'rmsprop':
        optimizer = optimizers.rmsprop(lr=params['optimizer']['rms_learning_rate'],
                                       rho=params['optimizer']['rms_rho'],
                                       epsilon=params['optimizer']['rms_epsilon'])
    elif params['optimizer']['type'] == 'adam':
        optimizer = optimizers.adam(lr=params['optimizer']['adam_learning_rate'],
                                    beta_1=params['optimizer']['adam_beta_1'], 
                                    beta_2=params['optimizer']['adam_beta_2'],
                                    epsilon=params['optimizer']['adam_epsilon'],
                                    amsgrad=params['optimizer']['adam_amsgrad'])
    elif params['optimizer']['type'] == 'nadam':
        optimizer = optimizers.nadam(lr=params['optimizer']['nadam_learning_rate'],
                                    beta_1=params['optimizer']['nadam_beta_1'], 
                                    beta_2=params['optimizer']['nadam_beta_2'],
                                    epsilon=params['optimizer']['nadam_epsilon'])
    elif params['optimizer']['type'] == 'sgd':
        optimizer = optimizers.sgd(lr=params['optimizer']['sgd_learning_rate'],
                                   momentum=params['optimizer']['sgd_momentum'],
                                   nesterov=params['optimizer']['sgd_nesterov'])
    
    if dmc_parameters["use_embedding_layer"]:
        input = Input(shape=(dmc_parameters["max_seq_len"],))
        model = Embedding(input_dim=dmc_parameters["one_hot_vector_len"], output_dim=params['embedding_layer_output'],input_length=dmc_parameters["max_seq_len"])(input)
        model = Dropout(rate = params['embedding_dropout'])(model)
    else:
        input = Input(shape=(dmc_parameters["max_seq_len"],dmc_parameters["one_hot_vector_len"]))
        model = input
    if params['bi_lstm1_units'] > 0:
        model = Bidirectional(CuDNNLSTM(units=params['bi_lstm1_units'], return_sequences=True))(model)
    if params['dense1_units'] > 0:
        model = TimeDistributed(Dense(params['dense1_units'], activation="relu"))(model)
    if params['bi_lstm2_units'] > 0:
        model = Bidirectional(CuDNNLSTM(units=params['bi_lstm2_units'], return_sequences=True))(model)
    if params['dense2_units'] > 0:
        model = TimeDistributed(Dense(params['dense2_units'], activation="relu"))(model)
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
    history = model.fit(X_tr, np.array(y_tr),batch_size=params['batch_size'],epochs=dmc_parameters["epochs"],validation_data=(X_vl,np.array(y_vl)),verbose=False,shuffle=True,callbacks=[es])
    loss, acc, prox = model.evaluate(x=X_vl, y=np.array(y_vl),batch_size=params['batch_size'], verbose=False)
    validation_labels = deepMirCut.pred2label(y_vl,dmc_parameters)
    validation_pred = model.predict(X_vl, verbose=False)
    pred_labels = deepMirCut.pred2label(validation_pred,dmc_parameters)
    fScore = f1_score(validation_labels, pred_labels)
    return loss, acc, prox, fScore

def f(params):
    loss, acc, prox,fScore = hyperopt_train_test(params)
    return {'loss': -fScore, 'real_loss': loss, 'accuracy':acc, 'proximity': prox, 'fScore': fScore , 'status': STATUS_OK}

def run_trials(pickle_file="tuning.pickle",trials_step=25):

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

    print_to_pickle_file(trials,pickle_file="tuning.pickle")
    print_trials(search_space,trials,output_file="tuning_trials.txt")
    print_best_trial(search_space,trials,output_file="tuning_best_trial.txt")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("USAGE: %s <train file> <validation file> <evaluations>\n"%(sys.argv[0]))
        exit()
    
    new_evals = int(sys.argv[3])
    trial_step = 20

    opts = {}
    opts['--verbose'] = True
    opts['--epochs'] = 20
    opts['--patience'] = 3
    opts['--use_embedding_layer'] = "True"
    opts['--use_crf_layer'] = "False"
    dmc_parameters = deepMirCut.load_parameters(opts)
    dmc_parameters['train_file'] = sys.argv[1]
    dmc_parameters['validation_file'] = sys.argv[2]

    trainSet = deepMirCut.readDataset(dmc_parameters["train_file"],dmc_parameters)
    new_trainSet = deepMirCut.dropLongSequences(trainSet,dmc_parameters)
    X_tr,y_tr = deepMirCut.prepareData(new_trainSet,dmc_parameters)
    validationSet = deepMirCut.readDataset(dmc_parameters["validation_file"],dmc_parameters)
    new_validationSet = deepMirCut.dropLongSequences(validationSet,dmc_parameters)
    X_vl,y_vl = deepMirCut.prepareData(new_validationSet,dmc_parameters)

    while (trial_step < new_evals):
        run_trials(pickle_file="tuning.pickle",trials_step=trial_step)
        new_evals -= trial_step

    run_trials(pickle_file="tuning.pickle",trials_step=new_evals)    
