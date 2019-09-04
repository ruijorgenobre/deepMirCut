import deepMirCut
import sys, os
import numpy as np
import pandas as pd
from keras.models import Model, Input
from keras.layers import CuDNNLSTM, LSTM, Embedding, Dense, TimeDistributed, Dropout, Bidirectional
from keras_contrib.layers import CRF
from keras_contrib import losses
from keras_contrib import metrics
from keras_contrib.utils import save_load_utils
from keras.callbacks import EarlyStopping
from seqeval.metrics import precision_score, recall_score, f1_score, classification_report
from hyperopt import Trials, STATUS_OK, tpe
from hyperas import optim
from hyperas.distributions import choice, uniform
import getopt

def data():
    import deepMirCut
    opts = {}
    opts['--verbose'] = True
    opts['--epochs'] = 20
    opts['--patience'] = 3
    opts['--use_embedding_layer'] = "True"
    opts['--use_crf_layer'] = "False"
    parameters = deepMirCut.load_parameters(opts)
    parameters['train_file'] = "Metazoa_trainSet_wFolds.txt"
    parameters['validation_file'] = "Metazoa_validateSet_wFolds.txt"    
    
    trainSet = deepMirCut.readDataset(parameters["train_file"],parameters)
    new_trainSet = deepMirCut.dropLongSequences(trainSet,parameters)
    X_tr,y_tr = deepMirCut.prepareData(new_trainSet,parameters)
    validationSet = deepMirCut.readDataset(parameters["validation_file"],parameters)
    new_validationSet = deepMirCut.dropLongSequences(validationSet,parameters)
    X_vl,y_vl = deepMirCut.prepareData(new_validationSet,parameters)
    return X_tr,y_tr,X_vl,y_vl

def model(X_tr,y_tr,X_vl,y_vl):
    batch_size = {{choice([16,32,64,128])}}
    embedding_dropout = {{uniform(0, 1)}}
    embedding_layer_output = {{choice([32,64,96,128,160])}}
    bi_lstm1_units = {{choice([0,32,64,96,128,160])}}
    dense1_units = {{choice([0,32,64,96,128,160])}}
    bi_lstm2_units = {{choice([0,32,64,96,128,160])}}
    dense2_units = {{choice([0,32,64,96,128,160])}}
    optimizer = {{choice(['rmsprop', 'adam', 'sgd'])}}
    
    if parameters["use_embedding_layer"]:
        input = Input(shape=(parameters["max_seq_len"],))
        model = Embedding(input_dim=parameters["one_hot_vector_len"], output_dim=embedding_layer_output,input_length=parameters["max_seq_len"])(input)
        model = Dropout(rate = embedding_dropout)(model)
    else:
        input = Input(shape=(parameters["max_seq_len"],parameters["one_hot_vector_len"]))
        model = input
    if bi_lstm1_units > 0:
        model = Bidirectional(CuDNNLSTM(units=bi_lstm1_units, return_sequences=True))(model)
    if dense1_units > 0:
        model = TimeDistributed(Dense(dense1_units, activation="relu"))(model)
    if bi_lstm2_units > 0:
        model = Bidirectional(CuDNNLSTM(units=bi_lstm2_units, return_sequences=True))(model)
    if dense2_units > 0:
        model = TimeDistributed(Dense(dense2_units, activation="relu"))(model)
    if parameters["use_crf_layer"]:
        crf = CRF(parameters["num_tags"])  # CRF layer
        out = crf(model)  # output
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss=losses.crf_loss, metrics=[metrics.crf_accuracy])
    else:
        out = TimeDistributed(Dense(parameters["num_tags"], activation="softmax"))(model)
        model = Model(input, out)
        model.compile(optimizer=optimizer, loss="categorical_crossentropy", metrics=["accuracy"])
    model.summary()
    es = EarlyStopping(monitor='val_loss', min_delta=0, patience=parameters["patience"], verbose=False, mode='min', restore_best_weights=True)
    history = model.fit(X_tr, np.array(y_tr),batch_size=batch_size,epochs=parameters["epochs"],validation_data=(X_vl,np.array(y_vl)),verbose=False,shuffle=True,callbacks=[es])
    loss, acc = model.evaluate(x=X_vl, y=np.array(y_vl),batch_size=batch_size, verbose=False)
    return {'loss': loss, 'status': STATUS_OK, 'model': model}

if __name__ == "__main__":
    max_evals = 100
    
    best_run, best_model = optim.minimize(model=model,
                                          data=data,
                                          algo=tpe.suggest,
                                          max_evals=max_evals,
                                          trials=Trials())
    X_train, Y_train, X_test, Y_test = data()
    print("Evalutation of best performing model:")
    print(best_model.evaluate(X_test, np.array(Y_test)))
    print(best_run)
    f = open("best_run.txt","w")
    f.write("%s\n"%(best_run))
    best_model.save('best_model.h5')
