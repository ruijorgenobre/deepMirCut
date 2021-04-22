#!/bin/bash

#Test set was evaluated with deepMirCut/testModel.py using the ensemble with the -d flag set to output decision values. (You may need to change the path depending on where you have deepMirCut/testModel.py installed.)
python3 ~/deepMirCut/testModel.py Metazoa_testSet_wFolds.txt -m seqBPRNA.model -o seqBPRNA -d --input_setting 2

#Scripts/graph_classificationDV_boxPlot.py takes the test set (Metazoa_testSet_wFolds.txt), the decision values returned by testModel.py (seqBPRNA_classification_DVs.txt), and an output prefix (seqBPRNA) as input and output several box plots showing decision values for predicted cutsites
python3 Scripts/graph_classificationDV_boxPlot.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA

#Scripts/graph_deltas_histogram.py the test set (Metazoa_testSet_wFolds.txt), the decision values returned by testModel.py (seqBPRNA_classification_DVs.txt), and an output prefix (seqBPRNA) as input and outputs a histogram for distances of predicted cutsites relative to actual cutsites.
python3 Scripts/graph_deltas_histogram.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA
