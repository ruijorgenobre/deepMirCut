# Scripts used in Best Model Testing

The Best Model was analyzed by running scripts in the following order:

Test set was evaluated with **deepMirCut/testModel&#46;py** using the ensemble with the -d flag set so that it would output decision values. (You may need to change the path depending on where you have deepMirCut/testModel.py installed.)
```sh
$ python3 ~/deepMirCut/testModel.py Metazoa_testSet_wFolds.txt -m seqBPRNA.model -o seqBPRNA -d --input_setting 2
```

**Scripts/graph_classificationDV_boxPlot&#46;py** takes the test set (Metazoa_testSet_wFolds.txt), the decision values returned by testModel&#46;py (seqBPRNA_classification_DVs.txt), and an output prefix (seqBPRNA) as input and outputs several box plots showing decision values for predicted cutsites
```sh
$ python3 Scripts/graph_classificationDV_boxPlot.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA
```

**Scripts/graph_deltas_histogram&#46;py** takes the test set (Metazoa_testSet_wFolds.txt), the decision values returned by testModel.py (seqBPRNA_classification_DVs.txt), and an output prefix (seqBPRNA) as input and outputs a histogram for distances of predicted cutsites relative to actual cutsites.
```sh
$ python3 Scripts/graph_deltas_histogram.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA
```
