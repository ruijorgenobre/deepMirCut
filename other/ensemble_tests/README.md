
# Scripts used in Ensemble testing

Ensembles were analyzed by running scripts in the following order:

**deepMirCut/testEnsembles\_average\_models\_pmf&#46;py** tests pmf for a set of ensembles increasing in size from 1 to 20.  It outputs a file displaying the pmf for each ensemble (seqBPRNA\_unweighted\_pmf.txt) and the pmf for the individual unweighted models (seqBPRNA\_unweighted\_pmf\_single\_models.txt)
```sh
$ python3 ~/deepMirCut/testEnsembles_average_models_pmf.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2
```

**deepMirCut/testEnsembles\_average\_models\_pse&#46;py** tests pse for a set of ensembles increasing in size from 1 to 20.  It outputs a file displaying the pse for the ensemble (seqBPRNA\_unweighted\_pse.txt) and the pse for the individual models (seqBPRNA\_unweighted\_pse\_single\_models.txt)
```sh
$ python3 ~/Scripts/deepMirCut/testEnsembles_average_models_pse.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2
```

**Scripts/graph\_unweighted\_ensemble\_pmf&#46;py** graphs pmf for ensembles increasing in size from 1 to 20
```sh
$ python3 Scripts/graph_unweighted_ensemble_pmf.py seqBPRNA_unweighted_pmf.txt seqBPRNA_unweighted_pmf_single_models.txt
```

**Scripts/graph\_unweighted\_ensemble\_pse&#46;py** graphs pse for ensembles increasing in size from 1 to 20
```sh
python3 Scripts/graph_unweighted_ensemble_pse.py seqBPRNA_unweighted_pse.txt seqBPRNA_unweighted_pse_single_models.txt
```

**deepMirCut/prepareEnsemble\_average\_models&#46;py** is used to prepare an ensemble of a size indicated by the -n option (in this case 12)
```sh
$ python3 ~/deepMirCut/prepareEnsemble_average_models.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2 -n 12
```

**The ensemble was tested using deepMirCut/testModel&#46;py**
```sh
$ python3 ~/deepMirCut/testModel.py Metazoa_testSet_wFolds.txt -L seqBPRNA_ensemble_list.txt -o seqBPRNA -d --input_setting 2
```

**Scripts/graph\_deltas\_histogram&#46;py** takes the test set (Metazoa\_testSet\_wFolds.txt), the decision values returned by the ensemble (seqBPRNA\_classification\_DVs.txt), and an output prefix (seqBPRNA) as input and uses this data to plot a histogram of predicted cutsites relative to annotated cutsites. 
```sh
$ python3 Scripts/graph_deltas_histogram.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA
```
