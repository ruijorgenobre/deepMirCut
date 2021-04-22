#!/bin/bash

#deepMirCut/testEnsembles_average_models_pmf.py tests pmf for a set of ensembles increasing in size from 1 to 20.  It outputs a file displaying the pmf for the ensemble (seqBPRNA_unweighted_pmf.txt) and the pmf for the individual unweighted models (seqBPRNA_unweighted_pmf_single_models.txt)
python3 ~/deepMirCut/testEnsembles_average_models_pmf.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2

#deepMirCut/testEnsembles_average_models_pse.py tests pse for a set of ensembles increasing in size from 1 to 20.  It outputs a file displaying the pse for the ensemble (seqBPRNA_unweighted_pse.txt) and the pse for the individual models (seqBPRNA_unweighted_pse_single_models.txt)
python3 ~/Scripts/deepMirCut/testEnsembles_average_models_pse.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2

#Scripts/graph_unweighted_ensemble_pmf.py graphs pmf for ensembles increasing in size from 1 to 20
python3 Scripts/graph_unweighted_ensemble_pmf.py seqBPRNA_unweighted_pmf.txt seqBPRNA_unweighted_pmf_single_models.txt

#Scripts/graph_unweighted_ensemble_pse.py graphs pse for ensembles increasing in size from 1 to 20
python3 Scripts/graph_unweighted_ensemble_pse.py seqBPRNA_unweighted_pse.txt seqBPRNA_unweighted_pse_single_models.txt

#deepMirCut/prepareEnsemble_average_models.py is used to prepare an ensemble of size 12
python3 ~/deepMirCut/prepareEnsemble_average_models.py Metazoa_validationSet_mult_wFolds.txt -L seqBPRNA_modelList.txt -o seqBPRNA -d --input_setting 2 -n 12

#the ensemble was tested using deepMirCut/testModel.py
python3 ~/deepMirCut/testModel.py Metazoa_testSet_wFolds.txt -L seqBPRNA_ensemble_list.txt -o seqBPRNA -d --input_setting 2

#Scripts/graph_deltas_histogram.py takes the test set (Metazoa_testSet_wFolds.txt), the decision values returned by the ensemble (seqBPRNA_classification_DVs.txt), and an output prefix (seqBPRNA) as input and uses this data to plot a histogram of offsets from the cutsite. 
python3 Scripts/graph_deltas_histogram.py Metazoa_testSet_wFolds.txt seqBPRNA_classification_DVs.txt seqBPRNA
