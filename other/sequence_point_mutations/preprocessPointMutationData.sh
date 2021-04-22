#!/bin/bash
#createMutationSet.pl creates takes the test set (Metazoa_testSet_wFolds.txt) as input, and outputs several files with mutations for nucleotides surrounding each cutsite (Metazoa_testSet_wFolds.txt_dicer3pMutations.txt, Metazoa_testSet_wFolds.txt_dicer5pMutations.txt, Metazoa_testSet_wFolds.txt_drosha3pMutations.txt, Metazoa_testSet_wFolds.txt_drosha5pMutations.txt).
perl Scripts/createMutationSet.pl Metazoa_testSet_wFolds.txt

#Mutation sets were evaluated using deepMirCut/getCutSiteScores.py.  The script outputs a file (Metazoa_testSet_wFolds.txt_cutsite_scores.txt) with decision values for each cutsite
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer5pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha5pMutations.txt -m ~/deepMirCut/seqBPRNA.model

