#!/bin/bash

#createMutationSet_wBPRNA.pl creates takes the test set (Metazoa_testSet_wFolds.txt) as input, and outputs several files with mutations for bpRNA surrounding each cutsite (Metazoa_testSet_wFolds.txt_dicer3pMutations.txt, Metazoa_testSet_wFolds.txt_dicer5pMutations.txt, Metazoa_testSet_wFolds.txt_drosha3pMutations.txt, Metazoa_testSet_wFolds.txt_drosha5pMutations.txt).
perl Scripts/createMutationSet_wBPRNA.pl Metazoa_testSet_wFolds.txt

#Mutation sets were evaluated using deepMirCut/getCutSiteScores.py.  The script outputs a file (Metazoa_testSet_wFolds.txt_cutsite_scores.txt) with decision values for each cutsite
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer5pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha5pMutations.txt -m ~/deepMirCut/seqBPRNA.model


#set of tests to check that mutations were done properly.  (mature.fa comes from mirbase)
#perl Scripts/checkMutationSet.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_drosha5pMutations.txt mature.fa
#perl Scripts/checkMutationSet.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_dicer5pMutations.txt mature.fa
#perl Scripts/checkMutationSet.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_dicer3pMutations.txt mature.fa
#perl Scripts/checkMutationSet.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_drosha3pMutations.txt mature.fa
