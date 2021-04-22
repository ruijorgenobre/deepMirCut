#!/bin/bash

#generatePointMutationFile_wReverse.pl takes as input the original test set file (Metazoa_testSet_wFolds.txt), the scores found for the test set (Metazoa_testSet_wFolds.txt_cutsite_scores.txt), the scores for the mutated test set (Metazoa_testSet_wFolds.txt_dicer3pMutations.txt_cutsite_scores.txt), and the cutsite to compare (DC3).  The script outputs a point mutation file (DC3_pointMutations.txt) with the change in decision values for each mutation. 
perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_dicer3pMutations.txt_cutsite_scores.txt DC3
perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_dicer5pMutations.txt_cutsite_scores.txt DC5
perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_drosha3pMutations.txt_cutsite_scores.txt DR3
perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_drosha5pMutations.txt_cutsite_scores.txt DR5

#getTrainingFrequencies.pl takes the training set as input (Metazoa_trainSet_wSimilar_mult_wFolds.txt), and returns a set of files with the frequencies of occurance for each nucleotide at positions surrounding around each cut site (DC3_trainPerc.txt, DC5_trainPerc.txt, DR3_trainPerc.txt, DR5_trainPerc.txt).
perl Scripts/getTrainingFrequencies.pl Metazoa_trainSet_wSimilar_mult_wFolds.txt

#analyzePointMutations_tScores.pl takes a point mutation file (DC3_pointMutations.txt), and a file giving the training percentage of each nucleodite for each position (DC3_trainPerc.txt), and a minimum training frequency (in this case 5 percent).  Then it returns a file containing the mean difference for each type of mutation (DC3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt), and the statistical significance for that mutation (C3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt)
perl Scripts/analyzePointMutations_tScores.pl DC3_pointMutations.txt DC3_trainPerc.txt 5
perl Scripts/analyzePointMutations_tScores.pl DC5_pointMutations.txt DC5_trainPerc.txt 5
perl Scripts/analyzePointMutations_tScores.pl DR3_pointMutations.txt DR3_trainPerc.txt 5
perl Scripts/analyzePointMutations_tScores.pl DR5_pointMutations.txt DR5_trainPerc.txt 5

#mutationHeatMap.pl takes as input a file for the mean difference in decision value of the point mutation (DC3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt), a value for the maximum difference (0.154779984977795), and a file for the significance of each point mutation (DC3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt), and then outputs a heatmap.
perl Scripts/mutationHeatMap.pl DC3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DC3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
perl Scripts/mutationHeatMap.pl DC5_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DC5_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
perl Scripts/mutationHeatMap.pl DR3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DR3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
perl Scripts/mutationHeatMap.pl DR5_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DR5_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt