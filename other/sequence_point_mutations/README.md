# Scripts for RNA Sequence Point Mutations

Point mutations were performed by running the following scripts in order:

**createMutationSet&#46;pl** takes the test set (Metazoa\_testSet\_wFolds.txt) as input, and outputs several files with mutations for nucleotides surrounding each cutsite (Metazoa\_testSet\_wFolds.txt\_dicer3pMutations.txt, Metazoa\_testSet\_wFolds.txt\_dicer5pMutations.txt, Metazoa\_testSet\_wFolds.txt\_drosha3pMutations.txt, Metazoa\_testSet\_wFolds.txt\_drosha5pMutations.txt).
```sh
$ perl Scripts/createMutationSet.pl Metazoa_testSet_wFolds.txt
```
Mutation sets were evaluated using **deepMirCut/getCutSiteScores&#46;py**.  The script outputs a file (Metazoa\_testSet\_wFolds.txt\_cutsite_scores.txt) with decision values for each cutsite
```sh
$ python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt -m ~/deepMirCut/seqBPRNA.model
$ python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
$ python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_dicer5pMutations.txt -m ~/deepMirCut/seqBPRNA.model
$ python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha3pMutations.txt -m ~/deepMirCut/seqBPRNA.model
$ python3 ~/deepMirCut/getCutSiteScores.py Metazoa_testSet_wFolds.txt_drosha5pMutations.txt -m ~/deepMirCut/seqBPRNA.model
```

**generatePointMutationFile\_wReverse&#46;pl** takes as input the original test set file (Metazoa\_testSet\_wFolds.txt), the scores found for the test set (Metazoa\_testSet\_wFolds.txt\_cutsite\_scores.txt), the scores for the mutated test set (Metazoa\_testSet\_wFolds.txt\_dicer3pMutations.txt\_cutsite_scores.txt), and the cutsite to perform point mutations on (DC3).  The script outputs a point mutation file (DC3\_pointMutations.txt) with the change in decision values for each mutation.
```sh
$ perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_dicer3pMutations.txt_cutsite_scores.txt DC3
$ perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_dicer5pMutations.txt_cutsite_scores.txt DC5
$ perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_drosha3pMutations.txt_cutsite_scores.txt DR3
$ perl Scripts/generatePointMutationFile_wReverse.pl Metazoa_testSet_wFolds.txt Metazoa_testSet_wFolds.txt_cutsite_scores.txt Metazoa_testSet_wFolds.txt_drosha5pMutations.txt_cutsite_scores.txt DR5
```

**getTrainingFrequencies&#46;pl** takes the training set as input (Metazoa\_trainSet\_wSimilar\_mult\_wFolds.txt), and returns a set of files with the frequencies of occurrence for each nucleotide at positions surrounding around each cut site (DC3\_trainPerc.txt, DC5\_trainPerc.txt, DR3\_trainPerc.txt, DR5\_trainPerc.txt) as output.
```sh
$ perl Scripts/getTrainingFrequencies.pl Metazoa_trainSet_wSimilar_mult_wFolds.txt
```

**analyzePointMutations\_tScores&#46;pl** takes a point mutation file (DC3\_pointMutations.txt), and a file giving the training percentage of each nucleotide for each position (DC3\_trainPerc.txt), and a minimum training frequency (in this case 5 percent).  Then it returns a file containing the mean difference for each type of mutation (DC3\_pointMutations.txt\_stats\_wTScore\_minTrainFreq5.txt), and the statistical significance for that mutation (C3\_pointMutations.txt\_stats\_wTScore\_minTrainFreq5\_pvals.txt)
```sh
$ perl Scripts/analyzePointMutations_tScores.pl DC3_pointMutations.txt DC3_trainPerc.txt 5
$ perl Scripts/analyzePointMutations_tScores.pl DC5_pointMutations.txt DC5_trainPerc.txt 5
$ perl Scripts/analyzePointMutations_tScores.pl DR3_pointMutations.txt DR3_trainPerc.txt 5
$ perl Scripts/analyzePointMutations_tScores.pl DR5_pointMutations.txt DR5_trainPerc.txt 5
```

**mutationHeatMap&#46;pl** takes as input a file for the mean difference in decision value of the point mutation (DC3\_pointMutations.txt\_stats\_wTScore\_minTrainFreq5.txt), a value for the maximum difference (0.154779984977795), and a file for the significance of each point mutation (DC3\_pointMutations.txt\_stats\_wTScore\_minTrainFreq5\_pvals.txt), and then outputs a heatmap.
```sh
$ perl Scripts/mutationHeatMap.pl DC3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DC3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
$ perl Scripts/mutationHeatMap.pl DC5_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DC5_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
$ perl Scripts/mutationHeatMap.pl DR3_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DR3_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
$ perl Scripts/mutationHeatMap.pl DR5_pointMutations.txt_stats_wTScore_minTrainFreq5.txt 0.154779984977795 DR5_pointMutations.txt_stats_wTScore_minTrainFreq5_pvals.txt
```

**getCutsiteSeqs&#46;pl** takes the test set as input (Metazoa\_testSet\_wFolds.txt), and generates a set of fasta files for sequences surrounding each cutsite (Metazoa\_testSet\_wFolds.txt\_DR5.fa, Metazoa\_testSet\_wFolds.txt\_DC5.fa, Metazoa\_testSet\_wFolds.txt\_DC3.fa, Metazoa\_testSet\_wFolds.txt\_DR3.fa).
```sh
$ perl Scripts/getCutsiteSeqs.pl Metazoa_testSet_wFolds.txt
```

[weblogo](https://weblogo.berkeley.edu/) was use to create logos
```sh
$ weblogo -f Metazoa_testSet_wFolds.txt_DR5.fa -o Metazoa_testSet_wFolds.txt_DR5_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
$ weblogo -f Metazoa_testSet_wFolds.txt_DC5.fa -o Metazoa_testSet_wFolds.txt_DC5_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
$ weblogo -f Metazoa_testSet_wFolds.txt_DC3.fa -o Metazoa_testSet_wFolds.txt_DC3_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
$ weblogo -f Metazoa_testSet_wFolds.txt_DR3.fa -o Metazoa_testSet_wFolds.txt_DR3_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
```