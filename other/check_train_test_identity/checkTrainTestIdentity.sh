#!/bin/bash

#prepareSeqsForIdentityTests.pl takes as input the training set (Metazoa_trainSet.txt), the trainSet with similar miRs that have less than 80 percent identity with test and validation sets added back in (Metazoa_trainSet_wSimilar.txt), the training set with all similar mirs added back in (Metazoa_trainSet_wSimilar.txt.all), the validation set (Metazoa_validationSet.txt), the test set (Metazoa_testSet.txt), and the hairpin.fa file from miRBase.  The Script produces fastas of each sequence without the extended buffer regions. 
perl Scripts/prepareSeqsForIdentityTests.pl Metazoa_trainSet.txt Metazoa_trainSet_wSimilar.txt Metazoa_trainSet_wSimilar.txt.all Metazoa_validationSet.txt Metazoa_testSet.txt hairpin.fa

#compareMaxNeedleIdentity.pl takes as input two fasta files, and an output_prefix.  It compares every entry of the first fasta file with all entries of the second fatsta file and prints the maximum identity in a text file named using the output prefix.
perl Scripts/compareMaxNeedleIdentity.pl testSet_precursors.fa trainSet_wSimilar_precursors.fa testSet_vs_trainSetWSimilar
perl Scripts/compareMaxNeedleIdentity.pl testSet_precursors.fa trainSet_precursors.fa testSet_vs_trainSet
perl Scripts/compareMaxNeedleIdentity.pl validationSet_precursors.fa trainSet_wSimilar_precursors.fa validationSet_vs_trainSetWSimilar
perl Scripts/compareMaxNeedleIdentity.pl validationSet_precursors.fa trainSet_precursors.fa validationSet_vs_trainSet
perl Scripts/compareMaxNeedleIdentity.pl testSet_precursors.fa trainSet_wSimilar_precursors_all.fa testSet_vs_trainSetWSimilarWithoutFiltering
perl Scripts/compareMaxNeedleIdentity.pl validationSet_precursors.fa trainSet_wSimilar_precursors_all.fa validationSet_vs_trainSetWSimilarIncludingWithoutFiltering
