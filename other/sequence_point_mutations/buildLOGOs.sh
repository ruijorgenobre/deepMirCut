#!/bin/bash

#getCutsiteSeqs.pl takes the test set as input (Metazoa_testSet_wFolds.txt), and generates a set of fasta files for sequences surrounding each cutsite (Metazoa_testSet_wFolds.txt_DR5.fa, Metazoa_testSet_wFolds.txt_DC5.fa, Metazoa_testSet_wFolds.txt_DC3.fa, Metazoa_testSet_wFolds.txt_DR3.fa).
perl Scripts/getCutsiteSeqs.pl Metazoa_testSet_wFolds.txt

#weblogo was use to create logos
weblogo -f Metazoa_testSet_wFolds.txt_DR5.fa -o Metazoa_testSet_wFolds.txt_DR5_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
weblogo -f Metazoa_testSet_wFolds.txt_DC5.fa -o Metazoa_testSet_wFolds.txt_DC5_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
weblogo -f Metazoa_testSet_wFolds.txt_DC3.fa -o Metazoa_testSet_wFolds.txt_DC3_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600
weblogo -f Metazoa_testSet_wFolds.txt_DR3.fa -o Metazoa_testSet_wFolds.txt_DR3_logo.eps -F eps -c classic -S 0.2 -Y Yes --ticmarks 0.1 --resolution 600

