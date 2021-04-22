#!/bin/bash

#run the testModel.py on the test sets using the ensemble with the -d flag set to output decision values.  You may need to change the path depending on where you have deepMirCut/testModel.py installed.
python3 ~/deepMirCut/testModel.py oneProd_test_3pAnnot_5pUnannot.txt -L seqBPRNA_ensemble_list.txt -o 3pAnnot_5pUnannot -d --input_setting 2
python3 ~/deepMirCut/testModel.py oneProd_test_5pAnnot_3pUnannot.txt -L seqBPRNA_ensemble_list.txt -o 5pAnnot_3pUnannot -d --input_setting 2

#graph_deltas_histogram_unannot.py takes the test sets (oneProd_test_3pAnnot_5pUnannot.txt, oneProd_test_5pAnnot_3pUnannot.txt) and the decision values returned by deepMirCut (3pAnnot_5pUnannot_classification_DVs.txt, 5pAnnot_3pUnannot_classification_DVs.txt) and an output prefix (unannotated) as input.  It produces a graph (unannotated_hist.svg) as output.
python3 Scripts/graph_deltas_histogram_unannot.py oneProd_test_3pAnnot_5pUnannot.txt oneProd_test_5pAnnot_3pUnannot.txt 3pAnnot_5pUnannot_classification_DVs.txt 5pAnnot_3pUnannot_classification_DVs.txt unannotated

#getOneProductTestSetReads.pl takes the deepMirCut dataset (mir_dataset.txt), the test file (oneProd_test_3pAnnot_5pUnannot.txt or oneProd_test_5pAnnot_3pUnannot.txt), a fasta of the human genome (hg38.fa), and a bam list containing the locations of bam files with reads mapped to the genome (bamList.txt).  It outputs a file with a list of read counts 
perl Scripts/getOneProductTestSetReads.pl mir_dataset.txt oneProd_test_3pAnnot_5pUnannot.txt hg38.fa bamList.txt
perl Scripts/getOneProductTestSetReads.pl mir_dataset.txt oneProd_test_5pAnnot_3pUnannot.txt hg38.fa bamList.txt

#graphOneProductReads.py takes a file with a list of read counts, the test set file, a file with decision values output by deepMirCut, and a side (5p or 3p).  It outputs a readcount graph as output.
python3 Scripts/graphOneProductReads.py oneProd_test_3pAnnot_5pUnannot.txt.counts oneProd_test_3pAnnot_5pUnannot.txt 3pAnnot_5pUnannot_classification_DVs.txt 5p
python3 Scripts/graphOneProductReads.py oneProd_test_5pAnnot_3pUnannot.txt.counts oneProd_test_5pAnnot_3pUnannot.txt 5pAnnot_3pUnannot_classification_DVs.txt 3p
