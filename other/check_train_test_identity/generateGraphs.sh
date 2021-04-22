#!/bin/bash
#graphMaxIdentity.py graphs a histogram of mirs falling withing various identity thresholds
python3 Scripts/graphMaxIdentity.py testSet_vs_trainSetWSimilar.txt "Test Set vs Train Set with Similar miRs"
python3 Scripts/graphMaxIdentity.py validationSet_vs_trainSetWSimilar.txt "Validation Set vs Train Set with Similar miRs"
python3 Scripts/graphMaxIdentity.py testSet_vs_trainSet.txt "Test Set vs Train Set"
python3 Scripts/graphMaxIdentity.py validationSet_vs_trainSet.txt "Validation Set vs Train Set"
python3 Scripts/graphMaxIdentity.py validationSet_vs_trainSetWSimilarIncludingWithoutFiltering.txt "Validation Set vs Train Set with Similar - no filtering"
python3 Scripts/graphMaxIdentity.py testSet_vs_trainSetWSimilarWithoutFiltering.txt "Validation Set vs Train Set with Similar - no filtering"
