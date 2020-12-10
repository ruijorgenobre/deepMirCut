#!/bin/bash

dir=$(realpath "$(dirname "${BASH_SOURCE[0]}")")
prepareData=$dir"/prepareData.pl"
echo "deep-mir-cut's prepareData.pl script requires bpRNA to operate (https://github.com/hendrixlab/bpRNA)."
echo ""
while : ; do
    echo -n "Plese enter the full path for bpRNA.pl: "
    read bpRNA
    if test -e "$bpRNA"; then
	sed -i '3s|.*|do "'$bpRNA'";|g' $prepareData
	break;
    else
	echo "file not found: $bpRNA"
    fi
done


seqBPRNAList=$dir'/seqBPRNA_ensemble_list.txt' 
sed -E -i 's|^.*ensemble_models\/model_([0-9]*)_([0-9]*).model|'$dir'\/ensemble_models\/model_\1_\2.model|' $seqBPRNAList


