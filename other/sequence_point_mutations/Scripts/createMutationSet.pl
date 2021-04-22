#!/usr/bin/perl -w
use strict;

my $parameters = {
    downstreamDist => 5,
    upstreamDist => 5,
};

my $usage = "USAGE:\n$0 <validation set>\n";

my $validationSetFile = $ARGV[0] or die $usage;

my $drosha5pMutationFile = $validationSetFile . "_drosha5pMutations.txt";
my $dicer5pMutationFile = $validationSetFile . "_dicer5pMutations.txt";
my $dicer3pMutationFile = $validationSetFile . "_dicer3pMutations.txt";
my $drosha3pMutationFile = $validationSetFile . "_drosha3pMutations.txt";

my $validationSet = readValidationSetFile($validationSetFile);
generateMutations($validationSet,$drosha5pMutationFile,$dicer5pMutationFile,$dicer3pMutationFile,$drosha3pMutationFile,$parameters);

sub generateMutations {
    my($validationSet,$drosha5pMutationFile,$dicer5pMutationFile,$dicer3pMutationFile,$drosha3pMutationFile,$parameters) = @_;
    my @mutationList = ('A','C','G','U');
    open(DR5MF,">$drosha5pMutationFile") or die "failed to open $drosha5pMutationFile for writing\n";
    open(DC5MF,">$dicer5pMutationFile") or die "failed to open $dicer5pMutationFile for writing\n";
    open(DC3MF,">$dicer3pMutationFile") or die "failed to open $dicer3pMutationFile for writing\n";
    open(DR3MF,">$drosha3pMutationFile") or die "failed to open $drosha3pMutationFile for writing\n";
    foreach my $mirInfo (@{$validationSet}) {
	my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$mirInfo};
	my $seqLen = length($seq);
	my @seqChars = split(//,$seq);
	my($drosha5pStart,$drosha5pStop) = readCut($drosha5p);
	my($dicer5pStart,$dicer5pStop) = readCut($dicer5p);
	my($dicer3pStart,$dicer3pStop) = readCut($dicer3p);
	my($drosha3pStart,$drosha3pStop) = readCut($drosha3p);
	for (my $itr = -$parameters->{upstreamDist} + 1; $itr <= $parameters->{downstreamDist}; $itr++) {
	    #cut is between where itr=0 and itr=1
	    foreach my $mutation (@mutationList) {
		if ($drosha5pStart + $itr >= 0 && $drosha5pStart + $itr < $seqLen && $seqChars[$drosha5pStart + $itr] ne $mutation && checkNucleotideSet($seqChars[$drosha5pStart + $itr],\@mutationList)) {
		    my @seqArr = split(//,$seq);
		    $seqArr[$drosha5pStart + $itr] = $mutation;
		    my $newSeq = join("",@seqArr);
		    my $newId = $id . "_drosha5p_$itr$mutation";
		    print DR5MF "$newId\t$name\t$mirbaseId\t$product5p\t$product3p\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$newSeq\t$fold\t$bpRNAOutput\n";
		}
		if ($dicer5pStart + $itr >= 0 && $dicer5pStart + $itr < $seqLen && $seqChars[$dicer5pStart + $itr] ne $mutation && checkNucleotideSet($seqChars[$dicer5pStart + $itr],\@mutationList)) {
		    my @seqArr = split(//,$seq);
		    $seqArr[$dicer5pStart + $itr] = $mutation;
		    my $newSeq = join("",@seqArr);
		    my $newId = $id . "_dicer5p_$itr$mutation";
		    print DC5MF "$newId\t$name\t$mirbaseId\t$product5p\t$product3p\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$newSeq\t$fold\t$bpRNAOutput\n";
		}
		if ($dicer3pStart + $itr >= 0 && $dicer3pStart + $itr < $seqLen && $seqChars[$dicer3pStart + $itr] ne $mutation && checkNucleotideSet($seqChars[$dicer3pStart + $itr],\@mutationList)) {
		    my @seqArr = split(//,$seq);
		    $seqArr[$dicer3pStart + $itr] = $mutation;
		    my $newSeq = join("",@seqArr);
		    my $newId = $id . "_dicer3p_$itr$mutation";
		    print DC3MF "$newId\t$name\t$mirbaseId\t$product5p\t$product3p\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$newSeq\t$fold\t$bpRNAOutput\n";
		}
		if ($drosha3pStart + $itr >= 0 && $drosha3pStart + $itr < $seqLen && $seqChars[$drosha3pStart + $itr] ne $mutation && checkNucleotideSet($seqChars[$drosha3pStart + $itr],\@mutationList)) {
		    my @seqArr = split(//,$seq);
		    $seqArr[$drosha3pStart + $itr] = $mutation;
		    my $newSeq = join("",@seqArr);
		    my $newId = $id . "_drosha3p_$itr$mutation";
		    print DR3MF "$newId\t$name\t$mirbaseId\t$product5p\t$product3p\t$drosha5p\t$dicer5p\t$dicer3p\t$drosha3p\t$hpStart\t$hpStop\t$newSeq\t$fold\t$bpRNAOutput\n";
		}
	    }
	}
    }
    close(DR5MF);
    close(DC5MF);
    close(DC3MF);
    close(DR3MF);
}

sub readValidationSetFile {
    my($validationSetFile) = @_;
    my @validationSet;
    open(VSF,$validationSetFile) or die "failed to open $validationSetFile for reading\n";
    while(<VSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = split(/\t/);
	    push(@validationSet,[$id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput]);
	}
    }
    close(VSF);
    return \@validationSet;
}

sub readCut {
    my($cut) = @_;
    my($cutStart,$cutStop) = split(/,/,$cut);
    $cutStart--;  #convert to zero based
    $cutStop--;   #convert to zero based
    return($cutStart,$cutStop);
}

sub checkNucleotideSet {
    my($nucleotide,$nucleotideSet) = @_;
    for (my $itr = 0; $itr < @{$nucleotideSet}; $itr++) {
	if ($nucleotide eq $nucleotideSet->[$itr]) {
	    return 1;
	}
    }
    return 0;
}
