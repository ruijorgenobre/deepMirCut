#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <train set> <train set with similar> <train set with similar2> <validation set> <test set> <hairpins.fa>\n";

my $trainSetFile = $ARGV[0] or die $usage;
my $trainSetWithSimilarFile = $ARGV[1] or die $usage;
my $trainSetWithSimilarFile2 = $ARGV[2] or die $usage;
my $validationSetFile = $ARGV[3] or die $usage;
my $testSetFile = $ARGV[4] or die $usage;
my $hairpinsFasta = $ARGV[5] or die $usage;

my $trainSetOutputFasta = "trainSet_precursors.fa";
my $trainSetWithSimilarOutputFasta = "trainSet_wSimilar_precursors.fa";
my $trainSetWithSimilarOutputFasta2 = "trainSet_wSimilar_precursors_all.fa";
my $validationSetOutputFasta = "validationSet_precursors.fa";
my $testSetOutputFasta = "testSet_precursors.fa";

my $trainSetNames = getNamesFromDataSet($trainSetFile);
my $trainSetWithSimilarNames = getNamesFromDataSet($trainSetWithSimilarFile);
my $trainSetWithSimilarNames2 = getNamesFromDataSet($trainSetWithSimilarFile2);
my $validationSetNames = getNamesFromDataSet($validationSetFile);
my $testSetNames = getNamesFromDataSet($testSetFile);
my $hairpinSeqs = readFastaFile($hairpinsFasta);

printOutputFasta($trainSetNames,$hairpinSeqs,$trainSetOutputFasta);
printOutputFasta($trainSetWithSimilarNames,$hairpinSeqs,$trainSetWithSimilarOutputFasta);
printOutputFasta($trainSetWithSimilarNames2,$hairpinSeqs,$trainSetWithSimilarOutputFasta2);
printOutputFasta($validationSetNames,$hairpinSeqs,$validationSetOutputFasta);
printOutputFasta($testSetNames,$hairpinSeqs,$testSetOutputFasta);

sub printOutputFasta {
    my($names,$hairpinSeqs,$outputFasta) = @_;
    open(OPTF,">$outputFasta") or die "failed to open $outputFasta for writing\n";
    foreach my $name (@{$names}) {
	my $sequence = $hairpinSeqs->{$name};
	$sequence =~ tr/uU/tT/;
	print OPTF ">$name\n$sequence\n";
    }
    close(OPTF);
}

sub getNamesFromDataSet {
    my($dataSetFile) = @_;
    my %namesHash;
    open(DS,$dataSetFile) or die "failed to open $dataSetFile\n";
    while (<DS>) {
	chomp;
	unless (/^\#/) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$sequence) = split(/\t/);
	    $namesHash{$name} = 1;
	}
    }
    close(DS);
    my @names = keys %namesHash;
    return \@names;
}

sub readClusterFile {
    my($CDHitClusterFile) = @_;
    my %clusterHash;
    my $currentCluster;
    open(CDH,$CDHitClusterFile) or die "failed to open $CDHitClusterFile\n";
    while (<CDH>) {
	chomp;
	if ( /^>/ ) {
	    s/>//;
	    my $clusterName = $_;
	    $currentCluster = $clusterName;
	} else {
	    my($num,$aa,$name) = split(/\s+/);
	    ($name) = $name =~ /^>(.*?)\.\.\./;
	    $clusterHash{$currentCluster}{$name} = 1;
	}
    }
    close(CDH);
    my %cluster;
    foreach my $clusterName (%clusterHash) {
	foreach my $name (keys %{$clusterHash{$clusterName}}) {
	    foreach my $name2 (keys %{$clusterHash{$clusterName}}) {
		if ($name ne $name2) {
		    push(@{$cluster{$name}},$name2);
		}
	    }
	}
    }
    return \%cluster;
}

sub addClusterToNames {
    my($names,$cluster) = @_;
    my %newNamesHash;
    foreach my $name (@{$names}) {
	$newNamesHash{$name} = 1;
	if ($cluster->{$name}) {
	    foreach my $name2 (@{$cluster->{$name}}) {
		$newNamesHash{$name2} = 1;
	    }
	}
    }
    my @newNames = keys %newNamesHash;
    return \@newNames;
}

sub readFastaFile {
    my($fasta) = @_;
    my %sequences;
    open(FA,$fasta) or die "failed to open $fasta for reading\n";
    my $currTag;
    while (<FA>) {
	chomp;
	if ( /^>/ ) {
	    ($currTag) = split;
	    $currTag =~ s/>//;
	} else {
	    $sequences{$currTag} .= $_;
	}
    }
    close(FA);
    return \%sequences;
}
