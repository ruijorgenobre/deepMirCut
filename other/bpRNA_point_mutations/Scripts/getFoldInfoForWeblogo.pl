#!/usr/bin/perl -w 
use strict;

my $parameters = {
    startDist => 10,
    stopDist => 10
};

my $usage = "USAGE:\n$0 <dataset file>\n";

my $dataSetFile = $ARGV[0] or die $usage;

my $outputFileDR5 = $dataSetFile . "_DR5.fa";
my $outputFileDC5 = $dataSetFile . "_DC5.fa";
my $outputFileDC3 = $dataSetFile . "_DC3.fa";
my $outputFileDR3 = $dataSetFile . "_DR3.fa";

my $cutSiteSeqs = getSeqsFromDataset($dataSetFile,$parameters);
printOutputFile($cutSiteSeqs->{"DR5"},$outputFileDR5);
printOutputFile($cutSiteSeqs->{"DC5"},$outputFileDC5);
printOutputFile($cutSiteSeqs->{"DC3"},$outputFileDC3);
printOutputFile($cutSiteSeqs->{"DR3"},$outputFileDR3);

sub printOutputFile {
    my($seqData,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    for my $seqInfo (@{$seqData}) {
	my($id,$name,$seq) = @{$seqInfo};
	print OPTF ">$id\n";
	print OPTF "$seq\n";
    }
    close(OPTF);
}

sub getSeqsFromDataset {
    my($dataSetFile,$parameters) = @_;
    my %cutSiteSeqs;
    open(DSF,$dataSetFile) or die "failed to open $dataSetFile\n";
    while (<DSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$miRBaseId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$seq,$fold,$bpRNA) = split(/\t/);
	    #Drosha 5p cutsite
	    my($DR5LeftN,$DR5SeqStart,$DR5SeqStop,$DR5RightN) = getSeqSpan($DR5,$seq,$parameters);
	    my $DR5Seq = 'N' x $DR5LeftN . getBPRNAString($fold,$bpRNA,$DR5SeqStart,$DR5SeqStop-$DR5SeqStart+1) . 'N' x $DR5RightN;
	    #if ($DR5LeftN || $DR5RightN) {
		#print "DR5\tright\t$DR5LeftN\t$DR5RightN\t$DR5Seq\t$prod5p\n";
	    #}
	    push(@{$cutSiteSeqs{"DR5"}},[$id,$name,$DR5Seq]);
	    #Dicer 5p cutsite
	    my($DC5LeftN,$DC5SeqStart,$DC5SeqStop,$DC5RightN) = getSeqSpan($DC5,$seq,$parameters);
	    my $DC5Seq = 'N' x $DC5LeftN . getBPRNAString($fold,$bpRNA,$DC5SeqStart,$DC5SeqStop-$DC5SeqStart+1) . 'N' x $DC5RightN;
	    #if ($DC5LeftN || $DC5RightN) {
		#print "DC5\tleft\t$DC5LeftN\t$DC5RightN\t$DC5Seq\t$prod5p\n";
	    #}
	    push(@{$cutSiteSeqs{"DC5"}},[$id,$name,$DC5Seq]);
	    #Dicer 3p cutsite
	    my($DC3LeftN,$DC3SeqStart,$DC3SeqStop,$DC3RightN) = getSeqSpan($DC3,$seq,$parameters);
	    my $DC3Seq = 'N' x $DC3LeftN . getBPRNAString($fold,$bpRNA,$DC3SeqStart,$DC3SeqStop-$DC3SeqStart+1) . 'N' x $DC3RightN;
	    #if ($DC3LeftN || $DC3RightN) {
		#print "DC3\tright\t$DC3LeftN\t$DC3RightN\t$DC5Seq\t$prod3p\n";
	    #}
	    push(@{$cutSiteSeqs{"DC3"}},[$id,$name,$DC3Seq]);
	    #Drosha 3p cutsite
	    my($DR3LeftN,$DR3SeqStart,$DR3SeqStop,$DR3RightN) = getSeqSpan($DR3,$seq,$parameters);
	    my $DR3Seq = 'N' x $DR3LeftN . getBPRNAString($fold,$bpRNA,$DR3SeqStart,$DR3SeqStop-$DR3SeqStart+1) . 'N' x $DR3RightN;
	    #if ($DR3LeftN || $DR3RightN) {
		#print "DR3\tleft\t$DR3LeftN\t$DR3RightN\t$DR3Seq\t$prod3p\n";
	    #}
	    push(@{$cutSiteSeqs{"DR3"}},[$id,$name,$DR3Seq]);
	}
    }
    close(DSF);
    return \%cutSiteSeqs;
}

sub getSeqSpan {
    my($cutsite,$seq,$parameters) = @_;
    my $seqLen = length($seq);
    my($cut5p,$cut3p) = readCutSite($cutsite);
    my($seqStart,$seqStop) = ($cut5p - ($parameters->{startDist} - 1),$cut3p + ($parameters->{stopDist} - 1));
    my($leftN,$rightN) = (0,0);
    if ($seqStart < 0) {
	$leftN = -$seqStart;
	$seqStart = 0;
    }
    if ($seqStop >= $seqLen) {
	$rightN = $seqStop - ($seqLen - 1);
	$seqStop = $seqLen - 1;
    }
    return($leftN,$seqStart,$seqStop,$rightN);
}

sub readCutSite {
    my($cutsite) = @_;
    my($cut5p,$cut3p) = split(/,/,$cutsite);
    $cut5p -= 1; #convert to 0-based
    $cut3p -= 1; #convert to 0-based
    return($cut5p,$cut3p);
}

sub getBPRNAString {
    my($fold,$bpRNA,$start,$len) = @_;
    my $foldStr = substr($fold,$start,$len);
    my $bpRNAStr = substr($bpRNA,$start,$len);
    my @foldChars = split(//,$foldStr);
    my @bpRNAChars = split(//,$bpRNAStr);
    my @outChars;
    for (my $itr = 0; $itr < @foldChars; $itr++) {
	if ($foldChars[$itr] eq '(') {
	    push(@outChars, 'L')
	} elsif ($foldChars[$itr] eq ')') {
	    push(@outChars, 'R')
	} else {
	    push(@outChars, $bpRNAChars[$itr])
	}
    }
    my $outStr = join('',@outChars);
    return $outStr;
}
