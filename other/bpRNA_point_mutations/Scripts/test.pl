#!/usr/bin/perl -w
use strict;

my $parameters = {
    downstreamDist => 10,
    upstreamDist => 10,
};

my $usage = "USAGE:\n$0 <validation File> <cutsite [DR5|DC5|DC3|DR3]>\n";

my $validationSetFile = $ARGV[0] or die $usage;
my $cutsite = $ARGV[1] or die $usage;


my $validationSet = readValidationSetFile($validationSetFile);
my $unmutatedBases = getUnmutated($validationSet,$cutsite,$parameters);


exit();


sub getUnmutated {
    my($validationSet,$cutsite,$parameters) = @_;
    my %unmutatedBases;
    foreach my $baseId (keys %{$validationSet}) {
	my($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$validationSet->{$baseId}};
	my $tBPRNA = combineDotBracketAndBPRNA($fold,$bpRNAOutput);
	my @tBPRNAChars = split('',$tBPRNA);
	my @tSeqChars = split('',$seq);
	for (my $mutationPos = -$parameters->{"downstreamDist"} + 1; $mutationPos <= $parameters->{"upstreamDist"}; $mutationPos++) {
	    
	    #double check for off-by-one error;
	    
	    my($gPos) = -1;
	    if ($cutsite eq "DR5") {
		($gPos) = readCut($drosha5p);
		$gPos += $mutationPos;
	    } elsif ($cutsite eq "DC5") {
		($gPos) = readCut($dicer5p);
		$gPos += $mutationPos;
	    } elsif ($cutsite eq "DC3") {
		($gPos) = readCut($dicer3p);
		$gPos += $mutationPos;
	    } elsif ($cutsite eq "DR3") {
		($gPos) = readCut($drosha3p);
		$gPos += $mutationPos;
	    } else {
		die "Could not determine cutsite for $baseId $cutsite\n";
	    }
	    #$unmutatedBases{$baseId}{$mutationPos} = $tBPRNAChars[$gPos];
	    $unmutatedBases{$baseId}{$mutationPos} = $tSeqChars[$gPos];
	}
    }
    foreach my $baseId (sort {$a cmp $b} keys %unmutatedBases) {
	my $str = "";
	foreach my $pos (sort {$a <=> $b} keys %{$unmutatedBases{$baseId}}) {
	    $str .= $unmutatedBases{$baseId}{$pos};
	}
	print "$baseId\t$str\n";
    }


    return \%unmutatedBases;
}

sub readCut {
    my($cut) = @_;
    my($cutStart,$cutStop) = split(/,/,$cut);
    $cutStart--;  #convert to zero based
    $cutStop--;   #convert to zero based
    return($cutStart,$cutStop);
}

sub combineDotBracketAndBPRNA {
    my($fold,$bpRNAOutput) = @_;
    if (length($fold) != length($bpRNAOutput)) {
	die "Error: length of fold and $bpRNAOutput different";
    }
    my @foldChars = split('',$fold);
    my @bpRNAChars = split('',$bpRNAOutput);
    for (my $itr = 0; $itr < @foldChars; $itr++) {
	if ($foldChars[$itr] ne '(' && $foldChars[$itr] ne ')') {
	    $foldChars[$itr] = $bpRNAChars[$itr];
	}
    }
    my $modifiedBPRNA = join('',@foldChars);
    return $modifiedBPRNA;
}

sub readScoresFile {
    my($scoresFile,$cutsite) = @_;
    my %scores;
    open(SCRS,$scoresFile) or die "failed to open  $scoresFile\n";
    while (<SCRS>) {
	chomp;
	unless ( /^#/ ) {
	    my($baseId,$name,$DR5,$DC5,$DC3,$DR3) = split(/\t/);
	    if ($cutsite eq "DR5") {
		$scores{$baseId} = $DR5;
	    }
	    if ($cutsite eq "DC5") {
		$scores{$baseId} = $DC5;
	    }
	    if ($cutsite eq "DC3") {
		$scores{$baseId} = $DC3;
	    }
	    if ($cutsite eq "DR3") {
		$scores{$baseId} = $DR3;
	    }
	}
    }
    close(SCRS);
    return \%scores;
}

sub readMutationScoresFile {
    my($mutationScoresFile,$cutsite) = @_;
    my %mutationScores;
    open(CSDV,$mutationScoresFile) or die "failed to open $mutationScoresFile\n";
    while (<CSDV>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$DR5,$DC5,$DC3,$DR3) = split(/\t/);
	    my($baseId,$pos,$mutation,$relevantCut) = getMutationInfo($id);
	    if ($cutsite eq "DR5") {
		if ($relevantCut ne "drosha5p") {
		    print "Warning: $cutsite does not match $relevantCut\n";
		}
		push(@{$mutationScores{$baseId}},[$pos,$mutation,$DR5]);
	    }
	    if ($cutsite eq "DC5") {
		if ($relevantCut ne "dicer5p") {
		    print "Warning: $cutsite does not match $relevantCut\n";
		}
		push(@{$mutationScores{$baseId}},[$pos,$mutation,$DC5]);
	    }
	    if ($cutsite eq "DC3") {
		if ($relevantCut ne "dicer3p") {
		    print "Warning: $cutsite does not match $relevantCut\n";
		}
		push(@{$mutationScores{$baseId}},[$pos,$mutation,$DC3]);
	    }
	    if ($cutsite eq "DR3") {
		if ($relevantCut ne "drosha3p") {
		    print "Warning: $cutsite does not match $relevantCut\n";
		}
		push(@{$mutationScores{$baseId}},[$pos,$mutation,$DR3]);
	    }
	}
    }
    close(CSDV);
    return \%mutationScores;
}

sub getMutationInfo {
    my($id) = @_;
    my(@idParts) = split(/\_/,$id);
    my($mutationInfo) = pop(@idParts);
    my($relevantCut) = pop(@idParts);
    my($baseId) = join("\_",@idParts);
    my($pos,$mutation) = $mutationInfo =~ /(^-?\d+)(.*)/;
    return($baseId,$pos,$mutation,$relevantCut);
}


sub readValidationSetFile {
    my($validationSetFile) = @_;
    my %validationSet;
    my %loadedIds;
    open(VSF,$validationSetFile) or die "failed to open $validationSetFile for reading\n";
    while(<VSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = split(/\t/);
	    if ($loadedIds{$id}) {
		die "$id already loaded in validation set\n";
	    }
	    $loadedIds{$id} = 1;
	    @{$validationSet{$id}} = ($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput);
	}
    }
    close(VSF);
    return \%validationSet;
}
