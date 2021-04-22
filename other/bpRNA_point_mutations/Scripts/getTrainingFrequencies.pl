#!/usr/bin/perl -w
use strict;

my $parameters = {
    downstreamDist => 10,
    upstreamDist => 10,
};

my $usage = "USAGE:\n$0 <train file>\n";

my $trainFile = $ARGV[0] or die $usage;

my $trainSet = readTrainSetFile($trainFile);

my $outputCountsFile = "train_counts.txt";
my $outputPercentsFile = "train_frequencies.txt";
my $outputGridFileDR5 = "DR5_trainPerc.txt";
my $outputGridFileDC5 = "DC5_trainPerc.txt";
my $outputGridFileDC3 = "DC3_trainPerc.txt";
my $outputGridFileDR3 = "DR3_trainPerc.txt";

my($unmutatedDR5,$unmutatedDR5nucs) = getUnmutated($trainSet,"DR5",$parameters);
my($unmutatedDC5,$unmutatedDC5nucs) = getUnmutated($trainSet,"DC5",$parameters);
my($unmutatedDC3,$unmutatedDC3nucs) = getUnmutated($trainSet,"DC3",$parameters);
my($unmutatedDR3,$unmutatedDR3nucs) = getUnmutated($trainSet,"DR3",$parameters);
my($countsDR5,$percentDR5) = getFldCounts($unmutatedDR5,$parameters);
my($countsDC5,$percentDC5) = getFldCounts($unmutatedDC5,$parameters);
my($countsDC3,$percentDC3) = getFldCounts($unmutatedDC3,$parameters);
my($countsDR3,$percentDR3) = getFldCounts($unmutatedDR3,$parameters);

printOutputFile($outputCountsFile,$countsDR5,$countsDC5,$countsDC3,$countsDR3,$parameters);
printOutputFile($outputPercentsFile,$percentDR5,$percentDC5,$percentDC3,$percentDR3,$parameters);
printGridOutputFile($outputGridFileDR5,$percentDR5);
printGridOutputFile($outputGridFileDC5,$percentDC5);
printGridOutputFile($outputGridFileDC3,$percentDC3);
printGridOutputFile($outputGridFileDR3,$percentDR3);

sub printGridOutputFile {
    my($gridOutputFile,$percents) = @_;
    open(GOPTF,">$gridOutputFile") or die "failed to open $gridOutputFile for writing\n";
    my @bases = ('(',')','B','I','H');
    for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 	
	for (my $bIdx = 0; $bIdx < @bases; $bIdx++) {
	    my $perc = $percents->{$bases[$bIdx]}{$pos} * 100;
	    print GOPTF "$pos\t$bIdx\t$perc\n";
	}    
    }
    close(GOPTF);
}

sub printOutputFile {
    my($outputFile,$countsDR5,$countsDC5,$countsDC3,$countsDR3,$parameters) = @_;
    my @bases = ('(',')','B','I','H');
    my $positions = "";
    for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 	
	$positions .= ($pos <= 0) ? "\t" . ($pos - 1) : "\t$pos";
    }
    $positions .= "\n";
    my $linesDR5 = "Drosha 5p-arm\n$positions";
    foreach my $bpRNAChar (@bases) {
 	$linesDR5 .= $bpRNAChar;
 	for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 
	    if ($countsDR5->{$bpRNAChar}{$pos} =~ /\D/) {
		$linesDR5 .= sprintf("\t%.3f",$countsDR5->{$bpRNAChar}{$pos});
	    } else {
		$linesDR5 .=  "\t$countsDR5->{$bpRNAChar}{$pos}";
	    }
 	}
	$linesDR5 .= "\n";
    }
    my $linesDC5 = "Dicer 5p-arm\n$positions";
    foreach my $bpRNAChar (@bases) {
 	$linesDC5 .= $bpRNAChar;
 	for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 
	    if ($countsDC5->{$bpRNAChar}{$pos} =~ /\D/) {
		$linesDC5 .= sprintf("\t%.3f",$countsDC5->{$bpRNAChar}{$pos});
	    } else {
		$linesDC5 .= "\t$countsDC5->{$bpRNAChar}{$pos}";
	    }
 	}
	$linesDC5 .= "\n";
    }
    my $linesDC3 = "Dicer 3p-arm\n$positions";
    foreach my $bpRNAChar (@bases) {
 	$linesDC3 .= $bpRNAChar;
 	for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 
	    if ($countsDC3->{$bpRNAChar}{$pos} =~ /\D/) {
		$linesDC3 .= sprintf("\t%.3f",$countsDC3->{$bpRNAChar}{$pos});
	    } else {
		$linesDC3 .= "\t$countsDC3->{$bpRNAChar}{$pos}";
	    }
 	}
	$linesDC3 .= "\n";
    }
    my $linesDR3 = "Drosha 3p-arm\n$positions";
    foreach my $bpRNAChar (@bases) {
 	$linesDR3 .= $bpRNAChar;
 	for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) {
	    if ($countsDR3->{$bpRNAChar}{$pos} =~ /\D/) {
		$linesDR3 .= sprintf("\t%.3f",$countsDR3->{$bpRNAChar}{$pos});
	    } else {
		$linesDR3 .= "\t$countsDR3->{$bpRNAChar}{$pos}";
	    }
 	}
	$linesDR3 .= "\n";
    }
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n"; 
    print OPTF $linesDR5 . "\n";
    print OPTF $linesDC5 . "\n";
    print OPTF $linesDC3 . "\n";
    print OPTF $linesDR3 . "\n";
    close(OPTF);
}

sub getFldCounts {
    my($unmutatedFlds,$parameters) = @_;
    my $count = 0;
    my %unmutatedCount;
    my %unmutatedPercent;
    foreach my $baseId (keys %{$unmutatedFlds}) {
	$count += 1;
	for my $pos (sort {$a <=> $b} keys %{$unmutatedFlds->{$baseId}}) {
	    my $bpRNAChar = $unmutatedFlds->{$baseId}{$pos};
	    if ($unmutatedCount{$bpRNAChar}{$pos}) {
		$unmutatedCount{$bpRNAChar}{$pos} += 1;
	    } else {
		$unmutatedCount{$bpRNAChar}{$pos} = 1;
	    }
	}
    }
    foreach my $bpRNAChar (keys %unmutatedCount) {
	foreach my $pos (keys %{$unmutatedCount{$bpRNAChar}}) {
	    $unmutatedPercent{$bpRNAChar}{$pos} = $unmutatedCount{$bpRNAChar}{$pos} / $count;
	}
    }    
    return (\%unmutatedCount,\%unmutatedPercent);
}

sub getUnmutated {
    my($trainSet,$cutsite,$parameters) = @_;
    my %unmutatedFlds;
    my %unmutatedNuc;
    foreach my $baseId (keys %{$trainSet}) {
	my($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$trainSet->{$baseId}};
	my $tBPRNA = combineDotBracketAndBPRNA($fold,$bpRNAOutput);
	my $bpRNALen = length($tBPRNA);
	my @tBPRNAChars = split('',$tBPRNA);
	my @nucChars = split('',$seq);
	for (my $pos = -$parameters->{"downstreamDist"} + 1; $pos <= $parameters->{"upstreamDist"}; $pos++) { 
	    my($gPos,$gStop) = (-1,-1);
	    if ($cutsite eq "DR5") {
		($gPos,$gStop) = readCut($drosha5p);
		$gPos += $pos;
	    } elsif ($cutsite eq "DC5") {
		($gPos,$gStop) = readCut($dicer5p);
		$gPos += $pos;
	    } elsif ($cutsite eq "DC3") {
		($gPos,$gStop) = readCut($dicer3p);
		$gPos += $pos;
	    } elsif ($cutsite eq "DR3") {
		($gPos,$gStop) = readCut($drosha3p);
		$gPos += $pos;
	    } else {
		die "Could not determine cutsite for $baseId $cutsite\n";
	    }
	    if ($gPos >= 0 && $gPos < $bpRNALen) {
		$unmutatedFlds{$baseId}{$pos} = $tBPRNAChars[$gPos];
		$unmutatedNuc{$baseId}{$pos} = $nucChars[$gPos];
	    }
	}
    }
    return(\%unmutatedFlds,\%unmutatedNuc);
}

sub readTrainSetFile {
    my($trainSetFile) = @_;
    my %trainSet;
    my %loadedIds;
    open(VSF,$trainSetFile) or die "failed to open $trainSetFile for reading\n";
    while(<VSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = split(/\t/);
	    if ($loadedIds{$id}) {
		die "$id already loaded in train set\n";
	    }
	    $loadedIds{$id} = 1;
	    @{$trainSet{$id}} = ($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput);
	}
    }
    close(VSF);
    return \%trainSet;
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

sub min {
    my($x,$y) = @_;
    if($x < $y) {
        return $x;
    } else {
        return $y;
    }
}

sub max {
    my($x,$y) = @_;
    if($x > $y) {
        return $x;
    } else {
        return $y;
    }
}




# foreach my $baseId (keys %{$unmutatedDR5}) {
#     my $leftStr5p = "";
#     for my $pos (sort {$a <=> $b} keys %{$unmutatedDR5->{$baseId}}) {
# 	$leftStr5p .= $unmutatedDR5nucs->{$baseId}{$pos};
#     }
#     my $rightStr5p = "";
#     for my $pos (sort {$a <=> $b} keys %{$unmutatedDC5->{$baseId}}) {
# 	$rightStr5p .= $unmutatedDC5nucs->{$baseId}{$pos};
#     }
#     my $leftStr3p = "";
#     for my $pos (sort {$a <=> $b} keys %{$unmutatedDC3->{$baseId}}) {
# 	$leftStr3p .= $unmutatedDC3nucs->{$baseId}{$pos};
#     }
#     my $rightStr3p = "";
#     for my $pos (sort {$a <=> $b} keys %{$unmutatedDR3->{$baseId}}) {
# 	$rightStr3p .= $unmutatedDR3nucs->{$baseId}{$pos};
#     }
#     my($name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$trainSet->{$baseId}};
#     print "$name\t$leftStr5p\t$rightStr5p\t$leftStr3p\t$rightStr3p\n";
# }
