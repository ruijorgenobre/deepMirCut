#!/usr/bin/perl -w
use strict;

my $parameters = {
    downstreamDist => 10,
    upstreamDist => 10,
};

my $usage = "USAGE:\n$0 <validation set file> <mutation set file> <mirbase product fasta>\n";

my $validationSetFile = $ARGV[0] or die $usage;
my $mutationSetFile = $ARGV[1] or die $usage;
my $miRBaseProductFasta = $ARGV[2] or die $usage;

my $validationSet = readValidationSetFile($validationSetFile);
my($mutationSet,$site) = readMutationSetFile($mutationSetFile);
my $prodTestSeqs = loadMiRBaseFasta($miRBaseProductFasta);

my $outputFile = $mutationSetFile . "_counts.txt";


print $site . "\n";
if ($site eq "drosha5p") {
    my($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts) = checkDrosha5pMutations($validationSet,$mutationSet,$prodTestSeqs,$parameters);
    printOutputFile($outputFile,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
} elsif ($site eq "dicer5p") {
    my($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts) = checkDicer5pMutations($validationSet,$mutationSet,$prodTestSeqs,$parameters);
    printOutputFile($outputFile,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
} elsif ($site eq "dicer3p") {
    my($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts) = checkDicer3pMutations($validationSet,$mutationSet,$prodTestSeqs,$parameters);
    printOutputFile($outputFile,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
} elsif ($site eq "drosha3p") {
    my($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts) = checkDrosha3pMutations($validationSet,$mutationSet,$prodTestSeqs,$parameters);
    printOutputFile($outputFile,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
} else {
    die "$site not recognized\n";
}


sub printOutputFile {
    my($outputFile,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    my @mutationList = ('(',')','B','I','H');
    my @sortedPositions = sort {$a <=> $b} keys %{$mutationFrequencies};
    my @printablePostitions;
    foreach my $pos (@sortedPositions) {
	if ($pos <= 0) {
	    push(@printablePostitions,$pos-1);
	} else {
	    push(@printablePostitions,$pos)
	}
    }
    print OPTF "\t".join("\t",@printablePostitions)."\n";
    foreach my $mutationTo (@mutationList) {
	print OPTF $mutationTo;
	foreach my $pos (@sortedPositions) {
	    print OPTF "\t".$mutationFrequencies->{$pos}{$mutationTo};
	}
	print OPTF "\n";
    }
    print OPTF "\n\n\n";
    print OPTF "\t".join("\t",@printablePostitions)."\n";
    foreach my $mutationTo (@mutationList) {
	foreach my $mutationFrom (@mutationList) {
	    unless ($mutationFrom eq $mutationTo) {
		print OPTF "$mutationFrom->$mutationTo";
		foreach my $pos (@sortedPositions) {
		    print OPTF "\t".$extMutationFrequencies->{$pos}{$mutationFrom}{$mutationTo};
		}
		print OPTF "\n";
	    }
	}
    }
    print OPTF "\n\n";
    print OPTF "Unmutated Character Counts:\n";
    foreach my $unMutatedChar (keys %{$nonMutationCharCounts}) {
	print OPTF "$unMutatedChar\t$nonMutationCharCounts->{$unMutatedChar}\n";
    }
    close(OPTF);
}

sub checkDrosha5pMutations {
    my($validationSet,$mutationSet,$prodTestSeqs,$parameters) = @_;
    my $nonMutationCharCounts = {};
    my $mutationFrequencies = {};
    my $extMutationFrequencies = {};
    foreach my $validationSetEntry (@{$validationSet}) {
	my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$validationSetEntry};
	my($drosha5pCutL,$drosha5pCutR) = readCut($drosha5p);
	my($dicer5pCutL,$dicer5pCutR) = readCut($dicer5p);
	my $prod5pTestSeq = substr($seq,$drosha5pCutR,$dicer5pCutL-$drosha5pCutR+1);
	if ($prod5pTestSeq ne $prodTestSeqs->{$product5p}) {
	    die "Error: product sequences not equal for $product5p: $prod5pTestSeq vs $prodTestSeqs->{$product5p}\n"; 
	}
	my $downStreamStart = ($drosha5pCutL - ($parameters->{downstreamDist} - 1) > 0) ? $drosha5pCutL - $parameters->{downstreamDist} + 1 : 0;
	my $upStreamStop = ($drosha5pCutR + ($parameters->{downstreamDist} - 1) < length($seq) ) ? $drosha5pCutR + ($parameters->{downstreamDist} - 1) : length($seq) - 1;
	if ($upStreamStop - $downStreamStart + 1 < 20) {
	    print "$id\t$downStreamStart,$upStreamStop\n";
	}
	#print $id . "\n";
	checkMutations($drosha5pCutL,$downStreamStart,$upStreamStop,$fold,$bpRNAOutput,$mutationSet->{$id},$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts,$parameters);
    }
    return($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
}

sub checkDicer5pMutations {
    my($validationSet,$mutationSet,$prodTestSeqs,$parameters) = @_;
    my $nonMutationCharCounts = {};
    my $mutationFrequencies = {};
    my $extMutationFrequencies = {};
    foreach my $validationSetEntry (@{$validationSet}) {
	my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$validationSetEntry};
	my($drosha5pCutL,$drosha5pCutR) = readCut($drosha5p);
	my($dicer5pCutL,$dicer5pCutR) = readCut($dicer5p);
	my $prod5pTestSeq = substr($seq,$drosha5pCutR,$dicer5pCutL-$drosha5pCutR+1);
	if ($prod5pTestSeq ne $prodTestSeqs->{$product5p}) {
	    die "Error: product sequences not equal for $product5p: $prod5pTestSeq vs $prodTestSeqs->{$product5p}\n"; 
	}
	my $downStreamStart = ($dicer5pCutL - ($parameters->{downstreamDist} - 1) > 0) ? $dicer5pCutL - $parameters->{downstreamDist} + 1 : 0;
	my $upStreamStop = ($dicer5pCutR + ($parameters->{downstreamDist} - 1) < length($seq) ) ? $dicer5pCutR + ($parameters->{downstreamDist} - 1) : length($seq) - 1;
	if ($upStreamStop - $downStreamStart + 1 < 20) {
	    print "$id\t$downStreamStart,$upStreamStop\n";
	}
	#print $id . "\n";
	checkMutations($dicer5pCutL,$downStreamStart,$upStreamStop,$fold,$bpRNAOutput,$mutationSet->{$id},$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts,$parameters);
    }
    return($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
}

sub checkDicer3pMutations {
    my($validationSet,$mutationSet,$prodTestSeqs,$parameters) = @_;
    my $nonMutationCharCounts = {};
    my $mutationFrequencies = {};
    my $extMutationFrequencies = {};
    foreach my $validationSetEntry (@{$validationSet}) {
	my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$validationSetEntry};
	my($dicer3pCutL,$dicer3pCutR) = readCut($dicer3p);
	my($drosha3pCutL,$drosha3pCutR) = readCut($drosha3p);
	my $prod3pTestSeq = substr($seq,$dicer3pCutR,$drosha3pCutL-$dicer3pCutR+1);
	if ($prod3pTestSeq ne $prodTestSeqs->{$product3p}) {
	    die "Error: product sequences not equal for $product3p: $prod3pTestSeq vs $prodTestSeqs->{$product3p}\n"; 
	}
	my $downStreamStart = ($dicer3pCutL - ($parameters->{downstreamDist} - 1) > 0) ? $dicer3pCutL - $parameters->{downstreamDist} + 1 : 0;
	my $upStreamStop = ($dicer3pCutR + ($parameters->{downstreamDist} - 1) < length($seq) ) ? $dicer3pCutR + ($parameters->{downstreamDist} - 1) : length($seq) - 1;
	if ($upStreamStop - $downStreamStart + 1 < 20) {
	    print "$id\t$downStreamStart,$upStreamStop\n";
	}
	#print $id . "\n";
	checkMutations($dicer3pCutL,$downStreamStart,$upStreamStop,$fold,$bpRNAOutput,$mutationSet->{$id},$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts,$parameters);
    }
    return($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
}

sub checkDrosha3pMutations {
    my($validationSet,$mutationSet,$prodTestSeqs,$parameters) = @_;
    my $nonMutationCharCounts = {};
    my $mutationFrequencies = {};
    my $extMutationFrequencies = {};
    foreach my $validationSetEntry (@{$validationSet}) {
	my($id,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = @{$validationSetEntry};
	my($dicer3pCutL,$dicer3pCutR) = readCut($dicer3p);
	my($drosha3pCutL,$drosha3pCutR) = readCut($drosha3p);
	my $prod3pTestSeq = substr($seq,$dicer3pCutR,$drosha3pCutL-$dicer3pCutR+1);
	if ($prod3pTestSeq ne $prodTestSeqs->{$product3p}) {
	    die "Error: product sequences not equal for $product3p: $prod3pTestSeq vs $prodTestSeqs->{$product3p}\n"; 
	}
	my $downStreamStart = ($drosha3pCutL - ($parameters->{downstreamDist} - 1) > 0) ? $drosha3pCutL - $parameters->{downstreamDist} + 1 : 0;
	my $upStreamStop = ($drosha3pCutR + ($parameters->{downstreamDist} - 1) < length($seq) ) ? $drosha3pCutR + ($parameters->{downstreamDist} - 1) : length($seq) - 1;
	if ($upStreamStop - $downStreamStart + 1 < 20) {
	    print "$id\t$downStreamStart,$upStreamStop\n";
	}
	#print $id . "\n";
	checkMutations($drosha3pCutL,$downStreamStart,$upStreamStop,$fold,$bpRNAOutput,$mutationSet->{$id},$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts,$parameters);
    }
    return($mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts);
}

sub checkMutations {
    my($cutsite,$downStreamStart,$upStreamStop,$fold,$bpRNAOutput,$mutationSet,$mutationFrequencies,$extMutationFrequencies,$nonMutationCharCounts,$parameters) = @_;
    my $modifiedBPRNA = combineDotBracketAndBPRNA($fold,$bpRNAOutput);
    my @mutationList = ('(',')','B','I','H');
    my %mutationListHash = map { $_ => 1 } @mutationList;
    my %mutationTestMatrix;
    my @mutationSegment = split('',substr($modifiedBPRNA,$downStreamStart,$upStreamStop-$downStreamStart+1));
    for (my $itr = 0; $itr < @mutationSegment; $itr++) {
	if ($mutationListHash{$mutationSegment[$itr]}) {
	    foreach my $mutchar (@mutationList) {
		if ($mutchar ne $mutationSegment[$itr]) {
		    #print "$mutationSegment[$itr],$itr,$mutchar\n";
		    $mutationTestMatrix{$itr}{$mutchar} = 0;
		}
	    }
	} else {
	    if ($nonMutationCharCounts->{$mutationSegment[$itr]}) {
		$nonMutationCharCounts->{$mutationSegment[$itr]} += 1;
		#print $mutationSegment[$itr] . "\n";
	    } else {
		$nonMutationCharCounts->{$mutationSegment[$itr]} = 1;
	    }
	}
    }
    foreach my $mutationSetEntry (@{$mutationSet}) {
	my($mutationId,$mname,$mmirbaseId,$mproduct5p,$mproduct3p,$mdrosha5p,$mdicer5p,$mdicer3p,$mdrosha3p,$mhpStart,$mhpStop,$mseq,$mfold,$mbpRNAOutput) = @{$mutationSetEntry};
	my $mutationSeq = combineDotBracketAndBPRNA($mfold,$mbpRNAOutput);
	my($mid,$mType,$mutationStr) = $mutationId =~ /(.*)\_(.*)\_(.*)$/;
	my(@mutationChars) = split('',$mutationStr);
	my $mutation = pop(@mutationChars);
	my $mutationPos = int(join('',@mutationChars));
	#print $mutationId . "\t".$mid."\t".$mutationPos."\t".$mutation . "\n";
	my($mutationFrom,$mutationTo) = checkMutationSite($modifiedBPRNA,$mutationSeq,$cutsite + $mutationPos);
	$mutationFrequencies->{$mutationPos}{$mutationTo} += 1;
	$extMutationFrequencies->{$mutationPos}{$mutationFrom}{$mutationTo} += 1;
	$mutationTestMatrix{$mutationPos + $parameters->{downstreamDist} - 1}{$mutationTo} += 1;
    }
    foreach my $mutationPos (keys %mutationTestMatrix) {
	foreach my $mutationChar (keys %{$mutationTestMatrix{$mutationPos}}) {
	    if ($mutationTestMatrix{$mutationPos}{$mutationChar} != 1) {
		die "$mutationPos was not mutated to $mutationChar in $modifiedBPRNA\n";
	    }
	}
    }
}

sub checkMutationSite {
    my($modifiedBPRNA,$mutationSeq,$mutationSite) = @_;
    my $modifiedBPRNALeft = substr($modifiedBPRNA,0,$mutationSite);
    my $mutationSeqLeft = substr($mutationSeq,0,$mutationSite);
    my $modifiedBPRNARight = substr($modifiedBPRNA,$mutationSite+1,length($modifiedBPRNA)-$mutationSite);
    my $mutationSeqRight = substr($mutationSeq,$mutationSite+1,length($modifiedBPRNA)-$mutationSite);
    my $BPRNAChar = substr($modifiedBPRNA,$mutationSite,1);
    my $mutationChar = substr($mutationSeq,$mutationSite,1);
    #print "$modifiedBPRNA\n$modifiedBPRNALeft\t$BPRNAChar\t$modifiedBPRNARight\n";
    #print "$mutationSeq\n$mutationSeqLeft\t$mutationChar\t$mutationSeqRight\n";
    if ($modifiedBPRNALeft ne $mutationSeqLeft) {
	die "sequences left of the mutation are not the same: $modifiedBPRNALeft vs $mutationSeqLeft\n";
    }
    if ($modifiedBPRNARight ne $mutationSeqRight) {
	die "sequences right of the mutation are not the same: $modifiedBPRNARight vs $mutationSeqRight\n";
    }
    if ($BPRNAChar eq $mutationChar) {
	die "mutation is the same: $BPRNAChar vs $mutationChar at $mutationSite\n$modifiedBPRNA\n";
    }
    return($BPRNAChar,$mutationChar);
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

sub readMutationSetFile {
    my($mutationSetFile) = @_;
    my %mutationSet;
    my $site;
    my $FIRST = 1;
    open(VSF,$mutationSetFile) or die "failed to open $mutationSetFile for reading\n";
    while(<VSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($mutationId,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput) = split(/\t/);
	    my($id,$mType,$mutation) = $mutationId =~ /(.*)\_(.*)\_(.*)$/;
	    if ($FIRST) {
		$FIRST = 0;
		$site = $mType;
	    } elsif ($site ne $mType) {
		die "Error: more than one mutation type in file: $site vs $mType\n";
	    }
	    push(@{$mutationSet{$id}},[$mutationId,$name,$mirbaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$seq,$fold,$bpRNAOutput]);
	}
    }
    close(VSF);
    return(\%mutationSet,$site);
}

sub readCut {
    my($cut) = @_;
    my($cutStart,$cutStop) = split(/,/,$cut);
    $cutStart--;  #convert to zero based
    $cutStop--;   #convert to zero based
    return($cutStart,$cutStop);
}

sub loadMiRBaseFasta {
    my($miRBaseHPFasta) = @_;
    my %miRBaseSeqs;
    open(MBHPF,$miRBaseHPFasta) or die "failed to open $miRBaseHPFasta\n";
    my $seqId = "";
    while (<MBHPF>) {
	chomp;
	my $line = $_;
	if ( $line =~ /^>/ ) {
	    $line =~ s/>//;
	    ($seqId) = split(/\s+/,$line);
	    unless ($miRBaseSeqs{$seqId}) {
		$miRBaseSeqs{$seqId} = "";
	    } else {
		print "Warning: $seqId already used as a key for another sequence in $miRBaseHPFasta";
	    }
	} else {
	    $miRBaseSeqs{$seqId} .= $line;
	}
    }
    close(MBHPF);
    return \%miRBaseSeqs;
}
