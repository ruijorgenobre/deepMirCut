#!/usr/bin/perl
use Bio::DB::Sam;
use strict;

my $usage = "USAGE: <mir_dataset> <oneProd_testSet.txt>\n";

my $mirDataSetFile = $ARGV[0] or die $usage;
my $oneProductTestSetFile = $ARGV[1] or die $usage;
my $genomeFile = $ARGV[2] or die $usage;
my $bamListFile = $ARGV[3] or die $usage;
my $outputFile = $oneProductTestSetFile . ".counts";

my $speciesPrefix = "hsa";


my $mirData = loadMiRDataSetSpecies($mirDataSetFile,$speciesPrefix);
my $testSetData = loadOneProductTestSet($oneProductTestSetFile);
my $testSetGenomicLocations = getTestSetGenomicLocations($testSetData,$mirData,$genomeFile);
my $testSetReads = getTestSetReads($testSetData,$testSetGenomicLocations,$bamListFile);
my $combinedTestSetReads = combineTestSetReads($testSetReads);
my $logScaleTestSetReads = convertToLogScale($combinedTestSetReads,10);
printOutputFile($outputFile,$combinedTestSetReads,$testSetGenomicLocations);

sub printOutputFile {
    my($outputFile,$combinedTestSetReads,$testSetGenomicLocations) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    foreach my $id (keys %{$combinedTestSetReads}) {
	my $location = $testSetGenomicLocations->{$id};
	print OPTF "$id\t$location\t".join(",",@{$combinedTestSetReads->{$id}}) . "\n";
    }
    close(OPTF);
}

sub convertToLogScale {
    my($combinedTestSetReads,$base) = @_;
    if ($base < 2) {
	die "base must be greater than 1";
    }
    my $logScaleTestSetReads = {};
    foreach my $id (keys %{$combinedTestSetReads}) {
	my @logScaleReads;
	for (my $itr = 0; $itr < @{$combinedTestSetReads->{$id}}; $itr++) {
	    $logScaleReads[$itr] = (${$combinedTestSetReads->{$id}}[$itr] > 0) ? log(${$combinedTestSetReads->{$id}}[$itr] + 1) / log($base) : 0;
	}
	@{$logScaleTestSetReads->{$id}} = @logScaleReads;
    }    
    return $logScaleTestSetReads;
}

sub combineTestSetReads {
    my($testSetReads) = @_;
    my $combinedTestSetReads = {};
    foreach my $id (keys %{$testSetReads}) {
	my @combinedReads;
	foreach my $sample (keys %{$testSetReads->{$id}}) {
	    for (my $itr = 0; $itr < @{$testSetReads->{$id}{$sample}}; $itr++) {
		$combinedReads[$itr] += ${$testSetReads->{$id}{$sample}}[$itr];
	    }
	}
	@{$combinedTestSetReads->{$id}} = @combinedReads; 
    }
    return $combinedTestSetReads;
}

sub getTestSetReads {
    my($testSetData,$testSetGenomicLocations,$bamListFile) = @_;
    my $testSetReads = {};
    my $bamList = loadBamList($bamListFile);
    foreach my $testSetEntry (@{$testSetData}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$newSequence) = @{$testSetEntry};
	my $location = $testSetGenomicLocations->{$id};
	getReadCounts($testSetReads,$id,$location,$bamList);
    }
    return $testSetReads;
}

sub getReadCounts {
    my($testSetReads,$id,$location,$bamList) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    foreach my $bamData (@{$bamList}) {
	my $chromLen = $stop - $start + 1;
	my @readCounts = (0) x $chromLen;
	my $readCountsLen = @readCounts;
	my($bam,$bamHeader,$bamIndex,$bamFile,$sample,$totalMapped) = @{$bamData};
	my($tid,$chromStart,$chromStop) = $bamHeader->parse_region("$chrom:$start-$stop"); #converts start and end to zero based coordinates
	my $callBack = sub {
	    my($alignment,$data) = @_;
	    my($readCounts,$chromStart,$chromStop,$chromLen,$strand) = @{$data};
	    my $rStrand = mapBamStrand($alignment->strand);
	    if ($strand eq $rStrand) {
		my $id = $alignment->qname;
		my $rStart = $alignment->start - 1;  #converted to 0-based coordinates
		my $rStop = $alignment->end - 1;  #converted to 0-based coordinates
		my $relStart = $rStart - $chromStart;
		my $relStop = $rStop - $chromStart;
		my $seq = $alignment->qseq;
		#$seq = reverseComplement($seq) if ($strand eq '-');
		#print "$seq\t$chromStart\t$chromStop\t$rStart..$rStop\t$relStart..$relStop\n";
		my $hitCount = $alignment->get_tag_values('NH');
		my $count;
		if ($id =~ /.*_x(\d+)$/) {
		    ($count) = $id =~ /.*_x(\d+)$/;
		} else {
		    $count = 1;
		}
		if ($relStop < $relStart) {
		    print "error: relStop = $relStop, relStart = $relStart\n";
		    exit()
		}		
		for (my $itr = $relStart; $itr <= $relStop; $itr++) {
		    if ($itr >= 0 && $itr < $chromLen) {
			$readCounts[$itr] += $count;
		    }
		}
	    }
	};
	my $callBackData = [\@readCounts,$chromStart,$chromStop,$chromLen,$strand];
	my $code = $bamIndex->fetch($bam,$tid,$chromStart,$chromStop,$callBack,$callBackData);
	unless (@readCounts == $readCountsLen) {
	    die "count different for $id in getReadCounts()\n";
	}
	if ($strand eq "-") {
	    @{$testSetReads->{$id}{$sample}} = reverse(@readCounts);
	} else {
	    @{$testSetReads->{$id}{$sample}} = @readCounts;
	}
    }
}

sub loadBamList {
    my($bamListFile) = @_;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	chomp;
	my($label,$bamFile) = split;
	my $bam = loadBamFile($bamFile);
	my $bamHeader = $bam->header;
	my $bamIndex = loadBamIndex($bamFile);
	my $totalMapped = getBamTotal($bamFile);
	push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,$label,$totalMapped]);
    }
    return \@bamList;
}

sub getTestSetGenomicLocations {
    my($testSetData,$mirData,$genomeFile) = @_;
    my $maxBuffer = 300;
    my %testSetGenomicLocations;
    foreach my $testSetInfo (@{$testSetData}) {
	my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$newSequence) = @{$testSetInfo};
	if (@{$mirData->{"$miRBaseId\_$name"}} > 1) {
	    print "Error: more than one entry for $miRBaseId and $name\n";
	    exit();
	}
	my($mdSpecies,$mdMiRBaseId,$mdGenomicLocation,$mdAltGenomicLocation,$mdProduct5p,$mdProduct3p,$mdDrosha5p,$mdDicer5p,$mdDicer3p,$mdDrosha3p,$mdHpStart,$mdHpStop,$mdLeftBuffer,$mdRightBuffer,$mdSeqTag,$mdSequence) = @{${$mirData->{"$miRBaseId\_$name"}}[0]};
	my($chrom,$start,$end,$strand) = parseLocation($mdAltGenomicLocation);
	my($newStart,$newStop) = (-1,-1);
	if ($strand eq "+") {
	    $start += ($maxBuffer - $mdLeftBuffer);
	    $end -= ($maxBuffer - $mdRightBuffer);
	    $newStart = $start + $mdHpStart - $hpStart;
	    $newStop = $newStart + length($newSequence) - 1;
	} elsif ($strand eq "-") {
	    $start += ($maxBuffer - $mdRightBuffer);
	    $end -= ($maxBuffer - $mdLeftBuffer);
	    $newStop = $end - $mdHpStart + $hpStart;
	    $newStart = $newStop - length($newSequence) + 1;
	} else {
	    print "error: strand $strand not recognized\n";
	    exit();
	}
	my $tsGenomicLocation = "$chrom:$newStart..$newStop:$strand";
	$testSetGenomicLocations{$id} = $tsGenomicLocation;
	my $seq = extractSequenceUsingIndex($genomeFile,$tsGenomicLocation);
	#print "strand:".$strand . "\n";
	#print $newSequence . "\n";
	#print $seq . "\n";
	if ($newSequence ne $seq) {
	    print "Error: sequences not equal\n";
	    exit();
	}
    }
    return \%testSetGenomicLocations;
}

sub loadOneProductTestSet {
    my($oneProductTestSetFile) = @_;
    my @testSetData;
    open(TSDATA,$oneProductTestSetFile) or die "failed to open $oneProductTestSetFile for reading\n";
    while (<TSDATA>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$newSequence) = split(/\t/);
	    push(@testSetData,[$id,$name,$miRBaseId,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$newSequence]);
	}
    }
    close(TSDATA);
    return(\@testSetData);
}

sub loadMiRDataSetSpecies {
    my($mirDataSetFile,$speciesPrefix) = @_;
    my %mirData;
    open(MDSF,$mirDataSetFile) or die "failed to open $mirDataSetFile\n";
    while (<MDSF>) {
	chomp;
	unless ( /^#/ ) {
	    my($species,$name,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = split(/\t/);
	    if ($species eq $speciesPrefix) {
		push(@{$mirData{"$miRBaseId\_$name"}},[$species,$miRBaseId,$genomicLocation,$altGenomicLocation,$product5p,$product3p,$drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence]);
	    }
	}
    }
    close(MDSF);
    return \%mirData;
}

sub parseLocation {
    my($location)=@_;
    my($chrom,$start,$end,$strand);    
    if($location =~ /(.*)\:(-?\d+)\-(-?\d+)\:(.*)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;
    } elsif($location =~ /(.*)\:(-?\d+)\-(-?\d+)/) {
        $chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)\:(.*)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand=$4;   	
    } elsif($location =~ /(.*)\:(-?\d+)\.\.(-?\d+)/) {
	$chrom=$1;
        $start=$2;
        $end=$3;
        $strand="+";   	
    }
    return ($chrom,$start,$end,$strand);
}    

sub reverseComplement {
# Returns the reverse complement of the input sequence.
    my($seq)=@_;
    $seq =~ tr/acgturykmbdhvACGTURYKMBDHV/tgcaayrmkvhdbTGCAAYRMKVHDB/;
    $seq=reverse($seq);
    return $seq;
}

sub extractSequenceUsingIndex {
    my($fastaFile,$location) = @_;
    my($chrom,$start,$stop,$strand) = parseLocation($location);
    my $cmd = "samtools faidx $fastaFile \'$chrom:$start-$stop\'";
    my @fastaOutput = `$cmd`;
    my $locationTag;
    my $sequence = "";
    for my $line (@fastaOutput) {
        chomp($line);
        if ($line =~ /^>/) {
            ($locationTag) = $line =~ /^>(.*)$/;
        } else {
            $sequence .= $line;
        }
    }
    if ($strand eq '-') {
        $sequence = reverseComplement($sequence);
    }
    $sequence =~ tr/atTucg/AUUUCG/;
    if ($sequence eq "") {
        die "Error: failed to find sequence for $location\n";
    }
    return $sequence;
}


sub loadBamFile {
    my($bamFile) = @_;
    my $bam = Bio::DB::Bam->open( $bamFile );
    return $bam;
}

sub loadBamIndex {
    my($bamFile) = @_;
    my $reIndex;  #changed to 1 if the index file doesn't exist
    my $bamIndex =  Bio::DB::Bam->index($bamFile,$reIndex);
    die "failed to load index for $bamFile\n" if ($reIndex);
    return $bamIndex;
}

sub loadBamList {
    my($bamListFile) = @_;
    open(BLF,$bamListFile) or die "could not open $bamListFile\n";
    my @bamList;
    while(<BLF>) {
	chomp;
	my($label,$bamFile) = split;
	my $bam = loadBamFile($bamFile);
	my $bamHeader = $bam->header;
	my $bamIndex = loadBamIndex($bamFile);
	my $totalMapped = getBamTotal($bamFile);
	push(@bamList,[$bam,$bamHeader,$bamIndex,$bamFile,$label,$totalMapped]);
    }
    return \@bamList;
}

sub getBamTotal {
    my($bamFile) = @_;
    my $bamFlagstatFile = $bamFile . ".flagstat";
    unless(-e $bamFlagstatFile) {
	# if flagstat file hasn't been created, make one.
	system("samtools flagstat $bamFile > $bamFlagstatFile");
    }    
    open(BFS,$bamFlagstatFile) or die "Could not open $bamFlagstatFile.\n";
    while(<BFS>) {
	if(/(\d+) \+ \d+ mapped/) {
	    close(BFS);
	    return $1;
	}
    }
    close(BFS);
    # if we made it here, something went wrong with parsing the flagstat file
    die "Error. Could not parse the flagstat file here: $bamFlagstatFile. Older version of samtools perhaps?\n"
}

sub mapBamStrand {
    my $strand = shift;
    if($strand eq "-1") {
	return "-";
    }
    return "+";
}
