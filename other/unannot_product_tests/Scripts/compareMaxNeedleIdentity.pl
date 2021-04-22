#!/usr/bin/perl -w

my $usage = "USAGE:\n$0 <fasta1> <fasta2> <output prefix>\n";

my $fasta1 = $ARGV[0] or die $usage;
my $fasta2 = $ARGV[1] or die $usage;
my $outputPrefix = $ARGV[2] or die $usage;

my $outputFile = $outputPrefix . ".txt";

my $sequences1 = readFastaFile($fasta1);
my $maxIdentities = getMaxIdentities($sequences1,$fasta2);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
foreach my $identityInfo (@{$maxIdentities}) {
    my($id1,$id2,$maxIdentity) = @{$identityInfo};
    print OPTF "$id1\t$id2\t$maxIdentity\n";
}
close(OPTF);

sub getMaxIdentities {
    my($sequences1,$fasta2) = @_;
    my @maxIdentities;
    my $hasEqualIds = 0;
    foreach my $id (keys %{$sequences1}) {
	my $sequence = "asis:" . $sequences1->{$id};
	my $alignments = performNeedleAlignments($sequence,$fasta2);
	my $maxIdentity = -1;
	my $maxId2;
	foreach my $id2 (keys %{$alignments}) {
	    my($identity,$similarity,$gaps,$score) = @{$alignments->{$id2}};
	    if ($id eq $id2) {
		$hasEqualIds = 1;
	    }
	    if ($identity > $maxIdentity && $id ne $id2) {
		$maxIdentity = $identity;
		$maxId2 = $id2;
	    }
	}
	push(@maxIdentities,[$id,$maxId2,$maxIdentity]);
	#print "$id\t$maxId2\t$maxIdentity\n";
    }
    if ($hasEqualIds) {
	print "warning: both fasta files have ids which are the same.\n";
    }
    return \@maxIdentities;
}

sub performNeedleAlignments {
    my($sequenceA,$sequenceB) = @_;
    my %alignments;
    my $cmd = "needle -asequence \"$sequenceA\" -bsequence \"$sequenceB\" -auto -stdout";
    open(NDLE, "$cmd |");
    while (<NDLE>) {
	chomp;
	my $line = $_;
	if ( $line =~ /^#\s+2:/ ) {
	    my($id) = $line =~ /^#\s+2:\s+(.*)/;
	    my $identFound = 0;
	    while (<NDLE>) {
		chomp;
		$line = $_;
		if ( $line =~ /^#\s+2:/ ) {
		    die "needle output improperly parsed ($id)\n";
		} elsif ( $line =~ /^#\s+Identity:/ ) {
		    $identFound = 1;
		    last;
		}
	    }
	    unless ($identFound) {
		die "failed to find line for identity\n";
	    }
	    my($identity) = $line =~ /^#\s+Identity:\s+.+?\s+\(\s*(.*)%\)/;
	    my $simFound = 0;
	    while (<NDLE>) {
		chomp;
		$line = $_;
		if ( $line =~ /^#\s+2:/ ) {
		    die "needle output improperly parsed ($id)\n";
		} elsif ( $line =~ /^#\s+Similarity:/ ) {
		    $simFound = 1;
		    last;
		}
	    }
	    unless ($simFound) {
		die "failed to find line for similarity\n";
	    }
	    my($similarity) = $line =~ /^#\s+Similarity:\s+.+?\s+\(\s*(.*)%\)/;
	    my $gapsFound = 0;
	    while (<NDLE>) {
		chomp;
		$line = $_;
		if ( $line =~ /^#\s+2:/ ) {
		    die "needle output improperly parsed ($id)\n";
		} elsif ( $line =~ /^#\s+Gaps:/ ) {
		    $gapsFound = 1;
		    last;
		}
	    }
	    unless ($gapsFound) {
		die "failed to find line for gaps\n";
	    }
	    my($gaps) = $line =~ /^#\s+Gaps:\s+.+?\s+\(\s*(.*)%\)/;
	    my $scoreFound = 0;
	    while (<NDLE>) {
		chomp;
		$line = $_;
		if ( $line =~ /^#\s+2:/ ) {
		    die "needle output improperly parsed ($id)\n";
		} elsif ( $line =~ /^#\s+Score:/ ) {
		    $scoreFound = 1;
		    last;
		}
	    }
	    unless ($scoreFound) {
		die "failed to find line for score\n";
	    }
	    my($score) = $_ =~ /^#\s+Score:\s+(.*)/;
	    @{$alignments{$id}} = ($identity,$similarity,$gaps,$score);
	}
    }
    close(NDLE);
    return \%alignments;
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
