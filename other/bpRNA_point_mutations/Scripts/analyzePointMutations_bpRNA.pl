#!/usr/bin/perl -w
my $usage = "USAGE:\n$0 <point mutation file>\n";

my $pointMutationFile = $ARGV[0] or die "$usage\n";

my $outputFile = $pointMutationFile . "\_stats.txt";

my @bases = ('(',')','B','I','H');

my $pointMutationDeltas = readPointMutationFile($pointMutationFile);
my($calculations,$maxMean) = runCalculations($pointMutationDeltas,\@bases);

open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
print OPTF "#max=$maxMean\n";
foreach $calcInfo (@{$calculations}) {
    print OPTF join("\t",@{$calcInfo}) . "\n";
}
close(OPTF);

sub runCalculations {
    my($pointMutationDeltas,$bases) = @_;
    my @calculations;
    my $maxMean = -~0;
    foreach my $pos (sort {$a <=> $b} keys %{$pointMutationDeltas}) {
	for (my $bIdx = 0; $bIdx < @{$bases}; $bIdx++) {
	    my $base = $bases->[$bIdx];
	    my($mean,$variance,$count) = calculateMean($pointMutationDeltas->{$pos}{$base});
	    if (abs($mean) > $maxMean) {
		$maxMean = abs($mean);
	    }
	    push(@calculations,[$pos,$bIdx,$mean,$variance,$count]);
	}
    }
    return(\@calculations,$maxMean);
}

sub calculateMean {
    my($values) = @_;
    my $sum = 0;
    my $sumsq = 0;
    my $count = 0;
    if ($values) {
	foreach my $val (@{$values}) {
	    $sum += $val;
	    $sumsq += $val * $val;
	    $count += 1;
	}
	my $mean = $sum / $count;
	my $variance = ($sumsq - $sum*$sum/$count) / ($count - 1);
	return($mean,$variance,$count);
    }
    return(0,0,0);
}

sub readPointMutationFile {
    my($pointMutationFile) = @_;
    my %pointMutationDeltas;
    open(PMF,$pointMutationFile) or die "failed to open $pointMutationFile\n";
    while (<PMF>) {
	chomp;
	unless ( /^#/ ) {
	    my($id,$position,$base,$scoreS,$scoreT,$delta) = split(/\t/);
	    push(@{$pointMutationDeltas{$position}{$base}},$delta);
	}
    }
    close(PMF);
    return \%pointMutationDeltas;
}
