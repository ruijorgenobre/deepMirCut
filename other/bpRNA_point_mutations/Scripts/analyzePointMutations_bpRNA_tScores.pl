#!/usr/bin/perl -w
use Statistics::DependantTTest;

my $usage = "USAGE:\n$0 <point mutation file>\n";

my $pointMutationFile = $ARGV[0] or die "$usage\n";

my $tScoreOutputFile = $pointMutationFile . "\_stats_wTScore.txt";
my $pValOutputFile = $pointMutationFile . "\_stats_wTScore_pvals.txt";

my @bases = ('(',')','B','I','H');

my $pointMutationDeltas = readPointMutationFile($pointMutationFile);
my($calculations,$maxTScore) = runCalculations($pointMutationDeltas,\@bases);

open(OPTF,">$tScoreOutputFile") or die "failed to open $tScoreOutputFile for writing\n";
print OPTF "#max=$maxTScore\n";
foreach $calcInfo (@{$calculations}) {
    my($pos,$bIdx,$cellMean,$cellVariance,$cellCount,$tScore,$pValue,$qValue) = @{$calcInfo};
    print OPTF "$pos\t$bIdx\t$cellMean\n";
}
close(OPTF);

open(OPTF,">$pValOutputFile") or die "failed to open $pValOutputFile for writing\n";
foreach $calcInfo (@{$calculations}) {
    my($pos,$bIdx,$cellMean,$cellVariance,$cellCount,$tScore,$pValue,$qValue) = @{$calcInfo};
    print OPTF "$pos\t$bIdx\t$pValue\t$qValue\n";
}
close(OPTF);


sub runCalculations {
    my($pointMutationDeltas,$bases) = @_;
    my @calculations;
    my $maxTScore = -~0;
    my $maxMean = -~0;
    my $gMean = 0;  #mean if the null hypothesis were true
    my $count = 0;
    foreach my $pos (sort {$a <=> $b} keys %{$pointMutationDeltas}) {
	$count += @{$bases};
    }
    print "cellCount = $count\n";
    foreach my $pos (sort {$a <=> $b} keys %{$pointMutationDeltas}) {
	for (my $bIdx = 0; $bIdx < @{$bases}; $bIdx++) {
	    my $base = $bases->[$bIdx];
	    my($cellMean,$cellVariance,$cellCount) = calculateMean($pointMutationDeltas->{$pos}{$base});
	    my $cellStd = sqrt($cellVariance);
	    my $denom = $cellStd / sqrt($cellCount);
	    my @before_values = ($gMean) x @{$pointMutationDeltas->{$pos}{$base}};
	    my $t_test = new Statistics::DependantTTest;
	    $t_test->load_data('before',@before_values);
	    $t_test->load_data('after',@{$pointMutationDeltas->{$pos}{$base}});
	    my ($tValue,$deg_freedom) = $t_test->perform_t_test('after','before');
	    my $cdf = Statistics::Distributions::tprob($deg_freedom,$tValue);
	    my $pValue = 2 * min($cdf,1-$cdf);
	    my $qValue = $pValue * $count; #Bonferroni corrected p-value
	    if (abs($tValue) > $maxTScore) {
		$maxTScore = abs($tValue);
	    }
	    if (abs($cellMean) > $maxMean) {
		$maxMean = abs($cellMean);
	    }
	    #print "$tValue\t$deg_freedom\t$pValue\n";
	    #my $tValue2 = ($cellMean - $gMean) / $denom;
	    #print "$tValue\t$tValue2\n";
	    push(@calculations,[$pos,$bIdx,$cellMean,$cellVariance,$cellCount,$tValue,$pValue,$qValue]);
	}
    }
    return(\@calculations,$maxMean);
}

sub runCalculations_old {
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
	my $testVarSum = 0;
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

sub min {
    my($x,$y) = @_;
    if($x < $y) {
        return $x;
    } else {
        return $y;
    }
}
