#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n$0 <oneProductTestSet> <miRPreprocess hairpins> <miRPreprocess products> <min count>\n";

my $oneProductSetFile = $ARGV[0] or die $usage;
my $miRPreprocessHairpins = $ARGV[1] or die $usage;
my $miRPreprocessProducts = $ARGV[2] or die $usage;
my $minTotal = (@ARGV > 3) ? $ARGV[3] : 0;

my($prodSet5p,$prodSet3p) = readOneProductSetFile($oneProductSetFile);
my $sequences = loadPrecursorSequences($miRPreprocessHairpins);
my $products = loadProductFile($miRPreprocessProducts);

my $outputFile5p = "oneProd_test_5pAnnot_3pUnannot.txt";
my $outputFile3p = "oneProd_test_3pAnnot_5pUnannot.txt";
my $outputFasta5p = "oneProd_3pUnannot.fa";
my $outputFasta3p = "oneProd_5pUnannot.fa";

#checkProductsAgainstHPSeqeunces($products,$sequences);
#checkProductsAgainstTestSet($products,$prodSet5p,"5p");
#checkProductsAgainstTestSet($products,$prodSet3p,"3p");

my $newProdSet5p = addCutsToSet($products,$prodSet5p,"3p",$minTotal);
my $newProdSet3p = addCutsToSet($products,$prodSet3p,"5p",$minTotal);

printOutputFiles($newProdSet5p,"3p",$outputFile5p,$outputFasta5p);
printOutputFiles($newProdSet3p,"5p",$outputFile3p,$outputFasta3p);

sub printOutputFiles {
    my($prodSet,$newSide,$outputFile,$outputFasta) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    open(OFASTA,">$outputFasta") or die "failed to open $outputFasta for writing\n";
    foreach my $tag (keys %{$prodSet}) {
	my($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA) = @{$prodSet->{$tag}};
	print OPTF "$id\t$tag\t$mirId\t$prod5p\t$prod3p\t$DR5\t$DC5\t$DC3\t$DR3\t$hpStart\t$hpStop\t$sequence\t$fold\t$bpRNA\n";
	if ($newSide eq "5p") {
	    print OFASTA ">$tag-$newSide\n". getSeqBetweenCuts($sequence,$DR5,$DC5) . "\n";
	} elsif ($newSide eq "3p") {
	    print OFASTA ">$tag-$newSide\n". getSeqBetweenCuts($sequence,$DC3,$DR3) . "\n";
	} else {
	    die "side = $newSide not understood\n";
	}
    }
    close(OPTF);
    close(OFASTA);
}

sub addCutsToSet {
    my($products,$prodSet,$newSide,$minTotal) = @_;
    my %newSet;
    foreach my $tag (keys %{$prodSet}) {
	my($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA) = @{$prodSet->{$tag}};
	if ($products->{$tag}{$newSide}) {
	    my($prodStart,$prodStop,$prodSeq,$total) = @{$products->{$tag}{$newSide}};
	    if ($total >= $minTotal) {
		if ($newSide eq "5p") {
		    my($newDR5,$newDC5) = convertStartStopIntoCuts($hpStart,$prodStart,$prodStop);
		    @{$newSet{$tag}} = ($id,$mirId,$prod5p,$prod3p,$newDR5,$newDC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA);
		} elsif ($newSide eq "3p") {
		    my($newDC3,$newDR3) = convertStartStopIntoCuts($hpStart,$prodStart,$prodStop);
		    @{$newSet{$tag}} = ($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$newDC3,$newDR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA);
		} else {
		    die "side = $newSide not understood\n";
		}
	    }
	}
    }
    return \%newSet;
}

sub getSeqBetweenCuts {
    my($sequence,$cut1,$cut2) = @_;
    if (($cut1 =~ /-/) || ($cut2 =~ /-/)) {
	die "Cant get seqeunce. $cut1 or $cut2 are not numerical."
    }
    my($lcut1,$rcut1) = split(",",$cut1);
    my($lcut2,$rcut2) = split(",",$cut2);
    $rcut1 -= 1;
    $lcut2 -= 1;
    return substr($sequence,$rcut1,$lcut2-$rcut1+1);
}

sub convertStartStopIntoCuts {
    my($hpStart,$prodStart,$prodStop) = @_;
    my $cut1 = "" . ($hpStart+$prodStart-1) . "," . ($hpStart+$prodStart); 
    my $cut2 = "" . ($hpStart+$prodStop) . "," . ($hpStart+$prodStop+1);
    return($cut1,$cut2);
}

sub checkProductsAgainstTestSet {
    my($products,$prodSet,$side) = @_;
    foreach my $tag (keys %{$prodSet}) {
	my($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA) = @{$prodSet->{$tag}};
	if ($products->{$tag}{$side}) {
	    my($prodStart,$prodStop,$prodSeq,$total) = @{$products->{$tag}{$side}};
	    if ($total > 1) {
		if ($side eq "5p") {
		    my $testSeq = getSeqBetweenCuts($sequence,$DR5,$DC5);
		    my($newDR5,$newDC5) = convertStartStopIntoCuts($hpStart,$prodStart,$prodStop);
		    my $testSeq2 = getSeqBetweenCuts($sequence,$newDR5,$newDC5);
		    print "$DR5\t$newDR5\t$DC5\t$newDC5\n";
		    print "$tag $side $total\n$sequence\n$testSeq\n$testSeq2\n$prodSeq\n\n";
		} elsif ($side eq "3p") {
		    my $testSeq = getSeqBetweenCuts($sequence,$DC3,$DR3);
		    my($newDC3,$newDR3) = convertStartStopIntoCuts($hpStart,$prodStart,$prodStop);
		    my $testSeq2 = getSeqBetweenCuts($sequence,$newDC3,$newDR3);
		    print "$DC3\t$newDC3\t$DR3\t$newDR3\n";
		    print "$tag $side $total\n$sequence\n$testSeq\n$testSeq2\n$prodSeq\n\n";
		} else {
		    die "side = $side not understood\n";
		}
	    }
	} else {
	    print "$tag $side not in miRPreprocessed products\n";
	}
    }
}

sub checkProductsAgainstHPSeqeunces {
    my($products,$sequences) = @_;
    foreach my $tag (keys %{$products}) {
	foreach my $side (keys %{$products->{$tag}}) {
	    my($prodStart,$prodStop,$prodSeq,$total) = @{$products->{$tag}{$side}};
	    my $testSeq = substr($sequences->{$tag},$prodStart,$prodStop-$prodStart+1);
	    if (length($prodSeq) < $prodStop-$prodStart+1) {
		print "Warning: $tag $side goes off end of sequence in product file.\n";
	    }
	    if ($prodSeq ne $testSeq) {
		print "Sequences for $tag $side between hairpins and products file are different:\n$prodSeq\n$testSeq\n\n";
	    }
	}
    }

}

sub readOneProductSetFile {
    my($oneProductSetFile) = @_;
    my %prodSet5p;
    my %prodSet3p;
    open(OPSF,$oneProductSetFile) or die "failed to open $oneProductSetFile\n";
    while (<OPSF>) {
	chomp;
	my($id,$name,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA) = split(/\t/);
	if ($prod5p ne '-') {
	    if ($prodSet5p{$name} || $prodSet3p{$name}) {
		die "$name appears more than once in dataset\n";
	    }
	    @{$prodSet5p{$name}} = ($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA);
	} elsif ($prod3p ne '-') {
	    if ($prodSet3p{$name} || $prodSet5p{$name}) {
		die "$name appears more than once in dataset\n";
	    }
	    @{$prodSet3p{$name}} = ($id,$mirId,$prod5p,$prod3p,$DR5,$DC5,$DC3,$DR3,$hpStart,$hpStop,$sequence,$fold,$bpRNA);
	} else {
	    die "There are no products for $name\n";
	}
    }
    close(OPSF);
    return(\%prodSet5p,\%prodSet3p);
}

sub loadPrecursorSequences {
    my($miRPreprocessHairpins) = @_;
    my %sequences;
    open(MPH,$miRPreprocessHairpins) or die "failed to open $miRPreprocessHairpins\n";
    while (<MPH>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$chrom,$start,$stop,$strand,$leftCenter,$rightCenter,$totalSense,$totalAntisense,$mfe,$sequence) = split(/\t/);
	    $sequences{$tag} = $sequence;
	}
    }
    close(MPH);
    return \%sequences;
}


sub loadProductFile {
    my($miRPreprocessProducts) = @_;
    my %products;
    open(MPP,$miRPreprocessProducts) or die "failed to open $miRPreprocessProducts\n";
    while (<MPP>) {
	chomp;
	unless ( /^#/ ) {
	    my($tag,$side,$type,$total,$totalMostAbundant,$adjusted,$adjTotalMostAbundant,$start,$stop,$strand,$sequence) = split(/\t/);
	    if ($type eq "miR") {
		@{$products{$tag}{$side}} = ($start,$stop,$sequence,$total);
	    } 
	}
    }
    close(MPP);
    return \%products;
}
