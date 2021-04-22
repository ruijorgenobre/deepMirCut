#!/usr/bin/perl -w
use strict;

my $usage = "USAGE:\n <mirdataset file> <hairpins fasta> <miRBase gff> <species prefix>\n";

my $mirDatasetFile = $ARGV[0] or die $usage;
my $hairpinFasta = $ARGV[1] or die $usage;
my $miRBaseGff = $ARGV[2] or die $usage;
my $queryPrefix = $ARGV[3] or die $usage;

my $outputFile = "$queryPrefix\_OneProduct.fa";

my $sequences = readFastaFile($hairpinFasta);
my $mirSpecies = getOneProdSpeciesPrefixes($mirDatasetFile);
my($hairpins,$products) = readMirbaseGff3_noRepeatPrecursors($miRBaseGff,$sequences);
my $gffOneProductSet = getGffOneProductMirs($hairpins,$products);
splitFastaBySpecies($sequences,$mirSpecies,$gffOneProductSet,$queryPrefix,$outputFile);

sub splitFastaBySpecies {
    my($sequences,$mirSpecies,$gffOneProductSet,$queryPrefix,$outputFile) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    my %collectedSeqs;
    foreach my $name (keys %{$sequences}) {
	if ($mirSpecies->{$name}) {
	    my $speciesPrefix = $mirSpecies->{$name};
	    if ($speciesPrefix eq $queryPrefix && $gffOneProductSet->{$name}) {
		$collectedSeqs{$name} = $sequences->{$name};
	    }
	}
    }
    foreach my $name (sort(keys %collectedSeqs)) {
	print OPTF ">$name\n$sequences->{$name}\n";
    }
    close(OPTF);
}

sub getGffOneProductMirs {
    my($hairpins,$products) = @_;
    my %gffOneProductSet;
    foreach my $chrom (keys %{$hairpins}) {
	foreach my $hairpinInfo (@{$hairpins->{$chrom}}) {
	    my($start,$stop,$strand,$id,$name) = @{$hairpinInfo};
	    if ($products->{$id}) {
		my $count = @{$products->{$id}};
		if ($count == 1) {
		    $gffOneProductSet{$name} = 1;
		}
	    }
	}
    }
    return \%gffOneProductSet;
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

sub getOneProdSpeciesPrefixes {
    my($mirDatasetFile) = @_;
    my %mirSpecies;

    open(MDF,$mirDatasetFile) or die "failed to open $mirDatasetFile\n";
    while (<MDF>) {
	chomp;
	unless ( /^\#/ ) {
	    my($speciesPrefix,$name,$id,$genomicLocation,$altGenomicLocation,$product5p,$product3p,
	       $drosha5p,$dicer5p,$dicer3p,$drosha3p,$hpStart,$hpStop,$leftBuffer,$rightBuffer,$seqTag,$sequence) = split(/\t/);
	    if (($drosha5p eq '-' || $dicer5p eq '-' || $dicer3p eq '-' || $drosha3p eq '-') &&
		not($drosha5p eq '-' && $dicer5p eq '-' && $dicer3p eq '-' && $drosha3p eq '-')) {
		$mirSpecies{$name} = $speciesPrefix;
	    }
	}
    }
    close(MDF);
    return \%mirSpecies;
}

sub readOrganismFile {
    my($organismFile) = @_;
    my %organismLineage;
    open(ORGF,$organismFile) or die "failed to open $organismFile\n";
    while (<ORGF>) {
	chomp;
	unless ( /^\#/ ) {
	    my($organism,$division,$name,$tree) = split(/\t/);
	    my @lineage = split(';',$tree); 
	    @{$organismLineage{$organism}} = @lineage;
	}
    }
    close(ORGF);
    return \%organismLineage;
}

sub readMirbaseGff3_noRepeatPrecursors {
    my($mirbaseGff3,$sequences) = @_;
    my %hairpins;
    my %products;
    my %seqHash;
    open(MBGFF3,$mirbaseGff3) or die "could not open $mirbaseGff3\n";
    my $lastHairpinId = 0;
    while(<MBGFF3>) {
	unless(/^\#/) {
	    chomp;
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$info) = split(/\t/);
	    my %info;
	    my @terms = split(/;/,$info);
	    foreach my $term (@terms) {
		my($key,$value) = $term =~ /(.*)=(.*)/;
		$info{$key} = $value;
	    }
	    if($type eq "miRNA_primary_transcript") {
		# hairpin region from the gff file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n"; 
		my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $seq = $sequences->{$name};
		if (!($id =~ /\_/ ) && !($seqHash{$seq})) {
		    push(@{$hairpins{$chrom}},[$start,$stop,$strand,$id,$name]);
		    $seqHash{$seq} = 1;
		} else {
		    print "excluding $id because it's a repeat\n";
		}
		$lastHairpinId = $id;
	    }
	    if($type eq "miRNA") {
		# mature product from the hairpin file
		my $id = $info{ID} or die "No ID found for the line:\n$_\n";
                my $name = $info{Name} or die "No Name found for the line:\n$_\n";
		my $parentId = $info{Derives_from} or die "No Derives_from found for the line:\n$_\n";
		if ($lastHairpinId eq $parentId) {
		    push(@{$products{$parentId}},[$chrom,$start,$stop,$strand,$id,$name]);
		} else {
		    print "excluding $id because it's belongs to $lastHairpinId\n"
		}
	    }
	}
    }
    close(MBGFF3);
    return(\%hairpins,\%products);
}
