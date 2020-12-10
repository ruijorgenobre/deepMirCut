#!/usr/bin/perl -w
use Graph;
do "./Scripts/bpRNA.pl";
use strict;

my $usage = "$0 <input fasta> <output file>\n";

my $inputFasta = $ARGV[0] or die $usage;
my $outputFile = $ARGV[1] or die $usage;
my $foldFile = $inputFasta . ".folds";


my $out = `RNAfold -i $inputFasta --noLP --noPS > $foldFile`;
print $out;
my($names,$folds,$seqs,$ids) = readFoldFile($foldFile);
my $context = getContextFromBPRNA($folds,$seqs);
printOutputFile($outputFile,$ids,$names,$folds,$seqs);

sub printOutputFile {
    my($outputFile,$ids,$names,$folds,$seqs) = @_;
    open(OPTF,">$outputFile") or die "failed to open $outputFile for writing\n";
    for my $id (@{$ids}) {
	print OPTF "$id\t$names->{$id}\t$seqs->{$id}\t$folds->{$id}\t$context->{$id}\n";
    }
    close(OPTF);
}

sub readFoldFile {
    my($foldFile) = @_;
    my $id_delim = "ex";
    my @newIds;
    my %names;
    my %folds;
    my %seqs;
    my $counter = 0;
    open(FF,$foldFile) or die "failed to open $foldFile\n";
    while (<FF>) {
	chomp;
	$_ =~ s/>//;
	my($id) = split;
	$counter++;
	my $new_id = "$id_delim$counter";
	push(@newIds,$new_id);
	$names{$new_id} = $id;
	my $sequence = <FF>;
	chomp($sequence);
	$seqs{$new_id} = $sequence;
	my $foldLine = <FF>;
	chomp($foldLine);
	my($fold,$mfe) = split(/\s+/,$foldLine);
	$folds{$new_id} = $fold;
    }
    close(FF);
    return(\%names,\%folds,\%seqs,\@newIds);
}

sub readInputFasta {
    my($inputFasta) = @_;
    open(INFA,$inputFasta) or die "failed to open $inputFasta\n";
    my %sequences;
    my $currTag;
    while (<INFA>) {
	chomp;
	if ( /^>/ ) {
	    ($currTag) = split;
	    $currTag =~ s/^>//;
	    $sequences{$currTag} = "";
	} else {
	    $sequences{$currTag} .= $_;
	}
    }
    close(INFA);
    return \%sequences;
}

sub getContextFromBPRNA {
    my($folds,$seqs) = @_;
    my %context;
    foreach my $id (sort keys %{$folds}) {
	my $sequence = $seqs->{$id};
	my $fold = $folds->{$id};
	$context{$id} = dotBracketToStructureArray($sequence,$fold);
    }
    return \%context;
}
