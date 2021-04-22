#!/bin/perl -w 

$usage = "Usage:\n$0 <input fasta> <output fasta>\n";
$file = $ARGV[0] or die $usage;
$outFile = $ARGV[1] or die $usage;

open(F,$file);
open(OF,">$outFile");
while(<F>) {
    chomp;
    if(/>/) {
        print OF "$_\n";
    } else {
        s/U/T/g;
        print OF "$_\n";
    }
}
close(F);
close(OF)
