#!/usr/bin/perl -w
use strict;
 
my $filename = "KEGG_modules_downloaded_140114.txt";
open(IN, "<", $filename) or die "Can't read file $filename: $!\n";
 
my $m;
my %modules;

while (defined (my $line = <IN>)) {   
    chomp $line;
    my @tmp = split ' ', $line;
    if ($line =~ /^[D|E].*(M\d{5})\s+(.*)/){ # if line contains a module ID
        $m = "$1\t$2";
    } elsif ($tmp[1] =~ /K\d{5}/){ # if line contains a gene ID
        push(@{$modules{$m}}, $tmp[1]);
    } else {
        print "$line\n"
    }
}
close IN;
open(OUT, ">", "KEGG_modules.tab") or die "Cannot create file.\n";
 
foreach my $key (keys %modules){
    print OUT "$key\t";
    foreach my $gene (@{$modules{$key}}){
        print OUT "$gene;";
    }
    print OUT "\n";
}
