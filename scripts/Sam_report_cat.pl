#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

my $length = $ARGV[0]; # Cutoff 10000
my $align = $ARGV[1]; # 80%
my $out = $ARGV[2];
my @file = glob("*.report");
open(OUT, ">$out");
my %pos;
my %anno;

foreach my $file (@file) {
    my $name = basename($file);
    my @name = split("\\.", $name);
    my $sample = $name[0];
    open(F, $file);
    while (my $l = <F>) {
        chomp $l;
        my @a = split(",", $l);
        if ($a[3] >= $length && $a[4] >= $align) {
            print OUT "$name,$l\n";
            $pos{$sample}{$a[0]} = $a[4];
            $anno{$a[0]} = $a[6];
        }
    }
    close F;
}

my @pos = sort keys %pos;
my @contig = sort keys %anno;
print "ID,";
foreach my $pos (@pos) {
    print "$pos,";
}
print "Positive sample,Anno\n";
foreach my $contig (@contig) {
    print "$contig,";
    my $sum = 0;
    foreach my $pos (@pos) {
        if (exists $pos{$pos}{$contig}) {
            print "$pos{$pos}{$contig},";
            $sum++;
        } else {
            print "0,";
        }
    }
    print "$sum,$anno{$contig}\n";
}
