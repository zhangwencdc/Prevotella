#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#Time_All_sstr_data.tsv
my $id=$ARGV[1];

my %people;
open(F,$id);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];
}
close F;

open(OUT,">$file.people");
open(FILE,$file);my %sp;my %dp;my %species;my %compare_same;my %compare_diff;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
	unless(substr($a[0],0,1) eq "S"){$a[0]="S".$a[0];}
	unless(substr($a[1],0,1) eq "S"){$a[1]="S".$a[1];}
	print OUT "$l\t";
	print OUT "$people{$a[0]}\t$people{$a[1]}\t";$species{$a[2]}++;
	if($l=~/shared_strain/){
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP";$sp{$a[2]}++;$compare_same{$a[2]}++;}else{print OUT "DiffP";$dp{$a[2]}++;$compare_diff{$a[2]}++;}
	}else{
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP";$compare_same{$a[2]}++;}else{print OUT "DiffP";$compare_diff{$a[2]}++;}
	}
	print OUT "\n";
}
close FILE;
open(O,">$file.people.stat");
my @species=sort keys %species;
print O "Species,Compare within SampePeople,Share within Samepeople,Compare Between Different People,Share Between Different People\n";
foreach my $species (@species) {
	print O "$species,$compare_same{$species},$sp{$species},$compare_diff{$species},$dp{$species}\n";
}