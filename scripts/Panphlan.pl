#!/usr/bin/perl
use strict;
use warnings;
#conda activate biobakery
my $id=$ARGV[0];#All_sampleid.txt
my $p=$ARGV[1]."/".$ARGV[1];#Prevotella_copri
open(F,$id);my %fq1;my %fq2;
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	$fq1{$a[0]}=$a[1];
	$fq2{$a[0]}=$a[2];

}
close F;

my @id=sort keys %fq1;
my $pf=$p."_pangenome.tsv";
my $outdir=$ARGV[1]."map_results";
system "mkdir $outdir\n";
foreach my $id (@id) {
	my $fq1=$fq1{$id};
	my $fq2=$fq2{$id};
	my $o=$ARGV[1]."map_results/".$id."_erectale.csv";
	system "panphlan_map.py -p $pf --indexes $p -i $fq1 -o $o\n";

}

system "panphlan_profiling.py -i $outdir/  --o_matrix $ARGV[1].result_profile_erectale.tsv -p $pf --add_ref --o_covplot erectale_covplot\n";