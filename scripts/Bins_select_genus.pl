#!/usr/bin/perl
use strict;
use warnings;
#Usage:基于gtdbtk预测结果提取指定菌种的基因组

my $file=$ARGV[0];#gtdbtk.bac120.summary.tsv gtdbtk预测结果

my $target=$ARGV[1];#Prevotella

#my @a=split"_",$target;
my $species="g__".$target;

my $dir=$ARGV[2];##Bins存放路径


my %num;
open(F,$file);
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	unless($a[1]=~/$species/){next;}
	my @b = $a[0] =~ /^(.+)\.(\d+)\.(\w+)$/;
	my @b = split(/\./, $a[0]);
	#print "$b[0]\n";
	$num{$b[0]}.=",".$a[0].".fna";
}
close F;

my @sample=sort keys %num;

foreach my $sp (@sample) {
	my $b=$num{$sp};
	my @b=split",",$b;
	my $out=$target."_".$sp.".fasta";
	foreach my $b (@b) {
		if($b=~/[0-9a-zA-Z]/){
			my $file=$dir."/".$b;
			system "cat $file >>$out\n";
		}
	}
}