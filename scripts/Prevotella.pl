#!/usr/bin/perl
use strict;
use warnings;

my $genomelist="/home/zhangwen/Data/GTDB/genome.list.detail";

my $genus=$ARGV[0];

my $filelist="/home/zhangwen/project/2022Time/WGS_Data/Time_New/list";

#列出复核条件的基因组

open(G,$genomelist);
while(1){
	my $l=<G>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my @b=split"_",$a[2];
	unless($b[0] eq $genus){next;}
	if($b[1] eq "sp."){next;}
	my @c=split"/",$a[0];
	my $name=pop @c;
	print "cp $a[0] $name\n";
	open(F,$filelist);
	while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		my @d=split"\t",$line;
		my $n=$d[0]."_".$name;
		print "perl /home/zhangwen/bin/Target_bowtie-dell-Micro-v2.pl -1 $d[1] -2 $d[2] -o $n -T $name\n";
	}
	close F;
}
close G;