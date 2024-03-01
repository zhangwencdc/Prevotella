#!/usr/bin/perl
use strict;
use warnings;

my @file=glob "/home/zhangwen/project/2022Time/Prevotella/Fq/*_genomic.fna.fq1.filter.fq";
my %name;
foreach my $fq1 (@file) {
		my @a=split"_",$fq1;
		my $name=substr($a[0],46);
		print "$name";
		if(exists $name{$name}){next;}
		$name{$name}++;
		system "cat /home/zhangwen/project/2022Time/Prevotella/Fq/$name*_genomic.fna.fq1.filter.fq >tmp.fq1.fq\n";
		system "seqkit rmdup tmp.fq1.fq -o $name.Prevotella.fq1.fastq\n";
		system "cat /home/zhangwen/project/2022Time/Prevotella/Fq/$name*_genomic.fna.fq2.filter.fq >tmp.fq2.fq\n";
		system "seqkit rmdup tmp.fq2.fq -o $name.Prevotella.fq2.fastq\n";
		system "rm -rf tmp.fq1.fq tmp.fq2.fq\n";
}
