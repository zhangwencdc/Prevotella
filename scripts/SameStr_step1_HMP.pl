#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#conda activate biobakery

my $raw="/home/zhangwen/Data/HMP_Data/Fastq/";
my @file=glob "$raw/*.denovo_duplicates_marked.trimmed.1.fastq";
system "mkdir sams bowtie2 profiles consensus_markers";
foreach my $file (@file) {
        #print "$file\n";

        #print "$name\n";
        my $f=basename($file);
                #$f=substr($f,0,length($f)-11);
                #my $f2=substr($file,0,length($file)-11)."_2.clean.fq";
                #my @name=split"_",$f;
                my $name=substr($f,0,length($f)-41);
                my $fq=$name.".fq";
		system "metaphlan $file --bowtie2db /home/zhangwen/Data/Metaphlan/ --input_type fastq --nproc 30 --legacy-output -t rel_ab --bowtie2out $name.bowtieout --samout $name.sam.bz2 -o $name.profile.txt\n";
 
}
