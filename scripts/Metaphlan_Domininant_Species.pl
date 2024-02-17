#!/usr/bin/perl
use strict;
#use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#Usage：合并Metaphlan结果，基于指定的菌属，寻找优势菌种
my ($input,$outdir,$genus,$species,$tag,$HELP,$info);
GetOptions(
                "I:s"=>\$input,#WGS_merge.metaphlan:Metaphlan合并结果 merge_metaphlan_tables.py *.profile.tsv >WGS_merge.metaphlan
                "O:s"=>\$outdir,
                "G:s"=>\$genus,
              
                "Tag:s"=>\$tag,
				 "D|Infor:s"=>\$info,#样本 编号表
                "help"=>\$HELP
);
die `pod2text $0` if ($HELP );
die `pod2text $0` if (!$genus && !$species );
die `pod2text $0` if (!$info );
if(!$outdir){$outdir="./";}
if(!$input){$input="./";}
if(!$tag){$tag="Test";}



open(F,$info);
my %people;my %time;my %p;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];$p{$a[1]}++;
	$time{$a[0]}=$a[2];
}
close F;


if(defined $genus){
	open(OUT,">$outdir/$tag.$genus.dominant");
	print OUT "Sample\tPeople\tTime\tPercentage\tDorminant Species\n";
	
	open(F,$input);
	my $header=<F>;
	chomp $header;
	my $header=<F>;
	chomp $header;
	my @header=split"\t",$header;my $type=0; my %dormi;my %g;
	while(1){
			my $l=<F>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			my @b= split(/\|/, $a[0]);
			my $b=pop @b;#print "$b\n";
			if($b=~/g__/){
				if($b=~/$genus/){
					$type=1;
					my $n=@a;
					foreach  (2..($n-1)) {
						my $sample=$header[$_];
						$sample=substr($sample,0,length($sample)-8);#剔除文件末尾的.profile
						#print "$_,$header[$_],$sample\n$header\n";
						$g{$sample}=$a[$_];
						#print OUT "$sample\t$people{$sample}\t$time{$sample}\t$a[$_]\n";
					}
				
				}else{
					$type=0;
				}
			}else{
				if($b=~/s__/ && $b=~/$genus/){
					if($type==1){
						my $n=@a;
						foreach  (2..($n-1)) {
							my $sample=$header[$_];
							$sample=substr($sample,0,length($sample)-8);#剔除文件末尾的.profile
							#print O "$sample\t$b\t$people{$sample}\t$time{$sample}\t$a[$_]\n";
							if($a[$_]>0.8*$g{$sample}){$dormi{$sample}=$b;}
						}
					}
				}
			}
	}
	close F;
	
	
	my @sample=sort keys %g;
	foreach my $sample (@sample) {
		print OUT "$sample\t$people{$sample}\t$time{$sample}\t$g{$sample}\t$dormi{$sample}\n";
	}
}