#!/usr/bin/perl
use strict;
#use warnings;

my $bloomr="D:\\E盘\\Project\\科室成员志愿者\\Ecoli_paper\\Ecoli\\StrainPhlan\\Bloom\\Bloom.r"; #画图R语言路径

###find Bloom点
my $file=$ARGV[0]; ####Ecoli.percentage.txt


open(G,$file);
my %group;my %g;my %time;my %sum;my %avg;
while(1){
	my $line=<G>;
	unless($line){last;}
	chomp $line;
	if($line=~/^Sample/){next;}
	my @a=split"\t",$line;
	$group{$a[0]}=$a[1];$g{$a[1]}++;$time{$a[0]}=$a[2];$sum{$a[1]}+=$a[4];
}
close G;


open(FILE,$file);
my $line=<FILE>;
chomp $line;
open(OUT,">$file.bloom");
print OUT "SampleID,People,Taxon,Percentage,Avg for people\n";
my %bloom;my %taxon;my %previous;my %disappear;
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	print OUT "$line\t";
	my $avg=$sum{$a[1]}/$g{$a[1]};
	if($a[4]>=$avg*5){print OUT "Bloom";$bloom{$a[1]}++;}   #大于平均丰度5倍，被定义成bloom
	if($a[4]<$previous{$a[1]}*0.2 && $a[4]>0){print OUT "Disappear";$disappear{$a[1]}++;}#低于前一个时间点丰度1/5倍，被定义成Disappear
	$previous{$a[1]}=$a[4];
	print OUT "\n";
}
close OUT;
close FILE;

my @bloom=sort keys %bloom;
foreach my $bloom (@bloom) {
	print "$bloom,$bloom{$bloom},bloom\n";
}
my @disappear=sort keys %disappear;
foreach my $disappear (@disappear) {
	print "$disappear,$disappear{$disappear},disappear\n";
}
#my @group=sort keys %g;
#open OF, ">$file.all.filter";
#	my $id=0;
#	print OF "ID,SampleID,People,Time,Percentage\n";
#foreach my $group (@group) {
##	my $id=0;
##	open OF, ">$file.$group.filter";
##	print OF "ID,SampleID,People,Time,Percentage\n";
#	print  "$file.$group.filter\n";
#	open(FILE,$file);
#	my $line=<FILE>;
#	chomp $line;
#	while(1){
#	my $line=<FILE>;
#	unless($line){last;}
#	chomp $line;
#	my @a=split"\t",$line;
#	unless($a[1] eq $group){next;}
#	$id++;
#	my $time;
#	if(length($a[2])==2){$time=substr($a[2],0,1)."0".substr($a[2],1,1);}else{$time=$a[2];}
#		print OF "$id,$a[4],$a[1],$time,$a[4]\n";
#	}
#	system "Rscript $bloomr -i $file.$group.filter -o $file.$group.filter\n";
#}