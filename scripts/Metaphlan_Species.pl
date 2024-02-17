#!/usr/bin/perl
use strict;
#use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#Usage：合并Metaphlan结果，基于指定的菌种，生成趋势图
my ($input,$outdir,$genus,$species,$tag,$HELP,$info);
GetOptions(
                "I:s"=>\$input,#WGS_merge.metaphlan:Metaphlan合并结果 merge_metaphlan_tables.py *.profile.tsv >WGS_merge.metaphlan
                "O:s"=>\$outdir,
                "G:s"=>\$genus,
                "S:s"=>\$species,
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
	open(OUT,">$outdir/$tag.$genus.metaphlancat");
	print OUT "Sample\tPeople\tTime\tPercentage\n";
	open(O,">$outdir/$tag.$genus.species");
	print O "Sample\tSpecies\tPeople\tTime\tPercentage\n";my %sp;
	open(F,$input);
	my $header=<F>;
	chomp $header;
	my $header=<F>;
	chomp $header;
	my @header=split"\t",$header;my $type=0;
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
						print OUT "$sample\t$people{$sample}\t$time{$sample}\t$a[$_]\n";
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
							print O "$sample\t$b\t$people{$sample}\t$time{$sample}\t$a[$_]\n";
							if($a[$_]>0.01){$sp{$b}++;}
						}
					}
				}
			}
	}
	close F;
	
	
	system "Rscript $Bin/Kraken_people.r -i $outdir/$tag.$genus.metaphlancat -o $outdir/$tag.$genus.genus\n";
	close OUT;
	my @sp=sort keys %sp;
	foreach my $sp (@sp) {
		if($sp{$sp}>=10){
			open(F,"$outdir/$tag.$genus.species");
			open(OT,">$outdir/tmp");
			print OT "Sample\tSpecies\tPeople\tTime\tPercentage\n";
			while(1){
			my $l=<F>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			unless($a[1] eq $sp){next;}
			print OT "$l\n";
			}
			close F;
			system "Rscript $Bin/Kraken_people.r -i $outdir/tmp -o $outdir/$tag.$sp.species\n";
			system "rm -rf $outdir/tmp\n";close OT;
		}
	}
	my @people=sort keys %p;
	foreach my $p (@people) {
			unless($p{$p}>=3){next;}
			open(F,"$outdir/$tag.$genus.species");
			open(OT,">$outdir/tmp");
			print OT "ID\tSample\tSpecies\tPeople\tTime\tPercentage\n";my $tid=0;
			while(1){
			my $l=<F>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			unless($a[2] eq $p){next;}
			$tid++;
			print OT "$tid\t$l\n";
			}
			close F;
			system "Rscript $Bin/Kraken_species.r -i $outdir/tmp -o $outdir/$tag.$p\n";
			close OT;
	}
}elsif(defined $species){
	open(OUT,">$outdir/$tag.$species.metaphlancat");
	print OUT "Sample\tPeople\tTime\tPercentage\n";
	open(F,$input);
	my $header=<F>;
	chomp $header;
	my @header=split"\t",$header;my $type=0;
	while(1){
			my $l=<F>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			my @b=split"|",$a[0];
			my $b=pop @b;#print "$b\n";
			if($b=~/s__/){
				if($b=~/$species/){
					$type=1;
					my $n=@a;
					foreach  (2..($n-1)) {
						my $sample=$header[$_];
						$sample=substr($sample,0,length($sample)-8);#剔除文件末尾的.profile
						print OUT "$sample\t$people{$sample}\t$time{$sample}\t$a[$_]\n";
					}
				
				}else{
					$type=0;
				}
			}
	}
	close F;
	close OUT;
	system "Rscript $Bin/Kraken_people.r -i $outdir/$tag.$species.metaphlancat -o $outdir/$tag.$species.species\n";
}

