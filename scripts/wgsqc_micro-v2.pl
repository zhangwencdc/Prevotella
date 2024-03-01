#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;  
use Cwd;
my $path = getcwd;
#conda activate biobakery
##Version2:1、增加了Bin分箱功能;2、增加了Metaphlan功能，支持后续Strainphlain和SameStr分析；

my ($Keyname,$Ref,$input,$database,$type,$fa,$r,$level,$Outdir,$verbose,$db,$snp,$fq1,$fq2,$assemble,$host,$dores,$dovf,$virus);
my ($Verbose,$Help);


GetOptions(
        "R:s"=>\$Ref, #参考基因组
        "L|level:s"=>\$level,  #默认Level: S  ,S代表Species，也可以设置为G，代表Genus
        "Out|O:s"=>\$Outdir,
        "verbose"=>\$Verbose,
        "K|Tag|Key:s" =>\$Keyname,
		"r|rlen:s" =>\$r,
		"database|DB:s"=>\$db,###比对数据库
		"1|Fq1:s"=>\$fq1,###必选 fq1
		"2|Fq2:s"=>\$fq2,###fq2
		"Res|res:s"=>\$dores,####查找耐药基因
		"VF|vf:s"=>\$dovf,####查找毒力基因
		"Assemble"=>\$assemble,###组装
		"H|Host|host"=>\$host,###鉴定Host
		"V|Virus|virus"=>\$virus,###鉴定virus
        "help"=>\$Help
);
####
$Keyname ||= "Input";
$Outdir ||= ".";
$type ||=1;
$level ||="S";
$r ||=150;

die `pod2text $0` if ($Help);

if(substr($fq1,length($fq1)-3,3)=~/.gz/){system "gunzip $fq1\n";$fq1=substr($fq1,0,length($fq1)-3);}
if(defined $fq2){if(substr($fq2,length($fq2)-3,3)=~/.gz/){system "gunzip $fq2\n";$fq2=substr($fq2,0,length($fq2)-3);}}
###FastV
if(defined $virus){
	system "mkdir Virus\n";
	print "Fastv searching\n";
	if(defined $fq2){
		system "fastv -i $fq1 -I $fq2 -o Fastv/$Keyname.fq1.fastv.fq -O Virus/$Keyname.fq2.fastv.fq -c /home/zhangwen/Data/Fastv/viral.kc.fasta -h Virus/$Keyname.fastv.html -j Virus/$Keyname.fastv.json\n";
	}else{
		system "fastv -i $fq1 -o Virus/$Keyname.fq1.fastv.fq -c /home/zhangwen/Data/Fastv/viral.kc.fasta -h Virus/$Keyname.fastv.html -j Virus/$Keyname.fastv.json\n";
	}
}

###Kraken2+Bracken
print "Kraken2 searching\n";
system "mkdir Kraken2\n";
if(defined $fq2){
	system "kraken2 --db /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ --report Kraken2/$Keyname.report --output Kraken2/$Keyname.kraken --paired $fq1 $fq2\n";
	system "/home/zhangwen/bin/Bracken/bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i Kraken2/$Keyname.report -o Kraken2/$Keyname.report.bracken  -l $level -r $r\n";
	system "perl /home/zhangwen/bin/kraken2-translate.pl Kraken2/$Keyname.report >Kraken2/$Keyname.report.txt\n";
	system "ktImportText Kraken2/$Keyname.report.txt -o Kraken2/$Keyname.report.html";
}else{
	system "kraken2 --db /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE/ --report Kraken2/$Keyname.report --output Kraken2/$Keyname.kraken $fq1 \n";
	system "/home/zhangwen/bin/Bracken/bracken -d /home/zhangwen/Data/Kraken2/minikraken2_v2_8GB_201904_UPDATE -i Kraken2/$Keyname.report -o Kraken2/$Keyname.report.bracken  -l $level -r $r\n";
	system "perl /home/zhangwen/bin/kraken2-translate.pl Kraken2/$Keyname.report >Kraken2/$Keyname.report.txt\n";
	system "ktImportText Kraken2/$Keyname.report.txt -o Kraken2/$Keyname.report.html";
}
##Metaphlain
system "mkdir sams bowtie2 profiles consensus_markers\n";
  
system " metaphlan $fq1,$fq2 --input fastq -s sams/$Keyname.sam.bz2 --bowtie2out bowtie2/$Keyname.bowtie2.bz2 -o profiles/$Keyname.profile.tsv\n";
                system "sample2markers.py -i sams/$Keyname.sam.bz2 -o consensus_markers -n 8\n";

###VFF

if(defined $dovf){
	system "mkdir VFF_ARGs\n";
print "Virulent factor searching\n";
if(defined $fq2){
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro.pl -1 $fq1 -2 $fq2 -T /home/zhangwen/Data/VFF/VFDB_setB_nt.fas -o VFF_ARGs/$Keyname.vff.csv\n ";
}else{
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro.pl -1 $fq1 -T /home/zhangwen/Data/VFF/VFDB_setB_nt.fas -o VFF_ARGs/$Keyname.vff.csv\n ";
}
}
##Resistance
if(defined $dores){
	system "mkdir VFF_ARGs\n";
print "Resistance genes searching\n";
if(defined $fq2){
	system "perl /home/zhangwen/bin/Resistance/Raw_Target_Resistance-micro.pl -1 $fq1 -2 $fq2 -Key VFF_ARGs/$Keyname.res \n ";
}else{
	system "perl /home/zhangwen/bin/Resistance/Raw_Target_Resistance-micro.pl -1 $fq1 -Key VFF_ARGs/$Keyname.res \n ";
}
}


if(defined $Ref){

	print "Target searching\n";
	if(defined $fq2){
	system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro.pl -1 $fq1 -2 $fq2 -T $Ref -o VFF_ARGs/$Keyname.target.csv\n ";
	}else{
		system "perl /home/zhangwen/bin/Target_bowtie-dell-Micro.pl -1 $fq1 -T $Ref -o VFF_ARGs/$Keyname.target.csv\n ";
	}
}
###host
if(defined $host){
		system "mkdir Host\n";
	if(defined $fq2){
	system "perl /home/zhangwen/bin/Host.pl -1 $fq1 -2 $fq2 -o Host/$Keyname.host\n";
	}else{
		system "perl /home/zhangwen/bin/Host.pl -1 $fq1  -o Host/$Keyname.host\n";
	}
}
####
if(defined $assemble){
	system "mkdir Assembled\n";
	
	if(defined $fq2){
		
		system "spades.py -1 $fq1 -2 $fq2 -o Assembled/$Keyname --meta \n";
	}else{
		system "spades.py -1 $fq1 -o Assembled/$Keyname --meta\n";
	}

	my $genome="Assembled/".$Keyname.".assembled.fasta";
	system "cp Assembled/$Keyname/scaffolds.fasta $genome\n";
	system "bowtie2-build  $genome $genome\n";
   system "bowtie2  -x $genome -1 $fq1 -2 $fq2 | samtools sort --threads 30 -o Assembled/$Keyname.sort.bam - \n";
   system "jgi_summarize_bam_contig_depths --outputDepth Assembled/$Keyname.depth.txt Assembled/$Keyname.sort.bam\n";
   system "metabat2 -i $genome -a Assembled/$Keyname.depth.txt -o Assembled/$Keyname\n";
}

#system "multiqc $file -c /home/zhangwen/bin/multiqc_config.yaml\n";
print "Success Finished\n";
###############################
sub norRNA{
	my $fq=shift @_;
	my $kraken=shift @_;
	open(FILE,$kraken);
	my %rrna;
	while(1){
		my $line=<FILE>;
		unless ($line) {last;
		}
		chomp $line;
		my @a=split"\t",$line;
		if($a[0]eq "C"  ){$rrna{$a[1]}++;}
	}
	close FILE;
	my $result=$fq.".norRNA";
	open(OUT,">$result");
	open(FQ,$fq);
	while(1){
			my $line=<FQ>;
			unless($line){last;}
			chomp $line;
			my $seq=<FQ>;
			my $a=<FQ>;
			my $b=<FQ>;
			
			my @name=split" ",$line;
			my $name=substr($name[0],1);
			if(exists $rrna{$name}){next;}
			
			print OUT "$line\n$seq$a$b";
	}
	close FQ;
	return $result;
}

sub nohuman{
	my $fq=shift @_;
	my $kraken=shift @_;
	open(FILE,$kraken);
	my %human;
	while(1){
		my $line=<FILE>;
		unless ($line) {last;
		}
		chomp $line;
		my @a=split"\t",$line;
		if($a[2]eq "9605" || $a[2]eq "9606" ){$human{$a[1]}++;}
	}
	close FILE;
	my $result=$fq.".nohuman";
	open(OUT,">$result");
	open(FQ,$fq);
	while(1){
			my $line=<FQ>;
			unless($line){last;}
			chomp $line;
			my $seq=<FQ>;
			my $a=<FQ>;
			my $b=<FQ>;
			
			my @name=split" ",$line;
			my $name=substr($name[0],1);
			if(exists $human{$name}){next;}
			
			print OUT "$line\n$seq$a$b";
	}
	close FQ;
	return $result;
}