#!/usr/bin/perl
use strict;
#use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Data::Dumper;
use File::Path;  
use Cwd;
my $path = getcwd;
#docker run -it wgs:v1
##

my ($Keyname,$Ref,$input,$database,$type,$fa,$r,$level,$Outdir,$verbose,$db,$snp,$fq1,$fq2,$assemble,$host,$dores,$dovf,$virus,$nohuman);
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
		"NH|nohuman"=>\$nohuman,#筛除human序列
        "help"=>\$Help
);
####
unless(-e $Outdir){system "mkdir $Outdir\n";}
$Keyname ||= "Input";
$Outdir ||= ".";
$type ||=1;
$level ||="S";
$r ||=150;

die `pod2text $0` if ($Help);
die `pod2text $0` if (!defined $fq1);
die `pod2text $0` if (!defined $db);

if(substr($fq1,length($fq1)-3,3)=~/.gz/){system "gunzip $fq1\n";$fq1=substr($fq1,0,length($fq1)-3);}
if(defined $fq2){if(substr($fq2,length($fq2)-3,3)=~/.gz/){system "gunzip $fq2\n";$fq2=substr($fq2,0,length($fq2)-3);}}

unless(-e "$Outdir/01.Stat"){system "mkdir $Outdir/01.Stat\n";}
unless(-e "$Outdir/WGS_Report"){system "mkdir $Outdir/WGS_Report\n";}
print "Stat running\n";
if(defined $fq2){
	system "seqkit stats -a $fq1 $fq2 >$Outdir/01.Stat/$Keyname.Stat_Report.txt\n";
}else{
	system "seqkit stats -a $fq1 >$Outdir/01.Stat/$Keyname.Stat_Report.txt\n";
}
system "cp $Outdir/01.Stat/$Keyname.Stat_Report.txt $Outdir/WGS_Report/$Keyname.Stat_Report.txt\n";

unless(-e "$Outdir/02.Structure"){system "mkdir $Outdir/02.Structure\n";}
#system "source /root/miniconda3/etc/profile.d/conda.sh\n";
#system "conda activate kraken2\n";
###FastV
if(defined $virus){
	
	system "mkdir $Outdir/02.Structure/Virus\n";
	print "Fastv searching\n";
	if(defined $fq2){
		system "/root/miniconda3/envs/kraken/bin/fastv -i $fq1 -I $fq2 -o $Outdir/02.Structure/Virus/$Keyname.fq1.fastv.fq -O $Outdir/02.Structure/Virus/$Keyname.fq2.fastv.fq -c $db/Fastv/viral.kc.fasta -h $Outdir/02.Structure/Virus/$Keyname.fastv.html -j $Outdir/02.Structure/Virus/$Keyname.fastv.json\n";
		system "rm -rf $Outdir/02.Structure/Virus/$Keyname.fq1.fastv.fq $Outdir/02.Structure/Virus/$Keyname.fq2.fastv.fq\n";
	}else{
		system "/root/miniconda3/envs/kraken/bin/fastv -i $fq1 -o $Outdir/02.Structure/Virus/$Keyname.fq1.fastv.fq -c $db/Fastv/viral.kc.fasta -h $Outdir/02.Structure/Virus/$Keyname.fastv.html -j $Outdir/02.Structure/Virus/$Keyname.fastv.json\n";
		system "rm -rf $Outdir/02.Structure/Virus/$Keyname.fq1.fastv.fq\n";
	}
	system "cp $Outdir/02.Structure/Virus/$Keyname.fastv.html $Outdir/WGS_Report/$Keyname.fastv.html\n";
}

###Kraken2+Bracken
print "Kraken2 searching\n";
system "mkdir $Outdir/02.Structure/Kraken2\n";
if(defined $fq2){
	system "/root/miniconda3/envs/kraken/bin/kraken2 --db $db/Kraken2 --report $Outdir/02.Structure/Kraken2/$Keyname.report --output $Outdir/02.Structure/Kraken2/$Keyname.kraken --paired $fq1 $fq2\n";
	system "/root/miniconda3/envs/kraken/bin/bracken -d $db/Kraken2 -i $Outdir/02.Structure/Kraken2/$Keyname.report -o $Outdir/02.Structure/Kraken2/$Keyname.report.bracken  -l $level -r $r\n";
	system "perl $Bin/kraken2-translate.pl $Outdir/02.Structure/Kraken2/$Keyname.report >$Outdir/02.Structure/Kraken2/$Keyname.report.txt\n";
	system "/root/miniconda3/envs/kraken/bin/ktImportText $Outdir/02.Structure/Kraken2/$Keyname.report.txt -o $Outdir/02.Structure/Kraken2/$Keyname.report.html";
}else{
	system "/root/miniconda3/envs/kraken/bin/kraken2 --db $db/Kraken2 --report $Outdir/02.Structure/Kraken2/$Keyname.report --output $Outdir/02.Structure/Kraken2/$Keyname.kraken $fq1 \n";
	system "/root/miniconda3/envs/kraken/bin/bracken -d $db/Kraken2 -i $Outdir/02.Structure/Kraken2/$Keyname.report -o $Outdir/02.Structure/Kraken2/$Keyname.report.bracken  -l $level -r $r\n";
	system "perl $Bin/kraken2-translate.pl $Outdir/02.Structure/Kraken2/$Keyname.report >$Outdir/02.Structure/Kraken2/$Keyname.report.txt\n";
	system "/root/miniconda3/envs/kraken/bin/ktImportText $Outdir/02.Structure/Kraken2/$Keyname.report.txt -o $Outdir/02.Structure/Kraken2/$Keyname.report.html";
}
system "cp $Outdir/02.Structure/Kraken2/$Keyname.report.html $Outdir/WGS_Report/$Keyname.kraken2.html\n";
system "cp $Outdir/02.Structure/Kraken2/$Keyname.report.bracken  $Outdir/WGS_Report/$Keyname.bracken.report.txt\n ";

###Metaphlain 经测试，无法和其他程序合并运行，需要单独运行conda 
##system "conda activate metaphlan\n";
#system "mkdir $Outdir/02.Structure/Metaphlan\n";
#system "mkdir $Outdir/02.Structure/Metaphlan/sams $Outdir/02.Structure/Metaphlan/bowtie2 $Outdir/02.Structure/Metaphlan/profiles $Outdir/02.Structure/Metaphlan/consensus_markers\n";
#$ENV{'PATH'} .= ":/root/miniconda3/envs/metaphlan/lib/python3.6/site-packages/metaphlan/utils";
#export PATH="$PATH:/root/miniconda3/envs/metaphlan/lib/python3.6/site-packages/metaphlan/utils"
#if(defined $fq2)  {
#system " /root/miniconda3/envs/metaphlan/bin/metaphlan $fq1,$fq2 --input fastq -s $Outdir/02.Structure/Metaphlan/sams/$Keyname.sam.bz2 --bowtie2out $Outdir/02.Structure/Metaphlan/bowtie2/$Keyname.bowtie2.bz2 -o $Outdir/02.Structure/Metaphlan/profiles/$Keyname.profile.tsv --bowtie2db $db/Metaphlan --bowtie2_exe /root/miniconda3/envs/metaphlan/bin/bowtie2\n";
#                system "/root/miniconda3/envs/metaphlan/bin/sample2markers.py -i $Outdir/02.Structure/Metaphlan/sams/$Keyname.sam.bz2 -o $Outdir/02.Structure/Metaphlan/consensus_markers -n 8\n";
#}else{
#	system " /root/miniconda3/envs/metaphlan/bin/metaphlan $fq1 --input fastq -s $Outdir/02.Structure/Metaphlan/sams/$Keyname.sam.bz2 --bowtie2out $Outdir/02.Structure/Metaphlan/bowtie2/$Keyname.bowtie2.bz2 -o $Outdir/02.Structure/Metaphlan/profiles/$Keyname.profile.tsv --bowtie2db $db/Metaphlan --bowtie2_exe /root/miniconda3/envs/metaphlan/bin/bowtie2\n";
#                system "/root/miniconda3/envs/metaphlan/bin/sample2markers.py -i $Outdir/02.Structure/Metaphlan/sams/$Keyname.sam.bz2 -o $Outdir/02.Structure/Metaphlan/consensus_markers -n 8\n";
#}
#system "$Outdir/02.Structure/Metaphlan/profiles/$Keyname.profile.tsv $Outdir/WGS_Report/$Keyname.metaphlan.report.txt\n";

###VFF
unless(-e "$Outdir/03.Target"){system "mkdir $Outdir/03.Target\n";}

if(defined $dovf){
	system "mkdir $Outdir/03.Target/VFF_ARGs\n";
	print "Virulent factor searching\n";
	if(defined $fq2){
		system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -2 $fq2 -T $db/VFF/VFDB_setB_nt.fas -o $Outdir/03.Target/VFF_ARGs/$Keyname.vff.csv\n ";
	}else{
		system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -T $db/VFF/VFDB_setB_nt.fas -o $Outdir/03.Target/VFF_ARGs/$Keyname.vff.csv\n ";
	}
	system "cp $Outdir/03.Target/VFF_ARGs/$Keyname.vff.csv.bowtie.report $Outdir/WGS_Report/$Keyname.vff.csv\n";
}

##Resistance
if(defined $dores){
	unless(-e "$Outdir/03.Target/VFF_ARGs"){system "mkdir $Outdir/03.Target/VFF_ARGs\n";}
	print "Resistance genes searching\n";
	if(defined $fq2){
		system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -2 $fq2 -T $db/ARG/res.fas -o $Outdir/03.Target/VFF_ARGs/$Keyname.ARG.csv\n ";
	}else{
		system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -T $db/ARG/res.fas -o $Outdir/03.Target/VFF_ARGs/$Keyname.ARG.csv\n ";
	}
	system "cp $Outdir/03.Target/VFF_ARGs/$Keyname.ARG.csv.bowtie.report $Outdir/WGS_Report/$Keyname.ARG.csv\n";
}


if(defined $Ref){
	system "mkdir $Outdir/03.Target/Reference\n";
	print "Target searching\n";
	if(defined $fq2){
	system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -2 $fq2 -T $Ref -o $Outdir/03.Target/Reference/$Keyname.reference.csv\n ";
	}else{
		system "perl $Bin/Target_bowtie-docker.pl -1 $fq1 -T $Ref -o $Outdir/03.Target/Reference/$Keyname.reference.csv\n ";
	}
	system "cp $Outdir/03.Target/Reference/$Keyname.reference.csv.bowtie.report $Outdir/WGS_Report/$Keyname.reference.csv\n";
}
###host
if(defined $host){
		system "mkdir $Outdir/03.Target/Host\n";
	if(defined $fq2){
	system "perl $Bin/Host.pl -1 $fq1 -2 $fq2 -o $Outdir/03.Target/Host/$Keyname.host -T $db/Host/Mitochondrion.fasta \n";
	}else{
		system "perl $Bin/Host.pl -1 $fq1  -o $Outdir/03.Target/Host/$Keyname.host -T $db/Host/Mitochondrion.fasta\n";
	}
	system "cp $Outdir/03.Target/Host/$Keyname.host.bowtie.report $Outdir/WGS_Report/$Keyname.host.csv\n";
}

####组装分箱
if(defined $assemble){
	#system "conda deactivate\n";
	system "mkdir $Outdir/04.Assembled\n";
	
	if(defined $fq2){
		
		system "spades.py -1 $fq1 -2 $fq2 -o $Outdir/04.Assembled/$Keyname --meta \n";
	}else{
		system "spades.py -s $fq1 -o $Outdir/04.Assembled/$Keyname --meta\n";
	}

	my $genome=$Outdir."/04.Assembled/".$Keyname.".assembled.fasta";
	system "cp $Outdir/04.Assembled/$Keyname/scaffolds.fasta $genome\n";
	system "cp $Outdir/04.Assembled/$Keyname/scaffolds.fasta $Outdir/WGS_Report/$Keyname.assembled.fasta\n";
#	system "conda activate metaphlan\n";
	system "/root/miniconda3/envs/metaphlan/bin/bowtie2-build  $genome $genome\n";
   system "/root/miniconda3/envs/metaphlan/bin/bowtie2  -x $genome -1 $fq1 -2 $fq2 | /root/miniconda3/envs/kraken/bin/samtools sort -o $Outdir/04.Assembled/$Keyname.sort.bam  \n";
#   system "conda deactivate\n";
   system "jgi_summarize_bam_contig_depths --outputDepth $Outdir/04.Assembled/$Keyname.depth.txt $Outdir/04.Assembled/$Keyname.sort.bam\n";
   system "metabat2 -i $genome -a $Outdir/04.Assembled/$Keyname.depth.txt -o $Outdir/04.Assembled/$Keyname\n";
  
}


###去除人类序列
if (defined $nohuman) {
    my $target = "$db/Human/GCF_000001405.40_GRCh38.p14_genomic.fna";
    unless (-e "$target.1.bt2") {
        system "/root/miniconda3/envs/metaphlan/bin/bowtie2-build $target $target\n";
    }
    if (defined $fq2) {
        system "/root/miniconda3/envs/metaphlan/bin/bowtie2 -1 $fq1 -2 $fq2 -x $target -S $Outdir/03.Target/$Keyname.nohuman.sam --no-unal\n";
    } else {
        system "/root/miniconda3/envs/metaphlan/bin/bowtie2 -U $fq1 -x $target -S $Outdir/03.Target/$Keyname.nohuman.sam --no-unal\n";
    }

    open(FILE, "$Outdir/03.Target/$Keyname.nohuman.sam");
    my %read;
    while (my $line = <FILE>) {
        chomp $line;
        next if substr($line, 0, 1) eq "@";
        my @a = split " ", $line;
        next unless $a[2] =~ /[0-9a-zA-Z]/;
        my $tmp;
        foreach my $a (@a) {
            if ($a =~ /^MD:Z:/) {
                $tmp = $a;
            }
        }
        my $len = length($a[9]);
        my $l = length($tmp);
        my $num = "";
        my $match = 0;
        foreach (0..($l-1)) {
            my $site = substr($tmp, $_, 1);
            if ($site =~ /[0-9]/) {
                $num .= $site;
            } else {
                $match += $num;
                $num = "";
            }
        }
        $match += $num;
        next unless $len > 0;
        next unless $match > 0;
        my $maper = $match / $len * 100;
        next if $match < 70 && $maper < 80;
        my @b = split "/", $a[0];
        $read{$b[0]}++;
    }
    close FILE;

    my @key = keys %read;
    my %result;
    my $num = 0;
    foreach my $key (@key) {
        next if $read{$key} < 1;
        $result{$key}++;
        $num++;
    }

    open(my $O1, ">", "$Outdir/WGS_Report/$Keyname.nohuman.R1.fq");
    open(FILE, $fq1);
    while (my $line = <FILE>) {
        chomp $line;
        my $seq = <FILE>; chomp $seq;
        my $a = <FILE>; chomp $a;
        my $b = <FILE>; chomp $b;
        $line = substr($line, 1);
        my @name = split " ", $line;
        my @b = split "/", $name[0];
        my $name = $b[0];
        unless (exists $read{$name}) {
            print $O1 "@$name\n$seq\n$a\n$b\n";
        }
    }
    close FILE; close $O1;

    if (defined $fq2) {
        open(my $O2, ">", "$Outdir/WGS_Report/$Keyname.nohuman.R2.fq");
        open(FILE, $fq2);
        while (my $line = <FILE>) {
            chomp $line;
            my $seq = <FILE>; chomp $seq;
            my $a = <FILE>; chomp $a;
            my $b = <FILE>; chomp $b;
            $line = substr($line, 1);
            my @name = split " ", $line;
            my @b = split "/", $name[0];
            my $name = $b[0];
            unless (exists $result{$name}) {
                print $O2 "@$name\n$seq\n$a\n$b\n";
            }
        }
        close FILE; close $O2;
    }
}
#system "multiqc $file -c $Bin/multiqc_config.yaml\n";
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