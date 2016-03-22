#!/usr/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my ($help,$strand,$chr,$exon13,$exon14,$exon15,$exon16,$exon17,$exon20,$exon20b,$exon21,$exon22);
$chr="chr13";
$strand="+";
$exon13="NA-NA";
$exon14="NA-NA";
$exon15="NA-NA";
$exon16="NA-NA";
$exon17="NA-NA";
$exon20="NA-NA";
$exon20b="NA-NA";
$exon21="NA-NA";
$exon22="NA-NA";
GetOptions(
    'h|help' 	=>	\$help,
	's|strand=s' 	=>	\$strand,
	'c|chr=s'		=>	\$chr,
	'exon13=s'	=>	\$exon13,
	'exon14=s'	=>	\$exon14,
	'exon15=s'	=>	\$exon15,
	'exon16=s'	=>	\$exon16,
	'exon17=s'	=>	\$exon17,
	'exon20=s'	=>	\$exon20,
	'exon20b=s'	=>	\$exon20b,
	'exon21=s'	=>	\$exon21,
	'exon22=s'	=>	\$exon22
)||usage(); 
usage () if defined $help;

$ARGV[0]='-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

# Following get each position of each exon
my %Exon_pos;
$Exon_pos{13}=$exon13;
$Exon_pos{14}=$exon14;
$Exon_pos{15}=$exon15;
$Exon_pos{16}=$exon16;
$Exon_pos{17}=$exon17;
$Exon_pos{20}=$exon20;
$Exon_pos{21}=$exon21;
$Exon_pos{22}=$exon22;
$Exon_pos{"20b"}=$exon20b;

# Following get start and end of each exon
my %Exon_s_e;
my @Number=(13..17,20..22,"20b");
my @Exon_s_e_2;
my %PosToExon;
my %Count2;
foreach my $number (@Number){
	my $temp_exon=$Exon_pos{$number};
	$temp_exon=~/(.*)-(.*)/;
	$Exon_s_e{$number."_s"}=$1;
	$Exon_s_e{$number."_e"}=$2;
	push @Exon_s_e_2,($1,$2);
	$PosToExon{$1}=$number;
	$PosToExon{$2}=$number;
	$Count2{"OthersTo".$number}=0;
	$Count2{$number."ToOthers"}=0;
}

my @Junc_type=("13To14", "13To15", "13To16", "13To17",
			   "14To15", "14To16", "14To17",
			   "15To16", "15To17",
			   "16To17",
			   "20To21", "20To22",
			   "20bTo21", "20bTo22",
			   "21To22");

# Following get link position of each junction type
my %Junc;
my %Count;
if($strand eq "+"){
	foreach my $junc_type (@Junc_type){
		$junc_type =~ /(.*)To(.*)/;
		my $exon_1=$1;
		my $exon_2=$2;
		my $temp_s=$Exon_s_e{$exon_1."_e"};
		my $temp_e=$Exon_s_e{$exon_2."_s"};
		$Junc{$junc_type}=$temp_s.":".$temp_e;
		$Count{$junc_type}=0;
		$Count2{$junc_type}=0;
	}
}else{
	foreach my $junc_type (@Junc_type){
		$junc_type =~ /(.*)To(.*)/;
		my $exon_1=$2;
		my $exon_2=$1;
		my $temp_s=$Exon_s_e{$exon_1."_e"};
		my $temp_e=$Exon_s_e{$exon_2."_s"};
		$Junc{$junc_type}=$temp_s.":".$temp_e;
		$Count{$junc_type}=0;
		$Count2{$junc_type}=0;
	}
}

# Following read file  
while(<IN>){
	chomp;
	next if /^#/;
	my ($readCount, $junc_pos)=(split "\t")[4,12];
	foreach my $temp_junc_type (@Junc_type){
		next if $junc_pos ne $Junc{$temp_junc_type};
		if($junc_pos eq $Junc{$temp_junc_type}){
			$Count{$temp_junc_type}=$readCount;
			last;
		}
	}
	
	# Following considering other type junction
	if($strand eq "+"){
		foreach my $temp_pos (@Exon_s_e_2){
			next if $junc_pos !~ $temp_pos;
			my ($pos_1, $pos_2)=(split ":", $junc_pos);
			my ($number_1, $number_2);
			if(exists $PosToExon{$pos_1}){
				$number_1=$PosToExon{$pos_1};
			}else{
				$number_1="Others";
			}
			if(exists $PosToExon{$pos_2}){
				$number_2=$PosToExon{$pos_2};
			}else{
				$number_2="Others";
			}
			if($number_1 ne "Others" && $number_2 ne "Others"){
				$Count2{$number_1."To".$number_2}=$readCount;
			}else{
				$Count2{$number_1."To".$number_2}+=$readCount;
			}
		}
	}else{
		foreach my $temp_pos (@Exon_s_e_2){
			next if $junc_pos !~ $temp_pos;
			my ($pos_1, $pos_2)=(split ":", $junc_pos);
			my ($number_1, $number_2);
			if(exists $PosToExon{$pos_1}){
				$number_1=$PosToExon{$pos_1};
			}else{
				$number_1="Others";
			}
			if(exists $PosToExon{$pos_2}){
				$number_2=$PosToExon{$pos_2};
			}else{
				$number_2="Others";
			}
			if($number_1 ne "Others" && $number_2 ne "Others"){
				$Count2{$number_2."To".$number_1}=$readCount;
			}else{
				$Count2{$number_2."To".$number_1}+=$readCount;
			}
		}
	}
}

# Following calculate the percentage of each junction type
my $sum_13ToAll=$Count{"13To14"} + $Count{"13To15"} + $Count{"13To16"} + $Count{"13To17"};
my %Per;
($Per{"13To14"}, $Per{"13To15"}, $Per{"13To16"}, $Per{"13To17"})= (0) x 4;
if($sum_13ToAll>0){
	$Per{"13To14"}=sprintf("%.6f", $Count{"13To14"}/$sum_13ToAll);
	$Per{"13To15"}=sprintf("%.6f", $Count{"13To15"}/$sum_13ToAll);
	$Per{"13To16"}=sprintf("%.6f", $Count{"13To16"}/$sum_13ToAll);
	$Per{"13To17"}=sprintf("%.6f", $Count{"13To17"}/$sum_13ToAll);
}


my $sum_14ToAll=$Count{"14To15"} + $Count{"14To16"} + $Count{"14To17"};
($Per{"14To15"}, $Per{"14To16"}, $Per{"14To17"})=(0) x 3;
if($sum_14ToAll > 0){
	$Per{"14To15"}=sprintf("%.6f", $Count{"14To15"}/$sum_14ToAll);
	$Per{"14To16"}=sprintf("%.6f", $Count{"14To16"}/$sum_14ToAll);
	$Per{"14To17"}=sprintf("%.6f", $Count{"14To17"}/$sum_14ToAll);
}


my $sum_15ToAll=$Count{"15To16"} + $Count{"15To17"};
($Per{"15To16"}, $Per{"15To17"})=(0) x 2;
if($sum_15ToAll > 0){
	$Per{"15To16"}=sprintf("%.6f", $Count{"15To16"}/$sum_15ToAll);
	$Per{"15To17"}=sprintf("%.6f", $Count{"15To17"}/$sum_15ToAll);
}

my $sum_16ToAll=$Count{"16To17"};
$Per{"16To17"}=0;
if($sum_16ToAll > 0){
	$Per{"16To17"}=sprintf("%.6f", $Count{"16To17"}/$sum_16ToAll);
}

=pod
my $sum_20ToAll=$Count{"20To21"} + $Count{"20To22"} + $Count{"20bTo21"} + $Count{"20bTo22"};
($Per{"20To21"}, $Per{"20To22"}, $Per{"20bTo21"}, $Per{"20bTo22"})= (0) x 4;
if($sum_20ToAll>0){
	$Per{"20To21"}=sprintf("%.6f", $Count{"20To21"}/$sum_20ToAll);
	$Per{"20To22"}=sprintf("%.6f", $Count{"20To22"}/$sum_20ToAll);
	$Per{"20bTo21"}=sprintf("%.6f", $Count{"20bTo21"}/$sum_20ToAll);
	$Per{"20bTo22"}=sprintf("%.6f", $Count{"20bTo22"}/$sum_20ToAll);
}
=cut

my $sum_20ToAll=$Count{"20To21"} + $Count{"20To22"};
($Per{"20To21"}, $Per{"20To22"})= (0) x 2;
if($sum_20ToAll>0){
	$Per{"20To21"}=sprintf("%.6f", $Count{"20To21"}/$sum_20ToAll);
	$Per{"20To22"}=sprintf("%.6f", $Count{"20To22"}/$sum_20ToAll);
}

my $sum_20bToAll=$Count{"20bTo21"} + $Count{"20bTo22"};
($Per{"20bTo21"}, $Per{"20bTo22"})= (0) x 2;
if($sum_20bToAll>0){
	$Per{"20bTo21"}=sprintf("%.6f", $Count{"20bTo21"}/$sum_20bToAll);
	$Per{"20bTo22"}=sprintf("%.6f", $Count{"20bTo22"}/$sum_20bToAll);
}

my $sum_21ToAll=$Count{"21To22"};
$Per{"21To22"}=0;
if($sum_21ToAll > 0){
	$Per{"21To22"}=sprintf("%.6f", $Count{"21To22"}/$sum_21ToAll);
}

foreach my $junc_type (@Junc_type){
	if($Junc{$junc_type}=~/NA/){
		$Count{$junc_type}="NA";
		$Per{$junc_type}="NA";
	}
}

say join "\t",("Exon_Type","exon13-exon14","exon13-exon15","exon13-exon16","exon13-exon17",
			   "exon14-exon15","exon14-exon16","exon14-exon17",
			   "exon15-exon16","exon15-exon17",
			   "exon16-exon17",
			   "exon20-exon21","exon20-exon22",
			   "exon20b-exon21","exon20b-exon22",
			   "exon21-exon22");
say join "\t",("Reads_Count",$Count{"13To14"}, $Count{"13To15"}, $Count{"13To16"}, $Count{"13To17"},
			   $Count{"14To15"}, $Count{"14To16"}, $Count{"14To17"},
			   $Count{"15To16"}, $Count{"15To17"},
			   $Count{"16To17"},
			   $Count{"20To21"}, $Count{"20To22"},
			   $Count{"20bTo21"}, $Count{"20bTo22"},
			   $Count{"21To22"});
say join "\t",("Percent",$Per{"13To14"}, $Per{"13To15"}, $Per{"13To16"}, $Per{"13To17"},
			   $Per{"14To15"}, $Per{"14To16"}, $Per{"14To17"},
			   $Per{"15To16"}, $Per{"15To17"},
			   $Per{"16To17"},
			   $Per{"20To21"}, $Per{"20To22"},
			   $Per{"20bTo21"}, $Per{"20bTo22"},
			   $Per{"21To22"});


# Following calculate the percentage of each junction type (considering others)
my $sum2_13ToAll=$Count2{"13To14"} + $Count2{"13To15"} + $Count2{"13To16"} + $Count2{"13To17"} + $Count2{"13ToOthers"};
my %Per2;
($Per2{"13To14"}, $Per2{"13To15"}, $Per2{"13To16"}, $Per2{"13To17"}, $Per2{"13ToOthers"})= (0) x 5;
if($sum2_13ToAll>0){
	$Per2{"13To14"}=sprintf("%.6f", $Count2{"13To14"}/$sum2_13ToAll);
	$Per2{"13To15"}=sprintf("%.6f", $Count2{"13To15"}/$sum2_13ToAll);
	$Per2{"13To16"}=sprintf("%.6f", $Count2{"13To16"}/$sum2_13ToAll);
	$Per2{"13To17"}=sprintf("%.6f", $Count2{"13To17"}/$sum2_13ToAll);
	$Per2{"13ToOthers"}=sprintf("%.6f", $Count2{"13ToOthers"}/$sum2_13ToAll);
}


my $sum2_14ToAll=$Count2{"14To15"} + $Count2{"14To16"} + $Count2{"14To17"} + $Count2{"14ToOthers"};
($Per2{"14To15"}, $Per2{"14To16"}, $Per2{"14To17"}, $Per2{"14ToOthers"})=(0) x 4;
if($sum2_14ToAll > 0){
	$Per2{"14To15"}=sprintf("%.6f", $Count2{"14To15"}/$sum2_14ToAll);
	$Per2{"14To16"}=sprintf("%.6f", $Count2{"14To16"}/$sum2_14ToAll);
	$Per2{"14To17"}=sprintf("%.6f", $Count2{"14To17"}/$sum2_14ToAll);
	$Per2{"14ToOthers"}=sprintf("%.6f", $Count2{"14ToOthers"}/$sum2_14ToAll);
}


my $sum2_15ToAll=$Count2{"15To16"} + $Count2{"15To17"} + $Count2{"15ToOthers"};
($Per2{"15To16"}, $Per2{"15To17"}, $Per2{"15ToOthers"})=(0) x 3;
if($sum2_15ToAll > 0){
	$Per2{"15To16"}=sprintf("%.6f", $Count2{"15To16"}/$sum2_15ToAll);
	$Per2{"15To17"}=sprintf("%.6f", $Count2{"15To17"}/$sum2_15ToAll);
	$Per2{"15ToOthers"}=sprintf("%.6f", $Count2{"15ToOthers"}/$sum2_15ToAll);
}

my $sum2_16ToAll=$Count2{"16To17"} + $Count2{"16ToOthers"};
($Per2{"16To17"}, $Per2{"16ToOthers"})=(0) x 2;
if($sum2_16ToAll > 0){
	$Per2{"16To17"}=sprintf("%.6f", $Count2{"16To17"}/$sum2_16ToAll);
	$Per2{"16ToOthers"}=sprintf("%.6f", $Count2{"16ToOthers"}/$sum2_16ToAll);
}

if($Count2{"17ToOthers"} > 0){
	$Per2{"17ToOthers"}=1;
}else{
	$Per2{"17ToOthers"}=0;
}

=pod
my $sum2_20ToAll=$Count2{"20To21"} + $Count2{"20To22"} + $Count2{"20bTo21"} + $Count2{"20bTo22"} + $Count2{"20ToOthers"} + $Count2{"20bToOthers"};
($Per2{"20To21"}, $Per2{"20To22"}, $Per2{"20bTo21"}, $Per2{"20bTo22"}, $Per2{"20ToOthers"}, $Per2{"20bToOthers"})= (0) x 6;
if($sum2_20ToAll>0){
	$Per2{"20To21"}=sprintf("%.6f", $Count2{"20To21"}/$sum2_20ToAll);
	$Per2{"20To22"}=sprintf("%.6f", $Count2{"20To22"}/$sum2_20ToAll);
	$Per2{"20bTo21"}=sprintf("%.6f", $Count2{"20bTo21"}/$sum2_20ToAll);
	$Per2{"20bTo22"}=sprintf("%.6f", $Count2{"20bTo22"}/$sum2_20ToAll);
	$Per2{"20ToOthers"}=sprintf("%.6f", $Count2{"20ToOthers"}/$sum2_20ToAll);
	$Per2{"20bToOthers"}=sprintf("%.6f", $Count2{"20bToOthers"}/$sum2_20ToAll);
}
=cut

my $sum2_20ToAll=$Count2{"20To21"} + $Count2{"20To22"} + $Count2{"20ToOthers"};
($Per2{"20To21"}, $Per2{"20To22"}, $Per2{"20ToOthers"})= (0) x 3;
if($sum2_20ToAll>0){
	$Per2{"20To21"}=sprintf("%.6f", $Count2{"20To21"}/$sum2_20ToAll);
	$Per2{"20To22"}=sprintf("%.6f", $Count2{"20To22"}/$sum2_20ToAll);
	$Per2{"20ToOthers"}=sprintf("%.6f", $Count2{"20ToOthers"}/$sum2_20ToAll);
}

my $sum2_20bToAll=$Count2{"20bTo21"} + $Count2{"20bTo22"} + $Count2{"20bToOthers"};
($Per2{"20bTo21"}, $Per2{"20bTo22"}, $Per2{"20bToOthers"})= (0) x 3;
if($sum2_20bToAll>0){
	$Per2{"20bTo21"}=sprintf("%.6f", $Count2{"20bTo21"}/$sum2_20bToAll);
	$Per2{"20bTo22"}=sprintf("%.6f", $Count2{"20bTo22"}/$sum2_20bToAll);
	$Per2{"20bToOthers"}=sprintf("%.6f", $Count2{"20bToOthers"}/$sum2_20bToAll);
}

my $sum2_21ToAll=$Count2{"21To22"} + $Count2{"21ToOthers"};
($Per2{"21To22"}, , $Per2{"21ToOthers"})=(0) x 2;
if($sum2_21ToAll > 0){
	$Per2{"21To22"}=sprintf("%.6f", $Count2{"21To22"}/$sum2_21ToAll);
	$Per2{"21ToOthers"}=sprintf("%.6f", $Count2{"21ToOthers"}/$sum2_21ToAll);
}

if($Count2{"22ToOthers"} > 0){
	$Per2{"22ToOthers"}=1;
}else{
	$Per2{"22ToOthers"}=0;
}

foreach my $junc_type (@Junc_type){
	if($Junc{$junc_type}=~/NA/){
		$Count2{$junc_type}="NA";
		$Per2{$junc_type}="NA";
	}
}

foreach my $number (@Number){
	my $temp_exon=$Exon_pos{$number};
	if($temp_exon eq "NA-NA"){
		$Count2{"OthersTo".$number}="NA";
		$Count2{$number."ToOthers"}="NA";
		$Per2{"OthersTo".$number}="NA";
		$Per2{$number."ToOthers"}="NA";
	}else{
		if($Count2{"OthersTo".$number}>0){
			$Per2{"OthersTo".$number}=1;
		}else{
			$Per2{"OthersTo".$number}=0;
		}
	}
}

say STDERR join "\t",("Exon_Type","Others-exon13","exon13-exon14","exon13-exon15","exon13-exon16","exon13-exon17","exon13-Others",
			   "Others-exon14","exon14-exon15","exon14-exon16","exon14-exon17","exon14-Others",
			   "Others-exon15","exon15-exon16","exon15-exon17","exon15-Others",
			   "Others-exon16","exon16-exon17","exon16-Others",
			   "Others-exon17","exon17-Others",
			   "Others-exon20","exon20-exon21","exon20-exon22","exon20-Others",
			   "Others-exon20b","exon20b-exon21","exon20b-exon22","exon20b-Others",
			   "Others-exon21","exon21-exon22","exon21-Others",
			   "Others-exon22","exon22-Others");
say STDERR join "\t",("Reads_Count",$Count2{"OthersTo13"},$Count2{"13To14"}, $Count2{"13To15"}, $Count2{"13To16"}, $Count2{"13To17"},$Count2{"13ToOthers"},
			   $Count2{"OthersTo14"},$Count2{"14To15"}, $Count2{"14To16"}, $Count2{"14To17"},$Count2{"14ToOthers"},
			   $Count2{"OthersTo15"},$Count2{"15To16"}, $Count2{"15To17"},$Count2{"15ToOthers"},
			   $Count2{"OthersTo16"},$Count2{"16To17"},$Count2{"16ToOthers"},
			   $Count2{"OthersTo17"},$Count2{"17ToOthers"},
			   $Count2{"OthersTo20"},$Count2{"20To21"}, $Count2{"20To22"},$Count2{"20ToOthers"},
			   $Count2{"OthersTo20b"},$Count2{"20bTo21"}, $Count2{"20bTo22"},$Count2{"20bToOthers"},
			   $Count2{"OthersTo21"},$Count2{"21To22"},$Count2{"21ToOthers"},
			   $Count2{"OthersTo22"},$Count2{"22ToOthers"});
say STDERR join "\t",("Percent",$Per2{"OthersTo13"},$Per2{"13To14"}, $Per2{"13To15"}, $Per2{"13To16"}, $Per2{"13To17"},$Per2{"13ToOthers"},
			   $Per2{"OthersTo14"},$Per2{"14To15"}, $Per2{"14To16"}, $Per2{"14To17"},$Per2{"14ToOthers"},
			   $Per2{"OthersTo15"},$Per2{"15To16"}, $Per2{"15To17"},$Per2{"15ToOthers"},
			   $Per2{"OthersTo16"},$Per2{"16To17"},$Per2{"16ToOthers"},
			   $Per2{"OthersTo17"},$Per2{"17ToOthers"},
			   $Per2{"OthersTo20"},$Per2{"20To21"}, $Per2{"20To22"},$Per2{"20ToOthers"},
			   $Per2{"OthersTo20b"},$Per2{"20bTo21"}, $Per2{"20bTo22"},$Per2{"20bToOthers"},
			   $Per2{"OthersTo21"},$Per2{"21To22"},$Per2{"21ToOthers"},
			   $Per2{"OthersTo22"},$Per2{"22ToOthers"});


sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to calculate percentage of different junction types in two variable region, junction type must be added in junc.bed as "A:B" format
Usage: perl $scriptName -c chr13 -s - --exon13 A-B --exon 14 C-D *junc.bed >output.file
Output to STDOUT

-d 	--dataset 	name of the dataset 
-c 	--chr 		chr of CAMK2D [default: chr13]
-s 	--strand 	strand of CAMK2D [default: +]
	--exon13 	position of exon 13 [default: NA-NA]
	--exon14 	position of exon 14 [default: NA-NA]
	--exon15 	position of exon 15 [default: NA-NA]
	--exon16 	position of exon 16 [default: NA-NA]
	--exon17 	position of exon 17 [default: NA-NA]
	--exon20 	position of exon 20 [default: NA-NA]
	--exon21 	position of exon 21 [default: NA-NA]
	--exon22 	position of exon 22 [default: NA-NA]
	--exon20b 	position of exon 20b [default: NA-NA]
HELP
	exit(-1);
}