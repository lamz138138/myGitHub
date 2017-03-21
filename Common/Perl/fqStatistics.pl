#!/usr/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my($help);
GetOptions(
	    'h|help'	=>	\$help
	    )||usage(); 
usage () if defined $help;
$ARGV[0]='-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

my $maxLength=0;
my $minLength=0;
my $meanLength=0;
my $n50Length=0;
my $totalBases=0;
my $readNumber=0;

my %hash;
while(<IN>){
	chomp;
	my $readName=$_;
	chomp(my $sequence=<IN>);
	chomp(my $readName2=<IN>);
	chomp(my $quality=<IN>);
	$readNumber++;
	my $length=length($sequence);
	$totalBases+=$length;
	$minLength=$length if ($minLength!=0 && $minLength>$length) || $minLength==0;
	$maxLength=$length if ($maxLength!=0 && $maxLength<$length) || $maxLength==0;
	if(exists $hash{$length}){
		$hash{$length}+=1;
	}else{
		$hash{$length}=1;
	}
}

$meanLength=sprintf("%.4f", $totalBases/$readNumber);
my $sumLength=0;
foreach my $length (sort {$a<=>$b} keys %hash){
	my $lengthNumber=$hash{$length};
	for(my $i=0; $i<$lengthNumber; $i++){	
		$sumLength+=$length;
		if($sumLength>$totalBases/2){
			$n50Length=$length;
			last;
		}
	}
	last if $n50Length!=0;
}

say join "\t",("Max Length:", $maxLength);
say join "\t",("Min Length:", $minLength);
say join "\t",("Mean Length:", $meanLength);
say join "\t",("N50 Length:", $n50Length);
say join "\t",("Total Bases:", $totalBases);
say join "\t",("Reads Count:", $readNumber);


sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to get statistics for fastq
    Usage: perl $scriptName *fq >output
	If Input_gpe.file not specified, input from STDIN
	Output to STDOUT

	-h --help		Print this help information
HELP
	exit(-1);
}
