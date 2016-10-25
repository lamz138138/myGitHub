#!/usr/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my ($help, $sourceName);
$sourceName="source";
GetOptions(
	's|sourceName=s'	=>	\$sourceName,
    'h|help' 			=>	\$help
)||usage(); 
usage () if defined $help;

$ARGV[0]='-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";


say "##gff-version 3";
while(<IN>){
    chomp;
	my ($chr, $start, $end, $name, $score, $strand, $thickStart, $thickEnd, $itemGgb, $blockCount, $blockSizes, $blockStarts)=split "\t";
	my @BlockSizes=split ",",$blockSizes;
	my @BlockStarts=split ",",$blockStarts;
	my $mRNAName=$chr.":".$start."-".$end;
	my $outputAttributes="ID=".$mRNAName;
	say join "\t",($chr, $sourceName, "mRNA", $start+1, $end, $score, $strand, ".", $outputAttributes);
	for(my $i=0; $i<$blockCount; $i++){
		my $tempStart=$start+$BlockStarts[$i];
		my $tempEnd=$tempStart+$BlockSizes[$i];
		$outputAttributes="ID=".$chr.":".$tempStart."-".$tempEnd.";Parent=".$mRNAName;
		say join "\t",($chr, $sourceName, "exon", $tempStart+1, $tempEnd, $score, $strand, ".", $outputAttributes);
	}
}


sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to change bed12 to gff3 format
Usage: perl $scriptName *.bed >*.gff3
	if INPUT not specified, input from STDIN, output to STDOUT

	-s --sourceName 	The name of column source [default: source]
	-h --help			Print this help information
HELP
	exit(-1);
}