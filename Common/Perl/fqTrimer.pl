#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my ($help, $minLength, $headBases, $tailBases);
GetOptions(
	'h|help'		=>	\$help,
	'm|minLength=i'	=>	\$minLength,
	'e|headBases=i'	=>	\$headBases,
	't|tailBases=i'	=>	\$tailBases
)||usage(); 
usage () if defined $help;
$ARGV[0]='-' unless defined $ARGV[0];
open IN_1,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

my ($name_1, $sequence, $name_2, $quality, $sequenceLength, $newSequence, $newQuality);

while(<IN_1>){
    chomp;
    $name_1=$_;
    chomp ($sequence=<IN_1>);
    chomp ($name_2=<IN_1>);
    chomp ($quality=<IN_1>);
	$sequenceLength=length $sequence;
	$newSequence=$sequence;
	$newQuality=$quality;
	if(!defined $minLength){
		if(defined $tailBases){
			&trimTailBases;
		}
		if(defined $headBases){
			&trimHeadBases;
		}
		&output;
	}else{
		if($sequenceLength > $minLength){
			if(defined $tailBases){
				&trimTailBases;
				$sequenceLength=length $newSequence;
				if(defined $headBases && $sequenceLength > $minLength){
					&trimHeadBases;
					$sequenceLength=length $newSequence;
				}
			}elsif(defined $headBases){
				&trimHeadBases;
				$sequenceLength=length $newSequence;
			}
			if($sequenceLength > $minLength){
				&output;
			}
		}
	}
}

sub trimTailBases{
	$newSequence = substr($newSequence, 0, $sequenceLength-$tailBases);
	$newQuality = substr($newQuality, 0, $sequenceLength-$tailBases);
}

sub trimHeadBases{
	$newSequence = substr($newSequence, $headBases);
	$newQuality = substr($newQuality, $headBases);
}

sub output{
	say $name_1;
	say $newSequence;
	say $name_2;
	say $newQuality;
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to trim bases in fastq sequence
Usage: perl $scriptName *.fastq >trimedSequence

	-m 	--minLength	INT	Min length of output sequence (Default: NULL)
	-e 	--headBases	INT	Bases to be trimmed in the head of sequence (Default: NULL)
	-t 	--tailBases	INT	Bases to be trimmed in the tail of sequence (Default: NULL)
	-h 	--help           	print this help information screen
HELP
    exit(-1);
}