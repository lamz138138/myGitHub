#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my ($help, $selectReads);
GetOptions(
	'h|help'		=>	\$help,
	's|selectReads=s'	=>	\$selectReads
)||usage(); 
usage () if defined $help;
$ARGV[0]='-' unless defined $ARGV[0];
open IN_1,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";

my %readsName;

if(defined $selectReads){
	open READNAME, $selectReads;
	while(<READNAME>){
		chomp;
		$readsName{$_}=0;
	}
}

my $readName;
while(<IN_1>){
	chomp;
	if(defined $selectReads){
		if(/^>/){
			$readName=$_;
			$readName=~s/\s.*//;
			$readName=~s/^>//;
			say $_ if exists $readsName{$readName};
		}else{
			say $_ if exists $readsName{$readName};
		}
	}
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to process fasta
Usage: perl $scriptName *.fasta >output

	-s 	--selectReads		file contain name of selecting read
	-h 	--help           	print this help information screen
HELP
    exit(-1);
}