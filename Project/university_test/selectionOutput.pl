#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my($help,$dilimiter1,$dilimiter2,$colmun11,$colmun12);
$dilimiter1="\t";
$dilimiter2="\t";
$colmun11=1;
$colmun12=1;
GetOptions(
                'h|help'=>	\$help,
                'dilimiter1=s'=>   \$dilimiter1,
                'dilimiter2=s'=>   \$dilimiter2,
                'colmun11=s'=>     \$colmun11,
                'colmun12=s'=>     \$colmun12
                )||usage(); 
usage () if defined $help;
$ARGV[0]='-' unless defined $ARGV[0];
$ARGV[1]='-' unless defined $ARGV[1];
open IN_1,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";
open IN_2,"$ARGV[1]" or die "Can't open $ARGV[1]:$!";

my (%hash,@samerow,@file1Specific,@file2Specific,$hash);
while(<IN_1>){
    chomp;
    my @file1_data=split $dilimiter1,$_;
    #if(defined $file1_data[$colmun11-1])
    #{
    #	$file1_data[$colmun11-1]=~s/(^\s*|\s*$)//g;
    #    $hash{$file1_data[$colmun11-1]}=0;
    #}
    my @colmuns=split ",",$colmun11;
    my $key=$file1_data[$colmuns[0]-1];
    $key=~s/(^\s*|\s*$)//g;
    for(my $i=1; $i<@colmuns; $i++){
	my $tempKey=$file1_data[$colmuns[$i]-1];
	$tempKey=~s/(^\s*|\s*$)//g;
	$key.=",".$tempKey;
    }
    $hash{$key}=0 if $key ne "";
}

while(<IN_2>){
    chomp;
    my @file2_data=split $dilimiter2,$_;
    my @colmuns=split ",",$colmun12;
    my $key=$file2_data[$colmuns[0]-1];
    $key=~s/(^\s*|\s*$)//g;
    for(my $i=1; $i<@colmuns; $i++){
	my $tempKey=$file2_data[$colmuns[$i]-1];
	$tempKey=~s/(^\s*|\s*$)//g;
	$key.=",".$tempKey;
    }
    if($key ne ""){
	if(exists $hash{$key}){
	    say $_;
	}else{
	    say STDERR $_;
	}
    }
    #if(defined $file2_data[$colmun12-1]){
    #    $file2_data[$colmun12-1]=~s/(^\s*|\s*$)//g;
    #	if(exists $hash{$file2_data[$colmun12-1]}){
    #    	say $_;
    #	}else{
    #    	say STDERR $_;
    #	}
    #}
}  

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script used to output file2 data crosse with file1, file2 uniq data output to stderr
Usage: perl $scriptName file1 file2 >inFile1AndFile2 2>uniqInFile2
    if INPUT not specified, input from STDIN
    output to STDOUT

    -h --help           print this help information screen
    --dilimiter1    	dilimiter of file 1 (default is "\\t")
    --dilimiter2    	dilimiter of file 2 (default is "\\t")
    --colmun11      	colmun1s used in file 1 (default is colmunl 1 )
    --colmun12      	colmun1s used in file 2 (default is colmunl 1 )
HELP
    exit(-1);
}
