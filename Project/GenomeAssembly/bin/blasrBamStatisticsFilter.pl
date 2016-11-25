#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my($help, $longest, $coverage, $identity, $accuracy);
$coverage=0.75;
$identity=0.6;
$accuracy=0.0;
GetOptions(
    'l|longest'=>   \$longest,
    'c|coverage=f'=>    \$coverage,
    'i|identity=f'=>    \$identity,
    'a|accuracy=f'=>    \$accuracy,
    'h|help'=>	\$help
    )||usage(); 
usage () if defined $help;
$ARGV[0]='-' unless defined $ARGV[0];
open IN,"samtools view -h $ARGV[0] |" or die "Can't open $ARGV[0]:$!";

my @outputString; #This use to retain longest alignment
my $readName_1; #This use to retain before read name
my $alignLength_1;
my %hash; #This use to retain coverage and its count
while(<IN>){
    chomp;
    if(/^@/){
        say $_;
        next;
    }else{
        my ($readName_2, $ciga)=(split "\t")[0,5];
        my $myCoverage=&getTagValue("VC", $_);
        my $myIdentity=&getTagValue("VI", $_);
        my $mySimilarity=&getTagValue("VS", $_);
        my $myAccuracy=&getTagValue("VA", $_);
        my $myAS=&getTagValue("AS", $_);    #Alignment score
        my $myXL=&getTagValue("XL", $_);    #Alignment length
        if(!defined $readName_1 || $readName_2 ne $readName_1){
            &output if @outputString;
            if(defined $readName_1){
                %hash=();
            }else{
                $readName_1=$readName_2;
            }
            if(defined $longest){
                &valueToArray($_, $myXL, $myAS);
            }else{
                if($myCoverage>=$coverage && $myIdentity>=$identity && $myAccuracy>=$accuracy){
                   &valueToArray($_, $myXL, $myAS);
                }
            }
            $readName_1=$readName_2;
            $alignLength_1=$myXL;
        }else{
            if(defined $longest){
                if($myXL > $alignLength_1){
                    @outputString=();
                    &valueToArray($_, $myXL, $myAS);
                    $alignLength_1=$myXL;
                }elsif($myXL == $alignLength_1){
                    &valueToArray($_, $myXL, $myAS);
                }
            }else{
                if($myCoverage>=$coverage && $myIdentity>=$identity && $myAccuracy>=$accuracy){
                   &valueToArray($_, $myXL, $myAS);
                }
            }
        }
        &valueAssign($myCoverage, $myIdentity, $mySimilarity, $myAccuracy);
    }
}
&output if @outputString;

sub valueToArray{
    my ($myAlignValue, $myXL, $myAS)=@_;
    $hash{"alignLengths"}{$myXL}+=1;
    $hash{"alignScores"}{$myAS}+=1;
    push @outputString, $myAlignValue;
}

sub valueAssign{
    my ($coverage, $identity, $similarity, $accuracy)=@_;
    $hash{"coverage"}{$coverage}++;
    $hash{"identity"}{$identity}++;
    $hash{"similarity"}{$similarity}++;
    $hash{"accuracy"}{$accuracy}++;
}

sub valueSort{
    my @topCoverage=sort {$b<=>$a} keys %{$hash{"coverage"}};
    my @topIdentity=sort {$b<=>$a} keys %{$hash{"identity"}};
    my @topSimilarity=sort {$b<=>$a} keys %{$hash{"similarity"}};
    my @topAccuracy=sort {$b<=>$a} keys %{$hash{"accuracy"}};
    my @topAlignLengths=sort {$b<=>$a} keys %{$hash{"alignLengths"}};
    my @topAlignScores=sort {$b<=>$a} keys %{$hash{"alignScores"}};
    for(my $i=1; $i<3; $i++){
        if(!defined $topCoverage[$i]){
            $topCoverage[$i]="NA";
            $hash{"coverage"}{$topCoverage[$i]}="NA";
        }
        if(!defined $topIdentity[$i]){
            $topIdentity[$i]="NA";
            $hash{"identity"}{$topIdentity[$i]}="NA";
        }
        if(!defined $topSimilarity[$i]){
            $topSimilarity[$i]="NA";
            $hash{"similarity"}{$topSimilarity[$i]}="NA";
        }
        if(!defined $topAccuracy[$i]){
            $topAccuracy[$i]="NA";
            $hash{"accuracy"}{$topAccuracy[$i]}="NA";
        }
        if(!defined $topAlignLengths[$i]){
            $topAlignLengths[$i]="NA";
            $hash{"alignLengths"}{$topAlignLengths[$i]}="NA";
        }
        if(!defined $topAlignScores[$i]){
            $topAlignScores[$i]="NA";
            $hash{"alignScores"}{$topAlignScores[$i]}="NA";
        }
        
    }
    my $mySC="SC:Z:".$topCoverage[0]."-".$hash{"coverage"}{$topCoverage[0]}.",".$topCoverage[1]."-".$hash{"coverage"}{$topCoverage[1]}.",".$topCoverage[2]."-".$hash{"coverage"}{$topCoverage[2]}; #SC: string of coverage
    my $mySI="SI:Z:".$topIdentity[0]."-".$hash{"identity"}{$topIdentity[0]}.",".$topIdentity[1]."-".$hash{"identity"}{$topIdentity[1]}.",".$topIdentity[2]."-".$hash{"identity"}{$topIdentity[2]}; #SI: string of identity
    my $mySS="SS:Z:".$topSimilarity[0]."-".$hash{"similarity"}{$topSimilarity[0]}.",".$topSimilarity[1]."-".$hash{"similarity"}{$topSimilarity[1]}.",".$topSimilarity[2]."-".$hash{"similarity"}{$topSimilarity[2]}; #SS: string of similarity
    my $mySA="SA:Z:".$topAccuracy[0]."-".$hash{"accuracy"}{$topAccuracy[0]}.",".$topAccuracy[1]."-".$hash{"accuracy"}{$topAccuracy[1]}.",".$topAccuracy[2]."-".$hash{"accuracy"}{$topAccuracy[2]}; #SA: string of accuracy
    my $mySL="SL:Z:".$topAlignLengths[0]."-".$hash{"alignLengths"}{$topAlignLengths[0]}.",".$topAlignLengths[1]."-".$hash{"alignLengths"}{$topAlignLengths[1]}.",".$topAlignLengths[2]."-".$hash{"alignLengths"}{$topAlignLengths[2]}; #SL: string of alignLengths
    my $myLS="LS:Z:".$topAlignScores[0]."-".$hash{"alignScores"}{$topAlignScores[0]}.",".$topAlignScores[1]."-".$hash{"alignScores"}{$topAlignScores[1]}.",".$topAlignScores[2]."-".$hash{"alignScores"}{$topAlignScores[2]}; #LS: string of alignScores
    $topAlignLengths[1]="0" if $topAlignLengths[1] eq "NA";
    my $myLD="LD:i:".($topAlignLengths[0]-$topAlignLengths[1]); #LD: length different
    return ($mySC, $mySI, $mySS, $mySA, $mySL, $myLS, $myLD);
}

sub getTagValue{
    my ($tag, $alignment)=@_;
    my $tagValue=$alignment;
    $tagValue=~s/.*$tag:.://;
    $tagValue=~s/\s.*//;
    return $tagValue;
}

sub output{
    my ($mySC, $mySI, $mySS, $mySA, $mySL, $myLS, $myLD)=&valueSort();
    my $myHC=scalar @outputString; #HC: hits counts
    $myHC="HC:i:".$myHC;
    foreach my $alignment (@outputString){
        $alignment=~s/\tSC:.*//;
        say join "\t",($alignment, $mySC, $mySI, $mySS, $mySA, $mySL, $myLS, $myLD, $myHC);
    }
    @outputString=();
}

sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script use to get the longest alignment (the highest coverage) for pacbio reads aligned by blasr, bam must be sorted by name
Usage: perl $scriptName *sorted.bam >OUTPUT
    If INPUT not specified, input from STDIN, output to STDOUT
    AS:	Alignment score
    XL:	Alignment length
    NM:	Number of mismatch
    XQ:	Reads length
    VC:	Value of coverage
    VI:	Value of identity
    VS:	Value of similarity
    VA:	Value of accuracy
    SC:	Top three coverage
    SI:	Top three identity
    SS:	Top three similarity
    SA:	Top three accuracy
    SL:	Top three align length
    LS:	Top three align score
    LD:	Length different between longest and second one

    -c --coverage   FLOAT   Report alignments only if their coverage isn't less than this [defalut: 0.75]
    -i --identiy    FLOAT   Report alignments only if their identiy isn't less than this [defalut: 0.6]
    -a --accuracy   FLOAT   Report alignments only if their accuracy isn't less than this [defalut: 0.0]
    -l --longest    BOOL    Report the longest alignments, no matter how is the value of coverage and identity
    -h --help               Print this help information screen
HELP
    exit(-1);
}