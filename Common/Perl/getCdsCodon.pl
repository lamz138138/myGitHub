#!/usr/bin/env perl
use strict;
use 5.010;
use Getopt::Long;
use File::Basename;
use lib dirname $0;

my ($help, $bin, $chr_input, $pos_left, $pos_right);
$pos_left=1;
$pos_right=100;
GetOptions(
	'b|bin'				=>	\$bin,
	'c|chr_input=s'		=>	\$chr_input,
	's|pos_left=i'		=>	\$pos_left,
	'e|pos_right=i'		=>	\$pos_right,
    'h|help' 			=>	\$help
)||usage(); 
usage () if defined $help;

$ARGV[0]='-' unless defined $ARGV[0];
open IN,"$ARGV[0]" or die "Can't open $ARGV[0]:$!";


while(<IN>){
	chomp;
	my @input=split "\t";
	shift @input if defined $bin;
	my ($name_1, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds, $id, $name_2, $cdsStartStat, $cdsEndStat, $exonFrames)=@input[0..14];
	next if $cdsStartStat ne "cmpl" || $cdsEndStat ne "cmpl";
	next if $chr ne $chr_input;
	my $gpeData=join "\t",($name_1, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds, $id, $name_2, $cdsStartStat, $cdsEndStat, $exonFrames);
	my $temp_pos_left=$pos_left;
	my $temp_pos_right=$pos_right;
	if($pos_left >= $cdsEnd || $pos_right <= $cdsStart){
		&outputWarnning(join "\t",($gpeData, $chr_input.":".$pos_left."-".$pos_right));
	}else{
		my @ExonStarts=split ",",$exonStarts;
		my @ExonEnds=split ",",$exonEnds;
		my @ExonFrames=split ",",$exonFrames;
		my ($index_left, $index_right, $index_cdsStart, $index_cdsEnd);
		for(my $i=0; $i<$exonCount; $i++){
			if(!defined $index_left){
				if($pos_left<$ExonStarts[$i]){
					$index_left=$i;
					$pos_left=$ExonStarts[$i];
				}elsif($pos_left >= $ExonStarts[$i] && $pos_left <= $ExonEnds[$i]){
					$index_left=$i
				}
			}
			if(!defined $index_right){
				if($pos_right > $ExonStarts[$i] && $pos_right <= $ExonEnds[$i]){
					$index_right=$i;
				}elsif($i<$exonCount-1 && $pos_right > $ExonEnds[$i] && $pos_right <= $ExonStarts[$i+1]){
					$index_right=$i;
					$pos_right=$ExonEnds[$i];
				}elsif($i==$exonCount-1 && $pos_right > $ExonEnds[$i]){
					$index_right=$i;
					$pos_right=$ExonEnds[$i];
				}
			}
			$index_cdsStart=$i if $cdsStart >= $ExonStarts[$i] && $cdsStart <= $ExonEnds[$i];
			$index_cdsEnd=$i if $cdsEnd >= $ExonStarts[$i] && $cdsEnd <= $ExonEnds[$i];
			last if defined $index_left && defined $index_right && defined $index_cdsStart && defined $index_cdsEnd;
		}
		if($pos_left > $pos_right){
			say STDERR "Warnning: this region is in intron:";
			say STDERR join "\t",($gpeData, $chr.":".$temp_pos_left."-".$temp_pos_right);
			last;
		}
		for(my $i=0; $i<$exonCount; $i++){
			$index_left=$i if $pos_left >= $ExonStarts[$i] && $pos_left <= $ExonEnds[$i];
			$index_right=$i if $pos_right >= $ExonStarts[$i] && $pos_right <= $ExonEnds[$i];
			$index_cdsStart=$i if $cdsStart >= $ExonStarts[$i] && $cdsStart <= $ExonEnds[$i];
			$index_cdsEnd=$i if $cdsEnd >= $ExonStarts[$i] && $cdsEnd <= $ExonEnds[$i];
			last if defined $index_left && defined $index_right && defined $index_cdsStart && defined $index_cdsEnd;
		}
		if($pos_left < $cdsStart){
			$pos_left=$cdsStart;
			$index_left=$index_cdsStart;
		}
		if($pos_right > $cdsEnd){
			$pos_right=$cdsEnd;
			$index_right=$index_cdsEnd;
		}
		if($strand eq "+"){
			# leftBase  标记坐标需移动多少bp，&getPos计算移动后的坐标
			if($pos_left != $cdsStart){
				my $leftBase;
				if($index_left == $index_cdsStart){
					$leftBase=($pos_left-$cdsStart)%3;
				}else{
					$leftBase=($pos_left-$ExonStarts[$index_left]-3+$ExonFrames[$index_left])%3;
				}
				$leftBase=3-$leftBase if $leftBase >0;
				($pos_left, $index_left)=&getPos($strand, $exonStarts, $exonEnds, $leftBase, $index_left, $pos_left);
			}
			&getCdsRegion($gpeData, $index_left, $index_right, $pos_left, $pos_right);
		}else{
			if($pos_right != $cdsEnd){
				my $leftBase;
				if($index_right == $index_cdsEnd){
					$leftBase=($cdsEnd-$pos_right)%3;
				}else{
					$leftBase=($ExonEnds[$index_right]-$pos_right-3+$ExonFrames[$index_right])%3;
				}
				$leftBase=3-$leftBase if $leftBase >0;
				($pos_right, $index_right)=&getPos($strand, $exonStarts, $exonEnds, $leftBase, $index_right, $pos_right);
			}
			&getCdsRegion($gpeData, $index_left, $index_right, $pos_left, $pos_right);
		}
	}
}

sub outputWarnning{
	my ($stringData)=@_;
	say STDERR "Warnning: this region is outside of cds region:";
	say STDERR $stringData;
	last;
}

sub getPos{
	my ($strand, $exonStarts, $exonEnds, $leftBase, $index, $pos)=@_;
	my @ExonStarts=split ",",$exonStarts;
	my @ExonEnds=split ",",$exonEnds;
	my $tempIndex=$index;
	if($strand eq "+"){
		my $tempLength=$ExonEnds[$index]-$pos;
		while($tempLength<=$leftBase){
			$leftBase-=$tempLength;
			$index++;
			$tempLength=$ExonEnds[$index]-$ExonStarts[$index];
		}
		if($tempIndex==$index){
			$pos+=$leftBase;
		}else{
			$pos=$ExonStarts[$index]+$leftBase;
		}
	}else{
		my $tempLength=$pos-$ExonStarts[$index];
		while($tempLength<=$leftBase){
			$leftBase-=$tempLength;
			$index--;
			$tempLength=$ExonEnds[$index]-$ExonStarts[$index];
		}if($tempIndex==$index){
			$pos-=$leftBase;
		}else{
			$pos=$ExonEnds[$index]-$leftBase;
		}
	}
	return ($pos, $index);
}

sub getCdsRegion{
	my ($gpeData, $index_left, $index_right, $pos_left, $pos_right)=@_;
	my ($name_1, $chr, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds, $id, $name_2, $cdsStartStat, $cdsEndStat, $exonFrames)=split "\t", $gpeData;
	my @ExonStarts=split ",",$exonStarts;
	my @ExonEnds=split ",",$exonEnds;
	my @ExonFrames=split ",",$exonFrames;
	my $length;
	if($index_left != $index_right){
		$length=$ExonEnds[$index_left]-$pos_left;
		$length+=$pos_right-$ExonStarts[$index_right];
	}else{
		$length=$pos_right-$pos_left;
	}
	for(my $i=$index_left+1; $i<$index_right; $i++){
		$length+=$ExonEnds[$i]-$ExonStarts[$i];
	}
	my $leftBase=$length%3;
	if($strand eq "+"){
		my $tempLength=$pos_right-$ExonStarts[$index_right];
		my $tempIndex=$index_right;
		while($tempLength<=$leftBase){
			$leftBase-=$tempLength;
			$index_right--;
			$tempLength=$ExonEnds[$index_right]-$ExonStarts[$index_right];
		}
		if($tempIndex==$index_right){
			$pos_right-=$leftBase;
		}else{
			$pos_right=$ExonEnds[$index_right]-$leftBase;
		}
	}else{
		my $tempLength=$ExonEnds[$index_left]-$pos_left;
		my $tempIndex=$index_left;
		while($tempLength<=$leftBase){
			$leftBase-=$tempLength;
			$index_left++;
			$tempLength=$ExonEnds[$index_left]-$ExonStarts[$index_left];
		}
		if($tempIndex==$index_left){
			$pos_left+=$leftBase;
		}else{
			$pos_left=$ExonStarts[$index_left]+$leftBase;
		}
	}
	my $output_starts=$pos_left;
	my $output_ends;
	my $output_frames;
	if($index_left != $index_right){
		$output_ends=$ExonEnds[$index_left];
		if($strand eq "+"){
			$output_frames="0";
		}else{
			$output_frames=$ExonFrames[$index_left];
		}
		for(my $i=$index_left+1; $i<$index_right; $i++){
			$output_starts.=",".$ExonStarts[$i];
			$output_ends.=",".$ExonEnds[$i];
			$output_frames.=",".$ExonFrames[$i];
		}
		$output_starts.=",".$ExonStarts[$index_right];
		$output_ends.=",".$pos_right;
		if($strand eq "+"){
			$output_frames.=",".$ExonFrames[$index_right];
		}else{
			$output_frames.=",0";
		}
	}else{
		$output_ends=$pos_right;
		$output_frames="0";
	}
	if($pos_right-$pos_left>2){
		say join "\t",($name_1, $chr, $strand, $pos_left, $pos_right, $pos_left, $pos_right, $index_right-$index_left+1, $output_starts, $output_ends, $id, $name_2, $cdsStartStat, $cdsEndStat, $output_frames);
	}else{
		say STDERR ("Warnning: this region is too short:");
		say STDERR join "\t",@_;
	}
}


sub usage{
    my $scriptName=basename $0;
print <<HELP;
This script was used to get cds codon in specific region
Usage: perl $scriptName -s 1 -e 100 *gpe >output
	if gpe is not specified, input from STDIN, output to STDOUT

	-b 	--bin		Have bin column
	-s 	--pos_left	The start position [default: 1]
	-e 	--pos_right 	The end position [default: 100]
	-h 	--help 		Print this help information
HELP
	exit(-1);
}