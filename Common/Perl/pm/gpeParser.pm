#!/usr/bin/perl -w
package gpeParser;

use strict;
use 5.010;
require Exporter;
use List::Util qw/sum/;
use DBI;

#our @ISA = qw/Exporter/;
#our @EXPORT = qw/getCDSLength/;
1;

######################################################
#Sub procedures for geting index
######################################################
sub CDSStartIndex(){
    my ($CDSStart, $exonStarts, $exonEnds)=@_;
    for(my $exonSIndex = 0; $exonSIndex <= scalar(@$exonStarts)-1; $exonSIndex++ ){
        if($CDSStart <= $exonEnds->[$exonSIndex]){#arrive at the CDSStartIndex
            if ($CDSStart >= $exonStarts->[$exonSIndex]){
                return $exonSIndex;
            }else{
                return;
            }
        }
    }
}

sub CDSEndIndex(){
    my ($CDSEnd,$exonStarts,$exonEnds)=@_;
    for(my $exonEIndex = scalar(@$exonStarts)-1; $exonEIndex >=0; $exonEIndex-- ){
        if($CDSEnd >= $exonStarts->[$exonEIndex]){#arrive at the CDSEndIndex
            if ($CDSEnd <= $exonEnds->[$exonEIndex]){
                return $exonEIndex;
            }else{
                return;
            }
        }
    }
}

sub CDSIndex(){
    my ($CDSStart,$CDSEnd,$exonStarts,$exonEnds)=@_;
    my ($exonSIndex,$exonEIndex);
    for($exonSIndex=0;$exonSIndex < scalar(@$exonStarts);$exonSIndex++ ){
        if($CDSStart <= $exonEnds->[$exonSIndex]){#arrive at the CDSStartIndex
            return if ($CDSStart < $exonStarts->[$exonSIndex]);
            last;
        }
    }
    for($exonEIndex=scalar(@$exonStarts)-1;$exonEIndex >=0;$exonEIndex-- ){
        if($CDSEnd >= $exonStarts->[$exonEIndex]){#arrive at the CDSEndIndex
            return if ($CDSEnd > $exonEnds->[$exonEIndex]);
            last;
        }
    }
    return ($exonSIndex, $exonEIndex);
}



######################################################
#Sub procedures for geting length
######################################################
sub getExonsLength(){#debugged
    my($exonStarts, $exonEnds) = @_;
    my @exonLengths = map{$exonEnds->[$_]-$exonStarts->[$_]}(0..scalar(@$exonStarts)-1);
    return sum(@exonLengths);
}

sub getCDSLength(){#debugged
    my ($CDSStart, $CDSEnd, $exonStarts, $exonEnds)=@_;
    &getExonsLength( &getCDSExons($CDSStart, $CDSEnd, $exonStarts, $exonEnds) );
}

######################################################
#Sub procedures for geting blocks
######################################################
sub getCDSExons(){#debugged
    my ($CDSStart, $CDSEnd, $exonStarts, $exonEnds)=@_;
    my ($exonSIndex,$exonEIndex);
    for($exonSIndex=0; $exonSIndex < scalar(@$exonStarts); $exonSIndex++ ){
        if($CDSStart <= $exonEnds->[$exonSIndex]){#arrive at the CDSStartIndex
            return if ($CDSStart < $exonStarts->[$exonSIndex]);
            last;
        }
    }
    for($exonEIndex=scalar(@$exonStarts)-1;$exonEIndex >=0;$exonEIndex-- ){
        if($CDSEnd >= $exonStarts->[$exonEIndex]){#arrive at the CDSEndIndex
            return if ($CDSEnd > $exonEnds->[$exonEIndex]);
            last;
        }
    }
    my @CDSExonStarts=@{$exonStarts}[$exonSIndex..$exonEIndex];
    $CDSExonStarts[0]=$CDSStart;
    my @CDSExonEnds=@{$exonEnds}[$exonSIndex..$exonEIndex];
    $CDSExonEnds[-1]=$CDSEnd;   
    return (\@CDSExonStarts,\@CDSExonEnds);
}

sub getLeftUTRExons(){ # debugged
    my ($CDSStart, $exonStarts, $exonEnds) = @_;
    my $exonCount = @$exonStarts;
    for( my $exonSIndex = 0; $exonSIndex < $exonCount; $exonSIndex++){
	if($exonStarts->[$exonSIndex] <= $CDSStart && $CDSStart < $exonEnds->[$exonSIndex]){
	    my (@leftUTRExonStarts, @leftUTRExonEnds);
	    if($exonStarts->[$exonSIndex] == $CDSStart){
		$exonSIndex--;
		@leftUTRExonStarts = @{$exonStarts}[0..$exonSIndex];
		@leftUTRExonEnds = @{$exonEnds}[0..$exonSIndex];
	    }else{
		@leftUTRExonStarts = @{$exonStarts}[0..$exonSIndex];
		@leftUTRExonEnds = @{$exonEnds}[0..$exonSIndex];
		$leftUTRExonEnds[-1] = $CDSStart;
	    }
	    return (\@leftUTRExonStarts, \@leftUTRExonEnds);
	}elsif($CDSStart < $exonStarts->[$exonSIndex]){
	    return;
	}
    }
}

sub getRightUTRExons(){ #debugged
    my ($CDSEnd, $exonStarts, $exonEnds) = @_;
    my $maxIndex = @$exonStarts -1;
    for( my $i = $maxIndex; $i >= 0; $i--){
	if($exonStarts->[$i] < $CDSEnd && $CDSEnd <= $exonEnds->[$i]){
	    my (@UTRExonStarts, @UTRExonEnds);
	    if($exonEnds->[$i] == $CDSEnd){
		$i++;
		@UTRExonStarts = @{$exonStarts}[$i..$maxIndex];
		@UTRExonEnds = @{$exonEnds}[$i..$maxIndex];
	    }else{
		@UTRExonStarts = @{$exonStarts}[$i..$maxIndex];
		@UTRExonEnds = @{$exonEnds}[$i..$maxIndex];
		$UTRExonStarts[0] = $CDSEnd;
	    }
	    return (\@UTRExonStarts, \@UTRExonEnds);
	}elsif($CDSEnd > $exonEnds->[$i]){
	    return;
	}
    }
}

sub getUTRExons(){ # not debugged
    my ($CDSStart, $CDSEnd, $exonStarts, $exonEnds) = @_;
    my ($exonSIndex, $exonEIndex);
    for($exonSIndex =0; $exonSIndex < scalar(@$exonStarts); $exonSIndex++ ){
        if($CDSStart <= $exonEnds->[$exonSIndex]){#arrive at the CDSStartIndex
            return if ($CDSStart < $exonStarts->[$exonSIndex]);
            last;
        }
    }
    for($exonEIndex = scalar(@$exonStarts)-1; $exonEIndex >= 0; $exonEIndex-- ){
        if($CDSEnd >= $exonStarts->[$exonEIndex]){#arrive at the CDSEndIndex
            return if ($CDSEnd > $exonEnds->[$exonEIndex]);
            last;
        }
    }
    my @UTRExonStarts = @{$exonStarts}[0..$exonSIndex];
    my @UTRExonEnds = @{$exonEnds}[0..$exonSIndex];
    $UTRExonEnds[-1] = $CDSStart;
    
    push @UTRExonStarts, @{$exonStarts}[$exonEIndex..-1];
    push @UTRExonEnds, @{$exonEnds}[$exonEIndex..-1];
    $UTRExonStarts[$exonEIndex] = $CDSEnd;
    
    return (\@UTRExonStarts, \@UTRExonEnds);
}

sub getIntrons{
    my ($exonStarts, $exonEnds) = @_;
    my @exonStarts = @$exonStarts;
    my @exonEnds = @$exonEnds;
    shift @exonStarts;	#change to be intron ends
    pop @exonEnds;	#change to be intron starts
    return (\@exonEnds, \@exonStarts);
}

sub getNewExons{
    my ($oldTransExonStarts,$oldTransExonEnds,$newTransExonStarts,$newTransExonEnds)=@_;
    my @newExons=();
    for (my $newTransExonIndex = 0 ; $newTransExonIndex < scalar(@$newTransExonStarts) ; $newTransExonIndex++){
        my ($newTransExonStart,$newTransExonEnd) = ($newTransExonStarts->[$newTransExonIndex],$newTransExonEnds->[$newTransExonIndex]);
        if( $newTransExonEnd < $oldTransExonStarts->[0]){
            push @newExons,[$newTransExonStart,$newTransExonEnd];
        }else{
            for(my $oldTransExonIndex=0 ; $oldTransExonIndex < scalar(@$oldTransExonStarts) ; $oldTransExonIndex++ ){
                if( $newTransExonStart > $oldTransExonEnds->[$oldTransExonIndex] &&
                    (  $oldTransExonIndex == scalar(@$oldTransExonStarts)-1 ||#to the last exon of old transcript
                        $newTransExonEnd < $oldTransExonStarts->[$oldTransExonIndex+1]
                    )
                   ){
                    push @newExons,[$newTransExonStart,$newTransExonEnd];
                    last;
                }
            }
        }
    }
    return \@newExons;
}

sub getRefineExons{
    my ($oldTransExonStarts,$oldTransExonEnds,$newTransExonStarts,$newTransExonEnds)=@_;
    my @refinedExons=();
    for (my $newTransExonIndex = 0 ; $newTransExonIndex < scalar(@$newTransExonStarts) ; $newTransExonIndex++){
        my ($newTransExonStart,$newTransExonEnd) = ($newTransExonStarts->[$newTransExonIndex],$newTransExonEnds->[$newTransExonIndex]);
        next if scalar(@$oldTransExonStarts) == 1;#only one exon in old transcript
        if( $newTransExonStart < $oldTransExonEnds->[0] &&
            $oldTransExonStarts->[0] < $newTransExonEnd &&
            $newTransExonEnd < $oldTransExonStarts->[1] &&
            $newTransExonEnd != $oldTransExonEnds->[0]){#refined exon is the first exon of old transcript
            push @refinedExons,[$oldTransExonStarts->[0],$oldTransExonEnds->[0],$newTransExonStart,$newTransExonEnd];
        }elsif(
                $oldTransExonEnds->[-2] < $newTransExonStart &&
                $newTransExonStart < $oldTransExonEnds->[-1] &&
                $newTransExonStart != $oldTransExonStarts->[-1] &&
                $newTransExonEnd > $oldTransExonStarts->[-1]
                 ){#refined exon is the last exon of old transcript
            push @refinedExons,[$oldTransExonStarts->[-1],$oldTransExonEnds->[-1],$newTransExonStart,$newTransExonEnd];
        }else{#refined exon is the inner exon of old transcript
            #start to transverse from the second exon to last but one exon of old transcript 
            for(my $oldTransExonIndex=1 ; $oldTransExonIndex < scalar(@$oldTransExonStarts)-1 ; $oldTransExonIndex++ ){
                my ($oldTransExonStart,$oldTransExonEnd)=($oldTransExonStarts->[$oldTransExonIndex],$oldTransExonEnds->[$oldTransExonIndex]);
                if($oldTransExonEnds->[$oldTransExonIndex-1] < $newTransExonStart &&
                   $newTransExonStart < $oldTransExonEnd &&
                   $oldTransExonStart < $newTransExonEnd &&
                   $newTransExonEnd < $oldTransExonStarts->[$oldTransExonIndex+1] &&
                   !($oldTransExonStart == $newTransExonStart && $oldTransExonEnd == $newTransExonEnd)  ){
                    next if ($newTransExonIndex == 0 && $newTransExonStart > $oldTransExonStart);
                    next if ($newTransExonIndex == scalar(@$newTransExonStarts)-1 && $newTransExonEnd < $oldTransExonEnd);
                    push @refinedExons,[$oldTransExonStart,$oldTransExonEnd,$newTransExonStart,$newTransExonEnd];
                    last;
                }
            }
        }
    }
    return \@refinedExons;
}

####    Argument            Type        Description
#       blockStarts         array ref   a reference to an array with block starts
#       blockEnds           array ref   a reference to an array with block ends
####    Return              Type        Description
#       sewedBlockStarts    array ref   a reference to an array with merged block starts
#       sewedBlockEnds      array ref   a reference to an array with merged block ends
sub getSewedExon{
    my ($blockStarts, $blockEnds) = @_;
    my (@sewedBlockStarts, @sewedBlockEnds);
    my $i;
    for ($i=0; $i < @$blockStarts - 1; $i++){
        if($blockEnds->[$i] == $blockStarts->[$i+1]){
            push @sewedBlockStarts, $blockStarts->[$i];
            $i++ while($i < @$blockStarts - 1 && $blockEnds->[$i] == $blockStarts->[$i+1]);
            push @sewedBlockEnds, $blockEnds->[$i];
            return (\@sewedBlockStarts, \@sewedBlockEnds) if $i == @$blockStarts - 1;
        }else{
            push @sewedBlockStarts, $blockStarts->[$i];
            push @sewedBlockEnds, $blockEnds->[$i];
        }
    }
    push @sewedBlockStarts, $blockStarts->[-1];
    push @sewedBlockEnds, $blockEnds->[-1];
    return (\@sewedBlockStarts, \@sewedBlockEnds);
}

####    Argument        Type    Description
#       blockStarts1    string  block starts of the 1st transript separated by comma
#       blockEnd1       string  block end of the 1st transript separated by comma
#       blockStarts2    string  block starts of the 2nd transript separated by comma
#       blockEnd1       string  block ends of the 2nd transript separated by comma
####    Return          Type    Description
#       mergedTrans     array	an array with refs of mergedStarts array and mergedEnds array
sub getMergedTrans(){
    my ($blockStarts1, $blockEnds1, $blockStarts2, $blockEnds2) = @_;
    my @blockStarts1 = @{$blockStarts1};
    my @blockEnds1 = @{$blockEnds1};
    my @blockStarts2 = @{$blockStarts2};
    my @blockEnds2 = @{$blockEnds2};
    my (@mergedStarts, @mergedEnds);
    my %hash_mergedTran;
    my ($i,$j);
    my (@a, @b);
    $j=0; #j use to indicate the variable of @blockStarts1 in use
    foreach (@blockStarts1){
		if($j>0){
			if($blockStarts1[$j]<$blockEnds1[$j-1]){
				if($blockEnds1[$j]<=$blockEnds1[$j-1]){
					$blockStarts1[$j]=$blockStarts1[$j-1];
					$blockEnds1[$j]=$blockEnds1[$j-1];
					$j++;
					next;
				}else{
					$blockStarts1[$j]=$blockStarts1[$j-1];
				}
			}elsif($blockStarts1[$j]==$blockEnds1[$j-1]){
				$blockStarts1[$j]=$blockStarts1[$j-1];
			}else{
				$hash_mergedTran{$blockStarts1[$j-1]}=$blockEnds1[$j-1];
			}
		}
		my $mergedStart=$blockStarts1[$j];
		$i=0; #i use to indicate the variable of @blockStarts2 in use
		if(@blockStarts2){
			foreach (@blockStarts2){
				if($_<$mergedStart){
					if($blockEnds2[$i]<$mergedStart){
						$hash_mergedTran{$blockStarts2[$i]}=$blockEnds2[$i];
					}elsif($blockEnds2[$i]==$mergedStart){
						$blockStarts1[$j]=$blockStarts2[$i];
					}else{
						if($blockEnds2[$i]<=$blockEnds1[$j]){
							$blockStarts1[$j]=$blockStarts2[$i];
						}else{
							$blockStarts1[$j]=$blockStarts2[$i];
							$blockEnds1[$j]=$blockEnds2[$i];
						}
					}
				}elsif($_>=$mergedStart){
					if($_>$blockEnds1[$j]){
						last;
					}elsif($blockEnds2[$i]>=$blockEnds1[$j]){
						$blockEnds1[$j]=$blockEnds2[$i];
					}
				}
				$i++;											
			}
			if($i>0){
				for(0..$i-1){
					shift @blockStarts2;
					shift @blockEnds2;
				}
			}
		}else{
			if($j<scalar @blockStarts1){
				$hash_mergedTran{$blockStarts1[$j]}=$blockEnds1[$j];
			}
		}
		$j++;	
    }
	$j=scalar @blockStarts1;
	$hash_mergedTran{$blockStarts1[$j-1]}=$blockEnds1[$j-1];
    if(@blockStarts2){
		my $i=scalar @blockStarts2;
		foreach(0..$i-1){
			$hash_mergedTran{$blockStarts2[$_]}=$blockEnds2[$_];
		}
    }
    @a= sort {$a <=> $b } keys %hash_mergedTran;
    @b= sort {$a <=> $b } values %hash_mergedTran;
    return (\@a,\@b);    
}


sub getOverlap(){
    my ($blockStarts1, $blockEnds1, $blockStarts2, $blockEnds2) = @_;
    my @blockStarts1 = @$blockStarts1;
    my @blockEnds1 = @$blockEnds1;
    my @blockStarts2 = @$blockStarts2;
    my @blockEnds2 = @$blockEnds2;
    my $j = 0;
    my (@overlapStarts, @overlapEnds);
I:  for(my $i = 0; $i <= $#blockStarts1; $i++){
	my ($blockStart1, $blockEnd1) = ($blockStarts1[$i], $blockEnds1[$i]);
	for(;$j <= $#blockStarts2; $j++){
	    my ($blockStart2, $blockEnd2) = ($blockStarts2[$j], $blockEnds2[$j]);
	    if($blockStart2 >= $blockEnd1){ # 2 is at right of 1
		next I;
	    }elsif($blockEnd2 <= $blockStart1){
		next;
	    }else{ # intersect
		if($blockStart2 <= $blockStart1){
		    if($blockEnd2 <= $blockEnd1){
			push @overlapStarts, $blockStarts1[$i];
			push @overlapEnds, $blockEnds2[$j];
		    }else{
			push @overlapStarts, $blockStarts1[$i];
			push @overlapEnds, $blockEnds1[$i];
			next I;
		    }
		}else{
		    if($blockEnd2 <= $blockEnd1){
			push @overlapStarts, $blockStarts2[$j];
			push @overlapEnds, $blockEnds2[$j];
		    }else{
			push @overlapStarts, $blockStarts2[$j];
			push @overlapEnds, $blockEnds1[$i];
			next I;
		    }
		}
	    }
	}
    }
    return (\@overlapStarts, \@overlapEnds);
}

sub getOverlapLength(){
    &getExonsLength(&getOverlap(@_));
}
sub getSubtract(){
    my ($blockStarts1, $blockEnds1, $blockStarts2, $blockEnds2) = @_;
    my @blockStarts1 = @$blockStarts1;
    my @blockEnds1 = @$blockEnds1;
    my @blockStarts2 = @$blockStarts2;
    my @blockEnds2 = @$blockEnds2;
    my $j = 0;
    my $offset =0;
I:  for(my $i = 0; $i <= $#blockStarts1; $i++){
	my ($blockStart1, $blockEnd1) = ($blockStarts1[$i+$offset], $blockEnds1[$i+$offset]);
	for(;$j <= $#blockStarts2; $j++){
	    my ($blockStart2, $blockEnd2) = ($blockStarts2[$j], $blockEnds2[$j]);
	    if($blockStart2 >= $blockEnd1){ # 2 is at right of 1
		last;
	    }elsif($blockEnd2 <= $blockStart1){
		next;
	    }else{ # intersect
		if($blockStart2 <= $blockStart1){
		    if($blockEnd2 < $blockEnd1){
			$blockStarts1[$i+$offset] = $blockEnd2;
			#$i--; # -- in advance and then ++ (by $i++) to stay at current block1
		    }else{
			splice @blockStarts1, $i, 1;
			splice @blockEnds1, $i, 1;
			$offset--;
			next I;
		    }
		}else{
		    if($blockEnd2 < $blockEnd1){
			splice @blockStarts1, $i+1, 0, $blockEnd2;
			splice @blockEnds1, $i, 0, $blockStart2;
			#$i--;  # -- in advance and then ++ (by $i++) to stay at current block1
			$offset++;
		    }else{
			$blockEnds1[$i+$offset] = $blockStart2;
			next I;
		    }
		}
	    }
	    #last;
	}
    }
    return (\@blockStarts1, \@blockEnds1);
}


sub getFwdExonsByRelCoor(){
    my ($exonStarts, $exonEnds, $relStart, $relEnd) = @_;
    my $exonCount = @$exonStarts;
    my ($start, $end, $startI) = (0);
    my (@subExonStarts, @subExonEnds);
    for (my $i = 0; $i < $exonCount; $i++){
	$start += $exonEnds->[$i] - $exonStarts->[$i];
	if($start > $relStart){
	    $startI = $i;
	    $end = $start - ($exonEnds->[$i] - $exonStarts->[$i]);
	    $start = $exonEnds->[$i] - ($start - $relStart);
	    for (; $i < $exonCount; $i++){
		$end += $exonEnds->[$i] - $exonStarts->[$i];
		if($end >= $relEnd){
		    $end = $exonEnds->[$i] - ($end - $relEnd);
		    @subExonStarts = @$exonStarts[$startI..$i];
		    @subExonEnds = @$exonEnds[$startI..$i];
		    $subExonStarts[0] = $start;
		    $subExonEnds[-1] = $end;
		    return (\@subExonStarts, \@subExonEnds);
		}
	    }
	}
    }    
}

sub getRevExonsByRelCoor(){
    my ($exonStarts, $exonEnds, $relStart, $relEnd) = @_;
    my $exonCount = @$exonStarts;
    my ($end, $start, $endI) = (0);
    my (@subExonStarts, @subExonEnds);
    for (my $i = $exonCount -1; $i >= 0; $i--){
	$end += $exonEnds->[$i] - $exonStarts->[$i];
	if($end > $relStart){
	    $endI = $i;
	    $start = $end - ($exonEnds->[$i] - $exonStarts->[$i]);
	    $end = $exonStarts->[$i] + ($end - $relStart);
	    for (; $i >= 0; $i--){
		$start += $exonEnds->[$i] - $exonStarts->[$i];
		if($start >= $relEnd){
		    $start = $exonStarts->[$i] + ($start - $relEnd);
		    @subExonStarts = @$exonStarts[$i..$endI];
		    @subExonEnds = @$exonEnds[$i..$endI];
		    $subExonStarts[0] = $start;
		    $subExonEnds[-1] = $end;
		    return (\@subExonStarts, \@subExonEnds);
		}
	    }
	}
    }    
}

######################################################
#Sub procedures for frame
######################################################

####    Argument    Type    Description
#       CDSStart    scalar  CDS start
#       CDSEnd      scalar  CDS end
#       strand      scalar  strand
#       exonStarts  scalar  exon starts (seperated by comma)
#       exonEnds    scalar  exon ends (seperated by comma)

####    Return      Type    Description
#       exonFrames  array   exon frames
sub getExonFrames{
    my ($cds_s, $cds_e, $strand, $exonSs, $exonEs) = @_;
    my @exons_s=split ",", $exonSs;
    my @exons_e=split ",", $exonEs;
    my $exon_count=scalar @exons_s;
    my ($exon_Frame,@exon_Frames);
    my $i=0;
    if($strand=~/\+/){
        while($i<$exon_count){
            if($exons_s[$i] < $cds_s){ 
                $exon_Frame = ($exons_e[$i] <= $cds_s) ? -1: 0;
            }elsif($exons_s[$i] == $cds_s){
                $exon_Frame=0;
            }else{
                if($exons_s[$i] < $cds_e){
                    if($exons_s[$i-1]<$cds_s){
                        $exon_Frame = ($exons_e[$i-1]-$cds_s)%3;
                    }else{
			$exon_Frame = ($exons_e[$i-1]-$exons_s[$i-1] + $exon_Frames[$i-1]-3)%3;
                    }
                }else{
                    $exon_Frame = -1;
                }
            }
            push @exon_Frames, $exon_Frame;
            $i++;
        }
    }else{
        $i = $exon_count - 1;
        while($i >= 0){
            if($exons_e[$i] > $cds_e){  
                $exon_Frame = $exons_s[$i] >= $cds_e ? -1 : 0;
            }elsif($exons_e[$i]==$cds_e){ 
                $exon_Frame = ($cds_s == $cds_e) ? -1 : 0;
            }else{
                if($exons_e[$i] > $cds_s){ 
                    if($exons_e[$i+1]>$cds_e){
                        $exon_Frame = ($cds_e - $exons_s[$i+1])%3;
                    }else{
			$exon_Frame = ($exons_e[$i+1] - $exons_s[$i+1] + $exon_Frames[$exon_count-$i-2] - 3) % 3;
                    }
                }else{
                    $exon_Frame = -1;
                }
            }
            push @exon_Frames, $exon_Frame;  
            $i--;
        }
        @exon_Frames = reverse @exon_Frames;
    }
    return \@exon_Frames;
}

####    Argument    Type    Description
#       gpeRecord   string  a line of gpe file

####    Return      Type    Description
#       exonFrames  array   exon frames
sub getExonFrames2{
    my ($gpeRecord) = @_;
    my ($strand, $cds_s, $cds_e, $exonSs, $exonEs) = (split "\t", $gpeRecord)[2,5,6,8,9];
    return getExonFrames($cds_s, $cds_e, $strand, $exonSs, $exonEs);
}

######################################################
#Sub procedures for others
######################################################

sub getSizes{
    my ($blockStarts, $blockEnds) = @_;
    map{$blockEnds->[$_] - $blockStarts->[$_]}0..(@$blockStarts-1);
}

####    Argument    Type    Description
#       fileHandle  FH      A file handle of gpe file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, name] as values
sub buildHash{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($name, $chr, $strand, $start, $end) = @fields[0..4];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, $name];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, $name] ];
        }        
    }
    while(my ($chr, $chrV) = each %hash){
        for my $strand (keys %$chrV){
            my @sortedStrandV = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$chrV->{$strand}};
            $chrV->{$strand} = \@sortedStrandV;
        }
    }
    return %hash;
}


####    Argument    Type    Description
#       fileHandle  FH      A file handle of gpe file

####    Return      Type    Description
#       hash        hash    A hash with chr, strand as key and with array containing [start, end, line] as values
sub buildHash2{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($chr, $strand, $start, $end) = @fields[1..4];
        if(exists $hash{$chr}{$strand}){
            push @{$hash{$chr}{$strand}}, [$start, $end, (join "\t", @fields)];
        }else{
            $hash{$chr}{$strand} = [ [$start, $end, (join "\t", @fields)] ];
        }        
    }
    while(my ($chr, $chrV) = each %hash){
        for my $strand (keys %$chrV){
            my @sortedStrandV = sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1]}@{$chrV->{$strand}};
            $chrV->{$strand} = \@sortedStrandV;
        }
    }
    return %hash;
}

sub buildHash3{
    my ($fh, $bin) = @_;
    my %hash;
    while(<$fh>){
        chomp;
        my @fields = split "\t";
        shift @fields if defined $bin;
        my ($name, $chr, $start, $end) = @fields[0, 1, 3, 4];
	if(exists $hash{"$chr"}){
	    push @{$hash{"$chr"}}, [$start, $end, (join "\t", @fields)];
	}else{
	    $hash{"$chr"} = [ [$start, $end, (join "\t", @fields)] ];
	}        
    }
    for my $chr(keys %hash){
	my @sorted = sort{$a->[0]<=>$b->[0]}(@{$hash{$chr}});
	$hash{"$chr"} = \@sorted;
    }
    return %hash;
}

sub getTransName{#no finished
    my ($start, $end, $hash, $chr, $strand) = @_;
    my @trans = @{$hash->{$chr}{$strand}};
    my ($left, $right) = (0, $#trans);
    my $mid = int( ($left + $right)/2 );;
    while($mid != $left && $start != $trans[$mid]->[0]){        
        if( $start < $trans[$mid]->[0]){
            $right = $mid;
        }elsif( $start > $trans[$mid]->[0]){
            $left = $mid;
        }
        $mid = int( ($left + $right)/2 );
    }
    # To this statement, $tmp[$mid] <= $start< $tmp[$right]
    my @transName;
    for( my $i = $mid; $i >= 0; $i--){
        if($start < $trans[$i]->[1] && $end > $trans[$i]->[0]){
            push @transName, $trans[$i]->[2];
        }
    }
    @transName = reverse @transName;
}

sub isTransRead{
    my ($readS, $readE, $exonStarts, $exonEnds, $slop) = @_;
    $slop = 0 unless defined $slop;
    my ($isStartInIntron,$isEndInIntron) = (0, 0);
    my ($exonSIndex, $exonEIndex);
    for($exonSIndex = 0; $exonSIndex < scalar(@$exonStarts); $exonSIndex++){#find exonSIndex of read
         if($readS <= $exonEnds->[$exonSIndex]){
            $isStartInIntron = 1 if($readS < $exonStarts->[$exonSIndex] - $slop);#need left extending offset
            last;
         }
    }
    for($exonEIndex = $exonSIndex; $exonEIndex < scalar(@$exonStarts); $exonEIndex++){#find exonEIndex of read
         if($readE <= $exonEnds->[$exonEIndex] + $slop){#need right extending offset
             $isEndInIntron = 1 if($readE < $exonStarts->[$exonEIndex]);
             last;
         }
    }
    return ($isStartInIntron == 0 && $isEndInIntron == 0); #read is two-side in exon
}


sub isCanonicalSite{
    my ($chr, $start, $strand, $type, $seqHash) = @_;
    return unless exists $seqHash->{$chr};
    if( $strand eq '+'){
        if($type eq 'donor'){
            return 1 if (substr $seqHash->{$chr}, $start, 2) =~/(G[TC])|(AT)/g;
        }else{ # acceptor
            return 1 if (substr $seqHash->{$chr}, $start, 2) =~/A[GC]/g;
        }
    }else{
        if($type eq 'donor'){
            return 1 if (substr $seqHash->{$chr}, $start, 2) =~/([AG]C)|(AT)/g;
        }else{ # acceptor
            return 1 if (substr $seqHash->{$chr}, $start, 2) =~/[CG]T/g;
        }
    }
}

sub mysqlCreateAndLoad{
    my ($gpeFile, $bin) = @_;
    my $dsn = "DBI:mysql:perl;mysql_socket=/rd1/mysql/mysql.sock;mysql_local_infile=1";
    my $user = 'perl';
    my $dbh = DBI->connect($dsn, $user, "", { RaiseError => 1 }) || die "Could not connect to database: $DBI::errstr";

    $dbh->do('DROP TABLE IF EXISTS gpe');
    my $mysqlStat = 'CREATE TABLE gpe (';
    $mysqlStat .= 'bin INT,' if defined $bin;
    $mysqlStat .= 'name VARCHAR(50),
                    chr VARCHAR(10),
                    strand VARCHAR(2),
                    start INT, end INT,
                    cdsStart INT, cdsEnd INT,
                    exonCount INT,
                    exonStarts TEXT, exonEnds TEXT,
                    score INT,
                    name2 VARCHAR(50),
                    cdsStartStat VARCHAR(10), cdsEndStat VARCHAR(10),
                    exonFrames TEXT,
                    INDEX(name),
                    INDEX(chr),
                    INDEX(strand),
                    INDEX(start),
                    INDEX(end)
                    ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
                ';
    $dbh->do($mysqlStat);
    $mysqlStat = "load data local infile '$gpeFile' into table gpe(";
    $mysqlStat .= 'bin,' if defined $bin;
    $mysqlStat .= 'name, chr, strand, start, end, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat,
                    exonFrames)';
    $dbh->do($mysqlStat);
    $dbh;
}

sub mysqlSelect{
    my ($dbh, $query) = @_;
    my $sth = $dbh->prepare($query);
    die "Error:" . $sth->errstr . "\n" unless ($sth->execute());
    $sth;
}

sub mysqlDrop{
    my ($dbh) = @_;
    $dbh->do('DROP TABLE gpe');
}


sub getConsensusIntronN(){
    my ($exonStarts1, $exonEnds1, $exonStarts2, $exonEnds2, $offset) = @_;
    $offset = 0 unless defined $offset;
    my ($intronStarts1, $intronEnds1) = getIntrons($exonStarts1, $exonEnds1);
    my ($intronStarts2, $intronEnds2) = getIntrons($exonStarts2, $exonEnds2);
    my ($j, $consensusN) = (0, 0);
I:  for(my $i = 0; $i < @$intronStarts1; $i++){
	for(; $j < @$intronStarts2; $j++){
	    next I if($intronStarts2->[$j] > $intronEnds1->[$i]);
	    if($intronStarts1->[$i] - $offset <= $intronStarts2->[$j] && $intronStarts2->[$j] <= $intronStarts1->[$i] + $offset &&
		$intronEnds1->[$i] - $offset <= $intronEnds2->[$j] && $intronEnds2->[$j] <= $intronEnds1->[$i] + $offset){
		$consensusN++;
	    }
	}
    }
    return $consensusN;
}

