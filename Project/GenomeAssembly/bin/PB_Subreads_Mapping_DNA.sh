#!/bin/bash

usage(){
	cat << EOF
Description:
	This script was used to get CCS reads and map reads to reference

Usage:
	PB_Subreads_Mapping_DNA.sh -r reference -w workDir
Debug:
	The split module (-m) need to be tested
Options:
	-r      STR     (Optional) The path to the reference [Default: /rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis/references/rheMac2] 
	-w      STR     (Optional) The path of the work directory [Default: .]
	-d	STR	(Optional) The path contain the *bax.h5 files [Defalut: ./rawData]
	-s	STR	(Optional) The reference specie [Default: rheMac2]
	-t	INT	(Optional) Number of thread [Default: 10]
	-c	INT	(Optional) Split the work to this number of chunk (-s must be assigned) [default: 10]
	-m		(Optional) Mapping reads to reference
	-s		(Optional) Split the work to small chunk
	-v		(Optional) Call variants (-m must be assgined)
	-h              Print this information
EOF
	exit 0
}

#[ $1 ] || usage
[ $# -eq 0 ] && usage

workPath="$PWD"
referencePath="/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis/references/rheMac2"
dataPath="$PWD/rawData"
specie="rheMac2"
threadNum="10"
chunkNum="10"
while getopts "hd:r:s:t:w:mvc:s" OPTION
  do
    case $OPTION in 
	h) usage;;
	d) dataPath=$OPTARG;;
	t) threadNum=$OPTARG;;
	r) referencePath=$OPTARG;;
	s) specie=$OPTARG;;
	m) mappingRef=1;;
	v) variantsCall=1;;
	c) chunkNum=$OPTARG;;
	s) splitChunk=1;;
	w) workPath=$OPTARG;;
    esac
  done
shift $((OPTIND - 1))

SEYMOUR_HOME=/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936
source $SEYMOUR_HOME/etc/setup.sh

mkdir -p "$workPath"/result
mkdir -p "$workPath"/data
mkdir -p "$workPath"/tmpdir

# 1. P_Fetch Module
# 1) toFofn
find "$dataPath" -name *bax.h5 >"$workPath"/input.fofn
# 2) overviewRpt && adapterRpt (Background) 
# a) overviewRpt
overview_report.py --debug "$workPath"/input.fofn "$workPath"/results/overview.json &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o overview.html -rules "$workPath"/results/.rules_overview.xml &
# b) adapterRpt
pbreport.py adapter "$workPath"/results filter_reports_adapters.json "$workPath"/input.fofn &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o filter_reports_adapters.html -rules "$workPath"/results/.rules_filter_reports_adapters.xml &

# 2. P_Filter Module
# 1) filter.plsFofn.Scatter (can split to many chunk)
if [ -z "$splitChunk" ]; then
  ln -s "$workPath"/input.fofn > "$workPath"/input.chunk.fofn
else
  TOTAL_LINES=`cat "$workPath"/input.fofn | wc -l`
  for((i=0;i<$TOTAL_LINES;i++))
    do
      awk -v chunkNum=$chunkNum -v i=$i "($TOTAL_LINES-NR+1)%$chunkNum==$i" "$workPath"/input.fofn > "$workPath"/input.chunk"$i"of"$chunkNum".fofn
    done
fi
# 2) filter_*of* (run for each chunk)
ls "$workPath"/data/filtered_regions.chunk*.fofn | while read file
  do
    sample=$( echo $file | sed 's#filtered_regions.##;s#.fofn##' )
    filter_plsh5.py --debug --filter='MinReadScore=0.7500,MinSRL=50,MinRL=50' --trim='True' --outputDir="$workPath"/data/filtered_regions --outputSummary="$workPath"/data/filtered_summary."$sample".csv --outputFofn="$workPath"/data/filtered_regions."$sample".fofn "$workPath"/input."$sample".fofn
  done
# 3) subreads_*of* (run for each chunk)
ls "$workPath"/data/filtered_regions.chunk*.fofn | while read file
  do
    sample=$( echo $file | sed 's#filtered_regions.##;s#.fofn##' )
    pls2fasta "$workPath"/input."$sample".fofn "$workPath"/data/filtered_subreads."$sample".fasta -trimbyregion -regiontable "$workPath"/data/filtered_regions."$sample".fofn
    pls2fasta "$workPath"/input."$sample".fofn "$workPath"/data/filtered_subreads."$sample".fastq -trimbyregion -regiontable "$workPath"/data/filtered_regions."$sample".fofn -fastq
  done
# 4) filter.rgnfofn.gather (merge chunk)
cat "$workPath"/data/filtered_regions.chunk*.fofn > "$workPath"/data/filtered_regions.fofn
syncpermoviefofn.py --debug "$workPath"/data/filtered_regions.fofn "$workPath"/input.fofn > "$workPath"/data/filtered_regions.fofn.tmp
mv "$workPath"/data/filtered_regions.fofn.tmp "$workPath"/data/filtered_regions.fofn
# 5) filter.summary.gather (merge chunk)
cat "$workPath"/data/filtered_summary.chunk*.csv > "$workPath"/data/filtered_summary.csv
# 6) subreads.subreadfastq.gather (merge chunk, Background)
cat "$workPath"/data/filtered_subreads.chunk*.fastq > "$workPath"/data/filtered_subreads.fastq &
# 7) subreads.subreads.gather (merge chunk, Background)
cat "$workPath"/data/filtered_subreads.chunk*.fasta > "$workPath"/data/filtered_subreads.fasta &
# 8) subreadsummary
filter_subread_summary.py "$workPath"/data/filtered_regions.fofn --output="$workPath"/data/filtered_subread_summary.csv --debug

# 3. P_FilterReports
# 1) subreadrpt && statsrpt && loadingrpt (background)
# a) subreadrpt
filter_subread.py --debug --report="$workPath"/results/filter_reports_filter_subread_stats.json --output="$workPath"/results "$workPath"/data/filtered_subread_summary.csv &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o filter_reports_filter_subread_stats.html -rules "$workPath"/results/.rules_filter_reports_filter_subread_stats.xml &
# b) statsrpt
filter_stats.py --debug --output="$workPath"/results --report="$workPath"/results/filter_reports_filter_stats.json "$workPath"/data/filtered_summary.csv &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o filter_reports_filter_stats.html -rules "$workPath"/results/.rules_filter_reports_filter_stats.xml &
# c) loadingrpt
pbreport.py loading --debug "$workPath"/results filter_reports_loading.json "$workPath"/data/filtered_summary.csv &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o filter_reports_loading.html -rules "$workPath"/results/.rules_filter_reports_loading.xml &

if [ ! -z $mappingRef ]; then
  # 4. p_mapping
  # 1) align (referencePath in --algorithmoptions may have some problem)
  if [ -z $splitChunk ]; then
    # a) don't split fofn
    #pbalign "$workPath/input.fofn" "$referencePath" "$workPath/data/aligned_reads.cmp.h5" --seed=1 --minaccuracy=0.75 --minlength=50 --concordant --algorithmoptions="-usequality" --algorithmoptions=' -minmatch 12 -bestn 10 -minpctidentity 70.0 -sa /rd1/user/zhongxm/project/genome_assembly/tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/lambda/sequence/lambda.fasta.sa' --hitpolicy=randombest --tmpdir="$workPath"/tmpdir --nproc="$threadNum" --regiontable="$workPath"/data/filtered_regions.fofn
    pbalign "$workPath/input.fofn" "$referencePath" "$workPath/data/aligned_reads.cmp.h5" --seed=1 --minaccuracy=0.75 --minlength=50 --concordant --algorithmoptions="-usequality" --algorithmoptions=' -minmatch 12 -bestn 10 -minpctidentity 70.0 -sa $referencePath/sequence/*.fasta.sa' --hitpolicy=randombest --tmpdir="$workPath"/tmpdir --nproc="$threadNum" --regiontable="$workPath"/data/filtered_regions.fofn
    loadchemistry.py "$workPath"/input.fofn "$workPath"/data/aligned_reads.cmp.h5
    loadpulses "$workPath"/input.fofn "$workPath"/data/aligned_reads.cmp.h5 -metrics deletionqv,ipd,insertionqv,pulsewidth,qualityvalue,mergeqv,substitutionqv,deletiontag -byread
  else
    # b) split fofn (we can split reads_of_insert.fofn to different chunk and align, then merge these files)
    # i) align.plsFofn.Scatter 
    total_lines=`cat "$workPath"/input.fofn | wc -l`
    for((i=0;i<$total_lines;i++))
      do
        awk -v chunkNum=$chunkNum -v i=$i "($total_lines-nr+1)%$chunkNum==$i" "$workPath"/reads_of_insert.fofn > "$workPath"/input.chunk"$i"of"$chunkNum".fofn
      done
    # ii) align_*of*
    ls "$workPath"/input.chunk*of*.fofn | while read file
      do
        i=$( echo $file | sed 's#.*input.##;s#of.*fofn##' )
        pbalign $file "$referencePath" "$workPath/data/aligned_reads.chunk"$i"of"$chunkNum".cmp.h5" --seed=1 --minaccuracy=0.75 --minlength=50 --concordant --algorithmoptions="-usequality" --algorithmoptions=' -minmatch 12 -bestn 10 -minpctidentity 70.0 -sareferencePath/sequence/*.fasta.sa' --hitpolicy=randombest --tmpdir="$workPath"/tmpdir --nproc="$threadNum" --regiontable="$workPath"/data/filtered_regions.chunk"$i"of"$chunkNum".fofn
        loadchemistry.py $file "$workPath/data/aligned_reads.chunk"$i"of"$chunkNum".cmp.h5"
	loadpulses $file "$workPath/data/aligned_reads.chunk"$i"of"$chunkNum".cmp.h5" -metrics deletionqv,ipd,insertionqv,pulsewidth,qualityvalue,mergeqv,substitutionqv,deletiontag -byread
      done
    # iii) align.cmpH5.Gather 
    assertcmph5nonempty.py --debug "$workParth"/data/aligned_reads.chunk*of"$chunkNum".cmp.h5
    cmph5tools.py  merge --outfile="$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/aligned_reads.chunk*of"$chunkNum".cmp.h5
fi
  # 2) sort
  cmph5tools.py -vv sort --deep --inplace "$workPath"/data/aligned_reads.cmp.h5
  # 3) repack (this part could be deleted)
  h5repack -f gzip=1 "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/aligned_reads.cmp.h5_tmp
  mv "$workPath"/data/aligned_reads.cmp.h5_tmp "$workPath"/data/aligned_reads.cmp.h5
  # 4) unmapped && sambam (background)
  # a) unmapped 
  extractunmappedsubreads.py "$workPath"/data/filtered_subreads.fasta "$workPath"/data/aligned_reads.cmp.h5  > "$workPath"/data/unmappedsubreads.fasta &
  # b) sambam 
  (pbsamtools --bam --outfile "$workPath/data/aligned_reads.sam" --refrepos "$referencePath" --readgroup "movie" "$workPath/data/aligned_reads.cmp.h5" && rm "$workPath"/data/aligned_reads.sam) &
  # 5) covgff
  summarize_coverage.py --numregions=500 "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/alignment_summary.gff
  # 6) gff2bed (background)
  gfftobed.py --name=meancoverage --description="mean coverage of genome in fixed interval regions" coverage "$workPath"/data/alignment_summary.gff > "$workPath"/data/coverage.bed &

  # 5. p_mappingreports
  # 1) statsjsonreport && coveragejsonreport (background)
  # a) statsjsonreport
  mapping_stats.py --debug --output="$workPath"/results --mode external --filtersummary="$workPath"/data/filtered_summary.csv "$workPath"/input.fofn "$workPath"/data/filtered_regions.fofn "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/results/mapping_stats_report.json &
  smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o mapping_stats_report.html -rules "$workPath"/results/.rules_mapping_stats_report.xml &
  # b) coveragejsonreport 
  # referencePath may have some problem
  #pbreport.py coverage --debug "$workPath"/results mapping_coverage_report.json '/rd1/user/zhongxm/project/genome_assembly/tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2' "$workPath"/data/alignment_summary.gff
  pbreport.py coverage --debug "$workPath"/results mapping_coverage_report.json '$referencePath' "$workPath"/data/alignment_summary.gff &
  smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o mapping_coverage_report.html -rules "$workPath"/results/.rules_mapping_coverage_report.xml &
  
  if [ ! -z $variantsCall ]; then
    # 6. p_genomicconsensus
    # 1) writecontiglist
    if [ -z $splitChunk ]; then
      echo "ref000001" > "$workPath"/data/contig_ids.txt
    else
      [ -e "$workPath"/data/contig_ids.txt ] && rm "$workPath"/data/contig_ids.txt
      for((i=1;i<=$chunkNum;i++))
        do
	  echo "ref0"$chunkNum >>"$workPath"/data/contig_ids.txt
	done
    fi
    # 2) callvariantswithconsensus
    if [ -z $splitChunk ]; then
      #variantcaller.py --skipunrecognizedcontigs -w "$workPath"/data/contig_ids.txt  -x 5 -q 40 -p "$seymour_home"/analysis/etc/algorithm_parameters/2015-11 -v -j31 --algorithm=quiver "$workPath"/data/aligned_reads.cmp.h5 -r '/rd1/user/zhongxm/project/genome_assembly/tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2/sequence/rheMac2.fasta' -o "$workPath"/data/variants.gff -o "$workPath"/data/consensus.fasta.gz -o "$workPath"/data/consensus.fastq.gz
      variantcaller.py --skipunrecognizedcontigs -w "$workPath"/data/contig_ids.txt  -x 5 -q 40 -p "$seymour_home"/analysis/etc/algorithm_parameters/2015-11 -v -j "$threadNum" --algorithm=quiver "$workPath"/data/aligned_reads.cmp.h5 -r "$referencePath/sequence/rheMac2.fasta" -o "$workPath"/data/variants.gff -o "$workPath"/data/consensus.fasta.gz -o "$workPath"/data/consensus.fastq.gz
    else
      # i) callvariantswithconsensus.contig_list.scatter
      total_lines=`cat "$workPath"/data/contig_ids.txt | wc -l`
      for((i=0;i<$total_lines;i++))
        do
          awk -v chunkNum=$chunkNum -v i=$i "($total_lines-nr+1)%$chunkNum==$i" "$workPath"/contig_ids.txt > "$workPath"/contig_ids.chunk"$i"of"$chunkNum".txt
        done
      # ii) callvariantswithconsensus_*of*
      ls "$workPath"/data/contig_ids.chunk*of*.txt | while read file
	do
	  sample=$( echo $file | sed 's#.*contig_ids.##;s#.txt##' )
	  variantcaller.py --skipunrecognizedcontigs -w $file  -x 5 -q 40 -p "$seymour_home"/analysis/etc/algorithm_parameters/2015-11 -v -j "$threadNum" --algorithm=quiver "$workPath"/data/aligned_reads.cmp.h5 -r "$referencePath/sequence/rheMac2.fasta" -o "$workPath"/data/variants."$sample".gff -o "$workPath"/data/consensus."$sample".fasta.gz -o "$workPath"/data/consensus."$sample".fastq.gz
	done
      # iii) callvariantswithconsensus.consensusfasta.gather && callvariantswithconsensus.consensusfastq.gather (Background)
      (gunzip -c "$workPath"/data/consensus.chunk*of*.fasta.gz | gzip -c > "$workPath"/data/consensus.fasta.gz && gunzip -t "$workPath"/data/consensus.fasta.gz) & 
      (gunzip -c "$workPath"/data/consensus.chunk*of*.fastq.gz | gzip -c > "$workPath"/data/consensus.fastq.gz && gunzip -t "$workPath"/data/consensus.fastq.gz ) &
      # iv) callVariantsWithConsensus.variantsGff.Gather
      cat "$workPath"/data/variants.chunk*of*.gff > "$workPath"/data/variants.gff
      sed -e 1b -e /gff-version/d "$workPath"/data/variants.gff > "$workPath"/tmpdir/variants.gff.tmp
      sed -n '/##/ p' "$workPath"/tmpdir/variants.gff.tmp >> "$workPath"/tmpdir/variants.gff.tmp2
      sed '/##/d' "$workPath"/tmpdir/variants.gff.tmp >> "$workPath"/tmpdir/variants.gff.tmp2
      mv "$workPath"/tmpdir/variants.gff.tmp2 "$workPath"/data/variants.gff
    fi
    # 3) makeVcf && makeBed (Background)
    # a) makeVcf 
    gffToVcf.py --globalReference="$specie" "$workPath"/data/variants.gff > "$workPath"/data/variants.vcf &
    # b) makeBed 
    gffToBed.py --name=variants --description='PacBio: snps, insertions, and deletions derived from consensus calls against reference' variants "$workPath"/data/variants.gff > "$workPath"/data/variants.bed &
    # 4) enrichAlnSummary
    summarizeConsensus.py --variantsGff "$workPath"/data/variants.gff "$workPath"/data/alignment_summary.gff --output "$workPath"/tmpdir/tmpY41til.gff
    mv "$workPath"/tmpdir/tmpY41til.gff "$workPath"/data/alignment_summary.gff
    # 5) zipVariants (Background and this part could be deleted)
    gzip -f "$workPath"/data/variants.gff &

    # 7. P_ConsensusReports
    # 1) topVariantsReport && variantsJsonReport (Background)
    # a) topVariantsReport 
    #pbreport.py topvariants --debug "$workPath"/results top_variants_report.json "$workPath"/data/variants.gff '/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2'
    pbreport.py topvariants --debug "$workPath"/results top_variants_report.json "$workPath"/data/variants.gff "$referencePath" &
    smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o top_variants_report.html -rules "$workPath"/results/.rules_top_variants_report.xml &
    # b) variantsJsonReport
    #pbreport.py variants --debug "$workPath"/results variants_report.json '/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2' "$workPath"/data/alignment_summary.gff "$workPath"/data/variants.gff
    pbreport.py variants --debug "$workPath"/results variants_report.json "$referencePath" "$workPath"/data/alignment_summary.gff "$workPath"/data/variants.gff &
    smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o variants_report.html -rules "$workPath"/results/.rules_variants_report.xml &
  fi
fi
