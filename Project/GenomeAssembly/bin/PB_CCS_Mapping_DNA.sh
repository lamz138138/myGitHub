#!/bin/bash

usage(){
	cat << EOF
Description:
	This script was used to get CCS reads and map reads to reference

Usage:
	PB_CCS_Mapping_DNA.sh -r reference -w workDir -d rawData

Debug:
	The split module (-m) need to be tested

Options:
	-r      STR     (Optional) The path to the reference [Default: /rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis/references/rheMac2] 
	-w      STR     (Optional) The path of the work directory [Default: .]
	-d	STR	(Optional) The path contain the *bax.h5 files [Defalut: ./rawData]
	-t	INT	(Optional) Number of thread [Default: 10]
	-c	INT	(Optional) Split the work to this number of chunk (-s must be assigned) [default: 10]
	-m		(Optional) Mapping reads to reference
	-s		(Optional) Split the work to small chunk
	-h              Print this information
EOF
	exit 0
}

#[ $1 ] || usage
[ $# -eq 0 ] && usage

workPath="$PWD"
referencePath="/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis/references/rheMac2"
dataPath="$PWD/rawData"
threadNum="10"
chunkNum="10"
while getopts "hc:d:mr:t:w:s" OPTION
  do
    case $OPTION in 
	h) usage;;
	d) dataPath=$OPTARG;;
	m) mappingRef=1;;
	t) threadNum=$OPTARG;;
	r) referencePath=$OPTARG;;
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

# 2. P_CCS Module
# 1) generateCCS
if [ -z "$splitChunk" ]; then
  ConsensusTools.sh CircularConsensus --minFullPasses 1 --minPredictedAccuracy 90 --parameters "SEYMOUR_HOME"/analysis/etc/algorithm_parameters/2015-11 --numThreads "$threadNum" --fofn "$workPath"/input.fofn -o "$workPath"/data
else
  # a) generateCCS.inputPlsFofn.Scatter
  TOTAL_LINES=`cat "$workPath"/input.fofn | wc -l`
  for((i=0;i<$TOTAL_LINES;i++))
    do
      awk -v chunkNum=$chunkNum -v i=$i "($TOTAL_LINES-NR+1)%$chunkNum==$i" "$workPath"/input.fofn > "$workPath"/input.chunk"$i"of"$chunkNum".fofn
    done
  # b) generateCCS_*of*
  ls "$workPath"/input.chunk*of*.fofn | while read file
    do
      ConsensusTools.sh CircularConsensus --minFullPasses 1 --minPredictedAccuracy 90 --parameters "SEYMOUR_HOME"/analysis/etc/algorithm_parameters/2015-11 --numThreads "$threadNum" --fofn $file -o "$workPath"/data
    done
  # c) generateCCS.ccsSentinel.Gather
fi
# 2) gatherFastx
cat "$workPath"/data/*.ccs.fasta >"$workPath"/reads_of_insert.fasta
cat "$workPath"/data/*.ccs.fastq >"$worPath"/reads_of_insert.fastq
# 3) toReadsOfInsertFofn
ls "$workPath"/data/*.ccs.h5 >"$workPath"/reads_of_insert.fofn
# 4) readsOfInsertJsonReport (Background)
reads_of_insert_report.py --debug --output-dir "$workPath"/results "$workPath"/reads_of_insert.fofn "$workPath"/results/reads_of_insert_report.json &
smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o reads_of_insert_report.html -rules "$workPath"/results/.rules_reads_of_insert_report.xml &

if [ ! -z $mappingRef ]; then
  # 3. P_Mapping
  # 1) alignCCS (referencePath in --algorithmOptions may have some problem)
  if [ -z $splitChunk ]; then
    # a) don't split fofn
    #pbalign "$workPath/reads_of_insert.fofn" "$referencePath" "$workPath/data/aligned_reads.cmp.h5" --seed=1 --minAccuracy=0.75 --minLength=50 --useccs=useccsdenovo --algorithmOptions=' -minMatch 12 -bestn 1 -minPctIdentity 70.0 -sa /rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2/sequence/rheMac2.fasta.sa' --hitPolicy=randombest --tmpDir="$workPath"/tmpdir --nproc="$threadNum"
    pbalign "$workPath/reads_of_insert.fofn" "$referencePath" "$workPath/data/aligned_reads.cmp.h5" --seed=1 --minAccuracy=0.75 --minLength=50 --useccs=useccsdenovo --algorithmOptions=" -minMatch 12 -bestn 1 -minPctIdentity 70.0 -sa $referencePath/sequence/*.fasta.sa" --hitPolicy=randombest --tmpDir="$workPath"/tmpdir --nproc="$threadNum"
    loadChemistry.py "$workPath"/reads_of_insert.fofn "$workPath"/data/aligned_reads.cmp.h5
    loadPulses "$workPath"/reads_of_insert.fofn "$workPath"/data/aligned_reads.cmp.h5 -metrics QualityValue -byread
  else
    # b) split fofn (we can split reads_of_insert.fofn to different chunk and align, then merge these files)
    # i) alignCCS.plsFofn.Scatter
    TOTAL_LINES=`cat "$workPath"/reads_of_insert.fofn | wc -l`
    for((i=0;i<$TOTAL_LINES;i++))
      do
        awk -v chunkNum=$chunkNum -v i=$i "($TOTAL_LINES-NR+1)%$chunkNum==$i" "$workPath"/reads_of_insert.fofn > "$workPath"/reads_of_insert.chunk"$i"of"$chunkNum".fofn
      done
    # ii) alignCCS_*of*
    ls "$workPath"/reads_of_insert.chunk*of*.fofn | while read file
      do
        sample=$( echo $file | sed 's#.*reads_of_insert.##;s#.fofn##' )
        pbalign $file "$referencePath" "$workPath/data/aligned_reads."$sample".cmp.h5" --seed=1 --minAccuracy=0.75 --minLength=50 --useccs=useccsdenovo --algorithmOptions=' -minMatch 12 -bestn 1 -minPctIdentity 70.0 -sa $referencePath/sequence/*.fasta.sa' --hitPolicy=randombest --tmpDir=/data7/lcy/zhongxm/tools/smrtanalysis_2.3/tmpdir --nproc="$threadNum"
        loadChemistry.py "$workPath"/reads_of_insert."$sample".fofn "$workPath"/data/aligned_reads."$sample".cmp.h5
	loadPulses "$workPath"/reads_of_insert."$sample".fofn "$workPath"/data/aligned_reads."$sample".cmp.h5 -metrics QualityValue -byread
      done
    # iii) alignCCS.cmpH5.Gather 
    assertCmpH5NonEmpty.py --debug "$workParh"/data/aligned_reads.chunk*of"$chunkNum".cmp.h5
    cmph5tools.py  merge --outFile="$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/aligned_reads.chunk*of"$chunkNum".cmp.h5
  fi
  # 2) sort
  cmph5tools.py -vv sort --deep --inPlace "$workPath"/data/aligned_reads.cmp.h5
  # 3) repack (This part could be deleted)
  h5repack -f GZIP=1 "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/aligned_reads.cmp.h5_TMP
  mv "$workPath"/data/aligned_reads.cmp.h5_TMP "$workPath"/data/aligned_reads.cmp.h5
  # 4) unmapped && samBam (Background)
  # a) unmapped
  extractUnmappedSubreads.py "$workPath"/data/reads_of_insert.fasta "$workPath"/data/aligned_reads.cmp.h5  > "$workPath"/data/unmappedSubreads.fasta &
  # b) samBam
  (pbsamtools --bam --outfile "$workPath/data/aligned_reads.sam" --refrepos "$referencePath" --readGroup "movie" "$workPath/data/aligned_reads.cmp.h5" && rm "$workPath"/data/aligned_reads.sam) &
  # 5) covGFF
  summarize_coverage.py --numRegions=500 "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/data/alignment_summary.gff
  # 6) gff2Bed
  gffToBed.py --name=meanCoverage --description="Mean coverage of genome in fixed interval regions" coverage "$workPath"/data/alignment_summary.gff > "$workPath"/data/coverage.bed

  # 4. P_MappingReports
  # 1) statsJsonReport && coverageJsonReport (Background)
  # a) statsJsonReport 
  mapping_stats.py --debug "$workPath"/results --mode external "$workPath"/input.fofn "$workPath"/input.fofn "$workPath"/data/aligned_reads.cmp.h5 "$workPath"/results/mapping_stats_report.json &
  smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o mapping_stats_report.html -rules "$workPath"/results/.rules_mapping_stats_report.xml &
  # b) coverageJsonReport
  # referencePath may have some problem
  #pbreport.py coverage --debug "$workPath"/results mapping_coverage_report.json '/rd1/user/zhongxm/Project/Genome_Assembly/Tools/smrtanalysis_2.3/install/smrtanalysis_2.3.0.140936/common/references/rheMac2' "$workPath"/data/alignment_summary.gff
  pbreport.py coverage --debug "$workPath"/results mapping_coverage_report.json "$referencePath" "$workPath"/data/alignment_summary.gff &
  smrtreporter -basedir "$workPath"/results -headinclude "$workPath"/results/.martin_header.html --html -o mapping_coverage_report.html -rules "$workPath"/results/.rules_mapping_coverage_report.xml &
fi
