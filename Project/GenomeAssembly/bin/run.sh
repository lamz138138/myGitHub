#!/bin/bash

export workPath="/data7/lcy/zhongxm/PacBio"
export toolsPath="/data7/lcy/zhongxm/tools"

#1. ConsensusTools (get ROI from h5 files) (It is recommend to run in SMRT Portal)
export consensusPath=/data7/lcy/zhongxm/tools/smrtanalysis/current/analysis/bin
cd $workPath/roi/demo
  # 1) run in smp01
    mkdir smp01 && cd smp01
    ln ../data .
    ls $PWD/data/*bax.h5 >input.fofn
    $consensusPath/ConsensusTools.sh CircularConsensus --logFile=consensus.log -n 20 --fofn input.fofn 1> consensus.err 2>&1
    cd ..
  # 2) run with smrtanalysis
    export SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis
    bash
    source $SMRT_ROOT/current/etc/setup.sh
    # a) RS_Subreads
      mkdir subReads && cd subReads
      mkdir log
      mkdir output
      ln -s ../data .
      ls $PWD/data/*bas.h5 >input.fofn
      fofnToSmrtpipeInput.py input.fofn >input.xml 2>log/fofnToSmrtpipeInput.err
      #cp $SMRT_ROOT/current/common/protocols/RS_Subreads.1.xml . 
      smrtpipe.py --distribute --output=output --params=params.xml xml:input.xml >log/RS_Subreads.log 2>log/RS_Subreads.err
      cd ..
    # b) 
    exit
    exit
  # 3) pacbio test data
    cd /data7/lcy/zhongxm/PacBio/data/RS2_test
    /data7/lcy/zhongxm/tools/rar/unrar x part1.rar
    cd -
    mkdir roi/RS2_test
    cd roi/RS2_test
    export SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis
    bash
    source $SMRT_ROOT/current/etc/setup.sh
    # a) RS_Subreads
      mkdir subReads && cd subReads
      mkdir log
      mkdir output
      ln -s ../../../data/RS2_test/B05_1 data
      ls $PWD/data/Analysis_Results/*bas.h5 >input.fofn
      fofnToSmrtpipeInput.py input.fofn >input.xml 2>log/fofnToSmrtpipeInput.err
      #smrtpipe.py --distribute --output=output --params=params_subreads.xml xml:input.xml >log/RS_Subreads.log 2>log/RS_Subreads.err
      cd ..
    # b) RS_CCS
      #mkdir ccsReads && cd ccsReads
      #mkdir log
      #mkdir output
      #ln -s ../../../data/RS2_test/B05_1 data
      #ls $PWD/data/Analysis_Results/*bas.h5 >input.fofn
      #fofnToSmrtpipeInput.py input.fofn >input.xml 2>log/fofnToSmrtpipeInput.err
      #smrtpipe.py --distribute --output=output --params=params_ccs.xml xml:input.xml >log/RS_Subreads.log 2>log/RS_Subreads.err
      #cd ..
    exit
    exit

#2. FALCON-integrate (should run in calculation node)
  export PATH=/data7/lcy/zhongxm/PacBio/bin:$PATH
  # 1) demo
    cd $workPath/falcon/demo
    find $PWD/data -name "*.fasta" > input.fofn
    #cd -
    #cp qsubScript/template.qsub qsubScript/falcon.demo.qsub
    #sed -i "s/jobName/falcon.demo/" qsubScript/falcon.demo.qsub
    #sed -i '/script/d' qsubScript/falcon.demo.qsub
    #echo ". $toolsPath/FALCON-integrate/fc_env/bin/activate" >> qsubScript/falcon.demo.qsub
    #echo "cd $workPath/falcon/demo" >> qsubScript/falcon.demo.qsub 
    #echo "fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err" >> qsubScript/falcon.demo.qsub
    #echo "deactivate" >>qsubScript/falcon.demo.qsub
    #qsub qsubScript/falcon.demo.qsub >log/falcon.demo.log 2>log/falcon.demo.err
    . $toolsPath/FALCON-integrate/fc_env/bin/activate
    fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
    deactivate    

  # 2) demo2 (run in smp01)
    mkdir $workPath/falcon/demo2
    cd $workPath/falcon/demo2
    ln -s ../demo/data .
    find $PWD/data/ -name "*.fasta" | grep -v v1 > input.fofn
    . $toolsPath/FALCON-integrate/fc_env/bin/activate
    fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
    deactivate    
    
  # 3) simulate
    export gSize=48000000
    export targetCov=25
    export minLength=4000
    . $toolsPath/FALCON-integrate/fc_env/bin/activate
    for type in 20XCLR+80XCCS 30XCLR+70XCCS 100XCLR 100XCLR_v2
      do
	mkdir -p $workPath/simulate/$type/falcon
	cd $workPath/simulate/$type/falcon
	mkdir data
	ls /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/$type/*fq | while read file
	  do
	    sample=$( basename $file .fq )
	    fq2fa.pl $file | falcon_name_fastq.pl -i - -t data/"$sample"_old2new.txt -o "$sample".fasta
	    mv "$sample".fasta data
	  done
	find $PWD/data/ -name "*.fasta" > input.fofn
	cp $workPath/falcon/demo/fc_run.cfg .
	fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
	# When there are more then 25X, assembly again with 25X dataset
        #mv 2-asm-falcon 2-asm-falcon_v1
	#minLenMeetTargetCoverage=$(grep -v '^#' 1-preads_ovl/preads4falcon.fasta.faCount | grep -v total | sort -k2,2nr | awk -v gSize=$gSize -v tCov=$targetCov '{t+=$2;if(t>=gSize*tCov)\
	#    {print $2;exit}}')
	#minLenMeetTargetCoverage=${minLenMeetTargetCoverage:-$minLength}
	#sed -i "/length_cutoff_pr =/ s/length_cutoff_pr =.*/length_cutoff_pr = $minLenMeetTargetCoverage/" fc_run.cfg
	#fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
      done
    deactivate    

    # Following is some mark in 20XCLR+80XCCS
    #fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
    #faCount 1-preads_ovl/preads4falcon.fasta >1-preads_ovl/preads4falcon.fasta.faCount
    #mv 2-asm-falcon 2-asm-falcon_v1
    #cp fc_run.cfg fc_run.v1.cfg

  # 4) simulate2
    export gSize=48000000
    export targetCov=25
    export minLength=4000
    . $toolsPath/FALCON-integrate/fc_env/bin/activate
    for type in 30XCLR_RS+70XCLR_Sequel 100XCLR_RS 100XCLR_Sequel 30XCLR_RS+70XCCS_Sequel 30XCLR_RS+70XCCS_RS
      do
	mkdir -p $workPath/simulate2/$type/falcon
	cd $workPath/simulate2/$type/falcon
	find /data7/lcy/zhongxm/PacBio/simulate2/data/$type/ -name "*.fasta" > input.fofn
	cp $workPath/falcon/fc_run.cfg .
	fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
	# When there are more then 25X, assembly again with 25X dataset
        #mv 2-asm-falcon 2-asm-falcon_v1
	#minLenMeetTargetCoverage=$(grep -v '^#' 1-preads_ovl/preads4falcon.fasta.faCount | grep -v total | sort -k2,2nr | awk -v gSize=$gSize -v tCov=$targetCov '{t+=$2;if(t>=gSize*tCov)\
	#    {print $2;exit}}')
	#minLenMeetTargetCoverage=${minLenMeetTargetCoverage:-$minLength}
	#sed -i "/length_cutoff_pr =/ s/length_cutoff_pr =.*/length_cutoff_pr = $minLenMeetTargetCoverage/" fc_run.cfg
	#fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
      done
    deactivate    

  # 5) simulate3 and simulat4
    export gSize=48000000
    export targetCov=25
    export minLength=10000
    . $toolsPath/FALCON-integrate/fc_env/bin/activate
    for dir in simulat3 simulate4
      do
	for type in 30XCLR_RS+30XCLR_Sequel 30XCLR_RS+40XCLR_Sequel 30XCLR_RS+50XCLR_Sequel 30XCLR_RS+60XCLR_Sequel
	  do
	    mkdir -p $workPath/"$dir"/$type/falcon
	    cd $workPath/"$dir"/$type/falcon
	    find /data7/lcy/zhongxm/PacBio/"$dir"/data/$type/ -name "*.fasta" > input.fofn
	    cp $workPath/falcon/fc_run.cfg .
	    fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
	    # When there are more then 25X, assembly again with 25X dataset
	    #mv 2-asm-falcon 2-asm-falcon_v1
	    #minLenMeetTargetCoverage=$(grep -v '^#' 1-preads_ovl/preads4falcon.fasta.faCount | grep -v total | sort -k2,2nr | \
	    #	awk -v gSize=$gSize -v tCov=$targetCov '{t+=$2;if(t>=gSize*tCov) {print $2;exit}}')
	    #minLenMeetTargetCoverage=${minLenMeetTargetCoverage:-$minLength}
	    #sed -i "/length_cutoff_pr =/ s/length_cutoff_pr =.*/length_cutoff_pr = $minLenMeetTargetCoverage/" fc_run.cfg
	    #fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err
	  done
      done
    deactivate    

  # 6) 90G data
    mkdir /data7/lcy/zhongxm/PacBio/falcon/Stage_1 && cd /data7/lcy/zhongxm/PacBio/falcon/Stage_1
    mkdir data && cd data
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016469/data/filtered_subreads.fasta CCS_Subreads_60G.fa
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016470/data/filtered_subreads.fasta CLR_Subreads_17G.fa
    perl -e '
    	my $myHandle="NA";
    	open IN,$ARGV[0];
	while(<IN>){
		if(/^>/){
			chomp;
			my $fileName=(split "/",$_)[0];
			if($fileName ne $myHandle){
				close OUT if $myHandle ne "NA";
				my $myHandle=$fileName;
				open OUT, ">$myHandle".".fa";
			}
			print OUT $_;
			print OUT "\n";
		}else{
			print OUT $_;
		}
	}' CLR_Subreads_17G.fa
    cd ..
    find "$PWD"/data/ -name "*.fasta" > input.fofn
    cp /data7/lcy/zhongxm/PacBio/simulate4/30XCLR_RS+70XCLR_Sequel/falcon/fc_run.cfg .
    fc_run.py fc_run.cfg >fc_run.log 2>fc_run.err


#3. PBcR
  export pbcrPath=/data7/lcy/zhongxm/tools/wgs-8.3rc2/Linux-amd64/bin
  #export pbcrPath=/data7/lcy/zhongxm/tools/smrtanalysis/current/analysis/bin
  export PATH=/data7/lcy/zhongxm/tools/amos/bin:/data8/lcy01/tools/bowtie2-2.2.3:$PATH
  #export LD_LIBRARY_PATH=/data7/lcy/zhongxm/tools/gcc/destDir/lib:/data7/lcy/zhongxm/tools/gcc/destDir/lib64:$LD_LIBRARYPATH
  cd /data7/lcy/zhongxm/tools/smrtanalysis/current/analysis/lib
  mv libstdc++.so.6 libstdc++.so.6_v1
  ln -s /data7/lcy/zhongxm/tools/gcc/destDir/lib64/libstdc++.so.6 .
  cd -
  # 1) demo
    cd $workPath/pbcr/demo
    tar -zxvf sampleData.tar.gz
    cd sampleData
    mkdir log
    java -jar convertFastaAndQualToFastq.jar pacbio.filtered_subreads.fasta >pacbio.filtered_subreads.fastq 2>log/convertFastaAndQualToFastq.err
    # a) in smp node
    $pbcrPath/PBcR -threads 20 -length 500 -partitions 200 -l lambdaIll -s pacbio.spec -fastq pacbio.filtered_subreads.fastq genomeSize=50000 illumina.frg > log/pbcr.log 2> log/pbcr.err
    # b) in SGE
    $pbcrPath/PBcR -threads 20 -length 500 -partitions 200 -l lambdaIll -s pacbio.SGE.spec -fastq pacbio.filtered_subreads.fastq genomeSize=50000 illumina.frg >log/pbcr.log 2>log/pbcr.err
  # 2) simulate
    # a) 20XCLR+80XCCS 30XCLR+70XCCS 30XCLR+100XNGS
    for type in 20XCLR+80XCCS 30XCLR+70XCCS 30XCLR+100XNGS
      do
	mkdir -p $workPath/simulate/$type/pbcr
	cd $workPath/simulate/$type/pbcr
	mkdir data
	cd data
	if [ "$type" == "30XCLR+100XNGS" ]; then
	  ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100X150bps/chr19_art1.fq chr19_illumina_1.fq
	  ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100X150bps/chr19_art2.fq chr19_illumina_2.fq
	  ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/30XCLR+70XCCS/chr19_CLR.fq .
	else
	  ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/$type/chr19_CCS.fq .
	  ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/$type/chr19_CLR.fq .
	fi
	cd ..
	cp $workPath/pbcr/demo/sampleData/pacbio.SGE.spec .
	if [ "$type" == "30XCLR+100XNGS" ]; then
	  $pbcrPath/fastqToCA -libraryname illumina -technology illumina -insertsize 200 30 -mates data/chr19_illumina_1.fq,data/chr19_illumina_2.fq >illuminar.frg 2>fastqToCA.err
	  $pbcrPath/PBcR -threads 20 -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq data/chr19_CLR.fq genomeSize=48000000 illuminar.frg >pbcr.log 2>pbcr.err
	else
	  $pbcrPath/fastqToCA -libraryname CCS -technology pacbio-ccs -reads data/chr19_CCS.fq >CCS.frg 2>fastqToCA.err
	  $pbcrPath/PBcR -threads 20 -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq data/chr19_CLR.fq genomeSize=48000000 CCS.frg >pbcr.log 2>pbcr.err
	fi
      done
    # b) 100XCLR
    export PATH=/data7/lcy/zhongxm/tools/java/bin:/data7/lcy/zhongxm/tools/wgs-8.3rc2/Linux-amd64/bin:$PATH  #java version 1.8 or newer is required for MHAP, enabling MHAP overlapper to align long reads to long reads
    mkdir -p $workPath/simulate/100XCLR/pbcr
    cd $workPath/simulate/100XCLR/pbcr
    mkdir data
    cd data
    ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100XCLR/chr19_CLR.fq .
    cd ..
    cp $workPath/pbcr/demo/sampleData/pacbio.SGE.spec .
    #echo "blasr = -minReadLength 200 -maxScore -1000 -maxLCPLength 16 -bestn 24 -nCandidates 24" >>pacbio.SGE.spec
    echo "ovlConcurrency = 8" >>pacbio.SGE.spec
    # the sensitive parameter were added
    $pbcrPath/PBcR -length 500 -l pbcr -sensitive -s pacbio.SGE.spec -fastq data/chr19_CLR.fq genomeSize=48000000 "blasr=-minReadLength 200 -maxScore -1000 -maxLCPLength 16 -bestn 24 -nCandidates 24" >pbcr.log 2>pbcr.err
    export PATH=$( echo $PATH | sed "s#/data7/lcy/zhongxm/tools/java/bin:/data7/lcy/zhongxm/tools/wgs-8.3rc2/Linux-amd64/bin:##" ):$PATH
    # c) 20XCLR+80XCCS only correction
    mkdir -p $workPath/simulate/20XCLR+80XCCS/pbcr_cor
    cd $workPath/simulate/20XCLR+80XCCS/pbcr_cor
    mkdir data
    cd data
    ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/20XCLR+80XCCS/chr19_CCS.fq .
    ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/20XCLR+80XCCS/chr19_CLR.fq .
    cd ..
    cp $workPath/pbcr/demo/sampleData/pacbio.SGE.spec .
    $pbcrPath/fastqToCA -libraryname CCS -technology pacbio-ccs -reads data/chr19_CCS.fq >CCS.frg 2>fastqToCA.err
    $pbcrPath/PBcR QV=51 assemble=0 -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq data/chr19_CLR.fq genomeSize=48000000 CCS.frg >pbcr.log 2>pbcr.err
  # 3) simulate2
    for type in 30XCLR_RS+70XCLR_Sequel 100XCLR_RS 100XCLR_Sequel 30XCLR_RS+70XCCS_Sequel 30XCLR_RS+70XCCS_RS 30XCLR_RS+100XNGS
      do
	mkdir -p $workPath/simulate2/$type/pbcr
	cd $workPath/simulate2/$type/pbcr
	mkdir data
	cd data
	ls /data7/lcy/zhongxm/PacBio/simulate2/data/"$type"/*fq | while read file
	  do
	    ln -s $file .
	  done
	cd ..
	cp $workPath/pbcr/pacbio.SGE.spec .
	if [ "$type" == "30XCLR_RS+100XNGS" ]; then
	  $pbcrPath/fastqToCA -libraryname illumina -technology illumina -insertsize 200 30 -mates data/chr19_illumina_1.fq,data/chr19_illumina_2.fq >illuminar.frg 2>fastqToCA.err
	  $pbcrPath/PBcR QV=51 -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq data/chr19_CLR_30X.fq genomeSize=48000000 illuminar.frg >pbcr.log 2>pbcr.err
	elif [ "$type" == "30XCLR_RS+70XCCS_Sequel" ] || [ "$type" == "30XCLR_RS+70XCCS_RS" ]
	  $pbcrPath/fastqToCA -libraryname CCS -technology pacbio-ccs -reads data/chr19_CCS_70X.fq >CCS.frg 2>fastqToCA.err
	  $pbcrPath/PBcR QV=51 -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq data/chr19_CLR_30X.fq genomeSize=48000000 CCS.frg >pbcr.log 2>pbcr.err
	else
	  cat data/*fq >data/all.fq
          export PATH=/data7/lcy/zhongxm/tools/java/bin:/data7/lcy/zhongxm/tools/wgs-8.3rc2/Linux-amd64/bin:$PATH  #java version 1.8 or newer is required for MHAP, enabling MHAP overlapper to align long reads to long reads
          echo "ovlConcurrency = 8" >>pacbio.SGE.spec
          # the sensitive parameter were added
          $pbcrPath/PBcR -length 500 -l pbcr -sensitive -s pacbio.SGE.spec -fastq data/chr19_CLR.fq genomeSize=48000000 "blasr=-minReadLength 200 -maxScore -1000 -maxLCPLength 16 -bestn 24 -nCandidates 24" >pbcr.log 2>pbcr.err
          export PATH=$( echo $PATH | sed "s#/data7/lcy/zhongxm/tools/java/bin:/data7/lcy/zhongxm/tools/wgs-8.3rc2/Linux-amd64/bin:##" ):$PATH
	fi
      done
  # 4) RS2_test
    cd /data9/lcy01/zhongxm/PacBio/RS2_test/pbcr
    ln -s /data9/lcy01/zhongxm/NGS/data/processed/fqTrimer/500_1.fq .
    ln -s /data9/lcy01/zhongxm/NGS/data/processed/fqTrimer/500_2.fq .
    ln -s /data7/lcy/zhongxm/PacBio/blasr/B05_1/filtered_subreads.fastq .
    cp $workPath/pbcr/demo/sampleData/pacbio.spec 
    $pbcrPath/fastqToCA -libraryname illumina -technology illumina -insertsize 150 40 -mates 500_1.fq,500_2.fq >illuminar.frg 2>fastqToCA.err
    $pbcrPath/PBcR QV=51 assemble=0 -threads 20 merylMemory=24000 -length 500 -partitions 200 -l pbcr -s pacbio.spec -fastq filtered_subreads.fastq genomeSize=3000000000 illuminar.frg \
    	>pbcr.log 2>pbcr.err
    mkdir v2
    cp /data7/lcy/zhongxm/PacBio/pbcr/pacbio.SGE.spec .
    $pbcrPath/PBcR QV=51 assemble=0 -l pbcr -s pacbio.SGE.spec -fastq filtered_subreads.fastq genomeSize=3000000000 illuminar.frg >pbcr.log 2>pbcr.err
  # 5) Stage_1 CLR
    cd /data7/lcy/zhongxm/PacBio/pbcr/Stage_1
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016468/data/reads_of_insert.fastq CCS.fq
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016470/data/filtered_subreads.fastq CLR_1.fq
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016457/data/filtered_subreads.fastq CLR_2.fq
    cp $workPath/pbcr/pacbio.SGE.spec .
    $pbcrPath/fastqToCA -libraryname CCS -technology pacbio-ccs -reads CCS.fq >CCS.frg 2>fastqToCA.err
    #echo "$PWD/CLR_1.fq" >fastq.file
    #echo "$PWD/CLR_2.fq" >>fastq.file
    cat CLR_1.fq CLR_2.fq >CLR.fq
    $pbcrPath/PBcR -length 500 -partitions 200 -l pbcr -s pacbio.SGE.spec -fastq CLR.fq genomeSize=3000000000 CCS.frg >pbcr.log 2>pbcr.err

  cd /data7/lcy/zhongxm/tools/smrtanalysis/current/analysis/lib
  mv libstdc++.so.6_v1 libstdc++.so.6
  cd -

#4. Canu
  export canuPath=/data7/lcy/zhongxm/tools/canu/Linux-amd64/bin
  export PATH=/data7/lcy/zhongxm/tools/gnuplot/destDir/bin:$PATH
  #export JAVA_HOME=/data7/lcy/zhongxm/tools/java  # canu find java in JAVA_HOME first, but it have some problems, so use java parameter
  #export PERL5LIB=/data7/lcy/zhongxm/tools/perl/lib:$PERL5LIB
  # 1) demo
    cd $workPath/canu/demo
    # a) in smp node
    $canuPath/canu -p ecoli -d ecoli-auto \
     	genomeSize=4.8m \
    	useGrid=0 \
	java=/data7/lcy/zhongxm/tools/java/bin/java \
      	-pacbio-raw p6.25x.fastq \
    	>canu.log 2>canu.err
    # b) in SGE
    $canuPath/canu -p ecoli -d ecoli-auto \
     	genomeSize=4.8m \
	gridEngineThreadsOption="-pe zhongxm THREADS" \
	gridEngineMemoryOption="-l mf=MEMORY" \
	gridOptionsJobName="jobName" \
	java=/data7/lcy/zhongxm/tools/java/bin/java \
	gridOptions="-V -q smp.q -S /bin/bash" \
      	-pacbio-raw p6.25x.fastq \
	>canu.log 2>canu.err
  # 2) simulate
    # a) 20XCLR+80XCCS 30XCLR+70XCCS
    for type in 20XCLR+80XCCS 30XCLR+70XCCS
      do
	mkdir -p $workPath/simulate/$type/canu
	cd $workPath/simulate/$type/canu
	mkdir data
	cd data
	ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/$type/chr19_CLR.fq .
	ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/$type/chr19_CCS.fq .
	cd ..
	$canuPath/canu -p canu -d canu \
	    genomeSize=48m \
	    gridEngineThreadsOption="-pe zhongxm THREADS" \
	    gridEngineMemoryOption="-l mf=MEMORY" \
	    gridOptionsJobName="canuJob" \
	    java=/data7/lcy/zhongxm/tools/java/bin/java \
	    gridOptions="-V -q smp.q -S /bin/bash" \
	    -pacbio-raw data/chr19_CLR.fq \
	    -pacbio-raw data/chr19_CCS.fq \
	    >canu.log 2>canu.err
      done
    # b) 100XCLR
      # i)
      mkdir -p $workPath/simulate/100XCLR/canu
      cd $workPath/simulate/100XCLR/canu
      mkdir data
      cd data
      ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100XCLR/chr19_CLR.fq .
      cd ..
      # corErrorRate was added
      $canuPath/canu -p canu -d canu \
	  genomeSize=48m \
	  gridEngineThreadsOption="-pe zhongxm THREADS" \
	  gridEngineMemoryOption="-l mf=MEMORY" \
	  gridOptionsJobName="canuJob" \
	  java=/data7/lcy/zhongxm/tools/java/bin/java \
	  gridOptions="-V -q smp.q -S /bin/bash" \
	  corErrorRate=0.50 \
	  -pacbio-raw data/chr19_CLR.fq \
	  >canu.log 2>canu.err
      # ii)
      mkdir -p $workPath/simulate/100XCLR_v2/canu
      cd $workPath/simulate/100XCLR_v2/canu
      mkdir data
      cd data
      ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100XCLR_v2/chr19_30X.CLR_0001.fastq chr19_CLR_30X.fq
      ln -s /data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/100XCLR_v2/chr19_70X.CLR_0001.fastq chr19_CLR_70X.fq
      cd ..
      $canuPath/canu -p canu -d canu \
	  genomeSize=48m \
	  gridEngineThreadsOption="-pe zhongxm THREADS" \
	  gridEngineMemoryOption="-l mf=MEMORY" \
	  gridOptionsJobName="canuJob" \
	  java=/data7/lcy/zhongxm/tools/java/bin/java \
	  gridOptions="-V -q smp.q -S /bin/bash" \
	  -pacbio-raw data/chr19_CLR_30X.fq \
	  -pacbio-raw data/chr19_CLR_70X.fq \
	  >canu.log 2>canu.err
  # 3) simulate2
    for type in 30XCLR_RS+70XCLR_Sequel 100XCLR_RS 100XCLR_Sequel 30XCLR_RS+70XCCS_Sequel 30XCLR_RS+70XCCS_RS
      do
	mkdir -p $workPath/simulate2/$type/canu
	cd $workPath/simulate2/$type/canu
	mkdir data
	cd data
	ls /data7/lcy/zhongxm/PacBio/simulate2/data/"$type"/*fq | while read file
	  do
	    ln -s $file .
	  done
	cd ..
	$canuPath/canu -p canu -d canu \
	    genomeSize=48m \
	    gridEngineThreadsOption="-pe zhongxm THREADS" \
	    gridEngineMemoryOption="-l mf=MEMORY" \
	    gridOptionsJobName="canuJob" \
	    java=/data7/lcy/zhongxm/tools/java/bin/java \
	    gridOptions="-V -q smp.q -S /bin/bash" \
	    -pacbio-raw data/*30X.fq \
	    -pacbio-raw data/*70X.fq \
	    >canu.log 2>canu.err
      done
  # 4) Stage_1 CLR
    cd /data7/lcy/zhongxm/PacBio/canu/Stage_1/data
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016469/data/filtered_subreads.fastq CCS.fq
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016470/data/filtered_subreads.fastq CLR_1.fq
    ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016457/data/filtered_subreads.fastq CLR_2.fq
    $canuPath/canu -p canu -d canu \
	genomeSize=3g \
	gridEngineThreadsOption="-pe zhongxm THREADS" \
	gridEngineMemoryOption="-l mf=MEMORY" \
	gridOptionsJobName="canuJob" \
	java=/data7/lcy/zhongxm/tools/java/bin/java \
	gridOptions="-V -q zhongxmQ -S /bin/bash" \
	gnuplotImageFormat=svg \
	-pacbio-raw data/CCS.fq \
	-pacbio-raw data/CLR_1.fq \
	-pacbio-raw data/CLR_2.fq \
	>canu.log 2>canu.err

#5. quast
  export quastPath="/data7/lcy/zhongxm/tools/quast"
  export referencePath="/data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/chr19.fa"
  # 1) simulate
    for type in 20XCLR+80XCCS 30XCLR+70XCCS 30XCLR+100XNGS 100XCLR
      do
	for tool in pbcr falcon canu
	  do
	    mkdir -p $workPath/simulate/$type/mummer/$tool
	    cd $workPath/simulate/$type/mummer/$tool
	    if [ "$tool" == "pbcr" ]; then
              ln -s ../../pbcr/pbcr/9-terminator/asm.scf.fasta scf.fasta
	    elsif [ "$tool" == "falcon" ]; then
              ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
	    else
              ln -s ../../canu/canu/canu.contigs.fasta  scf.fasta
	    fi
	    $quastPath/quast.py scf.fasta \
	      -R $referencePath \
	      -o output \
	      --threads 25 \
	      -s \
	      -e \
	      >quast.log 2>quast.err
	  done
      done

#6. MUMMer
  export mummerPath="/data7/lcy/zhongxm/tools/MUMmer"
  export PATH=/data7/lcy/zhongxm/tools/gnuplot/destDir/bin:$PATH
  export referencePath="/data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/chr19.fa"
  # 1) simulate
    # a) 20XCLR+80XCCS
    # i) pbcr
    for type in 20XCLR+80XCCS 30XCLR+70XCCS 30XCLR+100XNGS 100XCLR
      do
	for tool in pbcr falcon canu
	  do
	    mkdir -p $workPath/simulate/$type/mummer/$tool
	    cd $workPath/simulate/$type/mummer/$tool
	    if [ "$tool" == "pbcr" ]; then
              ln -s ../../pbcr/pbcr/9-terminator/asm.scf.fasta scf.fasta
	    elsif [ "$tool" == "falcon" ]; then
              ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
	    else
              ln -s ../../canu/canu/canu.contigs.fasta  scf.fasta
	    fi
	    $mummerPath/nucmer -p nucmer $referencePath scf.fasta >nucmer.log 2>nucmer.err
	    $mummerPath/mummerplot  -f \
	      -p mummer_output \
	      -s large \
	      -R $referencePath \
	      -Q scf.fasta \
	      -t postscript \
	      nucmer.delta \
	      >mummerplot.log 2>mummerplot.err
	    sed -i "/set output/ s/ps/svg/" mummer_output.gp
	    sed -i "/set terminal/ s/.*/set terminal svg/" mummer_output.gp
	    sed -i "/set mouse/ s/^set/#set/" mummer_output.gp
	    sed -i "/set size/ s/.*/set size 1024,1024/" mummer_output.gp
	    gnuplot mummer_output.gp 
	  done
      done

#7. Statistics after correction
  export PATH=/data7/lcy/zhongxm/PacBio/bin:$PATH
  # 1) simulate
    export gSize=48000000
    for type in 20XCLR+80XCCS 30XCLR+70XCCS 30XCLR+100XNGS 100XCLR
      do
	for tool in pbcr falcon canu
	  do
	    mkdir -p $workPath/simulate/$type/statistics_correct/$tool
	    cd $workPath/simulate/$type/statistics_correct/$tool
	    if [ "$tool" == "pbcr" ]; then
	      ln -s ../../pbcr/pbcr.fasta correct.fasta
	    elsif [ "$tool" == "falcon" ]; then
	      ln -s ../../falcon/1-preads_ovl/preads4falcon.fasta correct.fasta
	    else
	      cp ../../canu/canu/canu.trimmedReads.fasta.gz correct.fasta.gz
	      gunzip correct.fasta.gz
	    fi
	    faCount correct.fasta >correct.faCount
	    totalLength=$( tail -n 1 correct.faCount | cut -f2 )
	    echo "Coverage: " >coverage.tsv
	    echo "scale=6; $totalLength/$gSize" | bc >>coverage.tsv
	    totalReads=$( wc -l correct.faCount | cut -f1 -d " " )
	    totalReads=$[ $totalReads - 2 ]
	    echo "Total reads: $totalReads" >>coverage.tsv
	    readsLengthDistribution.R -file=correct.faCount -p=correctReadsDis.pdf -xlab=length -ylab=count -width=10 -height=10 >readsLengthDis.log 2>&1
	    if [ "$tool" == "canu" ]; then
	      rm -f correct.fasta
	    fi
	  done
      done

#8. Blast (map assembly contig to reference)
  export referencePath="/data7/lcy/zhongxm/PacBio/simulate/chr19ForXiehe/chr19.fa"
  export blastPath="/data7/lcy/zhongxm/tools/ncbi-blast-2.4.0+/bin"
  # 1) simulat2
  outputDir="/data7/lcy/zhongxm/PacBio/simulate2"
  mkdir -p "$outputDir"/blastDb/log
  $blastPath/makeblastdb -in $referencePath -dbtype=nucl -out $outputDir/blastDb/blastDb >"$outputDir"/blastDb/log/blastDb.log 2>"$outputDir"/blastDb/log/blastDb.err
  for type in 30XCLR_RS+70XCLR_Sequel 100XCLR_RS 100XCLR_Sequel 30XCLR_RS+70XCCS_Sequel 30XCLR_RS+70XCCS_RS
    do
      for tool in pbcr falocn pbcr
        do
	  mkdir -p $workPath/simulate2/$type/blast/$tool
	  cd $workPath/simulate2/$type/blast/$tool
	  if [ "$tool" == "pbcr" ]; then
	    ln -s ../../pbcr/pbcr/9-terminator/asm.scf.fasta scf.fasta
	  elsif [ "$tool" == "falcon" ]; then
	    ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
	  else
	    ln -s ../../canu/canu/canu.contigs.fasta  scf.fasta
	  fi
          $blastPath/blastn -query scf.fasta -db $outputDir/blastDb/blastDb -perc_identity 99 -evalue 100 -word_size 50 -out scf.blast.out -outfmt 7 >blast.log 2>blast.err
          #mkdir -p $workPath/simulate2/$type/blast/$tool/blastOut/scaffSplit
          #cd $workPath/simulate2/$type/blast/$tool/blastOut/scaffSplit
          #ln -s ../../scf.fasta scaff.fa
          #perl /data7/lcy/zhongxm/tools/fasta_splitter.pl -n-parts-total 30 scaff.fa >"$outputDir"/log/faSplit.log 2>"$outputDir"/log/faSplit.err
          #cd -
          #export LD_LIBRARY_PATH=/data7/lcy/zhongxm/tools/zlib/destDir/lib:$LD_LIBRARY_PATH
          #ls blastOut/scaffSplit/scaff.part* | while read file
          #  do
          #    sample=$( basename $file .fa | sed 's/.*-//;s/^0//' )
          #    $blastPath/blastn -query $file -db $outputDir/blastDb/blastDb -perc_identity 99 -evalue 100 -word_size 50 -out $outputDir/blastOut/scaff.$sample.blast.out -outfmt 7 \
          #       >"$outputDir"/log/blast."$sample".log 2>"$outputDir"/log/blast."$sample".err
          #  done
          #cat $outputDir/blastOut/scaffSplit/* >$outputDir/blastOut/scaffSplit/assembly.blast.out
        done

#9. Quiver (It is recommend to run in SMRT Portal)
  export SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis
  bash
  source $SMRT_ROOT/current/etc/setup.sh
  # 1) simulat2
  outputDir="/data7/lcy/zhongxm/PacBio/simulate2/"
  mkdir -p "$outputDir"/100XCLR_RS/quiver/falcon
  cd "$outputDir"/100XCLR_RS/quiver/falcon 
  ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
  ln -s /data7/lcy/zhongxm/PacBio/data/RS2_test/B05_1 data
  ls $PWD/data/Analysis_Results/*bas.h5 >input.fofn
  mkdir log
  #sawriter scf.fasta.sa scf.fasta >log/sawriter.log 2>log/sawriter.err
  #blasr data/Analysis_Results/*bax.h5 scf.fasta -sa scf.fasta.sa -m 2 -out blasr.xml -maxAnchorsPerPosition 1000 -nproc 20 >log/blasr.log 2>log/blasr.err
  pbalign --nproc 20 --forQuiver input.fofn scf.fasta pbalign.cmp.h5 >log/pbalign.log 2>log/pbalign.err
  samtools faidx scf.fasta
  quiver pbalign.cmp.h5 --referenceFilename scf.fasta -o consensus.fasta -j 20 >log/quiver.log 2>log/quiver.err  

  exit
  exit

#10. RepeatMasker
  export repeatMaskerPath=/data7/lcy/zhongxm/tools/repeatmasker/RepeatMasker
  # 1) simulat2
  outputDir="/data7/lcy/zhongxm/PacBio/simulate2/"
  mkdir -p "$outputDir"/100XCLR_RS/repeatMasker/falcon
  cd "$outputDir"/100XCLR_RS/repeatMasker/falcon
  mkdir log
  ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
  mkdir output
  $repeatMaskerPath/RepeatMasker -s -pa 20 -species "Rhesus monkey" -dir output scf.fasta >log/repeatMasker.log 2>log/repeatMasker.err

#11. Tandem repeats finder
  export tandemRepeatsFinderPath=/data7/lcy/zhongxm/tools/tandemRepeatsFinder
  # 1) simulat2
  outputDir="/data7/lcy/zhongxm/PacBio/simulate2/"
  mkdir -p "$outputDir"/100XCLR_RS/tandemRepeatsFinder/falcon
  cd "$outputDir"/100XCLR_RS/tandemRepeatsFinder/falcon
  mkdir log
  ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
  $tandemRepeatsFinderPath/trf scf.fasta 2 7 7 80 10 50 500 -d >log/trf.log 2>log/trf.err

#12. Segmental Duplication
  export repeatMaskerPath=/data7/lcy/zhongxm/tools/repeatmasker/RepeatMasker
  export blastPath="/data7/lcy/zhongxm/tools/ncbi-blast-2.4.0+/bin"
  # 1) simulat2
  outputDir="/data7/lcy/zhongxm/PacBio/simulate2/"
  mkdir -p "$outputDir"/100XCLR_RS/segmentalDuplication/falcon
  cd "$outputDir"/100XCLR_RS/segmentalDuplication/falcon
  mkdir log
  ln -s ../../falcon/2-asm-falcon/p_ctg.fa scf.fasta
  # Following mask the sequence
  mkdir output
  $repeatMaskerPath/RepeatMasker -s -pa 20 -xsmall -species "Rhesus monkey" -dir output scf.fasta >log/repeatMasker.log 2>log/repeatMasker.err
  # Following blast
  $blastPath/makeblastdb -in output/scf.fasta.masked -dbtype=nucl -out output/blastDb/blastDb >log/blastDb.log 2>log/blastDb.err
  $blastPath/blastn -query scf.fasta -db output/blastDb/blastDb  -lcase_masking -out scf.blast.out -outfmt "17 std score nident gps"  >log/blast.log 2>log/blast.err

#13. SMRT Analysis
export SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis_2.3
export SMRT_USER=lcy01
export SMRT_GROUP=mobile
$SMRT_ROOT/admin/bin/smrtportald-initd start # This would start mysql and tomcat
#$SMRT_ROOT/admin/bin/mysqld start
#$SMRT_ROOT/admin/bin/tomcatd start
#$SMRT_ROOT/admin/bin/kodosd start

#14. Bwa
ssh c0127
# 1) Suffix Array
mkdir -p index/log && cd index
ln -s ~/ref/rheMac2.fa .
bwa index -a bwtsw rheMac2.fa >log/index.log 2>log/index.err 
cd ..
# 2) Mapping
mkdir -p B05_1/log && cd B05_1
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016443/data/filtered_subreads.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016443/data/filtered_subreads.fastq .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016444/data/reads_of_insert.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016444/data/reads_of_insert.fastq .
mkdir -p G12_1/log && cd G12_1
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016439/data/filtered_subreads.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016439/data/filtered_subreads.fastq .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016445/data/reads_of_insert.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016445/data/reads_of_insert.fastq .
for myData in B05_1 G12_1 H12_1
  do
    cd $myData
    for type in filtered_subreads reads_of_insert
      do
        sample="$myData""_""$type"
	bwa mem -t 4 -x pacbio ../index/rheMac2.fa "$type".fastq >bwa_"$sample".sam 2>log/bwa_"$sample".err
        samtools view -bhS bwa_"$sample".sam -o bwa_"$sample".bam 
        samtools sort bwa_"$sample".bam bwa_"$sample".sorted
      done
  done
# get case for sequence alignment
cd B05_1
sed -n 56457,56509p reads_of_insert.fasta >7773.fa
sed -n 876544,877387p reads_of_insert.fasta >>7773.fa

#15. Blasr
ssh smp03
export SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis
bash
source $SMRT_ROOT/current/etc/setup.sh
# 1) Suffix Array
mkdir -p sa/log && cd sa
ln -s ~/ref/rheMac2.fa scf.fasta
sawriter scf.fasta.sa scf.fasta >log/sawriter.log 2>log/sawriter.err
cd ..
# 2) Mapping
# a) Single Data
mkdir -p B05_1/log && cd B05_1
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016443/data/filtered_subreads.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016443/data/filtered_subreads.fastq .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016444/data/reads_of_insert.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016444/data/reads_of_insert.fastq .
mkdir -p G12_1/log && cd G12_1
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016439/data/filtered_subreads.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016439/data/filtered_subreads.fastq .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016445/data/reads_of_insert.fasta .
ln -s /data7/lcy/zhongxm/tools/smrtanalysis/userdata/jobs/016/016445/data/reads_of_insert.fastq .
for myData in B05_1 G12_1 H12_1
  do
    cd $myData
    for type in filtered_subreads reads_of_insert sub_cor
      do
        sample="$myData""_""$type"
        blasr "$type".fasta ../sa/scf.fasta -sa ../sa/scf.fasta.sa -sam -clipping soft -out "$sample".sam -unaligned "$sample".unaligned -nproc 8 \
    		>log/blasr."$sample".log 2>log/blasr."$sample".err
        samtools view -bhS "$sample".sam -o "$sample".bam
        samtools sort "$sample".bam "$sample".sorted
	samtools view -bh -q 30 "$sample".sorted.bam -o "$sample".q30.bam
	samtools sort "$sample".q30.bam "$sample".q30.sorted
      done
  done
# b) RSII 60G data
mkdir -p RSII/Stage_1/{CCS_60G,CLR_17G}
mkdir -p RSII/Stage_1/CCS_60G/{Subreads,CCS}/log
cd RSII/Stage_1/CCS_60G/CCS
ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016468/data/reads_of_insert.fasta CCS
cd RSII/Stage_1/CCS_60G/Subreads
ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016469/data/filtered_subreads.fasta .
cd RSII/Stage_1/CLR_17G/Subreads
ln -s /data9/lcy01/zhongxm/tools/smrtanalysis/userdata/jobs/016/016470/data/filtered_subreads.fasta .
for type in Subreads CCS
  do
    if [ "$type" == "Subreads" ]; then
      sample="subreads"
    else
      sample="roi"
    fi
    cd $type
    blasr "$type".fasta ../../../../sa/scf.fasta -sa ../../../../sa/scf.fasta.sa -sam -clipping soft -out "$sample".sam -unaligned "$sample".unaligned -nproc 8 \
	    >log/blasr."$sample".log 2>log/blasr."$sample".err
    samtools view -bhS "$sample".sam -o "$sample".bam
    samtools sort "$sample".bam "$sample".sorted
  done

#16. SSPACE
export sspacePath=/data7/lcy/zhongxm/tools/sspace_basic
1) demo
ssh c0220
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/SRR001665/SRR001665_2.fastq.gz
ln -s /data7/lcy/zhongxm/tools/sspace_basic/example/contigs_abyss.fasta .
ln -s /data7/lcy/zhongxm/tools/sspace_basic/example/libraries.txt .
$sspacePath/SSPACE_Basic_v2.0.pl -l libraries.txt -s contigs_abyss.fasta -k 5 -a 0.7 -x 0 -b demo_scaffold -T 30 -v 1 -p > sspace.log 2> sspace.err
