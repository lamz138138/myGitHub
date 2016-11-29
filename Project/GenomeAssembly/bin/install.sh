#!/bin/bash

#1. FALCON-integrate
git clone git://github.com/PacificBiosciences/FALCON-integrate.git
cd FALCON-integrate
git checkout 0.4.0
make init
make virtualenv
make check
make -j install
make test

#2. Python2.7.12
wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
tar -zxvf Python-2.7.12.tgz
cd Python-2.7.12
export LDFLAGS="-L/data7/lcy/zhongxm/tools/zlib/destDir/lib"
./configure --prefix=$PWD
make 
cd -
ln -s Python-2.7.12 Python

#3. Smrtanalysis
wget https://s3.amazonaws.com/files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis_2.3.0.140936.run
wget https://s3.amazonaws.com/files.pacb.com/software/smrtanalysis/2.3.0/smrtanalysis-patch_2.3.0.140936.p4.run
SMRT_ROOT=/data7/lcy/zhongxm/tools/smrtanalysis_2.3
SMRT_USER=lcy01
SMRT_GROUP=mobile
mkdir $SMRT_ROOT
bash smrtanalysis_2.3.0.140936.run -p smrtanalysis-patch_2.3.0.140936.p5.run --rootdir $SMRT_ROOT
ln -s smrtanalysis_2.3 smrtanalysis 

#4. Pbcc
git clone https://github.com/PacificBiosciences/pbccs
cd pbccs
git submodule update --init # delete --remote
mkdir build
cd build
cmake ..
make -j 10

#5. LACHESIS
export PATH=/data7/lcy/zhongxm/tools/boost/destDir:/data7/lcy/zhongxm/tools/gcc/destDir/bin:$PATH
export LD_LIBRARY_PATH=/data7/lcy/zhongxm/tools/boost/destDir/stage/lib:/data7/lcy/zhongxm/tools/zlib/destDir/lib:/data7/lcy/zhongxm/tools/gcc/destDir/lib64:/data7/lcy/zhongxm/tools/gcc/gmp/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpfr/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpc/destDir/lib:$LD_LIBRARY_PATH

  # 1) gcc
    mkdir gcc
    # a) gmp
      mkdir gmp && cd gmp
      wget https://ftp.gnu.org/gnu/gmp/gmp-6.1.1.tar.bz2
      tar -xvjf gmp-6.1.1.tar.bz2
      mv gmp-6.1.1 srcDir
      mkdir destDir
      mkdir objDir && cd objDir
      ../srcDir/configure --prefix=${PWD/objDir/destDir}
      make -j 15
      make install
      cd ../../
    # b) mpfr
      mkdir mpfr && cd mpfr
      wget http://ftp.gnu.org/gnu/mpfr/mpfr-3.1.4.tar.gz
      tar -zxvf mpfr-3.1.4.tar.gz
      mv mpfr-3.1.4 srcDir
      mkdir destDir
      mkdir objDir && cd objDir
      ../srcDir/configure --prefix=${PWD/objDir/destDir} --with-gmp=${PWD/mpfr\/objDir/gmp\/destDir}
      make -j 15
      make install
      cd ../../
    # c) mpc
      mkdir mpc && cd mpc
      wget ftp://ftp.gnu.org/gnu/mpc/mpc-1.0.3.tar.gz
      tar -zxvf mpc-1.0.3.tar.gz
      mv mpc-1.0.3 srcDir
      mkdir destDir
      mkdir objDir && cd objDir
      ../srcDir/configure --prefix=${PWD/objDir/destDir} --with-gmp=${PWD/mpc\/objDir/gmp\/destDir} --with-mpfr=${PWD/mpc\/objDir/mpfr\/destDir}
      make -j 15
      make install
      cd ../../
    # d) gcc 
      export LD_LIBRARY_PATH=$PWD/gmp/destDir/lib:$PWD/mpfr/destDir/lib:$PWD/mpc/destDir/lib:$LD_LIBRARY_PATH
      wget http://mirror.hust.edu.cn/gnu/gcc/gcc-5.4.0/gcc-5.4.0.tar.gz
      tar -zxvf gcc-5.4.0.tar.gz
      mv gcc-5.4.0 srcDir
      mkdir destDir
      mkdir objDir && cd objDir
      ../srcDir/configure --prefix=${PWD/objDir/destDir} --with-gmp=${PWD/objDir/gmp\/destDir} --with-mpfr=${PWD/objDir/mpfr\/destDir} --with-mpc=${PWD/objDir/mpc\/destDir} --disable-multilib
      make -j 15
      make install
      cd ..
    cd ..
    export LD_LIBRARY_PATH=$PWD/gcc/destDir/lib64:$LD_LIBRARY_PATH
    export PATH=$PWD/gcc/destDir/bin:$PATH
  # 2) zlib
    mkdir zlib && cd zlib
    wget http://zlib.net/zlib-1.2.8.tar.gz
    tar -zxvf zlib-1.2.8.tar.gz
    mkdir destDir
    mv zlib-1.2.8 srcDir && cd srcDir
    ./configure --prefix=${PWD/srcDir/destDir}
    make -j 15
    make install
    cd ../../
    export LD_LIBRARY_PATH=$PWD/zlib/destDir/lib:$LD_LIBRARY_PATH # may should revised /etc/ld.so.conf.d/zlib.conf
  # 3) boost
    mkdir boost && cd boost
    wget https://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz
    tar -zxvf boost_1_61_0.tar.gz
    mv boost_1_61_0 destDir
    cd destDir
    ./bootstrap.sh --prefix=$PWD
    ./b2 # may should sudo 
    #./b2 install --prefix=$PWD
    export PATH=$PWD:$PATH
    export LD_LIBRARY_PATH=$PWD/stage/lib:$LD_LIBRARY_PATH 
    cd ../../     
  # 4) samtools
  # 5) BWA
  # 6) BLAST
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.4.0+-x64-linux.tar.gz
    tar -zxvf ncbi-blast-2.4.0+-x64-linux.tar.gz
  # 7) Bedtools
  # 8) Lachesis
    mkdir lachesis && cd lachesis
    wget https://github.com/shendurelab/LACHESIS/archive/master.zip 
    unzip master
    cd LACHESIS-master
    export LACHESIS_BOOST_DIR=${PWD/lachesis\/LACHESIS-master/boost\/destDir}
    export LACHESIS_SAMTOOLS_DIR=/data8/lcy01/tools/samtools-0.1.19
    make -j 15

#4. Picard
wget https://github.com/broadinstitute/picard/releases/download/2.5.0/picard-tools-2.5.0.zip
unzip picard-tools-2.5.0.zip
ln -s picard-tools-2.5.0 picard

#5. fasta-split
wget http://kirill-kryukov.com/study/tools/fasta-splitter/files/fasta-splitter-0.1.1.zip

#6. virtualenv
wget https://pypi.python.org/packages/5c/79/5dae7494b9f5ed061cff9a8ab8d6e1f02db352f3facf907d9eb614fb80e9/virtualenv-15.0.2.tar.gz --no-check-certificate
tar -zxvf virtualenv-15.0.2.tar.gz
ln -s virtualenv-15.0.2 virtualenv
mkdir myVirtualenv
unset PYTHONPATH
/usr/bin/python virtualenv/virtualenv.py -p python2.7 /data7/lcy/zhongxm/tools/myVirtualenv
unset LD_LIBRARY_PATH
unset PYTHONPATH
source /data7/lcy/zhongxm/tools/myVirtualenv/bin/activate
pip install matplotlib
#export PYTHONPATH="/data7/lcy/zhongxm/tools/Python-2.7.12/Lib"
#python virtualenv/virtualenv.py myVirtualenv
#export PYTHONPATH=".:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/pypeflow-0.1.1-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/rdfextras-0.4-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/rdflib-3.4.0-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/pyparsing-1.5.7-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/html5lib-1.0b8-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/isodate-0.5.4-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/six-1.10.0-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages/falcon_kit-0.4.0-py2.7-linux-x86_64.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/FALCON/.eggs/networkx-1.11-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/FALCON/.eggs/decorator-4.0.10-py2.7.egg:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python27.zip:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/plat-linux2:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/lib-tk:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/lib-old:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/lib-dynload:/usr/local/lib/python2.7:/data7/lcy/zhongxm/tools/FALCON-integrate/fc_env/lib/python2.7/site-packages:/data7/lcy/zhongxm/tools/Python-2.7.12/Lib"

#7. manta
#wget https://github.com/Illumina/manta/releases/download/v0.29.6/manta-0.29.6.centos5_x86_64.tar.bz2
#tar -xjf manta-0.29.6.centos5_x86_64.tar.bz2
#ln -s manta-0.29.6.centos5_x86_64 manta

#8. quast
#PATH=/data7/lcy/zhongxm/tools/gcc/destDir/bin:$PATH
#LD_LIBRARY_PATH=/data7/lcy/zhongxm/tools/boost/destDir/stage/lib:/data7/lcy/zhongxm/tools/zlib/destDir/lib:/data7/lcy/zhongxm/tools/gcc/destDir/lib64:/data7/lcy/zhongxm/tools/gcc/gmp/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpfr/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpc/destDir/lib:$LD_LIBRARY_PATH
git clone git://github.com/matplotlib/matplotlib.git
wget https://github.com/Illumina/manta/releases/download/v0.29.6/manta-0.29.6.centos5_x86_64.tar.bz2
wget https://downloads.sourceforge.net/project/quast/quast-4.1.tar.gz
ssh smp01
#unset LD_LIBRARY_PATH
#unset PYTHONPATH
#source /data7/lcy/zhongxm/tools/myVirtualenv/bin/activate
#pip install matplotlib
cp -r myVirtualenv/lib/python2.7/site-packages/matplotlib* Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/numpy* Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/pyparsing* Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/cycler* Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/pytz* Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/six** Python/lib/python2.7/site-packages/
cp -r myVirtualenv/lib/python2.7/site-packages/dateutil* Python/lib/python2.7/site-packages/
tar -xzf quast-4.1.tar.gz
cd quast-4.1
mkdir libs/manta/build
cp ../manta-0.29.6.centos5_x86_64.tar.bz2 libs/manta/build/manta.tar.bz2
sed -i '479,488d' libs/reads_analyzer.py
sed -i '489,490d' libs/reads_analyzer.py
#edit reads_analyzer.py by oundenting line 479,488
python quast.py --test
python quast.py --test-sv
cd -
ln -s quast-4.1 quast

#9. MUMmer
wget http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
tar -zxvf MUMmer3.23.tar.g
cd MUMmer3.23
make check
make install
cd -
ln -s MUMmer3.23 MUMer

#10. gunplot
mkdir gnuplot && cd gnuplot
wget https://sourceforge.net/projects/gnuplot/files/gnuplot/5.0.4/gnuplot-5.0.4.tar.gz --no-check-certificate
tar -zxvf gnuplot-5.0.4.tar.gz
mv gnuplot-5.0.4 srcDir
mkdir objDir destDir
cd objDir
../srcDir/configure --prefix=${PWD/objDir/destDir} --disable-mouse
make
make install

#11. blat
mkdir blat
cd blat
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/FOOTER.txt
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer
chmod u+x blat gfClient gfServer
cd ..

#12. perl
wget http://www.cpan.org/src/5.0/perl-5.24.0.tar.gz
tar -xvzf perl-5.24.0.tar.gz
cd perl-5.24.0
./Configure -des -Dprefix=$PWD >configure.log 2>&1
make >make.log 2>&1
make test
make install >makeInstall.log 2>&1
cd ..
ln -s perl-5.24.0 perl

#13. Qt4
mkdir qt4
cd qt4
wget https://download.qt.io/archive/qt/4.8/4.8.6/qt-everywhere-opensource-src-4.8.6.tar.gz
tar -zxvf qt-everywhere-opensource-src-4.8.6.tar.gz
mv qt-everywhere-opensource-src-4.8.6 srcDir
mkdir objDir destDir
cd objDir
../srcDir/configure -prefix ${PWD/objDir/destDir}
gmake >make.log 2>&1
gmake install >makeInstall.log 2>&1
cd ../../

#14. amos
export PATH=/data7/lcy/zhongxm/tools/MUMmer:/data7/lcy/zhongxm/tools/blat:$PATH
export PERL5LIB=/data7/lcy/zhongxm/tools/perl/lib:$PERL5LIB
tar -zxvf amos-3.1.0.tar.gz
cd amos-3.1.0
sed -i '1i #include <getopt.h>' src/Align/find-tandem.cc
./configure --with-Boost-dir=/data7/lcy/zhongxm/tools/boost/destDir --with-qmake-qt4=/data7/lcy/zhongxm/tools/qt4/destDir/bin/qmake >configure.log 2>&1
make >make.log 2>&1
make install >makeInstall.log 2>&1
cd ..
ln -s amos-3.1.0 amos

#15. PBcR
#export PATH=/data7/lcy/zhongxm/tools/boost/destDir:/data7/lcy/zhongxm/tools/gcc/destDir/bin:$PATH
#export LD_LIBRARY_PATH=/data7/lcy/zhongxm/tools/boost/destDir/stage/lib:/data7/lcy/zhongxm/tools/zlib/destDir/lib:/data7/lcy/zhongxm/tools/gcc/destDir/lib64:/data7/lcy/zhongxm/tools/gcc/gmp/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpfr/destDir/lib:/data7/lcy/zhongxm/tools/gcc/mpc/destDir/lib:$LD_LIBRARY_PATH
export PATH=/data7/lcy/zhongxm/tools/blat:/data7/lcy/zhongxm/tools/MUMmer:/data7/lcy/zhongxm/tools/amos:$PATH
bzip2 -dc wgs-8.3rc2.tar.bz2 | tar -xvf -
cd wgs-8.3rc2
cd kmer
make install >makeInstall.log 2>&1
cd ..
cd src 
make >make.log 2>&1
cd ..
sed -i 's#qsub#qsub -V -S /bin/bash#;s#pe threads#pe zhongxm#;s#mem=2GB#mf=2g#' Linux-amd64/bin/PBcR


#17. canu
#wget https://github.com/marbl/canu/releases/download/v1.3/canu-1.3.Linux-amd64.tar.bz2
#bzip2 -dc canu-1.3*.tar.bz2 | tar -xf -
#export PATH=/data7/lcy/zhongxm/tools/java/bin:/data7/lcy/zhongxm/tools/gnuplot/destDir/bin:$PATH
wget https://github.com/marbl/canu/archive/v1.3.tar.gz
tar -zxvf v1.3.tar.gz
cd canu-1.3/src
make -j 15 >make.log 2>&1
cd ../../
cd /data7/lcy/zhongxm/tools/canu-1.3/Linux-amd64/bin/lib/canu
sed -i "s#set terminal png#set terminal svg#;s#\.png#\.svg#" CorrectReads.pm
sed -i "s#set terminal png#set terminal svg#;s#\.png#\.svg#" Gatekeeper.pm
sed -i "s#set terminal png#set terminal svg#;s#\.png#\.svg#" Meryl.pm
sed -i "s/png/svg/g" HTML.pm
#sed -i "50 s/^use/#use/" Defaults.pm
#sed -i "320,330d" Defaults.pm
#sed -i "320i return(10);" Defaults.pm


#18. RepeatMasker
mkdir repeatmasker && cd repeatmasker
mkdir bin
 # 1) Cross_Match
   # a) Consed
     mkdir consed && cd consed
     tar -zxvf consed_linux.tar.gz
   # b) Phrap
     mkdir phrap && cd phrap
     uudecode temp.mail
     tar -zxvf distrib.tar.gz
     make
 # 2) RMBlast
   mkdir rmblast && cd rmblast
   # a) RMBlast Binaries
     wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/2.2.28/ncbi-rmblastn-2.2.28-x64-linux.tar.gz
   # b) BLAST+ Binaries
     wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/ncbi-blast-2.2.28+-x64-linux.tar.gz
   tar zxvf ncbi-rmblastn-2.2.28-x64-linux.tar.gz
   tar zxvf ncbi-blast-2.2.28+-x64-linux.tar.gz
   cp -R ncbi-rmblastn-2.2.28/* ncbi-blast-2.2.28+/
   rm -rf ncbi-rmblastn-2.2.28
   mv ncbi-blast-2.2.28+/* .
   rm -rf rmblast-2.2.28
 # 3) HMMER
   mkdir hmmer && cd hmmer
   wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
   tar -zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
   mv hmmer-3.1b2-linux-intel-x86_64/* .
   rm -rf hmmer-3.1b2-linux-intel-x86_64
   ./configure --prefix $PWD
   make
   make check
   make install
 # 4) ABBlast/WUBlast (didn't install for money)
 # 5) Tandem Repeat Finder
   mkdir tandemRepeatFinder && cd tandemRepeatFinder
   chmod u+x trf404.linux64
   ln -s trf404.linux64 trf
 # 6) Repeat Database
   # a) Dfam
     mkdir dfam && cd dfam
     wget http://www.dfam.org/web_download/Current_Release/Dfam.hmm.gz
     gunzip Dfam.hmm.gz
     #mv Dfam.hmm RepeatMasker/Libraries
   # b) RepBase
 # 7) RepeatMasker
   export PATH=$PWD/tandemRepeatFinder:$PWD/consed:$PWD/hmmer/bin:$PWD/rmblast/bin:$PATH
   wget http://www.repeatmasker.org/RepeatMasker-open-4-0-6.tar.gz
   tar -zxvf RepeatMasker-open-4-0-6.tar.gz
   cp repeatmaskerlibraries-20160829.tar.gz RepeatMasker
   cd RepeatMasker
   tar -zxvf repeatmaskerlibraries-20160829.tar.gz
   rm repeatmaskerlibraries-20160829.tar.gz
   perl ./configure

#19. tandemRepeatsFinder
mkdir tandemRepeatsFinder && cd tandemRepeatsFinder
chmod u+x trf409.legacylinux64
ln -s trf409.legacylinux64 trf

#20. mega
mkdir mega && cd mega
tar -zxvf megacc-7.0.20-1.x86_64.tar.gz
