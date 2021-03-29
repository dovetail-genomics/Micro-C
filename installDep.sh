#!/usr/bin/env bash

################################################################################
## These are the step to install all the dependencies
## used in the examples of this repo
## on Ubuntu 18.04 with the lates version as of 04/11/2020
################################################################################
sudo apt-get -y update && sudo -y apt-get upgrade

## Install GCC make and python and pip
sudo apt-get install -y build-essential python3 python3-pip

## This will put python3 as the default binairies call by python
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

################################################################################
## Install bedtools
################################################################################
sudo apt-get -y install bedtools

################################################################################
## Install deeptools and pairtools
################################################################################
## these are required for pyBigWig
sudo apt-get install -y zlib1g-dev libcurl4

pip3 install \
    pysam \
    tabulate \
    numpy \
    scipy \
    py2bit \
    matplotlib \
    pyBigWig \
    deeptools \
    pandas
    
   pip3 install pairtools

################################################################################
## Install bwa
################################################################################
git clone https://github.com/lh3/bwa.git
cd bwa; make -j $(nproc)
ls
./bwa
sudo cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin
cd

################################################################################
## Install samtools
################################################################################
# First install the following dependencies:
sudo apt-get install -y libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev

## Then download, compile and assemble samtools
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar xvf samtools-1.11.tar.bz2 

cd samtools-1.11/
./configure 
make -j $(nproc)
sudo make install

## install htslib
cd htslib-1.11
make
sudo make install
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
cd

###############
##install liblz4
##############
sudo apt-get install liblz4-tool

################################################################################
## Install preseq
################################################################################
wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz
tar xvf preseq-3.1.2.tar.gz
cd preseq-3*/
./configure --enable-hts CPPFLAGS='-I /home/ubuntu/samtools-1.11/htslib-1.11/' LDFLAGS='-L/home/ubuntu/samtools-1.11/htslib-1.11/'
make -j $(nproc)
sudo make install
cd
