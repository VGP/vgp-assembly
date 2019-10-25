FROM ubuntu:16.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ARG git_commit

# update and install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive \
    apt-get -y install time git make wget autoconf gcc g++ vim tcl8.6-dev autoconf autopoint python-dev \
                       perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev \
                       gfortran build-essential fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude \
                       libreadline-dev openjdk-8-jdk libpcre3 libpcre3-dev libbamtools-dev libsuitesparse-dev \
                       liblpsolve55-dev libboost-iostreams-dev
#&& \
#    apt-get clean && \
#    apt-get purge && \
#    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

### modules
# get modules
WORKDIR /opt
RUN git clone https://github.com/cea-hpc/modules.git && \
    cd /opt/modules && \
    git fetch && \
    git checkout $git_commit && \
    ./configure --with-tclx-ver=8.6 --with-tcl-ver=8.6 && \
    make && \
    make install

# user configuration
RUN echo "source /usr/local/Modules/init/bash" >>/root/.bashrc && \
    echo "module use /root/modules/" >>/root/.bashrc && \
    mkdir /root/modules && \
    mkdir /root/bin

### cmake build tool
WORKDIR /opt/cmake_install
RUN mkdir /opt/cmake && \
    wget https://cmake.org/files/v3.13/cmake-3.13.5-Linux-x86_64.sh && \
    sh cmake-3.13.5-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && \
    ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake && \
    rm -r /opt/cmake_install


### minimap2
# 2.11 
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.11.tar.gz && \
    tar -xvf v2.11.tar.gz && \
    rm -r /opt/minimap2/v2.11.tar.gz && \
    cd minimap2-2.11 && \
    make && \
    mkdir /root/bin/minimap2_2.11 && \
    cp minimap2 /root/bin/minimap2_2.11 && \
    mkdir /root/modules/minimap2 && \
    echo "#%Module" >>/root/modules/minimap2/2.11 && \
    echo "append-path PATH /root/bin/minimap2_2.11" >>/root/modules/minimap2/2.11 && \
    echo "#%Module" >>/root/modules/minimap2/.modulerc && \
    echo "module-version /2.11 default" >>/root/modules/minimap2/.modulerc

### samtools
# 1.9
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar xvf samtools-1.9.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.9.tar.bz2 && \
    cd samtools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses && \
    make && \
    mkdir /root/bin/samtools_1.9 && \
    cp samtools /root/bin/samtools_1.9 && \
    mkdir /root/modules/samtools && \
    echo "#%Module" >>/root/modules/samtools/1.9 && \
    echo "append-path PATH /root/bin/samtools_1.9" >>/root/modules/samtools/1.9 && \
    echo "#%Module" >>/root/modules/samtools/.modulerc && \
    echo "module-version /1.9 default" >>/root/modules/samtools/.modulerc 

### bedtools
## 1.9
WORKDIR /opt/bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && \
    rm -r /opt/bedtools/bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    mkdir /root/bin/bedtools_2.28.0 && \
    mv /opt/bedtools/bedtools2/bin/* /root/bin/bedtools_2.28.0 && \
    mkdir /root/modules/bedtools && \
    echo "#%Module" >>/root/modules/bedtools/2.28.0 && \
    echo "append-path PATH /root/bin/bedtools_2.28.0" >>/root/modules/bedtools/2.28.0 && \
    echo "#%Module" >>/root/modules/bedtools/.modulerc && \
    echo "module-version /2.28.0 default" >>/root/modules/bedtools/.modulerc

### bcftools
## 1.9
WORKDIR /opt/bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar xvf bcftools-1.9.tar.bz2 && \
    rm -r /opt/bcftools/bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses && \
    make && \
    mkdir /root/bin/bcftools_1.9 && \
    cp bcftools /root/bin/bcftools_1.9 && \
    mkdir /root/modules/bcftools && \
    echo "#%Module" >>/root/modules/bcftools/1.9 && \
    echo "append-path PATH /root/bin/bcftools_1.9" >>/root/modules/bcftools/1.9 && \
    echo "#%Module" >>/root/modules/bcftools/.modulerc && \
    echo "module-version /1.9 default" >>/root/modules/bcftools/.modulerc

### bamtools
## 2.5.1
WORKDIR /opt/bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.1.tar.gz && \
    tar -xvf v2.5.1.tar.gz && \
    rm -r /opt/bamtools/v2.5.1.tar.gz && \
    cd bamtools-2.5.1 && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools/bamtools-2.5.1/build ../ && \
    make install && \
    mkdir /root/bin/bamtools_2.5.1 && \
    cp bin/bamtools /root/bin/bamtools_2.5.1 && \
    mkdir /root/modules/bamtools && \
    echo "#%Module" >>/root/modules/bamtools/2.5.1 && \
    echo "append-path PATH /root/bin/bamtools_2.5.1" >>/root/modules/bamtools/2.5.1 && \
    echo "#%Module" >>/root/modules/bamtools/.modulerc && \
    echo "module-version /2.5.1 default" >>/root/modules/bamtools/.modulerc

### blast
# 2.2.30+
WORKDIR /opt/blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    tar -xvf ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    rm -r /opt/blast/ncbi-blast-2.2.30+-x64-linux.tar.gz && \
    mkdir /root/bin/blast_2.2.30+ && \
    mv /opt/blast/ncbi-blast-2.2.30+/bin/* /root/bin/blast_2.2.30+/ && \
    mkdir /root/modules/blast && \
    echo "#%Module" >>/root/modules/blast/2.2.30+ && \
    echo "append-path PATH /root/bin/blast_2.2.30+" >>/root/modules/blast/2.2.30+ && \
    echo "#%Module" >>/root/modules/blast/.modulerc && \
    echo "module-version /2.2.30+ default" >>/root/modules/blast/.modulerc

### hmmer
# 3.2.1
WORKDIR /opt/hmmer
RUN wget https://github.com/EddyRivasLab/hmmer/archive/hmmer-3.2.1.tar.gz && \
    tar -xvf hmmer-3.2.1.tar.gz && \
    rm -r /opt/hmmer/hmmer-3.2.1.tar.gz && \
    cd hmmer-hmmer-3.2.1 && \
    git clone https://github.com/EddyRivasLab/easel.git && \
    cd easel && \
    git checkout 924d7efff9a765d1d5807904f9686e2513a1e4f2 && \
    cd .. && \
    autoconf && \
    ./configure --prefix /opt/hmmer/hmmer-hmmer-3.2.1/install/ && \
    make && \
    make install && \
    mkdir /root/bin/hmmer_3.2.1 && \
    cp install/bin/* /root/bin/hmmer_3.2.1 && \
    mkdir /root/modules/hmmer && \
    echo "#%Module" >>/root/modules/hmmer/3.2.1 && \
    echo "append-path PATH /root/bin/hmmer_3.2.1" >>/root/modules/hmmer/3.2.1 && \
    echo "#%Module" >>/root/modules/hmmer/.modulerc && \
    echo "module-version /3.2.1 default" >>/root/modules/hmmer/.modulerc

### R
# 3.2.1
WORKDIR /opt/R
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.1.tar.gz && \
    tar xvf R-3.2.1.tar.gz && \
    rm -r /opt/R/R-3.2.1.tar.gz && \
    cd R-3.2.1 && \
    ./configure --enable-utf8 && \
    make && \
    mkdir /root/bin/R_3.2.1 && \
    cp -r bin/* /root/bin/R_3.2.1 && \
    mkdir /root/modules/R && \
    echo "#%Module" >>/root/modules/R/3.2.1 && \
    echo "append-path PATH /root/bin/R_3.2.1" >>/root/modules/R/3.2.1 && \
    echo "#%Module" >>/root/modules/R/.modulerc && \
    echo "module-version /3.2.1 default" >>/root/modules/R/.modulerc

### python
# 3.6.0
WORKDIR /opt/python3
RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
    tar xvf Python-3.6.0.tgz && \
    rm -r /opt/python3/Python-3.6.0.tgz && \
    cd Python-3.6.0 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_3.6.0 && \
    cp python /root/bin/python_3.6.0 && \
    mkdir /root/modules/python && \
    echo "#%Module" >>/root/modules/python/3.6 && \
    echo "append-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6 && \
    echo "#%Module" >>/root/modules/python/.modulerc && \
    echo "module-version /3.6 default" >>/root/modules/python/.modulerc

### htslib
# 1.9
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xvf htslib-1.9.tar.bz2 && \
    rm -r /opt/htslib/htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure && \
    make