FROM ubuntu:16.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

ARG modules_git_commit=c999abdc958fe45c735dca093555b5413aecb678

# update and install dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive \
    apt-get -y install time git make wget autoconf gcc g++ vim tcl8.6-dev autoconf autopoint python-dev \
                       perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev \
                       gfortran build-essential fort77 xorg-dev libblas-dev gcc-multilib gobjc++ aptitude \
                       libreadline-dev openjdk-8-jdk libpcre3 libpcre3-dev libbamtools-dev libsuitesparse-dev \
                       liblpsolve55-dev libboost-iostreams-dev libhdf5-dev libboost-all-dev libpng-dev \
                       graphviz sudo
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
    git checkout $modules_git_commit && \
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
# 2.17
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17.tar.bz2 && \
	tar xvf minimap2-2.17.tar.bz2 && \
	cd minimap2-2.17 && \
	make && \
	mkdir /root/bin/minimap2_2.17  && \
    mkdir /root/modules/minimap2 && \
	cp -r minimap2 /root/bin/minimap2_2.17 && \
	echo "#%Module" >>/root/modules/minimap2/2.17 && \
    echo "append-path PATH /root/bin/minimap2_2.17" >>/root/modules/minimap2/2.17 && \
    echo "#%Module" >>/root/modules/minimap2/.modulerc && \
    echo "module-version /2.17 default" >>/root/modules/minimap2/.modulerc
# 2.11 
WORKDIR /opt/minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.11.tar.gz && \
    tar -xvf v2.11.tar.gz && \
    rm -r /opt/minimap2/v2.11.tar.gz && \
    cd minimap2-2.11 && \
    make && \
    mkdir /root/bin/minimap2_2.11 && \
    cp minimap2 /root/bin/minimap2_2.11 && \
    echo "#%Module" >>/root/modules/minimap2/2.11 && \
    echo "append-path PATH /root/bin/minimap2_2.11" >>/root/modules/minimap2/2.11

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

# 1.4.1
WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 && \
    tar xvf samtools-1.4.1.tar.bz2 && \
    rm -r /opt/samtools/samtools-1.4.1.tar.bz2 && \
    cd samtools-1.4.1/ && \
    autoheader && \
    autoconf -Wno-header && \
    ./configure --without-curses && \
    make && \
    mkdir /root/bin/samtools_1.4.1 && \
    cp samtools /root/bin/samtools_1.4.1 && \
    echo "#%Module" >>/root/modules/samtools/1.4.1 && \
    echo "append-path PATH /root/bin/samtools_1.4.1" >>/root/modules/samtools/1.4.1

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
# 2.7
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/2.7/Python-2.7.tgz && \
    tar xvf Python-2.7.tgz && \
    rm /opt/python/Python-2.7.tgz && \
    cd Python-2.7 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_2.7 && \
    cp python /root/bin/python_2.7 && \
    mkdir /root/modules/python && \
    echo "#%Module" >>/root/modules/python/2.7 && \
    echo "prepend-path PATH /root/bin/python_2.7" >>/root/modules/python/2.7 && \
    echo "#%Module" >>/root/modules/python/.modulerc && \
    echo "module-version /2.7 default" >>/root/modules/python/.modulerc
# 3.6.0
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/3.6.0/Python-3.6.0.tgz && \
    tar xvf Python-3.6.0.tgz && \
    rm /opt/python/Python-3.6.0.tgz && \
    cd Python-3.6.0 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_3.6.0 && \
    cp python /root/bin/python_3.6.0 && \
    echo "#%Module" >>/root/modules/python/3.6 && \
    echo "prepend-path PATH /root/bin/python_3.6.0" >>/root/modules/python/3.6
# 3.5.0
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/3.5.0/Python-3.5.0.tgz && \
    tar xvf Python-3.5.0.tgz && \
    rm /opt/python/Python-3.5.0.tgz && \
    cd Python-3.5.0 && \
    ./configure && \
    make && \
    mkdir /root/bin/python_3.5.0 && \
    cp python /root/bin/python_3.5.0 && \
    echo "#%Module" >>/root/modules/python/3.5 && \
    echo "prepend-path PATH /root/bin/python_3.5.0" >>/root/modules/python/3.5

### htslib
# 1.9
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xvf htslib-1.9.tar.bz2 && \
    rm -r /opt/htslib/htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    ./configure && \
    make

### bwa
# 0.7.17
WORKDIR /opt/bwa
RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 && \
	tar xvf bwa-0.7.17.tar.bz2 && \
	rm bwa-0.7.17.tar.bz2 && \
	cd bwa-0.7.17 && \
	make && \
	mkdir /root/bin/bwa_0.7.17  && \
	cp bwa qualfa2fq.pl xa2multi.pl bwakit/* /root/bin/bwa_0.7.17 && \
	mkdir /root/modules/bwa && \
	echo "#%Module" >>/root/modules/bwa/0.7.17 && \
    echo "append-path PATH /root/bin/bwa_0.7.17" >>/root/modules/bwa/0.7.17 && \
    echo "#%Module" >>/root/modules/bwa/.modulerc && \
    echo "module-version /0.7.17 default" >>/root/modules/bwa/.modulerc
 # 0.7.13
 WORKDIR /opt/bwa
 RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.13.tar.bz2 && \
    tar xvf bwa-0.7.13.tar.bz2 && \
    rm bwa-0.7.13.tar.bz2 && \
    cd bwa-0.7.13 && \
    make && \
    mkdir /root/bin/bwa_0.7.13  && \
    cp bwa qualfa2fq.pl xa2multi.pl bwakit/* /root/bin/bwa_0.7.13 && \
    echo "#%Module" >>/root/modules/bwa/0.7.13 && \
     echo "append-path PATH /root/bin/bwa_0.7.13" >>/root/modules/bwa/0.7.13 && \
     echo "#%Module" >>/root/modules/bwa/.modulerc


### mash
# 2.2
WORKDIR /opt/mash
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
	tar xvf mash-Linux64-v2.2.tar && \
	rm mash-Linux64-v2.2.tar && \
	cd mash-Linux64-v2.2 && \
	mkdir /root/bin/mash_v2.2  && \
	cp -r . /root/bin/mash_v2.2 && \
	mkdir /root/modules/mash && \
	echo "#%Module" >>/root/modules/mash/2.2 && \
    echo "append-path PATH /root/bin/mash_v2.2" >>/root/modules/mash/2.2 && \
    echo "#%Module" >>/root/modules/mash/.modulerc && \
    echo "module-version /2.2 default" >>/root/modules/mash/.modulerc

### mummer
# 3.9.4alpha
WORKDIR /opt/mummer
RUN wget https://github.com/mummer4/mummer/releases/download/v3.9.4alpha/mummer-3.9.4alpha.tar.gz && \
	tar xvf mummer-3.9.4alpha.tar.gz && \
	rm mummer-3.9.4alpha.tar.gz && \
	cd mummer-3.9.4alpha && \
	./configure && \
	make && \
	make install && \
	ldconfig && \
	mkdir /root/bin/mummer_3.9.4alpha  && \
	ls -la | grep -v drwx | grep -v config | grep rwx | awk '{print $9}' | xargs -n1 -I{} cp {} /root/bin/mummer_3.9.4alpha/ && \
	mkdir /root/modules/mummer && \
	echo "#%Module" >>/root/modules/mummer/3.9.4 && \
    echo "append-path PATH /root/bin/mummer_3.9.4alpha" >>/root/modules/mummer/3.9.4 && \
    echo "#%Module" >>/root/modules/mummer/.modulerc && \
    echo "module-version /3.9.4 default" >>/root/modules/mummer/.modulerc

### ngmlr
# 0.2.6
WORKDIR /opt/ngmlr
RUN wget https://github.com/philres/ngmlr/releases/download/v0.2.6/ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
	tar xvf ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
	rm ngmlr-0.2.6-beta-linux-x86_64.tar.gz && \
	cd ngmlr-0.2.6 && \
	mkdir /root/bin/ngmlr_0.2.6  && \
	cp -r . /root/bin/ngmlr_0.2.6 && \
	mkdir /root/modules/ngmlr && \
	echo "#%Module" >>/root/modules/ngmlr/0.2.6 && \
    echo "append-path PATH /root/bin/ngmlr_0.2.6" >>/root/modules/ngmlr/0.2.6 && \
    echo "#%Module" >>/root/modules/ngmlr/.modulerc && \
    echo "module-version /0.2.6 default" >>/root/modules/ngmlr/.modulerc


### gnuplot
# 5.2.7
WORKDIR /opt/gnuplot
RUN wget https://sourceforge.net/projects/gnuplot/files/gnuplot/5.2.7/gnuplot-5.2.7.tar.gz && \
	tar xvf gnuplot-5.2.7.tar.gz && \
	rm gnuplot-5.2.7.tar.gz && \
	cd gnuplot-5.2.7 && \
	./configure --prefix=/opt/gnuplot/gnuplot-5.2.7 && \
	make && \
	make install && \
	mkdir /root/bin/gnuplot_5.2.7  && \
	cp -r bin/* /root/bin/gnuplot_5.2.7 && \
	mkdir /root/modules/gnuplot && \
	echo "#%Module" >>/root/modules/gnuplot/5.2.7 && \
    echo "append-path PATH /root/bin/gnuplot_5.2.7" >>/root/modules/gnuplot/5.2.7 && \
    echo "#%Module" >>/root/modules/gnuplot/.modulerc && \
    echo "module-version /5.2.7 default" >>/root/modules/gnuplot/.modulerc

### perl
# 5.18.2
WORKDIR /opt/perl
RUN wget https://www.cpan.org/src/5.0/perl-5.18.2.tar.gz && \
    tar -xvf  perl-5.18.2.tar.gz && \
    rm -r /opt/perl/perl-5.18.2.tar.gz && \
    cd perl-5.18.2 && \
    ./Configure -des -Dprefix=/opt/perl/perl-5.18.2 && \
    make && \
    make install && \
    mkdir /root/bin/perl_5.18.2 && \
    cp perl /root/bin/perl_5.18.2 && \
    mkdir /root/modules/perl && \
    echo "#%Module" >>/root/modules/perl/5.18.2 && \
    echo "prepend-path PATH /root/bin/perl_5.18.2" >>/root/modules/perl/5.18.2

### freebayes
# 1.3.1
WORKDIR /opt/freebayes
RUN git clone --recursive git://github.com/ekg/freebayes.git && \
    cd freebayes && \
    make && \
    mkdir /root/bin/freebayes_1.3.1  && \
    cp bin/freebayes /root/bin/freebayes_1.3.1 && \
    mkdir /root/modules/freebayes && \
     echo "#%Module" >>/root/modules/freebayes/1.3.1 && \
     echo "append-path PATH /root/bin/freebayes_1.3.1" >>/root/modules/freebayes/1.3.1 && \
     echo "#%Module" >>/root/modules/freebayes/.modulerc && \
     echo "module-version /1.3.1 default" >>/root/modules/freebayes/.modulerc
# 1.2.0
WORKDIR /opt/freebayes
RUN cd freebayes && \
    git checkout v1.2.0 && \
    make && \
    mkdir /root/bin/freebayes_1.2.0  && \
    cp bin/freebayes /root/bin/freebayes_1.2.0 && \
    echo "#%Module" >>/root/modules/freebayes/1.2.0 && \
    echo "append-path PATH /root/bin/freebayes_1.2.0" >>/root/modules/freebayes/1.2.0

### picard
# 2.9.2
WORKDIR /opt/picard
RUN git clone https://github.com/broadinstitute/picard.git && \
    cd picard/ && \
    git checkout 2.9.2 && \
    ./gradlew shadowJar && \
    mv build/libs/picard-2.9.2-SNAPSHOT-all.jar build/libs/picard.jar && \
    mkdir /root/bin/picard_2.9.2 && \
    cp build/libs/picard.jar /root/bin/picard_2.9.2 && \
    mkdir /root/modules/picard && \
    echo "#%Module" >>/root/modules/picard/2.9.2 && \
    echo "setenv PICARD_PATH /root/bin/picard_2.9.2" >> /root/modules/picard/2.9.2
# 2.20.6
WORKDIR /opt/picard
RUN cd picard/ && \
    git checkout 2.20.6 && \
    ./gradlew shadowJar && \
    mv build/libs/picard-2.20.6-SNAPSHOT-all.jar build/libs/picard.jar && \
    mkdir /root/bin/picard_2.20.6 && \
    cp build/libs/picard.jar /root/bin/picard_2.20.6 && \
    echo "#%Module" >>/root/modules/picard/2.20.6 && \
    echo "setenv PICARD_PATH /root/bin/picard_2.20.6" >> /root/modules/picard/2.20.6
